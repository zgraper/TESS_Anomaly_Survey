"""Counterexamples from the independent v0.3.0-beta audit."""
import hashlib
import json
import sqlite3
import threading

import numpy as np
import pytest

from tess_anomaly_survey import db, streaming
from tess_anomaly_survey import vetting
from tess_anomaly_survey.experiments import ExperimentStore, SurveyExecutionManager
from tess_anomaly_survey.web import SurveyDatabaseAdapter, resolve_light_curve_series
from test_streaming import create_synthetic_fits


def test_same_tic_two_runs_are_isolated_and_source_immutable(tmp_path):
    path = tmp_path / 'source.sqlite'
    cx = db.connect(path)
    for index, run in enumerate(('a', 'b'), 1):
        cx.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?)', (run, 'now', 'hash', 'fixed_box_2h', .01, '{}', 'cat.csv'))
        cx.execute('INSERT INTO candidates(id,run_id,tic,sector,source_file,center_btjd,score,width_hr,depth_ppm,n_support_windows,review_status,note) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
                   (index, run, 'TIC:1', 1, 'synthetic.npz', 1000., 100., 2., 1000., 1, 'unreviewed', ''))
    cx.commit()
    cx.close()
    before = path.read_bytes()
    with pytest.raises(ValueError, match='explicit run_id'):
        SurveyDatabaseAdapter(path)
    other_source = tmp_path / 'other.sqlite'
    other_source.write_bytes(before)
    with pytest.raises(ValueError, match='Scientific database'):
        SurveyDatabaseAdapter(path, run_id='a', review_db_path=other_source)
    assert other_source.read_bytes() == before
    a, b = [SurveyDatabaseAdapter(path, run_id=r) for r in ('a', 'b')]
    assert a.query_stars()['stars'][0]['total_events'] == 1
    a.record_event_review(1, 'interesting')
    a.record_star_review('TIC:1', 'interesting')
    assert b.get_star_details('TIC:1')['events'][0]['review_status'] == 'unreviewed'
    assert b.get_star_details('TIC:1')['star_status'] == 'unreviewed'
    with pytest.raises(ValueError, match='selected run'):
        a.record_event_review(2, 'artifact')
    with pytest.raises(ValueError, match='another source or run'):
        SurveyDatabaseAdapter(path, run_id='b', review_db_path=a.review_db_path)
    detail = a.get_star_details('TIC:1')
    assert detail['events'][0]['scientific_classification'] == 'not_run'
    assert detail['events'][0]['has_shoulder_flare'] is None
    assert detail['events'][0]['coincidence_risk_level'] == 'unknown'
    store = ExperimentStore(tmp_path / 'exports')
    for native_id, adapter, expected_status in [('a', a, 'interesting'), ('b', b, 'unreviewed')]:
        experiment = store.create_experiment(native_id, '', '', ['1'], [1], is_read_only=True)
        run_id = experiment['runs'][0]['run_id']
        store.update_run(run_id, output_db_path=str(path), source_run_id=native_id,
                         review_db_path=adapter.review_db_path)
        exported = json.loads(store.export_run_catalog(experiment['id'], run_id, 'json')[0])
        assert len(exported['candidates']) == 1
        assert exported['candidates'][0]['run_id'] == native_id
        assert exported['candidates'][0]['human_event_review'] == expected_status
        assert exported['provenance']['source_run_id'] == native_id
    assert path.read_bytes() == before

    # Legacy CLI exports must use the same scoped overlays without opening the
    # scientific database writable or silently merging runs.
    import csv
    from types import SimpleNamespace
    from tess_anomaly_survey.cli import do_export
    output = tmp_path / 'reviewed.csv'
    args = SimpleNamespace(db=str(path), run_id='a', status='interesting',
                           output=str(output), limit=100, max_per_star=0)
    do_export(args)
    with output.open(newline='', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 1 and rows[0]['run_id'] == 'a'
    assert rows[0]['review_status'] == 'interesting'
    args.run_id = None
    with pytest.raises(ValueError, match='explicit run_id'):
        do_export(args)
    assert path.read_bytes() == before


def test_historical_telemetry_is_native_run_scoped(tmp_path, monkeypatch):
    path = tmp_path / 'streaming_pilot.sqlite'
    with sqlite3.connect(path) as cx:
        cx.executescript('''
            CREATE TABLE runs(run_id TEXT);
            CREATE TABLE candidates(run_id TEXT,scientific_classification TEXT);
            CREATE TABLE observations(run_id TEXT,status TEXT);
            INSERT INTO runs VALUES ('a'),('b');
            INSERT INTO candidates VALUES ('a','verified_localized_dimming'),('b','not_run'),('b','not_run');
            INSERT INTO observations VALUES ('a','completed'),('b','completed'),('b','no_anomalies'),('b','failed');
        ''')
    before = path.read_bytes()
    monkeypatch.setenv('TESS_SURVEY_DIR', str(tmp_path))
    store = ExperimentStore(tmp_path / 'workspace')
    runs = {r['source_run_id']: r for r in store.get_experiment('exp_streaming_pilot')['runs']}
    assert (runs['a']['completed_obs'], runs['a']['candidate_count'], runs['a']['verified_count']) == (1, 1, 1)
    assert (runs['b']['completed_obs'], runs['b']['candidate_count'], runs['b']['verified_count']) == (2, 2, 0)
    assert runs['b']['total_targets'] == 3 and runs['b']['failed_obs'] == 1
    assert path.read_bytes() == before


def test_launch_rejects_readonly_and_foreign_runs_before_output(tmp_path):
    store = ExperimentStore(tmp_path / 'workspace')
    a = store.create_experiment(name='a', description='', scientific_objective='', targets=['1'], sectors=[1], is_read_only=True)
    b = store.create_experiment(name='b', description='', scientific_objective='', targets=['2'], sectors=[1])
    manager = SurveyExecutionManager(store)
    with pytest.raises(ValueError, match='read-only'):
        manager.launch_survey(a['id'], a['runs'][0]['run_id'])
    with pytest.raises(ValueError, match='belong'):
        manager.launch_survey(a['id'], b['runs'][0]['run_id'])


def test_resume_rejects_changed_configuration(tmp_path):
    path = str(tmp_path / 'run.sqlite')
    streaming.stream_pipeline([], path, 'r', check_dev_universe=False)
    with pytest.raises(ValueError, match='configuration'):
        streaming.stream_pipeline([], path, 'r', threshold_t=90., check_dev_universe=False)
    with sqlite3.connect(path) as cx:
        value, fingerprint = cx.execute('SELECT threshold_t,config_hash FROM runs').fetchone()
    assert value == streaming.FROZEN_THRESHOLD_T
    assert fingerprint != streaming.FROZEN_CONFIG_HASH


def test_product_cache_and_fits_identity(tmp_path):
    first, second = tmp_path / 'a.fits', tmp_path / 'b.fits'
    create_synthetic_fits(str(first), tic=1)
    create_synthetic_fits(str(second), tic=1, inject_dip=True)
    a = resolve_light_curve_series('1', 1, source_file=str(first), allow_download=False)
    b = resolve_light_curve_series('1', 1, source_file=str(second), allow_download=False)
    assert not np.array_equal(a[1], b[1])
    with pytest.raises(ValueError, match='TIC'):
        resolve_light_curve_series('2', 1, source_file=str(first), allow_download=False)
    with pytest.raises(ValueError, match='sector'):
        streaming.validate_fits_identity(first, '1', 2)
    with pytest.raises(FileNotFoundError, match='substitution'):
        resolve_light_curve_series('1', 1, source_file='missing-synthetic.npz')


def test_cancel_commit_resume_duplicates_and_restart(tmp_path, monkeypatch):
    path = str(tmp_path / 'cancel.sqlite')
    cancel = threading.Event()
    calls = []
    def resolve(tic, sector=None, **kwargs):
        filename = f'tess-s0001-{tic}-lc.fits'
        return dict(sector=1, product_filename=filename, download_url='https://example.invalid/' + filename)
    def worker(item, *args):
        tic = str(item['tic'])
        calls.append(tic)
        if len(calls) == 1:
            cancel.set()
        return dict(success=True, tic=tic, sector=1, obs_id=tic,
                    product_filename=item['product_filename'], data_uri='mast:' + tic,
                    source_url='public-test-fixture', file_sha256='test', prvversion='test',
                    calib_level=2, file_bytes=0, download_sec=0., wall_sec=0.,
                    result=dict(camera=1, ccd=1, tmag=10., teff=5000., status='completed',
                                n_cadences=1000, n_usable_cadences=1000, n_windows=1,
                                n_alert_windows=0, n_events=0, n_verified_events=0,
                                candidates=[], process_sec=0.))
    monkeypatch.setattr(streaming, 'resolve_mast_product', resolve)
    monkeypatch.setattr(streaming, '_worker_process_task', worker)
    targets = [{'tic': '1'}, {'tic': '1'}, {'tic': '1', 'sector': 1}, {'tic': '2'}, {'tic': '3'}]
    result = streaming.stream_pipeline(targets, path, 'r', max_workers=1, check_dev_universe=False, cancel_event=cancel)
    assert result['cancelled'] and calls == ['1']
    with sqlite3.connect(path) as cx:
        assert cx.execute('SELECT COUNT(*) FROM observations').fetchone()[0] == 1
    cancel.clear()
    result = streaming.stream_pipeline(targets, path, 'r', max_workers=1, check_dev_universe=False, cancel_event=cancel)
    assert calls == ['1', '2', '3']
    with sqlite3.connect(path) as cx:
        assert cx.execute('SELECT COUNT(*) FROM observations').fetchone()[0] == 3
    result = streaming.stream_pipeline(targets, path, 'r', max_workers=1, check_dev_universe=False)
    assert result['n_processed'] == 0
    store = ExperimentStore(tmp_path / 'workspace')
    exp = store.create_experiment('restart', '', '', ['1'], [1])
    run_id = exp['runs'][0]['run_id']
    store.update_run(run_id, status='running')
    SurveyExecutionManager(store)
    assert store.get_run(run_id)['status'] == 'interrupted'


def test_versioned_coherence_comparison():
    time = np.arange(720) * 2 / 1440
    good = np.ones(720, dtype=bool)
    isolated = np.zeros(720)
    isolated[[340, 360, 380]] = -.01
    kwargs = dict(center_btjd=time[360], mf_depth_ppm=10000.)
    assert vetting.verify_event_photometry(time, isolated, good, **kwargs)['is_verified']
    result = vetting.verify_event_photometry_v2(time, isolated, good, **kwargs)
    assert not result['is_verified'] and result['verifier_version'] == 'coherence_v2'
    contiguous = np.zeros(720)
    contiguous[340:381] = -.01
    for method in (vetting.verify_event_photometry, vetting.verify_event_photometry_v2):
        assert method(time, contiguous, good, **kwargs)['is_verified']
    # Masking must not join physically disjoint samples into a coherent run.
    good[341:380] = False
    assert not vetting.verify_event_photometry_v2(time, contiguous, good, **kwargs)['is_verified']


def test_partial_failure_recovery_replaces_failed_observation(tmp_path):
    first, second = tmp_path / 'first.fits', tmp_path / 'second.fits'
    create_synthetic_fits(str(first), tic=1)
    create_synthetic_fits(str(second), tic=2)
    targets = [dict(tic=str(tic), sector=1, filename=f'tess-s0001-{tic}_lc.fits', source_url=str(first)) for tic in (1, 2)]
    path = str(tmp_path / 'partial.sqlite')
    kwargs = dict(check_dev_universe=False, threshold_t=1e9, max_workers=1, rate_limit_delay_sec=0)
    result = streaming.stream_pipeline(targets, path, 'recovery', **kwargs)
    assert result['n_successful'] == 1 and result['n_failed'] == 1
    targets[1]['source_url'] = str(second)
    result = streaming.stream_pipeline(targets, path, 'recovery', **kwargs)
    assert result['n_successful'] == 1 and result['n_skipped'] == 1 and result['n_failed'] == 0
    with sqlite3.connect(path) as cx:
        assert cx.execute('SELECT COUNT(*) FROM observations').fetchone()[0] == 2
        assert cx.execute("SELECT COUNT(*) FROM observations WHERE status='failed'").fetchone()[0] == 0
        assert cx.execute('SELECT COUNT(*) FROM failures').fetchone()[0] == 1
