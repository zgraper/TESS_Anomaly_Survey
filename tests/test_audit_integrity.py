"""Counterexamples from the independent v0.3.0-beta audit."""
import hashlib
import sqlite3

import numpy as np
import pytest

from tess_anomaly_survey import db, streaming
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
