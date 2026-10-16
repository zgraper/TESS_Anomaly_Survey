"""Tests for periodicity detection, harmonic analysis, system persistence, and CLI integration."""
import sqlite3
import numpy as np
import pytest
from pathlib import Path

from tess_anomaly_survey import periodicity, db, cli, web


def test_periodicity_synthetic_4_dips():
    """Test standard 4-dip transit system with P=3.5d."""
    t0 = 100.0
    p_true = 3.5
    events = [
        {'center_btjd': t0 + k * p_true, 'score': 15.0 + k * 0.2, 'width_hr': 2.5, 'depth_ppm': 1200.0}
        for k in range(4)
    ]
    res = periodicity.find_periodicity(events, min_events=3, max_period_days=10.0)
    assert res.is_periodic is True
    assert abs(res.period_days - p_true) < 0.01
    assert abs(res.epoch_btjd - t0) < 0.05
    assert res.n_events == 4
    assert res.max_residual_days < 0.05
    assert res.period_snr > 20.0


def test_strict_minimum_three_events():
    """Ensure N < 3 strictly rejects periodic classification."""
    # N = 2 dips
    events_2 = [
        {'center_btjd': 100.0, 'score': 12.0, 'width_hr': 2.0, 'depth_ppm': 1000.0},
        {'center_btjd': 103.5, 'score': 14.0, 'width_hr': 2.0, 'depth_ppm': 1000.0},
    ]
    res_2 = periodicity.find_periodicity(events_2, min_events=3)
    assert res_2.is_periodic is False
    assert res_2.n_events == 2
    assert "Insufficient" in res_2.note

    # N = 1 dip
    events_1 = [{'center_btjd': 100.0, 'score': 10.0, 'width_hr': 2.0, 'depth_ppm': 1000.0}]
    res_1 = periodicity.find_periodicity(events_1, min_events=3)
    assert res_1.is_periodic is False
    assert res_1.n_events == 1
    assert res_1.classification == "single_event"


def test_non_periodic_random_dips():
    """Ensure irregularly spaced dips are rejected from periodic classification."""
    events = [
        {'center_btjd': 100.0, 'score': 12.0, 'width_hr': 2.0, 'depth_ppm': 1000.0},
        {'center_btjd': 101.37, 'score': 10.0, 'width_hr': 2.0, 'depth_ppm': 1000.0},
        {'center_btjd': 106.82, 'score': 14.0, 'width_hr': 2.0, 'depth_ppm': 1000.0},
        {'center_btjd': 115.19, 'score': 11.0, 'width_hr': 2.0, 'depth_ppm': 1000.0},
    ]
    res = periodicity.find_periodicity(events, min_events=3, max_period_days=10.0)
    assert res.is_periodic is False
    assert res.n_events == 4


def test_odd_even_eclipsing_binary_depth_disparity():
    """Test resolution of half-period alias for eclipsing binaries with deep primary and shallow secondary."""
    # Dip recurrence occurs every 1.75 days.
    # Even dips (primary): 5000 ppm. Odd dips (secondary): 1000 ppm.
    t0 = 50.0
    half_p = 1.75
    events = []
    for k in range(6):
        depth = 5000.0 if (k % 2 == 0) else 1000.0
        events.append({
            'center_btjd': t0 + k * half_p,
            'score': 20.0 if (k % 2 == 0) else 8.0,
            'width_hr': 3.0,
            'depth_ppm': depth
        })
    res = periodicity.find_periodicity(events, min_events=3, max_period_days=10.0)
    assert res.is_periodic is True
    # The dip recurrence period is preserved without destructive doubling, and alternate period 2*P is reported
    assert abs(res.period_days - 1.75) < 0.02
    assert abs(res.alternate_period_days - 3.50) < 0.02
    assert res.odd_even_sigma > 3.0


def test_momentum_dump_warning():
    """Test warning note when period is near the ~2.50d TESS momentum dump harmonic."""
    t0 = 20.0
    p = 2.50
    events = [
        {'center_btjd': t0 + k * p, 'score': 10.0, 'width_hr': 2.0, 'depth_ppm': 800.0}
        for k in range(4)
    ]
    res = periodicity.find_periodicity(events, min_events=3, max_period_days=10.0)
    assert res.is_periodic is True
    assert "momentum dump" in res.note.lower()


def test_compute_phase():
    """Test orbital phase folding metric."""
    times = np.array([10.0, 11.0, 12.0, 13.0, 14.0])
    epoch = 10.0
    period = 4.0
    phase = periodicity.compute_phase(times, epoch, period)
    # 10.0 -> 0.0
    # 11.0 -> 0.25
    # 12.0 -> 0.5 (or -0.5)
    # 13.0 -> -0.25
    # 14.0 -> 0.0
    assert abs(phase[0] - 0.0) < 1e-6
    assert abs(phase[1] - 0.25) < 1e-6
    assert abs(phase[3] - (-0.25)) < 1e-6
    assert abs(phase[4] - 0.0) < 1e-6
    assert np.all(phase >= -0.5) and np.all(phase <= 0.5)


def test_db_candidate_systems(tmp_path):
    """Test db schema migration and system storage functions."""
    db_file = tmp_path / "test_periodicity.sqlite"
    cx = db.connect(db_file)

    # Insert test run and candidate events
    with cx:
        cx.execute("INSERT INTO runs VALUES ('run1', '2026-01-01', 'hash1', 'fixed_box_2h', 0.01, '{}', 'cat.csv')")
        cur1 = cx.execute('''INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd, score, width_hr, depth_ppm)
            VALUES ('run1', 'TIC_100', 1, 's1.npz', 10.0, 15.0, 2.0, 1000.0)''')
        c1 = cur1.lastrowid
        cur2 = cx.execute('''INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd, score, width_hr, depth_ppm)
            VALUES ('run1', 'TIC_100', 1, 's1.npz', 13.5, 14.0, 2.0, 1000.0)''')
        c2 = cur2.lastrowid
        cur3 = cx.execute('''INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd, score, width_hr, depth_ppm)
            VALUES ('run1', 'TIC_100', 1, 's1.npz', 17.0, 16.0, 2.0, 1000.0)''')
        c3 = cur3.lastrowid

    # Create PeriodCandidate result
    result = periodicity.PeriodCandidate(
        period_days=3.5, epoch_btjd=10.0, n_events=3, is_periodic=True,
        odd_even_sigma=0.5, max_residual_days=0.01, period_snr=26.0, note="Synthetic test system"
    )

    sys_id = db.save_candidate_system(cx, 'run1', 'TIC_100', result, [c1, c2, c3])
    assert sys_id > 0

    # Retrieve system
    sys_info = db.get_system(cx, sys_id)
    assert sys_info is not None
    assert sys_info['tic'] == 'TIC_100'
    assert abs(sys_info['period_days'] - 3.5) < 1e-4
    assert sys_info['is_periodic'] is True
    assert sys_info['n_events'] == 3

    # Check candidates have system_id linked
    rows = cx.execute('SELECT system_id FROM candidates WHERE id IN (?,?,?)', (c1, c2, c3)).fetchall()
    assert all(r[0] == sys_id for r in rows)

    # Check star summary includes system
    summary = db.get_star_summary(cx, 'TIC_100')
    assert summary['system'] is not None
    assert summary['system']['id'] == sys_id
    assert summary['system']['is_periodic'] is True

    # Check get_systems_for_run
    systems = db.get_systems_for_run(cx, 'run1')
    assert len(systems) == 1
    assert systems[0]['id'] == sys_id

    # Check stats count
    stats = db.get_survey_review_stats(cx)
    assert stats['stars']['periodic_systems'] == 1

    cx.close()


def test_cli_analyze_periods(tmp_path):
    """Test CLI analyze-periods command populating candidate_systems."""
    db_file = tmp_path / "cli_test.sqlite"
    cx = db.connect(db_file)
    with cx:
        cx.execute("INSERT INTO runs VALUES ('runA', '2026-01-01', 'h1', 'fixed_box_2h', 0.01, '{}', 'c.csv')")
        # Star 1: 4 dips, P=4.0d
        for k in range(4):
            cx.execute('''INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd, score, width_hr, depth_ppm)
                VALUES ('runA', 'TIC_P', 1, 'f1.npz', ?, 12.0, 2.0, 1000.0)''', (100.0 + k * 4.0,))
        # Star 2: 2 dips (should be skipped due to min_events=3)
        for k in range(2):
            cx.execute('''INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd, score, width_hr, depth_ppm)
                VALUES ('runA', 'TIC_SHORT', 1, 'f2.npz', ?, 12.0, 2.0, 1000.0)''', (200.0 + k * 3.0,))
    cx.close()

    # Run CLI command
    cli.main(['analyze-periods', '--db', str(db_file), '--min-events', '3', '--max-period', '10.0'])

    cx = db.connect(db_file)
    p_sys = cx.execute("SELECT period_days, is_periodic FROM candidate_systems WHERE tic='TIC_P'").fetchone()
    assert p_sys is not None
    assert abs(p_sys[0] - 4.0) < 0.05
    assert p_sys[1] == 1

    # Star 2 should not have a candidate system record
    short_sys = cx.execute("SELECT 1 FROM candidate_systems WHERE tic='TIC_SHORT'").fetchone()
    assert short_sys is None
    cx.close()


def test_web_render_phase_plot():
    """Test phase plot PNG rendering function."""
    # When is_periodic is False
    png_bytes = web.render_phase_plot_png(
        source_file="nonexistent.npz", center_btjd=10.0, width_hr=2.0, score=10.0, depth_ppm=1000.0,
        tic="TIC_123", epoch_btjd=10.0, period_days=0.0, is_periodic=False, odd_even_sigma=0.0, root_path=None
    )
    assert png_bytes.startswith(b'\x89PNG\r\n\x1a\n')

    # When no root_path specified
    png_bytes_periodic = web.render_phase_plot_png(
        source_file="nonexistent.npz", center_btjd=10.0, width_hr=2.0, score=10.0, depth_ppm=1000.0,
        tic="TIC_123", epoch_btjd=10.0, period_days=3.5, is_periodic=True, odd_even_sigma=0.0, root_path=None
    )
    assert png_bytes_periodic.startswith(b'\x89PNG\r\n\x1a\n')


def test_web_phase_endpoint_live(tmp_path):
    """Test /api/plot/phase/<cid> and /api/candidates periodicity endpoints over HTTP."""
    import json
    import threading
    import urllib.request
    from http.server import HTTPServer

    db_file = tmp_path / "web_periodicity.sqlite"
    cx = db.connect(db_file)
    with cx:
        cx.execute("INSERT INTO runs VALUES ('r1', '2026-01-01', 'h', 'fixed_box_2h', 0.01, '{}', 'c.csv')")
        cur = cx.execute('''INSERT INTO candidates(id, run_id, tic, sector, source_file, center_btjd,
                            score, width_hr, depth_ppm)
                            VALUES (1, 'r1', 'TIC_WEB', 1, 'f.npz', 10.0, 15.0, 2.0, 1000.0)''')
    res = periodicity.PeriodCandidate(
        period_days=3.5, epoch_btjd=10.0, n_events=3, is_periodic=True,
        odd_even_sigma=1.2, max_residual_days=0.01, period_snr=25.0,
        period_uncertainty_days=0.001, epoch_uncertainty_days=0.002,
        alternate_period_days=7.0, classification="periodic", note="Live test"
    )
    db.save_candidate_system(cx, 'r1', 'TIC_WEB', res, [1])
    cx.close()

    class CustomHandler(web.AnomalyReviewHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(db_file, tmp_path, *args, **kwargs)

    server = HTTPServer(("127.0.0.1", 0), CustomHandler)
    host, port = server.server_address
    t = threading.Thread(target=server.serve_forever, daemon=True)
    t.start()
    base_url = f"http://{host}:{port}"

    try:
        # Test candidate list includes period_days, is_periodic, uncertainty, and alt period
        with urllib.request.urlopen(f"{base_url}/api/candidates") as resp:
            assert resp.status == 200
            data = json.loads(resp.read().decode('utf-8'))
            assert len(data['candidates']) == 1
            cand = data['candidates'][0]
            assert cand['is_periodic'] is True
            assert abs(cand['period_days'] - 3.5) < 1e-4
            assert cand['period_uncertainty_days'] == 0.001
            assert cand['alternate_period_days'] == 7.0
            assert cand['system_classification'] == "periodic"

        # Test candidate detail includes system info
        with urllib.request.urlopen(f"{base_url}/api/candidate/1") as resp:
            assert resp.status == 200
            data = json.loads(resp.read().decode('utf-8'))
            sys_info = data['star']['system']
            assert sys_info is not None
            assert sys_info['is_periodic'] is True
            assert abs(sys_info['period_days'] - 3.5) < 1e-4
            assert sys_info['period_uncertainty_days'] == 0.001
            assert sys_info['alternate_period_days'] == 7.0

        # Test phase plot endpoint
        with urllib.request.urlopen(f"{base_url}/api/plot/phase/1") as resp:
            assert resp.status == 200
            assert resp.headers.get_content_type() == "image/png"
            png_bytes = resp.read()
            assert png_bytes.startswith(b'\x89PNG\r\n\x1a\n')

        # Test phase plot with ?alt=1 query param
        with urllib.request.urlopen(f"{base_url}/api/plot/phase/1?alt=1") as resp:
            assert resp.status == 200
            assert resp.headers.get_content_type() == "image/png"
            png_bytes_alt = resp.read()
            assert png_bytes_alt.startswith(b'\x89PNG\r\n\x1a\n')
    finally:
        server.shutdown()
        server.server_close()


def test_missing_transits_data_gap():
    """Test period recovery when observations have missing transits due to downlink gaps."""
    # True period is 2.80 days, but cycles 2, 5, 6 fall into data gaps
    t0 = 50.0
    p_true = 2.80
    observed_cycles = [0, 1, 3, 4, 7]
    events = [
        {'center_btjd': t0 + k * p_true, 'score': 14.0, 'width_hr': 2.5, 'depth_ppm': 1500.0}
        for k in observed_cycles
    ]
    res = periodicity.find_periodicity(events, min_events=3, max_period_days=10.0)
    assert res.is_periodic is True
    assert abs(res.period_days - p_true) < 0.02
    assert abs(res.epoch_btjd - t0) < 0.05
    assert res.n_events == 5


def test_variable_eclipse_depths_with_starspots():
    """Test candidate with depth variations (e.g. active starspots crossing transit)."""
    t0 = 80.0
    p_true = 3.20
    # Modulated depths between 8,000 and 12,000 ppm
    depths = [10000.0, 11500.0, 8900.0, 10200.0, 12000.0]
    events = [
        {'center_btjd': t0 + k * p_true, 'score': 16.0, 'width_hr': 2.0, 'depth_ppm': depths[k]}
        for k in range(5)
    ]
    res = periodicity.find_periodicity(events, min_events=3, max_period_days=10.0)
    assert res.is_periodic is True
    assert abs(res.period_days - p_true) < 0.02
    assert res.n_events == 5


def test_irregular_boyajian_dust_dips():
    """Test chaotic, irregularly spaced dips (e.g. circumstellar dust or flare overshoot)."""
    # Randomly spaced dips without an underlying periodic ephemeris
    events = [
        {'center_btjd': 100.0, 'score': 25.0, 'width_hr': 4.0, 'depth_ppm': 5000.0},
        {'center_btjd': 102.14, 'score': 18.0, 'width_hr': 3.0, 'depth_ppm': 3000.0},
        {'center_btjd': 107.81, 'score': 30.0, 'width_hr': 6.0, 'depth_ppm': 8000.0},
        {'center_btjd': 115.33, 'score': 22.0, 'width_hr': 3.5, 'depth_ppm': 4500.0},
        {'center_btjd': 117.90, 'score': 15.0, 'width_hr': 2.0, 'depth_ppm': 2000.0},
    ]
    res = periodicity.find_periodicity(events, min_events=3, max_period_days=10.0)
    assert res.is_periodic is False
    assert res.period_days == 0.0


def test_momentum_dump_subharmonic_warning():
    """Test warning note when period is near the ~1.25d TESS momentum dump subharmonic."""
    t0 = 10.0
    p = 1.25
    events = [
        {'center_btjd': t0 + k * p, 'score': 12.0, 'width_hr': 2.0, 'depth_ppm': 700.0}
        for k in range(4)
    ]
    res = periodicity.find_periodicity(events, min_events=3, max_period_days=10.0)
    assert res.is_periodic is True
    assert "momentum dump" in res.note.lower()


def test_cluster_candidate_events():
    """Verify that intra-dip candidate rows are properly clustered to the highest SNR event."""
    raw_events = [
        # Dip 1 (at ~10.0 BTJD, 3 overlapping rows)
        {'center_btjd': 10.00, 'score': 20.0, 'width_hr': 2.0, 'depth_ppm': 5000.0},
        {'center_btjd': 10.02, 'score': 28.0, 'width_hr': 2.0, 'depth_ppm': 5200.0},
        {'center_btjd': 9.98, 'score': 15.0, 'width_hr': 2.0, 'depth_ppm': 4800.0},
        # Dip 2 (at ~12.5 BTJD, 2 overlapping rows)
        {'center_btjd': 12.50, 'score': 32.0, 'width_hr': 2.0, 'depth_ppm': 5100.0},
        {'center_btjd': 12.51, 'score': 18.0, 'width_hr': 2.0, 'depth_ppm': 4900.0},
        # Dip 3 (at ~15.0 BTJD, 2 overlapping rows)
        {'center_btjd': 15.00, 'score': 30.0, 'width_hr': 2.0, 'depth_ppm': 5050.0},
        {'center_btjd': 14.99, 'score': 12.0, 'width_hr': 2.0, 'depth_ppm': 4700.0},
    ]
    clustered = periodicity.cluster_candidate_events(raw_events, merge_radius_days=0.1)
    assert len(clustered) == 3
    # Check that highest score was selected for each cluster
    assert clustered[0]['score'] == 28.0
    assert abs(clustered[0]['center_btjd'] - 10.02) < 1e-4
    assert clustered[1]['score'] == 32.0
    assert abs(clustered[1]['center_btjd'] - 12.50) < 1e-4
    assert clustered[2]['score'] == 30.0
    assert abs(clustered[2]['center_btjd'] - 15.00) < 1e-4

    # Now verify find_periodicity handles the raw unclustered events seamlessly
    res = periodicity.find_periodicity(raw_events, min_events=3, max_period_days=10.0)
    assert res.is_periodic is True
    assert abs(res.period_days - 2.50) < 0.05
    assert res.n_events == 3


def test_formal_uncertainty_and_alternate_period():
    """Verify least-squares ephemeris uncertainties and EB alternate period preservation."""
    t0 = 20.0
    p = 2.0
    # Simulate EB with primary (depth 10,000 ppm) and secondary (depth 2,000 ppm) eclipses
    events = []
    for k in range(8):
        depth = 10000.0 if k % 2 == 0 else 2000.0
        events.append({
            'center_btjd': t0 + k * p,
            'score': 25.0 if k % 2 == 0 else 12.0,
            'width_hr': 2.0,
            'depth_ppm': depth
        })
    res = periodicity.find_periodicity(events, min_events=3, max_period_days=10.0)
    assert res.is_periodic is True
    assert res.period_uncertainty_days is not None
    assert res.period_uncertainty_days > 0.0
    assert res.epoch_uncertainty_days is not None
    assert res.epoch_uncertainty_days > 0.0
    assert res.alternate_period_days is not None
    # For alternating depths, alternate period should be 2*P (4.0 d)
    assert abs(res.alternate_period_days - 4.0) < 0.05
    assert len(res.competing_solutions) >= 1
    assert res.classification == "periodic"


def test_database_periodicity_migration_and_persistence(tmp_path):
    """Test safe schema migration, column additions, and persistence for candidate_systems."""
    import sqlite3
    db_file = tmp_path / "test_migration.sqlite"
    # Create legacy-style database without Phase 3.1 columns
    legacy_schema = '''
    CREATE TABLE runs (id TEXT PRIMARY KEY, created_utc TEXT, calibration_sha256 TEXT,
                       method TEXT, target_alert_fraction REAL, settings_json TEXT, catalog_path TEXT);
    CREATE TABLE candidates (id INTEGER PRIMARY KEY, run_id TEXT, tic TEXT, sector INTEGER,
                             source_file TEXT, center_btjd REAL, score REAL, width_hr REAL, depth_ppm REAL,
                             n_support_windows INTEGER DEFAULT 1, review_status TEXT DEFAULT 'unreviewed', note TEXT DEFAULT '');
    CREATE TABLE candidate_systems (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        run_id TEXT NOT NULL,
        tic TEXT NOT NULL,
        n_events INTEGER NOT NULL,
        period_days REAL,
        epoch_btjd REAL,
        period_snr REAL,
        is_periodic INTEGER NOT NULL DEFAULT 0,
        odd_even_sigma REAL,
        classification TEXT NOT NULL DEFAULT 'unreviewed',
        note TEXT NOT NULL DEFAULT '',
        created_utc TEXT NOT NULL
    );
    '''
    raw_con = sqlite3.connect(str(db_file))
    raw_con.executescript(legacy_schema)
    raw_con.execute("INSERT INTO runs VALUES ('r1', '2026-01-01', 'h', 'fixed_box_2h', 0.01, '{}', 'c.csv')")
    raw_con.execute("INSERT INTO candidates(id, run_id, tic, sector, source_file, center_btjd, score, width_hr, depth_ppm) VALUES (1, 'r1', 'TIC_MIG', 1, 's.npz', 10.0, 15.0, 2.0, 1000.0)")
    raw_con.close()

    # Now connect via db.connect, which should detect and safely add the missing columns
    cx = db.connect(db_file)
    cols = [r[1] for r in cx.execute("PRAGMA table_info(candidate_systems)").fetchall()]
    assert "period_uncertainty_days" in cols
    assert "epoch_uncertainty_days" in cols
    assert "alternate_period_days" in cols
    assert "competing_solutions_json" in cols

    # Save a PeriodCandidate with uncertainties and alternate period
    cand_sys = periodicity.PeriodCandidate(
        period_days=2.3682, epoch_btjd=1411.6123, n_events=12, is_periodic=True,
        odd_even_sigma=298.8, max_residual_days=0.0069, period_snr=6460.4,
        period_uncertainty_days=0.0002, epoch_uncertainty_days=0.0001,
        alternate_period_days=4.7364, classification="periodic",
        competing_solutions=[{'period_days': 4.7364, 'note': 'EB 2*P'}],
        note="EB primary/secondary"
    )
    sys_id = db.save_candidate_system(cx, 'r1', 'TIC_MIG', cand_sys, [1])
    assert sys_id > 0

    # Retrieve via get_system
    sys_record = db.get_system(cx, sys_id)
    assert sys_record is not None
    assert abs(sys_record['period_days'] - 2.3682) < 1e-4
    assert abs(sys_record['period_uncertainty_days'] - 0.0002) < 1e-4
    assert abs(sys_record['alternate_period_days'] - 4.7364) < 1e-4
    assert sys_record['classification'] == "periodic"
    assert len(sys_record['competing_solutions']) == 1

    # Retrieve via get_star_summary
    star_sum = db.get_star_summary(cx, 'TIC_MIG')
    assert star_sum['system'] is not None
    assert abs(star_sum['system']['period_days'] - 2.3682) < 1e-4
    assert abs(star_sum['system']['period_uncertainty_days'] - 0.0002) < 1e-4
    assert abs(star_sum['system']['alternate_period_days'] - 4.7364) < 1e-4
    cx.close()


def test_adaptive_clustering_broad_transits():
    """Verify adaptive clustering safely merges intra-dip detections spanning > 1 hour on broad transits."""
    # 4-hour transit with 2 triggers separated by 70 min (1.17h = 0.0486d)
    events_4h = []
    for k in range(5):
        t0 = 10.0 + k * 4.2
        events_4h.append({'center_btjd': t0, 'score': 22.0, 'width_hr': 4.0, 'depth_ppm': 8000.0})
        events_4h.append({'center_btjd': t0 + 0.0486, 'score': 25.0, 'width_hr': 4.0, 'depth_ppm': 8500.0})

    dips_4h = periodicity.cluster_candidate_events(events_4h)
    assert len(dips_4h) == 5, f"Expected 5 dips, got {len(dips_4h)} (should not split 4h transit)"
    res_4h = periodicity.find_periodicity(events_4h)
    assert res_4h.is_periodic
    assert abs(res_4h.period_days - 4.2) < 0.01

    # 6-hour transit with 3 triggers spanning 160 min (2.67h = 0.1111d)
    events_6h = []
    for k in range(4):
        t0 = 10.0 + k * 5.5
        events_6h.append({'center_btjd': t0, 'score': 18.0, 'width_hr': 6.0, 'depth_ppm': 12000.0})
        events_6h.append({'center_btjd': t0 + 0.0556, 'score': 24.0, 'width_hr': 6.0, 'depth_ppm': 13000.0})
        events_6h.append({'center_btjd': t0 + 0.1111, 'score': 19.0, 'width_hr': 6.0, 'depth_ppm': 12500.0})

    dips_6h = periodicity.cluster_candidate_events(events_6h)
    assert len(dips_6h) == 4, f"Expected 4 dips, got {len(dips_6h)} (should not split 6h transit)"
    res_6h = periodicity.find_periodicity(events_6h)
    assert res_6h.is_periodic
    assert abs(res_6h.period_days - 5.5) < 0.01


def test_adaptive_clustering_close_distinct_dips():
    """Verify adaptive clustering keeps closely spaced distinct dips separate."""
    # Eccentric binary with narrow eclipses (w=1.2h) separated by 2.2h (0.0917d)
    events_ecc = []
    for k in range(5):
        t0 = 10.0 + k * 3.5
        events_ecc.append({'center_btjd': t0, 'score': 30.0, 'width_hr': 1.2, 'depth_ppm': 10000.0})
        events_ecc.append({'center_btjd': t0 + 0.0917, 'score': 15.0, 'width_hr': 1.2, 'depth_ppm': 3000.0})

    dips_ecc = periodicity.cluster_candidate_events(events_ecc)
    assert len(dips_ecc) == 10, f"Expected 10 dips (5 primary + 5 secondary), got {len(dips_ecc)}"

    # Dual-planet transit (w=1.5h) separated by 1.8h (0.0750d)
    events_mp = []
    for k in range(4):
        t0 = 5.0 + k * 4.0
        events_mp.append({'center_btjd': t0, 'score': 18.0, 'width_hr': 1.5, 'depth_ppm': 1200.0})
        events_mp.append({'center_btjd': t0 + 0.0750, 'score': 16.0, 'width_hr': 1.5, 'depth_ppm': 900.0})

    dips_mp = periodicity.cluster_candidate_events(events_mp)
    assert len(dips_mp) == 8, f"Expected 8 distinct dips, got {len(dips_mp)}"


def test_adaptive_clustering_order_independence():
    """Verify adaptive clustering output is independent of candidate input ordering."""
    import random
    events = [
        {'center_btjd': 100.0, 'score': 15.0, 'width_hr': 2.0},
        {'center_btjd': 100.02, 'score': 25.0, 'width_hr': 2.0}, # intra-dip duplicate
        {'center_btjd': 103.5, 'score': 20.0, 'width_hr': 2.0},
        {'center_btjd': 107.0, 'score': 18.0, 'width_hr': 2.0},
        {'center_btjd': 107.015, 'score': 30.0, 'width_hr': 2.0}, # intra-dip duplicate
    ]
    u1 = periodicity.cluster_candidate_events(events)

    # Scramble order
    rng = random.Random(42)
    scrambled = list(events)
    rng.shuffle(scrambled)
    u2 = periodicity.cluster_candidate_events(scrambled)

    assert len(u1) == len(u2) == 3
    for d1, d2 in zip(u1, u2):
        assert abs(d1['center_btjd'] - d2['center_btjd']) < 1e-6
        assert d1['score'] == d2['score']



