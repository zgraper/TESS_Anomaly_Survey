import sqlite3
import pytest
from tess_anomaly_survey import db

def test_duplicate_merging_same_physical_event():
    """Demonstrate that under existing logic, two detections of the same 2-hour dip

    differing by 40 minutes (well within the 2-hour event duration) fail to merge
    and erroneously create two separate candidate records.
    """
    con = db.connect(":memory:")
    run_id = "test_run"
    con.execute("INSERT INTO runs VALUES (?,?,?,?,?,?,?)",
                (run_id, "2026-10-13T00:00:00Z", "hash", "fixed_box_2h", 0.01, "{}", "cat.csv"))
    cadence = 2.0 / 1440.0  # 2 minutes in days (~0.001389 d)
    row = {"tic": "12345678", "sector": 1}
    source = "test_lc.npz"

    # Detection 1 from Window 1: 2-hour event centered at 1000.485 BTJD
    det1 = {
        "center_btjd": 1000.485,
        "width_hr": 2.0,
        "depth_ppm": 1500.0,
        "score": 8.5
    }
    db.upsert_candidate(con, run_id, row, source, det1, cadence)

    # Detection 2 from Window 2 (overlapping): same 2-hour event centered at 1000.515 BTJD
    # Time separation is 0.030 days = 43.2 minutes, well within the 2-hour (120 min) transit!
    det2 = {
        "center_btjd": 1000.515,
        "width_hr": 2.0,
        "depth_ppm": 1600.0,
        "score": 9.2
    }
    db.upsert_candidate(con, run_id, row, source, det2, cadence)

    candidates = con.execute("SELECT id, center_btjd, width_hr, score, n_support_windows FROM candidates").fetchall()
    
    # We expect these two overlapping detections of the same physical event to be merged into ONE candidate
    # with 2 supporting windows.
    assert len(candidates) == 1, f"Expected 1 merged candidate for the same event, got {len(candidates)}: {candidates}"
    assert candidates[0][4] == 2, f"Expected n_support_windows=2, got {candidates[0][4]}"
    assert candidates[0][3] == 9.2, f"Expected higher score (9.2) to be retained, got {candidates[0][3]}"

    # Verify that BOTH original detection records are preserved in candidate_detections table
    detections = con.execute("SELECT center_btjd, score, width_hr FROM candidate_detections WHERE candidate_id=?",
                             (candidates[0][0],)).fetchall()
    assert len(detections) == 2, f"Expected 2 detection records in candidate_detections, got {len(detections)}"
    assert detections[0] == (1000.485, 8.5, 2.0)
    assert detections[1] == (1000.515, 9.2, 2.0)


def test_distinct_events_not_excessively_merged():
    """Verify that two distinct transit dips (separated by more than duration) are NOT merged."""
    con = db.connect(":memory:")
    run_id = "test_run"
    con.execute("INSERT INTO runs VALUES (?,?,?,?,?,?,?)",
                (run_id, "2026-10-13T00:00:00Z", "hash", "fixed_box_2h", 0.01, "{}", "cat.csv"))
    cadence = 2.0 / 1440.0
    row = {"tic": "12345678", "sector": 1}
    source = "test_lc.npz"

    # Event 1: 1-hour dip at 1000.0 BTJD
    det1 = {"center_btjd": 1000.0, "width_hr": 1.0, "depth_ppm": 2000.0, "score": 10.0}
    db.upsert_candidate(con, run_id, row, source, det1, cadence)

    # Event 2: 1-hour dip at 1000.25 BTJD (6 hours later - clearly distinct event!)
    det2 = {"center_btjd": 1000.25, "width_hr": 1.0, "depth_ppm": 1800.0, "score": 9.5}
    db.upsert_candidate(con, run_id, row, source, det2, cadence)

    candidates = con.execute("SELECT id, center_btjd FROM candidates").fetchall()
    assert len(candidates) == 2, f"Expected 2 separate candidates for distinct events, got {len(candidates)}"


def test_backup_created_on_schema_migration(tmp_path):
    """Verify that an existing database without candidate_detections triggers a backup before migration."""
    db_file = tmp_path / "legacy.sqlite"
    # Create legacy database with v0.1 schema (without candidate_detections)
    legacy_con = sqlite3.connect(str(db_file))
    legacy_con.executescript("""
        CREATE TABLE runs (id TEXT PRIMARY KEY, created_utc TEXT, calibration_sha256 TEXT,
                           method TEXT, target_alert_fraction REAL, settings_json TEXT, catalog_path TEXT);
        CREATE TABLE processed_files (run_id TEXT, source_file TEXT, tic TEXT, sector INTEGER,
                                      n_windows INTEGER, n_alert_windows INTEGER, status TEXT, error TEXT,
                                      PRIMARY KEY(run_id, source_file));
        CREATE TABLE windows (run_id TEXT, source_file TEXT, tic TEXT, start_index INTEGER, start_btjd REAL,
                              noise REAL, max_score REAL, cutoff REAL, noise_bin INTEGER, flagged INTEGER,
                              PRIMARY KEY(run_id, source_file, start_index));
        CREATE TABLE candidates (id INTEGER PRIMARY KEY AUTOINCREMENT, run_id TEXT, tic TEXT, sector INTEGER,
                                 source_file TEXT, center_btjd REAL, score REAL, width_hr REAL, depth_ppm REAL,
                                 n_support_windows INTEGER DEFAULT 1, review_status TEXT DEFAULT 'unreviewed', note TEXT DEFAULT '');
    """)
    legacy_con.execute("INSERT INTO runs VALUES ('r1', '2026-10-12', 'h', 'm', 0.01, '{}', 'c')")
    legacy_con.commit()
    legacy_con.close()

    # Now open with db.connect
    migrated_con = db.connect(db_file)
    # Check that candidate_detections table was created
    table = migrated_con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='candidate_detections'").fetchone()
    assert table is not None
    migrated_con.close()

    # Verify a backup file was generated in the directory
    backups = list(tmp_path.glob("legacy.sqlite.bak_*"))
    assert len(backups) == 1, f"Expected 1 backup file, found {len(backups)}: {backups}"


def test_repeat_scans_do_not_duplicate_detections():
    """Verify that scanning the exact same detection twice does not add duplicates or inflate support count."""
    con = db.connect(":memory:")
    run_id = "test_run"
    con.execute("INSERT INTO runs VALUES (?,?,?,?,?,?,?)",
                (run_id, "2026-10-13T00:00:00Z", "hash", "fixed_box_2h", 0.01, "{}", "cat.csv"))
    cadence = 2.0 / 1440.0
    row = {"tic": "99999", "sector": 1}
    source = "lc.npz"

    det = {"center_btjd": 1000.5, "width_hr": 2.0, "depth_ppm": 1200.0, "score": 8.0}
    db.upsert_candidate(con, run_id, row, source, det, cadence)
    # Re-run same scan / detection
    db.upsert_candidate(con, run_id, row, source, det, cadence)

    cand = con.execute("SELECT id, n_support_windows FROM candidates").fetchall()
    assert len(cand) == 1
    assert cand[0][1] == 1, f"Expected n_support_windows=1, got {cand[0][1]}"

    dets = con.execute("SELECT id FROM candidate_detections WHERE candidate_id=?", (cand[0][0],)).fetchall()
    assert len(dets) == 1, f"Expected 1 detection row, got {len(dets)}"


def test_processing_order_independence_chained_detections():
    """Verify that detections chaining an event are merged identically regardless of arrival order."""
    cadence = 2.0 / 1440.0
    row = {"tic": "88888", "sector": 1}
    source = "lc.npz"

    # Three overlapping detections: left (1000.00), middle (1000.03), right (1000.06)
    # Left and Right differ by 0.06d (86 min > 60 min horizon), but Middle bridges them!
    det_a = {"center_btjd": 1000.00, "width_hr": 2.0, "depth_ppm": 1000.0, "score": 8.0}
    det_b = {"center_btjd": 1000.03, "width_hr": 2.0, "depth_ppm": 1200.0, "score": 10.0}
    det_c = {"center_btjd": 1000.06, "width_hr": 2.0, "depth_ppm": 1100.0, "score": 9.0}

    # Order 1: sequential [A, B, C]
    con1 = db.connect(":memory:")
    con1.execute("INSERT INTO runs VALUES ('r1', 'now', 'h', 'm', 0.01, '{}', 'c')")
    for d in [det_a, det_b, det_c]:
        db.upsert_candidate(con1, "r1", row, source, d, cadence)
    cand1 = con1.execute("SELECT center_btjd, score, n_support_windows FROM candidates").fetchall()
    dets1 = con1.execute("SELECT count(*) FROM candidate_detections").fetchone()[0]

    # Order 2: non-sequential [A, C, B] where Middle arrives last and bridges Left and Right!
    con2 = db.connect(":memory:")
    con2.execute("INSERT INTO runs VALUES ('r2', 'now', 'h', 'm', 0.01, '{}', 'c')")
    for d in [det_a, det_c, det_b]:
        db.upsert_candidate(con2, "r2", row, source, d, cadence)
    cand2 = con2.execute("SELECT center_btjd, score, n_support_windows FROM candidates").fetchall()
    dets2 = con2.execute("SELECT count(*) FROM candidate_detections").fetchone()[0]

    assert len(cand1) == 1, f"Expected 1 merged candidate in Order 1, got {len(cand1)}"
    assert len(cand2) == 1, f"Expected 1 merged candidate in Order 2, got {len(cand2)}"
    assert cand1[0][1] == cand2[0][1] == 10.0  # Peak score preserved
    assert cand1[0][2] == cand2[0][2] == 3     # 3 supporting windows in both
    assert dets1 == dets2 == 3                 # All 3 individual detections preserved in both


def test_safe_wal_mode_backup(tmp_path):
    """Verify that backup_database works correctly on a database in WAL mode with active transactions."""
    db_file = tmp_path / "wal_test.sqlite"
    con = sqlite3.connect(str(db_file))
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("CREATE TABLE test_data (id INT PRIMARY KEY, val TEXT)")
    con.execute("INSERT INTO test_data VALUES (1, 'alpha'), (2, 'beta')")
    con.commit()

    # Create safe backup while connection is open
    bak = db.backup_database(db_file)
    assert bak.exists()

    # Verify backup is a valid standalone SQLite database containing the data
    bak_con = sqlite3.connect(str(bak))
    rows = bak_con.execute("SELECT id, val FROM test_data ORDER BY id").fetchall()
    assert rows == [(1, "alpha"), (2, "beta")]
    bak_con.close()
    con.close()


