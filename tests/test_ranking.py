import csv
from pathlib import Path
import sqlite3
import pytest
from tess_anomaly_survey import db
from tess_anomaly_survey.cli import main

def test_star_diverse_export(tmp_path):
    db_file = tmp_path / "diverse.sqlite"
    con = db.connect(db_file)
    run_id = "run_ranking"
    con.execute("INSERT INTO runs VALUES (?,?,?,?,?,?,?)",
                (run_id, "2026-10-13T00:00:00Z", "hash", "fixed_box_2h", 0.01, "{}", "cat.csv"))
    
    # Insert 5 candidates for TIC 111 (all high scores: 50, 49, 48, 47, 46)
    for i, score in enumerate([50.0, 49.0, 48.0, 47.0, 46.0]):
        con.execute("""INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd,
                       score, width_hr, depth_ppm, n_support_windows)
                       VALUES (?,?,?,?,?,?,?,?,?)""",
                    (run_id, "111", 1, "s1.npz", 1000.0 + i, score, 2.0, 1000.0, 1))

    # Insert 2 candidates for TIC 222 (moderate scores: 30, 29)
    for i, score in enumerate([30.0, 29.0]):
        con.execute("""INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd,
                       score, width_hr, depth_ppm, n_support_windows)
                       VALUES (?,?,?,?,?,?,?,?,?)""",
                    (run_id, "222", 1, "s2.npz", 2000.0 + i, score, 2.0, 1000.0, 1))

    # Insert 2 candidates for TIC 333 (moderate scores: 20, 19)
    for i, score in enumerate([20.0, 19.0]):
        con.execute("""INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd,
                       score, width_hr, depth_ppm, n_support_windows)
                       VALUES (?,?,?,?,?,?,?,?,?)""",
                    (run_id, "333", 1, "s3.npz", 3000.0 + i, score, 2.0, 1000.0, 1))
    con.commit()
    con.close()

    # Test 1: Default export (--max-per-star 2, limit 5)
    # Without star-diversity, TIC 111 would monopolize all 5 slots!
    # With default max-per-star=2, TIC 111 gets 2, TIC 222 gets 2, TIC 333 gets 1
    out_default = tmp_path / "export_default.csv"
    main(["export", "--db", str(db_file), "--output", str(out_default), "--limit", "5"])
    rows_default = list(csv.DictReader(out_default.open(encoding="utf-8")))
    assert len(rows_default) == 5
    tics_default = [r["tic"] for r in rows_default]
    assert tics_default.count("111") == 2
    assert tics_default.count("222") == 2
    assert tics_default.count("333") == 1

    # Test 2: Unrestricted export (--max-per-star 0, limit 5)
    # TIC 111 monopolizes all 5 slots because its scores are 50, 49, 48, 47, 46
    out_unrestricted = tmp_path / "export_unrestricted.csv"
    main(["export", "--db", str(db_file), "--output", str(out_unrestricted), "--limit", "5", "--max-per-star", "0"])
    rows_unrestricted = list(csv.DictReader(out_unrestricted.open(encoding="utf-8")))
    assert len(rows_unrestricted) == 5
    tics_unrestricted = [r["tic"] for r in rows_unrestricted]
    assert tics_unrestricted.count("111") == 5

    # Test 3: Strict 1-per-star export (--max-per-star 1, limit 3)
    out_strict = tmp_path / "export_strict.csv"
    main(["export", "--db", str(db_file), "--output", str(out_strict), "--limit", "3", "--max-per-star", "1"])
    rows_strict = list(csv.DictReader(out_strict.open(encoding="utf-8")))
    assert len(rows_strict) == 3
    tics_strict = [r["tic"] for r in rows_strict]
    assert tics_strict == ["111", "222", "333"]


def test_normalized_ranking(tmp_path):
    """Test star-normalized ranking promotes high-significance outliers on quiet stars."""
    db_file = tmp_path / "norm.sqlite"
    con = db.connect(db_file)
    run_id = "run_norm"
    con.execute("INSERT INTO runs VALUES (?,?,?,?,?,?,?)",
                (run_id, "2026-10-14T00:00:00Z", "hash", "box_2h", 0.01, "{}", "cat.csv"))

    # Active noisy star: 10 candidates with elevated scores (mean ~ 80, std ~ 5)
    for i in range(10):
        con.execute("""INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd,
                       score, width_hr, depth_ppm, n_support_windows)
                       VALUES (?,?,?,?,?,?,?,?,?)""",
                    (run_id, "noisy_star", 1, "s1.npz", 1000.0 + i, 80.0 + (i % 3) * 2.0, 2.0, 5000.0, 1))

    # Quiet star: 5 candidates, mostly low scores (~10), but ONE prominent dip of score 35 (z-score >> 5)
    for i, sc in enumerate([10.0, 10.5, 9.8, 10.2, 35.0]):
        con.execute("""INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd,
                       score, width_hr, depth_ppm, n_support_windows)
                       VALUES (?,?,?,?,?,?,?,?,?)""",
                    (run_id, "quiet_star", 1, "s2.npz", 2000.0 + i, sc, 2.0, 1000.0, 1))
    con.commit()
    con.close()

    # Under raw ranking, noisy_star candidates (scores 80-84) will beat quiet_star dip (score 35)
    out_raw = tmp_path / "export_raw.csv"
    main(["export", "--db", str(db_file), "--output", str(out_raw), "--ranking", "raw", "--limit", "2", "--max-per-star", "1"])
    rows_raw = list(csv.DictReader(out_raw.open(encoding="utf-8")))
    assert rows_raw[0]["tic"] == "noisy_star"

    # Under normalized ranking, quiet_star's prominent 35 score is a massive outlier on its star!
    out_norm = tmp_path / "export_norm.csv"
    main(["export", "--db", str(db_file), "--output", str(out_norm), "--ranking", "normalized", "--limit", "2", "--max-per-star", "1"])
    rows_norm = list(csv.DictReader(out_norm.open(encoding="utf-8")))
    assert rows_norm[0]["tic"] == "quiet_star"
    assert float(rows_norm[0]["normalized_score"]) > float(rows_norm[1]["normalized_score"])


def test_group_by_star_export(tmp_path):
    """Test group-by-star export aggregation."""
    db_file = tmp_path / "star_grp.sqlite"
    con = db.connect(db_file)
    run_id = "run_star_grp"
    con.execute("INSERT INTO runs VALUES (?,?,?,?,?,?,?)",
                (run_id, "2026-10-14T00:00:00Z", "hash", "box_2h", 0.01, "{}", "cat.csv"))

    for i in range(4):
        con.execute("""INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd,
                       score, width_hr, depth_ppm, n_support_windows)
                       VALUES (?,?,?,?,?,?,?,?,?)""",
                    (run_id, "starA", 1, "s1.npz", 1000.0 + i, 20.0 + i, 2.0, 1000.0, 1))

    for i in range(2):
        con.execute("""INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd,
                       score, width_hr, depth_ppm, n_support_windows)
                       VALUES (?,?,?,?,?,?,?,?,?)""",
                    (run_id, "starB", 1, "s2.npz", 2000.0 + i, 40.0 + i, 3.0, 2000.0, 1))
    con.commit()
    con.close()

    out_star = tmp_path / "export_stars.csv"
    main(["export", "--db", str(db_file), "--output", str(out_star), "--group-by-star"])
    rows = list(csv.DictReader(out_star.open(encoding="utf-8")))
    assert len(rows) == 2
    # starB has higher max_score (41 vs 23)
    assert rows[0]["tic"] == "starB"
    assert int(rows[0]["candidate_count"]) == 2
    assert float(rows[0]["max_score"]) == 41.0
    assert rows[1]["tic"] == "starA"
    assert int(rows[1]["candidate_count"]) == 4


def test_group_by_system_export(tmp_path):
    """Test export of candidate systems grouped by periodic system."""
    db_file = tmp_path / "sys_grp.sqlite"
    con = db.connect(db_file)
    run_id = "run_sys_grp"
    con.execute("INSERT INTO runs VALUES (?,?,?,?,?,?,?)",
                (run_id, "2026-10-14T00:00:00Z", "hash", "box_2h", 0.01, "{}", "cat.csv"))

    con.execute("""INSERT INTO candidate_systems(run_id, tic, n_events, period_days,
                   epoch_btjd, period_snr, is_periodic, odd_even_sigma, classification,
                   note, period_uncertainty_days, epoch_uncertainty_days,
                   alternate_period_days, created_utc)
                   VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                (run_id, "TIC_P1", 5, 3.25, 1400.1, 85.0, 1, 1.2, "periodic",
                 "Candidate periodic", 0.0003, 0.001, 1.625, "2026-10-14T00:00:00Z"))

    con.execute("""INSERT INTO candidate_systems(run_id, tic, n_events, period_days,
                   epoch_btjd, period_snr, is_periodic, odd_even_sigma, classification,
                   note, period_uncertainty_days, epoch_uncertainty_days,
                   alternate_period_days, created_utc)
                   VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                (run_id, "TIC_NP", 2, 0.0, 1400.1, 15.0, 0, 0.0, "no_periodicity_detected",
                 "Insufficient events", 0.0, 0.0, 0.0, "2026-10-14T00:00:00Z"))
    con.commit()
    con.close()

    out_sys = tmp_path / "export_systems.csv"
    main(["export", "--db", str(db_file), "--output", str(out_sys), "--group-by-system"])
    rows = list(csv.DictReader(out_sys.open(encoding="utf-8")))
    assert len(rows) == 2
    assert rows[0]["tic"] == "TIC_P1"
    assert rows[0]["is_periodic"] == "1"
    assert float(rows[0]["period_days"]) == 3.25
    assert float(rows[0]["period_uncertainty_days"]) == 0.0003
    assert rows[1]["tic"] == "TIC_NP"
    assert rows[1]["is_periodic"] == "0"

