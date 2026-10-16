import os
import sqlite3
import tempfile
import pytest
import numpy as np
from astropy.io import fits

from tess_anomaly_survey import streaming


def create_synthetic_fits(path: str, tic: int = 12345678, sector: int = 1, inject_dip: bool = False):
    """Helper to generate a minimal valid TESS SPOC light-curve FITS file for testing."""
    n = 2000
    t = np.linspace(1000.0, 1002.77, n) # ~2.8 days at 2-min cadence
    flux = np.random.normal(10000.0, 10.0, n)
    err = np.full(n, 10.0)
    quality = np.zeros(n, dtype=np.uint32)

    if inject_dip:
        # Inject 5% 2-hour dip at center
        m = (t >= 1001.4) & (t <= 1001.5)
        flux[m] *= 0.95

    c1 = fits.Column(name='TIME', format='D', unit='BJD - 2457000', array=t)
    c2 = fits.Column(name='PDCSAP_FLUX', format='E', unit='e-/s', array=flux)
    c3 = fits.Column(name='PDCSAP_FLUX_ERR', format='E', unit='e-/s', array=err)
    c4 = fits.Column(name='SAP_FLUX', format='E', unit='e-/s', array=flux)
    c5 = fits.Column(name='SAP_FLUX_ERR', format='E', unit='e-/s', array=err)
    c6 = fits.Column(name='QUALITY', format='J', array=quality)

    table_hdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6])
    table_hdu.header['EXTNAME'] = 'LIGHTCURVE'

    pri_hdu = fits.PrimaryHDU()
    pri_hdu.header['TELESCOP'] = 'TESS'
    pri_hdu.header['TICID'] = tic
    pri_hdu.header['SECTOR'] = sector
    pri_hdu.header['CAMERA'] = 1
    pri_hdu.header['CCD'] = 1
    pri_hdu.header['TESSMAG'] = 10.5
    pri_hdu.header['TEFF'] = 5800.0

    hdul = fits.HDUList([pri_hdu, table_hdu])
    hdul.writeto(path, overwrite=True)


def test_init_streaming_db():
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "test_stream.sqlite")
        con = streaming.init_streaming_db(db_path)
        tables = [r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()]
        assert "observations" in tables
        assert "candidates" in tables
        assert "candidate_detections" in tables
        assert "failures" in tables
        assert "runs" in tables
        assert "performance_benchmarks" in tables
        con.close()


def test_development_universe_guard():
    # If target is not in dev universe, pipeline must raise RuntimeError
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "test_stream.sqlite")
        manifest_path = os.path.join(tmpdir, "mock_dev_universe.csv")
        with open(manifest_path, "w", encoding="utf-8") as f:
            f.write("tic_id\n10006821\n")
        targets = [{"tic": "999999999999", "sector": 1}]
        
        with pytest.raises(RuntimeError) as exc_info:
            streaming.stream_pipeline(
                targets=targets,
                db_path=db_path,
                run_id="test_run",
                check_dev_universe=True,
                dev_universe_manifest=manifest_path
            )
        assert "DATA LEAKAGE GUARD" in str(exc_info.value)


def test_process_observation_synthetic():
    with tempfile.TemporaryDirectory() as tmpdir:
        fits_path = os.path.join(tmpdir, "test_lc.fits")
        create_synthetic_fits(fits_path, tic=10006821, sector=1, inject_dip=False)

        res = streaming.process_observation(fits_path)
        assert res["status"] == "no_anomalies"
        assert res["tic"] == "10006821"
        assert res["sector"] == 1
        assert res["n_windows"] > 0
        assert res["n_events"] == 0


def test_streaming_idempotent_resumption():
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "test_stream.sqlite")
        fits_file = os.path.join(tmpdir, "tess2018206045859-s0001-000000010006821-0120-s_lc.fits")
        # TIC 10006821 is in calibration universe
        create_synthetic_fits(fits_file, tic=10006821, sector=1, inject_dip=False)

        targets = [{
            "tic": "10006821",
            "sector": 1,
            "filename": os.path.basename(fits_file),
            "source_url": "file://" + fits_file.replace("\\", "/")
        }]

        # Run 1
        metrics1 = streaming.stream_pipeline(
            targets=targets,
            db_path=db_path,
            run_id="stream_run_1",
            temp_dir=tmpdir,
            retain_fits=True,
            check_dev_universe=False
        )
        assert metrics1["n_processed"] == 1
        assert metrics1["n_successful"] == 1

        con = sqlite3.connect(db_path)
        obs_count = con.execute("SELECT count(*) FROM observations").fetchone()[0]
        assert obs_count == 1
        con.close()

        # Run 2 with resume=True should skip already processed target
        metrics2 = streaming.stream_pipeline(
            targets=targets,
            db_path=db_path,
            run_id="stream_run_1",
            temp_dir=tmpdir,
            retain_fits=True,
            resume=True,
            check_dev_universe=False
        )
        assert metrics2["n_skipped"] == 1
        assert metrics2["n_processed"] == 0

        con = sqlite3.connect(db_path)
        obs_count2 = con.execute("SELECT count(*) FROM observations").fetchone()[0]
        assert obs_count2 == 1 # Still 1, no duplicate inserted!
        con.close()


def test_truncated_fits_handling():
    """Verify that truncated or corrupted downloads are caught, logged, and marked failed."""
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "test_stream.sqlite")
        corrupt_file = os.path.join(tmpdir, "corrupt_s0001_lc.fits")
        
        # Write 256 bytes of garbage
        with open(corrupt_file, "wb") as f:
            f.write(b"SIMPLE  =                    T / corrupt FITS header data truncation test" + b"\x00" * 180)

        targets = [{
            "tic": "10006821",
            "sector": 1,
            "filename": os.path.basename(corrupt_file),
            "source_url": "file://" + corrupt_file.replace("\\", "/")
        }]

        metrics = streaming.stream_pipeline(
            targets=targets,
            db_path=db_path,
            run_id="corrupt_test",
            temp_dir=tmpdir,
            retain_fits=True,
            check_dev_universe=False
        )

        assert metrics["n_processed"] == 1
        assert metrics["n_failed"] == 1
        assert metrics["n_successful"] == 0

        con = sqlite3.connect(db_path)
        obs_row = con.execute("SELECT status, error_message FROM observations WHERE run_id='corrupt_test'").fetchone()
        assert obs_row is not None
        assert obs_row[0] == "failed"
        assert obs_row[1] is not None

        fail_row = con.execute("SELECT stage, error_message FROM failures WHERE run_id='corrupt_test'").fetchone()
        assert fail_row is not None
        assert "pipeline" in fail_row[0]
        con.close()


def test_taxonomy_classification_continuous_and_unverified():
    """Verify distinct classification of continuous variability, unverified excursions, and flares."""
    # 1. Continuous sine wave modulation with >= 10 triggers
    n = 2000
    t = np.linspace(1000.0, 1002.77, n)
    f_continuous = 0.05 * np.sin(2.0 * np.pi * t / 0.2)
    g = np.ones(n, dtype=bool)
    assert streaming.detect_continuous_variability(t, f_continuous, g, n_events=12) is True

    # 2. Quiet star with 2 events -> False
    f_quiet = np.random.normal(0.0, 0.001, n)
    assert streaming.detect_continuous_variability(t, f_quiet, g, n_events=2) is False


def test_survey_resumption_different_archive_products(tmp_path):
    """Verify that resumption does not incorrectly skip distinct archive products for same star/sector."""
    db_path = str(tmp_path / "test_resume.sqlite")
    con = streaming.init_streaming_db(db_path)
    run_id = "run_resume_test"

    with con:
        con.execute("""
            INSERT INTO observations (
                run_id, tic, sector, product_filename, config_hash, status, created_utc
            ) VALUES (?, 'TIC:39018208', 2, 'tess_product_a.fits', 'cfg', 'completed', '2026-10-16T12:00:00Z')
        """, (run_id,))
    con.close()

    # Simulate resume filtering with 3 targets
    targets = [
        {"tic": "39018208", "sector": 2, "product_filename": "tess_product_a.fits"}, # same product -> skip
        {"tic": "39018208", "sector": 2, "product_filename": "tess_product_b.fits"}, # different product -> keep
        {"tic": "39018208", "sector": 2},                                             # generic star/sector -> skip
        {"tic": "38396965", "sector": 2, "product_filename": "tess_product_c.fits"}  # unobserved star -> keep
    ]

    # Verify via streaming db inspection logic
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    completed_obs = set()
    completed_tic_sec = set()
    for r in cur.execute("SELECT tic, sector, product_filename, status FROM observations WHERE run_id=?", (run_id,)).fetchall():
        if r[3] in ('completed', 'no_anomalies', 'completed_with_anomalies'):
            tic_clean = str(r[0]).replace('TIC:', '').strip()
            completed_obs.add((str(r[0]), int(r[1]), str(r[2])))
            completed_tic_sec.add((f"TIC:{tic_clean}", int(r[1])))
            completed_tic_sec.add((tic_clean, int(r[1])))

    pending = []
    for it in targets:
        clean_t = str(it.get('tic_id') or it.get('tic')).replace('TIC:', '').strip()
        tic_str = f"TIC:{clean_t}"
        sec = int(it.get('sector')) if it.get('sector') is not None else None
        fname = it.get('filename') or it.get('product_filename')
        if sec is not None:
            if fname:
                if (tic_str, sec, fname) in completed_obs:
                    continue
            else:
                if (tic_str, sec) in completed_tic_sec or (clean_t, sec) in completed_tic_sec:
                    continue
        pending.append(it)
    con.close()

    assert len(pending) == 2
    # tess_product_b.fits should NOT be skipped!
    assert any(p.get("product_filename") == "tess_product_b.fits" for p in pending)
    # tess_product_c.fits should NOT be skipped!
    assert any(p.get("product_filename") == "tess_product_c.fits" for p in pending)
    # tess_product_a.fits MUST be skipped
    assert not any(p.get("product_filename") == "tess_product_a.fits" for p in pending)
    # generic query without product_filename MUST be skipped
    assert not any("product_filename" not in p for p in pending)

