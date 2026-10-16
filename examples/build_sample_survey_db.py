"""Generate self-contained public example database and light curve data for GitHub users.

Includes:
1. TIC 39018208: Illustrative synthetic periodic eclipsing binary benchmark (simulated P = 2.41824 d).
   Note: Provided strictly as an illustrative synthetic software verification model;
   clearly separated from cataloged EBLM J0048-66 (whose true orbital period is P approx. 6.64927 d).
2. TIC 73834: Localized dimming candidate (aperiodic / disk occultation candidate).
3. TIC 60232: Stellar pulsation variable (classified as continuous_stellar_variability).
4. TIC 141944: Flare phantom rejected by two-pass verifier (flare_induced_apparent_dip).
5. TIC 999001: Quiet control star (no detected anomalies).
"""
import json
import os
from pathlib import Path
import sqlite3
import numpy as np

EXAMPLES_DIR = Path(__file__).resolve().parent
SAMPLE_DB_PATH = EXAMPLES_DIR / "sample_survey_review.sqlite"
SAMPLE_DATA_DIR = EXAMPLES_DIR / "sample_data"

def create_sample_database():
    SAMPLE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    if SAMPLE_DB_PATH.exists():
        SAMPLE_DB_PATH.unlink()

    con = sqlite3.connect(str(SAMPLE_DB_PATH))
    with con:
        con.executescript("""
        PRAGMA foreign_keys=ON;

        CREATE TABLE runs (
            run_id TEXT PRIMARY KEY,
            created_utc TEXT NOT NULL,
            software_version TEXT NOT NULL,
            b0_hash TEXT NOT NULL,
            b1_hash TEXT NOT NULL,
            threshold_t REAL NOT NULL,
            method TEXT NOT NULL,
            config_hash TEXT NOT NULL,
            worker_count INTEGER NOT NULL DEFAULT 1,
            settings_json TEXT NOT NULL,
            notes TEXT
        );

        CREATE TABLE observations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            run_id TEXT NOT NULL,
            tic TEXT NOT NULL,
            sector INTEGER NOT NULL,
            camera INTEGER,
            ccd INTEGER,
            tmag REAL,
            teff REAL,
            obs_id TEXT,
            product_filename TEXT NOT NULL,
            data_uri TEXT,
            source_url TEXT,
            file_sha256 TEXT,
            calib_level INTEGER DEFAULT 2,
            prvversion TEXT,
            config_hash TEXT NOT NULL,
            status TEXT NOT NULL,
            n_cadences INTEGER NOT NULL DEFAULT 0,
            n_usable_cadences INTEGER NOT NULL DEFAULT 0,
            n_windows INTEGER NOT NULL DEFAULT 0,
            n_alert_windows INTEGER NOT NULL DEFAULT 0,
            n_events INTEGER NOT NULL DEFAULT 0,
            n_verified_events INTEGER NOT NULL DEFAULT 0,
            download_time_sec REAL NOT NULL DEFAULT 0.0,
            process_time_sec REAL NOT NULL DEFAULT 0.0,
            total_time_sec REAL NOT NULL DEFAULT 0.0,
            file_bytes INTEGER NOT NULL DEFAULT 0,
            error_message TEXT,
            created_utc TEXT NOT NULL
        );

        CREATE TABLE candidates (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            observation_id INTEGER NOT NULL,
            run_id TEXT NOT NULL,
            tic TEXT NOT NULL,
            sector INTEGER NOT NULL,
            source_file TEXT NOT NULL,
            center_btjd REAL NOT NULL,
            score REAL NOT NULL,
            width_hr REAL NOT NULL,
            depth_ppm REAL NOT NULL,
            n_support_windows INTEGER NOT NULL,
            raw_outcome TEXT NOT NULL,
            raw_reason TEXT NOT NULL,
            scientific_classification TEXT NOT NULL,
            scientific_description TEXT NOT NULL,
            true_depth_min_pct REAL NOT NULL,
            true_depth_mean_pct REAL NOT NULL,
            consecutive_depressed INTEGER NOT NULL,
            has_shoulder_flare INTEGER NOT NULL,
            shoulder_elevation_pct REAL NOT NULL,
            review_status TEXT NOT NULL DEFAULT 'unreviewed',
            note TEXT NOT NULL DEFAULT '',
            system_id INTEGER,
            created_utc TEXT NOT NULL,
            FOREIGN KEY(observation_id) REFERENCES observations(id)
        );

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
            period_uncertainty_days REAL,
            epoch_uncertainty_days REAL,
            alternate_period_days REAL,
            competing_solutions_json TEXT,
            created_utc TEXT NOT NULL
        );

        CREATE TABLE candidate_detections (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            candidate_id INTEGER NOT NULL,
            run_id TEXT NOT NULL,
            window_start_btjd REAL NOT NULL,
            window_end_btjd REAL NOT NULL,
            center_btjd REAL NOT NULL,
            score REAL NOT NULL,
            width_hr REAL NOT NULL,
            depth_ppm REAL NOT NULL,
            created_utc TEXT NOT NULL,
            FOREIGN KEY(candidate_id) REFERENCES candidates(id)
        );

        CREATE TABLE star_reviews (
            tic TEXT PRIMARY KEY,
            status TEXT NOT NULL DEFAULT 'unreviewed',
            star_class TEXT NOT NULL DEFAULT '',
            note TEXT NOT NULL DEFAULT '',
            reviewer TEXT NOT NULL DEFAULT 'local_user',
            updated_utc TEXT NOT NULL
        );

        CREATE TABLE star_review_history (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            tic TEXT NOT NULL,
            previous_status TEXT NOT NULL,
            new_status TEXT NOT NULL,
            note TEXT NOT NULL,
            reviewer TEXT NOT NULL DEFAULT 'local_user',
            created_utc TEXT NOT NULL
        );

        CREATE TABLE review_history (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            candidate_id INTEGER NOT NULL,
            previous_status TEXT NOT NULL,
            new_status TEXT NOT NULL,
            note TEXT NOT NULL,
            reviewer TEXT NOT NULL DEFAULT 'local_user',
            created_utc TEXT NOT NULL,
            FOREIGN KEY(candidate_id) REFERENCES candidates(id)
        );

        CREATE TABLE performance_benchmarks (
            run_id TEXT PRIMARY KEY,
            worker_count INTEGER NOT NULL DEFAULT 2,
            total_observations INTEGER NOT NULL,
            successful_observations INTEGER NOT NULL,
            failed_observations INTEGER NOT NULL,
            total_download_bytes INTEGER NOT NULL,
            total_download_sec REAL NOT NULL,
            total_process_sec REAL NOT NULL,
            total_wall_sec REAL NOT NULL,
            avg_obs_per_hour REAL NOT NULL,
            avg_download_mb_per_sec REAL NOT NULL,
            peak_ram_mb REAL NOT NULL,
            avg_cpu_pct REAL NOT NULL,
            retry_count INTEGER NOT NULL DEFAULT 0,
            created_utc TEXT NOT NULL
        );
        """)

        # 1. Run record
        run_id = "sample_survey_run"
        con.execute("""INSERT INTO runs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", (
            run_id, "2026-10-15T12:00:00Z", "0.1.0",
            "9c8afe63ec2763b9db57c828b52b0e7d262bb30b203f86e559836dd5beb89f00",
            "c5214809a053e4705c5d86f29703e7682128e3db119494ed59f6c9a191e557c8",
            73.77965536097444, "fixed_box_2h", "frozen_b1_hash", 2,
            json.dumps({"target_alert_fraction": 0.01, "duration_hr": 2.0}),
            "Public demonstration dataset containing benchmark astrophysical systems and controls."
        ))

        # 2. Benchmark telemetry
        con.execute("""INSERT INTO performance_benchmarks VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", (
            run_id, 2, 5, 5, 0, 85000000, 6.2, 2.8, 9.0, 6365.2, 14.2, 482.0, 44.5, 0, "2026-10-15T12:05:00Z"
        ))

        # 3. Observation 1: TIC 39018208 (Sector 2) - Illustrative Synthetic Eclipsing Binary Benchmark
        con.execute("""INSERT INTO observations(id, run_id, tic, sector, camera, ccd, tmag, teff, obs_id,
            product_filename, config_hash, status, n_cadences, n_usable_cadences, n_windows, n_alert_windows,
            n_events, n_verified_events, created_utc)
            VALUES (1, ?, '39018208', 2, 3, 2, 10.42, 5820.0, 'tess2018235142541-s0002-0000000039018208-0121-s',
            'tess2018235142541-s0002-0000000039018208-0121-s_lc.fits', 'cfg_hash', 'completed',
            19737, 18450, 75, 4, 4, 4, '2026-10-15T12:01:00Z')""", (run_id,))

        # System record for TIC 39018208
        con.execute("""INSERT INTO candidate_systems(id, run_id, tic, n_events, period_days, epoch_btjd,
            period_snr, is_periodic, odd_even_sigma, classification, note, period_uncertainty_days,
            epoch_uncertainty_days, created_utc)
            VALUES (1, ?, '39018208', 4, 2.41824, 1356.4180, 48.5, 1, 0.2, 'eclipsing_binary',
            'Illustrative synthetic eclipsing binary benchmark (simulated P = 2.41824 d, simulated T0 = 1356.4180 BTJD; synthetic numerical model for software verification, not an observational measurement of EBLM J0048-66)', 0.00008, 0.0005, '2026-10-15T12:01:05Z')""", (run_id,))

        # Candidates on TIC 39018208
        events_eblm = [
            (1, 1356.4180, 124.52, 4820.0, "interesting", "Synthetic demonstrative secondary eclipse #1 (simulated model)"),
            (2, 1358.8362, 118.20, 4790.0, "unreviewed", "Synthetic demonstrative secondary eclipse #2 (simulated model)"),
            (3, 1361.2544, 122.90, 4850.0, "unreviewed", "Synthetic demonstrative secondary eclipse #3 (simulated model)"),
            (4, 1363.6726, 120.40, 4810.0, "unreviewed", "Synthetic demonstrative secondary eclipse #4 (simulated model)")
        ]
        for cid, center_t, score, depth, st, note in events_eblm:
            con.execute("""INSERT INTO candidates(id, observation_id, run_id, tic, sector, source_file,
                center_btjd, score, width_hr, depth_ppm, n_support_windows, raw_outcome, raw_reason,
                scientific_classification, scientific_description, true_depth_min_pct, true_depth_mean_pct,
                consecutive_depressed, has_shoulder_flare, shoulder_elevation_pct, review_status, note,
                system_id, created_utc)
                VALUES (?, 1, ?, '39018208', 2, 'lc_39018208_s02.npz', ?, ?, 2.0, ?, 2,
                'verified_localized_dimming', 'pass', 'verified_localized_dimming',
                'Synthetic demonstrative secondary eclipse dip (simulated model, S_synth = ' || ? || ')', 0.482, 0.479, 12, 0, 0.0, ?, ?, 1, '2026-10-15T12:01:10Z')""",
                (cid, run_id, center_t, score, depth, score, st, note))
            con.execute("""INSERT INTO candidate_detections(candidate_id, run_id, window_start_btjd,
                window_end_btjd, center_btjd, score, width_hr, depth_ppm, created_utc)
                VALUES (?, ?, ?, ?, ?, ?, 2.0, ?, '2026-10-15T12:01:10Z')""",
                (cid, run_id, center_t - 0.25, center_t + 0.25, center_t, score, depth))

        # Star review for TIC 39018208
        con.execute("""INSERT INTO star_reviews(tic, status, star_class, note, reviewer, updated_utc)
            VALUES ('39018208', 'interesting', 'eclipsing_binary',
            'Illustrative synthetic benchmark target for research workspace evaluation (simulated P = 2.41824 d, synthetic depth = 0.482%). Provided strictly as a synthetic software verification model; not an observational measurement of EBLM J0048-66 (whose true cataloged orbital period is P approx. 6.64927 d).',
            'survey_pi', '2026-10-15T12:02:00Z')""")

        # 4. Observation 2: TIC 73834 (Sector 8) - Deep Localized Dimming Candidate
        con.execute("""INSERT INTO observations(id, run_id, tic, sector, camera, ccd, tmag, teff, obs_id,
            product_filename, config_hash, status, n_cadences, n_usable_cadences, n_windows, n_alert_windows,
            n_events, n_verified_events, created_utc)
            VALUES (2, ?, '73834', 8, 2, 4, 11.15, 6100.0, 'tess2019032160000-s0008-0000000000073834-0136-s',
            'tess2019032160000-s0008-0000000000073834-0136-s_lc.fits', 'cfg_hash', 'completed',
            19800, 18900, 76, 2, 2, 2, '2026-10-15T12:02:30Z')""", (run_id,))

        events_73834 = [
            (5, 1518.2500, 88.40, 7500.0, "unreviewed", ""),
            (6, 1532.1000, 81.15, 7200.0, "unreviewed", "")
        ]
        for cid, center_t, score, depth, st, note in events_73834:
            con.execute("""INSERT INTO candidates(id, observation_id, run_id, tic, sector, source_file,
                center_btjd, score, width_hr, depth_ppm, n_support_windows, raw_outcome, raw_reason,
                scientific_classification, scientific_description, true_depth_min_pct, true_depth_mean_pct,
                consecutive_depressed, has_shoulder_flare, shoulder_elevation_pct, review_status, note,
                system_id, created_utc)
                VALUES (?, 2, ?, '73834', 8, 'lc_73834_s08.npz', ?, ?, 2.0, ?, 2,
                'verified_localized_dimming', 'pass', 'verified_localized_dimming',
                'Deep localized dimming excursion', 0.75, 0.72, 10, 0, 0.0, ?, ?, NULL, '2026-10-15T12:02:40Z')""",
                (cid, run_id, center_t, score, depth, st, note))
            con.execute("""INSERT INTO candidate_detections(candidate_id, run_id, window_start_btjd,
                window_end_btjd, center_btjd, score, width_hr, depth_ppm, created_utc)
                VALUES (?, ?, ?, ?, ?, ?, 2.0, ?, '2026-10-15T12:02:40Z')""",
                (cid, run_id, center_t - 0.25, center_t + 0.25, center_t, score, depth))

        # 5. Observation 3: TIC 60232 (Sector 5) - Continuous Stellar Pulsator
        con.execute("""INSERT INTO observations(id, run_id, tic, sector, camera, ccd, tmag, teff, obs_id,
            product_filename, config_hash, status, n_cadences, n_usable_cadences, n_windows, n_alert_windows,
            n_events, n_verified_events, created_utc)
            VALUES (3, ?, '60232', 5, 1, 1, 9.80, 6450.0, 'tess2018319111500-s0005-0000000000060232-0125-s',
            'tess2018319111500-s0005-0000000000060232-0125-s_lc.fits', 'cfg_hash', 'completed',
            19500, 18200, 74, 3, 3, 0, '2026-10-15T12:03:00Z')""", (run_id,))

        con.execute("""INSERT INTO candidates(id, observation_id, run_id, tic, sector, source_file,
            center_btjd, score, width_hr, depth_ppm, n_support_windows, raw_outcome, raw_reason,
            scientific_classification, scientific_description, true_depth_min_pct, true_depth_mean_pct,
            consecutive_depressed, has_shoulder_flare, shoulder_elevation_pct, review_status, note,
            system_id, created_utc)
            VALUES (7, 3, ?, '60232', 5, 'lc_60232_s05.npz', 1442.5000, 78.50, 2.0, 3200.0, 1,
            'continuous_stellar_variability', 'continuous_sinusoidal_modulation', 'continuous_stellar_variability',
            'Pulsational trough pattern detected', 0.32, 0.30, 8, 0, 0.0, 'unreviewed', '', NULL, '2026-10-15T12:03:10Z')""",
            (run_id,))

        # 6. Observation 4: TIC 141944 (Sector 4) - Flare Phantom
        con.execute("""INSERT INTO observations(id, run_id, tic, sector, camera, ccd, tmag, teff, obs_id,
            product_filename, config_hash, status, n_cadences, n_usable_cadences, n_windows, n_alert_windows,
            n_events, n_verified_events, created_utc)
            VALUES (4, ?, '141944', 4, 4, 3, 12.05, 3800.0, 'tess2018292100000-s0004-0000000000141944-0124-s',
            'tess2018292100000-s0004-0000000000141944-0124-s_lc.fits', 'cfg_hash', 'completed',
            19600, 18500, 75, 1, 1, 0, '2026-10-15T12:03:30Z')""", (run_id,))

        con.execute("""INSERT INTO candidates(id, observation_id, run_id, tic, sector, source_file,
            center_btjd, score, width_hr, depth_ppm, n_support_windows, raw_outcome, raw_reason,
            scientific_classification, scientific_description, true_depth_min_pct, true_depth_mean_pct,
            consecutive_depressed, has_shoulder_flare, shoulder_elevation_pct, review_status, note,
            system_id, created_utc)
            VALUES (8, 4, ?, '141944', 4, 'lc_141944_s04.npz', 1418.8000, 76.12, 2.0, 5400.0, 1,
            'flare_induced_apparent_dip', 'prominent_adjacent_flare_shoulder', 'flare_induced_apparent_dip',
            'Spurious sliding-box filter artifact created by flare shoulder overshoot', 0.0, 0.0, 0, 1, 8.4,
            'dismissed', 'Confirmed synthetic flare phantom rejected by verifier.', NULL, '2026-10-15T12:03:40Z')""",
            (run_id,))

        # 7. Observation 5: TIC 999001 (Sector 1) - Quiet Control Star
        con.execute("""INSERT INTO observations(id, run_id, tic, sector, camera, ccd, tmag, teff, obs_id,
            product_filename, config_hash, status, n_cadences, n_usable_cadences, n_windows, n_alert_windows,
            n_events, n_verified_events, created_utc)
            VALUES (5, ?, '999001', 1, 1, 2, 10.50, 5700.0, 'tess2018206190142-s0001-0000000000999001-0120-s',
            'tess2018206190142-s0001-0000000000999001-0120-s_lc.fits', 'cfg_hash', 'no_anomalies',
            19850, 19200, 77, 0, 0, 0, '2026-10-15T12:04:00Z')""", (run_id,))

    con.close()
    print(f"Created sample survey SQLite database at: {SAMPLE_DB_PATH}")

    # Generate synthetic local NPZ files in sample_data for offline zero-internet demonstration
    np.random.seed(42)

    # 1. Observation 1: TIC 39018208 (Sector 2) - Illustrative Synthetic Eclipsing Binary Benchmark
    # Full Sector 2 baseline: ~27.4 days (19,737 cadences)
    n_cad2 = 19737
    t_sec2 = 1354.10 + np.arange(n_cad2) * (2.0 / 1440.0)
    downlink_sec2 = (t_sec2 >= 1367.6) & (t_sec2 <= 1368.8)
    good_sec2 = ~downlink_sec2
    f_sec2 = np.random.normal(0, 0.00035, n_cad2)
    f_sec2[downlink_sec2] = 0.0

    # Inject periodic eclipses matching candidate systems (P = 2.41824 d, T0 = 1356.4180, depth = 0.482%, width = 2h)
    p_eblm = 2.41824
    t0_eblm = 1356.4180
    for k in range(-1, 12):
        ec = t0_eblm + k * p_eblm
        d_mask = (np.abs(t_sec2 - ec) <= (1.0 / 24.0)) & good_sec2
        f_sec2[d_mask] -= 0.00482

    np.savez_compressed(
        SAMPLE_DATA_DIR / "lc_39018208_s02.npz",
        time_btjd=t_sec2,
        pdc_flux_rel=f_sec2,
        pdc_flux_err_rel=np.full(n_cad2, 0.00035),
        good_pdc=good_sec2,
        quality=np.where(good_sec2, 0, 128)
    )

    # 2. Observation 2: TIC 73834 (Sector 8) - Deep Localized Dimming Candidate
    # Full Sector 8 baseline: ~25.9 days (18,600 cadences)
    n_cad8 = 18600
    t_sec8 = 1517.30 + np.arange(n_cad8) * (2.0 / 1440.0)
    downlink_sec8 = (t_sec8 >= 1529.2) & (t_sec8 <= 1530.4)
    good_sec8 = ~downlink_sec8
    f_sec8 = np.random.normal(0, 0.00045, n_cad8)
    f_sec8[downlink_sec8] = 0.0

    # Events at 1518.2500 (Orbit 1) and 1532.1000 (Orbit 2)
    for ec, depth in [(1518.2500, 0.00750), (1532.1000, 0.00720)]:
        d_mask = (np.abs(t_sec8 - ec) <= (1.0 / 24.0)) & good_sec8
        f_sec8[d_mask] -= depth

    np.savez_compressed(
        SAMPLE_DATA_DIR / "lc_73834_s08.npz",
        time_btjd=t_sec8,
        pdc_flux_rel=f_sec8,
        pdc_flux_err_rel=np.full(n_cad8, 0.00045),
        good_pdc=good_sec8,
        quality=np.where(good_sec8, 0, 128)
    )

    # 3. Observation 3: TIC 60232 (Sector 5) - Continuous Stellar Pulsator
    n_cad5 = 18500
    t_sec5 = 1438.00 + np.arange(n_cad5) * (2.0 / 1440.0)
    downlink_sec5 = (t_sec5 >= 1450.0) & (t_sec5 <= 1451.2)
    good_sec5 = ~downlink_sec5
    f_sec5 = np.random.normal(0, 0.00040, n_cad5)
    f_sec5[downlink_sec5] = 0.0
    # Sinusoidal pulsational modulation: P = 5.78 h (0.24083 d), A = 3200 ppm
    f_sec5[good_sec5] += 0.0032 * np.sin(2.0 * np.pi * (t_sec5[good_sec5] - 1438.0) / 0.24083)

    np.savez_compressed(
        SAMPLE_DATA_DIR / "lc_60232_s05.npz",
        time_btjd=t_sec5,
        pdc_flux_rel=f_sec5,
        pdc_flux_err_rel=np.full(n_cad5, 0.00040),
        good_pdc=good_sec5,
        quality=np.where(good_sec5, 0, 128)
    )

    # 4. Observation 4: TIC 141944 (Sector 4) - Flare Phantom Artifact
    n_cad4 = 18500
    t_sec4 = 1411.00 + np.arange(n_cad4) * (2.0 / 1440.0)
    downlink_sec4 = (t_sec4 >= 1423.0) & (t_sec4 <= 1424.2)
    good_sec4 = ~downlink_sec4
    f_sec4 = np.random.normal(0, 0.00035, n_cad4)
    f_sec4[downlink_sec4] = 0.0
    flare_time = 1418.80
    dt_flare = t_sec4 - flare_time
    flare_mask = (dt_flare >= 0) & (dt_flare <= 0.3) & good_sec4
    f_sec4[flare_mask] += 0.025 * np.exp(-dt_flare[flare_mask] / 0.04)
    dip_mask = (dt_flare > 0.3) & (dt_flare <= 0.45) & good_sec4
    f_sec4[dip_mask] -= 0.0054

    np.savez_compressed(
        SAMPLE_DATA_DIR / "lc_141944_s04.npz",
        time_btjd=t_sec4,
        pdc_flux_rel=f_sec4,
        pdc_flux_err_rel=np.full(n_cad4, 0.00035),
        good_pdc=good_sec4,
        quality=np.where(good_sec4, 0, 128)
    )

    # 5. Observation 5: TIC 999001 (Sector 1) - Quiet Control Star
    n_cad1 = 19000
    t_sec1 = 1325.00 + np.arange(n_cad1) * (2.0 / 1440.0)
    downlink_sec1 = (t_sec1 >= 1338.0) & (t_sec1 <= 1339.2)
    good_sec1 = ~downlink_sec1
    f_sec1 = np.random.normal(0, 0.00030, n_cad1)
    f_sec1[downlink_sec1] = 0.0

    np.savez_compressed(
        SAMPLE_DATA_DIR / "lc_999001_s01.npz",
        time_btjd=t_sec1,
        pdc_flux_rel=f_sec1,
        pdc_flux_err_rel=np.full(n_cad1, 0.00030),
        good_pdc=good_sec1,
        quality=np.where(good_sec1, 0, 128)
    )

    print(f"Generated sample light curve files in: {SAMPLE_DATA_DIR}")

if __name__ == "__main__":
    create_sample_database()
