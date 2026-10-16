"""Synthetic signal injection and recovery benchmarking for TESS Anomaly Survey.

Injects synthetic transit and asymmetric dimming signals into quality-controlled TESS
observations, evaluating recovery fraction, timing accuracy, matched unmodified controls,
and false-alarm incidence under calibrated survey thresholds.
"""
from __future__ import annotations
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import numpy as np

from . import core
import experiment_b1 as b1


def wilson_score_interval(successes: int, trials: int, confidence: float = 0.95) -> tuple[float, float]:
    """Compute Wilson score confidence interval for a binomial proportion."""
    if trials == 0:
        return 0.0, 0.0
    p = successes / trials
    z = 1.95996  # 95% confidence standard normal quantile
    denom = 1.0 + z**2 / trials
    centre = (p + z**2 / (2.0 * trials)) / denom
    margin = (z * math.sqrt((p * (1.0 - p) + z**2 / (4.0 * trials)) / trials)) / denom
    lower = max(0.0, centre - margin)
    upper = min(1.0, centre + margin)
    return lower, upper


def generate_signal_profile(
    times: np.ndarray,
    t0: float,
    duration_hr: float,
    depth_ppm: float,
    shape: str = "box"
) -> np.ndarray:
    """Generate fractional flux diminution profile (0.0 = no dimming, delta = maximum dip).

    Args:
        times: array of observation timestamps in BTJD
        t0: injection center / onset epoch in BTJD
        duration_hr: event duration in hours
        depth_ppm: event depth in parts per million
        shape: 'box' (flat-bottom transit) or 'asymmetric' (steep ingress, exponential egress)

    Returns:
        Array of fractional depth diminution values (positive dimming).
    """
    duration_days = duration_hr / 24.0
    delta = depth_ppm * 1e-6
    dimming = np.zeros_like(times, dtype=float)

    if shape == "box":
        dt = np.abs(times - t0)
        # Trapezoid with 15% ingress/egress ramps for realistic smoothing
        half_w = duration_days / 2.0
        ramp = 0.15 * duration_days
        in_flat = dt <= (half_w - ramp)
        in_ramp = (dt > (half_w - ramp)) & (dt <= half_w)
        dimming[in_flat] = delta
        if np.any(in_ramp) and ramp > 0:
            dimming[in_ramp] = delta * (half_w - dt[in_ramp]) / ramp
    elif shape == "asymmetric":
        # Cometary / dust-clump proxy: instantaneous ingress, exponential recovery
        # Scale decay timescale so that 90% of recovery occurs within duration_days
        tau = max(1e-4, duration_days / 2.5)
        dt = times - t0
        egress = dt >= 0.0
        dimming[egress] = delta * np.exp(-dt[egress] / tau)
        # 10% ramp before t0 for continuity
        ingress_ramp = (dt < 0.0) & (dt >= -0.1 * duration_days)
        if np.any(ingress_ramp):
            dimming[ingress_ramp] = delta * (1.0 + dt[ingress_ramp] / (0.1 * duration_days))
    else:
        raise ValueError(f"Unknown signal shape: {shape}")

    return dimming


def run_benchmark(
    catalog_path: str | Path,
    root: str | Path,
    calibration_path: str | Path,
    depths_ppm: list[float] | None = None,
    durations_hr: list[float] | None = None,
    shapes: list[str] | None = None,
    n_injections_per_cell: int = 20,
    seed: int = 2026,
    output_path: str | Path | None = None,
    progress_callback: callable | None = None
) -> dict:
    """Run systematic injection-recovery benchmark across depths, durations, and shapes.

    Args:
        catalog_path: path to validation or benchmark light curve catalog CSV
        root: root directory containing processed TESS NPZs
        calibration_path: path to calibration JSON file
        depths_ppm: list of injection depths in ppm (default: [200, 500, 1000, 2000, 5000, 10000])
        durations_hr: list of durations in hours (default: [0.5, 1.0, 2.0, 4.0, 8.0])
        shapes: list of signal families (default: ['box', 'asymmetric'])
        n_injections_per_cell: number of injection trials per grid cell
        seed: random seed for reproducibility
        output_path: optional JSON output path
        progress_callback: optional progress reporting callback

    Returns:
        Structured dictionary of benchmark results and stratified summaries.
    """
    cal = core.read_json(calibration_path)
    cal_method = cal['method']
    cal_settings = cal['settings']

    if depths_ppm is None:
        depths_ppm = [200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0]
    if durations_hr is None:
        durations_hr = [0.5, 1.0, 2.0, 4.0, 8.0]
    if shapes is None:
        shapes = ["box", "asymmetric"]

    # Read eligible light curves
    rows = list(core.eligible_catalog_rows(catalog_path, root))
    if not rows:
        raise RuntimeError("No eligible catalog rows found for benchmark")

    rng = np.random.default_rng(seed)

    # Pre-select clean observation files
    clean_files = []
    for r in rows:
        p = Path(root) / r['processed_file']
        clean_files.append((r, p))

    # Pre-harvest valid windows to inject into
    sampled_windows = []
    cache = {}
    for r, p in clean_files:
        try:
            for w in core.iter_windows(p, **cal_settings):
                # Verify that window has high good quality fraction
                good_frac = float(np.mean(w['good']))
                if good_frac >= 0.98:
                    sampled_windows.append((r, p, w))
                if len(sampled_windows) >= max(300, n_injections_per_cell * 15):
                    break
        except Exception:
            continue
        if len(sampled_windows) >= max(300, n_injections_per_cell * 15):
            break

    if not sampled_windows:
        raise RuntimeError("No valid windows found with required quality criteria")

    grid_results = []
    total_cells = len(shapes) * len(durations_hr) * len(depths_ppm)
    cell_idx = 0

    for shape in shapes:
        for duration_hr in durations_hr:
            for depth_ppm in depths_ppm:
                cell_idx += 1
                n_trials = n_injections_per_cell
                recovered_count = 0
                unrelated_alert_count = 0
                control_alert_count = 0
                timing_errors_min = []
                recovered_scores = []
                control_scores = []

                # Select n_trials windows with replacement
                window_indices = rng.choice(len(sampled_windows), size=n_trials, replace=True)

                for trial_i, win_idx in enumerate(window_indices):
                    r, p, w = sampled_windows[win_idx]
                    t = w['time']
                    f = w['flux']
                    good = w['good']
                    cadence = w['cadence_days']

                    # 1. Evaluate matched unmodified control window
                    _, noise_ctrl, top_score_ctrl = core.score_window(w, cal_method, cache)
                    thresh_ctrl, _ = core.threshold_for(cal, noise_ctrl)
                    is_ctrl_alert = bool(top_score_ctrl > thresh_ctrl)
                    control_scores.append(float(top_score_ctrl))
                    if is_ctrl_alert:
                        control_alert_count += 1

                    # 2. Synthetic Injection
                    # Place injection center in the central 50% of the window to ensure full duration fits
                    t_span = t[-1] - t[0]
                    duration_days = duration_hr / 24.0
                    t_min_inj = t[0] + 0.25 * t_span
                    t_max_inj = t[-1] - 0.25 * t_span

                    if t_max_inj <= t_min_inj:
                        t0_inj = float(t[len(t) // 2])
                    else:
                        t0_inj = float(rng.uniform(t_min_inj, t_max_inj))

                    dimming = generate_signal_profile(t, t0_inj, duration_hr, depth_ppm, shape=shape)
                    # In relative flux (where baseline is ~0.0), a dip of depth delta is f - dimming
                    f_inj = f - dimming

                    # Create injected window dict
                    w_inj = dict(w, flux=f_inj)

                    # Exact survey scoring on injected series
                    raw_inj, noise_inj, top_score_inj = core.score_window(w_inj, cal_method, cache)
                    thresh_inj, _ = core.threshold_for(cal, noise_inj)

                    # Run exact survey deduplication
                    candidates = b1.deduplicate(raw_inj, cadence, thresh_inj)

                    # Evaluate correct-time recovery vs unrelated alert
                    # Candidate match window: max(1.5 * cadence, 0.5 * duration)
                    match_tolerance_days = max(1.5 * cadence, 0.5 * duration_days)
                    hit = False
                    unrelated = False
                    best_cand_score = 0.0
                    best_time_err = 0.0

                    for cand in candidates:
                        c_time = float(cand['center_btjd'])
                        c_score = float(cand['score'])
                        dt_days = abs(c_time - t0_inj)

                        if dt_days <= match_tolerance_days:
                            hit = True
                            best_cand_score = max(best_cand_score, c_score)
                            best_time_err = dt_days * 1440.0  # in minutes
                        else:
                            unrelated = True

                    if hit:
                        recovered_count += 1
                        recovered_scores.append(best_cand_score)
                        timing_errors_min.append(best_time_err)
                    if unrelated:
                        unrelated_alert_count += 1

                # Cell summary statistics
                rec_rate = recovered_count / n_trials
                rec_ci_low, rec_ci_high = wilson_score_interval(recovered_count, n_trials)
                ctrl_rate = control_alert_count / n_trials
                ctrl_ci_low, ctrl_ci_high = wilson_score_interval(control_alert_count, n_trials)
                unrel_rate = unrelated_alert_count / n_trials

                cell_summary = {
                    'shape': shape,
                    'duration_hr': duration_hr,
                    'depth_ppm': depth_ppm,
                    'trials': n_trials,
                    'recovered': recovered_count,
                    'recovery_rate': round(rec_rate, 4),
                    'recovery_ci_95': [round(rec_ci_low, 4), round(rec_ci_high, 4)],
                    'control_alerts': control_alert_count,
                    'control_false_alarm_rate': round(ctrl_rate, 4),
                    'control_ci_95': [round(ctrl_ci_low, 4), round(ctrl_ci_high, 4)],
                    'unrelated_alerts': unrelated_alert_count,
                    'unrelated_alert_rate': round(unrel_rate, 4),
                    'mean_recovered_score': round(float(np.mean(recovered_scores)), 2) if recovered_scores else 0.0,
                    'mean_control_score': round(float(np.mean(control_scores)), 2) if control_scores else 0.0,
                    'mean_timing_error_min': round(float(np.mean(timing_errors_min)), 2) if timing_errors_min else 0.0,
                }
                grid_results.append(cell_summary)

                if progress_callback:
                    progress_callback(cell_idx, total_cells, cell_summary)

    benchmark_manifest = {
        'format_version': 'SURVEY_BENCHMARK_V1',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'calibration_path': str(Path(calibration_path).resolve()),
        'catalog_path': str(Path(catalog_path).resolve()),
        'method': cal_method,
        'target_alert_fraction': cal.get('target_alert_fraction'),
        'seed': seed,
        'n_injections_per_cell': n_injections_per_cell,
        'depths_ppm': depths_ppm,
        'durations_hr': durations_hr,
        'shapes': shapes,
        'total_evaluations': len(grid_results) * n_injections_per_cell,
        'grid_results': grid_results,
        'notes': 'Measured correct-time recovery fraction and matched control false alarms across synthetic signal parameter grid under exact survey scoring and windowing.'
    }

    if output_path:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        core.write_json(out_p, benchmark_manifest)
        # Also write summary CSV alongside JSON
        csv_p = out_p.with_suffix('.csv')
        export_benchmark_csv(grid_results, csv_p)

    return benchmark_manifest


def export_benchmark_csv(grid_results: list[dict], output_csv: Path | str) -> None:
    """Write benchmark grid results to machine-readable CSV format."""
    cols = [
        'shape', 'duration_hr', 'depth_ppm', 'trials', 'recovered',
        'recovery_rate', 'recovery_ci_low', 'recovery_ci_high',
        'control_alerts', 'control_false_alarm_rate', 'control_ci_low', 'control_ci_high',
        'unrelated_alerts', 'unrelated_alert_rate', 'mean_recovered_score',
        'mean_control_score', 'mean_timing_error_min'
    ]
    with Path(output_csv).open('w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerow(cols)
        for r in grid_results:
            writer.writerow([
                r['shape'], r['duration_hr'], r['depth_ppm'], r['trials'], r['recovered'],
                r['recovery_rate'], r['recovery_ci_95'][0], r['recovery_ci_95'][1],
                r['control_alerts'], r['control_false_alarm_rate'], r['control_ci_95'][0], r['control_ci_95'][1],
                r['unrelated_alerts'], r['unrelated_alert_rate'], r['mean_recovered_score'],
                r['mean_control_score'], r['mean_timing_error_min']
            ])
