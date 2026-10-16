"""Periodicity detection, harmonic analysis, and ephemeris fitting for astronomical events.

Uses epoch-difference clustering, dual-tolerance linear ephemeris fitting,
competing period preservation, and formal least-squares uncertainty calculation.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import math
import numpy as np


@dataclass
class PeriodCandidate:
    period_days: float
    epoch_btjd: float
    n_events: int
    is_periodic: bool
    odd_even_sigma: float
    max_residual_days: float
    period_snr: float
    period_uncertainty_days: float = 0.0
    epoch_uncertainty_days: float = 0.0
    alternate_period_days: float = 0.0
    classification: str = "no_periodicity_detected"  # 'periodic', 'no_periodicity_detected', 'confirmed_nonperiodic', 'single_event'
    competing_solutions: list = field(default_factory=list)
    note: str = ""


def cluster_candidate_events(events: list[dict], merge_radius_days: float = None) -> list[dict]:
    """Merge intra-dip duplicate candidate detections into unique physical dips,
    preserving support window provenance and selecting the peak SNR detection.
    When merge_radius_days is None, uses adaptive duration-matched clustering:
    rad = max(3 * cadence, 0.55 * max_width), which prevents artificial merging
    of closely spaced distinct eclipses while safely combining intra-dip duplicates.
    """
    if not events:
        return []
    sorted_ev = sorted(events, key=lambda e: float(e['center_btjd']))
    cadence_days = 2.0 / 1440.0

    clusters = []
    for e in sorted_ev:
        t = float(e['center_btjd'])
        if not clusters:
            clusters.append([e])
            continue

        prev = clusters[-1][-1]
        dt = t - float(prev['center_btjd'])

        if merge_radius_days is not None:
            rad = merge_radius_days
        else:
            w_d = float(e.get('width_hr', 2.0)) / 24.0
            prev_w_d = float(prev.get('width_hr', 2.0)) / 24.0
            rad = max(3.0 * cadence_days, 0.55 * max(w_d, prev_w_d))

        if dt <= rad:
            clusters[-1].append(e)
        else:
            clusters.append([e])

    unique_dips = []
    for c in clusters:
        best = max(c, key=lambda x: float(x.get('score', 0.0)))
        dip = dict(best)
        dip['n_support_windows'] = sum(int(x.get('n_support_windows', 1)) for x in c)
        dip['duplicate_count'] = len(c)
        unique_dips.append(dip)
    return unique_dips


def find_periodicity(events: list[dict], cadence_days: float = 2.0 / 1440.0,
                     min_events: int = 3, max_period_days: float = 20.0,
                     min_period_days: float = 0.2) -> PeriodCandidate:
    """Analyze a collection of candidate dimming events on a star for periodicity.

    Args:
        events: list of dicts with at least 'center_btjd', 'score', 'width_hr', 'depth_ppm'
        cadence_days: measurement cadence in days (default: 2 min = ~0.001389 d)
        min_events: minimum events required for periodic classification (strictly >= 3)
        max_period_days: upper bound on search period
        min_period_days: lower bound on search period (near stellar Roche limit)

    Returns:
        PeriodCandidate dataclass with detection metrics, uncertainties, and competing periods.
    """
    if not events:
        return PeriodCandidate(
            period_days=0.0, epoch_btjd=0.0, n_events=0, is_periodic=False,
            odd_even_sigma=0.0, max_residual_days=0.0, period_snr=0.0,
            classification="no_periodicity_detected", note="No candidate events found"
        )

    # 1. Cluster intra-dip duplicate candidate rows
    unique_dips = cluster_candidate_events(events)
    n_unique = len(unique_dips)

    if n_unique < min_events:
        class_label = "single_event" if n_unique == 1 else "no_periodicity_detected"
        note_str = (f"Single isolated dip detected ({len(events)} raw detections; requires >= {min_events} distinct dips)"
                    if n_unique == 1 else
                    f"Insufficient distinct dips ({n_unique} unique dips from {len(events)} detections; requires >= {min_events})")
        return PeriodCandidate(
            period_days=0.0,
            epoch_btjd=float(unique_dips[0]['center_btjd']),
            n_events=n_unique,
            is_periodic=False,
            odd_even_sigma=0.0,
            max_residual_days=0.0,
            period_snr=float(max(float(e.get('score', 0.0)) for e in unique_dips)),
            classification=class_label,
            note=note_str
        )

    # Sort unique events chronologically
    sorted_events = sorted(unique_dips, key=lambda e: float(e['center_btjd']))
    times = np.array([float(e['center_btjd']) for e in sorted_events], dtype=float)
    scores = np.array([float(e.get('score', 1.0)) for e in sorted_events], dtype=float)
    depths = np.array([float(e.get('depth_ppm', 0.0)) for e in sorted_events], dtype=float)
    widths = np.array([float(e.get('width_hr', 2.0)) for e in sorted_events], dtype=float)

    mean_duration_days = float(np.mean(widths)) / 24.0
    # Dual timing criteria: RMS tolerance <= 12% duration, Max tolerance <= 25% duration
    tol_rms_days = max(2.0 * cadence_days, 0.12 * mean_duration_days)
    tol_max_days = max(3.0 * cadence_days, 0.25 * mean_duration_days)

    # 2. Pairwise time differences
    diffs = []
    n = len(times)
    for i in range(n):
        for j in range(i + 1, n):
            d = times[j] - times[i]
            if d >= min_period_days:
                diffs.append(d)

    if not diffs:
        return PeriodCandidate(
            period_days=0.0, epoch_btjd=float(times[0]), n_events=n,
            is_periodic=False, odd_even_sigma=0.0, max_residual_days=0.0,
            period_snr=0.0, classification="no_periodicity_detected",
            note="No valid pairwise differences above min_period"
        )

    # 3. Candidate trial periods from subharmonic divisors of differences
    trial_periods = set()
    for d in diffs:
        for m in range(1, 8):
            p = d / m
            if min_period_days <= p <= max_period_days:
                trial_periods.add(round(p, 3))

    sorted_trials = sorted(trial_periods)
    thinned_trials = []
    last_p = -1.0
    for p in sorted_trials:
        if p - last_p >= 0.005:
            thinned_trials.append(p)
            last_p = p

    # Dense grid search around trial periods
    refined_periods = []
    for p in thinned_trials:
        refined_periods.extend(np.linspace(p * 0.985, p * 1.015, 11))

    if not refined_periods:
        refined_periods = np.linspace(min_period_days, max_period_days, 1000)

    # 4. Fit linear ephemeris t_i ~ T_0 + n_i * P and collect valid candidate periods
    t0_candidate = float(times[0])
    valid_candidates = []
    min_observed_max_res = float('inf')

    for p in refined_periods:
        if p < min_period_days or p > max_period_days:
            continue
        exact_n = (times - t0_candidate) / p
        int_n = np.round(exact_n).astype(int)

        # Quick vectorized check
        deviations = np.abs(exact_n - int_n) * p
        max_dev = float(np.max(deviations))
        if max_dev < min_observed_max_res:
            min_observed_max_res = max_dev

        if max_dev > tol_max_days:
            continue

        if len(set(int_n)) != len(times):
            continue

        steps = [int(s) for s in np.diff(int_n)]
        g = math.gcd(*steps) if steps else 1
        if g > 1:
            p_fund = p * g
            int_n = int_n // g
        else:
            p_fund = p

        if p_fund < min_period_days or p_fund > max_period_days:
            continue

        residuals = np.abs(times - (t0_candidate + int_n * p_fund))
        max_res = float(np.max(residuals))
        rms_res = float(np.sqrt(np.mean(residuals ** 2)))

        if max_res <= tol_max_days and rms_res <= tol_rms_days:
            span = int(int_n[-1] - int_n[0])
            completeness = len(times) / (span + 1)
            valid_candidates.append({
                'period': p_fund,
                'epoch': t0_candidate,
                'epochs': int_n,
                'max_res': max_res,
                'rms_res': rms_res,
                'completeness': completeness
            })

    if not valid_candidates:
        # Distinguish between confirmed_nonperiodic and no_periodicity_detected
        is_confirmed_nonperiodic = (min_observed_max_res > 0.5 * mean_duration_days) and (n >= 4)
        class_label = "confirmed_nonperiodic" if is_confirmed_nonperiodic else "no_periodicity_detected"
        note_str = ("Confirmed non-periodic event distribution (residuals exceed physical duration limits >60 min)"
                    if is_confirmed_nonperiodic else
                    "No periodic ephemeris satisfies timing tolerance")
        return PeriodCandidate(
            period_days=0.0,
            epoch_btjd=float(times[0]),
            n_events=n,
            is_periodic=False,
            odd_even_sigma=0.0,
            max_residual_days=min_observed_max_res if min_observed_max_res != float('inf') else 0.0,
            period_snr=0.0,
            classification=class_label,
            note=note_str
        )

    # Sort candidates preferring higher completeness (fundamental period) then smaller residuals
    valid_candidates.sort(key=lambda c: (-c['completeness'], c['rms_res']))
    best = valid_candidates[0]
    best_period = best['period']
    best_epochs_int = best['epochs']

    # 5. Refine ephemeris via least-squares linear fit
    A = np.vstack([np.ones_like(best_epochs_int), best_epochs_int]).T
    fitted_t0, fitted_p = np.linalg.lstsq(A, times, rcond=None)[0]
    final_res = np.abs(times - (fitted_t0 + best_epochs_int * fitted_p))
    final_max_res = float(np.max(final_res))
    final_rms_res = float(np.sqrt(np.mean(final_res ** 2)))

    if final_max_res > tol_max_days or final_rms_res > tol_rms_days:
        return PeriodCandidate(
            period_days=float(best_period),
            epoch_btjd=float(best['epoch']),
            n_events=n,
            is_periodic=False,
            odd_even_sigma=0.0,
            max_residual_days=final_max_res,
            period_snr=0.0,
            classification="no_periodicity_detected",
            note="Linear ephemeris residuals exceed tolerance after least-squares refinement"
        )

    # Formal least-squares uncertainties
    denom = float(np.sum((best_epochs_int - np.mean(best_epochs_int)) ** 2))
    sigma_p = float(final_rms_res / math.sqrt(denom)) if denom > 0 else 0.0
    sigma_t0 = float(final_rms_res * math.sqrt(np.sum(best_epochs_int ** 2) / (n * denom))) if denom > 0 else 0.0

    # 6. Check odd/even eclipse disparity (eclipsing binary check)
    odd_mask = (best_epochs_int % 2 == 1)
    even_mask = (best_epochs_int % 2 == 0)

    odd_even_sigma = 0.0
    alternate_period = 0.0
    if np.any(odd_mask) and np.any(even_mask):
        odd_depths = depths[odd_mask]
        even_depths = depths[even_mask]
        mean_odd = float(np.mean(odd_depths))
        mean_even = float(np.mean(even_depths))
        var_odd = float(np.var(odd_depths, ddof=1)) if len(odd_depths) > 1 else (0.1 * mean_odd) ** 2
        var_even = float(np.var(even_depths, ddof=1)) if len(even_depths) > 1 else (0.1 * mean_even) ** 2
        denom_var = math.sqrt(max(1.0, var_odd / len(odd_depths) + var_even / len(even_depths)))
        odd_even_sigma = abs(mean_odd - mean_even) / denom_var

        if odd_even_sigma > 3.0:
            # Preserving competing solutions without destructive overwriting:
            # Primary period is the dip recurrence interval; alternate is the full binary period.
            alternate_period = float(fitted_p * 2.0)
        else:
            alternate_period = float(fitted_p / 2.0)

    # 7. Compile competing solutions
    competing = []
    seen_periods = set()
    for sol in valid_candidates:
        round_p = round(sol['period'], 3)
        if round_p not in seen_periods:
            seen_periods.add(round_p)
            competing.append({
                'period_days': float(sol['period']),
                'rms_res_days': float(sol['rms_res']),
                'max_res_days': float(sol['max_res']),
                'completeness': float(sol['completeness'])
            })
            if len(competing) >= 3:
                break

    # 8. Check momentum dump harmonic
    note = "Candidate periodic system"
    if odd_even_sigma > 3.0:
        note += f" (Eclipsing binary candidate with alternate P_binary={alternate_period:.4f}d, Z={odd_even_sigma:.1f} sigma)"
    if (abs(fitted_p - 2.50) < 0.06 or abs(fitted_p - 1.25) < 0.04 or
        (alternate_period > 0 and abs(alternate_period - 2.50) < 0.06)):
        note += " (Warning: near TESS 2.5d momentum dump cycle)"

    period_snr = float(np.sqrt(np.sum(scores ** 2)))

    return PeriodCandidate(
        period_days=float(fitted_p),
        period_uncertainty_days=sigma_p,
        epoch_btjd=float(fitted_t0),
        epoch_uncertainty_days=sigma_t0,
        n_events=n,
        is_periodic=True,
        odd_even_sigma=float(odd_even_sigma),
        max_residual_days=float(final_max_res),
        period_snr=period_snr,
        alternate_period_days=float(alternate_period),
        classification="periodic",
        competing_solutions=competing,
        note=note
    )


def compute_phase(time_btjd: np.ndarray, epoch_btjd: float, period_days: float) -> np.ndarray:
    """Compute phase in range [-0.5, 0.5] centered on transit epoch."""
    if period_days <= 0:
        return np.zeros_like(time_btjd)
    phase = ((time_btjd - epoch_btjd) / period_days) % 1.0
    return np.where(phase > 0.5, phase - 1.0, phase)
