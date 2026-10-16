"""Astrophysical candidate vetting, instrumental systematics decontamination, and shortlist prioritization.

Provides automated tools to identify simultaneous cross-star spacecraft systematics
(e.g., momentum-wheel dumps, orbit crossing glitches), compute decontaminated event
significance, categorize candidate stars into scientific vetting tiers, and rank
anomalies for human inspection and follow-up.
"""
from __future__ import annotations
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
import numpy as np

from . import periodicity


def detect_simultaneous_clusters(
    events: List[Dict[str, Any]],
    time_tolerance_days: float = 0.04,
    min_co_occurring_stars: int = 3
) -> Tuple[List[Dict[str, Any]], Set[Tuple[str, float]]]:
    """Identify simultaneous cross-star event clusters representing instrumental systematics.

    TESS experiences periodic spacecraft momentum-wheel desaturations, reaction-wheel torque dumps,
    and pointing jitter that trigger simultaneous anomaly detections across multiple stars on the
    same camera/CCD.

    Args:
        events: List of event dicts with keys 'tic' and 'center_btjd' (and optionally 'score', 'sector').
        time_tolerance_days: Max timestamp difference to group events (default 0.04d ~ 58 min).
        min_co_occurring_stars: Minimum distinct TIC stars required to classify cluster as systematic.

    Returns:
        simultaneous_clusters: List of cluster dicts detailing mean_btjd, star_count, and tics.
        flagged_events: Set of (tic, center_btjd) tuples falling into simultaneous clusters.
    """
    if not events:
        return [], set()

    sorted_events = sorted(events, key=lambda x: float(x['center_btjd']))
    time_clusters: List[List[Dict[str, Any]]] = []
    curr: List[Dict[str, Any]] = [sorted_events[0]]

    for ev in sorted_events[1:]:
        t = float(ev['center_btjd'])
        t_prev = float(curr[-1]['center_btjd'])
        if t - t_prev <= time_tolerance_days:
            curr.append(ev)
        else:
            time_clusters.append(curr)
            curr = [ev]
    if curr:
        time_clusters.append(curr)

    simultaneous_clusters = []
    flagged_events = set()

    for cl in time_clusters:
        distinct_tics = list(set(e['tic'] for e in cl))
        if len(distinct_tics) >= min_co_occurring_stars:
            mean_t = float(np.mean([float(e['center_btjd']) for e in cl]))
            sector = cl[0].get('sector')
            simultaneous_clusters.append({
                'mean_btjd': mean_t,
                'sector': sector,
                'star_count': len(distinct_tics),
                'tics': distinct_tics,
                'n_events': len(cl)
            })
            for e in cl:
                flagged_events.add((e['tic'], float(e['center_btjd'])))

    return simultaneous_clusters, flagged_events


def compute_star_vetting_metrics(
    cands: List[Dict[str, Any]],
    flagged_systematics: Set[Tuple[str, float]],
    sys_info: Optional[Dict[str, Any]] = None,
    cat_match: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """Compute comprehensive astrophysical vetting and decontamination metrics for a single star.

    Args:
        cands: Raw candidate records for this star.
        flagged_systematics: Set of (tic, center_btjd) tuples flagged as cross-star systematics.
        sys_info: Periodicity analysis results from candidate_systems table (if available).
        cat_match: Existing catalog identification dict (if available).

    Returns:
        Structured dictionary with clean event counts, scores, depths, systematics fraction, and tier.
    """
    if not cands:
        raise ValueError("Cannot compute vetting metrics for empty candidate list")

    tic = cands[0]['tic']
    sectors = sorted(list(set(c['sector'] for c in cands if c.get('sector') is not None)))
    cands_sorted = sorted(cands, key=lambda x: float(x['center_btjd']))

    # Duration-matched adaptive clustering to find unique physical dips
    events_in = [
        {
            'center_btjd': float(c['center_btjd']),
            'score': float(c['score']),
            'width_hr': float(c.get('width_hr', 2.0)),
            'depth_ppm': float(c.get('depth_ppm', 0.0)),
            'id': c.get('id')
        }
        for c in cands_sorted
    ]
    unique_dips = periodicity.cluster_candidate_events(events_in)
    n_dips = len(unique_dips)

    # Classify clean vs systematic dips
    clean_dips = []
    systematic_dips = []
    for d in unique_dips:
        t0 = float(d['center_btjd'])
        if (tic, t0) in flagged_systematics:
            systematic_dips.append(d)
        else:
            clean_dips.append(d)

    n_clean = len(clean_dips)
    n_sys = len(systematic_dips)
    f_sys = n_sys / n_dips if n_dips > 0 else 0.0

    scores = [d['score'] for d in unique_dips]
    clean_scores = [d['score'] for d in clean_dips] if clean_dips else [0.0]
    depths = [d.get('depth_ppm', 0.0) for d in unique_dips]
    clean_depths = [d.get('depth_ppm', 0.0) for d in clean_dips] if clean_dips else [0.0]
    widths = [d.get('width_hr', 2.0) for d in unique_dips]

    max_score = float(max(scores)) if scores else 0.0
    max_clean_score = float(max(clean_scores)) if clean_scores else 0.0
    max_depth = float(max(depths)) if depths else 0.0
    med_depth = float(np.median(depths)) if depths else 0.0
    max_clean_depth = float(max(clean_depths)) if clean_depths else 0.0
    med_width = float(np.median(widths)) if widths else 2.0

    # Periodicity status
    is_periodic = bool(sys_info.get('is_periodic', False)) if sys_info else False
    p_days = float(sys_info.get('period_days', 0.0)) if sys_info else 0.0
    p_class = sys_info.get('classification', 'single_event' if n_dips == 1 else 'untested') if sys_info else ('single_event' if n_dips == 1 else 'untested')

    # Catalog info
    cat_match = cat_match or {}
    cat_id = cat_match.get('catalog_identifier', 'Uncataloged')
    astro_nature = cat_match.get('astrophysical_nature', 'Pending Vetting')
    match_cat = cat_match.get('match_category', 'None')

    # Tier assignment
    if is_periodic:
        if cat_id in ('Uncataloged', 'None', '') or 'Unverified' in match_cat:
            tier = "Tier 4: Unverified Periodic System"
        else:
            tier = "Tier 5: Confirmed Catalog Object"
    elif f_sys >= 0.70 and n_dips >= 4:
        tier = "Tier 5: Instrumental Systematics"
    elif max_clean_depth >= 100000.0 or max_clean_score >= 500.0:
        tier = "Tier 1: Extreme Non-Periodic Anomaly"
    elif n_clean >= 2:
        tier = "Tier 2: Clean Multi-Dip Non-Periodic"
    elif n_clean == 1:
        tier = "Tier 3: Clean Single-Dip Candidate"
    else:
        tier = "Tier 5: Low-Confidence / Systematic"

    # Priority score Q
    # Q = log10(max(max_clean_score, 1)) * sqrt(max(n_clean, 1)) * (1 - f_sys)
    q_score = 0.0
    if max_clean_score > 0.0 and tier not in ("Tier 5: Instrumental Systematics", "Tier 5: Confirmed Catalog Object"):
        q_score = math.log10(max(max_clean_score, 1.0)) * math.sqrt(max(n_clean, 1)) * (1.0 - 0.75 * f_sys)
        # Boost unverified periodic systems for targeted audit
        if tier == "Tier 4: Unverified Periodic System":
            q_score *= 1.25

    return {
        'tic': tic,
        'sectors': sectors,
        'n_cands': len(cands),
        'n_dips': n_dips,
        'n_clean_dips': n_clean,
        'n_sys_dips': n_sys,
        'systematic_fraction': round(f_sys, 3),
        'max_score': round(max_score, 1),
        'max_clean_score': round(max_clean_score, 1),
        'max_depth_ppm': round(max_depth, 0),
        'median_depth_ppm': round(med_depth, 0),
        'max_clean_depth_ppm': round(max_clean_depth, 0),
        'median_width_hr': round(med_width, 2),
        'is_periodic': is_periodic,
        'period_days': round(p_days, 5),
        'period_class': p_class,
        'catalog_identifier': cat_id,
        'astrophysical_nature': astro_nature,
        'match_category': match_cat,
        'tier': tier,
        'priority_score': round(q_score, 2),
        'clean_dips': clean_dips,
        'systematic_dips': systematic_dips
    }


def rank_vetted_shortlist(
    star_summaries: List[Dict[str, Any]],
    target_count: int = 24
) -> List[Dict[str, Any]]:
    """Rank vetted candidate stars and produce a balanced research shortlist.

    Balanced selection preserves diversity:
    - Both Tier 4 Unverified Periodic Candidates (TIC 294261093, TIC 4629065) for dedicated vetting
    - Top Tier 1 Extreme Non-Periodic Anomalies
    - Top Tier 2 Clean Multi-Dip Candidates
    - Top Tier 3 Clean Single-Dip Candidates
    """
    by_tier = defaultdict(list)
    for s in star_summaries:
        by_tier[s['tier']].append(s)

    # Sort within each tier by priority score descending
    for tier in by_tier:
        by_tier[tier].sort(key=lambda x: -x['priority_score'])

    shortlist = []
    # 1. Tier 1 Extreme Non-Periodic (top 8)
    shortlist.extend(by_tier.get("Tier 1: Extreme Non-Periodic Anomaly", [])[:8])

    # 2. Tier 2 Clean Multi-Dip (top 8)
    shortlist.extend(by_tier.get("Tier 2: Clean Multi-Dip Non-Periodic", [])[:8])

    # 3. Tier 3 Clean Single-Dip (top 6)
    shortlist.extend(by_tier.get("Tier 3: Clean Single-Dip Candidate", [])[:6])

    # Sort main discovery candidates by priority score descending
    shortlist.sort(key=lambda x: -x['priority_score'])

    # 4. Append Tier 4 Unverified Periodic Systems for explicit audit
    tier4 = by_tier.get("Tier 4: Unverified Periodic System", [])
    shortlist.extend(tier4)

    return shortlist[:target_count]


def verify_event_photometry_v2(time, flux, good, center_btjd, width_hr=2.0,
                               mf_depth_ppm=0.0, min_consecutive_cadences=2):
    """Version 2 coherence gate; legacy results remain reproducible with v1.

    Reuse v1 continuum/flare measurements but require consecutive original
    cadences AND timestamp adjacency. Masked points and gaps break a run.
    This is a post-detection heuristic, not a recalibrated detector or FAP.
    """
    result = verify_event_photometry(time, flux, good, center_btjd, width_hr,
                                     mf_depth_ppm, min_consecutive_cadences)
    result['verifier_version'] = 'coherence_v2'
    if 'robust_continuum' not in result:
        result['has_flare'] = None
        return result
    time, flux, good = np.asarray(time), np.asarray(flux), np.asarray(good, dtype=bool)
    differences = np.diff(time)
    positive = differences[np.isfinite(differences) & (differences > 0)]
    cadence = float(np.median(positive)) if positive.size else 0.0
    depressed = (good & np.isfinite(time) & np.isfinite(flux)
                 & (np.abs(time - center_btjd) <= width_hr / 48.0)
                 & (flux < result['robust_continuum'] - 1.5 * result['local_noise']))
    longest = current = 0
    for index, value in enumerate(depressed):
        adjacent = index > 0 and 0 < time[index] - time[index - 1] <= 1.5 * cadence
        current = (current + 1 if adjacent else 1) if value else 0
        longest = max(longest, current)
    result['consecutive_depressed'] = longest
    if result['is_verified'] and longest < min_consecutive_cadences:
        result.update(outcome='probable_phantom', is_verified=False,
                      reason='insufficient_timestamp_coherent_cadences_v2')
    return result


def verify_event_photometry(
    time: np.ndarray,
    flux: np.ndarray,
    good: np.ndarray,
    center_btjd: float,
    width_hr: float = 2.0,
    mf_depth_ppm: float = 0.0,
    min_consecutive_cadences: int = 2
) -> Dict[str, Any]:
    """Robust post-detection verification stage for candidate dimming events.

    Rejects phantom dips caused by flare-contaminated shoulders, positive outliers,
    or baseline elevation without confirmed negative flux excursions. Requires
    coherent negative flux across multiple consecutive cadences to reject isolated
    single-cadence downward spikes.

    Args:
        time: Array of BTJD timestamps.
        flux: Array of relative flux measurements (quiescent median ~ 0.0).
        good: Boolean validity mask (True for good cadences).
        center_btjd: Midpoint timestamp of candidate event.
        width_hr: Event duration in hours (default 2.0).
        mf_depth_ppm: Original matched-filter depth in parts per million.
        min_consecutive_cadences: Minimum consecutive depressed cadences required (default 2).

    Returns:
        Dict with 'outcome' ('verified_dimming' | 'probable_phantom' | 'inconclusive'),
        'is_verified' (bool), 'reason', 'true_depth_min_pct', 'true_depth_mean_pct',
        'consecutive_depressed', 'shoulder_elev_pct', 'has_flare', and baseline metrics.
    """
    t0 = float(center_btjd)
    half_w = (width_hr / 24.0) / 2.0

    # 1. Local slice +/- 12 hours
    m_local = (time >= t0 - 0.5) & (time <= t0 + 0.5) & good
    if np.sum(m_local) < 10:
        return {
            'outcome': 'inconclusive',
            'is_verified': False,
            'reason': 'insufficient_local_cadences',
            'true_depth_min_pct': 0.0,
            'true_depth_mean_pct': 0.0,
            'consecutive_depressed': 0,
            'shoulder_elev_pct': 0.0,
            'has_flare': False
        }

    t_loc = time[m_local]
    f_loc = flux[m_local]

    # 2. Robust Quiescent Continuum (excluding event core and positive flare outliers)
    m_out = abs(t_loc - t0) > (half_w * 2.0)
    if np.sum(m_out) >= 10:
        f_out = f_loc[m_out]
        med_out = float(np.median(f_out))
        mad_out = float(np.median(abs(f_out - med_out)))
        sig_out = 1.4826 * mad_out if mad_out > 0 else float(np.std(f_out))
        quiet_out = f_out[f_out < med_out + 3.0 * sig_out]
        robust_cont = float(np.median(quiet_out)) if len(quiet_out) > 5 else med_out
        noise = 1.4826 * float(np.median(abs(quiet_out - robust_cont))) if len(quiet_out) > 5 else sig_out
    else:
        robust_cont = float(np.median(f_loc))
        noise = float(np.std(f_loc))

    noise = max(noise, 1e-5)

    # 3. Core Box and Shoulder Measurements
    m_core = (t_loc >= t0 - half_w) & (t_loc <= t0 + half_w)
    m_left = (t_loc >= t0 - 2.0 * half_w) & (t_loc < t0 - half_w)
    m_right = (t_loc > t0 + half_w) & (t_loc <= t0 + 2.0 * half_w)

    n_core = int(np.sum(m_core))
    if n_core < 2:
        return {
            'outcome': 'inconclusive',
            'is_verified': False,
            'reason': 'no_core_data',
            'true_depth_min_pct': 0.0,
            'true_depth_mean_pct': 0.0,
            'consecutive_depressed': 0,
            'shoulder_elev_pct': 0.0,
            'has_flare': False
        }

    f_core = f_loc[m_core]
    f_core_mean = float(np.mean(f_core))
    f_core_min = float(np.min(f_core))

    f_left_mean = float(np.mean(f_loc[m_left])) if np.sum(m_left) > 0 else robust_cont
    f_right_mean = float(np.mean(f_loc[m_right])) if np.sum(m_right) > 0 else robust_cont
    f_left_max = float(np.max(f_loc[m_left])) if np.sum(m_left) > 0 else robust_cont
    f_right_max = float(np.max(f_loc[m_right])) if np.sum(m_right) > 0 else robust_cont

    shoulder_base = 0.5 * (f_left_mean + f_right_mean)
    shoulder_elev = shoulder_base - robust_cont
    shoulder_flare_peak = max(f_left_max, f_right_max) - robust_cont
    has_flare = (shoulder_flare_peak > 4.0 * noise) and (shoulder_flare_peak > 0.008)

    true_depth_min = robust_cont - f_core_min
    true_depth_mean = robust_cont - f_core_mean
    mf_depth_frac = mf_depth_ppm / 1e6 if mf_depth_ppm > 0 else 0.005

    # 4. Multi-cadence coherence calculation
    depressed_mask = (f_core < robust_cont - 1.5 * noise).astype(int)
    max_consecutive = 0
    curr_run = 0
    for v in depressed_mask:
        if v == 1:
            curr_run += 1
            if curr_run > max_consecutive:
                max_consecutive = curr_run
        else:
            curr_run = 0

    n_depressed = int(np.sum(depressed_mask))

    # 5. Determine Three-State Outcome
    outcome = 'verified_dimming'
    is_verified = True
    reason = 'genuine_dimming_confirmed'

    if true_depth_min < 2.0 * noise:
        outcome = 'probable_phantom'
        is_verified = False
        reason = 'no_significant_negative_excursion'
    elif max_consecutive < min_consecutive_cadences and n_core >= 3 and n_depressed <= 1:
        outcome = 'probable_phantom'
        is_verified = False
        reason = 'single_cadence_negative_outlier'
    elif (shoulder_elev > 0.5 * mf_depth_frac) and (true_depth_min < 0.5 * mf_depth_frac):
        outcome = 'probable_phantom'
        is_verified = False
        reason = 'phantom_shoulder_elevation'
    elif has_flare and (true_depth_mean < 0.25 * shoulder_flare_peak) and (max_consecutive < min_consecutive_cadences):
        outcome = 'probable_phantom'
        is_verified = False
        reason = 'flare_shoulder_artifact'

    return {
        'outcome': outcome,
        'is_verified': is_verified,
        'reason': reason,
        'true_depth_min_pct': true_depth_min * 100.0,
        'true_depth_mean_pct': true_depth_mean * 100.0,
        'consecutive_depressed': max_consecutive,
        'n_core_depressed': n_depressed,
        'n_core_total': n_core,
        'shoulder_elev_pct': shoulder_elev * 100.0,
        'shoulder_flare_peak_pct': shoulder_flare_peak * 100.0,
        'has_flare': has_flare,
        'robust_continuum': robust_cont,
        'local_noise': noise
    }

