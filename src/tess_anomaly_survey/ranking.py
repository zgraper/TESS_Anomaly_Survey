"""Candidate ranking, star diversity filtering, and multi-tier candidate grouping.

Provides deterministic candidate ranking by raw matched-filter score or star-normalized
significance, and supports star-level and periodic system-level export grouping.
"""
from __future__ import annotations
from collections import defaultdict
import math
import numpy as np
import sqlite3


def compute_star_statistics(candidates: list[dict]) -> dict[str, dict[str, float]]:
    """Compute robust central tendency and scale of candidate scores per TIC star.

    For stars with multiple candidates, computes median and robust MAD scale.
    For stars with single candidates, uses global candidate scale.
    """
    by_star = defaultdict(list)
    for c in candidates:
        by_star[c['tic']].append(float(c['score']))

    all_scores = [float(c['score']) for c in candidates]
    global_med = float(np.median(all_scores)) if all_scores else 0.0
    global_mad = float(np.median(np.abs(np.array(all_scores) - global_med))) * 1.4826 if all_scores else 1.0
    global_scale = max(global_mad, 1.0)

    stats = {}
    for tic, scores in by_star.items():
        arr = np.array(scores, dtype=float)
        n = len(arr)
        med = float(np.median(arr))
        if n >= 3:
            mad = float(np.median(np.abs(arr - med))) * 1.4826
            scale = mad if mad >= 1.0 else (float(np.std(arr)) if float(np.std(arr)) >= 1.0 else 1.0)
        elif n == 2:
            scale = max(float(np.std(arr)), 1.0)
        else:
            scale = global_scale
        stats[tic] = {'median': med, 'scale': scale, 'count': n, 'mean': float(np.mean(arr))}
    return stats


def compute_candidate_ranking(
    candidates: list[dict],
    method: str = "raw",
    max_per_star: int = 2,
    limit: int = 100
) -> list[dict]:
    """Rank candidates deterministically with optional star diversity and score normalization.

    Args:
        candidates: list of candidate dicts
        method: 'raw' (matched-filter score) or 'normalized' (star-normalized significance)
        max_per_star: maximum candidates to export per TIC star (0 for unrestricted)
        limit: maximum total candidates to return

    Returns:
        List of candidate dicts with added 'normalized_score' and deterministic ordering.
    """
    if not candidates:
        return []

    star_stats = compute_star_statistics(candidates)

    annotated = []
    for c in candidates:
        item = dict(c)
        tic = item['tic']
        score = float(item['score'])
        st = star_stats.get(tic, {'median': score, 'scale': 1.0, 'count': 1})
        if st['count'] >= 2:
            norm_score = (score - st['median']) / st['scale']
        else:
            norm_score = (score - st['median']) / st['scale']
        item['normalized_score'] = round(float(norm_score), 3)
        annotated.append(item)

    if method == "normalized":
        # Deterministic sorting: highest normalized score, tie-break by raw score, center_btjd, id
        sorted_cands = sorted(
            annotated,
            key=lambda x: (
                -float(x['normalized_score']),
                -float(x['score']),
                float(x.get('center_btjd', 0.0)),
                int(x.get('id', 0))
            )
        )
    else:
        # Default 'raw': highest score, tie-break by center_btjd, id
        sorted_cands = sorted(
            annotated,
            key=lambda x: (
                -float(x['score']),
                -float(x.get('normalized_score', 0.0)),
                float(x.get('center_btjd', 0.0)),
                int(x.get('id', 0))
            )
        )

    # Apply star diversity filter deterministically
    star_counts = defaultdict(int)
    ranked = []
    for c in sorted_cands:
        tic = c['tic']
        if max_per_star > 0 and star_counts[tic] >= max_per_star:
            continue
        ranked.append(c)
        star_counts[tic] += 1
        if len(ranked) >= limit:
            break

    return ranked


def group_candidates_by_star(
    candidates: list[dict],
    method: str = "raw",
    limit: int = 100
) -> list[dict]:
    """Aggregate candidates to star-level records.

    Args:
        candidates: list of candidate dicts
        method: 'raw' (rank stars by peak raw score) or 'normalized' (rank stars by peak normalized score)
        limit: maximum number of stars to return

    Returns:
        List of star-level summary dicts with deterministic ordering.
    """
    if not candidates:
        return []

    star_stats = compute_star_statistics(candidates)
    by_star = defaultdict(list)
    for c in candidates:
        by_star[c['tic']].append(c)

    star_records = []
    for tic, cands in by_star.items():
        st = star_stats[tic]
        scores = [float(x['score']) for x in cands]
        depths = [float(x.get('depth_ppm', 0.0)) for x in cands]
        widths = [float(x.get('width_hr', 0.0)) for x in cands]
        centers = [float(x.get('center_btjd', 0.0)) for x in cands]
        sectors = [int(x.get('sector', 0)) for x in cands if x.get('sector') is not None]

        # Top candidate by score
        top_cand = max(cands, key=lambda x: (float(x['score']), -float(x.get('center_btjd', 0.0))))
        max_norm = round((max(scores) - st['median']) / st['scale'], 3)

        # Aggregate review status
        statuses = set(x.get('review_status', 'unreviewed') for x in cands)
        if len(statuses) == 1:
            agg_status = list(statuses)[0]
        elif 'interesting' in statuses:
            agg_status = 'interesting'
        elif 'unreviewed' in statuses:
            agg_status = 'partially_reviewed'
        else:
            agg_status = 'reviewed'

        rec = {
            'tic': tic,
            'sector': sectors[0] if sectors else 0,
            'candidate_count': len(cands),
            'max_score': round(float(max(scores)), 2),
            'mean_score': round(float(np.mean(scores)), 2),
            'max_normalized_score': max_norm,
            'first_center_btjd': round(float(min(centers)), 4),
            'last_center_btjd': round(float(max(centers)), 4),
            'median_depth_ppm': round(float(np.median(depths)), 1),
            'median_width_hr': round(float(np.median(widths)), 2),
            'top_candidate_id': top_cand.get('id'),
            'review_status': agg_status,
        }
        star_records.append(rec)

    if method == "normalized":
        sorted_stars = sorted(
            star_records,
            key=lambda x: (-float(x['max_normalized_score']), -float(x['max_score']), str(x['tic']))
        )
    else:
        sorted_stars = sorted(
            star_records,
            key=lambda x: (-float(x['max_score']), -float(x['max_normalized_score']), str(x['tic']))
        )

    return sorted_stars[:limit]


def export_periodic_systems(
    cx: sqlite3.Connection,
    run_id: str | None = None,
    method: str = "raw",
    limit: int = 100
) -> list[dict]:
    """Retrieve and rank periodic systems from candidate_systems.

    Args:
        cx: SQLite connection
        run_id: optional run ID filter
        method: 'raw' (rank by period_snr / composite stat) or 'normalized'
        limit: maximum systems to return

    Returns:
        List of system dicts with deterministic ordering.
    """
    cond = "WHERE run_id=?" if run_id else ""
    vals = [run_id] if run_id else []

    cur = cx.execute(f"""
        SELECT id, run_id, tic, n_events, period_days, epoch_btjd,
               period_snr, is_periodic, odd_even_sigma, classification,
               note, period_uncertainty_days, epoch_uncertainty_days,
               alternate_period_days, created_utc
        FROM candidate_systems {cond}
    """, vals)

    systems = []
    for r in cur:
        snr = float(r[6]) if r[6] is not None else 0.0
        p_days = float(r[4]) if r[4] is not None else 0.0
        is_per = int(r[7])
        systems.append({
            'system_id': r[0],
            'run_id': r[1],
            'tic': r[2],
            'n_events': r[3],
            'period_days': round(p_days, 5) if p_days else None,
            'epoch_btjd': round(float(r[5]), 5) if r[5] is not None else None,
            'period_snr': round(snr, 2),
            'is_periodic': is_per,
            'odd_even_sigma': round(float(r[8]), 2) if r[8] is not None else None,
            'classification': r[9],
            'note': r[10] or '',
            'period_uncertainty_days': round(float(r[11]), 6) if r[11] is not None else None,
            'epoch_uncertainty_days': round(float(r[12]), 6) if r[12] is not None else None,
            'alternate_period_days': round(float(r[13]), 5) if r[13] is not None else None,
            'created_utc': r[14],
        })

    if not systems:
        return []

    # Deterministic sorting: periodic systems first, then by period_snr DESC, period_days ASC, tic ASC
    sorted_sys = sorted(
        systems,
        key=lambda x: (
            -int(x['is_periodic']),
            -float(x['period_snr']),
            float(x['period_days'] or 9999.0),
            str(x['tic'])
        )
    )
    return sorted_sys[:limit]
