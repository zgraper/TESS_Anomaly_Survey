"""Frozen B1 detector adapter, survey windowing, and empirical calibration.

Survey calibration is independent from the B2.2 frozen experimental operating
point. Nothing here modifies the B1 detector's candidate-score function.
"""
from __future__ import annotations
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import inspect
import json
import math
import numpy as np
import experiment_b0 as b0
import experiment_b1 as b1

ORIGINAL_B0_SHA256 = "9c8afe63ec2763b9db57c828b52b0e7d262bb30b203f86e559836dd5beb89f00"
ORIGINAL_B1_SHA256 = "c5214809a053e4705c5d86f29703e7682128e3db119494ed59f6c9a191e557c8"
METHODS = tuple(b1.METHODS)
TARGETS = (0.05, 0.01, 0.005, 0.001)


def sha256_file(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda: stream.read(1024*1024),b''):
            h.update(data)
    return h.hexdigest()


def frozen_check():
    for module, expected in ((b0, ORIGINAL_B0_SHA256),(b1,ORIGINAL_B1_SHA256)):
        actual=sha256_file(Path(inspect.getfile(module)))
        if actual != expected:
            raise RuntimeError(f"Frozen {module.__name__} source hash changed: {actual}; expected {expected}. Restore originals; don't silently alter published detectors.")


def now():
    return datetime.now(timezone.utc).isoformat()


def read_catalog(path):
    with Path(path).open('r',newline='',encoding='utf-8-sig') as fh:
        rows=list(csv.DictReader(fh))
    if not rows or not {'tic_id','processed_file','sector'}.issubset(rows[0]):
        raise ValueError("Catalog must have tic_id, sector, processed_file columns")
    seen=set()
    for r in rows:
        if not r.get('processed_file'):
            continue
        # Require an ID for cross-catalog star-disjointness; skip unknown IDs.
        raw=str(r.get('tic_id','')).strip()
        if not raw or raw.lower() in {'unknown','nan','none'}:
            continue
        key=b0.star_key(r)
        ident=(key,str(r['processed_file']).replace('\\','/'))
        if ident not in seen:
            seen.add(ident)
            yield r


def series(path):
    """Full source observation, using exactly the preprocessed B0 NPZ field names."""
    with np.load(path,allow_pickle=False) as npz:
        t=np.asarray(npz['time_btjd'],dtype=float)
        f=np.asarray(npz['pdc_flux_rel'],dtype=float)
        good=np.asarray(npz['good_pdc'],dtype=bool)
    if len(t)!=len(f) or len(f)!=len(good):
        raise ValueError("Misaligned time, flux and mask")
    good &= np.isfinite(t)&np.isfinite(f)&(1+f>0)
    gaps=np.diff(t)
    valid=gaps[np.isfinite(gaps)&(gaps>0)&(gaps<.03)]
    if len(valid)<10:
        raise ValueError("Insufficient cadence information")
    cadence=float(np.median(valid))
    if not (0<cadence<.03):
        raise ValueError("Invalid cadence")
    return t,f,good,cadence


def iter_windows(path, length=512, stride=256,min_good=.97,max_windows=0):
    """Keep temporal gaps; do not compress missing cadences or interpolate.

    Zero max_windows means scan all windows. When capped, evenly-spaced
    *time* windows are selected without referencing flux scores.
    """
    t,f,good,cadence=series(path)
    adj=(np.isfinite(np.diff(t)) & (np.diff(t)>0) &
         (np.diff(t)<=1.55*cadence))
    boundaries=np.r_[0,np.flatnonzero(~adj)+1,len(t)]
    starts=[]
    for a,b in zip(boundaries[:-1],boundaries[1:]):
        if b-a<length:
            continue
        these=list(range(int(a),int(b)-length+1,stride))
        end=int(b)-length
        if not these or these[-1]!=end:
            these.append(end)
        starts.extend(these)
    if max_windows>0 and len(starts)>max_windows:
        idx=np.linspace(0,len(starts)-1,max_windows,dtype=int)
        starts=[starts[int(i)] for i in idx]
    for st in starts:
        sl=slice(st,st+length)
        if np.mean(good[sl])<min_good:
            continue
        yield {'time':t[sl],'flux':f[sl],'good':good[sl],
               'cadence_days':cadence,'start':st,
               'window_start_btjd':float(t[st]),'window_end_btjd':float(t[st+length-1])}


def template_cache(method,cadence_days):
    if method!='planet_template':
        return None
    try:
        import batman
    except ImportError as e:
        raise RuntimeError("planet_template requires pip install batman-package") from e
    return b1.templates(cadence_days, batman)


def score_window(window,method,cache_by_cadence):
    cd=window['cadence_days']
    key=(method,round(cd,8))
    if key not in cache_by_cadence:
        cache_by_cadence[key]=template_cache(method,cd)
    raw, noise=b1.candidate_scores(window['time'],window['flux'],window['good'],
                                   cd,method,cache_by_cadence[key])
    return raw,float(noise),float(raw[0]['score'] if raw else 0.0)


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False),encoding='utf-8')


def empirical_quantile(scores,q):
    """Conservative nearest-rank cutoff; alerts require score > threshold."""
    arr=np.sort(np.asarray(scores,dtype=float))
    if not len(arr):
        raise ValueError('No calibration scores')
    index=min(len(arr)-1,max(0,int(math.ceil(q*len(arr))-1)))
    return float(arr[index])


def compute_calibration_uncertainty(per_tic, target, grouping, cutpoints, n_boot=200, seed=2026):
    """Star-level block bootstrap resampling to estimate threshold uncertainty.

    Resamples distinct TIC stars with replacement to preserve intra-star cadence
    and window correlations.
    """
    rng = np.random.default_rng(seed)
    star_keys = list(per_tic.keys())
    n_stars = len(star_keys)
    if n_stars < 2 or n_boot <= 0:
        return {}

    pooled_boots = []
    tertile_boots = {0: [], 1: [], 2: []}

    for _ in range(n_boot):
        sampled_keys = rng.choice(star_keys, size=n_stars, replace=True)
        boot_windows = []
        for k in sampled_keys:
            boot_windows.extend(per_tic[k])

        if not boot_windows:
            continue

        pooled_scores = [w['max_score'] for w in boot_windows]
        pooled_boots.append(empirical_quantile(pooled_scores, 1.0 - target))

        if grouping == 'noise_tertiles':
            noises = np.array([w['noise'] for w in boot_windows])
            tbins = np.searchsorted(cutpoints, noises, side='right')
            for g in range(3):
                g_scores = [boot_windows[i]['max_score'] for i in np.flatnonzero(tbins == g)]
                if g_scores:
                    tertile_boots[g].append(empirical_quantile(g_scores, 1.0 - target))

    def _summarize(boots):
        if not boots:
            return {'std_err': 0.0, 'ci_95': [0.0, 0.0], 'relative_stderr': 0.0}
        arr = np.sort(np.array(boots, dtype=float))
        se = float(np.std(arr))
        low = float(np.percentile(arr, 2.5))
        high = float(np.percentile(arr, 97.5))
        med = float(np.median(arr))
        rse = float(se / med) if med > 0 else 0.0
        return {
            'std_err': round(se, 4),
            'ci_95': [round(low, 4), round(high, 4)],
            'relative_stderr': round(rse, 4)
        }

    res = {
        'bootstrap_replicates': n_boot,
        'star_level_block_bootstrap': True,
        'pooled': _summarize(pooled_boots)
    }
    if grouping == 'noise_tertiles':
        res['noise_tertiles'] = {str(g): _summarize(tertile_boots[g]) for g in range(3)}
    return res


def build_calibration(rows, root, method, settings, target, grouping, min_stars, min_windows, allow_small, bootstrap=200, seed=2026):
    if not (0 < target < 1):
        raise ValueError('Alert target must be fraction strictly between zero and one')
    if grouping not in ('pooled', 'noise_tertiles'):
        raise ValueError('Unknown grouping')
    window_rows = []
    per_tic = defaultdict(list)
    failed = []
    cache = {}
    for row in rows:
        tic = b0.star_key(row)
        path = Path(root) / row['processed_file']
        try:
            for window in iter_windows(path, **settings):
                _, noise, top = score_window(window, method, cache)
                sample = {'tic': tic, 'max_score': top, 'noise': noise}
                window_rows.append(sample)
                per_tic[tic].append(sample)
        except Exception as e:
            failed.append({'tic': tic, 'file': str(path), 'error': str(e)})
    star_count = len(per_tic)
    if not window_rows:
        raise RuntimeError('No valid calibration windows')

    expected_tail = len(window_rows) * target
    effective_tail = expected_tail / 3.0 if grouping == 'noise_tertiles' else expected_tail

    if not allow_small:
        if star_count < min_stars or len(window_rows) < min_windows:
            raise RuntimeError(
                f'Calibration insufficient: {star_count} distinct stars and {len(window_rows)} windows; '
                f'required at least {min_stars} stars and {min_windows} windows. '
                f'Use --allow-small-calibration ONLY for tests, not a claimed {target*100:.2f}% rate.'
            )
        if effective_tail < 30.0:
            min_req = int(math.ceil(30.0 / (target / (3.0 if grouping == 'noise_tertiles' else 1.0))))
            raise RuntimeError(
                f'Calibration sample size inadequate for target alert fraction {target}: '
                f'expected tail events N*alpha = {effective_tail:.1f} < 30. '
                f'Requires at least {min_req} windows to reliably constrain this tail quantile. '
                f'Use --allow-small-calibration only for exploratory tests.'
            )

    allnoise = np.array([r['noise'] for r in window_rows])
    cuts = [float(np.quantile(allnoise, 1/3)), float(np.quantile(allnoise, 2/3))]
    bins = np.searchsorted(cuts, allnoise, side='right')
    thresholds = {'pooled': empirical_quantile([x['max_score'] for x in window_rows], 1 - target)}
    if grouping == 'noise_tertiles':
        for g in range(3):
            inds = np.flatnonzero(bins == g)
            starg = {window_rows[int(i)]['tic'] for i in inds}
            if not allow_small and (len(inds) < max(100, int(min_windows/3)) or len(starg) < max(100, int(min_stars/3))):
                raise RuntimeError(f'Tertile {g} has insufficient calibration diversity: {len(starg)} stars, {len(inds)} windows')
            if not len(inds):
                raise RuntimeError(f'Noise tertile {g} is empty')
            thresholds[str(g)] = empirical_quantile([window_rows[int(i)]['max_score'] for i in inds], 1 - target)

    uncertainty = {}
    if bootstrap > 0 and len(per_tic) >= 2:
        uncertainty = compute_calibration_uncertainty(
            per_tic, target, grouping, cuts, n_boot=bootstrap, seed=seed
        )

    # Compute tail-survival diagnostics for extreme quantiles
    all_scores = np.sort([r['max_score'] for r in window_rows])
    diag_quantiles = [0.90, 0.95, 0.98, 0.99, 0.995, 0.999]
    tail_survival = []
    for q in diag_quantiles:
        if q <= (1.0 - 1.0 / len(all_scores)):
            val = float(np.quantile(all_scores, q))
            tail_survival.append({
                'quantile': q,
                'target_alert_fraction': round(1.0 - q, 4),
                'score_threshold': round(val, 4),
                'n_exceedances': int(np.sum(all_scores > val))
            })

    diagnostics = {
        'expected_tail_events': round(expected_tail, 1),
        'effective_tail_events_per_bin': round(effective_tail, 1),
        'tail_stability_status': 'adequate' if effective_tail >= 30.0 else 'small_sample_exploratory',
        'tail_survival_curve': tail_survival,
    }

    return {'format_version': 'SURVEY_CALIBRATION_V1', 'created_utc': now(),
            'method': method, 'target_alert_fraction': target, 'grouping': grouping,
            'cutpoints': cuts, 'thresholds': thresholds, 'uncertainty': uncertainty,
            'diagnostics': diagnostics, 'settings': settings,
            'n_stars': star_count, 'n_windows': len(window_rows),
            'stars': sorted(per_tic), 'n_skipped_files': len(failed),
            'failed_example': failed[:30],
            'b0_sha256': ORIGINAL_B0_SHA256, 'b1_sha256': ORIGINAL_B1_SHA256,
            'calibration_window_definition': 'one maximum raw positive local score per unmodified, quality-controlled 512-cadence window; overlapping windows are correlated',
            'small_calibration_test_only': bool(allow_small and (star_count < min_stars or len(window_rows) < min_windows or effective_tail < 30.0)),
            'notes': 'Target is empirical control-window alert fraction on the calibration dataset, NOT a guaranteed survey rate or astrophysical false-discovery rate. Validate on star-disjoint observations.'}


def threshold_for(cal,noise):
    group=int(np.searchsorted(cal['cutpoints'],noise,side='right'))
    t=cal['thresholds'][str(group) if cal['grouping']=='noise_tertiles' else 'pooled']
    return float(t),group


def eligible_catalog_rows(cat,root):
    for row in read_catalog(cat):
        path=Path(root)/row['processed_file']
        if path.is_file():
            yield row
