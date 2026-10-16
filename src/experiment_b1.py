"""TESS Experiment B1: star-disjoint, calibrated *localized* injection recovery.

Development ONLY. Inputs are existing sectors 1--10 train-split TESS NPZs.
No B external holdout is read. Requires experiment_b0.py in same directory.

Compared algorithms:
  fixed_box_2h      : single 2-hour box / sideband template
  multiscale_box    : box templates 0.5/1/2/4/8 hours
  planet_template   : limb-darkened batman transit templates 0.5/1/2/4/8 hours

Threshold: 95th percentile of MAX score per presumed-background
calibration STAR. Noise tertiles permitted only if the calibration sample
contains at least --min-bin-stars unique stars in EVERY tertile; otherwise
pooled thresholds are the prespecified fallback. All event scoring uses its
OWN measured noise/bin and normalization; paired clean flux is inaccessible
to the detector. Positive scores must cross the threshold AND be localized
within half the realized event duration. All metrics are development-only.

Quick:
  python experiment_b1.py --self-test
  python experiment_b1.py --root "G:\\My Drive\\TESS" --run-name b1_smoke_01 --calibration-stars 30 --development-stars 12 --plot-count 4
Main:
  python experiment_b1.py --root "G:\\My Drive\\TESS" --run-name b1_01
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

try:
    import experiment_b0 as b0
except ImportError as e:
    raise SystemExit("Place experiment_b1.py alongside experiment_b0.py in C:\\Users\\zgrap\\Projects") from e

VERSION = "B1_STAR_DISJOINT_LOCALIZED_DEVELOPMENT_V1"
METHODS = ("fixed_box_2h", "multiscale_box", "planet_template")
KINDS = ("planet_batman", "dusty_tail_proxy")
WIDTHS_HR = (0.5, 1.0, 2.0, 4.0, 8.0)
DEPTHS = (200, 500, 1000, 2000, 5000, 10000)
DURATIONS = WIDTHS_HR


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_csv(path, records, cols=None):
    path = Path(path)
    if not records and cols is None:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=cols or list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False), encoding="utf-8")


def robust_noise(time, flux, good, cadence_days):
    good = np.asarray(good, bool) & np.isfinite(flux) & np.isfinite(time)
    consecutive = (good[:-1] & good[1:] &
                   (np.diff(time) > 0) & (np.diff(time) < 1.55 * cadence_days))
    diff = np.diff(flux)[consecutive]
    if diff.size < 10:
        raise ValueError("Too few adjacent good flux differences for noise estimate")
    sigma = 1.4826 * np.median(np.abs(diff - np.median(diff))) / math.sqrt(2)
    return float(max(sigma, 1e-6))


def templates(cadence_days, batman):
    """Cached dimensionless per-duration shapes; physically modeled planet."""
    rng = np.random.default_rng(95014)
    t = np.arange(512, dtype=float) * cadence_days
    result = {}
    for width_hr in WIDTHS_HR:
        width = max(3, int(round(width_hr / (24 * cadence_days))))
        if width >= 340:
            raise ValueError("Window too short for widest candidate and sidebands")
        transmission, _ = b0.build_transit(t, float(t[256]), 2000,
                                            float(width_hr), float(cadence_days), rng, batman)
        attenuation = 1 - transmission
        # Resample the physically generated profile to this detector's width.
        # A fixed template library is learned from no evaluated observations.
        center = 256
        x_source = (np.arange(512) - center) * cadence_days * 24
        x_target = (np.arange(width) - (width - 1) / 2) * cadence_days * 24
        shape = np.interp(x_target, x_source, attenuation, left=0, right=0)
        if not np.isfinite(shape).all() or shape.max() <= 0:
            raise RuntimeError(f"Zero batman template for width {width_hr}h")
        shape = np.clip(shape / shape.max(), 0, 1)
        result[width_hr] = shape
    return result


def candidate_scores(time, flux, good, cadence_days, method, shape_cache=None):
    """Return every raw positive template candidate; use own-observation noise.

    Baseline is determined from sidebands outside candidate event. For all
    methods the identical projection and sideband logic is used; only shapes
    and width catalogs differ. We never access the paired control observation.
    """
    if method not in METHODS:
        raise ValueError(f"Unknown method {method}")
    n = len(time)
    valid = np.asarray(good, bool) & np.isfinite(flux) & np.isfinite(time)
    if np.mean(valid) < .90:
        raise ValueError("Too many invalid points")
    noise = robust_noise(time, flux, valid, cadence_days)
    f = np.where(valid, flux, 0.0).astype(float)
    weights = valid.astype(float)
    cs = np.r_[0.0, np.cumsum(f)]
    cn = np.r_[0.0, np.cumsum(weights)]
    out = []
    widths = (2.0,) if method == "fixed_box_2h" else WIDTHS_HR
    for width_hr in widths:
        width = max(3, int(round(width_hr / (24 * cadence_days))))
        side = min(60, max(8, width // 2), max(8, (n - width) // 4))
        if side * 2 + width > n:
            continue
        starts = np.arange(side, n - width - side + 1, max(1, width // 10), dtype=int)
        if method == "planet_template":
            if shape_cache is None:
                raise RuntimeError("Transit templates not available")
            shape = shape_cache[width_hr]
        else:
            shape = np.ones(width, dtype=float)
        for st in starts:
            en = int(st + width)
            left_count = cn[st] - cn[st - side]
            right_count = cn[en + side] - cn[en]
            if left_count < .90 * side or right_count < .90 * side:
                continue
            v = valid[st:en]
            if v.mean() < .95:
                continue
            baseline = .5 * ((cs[st] - cs[st - side]) / left_count +
                             (cs[en + side] - cs[en]) / right_count)
            phi = shape * v
            norm = float(np.dot(phi, phi))
            if norm < 2:
                continue
            residual = baseline - f[st:en]
            best_depth = float(np.dot(residual, phi) / norm)
            if best_depth <= 0:
                continue
            snr = best_depth * math.sqrt(norm) / noise
            mid = float(np.interp(st + (width - 1) / 2, np.arange(n), time))
            out.append({"score": float(snr), "center_btjd": mid,
                        "start_btjd": float(time[st]), "end_btjd": float(time[en - 1]),
                        "width_hr": float(width_hr), "depth_ppm": float(best_depth * 1e6)})
    out.sort(key=lambda c: c["score"], reverse=True)
    return out, noise


def deduplicate(candidates, cadence_days, threshold=0.0):
    """NMS across widths. Select all distinct alerts above a fixed threshold.

    The suppression rule is defined WITHOUT knowing injected event position.
    Only output candidates above threshold; called with zero for unthresholded
    top-K inspection. Strongest candidate always preserved.
    """
    kept=[]
    for c in candidates:
        if c["score"] <= threshold:
            break
        if all(abs(c["center_btjd"] - x["center_btjd"]) >
               max(2 * cadence_days, .25 * min(c["width_hr"], x["width_hr"]) / 24)
               for x in kept):
            kept.append(c)
    return kept


def calibrate(calib, min_bin_stars, quantile):
    """Calibrate at star level using *one* control maximum per distinct TIC."""
    noisy = sorted(float(o["noise"]) for o in calib)
    cuts = [float(np.quantile(noisy, q)) for q in (1 / 3, 2 / 3)]
    def group(n):
        return int(np.searchsorted(cuts, n, side="right"))
    groups = [group(o["noise"]) for o in calib]
    counts = [groups.count(i) for i in range(3)]
    use_tertiles = min(counts) >= min_bin_stars and cuts[0] < cuts[1]
    thresholds = {}
    for method in METHODS:
        scores = np.asarray([o["maxima"][method] for o in calib], dtype=float)
        thresholds[method] = {"pooled":float(np.quantile(scores, quantile))}
        if use_tertiles:
            for g in range(3):
                thresholds[method][str(g)] = float(np.quantile(scores[np.asarray(groups)==g], quantile))
    return {"cutpoints":cuts,"counts":counts,"tertiles_active":bool(use_tertiles),
            "minimum_per_bin":min_bin_stars,"quantile":quantile,"thresholds":thresholds}


def applicable_threshold(cal, method, noise):
    group = int(np.searchsorted(cal["cutpoints"], noise, side="right"))
    mode = str(group) if cal["tertiles_active"] else "pooled"
    return cal["thresholds"][method][mode], group, mode


def excluded_b0(root):
    """Automatically exclude prior B0 pilot TICs wherever manifest exists."""
    excluded = set()
    paths = sorted((root / "experiments" / "experiment_b").glob("b0*/injection_manifest.csv"))
    for p in paths:
        with p.open(encoding="utf-8-sig", newline="") as f:
            for row in csv.DictReader(f):
                if str(row.get("tic_id", "")).strip():
                    excluded.add(b0.star_key(row))
                else:
                    source = (row.get("source_npz") or row.get("processed_file") or "").strip()
                    if source:
                        excluded.add("FILE:" + source.lower().replace("\\", "/"))
    return excluded, [str(p) for p in paths]


def select_rows(root, seed):
    rows = b0.choose_catalog_rows(b0.read_catalog(root), seed,
                                  allowed_splits=("train",))
    old, manifests = excluded_b0(root)
    eligible = [r for r in rows if b0.star_key(r) not in old]
    if not eligible:
        raise ValueError("No eligible training-split stars after B0 exclusions")
    return eligible, old, manifests


def choose_windows(root, rows, n_calib, n_dev, rng, length, min_good, max_files):
    """Obtain one usable contiguous window per distinct TIC; no dev calibration overlap."""
    selected = []
    rejected = []
    required = n_calib+n_dev
    for row in rows[:max_files]:
        if len(selected) >= required:
            break
        try:
            path = root / row["processed_file"]
            w = b0.read_npz_window(path, length, rng, min_good=min_good)
            # The full 8-hour planetary/dusty event must avoid flagged cadences
            # for any middle-of-window center chosen by inject_one(). Reject
            # windows with an inner bad cadence rather than silently resampling
            # only easy event durations/locations.
            if not np.all(w["good"][int(length*.20):int(length*.92)]):
                raise ValueError("Bad cadence in central long-event injection region")
            w["row"] = row
            selected.append(w)
        except (FileNotFoundError, ValueError, KeyError, OSError, RuntimeError) as e:
            rejected.append({"source_npz":row["processed_file"],"reason":str(e)})
    if len(selected) != required:
        raise RuntimeError(f"Need {required} stars but found {len(selected)} usable after "
                           f"{min(max_files,len(rows))} sources; increase --max-files or lower --min-good.")
    return selected[:n_calib], selected[n_calib:], rejected


def injection_requests(n, rng):
    """Balanced randomized blocks across shape, absolute depth and duration."""
    grid = [(k,d,h) for k in KINDS for d in DEPTHS for h in DURATIONS]
    output=[]
    while len(output)<n:
        block=list(grid)
        rng.shuffle(block)
        output.extend(block)
    return output[:n]


def inject_one(w, kind, depth, duration, rng, batman):
    n = len(w["time"])
    center_index = int(rng.integers(int(n * .46), int(n * .54)))
    center = float(w["time"][center_index])
    if kind == "planet_batman":
        transmission, phys = b0.build_transit(w["time"], center, depth, duration,
                                               w["cadence_days"], rng, batman)
    else:
        transmission, phys = b0.build_dusty_tail(w["time"], center, depth, duration,
                                                 w["cadence_days"])
    amount = 1 - transmission
    realized_ppm = float(amount.max() * 1e6)
    active = np.flatnonzero(amount > .05 * amount.max())
    if len(active) < 2:
        raise ValueError("Undersampled injected event")
    a, z = int(active[0]), int(active[-1])
    if a < 3 or z >= n - 3 or not w["good"][a:z+1].all():
        raise ValueError("Injected event intersects gap, bad flags, or window edge")
    actual_duration = (w["time"][z] - w["time"][a] + w["cadence_days"]) * 24
    control, _, _, injected, _, _ = b0.inject_flux(w, transmission)
    # b0.inject_flux produces independent normalizations for each observation.
    return {"center":center,"start":float(w["time"][a]),"end":float(w["time"][z]),
            "duration_hr":float(actual_duration),"realized_depth_ppm":realized_ppm,
            "transmission":transmission,"control":control,"injected":injected,
            "params":phys}


def score_window(w, flux, cache):
    maxima = {}
    all_scores = {}
    noise_value = None
    for method in METHODS:
        raw, noise = candidate_scores(w["time"], flux, w["good"],
                                       w["cadence_days"], method, cache)
        maxima[method] = raw[0]["score"] if raw else 0.0
        all_scores[method] = raw
        noise_value = noise
    return {"maxima":maxima,"candidates":all_scores,"noise":float(noise_value)}


def match_event(alerts, truth):
    tolerance = truth["duration_hr"] / 48   # days: 1/2 of event duration, from truth
    near = [a for a in alerts if abs(a["center_btjd"]-truth["center"]) <= tolerance]
    if not near:
        return False, None, None
    best = near[0]
    return True, float(abs(best["center_btjd"]-truth["center"]) * 24 * 60), best


def wilson(success, n, z=1.95996398454):
    if n == 0: return [None,None]
    p=success/n
    denom=1+z*z/n
    mid=(p+z*z/(2*n))/denom
    half=z/denom*math.sqrt((p*(1-p)+z*z/(4*n))/n)
    return [float(max(0,mid-half)),float(min(1,mid+half))]


def development_breakdown(events, predictions):
    by_id = {x["event_id"]: x for x in events}
    out = []
    for grouping in ("kind", "target_duration_hr", "target_depth_ppm"):
        for group in sorted(set(e[grouping] for e in events)):
            for method in METHODS:
                rows = [p for p in predictions if p["method"] == method and
                        by_id[p["event_id"]][grouping] == group]
                n = len(rows)
                detected = sum(r["localized_recovery"] for r in rows)
                alerted = sum(r["control_flagged"] for r in rows)
                out.append({"variable":grouping,"value":str(group),"method":method,
                            "n_stars":n,"localized_recovered":detected,"localized_recall":detected/n,
                            "background_alerted":alerted,"background_alert_rate":alerted/n,
                            "mean_unrelated_alerts":float(np.mean([r["n_unrelated_injected_alerts"]
                                                                      for r in rows]))})
    return out


def paired_comparison(predictions, seed):
    """Exploratory, star-paired development bootstrap; not independent final evidence."""
    ids = sorted(set(row["event_id"] for row in predictions))
    mapping = {(r["event_id"],r["method"]):int(r["localized_recovery"]) for r in predictions}
    rng = np.random.default_rng(seed+9917)
    out=[]
    for left,right in (("multiscale_box", "fixed_box_2h"),
                       ("planet_template", "multiscale_box"),
                       ("planet_template", "fixed_box_2h")):
        differences = np.asarray([mapping[(eid,left)]-mapping[(eid,right)] for eid in ids],dtype=float)
        boots = [float(np.mean(rng.choice(differences,size=len(differences),replace=True)))
                 for _ in range(1000)]
        out.append({"detector_a":left,"detector_b":right,"n_stars":len(ids),
                    "a_minus_b_recall_pp":100*float(np.mean(differences)),
                    "bootstrap_lower_pp":100*float(np.quantile(boots,.025)),
                    "bootstrap_upper_pp":100*float(np.quantile(boots,.975)),
                    "a_only":int(np.count_nonzero(differences==1)),
                    "b_only":int(np.count_nonzero(differences==-1))})
    return out


def format_report(metrics):
    s = ["# Experiment B1 — Star-disjoint calibrated localized development",
         "", f"- Generated (UTC): {metrics['timestamp_utc']}",
         f"- Distinct calibration stars: **{metrics['n_calibration_stars']}**",
         f"- Distinct development stars with successful injections: **{metrics['n_development_injections']}**",
         f"- Calibration: {'noise tertiles' if metrics['calibration']['tertiles_active'] else 'prespecified pooled fallback'}",
         f"- Target calibration background alert rate: **{metrics['target_fpr']:.1%} per 512-cadence window**",
         "- Input: Sectors 1–10, catalog split=train; external sectors not read",
         "- Injection depth chosen independently of stellar noise",
         "", "## Primary development outcomes", "",
         "| Detector | Correct-time recall | Clean background FPR | Mean unrelated alerts per injected window |",
         "|---|---:|---:|---:|"]
    for name,m in metrics["methods"].items():
        s.append(f"| {name} | {m['localized_recall']:.1%} ({m['localized_recovered']}/{m['n_injections']}) | "
                 f"{m['control_alert_rate']:.1%} ({m['control_alerted']}/{m['n_injections']}) | "
                 f"{m['mean_unrelated_alerts_per_injection']:.2f} |")
    s += ["", "## Exploratory paired method differences", "",
          "| Detector A - Detector B | Recall difference | Star bootstrap 95% interval |",
          "|---|---:|---:|"]
    for row in metrics["paired_method_comparisons"]:
        s.append(f"| {row['detector_a']} − {row['detector_b']} | "
                 f"{row['a_minus_b_recall_pp']:+.1f} pp | "
                 f"{row['bootstrap_lower_pp']:+.1f} to {row['bootstrap_upper_pp']:+.1f} pp |")
    s += ["", "## Detection by target duration (all signal families and depths)", "",
          "| Duration | Fixed 2-hour box | Multiscale box | Transit template |",
          "|---|---:|---:|---:|"]
    duration_rows=[r for r in metrics["breakdown"] if r["variable"]=="target_duration_hr"]
    for duration in sorted(set(r['value'] for r in duration_rows),key=float):
        vals={r['method']:r for r in duration_rows if r['value']==duration}
        s.append(f"| {duration} h | " + " | ".join(
            f"{vals[m]['localized_recall']:.1%} ({vals[m]['localized_recovered']}/{vals[m]['n_stars']})"
            for m in METHODS) + " |")
    s += ["", "## Known limitations", "",
          "- **Development only**, including any visually reviewed B0 stars excluded when their manifests are present.",
          "- Controls are presumed backgrounds, not confirmed astrophysical negatives.",
          "- The fixed two-hour box is a reference comparator, **not an exact implementation of A3.1's original filter**.",
          "- All three methods share local sideband contrast and noise scaling; the transit filter uses batman-derived hypothetical-star shape templates.",
          "- The transit template shares batman model assumptions with injected planetary events, favoring a matched generative family; dusty-tail results stress this assumption.",
          "- Thresholds are empirical, tail-sensitive and selected on finite calibration stars; repeated windows were not treated as independent.",
          "- No neural model is evaluated in this B1 pilot; any such comparison is a separately scoped next stage.",
          "- Dusty profiles are phenomenological and planetary host parameters are simulated, not TIC measurements.",
          "- This is neither a calibrated astrophysical false-discovery rate nor an exoplanet discovery.",
          "- Further method changes after reviewing B1 require an entirely separate, newly frozen external test for confirmatory claims.",
          "", "## Outputs", "", "See `b1_metrics.json`, `star_manifest.csv`, `calibration_scores.csv`, "
          "`development_events.csv`, `development_predictions.csv`, `b1_breakdown.csv`, "
          "`paired_comparisons.csv` and `candidate_alerts.csv`.", ""]
    return "\n".join(s)


def plot_examples(items, output, how_many):
    if how_many <= 0: return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    output.mkdir(parents=True, exist_ok=True)
    for obj in items[:how_many]:
        w=obj["w"]; truth=obj["truth"]
        t=w["time"]
        fig,ax=plt.subplots(figsize=(10.5,4.5),constrained_layout=True)
        valid=w["good"]
        ax.plot(t[valid], truth["control"][valid]*1e6, ".", ms=2.2, alpha=.60,label="Control")
        ax.plot(t[valid], truth["injected"][valid]*1e6, ".", ms=2.2, alpha=.60,label="Injected")
        ax.axvspan(truth["start"],truth["end"],alpha=.12,label="True event")
        ax.axvline(truth["center"],ls="--",lw=1.2,label="True center")
        for method, candidate in obj["top_candidates"].items():
            if candidate:
                ax.axvline(candidate["center_btjd"],ls=":",lw=1.1,label=method)
        ax.set(xlabel="TESS BTJD", ylabel="Relative PDCSAP flux (ppm)",
               title=f"B1 {obj['kind']} TIC {obj['tic']} {obj['depth']} ppm / {obj['duration']} h")
        ax.grid(alpha=.2)
        ax.legend(fontsize=7,ncol=2)
        fig.savefig(output/f"B1_{obj['index']:03d}.png",dpi=140)
        plt.close(fig)


def run(args):
    if args.calibration_stars < 10 or args.development_stars < 6:
        raise ValueError("At least 10 calibration and 6 development stars required")
    if not (0.0 < args.target_fpr < 0.25):
        raise ValueError("--target-fpr must be between 0 and 0.25")
    if args.min_bin_stars < 5:
        raise ValueError("--min-bin-stars must be >=5")
    if args.window_length != 512:
        raise ValueError("B1 v1 requires 512-cadence windows; duration templates are fixed to 512")
    try: import batman
    except ImportError as e: raise RuntimeError("batman not importable: use your astroenv-b0 Python 3.13 environment") from e
    root=args.root.expanduser()
    out=root/"experiments"/"experiment_b"/args.run_name
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite B1 results: {out}. Choose a new --run-name")
    available,excluded,excluded_paths=select_rows(root,args.seed)
    rng=np.random.default_rng(args.seed)
    calibrators,dev,skipped=choose_windows(root,available,args.calibration_stars,
                              args.development_stars,rng,args.window_length,args.min_good,args.max_files)
    cal_stars = {b0.star_key(w["row"]) for w in calibrators}
    dev_stars = {b0.star_key(w["row"]) for w in dev}
    assert not (cal_stars & dev_stars),"calibration/development TIC leakage"
    assert not ((cal_stars|dev_stars)&excluded),"B0 pilot TIC included"
    cache=templates(float(np.median([x["cadence_days"] for x in calibrators+dev])),batman)
    # Cache source templates for nominal 2-minute cadence. All selected source
    # cadences must be close: avoid inadvertently applying mismatched templates.
    cad0=float(np.median([x["cadence_days"] for x in calibrators+dev]))
    if any(abs(w["cadence_days"]/cad0-1)>.04 for w in calibrators+dev):
        raise ValueError("Mixed cadence populations >4%. Split cohorts by cadence before evaluation")
    out.mkdir(parents=True,exist_ok=False)
    star_manifest=[]; cal_records=[]; calibrations=[]
    for i,w in enumerate(calibrators,1):
        f,_,_=b0.renormalize_observation(w["flux"],w["flux_err"],w["good"])
        scores=score_window(w,f,cache)
        item={"star_key":b0.star_key(w["row"]),"noise":scores["noise"],"maxima":scores["maxima"]}
        calibrations.append(item)
        r={"role":"calibration","star_key":item["star_key"],"tic_id":w["row"]["tic_id"],
           "sector":w["row"]["sector"],"source_npz":w["row"]["processed_file"],
           "window_start":w["start"],"noise_ppm":round(scores["noise"]*1e6,2)}
        star_manifest.append(r)
        cal_records.append({**r,**{name+"_max":float(scores["maxima"][name]) for name in METHODS}})
        if i%50==0 or i==len(calibrators): print(f"Calibrating {i}/{len(calibrators)} stars",flush=True)
    cal=calibrate(calibrations,args.min_bin_stars,1-args.target_fpr)
    print(f"Calibrated: group sizes {cal['counts']}, mode "
          f"{'noise_tertiles' if cal['tertiles_active'] else 'pooled'}",flush=True)
    combos=injection_requests(len(dev),rng)
    events=[];predictions=[];alerts_rows=[]; plots=[];inject_failures=[]
    for i,(w,(kind,depth,duration)) in enumerate(zip(dev,combos),1):
        # Keep target fixed; if geometry/gap rejected, select a different
        # event center in the SAME star, not another calibration star.
        for attempt in range(30):
            try:
                truth=inject_one(w,kind,depth,duration,rng,batman)
                break
            except (ValueError,RuntimeError) as e:
                if attempt == 29:
                    raise RuntimeError(f"Cannot inject {kind}/{depth}ppm/{duration}h into "
                                       f"TIC {w['row']['tic_id']} after 30 tries: {e}") from e
        tic=b0.star_key(w["row"])
        control=score_window(w,truth["control"],cache)
        injected=score_window(w,truth["injected"],cache)
        star_manifest.append({"role":"development","star_key":tic,"tic_id":w["row"]["tic_id"],
          "sector":w["row"]["sector"],"source_npz":w["row"]["processed_file"],
          "window_start":w["start"],"noise_ppm":round(control["noise"]*1e6,2)})
        eid=f"B1_{i:04d}"
        events.append({"event_id":eid,"star_key":tic,"tic_id":w["row"]["tic_id"],
          "sector":w["row"]["sector"],"source_npz":w["row"]["processed_file"],
          "kind":kind,"target_depth_ppm":depth,"target_duration_hr":duration,
          "realized_depth_ppm":truth["realized_depth_ppm"],"realized_duration_hr":truth["duration_hr"],
          "center_btjd":truth["center"],"event_start_btjd":truth["start"],
          "event_end_btjd":truth["end"],"physical_params_json":json.dumps(truth["params"],sort_keys=True)})
        top_candidates={}
        for method in METHODS:
            tc,group_c,mode_c=applicable_threshold(cal,method,control["noise"])
            tp,group_p,mode_p=applicable_threshold(cal,method,injected["noise"])
            neg_alerts=deduplicate(control["candidates"][method],w["cadence_days"],tc)
            pos_alerts=deduplicate(injected["candidates"][method],w["cadence_days"],tp)
            localized,offset,candidate=match_event(pos_alerts,truth)
            unrelated=[a for a in pos_alerts if abs(a["center_btjd"]-truth["center"])>
                       truth["duration_hr"]/48]
            top_candidates[method]=injected["candidates"][method][0] if injected["candidates"][method] else None
            predictions.append({"event_id":eid,"method":method,"kind":kind,"target_depth_ppm":depth,
                "target_duration_hr":duration,"noise_ppm_injected":injected["noise"]*1e6,
                "noise_ppm_control":control["noise"]*1e6,"noise_group_control":group_c,
                "noise_group_injected":group_p,"mode_control":mode_c,"mode_injected":mode_p,
                "threshold_control":tc,"threshold_injected":tp,
                "max_score_control":control["maxima"][method],
                "max_score_injected":injected["maxima"][method],
                "control_flagged":bool(neg_alerts),"localized_recovery":bool(localized),
                "localization_error_min":offset,"n_control_alerts":len(neg_alerts),
                "n_injected_alerts":len(pos_alerts),"n_unrelated_injected_alerts":len(unrelated),
                "top_candidate_btjd":top_candidates[method]["center_btjd"] if top_candidates[method] else ""})
            for observation,candidates,threshold in (("control",control["candidates"][method],tc),
                                                      ("injected",injected["candidates"][method],tp)):
                active=deduplicate(candidates,w["cadence_days"],threshold)
                for rank,a in enumerate(active,1):
                    alerts_rows.append({"event_id":eid,"method":method,"observation":observation,
                        "alert_rank":rank,"threshold":threshold,"localized_truth":
                        (bool(abs(a["center_btjd"]-truth["center"])<=truth["duration_hr"]/48)
                         if observation=="injected" else ""),**a})
        if i <= args.plot_count:
            plots.append({"index":i,"w":w,"truth":truth,"kind":kind,
                          "tic":w["row"]["tic_id"],"depth":depth,"duration":duration,
                          "top_candidates":top_candidates})
        if i%10==0 or i==len(dev): print(f"Evaluating {i}/{len(dev)} injections",flush=True)
    summary={}
    for method in METHODS:
        rr=[r for r in predictions if r["method"]==method]
        k=sum(r["localized_recovery"] for r in rr)
        controls=sum(r["control_flagged"] for r in rr)
        summary[method]={"n_injections":len(rr),"localized_recovered":k,
             "localized_recall":k/len(rr),"recall_wilson_95":wilson(k,len(rr)),
             "control_alerted":controls,"control_alert_rate":controls/len(rr),
             "control_alert_wilson_95":wilson(controls,len(rr)),
             "mean_unrelated_alerts_per_injection":float(np.mean([r["n_unrelated_injected_alerts"] for r in rr])),
             "total_unrelated_alerts":sum(r["n_unrelated_injected_alerts"] for r in rr)}
    breakdown=development_breakdown(events,predictions)
    pairs=paired_comparison(predictions,args.seed)
    metrics={"version":VERSION,"timestamp_utc":utc(),"development_only":True,
        "external_holdout_opened":False,"calibration_split":"one control window per unique TIC; sectors 1-10 train only",
        "n_calibration_stars":len(calibrators),"n_development_injections":len(dev),
        "target_fpr":args.target_fpr,"methods":summary,"calibration":cal,
        "breakdown":breakdown,"paired_method_comparisons":pairs,
        "n_b0_pilot_stars_excluded":len(excluded),"b0_exclusion_manifests":excluded_paths,
        "n_sources_skipped":len(skipped),"seed":args.seed,
        "code_sha256":b0.sha256_file(Path(__file__)),
        "b0_code_sha256":b0.sha256_file(Path(b0.__file__)),
        "catalog_sha256":b0.sha256_file(root/"reports"/"dataset_catalog.csv"),
        "shape_template_kind":"batman quadratic limb-darkened circular hypothetical-star transit at median cadence",
        "scoring":"sideband-depth projection / point-to-point robust sigma; per observed window",
        "localization_rule":"alert center within half of observed >5% model attenuation duration",
        "injection_strength":"target ppm chosen from fixed grid, independent of measured stellar noise",
        "controls":"unmodified real TESS backgrounds are presumed rather than verified negatives",
        "data_scope":"no files outside original sectors 1-10 train split read"}
    write_json(out/"b1_metrics.json",metrics)
    write_csv(out/"star_manifest.csv",star_manifest)
    write_csv(out/"calibration_scores.csv",cal_records)
    write_csv(out/"development_events.csv",events)
    write_csv(out/"development_predictions.csv",predictions)
    write_csv(out/"b1_breakdown.csv",breakdown)
    write_csv(out/"paired_comparisons.csv",pairs)
    write_csv(out/"candidate_alerts.csv",alerts_rows,cols=["event_id","method","observation","alert_rank",
        "threshold","localized_truth","score","center_btjd","start_btjd","end_btjd","width_hr","depth_ppm"])
    write_csv(out/"skipped_sources.csv",skipped,cols=["source_npz","reason"])
    plot_examples(plots,out/"plots",args.plot_count)
    (out/"B1_REPORT.md").write_text(format_report(metrics),encoding="utf-8")
    print("Saved:",out,flush=True)
    return metrics


def fake_batman():
    """Synthetic TEST ONLY model for self-test if installed batman unavailable."""
    class TransitParams: pass
    class TransitModel:
        def __init__(self,p,times,**kwargs): self.t=np.asarray(times); self.t0=p.t0
        def light_curve(self,p):
            # Designed ONLY to check numerical integration/interface, not physical correctness.
            b = p.a * math.cos(math.radians(p.inc))
            ang=math.sqrt(max(0,(1+p.rp)**2-b*b))/math.sqrt(p.a*p.a-b*b)
            dur=p.per / math.pi * math.asin(min(1,ang))
            z=(self.t-self.t0)/max(dur,1e-6)
            return 1 - p.rp*p.rp*np.maximum(0,1-(2*z)**8)
    class Fake: pass
    Fake.TransitParams=TransitParams; Fake.TransitModel=TransitModel
    return Fake


def self_test():
    try:
        import batman
        has_batman=True
    except ImportError:
        batman=fake_batman();has_batman=False
    rng=np.random.default_rng(782)
    c=2/1440
    n=512
    t=1500+np.arange(n)*c
    good=np.ones(n,bool)
    base=rng.normal(0,.00035,n)
    e=np.full(n,.00035)
    w={"time":t,"flux":base,"flux_err":e,"good":good,
       "cadence_days":c,"start":0,"row":{"tic_id":"11","sector":"1","processed_file":"s.npz"}}
    cache=templates(c,batman)
    normal,_,_=b0.renormalize_observation(base,e,good)
    original=score_window(w,normal,cache)
    assert len(original["maxima"])==3 and all(x>=0 for x in original["maxima"].values())
    transmission,_=b0.build_dusty_tail(t,float(t[256]),10000,2,c)
    ctrl,_,_,inj,_,_=b0.inject_flux(w,transmission)
    assert np.allclose(ctrl,normal)
    scores=score_window(w,inj,cache)
    assert scores["maxima"]["multiscale_box"]>original["maxima"]["multiscale_box"]
    candidates=deduplicate(scores["candidates"]["multiscale_box"],c,threshold=0)
    truth={"center":float(t[256]),"duration_hr":2.0}
    assert match_event(candidates[:10],truth)[0]
    fake_cal=[]
    for i in range(60):
        fake_cal.append({"noise":i+.01,"maxima":{k: float(i%15) for k in METHODS}})
    cal=calibrate(fake_cal,15,.95)
    assert cal["tertiles_active"] and cal["counts"]==[20,20,20]
    assert not calibrate(fake_cal,30,.95)["tertiles_active"]
    assert applicable_threshold(cal,METHODS[0],1)[1]==0
    # no leakage from synthetic injected source into normalized control
    assert np.max(np.abs(normal-ctrl))<1e-12
    # Use actual batman if available; no fake physical validation claims.
    if has_batman:
        tr,_=b0.build_transit(t,float(t[256]),2000,2,c,rng,batman)
        assert 1600 < (1-tr.min())*1e6 < 2400
    print("PASS: B1 filtering, calibration/tertile fallback, localization, independent renormalization")
    print("PASS: actual batman transit tests" if has_batman else
          "SKIP: physical batman test (not installed here); fake interface tested for numerical plumbing only")
    return 0


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root",type=Path,default=Path(b0.DEFAULT_ROOT))
    p.add_argument("--run-name",default="b1_01")
    p.add_argument("--calibration-stars",type=int,default=330)
    p.add_argument("--development-stars",type=int,default=180)
    p.add_argument("--min-bin-stars",type=int,default=100)
    p.add_argument("--max-files",type=int,default=2500)
    p.add_argument("--min-good",type=float,default=.97)
    p.add_argument("--window-length",type=int,default=512)
    p.add_argument("--target-fpr",type=float,default=.05)
    p.add_argument("--seed",type=int,default=113)
    p.add_argument("--plot-count",type=int,default=12)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args(argv)
    try:
        return self_test() if args.self_test else (run(args) and 0)
    except (RuntimeError,ValueError,OSError,FileExistsError,AssertionError) as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}",file=sys.stderr)
        return 1


if __name__=="__main__":
    sys.exit(main())
