"""TESS Experiment B0: physically informed injection and localization smoke test.

Reads existing G:/My Drive/TESS/reports/dataset_catalog.csv and processed/*.npz.
Original data never modified. No GPU and no model training required.

Install: python -m pip install numpy matplotlib batman-package
Test:    python experiment_b0.py --self-test
Run:     python experiment_b0.py --root "G:\\My Drive\\TESS" --run-name b0_01

Planet transit injections: `batman` with quadratic limb darkening, a circular
Keplerian orbit and finite-cadence exposure integration. Dusty-tail profiles:
explicitly phenomenological asymmetric proxies, NOT a dust radiative model.
The baseline detector is illustrative, uncalibrated, and not an A3.1 replica.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import sys
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

VERSION = "B0_INJECTION_LOCALIZATION_V1"
DEFAULT_ROOT = r"G:\My Drive\TESS"
DEPTH_PPM = (200, 500, 1000, 2000, 5000, 10000)
DURATION_HR = (0.5, 1, 2, 4, 8)
SHAPES = ("planet_batman", "dusty_tail_proxy")
# Candidate box widths are in HOURS; no future-test outcomes used to select them.
FILTER_WIDTHS_HR = (0.5, 1, 2, 4, 8)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_catalog(root):
    file = root / "reports" / "dataset_catalog.csv"
    if not file.is_file():
        raise FileNotFoundError(f"Dataset catalog not found: {file}")
    with file.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or not {"processed_file", "sector", "tic_id", "split"}.issubset(rows[0]):
        raise ValueError("Catalog is empty or does not match preprocess_tess.py format")
    return rows


def star_key(row):
    tic = str(row.get("tic_id", "")).strip()
    if tic and tic.lower() not in ("nan", "none", "unknown"):
        try:
            return "TIC:" + str(int(tic))
        except ValueError:
            return "TIC:" + tic
    # Unknown target IDs cannot safely be linked across sectors. Document fallback.
    return "FILE:" + str(row["processed_file"]).lower().replace("\\", "/")


def choose_catalog_rows(rows, seed, allowed_splits=("train",)):
    eligible = []
    for row in rows:
        try:
            sector = int(row.get("sector", -1))
        except ValueError:
            continue
        if (sector in range(1, 11) and row.get("status") == "ok" and
            row.get("split") in allowed_splits and row.get("processed_file")):
            eligible.append(row)
    rng = random.Random(seed)
    groups = defaultdict(list)
    for r in eligible:
        groups[star_key(r)].append(r)
    keys = sorted(groups)
    rng.shuffle(keys)
    chosen = []
    for k in keys:
        # Exactly one source observation per TIC in this pilot.
        group = sorted(groups[k], key=lambda z: (int(z["sector"]), z["processed_file"]))
        chosen.append(rng.choice(group))
    return chosen


def read_npz_window(path, length, rng, min_good=0.97):
    """Sample an actual contiguous TIME window, without compressing gaps.

    QUALITY remains available; reject large timestamp gaps rather than stitching
    disjoint cadences into an artificial sequence. Small flagged intervals outside
    the injection are allowed and remain flagged in outputs.
    """
    with np.load(path, allow_pickle=False) as z:
        needed = {"time_btjd", "pdc_flux_rel", "pdc_flux_err_rel", "good_pdc", "quality"}
        if not needed.issubset(z.files):
            raise ValueError(f"Missing NPZ fields: {needed - set(z.files)}")
        t = z["time_btjd"].astype(np.float64)
        f = z["pdc_flux_rel"].astype(np.float64)
        e = z["pdc_flux_err_rel"].astype(np.float64)
        good = z["good_pdc"].astype(bool)
        quality = z["quality"].copy()
    good &= np.isfinite(t) & np.isfinite(f) & (1+f > 0)
    if len(t) < length or np.count_nonzero(good) < length * min_good:
        raise ValueError("Insufficient good cadences")
    dt0 = np.diff(t)
    pos = dt0[np.isfinite(dt0) & (dt0 > 0) & (dt0 < 0.03)]
    if not pos.size:
        raise ValueError("Cannot determine sampling cadence")
    cadence = float(np.median(pos))
    # Real mission time steps are not evenly spaced over gaps. Never cross one.
    adjacent = np.isfinite(dt0) & (dt0 > 0) & (dt0 <= 1.55*cadence)
    cuts = np.r_[0, np.flatnonzero(~adjacent)+1, len(t)]
    candidates = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        if b-a < length:
            continue
        # Avoid scanning hundreds of thousands of all possible windows.
        for start in np.linspace(a, b-length, min(12, max(1, (b-a)//length*3)), dtype=int):
            sub = good[start:start+length]
            if np.mean(sub) >= min_good:
                candidates.append(int(start))
    if not candidates:
        raise ValueError("No contiguous, sufficiently clean observation window")
    start = int(rng.choice(candidates))
    sl = slice(start, start+length)
    return {"time": t[sl], "flux": f[sl], "flux_err": e[sl],
            "good": good[sl], "quality": quality[sl], "start": start,
            "cadence_days": cadence}


def kepler_scaled_a(period_days, mass_solar, radius_solar):
    # Kepler's third law for negligible planetary mass.
    return 215.032 * mass_solar**(1/3) * (period_days/365.25)**(2/3) / radius_solar


def duration_circular_hr(period_days, mass_solar, radius_solar, radius_ratio, impact):
    a = kepler_scaled_a(period_days, mass_solar, radius_solar)
    if a <= 1 + radius_ratio or impact >= 1+radius_ratio or impact >= a:
        return 0.0
    numerator = (1+radius_ratio)**2 - impact**2
    sin_arg = math.sqrt(max(numerator, 0)) / math.sqrt(a*a-impact*impact)
    if sin_arg >= 1:
        return float("nan")
    return 24 * period_days / math.pi * math.asin(sin_arg)


def solve_period(duration_hr, mass_solar, radius_solar, radius_ratio, impact):
    # Bounds 0.25..250 days; ensure circular orbits consistent with Kepler's law.
    lo, hi = 0.25, 250.0
    d_lo = duration_circular_hr(lo, mass_solar, radius_solar, radius_ratio, impact)
    d_hi = duration_circular_hr(hi, mass_solar, radius_solar, radius_ratio, impact)
    if (not np.isfinite(d_lo) or not np.isfinite(d_hi) or
        not d_lo < duration_hr < d_hi):
        return None
    for _ in range(65):
        mid = (lo+hi)*0.5
        if duration_circular_hr(mid, mass_solar, radius_solar, radius_ratio, impact) < duration_hr:
            lo = mid
        else:
            hi = mid
    return (lo+hi)*0.5


def build_transit(times, center, requested_depth_ppm, duration_hr, cadence_days, rng, batman):
    """Actual batman planet limb-darkened light curve at real timestamps.

    Stellar mass and radius are SIMULATED assumptions, not measured TIC values.
    Eccentricity zero, random non-grazing impact, and P solved from duration using
    Kepler's third law. Adjust Rp to approximately match a requested depth bin.
    """
    for attempt in range(90):
        r_star = float(rng.uniform(0.20, 0.42) if duration_hr <= 0.5 else
                       rng.uniform(0.9, 1.5) if duration_hr >= 8 else
                       rng.uniform(0.22, 1.35))
        m_star = r_star**1.1   # Simplified mass-radius assumption; not catalog fitted.
        b = float(rng.uniform(0.0, 0.55))
        u1, u2 = 0.30, 0.20
        if u1 + u2 >= 1:
            raise RuntimeError("Invalid limb-darkening coefficients")
        rp = math.sqrt(requested_depth_ppm * 1e-6) * 0.93
        for _ in range(5):
            period = solve_period(duration_hr, m_star, r_star, rp, b)
            if period is None:
                break
            # Exclude a second transit inside this short, isolated-event window.
            half_span = max(center-float(times[0]), float(times[-1])-center)
            if period <= half_span + duration_hr/48 + cadence_days:
                break
            a = kepler_scaled_a(period, m_star, r_star)
            inc = math.degrees(math.acos(b/a))
            p = batman.TransitParams()
            p.t0, p.per, p.rp, p.a, p.inc = center, period, rp, a, inc
            p.ecc, p.w, p.u, p.limb_dark = 0.0, 90.0, [u1, u2], "quadratic"
            # Evaluate on real cadence timestamps; expose the finite integration.
            try:
                model = batman.TransitModel(p, times, supersample_factor=5,
                                           exp_time=float(cadence_days), nthreads=1)
                transmission = np.asarray(model.light_curve(p), dtype=np.float64)
            except Exception as exc:
                raise RuntimeError(f"batman transit generation failed: {exc}") from exc
            if transmission.shape != times.shape or not np.all(np.isfinite(transmission)):
                raise RuntimeError("batman produced invalid transit model")
            realized = (1 - float(np.min(transmission))) * 1e6
            if realized <= 0:
                break
            if abs(realized / requested_depth_ppm - 1) <= 0.08:
                dur = duration_circular_hr(period, m_star, r_star, rp, b)
                return transmission, {
                    "simulated_radius_solar":r_star, "simulated_mass_solar":m_star,
                    "period_days":period, "radius_ratio":rp, "impact_parameter":b,
                    "scaled_a":a, "inclination_deg":inc, "eccentricity":0.0,
                    "limb_darkening": "quadratic", "u1":u1, "u2":u2,
                    "exposure_sec":cadence_days*86400, "supersample_factor":5,
                    "modeled_duration_hr":dur,
                    "physical_parameters_note":"Hypothetical circular planet orbit; simulated solar-unit star, not TIC-derived",
                }
            rp *= math.sqrt(requested_depth_ppm / realized)
            if not (0.004 < rp < 0.25):
                break
    raise ValueError(f"Could not construct physical circular transit for {duration_hr}h / {requested_depth_ppm}ppm")


def build_dusty_tail(times, center, depth_ppm, duration_hr, cadence_days):
    """Phenomenological smooth-fast-ingress / exponential-slow-egress proxy.

    This is NOT a physical exocomet forward model. The nominal duration is the
    >5%-of-maximum-attenuation interval, approximately; metadata records realized.
    """
    nominal_days = duration_hr / 24
    # Approximate full width at 5% attenuation: 2.45*sigma + 2.996*tau.
    sigma = nominal_days / (2.45+2.996*2.3)
    tau = sigma*2.3
    sub = (np.arange(7)-3)*cadence_days/7
    dt = times[:,None] + sub[None,:] - center
    attenuation = np.where(dt<0, np.exp(-0.5*(dt/sigma)**2), np.exp(-np.maximum(dt,0)/tau))
    transmission = 1 - depth_ppm*1e-6 * attenuation.mean(axis=1)
    return transmission, {"ingress_sigma_hr": sigma*24,
                          "tail_efold_hr": tau*24,
                          "supersample_factor":7,
                          "physical_parameters_note":"Phenomenological dusty-tail proxy, not an optical-depth/scattering forward model"}


def renormalize_observation(flux_relative, uncertainty_rel, good):
    # Inference gets only the current observation. This never receives an
    # unmodified partner's baseline, reproducing the corrected A3.1 constraint.
    flux = np.asarray(flux_relative, dtype=np.float64).copy()
    raw = 1.0+flux
    use = good & np.isfinite(raw) & (raw > 0)
    if not np.any(use):
        raise ValueError("No valid flux for independent renormalization")
    median = float(np.median(raw[use]))
    if not np.isfinite(median) or median <= 0:
        raise ValueError("Invalid observation median")
    out = np.full(len(flux), np.nan, dtype=np.float64)
    out[use] = raw[use] / median - 1
    err = np.full(len(flux), np.nan, dtype=np.float64)
    err[use] = uncertainty_rel[use] / median
    return out, err, median


def inject_flux(window, transmission):
    if transmission.shape != window["time"].shape:
        raise ValueError("Transmission/time misaligned")
    if not np.all(np.isfinite(transmission)) or np.any(transmission <= 0) or np.any(transmission > 1+1e-9):
        raise ValueError("Invalid attenuating transmission")
    base = window["flux"]
    orig_err = window["flux_err"]
    good = window["good"]
    # Apply multiplicatively to flux; uncertainty from photon measurement gets
    # the same transformation (approximation, no new Poisson noise resampling).
    injected = (1+base)*transmission - 1
    err_injected = orig_err*transmission
    control, control_err, control_median = renormalize_observation(base, orig_err, good)
    altered, altered_err, altered_median = renormalize_observation(injected, err_injected, good)
    return (control, control_err, control_median,
            altered, altered_err, altered_median)


def best_local_candidates(time, flux, good, cadence_days, top_k=10):
    """Exploratory box contrast vs neighboring sidebands; not frozen A3 filter.

    Returns located alerts. No threshold selected or false-positive claim.
    Avoid flagged samples and discard candidate boxes with <95% good cadences.
    """
    f = np.asarray(flux, dtype=np.float64)
    valid = good & np.isfinite(f) & np.isfinite(time)
    if np.count_nonzero(valid) < 0.8*len(f):
        return []
    sig = np.diff(f)
    adj = valid[1:] & valid[:-1] & (np.diff(time)>0) & (np.diff(time)<1.6*cadence_days)
    d = sig[adj]
    # Robust local scale; floor avoids pathological near-zero division.
    robust_sigma = max(1e-5, float(1.4826*np.median(np.abs(d-np.median(d)))/math.sqrt(2))) if d.size>=10 else 0.001
    candidates=[]
    for duration in FILTER_WIDTHS_HR:
        width = max(3, int(round(duration/(24*cadence_days))))
        if width > int(len(f)*.65):
            continue
        # Use shorter neighboring sidebands for long events in 512-cadence windows.
        side = min(60, max(5, width//2), max(5, (len(f)-width)//4))
        for start in range(side, len(f)-side-width+1, max(1,width//8)):
            core = slice(start,start+width)
            left = slice(start-side,start)
            right= slice(start+width,start+width+side)
            if (np.mean(valid[core]) < 0.95 or np.mean(valid[left]) < 0.9 or np.mean(valid[right]) < 0.9):
                continue
            core_mean = float(np.mean(f[core][valid[core]]))
            side_mean = .5*(float(np.mean(f[left][valid[left]]))+float(np.mean(f[right][valid[right]])))
            depth = side_mean-core_mean
            if depth <= 0:
                continue
            score = depth * math.sqrt(width)/robust_sigma
            center_idx = start + (width-1)/2
            candidates.append({"score":score, "center_btjd": float(np.interp(center_idx,np.arange(len(time)),time)),
                               "start_btjd":float(time[start]), "end_btjd":float(time[start+width-1]),
                               "width_hr": duration})
    candidates.sort(key=lambda x:x["score"], reverse=True)
    # Nonmaximum suppression by center proximity. Preserve several distinct alerts.
    selected=[]
    for c in candidates:
        if all(abs(c["center_btjd"]-s["center_btjd"]) >
               max(cadence_days*2, 0.25*min(c["width_hr"],s["width_hr"])/24) for s in selected):
            selected.append(c)
            if len(selected)>=top_k:
                break
    return selected


def event_localized(candidates, center, duration_hr, limit):
    # B protocol: localization within half the true event duration.
    tol_days = (duration_hr/2)/24
    for rank,c in enumerate(candidates[:limit],1):
        if abs(c["center_btjd"]-center) <= tol_days:
            return True, rank, abs(c["center_btjd"]-center)*24*60
    return False, None, None


def make_plots(items, dest, plot_count):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dest.mkdir(parents=True,exist_ok=True)
    for x in items[:plot_count]:
        fig,axs=plt.subplots(2,1,figsize=(11,6.3),sharex=True,constrained_layout=True)
        t=x["time"]
        good=x["good"]
        axs[0].plot(t[good],x["control"][good]*1e6,".",markersize=2,alpha=.6,label="Original / control")
        axs[0].plot(t[good],x["injected"][good]*1e6,".",markersize=2,alpha=.6,label="Injected")
        axs[0].axvspan(x["start_btjd"],x["end_btjd"],alpha=.13,label="Event interval")
        axs[0].axvline(x["center_btjd"],ls="--",alpha=.65,label="True center")
        if x["alert"]:
            axs[0].axvline(x["alert"]["center_btjd"],c="firebrick",ls=":",label="Top filter candidate")
        axs[0].set_ylabel("Relative flux (ppm)")
        axs[0].legend(loc="best",fontsize=8,ncol=3)
        axs[1].plot(t,(1-x["transmission"])*1e6,lw=1.6,label="Injected attenuation")
        axs[1].set(xlabel="TESS BTJD",ylabel="Model attenuation (ppm)")
        for ax in axs: ax.grid(alpha=.2)
        fig.suptitle(f"B0 {x['kind']}  | TIC {x['tic_id']}  | sector {x['sector']}  | target {x['target_depth_ppm']} ppm")
        fig.savefig(dest / f"{x['injection_id']}.png",dpi=130)
        plt.close(fig)


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w",encoding="utf-8",newline="") as fh:
        writer = csv.DictWriter(fh,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(root, run_name, n_transits, n_dusty, seed, length, max_files, min_good, plot_count):
    if n_transits < 0 or n_dusty < 0 or not (n_transits+n_dusty):
        raise ValueError("Request at least one transit or dusty-tail injection")
    if length < 400:
        raise ValueError("B0 requires window-length >= 400 for 8-hour events")
    if not 0.80 <= min_good <= 1:
        raise ValueError("--min-good must be in [0.80, 1]")
    try:
        import batman
    except ImportError as exc:
        raise SystemExit("Install transit package in this Python environment: python -m pip install batman-package") from exc
    out = root/"experiments"/"experiment_b"/run_name
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite existing B0 run: {out}. Use a NEW --run-name.")
    rows = read_catalog(root)
    chosen = choose_catalog_rows(rows,seed)
    print(f"Catalog {len(rows):,} observations; eligible one-per-star development sources: {len(chosen):,}")
    if not chosen:
        raise ValueError("No eligible Sectors 1–10 train/ok NPZ files")
    out.mkdir(parents=True,exist_ok=False)
    rng = np.random.default_rng(seed)
    requests = ["planet_batman"]*n_transits+["dusty_tail_proxy"]*n_dusty
    rng.shuffle(requests)
    targets=[]; candidates=[]; rowsout=[]; alertsout=[]; failures=[]
    try:
        for i,row in enumerate(chosen[:max_files],1):
            if len(targets) >= len(requests):
                break
            rel = str(row["processed_file"])
            path = root/rel
            try:
                if not path.is_file():
                    raise FileNotFoundError(str(path))
                w=read_npz_window(path,length,rng,min_good)
                kind = requests[len(targets)]
                depth=int(rng.choice(DEPTH_PPM))
                duration=float(rng.choice(DURATION_HR))
                # Center is near middle; ensure time interval fits even for 8h event.
                cidx=int(rng.integers(int(length*.46),int(length*.54)))
                center=float(w["time"][cidx])
                if kind=="planet_batman":
                    model,phys=build_transit(w["time"],center,depth,duration,w["cadence_days"],rng,batman)
                else:
                    model,phys=build_dusty_tail(w["time"],center,depth,duration,w["cadence_days"])
                realized=float((1-model.min())*1e6)
                # Truth uses model attenuation >5% of peak. Captures finite-exposure edges.
                active=np.flatnonzero((1-model)>0.05*(1-model.min()))
                if len(active)<2:
                    raise ValueError("Event insufficiently sampled at TESS cadence")
                a,b=int(active[0]),int(active[-1])
                if a<3 or b>length-4 or np.any(~w["good"][a:b+1]):
                    raise ValueError("Injection intersects invalid cadences or edge")
                start,end=float(w["time"][a]),float(w["time"][b])
                actual_duration=(end-start+w["cadence_days"])*24
                ctrl,ec,mc,inj,ei,mi = inject_flux(w,model)
                found=best_local_candidates(w["time"],inj,w["good"],w["cadence_days"])
                found_clean=best_local_candidates(w["time"],ctrl,w["good"],w["cadence_days"])
                localized_1,rank1,dist1=event_localized(found,center,actual_duration,1)
                localized_5,rank5,dist5=event_localized(found,center,actual_duration,5)
                uid=f"B0_{len(targets)+1:03d}"
                meta={"injection_id":uid, "kind":kind, "source_npz":rel,
                      "source_sha256":sha256_file(path),"tic_id":row["tic_id"],"sector":int(row["sector"]),
                      "catalog_split":row["split"], "window_start_index":w["start"],
                      "window_length":length, "n_good":int(np.sum(w["good"])),
                      "cadence_sec":w["cadence_days"]*86400,
                      "target_depth_ppm":depth, "realized_model_depth_ppm":round(realized,3),
                      "target_duration_hr":duration, "realized_model_duration_hr":round(actual_duration,4),
                      "center_btjd":center, "start_btjd":start, "end_btjd":end,
                      "baseline_control":mc,"baseline_injected":mi,
                      "top_candidate_score":found[0]["score"] if found else None,
                      "top_candidate_btjd":found[0]["center_btjd"] if found else None,
                      "control_top_score":found_clean[0]["score"] if found_clean else None,
                      "top1_localized":bool(localized_1),"top5_localized":bool(localized_5),
                      "first_localized_rank":rank5,"localization_error_minutes":dist5,
                      "model_parameters_json":json.dumps(phys,sort_keys=True)}
                rowsout.append(meta)
                for j,alert in enumerate(found,1):
                    alertsout.append({"injection_id":uid,"rank":j,**alert})
                targets.append({"injection_id":uid,"kind":kind,"tic_id":row["tic_id"],
                                "sector":int(row["sector"]),"target_depth_ppm":depth,
                                "time":w["time"],"control":ctrl,"control_err":ec,
                                "injected":inj,"injected_err":ei,"good":w["good"],
                                "quality":w["quality"],"transmission":model,
                                "center_btjd":center,"start_btjd":start,"end_btjd":end,
                                "alert":found[0] if found else None})
                print(f"[{len(targets):2d}/{len(requests)}] {kind}: sector={row['sector']} TIC={row['tic_id']} "
                      f"depth={depth}ppm duration={duration:g}h top1_localized={localized_1}",flush=True)
            except (FileNotFoundError,ValueError,OSError,KeyError,RuntimeError,OverflowError) as exc:
                failures.append({"source_npz":rel,"reason":str(exc)})
                if len(failures)<=8:
                    print(f"  Skipping {rel}: {exc}",file=sys.stderr)
        if len(targets) != len(requests):
            raise RuntimeError(f"Only built {len(targets)}/{len(requests)} injections after {min(max_files,len(chosen))} files; "
                               "increase --max-files, reduce --min-good, or inspect failures.csv")
        n=len(targets)
        arrays={}
        for key in ("time","control","control_err","injected","injected_err","good","quality","transmission"):
            arrays[key]=np.stack([o[key] for o in targets])
        arrays["injection_ids"]=np.asarray([o["injection_id"] for o in targets])
        np.savez_compressed(out/"injected_windows.npz",**arrays)
        write_csv(out/"injection_manifest.csv",rowsout)
        write_csv(out/"ranked_candidates.csv",alertsout)
        write_csv(out/"skipped_sources.csv",failures)
        make_plots(targets,out/"plots",plot_count)
        summary={"version":VERSION,"timestamp_utc":utc_now(),"mode":"development_pilot_only",
                 "external_holdout_opened":False,"input_root":str(root),"seed":seed,
                 "requested":{"planet_batman":n_transits,"dusty_tail_proxy":n_dusty},
                 "completed":{kind:sum(r["kind"]==kind for r in rowsout) for kind in SHAPES},
                 "candidate_filter":"uncalibrated local box contrast / sidebands; not A3 matched filter",
                 "localization_tolerance":"half realized >5%-depth event duration",
                 "top1_localized":sum(r["top1_localized"] for r in rowsout),
                 "top5_localized":sum(r["top5_localized"] for r in rowsout),
                 "skipped_source_count":len(failures),
                 "modeling_notes":["batman quadratic limb darkening, circular Kepler orbit, 5x finite-exposure sampling",
                                   "Stellar M,R are hypothetical, not TIC-derived; target durations/depth bins are design inputs",
                                   "Dust tails are phenomenological and do not model optical scattering",
                                   "No calibration; top-1/top-5 localization are smoke-test diagnostics, NOT science detection rates"]}
        (out/"b0_metrics.json").write_text(json.dumps(summary,indent=2,allow_nan=False),encoding="utf-8")
        text=["# Experiment B0 — Injection and localization pilot","",
              f"- Generated (UTC): {summary['timestamp_utc']}",
              f"- Simulated transits: {summary['completed']['planet_batman']}",
              f"- Dusty-tail proxies: {summary['completed']['dusty_tail_proxy']}",
              f"- One different TIC per injection: {len(targets)}",
              f"- Top-1 candidate within event: {summary['top1_localized']}/{n}",
              f"- Top-5 contains correctly timed candidate: {summary['top5_localized']}/{n}",
              f"- Rejected source attempts: {len(failures)}", "",
              "## Scope and limits", "",
              "These figures are *uncalibrated* engineering checks, not validated sensitivity or false-positive estimates.",
              "The transit family is limb-darkened batman with circular Kepler parameters from *simulated* stellar mass/radius.",
              "Dusty tails are phenomenological and not a physical exocomet model.",
              "Background TESS variability and unverified anomalies remain in the source observations.",
              "All data are Sectors 1–10, split=train, status=ok; B's proposed external holdout is never accessed.",
              "Source NPZ data and previous A datasets are unchanged.",
              "", "## Next experiment gates", "",
              "1. Inspect plots for transit shape and correct time alignment.",
              "2. Check depth and duration realized-versus-requested, gaps and quality provenance.",
              "3. Freeze location matching, then build independent calibration and time-localized detector evaluation in B1/B2.",
              "4. Obtain fresh external sectors only before the final independent B evaluation."]
        (out/"B0_REPORT.md").write_text("\n".join(text)+"\n",encoding="utf-8")
        print(f"Saved B0 results in: {out}")
        return out
    except BaseException:
        # Preserve skipped source provenance and any partial results for debugging,
        # but never write success report for an incomplete run.
        write_csv(out/"skipped_sources.csv",failures)
        raise


def self_test():
    rng=np.random.default_rng(42)
    t=1500+np.arange(512)*2/1440
    good=np.ones(len(t),dtype=bool)
    normal=np.zeros(len(t),dtype=float)
    err=np.full(len(t),1e-4)
    w={"time":t,"flux":normal,"flux_err":err,"good":good,
       "quality":np.zeros(len(t),np.uint32),"start":0,"cadence_days":2/1440}
    template,proxy=build_dusty_tail(t,t[256],0.01*1e6,2.0,w["cadence_days"])
    assert np.all((template>0)&(template<=1)) and template.min() < .995
    before,eb,mb,after,ea,ma=inject_flux(w,template)
    assert np.allclose(before[good],0) and np.nanmin(after) < -0.009
    alerts=best_local_candidates(t,after,good,2/1440)
    ok,rank,off=event_localized(alerts,t[256],2.0,5)
    assert ok and rank is not None and off < 60, (rank,off)
    assert abs(duration_circular_hr(10,1.0,1.0,0.1,0.3) -
               duration_circular_hr(solve_period(duration_circular_hr(10,1,1,.1,.3),1,1,.1,.3),1,1,.1,.3))<1e-7
    # Validate NPZ loader and original timestamps & masks without batman dependency.
    with tempfile.TemporaryDirectory() as d:
        path=Path(d)/"synthetic.npz"
        q=np.zeros(1024,np.uint32)
        q[20]=1
        tw=1500+np.arange(1024)*2/1440
        tw[900:]+=0.15  # artificial gap, should not be crossed
        fw=rng.normal(0,.0002,len(tw))
        mask=(q==0)
        np.savez_compressed(path,time_btjd=tw,pdc_flux_rel=fw,
                            pdc_flux_err_rel=np.full(len(tw),.0002),
                            good_pdc=mask,quality=q)
        sample=read_npz_window(path,512,rng,min_good=.98)
        assert len(sample["time"])==512
        assert (np.diff(sample["time"]) < 1.55*sample["cadence_days"]).all()
        assert (sample["good"].mean()>=.98)
    try:
        import batman
    except ImportError:
        print("PASS: numerical, dusty-tail, normalization, localization, FITS-derived NPZ loader tests")
        print("SKIP: actual batman transit simulation (install batman-package to test it)")
        return 0
    model,params=build_transit(t,t[256],2000,2,t[1]-t[0],rng,batman)
    assert 1700 <= (1-model.min())*1e6 <= 2200
    assert .25<=params["period_days"]<=250
    print("PASS: B0 numerical / dusty-tail / NPZ / localization AND actual batman transit self-tests")
    return 0


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--root",type=Path,default=Path(DEFAULT_ROOT))
    p.add_argument("--run-name",default="b0_01")
    p.add_argument("--n-transits",type=int,default=40)
    p.add_argument("--n-dusty",type=int,default=10)
    p.add_argument("--window-length",type=int,default=512)
    p.add_argument("--max-files",type=int,default=1200)
    p.add_argument("--min-good",type=float,default=.97)
    p.add_argument("--plot-count",type=int,default=12)
    p.add_argument("--seed",type=int,default=73)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args(argv)
    try:
        if args.self_test:
            return self_test()
        run(args.root.expanduser(),args.run_name,args.n_transits,args.n_dusty,args.seed,
            args.window_length,args.max_files,args.min_good,args.plot_count)
        return 0
    except (ValueError,FileNotFoundError,RuntimeError,FileExistsError) as exc:
        print(f"ERROR: {exc}",file=sys.stderr)
        return 1


if __name__=="__main__":
    sys.exit(main())
