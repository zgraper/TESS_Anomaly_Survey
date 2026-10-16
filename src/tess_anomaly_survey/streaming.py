"""Streaming ingest and real-time processing of TESS light curves from MAST.

Processes NASA TESS SPOC 2-minute cadence light curves directly from MAST into
a dedicated Phase 7 SQLite results database, enforcing temporary-file cleanup,
bounded parallel acquisition, idempotent resumption, strict development-universe
guards, refined 5-state scientific taxonomy, and high-resolution performance telemetry.
"""
from __future__ import annotations

import concurrent.futures
import csv
import hashlib
import json
import logging
import math
import os
import queue
import re
import shutil
import sqlite3
import sys
import tempfile
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple

import numpy as np
import requests
try:
    import psutil
except ImportError:
    psutil = None
from astropy.io import fits

import experiment_b0 as b0
import experiment_b1 as b1
from . import __version__, core, db, periodicity, vetting

logger = logging.getLogger("tess_streaming")

# Lightkurve / SPOC standard default quality bitmask (excludes severe cadences)
TESS_DEFAULT_QUALITY_MASK = 17087

# Calibrated nominal 1% alert threshold for fixed_box_2h detector
FROZEN_THRESHOLD_T = 73.77965536097444

# Frozen pipeline configuration hash from Phase 6B
FROZEN_CONFIG_HASH = "dc32c430b5aae7d8af069a6359efd8c49e06736994c600f258f2852eb3269d3c"

MAST_DIRECT_DOWNLOAD_URL = "https://mast.stsci.edu/api/v0.1/Download/file/?uri="


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def init_streaming_db(db_path: str) -> sqlite3.Connection:
    """Create or connect to the dedicated Phase 7 results database with strict schemas."""
    os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)
    con = sqlite3.connect(db_path, timeout=60.0)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    con.execute("PRAGMA foreign_keys=ON")

    with con:
        con.execute("""
        CREATE TABLE IF NOT EXISTS runs (
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
        )""")

        con.execute("""
        CREATE TABLE IF NOT EXISTS observations (
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
            status TEXT NOT NULL, -- 'completed', 'no_anomalies', 'failed', 'skipped'
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
            created_utc TEXT NOT NULL,
            UNIQUE(run_id, tic, sector, product_filename)
        )""")

        con.execute("""
        CREATE TABLE IF NOT EXISTS candidates (
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
            created_utc TEXT NOT NULL,
            FOREIGN KEY(observation_id) REFERENCES observations(id)
        )""")

        con.execute("""
        CREATE TABLE IF NOT EXISTS candidate_detections (
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
        )""")

        con.execute("""
        CREATE TABLE IF NOT EXISTS failures (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            run_id TEXT NOT NULL,
            tic TEXT,
            sector INTEGER,
            stage TEXT NOT NULL,
            error_message TEXT NOT NULL,
            traceback TEXT,
            created_utc TEXT NOT NULL
        )""")

        con.execute("""
        CREATE TABLE IF NOT EXISTS performance_benchmarks (
            run_id TEXT PRIMARY KEY,
            worker_count INTEGER NOT NULL DEFAULT 1,
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
        )""")

        # Migration helper: ensure source_file, review_status, note exist if table was previously created
        cur = con.cursor()
        cols = [r[1] for r in cur.execute("PRAGMA table_info(candidates)").fetchall()]
        if "source_file" not in cols:
            cur.execute("ALTER TABLE candidates ADD COLUMN source_file TEXT NOT NULL DEFAULT ''")
        if "review_status" not in cols:
            cur.execute("ALTER TABLE candidates ADD COLUMN review_status TEXT NOT NULL DEFAULT 'unreviewed'")
        if "note" not in cols:
            cur.execute("ALTER TABLE candidates ADD COLUMN note TEXT NOT NULL DEFAULT ''")
        
        obs_cols = [r[1] for r in cur.execute("PRAGMA table_info(observations)").fetchall()]
        if "obs_id" not in obs_cols:
            cur.execute("ALTER TABLE observations ADD COLUMN obs_id TEXT")
        if "file_sha256" not in obs_cols:
            cur.execute("ALTER TABLE observations ADD COLUMN file_sha256 TEXT")
        if "calib_level" not in obs_cols:
            cur.execute("ALTER TABLE observations ADD COLUMN calib_level INTEGER DEFAULT 2")
        if "prvversion" not in obs_cols:
            cur.execute("ALTER TABLE observations ADD COLUMN prvversion TEXT")
        if "config_hash" not in obs_cols:
            cur.execute("ALTER TABLE observations ADD COLUMN config_hash TEXT NOT NULL DEFAULT ''")

    return con


def load_development_universe(manifest_path: str = "research/phase6b/excluded_development_stars.csv") -> Set[str]:
    """Load the verified 8,984-star development universe to prevent data leakage."""
    p = Path(manifest_path)
    if not p.is_file():
        raise FileNotFoundError(f"Development universe manifest missing: {manifest_path}")
    tics = set()
    with p.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            raw = str(r.get("tic_id", "")).replace("TIC:", "").strip()
            if raw:
                tics.add(raw)
    return tics


def check_development_guard(tic: str, dev_tics: Set[str]) -> None:
    """Enforce strict isolation between development universe and independent validation cohort."""
    t = str(tic).replace("TIC:", "").strip()
    if t not in dev_tics:
        raise RuntimeError(f"DATA LEAKAGE GUARD: TIC {tic} is not in the frozen development universe of 8,984 stars! Access denied.")


def resolve_mast_product(
    tic: str,
    sector: Optional[int] = None,
    astroquery_first: bool = True,
    session: Optional[requests.Session] = None
) -> Dict[str, Any]:
    """Query MAST to obtain the official SPOC 2-minute cadence light curve product metadata."""
    clean_tic = str(tic).replace("TIC:", "").strip()
    
    # 1. Try Astroquery if requested
    if astroquery_first:
        try:
            from astroquery.mast import Observations
            query_kwargs: Dict[str, Any] = {
                "provenance_name": "SPOC",
                "target_name": clean_tic,
            }
            if sector is not None:
                query_kwargs["sequence_number"] = int(sector)
            
            obs = Observations.query_criteria(**query_kwargs)
            if len(obs) > 0:
                prod = Observations.get_product_list(obs[0])
                spoc_lc = prod[prod["productSubGroupDescription"] == "LC"]
                if len(spoc_lc) > 0:
                    row = spoc_lc[0]
                    filename = str(row["productFilename"])
                    data_uri = str(row["dataURI"])
                    download_url = f"{MAST_DIRECT_DOWNLOAD_URL}{data_uri}"
                    sec = int(obs[0]["sequence_number"]) if "sequence_number" in obs[0].colnames else (sector or -1)
                    obs_id_val = str(obs[0]["obs_id"]) if "obs_id" in obs[0].colnames else ""
                    prv_ver = str(row.get("prvversion", ""))
                    cal_lvl = int(row.get("calib_level", 2))
                    return {
                        "tic": clean_tic,
                        "sector": sec,
                        "obs_id": obs_id_val,
                        "product_filename": filename,
                        "data_uri": data_uri,
                        "download_url": download_url,
                        "prvversion": prv_ver,
                        "calib_level": cal_lvl,
                        "method": "astroquery"
                    }
        except Exception as e:
            logger.debug(f"Astroquery resolution for TIC {clean_tic} failed ({e}), falling back to direct MAST API.")

    # 2. Direct MAST Mashup API query
    url = "https://mast.stsci.edu/api/v0/invoke"
    params = {
        "request": json.dumps({
            "service": "Mast.Caom.Filtered",
            "format": "json",
            "params": {
                "columns": "obsid,obs_id,sequence_number,target_name",
                "filters": [
                    {"paramName": "provenance_name", "values": ["SPOC"]},
                    {"paramName": "target_name", "values": [clean_tic]},
                    {"paramName": "dataproduct_type", "values": ["timeseries"]}
                ]
            }
        })
    }
    s = session or requests.Session()
    r = s.get(url, params=params, timeout=20)
    r.raise_for_status()
    data = r.json()
    records = data.get("data", [])
    if sector is not None:
        records = [rec for rec in records if int(rec.get("sequence_number", -1)) == int(sector)]
    
    if not records:
        raise FileNotFoundError(f"No SPOC observations found on MAST for TIC {clean_tic} (sector={sector})")

    obs_id = records[0]["obsid"]
    sec = int(records[0]["sequence_number"])

    # Query products for this obs_id
    p_params = {
        "request": json.dumps({
            "service": "Mast.Caom.Products",
            "format": "json",
            "params": {"obsid": obs_id}
        })
    }
    pr = s.get(url, params=p_params, timeout=20)
    pr.raise_for_status()
    prods = pr.json().get("data", [])
    lc_prods = [p for p in prods if p.get("productSubGroupDescription") == "LC" or str(p.get("productFilename", "")).endswith("_lc.fits")]
    if not lc_prods:
        raise FileNotFoundError(f"No LC product found for TIC {clean_tic} obs_id {obs_id}")

    target_prod = lc_prods[0]
    filename = target_prod["productFilename"]
    data_uri = target_prod["dataURI"]
    download_url = f"{MAST_DIRECT_DOWNLOAD_URL}{data_uri}"

    return {
        "tic": clean_tic,
        "sector": sec,
        "obs_id": str(obs_id),
        "product_filename": filename,
        "data_uri": data_uri,
        "download_url": download_url,
        "prvversion": str(target_prod.get("prvversion", "")),
        "calib_level": int(target_prod.get("calib_level", 2)),
        "method": "mashup_api"
    }


# Backwards-compatible alias
resolve_mast_spoc_product = resolve_mast_product


def download_light_curve(
    download_url: str,
    dest_path: str,
    session: Optional[requests.Session] = None,
    timeout: float = 45.0,
    max_retries: int = 3,
    backoff_factor: float = 1.5
) -> Tuple[int, float, str]:
    """Download FITS product from MAST with retries, exponential backoff, SHA-256 checksum, and size validation."""
    s = session or requests.Session()
    last_err = None
    t0 = time.perf_counter()

    # Handle local file paths and file:// URLs directly for offline/testing support
    if download_url.startswith("file://") or os.path.isfile(download_url):
        local_src = download_url.replace("file://", "")
        if os.name == 'nt' and local_src.startswith('/') and len(local_src) > 2 and local_src[2] == ':':
            local_src = local_src[1:] # Strip leading slash on Windows /C:/...
        if os.path.abspath(local_src) != os.path.abspath(dest_path):
            shutil.copyfile(local_src, dest_path)
        bytes_written = os.path.getsize(dest_path)
        h = hashlib.sha256()
        with open(dest_path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(chunk)
        elapsed = time.perf_counter() - t0
        return bytes_written, elapsed, h.hexdigest()

    for attempt in range(1, max_retries + 1):
        try:
            r = s.get(download_url, stream=True, timeout=timeout)
            r.raise_for_status()
            bytes_written = 0
            h = hashlib.sha256()
            with open(dest_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=65536):
                    if chunk:
                        f.write(chunk)
                        h.update(chunk)
                        bytes_written += len(chunk)
            
            if bytes_written < 10000:
                raise ValueError(f"Downloaded file too small ({bytes_written} bytes); possible truncated download")

            # Validate FITS header
            with fits.open(dest_path) as hdul:
                if len(hdul) < 2 or hdul[1].data is None:
                    raise ValueError("Downloaded file is not a valid TESS light curve FITS file")

            elapsed = time.perf_counter() - t0
            return bytes_written, elapsed, h.hexdigest()

        except Exception as e:
            last_err = e
            if os.path.exists(dest_path):
                try:
                    os.remove(dest_path)
                except OSError:
                    pass
            if attempt < max_retries:
                sleep_time = backoff_factor ** attempt
                logger.warning(f"Download attempt {attempt} failed ({e}). Retrying in {sleep_time:.1f}s...")
                time.sleep(sleep_time)

    raise RuntimeError(f"Failed to download {download_url} after {max_retries} attempts: {last_err}") from last_err


def extract_fits_light_curve(fits_path: str) -> Dict[str, Any]:
    """Extract and calibrate time, relative normalized flux, and quality mask directly from FITS."""
    with fits.open(fits_path) as hdul:
        h0 = hdul[0].header
        h1 = hdul[1].header
        table = hdul[1].data

        camera = int(h0.get("CAMERA", 0))
        ccd = int(h0.get("CCD", 0))
        v_tmag = h0.get("TESSMAG")
        tmag = float(v_tmag) if v_tmag is not None else math.nan
        v_teff = h0.get("TEFF")
        teff = float(v_teff) if v_teff is not None else math.nan
        sector = int(h0.get("SECTOR", 0))
        tic = str(h0.get("TICID", ""))

        time_raw = np.asarray(table["TIME"], dtype=float)
        pdc_raw = np.asarray(table["PDCSAP_FLUX"], dtype=float)
        quality = np.asarray(table["QUALITY"], dtype=np.uint32)

    good_time = np.isfinite(time_raw)
    bad_quality = (quality & TESS_DEFAULT_QUALITY_MASK) != 0
    good_base = good_time & ~bad_quality & np.isfinite(pdc_raw)

    n_cadences = len(time_raw)
    n_usable = int(np.sum(good_base))

    if n_usable < 512:
        return {
            "status": "insufficient_usable_cadences",
            "tic": tic, "sector": sector, "camera": camera, "ccd": ccd,
            "tmag": tmag, "teff": teff, "n_cadences": n_cadences, "n_usable": n_usable,
            "time": time_raw, "flux": np.zeros_like(time_raw), "good": good_base, "cadence": 2.0 / 1440.0
        }

    med_pdc = float(np.median(pdc_raw[good_base]))
    if not (np.isfinite(med_pdc) and med_pdc > 0):
        return {
            "status": "invalid_median_flux",
            "tic": tic, "sector": sector, "camera": camera, "ccd": ccd,
            "tmag": tmag, "teff": teff, "n_cadences": n_cadences, "n_usable": n_usable,
            "time": time_raw, "flux": np.zeros_like(time_raw), "good": good_base, "cadence": 2.0 / 1440.0
        }

    flux_norm = np.full_like(pdc_raw, np.nan, dtype=np.float32)
    flux_norm[good_base] = (pdc_raw[good_base] / med_pdc - 1.0).astype(np.float32)
    good = good_base & (1.0 + flux_norm > 0)

    gaps = np.diff(time_raw[good])
    valid_gaps = gaps[np.isfinite(gaps) & (gaps > 0) & (gaps < 0.03)]
    cadence = float(np.median(valid_gaps)) if len(valid_gaps) >= 10 else 2.0 / 1440.0

    return {
        "status": "ok",
        "tic": tic,
        "sector": sector,
        "camera": camera,
        "ccd": ccd,
        "tmag": tmag,
        "teff": teff,
        "n_cadences": n_cadences,
        "n_usable": int(np.sum(good)),
        "time": time_raw,
        "flux": flux_norm,
        "good": good,
        "cadence": cadence
    }


def iter_fits_windows(time_arr, flux_arr, good_arr, cadence, length=512, stride=256, min_good=0.97):
    """Slice contiguous 512-cadence time windows from FITS time-series without gap interpolation."""
    adj = (np.isfinite(np.diff(time_arr)) & (np.diff(time_arr) > 0) &
           (np.diff(time_arr) <= 1.55 * cadence))
    boundaries = np.r_[0, np.flatnonzero(~adj) + 1, len(time_arr)]
    starts = []
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        if b - a < length:
            continue
        these = list(range(int(a), int(b) - length + 1, stride))
        end = int(b) - length
        if not these or these[-1] != end:
            these.append(end)
        starts.extend(these)

    for st in starts:
        sl = slice(st, st + length)
        if np.mean(good_arr[sl]) < min_good:
            continue
        yield {
            "time": time_arr[sl],
            "flux": flux_arr[sl],
            "good": good_arr[sl],
            "cadence_days": cadence,
            "start": st,
            "window_start_btjd": float(time_arr[st]),
            "window_end_btjd": float(time_arr[st + length - 1])
        }


def detect_continuous_variability(time_arr: np.ndarray, flux_arr: np.ndarray, good_arr: np.ndarray, n_events: int) -> bool:
    """Evaluate whether star exhibits dense continuous periodic or rotational modulation."""
    if n_events < 8:
        return False
    # Check span and density of good flux
    valid_t = time_arr[good_arr]
    valid_f = flux_arr[good_arr]
    if len(valid_t) < 1000:
        return False

    # Check for dense recurring modulation using ratio of global peak-to-trough range to point-to-point scatter
    p95 = np.percentile(valid_f, 95)
    p05 = np.percentile(valid_f, 5)
    f_diff = np.diff(valid_f)
    sigma_diff = (1.4826 / np.sqrt(2)) * np.median(np.abs(f_diff - np.median(f_diff)))
    
    # Continuous stellar variability exhibits large smooth modulation (ratio > 5.0) with many triggers
    if (p95 - p05) > 5.0 * max(sigma_diff, 1e-9) and n_events >= 8:
        return True
    return False


def process_observation(
    fits_path: str,
    threshold_t: float = FROZEN_THRESHOLD_T,
    method: str = "fixed_box_2h"
) -> Dict[str, Any]:
    """Execute complete B1 detector, adaptive clustering, and refined verification pipeline on a FITS file."""
    t_start = time.perf_counter()
    extracted = extract_fits_light_curve(fits_path)
    filename = os.path.basename(fits_path)

    if extracted["status"] != "ok":
        return {
            "status": extracted["status"],
            "tic": extracted["tic"],
            "sector": extracted["sector"],
            "camera": extracted["camera"],
            "ccd": extracted["ccd"],
            "tmag": extracted["tmag"],
            "teff": extracted["teff"],
            "n_cadences": extracted["n_cadences"],
            "n_usable_cadences": extracted["n_usable"],
            "n_windows": 0,
            "n_alert_windows": 0,
            "n_events": 0,
            "n_verified_events": 0,
            "candidates": [],
            "process_sec": time.perf_counter() - t_start
        }

    t_arr = extracted["time"]
    f_arr = extracted["flux"]
    g_arr = extracted["good"]
    cadence = extracted["cadence"]

    window_count = 0
    raw_alerts = []
    
    # 1. Slide 512-cadence windows and evaluate fixed_box_2h matched filter
    for w in iter_fits_windows(t_arr, f_arr, g_arr, cadence, length=512, stride=256, min_good=0.97):
        window_count += 1
        raw_cands, noise = b1.candidate_scores(w["time"], w["flux"], w["good"], cadence, method, None)
        if raw_cands:
            top = raw_cands[0]
            if float(top["score"]) >= threshold_t:
                raw_alerts.append({
                    "center_btjd": float(top["center_btjd"]),
                    "score": float(top["score"]),
                    "width_hr": float(top["width_hr"]),
                    "depth_ppm": float(top["depth_ppm"]),
                    "window_start_btjd": w["window_start_btjd"],
                    "window_end_btjd": w["window_end_btjd"]
                })

    n_alert_windows = len(raw_alerts)

    # 2. Adaptive Duration Event Clustering
    if not raw_alerts:
        return {
            "status": "no_anomalies",
            "tic": extracted["tic"],
            "sector": extracted["sector"],
            "camera": extracted["camera"],
            "ccd": extracted["ccd"],
            "tmag": extracted["tmag"],
            "teff": extracted["teff"],
            "n_cadences": extracted["n_cadences"],
            "n_usable_cadences": extracted["n_usable"],
            "n_windows": window_count,
            "n_alert_windows": 0,
            "n_events": 0,
            "n_verified_events": 0,
            "candidates": [],
            "process_sec": time.perf_counter() - t_start
        }

    raw_alerts_sorted = sorted(raw_alerts, key=lambda x: x["center_btjd"])
    clusters = periodicity.cluster_candidate_events(raw_alerts_sorted)

    candidates = []
    n_verified = 0

    # Check star-level continuous variability
    is_continuous_var = detect_continuous_variability(t_arr, f_arr, g_arr, len(clusters))

    # 3. Two-pass flare-aware photometric verification & refined five-state taxonomy
    for cl in clusters:
        center = float(cl["center_btjd"])
        width = float(cl["width_hr"])
        score = float(cl["score"])
        depth = float(cl["depth_ppm"])
        n_supp = int(cl.get("n_support_windows", 1))

        ver_res = vetting.verify_event_photometry_v2(
            t_arr, f_arr, g_arr,
            center_btjd=center,
            width_hr=width,
            mf_depth_ppm=depth,
            min_consecutive_cadences=2
        )

        raw_out = ver_res["outcome"]
        raw_reason = ver_res["reason"]
        has_flare = ver_res["has_flare"]
        sh_elev = ver_res["shoulder_elev_pct"]
        true_min = ver_res["true_depth_min_pct"]
        mf_depth_pct = depth / 10000.0

        # Assign Refined 5-State Scientific Taxonomy
        # (Separating photometric verification from physical source attribution)
        if raw_out == "inconclusive":
            sci_class = "inconclusive"
            sci_desc = "Data gap or boundary proximity preventing robust baseline"
        elif raw_reason in ["phantom_shoulder_elevation", "flare_shoulder_artifact"] or (
            has_flare and sh_elev > 0.4 * mf_depth_pct and true_min < 0.5 * mf_depth_pct
        ):
            sci_class = "flare_induced_apparent_dip"
            sci_desc = f"Apparent dip induced by elevated flare shoulder (elevation={sh_elev:.2f}%)"
        elif is_continuous_var and raw_out == "probable_phantom":
            sci_class = "continuous_stellar_variability"
            sci_desc = "Wave trough of continuous periodic/rotational stellar modulation"
        elif raw_out == "verified_dimming":
            # Genuine localized dimming (including eclipsing binary eclipses and planetary transits)
            sci_class = "verified_localized_dimming"
            sci_desc = "Genuine localized dimming confirmed with multi-cadence coherence"
            n_verified += 1
        else:
            # Unverified non-flare, non-continuous excursion (Class 4 in 5-state taxonomy)
            sci_class = "probable_instrumental_artifact"
            sci_desc = "Unverified flux excursion without multi-cadence coherence (provisional instrumental / noise artifact)"

        candidates.append({
            "source_file": filename,
            "center_btjd": center,
            "score": score,
            "width_hr": width,
            "depth_ppm": depth,
            "n_support_windows": n_supp,
            "raw_outcome": raw_out,
            "raw_reason": raw_reason,
            "scientific_classification": sci_class,
            "scientific_description": sci_desc + " [post-detection coherence_v2]",
            "true_depth_min_pct": ver_res["true_depth_min_pct"],
            "true_depth_mean_pct": ver_res["true_depth_mean_pct"],
            "consecutive_depressed": ver_res["consecutive_depressed"],
            "has_shoulder_flare": None if has_flare is None else int(has_flare),
            "shoulder_elevation_pct": sh_elev,
            "supporting_detections": cl.get("supporting_detections", [])
        })

    return {
        "status": "completed",
        "tic": extracted["tic"],
        "sector": extracted["sector"],
        "camera": extracted["camera"],
        "ccd": extracted["ccd"],
        "tmag": extracted["tmag"],
        "teff": extracted["teff"],
        "n_cadences": extracted["n_cadences"],
        "n_usable_cadences": extracted["n_usable"],
        "n_windows": window_count,
        "n_alert_windows": n_alert_windows,
        "n_events": len(candidates),
        "n_verified_events": n_verified,
        "candidates": candidates,
        "process_sec": time.perf_counter() - t_start
    }


def validate_fits_identity(path, tic, sector=None):
    """Reject photometry inconsistent with the requested SPOC two-minute product."""
    with fits.open(path, memmap=False) as hdul:
        header = hdul[0].header
        if str(header.get("TICID", "")).strip() != str(tic).replace("TIC:", "").strip():
            raise ValueError("FITS TIC does not match requested target")
        if sector is not None and int(header.get("SECTOR", -1)) != int(sector):
            raise ValueError("FITS sector does not match requested observation")
        times = np.asarray(hdul[1].data["TIME"], dtype=float)
        delta = np.diff(times[np.isfinite(times)])
        cadence = float(np.median(delta[delta > 0])) if np.any(delta > 0) else float('nan')
        if not math.isfinite(cadence) or abs(cadence * 86400 - 120) > 1.2:
            raise ValueError("FITS cadence is not the supported two-minute cadence")


def _worker_process_task(
    item: Dict[str, Any],
    temp_cache: Path,
    session: requests.Session,
    worker_id: int,
    retain_fits: bool,
    check_dev_universe: bool,
    dev_tics: Set[str],
    rate_limit_delay_sec: float,
    threshold_t: float = FROZEN_THRESHOLD_T,
    method: str = "fixed_box_2h"
) -> Dict[str, Any]:
    """Worker task executing download and in-memory processing for a single observation."""
    tic = str(item.get("tic_id") or item.get("tic")).replace("TIC:", "").strip()
    sector = int(item.get("sector")) if item.get("sector") is not None else None
    t_start = time.perf_counter()
    dest_file = None

    if check_dev_universe and tic not in dev_tics:
        return {
            "success": False,
            "tic": tic, "sector": sector,
            "stage": "validation_guard",
            "error": f"DATA LEAKAGE GUARD: TIC {tic} is outside development universe!",
            "wall_sec": time.perf_counter() - t_start
        }

    try:
        if item.get('_resolution_error'):
            raise ValueError(item['_resolution_error'])
        # Resolve product
        dl_url = item.get("source_url") or item.get("download_url")
        fname = item.get("filename") or item.get("product_filename")
        if dl_url and fname:
            download_url = dl_url
            filename = fname
            sec = sector if sector is not None else int(re.search(r"-s(\d{4})-", filename).group(1))
            data_uri = f"mast:TESS/product/{filename}"
            obs_id = str(item.get("obs_id", ""))
            prv_ver = str(item.get("prvversion", ""))
            cal_lvl = int(item.get("calib_level", 2))
        else:
            resolved = resolve_mast_product(tic, sector, astroquery_first=True, session=session)
            download_url = resolved["download_url"]
            filename = resolved["product_filename"]
            sec = resolved["sector"]
            data_uri = resolved["data_uri"]
            obs_id = resolved["obs_id"]
            prv_ver = resolved["prvversion"]
            cal_lvl = resolved["calib_level"]

        if Path(filename).name != filename or '/' in filename or '\\' in filename:
            raise ValueError("Invalid archive product filename")
        # Task-unique destination file
        dest_file = temp_cache / f"w{worker_id}_{filename}"
        n_bytes, dl_sec, sha256_hash = download_light_curve(download_url, str(dest_file), session=session)

        validate_fits_identity(dest_file, tic, sec)
        # Process observation
        res = process_observation(str(dest_file), threshold_t=threshold_t, method=method)

        if rate_limit_delay_sec > 0:
            time.sleep(rate_limit_delay_sec)

        return {
            "success": True,
            "tic": tic,
            "sector": sec,
            "obs_id": obs_id,
            "product_filename": filename,
            "data_uri": data_uri,
            "source_url": download_url,
            "file_sha256": sha256_hash,
            "prvversion": prv_ver,
            "calib_level": cal_lvl,
            "file_bytes": n_bytes,
            "download_sec": dl_sec,
            "dest_file": dest_file,
            "result": res,
            "wall_sec": time.perf_counter() - t_start
        }

    except Exception as e:
        tb = traceback.format_exc()
        if dest_file and dest_file.exists() and not retain_fits:
            try:
                dest_file.unlink()
            except OSError:
                pass
        return {
            "success": False,
            "tic": tic,
            "sector": sector or -1,
            "stage": "pipeline",
            "error": str(e),
            "traceback": tb,
            "wall_sec": time.perf_counter() - t_start
        }


def execution_configuration(threshold_t=FROZEN_THRESHOLD_T, method="fixed_box_2h"):
    effective_config = {
        "length": 512, "stride": 256, "min_good": 0.97,
        "threshold_t": threshold_t, "method": method,
        "quality_mask": TESS_DEFAULT_QUALITY_MASK,
        "b0_hash": core.ORIGINAL_B0_SHA256, "b1_hash": core.ORIGINAL_B1_SHA256,
        "post_detection_method": "coherence_v2",
    }
    effective_config['scientific_source_hashes'] = {
        name: hashlib.sha256(Path(__file__).with_name(name + '.py').read_bytes().replace(b'\r\n', b'\n')).hexdigest()
        for name in ('streaming', 'vetting', 'periodicity', 'core')
    }
    config_json = json.dumps(effective_config, sort_keys=True, separators=(",", ":"))
    # This execution fingerprint covers post-processing too; retain the historical
    # calibration hash separately instead of falsely equating the two contracts.
    config_hash = hashlib.sha256(config_json.encode()).hexdigest()

    return effective_config, config_json, config_hash


def stream_pipeline(
    targets: List[Dict[str, Any]],
    db_path: str,
    run_id: str,
    temp_dir: Optional[str] = None,
    retain_fits: bool = False,
    resume: bool = True,
    check_dev_universe: bool = True,
    dev_universe_manifest: str = "research/phase6b/excluded_development_stars.csv",
    rate_limit_delay_sec: float = 0.02,
    max_observations: int = 0,
    max_workers: int = 2,
    log_every: int = 20,
    threshold_t: float = FROZEN_THRESHOLD_T,
    method: str = "fixed_box_2h",
    progress_callback: Optional[Callable[[Dict[str, Any]], None]] = None,
    cancel_event: Optional[Any] = None,
    notes: Optional[str] = None
) -> Dict[str, Any]:
    """Stream, process, and persist a cohort of TESS light curves with configurable parallelism."""
    core.frozen_check()

    if method != "fixed_box_2h":
        raise ValueError("Only the fixed two-hour detector is supported.")
    if not math.isfinite(threshold_t) or threshold_t <= 0:
        raise ValueError("Threshold must be finite and positive.")
    effective_config, config_json, config_hash = execution_configuration(threshold_t, method)

    dev_tics = set()
    if check_dev_universe:
        dev_tics = load_development_universe(dev_universe_manifest)
        chk_targets = targets[:max_observations] if max_observations > 0 else targets
        for it in chk_targets:
            t = str(it.get('tic_id') or it.get('tic')).replace('TIC:', '').strip()
            check_development_guard(t, dev_tics)

    con = init_streaming_db(db_path)
    cur = con.cursor()

    existing = cur.execute("SELECT settings_json FROM runs WHERE run_id=?", (run_id,)).fetchone()
    if existing:
        if not resume or json.loads(existing[0] or "{}") != effective_config:
            con.close()
            raise ValueError("Existing run configuration differs or lacks exact provenance; create a new run.")

    # Record run entry
    with con:
        con.execute("""
        INSERT OR IGNORE INTO runs (run_id, created_utc, software_version, b0_hash, b1_hash, threshold_t, method, config_hash, worker_count, settings_json, notes)
        VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
        (run_id, utc_now(), __version__, core.ORIGINAL_B0_SHA256, core.ORIGINAL_B1_SHA256,
         threshold_t, method, config_hash, max_workers,
         config_json,
         notes or "Streaming survey run"))

    temp_cache = Path(temp_dir) if temp_dir else Path(tempfile.mkdtemp(prefix="tess_stream_cache_"))
    temp_cache.mkdir(parents=True, exist_ok=True)

    completed_obs = set()
    completed_tic_sec = set()
    if resume:
        for r in cur.execute("SELECT tic, sector, product_filename, status FROM observations WHERE run_id=?", (run_id,)).fetchall():
            if r[3] in ('completed', 'no_anomalies', 'completed_with_anomalies'):
                tic_clean = str(r[0]).replace('TIC:', '').strip()
                completed_obs.add((str(r[0]), int(r[1]), str(r[2])))
                completed_tic_sec.add((f"TIC:{tic_clean}", int(r[1])))
                completed_tic_sec.add((tic_clean, int(r[1])))

    # Resolve identities before scheduling, then deduplicate exact archive products.
    pending_targets = []
    seen_products = set()
    seen_requests = set()
    n_skipped = 0
    for original in targets:
        if cancel_event and cancel_event.is_set():
            break
        item = dict(original)
        clean_t = str(item.get('tic_id') or item.get('tic')).replace('TIC:', '').strip()
        sec = int(item['sector']) if item.get('sector') is not None else None
        fname = item.get('filename') or item.get('product_filename')
        request_key = (clean_t, sec, fname)
        if request_key in seen_requests:
            n_skipped += 1
            continue
        seen_requests.add(request_key)
        if not fname or not (item.get('source_url') or item.get('download_url')):
            try:
                resolved = resolve_mast_product(clean_t, sec)
                if fname and fname != resolved['product_filename']:
                    raise ValueError('Requested product differs from resolved archive product')
                item.update(resolved)
                item['filename'] = resolved['product_filename']
                fname, sec = resolved['product_filename'], resolved['sector']
            except Exception as exc:
                item['_resolution_error'] = str(exc)
        key = (f'TIC:{clean_t}', sec, fname)
        if key in seen_products or (resume and key in completed_obs):
            n_skipped += 1
            continue
        seen_products.add(key)
        pending_targets.append(item)
        if max_observations > 0 and len(pending_targets) >= max_observations:
            break

    total_download_bytes = 0
    total_download_time = 0.0
    total_process_time = 0.0
    wall_start = time.perf_counter()

    n_processed = 0
    n_successful = 0
    n_failed = 0
    retry_count = 0

    proc = psutil.Process() if psutil else None
    peak_ram_mb = 0.0
    cpu_measurements = []

    # Thread pool and shared HTTP sessions
    sessions = [requests.Session() for _ in range(max(1, max_workers))]
    for s in sessions:
        ad = requests.adapters.HTTPAdapter(pool_connections=10, pool_maxsize=10, max_retries=3)
        s.mount("https://", ad)

    try:
        # Bounded thread pool execution
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            def completed_jobs():
                remaining = iter(enumerate(pending_targets))
                active = set()
                def fill():
                    while len(active) < max_workers and not (cancel_event and cancel_event.is_set()):
                        pair = next(remaining, None)
                        if pair is None:
                            break
                        idx, item = pair
                        active.add(executor.submit(
                            _worker_process_task, item, temp_cache, None, idx,
                            retain_fits, check_dev_universe, dev_tics,
                            rate_limit_delay_sec, threshold_t, method))
                fill()
                while active:
                    done, active = concurrent.futures.wait(active, return_when=concurrent.futures.FIRST_COMPLETED)
                    for future in done:
                        yield future
                    # Completed observations are committed before new work is submitted.
                    fill()
            for fut in completed_jobs():
                n_processed += 1
                out = fut.result()
                dest_file = out.get("dest_file")

                try:
                    if out["success"]:
                        res = out["result"]
                        n_bytes = out["file_bytes"]
                        dl_sec = out["download_sec"]
                        total_download_bytes += n_bytes
                        total_download_time += dl_sec
                        total_process_time += res["process_sec"]

                        tic_tag = f"TIC:{out['tic']}"
                        sec = out["sector"]

                        with con:
                            cur.execute("DELETE FROM observations WHERE run_id=? AND tic=? AND sector=? AND status='failed'", (run_id, tic_tag, sec))
                            cur.execute("""
                            INSERT INTO observations (
                                run_id, tic, sector, camera, ccd, tmag, teff, obs_id, product_filename,
                                data_uri, source_url, file_sha256, calib_level, prvversion, config_hash,
                                status, n_cadences, n_usable_cadences, n_windows, n_alert_windows, n_events, n_verified_events,
                                download_time_sec, process_time_sec, total_time_sec, file_bytes, error_message, created_utc
                            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                            (run_id, tic_tag, sec, res["camera"], res["ccd"], res["tmag"], res["teff"],
                             out["obs_id"], out["product_filename"], out["data_uri"], out["source_url"],
                             out["file_sha256"], out["calib_level"], out["prvversion"], config_hash,
                             res["status"], res["n_cadences"], res["n_usable_cadences"], res["n_windows"],
                             res["n_alert_windows"], res["n_events"], res["n_verified_events"],
                             dl_sec, res["process_sec"], out["wall_sec"], n_bytes, None, utc_now()))

                            obs_id_db = cur.lastrowid

                            for c in res["candidates"]:
                                cur.execute("""
                                INSERT INTO candidates (
                                    observation_id, run_id, tic, sector, source_file, center_btjd, score, width_hr, depth_ppm,
                                    n_support_windows, raw_outcome, raw_reason, scientific_classification,
                                    scientific_description, true_depth_min_pct, true_depth_mean_pct,
                                    consecutive_depressed, has_shoulder_flare, shoulder_elevation_pct,
                                    review_status, note, created_utc
                                ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                                (obs_id_db, run_id, tic_tag, sec, c["source_file"], c["center_btjd"], c["score"],
                                 c["width_hr"], c["depth_ppm"], c["n_support_windows"], c["raw_outcome"],
                                 c["raw_reason"], c["scientific_classification"], c["scientific_description"],
                                 c["true_depth_min_pct"], c["true_depth_mean_pct"], c["consecutive_depressed"],
                                 c["has_shoulder_flare"], c["shoulder_elevation_pct"], "unreviewed", "", utc_now()))

                                cand_id_db = cur.lastrowid

                                for det in c.get("supporting_detections", []):
                                    cur.execute("""
                                    INSERT INTO candidate_detections (
                                        candidate_id, run_id, window_start_btjd, window_end_btjd, center_btjd,
                                        score, width_hr, depth_ppm, created_utc
                                    ) VALUES (?,?,?,?,?,?,?,?,?)""",
                                    (cand_id_db, run_id, det["window_start_btjd"], det["window_end_btjd"],
                                     det["center_btjd"], det["score"], det["width_hr"], det["depth_ppm"], utc_now()))

                        n_successful += 1

                    else: # Failure case
                        n_failed += 1
                        tic_tag = f"TIC:{out['tic']}"
                        with con:
                            cur.execute("""
                            INSERT INTO failures (run_id, tic, sector, stage, error_message, traceback, created_utc)
                            VALUES (?,?,?,?,?,?,?)""",
                            (run_id, tic_tag, out.get("sector"), out.get("stage", "pipeline"),
                             str(out.get("error")), out.get("traceback"), utc_now()))

                            cur.execute("""
                            INSERT OR REPLACE INTO observations (
                                run_id, tic, sector, camera, ccd, tmag, teff, obs_id, product_filename,
                                data_uri, source_url, file_sha256, calib_level, prvversion, config_hash,
                                status, n_cadences, n_usable_cadences, n_windows, n_alert_windows, n_events, n_verified_events,
                                download_time_sec, process_time_sec, total_time_sec, file_bytes, error_message, created_utc
                            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                            (run_id, tic_tag, out.get("sector", -1), 0, 0, math.nan, math.nan, "",
                             dest_file.name if dest_file else "unknown", None, None, None, 2, "", config_hash,
                             "failed", 0, 0, 0, 0, 0, 0, 0.0, 0.0, out["wall_sec"], 0, str(out.get("error")), utc_now()))

                finally:
                    if dest_file and dest_file.exists() and not retain_fits:
                        try:
                            dest_file.unlink()
                        except OSError:
                            pass

                # Progress & Resource Tracking
                rss_mb = (proc.memory_info().rss / (1024 * 1024)) if proc else 0.0
                if rss_mb > peak_ram_mb:
                    peak_ram_mb = rss_mb
                if proc:
                    try:
                        cpu_measurements.append(proc.cpu_percent())
                    except Exception:
                        pass

                wall_el = time.perf_counter() - wall_start
                rate_per_hr = (n_processed / wall_el) * 3600.0 if wall_el > 0 else 0.0
                eta_sec = ((len(pending_targets) - n_processed) / (n_processed / wall_el)) if (n_processed > 0 and wall_el > 0 and len(pending_targets) > n_processed) else 0.0

                if progress_callback:
                    try:
                        progress_callback({
                            "run_id": run_id,
                            "processed": n_processed,
                            "total": len(pending_targets),
                            "successful": n_successful,
                            "failed": n_failed,
                            "latest_tic": out.get("tic"),
                            "latest_success": out.get("success"),
                            "rate_per_hr": rate_per_hr,
                            "elapsed_sec": wall_el,
                            "eta_sec": eta_sec,
                            "peak_ram_mb": peak_ram_mb,
                            "error": out.get("error") if not out.get("success") else None
                        })
                    except Exception:
                        pass

                if n_processed % log_every == 0 or n_processed == len(pending_targets):
                    print(f"[{utc_now()[:19]}] [Workers={max_workers}] Processed {n_processed}/{len(pending_targets)} obs | "
                          f"Success: {n_successful} | Fail: {n_failed} | Rate: {rate_per_hr:.1f} obs/hr | RAM: {rss_mb:.1f} MB")

    finally:
        total_wall_sec = time.perf_counter() - wall_start
        avg_rate_hr = (n_successful / total_wall_sec * 3600.0) if total_wall_sec > 0 else 0.0
        avg_mb_sec = (total_download_bytes / (1024 * 1024) / total_download_time) if total_download_time > 0 else 0.0
        avg_cpu = float(np.mean(cpu_measurements)) if cpu_measurements else 0.0

        with con:
            con.execute("""
            INSERT OR REPLACE INTO performance_benchmarks (
                run_id, worker_count, total_observations, successful_observations, failed_observations,
                total_download_bytes, total_download_sec, total_process_sec, total_wall_sec,
                avg_obs_per_hour, avg_download_mb_per_sec, peak_ram_mb, avg_cpu_pct, retry_count, created_utc
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (run_id, max_workers, n_processed, n_successful, n_failed, total_download_bytes,
             total_download_time, total_process_time, total_wall_sec,
             avg_rate_hr, avg_mb_sec, peak_ram_mb, avg_cpu, retry_count, utc_now()))

        con.close()

        if not temp_dir and not retain_fits and temp_cache.exists():
            try:
                shutil.rmtree(temp_cache, ignore_errors=True)
            except OSError:
                pass

    return {
        "run_id": run_id,
        "worker_count": max_workers,
        "n_processed": n_processed,
        "n_successful": n_successful,
        "n_failed": n_failed,
        "n_skipped": n_skipped,
        "cancelled": bool(cancel_event and cancel_event.is_set()),
        "total_download_mb": total_download_bytes / (1024 * 1024),
        "total_download_sec": total_download_time,
        "total_process_sec": total_process_time,
        "total_wall_sec": total_wall_sec,
        "avg_obs_per_hour": avg_rate_hr,
        "avg_download_mb_sec": avg_mb_sec,
        "peak_ram_mb": peak_ram_mb,
        "avg_cpu_pct": avg_cpu,
        "retry_count": retry_count
    }
