"""Star-centered research investigation workspace and web review interface.

Provides a unified star-first research platform for TESS photometric anomaly exploration,
supporting both legacy and streaming Phase 7 / Phase 6B survey databases, on-demand
NASA MAST light curve streaming, scope-isolated researcher reviews, safe undo history,
and interactive Plotly.js / Matplotlib visualizations.
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import html
from http import HTTPStatus
from http.server import HTTPServer, ThreadingHTTPServer, BaseHTTPRequestHandler
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import secrets
import sqlite3
import sys
import threading
import time
import urllib.parse
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from . import __version__, core, db, experiments, periodicity, streaming

VALID_STATUSES = ('unreviewed', 'interesting', 'known', 'artifact', 'dismissed')
VALID_STAR_CLASSES = (
    'candidate_transit',
    'eclipsing_binary',
    'pulsator',
    'disk_occultation',
    'stellar_variability',
    'flare_star',
    'instrumental_artifact',
    'unresolved'
)

# In-memory LRU light curve cache: (tic, sector) -> {time, flux_ppm, good, cadence, source}
_LIGHT_CURVE_CACHE: Dict[Tuple[str, int], Dict[str, Any]] = {}
_CACHE_LOCK = threading.Lock()
MAX_IN_MEMORY_LIGHT_CURVES = 64


def get_default_cache_dir() -> Path:
    """Return default local cache directory for on-demand downloaded FITS products."""
    base = Path.home() / ".tess_anomaly_survey" / "cache"
    base.mkdir(parents=True, exist_ok=True)
    return base


def sanitize_tic(tic: Any) -> str:
    """Standardize TIC identifier to clean string."""
    raw = str(tic).strip()
    match = re.fullmatch(r"(?:TIC:?\s*)?(\d+)", raw, flags=re.IGNORECASE)
    return match.group(1) if match else raw


def resolve_light_curve_series(
    tic: str,
    sector: Optional[int] = None,
    source_file: Optional[str] = None,
    product_filename: Optional[str] = None,
    data_uri: Optional[str] = None,
    root_path: Optional[str] = None,
    cache_dir: Optional[str] = None,
    allow_download: bool = True,
    expected_sha256: Optional[str] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, str]:
    """Retrieve time, relative flux in PPM, good mask, and cadence for a target star.

    Attempts resolution via:
    1. In-memory cache
    2. Local root directory (preprocessed NPZ or raw FITS)
    3. Explicit local cache directory
    4. On-demand streaming download from NASA MAST archive
    """
    clean_tic = sanitize_tic(tic)
    local = None
    if source_file:
        proposed = (Path(root_path) / source_file) if root_path else Path(source_file)
        if proposed.is_file():
            local = proposed.resolve()
    fingerprint = (str(local), local.stat().st_size, local.stat().st_mtime_ns) if local else None
    cache_key = (clean_tic, sector, source_file, product_filename, data_uri, fingerprint, expected_sha256)
    with _CACHE_LOCK:
        cached = _LIGHT_CURVE_CACHE.get(cache_key)
        if cached:
            return cached['time'], cached['flux_ppm'], cached['good'], cached['cadence'], cached['source']
    if local:
        if expected_sha256 and hashlib.sha256(local.read_bytes()).hexdigest() != expected_sha256:
            raise ValueError("Local product checksum differs from the recorded observation")
        if local.name.lower().endswith(('.fits', '.fit', '.fits.gz')):
            streaming.validate_fits_identity(local, clean_tic, sector)
        t, f, good, cad = _load_local_light_curve(local)
        return _cache_and_return(cache_key, t, f, good, cad, f"local:{local.name}")
    # Missing synthetic/local data must never be replaced with real observations.
    if source_file and str(source_file).lower().endswith('.npz'):
        raise FileNotFoundError("Synthetic/local NPZ source unavailable; no archive substitution is permitted.")
    if not product_filename and data_uri:
        product_filename = data_uri.rsplit('/', 1)[-1]
    if not product_filename:
        raise FileNotFoundError("Exact archive product identity is unavailable.")
    if Path(product_filename).name != product_filename or '/' in product_filename or '\\' in product_filename:
        raise ValueError("Invalid archive product filename.")
    expected_uri = f"mast:TESS/product/{product_filename}"
    if data_uri and data_uri != expected_uri:
        raise ValueError("Archive URI and product filename disagree.")
    cache_path = Path(cache_dir) if cache_dir else get_default_cache_dir()
    cache_path.mkdir(parents=True, exist_ok=True)
    dest = cache_path / product_filename
    if not dest.is_file():
        if not allow_download:
            raise FileNotFoundError("Recorded archive product is not cached.")
        url = "https://mast.stsci.edu/api/v0.1/Download/file?uri=" + urllib.parse.quote(expected_uri, safe='')
        streaming.download_light_curve(url, str(dest), timeout=35.0, max_retries=2)
    if expected_sha256 and hashlib.sha256(dest.read_bytes()).hexdigest() != expected_sha256:
        raise ValueError("Archive product checksum differs from the recorded observation")
    streaming.validate_fits_identity(dest, clean_tic, sector)
    t, f, good, cad = _load_local_light_curve(dest)
    return _cache_and_return(cache_key, t, f, good, cad, f"mast:{product_filename}")


def _load_local_light_curve(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Parse time series from NPZ or FITS file format into (time, flux_ppm, good, cadence)."""
    p_str = str(path)
    if p_str.endswith(".npz"):
        t, f, good, cad = core.series(path)
        return t, f * 1e6, good, cad

    if p_str.endswith(".fits") or p_str.endswith(".fit") or p_str.endswith(".fits.gz"):
        res = streaming.extract_fits_light_curve(p_str)
        if res.get("status") in ("insufficient_usable_cadences", "unusable_flux"):
            raise ValueError(f"FITS file has insufficient usable cadences: {res.get('status')}")
        t = res["time"]
        f = res["flux"] * 1e6  # Convert relative fraction to PPM
        good = res["good"]
        cad = res["cadence"]
        return t, f, good, cad

    raise ValueError(f"Unsupported light curve file format: {path.name}")


def _cache_and_return(key: Tuple[str, int], t: np.ndarray, f_ppm: np.ndarray, good: np.ndarray, cad: float, source: str):
    with _CACHE_LOCK:
        if len(_LIGHT_CURVE_CACHE) >= MAX_IN_MEMORY_LIGHT_CURVES:
            # Evict oldest entry
            oldest_k = next(iter(_LIGHT_CURVE_CACHE))
            del _LIGHT_CURVE_CACHE[oldest_k]
        _LIGHT_CURVE_CACHE[key] = {
            'time': t, 'flux_ppm': f_ppm, 'good': good, 'cadence': cad, 'source': source
        }
    return t, f_ppm, good, cad, source


# ---------------------------------------------------------------------------
# Database Schema Introspection & Unified Adapter
# ---------------------------------------------------------------------------

class SurveyDatabaseAdapter:
    """Unified adapter wrapping legacy, streaming, and overlay review databases."""

    def __init__(self, db_path: str, review_db_path: Optional[str] = None,
                 read_only: bool = False, run_id: Optional[str] = None):
        self.db_path = str(Path(db_path).resolve())
        self.read_only = read_only
        with sqlite3.connect(Path(self.db_path).as_uri() + "?mode=ro", uri=True) as cx:
            tables = {r[0] for r in cx.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            ids = set()
            for table in ('runs', 'candidates', 'observations'):
                if table in tables:
                    cols = {r[1] for r in cx.execute(f"PRAGMA table_info({table})")}
                    key = 'run_id' if 'run_id' in cols else ('id' if table == 'runs' else None)
                    if key:
                        ids.update(str(r[0]) for r in cx.execute(f"SELECT DISTINCT {key} FROM {table}") if r[0] is not None)
        if run_id is None and len(ids) > 1:
            raise ValueError("Multiple source runs: select an explicit run_id.")
        if run_id is not None and ids and str(run_id) not in ids:
            raise ValueError("Requested run is absent from the source database.")
        self.run_id = str(run_id) if run_id is not None else next(iter(ids), None)
        self._multiple_runs = len(ids) > 1
        self.data_kind = 'synthetic demonstration' if self.run_id == 'sample_survey_run' else 'observational / consult source provenance'
        suffix = hashlib.sha256((self.run_id or 'legacy').encode()).hexdigest()[:16]
        self.review_db_path = str(Path(review_db_path).resolve()) if review_db_path else str(
            Path(self.db_path).with_suffix(f".{suffix}.reviews.sqlite"))
        if self.review_db_path == self.db_path:
            raise ValueError("Review overlay must be separate from scientific results.")
        if Path(self.review_db_path).is_file():
            with sqlite3.connect(Path(self.review_db_path).as_uri() + '?mode=ro', uri=True) as cx:
                tables = {r[0] for r in cx.execute("SELECT name FROM sqlite_master WHERE type='table'")}
                if tables & {'candidates', 'observations', 'runs'}:
                    raise ValueError("Scientific database cannot be used as a review overlay.")
        self._init_overlay_db(self.review_db_path)
        with sqlite3.connect(self.review_db_path) as cx:
            cx.execute("CREATE TABLE IF NOT EXISTS source_binding (source_path TEXT, run_id TEXT)")
            binding = (self.db_path, self.run_id or 'legacy')
            row = cx.execute("SELECT source_path, run_id FROM source_binding").fetchone()
            if row and tuple(row) != binding:
                raise ValueError("Review overlay belongs to another source or run.")
            if not row:
                cx.execute("INSERT INTO source_binding VALUES (?,?)", binding)
        cx = self.get_connection()
        try:
            self.schema_info = self._inspect_schema(cx)
        finally:
            cx.close()

    def _init_overlay_db(self, path: str):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        con = sqlite3.connect(path)
        with con:
            con.executescript('''
            PRAGMA foreign_keys=ON;
            CREATE TABLE IF NOT EXISTS star_reviews (
                tic TEXT PRIMARY KEY,
                status TEXT NOT NULL DEFAULT 'unreviewed',
                star_class TEXT NOT NULL DEFAULT '',
                note TEXT NOT NULL DEFAULT '',
                reviewer TEXT NOT NULL DEFAULT 'local_user',
                updated_utc TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS star_review_history (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                tic TEXT NOT NULL,
                previous_status TEXT NOT NULL,
                new_status TEXT NOT NULL,
                note TEXT NOT NULL,
                reviewer TEXT NOT NULL DEFAULT 'local_user',
                created_utc TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS candidate_review_overlay (
                candidate_id INTEGER PRIMARY KEY,
                status TEXT NOT NULL,
                note TEXT NOT NULL DEFAULT '',
                reviewer TEXT NOT NULL DEFAULT 'local_user',
                updated_utc TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS review_history (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                candidate_id INTEGER NOT NULL,
                previous_status TEXT NOT NULL,
                new_status TEXT NOT NULL,
                note TEXT NOT NULL,
                reviewer TEXT NOT NULL DEFAULT 'local_user',
                created_utc TEXT NOT NULL
            );
            ''')
        con.close()

    def get_connection(self) -> sqlite3.Connection:
        # The source is never opened writable, including when annotations are enabled.
        cx = sqlite3.connect(Path(self.db_path).as_uri() + "?mode=ro", uri=True, timeout=30.0)
        cx.create_function("selected_run", 0, lambda: self.run_id)
        tables = {r[0] for r in cx.execute("SELECT name FROM main.sqlite_master WHERE type='table'")}
        for table in ('candidates', 'observations', 'candidate_systems', 'candidate_detections',
                      'processed_files', 'performance_benchmarks'):
            if table in tables:
                cols = {r[1] for r in cx.execute(f"PRAGMA main.table_info({table})")}
                if 'run_id' in cols and self.run_id is not None:
                    cx.execute(f"CREATE TEMP VIEW {table} AS SELECT * FROM main.{table} WHERE CAST(run_id AS TEXT)=selected_run()")
        if 'star_reviews' not in tables or self._multiple_runs:
            cx.execute("CREATE TEMP VIEW star_reviews AS SELECT NULL AS tic, NULL AS status, NULL AS star_class, NULL AS note WHERE 0")
        if 'star_review_history' not in tables or self._multiple_runs:
            cx.execute("CREATE TEMP VIEW star_review_history AS SELECT NULL AS id, NULL AS tic, NULL AS previous_status, NULL AS new_status, NULL AS note, NULL AS reviewer, NULL AS created_utc WHERE 0")
        cx.execute("ATTACH DATABASE ? AS overlay", (self.review_db_path,))
        return cx

    def _inspect_schema(self, cx: sqlite3.Connection) -> Dict[str, Any]:
        tables = {r[0] for r in cx.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
        cand_cols = {c[1] for c in cx.execute("PRAGMA table_info(candidates)").fetchall()} if 'candidates' in tables else set()
        obs_cols = {c[1] for c in cx.execute("PRAGMA table_info(observations)").fetchall()} if 'observations' in tables else set()

        return {
            'has_observations': 'observations' in tables,
            'has_processed_files': 'processed_files' in tables,
            'has_candidate_systems': 'candidate_systems' in tables,
            'has_detections': 'candidate_detections' in tables,
            'has_benchmarks': 'performance_benchmarks' in tables,
            'has_scientific_classification': 'scientific_classification' in cand_cols,
            'has_coincidence_col': 'coincidence_risk_level' in cand_cols or 'coincidence_star_count' in cand_cols,
            'has_overlay': self.review_db_path is not None,
            'candidates_cols': cand_cols,
            'obs_cols': obs_cols
        }

    def get_survey_stats(self) -> Dict[str, Any]:
        cx = self.get_connection()
        try:
            total_cands = cx.execute("SELECT COUNT(*) FROM candidates").fetchone()[0]

            # Star aggregates
            if self.review_db_path:
                cur = cx.execute('''
                    SELECT c.tic,
                           COUNT(c.id) as n_cands,
                           SUM(CASE WHEN COALESCE(o_rev.status, c.review_status, 'unreviewed') != 'unreviewed' THEN 1 ELSE 0 END) as n_rev,
                           COALESCE(sr_ov.status, sr.status, 'unreviewed') as star_st
                    FROM candidates c
                    LEFT JOIN star_reviews sr ON c.tic = sr.tic
                    LEFT JOIN overlay.candidate_review_overlay o_rev ON c.id = o_rev.candidate_id
                    LEFT JOIN overlay.star_reviews sr_ov ON c.tic = sr_ov.tic
                    GROUP BY c.tic
                ''')
            else:
                cur = cx.execute('''
                    SELECT c.tic,
                           COUNT(c.id) as n_cands,
                           SUM(CASE WHEN c.review_status != 'unreviewed' THEN 1 ELSE 0 END) as n_rev,
                           COALESCE(sr.status, 'unreviewed') as star_st
                    FROM candidates c
                    LEFT JOIN star_reviews sr ON c.tic = sr.tic
                    GROUP BY c.tic
                ''')

            star_rows = cur.fetchall()
            total_stars = len(star_rows)
            fully_reviewed_stars = sum(1 for r in star_rows if r[1] > 0 and r[2] == r[1])
            partially_reviewed_stars = sum(1 for r in star_rows if 0 < r[2] < r[1])
            unreviewed_stars = sum(1 for r in star_rows if r[2] == 0)

            # Verification breakdown if present
            verif_counts = defaultdict(int)
            if self.schema_info['has_scientific_classification']:
                for v_row in cx.execute("SELECT scientific_classification, COUNT(*) FROM candidates GROUP BY scientific_classification"):
                    verif_counts[v_row[0]] = v_row[1]
            else:
                verif_counts['unclassified'] = total_cands

            # Periodic systems
            periodic_stars = 0
            if self.schema_info['has_candidate_systems']:
                periodic_stars = cx.execute("SELECT COUNT(DISTINCT tic) FROM candidate_systems WHERE is_periodic=1").fetchone()[0]

            # Event review status counts
            event_counts = {s: 0 for s in VALID_STATUSES}
            if self.review_db_path:
                ev_cur = cx.execute('''
                    SELECT COALESCE(o_rev.status, c.review_status, 'unreviewed'), COUNT(*)
                    FROM candidates c
                    LEFT JOIN overlay.candidate_review_overlay o_rev ON c.id = o_rev.candidate_id
                    GROUP BY COALESCE(o_rev.status, c.review_status, 'unreviewed')
                ''')
            else:
                ev_cur = cx.execute('''
                    SELECT c.review_status, COUNT(*)
                    FROM candidates c
                    GROUP BY c.review_status
                ''')
            for r in ev_cur:
                st = r[0]
                if st in event_counts:
                    event_counts[st] = r[1]

            return {
                'events': {
                    'total': total_cands,
                    'reviewed': total_cands - event_counts['unreviewed'],
                    'unreviewed': event_counts['unreviewed'],
                    'breakdown': event_counts
                },
                'stars': {
                    'total': total_stars,
                    'fully_reviewed': fully_reviewed_stars,
                    'partially_reviewed': partially_reviewed_stars,
                    'unreviewed': unreviewed_stars,
                    'periodic_count': periodic_stars
                },
                'verifications': dict(verif_counts),
                'read_only': self.read_only,
                'has_overlay': self.review_db_path is not None,
                'data_kind': self.data_kind,
                'db_name': Path(self.db_path).name
            }
        finally:
            cx.close()

    def query_stars(
        self,
        page: int = 1,
        limit: int = 25,
        status_filter: str = 'all',
        verification_filter: str = 'all',
        system_filter: str = 'all',
        coincidence_filter: str = 'all',
        search_query: str = '',
        sort_by: str = 'score',
        sort_order: str = 'desc'
    ) -> Dict[str, Any]:
        """Star-centered catalog query with expandable constituent events and multi-domain filtering."""
        cx = self.get_connection()
        try:
            # Build unified candidate rows with overlay review values
            ov_cand_join = "LEFT JOIN overlay.candidate_review_overlay o_rev ON c.id = o_rev.candidate_id" if self.review_db_path else ""
            ov_star_join = "LEFT JOIN overlay.star_reviews sr_ov ON c.tic = sr_ov.tic" if self.review_db_path else ""

            verif_col = "c.scientific_classification" if self.schema_info['has_scientific_classification'] else "'not_run'"

            if self.schema_info['has_candidate_systems']:
                sys_cols = """
                    s.period_days as period_days,
                    s.is_periodic as is_periodic,
                    s.classification as system_classification,
                    s.epoch_btjd as epoch_btjd,
                    s.odd_even_sigma as odd_even_sigma,
                    s.period_uncertainty_days as period_uncertainty_days,
                    s.alternate_period_days as alternate_period_days
                """
                sys_joins = """
                    LEFT JOIN candidate_systems s ON c.system_id = s.id

                """
            else:
                sys_cols = """
                    NULL as period_days,
                    0 as is_periodic,
                    '' as system_classification,
                    NULL as epoch_btjd,
                    0.0 as odd_even_sigma,
                    NULL as period_uncertainty_days,
                    NULL as alternate_period_days
                """
                sys_joins = ""

            sql = f'''
                SELECT c.id, c.tic, c.sector, c.source_file, c.center_btjd, c.score, c.width_hr, c.depth_ppm,
                       COALESCE({'o_rev.status, ' if self.review_db_path else ''}c.review_status, 'unreviewed') as eff_event_status,
                       COALESCE({'o_rev.note, ' if self.review_db_path else ''}c.note, '') as eff_event_note,
                       {verif_col} as verif_class,
                       COALESCE({'sr_ov.status, ' if self.review_db_path else ''}sr.status, 'unreviewed') as eff_star_status,
                       COALESCE({'sr_ov.star_class, ' if self.review_db_path else ''}sr.star_class, '') as eff_star_class,
                       COALESCE({'sr_ov.note, ' if self.review_db_path else ''}sr.note, '') as eff_star_note,
                       {sys_cols}
                FROM candidates c
                LEFT JOIN star_reviews sr ON c.tic = sr.tic
                {sys_joins}
                {ov_cand_join}
                {ov_star_join}
            '''

            rows = cx.execute(sql).fetchall()

            # Group into Star dictionary
            stars_map: Dict[str, Dict[str, Any]] = {}
            for r in rows:
                (cid, tic, sec, src_file, center_t, score, width, depth,
                 ev_st, ev_note, verif, star_st, star_cls, star_note,
                 p_days, is_per, sys_cls, epoch_t, oe_sig, p_err, alt_p) = r

                if tic not in stars_map:
                    stars_map[tic] = {
                        'tic': str(tic),
                        'clean_tic': sanitize_tic(tic),
                        'sectors': set(),
                        'events': [],
                        'max_score': -1.0,
                        'star_status': star_st,
                        'star_class': star_cls,
                        'star_note': star_note,
                        'period_days': p_days,
                        'is_periodic': bool(is_per) if is_per is not None else False,
                        'system_classification': sys_cls or '',
                        'epoch_btjd': epoch_t,
                        'odd_even_sigma': oe_sig or 0.0,
                        'period_uncertainty_days': p_err,
                        'alternate_period_days': alt_p,
                        'primary_verification': verif or 'not_run',
                        'coincidence_risk': 'unknown',  # Will refine if systematic
                        'has_dimming': False
                    }

                st_obj = stars_map[tic]
                if sec is not None:
                    st_obj['sectors'].add(int(sec))
                if score > st_obj['max_score']:
                    st_obj['max_score'] = float(score)

                if verif == 'verified_localized_dimming':
                    st_obj['has_dimming'] = True
                    st_obj['primary_verification'] = 'verified_localized_dimming'

                st_obj['events'].append({
                    'id': cid,
                    'tic': str(tic),
                    'sector': sec,
                    'source_file': src_file,
                    'center_btjd': round(center_t, 4),
                    'score': round(float(score), 2),
                    'width_hr': round(float(width), 2),
                    'depth_ppm': round(float(depth), 1),
                    'scientific_classification': verif,
                    'review_status': ev_st,
                    'note': ev_note,
                    'period_days': p_days,
                    'is_periodic': bool(is_per) if is_per is not None else False,
                    'system_classification': sys_cls or '',
                    'period_uncertainty_days': p_err,
                    'alternate_period_days': alt_p
                })

            # Format star list
            star_list = []
            for tic, sobj in stars_map.items():
                sobj['sectors'] = sorted(list(sobj['sectors']))
                sobj['sector_display'] = ", ".join(str(s) for s in sobj['sectors']) if sobj['sectors'] else "—"
                sobj['total_events'] = len(sobj['events'])
                sobj['reviewed_events'] = sum(1 for e in sobj['events'] if e['review_status'] != 'unreviewed')

                # Calculate progress
                if sobj['reviewed_events'] == 0:
                    sobj['progress'] = 'unreviewed'
                elif sobj['reviewed_events'] == sobj['total_events']:
                    sobj['progress'] = 'reviewed'
                else:
                    sobj['progress'] = 'partially_reviewed'

                # Sort events within star by score DESC
                sobj['events'].sort(key=lambda e: e['score'], reverse=True)
                star_list.append(sobj)

            # Apply Filters
            filtered = []
            q_clean = search_query.strip().lower()

            for s in star_list:
                # 1. Search Query
                if q_clean:
                    match_tic = q_clean in s['tic'].lower()
                    match_note = q_clean in s['star_note'].lower() or any(q_clean in e['note'].lower() for e in s['events'])
                    match_class = q_clean in s['star_class'].lower()
                    match_sec = any(q_clean == str(sec) for sec in s['sectors'])
                    if not (match_tic or match_note or match_class or match_sec):
                        continue

                # 2. Status Filter
                if status_filter != 'all':
                    if status_filter in VALID_STATUSES:
                        if s['star_status'] != status_filter:
                            continue
                    elif status_filter == 'unreviewed_progress' and s['progress'] != 'unreviewed':
                        continue

                # 3. Verification Filter
                if verification_filter != 'all':
                    if verification_filter == 'not_run' and not any(e['scientific_classification'] in (None, 'not_run', 'unknown', 'unavailable') for e in s['events']):
                        continue
                    if verification_filter == 'dimming' and not any(e['scientific_classification'] == 'verified_localized_dimming' for e in s['events']):
                        continue
                    elif verification_filter == 'variability' and not any(e['scientific_classification'] == 'continuous_stellar_variability' for e in s['events']):
                        continue
                    elif verification_filter == 'flare' and not any(e['scientific_classification'] == 'flare_induced_apparent_dip' for e in s['events']):
                        continue
                    elif verification_filter == 'artifact' and not any(e['scientific_classification'] == 'probable_instrumental_artifact' for e in s['events']):
                        continue

                # 4. System Filter
                if system_filter == 'periodic' and not s['is_periodic']:
                    continue
                elif system_filter == 'isolated' and s['is_periodic']:
                    continue

                # 5. Coincidence Filter
                if coincidence_filter != 'all':
                    if coincidence_filter == 'none' and s['coincidence_risk'] != 'none':
                        continue
                    elif coincidence_filter == 'high' and s['coincidence_risk'] != 'high':
                        continue

                filtered.append(s)

            # Apply Sorting
            reverse = (sort_order.lower() == 'desc')
            if sort_by == 'score':
                filtered.sort(key=lambda s: s['max_score'], reverse=reverse)
            elif sort_by == 'events':
                filtered.sort(key=lambda s: s['total_events'], reverse=reverse)
            elif sort_by == 'sector':
                filtered.sort(key=lambda s: (s['sectors'][0] if s['sectors'] else 0), reverse=reverse)
            elif sort_by == 'tic':
                filtered.sort(key=lambda s: s['clean_tic'], reverse=reverse)
            elif sort_by == 'period':
                filtered.sort(key=lambda s: (s['period_days'] or 0.0), reverse=reverse)
            elif sort_by == 'status':
                filtered.sort(key=lambda s: s['star_status'], reverse=reverse)

            total_stars = len(filtered)
            offset = max(0, (page - 1) * limit)
            page_items = filtered[offset:offset + limit]

            return {
                'total_stars': total_stars,
                'page': page,
                'limit': limit,
                'total_pages': max(1, math.ceil(total_stars / limit)),
                'stars': page_items
            }
        finally:
            cx.close()

    def get_star_details(self, tic: str) -> Optional[Dict[str, Any]]:
        """Get complete star research details including physical parameters, events, and reviews."""
        clean_tic = sanitize_tic(tic)
        cx = self.get_connection()
        try:
            # Query candidate rows for this star
            ov_cand_join = "LEFT JOIN overlay.candidate_review_overlay o_rev ON c.id = o_rev.candidate_id" if self.review_db_path else ""
            ov_star_join = "LEFT JOIN overlay.star_reviews sr_ov ON c.tic = sr_ov.tic" if self.review_db_path else ""

            verif_col = "c.scientific_classification" if self.schema_info['has_scientific_classification'] else "'not_run'"

            if self.schema_info['has_candidate_systems']:
                sys_cols = """
                    s.period_days, s.is_periodic, s.classification as system_classification,
                    s.epoch_btjd, s.odd_even_sigma, s.period_uncertainty_days, s.epoch_uncertainty_days,
                    s.alternate_period_days
                """
                sys_joins = "LEFT JOIN candidate_systems s ON c.system_id = s.id"
            else:
                sys_cols = """
                    NULL as period_days, 0 as is_periodic, '' as system_classification,
                    NULL as epoch_btjd, 0.0 as odd_even_sigma, NULL as period_uncertainty_days,
                    NULL as epoch_uncertainty_days, NULL as alternate_period_days
                """
                sys_joins = ""

            cand_cols = self.schema_info.get('candidates_cols', set())
            extra_cols = [
                "c.scientific_description" if 'scientific_description' in cand_cols else "'' as scientific_description",
                "c.true_depth_min_pct" if 'true_depth_min_pct' in cand_cols else "NULL as true_depth_min_pct",
                "c.has_shoulder_flare" if 'has_shoulder_flare' in cand_cols else "NULL as has_shoulder_flare",
                "c.shoulder_elevation_pct" if 'shoulder_elevation_pct' in cand_cols else "0.0 as shoulder_elevation_pct",
                "c.consecutive_depressed" if 'consecutive_depressed' in cand_cols else "0 as consecutive_depressed",
                "c.coincidence_risk_level" if 'coincidence_risk_level' in cand_cols else ("c.coincidence_star_count" if 'coincidence_star_count' in cand_cols else "'unknown' as coincidence_risk_level")
            ]
            extra_cols_sql = ", ".join(extra_cols)

            sql = f'''
                SELECT c.id, c.run_id, c.tic, c.sector, c.source_file, c.center_btjd, c.score, c.width_hr, c.depth_ppm,
                       c.n_support_windows,
                       COALESCE({'o_rev.status, ' if self.review_db_path else ''}c.review_status, 'unreviewed') as eff_event_status,
                       COALESCE({'o_rev.note, ' if self.review_db_path else ''}c.note, '') as eff_event_note,
                       {verif_col} as verif_class,
                       COALESCE({'sr_ov.status, ' if self.review_db_path else ''}sr.status, 'unreviewed') as eff_star_status,
                       COALESCE({'sr_ov.star_class, ' if self.review_db_path else ''}sr.star_class, '') as eff_star_class,
                       COALESCE({'sr_ov.note, ' if self.review_db_path else ''}sr.note, '') as eff_star_note,
                       {sys_cols},
                       {extra_cols_sql}
                FROM candidates c
                LEFT JOIN star_reviews sr ON c.tic = sr.tic
                {sys_joins}
                {ov_cand_join}
                {ov_star_join}
                WHERE c.tic IN (?,?,?)
                ORDER BY c.score DESC
            '''

            rows = cx.execute(sql, (tic, clean_tic, f"TIC:{clean_tic}")).fetchall()
            if not rows:
                sr_row = None
                try:
                    if self.review_db_path:
                        sr_row = cx.execute(
                            "SELECT status, star_class, note FROM overlay.star_reviews WHERE tic=? OR tic=?",
                            (clean_tic, f"TIC:{clean_tic}")
                        ).fetchone()
                    if not sr_row:
                        sr_row = cx.execute(
                            "SELECT status, star_class, note FROM star_reviews WHERE tic=? OR tic=?",
                            (clean_tic, f"TIC:{clean_tic}")
                        ).fetchone()
                except Exception:
                    pass

                obs_meta = {}
                if self.schema_info.get('has_observations'):
                    try:
                        obs_row = cx.execute('''
                            SELECT camera, ccd, tmag, teff, obs_id, product_filename, data_uri, sector
                            FROM observations WHERE tic=? OR tic=? LIMIT 1
                        ''', (clean_tic, f"TIC:{clean_tic}")).fetchone()
                        if obs_row:
                            obs_meta = {
                                'camera': obs_row[0], 'ccd': obs_row[1], 'tmag': obs_row[2], 'teff': obs_row[3],
                                'obs_id': obs_row[4], 'product_filename': obs_row[5], 'data_uri': obs_row[6],
                                'sector': obs_row[7]
                            }
                    except Exception:
                        pass

                valid_targets = getattr(self, "valid_targets", set())
                is_valid = (tic in valid_targets) or (clean_tic in valid_targets) or (f"TIC:{clean_tic}" in valid_targets)
                if sr_row or obs_meta or is_valid:
                    st_status = sr_row[0] if sr_row else 'unreviewed'
                    st_class = sr_row[1] if sr_row else ''
                    st_note = sr_row[2] if sr_row else ''
                    return {
                        'tic': f"TIC:{clean_tic}",
                        'clean_tic': clean_tic,
                        'star_status': st_status,
                        'star_class': st_class,
                        'note': st_note,
                        'star_note': st_note,
                        'total_events': 0,
                        'reviewed_events': 0,
                        'unreviewed_events': 0,
                        'progress': 'unreviewed',
                        'has_localized_dimmings': False,
                        'system': None,
                        'events': [],
                        'observation': obs_meta,
                        'history': []
                    }
                return None

            first = rows[0]
            star_status = first[13]
            star_class = first[14]
            star_note = first[15]

            # Periodic system info
            system_info = None
            if first[16] is not None:
                system_info = {
                    'period_days': first[16],
                    'is_periodic': bool(first[17]),
                    'classification': first[18] or '',
                    'epoch_btjd': first[19],
                    'odd_even_sigma': first[20] or 0.0,
                    'period_uncertainty_days': first[21],
                    'epoch_uncertainty_days': first[22],
                    'alternate_period_days': first[23]
                }

            # Observation physical parameters if available
            obs_meta = {}
            if self.schema_info['has_observations']:
                obs_row = cx.execute('''
                    SELECT camera, ccd, tmag, teff, obs_id, product_filename, data_uri, sector
                    FROM observations WHERE tic=? OR tic=? LIMIT 1
                ''', (tic, clean_tic)).fetchone()
                if obs_row:
                    obs_meta = {
                        'camera': obs_row[0], 'ccd': obs_row[1], 'tmag': obs_row[2], 'teff': obs_row[3],
                        'obs_id': obs_row[4], 'product_filename': obs_row[5], 'data_uri': obs_row[6],
                        'sector': obs_row[7]
                    }

            events = []
            for r in rows:
                events.append({
                    'id': r[0],
                    'run_id': r[1],
                    'tic': r[2],
                    'sector': r[3],
                    'source_file': r[4],
                    'center_btjd': round(float(r[5]), 4),
                    'score': round(float(r[6]), 2),
                    'width_hr': round(float(r[7]), 2),
                    'depth_ppm': round(float(r[8]), 1),
                    'n_support_windows': r[9],
                    'review_status': r[10],
                    'note': r[11],
                    'scientific_classification': r[12],
                    'scientific_description': r[24] if len(r) > 24 and r[24] is not None else '',
                    'true_depth_min_pct': round(float(r[25]), 3) if len(r) > 25 and r[25] is not None else None,
                    'has_shoulder_flare': bool(r[26]) if len(r) > 26 and r[26] is not None else None,
                    'shoulder_elevation_pct': round(float(r[27]), 2) if len(r) > 27 and r[27] is not None else 0.0,
                    'consecutive_depressed': r[28] if len(r) > 28 and r[28] is not None else 0,
                    'coincidence_risk_level': str(r[29]) if len(r) > 29 and r[29] is not None else 'unknown'
                })

            if self.schema_info['has_observations'] and 'observation_id' in self.schema_info['candidates_cols']:
                for event in events:
                    cursor = cx.execute("SELECT o.* FROM observations o JOIN candidates c ON c.observation_id=o.id AND c.run_id=o.run_id WHERE c.id=?", (event['id'],))
                    observation = cursor.fetchone()
                    event['observation'] = dict(zip([col[0] for col in cursor.description], observation)) if observation else None
                obs_meta = events[0].get('observation') or {}

            # Star history
            hist_table = "overlay.star_review_history" if self.review_db_path else "star_review_history"
            hist_cur = cx.execute(f'''
                SELECT previous_status, new_status, note, reviewer, created_utc
                FROM {hist_table} WHERE tic=? OR tic=? ORDER BY id DESC LIMIT 20
            ''', (clean_tic, f"TIC:{clean_tic}"))
            history = [{'previous_status': h[0], 'new_status': h[1], 'note': h[2],
                        'reviewer': h[3], 'created_utc': h[4]} for h in hist_cur.fetchall()]

            for row in cx.execute("SELECT h.previous_status,h.new_status,h.note,h.reviewer,h.created_utc,h.candidate_id FROM overlay.review_history h JOIN candidates c ON c.id=h.candidate_id WHERE c.tic=? OR c.tic=? ORDER BY h.id DESC LIMIT 20", (clean_tic, f"TIC:{clean_tic}")):
                history.append(dict(previous_status=row[0], new_status=row[1], note=row[2], reviewer=row[3], created_utc=row[4], candidate_id=row[5]))
            history.sort(key=lambda h: h['created_utc'], reverse=True)

            total_events = len(events)
            reviewed_events = sum(1 for e in events if e['review_status'] != 'unreviewed')

            return {
                'tic': tic,
                'clean_tic': clean_tic,
                'data_kind': self.data_kind,
                'execution': experiments.read_source_run(self.db_path, self.run_id),
                'star_status': star_status,
                'star_class': star_class,
                'note': star_note,
                'system': system_info,
                'observation': obs_meta,
                'total_events': total_events,
                'reviewed_events': reviewed_events,
                'unreviewed_events': total_events - reviewed_events,
                'progress': 'reviewed' if reviewed_events == total_events else ('partially_reviewed' if reviewed_events > 0 else 'unreviewed'),
                'events': events,
                'history': history
            }
        finally:
            cx.close()

    def record_star_review(self, tic: str, status: str, star_class: str = "", note: str = "", reviewer: str = "local_user"):
        if status not in VALID_STATUSES:
            raise ValueError(f"Invalid status: {status}")
        now_utc = datetime.now(timezone.utc).isoformat()
        cx = self.get_connection()
        try:
            target_cx = cx
            table_prefix = "overlay." if self.review_db_path else ""

            row = target_cx.execute(f"SELECT status FROM {table_prefix}star_reviews WHERE tic=?", (tic,)).fetchone()
            prev_status = row[0] if row else "unreviewed"

            with target_cx:
                target_cx.execute(f'''
                    INSERT INTO {table_prefix}star_reviews(tic, status, star_class, note, reviewer, updated_utc)
                    VALUES (?,?,?,?,?,?)
                    ON CONFLICT(tic) DO UPDATE SET status=excluded.status, star_class=excluded.star_class,
                    note=excluded.note, reviewer=excluded.reviewer, updated_utc=excluded.updated_utc
                ''', (tic, status, star_class, note, reviewer, now_utc))

                target_cx.execute(f'''
                    INSERT INTO {table_prefix}star_review_history(tic, previous_status, new_status, note, reviewer, created_utc)
                    VALUES (?,?,?,?,?,?)
                ''', (tic, prev_status, status, note, reviewer, now_utc))
        finally:
            cx.close()

    def record_event_review(self, candidate_id: int, status: str, note: str = "", reviewer: str = "local_user"):
        if status not in VALID_STATUSES:
            raise ValueError(f"Invalid status: {status}")
        now_utc = datetime.now(timezone.utc).isoformat()
        cx = self.get_connection()
        try:
            with cx:
                if not cx.execute("SELECT 1 FROM candidates WHERE id=?", (candidate_id,)).fetchone():
                    raise ValueError("Candidate is not part of the selected run.")
                if self.review_db_path:
                    # Write to overlay tables
                    row = cx.execute("SELECT status FROM overlay.candidate_review_overlay WHERE candidate_id=?", (candidate_id,)).fetchone()
                    if not row:
                        orig = cx.execute("SELECT review_status FROM candidates WHERE id=?", (candidate_id,)).fetchone()
                        prev_status = orig[0] if orig else "unreviewed"
                    else:
                        prev_status = row[0]

                    cx.execute('''
                        INSERT INTO overlay.candidate_review_overlay(candidate_id, status, note, reviewer, updated_utc)
                        VALUES (?,?,?,?,?)
                        ON CONFLICT(candidate_id) DO UPDATE SET status=excluded.status, note=excluded.note,
                        reviewer=excluded.reviewer, updated_utc=excluded.updated_utc
                    ''', (candidate_id, status, note, reviewer, now_utc))

                    cx.execute('''
                        INSERT INTO overlay.review_history(candidate_id, previous_status, new_status, note, reviewer, created_utc)
                        VALUES (?,?,?,?,?,?)
                    ''', (candidate_id, prev_status, status, note, reviewer, now_utc))
                else:
                    db.record_review(cx, candidate_id, status, note=note, reviewer=reviewer)
        finally:
            cx.close()

    def undo_star_review(self, tic: str) -> Dict[str, Any]:
        cx = self.get_connection()
        try:
            table_prefix = "overlay." if self.review_db_path else ""
            with cx:
                latest = cx.execute(f'''
                    SELECT id, previous_status FROM {table_prefix}star_review_history
                    WHERE tic=? ORDER BY id DESC LIMIT 1
                ''', (tic,)).fetchone()
                if not latest:
                    raise ValueError(f"No star review history to undo for {tic}")
                hist_id, prev_status = latest
                cx.execute(f"UPDATE {table_prefix}star_reviews SET status=? WHERE tic=?", (prev_status, tic))
                cx.execute(f"DELETE FROM {table_prefix}star_review_history WHERE id=?", (hist_id,))
                return {'tic': tic, 'reverted_to': prev_status}
        finally:
            cx.close()

    def undo_event_review(self, candidate_id: int) -> Dict[str, Any]:
        cx = self.get_connection()
        try:
            with cx:
                if self.review_db_path:
                    latest = cx.execute('''
                        SELECT id, previous_status FROM overlay.review_history
                        WHERE candidate_id=? ORDER BY id DESC LIMIT 1
                    ''', (candidate_id,)).fetchone()
                    if not latest:
                        raise ValueError(f"No review history to undo for candidate {candidate_id}")
                    hist_id, prev_status = latest
                    cx.execute("UPDATE overlay.candidate_review_overlay SET status=? WHERE candidate_id=?", (prev_status, candidate_id))
                    cx.execute("DELETE FROM overlay.review_history WHERE id=?", (hist_id,))
                    return {'candidate_id': candidate_id, 'reverted_to': prev_status}
                else:
                    return db.undo_review(cx, candidate_id)
        finally:
            cx.close()

    def bulk_review_star_events(self, tic: str, status: str, note: str = "", reviewer: str = "local_user", only_unreviewed: bool = True) -> int:
        if status not in VALID_STATUSES:
            raise ValueError(f"Invalid status: {status}")
        cx = self.get_connection()
        try:
            clean_tic = sanitize_tic(tic)
            cands = cx.execute("SELECT id, review_status FROM candidates WHERE tic=? OR tic=?", (clean_tic, f"TIC:{clean_tic}")).fetchall()
            updated = 0
            for cid, orig_st in cands:
                # Check current effective status
                if self.review_db_path:
                    ov_st = cx.execute("SELECT status FROM overlay.candidate_review_overlay WHERE candidate_id=?", (cid,)).fetchone()
                    curr_st = ov_st[0] if ov_st else orig_st
                else:
                    curr_st = orig_st

                if only_unreviewed and curr_st != 'unreviewed':
                    continue

                self.record_event_review(cid, status, note=note, reviewer=reviewer)
                updated += 1

            self.record_star_review(tic, status, note=note, reviewer=reviewer)
            return updated
        finally:
            cx.close()


# ---------------------------------------------------------------------------
# Static Plot Rendering (Matplotlib Fallback)
# ---------------------------------------------------------------------------

def render_plot_png(source_file, center_btjd, width_hr, score, depth_ppm, tic,
                    root_path=None, cache_dir=None, mode="zoom", hours=12.0, other_centers=None, sector=None, product_filename=None, data_uri=None, synthetic=False, expected_sha256=None):
    """Generate in-memory static PNG plot with on-demand fallback support."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 4.2))
    if synthetic:
        fig.suptitle("SYNTHETIC DEMONSTRATION — not astronomical detections", color='darkorange')

    try:
        t, f_ppm, good, cadence, src_label = resolve_light_curve_series(
            tic=tic, sector=sector, product_filename=product_filename, data_uri=data_uri, expected_sha256=expected_sha256, source_file=source_file, root_path=root_path, cache_dir=cache_dir, allow_download=True
        )
    except Exception as e:
        ax.text(0.5, 0.5, f"Light Curve Data Unavailable:\n{str(e)}\n\n(Pass --root <path> or check internet access for MAST streaming)",
                ha='center', va='center', fontsize=11, color='darkorange')
        ax.set_axis_off()
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=120, bbox_inches='tight')
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()

    if mode == "zoom":
        delta_hr = abs(t - center_btjd) * 24.0
        mask = (delta_hr <= hours) & good
        if not np_any(mask):
            mask = (delta_hr <= hours)
        ax.scatter(t[mask], f_ppm[mask], s=12, alpha=0.75, color="#2563eb", label="TESS PDC Flux")
        half_w = (width_hr / 24.0) / 2.0
        ax.axvspan(center_btjd - half_w, center_btjd + half_w, color="#f59e0b", alpha=0.25,
                   label=f"Template width ({width_hr:.1f}h)")
        ax.axvline(center_btjd, color="#dc2626", linestyle="--", linewidth=1.5,
                   label=f"Center (BTJD {center_btjd:.4f})")
        ax.axhline(0, color="gray", linestyle=":", linewidth=0.8)
        ax.set_xlabel("TESS Time (BTJD)", fontsize=10)
        ax.set_ylabel("Relative Flux (ppm)", fontsize=10)
        ax.set_title(f"{tic} | Zoom View (±{hours:.0f}h) | SNR {score:.2f} | Depth {depth_ppm:.0f} ppm [{src_label}]",
                     fontsize=11, fontweight="bold")
        ax.legend(loc="best", fontsize=9, framealpha=0.9)
    else:
        step = max(1, len(t) // 6000)
        idx = slice(None, None, step)
        ax.plot(t[idx], f_ppm[idx], color="#475569", alpha=0.5, linewidth=0.7)
        ax.scatter(t[idx], f_ppm[idx], s=3, color="#0284c7", alpha=0.6, label="PDC Flux")
        ax.axvline(center_btjd, color="#dc2626", linestyle="--", linewidth=1.8,
                   label=f"Selected Candidate ({center_btjd:.4f})")
        if other_centers:
            for oc in other_centers:
                if abs(oc - center_btjd) > 1e-4:
                    ax.axvline(oc, color="#f59e0b", linestyle=":", linewidth=1.2, alpha=0.8)
        ax.axhline(0, color="gray", linestyle=":", linewidth=0.8)
        ax.set_xlabel("TESS Time (BTJD)", fontsize=10)
        ax.set_ylabel("Relative Flux (ppm)", fontsize=10)
        ax.set_title(f"{tic} | Full Observation | Span: {t[-1]-t[0]:.1f} days [{src_label}]",
                     fontsize=11, fontweight="bold")
        ax.legend(loc="upper right", fontsize=9, framealpha=0.9)

    ax.grid(True, linestyle="--", alpha=0.4)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=120)
    plt.close(fig)
    buf.seek(0)
    return buf.getvalue()


def render_phase_plot_png(source_file, center_btjd, width_hr, score, depth_ppm, tic,
                          epoch_btjd, period_days, is_periodic, odd_even_sigma,
                          period_err=None, epoch_err=None, is_alt=False,
                          root_path=None, cache_dir=None, sector=None, product_filename=None, data_uri=None, synthetic=False, expected_sha256=None):
    """Generate in-memory phase-folded static PNG plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 4.2))
    if synthetic:
        fig.suptitle("SYNTHETIC DEMONSTRATION — not astronomical detections", color='darkorange')

    if not is_periodic or period_days <= 0:
        ax.text(0.5, 0.5, "Candidate is not part of a confirmed periodic system.\n(Requires N >= 3 coherent events and acceptable timing residuals)",
                ha='center', va='center', fontsize=11, color='gray')
        ax.set_axis_off()
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=120, bbox_inches='tight')
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()

    try:
        t, f_ppm, good, cadence, src_label = resolve_light_curve_series(
            tic=tic, sector=sector, product_filename=product_filename, data_uri=data_uri, expected_sha256=expected_sha256, source_file=source_file, root_path=root_path, cache_dir=cache_dir, allow_download=True
        )
    except Exception as e:
        ax.text(0.5, 0.5, f"Light Curve Data Unavailable:\n{str(e)}",
                ha='center', va='center', fontsize=11, color='darkorange')
        ax.set_axis_off()
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=120, bbox_inches='tight')
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()

    t_valid = t[good]
    f_valid = f_ppm[good]
    phase = periodicity.compute_phase(t_valid, epoch_btjd, period_days)

    cycle = np.round((t_valid - epoch_btjd) / period_days).astype(int)
    odd_mask = (cycle % 2 == 1)
    even_mask = (cycle % 2 == 0)

    if np_any(even_mask):
        ax.scatter(phase[even_mask], f_valid[even_mask], s=7, color="#2563eb", alpha=0.5,
                   label="Even transits (primary)")
    if np_any(odd_mask):
        ax.scatter(phase[odd_mask], f_valid[odd_mask], s=7, color="#f59e0b", alpha=0.5,
                   label="Odd transits (secondary)")

    # Binned median curve
    sort_idx = np.argsort(phase)
    p_sorted = phase[sort_idx]
    f_sorted = f_valid[sort_idx]
    n_bins = 60
    bin_edges = np.linspace(-0.5, 0.5, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    binned_f, binned_p = [], []
    for i in range(n_bins):
        in_bin = (p_sorted >= bin_edges[i]) & (p_sorted < bin_edges[i+1])
        if np_any(in_bin):
            binned_f.append(np.median(f_sorted[in_bin]))
            binned_p.append(bin_centers[i])
    if binned_p:
        ax.plot(binned_p, binned_f, color="#dc2626", linewidth=2.0, label="Binned median")

    ax.axvline(0, color="gray", linestyle="--", linewidth=1.0)
    ax.axhline(0, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlim(-0.5, 0.5)
    ax.set_xlabel("Phase (orbital cycles)", fontsize=10)
    ax.set_ylabel("Relative Flux (ppm)", fontsize=10)

    p_str = f"P = {period_days:.4f}" + (f" ± {period_err:.4f}" if period_err is not None else "") + " d"
    if is_alt:
        p_str += " (Alternate Solution)"
    t0_str = f"T0 = {epoch_btjd:.4f}" + (f" ± {epoch_err:.4f}" if epoch_err is not None else "") + " BTJD"
    eb_note = f" | Disparity: {odd_even_sigma:.1f}σ" if odd_even_sigma and odd_even_sigma > 0 else ""
    ax.set_title(f"{tic} | Phase-Folded | {p_str} | {t0_str}{eb_note}", fontsize=11, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.grid(True, linestyle="--", alpha=0.4)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=120)
    plt.close(fig)
    buf.seek(0)
    return buf.getvalue()


def np_any(arr):
    try:
        return arr.any()
    except Exception:
        return any(arr)


# ---------------------------------------------------------------------------
# HTTP Handler
# ---------------------------------------------------------------------------

class AnomalyReviewHandler(BaseHTTPRequestHandler):
    def __init__(self, db_path, root_path=None, *args, cache_dir=None, review_db=None, read_only=False, **kwargs):
        self.db_path = db_path
        self.root_path = root_path
        self.cache_dir = cache_dir
        self.review_db = review_db
        self.read_only = read_only
        self._default_adapter = None
        self._default_adapter_args = (db_path, review_db, read_only)
        super().__init__(*args, **kwargs)

    @property
    def adapter(self) -> SurveyDatabaseAdapter:
        if hasattr(self.server, "experiment_store") and self.server.experiment_store is not None:
            exp_id, run_id = self.server.experiment_store.get_active_selection()
            selected = urllib.parse.parse_qs(urllib.parse.urlparse(self.path).query).get('workspace_run_id', [None])[0]
            if selected:
                selected_run = self.server.experiment_store.get_run(selected)
                if not selected_run:
                    raise ValueError("Requested workspace run does not exist")
                run_id, exp_id = selected, selected_run['experiment_id']
            if run_id:
                run = self.server.experiment_store.get_run(run_id)
                exp = self.server.experiment_store.get_experiment(exp_id) if exp_id else None
                if run and run.get("output_db_path"):
                    out_path = str(Path(run["output_db_path"]).resolve())
                    if not os.path.exists(out_path):
                        con = streaming.init_streaming_db(out_path)
                        con.close()
                    rev_path = str(Path(run["review_db_path"]).resolve()) if run.get("review_db_path") else None
                    is_ro = bool(exp and exp.get("is_read_only", 0) == 1)
                    key = (out_path, rev_path, is_ro, run.get("source_run_id"))
                    if not hasattr(self.server, "adapters_cache"):
                        self.server.adapters_cache = {}
                    if key not in self.server.adapters_cache:
                        self.server.adapters_cache[key] = SurveyDatabaseAdapter(
                            out_path, review_db_path=rev_path, read_only=is_ro, run_id=run.get("source_run_id")
                        )
                    ad = self.server.adapters_cache[key]
                    if exp and exp.get("targets"):
                        ad.valid_targets = {sanitize_tic(t) for t in exp["targets"]}
                    return ad
        if self._default_adapter is None and self._default_adapter_args[0]:
            self._default_adapter = SurveyDatabaseAdapter(*self._default_adapter_args)
        if self._default_adapter is not None:
            return self._default_adapter
        raise RuntimeError("No active survey database adapter available.")

    def _is_localhost(self) -> bool:
        client_ip = self.client_address[0]
        return client_ip in ("127.0.0.1", "::1", "localhost")

    def log_message(self, format, *args):
        # Silence console log spam for assets
        pass

    def send_json(self, data, status=HTTPStatus.OK):
        payload = json.dumps(data).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-cache, no-store, must-revalidate")
        self.end_headers()
        self.wfile.write(payload)

    def send_png(self, png_bytes, status=HTTPStatus.OK):
        self.send_response(status)
        self.send_header("Content-Type", "image/png")
        self.send_header("Content-Length", str(len(png_bytes)))
        self.send_header("Cache-Control", "public, max-age=30")
        self.end_headers()
        self.wfile.write(png_bytes)

    def send_html(self, html_str, status=HTTPStatus.OK):
        payload = html_str.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def _request_allowed(self, write=False):
        port = self.server.server_address[1]
        authorities = {f'127.0.0.1:{port}', f'localhost:{port}', f'[::1]:{port}'}
        host = self.headers.get('Host', '').lower()
        origin = self.headers.get('Origin')
        allowed = self._is_localhost() and host in authorities
        if origin is not None:
            allowed = allowed and origin == 'http://' + host
        if self.headers.get('Sec-Fetch-Site') in ('cross-site', 'same-site'):
            allowed = False
        if write:
            allowed = allowed and self.headers.get('Content-Type', '').split(';')[0] == 'application/json'
            # Browser writes require the token embedded in the same-origin page.
            # Non-browser JSON clients without Origin remain supported on loopback.
            if origin is not None or self.headers.get('Sec-Fetch-Site'):
                allowed = allowed and secrets.compare_digest(self.headers.get('X-CSRF-Token', ''), getattr(self.server, 'csrf_token', 'missing'))
        if not allowed:
            self.send_json({'error': 'Forbidden: invalid authority, origin, or CSRF token'}, HTTPStatus.FORBIDDEN)
        return allowed

    def do_GET(self):
        if not self._request_allowed():
            return
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        query = urllib.parse.parse_qs(parsed.query)

        if path == "/":
            if not hasattr(self.server, 'csrf_token'):
                self.server.csrf_token = secrets.token_urlsafe(32)
            self.send_html(INDEX_HTML.replace('__CSRF_TOKEN__', self.server.csrf_token))
            return

        if path == "/api/experiments":
            if hasattr(self.server, "experiment_store") and self.server.experiment_store is not None:
                exps = self.server.experiment_store.list_experiments()
                act_exp, act_run = self.server.experiment_store.get_active_selection()
                self.send_json({
                    "experiments": exps,
                    "active_experiment": self.server.experiment_store.get_experiment(act_exp) if act_exp else None,
                    "active_run": self.server.experiment_store.get_run(act_run) if act_run else None,
                    "presets": experiments.DETECTOR_PRESETS
                })
            else:
                self.send_json({
                    "experiments": [],
                    "active_experiment": None,
                    "active_run": None,
                    "presets": experiments.DETECTOR_PRESETS
                })
            return

        if path == "/api/experiments/active":
            if hasattr(self.server, "experiment_store") and self.server.experiment_store is not None:
                act_exp, act_run = self.server.experiment_store.get_active_selection()
                exp = self.server.experiment_store.get_experiment(act_exp) if act_exp else None
                run = self.server.experiment_store.get_run(act_run) if act_run else None
                is_act = self.server.execution_manager.is_run_active(act_run) if (hasattr(self.server, "execution_manager") and act_run) else False
                self.send_json({
                    "experiment": exp,
                    "run": run,
                    "is_active": is_act
                })
            else:
                self.send_json({"experiment": None, "run": None, "is_active": False})
            return

        if path.startswith("/api/runs/") and path.endswith("/progress"):
            parts = path.split("/")
            run_id = parts[3]
            if hasattr(self.server, "execution_manager") and hasattr(self.server, "experiment_store"):
                telem = self.server.execution_manager.get_live_telemetry(run_id)
                logs = self.server.experiment_store.get_run_logs(run_id, 30)
                is_act = self.server.execution_manager.is_run_active(run_id)
                self.send_json({
                    "telemetry": telem,
                    "run": self.server.experiment_store.get_run(run_id),
                    "logs": logs,
                    "is_active": is_act
                })
            else:
                self.send_json({"error": "Execution manager unavailable"}, HTTPStatus.NOT_FOUND)
            return

        if path.startswith("/api/runs/") and path.endswith("/export"):
            parts = path.split("/")
            run_id = parts[3]
            fmt = query.get("format", ["csv"])[0].lower()
            if hasattr(self.server, "experiment_store") and self.server.experiment_store is not None:
                run = self.server.experiment_store.get_run(run_id)
                if not run:
                    self.send_json({"error": "Run not found"}, HTTPStatus.NOT_FOUND)
                    return
                try:
                    content, mime = self.server.experiment_store.export_run_catalog(run["experiment_id"], run_id, fmt)
                    ext = "csv" if fmt == "csv" else "json"
                    fname = f"tess_catalog_{run['experiment_id']}_{run_id}.{ext}"
                    payload = content.encode("utf-8")
                    self.send_response(HTTPStatus.OK)
                    self.send_header("Content-Type", f"{mime}; charset=utf-8")
                    self.send_header("Content-Disposition", f'attachment; filename="{fname}"')
                    self.send_header("Content-Length", str(len(payload)))
                    self.end_headers()
                    self.wfile.write(payload)
                    return
                except Exception as e:
                    self.send_json({"error": str(e)}, HTTPStatus.INTERNAL_SERVER_ERROR)
                    return
            else:
                self.send_json({"error": "Experiment store unavailable"}, HTTPStatus.NOT_FOUND)
            return

        if path == "/api/stats":
            stats = self.adapter.get_survey_stats()
            stats['root_path'] = str(self.root_path) if self.root_path else None
            self.send_json(stats)
            return

        # Star-centered explore catalog endpoint
        if path == "/api/stars":
            page = int(query.get('page', ['1'])[0])
            limit = int(query.get('limit', ['25'])[0])
            status_f = query.get('status', ['all'])[0]
            verif_f = query.get('verification', ['all'])[0]
            sys_f = query.get('system_type', ['all'])[0]
            coinc_f = query.get('coincidence', ['all'])[0]
            search_q = query.get('search', [''])[0]
            sort_by = query.get('sort_by', ['score'])[0]
            sort_order = query.get('sort_order', ['desc'])[0]

            res = self.adapter.query_stars(
                page=page, limit=limit, status_filter=status_f,
                verification_filter=verif_f, system_filter=sys_f,
                coincidence_filter=coinc_f, search_query=search_q,
                sort_by=sort_by, sort_order=sort_order
            )
            self.send_json(res)
            return

        # Star investigation detail endpoint
        if path.startswith("/api/star/") and not path.endswith("/lightcurve") and not path.endswith("/review"):
            tic = urllib.parse.unquote(path[len("/api/star/"):])
            detail = self.adapter.get_star_details(tic)
            if not detail:
                self.send_json({'error': 'Star not found'}, HTTPStatus.NOT_FOUND)
                return
            self.send_json(detail)
            return

        # Interactive light curve data for Plotly.js
        if path.startswith("/api/star/") and path.endswith("/lightcurve"):
            parts = path.split("/")
            tic = urllib.parse.unquote(parts[3])
            detail = self.adapter.get_star_details(tic)
            if not detail:
                self.send_json({'error': 'Star not found'}, HTTPStatus.NOT_FOUND)
                return

            selected = detail['events']
            event_id = query.get('candidate_id', [None])[0]
            sector_query = query.get('sector', [None])[0]
            if event_id:
                selected = [e for e in selected if str(e['id']) == event_id]
            elif sector_query:
                selected = [e for e in selected if str(e['sector']) == sector_query]
            event = selected[0] if selected else None
            if (event_id or sector_query) and not event:
                self.send_json({'error': 'Candidate/product not in selected run'}, HTTPStatus.NOT_FOUND)
                return
            observation = (event.get('observation') if event else detail.get('observation')) or {}
            sector = event['sector'] if event else observation.get('sector')
            source_file = event['source_file'] if event else None
            obs_prod = observation.get('product_filename')
            detail['events'] = [e for e in detail['events'] if (e.get('observation') or {}).get('id') == observation.get('id')] if observation else selected[:1]

            try:
                t, f_ppm, good, cadence, src_label = resolve_light_curve_series(
                    tic=tic, sector=sector, source_file=source_file, product_filename=obs_prod, data_uri=observation.get("data_uri"), expected_sha256=observation.get("file_sha256"),
                    root_path=self.root_path, cache_dir=self.cache_dir, allow_download=True
                )

                if detail.get('data_kind') == 'synthetic demonstration':
                    src_label = 'SYNTHETIC DEMONSTRATION: ' + src_label
                # Downsample full observation for rapid browser rendering if point count is large
                # Retains full density around detected event centers!
                event_centers = [e['center_btjd'] for e in detail['events']]
                downsample = len(t) > 6000
                if downsample:
                    step = max(1, len(t) // 4000)
                    near_event = np.zeros(len(t), dtype=bool)
                    for ec in event_centers:
                        near_event |= (np.abs(t - ec) * 24.0 <= 18.0) # Full cadence within +/-18h of any event
                    keep_mask = good & (near_event | (np.arange(len(t)) % step == 0))
                else:
                    keep_mask = good

                t_out = t[keep_mask].tolist()
                f_out = f_ppm[keep_mask].tolist()

                self.send_json({
                    'tic': tic,
                    'sector': sector,
                    'cadence_days': cadence,
                    'source': src_label,
                    'time': t_out,
                    'flux_ppm': f_out,
                    'events': detail['events'],
                    'system': detail['system']
                })
            except Exception as e:
                self.send_json({'error': f"Failed to retrieve light curve: {str(e)}"}, HTTPStatus.INTERNAL_SERVER_ERROR)
            return

        # Backwards-compatible candidates endpoint
        if path == "/api/candidates":
            status_filter = query.get('status', ['all'])[0]
            search_query = query.get('search', [''])[0].strip()
            tic_filter = query.get('tic', [''])[0].strip()
            max_per_star = int(query.get('max_per_star', ['0'])[0])
            limit = int(query.get('limit', ['100'])[0])
            offset = int(query.get('offset', ['0'])[0])

            star_res = self.adapter.query_stars(
                page=1, limit=10000, status_filter='all', search_query=search_query
            )

            cands = []
            for s in star_res['stars']:
                if tic_filter and s['tic'] != tic_filter and s['clean_tic'] != sanitize_tic(tic_filter):
                    continue
                for e in s['events']:
                    if status_filter != 'all' and e['review_status'] != status_filter:
                        continue
                    cands.append(e)

            page_cands = cands[offset:offset + limit] if limit > 0 else cands[offset:]
            self.send_json({
                'total': len(cands),
                'candidates': page_cands,
                'offset': offset,
                'limit': limit
            })
            return

        # Candidate detail endpoint
        if path.startswith("/api/candidate/"):
            try:
                cid = int(path.split("/")[-1])
            except ValueError:
                self.send_json({'error': 'Invalid candidate ID'}, HTTPStatus.BAD_REQUEST)
                return

            cx = self.adapter.get_connection()
            try:
                row = cx.execute("SELECT tic FROM candidates WHERE id=?", (cid,)).fetchone()
                if not row:
                    self.send_json({'error': 'Candidate not found'}, HTTPStatus.NOT_FOUND)
                    return
                tic = row[0]
            finally:
                cx.close()

            detail = self.adapter.get_star_details(tic)
            target_cand = next((e for e in detail['events'] if e['id'] == cid), None)

            # Get candidate detections provenance if available
            detections = []
            cx = self.adapter.get_connection()
            try:
                if self.adapter.schema_info['has_detections']:
                    d_cur = cx.execute('''SELECT id, center_btjd, score, width_hr, depth_ppm, created_utc
                        FROM candidate_detections WHERE candidate_id=? ORDER BY score DESC''', (cid,))
                    detections = [{'id': d[0], 'center_btjd': d[1], 'score': d[2],
                                   'width_hr': d[3], 'depth_ppm': d[4], 'created_utc': d[5]} for d in d_cur.fetchall()]
            finally:
                cx.close()

            self.send_json({
                'candidate': target_cand,
                'history': detail['history'],
                'detections': detections,
                'star': detail
            })
            return

        # Static Matplotlib plot endpoint (fallback)
        if path.startswith("/api/plot/"):
            parts = path.split("/")
            if len(parts) >= 5:
                mode = parts[3]
                try:
                    cid = int(parts[4])
                except ValueError:
                    self.send_json({'error': 'Invalid ID'}, HTTPStatus.BAD_REQUEST)
                    return

                hours = float(query.get('hours', ['12'])[0])
                alt_param = query.get('alt', ['0'])[0].lower() in ('1', 'true', 'yes')

                cx = self.adapter.get_connection()
                try:
                    row = cx.execute('''SELECT c.source_file, c.center_btjd, c.width_hr, c.score, c.depth_ppm, c.tic,
                                               s.epoch_btjd, s.period_days, s.is_periodic, s.odd_even_sigma,
                                               s.period_uncertainty_days, s.epoch_uncertainty_days, s.alternate_period_days
                                        FROM candidates c
                                        LEFT JOIN candidate_systems s ON c.system_id = s.id
                                        WHERE c.id=?''', (cid,)).fetchone()
                    if not row:
                        self.send_json({'error': 'Candidate not found'}, HTTPStatus.NOT_FOUND)
                        return
                    (source_file, center_btjd, width_hr, score, depth_ppm, tic,
                     epoch_btjd, period_days, is_periodic, odd_even_sigma,
                     p_err, t0_err, alt_period) = row

                    other_centers = [r[0] for r in cx.execute('SELECT center_btjd FROM candidates WHERE tic=?', (tic,)).fetchall()]
                finally:
                    cx.close()

                event = next(e for e in self.adapter.get_star_details(tic)['events'] if e['id'] == cid)
                observation = event.get('observation') or {}
                identity = dict(sector=event['sector'], product_filename=observation.get('product_filename'), data_uri=observation.get('data_uri'), expected_sha256=observation.get('file_sha256'), synthetic=self.adapter.data_kind == 'synthetic demonstration')
                if mode == "phase":
                    fold_period = period_days if period_days is not None else 0.0
                    is_alt = False
                    if alt_param and alt_period:
                        fold_period = alt_period
                        is_alt = True
                    png = render_phase_plot_png(source_file, center_btjd, width_hr, score, depth_ppm, tic,
                                                epoch_btjd=epoch_btjd if epoch_btjd is not None else center_btjd,
                                                period_days=fold_period,
                                                is_periodic=bool(is_periodic),
                                                odd_even_sigma=odd_even_sigma if odd_even_sigma is not None else 0.0,
                                                period_err=p_err if not is_alt else None,
                                                epoch_err=t0_err,
                                                is_alt=is_alt,
                                                root_path=self.root_path,
                                                cache_dir=self.cache_dir, **identity)
                else:
                    png = render_plot_png(source_file, center_btjd, width_hr, score, depth_ppm, tic,
                                          root_path=self.root_path, cache_dir=self.cache_dir,
                                          mode=mode, hours=hours, other_centers=other_centers, **identity)
                self.send_png(png)
                return

        self.send_json({'error': 'Not Found'}, HTTPStatus.NOT_FOUND)

    def do_POST(self):
        parsed = urllib.parse.urlparse(self.path)

        if not self._request_allowed(write=True):
            return

        length = int(self.headers.get("Content-Length", 0))
        if length <= 0 or length > 524288:
            self.send_json({'error': 'Payload size invalid'}, HTTPStatus.BAD_REQUEST)
            return
        body_bytes = self.rfile.read(length)
        try:
            data = json.loads(body_bytes.decode('utf-8'))
        except json.JSONDecodeError:
            self.send_json({'error': 'Malformed JSON'}, HTTPStatus.BAD_REQUEST)
            return

        # 0a. Switch Active Experiment / Run
        if parsed.path == "/api/experiments/active":
            exp_id = data.get("experiment_id")
            run_id = data.get("run_id")
            if not exp_id:
                self.send_json({"error": "Missing experiment_id"}, HTTPStatus.BAD_REQUEST)
                return
            if hasattr(self.server, "experiment_store") and self.server.experiment_store is not None:
                if not run_id:
                    runs = self.server.experiment_store.list_runs(exp_id)
                    if runs:
                        run_id = runs[0]["run_id"]
                if not run_id:
                    self.send_json({"error": f"No runs found for experiment {exp_id}"}, HTTPStatus.BAD_REQUEST)
                    return
                self.server.experiment_store.set_active_selection(exp_id, run_id)
                exp = self.server.experiment_store.get_experiment(exp_id)
                run = self.server.experiment_store.get_run(run_id)
                self.send_json({"success": True, "experiment": exp, "active_run": run, "run": run})
            else:
                self.send_json({"error": "Experiment store unavailable"}, HTTPStatus.BAD_REQUEST)
            return

        # 0b. Pre-flight Validate Targets
        if parsed.path == "/api/experiments/validate-targets":
            raw_targets = data.get("targets_raw") or data.get("targets_text", "")
            sectors = data.get("sectors", [])
            val_res = experiments.validate_and_estimate_targets(raw_targets, sectors=sectors)
            val_res["success"] = val_res.get("is_safe", False)
            val_res["total_unique_targets"] = val_res.get("total_unique_stars", 0)
            val_res["invalid_count"] = len(val_res.get("invalid_tokens", []))
            self.send_json(val_res)
            return

        # 0c. Create New Experiment
        if parsed.path == "/api/experiments/create":
            name = str(data.get("name", "")).strip()
            if not name:
                self.send_json({"error": "Experiment name is required"}, HTTPStatus.BAD_REQUEST)
                return
            desc = str(data.get("description", "")).strip()
            obj = str(data.get("scientific_objective") or data.get("objective", "")).strip()
            targets = data.get("targets", [])
            if not targets and data.get("targets_text"):
                val_t = experiments.validate_and_estimate_targets(data["targets_text"], sectors=data.get("sectors", []))
                targets = val_t.get("valid_tics", [])
            sectors = data.get("sectors", [])
            preset = data.get("config_preset", "frozen_fixed_box_2h")
            threshold_t = data.get("threshold_t")
            duration_hr = data.get("duration_hr")
            is_ro = bool(data.get("is_read_only", False))
            launch_immediately = bool(data.get("launch_immediately", False))

            if not targets:
                self.send_json({"error": "At least one target TIC is required"}, HTTPStatus.BAD_REQUEST)
                return

            if hasattr(self.server, "experiment_store") and self.server.experiment_store is not None:
                try:
                    exp = self.server.experiment_store.create_experiment(
                        name=name,
                        description=desc,
                        scientific_objective=obj,
                        targets=targets,
                        sectors=sectors,
                        config_preset=preset,
                        threshold_t=threshold_t,
                        duration_hr=duration_hr,
                        is_read_only=is_ro
                    )
                    run = exp["runs"][0] if exp.get("runs") else None
                    if launch_immediately and run and hasattr(self.server, "execution_manager"):
                        self.server.execution_manager.launch_survey(
                            experiment_id=exp["id"],
                            run_id=run["run_id"],
                            max_workers=int(data.get("workers", 2)),
                            resume=True
                        )
                    self.send_json({"success": True, "experiment": exp, "run": run})
                except Exception as e:
                    self.send_json({"error": str(e)}, HTTPStatus.INTERNAL_SERVER_ERROR)
            else:
                self.send_json({"error": "Experiment store unavailable"}, HTTPStatus.BAD_REQUEST)
            return

        # 0d. Launch Survey Execution
        if parsed.path == "/api/experiments/launch":
            exp_id = data.get("experiment_id")
            run_id = data.get("run_id")
            workers = int(data.get("max_workers", 2))
            resume = bool(data.get("resume", True))

            if not exp_id or not run_id:
                self.send_json({"error": "Missing experiment_id or run_id"}, HTTPStatus.BAD_REQUEST)
                return

            if hasattr(self.server, "execution_manager") and self.server.execution_manager is not None:
                try:
                    telem = self.server.execution_manager.launch_survey(
                        experiment_id=exp_id,
                        run_id=run_id,
                        max_workers=workers,
                        resume=resume
                    )
                    self.send_json({"success": True, "telemetry": telem, "run": self.server.experiment_store.get_run(run_id)})
                except RuntimeError as e:
                    self.send_json({"error": str(e)}, HTTPStatus.CONFLICT)
                except Exception as e:
                    self.send_json({"error": str(e)}, HTTPStatus.INTERNAL_SERVER_ERROR)
            else:
                self.send_json({"error": "Execution manager unavailable"}, HTTPStatus.BAD_REQUEST)
            return

        # 0e. Abort Survey Execution
        if parsed.path == "/api/runs/abort":
            run_id = data.get("run_id")
            if not run_id:
                self.send_json({"error": "Missing run_id"}, HTTPStatus.BAD_REQUEST)
                return
            if hasattr(self.server, "execution_manager") and self.server.execution_manager is not None:
                aborted = self.server.execution_manager.abort_survey(run_id)
                self.send_json({"success": True, "aborted": aborted})
            else:
                self.send_json({"error": "Execution manager unavailable"}, HTTPStatus.BAD_REQUEST)
            return

        # 1. Event Review
        if parsed.path == "/api/review":
            candidate_id = data.get('candidate_id')
            new_status = data.get('status')
            note = str(data.get('note', '')).strip()[:2000]
            reviewer = str(data.get('reviewer', 'web_user')).strip()[:50]

            if not isinstance(candidate_id, int) or new_status not in VALID_STATUSES:
                self.send_json({'error': 'Invalid review payload parameters'}, HTTPStatus.BAD_REQUEST)
                return

            try:
                self.adapter.record_event_review(candidate_id, new_status, note=note, reviewer=reviewer)
                cx = self.adapter.get_connection()
                try:
                    tic_row = cx.execute("SELECT tic FROM candidates WHERE id=?", (candidate_id,)).fetchone()
                    star_summary = self.adapter.get_star_details(tic_row[0]) if tic_row else None
                finally:
                    cx.close()

                self.send_json({
                    'success': True,
                    'candidate_id': candidate_id,
                    'status': new_status,
                    'note': note,
                    'star': star_summary
                })
            except Exception as e:
                self.send_json({'error': str(e)}, HTTPStatus.INTERNAL_SERVER_ERROR)
            return

        # 2. Event Review Undo
        if parsed.path == "/api/review/undo":
            candidate_id = data.get('candidate_id')
            if not isinstance(candidate_id, int):
                self.send_json({'error': 'Invalid candidate ID'}, HTTPStatus.BAD_REQUEST)
                return
            try:
                res = self.adapter.undo_event_review(candidate_id)
                self.send_json({'success': True, 'undo': res})
            except Exception as e:
                self.send_json({'error': str(e)}, HTTPStatus.BAD_REQUEST)
            return

        # 3. Star-Level Review
        if parsed.path == "/api/star/review":
            tic = data.get('tic')
            status = data.get('status')
            note = str(data.get('note', '')).strip()[:2000]
            reviewer = str(data.get('reviewer', 'web_user')).strip()[:50]
            star_class = str(data.get('star_class', '')).strip()[:50]

            if not tic or status not in VALID_STATUSES:
                self.send_json({'error': 'Invalid star review parameters'}, HTTPStatus.BAD_REQUEST)
                return

            try:
                self.adapter.record_star_review(tic, status, star_class=star_class, note=note, reviewer=reviewer)
                summary = self.adapter.get_star_details(tic)
                self.send_json({'success': True, 'star': summary})
            except Exception as e:
                self.send_json({'error': str(e)}, HTTPStatus.INTERNAL_SERVER_ERROR)
            return

        # 4. Star Review Undo
        if parsed.path == "/api/star/undo":
            tic = data.get('tic')
            if not tic:
                self.send_json({'error': 'Invalid TIC'}, HTTPStatus.BAD_REQUEST)
                return
            try:
                res = self.adapter.undo_star_review(tic)
                summary = self.adapter.get_star_details(tic)
                self.send_json({'success': True, 'undo': res, 'star': summary})
            except Exception as e:
                self.send_json({'error': str(e)}, HTTPStatus.BAD_REQUEST)
            return

        # 5. Star Bulk Event Review
        if parsed.path == "/api/star/bulk-review":
            tic = data.get('tic')
            status = data.get('status')
            note = str(data.get('note', '')).strip()[:2000]
            reviewer = str(data.get('reviewer', 'web_user')).strip()[:50]
            only_unreviewed = bool(data.get('only_unreviewed', True))

            if not tic or status not in VALID_STATUSES:
                self.send_json({'error': 'Invalid bulk review parameters'}, HTTPStatus.BAD_REQUEST)
                return

            try:
                updated_count = self.adapter.bulk_review_star_events(
                    tic, status, note=note, reviewer=reviewer, only_unreviewed=only_unreviewed
                )
                summary = self.adapter.get_star_details(tic)
                self.send_json({'success': True, 'updated_count': updated_count, 'star': summary})
            except Exception as e:
                self.send_json({'error': str(e)}, HTTPStatus.INTERNAL_SERVER_ERROR)
            return

        self.send_json({'error': 'Not Found'}, HTTPStatus.NOT_FOUND)


def serve_review_ui(
    db_path=None,
    root_path=None,
    host="127.0.0.1",
    port=8000,
    cache_dir=None,
    review_db=None,
    read_only=False,
    experiments_dir="experiments"
):
    """Launch the star-centered local research investigation workspace."""
    cache_path = Path(cache_dir) if cache_dir else get_default_cache_dir()
    cache_path.mkdir(parents=True, exist_ok=True)

    store = experiments.ExperimentStore(workspace_dir=experiments_dir or "experiments")
    exec_mgr = experiments.SurveyExecutionManager(store)

    if db_path:
        db_path_resolved = str(Path(db_path).resolve())
        found_run = False
        for r in store.list_runs():
            if r.get("output_db_path") and Path(r["output_db_path"]).resolve() == Path(db_path_resolved):
                store.set_active_selection(r["experiment_id"], r["run_id"])
                found_run = True
                break
        if not found_run:
            adhoc_exp = store.create_experiment(
                name=("Synthetic Demonstration" if Path(db_path).name == 'sample_survey_review.sqlite' else f"Survey: {Path(db_path).name}"),
                description=("Synthetic demonstration; not real astronomical detections" if Path(db_path).name == 'sample_survey_review.sqlite' else f"Ad-hoc research workspace for {Path(db_path).name}"),
                scientific_objective="Exploratory candidate review",
                targets=[],
                sectors=[],
                is_read_only=True
            )
            run_id = adhoc_exp["runs"][0]["run_id"]
            with sqlite3.connect(Path(db_path_resolved).as_uri() + '?mode=ro', uri=True) as source:
                columns = {r[1] for r in source.execute('PRAGMA table_info(runs)')}
                key = 'run_id' if 'run_id' in columns else 'id'
                native_ids = [str(r[0]) for r in source.execute(f'SELECT {key} FROM runs')] if columns else [None]
            for index, native_id in enumerate(native_ids):
                entry_id = run_id if index == 0 else store.create_run(adhoc_exp['id'], 'Source: ' + str(native_id))['run_id']
                overlay = review_db if len(native_ids) == 1 else None
                store.update_run(entry_id, output_db_path=db_path_resolved, review_db_path=overlay,
                                 source_run_id=native_id, status='completed')
            store.set_active_selection(adhoc_exp["id"], run_id)

    class CustomHandler(AnomalyReviewHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(db_path, root_path, *args, cache_dir=cache_path,
                             review_db=review_db, read_only=read_only, **kwargs)

    server = ThreadingHTTPServer((host, port), CustomHandler)
    server.experiment_store = store
    server.execution_manager = exec_mgr
    server.adapters_cache = {}
    server.initial_db_path = db_path
    url = f"http://{host}:{port}"

    act_exp, act_run = store.get_active_selection()
    cur_exp = store.get_experiment(act_exp) if act_exp else None

    print("=" * 72)
    print("* TESS ANOMALY SURVEY - RESEARCH INVESTIGATION WORKSPACE")
    print("=" * 72)
    print(f"Server URL   : {url}")
    if cur_exp:
        print(f"Experiment   : {cur_exp['name']} [{act_exp}]")
        print(f"Active Run   : {act_run}")
        print(f"Mode         : {'READ-ONLY HISTORICAL SOURCE' if cur_exp.get('is_read_only') else 'WRITABLE EXPERIMENT WORKSPACE'}")
    elif db_path:
        print(f"Database     : {db_path} {'[READ-ONLY]' if read_only else ''}")
    print(f"Light Curves : {root_path or 'On-Demand MAST Streaming & Cache'}")
    print(f"Cache Dir    : {cache_path}")
    print(f"Experiments  : {store.workspace_dir}")
    print("=" * 72)
    print("Serving requests on localhost. Press Ctrl+C to terminate.")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down workspace server.")
    finally:
        server.server_close()


# ---------------------------------------------------------------------------
# Single-Page Application (HTML / CSS / JS)
# ---------------------------------------------------------------------------

INDEX_HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>TESS Anomaly Survey — Review Workspace & Star Investigation</title>
<!-- Plotly.js for interactive time series analysis -->
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  :root {
    --bg-main: #0b0f17;
    --bg-card: #151d2c;
    --bg-card-hover: #1e293b;
    --border: #26334d;
    --border-subtle: #1e293b;
    --accent: #38bdf8;
    --accent-hover: #0284c7;
    --text-primary: #f8fafc;
    --text-secondary: #94a3b8;
    --text-muted: #64748b;
    --status-unreviewed: #64748b;
    --status-interesting: #10b981;
    --status-known: #3b82f6;
    --status-artifact: #f59e0b;
    --status-dismissed: #ef4444;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  html, body {
    width: 100%;
    max-width: 100vw;
    overflow-x: hidden;
  }
  body {
    background-color: var(--bg-main);
    color: var(--text-primary);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    font-size: 13px;
    line-height: 1.5;
    height: 100vh;
    display: flex;
    flex-direction: column;
  }
  header {
    background: #090d16;
    border-bottom: 1px solid var(--border);
    padding: 10px 16px;
    display: flex;
    flex-wrap: wrap;
    justify-content: space-between;
    align-items: center;
    gap: 10px;
    flex-shrink: 0;
  }
  header h1 { font-size: 1.15rem; font-weight: 700; color: var(--accent); display: flex; align-items: center; gap: 8px; flex-shrink: 0; }
  .nav-tabs { display: flex; background: #0b0f17; border: 1px solid var(--border); border-radius: 8px; padding: 2px; flex-wrap: wrap; }
  .nav-tab {
    padding: 5px 12px;
    border-radius: 6px;
    font-size: 12px;
    font-weight: 600;
    color: var(--text-secondary);
    background: transparent;
    border: none;
    cursor: pointer;
    transition: all 0.15s;
    white-space: nowrap;
  }
  .nav-tab.active { background: #0284c7; color: #fff; box-shadow: 0 1px 3px rgba(0,0,0,0.3); }

  .badge { padding: 2px 7px; border-radius: 4px; font-weight: 700; text-transform: uppercase; font-size: 10px; display: inline-block; white-space: nowrap; }
  .badge.unreviewed { background: rgba(100,116,139,0.25); color: #cbd5e1; border: 1px solid #64748b; }
  .badge.interesting { background: rgba(16,185,129,0.25); color: #34d399; border: 1px solid #10b981; }
  .badge.known { background: rgba(59,130,246,0.25); color: #60a5fa; border: 1px solid #3b82f6; }
  .badge.artifact { background: rgba(245,158,11,0.25); color: #fbbf24; border: 1px solid #f59e0b; }
  .badge.dismissed { background: rgba(239,68,68,0.25); color: #f87171; border: 1px solid #ef4444; }

  .prog-badge { padding: 2px 6px; border-radius: 4px; font-weight: 600; font-size: 10px; text-transform: uppercase; white-space: nowrap; }
  .prog-badge.unreviewed { color: #94a3b8; border: 1px solid #475569; background: rgba(71,85,105,0.2); }
  .prog-badge.partially_reviewed { color: #f59e0b; border: 1px solid #f59e0b; background: rgba(245,158,11,0.15); }
  .prog-badge.reviewed { color: #10b981; border: 1px solid #10b981; background: rgba(16,185,129,0.15); }

  .main-content {
    flex: 1;
    display: flex;
    flex-direction: column;
    overflow-y: auto;
    overflow-x: hidden;
    padding: 14px 16px;
    gap: 14px;
    min-width: 0;
  }

  /* Explore table styling */
  .toolbar {
    background: var(--bg-card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 10px 14px;
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: space-between;
    gap: 10px;
    min-width: 0;
  }
  .toolbar-group { display: flex; flex-wrap: wrap; align-items: center; gap: 8px; min-width: 0; }
  input[type="text"], select, textarea {
    background: #0b0f17;
    border: 1px solid var(--border);
    color: var(--text-primary);
    border-radius: 6px;
    padding: 6px 10px;
    font-size: 12px;
    outline: none;
    max-width: 100%;
  }
  input[type="text"]:focus, select:focus, textarea:focus { border-color: var(--accent); }

  .btn {
    background: var(--accent-hover);
    color: #fff;
    border: none;
    border-radius: 6px;
    padding: 6px 14px;
    font-size: 12px;
    font-weight: 600;
    cursor: pointer;
    transition: background 0.15s;
    display: inline-flex;
    align-items: center;
    gap: 5px;
    white-space: nowrap;
  }
  .btn:hover { background: #0369a1; }
  .btn-secondary { background: #1e293b; border: 1px solid var(--border); color: var(--text-primary); }
  .btn-secondary:hover { background: #334155; }
  .btn-sm { padding: 3px 8px; font-size: 11px; }

  .table-container {
    background: var(--bg-card);
    border: 1px solid var(--border);
    border-radius: 10px;
    overflow-x: auto;
    -webkit-overflow-scrolling: touch;
    width: 100%;
    min-width: 0;
    display: flex;
    flex-direction: column;
  }
  table.data-table {
    width: 100%;
    min-width: 760px;
    border-collapse: collapse;
    text-align: left;
  }
  table.data-table th {
    background: #090d16;
    padding: 10px 14px;
    font-size: 11px;
    font-weight: 700;
    text-transform: uppercase;
    color: var(--text-secondary);
    border-bottom: 1px solid var(--border);
    white-space: nowrap;
  }
  table.data-table td {
    padding: 10px 14px;
    border-bottom: 1px solid var(--border-subtle);
    font-size: 12px;
  }
  table.data-table tr.star-row:hover { background: rgba(56, 189, 248, 0.04); cursor: pointer; }
  table.data-table tr.star-row.active { background: rgba(56, 189, 248, 0.08); }

  .events-subpanel {
    background: #090d16;
    padding: 12px 16px 16px 42px;
    border-bottom: 1px solid var(--border);
  }

  /* Workspace Layout */
  .workspace-grid {
    display: grid;
    grid-template-columns: minmax(0, 1.8fr) minmax(0, 1.2fr);
    gap: 16px;
    flex: 1;
    min-width: 0;
  }
  @media (max-width: 1100px) {
    .workspace-grid { grid-template-columns: minmax(0, 1fr); }
  }

  .card {
    background: var(--bg-card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 14px;
    display: flex;
    flex-direction: column;
    gap: 12px;
    min-width: 0;
    word-break: break-word;
  }
  .card-header {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding-bottom: 10px;
    border-bottom: 1px solid var(--border);
    font-size: 12px;
    font-weight: 700;
    color: var(--text-primary);
    flex-wrap: wrap;
    gap: 8px;
  }

  .plot-container {
    width: 100%;
    height: 240px;
    background: #090d16;
    border-radius: 8px;
    border: 1px solid var(--border);
    position: relative;
    overflow: hidden;
    min-width: 0;
  }

  @media (max-width: 768px) {
    header { padding: 8px 12px; }
    header h1 { font-size: 1rem; }
    .main-content { padding: 10px; gap: 10px; }
    .toolbar { padding: 8px 10px; }
    .toolbar-group { width: 100%; }
    .toolbar-group input[type="text"], .toolbar-group select { width: 100%; }
    .plot-container { height: 210px; }
    .modal-dialog { width: 96vw; margin: 8px; }
  }
  @media (max-width: 480px) {
    header h1 span { display: none; }
    .nav-tabs { width: 100%; justify-content: space-between; }
    .nav-tab { padding: 5px 8px; font-size: 11px; }
  }

  /* Provenance tiers */
  .prov-box {
    background: #090d16;
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 8px 12px;
    display: flex;
    flex-direction: column;
    gap: 3px;
  }
  .prov-title { font-size: 10px; text-transform: uppercase; font-weight: 700; color: var(--text-muted); }

  /* Toast Notification */
  #toast-container {
    position: fixed;
    bottom: 20px;
    right: 20px;
    display: flex;
    flex-direction: column;
    gap: 8px;
    z-index: 9999;
  }
  .toast {
    background: #047857;
    color: #fff;
    padding: 8px 16px;
    border-radius: 6px;
    box-shadow: 0 4px 12px rgba(0,0,0,0.5);
    font-size: 12px;
    font-weight: 600;
    display: flex;
    align-items: center;
    gap: 8px;
    animation: slideIn 0.2s ease-out;
  }
  @keyframes slideIn { from { transform: translateY(20px); opacity: 0; } to { transform: translateY(0); opacity: 1; } }

  /* Modal styling */
  .modal-backdrop {
    position: fixed; top: 0; left: 0; right: 0; bottom: 0;
    background: rgba(0,0,0,0.75); z-index: 1000;
    display: flex; align-items: center; justify-content: center;
  }
  .modal-dialog {
    background: #111827; border: 1px solid var(--border);
    border-radius: 12px; width: 720px; max-width: 95vw;
    max-height: 90vh; overflow-y: auto; display: flex; flex-direction: column;
    box-shadow: 0 20px 25px -5px rgba(0, 0, 0, 0.5);
  }
  .modal-header {
    padding: 14px 20px; border-bottom: 1px solid var(--border);
    display: flex; justify-content: space-between; align-items: center;
  }
  .modal-body { padding: 20px; display: flex; flex-direction: column; gap: 14px; }
  .modal-footer {
    padding: 14px 20px; border-top: 1px solid var(--border);
    display: flex; justify-content: flex-end; gap: 10px; background: #090d16;
  }
  .preset-card {
    border: 1px solid var(--border); border-radius: 8px; padding: 12px 14px;
    background: #151d2c; cursor: pointer; transition: border-color 0.15s;
  }
  .preset-card.selected { border-color: var(--accent); background: #1a2538; }
  .progress-bar-container {
    width: 100%; background: #1e293b; border-radius: 6px; height: 12px; overflow: hidden;
  }
  .progress-bar-fill {
    height: 100%; background: #0284c7; width: 0%; transition: width 0.3s ease;
  }
  .pulse-indicator {
    display: inline-block; width: 8px; height: 8px; border-radius: 50%;
    background: #10b981; margin-right: 6px; box-shadow: 0 0 8px #10b981;
    animation: pulse 1.5s infinite;
  }
  @keyframes pulse { 0% { opacity: 0.4; } 50% { opacity: 1; } 100% { opacity: 0.4; } }
</style>
</head>
<body>

  <!-- Top App Bar with Experiment Selection -->
  <header>
    <div style="display:flex; align-items:center; gap:12px; flex-wrap:wrap; min-width:0; max-width:100%;">
      <h1>★ TESS Anomaly Survey <span style="font-size:0.75rem; color:var(--text-secondary); font-weight:normal;">Research Workspace</span></h1>
      
      <!-- Experiment / Run Selector -->
      <div style="display:flex; align-items:center; gap:8px; flex-wrap:wrap; min-width:0; max-width:100%;">
        <span style="font-size:11px; font-weight:700; color:var(--text-secondary); letter-spacing:0.5px;">EXPERIMENT:</span>
        <select id="sel-experiment" onchange="onSwitchExperiment(this.value)" style="background:var(--bg-card); color:#fff; border:1px solid var(--border); border-radius:6px; padding:4px 10px; font-size:12px; font-weight:600; min-width:180px; max-width:280px;"></select>
        <select id="sel-run" aria-label="Run" onchange="onSwitchExperiment(AppState.activeExperiment.id, this.value)" style="max-width:280px; background:var(--bg-card); color:#fff;"></select>
        <span id="badge-exp-mode" class="badge unreviewed">CONNECTING</span>
        <button class="btn btn-sm" onclick="openNewExperimentModal()" style="background:var(--accent); color:#0b0f17; font-weight:700; padding:4px 10px; font-size:11px; border:none; border-radius:6px; cursor:pointer;">+ New Experiment</button>
      </div>
    </div>

    <div style="display:flex; align-items:center; gap:12px; flex-wrap:wrap; min-width:0; max-width:100%;">
      <div class="nav-tabs">
        <button id="tab-btn-explore" class="nav-tab active" onclick="switchView('explore')">Explore Stars</button>
        <button id="tab-btn-workspace" class="nav-tab" onclick="switchView('workspace')">Investigation Workspace</button>
        <button id="tab-btn-progress" class="nav-tab" onclick="switchView('progress')">Survey Progress & Runs</button>
      </div>
      <div style="display:flex; align-items:center; gap:8px; font-size:12px; color:var(--text-secondary); flex-wrap:wrap;">
        <span id="app-db-label" style="font-family:monospace; background:rgba(255,255,255,0.05); padding:3px 8px; border-radius:4px; border:1px solid var(--border);">DB: Connected</span>
        <span id="app-save-mode" class="badge" style="background:#064e3b; color:#34d399; border:1px solid #059669;">Verified Localhost</span>
      </div>
    </div>
  </header>

  <!-- Notification Toasts -->
  <div id="toast-container"></div>

  <!-- VIEW 1: EXPLORE STARS -->
  <main id="view-explore" class="main-content">
    <!-- Active Run Status & Scientific Caveat Notice -->
    <div id="explore-running-banner" style="display:none; background:rgba(2,132,199,0.12); border:1px solid rgba(2,132,199,0.4); border-radius:8px; padding:8px 14px; font-size:12px; align-items:center; justify-content:space-between;">
      <div style="display:flex; align-items:center;">
        <span class="pulse-indicator"></span>
        <span>A survey run is currently executing in background (<strong id="banner-run-progress">0/0 obs</strong>). Intermediate results are visible below.</span>
      </div>
      <button class="btn btn-sm" style="font-size:11px; padding:2px 8px;" onclick="switchView('progress')">Monitor Survey Execution →</button>
    </div>

    <div style="background:rgba(255,255,255,0.02); border:1px solid var(--border-subtle); border-radius:6px; padding:6px 12px; font-size:11px; color:var(--text-secondary); display:flex; justify-content:space-between; align-items:center;">
      <span><strong>Scientific Caveat:</strong> Matched-filter detector score S quantifies local statistical flux deficit prominence; high scores indicate photometric prominence, not necessarily astrophysical rarity or planet transit probability.</span>
      <span id="explore-active-run-label" style="font-family:monospace; color:var(--accent);"></span>
    </div>

    <div class="toolbar">
      <div class="toolbar-group">
        <input type="text" id="filter-search" placeholder="Search TIC, Notes, Sector..." style="min-width:220px;" onkeyup="if(event.key==='Enter') applyFilters()">
        <select id="filter-status" onchange="applyFilters()">
          <option value="all">All Star Review Statuses</option>
          <option value="unreviewed">Unreviewed Only</option>
          <option value="interesting">Interesting</option>
          <option value="known">Known Systems</option>
          <option value="artifact">Artifacts</option>
          <option value="dismissed">Dismissed</option>
        </select>
        <select id="filter-verification" onchange="applyFilters()">
          <option value="all">All Photometric Classes</option>
          <option value="dimming">Verified Localized Dimmings</option>
            <option value="not_run">Verification not run / unavailable</option>
          <option value="variability">Continuous Variability</option>
          <option value="flare">Flare-Induced Dips</option>
          <option value="artifact">Instrumental Artifacts</option>
        </select>
        <select id="filter-system" onchange="applyFilters()">
          <option value="all">All Star Types</option>
          <option value="periodic">Periodic Systems (N ≥ 3)</option>
          <option value="isolated">Isolated / Single Dips</option>
        </select>
        <select id="filter-sort" onchange="applyFilters()">
          <option value="score">Sort by Max Score (High → Low)</option>
          <option value="events">Sort by Event Count</option>
          <option value="sector">Sort by Sector</option>
          <option value="tic">Sort by TIC ID</option>
          <option value="period">Sort by Period</option>
        </select>
        <button class="btn" onclick="applyFilters()">Filter</button>
        <button class="btn btn-secondary" onclick="resetFilters()">Reset</button>
      </div>
      <div style="display:flex; align-items:center; gap:12px; font-size:12px; color:var(--text-secondary);">
        <span id="explore-count-label">Loading stars...</span>
      </div>
    </div>

    <div class="table-container">
      <table class="data-table">
        <thead>
          <tr>
            <th style="width:30px;"></th>
            <th>Target Star (TIC)</th>
            <th>Sector(s)</th>
            <th style="text-align:center;">Events</th>
            <th>Max Score (SNR)</th>
            <th>Photometric Verification</th>
            <th>Periodicity / System</th>
            <th>Review Status</th>
            <th style="text-align:right;">Actions</th>
          </tr>
        </thead>
        <tbody id="stars-table-body">
          <tr><td colspan="9" style="text-align:center; padding:30px; color:var(--text-secondary);">Loading anomaly library stars...</td></tr>
        </tbody>
      </table>
      <div style="display:flex; justify-content:space-between; align-items:center; padding:12px 16px; background:#090d16; border-top:1px solid var(--border);">
        <span id="explore-pagination-info" style="color:var(--text-secondary); font-size:12px;">Page 1</span>
        <div style="display:flex; gap:8px;">
          <button id="btn-prev-page" class="btn btn-secondary btn-sm" onclick="prevPage()" disabled>Previous</button>
          <button id="btn-next-page" class="btn btn-secondary btn-sm" onclick="nextPage()">Next</button>
        </div>
      </div>
    </div>
  </main>

  <!-- VIEW 2: INVESTIGATION WORKSPACE -->
  <main id="view-workspace" class="main-content" style="display:none;">
    <!-- Star Header Bar with Preserved Filter Navigation -->
    <div style="background:var(--bg-card); border:1px solid var(--border); border-radius:10px; padding:12px 18px; display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:12px;">
      <div style="display:flex; align-items:center; gap:14px;">
        <button class="btn btn-secondary btn-sm" onclick="switchView('explore')">← Back to Explore</button>
        <div>
          <div style="display:flex; align-items:center; gap:8px;">
            <span id="ws-tic-label" style="font-size:1.1rem; font-weight:700; color:#fff; font-family:monospace;">TIC —</span>
            <span id="ws-star-class-badge" class="badge unreviewed">UNREVIEWED</span>
            <span id="ws-benchmark-badge" class="badge" style="background:#064e3b; color:#34d399; border:1px solid #059669; display:none;">BENCHMARK</span>
          </div>
          <div id="ws-star-subinfo" style="font-size:11px; color:var(--text-secondary); display:flex; gap:10px; flex-wrap:wrap; margin-top:2px;">
            <span>Sector: <b id="ws-sector-val" style="color:#fff;">—</b></span>
            <span>Tmag: <b id="ws-tmag-val" style="color:#fff;">—</b></span>
            <span>Teff: <b id="ws-teff-val" style="color:#fff;">—</b></span>
            <span>Total Events: <b id="ws-events-count-val" style="color:var(--accent);">—</b></span>
          </div>
        </div>
      </div>

      <!-- Preserved Navigation across Filter List -->
      <div style="display:flex; align-items:center; gap:8px; flex-wrap:wrap;">
        <button id="ws-btn-prev-star" class="btn btn-secondary btn-sm" onclick="navigateStar(-1)">← Previous Star</button>
        <span id="ws-star-index-label" style="font-size:11px; color:var(--text-secondary); font-family:monospace;">Star 1 of 1</span>
        <button id="ws-btn-next-star" class="btn btn-secondary btn-sm" onclick="navigateStar(1)">Next Star →</button>
      </div>
    </div>

    <!-- Main Workspace Grid -->
    <div class="workspace-grid">
      <!-- Left Column: Synchronized Interactive Plots -->
      <div style="display:flex; flex-direction:column; gap:14px; min-width:0;">

        <!-- Plot 1: Full Observation View -->
        <div class="card">
          <div class="card-header">
            <span>FULL OBSERVATION LIGHT CURVE (~27.4 DAYS SPAN)</span>
            <div style="display:flex; gap:8px; align-items:center; flex-wrap:wrap;">
              <span id="ws-full-lc-source" style="font-size:10px; color:var(--text-muted); font-family:monospace;">Source: Local</span>
              <button class="btn btn-secondary btn-sm" onclick="togglePlotMode()">Switch to Static PNG</button>
            </div>
          </div>
          <div id="plot-full-container" class="plot-container"></div>
        </div>

        <!-- Plot 2: Event Zoom View -->
        <div class="card">
          <div class="card-header">
            <span>EVENT ZOOM VIEW (±12 HOURS)</span>
            <div style="display:flex; gap:8px; align-items:center; flex-wrap:wrap;">
              <!-- Event selector buttons -->
              <div id="ws-event-tabs" style="display:flex; gap:4px; flex-wrap:wrap;"></div>
              <select id="ws-zoom-span" onchange="updateEventZoom()" style="font-size:11px; padding:2px 6px;">
                <option value="12">±12 h</option>
                <option value="6">±6 h</option>
                <option value="24">±24 h</option>
              </select>
            </div>
          </div>
          <div id="plot-zoom-container" class="plot-container"></div>
        </div>

        <!-- Plot 3: Phase Folded Curve (if periodic) -->
        <div id="card-phase-fold" class="card" style="display:none;">
          <div class="card-header">
            <span>PHASE-FOLDED LIGHT CURVE</span>
            <span id="ws-phase-info" style="font-size:11px; color:var(--accent); font-family:monospace;">P = — d</span>
          </div>
          <div id="plot-phase-container" class="plot-container"></div>
        </div>

      </div>

      <!-- Right Column: Multi-Domain Provenance & Dual Review Panel -->
      <div style="display:flex; flex-direction:column; gap:14px;">

        <!-- 5-Tier Scientific Provenance Card -->
        <div class="card">
          <div class="card-header">
            <span>SCIENTIFIC OBSERVATIONAL EVIDENCE</span>
            <span class="badge" style="background:rgba(56,189,248,0.15); color:var(--accent); border:1px solid var(--accent);">5 TIERS</span>
          </div>
          
          <div class="prov-box">
            <span class="prov-title">1. Matched-Filter Detector Metric</span>
            <div style="display:flex; justify-content:space-between; font-size:12px;">
              <span>Score: <b id="prov-score" style="color:#fbbf24; font-family:monospace;">—</b></span>
              <span style="color:var(--text-muted); font-family:monospace;">Execution threshold shown below</span>
            </div>
            <span id="prov-kernel-info" style="font-size:11px; color:var(--text-secondary);">Select a target star to evaluate detector significance.</span>
          </div>

          <div class="prov-box">
            <span class="prov-title">2. Automated Photometric Verification (5-State Machine)</span>
            <div id="prov-verif-class" style="font-weight:700; color:#34d399; font-size:12px;">—</div>
            <div style="display:flex; justify-content:space-between; font-size:11px; color:var(--text-secondary); margin-top:2px;">
              <span>Depth measurements: <b id="prov-depth" style="color:#fff;">—</b></span>
              <span>Flare Shoulder: <b id="prov-flare" style="color:#34d399;">—</b></span>
            </div>
            <span id="prov-verif-desc" style="font-size:11px; color:var(--text-secondary); margin-top:2px;">Awaiting target selection.</span>
          </div>

          <div class="prov-box">
            <span class="prov-title">3. Cross-Star Coincidence Risk</span>
            <div id="prov-coincidence-level" style="font-weight:700; color:#34d399; font-size:12px;">—</div>
            <span id="prov-coincidence-desc" style="font-size:11px; color:var(--text-secondary);">No target star or event selected.</span>
          </div>

          <div class="prov-box">
            <span class="prov-title">4. Provisional Morphological Interpretation</span>
            <div id="prov-provisional-class" style="font-weight:700; color:var(--accent); font-size:12px;">—</div>
            <span id="prov-provisional-desc" style="font-size:11px; color:var(--text-secondary);">Select a star and candidate event to inspect scientific evidence.</span>
          </div>
        </div>

        <!-- Dual-Scope Review Panel -->
        <div class="card" style="border-color:#0284c7;">
          <div class="card-header" style="color:var(--accent);">
            <span>5. RESEARCHER REVIEW & AUDIT</span>
            <span class="badge" style="background:#064e3b; color:#34d399; border:1px solid #059669;">ISOLATED SCOPES</span>
          </div>

          <!-- Star-Level Review Form -->
          <div style="background:#090d16; border:1px solid var(--border); border-radius:8px; padding:12px; display:flex; flex-direction:column; gap:8px;">
            <div style="display:flex; justify-content:space-between; align-items:center;">
              <span style="font-weight:700; font-size:12px; color:#fff;">Star-Level Review</span>
              <span id="star-rev-status-badge" class="badge unreviewed">unreviewed</span>
            </div>
            <div style="display:grid; grid-template-columns:1fr 1fr; gap:8px;">
              <select id="star-rev-status">
                <option value="unreviewed">Status: unreviewed</option>
                <option value="interesting">Status: interesting</option>
                <option value="known">Status: known</option>
                <option value="artifact">Status: artifact</option>
                <option value="dismissed">Status: dismissed</option>
              </select>
              <select id="star-rev-class">
                <option value="">Class: (Select Taxonomy)</option>
                <option value="candidate_transit">candidate_transit</option>
                <option value="eclipsing_binary">eclipsing_binary</option>
                <option value="pulsator">pulsator</option>
                <option value="disk_occultation">disk_occultation</option>
                <option value="stellar_variability">stellar_variability</option>
                <option value="instrumental_artifact">instrumental_artifact</option>
                <option value="unresolved">unresolved</option>
              </select>
            </div>
            <textarea id="star-rev-note" placeholder="Star-level research notes, literature crossmatches, or astrophysical comments..." style="height:55px; resize:vertical;"></textarea>
            <div style="display:flex; justify-content:space-between; align-items:center;">
              <label style="display:flex; align-items:center; gap:6px; font-size:11px; color:var(--text-secondary); cursor:pointer;">
                <input type="checkbox" id="star-bulk-apply"> Bulk apply to unreviewed events
              </label>
              <button class="btn btn-sm" onclick="saveStarReview()">Save Star Review</button>
            </div>
          </div>

          <!-- Event-Level Review Form -->
          <div style="background:#090d16; border:1px solid var(--border); border-radius:8px; padding:12px; display:flex; flex-direction:column; gap:8px;">
            <div style="display:flex; justify-content:space-between; align-items:center;">
              <span style="font-weight:700; font-size:12px; color:#fff;">Event-Level Review <span id="event-rev-id-label" style="color:var(--accent);">#—</span></span>
              <span id="event-rev-status-badge" class="badge unreviewed">unreviewed</span>
            </div>
            <select id="event-rev-status">
              <option value="unreviewed">Status: unreviewed</option>
              <option value="interesting">Status: interesting</option>
              <option value="known">Status: known</option>
              <option value="artifact">Status: artifact</option>
              <option value="dismissed">Status: dismissed</option>
            </select>
            <textarea id="event-rev-note" placeholder="Event-specific notes (e.g. centroid offset, asymmetry, timing)..." style="height:48px; resize:vertical;"></textarea>
            <div style="display:flex; justify-content:flex-end;">
              <button class="btn btn-sm" onclick="saveEventReview()">Save Event Review</button>
            </div>
          </div>

          <!-- Safe Undo and Audit History -->
          <div style="background:#090d16; border:1px solid var(--border); border-radius:8px; padding:10px; display:flex; flex-direction:column; gap:6px;">
            <div style="display:flex; justify-content:space-between; align-items:center; font-size:11px; font-weight:700; color:var(--text-secondary);">
              <span>AUDIT TRAIL & SAFE UNDO</span>
              <div style="display:flex; gap:6px;">
                <button class="btn btn-secondary btn-sm" style="font-size:10px; padding:2px 6px;" onclick="undoStarReview()">↺ Undo Star</button>
                <button class="btn btn-secondary btn-sm" style="font-size:10px; padding:2px 6px;" onclick="undoEventReview()">↺ Undo Event</button>
              </div>
            </div>
            <div id="ws-audit-history-list" style="display:flex; flex-direction:column; gap:4px; max-height:80px; overflow-y:auto; font-size:11px; color:var(--text-secondary);">
              <span>No prior modifications recorded.</span>
            </div>
          </div>

        </div>

      </div>
    </div>
  </main>

  <!-- VIEW 3: SURVEY EXECUTION & TELEMETRY -->
  <main id="view-progress" class="main-content" style="display:none;">
    <div class="card">
      <div class="card-header" style="display:flex; justify-content:space-between; align-items:center;">
        <div style="display:flex; align-items:center; gap:10px;">
          <span id="progress-header-title" style="font-weight:700; font-size:13px; letter-spacing:0.5px;">SURVEY EXECUTION MONITOR</span>
          <span id="progress-status-badge" class="badge">HISTORICAL BENCHMARK</span>
        </div>
        <div style="display:flex; align-items:center; gap:8px;">
          <button id="btn-export-csv" class="btn btn-sm btn-secondary" onclick="exportCatalog('csv')">Export Catalog (CSV)</button>
          <button id="btn-export-json" class="btn btn-sm btn-secondary" onclick="exportCatalog('json')">Export Provenance (JSON)</button>
          <button id="btn-abort-survey" class="btn btn-sm" style="display:none; background:#ef4444; color:#fff;" onclick="abortSurveyRun()">Abort Survey</button>
          <button id="btn-resume-survey" class="btn btn-sm" style="display:none; background:var(--accent); color:#0b0f17; font-weight:700;" onclick="resumeSurveyRun()">Resume Incomplete Run</button>
        </div>
      </div>

      <!-- Live Execution Progress Bar (Shown when active run is running) -->
      <div id="live-progress-container" style="display:none; margin-top:8px; flex-direction:column; gap:8px;">
        <div style="display:flex; justify-content:space-between; align-items:center; font-size:12px; color:var(--text-secondary);">
          <span style="display:flex; align-items:center;"><span class="pulse-indicator"></span> Streaming Survey Running: <strong id="live-proc-count" style="margin-left:4px; color:#fff;">0</strong> / <strong id="live-total-count" style="margin-right:4px;">0</strong> observations</span>
          <span id="live-percent-label" style="font-weight:700; color:var(--accent); font-family:monospace;">0%</span>
        </div>
        <div class="progress-bar-container">
          <div id="live-progress-bar" class="progress-bar-fill"></div>
        </div>
        <div style="display:flex; justify-content:flex-end;">
          <button class="btn btn-sm btn-secondary" style="font-size:11px;" onclick="switchView('explore')">View Intermediate Results in Explore →</button>
        </div>
      </div>

      <!-- Metric Cards -->
      <div style="display:grid; grid-template-columns:repeat(auto-fit, minmax(180px, 1fr)); gap:14px; margin-top:10px;">
        <div class="prov-box">
          <span class="prov-title">Observations</span>
          <span id="tel-obs-count" style="font-size:1.5rem; font-weight:700; color:#fff; font-family:monospace;">—</span>
          <span id="tel-obs-sub" style="font-size:11px; color:#34d399;">Completed / Ingested</span>
        </div>
        <div class="prov-box">
          <span class="prov-title">Unique Stars</span>
          <span id="tel-stars-count" style="font-size:1.5rem; font-weight:700; color:#fff; font-family:monospace;">—</span>
          <span id="tel-stars-sub" style="font-size:11px; color:var(--text-secondary);">Targets</span>
        </div>
        <div class="prov-box">
          <span class="prov-title">Candidates Detected</span>
          <span id="tel-candidates-count" style="font-size:1.5rem; font-weight:700; color:#f59e0b; font-family:monospace;">—</span>
          <span id="tel-verified-sub" style="font-size:11px; color:var(--text-secondary);">Anomaly alerts</span>
        </div>
        <div class="prov-box">
          <span class="prov-title">Pipeline Throughput</span>
          <span id="tel-throughput" style="font-size:1.5rem; font-weight:700; color:var(--accent); font-family:monospace;">—</span>
          <span id="tel-throughput-sub" style="font-size:11px; color:var(--text-secondary);">Obs / hour</span>
        </div>
        <div class="prov-box">
          <span class="prov-title">Elapsed / ETA</span>
          <span id="tel-eta" style="font-size:1.5rem; font-weight:700; color:#cbd5e1; font-family:monospace;">—</span>
          <span id="tel-eta-sub" style="font-size:11px; color:var(--text-secondary);">Time remaining</span>
        </div>
      </div>

      <!-- Provenance & Scientific Principles -->
      <div style="margin-top:14px; background:#0b0f17; border:1px solid var(--border); border-radius:8px; padding:12px 16px;">
        <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:8px;">
          <span style="font-weight:700; font-size:11px; color:var(--accent); text-transform:uppercase; letter-spacing:0.5px;">Configuration Provenance</span>
          <span id="provenance-hash" style="font-family:monospace; font-size:10px; color:var(--text-muted);"></span>
        </div>
        <p id="provenance-objective" style="color:#cbd5e1; font-size:12px; margin-bottom:6px;"></p>
        <p id="provenance-config" style="color:var(--text-secondary); font-size:11px; margin-bottom:8px;"></p>
        <div style="background:rgba(2,132,199,0.08); border-left:3px solid var(--accent); padding:8px 12px; border-radius:4px; font-size:11px; color:#cbd5e1;">
          <strong>Scientific Principle:</strong> Matched-filter detector score S reflects localized statistical flux deficit significance under the 2h template. High scores indicate high photometric prominence, not necessarily astrophysical rarity or scientific novelty.
        </div>
      </div>

      <!-- Activity Logs Feed -->
      <div style="margin-top:14px;">
        <div style="font-weight:600; font-size:11px; color:var(--text-secondary); margin-bottom:6px; text-transform:uppercase;">Execution Activity Logs</div>
        <div id="live-logs-box" style="background:#090d16; border:1px solid var(--border); border-radius:6px; padding:10px; max-height:150px; overflow-y:auto; font-family:monospace; font-size:11px; color:#94a3b8; display:flex; flex-direction:column; gap:4px;">
          <span style="color:var(--text-muted);">No active logs recorded.</span>
        </div>
      </div>
    </div>
  </main>

  <!-- MODAL: NEW SURVEY EXPERIMENT -->
  <div id="modal-new-experiment" class="modal-backdrop" style="display:none;" onclick="if(event.target===this) closeNewExperimentModal();">
    <div class="modal-dialog">
      <div class="modal-header">
        <h3 style="font-size:14px; font-weight:700; color:var(--accent); display:flex; align-items:center; gap:8px;">
          <span>NEW SURVEY EXPERIMENT</span>
        </h3>
        <button onclick="closeNewExperimentModal()" style="background:transparent; border:none; color:var(--text-secondary); cursor:pointer; font-size:20px; line-height:1;">&times;</button>
      </div>
      <div class="modal-body">
        <!-- Step 1: Scientific Identity -->
        <div>
          <label style="display:block; font-weight:600; margin-bottom:4px; font-size:12px;">Experiment Name *</label>
          <input type="text" id="new-exp-name" placeholder="e.g. Sector 8 Disk Occultation Search" style="width:100%;">
        </div>
        <div>
          <label style="display:block; font-weight:600; margin-bottom:4px; font-size:12px;">Scientific Objective *</label>
          <input type="text" id="new-exp-objective" placeholder="e.g. Targeted search for localized negative flux excursions in young stars" style="width:100%;">
        </div>
        <div>
          <label style="display:block; font-weight:600; margin-bottom:4px; font-size:12px;">Description / Context</label>
          <textarea id="new-exp-desc" rows="2" placeholder="Cohort criteria, literature references, or experimental hypotheses..." style="width:100%;"></textarea>
        </div>

        <!-- Step 2: Target & Sector Selection -->
        <div style="border-top:1px solid var(--border); padding-top:12px;">
          <label style="display:block; font-weight:600; margin-bottom:4px; font-size:12px;">Target TICs (Identifiers or CSV) *</label>
          <textarea id="new-exp-targets" rows="3" placeholder="Enter TICs separated by commas, spaces, or newlines (e.g. 39018208, 73834, 60232) or paste CSV..." style="width:100%; font-family:monospace;"></textarea>
          
          <div style="margin-top:8px;">
            <label style="display:block; font-weight:600; margin-bottom:6px; font-size:12px;">Observing Sectors (Select at least one):</label>
            <div id="sector-checkboxes" style="display:flex; flex-wrap:wrap; gap:8px;">
              <!-- 1 to 13 checkboxes populated by JS -->
            </div>
          </div>

          <div style="margin-top:10px; display:flex; align-items:center; gap:10px;">
            <button type="button" class="btn btn-sm btn-secondary" onclick="runTargetValidation()">Validate Targets & Estimate</button>
            <span id="target-val-status" style="font-size:11px; color:var(--text-secondary);"></span>
          </div>

          <!-- Pre-flight Estimate Panel -->
          <div id="preflight-box" style="display:none; margin-top:10px; background:#0b0f17; border:1px solid var(--border); border-radius:6px; padding:10px; font-size:11px;">
            <div style="display:grid; grid-template-columns:repeat(3, 1fr); gap:8px; text-align:center;">
              <div><strong id="pf-stars" style="color:#fff; font-size:13px;">0</strong><br><span style="color:var(--text-muted);">Unique Stars</span></div>
              <div><strong id="pf-obs" style="color:var(--accent); font-size:13px;">0</strong><br><span style="color:var(--text-muted);">Total Observations</span></div>
              <div><strong id="pf-mb" style="color:#34d399; font-size:13px;">0</strong> MB<br><span style="color:var(--text-muted);">Expected Download</span></div>
            </div>
            <div id="pf-warnings" style="margin-top:8px; color:#fbbf24; font-size:11px;"></div>
          </div>
        </div>

        <!-- Step 3: Detection Configuration -->
        <div style="border-top:1px solid var(--border); padding-top:12px;">
          <label style="display:block; font-weight:600; margin-bottom:6px; font-size:12px;">Detection Configuration</label>
          <div style="display:flex; flex-direction:column; gap:8px;">
            <div id="preset-card-frozen" class="preset-card selected" onclick="selectPreset('frozen')">
              <div style="display:flex; justify-content:space-between; align-items:center;">
                <strong style="color:#fff;">Frozen Fixed Box 2h (Recommended)</strong>
                <span class="badge" style="background:#064e3b; color:#34d399;">PUBLISHED PRESET</span>
              </div>
              <div style="font-size:11px; color:var(--text-secondary); margin-top:4px;">
                Threshold T = 73.7797 | Duration = 2.0h | B0/B1 SHA-256 Verified
              </div>
              <div style="font-size:10px; color:var(--text-muted); margin-top:4px;">
                Calibrated exclusively for 2-hour localized negative flux excursions. Does not perform transit fitting.
              </div>
            </div>

            <div id="preset-card-custom" class="preset-card" onclick="selectPreset('custom')">
              <div style="display:flex; justify-content:space-between; align-items:center;">
                <strong style="color:#fff;">Experimental Custom Configuration</strong>
                <span class="badge" style="background:#78350f; color:#fbbf24;">EXPERIMENTAL</span>
              </div>
              <div id="custom-cfg-inputs" style="display:none; margin-top:8px; gap:10px; grid-template-columns:1fr 1fr;">
                <div>
                  <label style="font-size:11px; color:var(--text-secondary);">Custom Threshold T:</label>
                  <input type="number" id="custom-threshold" step="0.1" value="73.78" style="width:100%;">
                </div>
                <div>
                  <label style="font-size:11px; color:var(--text-secondary);">Custom Box Duration (hours):</label>
                  <input type="number" id="custom-duration" readonly title="Only the validated two-hour template is supported" step="0.5" value="2.0" style="width:100%;">
                </div>
              </div>
            </div>
          </div>
        </div>

        <!-- Step 4: Execution Settings -->
        <div style="border-top:1px solid var(--border); padding-top:12px; display:flex; justify-content:space-between; align-items:center;">
          <div>
            <label style="font-weight:600; font-size:12px;">Parallel Workers:</label>
            <select id="new-exp-workers" style="margin-left:6px; padding:3px 8px;">
              <option value="1">1 (Conservative)</option>
              <option value="2" selected>2 (Recommended)</option>
              <option value="4">4 (High Throughput)</option>
            </select>
          </div>
          <div style="display:flex; align-items:center; gap:6px;">
            <input type="checkbox" id="new-exp-autolaunch" checked>
            <label for="new-exp-autolaunch" style="font-size:12px; cursor:pointer; color:#fff;">Launch Survey Immediately</label>
          </div>
        </div>
      </div>
      <div class="modal-footer">
        <button class="btn btn-secondary" onclick="closeNewExperimentModal()">Cancel</button>
        <button id="btn-create-experiment" class="btn" style="background:var(--accent); color:#0b0f17; font-weight:700;" onclick="submitCreateExperiment()">Create & Launch Experiment</button>
      </div>
    </div>
  </div>

  <script>
    // Global Application State
    const AppState = {
      currentView: 'explore',
      stars: [],
      filteredStars: [],
      currentStarIndex: 0,
      activeStar: null,
      activeEvent: null,
      activeLightCurve: null,
      currentPage: 1,
      pageSize: 25,
      expandedRows: new Set(),
      useStaticPlots: false,
      experiments: [],
      activeExperiment: null,
      activeRun: null,
      pollingTimer: null,
      selectedPreset: 'frozen'
    };

    // Navigation and Routing
    function switchView(viewName) {
      AppState.currentView = viewName;
      document.getElementById('view-explore').style.display = (viewName === 'explore') ? 'flex' : 'none';
      document.getElementById('view-workspace').style.display = (viewName === 'workspace') ? 'flex' : 'none';
      document.getElementById('view-progress').style.display = (viewName === 'progress') ? 'flex' : 'none';

      document.getElementById('tab-btn-explore').className = 'nav-tab' + (viewName === 'explore' ? ' active' : '');
      document.getElementById('tab-btn-workspace').className = 'nav-tab' + (viewName === 'workspace' ? ' active' : '');
      document.getElementById('tab-btn-progress').className = 'nav-tab' + (viewName === 'progress' ? ' active' : '');

      if (viewName === 'explore') {
        history.replaceState(null, '', '#view=explore');
      } else if (viewName === 'workspace' && AppState.activeStar) {
        history.replaceState(null, '', `#view=workspace&tic=${AppState.activeStar.tic}`);
        renderWorkspacePlots();
      } else if (viewName === 'progress') {
        history.replaceState(null, '', '#view=progress');
        checkRunStatus();
      }
    }

    // Initialize Sector Checkboxes
    function initSectorCheckboxes() {
      const container = document.getElementById('sector-checkboxes');
      if (!container || container.children.length > 0) return;
      for (let s = 1; s <= 13; s++) {
        const lbl = document.createElement('label');
        lbl.style.cssText = 'font-size:11px; display:inline-flex; align-items:center; gap:3px; background:#1e293b; padding:2px 6px; border-radius:4px; cursor:pointer;';
        lbl.innerHTML = `<input type="checkbox" name="sector-cb" value="${s}" ${s === 8 ? 'checked' : ''}> Sec ${s}`;
        container.appendChild(lbl);
      }
    }

    const nativeFetch = window.fetch.bind(window);
    window.fetch = (resource, options = {}) => {
      const url = new URL(resource, window.location.href);
      if (url.origin === location.origin && url.pathname.startsWith('/api/')) {
        if (AppState.activeRun) url.searchParams.set('workspace_run_id', AppState.activeRun.run_id);
      }
      options.headers = { ...options.headers, 'X-CSRF-Token': '__CSRF_TOKEN__' };
      return nativeFetch(url, options);
    };

    // Load Experiments list and active experiment
    async function loadExperiments() {
      try {
        const resp = await fetch('/api/experiments');
        const data = await resp.json();
        AppState.experiments = data.experiments || [];
        AppState.activeExperiment = data.active_experiment || null;
        AppState.activeRun = data.active_run || null;

        const select = document.getElementById('sel-experiment');
        select.innerHTML = '';

        const userGroup = document.createElement('optgroup');
        userGroup.label = 'Active & User Experiments';
        const histGroup = document.createElement('optgroup');
        histGroup.label = 'Historical Research Benchmarks (Read-Only)';

        AppState.experiments.forEach(exp => {
          const opt = document.createElement('option');
          opt.value = exp.id;
          opt.innerText = exp.name + (exp.is_historical ? ' [Historical]' : '');
          if (AppState.activeExperiment && exp.id === AppState.activeExperiment.id) {
            opt.selected = true;
          }
          if (exp.is_historical) {
            histGroup.appendChild(opt);
          } else {
            userGroup.appendChild(opt);
          }
        });

        if (userGroup.children.length > 0) select.appendChild(userGroup);
        if (histGroup.children.length > 0) select.appendChild(histGroup);

        updateExperimentBadge();
        updateProvenanceCard();
        checkRunStatus();
      } catch (e) {
        console.error("Failed to load experiments", e);
      }
    }

    function updateExperimentBadge() {
      const runs = document.getElementById('sel-run');
      runs.replaceChildren();
      for (const run of (AppState.activeExperiment?.runs || [])) {
        const option = document.createElement('option');
        option.value = run.run_id;
        option.textContent = run.run_label || run.run_id;
        option.selected = run.run_id === AppState.activeRun?.run_id;
        runs.appendChild(option);
      }
      const badge = document.getElementById('badge-exp-mode');
      const exp = AppState.activeExperiment;
      if (!badge) return;
      if (!exp) {
        badge.className = 'badge unreviewed';
        badge.innerText = 'NO EXPERIMENT';
        return;
      }
      if (exp.is_historical || exp.is_read_only) {
        badge.className = 'badge known';
        badge.innerText = 'READ-ONLY HISTORICAL';
        badge.title = 'Source database is mounted read-only (?mode=ro). Reviews saved to isolated overlay.';
      } else {
        badge.className = 'badge interesting';
        badge.innerText = 'ACTIVE EXPERIMENT';
        badge.title = 'Live research workspace database.';
      }
      const expLabel = document.getElementById('explore-active-run-label');
      if (expLabel) {
        expLabel.innerText = `Exp: ${exp.name}` + (AppState.activeRun ? ` (${AppState.activeRun.run_label || AppState.activeRun.run_id})` : '');
      }
    }

    async function onSwitchExperiment(expId, runId = null) {
      try {
        const resp = await fetch('/api/experiments/active', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ experiment_id: expId, run_id: runId })
        });
        const res = await resp.json();
        if (res.success) {
          showToast(`Switched active experiment to: ${res.experiment.name}`);
          AppState.activeExperiment = res.experiment;
          AppState.activeRun = res.run;
          updateExperimentBadge();
          updateProvenanceCard();
          AppState.currentPage = 1;
          AppState.expandedRows.clear();
          AppState.activeStar = null;
          AppState.activeEvent = null;
          AppState.activeLightCurve = null;
          AppState.lightCurveRequestId = (AppState.lightCurveRequestId || 0) + 1;
          AppState.stars = [];
          AppState.filteredStars = [];
          document.getElementById('stars-table-body').replaceChildren();
          document.getElementById('explore-count-label').textContent = 'Loading selected run…';
          document.getElementById('app-db-label').textContent = 'Loading selected run…';
          updateWorkspaceUI();
          switchView('explore');
          await loadStats();
          await fetchStars();
          checkRunStatus();
        } else {
          alert("Failed to switch experiment: " + (res.error || 'Unknown'));
        }
      } catch (e) {
        alert("Network error: " + e);
      }
    }

    function openNewExperimentModal() {
      initSectorCheckboxes();
      document.getElementById('modal-new-experiment').style.display = 'flex';
      document.getElementById('preflight-box').style.display = 'none';
      document.getElementById('target-val-status').innerText = '';
    }

    function closeNewExperimentModal() {
      document.getElementById('modal-new-experiment').style.display = 'none';
    }

    function selectPreset(preset) {
      AppState.selectedPreset = preset;
      const frozenCard = document.getElementById('preset-card-frozen');
      const customCard = document.getElementById('preset-card-custom');
      const customInputs = document.getElementById('custom-cfg-inputs');

      if (preset === 'frozen') {
        frozenCard.className = 'preset-card selected';
        customCard.className = 'preset-card';
        customInputs.style.display = 'none';
      } else {
        frozenCard.className = 'preset-card';
        customCard.className = 'preset-card selected';
        customInputs.style.display = 'grid';
      }
    }

    async function runTargetValidation() {
      const targetText = document.getElementById('new-exp-targets').value.trim();
      const checkedSectors = Array.from(document.querySelectorAll('input[name="sector-cb"]:checked')).map(cb => parseInt(cb.value));

      const statusLabel = document.getElementById('target-val-status');
      const preflightBox = document.getElementById('preflight-box');

      if (!targetText) {
        statusLabel.innerHTML = '<span style="color:#ef4444;">Please provide at least one target TIC.</span>';
        return;
      }
      if (checkedSectors.length === 0) {
        statusLabel.innerHTML = '<span style="color:#ef4444;">Please select at least one sector.</span>';
        return;
      }

      statusLabel.innerHTML = '<span style="color:var(--text-secondary);">Validating targets & estimating...</span>';

      try {
        const resp = await fetch('/api/experiments/validate-targets', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ targets_text: targetText, sectors: checkedSectors })
        });
        const data = await resp.json();
        if (data.success) {
          preflightBox.style.display = 'block';
          document.getElementById('pf-stars').innerText = data.total_unique_targets.toLocaleString();
          document.getElementById('pf-obs').innerText = data.total_observations.toLocaleString();
          document.getElementById('pf-mb').innerText = data.estimated_download_mb.toLocaleString();
          
          let warnHtml = '';
          if (data.invalid_count > 0) {
            warnHtml += `⚠️ Skipped ${data.invalid_count} invalid non-numeric tokens. `;
          }
          warnHtml += `Est. Execution Time: ~${data.estimated_runtime_seconds}s (${data.observations_per_target} sectors per star).`;
          document.getElementById('pf-warnings').innerHTML = warnHtml;
          statusLabel.innerHTML = '<span style="color:#34d399;">✓ Pre-flight validation passed!</span>';
        } else {
          statusLabel.innerHTML = `<span style="color:#ef4444;">Validation failed: ${data.error}</span>`;
          preflightBox.style.display = 'none';
        }
      } catch (e) {
        statusLabel.innerHTML = `<span style="color:#ef4444;">Network error: ${e}</span>`;
      }
    }

    async function submitCreateExperiment() {
      const name = document.getElementById('new-exp-name').value.trim();
      const objective = document.getElementById('new-exp-objective').value.trim();
      const desc = document.getElementById('new-exp-desc').value.trim();
      const targetText = document.getElementById('new-exp-targets').value.trim();
      const checkedSectors = Array.from(document.querySelectorAll('input[name="sector-cb"]:checked')).map(cb => parseInt(cb.value));
      const workers = parseInt(document.getElementById('new-exp-workers').value) || 2;
      const autolaunch = document.getElementById('new-exp-autolaunch').checked;

      if (!name) {
        alert("Please provide an Experiment Name.");
        return;
      }
      if (!objective) {
        alert("Please provide a Scientific Objective.");
        return;
      }
      if (!targetText) {
        alert("Please enter target TICs.");
        return;
      }
      if (checkedSectors.length === 0) {
        alert("Please select at least one observing sector.");
        return;
      }

      let threshold_t = 73.77965536097444;
      let duration_hours = 2.0;
      if (AppState.selectedPreset === 'custom') {
        threshold_t = parseFloat(document.getElementById('custom-threshold').value) || threshold_t;
        duration_hours = parseFloat(document.getElementById('custom-duration').value) || duration_hours;
      }

      const btn = document.getElementById('btn-create-experiment');
      btn.disabled = true;
      btn.innerText = "Creating Experiment...";

      try {
        const payload = {
          name: name,
          objective: objective,
          description: desc,
          targets_text: targetText,
          sectors: checkedSectors,
          workers: workers,
          threshold_t: threshold_t,
          config_preset: AppState.selectedPreset === 'custom' ? 'experimental_custom' : 'frozen_fixed_box_2h',
          duration_hr: duration_hours,
          launch_immediately: autolaunch
        };

        const resp = await fetch('/api/experiments/create', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(payload)
        });
        const res = await resp.json();

        if (res.success) {
          showToast(`Created experiment: ${name}`);
          closeNewExperimentModal();
          AppState.activeStar = null;
          AppState.activeEvent = null;
          AppState.activeLightCurve = null;
          AppState.lightCurveRequestId = (AppState.lightCurveRequestId || 0) + 1;
          await loadExperiments();
          await loadStats();
          await fetchStars();
          if (autolaunch && res.run) {
            switchView('progress');
            checkRunStatus();
          }
        } else {
          alert("Failed to create experiment: " + (res.error || 'Unknown error'));
        }
      } catch (e) {
        alert("Network error: " + e);
      } finally {
        btn.disabled = false;
        btn.innerText = "Create & Launch Experiment";
      }
    }

    function checkRunStatus() {
      if (AppState.pollingTimer) {
        clearInterval(AppState.pollingTimer);
        AppState.pollingTimer = null;
      }
      if (!AppState.activeRun) {
        updateProgressUI(null);
        return;
      }
      pollProgressOnce();
      if (AppState.activeRun.status === 'running' || AppState.activeRun.status === 'pending') {
        AppState.pollingTimer = setInterval(pollProgressOnce, 2000);
      }
    }

    async function pollProgressOnce() {
      if (!AppState.activeRun) return;
      const runId = AppState.activeRun.run_id;
      try {
        const resp = await fetch(`/api/runs/${runId}/progress`);
        const data = await resp.json();
        if (runId !== AppState.activeRun?.run_id) return;
        if (data.run) {
          AppState.activeRun = data.run;
          updateProgressUI(data);
          if (data.run.status !== 'running' && data.run.status !== 'pending') {
            if (AppState.pollingTimer) {
              clearInterval(AppState.pollingTimer);
              AppState.pollingTimer = null;
            }
            loadStats();
            fetchStars();
          }
        }
      } catch (e) {
        console.error("Progress poll failed", e);
      }
    }

    function updateProgressUI(data) {
      const banner = document.getElementById('explore-running-banner');
      const bannerText = document.getElementById('banner-run-progress');
      const liveContainer = document.getElementById('live-progress-container');
      const btnAbort = document.getElementById('btn-abort-survey');
      const btnResume = document.getElementById('btn-resume-survey');
      const statusBadge = document.getElementById('progress-status-badge');
      const headerTitle = document.getElementById('progress-header-title');

      if (!data || !data.run) {
        if (banner) banner.style.display = 'none';
        if (liveContainer) liveContainer.style.display = 'none';
        if (btnAbort) btnAbort.style.display = 'none';
        if (btnResume) btnResume.style.display = 'none';
        if (statusBadge) {
          const isHist = AppState.activeExperiment && AppState.activeExperiment.is_historical;
          statusBadge.className = 'badge ' + (isHist ? 'known' : 'unreviewed');
          statusBadge.innerText = isHist ? 'HISTORICAL BENCHMARK' : 'NO ACTIVE RUN';
        }
        return;
      }

      const run = data.run;
      const isRunning = run.status === 'running' || run.status === 'pending';
      const isCompleted = run.status === 'completed';
      const isAborted = run.status === 'aborted';
      const isFailed = ['failed', 'interrupted', 'completed_with_errors'].includes(run.status);

      if (headerTitle) {
        headerTitle.innerText = `SURVEY RUN #${run.run_id} (${run.status.toUpperCase()})`;
      }

      if (statusBadge) {
        let badgeClass = 'unreviewed';
        if (isRunning) badgeClass = 'interesting';
        else if (isCompleted) badgeClass = 'known';
        else if (isAborted) badgeClass = 'dismissed';
        else if (isFailed) badgeClass = 'artifact';
        statusBadge.className = 'badge ' + badgeClass;
        statusBadge.innerText = run.status.toUpperCase();
      }

      if (banner) {
        banner.style.display = isRunning ? 'flex' : 'none';
        if (bannerText) {
          bannerText.innerText = `${run.processed_targets}/${run.total_targets} obs (${run.candidate_count} candidates)`;
        }
      }

      if (liveContainer) {
        liveContainer.style.display = isRunning ? 'flex' : 'none';
        const pct = run.total_targets > 0 ? Math.min(100, (run.processed_targets / run.total_targets * 100)).toFixed(1) : '0.0';
        document.getElementById('live-proc-count').innerText = run.processed_targets;
        document.getElementById('live-total-count').innerText = run.total_targets;
        document.getElementById('live-percent-label').innerText = `${pct}%`;
        document.getElementById('live-progress-bar').style.width = `${pct}%`;
      }

      if (btnAbort) {
        btnAbort.style.display = isRunning ? 'inline-block' : 'none';
      }
      if (btnResume) {
        btnResume.style.display = (isAborted || isFailed) ? 'inline-block' : 'none';
      }

      document.getElementById('tel-obs-count').innerText = `${run.processed_targets} / ${run.total_targets}`;
      document.getElementById('tel-candidates-count').innerText = run.candidate_count;
      if (data.unique_stars !== undefined) {
        document.getElementById('tel-stars-count').innerText = data.unique_stars;
      }
      if (run.throughput_obs_hr) {
        document.getElementById('tel-throughput').innerText = Math.round(run.throughput_obs_hr).toLocaleString();
      } else {
        document.getElementById('tel-throughput').innerText = isRunning ? 'Calculating...' : '—';
      }
      if (run.eta_sec !== null && run.eta_sec !== undefined && isRunning) {
        const m = Math.floor(run.eta_sec / 60);
        const s = Math.round(run.eta_sec % 60);
        document.getElementById('tel-eta').innerText = `${m}m ${s}s`;
      } else if (run.elapsed_sec) {
        const m = Math.floor(run.elapsed_sec / 60);
        const s = Math.round(run.elapsed_sec % 60);
        document.getElementById('tel-eta').innerText = `${m}m ${s}s`;
        document.getElementById('tel-eta-sub').innerText = 'Total elapsed time';
      } else {
        document.getElementById('tel-eta').innerText = '—';
      }

      const logBox = document.getElementById('live-logs-box');
      if (logBox && data.logs && data.logs.length > 0) {
        logBox.innerHTML = '';
        data.logs.forEach(line => {
          const d = document.createElement('div');
          d.innerText = `${line.created_utc} [${line.level}] ${line.message}`;
          logBox.appendChild(d);
        });
        logBox.scrollTop = logBox.scrollHeight;
      }
    }

    function updateProvenanceCard() {
      const exp = AppState.activeExperiment;
      const run = AppState.activeRun;
      const hashSpan = document.getElementById('provenance-hash');
      const objP = document.getElementById('provenance-objective');
      const cfgP = document.getElementById('provenance-config');

      if (!exp) return;
      const provenance = run?.scientific_provenance || (exp.is_read_only ? {} : exp);

      if (hashSpan) {
        hashSpan.innerText = provenance.config_hash ? `Hash: ${provenance.config_hash.substring(0, 16)}...` : "Hash: unavailable";
      }
      if (objP) {
        objP.textContent = `Objective: ${exp.scientific_objective || exp.description || 'General anomaly exploration.'}`;
      }
      if (cfgP) {
        let cfg = `Method: ${provenance.method || (exp.is_read_only ? 'unavailable' : 'fixed_box_2h')} | Threshold T: ${provenance.threshold_t ?? 'unavailable'} | Supported cadence: 120s`;
        if (run) {
          cfg += ` | Output DB: ${run.output_db_path || 'sqlite'}`;
        }
        cfgP.textContent = cfg;
      }
    }

    async function abortSurveyRun() {
      if (!AppState.activeRun) return;
      if (!confirm(`Are you sure you want to abort Run #${AppState.activeRun.run_id}? Ingested observations will be preserved.`)) {
        return;
      }
      try {
        const resp = await fetch('/api/runs/abort', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ run_id: AppState.activeRun.run_id })
        });
        const res = await resp.json();
        if (res.success) {
          showToast(`Run #${AppState.activeRun.run_id} abort signal sent.`);
          pollProgressOnce();
        } else {
          alert("Failed to abort run: " + (res.error || 'Unknown'));
        }
      } catch (e) {
        alert("Network error: " + e);
      }
    }

    async function resumeSurveyRun() {
      if (!AppState.activeRun) return;
      try {
        const resp = await fetch('/api/experiments/launch', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            experiment_id: AppState.activeExperiment.id,
            run_id: AppState.activeRun.run_id,
            resume: true
          })
        });
        const res = await resp.json();
        if (res.success) {
          showToast(`Survey run #${AppState.activeRun.run_id} resumed!`);
          AppState.activeRun = res.run;
          checkRunStatus();
        } else {
          alert("Failed to resume run: " + (res.error || 'Unknown'));
        }
      } catch (e) {
        alert("Network error: " + e);
      }
    }

    async function exportCatalog(format) {
      if (!AppState.activeRun) return;
      try {
        const response = await fetch(`/api/runs/${AppState.activeRun.run_id}/export?format=${format}`);
        if (!response.ok) throw new Error(await response.text());
        const blob = await response.blob();
        const href = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = href;
        link.download = `tess_catalog_${AppState.activeRun.run_id}.${format}`;
        document.body.appendChild(link);
        link.click();
        link.remove();
        setTimeout(() => URL.revokeObjectURL(href), 30000);
        showToast(`Exported ${format.toUpperCase()} for ${AppState.activeRun.run_label || AppState.activeRun.run_id}`);
      } catch (error) {
        alert('Export failed: ' + error.message);
      }
    }

    // Load initial data
    window.addEventListener('DOMContentLoaded', async () => {
      initSectorCheckboxes();
      await loadExperiments();
      await loadStats();
      await fetchStars();

      // Hash routing
      const hash = window.location.hash;
      if (hash.includes('workspace')) {
        switchView('workspace');
        const params = new URLSearchParams(hash.substring(1));
        const tic = params.get('tic');
        if (tic) loadStarDetails(tic);
      } else if (hash.includes('progress')) {
        switchView('progress');
      } else if (hash.includes('modal=new')) {
        openNewExperimentModal();
      }
    });

    function showToast(msg) {
      const container = document.getElementById('toast-container');
      const toast = document.createElement('div');
      toast.className = 'toast';
      toast.textContent = '✓ ' + msg;
      container.appendChild(toast);
      setTimeout(() => {
        toast.style.opacity = '0';
        toast.style.transition = 'opacity 0.3s';
        setTimeout(() => toast.remove(), 300);
      }, 2500);
    }

    // Load Survey Stats
    async function loadStats() {
      const runId = AppState.activeRun?.run_id;
      try {
        const resp = await fetch('/api/stats');
        const data = await resp.json();
        if (runId !== AppState.activeRun?.run_id) return;
        document.getElementById('app-db-label').innerText = (data.data_kind === 'synthetic demonstration' ? 'SYNTHETIC DEMO: ' : 'DB: ') + (data.db_name || 'Connected');
        document.getElementById('tel-obs-count').innerText = AppState.activeRun ? `${AppState.activeRun.completed_obs} / ${AppState.activeRun.total_targets}` : 'Unavailable';
        document.getElementById('tel-stars-count').innerText = (data.stars ? data.stars.total : '—').toLocaleString();
        const dimmings = document.getElementById('tel-dimmings-count');
        if (dimmings) dimmings.innerText = (data.verifications?.verified_localized_dimming || 0).toLocaleString();
      } catch (e) {
        console.error("Failed to load stats", e);
      }
    }

    // Fetch Stars for Explore Page
    async function fetchStars() {
      const runId = AppState.activeRun?.run_id;
      const requestId = (AppState.catalogRequestId || 0) + 1;
      AppState.catalogRequestId = requestId;
      const status = document.getElementById('filter-status').value;
      const verif = document.getElementById('filter-verification').value;
      const sysType = document.getElementById('filter-system').value;
      const sort = document.getElementById('filter-sort').value;
      const search = document.getElementById('filter-search').value;

      const url = `/api/stars?page=${AppState.currentPage}&limit=${AppState.pageSize}&status=${status}&verification=${verif}&system_type=${sysType}&sort_by=${sort}&search=${encodeURIComponent(search)}`;
      try {
        const resp = await fetch(url);
        const data = await resp.json();
        if (runId !== AppState.activeRun?.run_id || requestId !== AppState.catalogRequestId) return;
        AppState.stars = data.stars || [];
        AppState.filteredStars = AppState.stars;
        renderStarsTable(data);
      } catch (e) {
        console.error("Failed to fetch stars", e);
      }
    }

    function renderStarsTable(data) {
      const tbody = document.getElementById('stars-table-body');
      tbody.innerHTML = '';

      document.getElementById('explore-count-label').innerText = `Showing ${data.stars.length} of ${data.total_stars} stars`;
      document.getElementById('explore-pagination-info').innerText = `Page ${data.page} of ${data.total_pages}`;
      document.getElementById('btn-prev-page').disabled = (data.page <= 1);
      document.getElementById('btn-next-page').disabled = (data.page >= data.total_pages);

      if (data.stars.length === 0) {
        tbody.innerHTML = '<tr><td colspan="9" style="text-align:center; padding:30px; color:var(--text-secondary);">No stars matched the selected filters.</td></tr>';
        return;
      }

      data.stars.forEach((star, idx) => {
        const isExpanded = AppState.expandedRows.has(star.tic);
        const row = document.createElement('tr');
        row.className = 'star-row';
        row.onclick = () => toggleRow(star.tic);

        const verifBadge = `<span class="badge ${star.primary_verification === 'verified_localized_dimming' ? 'interesting' : 'unreviewed'}">${star.primary_verification}</span>`;
        const statusBadge = `<span class="badge ${star.star_status}">${star.star_status}</span>`;

        let periodStr = '<span style="color:var(--text-muted);">Isolated</span>';
        if (star.is_periodic && star.period_days) {
          periodStr = `<b style="color:var(--accent);">P = ${star.period_days.toFixed(4)} d</b><br><span style="font-size:10px; color:var(--text-muted);">${star.system_classification || 'Periodic System'}</span>`;
        }

        row.innerHTML = `
          <td style="text-align:center; color:var(--accent); font-weight:bold;">${isExpanded ? '▼' : '▶'}</td>
          <td style="font-family:monospace; font-weight:bold; color:#fff;">${star.tic}</td>
          <td>${star.sector_display}</td>
          <td style="text-align:center; font-weight:bold;">${star.total_events}</td>
          <td style="font-family:monospace; color:#fbbf24;">${star.max_score.toFixed(2)}</td>
          <td>${verifBadge}</td>
          <td>${periodStr}</td>
          <td>${statusBadge}</td>
          <td style="text-align:right;">
            <button class="btn btn-sm" onclick="event.stopPropagation(); openStarInvestigation('${star.tic}', ${idx})">Investigate →</button>
          </td>
        `;
        tbody.appendChild(row);

        // Expanded Events Sub-Table
        if (isExpanded) {
          const subRow = document.createElement('tr');
          subRow.innerHTML = `
            <td colspan="9" class="events-subpanel">
              <div style="font-weight:700; font-size:11px; color:var(--text-secondary); margin-bottom:6px;">CONSTITUENT DETECTED EVENTS (${star.events.length}):</div>
              <table style="width:100%; border-collapse:collapse; font-size:11px;">
                <thead>
                  <tr style="color:var(--text-muted); border-bottom:1px solid #1e293b; text-align:left;">
                    <th style="padding:4px 6px;">Event ID</th>
                    <th style="padding:4px 6px;">Center (BTJD)</th>
                    <th style="padding:4px 6px;">Score</th>
                    <th style="padding:4px 6px;">Duration</th>
                    <th style="padding:4px 6px;">Depth (PPM)</th>
                    <th style="padding:4px 6px;">Verification</th>
                    <th style="padding:4px 6px;">Review</th>
                    <th style="padding:4px 6px; text-align:right;">Action</th>
                  </tr>
                </thead>
                <tbody>
                  ${star.events.map(ev => `
                    <tr style="border-bottom:1px solid #1e293b; cursor:pointer;" onclick="event.stopPropagation(); openStarInvestigation('${star.tic}', ${idx}, ${ev.id})">
                      <td style="padding:5px 6px; font-family:monospace; color:var(--accent);">#${ev.id}</td>
                      <td style="padding:5px 6px; font-family:monospace;">${ev.center_btjd}</td>
                      <td style="padding:5px 6px; font-family:monospace; color:#fbbf24;">${ev.score}</td>
                      <td style="padding:5px 6px;">${ev.width_hr} h</td>
                      <td style="padding:5px 6px; font-family:monospace;">${ev.depth_ppm}</td>
                      <td style="padding:5px 6px;"><span class="badge ${ev.scientific_classification === 'verified_localized_dimming' ? 'interesting' : 'unreviewed'}">${ev.scientific_classification}</span></td>
                      <td style="padding:5px 6px;"><span class="badge ${ev.review_status}">${ev.review_status}</span></td>
                      <td style="padding:5px 6px; text-align:right;"><button class="btn btn-secondary btn-sm" style="font-size:10px; padding:2px 6px;">Inspect Event</button></td>
                    </tr>
                  `).join('')}
                </tbody>
              </table>
            </td>
          `;
          tbody.appendChild(subRow);
        }
      });
    }

    function toggleRow(tic) {
      if (AppState.expandedRows.has(tic)) {
        AppState.expandedRows.delete(tic);
      } else {
        AppState.expandedRows.add(tic);
      }
      fetchStars();
    }

    function applyFilters() {
      AppState.currentPage = 1;
      fetchStars();
    }

    function resetFilters() {
      document.getElementById('filter-search').value = '';
      document.getElementById('filter-status').value = 'all';
      document.getElementById('filter-verification').value = 'all';
      document.getElementById('filter-system').value = 'all';
      document.getElementById('filter-sort').value = 'score';
      applyFilters();
    }

    function prevPage() {
      if (AppState.currentPage > 1) {
        AppState.currentPage--;
        fetchStars();
      }
    }

    function nextPage() {
      AppState.currentPage++;
      fetchStars();
    }

    // Open Investigation Workspace for Star
    async function openStarInvestigation(tic, starIndex, eventId = null) {
      AppState.currentStarIndex = (starIndex !== undefined) ? starIndex : 0;
      await loadStarDetails(tic, eventId);
      switchView('workspace');
    }

    async function loadStarDetails(tic, preferredEventId = null) {
      const requestId = (AppState.starRequestId || 0) + 1;
      AppState.starRequestId = requestId;
      const runId = AppState.activeRun?.run_id;
      try {
        const resp = await fetch(`/api/star/${encodeURIComponent(tic)}`);
        const data = await resp.json();
        if (requestId !== AppState.starRequestId || runId !== AppState.activeRun?.run_id) return;
        if (!resp.ok) throw new Error(data.error || 'Star unavailable');
        AppState.activeStar = data;
        if (AppState.currentView === 'workspace') history.replaceState(null, '', `#view=workspace&tic=${encodeURIComponent(tic)}`);

        // Pick active event
        if (preferredEventId) {
          AppState.activeEvent = data.events.find(e => e.id === preferredEventId) || data.events[0];
        } else {
          AppState.activeEvent = data.events[0] || null;
        }

        updateWorkspaceUI();
        await fetchLightCurveData(tic, AppState.activeEvent ? AppState.activeEvent.sector : null);
      } catch (e) {
        console.error("Failed to load star details", e);
      }
    }

    function updateWorkspaceUI() {
      const star = AppState.activeStar;
      if (!star) {
        document.getElementById('ws-tic-label').innerText = 'TIC —';
        document.getElementById('ws-star-class-badge').className = 'badge unreviewed';
        document.getElementById('ws-star-class-badge').innerText = 'UNREVIEWED';
        document.getElementById('ws-benchmark-badge').style.display = 'none';
        document.getElementById('ws-sector-val').innerText = '—';
        document.getElementById('ws-tmag-val').innerText = '—';
        document.getElementById('ws-teff-val').innerText = '—';
        document.getElementById('ws-events-count-val').innerText = '—';
        document.getElementById('ws-star-index-label').innerText = 'No star selected';
        document.getElementById('ws-btn-prev-star').disabled = true;
        document.getElementById('ws-btn-next-star').disabled = true;
        document.getElementById('ws-event-tabs').innerHTML = '';

        document.getElementById('star-rev-status').value = 'unreviewed';
        document.getElementById('star-rev-class').value = '';
        document.getElementById('star-rev-note').value = '';
        document.getElementById('star-rev-status-badge').className = 'badge unreviewed';
        document.getElementById('star-rev-status-badge').innerText = 'unreviewed';

        document.getElementById('event-rev-id-label').innerText = '#—';
        document.getElementById('event-rev-status').value = 'unreviewed';
        document.getElementById('event-rev-note').value = '';
        document.getElementById('event-rev-status-badge').className = 'badge unreviewed';
        document.getElementById('event-rev-status-badge').innerText = 'unreviewed';

        document.getElementById('prov-score').innerText = '—';
        document.getElementById('prov-kernel-info').innerText = 'Select a target star to evaluate detector significance.';
        document.getElementById('prov-verif-class').innerText = '—';
        document.getElementById('prov-depth').innerText = '—';
        document.getElementById('prov-flare').innerText = '—';
        document.getElementById('prov-verif-desc').innerText = 'Awaiting target selection.';
        document.getElementById('prov-coincidence-level').innerText = '—';
        document.getElementById('prov-coincidence-desc').innerText = 'No target star or event selected.';
        document.getElementById('prov-provisional-class').innerText = '—';
        document.getElementById('prov-provisional-desc').innerText = 'Select a star and candidate event to inspect scientific evidence.';
        document.getElementById('ws-audit-history-list').innerHTML = '<span>No target star selected.</span>';
        document.getElementById('card-phase-fold').style.display = 'none';
        for (const id of ['plot-full-container', 'plot-zoom-container', 'plot-phase-container']) {
          const container = document.getElementById(id);
          if (typeof Plotly !== 'undefined') Plotly.purge(container);
          container.textContent = 'Select a target star.';
        }
        return;
      }

      document.getElementById('ws-tic-label').innerText = star.tic;
      document.getElementById('ws-star-class-badge').className = 'badge ' + star.star_status;
      document.getElementById('ws-star-class-badge').innerText = (star.star_class ? star.star_class.toUpperCase() : star.star_status.toUpperCase());
      document.getElementById('ws-sector-val').innerText = (star.events[0] ? star.events[0].sector : '—');
      document.getElementById('ws-tmag-val').innerText = (star.observation && star.observation.tmag ? star.observation.tmag.toFixed(2) : '—');
      document.getElementById('ws-teff-val').innerText = (star.observation && star.observation.teff ? star.observation.teff.toFixed(0) + ' K' : '—');
      document.getElementById('ws-events-count-val').innerText = star.total_events + ' detected';

      const demoBadge = document.getElementById('ws-benchmark-badge');
      demoBadge.textContent = 'SYNTHETIC DEMONSTRATION';
      demoBadge.style.display = star.data_kind === 'synthetic demonstration' ? 'inline-block' : 'none';
      // Navigation counter
      const totalInFilter = AppState.filteredStars.length || 1;
      document.getElementById('ws-star-index-label').innerText = `Star ${AppState.currentStarIndex + 1} of ${totalInFilter}`;
      document.getElementById('ws-btn-prev-star').disabled = (AppState.currentStarIndex <= 0);
      document.getElementById('ws-btn-next-star').disabled = (AppState.currentStarIndex >= totalInFilter - 1);

      // Event selector tabs
      const tabsDiv = document.getElementById('ws-event-tabs');
      tabsDiv.innerHTML = '';
      star.events.forEach(ev => {
        const btn = document.createElement('button');
        const isActive = AppState.activeEvent && AppState.activeEvent.id === ev.id;
        btn.className = 'btn btn-sm ' + (isActive ? '' : 'btn-secondary');
        btn.style.fontSize = '10px';
        btn.style.padding = '2px 6px';
        btn.innerText = `#${ev.id} (${ev.center_btjd})`;
        btn.onclick = () => selectActiveEvent(ev);
        tabsDiv.appendChild(btn);
      });

      // Update Star Review Form
      document.getElementById('star-rev-status').value = star.star_status;
      document.getElementById('star-rev-class').value = star.star_class || '';
      document.getElementById('star-rev-note').value = star.note || '';
      document.getElementById('star-rev-status-badge').className = 'badge ' + star.star_status;
      document.getElementById('star-rev-status-badge').innerText = star.star_status;

      // Update Event Review Form & Scientific Provenance
      if (AppState.activeEvent) {
        const ev = AppState.activeEvent;
        document.getElementById('event-rev-id-label').innerText = '#' + ev.id;
        document.getElementById('event-rev-status').value = ev.review_status;
        document.getElementById('event-rev-note').value = ev.note || '';
        document.getElementById('event-rev-status-badge').className = 'badge ' + ev.review_status;
        document.getElementById('event-rev-status-badge').innerText = ev.review_status;

        // Tier 1: Matched-Filter Detection Metric
        document.getElementById('prov-score').innerText = `S = ${ev.score.toFixed(2)}`;
        document.getElementById('prov-kernel-info').innerText = `Matched ${ev.width_hr || 2.0}h box kernel (B1 Detector, T=${star.execution?.threshold_t ?? "unavailable"})`;

        // Tier 2: Automated Photometric Verification (5-State Machine)
        const verifClass = ev.scientific_classification || 'not_run';
        document.getElementById('prov-verif-class').innerText = verifClass.replace(/_/g, ' ').toUpperCase();
        document.getElementById('prov-verif-class').style.color = (verifClass === 'verified_localized_dimming') ? '#34d399' : (verifClass.includes('flare') ? '#fbbf24' : '#60a5fa');
        document.getElementById('prov-depth').innerText = `Matched-filter: ${ev.depth_ppm} ppm; photometric minimum: ${ev.true_depth_min_pct === null ? "unavailable" : ev.true_depth_min_pct + "%"}`;

        if (ev.has_shoulder_flare) {
          const elev = ev.shoulder_elevation_pct ? ` (+${ev.shoulder_elevation_pct.toFixed(2)}%)` : '';
          document.getElementById('prov-flare').innerText = `Elevated Shoulder${elev}`;
          document.getElementById('prov-flare').style.color = '#fbbf24';
        } else {
          document.getElementById('prov-flare').innerText = ev.has_shoulder_flare === false ? 'None Detected' : 'Not run / unavailable';
          document.getElementById('prov-flare').style.color = '#34d399';
        }
        document.getElementById('prov-verif-desc').innerText = ev.scientific_description || 
          (verifClass === 'verified_localized_dimming' 
            ? 'Confirmed localized flux depression with robust continuum recovery.' 
            : 'Automated 5-state photometric verification classification.');

        // Tier 3: Spatial Coincidence Risk
        const coincRisk = ev.coincidence_risk_level || star.coincidence_risk || 'unknown';
        const coincElem = document.getElementById('prov-coincidence-level');
        const coincDesc = document.getElementById('prov-coincidence-desc');
        if (coincRisk === 'high' || (ev.coincidence_star_count && ev.coincidence_star_count > 1)) {
          coincElem.innerText = 'Level: HIGH (Simultaneous Triggers)';
          coincElem.style.color = '#ef4444';
          coincDesc.innerText = 'Multi-star simultaneous trigger detected across field. Potential instrumental systematic or scattered light.';
        } else {
          coincElem.innerText = coincRisk === 'none' ? 'Level: NONE' : 'Level: UNKNOWN / NOT RUN';
          coincElem.style.color = '#34d399';
          coincDesc.innerText = coincRisk === 'none' ? 'Recorded coincidence analysis found no match.' : 'No coincidence analysis result is recorded.';
        }

        // Tier 4: Provisional Morphological Interpretation
        const provClassElem = document.getElementById('prov-provisional-class');
        const provDescElem = document.getElementById('prov-provisional-desc');
        if (star.system && star.system.is_periodic && star.system.period_days) {
          provClassElem.innerText = `${star.system.system_classification || star.system.classification || 'Periodic Eclipsing System'}`;
          provClassElem.style.color = 'var(--accent)';
          provDescElem.innerText = `Multi-event periodic system (P = ${star.system.period_days.toFixed(4)} d, T0 = ${star.system.epoch_btjd.toFixed(4)} BTJD).`;
        } else if (verifClass === 'flare_induced_apparent_dip') {
          provClassElem.innerText = 'Flare-Induced Apparent Dip';
          provClassElem.style.color = '#fbbf24';
          provDescElem.innerText = 'Apparent dip artifact resulting from elevated flare shoulder.';
        } else if (verifClass === 'continuous_variability') {
          provClassElem.innerText = 'Continuous Stellar Variability';
          provClassElem.style.color = '#60a5fa';
          provDescElem.innerText = 'Continuous periodic or quasi-periodic stellar pulsation / modulation.';
        } else {
          provClassElem.innerText = star.total_events > 1 ? 'Multiple events; period unavailable' : 'Isolated Single-Dip Candidate';
          provClassElem.style.color = '#34d399';
          provDescElem.innerText = 'Empirical localized negative excursion. Requires multi-sector or radial-velocity confirmation.';
        }
      }

      // Audit Trail
      const histDiv = document.getElementById('ws-audit-history-list');
      if (star.history && star.history.length > 0) {
        histDiv.innerHTML = star.history.map(h => `
          <div style="display:flex; justify-content:space-between; border-bottom:1px solid #1e293b; padding:2px 0;">
            <span>${h.created_utc.substring(0, 19).replace('T', ' ')} by <b>${h.reviewer}</b></span>
            <span class="badge ${h.new_status}">${h.previous_status} → ${h.new_status}</span>
          </div>
        `).join('');
      } else {
        histDiv.innerHTML = '<span>No prior modifications recorded.</span>';
      }

      // Phase Folding Card
      if (star.system && star.system.is_periodic && star.system.period_days) {
        document.getElementById('card-phase-fold').style.display = 'flex';
        document.getElementById('ws-phase-info').innerText = `P = ${star.system.period_days.toFixed(4)} d | T0 = ${star.system.epoch_btjd.toFixed(4)} BTJD`;
      } else {
        document.getElementById('card-phase-fold').style.display = 'none';
      }
    }

    function selectActiveEvent(ev) {
      AppState.activeEvent = ev;
      updateWorkspaceUI();
      fetchLightCurveData(AppState.activeStar.tic, ev.sector);
    }

    function navigateStar(delta) {
      const newIndex = AppState.currentStarIndex + delta;
      if (newIndex >= 0 && newIndex < AppState.filteredStars.length) {
        AppState.currentStarIndex = newIndex;
        const targetStar = AppState.filteredStars[newIndex];
        loadStarDetails(targetStar.tic);
      }
    }

    // Light Curve Fetching & Plotly Rendering
    async function fetchLightCurveData(tic, sector) {
      const fullCont = document.getElementById('plot-full-container');
      const zoomCont = document.getElementById('plot-zoom-container');
      const phaseCont = document.getElementById('plot-phase-container');

      const loadingHtml = `
        <div style="display:flex; flex-direction:column; align-items:center; justify-content:center; height:100%; gap:8px; color:var(--text-secondary);">
          <span class="pulse-indicator" style="width:12px; height:12px; margin:0;"></span>
          <span style="font-size:12px;">Retrieving light curve for TIC ${tic} (checking cache & NASA MAST archive)...</span>
        </div>`;
      const requestId = (AppState.lightCurveRequestId || 0) + 1;
      AppState.lightCurveRequestId = requestId;
      AppState.activeLightCurve = null;
      for (const container of [fullCont, zoomCont, phaseCont]) {
        if (container && typeof Plotly !== 'undefined') Plotly.purge(container);
      }
      fullCont.innerHTML = loadingHtml;
      zoomCont.innerHTML = loadingHtml;

      try {
        const url = `/api/star/${encodeURIComponent(tic)}/lightcurve?candidate_id=${AppState.activeEvent ? AppState.activeEvent.id : ''}`;
        const resp = await fetch(url);
        const data = await resp.json();

        if (!resp.ok || data.error) {
          throw new Error(data.error || `HTTP ${resp.status}`);
        }

        if (requestId !== AppState.lightCurveRequestId) return;
        fullCont.innerHTML = '';
        zoomCont.innerHTML = '';
        AppState.activeLightCurve = data;
        document.getElementById('ws-full-lc-source').innerText = 'Source: ' + (data.source || 'NASA MAST Archive');
        renderWorkspacePlots();
      } catch (e) {
        if (requestId !== AppState.lightCurveRequestId) return;
        console.error("Failed to load light curve series", e);
        const errorHtml = `
          <div style="display:flex; flex-direction:column; align-items:center; justify-content:center; height:100%; gap:8px; padding:16px; text-align:center;">
            <div style="color:#ef4444; font-size:12px; font-weight:700;">Light Curve Unavailable</div>
            <div style="color:var(--text-secondary); font-size:11px; max-width:420px; line-height:1.4;">${e.message || 'Light curve file missing locally and could not be retrieved from NASA MAST.'}</div>
            <button class="btn btn-sm btn-secondary" style="margin-top:4px;" onclick="fetchLightCurveData('${tic}', ${sector || 'null'})">↻ Retry Retrieval</button>
          </div>`;
        fullCont.innerHTML = errorHtml;
        zoomCont.innerHTML = errorHtml;
        if (phaseCont) phaseCont.innerHTML = '';
      }
    }

    function renderWorkspacePlots() {
      if (AppState.useStaticPlots || typeof Plotly === 'undefined') {
        renderStaticPlots();
        return;
      }
      renderPlotlyFull();
      renderPlotlyZoom();
      if (AppState.activeStar && AppState.activeStar.system && AppState.activeStar.system.is_periodic) {
        renderPlotlyPhase();
      }
    }

    function togglePlotMode() {
      AppState.useStaticPlots = !AppState.useStaticPlots;
      renderWorkspacePlots();
    }

    function renderStaticPlots() {
      const cid = AppState.activeEvent ? AppState.activeEvent.id : 1;
      document.getElementById('plot-full-container').innerHTML = `<img src="/api/plot/full/${cid}?workspace_run_id=${AppState.activeRun.run_id}" style="width:100%; height:100%; object-fit:contain;">`;
      document.getElementById('plot-zoom-container').innerHTML = `<img src="/api/plot/zoom/${cid}?hours=12&workspace_run_id=${AppState.activeRun.run_id}" style="width:100%; height:100%; object-fit:contain;">`;
      if (AppState.activeStar && AppState.activeStar.system && AppState.activeStar.system.is_periodic) {
        document.getElementById('plot-phase-container').innerHTML = `<img src="/api/plot/phase/${cid}?workspace_run_id=${AppState.activeRun.run_id}" style="width:100%; height:100%; object-fit:contain;">`;
      }
    }

    function renderPlotlyFull() {
      const lc = AppState.activeLightCurve;
      if (!lc || !lc.time) return;

      const traces = [{
        x: lc.time,
        y: lc.flux_ppm,
        mode: 'lines+markers',
        type: 'scattergl',
        line: { color: '#475569', width: 0.8 },
        marker: { size: 3, color: '#0284c7', opacity: 0.7 },
        name: 'PDC Flux'
      }];

      // Event vertical markers
      const shapes = [];
      const annotations = [];
      lc.events.forEach(ev => {
        const isActive = AppState.activeEvent && AppState.activeEvent.id === ev.id;
        shapes.push({
          type: 'line',
          x0: ev.center_btjd,
          x1: ev.center_btjd,
          y0: 0,
          y1: 1,
          yref: 'paper',
          line: {
            color: isActive ? '#ef4444' : '#f59e0b',
            width: isActive ? 2 : 1.2,
            dash: isActive ? 'dash' : 'dot'
          }
        });
        annotations.push({
          x: ev.center_btjd,
          y: 1,
          yref: 'paper',
          text: `#${ev.id}`,
          showarrow: false,
          font: { size: 9, color: isActive ? '#ef4444' : '#f59e0b' }
        });
      });

      const layout = {
        paper_bgcolor: '#090d16',
        plot_bgcolor: '#090d16',
        margin: { l: 50, r: 20, t: 25, b: 35 },
        xaxis: { title: { text: 'TESS Time (BTJD)', font: { size: 11, color: '#94a3b8' } }, tickfont: { size: 10, color: '#64748b' }, gridcolor: '#1e293b' },
        yaxis: { title: { text: 'Relative Flux (ppm)', font: { size: 11, color: '#94a3b8' } }, tickfont: { size: 10, color: '#64748b' }, gridcolor: '#1e293b' },
        shapes: shapes,
        annotations: annotations,
        showlegend: false
      };

      Plotly.react('plot-full-container', traces, layout, { responsive: true, displayModeBar: false });
    }

    function renderPlotlyZoom() {
      const lc = AppState.activeLightCurve;
      const ev = AppState.activeEvent;
      if (!lc || !lc.time || !ev) return;

      const spanHours = parseFloat(document.getElementById('ws-zoom-span').value) || 12;
      const spanDays = spanHours / 24.0;
      const tMin = ev.center_btjd - spanDays;
      const tMax = ev.center_btjd + spanDays;

      // Filter points within zoom window
      const xZoom = [];
      const yZoom = [];
      for (let i = 0; i < lc.time.length; i++) {
        if (lc.time[i] >= tMin && lc.time[i] <= tMax) {
          xZoom.push(lc.time[i]);
          yZoom.push(lc.flux_ppm[i]);
        }
      }

      const halfW = (ev.width_hr / 24.0) / 2.0;
      const traces = [{
        x: xZoom,
        y: yZoom,
        mode: 'markers',
        type: 'scatter',
        marker: { size: 5, color: '#38bdf8', opacity: 0.8 },
        name: 'Flux Points'
      }];

      const layout = {
        paper_bgcolor: '#090d16',
        plot_bgcolor: '#090d16',
        margin: { l: 50, r: 20, t: 20, b: 35 },
        xaxis: { range: [tMin, tMax], title: { text: 'TESS Time (BTJD)', font: { size: 11, color: '#94a3b8' } }, tickfont: { size: 10, color: '#64748b' }, gridcolor: '#1e293b' },
        yaxis: { title: { text: 'Relative Flux (ppm)', font: { size: 11, color: '#94a3b8' } }, tickfont: { size: 10, color: '#64748b' }, gridcolor: '#1e293b' },
        shapes: [
          {
            type: 'rect',
            x0: ev.center_btjd - halfW,
            x1: ev.center_btjd + halfW,
            y0: 0,
            y1: 1,
            yref: 'paper',
            fillcolor: '#f59e0b',
            opacity: 0.15,
            line: { width: 0 }
          },
          {
            type: 'line',
            x0: ev.center_btjd,
            x1: ev.center_btjd,
            y0: 0,
            y1: 1,
            yref: 'paper',
            line: { color: '#ef4444', width: 1.5, dash: 'dash' }
          }
        ],
        showlegend: false
      };

      Plotly.react('plot-zoom-container', traces, layout, { responsive: true, displayModeBar: false });
    }

    function renderPlotlyPhase() {
      const lc = AppState.activeLightCurve;
      const sys = AppState.activeStar ? AppState.activeStar.system : null;
      if (!lc || !lc.time || !sys || !sys.period_days) return;

      const p = sys.period_days;
      const t0 = sys.epoch_btjd || lc.time[0];

      const phases = [];
      const fluxes = [];
      const colors = [];

      for (let i = 0; i < lc.time.length; i++) {
        const ph = ((lc.time[i] - t0) / p) % 1.0;
        const normPh = ph > 0.5 ? ph - 1.0 : (ph < -0.5 ? ph + 1.0 : ph);
        const cycle = Math.round((lc.time[i] - t0) / p);
        phases.push(normPh);
        fluxes.push(lc.flux_ppm[i]);
        colors.push(cycle % 2 === 0 ? '#38bdf8' : '#f59e0b');
      }

      const traces = [{
        x: phases,
        y: fluxes,
        mode: 'markers',
        type: 'scatter',
        marker: { size: 4, color: colors, opacity: 0.6 },
        name: 'Phase Folded'
      }];

      const layout = {
        paper_bgcolor: '#090d16',
        plot_bgcolor: '#090d16',
        margin: { l: 50, r: 20, t: 20, b: 35 },
        xaxis: { range: [-0.5, 0.5], title: { text: 'Phase (orbital cycles)', font: { size: 11, color: '#94a3b8' } }, tickfont: { size: 10, color: '#64748b' }, gridcolor: '#1e293b' },
        yaxis: { title: { text: 'Relative Flux (ppm)', font: { size: 11, color: '#94a3b8' } }, tickfont: { size: 10, color: '#64748b' }, gridcolor: '#1e293b' },
        shapes: [
          { type: 'line', x0: 0, x1: 0, y0: 0, y1: 1, yref: 'paper', line: { color: '#64748b', dash: 'dash', width: 1 } }
        ],
        showlegend: false
      };

      Plotly.react('plot-phase-container', traces, layout, { responsive: true, displayModeBar: false });
    }

    function updateEventZoom() {
      renderPlotlyZoom();
    }

    // Review Actions with Scope Isolation
    async function saveStarReview() {
      if (!AppState.activeStar) return;
      const status = document.getElementById('star-rev-status').value;
      const starClass = document.getElementById('star-rev-class').value;
      const note = document.getElementById('star-rev-note').value;
      const bulkApply = document.getElementById('star-bulk-apply').checked;

      const endpoint = bulkApply ? '/api/star/bulk-review' : '/api/star/review';
      const payload = {
        tic: AppState.activeStar.tic,
        status: status,
        star_class: starClass,
        note: note,
        only_unreviewed: true
      };

      try {
        const resp = await fetch(endpoint, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(payload)
        });
        const res = await resp.json();
        if (res.success) {
          showToast(bulkApply ? `Star classified & ${res.updated_count} events updated!` : 'Star review updated successfully!');
          await loadStarDetails(AppState.activeStar.tic, AppState.activeEvent ? AppState.activeEvent.id : null);
          fetchStars();
        } else {
          alert('Error saving star review: ' + (res.error || 'Unknown'));
        }
      } catch (e) {
        alert('Network error: ' + e);
      }
    }

    async function saveEventReview() {
      if (!AppState.activeEvent) return;
      const status = document.getElementById('event-rev-status').value;
      const note = document.getElementById('event-rev-note').value;

      try {
        const resp = await fetch('/api/review', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
            candidate_id: AppState.activeEvent.id,
            status: status,
            note: note
          })
        });
        const res = await resp.json();
        if (res.success) {
          showToast(`Event #${AppState.activeEvent.id} review saved!`);
          await loadStarDetails(AppState.activeStar.tic, AppState.activeEvent.id);
          fetchStars();
        } else {
          alert('Error saving event review: ' + (res.error || 'Unknown'));
        }
      } catch (e) {
        alert('Network error: ' + e);
      }
    }

    async function undoStarReview() {
      if (!AppState.activeStar) return;
      try {
        const resp = await fetch('/api/star/undo', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ tic: AppState.activeStar.tic })
        });
        const res = await resp.json();
        if (res.success) {
          showToast(`Reverted star review to '${res.undo.reverted_to}'`);
          await loadStarDetails(AppState.activeStar.tic, AppState.activeEvent ? AppState.activeEvent.id : null);
          fetchStars();
        } else {
          alert('Undo failed: ' + (res.error || 'Unknown'));
        }
      } catch (e) {
        alert('Network error: ' + e);
      }
    }

    async function undoEventReview() {
      if (!AppState.activeEvent) return;
      try {
        const resp = await fetch('/api/review/undo', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ candidate_id: AppState.activeEvent.id })
        });
        const res = await resp.json();
        if (res.success) {
          showToast(`Reverted event #${AppState.activeEvent.id} to '${res.undo.reverted_to}'`);
          await loadStarDetails(AppState.activeStar.tic, AppState.activeEvent.id);
          fetchStars();
        } else {
          alert('Undo failed: ' + (res.error || 'Unknown'));
        }
      } catch (e) {
        alert('Network error: ' + e);
      }
    }

    // Responsive Plot Auto-Resize Observer
    if (typeof ResizeObserver !== 'undefined') {
      const plotResizeObserver = new ResizeObserver(() => {
        ['plot-full-container', 'plot-zoom-container', 'plot-phase-container'].forEach(id => {
          const el = document.getElementById(id);
          if (el && typeof Plotly !== 'undefined' && el.data) {
            Plotly.Plots.resize(el);
          }
        });
      });
      ['plot-full-container', 'plot-zoom-container', 'plot-phase-container'].forEach(id => {
        const el = document.getElementById(id);
        if (el) plotResizeObserver.observe(el);
      });
    }
  </script>
</body>
</html>
"""
