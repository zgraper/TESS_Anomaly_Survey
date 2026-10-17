"""Persistent Experiment and Survey Run Management for TESS Anomaly Survey.

Provides:
1. ExperimentStore: SQLite-backed registry of research experiments and survey runs.
2. Historical source discovery: Seamless read-only mounting of Phase 7, Phase 6B,
   and legacy research databases without altering source files.
3. SurveyExecutionManager: Thread-safe background execution of streaming surveys
   with real-time progress callbacks, ETA estimation, cancellation, and error logging.
4. Target validation and pre-flight estimation.
5. Catalog and provenance export (CSV/JSON).
"""
from __future__ import annotations

import contextlib
import csv
import io
import json
import math
import os
from pathlib import Path
import re
import sqlite3
import threading
import time
from typing import Any, Dict, List, Optional, Set, Tuple, Union
import uuid

from . import __version__, core, streaming

# Detector presets and constants
FROZEN_THRESHOLD_T = 73.77965536097444
FROZEN_METHOD = "fixed_box_2h"
FROZEN_DURATION_HR = 2.0
FROZEN_CONFIG_HASH = streaming.FROZEN_CONFIG_HASH

DETECTOR_PRESETS = {
    "frozen_fixed_box_2h": {
        "name": "Frozen Fixed Box 2h (Recommended)",
        "threshold_t": FROZEN_THRESHOLD_T,
        "duration_hr": FROZEN_DURATION_HR,
        "method": FROZEN_METHOD,
        "is_frozen": True,
        "calibration_source": "Inherited empirical cutoff; not a catalog false-alarm probability",
        "b0_sha256": core.ORIGINAL_B0_SHA256,
        "b1_sha256": core.ORIGINAL_B1_SHA256,
        "limitations": (
            "Calibrated exclusively for 2-hour localized negative flux excursions. "
            "Does not perform transit fitting, stellar modeling, or flare decontamination."
        )
    },
    "experimental_custom": {
        "name": "Experimental Custom Configuration",
        "threshold_t": FROZEN_THRESHOLD_T,
        "duration_hr": FROZEN_DURATION_HR,
        "method": "fixed_box_2h",
        "is_frozen": False,
        "calibration_source": "User-specified experimental parameters",
        "b0_sha256": core.ORIGINAL_B0_SHA256,
        "b1_sha256": core.ORIGINAL_B1_SHA256,
        "limitations": "Experimental configuration — Not part of the published frozen methodology."
    }
}


def sanitize_tic(val: Any) -> str:
    """Extract clean numeric TIC identifier."""
    s = str(val or "").strip()
    s = re.sub(r"^TIC:?", "", s, flags=re.IGNORECASE).strip()
    return re.sub(r"\D", "", s)


def utc_now() -> str:
    """Return ISO 8601 UTC timestamp."""
    return streaming.utc_now()


# ---------------------------------------------------------------------------
# Pre-Flight Target Validation & Resource Estimator
# ---------------------------------------------------------------------------

def validate_and_estimate_targets(
    targets_raw: str,
    sectors: Optional[List[Union[int, str]]] = None
) -> Dict[str, Any]:
    """Parse, validate, and estimate computational requirements for a target list.
    
    Accepts:
    - Comma-, whitespace-, newline-, or semicolon-delimited TIC numbers.
    - CSV text with a header containing 'tic', 'tic_id', or 'target'.
    """
    lines = [ln.strip() for ln in targets_raw.replace("\r", "").split("\n") if ln.strip()]
    raw_tokens: List[str] = []

    # Check if input looks like a CSV with headers
    is_csv_header = False
    if lines and ("," in lines[0] or "\t" in lines[0]):
        header = [c.strip().lower() for c in re.split(r"[,;\t]", lines[0])]
        tic_col_idx = None
        for idx, col in enumerate(header):
            if col in ("tic", "tic_id", "ticid", "target"):
                tic_col_idx = idx
                break
        if tic_col_idx is not None:
            is_csv_header = True
            for line in lines[1:]:
                parts = re.split(r"[,;\t]", line)
                if len(parts) > tic_col_idx:
                    raw_tokens.append(parts[tic_col_idx].strip())

    if not is_csv_header:
        for line in lines:
            # Strip TIC: or TIC prefix followed by spaces/colons
            norm_line = re.sub(r"\bTIC[:\s]*", "", line, flags=re.IGNORECASE)
            # Split line on commas, semicolons, or whitespace
            tokens = re.split(r"[,;\s\t]+", norm_line)
            for tok in tokens:
                t = tok.strip()
                if t:
                    raw_tokens.append(t)

    valid_tics: List[str] = []
    invalid_tokens: List[Dict[str, str]] = []
    seen_tics: Set[str] = set()
    dup_count = 0

    for tok in raw_tokens:
        clean = sanitize_tic(tok)
        if not clean:
            invalid_tokens.append({"token": tok, "reason": "Non-numeric target identifier"})
            continue
        try:
            val = int(clean)
            if val < 1 or val > 2_000_000_000:
                invalid_tokens.append({"token": tok, "reason": "Out of valid TESS TIC ID range (1 to 2,000,000,000)"})
                continue
        except ValueError:
            invalid_tokens.append({"token": tok, "reason": "Malformed integer format"})
            continue

        if clean in seen_tics:
            dup_count += 1
            continue

        seen_tics.add(clean)
        valid_tics.append(clean)

    # Validate Sectors
    parsed_sectors: List[int] = []
    if sectors:
        for s in sectors:
            try:
                sec_int = int(str(s).strip())
                if 1 <= sec_int <= 99:
                    if sec_int not in parsed_sectors:
                        parsed_sectors.append(sec_int)
            except ValueError:
                pass
    parsed_sectors.sort()

    n_stars = len(valid_tics)
    n_sectors = max(1, len(parsed_sectors))
    estimated_observations = n_stars * n_sectors
    estimated_download_mb = round(estimated_observations * 2.1, 1) # ~2.1 MB per SPOC 2-min cadence FITS
    estimated_process_sec = round(estimated_observations * 0.55, 1) # ~0.55s per obs with 2 workers

    warnings: List[str] = []
    if n_stars == 0:
        warnings.append("No valid target TIC identifiers found in input.")
    elif n_stars > 50:
        warnings.append(
            f"Target list contains {n_stars} stars (~{estimated_observations} observations). "
            f"Interactive browser surveys are best suited for <= 50 targets. "
            f"Expected download footprint is {estimated_download_mb} MB."
        )

    return {
        "valid_tics": valid_tics,
        "invalid_tokens": invalid_tokens,
        "duplicate_count": dup_count,
        "total_unique_stars": n_stars,
        "selected_sectors": parsed_sectors,
        "estimated_observations": estimated_observations,
        "estimated_download_mb": estimated_download_mb,
        "estimated_duration_sec": estimated_process_sec,
        "warnings": warnings,
        "is_safe": n_stars > 0
    }


# ---------------------------------------------------------------------------
# ExperimentStore: SQLite Persistence for Experiments and Runs
# ---------------------------------------------------------------------------

def read_source_run(path, run_id=None):
    """Read recorded execution metadata; never infer it from a workspace preset."""
    if not Path(path).is_file():
        return None
    try:
        with sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro', uri=True) as cx:
            cx.row_factory = sqlite3.Row
            cols = {r[1] for r in cx.execute('PRAGMA table_info(runs)')}
            key = 'run_id' if 'run_id' in cols else 'id'
            rows = cx.execute(f'SELECT * FROM runs WHERE CAST({key} AS TEXT)=?', (str(run_id),)).fetchall() if run_id else cx.execute('SELECT * FROM runs LIMIT 2').fetchall()
            return dict(rows[0]) if len(rows) == 1 else None
    except sqlite3.Error:
        return None


class ExperimentStore:
    """Manages persistent research experiments, survey runs, and historical imports."""

    def __init__(self, workspace_dir: Union[str, Path] = "experiments", db_path: Optional[str] = None):
        self.workspace_dir = Path(workspace_dir).resolve()
        self.workspace_dir.mkdir(parents=True, exist_ok=True)
        self.overlays_dir = self.workspace_dir / "overlays"
        self.overlays_dir.mkdir(parents=True, exist_ok=True)

        self.db_path = str(Path(db_path).resolve()) if db_path else str(self.workspace_dir / "workspace.sqlite")
        self._lock = threading.Lock()
        self._init_db()
        self._register_historical_sources()

    @contextlib.contextmanager
    def _connection(self):
        cx = sqlite3.connect(self.db_path, timeout=30.0)
        cx.row_factory = sqlite3.Row
        try:
            with cx:
                yield cx
        finally:
            cx.close()

    def _init_db(self):
        with self._lock, self._connection() as cx:
            cx.executescript("""
            PRAGMA foreign_keys=ON;

            CREATE TABLE IF NOT EXISTS experiments (
                id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT NOT NULL DEFAULT '',
                scientific_objective TEXT NOT NULL DEFAULT '',
                is_historical INTEGER NOT NULL DEFAULT 0,
                is_read_only INTEGER NOT NULL DEFAULT 0,
                config_preset TEXT NOT NULL DEFAULT 'frozen_fixed_box_2h',
                threshold_t REAL NOT NULL,
                duration_hr REAL NOT NULL,
                config_hash TEXT NOT NULL,
                software_version TEXT NOT NULL,
                target_count INTEGER NOT NULL DEFAULT 0,
                sectors_json TEXT NOT NULL DEFAULT '[]',
                targets_json TEXT NOT NULL DEFAULT '[]',
                created_utc TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS runs (
                run_id TEXT PRIMARY KEY,
                experiment_id TEXT NOT NULL,
                run_label TEXT NOT NULL DEFAULT '',
                status TEXT NOT NULL DEFAULT 'pending', -- pending, running, completed, failed, aborted
                output_db_path TEXT NOT NULL,
                review_db_path TEXT,
                total_targets INTEGER NOT NULL DEFAULT 0,
                processed_targets INTEGER NOT NULL DEFAULT 0,
                completed_obs INTEGER NOT NULL DEFAULT 0,
                failed_obs INTEGER NOT NULL DEFAULT 0,
                candidate_count INTEGER NOT NULL DEFAULT 0,
                verified_count INTEGER NOT NULL DEFAULT 0,
                throughput_obs_hr REAL NOT NULL DEFAULT 0.0,
                elapsed_sec REAL NOT NULL DEFAULT 0.0,
                eta_sec REAL,
                peak_ram_mb REAL NOT NULL DEFAULT 0.0,
                avg_cpu_pct REAL NOT NULL DEFAULT 0.0,
                started_utc TEXT,
                ended_utc TEXT,
                error_message TEXT,
                settings_json TEXT NOT NULL DEFAULT '{}',
                created_utc TEXT NOT NULL,
                FOREIGN KEY(experiment_id) REFERENCES experiments(id)
            );

            CREATE TABLE IF NOT EXISTS run_logs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                run_id TEXT NOT NULL,
                level TEXT NOT NULL DEFAULT 'INFO',
                message TEXT NOT NULL,
                created_utc TEXT NOT NULL,
                FOREIGN KEY(run_id) REFERENCES runs(run_id)
            );

            CREATE TABLE IF NOT EXISTS workspace_state (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );
            """)

    def _register_historical_sources(self):
        """Auto-discover historical research databases and register as read-only experiments."""
        historical_definitions = [
            {
                "id": "exp_phase7_development",
                "name": "Phase 7 Development Survey (Sectors 1–10)",
                "description": "Empirical survey across 9,753 observations in Sectors 1–10.",
                "scientific_objective": (
                    "Empirical localized dimming survey using frozen fixed_box_2h detector "
                    "with cross-camera and cross-star coincidence auditing."
                ),
                "run_id": "run_phase7_9753",
                "candidate_paths": [
                    Path(os.environ.get("TESS_SURVEY_DIR", "")) / "development_survey_9753.sqlite",
                    Path.home() / "Projects" / "TESS_Survey_Results" / "development_survey_9753.sqlite",
                    Path.home() / ".tess_anomaly_survey" / "databases" / "development_survey_9753.sqlite",
                    Path("research/phase7/development_survey_9753.sqlite")
                ],
                "sectors": list(range(1, 11)),
                "target_count": 8984
            },
            {
                "id": "exp_phase6b_validation",
                "name": "Phase 6B Independent Validation Study",
                "description": "Preregistered independent validation study (4,000 targets).",
                "scientific_objective": (
                    "Independent validation of frozen threshold T = 73.78 on unseen in-sector "
                    "and cross-sector (Sector 11) targets."
                ),
                "run_id": "run_phase6b_4000",
                "candidate_paths": [
                    Path(os.environ.get("TESS_SURVEY_DIR", "")) / "phase6b_validation.sqlite",
                    Path.home() / "Projects" / "TESS_Survey_Results" / "phase6b_validation.sqlite",
                    Path.home() / ".tess_anomaly_survey" / "databases" / "phase6b_validation.sqlite",
                    Path("research/phase6b/phase6b_validation.sqlite")
                ],
                "sectors": list(range(1, 12)),
                "target_count": 4000
            },
            {
                "id": "exp_legacy_library",
                "name": "Phase 3/4 Archival Anomaly Library",
                "description": "Archival candidate and benchmark system database.",
                "scientific_objective": "Initial development benchmark and early anomaly library.",
                "run_id": "run_legacy_107",
                "candidate_paths": [
                    Path(os.environ.get("TESS_SURVEY_DIR", "")) / "anomaly_library.sqlite",
                    Path.home() / "Projects" / "TESS_Survey_Results" / "anomaly_library.sqlite",
                    Path.home() / ".tess_anomaly_survey" / "databases" / "anomaly_library.sqlite",
                    Path("anomaly_library.sqlite")
                ],
                "sectors": [1, 2],
                "target_count": 107
            },
            {
                "id": "exp_streaming_pilot",
                "name": "Phase 7 Streaming Pilot Survey",
                "description": "Pilot streaming run (Sectors 1–10) validating zero-disk architecture.",
                "scientific_objective": "Pilot validation of streaming pipeline and worker performance.",
                "run_id": "run_streaming_pilot",
                "candidate_paths": [
                    Path(os.environ.get("TESS_SURVEY_DIR", "")) / "streaming_pilot.sqlite",
                    Path.home() / "Projects" / "TESS_Survey_Results" / "streaming_pilot.sqlite",
                    Path.home() / ".tess_anomaly_survey" / "databases" / "streaming_pilot.sqlite",
                    Path("streaming_pilot.sqlite")
                ],
                "sectors": list(range(1, 11)),
                "target_count": 1008
            },
            {
                "id": "exp_sample_demonstration",
                "name": "Public Demonstration Survey",
                "description": "Self-contained synthetic demonstration with illustrative signals and control targets; not astrophysical ground truth.",
                "scientific_objective": (
                    "Demonstration survey featuring an illustrative synthetic eclipsing binary benchmark (TIC 39018208, simulated P = 2.41824 d), "
                    "localized dimming candidate (TIC 73834), pulsator (TIC 60232), flare artifact (TIC 141944), and quiet control star (TIC 999001)."
                ),
                "run_id": "run_sample_demo",
                "candidate_paths": [
                    Path("examples/sample_survey_review.sqlite"),
                    Path(__file__).resolve().parent.parent.parent / "examples" / "sample_survey_review.sqlite"
                ],
                "sectors": [1, 2, 4, 5, 8],
                "target_count": 5
            }
        ]

        with self._lock, self._connection() as cx:
            for item in historical_definitions:
                found_path: Optional[Path] = None
                for cp in item["candidate_paths"]:
                    if cp.exists():
                        found_path = cp.resolve()
                        break
                if not found_path:
                    continue

                # Ensure overlay DB exists
                overlay_path = str(found_path.with_name(f"{found_path.stem}.reviews.sqlite"))

                # Register historical experiment
                cx.execute("""
                INSERT OR REPLACE INTO experiments (
                    id, name, description, scientific_objective, is_historical, is_read_only,
                    config_preset, threshold_t, duration_hr, config_hash, software_version,
                    target_count, sectors_json, targets_json, created_utc
                ) VALUES (?, ?, ?, ?, 1, 1, 'frozen_fixed_box_2h', ?, ?, ?, ?, ?, ?, '[]', ?)
                """, (
                    item["id"], item["name"], item["description"], item["scientific_objective"],
                    FROZEN_THRESHOLD_T, FROZEN_DURATION_HR, FROZEN_CONFIG_HASH, __version__,
                    item["target_count"], json.dumps(item["sectors"]), utc_now()
                ))

                # Count stats from historical database
                n_cands = 0
                n_verified = 0
                throughput = 0.0
                try:
                    with sqlite3.connect(f"file:{found_path.as_posix()}?mode=ro", uri=True) as h_cx:
                        cand_tables = {r[0] for r in h_cx.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
                        if "candidates" in cand_tables:
                            n_cands = h_cx.execute("SELECT COUNT(*) FROM candidates").fetchone()[0]
                            cand_cols = {c[1] for c in h_cx.execute("PRAGMA table_info(candidates)").fetchall()}
                            if "scientific_classification" in cand_cols:
                                n_verified = h_cx.execute(
                                    "SELECT COUNT(*) FROM candidates WHERE scientific_classification='verified_localized_dimming'"
                                ).fetchone()[0]
                        if "performance_benchmarks" in cand_tables:
                            pb = h_cx.execute("SELECT avg_obs_per_hour FROM performance_benchmarks ORDER BY created_utc DESC LIMIT 1").fetchone()
                            if pb:
                                throughput = float(pb[0])
                except Exception:
                    pass

                # Register historical run
                cx.execute("""
                INSERT OR REPLACE INTO runs (
                    run_id, experiment_id, run_label, status, output_db_path, review_db_path,
                    total_targets, processed_targets, completed_obs, failed_obs,
                    candidate_count, verified_count, throughput_obs_hr, elapsed_sec,
                    started_utc, ended_utc, settings_json, created_utc
                ) VALUES (?, ?, 'Historical Run', 'completed', ?, ?, ?, ?, ?, 0, ?, ?, ?, 0.0, ?, ?, '{}', ?)
                """, (
                    item["run_id"], item["id"], str(found_path), overlay_path,
                    item["target_count"], item["target_count"], item["target_count"],
                    n_cands, n_verified, throughput, utc_now(), utc_now(), utc_now()
                ))

            # One workspace record per native source run; annotations live in the workspace.
            cols = {r[1] for r in cx.execute("PRAGMA table_info(runs)")}
            if 'source_run_id' not in cols:
                cx.execute("ALTER TABLE runs ADD COLUMN source_run_id TEXT")
            imported = cx.execute("SELECT r.* FROM runs r JOIN experiments e ON e.id=r.experiment_id WHERE e.is_historical=1").fetchall()
            for original in imported:
                if '_source_' in original['run_id']:
                    continue
                record = dict(original)
                with sqlite3.connect(Path(record['output_db_path']).as_uri() + '?mode=ro', uri=True) as source:
                    tables = {r[0] for r in source.execute("SELECT name FROM sqlite_master WHERE type='table'")}
                    native_ids = []
                    if 'runs' in tables:
                        rc = {r[1] for r in source.execute("PRAGMA table_info(runs)")}
                        key = 'run_id' if 'run_id' in rc else 'id'
                        native_ids = [str(r[0]) for r in source.execute(f"SELECT {key} FROM runs ORDER BY {key}")]
                    for index, native_id in enumerate(native_ids):
                        entry = dict(record)
                        if index:
                            entry['run_id'] = record['run_id'] + '_source_' + str(index)
                        entry['source_run_id'] = native_id
                        entry['run_label'] = 'Historical: ' + native_id
                        entry['review_db_path'] = str(self.workspace_dir / 'overlays' / (entry['run_id'] + '.reviews.sqlite'))
                        entry['candidate_count'] = source.execute("SELECT COUNT(*) FROM candidates WHERE CAST(run_id AS TEXT)=?", (native_id,)).fetchone()[0]
                        cand_cols = {r[1] for r in source.execute('PRAGMA table_info(candidates)')}
                        entry['verified_count'] = source.execute("SELECT COUNT(*) FROM candidates WHERE CAST(run_id AS TEXT)=? AND scientific_classification='verified_localized_dimming'", (native_id,)).fetchone()[0] if 'scientific_classification' in cand_cols else 0
                        # Historical totals describe recorded observations for this run,
                        # never the nominal size of a different cohort or native run.
                        entry.update(total_targets=0, processed_targets=0, completed_obs=0,
                                     failed_obs=0, throughput_obs_hr=0.0)
                        if 'observations' in tables:
                            obs_cols = {r[1] for r in source.execute('PRAGMA table_info(observations)')}
                            if {'run_id', 'status'} <= obs_cols:
                                counts = dict(source.execute('SELECT status,COUNT(*) FROM observations WHERE CAST(run_id AS TEXT)=? GROUP BY status', (native_id,)))
                                entry['completed_obs'] = sum(counts.get(s, 0) for s in ('completed', 'no_anomalies'))
                                entry['failed_obs'] = counts.get('failed', 0)
                                entry['total_targets'] = sum(counts.values())
                                entry['processed_targets'] = entry['completed_obs'] + entry['failed_obs']
                        elif 'processed_files' in tables:
                            entry['completed_obs'] = source.execute('SELECT COUNT(*) FROM processed_files WHERE CAST(run_id AS TEXT)=?', (native_id,)).fetchone()[0]
                            entry['total_targets'] = entry['processed_targets'] = entry['completed_obs']
                        keys = list(entry)
                        cx.execute('INSERT OR REPLACE INTO runs (' + ','.join(keys) + ') VALUES (' + ','.join('?' for _ in keys) + ')', [entry[k] for k in keys])

            # Set default active experiment if unset
            cur_active = cx.execute("SELECT value FROM workspace_state WHERE key='active_experiment'").fetchone()
            if not cur_active:
                first_exp = cx.execute("SELECT id FROM experiments ORDER BY is_historical ASC, created_utc DESC LIMIT 1").fetchone()
                if first_exp:
                    first_run = cx.execute("SELECT run_id FROM runs WHERE experiment_id=? LIMIT 1", (first_exp[0],)).fetchone()
                    if first_run:
                        cx.execute("INSERT OR REPLACE INTO workspace_state VALUES ('active_experiment', ?)", (first_exp[0],))
                        cx.execute("INSERT OR REPLACE INTO workspace_state VALUES ('active_run', ?)", (first_run[0],))

    def list_experiments(self) -> List[Dict[str, Any]]:
        with self._lock, self._connection() as cx:
            rows = cx.execute("""
                SELECT e.*, COUNT(r.run_id) as run_count
                FROM experiments e
                LEFT JOIN runs r ON e.id = r.experiment_id
                GROUP BY e.id
                ORDER BY e.is_historical ASC, e.created_utc DESC
            """).fetchall()

            res = []
            for r in rows:
                d = dict(r)
                d["sectors"] = json.loads(d["sectors_json"]) if d["sectors_json"] else []
                d["targets"] = json.loads(d["targets_json"]) if d["targets_json"] else []
                res.append(d)
            return res

    def get_experiment(self, exp_id: str) -> Optional[Dict[str, Any]]:
        with self._lock, self._connection() as cx:
            row = cx.execute("SELECT * FROM experiments WHERE id=?", (exp_id,)).fetchone()
            if not row:
                return None
            exp = dict(row)
            exp["sectors"] = json.loads(exp["sectors_json"]) if exp["sectors_json"] else []
            exp["targets"] = json.loads(exp["targets_json"]) if exp["targets_json"] else []

            run_rows = cx.execute("""
                SELECT * FROM runs WHERE experiment_id=? ORDER BY created_utc DESC
            """, (exp_id,)).fetchall()
            exp["runs"] = [dict(rr) for rr in run_rows]
            return exp

    def create_experiment(
        self,
        name: str,
        description: str,
        scientific_objective: str,
        targets: List[str],
        sectors: List[int],
        config_preset: str = "frozen_fixed_box_2h",
        threshold_t: Optional[float] = None,
        duration_hr: Optional[float] = None,
        is_read_only: bool = False
    ) -> Dict[str, Any]:
        """Create a new scientific experiment definition."""
        preset = DETECTOR_PRESETS.get(config_preset, DETECTOR_PRESETS["frozen_fixed_box_2h"])
        eff_threshold = float(threshold_t) if threshold_t is not None else float(preset["threshold_t"])
        eff_duration = float(duration_hr) if duration_hr is not None else float(preset["duration_hr"])
        if eff_duration != 2.0:
            raise ValueError("Only a two-hour template is supported; custom duration is disabled.")
        if not math.isfinite(eff_threshold) or eff_threshold <= 0:
            raise ValueError("Threshold must be finite and positive.")

        exp_slug = re.sub(r"[^a-zA-Z0-9_]+", "_", name.strip().lower())[:30]
        ts_id = time.strftime("%Y%m%d_%H%M%S")
        uid_hex = uuid.uuid4().hex[:6]
        exp_id = f"exp_{ts_id}_{uid_hex}_{exp_slug}" if exp_slug else f"exp_{ts_id}_{uid_hex}"

        # Software & config hash: strictly only assigned if threshold, duration, and preset match published frozen detector
        is_truly_frozen = (
            abs(eff_threshold - FROZEN_THRESHOLD_T) < 1e-6
            and abs(eff_duration - FROZEN_DURATION_HR) < 1e-6
            and config_preset == "frozen_fixed_box_2h"
        )
        cfg_hash = streaming.execution_configuration(eff_threshold)[2]

        with self._lock, self._connection() as cx:
            cx.execute("""
            INSERT INTO experiments (
                id, name, description, scientific_objective, is_historical, is_read_only,
                config_preset, threshold_t, duration_hr, config_hash, software_version,
                target_count, sectors_json, targets_json, created_utc
            ) VALUES (?, ?, ?, ?, 0, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                exp_id, name, description, scientific_objective,
                1 if is_read_only else 0,
                config_preset, eff_threshold, eff_duration, cfg_hash, __version__,
                len(targets), json.dumps(sectors), json.dumps(targets), utc_now()
            ))

            # Automatically create the initial run for this experiment
            run_id = f"run_{ts_id}_{uid_hex}_01"
            exp_dir = self.workspace_dir / exp_id
            exp_dir.mkdir(parents=True, exist_ok=True)
            output_db = exp_dir / f"{run_id}.sqlite"
            review_db = exp_dir / f"{run_id}.reviews.sqlite"

            cx.execute("""
            INSERT INTO runs (
                run_id, experiment_id, run_label, status, output_db_path, review_db_path,
                total_targets, processed_targets, completed_obs, failed_obs,
                candidate_count, verified_count, throughput_obs_hr, elapsed_sec,
                settings_json, created_utc
            ) VALUES (?, ?, 'Initial Survey Run', 'pending', ?, ?, ?, 0, 0, 0, 0, 0, 0.0, 0.0, ?, ?)
            """, (
                run_id, exp_id, str(output_db), str(review_db) if review_db else None,
                len(targets) * max(1, len(sectors)),
                json.dumps({"max_workers": 2, "threshold_t": eff_threshold, "duration_hr": eff_duration}),
                utc_now()
            ))

            # Update active selection to this new experiment
            cx.execute("INSERT OR REPLACE INTO workspace_state VALUES ('active_experiment', ?)", (exp_id,))
            cx.execute("INSERT OR REPLACE INTO workspace_state VALUES ('active_run', ?)", (run_id,))

        return self.get_experiment(exp_id) # type: ignore

    def list_runs(self, experiment_id: Optional[str] = None) -> List[Dict[str, Any]]:
        with self._lock, self._connection() as cx:
            if experiment_id:
                rows = cx.execute("SELECT * FROM runs WHERE experiment_id=? ORDER BY created_utc DESC", (experiment_id,)).fetchall()
            else:
                rows = cx.execute("SELECT * FROM runs ORDER BY created_utc DESC").fetchall()
            return [dict(r) for r in rows]

    def get_run(self, run_id: str) -> Optional[Dict[str, Any]]:
        with self._lock, self._connection() as cx:
            row = cx.execute("SELECT * FROM runs WHERE run_id=?", (run_id,)).fetchone()
            result = dict(row) if row else None
        if result:
            result['scientific_provenance'] = read_source_run(result['output_db_path'], result.get('source_run_id'))
        return result

    def create_run(
        self,
        experiment_id: str,
        run_label: str = "",
        settings: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """Create a new run within an existing experiment."""
        exp = self.get_experiment(experiment_id)
        if not exp:
            raise ValueError(f"Experiment {experiment_id} not found.")

        ts_id = time.strftime("%Y%m%d_%H%M%S")
        uid_hex = uuid.uuid4().hex[:6]
        run_id = f"run_{ts_id}_{uid_hex}_{len(exp.get('runs', [])) + 1:02d}"
        exp_dir = self.workspace_dir / experiment_id
        exp_dir.mkdir(parents=True, exist_ok=True)
        output_db = exp_dir / f"{run_id}.sqlite"
        review_db = exp_dir / f"{run_id}.reviews.sqlite"

        total_obs = exp.get("target_count", 0) * max(1, len(exp.get("sectors", [])))

        with self._lock, self._connection() as cx:
            cx.execute("""
            INSERT INTO runs (
                run_id, experiment_id, run_label, status, output_db_path, review_db_path,
                total_targets, processed_targets, completed_obs, failed_obs,
                candidate_count, verified_count, throughput_obs_hr, elapsed_sec,
                settings_json, created_utc
            ) VALUES (?, ?, ?, 'pending', ?, ?, ?, 0, 0, 0, 0, 0, 0.0, 0.0, ?, ?)
            """, (
                run_id, experiment_id, run_label or f"Run {ts_id}",
                str(output_db), str(review_db) if review_db else None,
                total_obs, json.dumps(settings or {}), utc_now()
            ))

        return self.get_run(run_id) # type: ignore

    def update_run(self, run_id: str, **kwargs):
        """Update run execution state and telemetry."""
        allowed = {
            "status", "processed_targets", "completed_obs", "failed_obs",
            "candidate_count", "verified_count", "throughput_obs_hr",
            "elapsed_sec", "eta_sec", "peak_ram_mb", "avg_cpu_pct",
            "started_utc", "ended_utc", "error_message", "settings_json",
            "output_db_path", "review_db_path", "total_targets", "source_run_id"
        }
        updates = {k: v for k, v in kwargs.items() if k in allowed}
        if not updates:
            return

        set_clause = ", ".join(f"{k} = ?" for k in updates.keys())
        values = list(updates.values()) + [run_id]

        with self._lock, self._connection() as cx:
            cx.execute(f"UPDATE runs SET {set_clause} WHERE run_id = ?", values)

    def add_run_log(self, run_id: str, message: str, level: str = "INFO"):
        with self._lock, self._connection() as cx:
            cx.execute(
                "INSERT INTO run_logs (run_id, level, message, created_utc) VALUES (?, ?, ?, ?)",
                (run_id, level, message, utc_now())
            )

    def get_run_logs(self, run_id: str, limit: int = 50) -> List[Dict[str, Any]]:
        with self._lock, self._connection() as cx:
            rows = cx.execute("""
                SELECT level, message, created_utc
                FROM run_logs
                WHERE run_id=?
                ORDER BY id DESC LIMIT ?
            """, (run_id, limit)).fetchall()
            return [dict(r) for r in reversed(rows)]

    def get_active_selection(self) -> Tuple[Optional[str], Optional[str]]:
        with self._lock, self._connection() as cx:
            exp_row = cx.execute("SELECT value FROM workspace_state WHERE key='active_experiment'").fetchone()
            run_row = cx.execute("SELECT value FROM workspace_state WHERE key='active_run'").fetchone()
            return (exp_row[0] if exp_row else None, run_row[0] if run_row else None)

    def set_active_selection(self, experiment_id: str, run_id: str):
        run = self.get_run(run_id)
        if not run or run["experiment_id"] != experiment_id:
            raise ValueError("Run does not belong to experiment")
        with self._lock, self._connection() as cx:
            cx.execute("INSERT OR REPLACE INTO workspace_state VALUES ('active_experiment', ?)", (experiment_id,))
            cx.execute("INSERT OR REPLACE INTO workspace_state VALUES ('active_run', ?)", (run_id,))

    def export_run_catalog(self, experiment_id: str, run_id: str, fmt: str = "csv") -> Tuple[str, str]:
        """Export experiment candidate catalog and configuration provenance as CSV or JSON.
        
        Returns (content_string, mime_type).
        """
        exp = self.get_experiment(experiment_id)
        run = self.get_run(run_id)
        if not exp or not run:
            raise ValueError("Experiment or Run not found")
        if run["experiment_id"] != experiment_id:
            raise ValueError("Run does not belong to experiment")

        output_db = run["output_db_path"]
        review_db = run.get("review_db_path")
        read_only = exp.get("is_read_only", 0) == 1

        # Introspect and extract candidate rows
        records: List[Dict[str, Any]] = []
        actual_provenance = {}
        if os.path.exists(output_db):
            from .web import SurveyDatabaseAdapter
            adapter = SurveyDatabaseAdapter(output_db, review_db_path=review_db, read_only=read_only, run_id=run.get("source_run_id"))
            cx = adapter.get_connection()
            try:
                cand_cols = {col[1] for col in cx.execute("PRAGMA table_info(candidates)").fetchall()}
                verif_col = "c.scientific_classification" if "scientific_classification" in cand_cols else "'not_run'"
                if "scientific_reason" in cand_cols:
                    reason_col = "c.scientific_reason"
                elif "raw_reason" in cand_cols:
                    reason_col = "c.raw_reason"
                else:
                    reason_col = "''"

                status_col = "c.review_status" if "review_status" in cand_cols else ("c.status" if "status" in cand_cols else "'unreviewed'")
                note_col = "c.note" if "note" in cand_cols else ("c.review_note" if "review_note" in cand_cols else "''")

                ov_cand_join = "LEFT JOIN overlay.candidate_review_overlay o_rev ON c.id = o_rev.candidate_id" if adapter.review_db_path else ""
                ov_star_join = "LEFT JOIN overlay.star_reviews sr_ov ON c.tic = sr_ov.tic" if adapter.review_db_path else ""

                sql = f"""
                SELECT c.id, c.tic, c.sector, c.center_btjd, c.score, c.width_hr, c.depth_ppm,
                       {verif_col} as automated_classification,
                       {reason_col} as automated_reason,
                       COALESCE({'o_rev.status, ' if adapter.review_db_path else ''}{status_col}, 'unreviewed') as human_event_review,
                       COALESCE({'o_rev.note, ' if adapter.review_db_path else ''}{note_col}, '') as human_event_note,
                       COALESCE({'sr_ov.status, ' if adapter.review_db_path else ''}sr.status, 'unreviewed') as human_star_review,
                       COALESCE({'sr_ov.star_class, ' if adapter.review_db_path else ''}sr.star_class, '') as human_star_class,
                       COALESCE({'sr_ov.note, ' if adapter.review_db_path else ''}sr.note, '') as human_star_note
                FROM candidates c
                LEFT JOIN star_reviews sr ON c.tic = sr.tic
                {ov_cand_join}
                {ov_star_join}
                ORDER BY c.score DESC
                """
                rows = cx.execute(sql).fetchall()
                for r in rows:
                    records.append({
                        "candidate_id": r[0],
                        "tic": str(r[1]).replace("TIC:", "").strip(),
                        "sector": r[2],
                        "center_btjd": round(float(r[3]), 5),
                        "detector_score": round(float(r[4]), 2),
                        "width_hr": round(float(r[5]), 2),
                        "depth_ppm": round(float(r[6]), 1),
                        "automated_classification": r[7],
                        "automated_reason": r[8] or "",
                        "human_event_review": r[9],
                        "human_event_note": r[10],
                        "human_star_review": r[11],
                        "human_star_class": r[12],
                        "human_star_note": r[13]
                    })
                for record in records:
                    cursor = cx.execute("SELECT * FROM candidates WHERE id=?", (record['candidate_id'],))
                    raw = dict(zip([c[0] for c in cursor.description], cursor.fetchone()))
                    for key in ('run_id', 'observation_id', 'source_file', 'has_shoulder_flare', 'true_depth_min_pct'):
                        record[key] = raw.get(key)
                    record['coincidence_risk_level'] = raw.get('coincidence_risk_level') or 'unknown'
                    record['automated_classification'] = record['automated_classification'] or 'not_run'
                    if raw.get('observation_id') and adapter.schema_info['has_observations']:
                        obs = cx.execute("SELECT product_filename,data_uri,file_sha256 FROM observations WHERE id=?", (raw['observation_id'],)).fetchone()
                        if obs:
                            record.update(zip(('product_filename','data_uri','file_sha256'), obs))
                if adapter.run_id:
                    run_cols = {c[1] for c in cx.execute("PRAGMA table_info(runs)")}
                    key = 'run_id' if 'run_id' in run_cols else 'id'
                    cursor = cx.execute(f"SELECT * FROM runs WHERE CAST({key} AS TEXT)=?", (adapter.run_id,))
                    result = cursor.fetchone()
                    if result:
                        actual_provenance = dict(zip([c[0] for c in cursor.description], result))
            finally:
                cx.close()

        provenance = {
            "experiment_id": exp["id"],
            "experiment_name": exp["name"],
            "scientific_objective": exp["scientific_objective"],
            "run_id": run["run_id"],
            "detector_preset": exp["config_preset"],
            "threshold_t": exp["threshold_t"],
            "template_duration_hr": exp["duration_hr"],
            "config_hash": exp["config_hash"],
            "software_version": exp["software_version"],
            "exported_utc": utc_now(),
            "total_candidates": len(records),
            "scientific_caveat": (
                "Scientific Notice: Detector score S quantifies statistical flux deficit prominence "
                "under the matched 2-hour box filter; high detector scores do NOT necessarily indicate "
                "astrophysical novelty or planetary rarity."
            )
        }

        if actual_provenance:
            for key in ('threshold_t', 'config_hash', 'software_version', 'method', 'settings_json', 'b0_hash', 'b1_hash'):
                provenance[key] = actual_provenance.get(key, 'unavailable')
            provenance['source_run_id'] = actual_provenance.get('run_id', actual_provenance.get('id'))
        provenance['data_kind'] = 'synthetic demonstration' if run.get('source_run_id') == 'sample_survey_run' or Path(output_db).name == 'sample_survey_review.sqlite' or exp['id'] == 'exp_demo' or 'synthetic' in exp.get('description', '').lower() else 'observational / consult source provenance'

        if fmt.lower() == "json":
            payload = {
                "provenance": provenance,
                "candidates": records
            }
            return (json.dumps(payload, indent=2), "application/json")

        # CSV format with commented provenance metadata header
        out = io.StringIO()
        out.write("# ==========================================================================\n")
        out.write("# TESS ANOMALY SURVEY - CANDIDATE CATALOG EXPORT\n")
        for k, v in provenance.items():
            out.write(f"# {k}: {v}\n")
        out.write("# ==========================================================================\n")

        fieldnames = [
            "candidate_id", "tic", "sector", "center_btjd", "detector_score", "width_hr", "depth_ppm",
            "automated_classification", "automated_reason",
            "human_event_review", "human_event_note",
            "human_star_review", "human_star_class", "human_star_note",
            "run_id", "observation_id", "source_file", "product_filename", "data_uri", "file_sha256",
            "has_shoulder_flare", "true_depth_min_pct", "coincidence_risk_level"
        ]
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)

        return (out.getvalue(), "text/csv")


# ---------------------------------------------------------------------------
# SurveyExecutionManager: Background Survey Execution & Telemetry
# ---------------------------------------------------------------------------

class SurveyExecutionManager:
    """Manages thread-safe survey streaming execution, progress, and cancellation."""

    def __init__(self, store: ExperimentStore):
        self.store = store
        self._lock = threading.Lock()
        # Active jobs: run_id -> {thread, cancel_event, telemetry, db_path}
        self._active_jobs: Dict[str, Dict[str, Any]] = {}
        with store._connection() as cx:
            cx.execute("UPDATE runs SET status='interrupted', error_message='Application stopped during execution; resume explicitly.' WHERE status='running'")

    def is_run_active(self, run_id: str) -> bool:
        with self._lock:
            job = self._active_jobs.get(run_id)
            return bool(job and job["thread"].is_alive())

    def get_live_telemetry(self, run_id: str) -> Optional[Dict[str, Any]]:
        with self._lock:
            job = self._active_jobs.get(run_id)
            if job and "telemetry" in job:
                return dict(job["telemetry"])
        run = self.store.get_run(run_id)
        return run

    def abort_survey(self, run_id: str) -> bool:
        with self._lock:
            job = self._active_jobs.get(run_id)
            if not job or not job["thread"].is_alive():
                return False
            job["cancel_event"].set()
        self.store.add_run_log(run_id, "User requested survey execution cancellation.", "WARN")
        return True

    def launch_survey(
        self,
        experiment_id: str,
        run_id: str,
        max_workers: int = 2,
        rate_limit_delay_sec: float = 0.02,
        retain_fits: bool = False,
        resume: bool = True
    ) -> Dict[str, Any]:
        """Launch a persistent background survey streaming run."""
        exp = self.store.get_experiment(experiment_id)
        if not exp:
            raise ValueError(f"Experiment {experiment_id} not found.")

        run = self.store.get_run(run_id)
        if not run:
            raise ValueError(f"Run {run_id} not found.")

        if run["experiment_id"] != experiment_id:
            raise ValueError("Run does not belong to the supplied experiment.")
        if exp.get("is_read_only") or exp.get("is_historical"):
            raise ValueError("Historical and read-only experiments cannot execute surveys.")
        if not 1 <= max_workers <= 16:
            raise ValueError("Worker count must be between 1 and 16.")
        if float(exp["duration_hr"]) != 2.0:
            raise ValueError("Only the fixed two-hour detector is supported.")

        out_db = run["output_db_path"]

        with self._lock:
            # Check if any active run is targeting the same output database
            for j_id, j in self._active_jobs.items():
                if j["thread"].is_alive() and j["db_path"] == out_db:
                    raise RuntimeError(f"Database {out_db} is already locked by active run {j_id}.")

            cancel_event = threading.Event()
            telemetry = {
                "run_id": run_id,
                "status": "running",
                "processed": run.get("processed_targets", 0),
                "total": run.get("total_targets", 0),
                "successful": run.get("completed_obs", 0),
                "failed": run.get("failed_obs", 0),
                "candidate_count": run.get("candidate_count", 0),
                "rate_per_hr": run.get("throughput_obs_hr", 0.0),
                "elapsed_sec": run.get("elapsed_sec", 0.0),
                "eta_sec": None,
                "peak_ram_mb": 0.0,
                "latest_tic": None,
                "latest_error": None
            }

            th = threading.Thread(
                target=self._execute_run_thread,
                args=(exp, run, max_workers, rate_limit_delay_sec, retain_fits, resume, cancel_event),
                name=f"SurveyWorker-{run_id}",
                daemon=True
            )

            self._active_jobs[run_id] = {
                "thread": th,
                "cancel_event": cancel_event,
                "telemetry": telemetry,
                "db_path": out_db
            }

            self.store.update_run(
                run_id,
                status="running",
                started_utc=utc_now(),
                error_message=None
            )
            self.store.add_run_log(run_id, f"Survey run launched with {max_workers} worker threads.")

            th.start()

        return telemetry

    def _execute_run_thread(
        self,
        exp: Dict[str, Any],
        run: Dict[str, Any],
        max_workers: int,
        rate_limit_delay_sec: float,
        retain_fits: bool,
        resume: bool,
        cancel_event: threading.Event
    ):
        run_id = run["run_id"]
        out_db = run["output_db_path"]
        targets_list = exp.get("targets", [])
        sectors_list = exp.get("sectors", [])

        # Build target cohort
        cohort: List[Dict[str, Any]] = []
        if sectors_list:
            for t in targets_list:
                for s in sectors_list:
                    cohort.append({"tic": sanitize_tic(t), "sector": int(s)})
        else:
            for t in targets_list:
                cohort.append({"tic": sanitize_tic(t)})

        total_cohort = len(cohort)
        self.store.update_run(run_id, total_targets=total_cohort)

        t_start = time.perf_counter()

        def on_progress(pinfo: Dict[str, Any]):
            with self._lock:
                job = self._active_jobs.get(run_id)
                if job:
                    job["telemetry"].update({
                        "processed": pinfo["processed"],
                        "total": pinfo["total"],
                        "successful": pinfo["successful"],
                        "failed": pinfo["failed"],
                        "rate_per_hr": round(pinfo["rate_per_hr"], 1),
                        "elapsed_sec": round(pinfo["elapsed_sec"], 1),
                        "eta_sec": round(pinfo["eta_sec"], 1) if pinfo.get("eta_sec") else None,
                        "peak_ram_mb": round(pinfo["peak_ram_mb"], 1),
                        "latest_tic": pinfo.get("latest_tic"),
                        "latest_error": pinfo.get("error")
                    })

            # Check candidates count in output DB
            c_cnt = 0
            v_cnt = 0
            try:
                with sqlite3.connect(f"file:{Path(out_db).as_posix()}?mode=ro", uri=True) as cx:
                    c_cnt = cx.execute("SELECT COUNT(*) FROM candidates").fetchone()[0]
                    v_cnt = cx.execute(
                        "SELECT COUNT(*) FROM candidates WHERE scientific_classification='verified_localized_dimming'"
                    ).fetchone()[0]
            except Exception:
                pass

            self.store.update_run(
                run_id,
                processed_targets=pinfo["processed"],
                completed_obs=pinfo["successful"],
                failed_obs=pinfo["failed"],
                candidate_count=c_cnt,
                verified_count=v_cnt,
                throughput_obs_hr=pinfo["rate_per_hr"],
                elapsed_sec=pinfo["elapsed_sec"],
                eta_sec=pinfo.get("eta_sec"),
                peak_ram_mb=pinfo["peak_ram_mb"]
            )

        try:
            threshold_t = float(exp.get("threshold_t", FROZEN_THRESHOLD_T))
            method = "fixed_box_2h"

            self.store.add_run_log(run_id, f"Processing {len(cohort)} observations against database: {Path(out_db).name}")

            # Run streaming pipeline
            metrics = streaming.stream_pipeline(
                targets=cohort,
                db_path=out_db,
                run_id=run_id,
                max_workers=max_workers,
                rate_limit_delay_sec=rate_limit_delay_sec,
                retain_fits=retain_fits,
                resume=resume,
                check_dev_universe=False, # General user experiments are unconstrained by dev guard
                threshold_t=threshold_t,
                method=method,
                progress_callback=on_progress,
                cancel_event=cancel_event,
                notes=f"Experiment: {exp['name']}"
            )

            # Count final candidates
            final_cands = 0
            final_verif = 0
            try:
                with sqlite3.connect(f"file:{Path(out_db).as_posix()}?mode=ro", uri=True) as cx:
                    final_cands = cx.execute("SELECT COUNT(*) FROM candidates").fetchone()[0]
                    final_verif = cx.execute(
                        "SELECT COUNT(*) FROM candidates WHERE scientific_classification='verified_localized_dimming'"
                    ).fetchone()[0]
            except Exception:
                pass

            final_status = "aborted" if cancel_event.is_set() else ("completed_with_errors" if metrics['n_failed'] else "completed")
            completed_count = 0
            with sqlite3.connect(Path(out_db).as_uri() + '?mode=ro', uri=True) as cx:
                completed_count = cx.execute("SELECT COUNT(*) FROM observations WHERE run_id=? AND status != 'failed'", (run_id,)).fetchone()[0]
            tot_elapsed = time.perf_counter() - t_start

            self.store.update_run(
                run_id,
                status=final_status,
                completed_obs=completed_count,
                processed_targets=completed_count + metrics['n_failed'],
                failed_obs=metrics['n_failed'],
                ended_utc=utc_now(),
                candidate_count=final_cands,
                verified_count=final_verif,
                elapsed_sec=tot_elapsed
            )

            with self._lock:
                job = self._active_jobs.get(run_id)
                if job:
                    job["telemetry"]["status"] = final_status
                    job["telemetry"]["candidate_count"] = final_cands

            msg = f"Survey {final_status}: {metrics.get('n_successful', 0)} obs completed, {final_cands} candidates found."
            self.store.add_run_log(run_id, msg, "INFO" if final_status == "completed" else "WARN")

        except Exception as e:
            tot_elapsed = time.perf_counter() - t_start
            err_msg = str(e)
            self.store.update_run(
                run_id,
                status="failed",
                ended_utc=utc_now(),
                error_message=err_msg,
                elapsed_sec=tot_elapsed
            )
            with self._lock:
                job = self._active_jobs.get(run_id)
                if job:
                    job["telemetry"]["status"] = "failed"
                    job["telemetry"]["latest_error"] = err_msg

            self.store.add_run_log(run_id, f"Execution failed: {err_msg}", "ERROR")
