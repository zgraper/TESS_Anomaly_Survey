from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sqlite3

SCHEMA = '''
PRAGMA foreign_keys=ON;
CREATE TABLE IF NOT EXISTS runs (
 id TEXT PRIMARY KEY, created_utc TEXT NOT NULL, calibration_sha256 TEXT NOT NULL,
 method TEXT NOT NULL, target_alert_fraction REAL NOT NULL,
 settings_json TEXT NOT NULL, catalog_path TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS processed_files (
 run_id TEXT NOT NULL, source_file TEXT NOT NULL, tic TEXT NOT NULL,
 sector INTEGER, n_windows INTEGER NOT NULL, n_alert_windows INTEGER NOT NULL,
 status TEXT NOT NULL, error TEXT,
 PRIMARY KEY(run_id, source_file), FOREIGN KEY(run_id) REFERENCES runs(id)
);
CREATE TABLE IF NOT EXISTS windows (
 run_id TEXT NOT NULL, source_file TEXT NOT NULL, tic TEXT NOT NULL,
 start_index INTEGER NOT NULL, start_btjd REAL NOT NULL,
 noise REAL NOT NULL, max_score REAL NOT NULL, cutoff REAL NOT NULL,
 noise_bin INTEGER NOT NULL, flagged INTEGER NOT NULL,
 PRIMARY KEY(run_id, source_file, start_index),
 FOREIGN KEY(run_id) REFERENCES runs(id)
);
CREATE TABLE IF NOT EXISTS candidates (
 id INTEGER PRIMARY KEY AUTOINCREMENT,
 run_id TEXT NOT NULL, tic TEXT NOT NULL, sector INTEGER,
 source_file TEXT NOT NULL, center_btjd REAL NOT NULL,
 score REAL NOT NULL, width_hr REAL NOT NULL, depth_ppm REAL NOT NULL,
 n_support_windows INTEGER NOT NULL DEFAULT 1,
 review_status TEXT NOT NULL DEFAULT 'unreviewed', note TEXT NOT NULL DEFAULT '',
 system_id INTEGER,
 FOREIGN KEY(run_id) REFERENCES runs(id),
 FOREIGN KEY(system_id) REFERENCES candidate_systems(id)
);
CREATE INDEX IF NOT EXISTS idx_candidates_rank ON candidates(run_id,score DESC);
CREATE INDEX IF NOT EXISTS idx_candidates_tic_center ON candidates(run_id,tic,center_btjd);

CREATE TABLE IF NOT EXISTS candidate_systems (
 id INTEGER PRIMARY KEY AUTOINCREMENT,
 run_id TEXT NOT NULL,
 tic TEXT NOT NULL,
 n_events INTEGER NOT NULL,
 period_days REAL,
 epoch_btjd REAL,
 period_snr REAL,
 is_periodic INTEGER NOT NULL DEFAULT 0,
 odd_even_sigma REAL,
 classification TEXT NOT NULL DEFAULT 'unreviewed',
 note TEXT NOT NULL DEFAULT '',
 period_uncertainty_days REAL,
 epoch_uncertainty_days REAL,
 alternate_period_days REAL,
 competing_solutions_json TEXT,
 created_utc TEXT NOT NULL,
 FOREIGN KEY(run_id) REFERENCES runs(id)
);
CREATE INDEX IF NOT EXISTS idx_systems_tic ON candidate_systems(run_id, tic);

CREATE TABLE IF NOT EXISTS candidate_detections (
 id INTEGER PRIMARY KEY AUTOINCREMENT,
 candidate_id INTEGER NOT NULL,
 run_id TEXT NOT NULL,
 source_file TEXT NOT NULL,
 sector INTEGER,
 center_btjd REAL NOT NULL,
 score REAL NOT NULL,
 width_hr REAL NOT NULL,
 depth_ppm REAL NOT NULL,
 created_utc TEXT NOT NULL,
 FOREIGN KEY(candidate_id) REFERENCES candidates(id) ON DELETE CASCADE,
 FOREIGN KEY(run_id) REFERENCES runs(id)
);
CREATE INDEX IF NOT EXISTS idx_detections_candidate ON candidate_detections(candidate_id);
CREATE INDEX IF NOT EXISTS idx_detections_run_tic ON candidate_detections(run_id, source_file);

CREATE TABLE IF NOT EXISTS review_history (
 id INTEGER PRIMARY KEY AUTOINCREMENT,
 candidate_id INTEGER NOT NULL,
 previous_status TEXT NOT NULL,
 new_status TEXT NOT NULL,
 note TEXT NOT NULL,
 reviewer TEXT NOT NULL DEFAULT 'local_user',
 created_utc TEXT NOT NULL,
 FOREIGN KEY(candidate_id) REFERENCES candidates(id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_review_history_candidate ON review_history(candidate_id);

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
CREATE INDEX IF NOT EXISTS idx_star_history_tic ON star_review_history(tic);
'''

def backup_database(db_path: Path) -> Path:
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    bak = db_path.with_name(f"{db_path.name}.bak_{ts}")
    source_con = sqlite3.connect(str(db_path))
    try:
        dest_con = sqlite3.connect(str(bak))
        try:
            source_con.backup(dest_con)
        finally:
            dest_con.close()
    finally:
        source_con.close()
    return bak


def connect(db_path):
    p = Path(db_path)
    is_existing = str(db_path) != ":memory:" and p.exists() and p.stat().st_size > 0
    if is_existing:
        probe = sqlite3.connect(str(p))
        try:
            has_det = probe.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='candidate_detections'"
            ).fetchone()
            has_sys = probe.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='candidate_systems'"
            ).fetchone()
        finally:
            probe.close()
        if not has_det or not has_sys:
            backup_database(p)

    p.parent.mkdir(parents=True, exist_ok=True)
    cx = sqlite3.connect(str(p))
    with cx:
        cx.executescript(SCHEMA)
        # Ensure system_id column exists on candidates table for older databases
        col_names = [r[1] for r in cx.execute("PRAGMA table_info(candidates)").fetchall()]
        if col_names and "system_id" not in col_names:
            cx.execute("ALTER TABLE candidates ADD COLUMN system_id INTEGER REFERENCES candidate_systems(id)")

        # Ensure candidate_systems has Phase 3.1 uncertainty and competing solution columns
        sys_cols = [r[1] for r in cx.execute("PRAGMA table_info(candidate_systems)").fetchall()]
        for col, col_type in [
            ("period_uncertainty_days", "REAL"),
            ("epoch_uncertainty_days", "REAL"),
            ("alternate_period_days", "REAL"),
            ("competing_solutions_json", "TEXT"),
        ]:
            if sys_cols and col not in sys_cols:
                cx.execute(f"ALTER TABLE candidate_systems ADD COLUMN {col} {col_type}")

        # Seed legacy candidates into candidate_detections so historical provenance is preserved
        cx.execute('''INSERT OR IGNORE INTO candidate_detections(
            candidate_id, run_id, source_file, sector, center_btjd, score, width_hr, depth_ppm, created_utc)
            SELECT id, run_id, source_file, sector, center_btjd, score, width_hr, depth_ppm, datetime('now')
            FROM candidates
            WHERE id NOT IN (SELECT candidate_id FROM candidate_detections)''')
    cx.execute('PRAGMA busy_timeout=30000')
    return cx


def record_review(cx, candidate_id, new_status, note="", reviewer="local_user"):
    valid_statuses = {'unreviewed', 'interesting', 'known', 'artifact', 'dismissed'}
    if new_status not in valid_statuses:
        raise ValueError(f"Invalid status: {new_status}. Must be one of {valid_statuses}")
    now_utc = datetime.now(timezone.utc).isoformat()
    with cx:
        row = cx.execute('SELECT review_status, note FROM candidates WHERE id=?', (candidate_id,)).fetchone()
        if not row:
            raise ValueError(f"Candidate ID {candidate_id} not found")
        prev_status, prev_note = row
        cx.execute('UPDATE candidates SET review_status=?, note=? WHERE id=?',
                   (new_status, note, candidate_id))
        cx.execute('''INSERT INTO review_history(candidate_id, previous_status, new_status, note, reviewer, created_utc)
            VALUES (?,?,?,?,?,?)''',
            (candidate_id, prev_status, new_status, note, reviewer, now_utc))


def get_review_history(cx, candidate_id):
    cur = cx.execute('''SELECT id, previous_status, new_status, note, reviewer, created_utc
        FROM review_history WHERE candidate_id=? ORDER BY id DESC''', (candidate_id,))
    return [{'id': r[0], 'previous_status': r[1], 'new_status': r[2], 'note': r[3],
             'reviewer': r[4], 'created_utc': r[5]} for r in cur.fetchall()]


def record_star_review(cx, tic, new_status, note="", reviewer="local_user", star_class=""):
    valid_statuses = {'unreviewed', 'interesting', 'known', 'artifact', 'dismissed'}
    if new_status not in valid_statuses:
        raise ValueError(f"Invalid status: {new_status}. Must be one of {valid_statuses}")
    now_utc = datetime.now(timezone.utc).isoformat()
    with cx:
        row = cx.execute('SELECT status, note, star_class FROM star_reviews WHERE tic=?', (tic,)).fetchone()
        prev_status = row[0] if row else 'unreviewed'
        cx.execute('''INSERT INTO star_reviews(tic, status, star_class, note, reviewer, updated_utc)
            VALUES (?,?,?,?,?,?)
            ON CONFLICT(tic) DO UPDATE SET status=excluded.status, star_class=excluded.star_class,
            note=excluded.note, reviewer=excluded.reviewer, updated_utc=excluded.updated_utc''',
            (tic, new_status, star_class, note, reviewer, now_utc))
        cx.execute('''INSERT INTO star_review_history(tic, previous_status, new_status, note, reviewer, created_utc)
            VALUES (?,?,?,?,?,?)''',
            (tic, prev_status, new_status, note, reviewer, now_utc))


def undo_review(cx, candidate_id: int):
    """Revert the most recent review action on an individual candidate event."""
    with cx:
        latest = cx.execute('''SELECT id, previous_status, note FROM review_history
            WHERE candidate_id=? ORDER BY id DESC LIMIT 1''', (candidate_id,)).fetchone()
        if not latest:
            raise ValueError(f"No review history to undo for candidate {candidate_id}")
        hist_id, prev_status, prev_note = latest
        cx.execute('UPDATE candidates SET review_status=? WHERE id=?', (prev_status, candidate_id))
        cx.execute('DELETE FROM review_history WHERE id=?', (hist_id,))
        return {'candidate_id': candidate_id, 'reverted_to': prev_status}


def undo_star_review(cx, tic: str):
    """Revert the most recent star-level review action."""
    with cx:
        latest = cx.execute('''SELECT id, previous_status, note FROM star_review_history
            WHERE tic=? ORDER BY id DESC LIMIT 1''', (tic,)).fetchone()
        if not latest:
            raise ValueError(f"No star review history to undo for star {tic}")
        hist_id, prev_status, prev_note = latest
        cx.execute('UPDATE star_reviews SET status=? WHERE tic=?', (prev_status, tic))
        cx.execute('DELETE FROM star_review_history WHERE id=?', (hist_id,))
        return {'tic': tic, 'reverted_to': prev_status}


def bulk_review_star_events(cx, tic, new_status, note="", reviewer="local_user", only_unreviewed=True):
    """Classify remaining candidate events for a star, maintaining full event audit provenance."""
    valid_statuses = {'unreviewed', 'interesting', 'known', 'artifact', 'dismissed'}
    if new_status not in valid_statuses:
        raise ValueError(f"Invalid status: {new_status}. Must be one of {valid_statuses}")

    cond = "WHERE tic=? AND review_status='unreviewed'" if only_unreviewed else "WHERE tic=?"
    cand_ids = [r[0] for r in cx.execute(f"SELECT id FROM candidates {cond}", (tic,)).fetchall()]
    with cx:
        for cid in cand_ids:
            record_review(cx, cid, new_status, note=note, reviewer=reviewer)
        record_star_review(cx, tic, new_status, note=note, reviewer=reviewer)
    return len(cand_ids)


def save_candidate_system(cx, run_id, tic, result, candidate_ids):
    """Store or update periodicity analysis result for a star and link candidate IDs."""
    now_utc = datetime.now(timezone.utc).isoformat()
    competing_json = json.dumps(getattr(result, 'competing_solutions', []))
    classification = getattr(result, 'classification', 'unreviewed')
    p_err = getattr(result, 'period_uncertainty_days', None)
    t0_err = getattr(result, 'epoch_uncertainty_days', None)
    alt_p = getattr(result, 'alternate_period_days', None)

    with cx:
        existing = cx.execute('SELECT id FROM candidate_systems WHERE run_id=? AND tic=?', (run_id, tic)).fetchone()
        if existing:
            sys_id = existing[0]
            cx.execute('''UPDATE candidate_systems SET n_events=?, period_days=?, epoch_btjd=?,
                period_snr=?, is_periodic=?, odd_even_sigma=?, classification=?, note=?,
                period_uncertainty_days=?, epoch_uncertainty_days=?, alternate_period_days=?, competing_solutions_json=?,
                created_utc=?
                WHERE id=?''',
                (result.n_events, result.period_days, result.epoch_btjd, result.period_snr,
                 int(result.is_periodic), result.odd_even_sigma, classification,
                 result.note, p_err, t0_err, alt_p, competing_json, now_utc, sys_id))
        else:
            cur = cx.execute('''INSERT INTO candidate_systems(
                run_id, tic, n_events, period_days, epoch_btjd, period_snr, is_periodic,
                odd_even_sigma, classification, note, period_uncertainty_days, epoch_uncertainty_days,
                alternate_period_days, competing_solutions_json, created_utc)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',
                (run_id, tic, result.n_events, result.period_days, result.epoch_btjd,
                 result.period_snr, int(result.is_periodic), result.odd_even_sigma,
                 classification, result.note, p_err, t0_err, alt_p, competing_json, now_utc))
            sys_id = cur.lastrowid

        if candidate_ids:
            placeholders = ','.join('?' for _ in candidate_ids)
            cx.execute(f'UPDATE candidates SET system_id=? WHERE id IN ({placeholders})',
                       [sys_id] + list(candidate_ids))
    return sys_id


def get_system(cx, system_id):
    """Retrieve candidate system information by system_id."""
    row = cx.execute('''SELECT id, run_id, tic, n_events, period_days, epoch_btjd,
        period_snr, is_periodic, odd_even_sigma, classification, note,
        period_uncertainty_days, epoch_uncertainty_days, alternate_period_days,
        competing_solutions_json, created_utc
        FROM candidate_systems WHERE id=?''', (system_id,)).fetchone()
    if not row:
        return None
    return {
        'id': row[0], 'run_id': row[1], 'tic': row[2], 'n_events': row[3],
        'period_days': row[4], 'epoch_btjd': row[5], 'period_snr': row[6],
        'is_periodic': bool(row[7]), 'odd_even_sigma': row[8],
        'classification': row[9], 'note': row[10],
        'period_uncertainty_days': row[11], 'epoch_uncertainty_days': row[12],
        'alternate_period_days': row[13],
        'competing_solutions': json.loads(row[14]) if row[14] else [],
        'created_utc': row[15]
    }


def get_systems_for_run(cx, run_id):
    """Retrieve all candidate systems detected for a given run."""
    cur = cx.execute('''SELECT id, run_id, tic, n_events, period_days, epoch_btjd,
        period_snr, is_periodic, odd_even_sigma, classification, note,
        period_uncertainty_days, epoch_uncertainty_days, alternate_period_days,
        competing_solutions_json, created_utc
        FROM candidate_systems WHERE run_id=? ORDER BY period_snr DESC''', (run_id,))
    return [{
        'id': row[0], 'run_id': row[1], 'tic': row[2], 'n_events': row[3],
        'period_days': row[4], 'epoch_btjd': row[5], 'period_snr': row[6],
        'is_periodic': bool(row[7]), 'odd_even_sigma': row[8],
        'classification': row[9], 'note': row[10],
        'period_uncertainty_days': row[11], 'epoch_uncertainty_days': row[12],
        'alternate_period_days': row[13],
        'competing_solutions': json.loads(row[14]) if row[14] else [],
        'created_utc': row[15]
    } for row in cur.fetchall()]


def get_star_summary(cx, tic):
    cand_rows = cx.execute('''SELECT id, center_btjd, score, width_hr, depth_ppm, review_status, note, system_id
        FROM candidates WHERE tic=? ORDER BY score DESC''', (tic,)).fetchall()
    total = len(cand_rows)
    unreviewed = sum(1 for r in cand_rows if r[5] == 'unreviewed')
    reviewed = total - unreviewed
    if total == 0:
        progress = 'unreviewed'
    elif unreviewed == 0:
        progress = 'reviewed'
    elif reviewed > 0:
        progress = 'partially_reviewed'
    else:
        progress = 'unreviewed'

    star_row = cx.execute('SELECT status, star_class, note, reviewer, updated_utc FROM star_reviews WHERE tic=?', (tic,)).fetchone()
    history_cur = cx.execute('''SELECT id, previous_status, new_status, note, reviewer, created_utc
        FROM star_review_history WHERE tic=? ORDER BY id DESC''', (tic,))

    sys_row = cx.execute('''SELECT id, n_events, period_days, epoch_btjd, period_snr,
        is_periodic, odd_even_sigma, classification, note,
        period_uncertainty_days, epoch_uncertainty_days, alternate_period_days,
        competing_solutions_json, created_utc
        FROM candidate_systems WHERE tic=? ORDER BY id DESC LIMIT 1''', (tic,)).fetchone()
    system_info = None
    if sys_row:
        system_info = {
            'id': sys_row[0], 'n_events': sys_row[1], 'period_days': sys_row[2],
            'epoch_btjd': sys_row[3], 'period_snr': sys_row[4], 'is_periodic': bool(sys_row[5]),
            'odd_even_sigma': sys_row[6], 'classification': sys_row[7], 'note': sys_row[8],
            'period_uncertainty_days': sys_row[9], 'epoch_uncertainty_days': sys_row[10],
            'alternate_period_days': sys_row[11],
            'competing_solutions': json.loads(sys_row[12]) if sys_row[12] else [],
            'created_utc': sys_row[13]
        }

    return {
        'tic': tic,
        'total_events': total,
        'reviewed_events': reviewed,
        'unreviewed_events': unreviewed,
        'progress': progress,
        'star_status': star_row[0] if star_row else 'unreviewed',
        'star_class': star_row[1] if star_row else '',
        'note': star_row[2] if star_row else '',
        'reviewer': star_row[3] if star_row else '',
        'updated_utc': star_row[4] if star_row else '',
        'system': system_info,
        'history': [{'id': r[0], 'previous_status': r[1], 'new_status': r[2], 'note': r[3],
                     'reviewer': r[4], 'created_utc': r[5]} for r in history_cur.fetchall()],
        'events': [{'id': r[0], 'center_btjd': r[1], 'score': r[2], 'width_hr': r[3],
                    'depth_ppm': r[4], 'review_status': r[5], 'note': r[6],
                    'system_id': r[7]} for r in cand_rows]
    }


def get_survey_review_stats(cx):
    event_counts = {}
    for s in ('unreviewed', 'interesting', 'known', 'artifact', 'dismissed'):
        event_counts[s] = cx.execute('SELECT COUNT(*) FROM candidates WHERE review_status=?', (s,)).fetchone()[0]
    total_events = cx.execute('SELECT COUNT(*) FROM candidates').fetchone()[0]
    reviewed_events = total_events - event_counts['unreviewed']

    # Star-level statistics
    per_star = cx.execute('''SELECT tic,
        SUM(CASE WHEN review_status='unreviewed' THEN 1 ELSE 0 END) as unrev,
        COUNT(*) as total
        FROM candidates GROUP BY tic''').fetchall()

    total_stars = len(per_star)
    fully_reviewed_stars = sum(1 for r in per_star if r[1] == 0 and r[2] > 0)
    partially_reviewed_stars = sum(1 for r in per_star if 0 < r[1] < r[2])
    unreviewed_stars = sum(1 for r in per_star if r[1] == r[2] and r[2] > 0)

    # Classification counts from star_reviews
    star_status_counts = {}
    for s in ('unreviewed', 'interesting', 'known', 'artifact', 'dismissed'):
        star_status_counts[s] = cx.execute('SELECT COUNT(*) FROM star_reviews WHERE status=?', (s,)).fetchone()[0]

    # Periodic systems count
    has_systems = cx.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='candidate_systems'"
    ).fetchone()
    periodic_count = cx.execute(
        "SELECT COUNT(*) FROM candidate_systems WHERE is_periodic=1"
    ).fetchone()[0] if has_systems else 0

    return {
        'events': {
            'total': total_events,
            'reviewed': reviewed_events,
            'unreviewed': event_counts['unreviewed'],
            'counts': event_counts
        },
        'stars': {
            'total': total_stars,
            'fully_reviewed': fully_reviewed_stars,
            'partially_reviewed': partially_reviewed_stars,
            'unreviewed': unreviewed_stars,
            'classification_counts': star_status_counts,
            'periodic_systems': periodic_count
        }
    }


def begin_run(cx,run_id,cal_hash,cal,settings,catalog_path,created_utc):
    row=cx.execute('SELECT calibration_sha256,settings_json,catalog_path FROM runs WHERE id=?',(run_id,)).fetchone()
    serialized=json.dumps(settings,sort_keys=True)
    if row:
        if row!=(cal_hash,serialized,str(catalog_path)):
            raise RuntimeError('A run ID already exists for a different configuration')
    else:
        cx.execute('INSERT INTO runs VALUES (?,?,?,?,?,?,?)',(
            run_id,created_utc,cal_hash,cal['method'],cal['target_alert_fraction'],serialized,str(catalog_path)))
        cx.commit()


def is_done(cx,run_id,source_file):
    return cx.execute('SELECT 1 FROM processed_files WHERE run_id=? AND source_file=?',(
        run_id,source_file)).fetchone() is not None


def upsert_candidate(cx, run_id, row, source_file, cand, cadence):
    """Merge near-coincident detections from overlapping windows into a unique event,
    preserving all individual detection records and window provenance.
    Handles repeat scans and order-independent bridging of chained detections.
    """
    tic = row['tic']
    center = float(cand['center_btjd'])
    width = float(cand['width_hr'])
    score = float(cand['score'])
    depth = float(cand['depth_ppm'])
    sector = int(row['sector'])
    now_utc = datetime.now(timezone.utc).isoformat()

    # Search horizon covers the maximum template duration (up to 8h)
    search_radius = max(3 * cadence, 0.5 * max(8.0, width) / 24.0)
    neighbors = cx.execute('''SELECT id, center_btjd, width_hr, score, n_support_windows
        FROM candidates WHERE run_id=? AND tic=? AND center_btjd BETWEEN ? AND ?''',
        (run_id, tic, center - search_radius, center + search_radius)).fetchall()

    # Match if center separation is within half the maximum duration
    matching = [r for r in neighbors if abs(float(r[1]) - center) <=
                max(3 * cadence, 0.5 * max(float(r[2]), width) / 24.0)]

    if matching:
        nearest = min(matching, key=lambda x: abs(x[1] - center))
        candidate_id = nearest[0]

        # Check for repeat scan of exact same window detection
        already_present = cx.execute('''SELECT 1 FROM candidate_detections
            WHERE candidate_id=? AND run_id=? AND source_file=?
            AND ABS(center_btjd - ?) <= 1e-6 AND ABS(width_hr - ?) <= 1e-6''',
            (candidate_id, run_id, source_file, center, width)).fetchone()

        if already_present:
            # Repeat scan / re-processing: do not duplicate detection or increment windows
            return

        # If this detection bridges multiple existing candidate clusters (order-independence):
        if len(matching) > 1:
            for other in matching:
                other_id = other[0]
                if other_id == candidate_id:
                    continue
                # Re-parent all detections from the other candidate
                cx.execute('UPDATE candidate_detections SET candidate_id=? WHERE candidate_id=?',
                           (candidate_id, other_id))
                # Accumulate support window count
                cx.execute('UPDATE candidates SET n_support_windows=n_support_windows + ? WHERE id=?',
                           (other[4], candidate_id))
                # If other candidate had a higher score, update peak info
                if float(other[3]) > score:
                    score = float(other[3])
                    center = float(other[1])
                    width = float(other[2])
                    depth = float(other[3]) # score was compared, depth corresponds to peak
                # Remove redundant candidate row
                cx.execute('DELETE FROM candidates WHERE id=?', (other_id,))

        # Record this individual detection
        cx.execute('''INSERT INTO candidate_detections(candidate_id, run_id, source_file, sector,
            center_btjd, score, width_hr, depth_ppm, created_utc)
            VALUES (?,?,?,?,?,?,?,?,?)''',
            (candidate_id, run_id, source_file, sector, float(cand['center_btjd']),
             float(cand['score']), float(cand['width_hr']), depth, now_utc))

        if score > float(nearest[3]):
            cx.execute('''UPDATE candidates SET center_btjd=?, score=?, width_hr=?, depth_ppm=?,
                source_file=?, sector=?, n_support_windows=n_support_windows+1 WHERE id=?''',
                (center, score, width, depth, source_file, sector, candidate_id))
        else:
            cx.execute('UPDATE candidates SET n_support_windows=n_support_windows+1 WHERE id=?',
                       (candidate_id,))
    else:
        cur = cx.execute('''INSERT INTO candidates(run_id, tic, sector, source_file, center_btjd,
            score, width_hr, depth_ppm, n_support_windows) VALUES (?,?,?,?,?,?,?,?,1)''',
            (run_id, tic, sector, source_file, center, score, width, depth))
        candidate_id = cur.lastrowid
        cx.execute('''INSERT INTO candidate_detections(candidate_id, run_id, source_file, sector,
            center_btjd, score, width_hr, depth_ppm, created_utc)
            VALUES (?,?,?,?,?,?,?,?,?)''',
            (candidate_id, run_id, source_file, sector, center, score, width, depth, now_utc))



def statistics(cx):
    runs=[]
    for run in cx.execute('SELECT id,method,target_alert_fraction,created_utc FROM runs ORDER BY created_utc'):
        rid,method,target,when=run
        n,flag=cx.execute('SELECT COUNT(*),COALESCE(SUM(flagged),0) FROM windows WHERE run_id=?',(rid,)).fetchone()
        cand=cx.execute('SELECT COUNT(*) FROM candidates WHERE run_id=?',(rid,)).fetchone()[0]
        files=cx.execute('SELECT COUNT(*) FROM processed_files WHERE run_id=? AND status=?',(rid,'ok')).fetchone()[0]
        bad=cx.execute('SELECT COUNT(*) FROM processed_files WHERE run_id=? AND status=?',(rid,'error')).fetchone()[0]
        runs.append({'run_id':rid,'method':method,'nominal_target':target,'files_processed':files,
                     'files_failed':bad,'windows':n,'flagged_windows':flag,
                     'observed_flagged_fraction':flag/n if n else None,
                     'unique_review_candidates':cand,'created_utc':when})
    return runs
