"""
Comprehensive automated tests for Experiment Management and Research Workspace.
Tests cover:
- ExperimentStore metadata persistence and historical database auto-discovery
- Pre-flight target validation and resource estimation
- Background SurveyExecutionManager lifecycle, progress tracking, and abort
- Exporting catalogs with full detector configuration provenance and scientific caveats
- Non-destructive overlay reviews preserving bit-for-bit SHA-256 source database integrity
- Web API endpoints for experiments, runs, telemetry, and localhost security guards
"""

import hashlib
import json
import sqlite3
import threading
import time
import urllib.request
import urllib.error
from pathlib import Path

import pytest

from tess_anomaly_survey.experiments import (
    ExperimentStore,
    SurveyExecutionManager,
    validate_and_estimate_targets,
    FROZEN_THRESHOLD_T,
)
from tess_anomaly_survey.web import (
    SurveyDatabaseAdapter,
    AnomalyReviewHandler,
)


def _compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def test_experiment_store_discovery_and_crud(tmp_path):
    exp_dir = tmp_path / "experiments"
    db_dir = tmp_path / "databases"
    db_dir.mkdir()

    # Create dummy historical sqlite DB
    hist_db = db_dir / "development_survey_9753.sqlite"
    with sqlite3.connect(hist_db) as cx:
        cx.execute("CREATE TABLE candidates (id INTEGER PRIMARY KEY, tic_id TEXT, sector INTEGER)")
        cx.execute("INSERT INTO candidates VALUES (1, 'TIC:39018208', 8)")

    store = ExperimentStore(workspace_dir=exp_dir)

    exps = store.list_experiments()
    assert len(exps) >= 1

    # Create custom user experiment
    created = store.create_experiment(
        name="Sector 8 Disk Search",
        description="High SNR search",
        scientific_objective="Inspect circumstellar occultations",
        targets=["39018208", "73834"],
        sectors=[8],
        threshold_t=73.7797,
        duration_hr=2.0
    )
    assert created['id']
    assert created['name'] == "Sector 8 Disk Search"
    assert created['is_historical'] == 0
    assert len(created['runs']) >= 1
    init_run = created['runs'][0]

    # Verify active experiment switching
    store.set_active_selection(created['id'], init_run['run_id'])
    act_exp, act_run = store.get_active_selection()
    assert act_exp == created['id']
    assert act_run == init_run['run_id']

    # Update run status and add log
    store.update_run(init_run['run_id'], status='running', candidate_count=3)
    store.add_run_log(init_run['run_id'], "Processed TIC 39018208", level="INFO")

    updated_run = store.get_run(init_run['run_id'])
    assert updated_run['status'] == 'running'
    assert updated_run['candidate_count'] == 3

    logs = store.get_run_logs(init_run['run_id'])
    assert len(logs) == 1
    assert "Processed TIC 39018208" in logs[0]['message']


def test_validate_and_estimate_targets():
    # Test valid messy input
    text = "TIC 39018208, 73834\nTIC:60232  invalid_abc  12345"
    res = validate_and_estimate_targets(text, sectors=[1, 8])
    assert res['is_safe'] is True
    assert res['total_unique_stars'] == 4
    assert len(res['invalid_tokens']) == 1
    assert res['estimated_observations'] == 8  # 4 * 2 sectors
    assert res['estimated_download_mb'] == round(8 * 2.1, 1)
    assert res['estimated_duration_sec'] > 0

    # Test CSV with header
    csv_text = "tic_id,sector\n111111,1\n222222,8\n333333,8"
    res_csv = validate_and_estimate_targets(csv_text, sectors=[8])
    assert res_csv['is_safe'] is True
    assert res_csv['total_unique_stars'] == 3
    assert res_csv['estimated_observations'] == 3

    # Test empty input error
    res_empty = validate_and_estimate_targets("", sectors=[1])
    assert res_empty['is_safe'] is False
    assert len(res_empty['warnings']) > 0


def test_survey_execution_manager_lifecycle(tmp_path):
    exp_dir = tmp_path / "experiments"
    store = ExperimentStore(workspace_dir=exp_dir)
    manager = SurveyExecutionManager(store)

    exp = store.create_experiment(
        name="Lifecycle Test",
        description="Testing lifecycle",
        scientific_objective="Test manager execution",
        targets=["12345"],
        sectors=[1]
    )
    run_id = exp['runs'][0]['run_id']

    # Test aborting an idle run returns False
    assert manager.abort_survey(run_id) is False

    # Launch survey
    manager.launch_survey(
        experiment_id=exp['id'],
        run_id=run_id,
        max_workers=1,
        resume=True
    )

    assert manager.is_run_active(run_id) is True
    telem = manager.get_live_telemetry(run_id)
    assert telem is not None
    assert telem['status'] in ('running', 'completed', 'failed')

    # Test abort survey and wait for thread
    manager.abort_survey(run_id)
    with manager._lock:
        job = manager._active_jobs.get(run_id)
    if job:
        job["thread"].join(timeout=3.0)

    final_run = store.get_run(run_id)
    assert final_run['status'] in ('aborted', 'completed', 'failed')


def test_export_run_catalog_provenance(tmp_path):
    exp_dir = tmp_path / "experiments"
    store = ExperimentStore(workspace_dir=exp_dir)

    exp = store.create_experiment(
        name="Catalog Provenance Test",
        description="Verify export formats",
        scientific_objective="Verify export formats and headers",
        targets=["39018208"],
        sectors=[8],
        threshold_t=73.7797
    )
    run_id = exp['runs'][0]['run_id']
    run = store.get_run(run_id)
    out_db = Path(run['output_db_path'])

    # Populate candidate table in out_db
    with sqlite3.connect(out_db) as cx:
        cx.execute("""
            CREATE TABLE candidates (
                id INTEGER PRIMARY KEY,
                tic TEXT,
                sector INTEGER,
                center_btjd REAL,
                score REAL,
                width_hr REAL,
                depth_ppm REAL,
                scientific_classification TEXT,
                scientific_reason TEXT,
                review_status TEXT,
                review_note TEXT,
                created_utc TEXT
            )
        """)
        cx.execute("""
            INSERT INTO candidates VALUES (
                101, '39018208', 8, 1517.2, 92.5, 2.0, 1500.0,
                'verified_localized_dimming', 'Coherent box dip',
                'interesting', 'Clear box transit', '2026-10-16T12:00:00Z'
            )
        """)

    # Export CSV
    csv_content, mime_csv = store.export_run_catalog(exp['id'], run_id, fmt='csv')
    assert mime_csv == 'text/csv'
    assert "TESS ANOMALY SURVEY - CANDIDATE CATALOG EXPORT" in csv_content
    assert "threshold_t: 73.7797" in csv_content
    assert "39018208" in csv_content
    assert "verified_localized_dimming" in csv_content

    # Export JSON
    json_content, mime_json = store.export_run_catalog(exp['id'], run_id, fmt='json')
    assert mime_json == 'application/json'
    data = json.loads(json_content)
    assert 'provenance' in data
    assert data['provenance']['threshold_t'] > 73.0
    assert len(data['candidates']) == 1
    assert data['candidates'][0]['tic'] == '39018208'
    assert data['candidates'][0]['human_event_review'] == 'interesting'


def test_historical_database_isolation_and_integrity(tmp_path):
    # Create sample historical database
    source_db = tmp_path / "historical_archive.sqlite"
    with sqlite3.connect(source_db) as cx:
        cx.execute("""
            CREATE TABLE candidates (
                id INTEGER PRIMARY KEY,
                run_id TEXT,
                tic TEXT,
                sector INTEGER,
                source_file TEXT,
                center_btjd REAL,
                score REAL,
                width_hr REAL,
                depth_ppm REAL,
                n_support_windows INTEGER,
                scientific_classification TEXT,
                review_status TEXT DEFAULT 'unreviewed',
                note TEXT DEFAULT ''
            )
        """)
        cx.execute("""
            INSERT INTO candidates VALUES (
                1, 'run_hist', '39018208', 8, 's0008_tic39018208.fits', 1517.2, 92.5, 2.0, 1500.0, 1,
                'verified_localized_dimming', 'unreviewed', ''
            )
        """)

    initial_sha = _compute_sha256(source_db)

    # Mount with SurveyDatabaseAdapter in read-only mode pointing to overlay
    overlay_db = tmp_path / "reviews_overlay.sqlite"
    adapter = SurveyDatabaseAdapter(str(source_db), review_db_path=str(overlay_db), read_only=True)

    # Save candidate and star reviews through the overlay
    adapter.record_event_review(1, 'interesting', 'Validated in overlay')
    adapter.record_star_review('39018208', 'interesting', star_class='candidate_transit', note='Overlay star note')

    # Verify overlay has the reviews reflected in get_star_details
    detail = adapter.get_star_details('39018208')
    assert detail is not None
    assert detail['star_status'] == 'interesting'
    assert detail['star_class'] == 'candidate_transit'
    assert len(detail['events']) == 1
    assert detail['events'][0]['review_status'] == 'interesting'
    assert detail['events'][0]['note'] == 'Validated in overlay'

    # Verify bit-for-bit source database SHA-256 did NOT change
    post_sha = _compute_sha256(source_db)
    assert post_sha == initial_sha, "Source database was mutated during review!"


def test_web_server_experiment_api(tmp_path):
    exp_dir = tmp_path / "experiments"
    db_path = tmp_path / "survey_test.sqlite"

    with sqlite3.connect(db_path) as cx:
        cx.execute("""
            CREATE TABLE candidates (
                id INTEGER PRIMARY KEY,
                tic TEXT,
                sector INTEGER,
                center_btjd REAL,
                score REAL,
                width_hr REAL,
                depth_ppm REAL,
                scientific_classification TEXT,
                scientific_reason TEXT,
                review_status TEXT,
                review_note TEXT,
                created_utc TEXT
            )
        """)
        cx.execute("""
            INSERT INTO candidates VALUES (
                1, '39018208', 8, 1517.2, 92.5, 2.0, 1500.0,
                'verified_localized_dimming', 'Coherent dip',
                'unreviewed', '', '2026-10-16T12:00:00Z'
            )
        """)

    from http.server import HTTPServer

    store = ExperimentStore(workspace_dir=exp_dir)
    exec_mgr = SurveyExecutionManager(store)

    class CustomHandler(AnomalyReviewHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(str(db_path), None, *args, **kwargs)

    server = HTTPServer(('127.0.0.1', 0), CustomHandler)
    server.experiment_store = store
    server.execution_manager = exec_mgr
    server.adapters_cache = {}
    server.initial_db_path = str(db_path)
    port = server.server_port

    server_thread = threading.Thread(target=server.serve_forever, daemon=True)
    server_thread.start()

    base_url = f"http://127.0.0.1:{port}"

    try:
        # 1. GET /api/experiments
        req = urllib.request.Request(f"{base_url}/api/experiments")
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode())
            assert 'experiments' in data
            assert len(data['experiments']) >= 1

        # 2. POST /api/experiments/validate-targets
        val_payload = json.dumps({"targets_text": "39018208, 73834", "sectors": [8]}).encode()
        req = urllib.request.Request(
            f"{base_url}/api/experiments/validate-targets",
            data=val_payload,
            headers={'Content-Type': 'application/json'}
        )
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode())
            assert data['success'] is True
            assert data['total_unique_targets'] == 2

        # 3. POST /api/experiments/create
        create_payload = json.dumps({
            "name": "Integration Test Survey",
            "objective": "API testing",
            "targets_text": "39018208",
            "sectors": [8],
            "launch_immediately": False
        }).encode()
        req = urllib.request.Request(
            f"{base_url}/api/experiments/create",
            data=create_payload,
            headers={'Content-Type': 'application/json'}
        )
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode())
            assert data['success'] is True
            exp_id = data['experiment']['id']
            run_id = data['run']['run_id']

        # 4. POST /api/experiments/active
        switch_payload = json.dumps({"experiment_id": exp_id, "run_id": run_id}).encode()
        req = urllib.request.Request(
            f"{base_url}/api/experiments/active",
            data=switch_payload,
            headers={'Content-Type': 'application/json'}
        )
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode())
            assert data['success'] is True
            assert data['experiment']['id'] == exp_id

        # 5. GET /api/runs/<run_id>/export?format=json
        req = urllib.request.Request(f"{base_url}/api/runs/{run_id}/export?format=json")
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode())
            assert 'provenance' in data
            assert 'candidates' in data

    finally:
        server.shutdown()
        server.server_close()


def test_localhost_security_guard(tmp_path):
    exp_dir = tmp_path / "experiments"
    db_path = tmp_path / "security_test.sqlite"
    with sqlite3.connect(db_path) as cx:
        cx.execute("CREATE TABLE candidates (id INTEGER PRIMARY KEY, tic TEXT)")

    from http.server import HTTPServer

    store = ExperimentStore(workspace_dir=exp_dir)
    exec_mgr = SurveyExecutionManager(store)

    class RemoteSimulatedHandler(AnomalyReviewHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(str(db_path), None, *args, **kwargs)

        def _is_localhost(self) -> bool:
            return False  # Simulate non-localhost remote IP

    server = HTTPServer(('127.0.0.1', 0), RemoteSimulatedHandler)
    server.experiment_store = store
    server.execution_manager = exec_mgr
    server.adapters_cache = {}
    server.initial_db_path = str(db_path)
    port = server.server_port

    server_thread = threading.Thread(target=server.serve_forever, daemon=True)
    server_thread.start()

    base_url = f"http://127.0.0.1:{port}"

    try:
        # State-mutating POST must be rejected with 403 Forbidden
        payload = json.dumps({"name": "Disallowed Remote Experiment"}).encode()
        req = urllib.request.Request(
            f"{base_url}/api/experiments/create",
            data=payload,
            headers={'Content-Type': 'application/json'}
        )
        with pytest.raises(urllib.error.HTTPError) as exc_info:
            urllib.request.urlopen(req)
        assert exc_info.value.code == 403
    finally:
        server.shutdown()
        server.server_close()


def test_host_and_origin_security_guards(tmp_path):
    exp_dir = tmp_path / "experiments"
    db_path = tmp_path / "sec_test.sqlite"
    with sqlite3.connect(db_path) as cx:
        cx.execute("CREATE TABLE candidates (id INTEGER PRIMARY KEY, tic TEXT)")

    from http.server import HTTPServer

    store = ExperimentStore(workspace_dir=exp_dir)
    exec_mgr = SurveyExecutionManager(store)

    class CustomHandler(AnomalyReviewHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(str(db_path), None, *args, **kwargs)

    server = HTTPServer(('127.0.0.1', 0), CustomHandler)
    server.experiment_store = store
    server.execution_manager = exec_mgr
    server.adapters_cache = {}
    server.initial_db_path = str(db_path)
    port = server.server_port

    server_thread = threading.Thread(target=server.serve_forever, daemon=True)
    server_thread.start()

    base_url = f"http://127.0.0.1:{port}"

    try:
        payload = json.dumps({"name": "Test"}).encode()

        # 1. External Origin header -> 403 Forbidden
        req_origin = urllib.request.Request(
            f"{base_url}/api/experiments/create",
            data=payload,
            headers={'Content-Type': 'application/json', 'Origin': 'https://evil-attacker.com'}
        )
        with pytest.raises(urllib.error.HTTPError) as exc:
            urllib.request.urlopen(req_origin)
        assert exc.value.code == 403

        # 2. DNS Rebinding Host header -> 403 Forbidden
        req_host = urllib.request.Request(
            f"{base_url}/api/experiments/create",
            data=payload,
            headers={'Content-Type': 'application/json', 'Host': 'evil-domain.com:8000'}
        )
        with pytest.raises(urllib.error.HTTPError) as exc:
            urllib.request.urlopen(req_host)
        assert exc.value.code == 403
    finally:
        server.shutdown()
        server.server_close()

