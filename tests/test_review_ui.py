import json
from http.server import HTTPServer
import threading
import time
import urllib.request
import urllib.error
from pathlib import Path
import pytest
import sqlite3

from tess_anomaly_survey import db
from tess_anomaly_survey.web import AnomalyReviewHandler, render_plot_png, SurveyDatabaseAdapter


@pytest.mark.parametrize('path', ['/', '/api/stats', '/api/stars'])
@pytest.mark.parametrize('headers', [
    {'Host': 'evil.example'}, {'Host': 'localhost:1'},
    {'Origin': 'http://localhost:1'}, {'Origin': 'null'},
    {'Sec-Fetch-Site': 'cross-site'},
])
def test_read_authority_security(review_server, path, headers):
    with pytest.raises(urllib.error.HTTPError) as error:
        urllib.request.urlopen(urllib.request.Request(review_server + path, headers=headers))
    assert error.value.code == 403


def test_browser_write_requires_csrf_token(review_server):
    payload = json.dumps({'candidate_id': 1, 'status': 'interesting'}).encode()
    request = urllib.request.Request(review_server + '/api/review', data=payload,
        headers={'Content-Type': 'application/json', 'Origin': review_server})
    with pytest.raises(urllib.error.HTTPError) as error:
        urllib.request.urlopen(request)
    assert error.value.code == 403

@pytest.fixture
def test_env(tmp_path):
    db_file = tmp_path / "review_test.sqlite"
    con = db.connect(db_file)
    run_id = "test_run"
    con.execute("INSERT INTO runs VALUES (?,?,?,?,?,?,?)",
                (run_id, "2026-10-13T00:00:00Z", "hash", "fixed_box_2h", 0.01, "{}", "cat.csv"))

    # Candidate 1: Star 1001
    con.execute("""INSERT INTO candidates(id, run_id, tic, sector, source_file, center_btjd,
                   score, width_hr, depth_ppm, n_support_windows, review_status, note)
                   VALUES (1, ?, 'TIC:1001', 1, 's1.npz', 1000.5, 12.5, 2.0, 1500.0, 2, 'unreviewed', '')""",
                (run_id,))
    con.execute("""INSERT INTO candidate_detections(candidate_id, run_id, source_file, sector,
                   center_btjd, score, width_hr, depth_ppm, created_utc)
                   VALUES (1, ?, 's1.npz', 1, 1000.5, 12.5, 2.0, 1500.0, '2026-10-13T00:00:00Z')""",
                (run_id,))

    # Candidate 2: Star 1001 (second candidate for same star)
    con.execute("""INSERT INTO candidates(id, run_id, tic, sector, source_file, center_btjd,
                   score, width_hr, depth_ppm, n_support_windows, review_status, note)
                   VALUES (2, ?, 'TIC:1001', 1, 's1.npz', 1003.8, 8.2, 1.0, 900.0, 1, 'unreviewed', '')""",
                (run_id,))

    # Candidate 3: Star 2002
    con.execute("""INSERT INTO candidates(id, run_id, tic, sector, source_file, center_btjd,
                   score, width_hr, depth_ppm, n_support_windows, review_status, note)
                   VALUES (3, ?, 'TIC:2002', 1, 's2.npz', 2000.2, 15.0, 2.0, 2500.0, 1, 'unreviewed', '')""",
                (run_id,))

    con.commit()
    con.close()
    return db_file, tmp_path


@pytest.fixture
def review_server(test_env):
    db_file, root_dir = test_env
    # Pick port 0 for random available port
    class CustomHandler(AnomalyReviewHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(db_file, root_dir, *args, **kwargs)

    server = HTTPServer(("127.0.0.1", 0), CustomHandler)
    host, port = server.server_address
    t = threading.Thread(target=server.serve_forever, daemon=True)
    t.start()
    base_url = f"http://{host}:{port}"
    yield base_url
    server.shutdown()
    server.server_close()


def test_web_index_and_stats(review_server):
    # Test GET /
    with urllib.request.urlopen(f"{review_server}/") as resp:
        assert resp.status == 200
        content = resp.read().decode("utf-8")
        assert "TESS Anomaly Survey" in content
        assert "Review Workspace" in content

    # Test GET /api/stats (separate event and star metrics)
    with urllib.request.urlopen(f"{review_server}/api/stats") as resp:
        assert resp.status == 200
        stats = json.loads(resp.read().decode("utf-8"))
        ev = stats["events"]
        st = stats["stars"]
        assert ev["total"] == 3
        assert ev["unreviewed"] == 3
        assert ev["reviewed"] == 0
        assert st["total"] == 2
        assert st["unreviewed"] == 2
        assert st["fully_reviewed"] == 0
        assert st["partially_reviewed"] == 0


def test_web_candidates_filtering(review_server):
    # Filter all
    with urllib.request.urlopen(f"{review_server}/api/candidates?status=all") as resp:
        data = json.loads(resp.read().decode("utf-8"))
        assert data["total"] == 3
        assert len(data["candidates"]) == 3

    # Filter by TIC
    with urllib.request.urlopen(f"{review_server}/api/candidates?tic=TIC:2002") as resp:
        data = json.loads(resp.read().decode("utf-8"))
        assert data["total"] == 1
        assert data["candidates"][0]["tic"] == "TIC:2002"

    # Search filter
    with urllib.request.urlopen(f"{review_server}/api/candidates?search=1001") as resp:
        data = json.loads(resp.read().decode("utf-8"))
        assert data["total"] == 2


def test_web_candidate_detail_and_star_summary(review_server):
    # Test candidate 1 detail
    with urllib.request.urlopen(f"{review_server}/api/candidate/1") as resp:
        data = json.loads(resp.read().decode("utf-8"))
        assert data["candidate"]["id"] == 1
        assert data["candidate"]["tic"] == "TIC:1001"
        assert len(data["detections"]) == 1
        # Star summary shows 2 total events
        assert data["star"]["tic"] == "TIC:1001"
        assert data["star"]["total_events"] == 2
        assert data["star"]["progress"] == "unreviewed"


def test_individual_event_review_leaves_sibling_unreviewed(review_server, test_env):
    """Verify that reviewing candidate 1:

    1. Updates candidate 1 to 'interesting'
    2. Removes candidate 1 from status=unreviewed list
    3. Leaves candidate 2 on the SAME star as 'unreviewed'
    4. Updates the star status to 'partially_reviewed'
    """
    db_file, _ = test_env
    # Submit review for Candidate 1
    req_data = json.dumps({
        "candidate_id": 1,
        "status": "interesting",
        "note": "Promising localized dip",
        "reviewer": "test_observer"
    }).encode("utf-8")

    req = urllib.request.Request(f"{review_server}/api/review", data=req_data,
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as resp:
        assert resp.status == 200
        result = json.loads(resp.read().decode("utf-8"))
        assert result["success"] is True
        assert result["status"] == "interesting"
        # Star summary included in response shows partial progress
        assert result["star"]["progress"] == "partially_reviewed"
        assert result["star"]["reviewed_events"] == 1
        assert result["star"]["unreviewed_events"] == 1

    # Verify that Candidate 1 NO LONGER appears under status=unreviewed
    with urllib.request.urlopen(f"{review_server}/api/candidates?status=unreviewed") as resp:
        unrev_data = json.loads(resp.read().decode("utf-8"))
        unrev_ids = [c["id"] for c in unrev_data["candidates"]]
        assert 1 not in unrev_ids, "Reviewed candidate 1 should NOT appear under unreviewed"
        assert 2 in unrev_ids, "Sibling candidate 2 must remain in unreviewed"
        assert 3 in unrev_ids, "Candidate 3 must remain in unreviewed"

    # Verify directly in SQLite database
    con = sqlite3.connect(str(db_file))
    row1 = con.execute("SELECT review_status FROM candidates WHERE id=1").fetchone()[0]
    row2 = con.execute("SELECT review_status FROM candidates WHERE id=2").fetchone()[0]
    assert row1 == "unreviewed"
    with sqlite3.connect(SurveyDatabaseAdapter(str(db_file)).review_db_path) as overlay:
        assert overlay.execute("SELECT status FROM candidate_review_overlay WHERE candidate_id=1").fetchone()[0] == "interesting"
    assert row2 == "unreviewed", "Reviewing event 1 must NOT automatically alter event 2"
    con.close()


def test_star_review_and_bulk_classification(review_server, test_env):
    """Verify:

    1. Star-level classification applies without modifying individual event statuses.
    2. Bulk-classification of remaining unreviewed events updates candidate 2 with audit provenance.
    """
    db_file, _ = test_env

    # Review Candidate 1 first
    con = sqlite3.connect(str(db_file))
    db.record_review(con, 1, "interesting", "Event 1 dip")
    con.close()

    # 1. Submit star-level classification
    star_req_data = json.dumps({
        "tic": "TIC:1001",
        "status": "known",
        "star_class": "eclipsing_binary",
        "note": "Known contact binary system"
    }).encode("utf-8")
    req = urllib.request.Request(f"{review_server}/api/star/review", data=star_req_data,
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as resp:
        assert resp.status == 200
        res = json.loads(resp.read().decode("utf-8"))
        assert res["success"] is True
        assert res["star"]["star_status"] == "known"
        assert res["star"]["star_class"] == "eclipsing_binary"

    # Verify candidate 2 is still unreviewed
    con = sqlite3.connect(str(db_file))
    assert con.execute("SELECT review_status FROM candidates WHERE id=2").fetchone()[0] == "unreviewed"
    con.close()

    # 2. Submit bulk classification for TIC:1001
    bulk_data = json.dumps({
        "tic": "TIC:1001",
        "status": "known",
        "note": "Classified remaining events as known EB",
        "reviewer": "test_bulk_user",
        "only_unreviewed": True
    }).encode("utf-8")
    req_bulk = urllib.request.Request(f"{review_server}/api/star/bulk-review", data=bulk_data,
                                      headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req_bulk) as resp:
        assert resp.status == 200
        bulk_res = json.loads(resp.read().decode("utf-8"))
        assert bulk_res["success"] is True
        assert bulk_res["updated_count"] == 1  # Only candidate 2 was unreviewed
        assert bulk_res["star"]["progress"] == "reviewed"
        assert bulk_res["star"]["unreviewed_events"] == 0

    # Verify candidate 2 was updated and has an individual audit trail in review_history
    con = sqlite3.connect(str(db_file))
    cand1_status = con.execute("SELECT review_status FROM candidates WHERE id=1").fetchone()[0]
    cand2_status = con.execute("SELECT review_status FROM candidates WHERE id=2").fetchone()[0]
    assert cand1_status == "interesting", "Bulk review must NOT overwrite previously reviewed candidate 1"
    assert cand2_status == "unreviewed"
    con.close()
    con = sqlite3.connect(SurveyDatabaseAdapter(str(db_file)).review_db_path)
    assert con.execute("SELECT status FROM candidate_review_overlay WHERE candidate_id=2").fetchone()[0] == "known"
    hist = con.execute("SELECT previous_status, new_status, reviewer FROM review_history WHERE candidate_id=2").fetchone()
    assert hist == ("unreviewed", "known", "test_bulk_user")
    con.close()

    # 3. Verify that candidate 1 remains fully visible under status=interesting
    # even though its parent star TIC:1001 is classified as known!
    with urllib.request.urlopen(f"{review_server}/api/candidates?status=interesting") as resp:
        int_data = json.loads(resp.read().decode("utf-8"))
        int_ids = [c["id"] for c in int_data["candidates"]]
        assert 1 in int_ids, "Candidate 1 must remain visible under status=interesting filter"



def test_web_security_and_validation(review_server):
    # Test invalid status returns 400
    bad_data = json.dumps({"candidate_id": 1, "status": "invalid_status"}).encode("utf-8")
    req = urllib.request.Request(f"{review_server}/api/review", data=bad_data,
                                 headers={"Content-Type": "application/json"})
    with pytest.raises(urllib.error.HTTPError) as exc_info:
        urllib.request.urlopen(req)
    assert exc_info.value.code == 400

    # Test non-existent candidate detail returns 404
    with pytest.raises(urllib.error.HTTPError) as exc_info:
        urllib.request.urlopen(f"{review_server}/api/candidate/99999")
    assert exc_info.value.code == 404


def test_plot_rendering_fallback(review_server):
    # Plot endpoint returns valid PNG even when raw NPZ is not on disk (informative fallback)
    with urllib.request.urlopen(f"{review_server}/api/plot/zoom/1?hours=8") as resp:
        assert resp.status == 200
        assert resp.headers.get("Content-Type") == "image/png"
        png_bytes = resp.read()
        assert png_bytes.startswith(b"\x89PNG")


def test_star_centered_explore_api(review_server):
    """Verify star-centered Explore API returns one row per star with expandable constituent events."""
    with urllib.request.urlopen(f"{review_server}/api/stars") as resp:
        assert resp.status == 200
        data = json.loads(resp.read().decode("utf-8"))
        assert data["total_stars"] == 2
        assert len(data["stars"]) == 2

        star1 = next(s for s in data["stars"] if s["clean_tic"] == "1001")
        assert star1["total_events"] == 2
        assert len(star1["events"]) == 2
        assert star1["max_score"] == 12.5
        assert star1["events"][0]["score"] == 12.5
        assert star1["events"][1]["score"] == 8.2

    # Test searching by TIC
    with urllib.request.urlopen(f"{review_server}/api/stars?search=2002") as resp:
        data = json.loads(resp.read().decode("utf-8"))
        assert data["total_stars"] == 1
        assert data["stars"][0]["clean_tic"] == "2002"

    # Test sorting by score
    with urllib.request.urlopen(f"{review_server}/api/stars?sort_by=score&sort_order=desc") as resp:
        data = json.loads(resp.read().decode("utf-8"))
        assert data["stars"][0]["clean_tic"] == "2002"  # Score 15.0 > 12.5


def test_undo_event_review_restores_prior_status(review_server, test_env):
    """Verify safe undo restores previous event review status in audit history."""
    db_file, _ = test_env

    # 1. Update candidate 1 to 'interesting'
    req1 = urllib.request.Request(
        f"{review_server}/api/review",
        data=json.dumps({"candidate_id": 1, "status": "interesting", "note": "Dip"}).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req1) as resp:
        assert json.loads(resp.read().decode("utf-8"))["status"] == "interesting"

    # 2. Undo the review
    undo_req = urllib.request.Request(
        f"{review_server}/api/review/undo",
        data=json.dumps({"candidate_id": 1}).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(undo_req) as resp:
        assert resp.status == 200
        res = json.loads(resp.read().decode("utf-8"))
        assert res["success"] is True
        assert res["undo"]["reverted_to"] == "unreviewed"

    # Verify candidate 1 is unreviewed again in database
    con = sqlite3.connect(str(db_file))
    assert con.execute("SELECT review_status FROM candidates WHERE id=1").fetchone()[0] == "unreviewed"
    con.close()


def test_undo_star_review_restores_prior_status(review_server, test_env):
    """Verify safe undo restores previous star review status."""
    db_file, _ = test_env

    # 1. Update star review
    req1 = urllib.request.Request(
        f"{review_server}/api/star/review",
        data=json.dumps({"tic": "TIC:1001", "status": "interesting", "star_class": "candidate_transit"}).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req1) as resp:
        assert json.loads(resp.read().decode("utf-8"))["success"] is True

    # 2. Undo star review
    undo_req = urllib.request.Request(
        f"{review_server}/api/star/undo",
        data=json.dumps({"tic": "TIC:1001"}).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(undo_req) as resp:
        assert resp.status == 200
        res = json.loads(resp.read().decode("utf-8"))
        assert res["success"] is True
        assert res["undo"]["reverted_to"] == "unreviewed"

    con = sqlite3.connect(str(db_file))
    assert con.execute("SELECT status FROM star_reviews WHERE tic='TIC:1001'").fetchone() is None
    con.close()
    con = sqlite3.connect(SurveyDatabaseAdapter(str(db_file)).review_db_path)
    assert con.execute("SELECT status FROM star_reviews WHERE tic='TIC:1001'").fetchone()[0] == "unreviewed"
    con.close()


def test_overlay_review_database_leaves_source_database_unmodified(test_env, tmp_path):
    """Verify non-destructive review storage: source database remains bit-for-bit unchanged."""
    import hashlib
    db_file, root_dir = test_env

    # Compute source DB SHA-256 before reviews
    h_before = hashlib.sha256(db_file.read_bytes()).hexdigest()

    overlay_db = tmp_path / "overlay_reviews.sqlite"

    class OverlayHandler(AnomalyReviewHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(db_file, root_dir, *args, review_db=str(overlay_db), read_only=True, **kwargs)

    server = HTTPServer(("127.0.0.1", 0), OverlayHandler)
    host, port = server.server_address
    t = threading.Thread(target=server.serve_forever, daemon=True)
    t.start()
    base_url = f"http://{host}:{port}"

    try:
        # Submit star and event review
        req_star = urllib.request.Request(
            f"{base_url}/api/star/review",
            data=json.dumps({"tic": "TIC:1001", "status": "known", "star_class": "eclipsing_binary"}).encode("utf-8"),
            headers={"Content-Type": "application/json"}
        )
        with urllib.request.urlopen(req_star) as resp:
            assert resp.status == 200

        req_ev = urllib.request.Request(
            f"{base_url}/api/review",
            data=json.dumps({"candidate_id": 1, "status": "known", "note": "Overlay test"}).encode("utf-8"),
            headers={"Content-Type": "application/json"}
        )
        with urllib.request.urlopen(req_ev) as resp:
            assert resp.status == 200

        # Verify overlay database stored the review
        con_ov = sqlite3.connect(str(overlay_db))
        ov_star = con_ov.execute("SELECT status, star_class FROM star_reviews WHERE tic='TIC:1001'").fetchone()
        assert ov_star == ("known", "eclipsing_binary")
        ov_cand = con_ov.execute("SELECT status FROM candidate_review_overlay WHERE candidate_id=1").fetchone()
        assert ov_cand == ("known",)
        con_ov.close()

        # Check API reflects the overlay reviews
        with urllib.request.urlopen(f"{base_url}/api/star/TIC:1001") as resp:
            detail = json.loads(resp.read().decode("utf-8"))
            assert detail["star_status"] == "known"
            assert detail["events"][0]["review_status"] == "known"
    finally:
        server.shutdown()
        server.server_close()

    # Verify source database file is bit-for-bit UNTOUCHED!
    h_after = hashlib.sha256(db_file.read_bytes()).hexdigest()
    assert h_after == h_before, "Source research database must NOT be modified when overlay review DB is active!"


def test_sample_dataset_loads_and_serves_investigation(tmp_path):
    """Verify that the bundled sample dataset examples/sample_survey_review.sqlite serves cleanly."""
    sample_db = Path("examples/sample_survey_review.sqlite").resolve()
    sample_root = Path("examples/sample_data").resolve()
    assert sample_db.is_file(), "Sample database must exist"

    class SampleHandler(AnomalyReviewHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(str(sample_db), str(sample_root), *args, **kwargs)

    server = HTTPServer(("127.0.0.1", 0), SampleHandler)
    host, port = server.server_address
    t = threading.Thread(target=server.serve_forever, daemon=True)
    t.start()
    base_url = f"http://{host}:{port}"

    try:
        # Check stats
        with urllib.request.urlopen(f"{base_url}/api/stats") as resp:
            stats = json.loads(resp.read().decode("utf-8"))
            assert stats["stars"]["total"] == 4
            assert stats["events"]["total"] == 8

        # Check EBLM J0048-66 (TIC 39018208)
        with urllib.request.urlopen(f"{base_url}/api/star/39018208") as resp:
            star = json.loads(resp.read().decode("utf-8"))
            assert star["clean_tic"] == "39018208"
            assert star["star_status"] == "interesting"
            assert star["system"]["is_periodic"] is True
            assert star["system"]["period_days"] == pytest.approx(2.41824, rel=1e-4)
            assert len(star["events"]) == 4

        # Check light curve endpoint returns time and flux
        with urllib.request.urlopen(f"{base_url}/api/star/39018208/lightcurve") as resp:
            lc = json.loads(resp.read().decode("utf-8"))
            assert "time" in lc and len(lc["time"]) > 100
            assert "flux_ppm" in lc and len(lc["flux_ppm"]) == len(lc["time"])
            assert len(lc["events"]) == 4
    finally:
        server.shutdown()
        server.server_close()


