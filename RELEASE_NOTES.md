# 0.3.0b2 — Audit-corrected bounded beta

Tag: `v0.3.0-beta.2`. Prepared 2026-10-16 from verified correction commit `ac15694`; publication authorized 2026-11-03 for `zgraper/TESS_Anomaly_Survey`. Package metadata, application version and citation use the canonical Python prerelease version `0.3.0b2`.

- Read-only scientific source access with run-bound review overlays and launch guards.
- Exact archive product resolution, FITS target/sector/cadence checks, isolated caches and scoped exports.
- Immutable execution fingerprints, bounded scheduling, product deduplication and explicit restart recovery.
- Corrected experiment/progress contract, plot cleanup and supported two-hour custom threshold controls.
- Separate `coherence_v2` post-detection gate; historical v1 classifications and frozen B0/B1 are preserved.
- Numeric MAST Products identifiers, dependency-free CSV export, same-origin/CSRF protection, truthful evidence states and documentation.

This prerelease preserves the original `v0.3.0-beta` tag. See [audit corrections](docs/AUDIT_CORRECTIONS.md) for acceptance results and remaining limitations.

All 17 fixes and three documented restrictions from R01–R20 are included. The correction suite passed 108 tests with one expected private-catalog module skip; release preparation changes only version/publication metadata and receives an installed-package workflow smoke test.

Supported scope: fixed two-hour detection with explicit experimental custom thresholds; phase folding at recorded event-timing periods; Windows/Python 3.14.4 acceptance. Linux/macOS, other Python versions, physical tablets and long workloads remain unverified. Interactive plots depend on CDN Plotly; offline static plotting requires local data and the optional plotting dependencies. `coherence_v2` is separately versioned post-detection logic, with no new recovery-rate or false-alarm claims. Historical classifications, datasets, frozen B0/B1 and threshold remain unchanged.

---

# Release Notes

## v0.3.0-beta — Star-Centered Research Workspace & Public Beta (2026-10-16)

TESS Anomaly Survey `v0.3.0-beta` introduces a complete browser-based astronomical research workspace, live streaming survey execution against NASA MAST, star-centered investigation views, and non-destructive review overlays.

### Highlights & Key Features

1. **Star-Centered Research Workspace (`tess-anomaly review-ui`)**
   - **Explore Stars Catalog**: Star-aggregated triage view (1 row per target star by default) expandable into individual candidate events, with comprehensive filtering (unreviewed, interesting, false alarm, multi-event, period-folded).
   - **Interactive Investigation Workspace**: Synchronized observation and event-zoom panels, with folding at recorded periods derived from candidate-event timings. Browser Lomb–Scargle analysis and manual period overrides are not implemented.
   - **Non-Destructive Review Overlays**: Scientific candidate classifications (e.g., `interesting`, `known_transit`, `stellar_variability`, `instrumental_artifact`) and observational notes are persisted in companion `.reviews.sqlite` databases, ensuring historical research databases remain strictly immutable.
   - **Bundled Public Demo**: Ready out-of-the-box with a 4-star demonstration survey (`examples/sample_survey_review.sqlite`) and synthetic light curves (`examples/sample_data/`).

2. **Streaming Survey Engine & NASA MAST Integration**
   - **Direct MAST Streaming**: Query and stream calibrated light curve products directly from NASA's Barbara A. Mikulski Archive for Space Telescopes (MAST) via `astroquery.mast`.
   - **Real-Time Telemetry**: Survey execution progress is displayed live in the browser via asynchronous HTTP polling, reporting target completion rates, cadence throughput, memory utilization, and anomaly yields.
   - **Resumable Multi-Product Processing**: Exact archive product deduplication `(tic, sector, product_filename)` allows researchers to safely resume interrupted surveys without skipping distinct high-cadence products (e.g., SPOC 20-second vs 2-minute cadence).

3. **Security Hardening & Portability**
   - **Anti-DNS Rebinding**: Strict `Host` header enforcement rejecting requests not directed to `127.0.0.1` or `localhost`.
   - **Anti-CSRF & Origin Verification**: Mutating `POST` endpoints enforce origin validation to protect researchers running local servers.
   - **Path Traversal Defenses**: Strict validation on light-curve file paths and output directory destinations.
   - **Platform Validation**: The correction acceptance run used Windows and Python 3.14.4. Linux/macOS and other Python versions remain unverified. Console logging uses ASCII-safe formatting.

4. **Scientific Integrity & Frozen Detector Preservation**
   - Detector implementations and baseline parameters ($T = 73.77965536097444$ for 2-hour localized box filtering) remain byte-frozen and verified against original SHA-256 hashes.
   - Research databases from prior phases (Phase 6B independent validation and Phase 7 development survey) remain intact and reproducible.

---

## v0.1.0 — Initial Engineering Prototype (2026-10-12)

Initial exploratory software release establishing the CLI batch scanning framework, calibration workflows, and frozen B0/B1 baseline detector algorithms.
