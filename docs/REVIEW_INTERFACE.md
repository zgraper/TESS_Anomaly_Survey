# Star-Centered Research Investigation Workspace Guide

The **TESS Anomaly Survey Research Workspace** is a local, browser-based platform designed for astrophysical researchers, archival vetters, and open-source contributors to inspect, triage, and annotate deep photometric dimming events detected in NASA TESS light curves.

The interface adheres to a **star-centered research philosophy** rather than an event-centered alert stream. The primary entity of scrutiny is the target star (`TIC`), with all detected events, physical context, cross-sector observations, and astrophysical parameters presented as a coherent whole.

---

## 1. Quick Start & Multi-Platform Installation

The application runs locally on standard Python (3.10–3.14) with zero remote web server requirements and binds exclusively to loopback (`127.0.0.1`) by default.

### Clone and Environment Setup

```bash
# Clone the repository
git clone https://github.com/zgraper/tess-anomaly-survey.git
cd tess-anomaly-survey

# Create and activate a virtual environment
# Windows (PowerShell):
python -m venv .venv
.\.venv\Scripts\Activate.ps1

# Linux / macOS (bash / zsh):
python3 -m venv .venv
source .venv/bin/activate

# Install the survey package and dependencies
pip install -e ".[test,plots]"
```

### Launch with Bundled Example Dataset (Zero Internet / Immediate Demo)

The repository includes a self-contained demonstration database containing confirmed benchmark eclipsing binaries (EBLM J0048-66), unclassified localized dimmings, flare phantom rejects, and quiet control stars:

```bash
# Windows:
python -m tess_anomaly_survey review-ui --db examples/sample_survey_review.sqlite --root examples/sample_data

# Linux / macOS:
python3 -m tess_anomaly_survey review-ui --db examples/sample_survey_review.sqlite --root examples/sample_data
```

Once launched, open your web browser to:
👉 **`http://127.0.0.1:8000/`**

---

## 2. Launch Workflows & Database Configurations

### A. Non-Destructive Inspection of Research Databases (`--read-only` & `--review-db`)
To inspect frozen research databases (e.g. `development_survey_9753.sqlite` or `phase6b_validation.sqlite`) without risking any modification to original research files, use `--read-only`:

```bash
python -m tess_anomaly_survey review-ui \
    --db "path/to/development_survey.sqlite" \
    --read-only \
    --review-db "my_expert_reviews.sqlite"
```
- **Read-Only Protection:** The underlying survey database is mounted with SQLite URI `?mode=ro`.
- **Overlay Storage:** All star classifications, event triage statuses, researcher notes, and audit undo histories are saved to `my_expert_reviews.sqlite`, leaving the source database 100% bit-for-bit unchanged.

### B. On-Demand NASA MAST Streaming & Local Cache (`--cache-dir`)
If you do not have local preprocessed `.npz` or `.fits` files on disk, the workspace automatically retrieves public SPOC 2-minute cadence light curves directly from NASA's Mikulski Archive for Space Telescopes (MAST):

```bash
python -m tess_anomaly_survey review-ui \
    --db examples/sample_survey_review.sqlite \
    --cache-dir ~/.tess_anomaly_survey/cache
```
- Downloaded FITS products are validated and stored in `--cache-dir`.
- Subsequent visits to the same star load instantaneously from local disk or memory cache.

### CLI Launch Options Reference

| Option | Default | Description |
| :--- | :--- | :--- |
| `--db` | *Required* | Path to SQLite survey database (`sample_survey_review.sqlite`, `development_survey_9753.sqlite`, etc.) |
| `--root` | `None` | Local directory containing preprocessed `.npz` or `.fits` files |
| `--host` | `127.0.0.1` | Network binding host (restricted to loopback for security) |
| `--port` | `8000` | Port to serve |
| `--cache-dir` | `~/.tess_anomaly_survey/cache` | Local directory for on-demand downloaded MAST light curve FITS products |
| `--review-db` | `None` | External overlay SQLite database path for non-destructive review storage |
| `--read-only` | `False` | Enforce strictly read-only mode on the source database |

---

## 3. Workspace Architecture & Features

### View 1: Star-Centered Explore Page
- **One Row Per Star:** Displays each target star (`TIC`), observing sector(s), total detected events, maximum matched-filter score ($S_{\max}$), primary photometric verification class, and system periodicity.
- **Accordion Event Expansion:** Click any star row to expand its constituent detected events sub-table, displaying individual center times (BTJD), SNR scores, durations, and depths.
- **Multi-Domain Filters:**
  - *Review Status:* All, Unreviewed, Interesting, Known Systems, Artifacts, Dismissed.
  - *Photometric Verification:* Localized Dimmings, Continuous Stellar Variability, Flare Phantoms, Artifacts.
  - *System Classification:* Periodic Systems ($N \ge 3$ coherent dips) vs. Isolated / Single Dips.
  - *Search & Sorting:* Instant search across TIC IDs, user notes, and sectors; sort by Score, Events, Sector, TIC, or Period.

### View 2: Light-Curve Investigation Workspace
- **Synchronized Navigation:** Preserves your active search and filter order. Use `← Previous Star` and `Next Star →` to triage your shortlist without losing context.
- **Full Observation Curve (Plotly.js):** Interactive overview of the entire ~27-day sector. Displays dashed vertical markers for all detected events on the star. Clicking any marker jumps the zoom view to that event.
- **Event Zoom Curve:** High-cadence interactive view centered on the active event ($\pm 6\text{h}$, $\pm 12\text{h}$, or $\pm 24\text{h}$ span), showing the shaded duration template window ($2.0\text{h}$) and measured depth.
- **Phase-Folded Curve:** When periodicity analysis detects coherent orbital cycles, displays the folded light curve with binned median trend and odd/even transit disparity markers.
- **Static Matplotlib Fallback:** Click `Switch to Static PNG` to toggle server-rendered static figures if running on a machine without WebGL/JavaScript acceleration.

### Multi-Domain Taxonomy Separation
To preserve scientific rigor, five conceptual tiers are strictly separated:
1. **Detector Metric:** Matched-filter score ($S \ge 73.78$) and duration kernel ($2.0\text{h}$).
2. **Photometric Verification:** 5-state machine outcome (`verified_localized_dimming`, `continuous_stellar_variability`, `flare_induced_apparent_dip`, etc.) and measured true depth.
3. **Coincidence Risk:** Empirical Monte Carlo cross-star coincidence risk flag (None / Clean vs. Multi-Star Simultaneous Glitch).
4. **Provisional Interpretation:** Astrophysical hypothesis (e.g. EBLM secondary eclipse, candidate transit, pulsation trough).
5. **Researcher Review:** Human expert judgment (star-level and event-level classifications).

### Scope Isolation & Safe Undo Mechanism
- **Event-Level Independence:** Reviewing an individual event updates only that specific event. Sibling events on the same star remain unmodified.
- **Star-Level Independence:** Reviewing the star classifies the astrophysical target without altering individual event review statuses.
- **Bulk Option:** An explicit checkbox allows bulk-applying a star classification to remaining unreviewed events, generating individual audit entries for each event.
- **Safe Undo Stack:** Click `↺ Undo Star` or `↺ Undo Event` to safely revert any accidental review modification.

---

## 4. Experiment Management & Background Survey Streaming

The research workspace allows astronomers to define, execute, monitor, and export independent photometric survey experiments directly from the browser:

### A. Experiment Creation Wizard & Pre-Flight Validation
- Click **`+ New Experiment`** in the navigation header to launch the experiment wizard.
- **Scientific Identity:** Specify an experiment name, scientific objective, and experimental description.
- **Target Ingest:** Input target TICs via comma-, space-, or newline-delimited text or paste raw CSV manifests with headers.
- **Pre-Flight Resource Estimation:** Click *Validate Targets & Estimate* to compute unique target count, total observations across chosen sectors, download volume (MB), and estimated wall-clock runtime. Non-numeric or out-of-range tokens are automatically detected and flagged.
- **Detector Configuration:** Select the published frozen detector preset ($T = 73.7797$, fixed 2.0h box kernel) with B0/B1 SHA-256 provenance, or select an experimental threshold. Duration is fixed at two hours; unsupported durations are rejected.

### B. Background Execution & Real-Time Telemetry
- Streaming surveys run in a background daemon thread with concurrency locks preventing simultaneous write collisions on the output SQLite database.
- **Live Progress Monitor:** Real-time progress bar, observation ingestion counts, pipeline throughput (observations/hour), elapsed time, and dynamic ETA estimation.
- **Intermediate Result Exploration:** As observations finish and are committed to the SQLite catalog, researchers can switch to the **Explore Stars** tab and immediately investigate newly detected candidates without waiting for the survey run to conclude.
- **Run Control:** Clean, safe cancellation via **`Abort Survey`**, with subsequent **`Resume Incomplete Run`** support that skips already completed targets using database state.

### C. Catalog & Provenance Export (CSV & JSON)
- **Export Catalog (CSV):** Generates a CSV catalog of detected candidates joined with human review annotations, preceded by commented provenance headers recording detector thresholds, kernel durations, software versions, and scientific caveats.
- **Export Provenance (JSON):** Produces machine-readable JSON containing full configuration snapshots, B0/B1 verification hashes, and candidate arrays.

---

## 5. Current vs. Planned Capabilities

| Capability | Current Status | Notes |
| :--- | :--- | :--- |
| **Star-Centered Explore Page** | ✅ Fully Implemented | Expandable accordion rows, multi-tier filters, sorting |
| **Interactive Plotly Curves** | ✅ Fully Implemented | Full observation, zoom, and phase-folded interactive plotting |
| **Matplotlib Static Fallback** | ✅ Fully Implemented | In-memory PNG rendering via `/api/plot/...` |
| **On-Demand MAST Streaming** | ✅ Fully Implemented | Automatic FITS acquisition and local disk caching |
| **Historical Benchmark Auto-Discovery** | ✅ Fully Implemented | Auto-detects and mounts historical research DBs read-only |
| **Non-Destructive Overlay DB** | ✅ Fully Implemented | `--read-only` and overlay reviews preserving source SHA-256 |
| **Experiment Creation & Wizard** | ✅ Fully Implemented | Pre-flight target validation and resource estimation |
| **Background Streaming Pipeline** | ✅ Fully Implemented | Concurrent execution with progress callbacks and database safety |
| **Survey Progress & Live Logs** | ✅ Fully Implemented | Telemetry cards, throughput, ETA, and audit log terminal |
| **Abort & Safe Resume** | ✅ Fully Implemented | Cancellation event handling and exact-product resumption with identical scientific configuration |
| **Catalog & Provenance Export** | ✅ Fully Implemented | Standardized CSV and JSON exports with scientific caveats |
| **Dual-Scope Review & Safe Undo** | ✅ Fully Implemented | Scope isolation between star and event reviews with audit history |
| **Localhost Security Guard** | ✅ Fully Implemented | Validates loopback Host/authority and Origin on reads and writes; browser writes require a CSRF token |


The process must remain open while a survey runs. Restart marks unfinished jobs interrupted; resume explicitly. Plotly requires its CDN script. Offline static plots require the `plots` extra and local data. The period view folds recorded event-timing fits; it does not provide a browser periodogram or custom-period fitting. Missing evidence is displayed as unknown/not run. Synthetic demonstration plots and exports are labeled and never acquire substitute MAST observations.
