# Independent audit correction register

Date: 2026-10-16. Branch: `codex/audit-corrections`.
Original audited snapshot/tag: `0681349afaa51ee557b61913f17a468754a3a16e` / `v0.3.0-beta`.
Tested implementation: `3543a9b` (subsequent changes add the acceptance record and a standalone browser-function regression only).

**Recommendation: GO for the bounded public software beta, subject to the owner's review and publication approval.** No unresolved release-blocking finding remains within the workflow tested below. Nothing has been pushed, tagged for publication, or released. This recommendation does not validate historical research conclusions or establish new scientific performance for `coherence_v2`.

## Issue register

Paths below are repository-relative. Test names identify executable regression cases; browser/installation checks are separately recorded below. “Mitigated” means the advertised capability is restricted to what is implemented, rather than adding an unvalidated feature.

| ID | Disposition | Implementation / restriction | Verification |
|---|---|---|---|
| R01 | Fixed | [`.gitattributes`](../.gitattributes) pins both frozen detector files to LF; integrity checks unchanged. | `test_survey.py::test_frozen_sources`; Windows `core.autocrlf=true` clone, wheel and sdist bytes match both canonical hashes. |
| R02 | Fixed | [`web.py`](../src/tess_anomaly_survey/web.py): experiment/run selection, progress/log response parsing, observed counts, native downloads, clearing stale catalog; [`experiments.py`](../src/tess_anomaly_survey/experiments.py): historical telemetry scoped to native run. | `test_experiments_workspace.py::test_web_server_experiment_api`; `test_audit_integrity.py::test_historical_telemetry_is_native_run_scoped`; real browser creation, 1/1 observation progress, two candidates, CSV/JSON downloads. |
| R03 | Fixed | `experiments.py::SurveyExecutionManager.launch_survey` rejects read-only/historical experiments and foreign run ownership before opening output. | `test_launch_rejects_readonly_and_foreign_runs_before_output`; original-database protection checks reject launches and preserve bytes. |
| R04 | Fixed | `web.py::SurveyDatabaseAdapter` requires explicit scope for multi-run sources and filters source tables; joins use system/observation identity. `experiments.py::export_run_catalog` exports native provenance. Legacy CLI exports use scoped read-only adapters. | `test_same_tic_two_runs_are_isolated_and_source_immutable` checks distinct candidates, reviews, provenance and JSON/CLI exports for the same TIC; historical telemetry test. |
| R05 | Fixed | `web.py::resolve_light_curve_series` binds cache/plots to recorded source/product/data URI/hash; candidate selection chooses the corresponding observation. `streaming.py::validate_fits_identity` checks target, sector and two-minute cadence. | `test_product_cache_and_fits_identity`; different archive-product resume test; failed identity/retry regression; real browser reports the exact recorded Sector 2 product. |
| R06 | Fixed | `streaming.py::execution_configuration` fingerprints effective settings and scientific module hashes; existing runs reject mismatches. New provenance does not reuse the historical calibration hash. | `test_resume_rejects_changed_configuration`; installed real run/export provenance inspected. |
| R07 | Fixed | `streaming.py::stream_pipeline` resolves product identity before scheduling, deduplicates manifest requests and resolved products, and skips completed products on resume. | `test_cancel_commit_resume_duplicates_and_restart`; `test_streaming.py::test_streaming_idempotent_resumption` and `test_survey_resumption_different_archive_products`. |
| R08 | Fixed | Bounded outstanding futures; cancellation stops new submissions and persists already-started successful work. Restart marks orphaned running jobs interrupted. Successful retry replaces failed observation state. | `test_cancel_commit_resume_duplicates_and_restart`, `test_partial_failure_recovery_replaces_failed_observation`; manager lifecycle test. Shutdown survival is not advertised. |
| R09 | Fixed | `web.py` and run exports use `not_run`, `unknown`, and nullable flare evidence when results are absent; browser includes missing-verification filter. | Same-TIC regression asserts missing verifier, flare and coincidence states; browser shows unknown/not-run coincidence on real data. |
| R10 | Fixed, separately versioned | [`vetting.py::verify_event_photometry_v2`](../src/tess_anomaly_survey/vetting.py) requires timestamp-coherent adjacent good cadences; streaming records `coherence_v2`. Legacy verifier and historical classifications remain intact. | `test_versioned_coherence_comparison`: isolated dips accepted by v1 are rejected by v2; contiguous injected dips accepted by both; masked gaps break coherence. No new population-level performance claim. |
| R11 | Fixed | `streaming.py::resolve_mast_product` passes numeric `obsid` to the direct MAST Products fallback. | Live fallback resolved TIC 39018208 Sector 2, numeric observation 60873935, product listed below. |
| R12 | Fixed | Source SQLite connections always use `mode=ro`; default run-bound review overlays, source-binding checks and candidate membership validation prevent cross-run review writes. CLI review/export honors this boundary. | Same-TIC isolation test, overlay/source SHA tests, original-database protection checks; browser review persists after server restart. |
| R13 | Mitigated through truthful controls | Custom thresholds work and are labeled experimental; duration is restricted to supported fixed 2h and other values are rejected. UI/exports show actual threshold and separate matched-filter depth from photometric minimum. | Workspace CRUD/API tests; actual real-run threshold and distinct depth measurements inspected in browser. No new detector duration introduced. |
| R14 | Fixed | Synthetic identity appears in demo description, catalog, plots and exports; missing synthetic NPZ cannot fall back to real MAST flux. | Sample dataset test; synthetic browser inspection and same-TIC switching; source resolver explicitly rejects absent synthetic files. |
| R15 | Fixed | [`README.md`](../README.md) states 512 samples ≈17.07h and stride 256 ≈8.53h, distinguishing endpoint span. | Arithmetic/source review; detector settings unchanged. |
| R16 | Fixed | `web.py::_request_allowed` applies loopback client, exact Host authority and same-origin checks to read/write routes; browser writes require CSRF token and JSON. | `test_review_ui.py::test_read_authority_security` parameter cases, `test_browser_write_requires_csrf_token`, workspace Host/Origin tests. |
| R17 | Fixed | [`cli.py::do_stream`](../src/tess_anomaly_survey/cli.py) uses standard-library CSV export, with run scope and overlay status. | Ordinary installed environment has no pandas; streaming CSV export acceptance succeeds. |
| R18 | Mitigated through truthful restrictions | README, release notes and review guide describe candidate-timing periodicity and folding at recorded periods. No browser Lomb–Scargle/manual-period claim. Interactive plots require CDN Plotly; offline static PNG requires local photometry and `plots` extra. | Documentation/source inspection, periodicity and static endpoint tests, browser phase panel on synthetic demo. No fully offline interactive guarantee. |
| R19 | Fixed | `web.py` purges Plotly containers, clears loading state and ignores stale responses across target/run transitions. | Repeated Explore→Investigate, event changes, next/previous targets, real↔synthetic experiment changes and revisits; full/zoom/phase panels render. |
| R20 | Mitigated with sanitized provenance | Packaged [`calibration_provenance.json`](../src/tess_anomaly_survey/calibration_provenance.json) records inherited threshold/settings/hashes and limits. Research docs remove missing private-archive references and qualify historical and platform claims. | Wheel contains provenance; Windows/Python 3.14.4 installation executed. Private catalog module skipped; Linux/macOS and other Python versions remain unverified. Full scientific reproducibility archive deferred as nonblocking. |

## Change log and commits

- `cbc2356`: Phase 1 scientific source protection, review/run/product identity; gate **88 passed, 1 skipped** before Phase 2.
- `efe46b1`: execution lifecycle, browser contract, versioned coherence method.
- `9a41015`: acceptance counters, selection, review history, browser downloads.
- `32bd30a`: source provenance, synthetic disclosure, failure recovery.
- `1f72e26`: clear investigation state and discard stale responses.
- `47cae37`: preserve legacy nonnumeric target labels.
- `4018977`: truthful custom thresholds and inherited calibration labels.
- `9a99598`: clear catalog during transitions; missing-verification filter.
- `3f508ec`: native historical telemetry, scoped overlay-aware legacy exports, remaining documentation/depth labels.
- `3543a9b`: reject late progress responses after selecting a different run. `tests/progress_race.cjs` executes the actual embedded function with delayed responses and verifies both stale-response rejection and current-response updates.

## Acceptance results

The corrected public tree was cloned with `core.autocrlf=true` on Windows. A source distribution and wheel were built from that clone, installed normally into a fresh virtual environment, and subsequently reinstalled after final corrections. Tests ran outside the source directory with pytest's `pythonpath` override disabled, so imports came from installed `site-packages`.

- Final wheel suite: **108 passed, 1 skipped, 2 warnings**, 18.96 seconds. Targeted final identity/workspace regressions: **16 passed**. Additional Node progress-race regression: **passed** against installed `web.py`.
- One module-level skip concerns `test_catalog_integrity.py`: pandas/private research catalogs are absent. This suite does not reproduce those catalogs or historical adjudication.
- Source-distribution installation and frozen-source validation passed. Dependency consistency checked with `pip check`.
- Direct MAST fallback succeeded for one real observation. No large survey was executed.
- Cancellation/restart/failure tests used small deterministic fixtures, including synthetic FITS identity mismatch and successful retry. Already-started work may finish after cancellation; no process-independent background execution is claimed.
- All nine baselined historical SQLite files remained byte-identical. Seven native historical/demo run contexts rejected direct source writes and survey launches while accepting review writes in disposable overlays.
- Same TIC in two native runs retains independent candidates, review decisions and exports. Distinct FITS products sharing TIC/sector have separate cache entries.
- Browser acceptance created `Audit acceptance — real MAST Sector 2`, run `run_20261016_143128_532c57_01`, for TIC 39018208, Sector 2, default threshold, two workers. Execution completed one observation and detected two candidates. Product: `tess2018234235059-s0002-0000000039018208-0121-s_lc.fits` (`mast:TESS/product/` prefix in its data URI).
- Inspected both events and exact-product flux. Saved an **acceptance-only** “interesting” review for event 1 with the note “Acceptance test review only; no astrophysical conclusion.” Review persisted through restart; source candidate status remained unreviewed. CSV and JSON downloaded through the browser and retained native run/product provenance and overlay review.
- Repeated navigation and real/synthetic switching rendered plots without leaving loading indicators. Synthetic example period remains explicitly illustrative. Security tests rejected invalid Host, mismatched Origin/port and missing browser CSRF token.

## Scientific and repository integrity

| Frozen artifact | SHA-256 |
|---|---|
| B0 | `9c8afe63ec2763b9db57c828b52b0e7d262bb30b203f86e559836dd5beb89f00` |
| B1 | `c5214809a053e4705c5d86f29703e7682128e3db119494ed59f6c9a191e557c8` |

Both hashes match the public working tree, Windows clone, wheel and sdist. Threshold remains `73.77965536097444`; historical configuration identity remains `dc32c430b5aae7d8af069a6359efd8c49e06736994c600f258f2852eb3269d3c`. New executions record their own complete configuration fingerprint, including the separately versioned post-detection implementation.

Original tag remains at `0681349`. The private repository remains clean at its original commit; public/private commit intersection is empty. No private history, decoding keys, confidential datasets or manuscript changes were added. Original demo database bytes are unchanged. The public repository has no push remote configured.

## Limits and retained evidence

No unresolved blocker was found for the tested Windows workflow. Linux/macOS, other supported Python versions, physical tablets, prolonged workloads and scientific recovery/FAP benchmarking remain unverified. Legacy runs whose recorded configuration does not match the complete current configuration must be preserved and followed by a new run, not silently resumed under changed settings. `coherence_v2` is a new post-detection heuristic requiring owner review; it does not retroactively change published results.

Detailed local evidence is retained outside the public repository: `installed-wheel-tests.log`, `corrections-build.log`, `final-wheel-install.log`, `sdist-install.log`, `corrected-mast-fallback.log`, `acceptance-integrity.json`, `acceptance-integrity-check.py`, browser CSV/JSON exports, and real/synthetic/progress screenshots. The companion local acceptance report links those files and records final checks. Historical paths and research databases are intentionally excluded from this public document.
