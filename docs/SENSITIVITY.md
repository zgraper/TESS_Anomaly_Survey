# Sensitivity, calibration and the 1% alert target

## Three different controls

1. **Detection threshold:** Higher empirical maximum-score quantile (e.g. 99% vs 95%) yields fewer *calibration-set* windows that trigger, with potentially lower recovery of genuine weak events.
2. **Review budget:** A ranked top-100 export only limits researcher workload; it does not change detection sensitivity or statistical rate.
3. **Artifact vetting:** A separate future science layer can remove common instrumental signatures or previously cataloged variable stars, but must undergo its own injection-recovery evaluation.

The command `calibrate --target-alert-fraction 0.01` computes a per-window 99th-percentile **maximum localized detector score** over presumed-background calibration windows. We use a conservative nearest-rank quantile and flag only maxima *strictly greater* than the cutoff. For the optional noise-tertile calibration, cutpoints and per-bin quantiles are measured only on calibration windows.

**This is not a guaranteed ceiling.** Validation windows from unseen stars can have heavier score tails, and overlapping windows within stars are correlated. The `validate` command reports actual flagged-window incidence and star-level-bootstrap uncertainty on a wholly separate TIC cohort. Because those controls are not confirmed astrophysical negatives, we should say 'presumed-background alert rate' rather than astrophysical 'false-positive rate'.

Suggested minimums (rough heuristics, not statistical guarantees):

| Target | Calibration windows | Distinct calibration TICs | Caveat |
| --- | ---: | ---: | --- |
| 5% | 10,000 | 1,000+ | Score tails remain sector-dependent |
| 1% | 10,000–50,000 | 2,000–5,000+ | Validate on independent TICs |
| 0.5% | 20,000–100,000 | 2,000–5,000+ | Tail estimation more fragile |
| 0.1% | 100,000+ | 5,000–10,000+ | Extreme quantiles and drift require special care |

A larger number of *overlapping windows from the same few stars* does not replace having many independent stars. The program rejects configurations below minimum sample counts by default, but minimum counts are not an assurance of calibration quality.

## How the survey differs from B2.2

B2.2 assessed one prepared 512-cadence window per star with fixed injected signals. The survey searches **many overlapping windows per complete source observation**. Even the same detector can have a different score distribution and alert burden. Consequently, B2.2's existing 5% thresholds must *not* be recycled as validated 1% thresholds for large-scale scanning.

## Quantify missed signals

After calibrating stricter cutoffs, perform a **new star-disjoint injection-recovery study** with the same survey scanner, track recall as a function of depth, duration and event morphology, and test alert incidence on distinct controls. Do not re-use the formerly untouched B2.2 cohort as independent validation after its results have been inspected.

## Scalability

- The program reads NPZ files one at a time, creates a window-level audit trail and commits per completed source file.
- `--max-windows-per-star N` can reduce load for a *pilot*, but changes the scanning population. Calibration and actual scanning must use the same value.
- The 512-cadence window, overlapping stride, noise estimator and frozen 2-hour box match B1's core scoring method, but full-observation coverage and post-scan deduplication are new survey behaviors.
- `--candidates-per-window 3` limits how many separate local peaks are inserted into the candidate library *per flagged window*. Window-level incidence is still computed from **all raw candidate scores**. This is a deliberate review bottleneck, not an astrophysical filter.
- The program reports cumulative window alert frequency during a scan. It can differ from validation; stop and investigate major drift, but do not recalibrate on the same survey data while claiming its independent performance.
