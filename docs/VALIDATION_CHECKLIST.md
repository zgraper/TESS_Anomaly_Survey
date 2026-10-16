# Survey readiness checklist

- [ ] Review repository files and choose license / author metadata.
- [ ] Create explicit star-disjoint calibration, validation and scan catalog CSVs. Exclude test stars reserved for future confirmatory claims.
- [ ] Confirm calibration controls are representative across sectors, stellar noise, magnitude and quality.
- [ ] Fit frozen-detector per-window score quantiles for intended alert fraction using enough distinct stars/windows.
- [ ] Run a **star-disjoint** validation survey; report point estimate and star-bootstrap uncertainty for window alert fraction.
- [ ] If sub-1% is not supported, revise calibration on a **new calibration partition** and validate against a different untouched set; do not repeatedly tune on the same validation cohort and describe it as independent.
- [ ] Run a new injection-recovery benchmark using the same survey windowing/threshold configuration to quantify sensitivity loss at stricter thresholds.
- [ ] Pilot scan 100 real, unreserved stars; inspect event merging and number of unrelated alerts per TIC, not just flagged-window fraction.
- [ ] Review candidate light-curve plots; identify quality events and known astrophysical phenomena.
- [ ] Only then expand to large-scale datasets, and preserve immutable detector/calibration versions with each scan.
- [ ] Archive tagged release, metadata and reproducibility outputs with DOI for publication.
