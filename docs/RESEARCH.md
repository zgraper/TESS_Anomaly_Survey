# Research integrity and frozen sources

The original B1 statistical detectors, reproduced **byte-for-byte** in the `src/` directory, have these SHA-256 checksums:

| Source | SHA-256 |
|---|---|
| `src/experiment_b0.py` | `9c8afe63ec2763b9db57c828b52b0e7d262bb30b203f86e559836dd5beb89f00` |
| `src/experiment_b1.py` | `c5214809a053e4705c5d86f29703e7682128e3db119494ed59f6c9a191e557c8` |

The survey refuses to run when these files are silently edited. Private research artifacts are not shipped in this public repository. See [`calibration_provenance.json`](../src/tess_anomaly_survey/calibration_provenance.json) for the sanitized calibration identity and limitations. Historical statistics below are reported context, not independently reproduced by this software test suite.

A **separate survey calibration artifact** is created for the new scanning window population. It records the detector, source hashes, target, window settings, cutpoints, per-group empirical thresholds, number of stars/windows and the exact TIC set excluded from independent validation. It does not modify B1.2 calibration or retroactively improve the published B2.2 scores.

The user-reported B2.2 final external findings were 144/600 correct-time recoveries (24.0%, fixed 2-hour box with noise-tertile thresholds) and 37/600 (6.2%) presumed-background control windows containing at least one alert. The more complex multiscale and planet-shaped detectors each recovered 103/600 (17.2%). These are **injection-recovery statistics** under the frozen B2.2 protocol, not a scientifically vetted sample of actual unexplained celestial events.

The external protocol was amended before scoring to represent Sectors 11, 12 and 13 equally (200 stars each), after the original QC sector composition had been inspected. This post-QC, pre-scoring amendment should remain disclosed in any manuscript.

A complete reproducibility package for a paper additionally needs versioned input data, environment specifications, event-level results and the exact paper figure-generation scripts. These large assets have not been recreated here.
