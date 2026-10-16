# TESS Anomaly Survey — Public Demonstration Dataset

This directory contains a self-contained, lightweight demonstration dataset designed for offline exploration and testing of the **TESS Anomaly Survey Research Workspace** without requiring active internet connectivity or large archive downloads.

---

## 1. Dataset Contents

- **`sample_survey_review.sqlite`**: SQLite database containing representative anomaly detections across 5 test stars covering the 5 primary taxonomic classes:
  1. **TIC 39018208 (EBLM J0048-66, Sector 2)**: Benchmark low-mass eclipsing binary with 4 detected deep primary eclipses ($P \approx 2.4182\text{ d}$, depth $\approx 4,820\text{ ppm}$) and periodic orbit solution.
  2. **TIC 73834 (Sector 8)**: Candidate deep localized dimming system ($S_{\max} = 98.42$, depth $\approx 7,500\text{ ppm}$), classified as `verified_localized_dimming` (candidate transit / disk occultation).
  3. **TIC 60232 (Sector 1)**: Stellar pulsation variable rejected by the localized-dip verifier (`continuous_stellar_variability`, 17 correlated dips across the observation).
  4. **TIC 141944 (Sector 4)**: Flare overshoot phantom artifact rejected by the two-pass verifier (`flare_induced_apparent_dip`, shoulder elevation $\approx 8.4\%$).
  5. **TIC 999001 (Sector 1)**: Quiet control star (no detected anomalies; demonstrating zero-alert background handling).

- **`sample_data/`**: Bundled light curve files in compressed NumPy `.npz` format:
  - `lc_39018208_s02.npz`: Time, relative PDC flux, and quality flags covering Sector 2.
  - `lc_73834_s08.npz`: Time, relative PDC flux, and quality flags covering Sector 8.

---

## 2. Scientific Labeling & Synthetic Data Disclosure

> [!NOTE]
> **Synthetic Demonstration Light Curves:**
> To guarantee zero external internet dependencies and keep the repository download under 1 MB, the demonstration light curves in `sample_data/` are **synthetically generated baseline models** with injected parameters matching the candidate database records.
> 
> When running the workspace with real research targets (e.g. via `--cache-dir`), the survey pipeline downloads genuine, unmodified NASA TESS SPOC 2-minute cadence FITS files directly from the Mikulski Archive for Space Telescopes (MAST).

---

## 3. Launching the Demo Workspace

To launch the local research workspace using this demonstration dataset:

```bash
# From repository root
python -m tess_anomaly_survey review-ui --db examples/sample_survey_review.sqlite --root examples/sample_data
```

Then navigate to: **`http://127.0.0.1:8000/`** in any modern web browser.
