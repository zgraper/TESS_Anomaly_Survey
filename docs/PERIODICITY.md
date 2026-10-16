# Periodicity Analysis & Repeating Astronomical Event Classification

**Module:** `tess_anomaly_survey.periodicity`  
**Version:** v0.2.0 (Phase 3)

---

## 1. Scientific Overview & Assumptions

In the TESS Anomaly Survey, anomalous stellar dimming detections can represent:
1. **Periodic Stellar Variability or Transiting Exoplanets:** Events that repeat at regular orbital period $P$ from reference epoch $T_0$:
   $$t_k \approx T_0 + k \cdot P, \quad k \in \mathbb{Z}$$
2. **Singular / Non-Periodic Events ("Mono-Transits"):** Long-period exoplanets ($P > 27\,\text{days}$), irregular circumstellar dust dips (e.g. Boyajian's Star analogs), or stellar flares with overshoot.
3. **Spacecraft Systematics:** Momentum dump thruster firing cycles (typically recurring at $\sim 2.50\,\text{days}$ on TESS), Earth/Moon scattered light, or reaction wheel desaturations.

### Strict Scientific Constraints
- **Minimum Event Rule ($N \ge 3$):** A star is **never** claimed as periodic without at least three distinct, credible candidate events ($N \ge 3$) populating distinct orbital epochs. A star with 1 or 2 detections is strictly flagged as `is_periodic = False` (singular or inconclusive).
- **Ephemeris Timing Tolerance:** All candidate events must fit the linear ephemeris within a tolerance $\text{tol} = \max(2 \times \text{cadence}, 0.15 \times \bar{w})$, where $\bar{w}$ is the mean transit duration.
- **Harmonic Alias Resolution:** Subharmonic period trial overtones ($P/2, P/3, \dots$) are resolved using integer step Greatest Common Divisor (GCD) clustering to extract the true fundamental orbital period.
- **Eclipsing Binary Discrimination ($Z_{\text{odd/even}}$):** If alternating dips exhibit a statistically significant depth disparity ($Z_{\text{odd/even}} > 3.0\,\sigma$), the detector doubles the trial period to the true binary orbital period $P_{\text{binary}} = 2 \times P_{\text{trial}}$ to account for primary and secondary eclipses.
- **Systematics Warning:** If the period falls within $2.50 \pm 0.05\,\text{days}$ or $1.25 \pm 0.04\,\text{days}$, an automated warning is attached flagging potential spacecraft momentum dump contamination.

---

## 2. CLI Usage

### Running Periodicity Analysis

Analyze all candidate events in an anomaly library:

```powershell
tess-anomaly analyze-periods --db .\anomaly_library.sqlite
```

#### Optional Arguments
- `--run-id <id>`: Filter candidate events to a specific survey run.
- `--min-events <int>`: Minimum candidate events per star required to run periodicity analysis (default: 3; must be $\ge 3$).
- `--max-period <float>`: Maximum search period in days (default: 20.0 d).
- `--min-period <float>`: Minimum search period in days (default: 0.2 d, corresponding to stellar Roche breakup limit).

Example output:
```text
TIC 10001 (run 8a1f): PERIODIC P=3.5241 d, T0=1325.4120, N=4, SNR=28.4 Candidate periodic system
Periodicity analysis complete: 42 stars analyzed (>= 3 events), 8 periodic candidate systems detected.
```

---

## 3. Database Schema

Periodicity detections are stored transactionally in `candidate_systems`:

```sql
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
    created_utc TEXT NOT NULL,
    FOREIGN KEY(run_id) REFERENCES runs(id)
);
CREATE INDEX IF NOT EXISTS idx_systems_tic ON candidate_systems(run_id, tic);
```

Each candidate in `candidates` is linked to its system group via `system_id INTEGER REFERENCES candidate_systems(id)`.

---

## 4. Web Review Interface Diagnostics

When inspecting a periodic candidate in the review interface (`tess-anomaly review-ui --db ... --root ...`):
1. **Periodic Badges:** In the candidate list, periodic candidates display purple badges showing their fitted period, e.g. `P=3.52d`.
2. **Ephemeris Card:** Displays fitted orbital period $P$, reference epoch $T_0$, folded SNR, and odd/even disparity $Z_{\text{odd/even}}$.
3. **Phase-Folded Diagnostic Plot:**
   - Cadences folded by phase: $\phi(t) = ((t - T_0) / P) \pmod 1 - 0.5 \in [-0.5, 0.5]$.
   - Even (primary) transits displayed as blue points.
   - Odd (secondary) transits displayed as orange points.
   - Robust median curve overlaid in red to verify transit morphology.
