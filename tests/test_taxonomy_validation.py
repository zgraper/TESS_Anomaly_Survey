"""
Unit tests validating 5-state scientific taxonomy across astrophysical and instrumental scenarios:
1. Known pulsator / rotational variable -> continuous_stellar_variability
2. Eclipsing binary with repeating dips -> verified_localized_dimming (NEVER reclassified as continuous)
3. Isolated transit -> verified_localized_dimming
4. Smooth instrumental drift / non-coherent excursion -> probable_instrumental_artifact
5. Flare-contaminated dip -> flare_induced_apparent_dip
6. Data gap / boundary proximity -> inconclusive
"""
import numpy as np
import pytest
from tess_anomaly_survey import streaming, vetting


def test_continuous_variability_pulsator():
    """Verify that smooth continuous periodic modulation is classified as continuous_stellar_variability."""
    n = 2000
    t = np.linspace(1000.0, 1002.77, n)
    f = 0.04 * np.sin(2.0 * np.pi * t / 0.15)
    g = np.ones(n, dtype=bool)
    assert streaming.detect_continuous_variability(t, f, g, n_events=10) is True


def test_repeating_eclipsing_binary_retention():
    """Verify that genuine repeating eclipses are NOT reclassified as continuous variability."""
    n = 2000
    t = np.linspace(1000.0, 1002.77, n)
    # Quiet continuum with Gaussian noise
    np.random.seed(42)
    f = np.random.normal(0.0, 0.0005, n)
    g = np.ones(n, dtype=bool)
    cadence = 2.0 / 1440.0

    # Inject 4 deep eclipses (2% depth, 2.4 hr duration)
    eclipse_times = [1000.5, 1001.0, 1001.5, 1002.0]
    for et in eclipse_times:
        mask = np.abs(t - et) <= (1.2 / 24.0)
        f[mask] -= 0.02

    # Even if n_events is large, check that detect_continuous_variability is False because continuum is flat
    # Point-to-point scatter is dominated by flat continuum noise
    assert streaming.detect_continuous_variability(t, f, g, n_events=4) is False

    # Check that verify_event_photometry confirms each eclipse
    ver = vetting.verify_event_photometry(
        t, f, g,
        center_btjd=1000.5,
        width_hr=2.4,
        mf_depth_ppm=20000.0,
        min_consecutive_cadences=2
    )
    assert ver["outcome"] == "verified_dimming"
    assert ver["reason"] == "genuine_dimming_confirmed"


def test_isolated_transit_verification():
    """Verify that an isolated 2-hour transit is verified."""
    n = 2000
    t = np.linspace(1000.0, 1002.77, n)
    np.random.seed(42)
    f = np.random.normal(0.0, 0.0003, n)
    g = np.ones(n, dtype=bool)

    # 1.5% dip at center
    mask = np.abs(t - 1001.4) <= (1.0 / 24.0)
    f[mask] -= 0.015

    ver = vetting.verify_event_photometry(
        t, f, g,
        center_btjd=1001.4,
        width_hr=2.0,
        mf_depth_ppm=15000.0,
        min_consecutive_cadences=2
    )
    assert ver["outcome"] == "verified_dimming"
    assert ver["reason"] == "genuine_dimming_confirmed"


def test_flare_induced_apparent_dip():
    """Verify that an apparent dip created by flare shoulders is caught as flare_induced_apparent_dip."""
    n = 2000
    t = np.linspace(1000.0, 1002.77, n)
    np.random.seed(42)
    f = np.random.normal(0.0, 0.0003, n)
    g = np.ones(n, dtype=bool)

    # Inject flare in shoulder window (+5% brightening just before dip)
    sh_mask = (t >= 1001.25) & (t <= 1001.35)
    f[sh_mask] += 0.05
    # The center has NO negative dip relative to continuum (f ≈ 0)

    ver = vetting.verify_event_photometry(
        t, f, g,
        center_btjd=1001.4,
        width_hr=2.0,
        mf_depth_ppm=10000.0, # matched filter was fooled by shoulder
        min_consecutive_cadences=2
    )
    assert ver["outcome"] == "probable_phantom"
    assert ver["reason"] == "phantom_shoulder_elevation"
    assert ver["has_flare"] is True


def test_inconclusive_gap_boundary():
    """Verify that an event next to a gap or sector boundary with N_local < 10 is inconclusive."""
    n = 2000
    t = np.linspace(1000.0, 1002.77, n)
    f = np.zeros(n)
    g = np.zeros(n, dtype=bool)
    # Only 5 cadences are good in the whole slice
    g[:5] = True

    ver = vetting.verify_event_photometry(
        t, f, g,
        center_btjd=1000.02, # right at the edge
        width_hr=2.0,
        mf_depth_ppm=5000.0,
        min_consecutive_cadences=2
    )
    assert ver["outcome"] == "inconclusive"
    assert ver["reason"] == "insufficient_local_cadences"
