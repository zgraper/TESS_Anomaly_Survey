import numpy as np
import pytest
from tess_anomaly_survey import benchmark


def test_wilson_score_interval():
    """Verify Wilson score confidence interval boundaries and symmetry."""
    low, high = benchmark.wilson_score_interval(0, 100)
    assert low == 0.0
    assert 0.0 < high < 0.05

    low, high = benchmark.wilson_score_interval(100, 100)
    assert 0.95 < low < 1.0
    assert high == pytest.approx(1.0)

    low, high = benchmark.wilson_score_interval(50, 100)
    assert abs((low + high) / 2.0 - 0.5) < 0.01


def test_signal_profiles():
    """Verify box transit and asymmetric dimming profile generation."""
    times = np.linspace(0.0, 1.0, 720)  # 1 day of 2-min observations
    t0 = 0.5
    duration_hr = 2.0  # ~0.0833 days
    depth_ppm = 5000.0

    # Box profile
    p_box = benchmark.generate_signal_profile(times, t0, duration_hr, depth_ppm, shape="box")
    assert np.max(p_box) == pytest.approx(depth_ppm * 1e-6, rel=1e-3)
    # Outside duration should be zero
    outside = np.abs(times - t0) > (duration_hr / 24.0 / 2.0)
    assert np.all(p_box[outside] == 0.0)

    # Asymmetric profile
    p_asym = benchmark.generate_signal_profile(times, t0, duration_hr, depth_ppm, shape="asymmetric")
    assert np.max(p_asym) == pytest.approx(depth_ppm * 1e-6, rel=0.03)
    # Before t0 - ramp should be zero
    well_before = times < (t0 - 0.2 * (duration_hr / 24.0))
    assert np.all(p_asym[well_before] == 0.0)
    # Peak is near t0
    peak_idx = np.argmax(p_asym)
    assert abs(times[peak_idx] - t0) < 0.01
