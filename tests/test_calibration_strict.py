import math
import numpy as np
import pytest
from tess_anomaly_survey import core


def test_sample_size_gating():
    """Verify build_calibration refuses insufficient calibration sample sizes unless explicitly overridden."""
    from unittest.mock import patch
    def mock_iter_windows(path, **kwargs):
        for i in range(5):
            yield {'start': i, 'time': np.linspace(0, 1, 512), 'flux': np.ones(512), 'good': np.ones(512, bool), 'cadence_days': 2/1440}

    def mock_score_window(w, method, cache):
        return [], 0.001, 10.0

    fake_rows = [{'tic': f'TIC_{i}', 'processed_file': f'fake_{i}.npz'} for i in range(5)]
    settings = {'length': 512, 'stride': 256, 'min_good': 0.97, 'max_windows': 0}

    # Should raise RuntimeError without allow_small due to insufficient sample size
    with patch('tess_anomaly_survey.core.iter_windows', mock_iter_windows):
        with patch('tess_anomaly_survey.core.score_window', mock_score_window):
            with pytest.raises(RuntimeError) as excinfo:
                core.build_calibration(
                    fake_rows, root=".", method="fixed_box_2h", settings=settings,
                    target=0.01, grouping="pooled", min_stars=1000, min_windows=10000,
                    allow_small=False
                )
    assert "Calibration insufficient" in str(excinfo.value) or "inadequate" in str(excinfo.value)


def test_compute_calibration_uncertainty_star_bootstrap():
    """Verify star-level block bootstrap computes SE and 95% CI preserving star clusters."""
    # Create 20 stars, each with 10 correlated windows
    per_tic = {}
    rng = np.random.default_rng(42)
    for i in range(20):
        star_baseline = rng.normal(50.0, 10.0)
        # Windows on the same star share baseline noise
        window_scores = star_baseline + rng.normal(0.0, 3.0, size=10)
        per_tic[f"TIC_{i}"] = [{'max_score': float(s), 'noise': 0.001} for s in window_scores]

    cuts = [0.0005, 0.0015]
    unc = core.compute_calibration_uncertainty(per_tic, target=0.05, grouping='pooled', cutpoints=cuts, n_boot=100, seed=42)

    assert unc['star_level_block_bootstrap'] is True
    assert unc['bootstrap_replicates'] == 100
    pooled = unc['pooled']
    assert pooled['std_err'] > 0.0
    assert pooled['ci_95'][0] < pooled['ci_95'][1]
    assert pooled['relative_stderr'] > 0.0


def test_calibration_tail_survival_curve(tmp_path):
    """Test tail survival diagnostic generation under allow_small_calibration."""
    from unittest.mock import patch

    # Mock iter_windows and score_window to simulate 500 windows across 20 stars
    def mock_iter_windows(path, **kwargs):
        for i in range(25):
            yield {'start': i, 'time': np.linspace(0, 1, 512), 'flux': np.ones(512), 'good': np.ones(512, bool), 'cadence_days': 2/1440}

    scores_iter = iter(np.linspace(10, 100, 500))
    def mock_score_window(w, method, cache):
        return [], 0.001, next(scores_iter)

    fake_rows = [{'tic': f'TIC_{i}', 'processed_file': f'fake_{i}.npz'} for i in range(20)]
    settings = {'length': 512, 'stride': 256, 'min_good': 0.97, 'max_windows': 0}

    with patch('tess_anomaly_survey.core.iter_windows', mock_iter_windows):
        with patch('tess_anomaly_survey.core.score_window', mock_score_window):
            cal = core.build_calibration(
                fake_rows, root=".", method="fixed_box_2h", settings=settings,
                target=0.01, grouping="pooled", min_stars=10, min_windows=100,
                allow_small=True, bootstrap=50
            )

    assert 'thresholds' in cal
    assert 'pooled' in cal['thresholds']
    assert 'diagnostics' in cal
    diag = cal['diagnostics']
    assert 'tail_survival_curve' in diag
    curve = diag['tail_survival_curve']
    assert len(curve) > 0
    # Thresholds should strictly increase as quantile increases
    quants = [pt['quantile'] for pt in curve]
    threshs = [pt['score_threshold'] for pt in curve]
    assert quants == sorted(quants)
    assert threshs == sorted(threshs)
