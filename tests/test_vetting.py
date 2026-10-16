"""Unit tests for astrophysical candidate vetting and systematics decontamination."""
import math
import pytest
from tess_anomaly_survey import vetting


def test_detect_simultaneous_clusters():
    """Verify cross-star simultaneous clustering identifies spacecraft systematics."""
    events = [
        # Simultaneous glitch at t=100.0 across 4 stars
        {'tic': 'star1', 'center_btjd': 100.0, 'sector': 1},
        {'tic': 'star2', 'center_btjd': 100.01, 'sector': 1},
        {'tic': 'star3', 'center_btjd': 100.02, 'sector': 1},
        {'tic': 'star4', 'center_btjd': 100.015, 'sector': 1},
        
        # Dual-star coincidental event at t=105.0 (< min_co_occurring_stars=3)
        {'tic': 'star1', 'center_btjd': 105.0, 'sector': 1},
        {'tic': 'star2', 'center_btjd': 105.01, 'sector': 1},
        
        # Isolated astrophysical events
        {'tic': 'star3', 'center_btjd': 110.5, 'sector': 1},
        {'tic': 'star4', 'center_btjd': 115.2, 'sector': 1},
    ]

    clusters, flagged = vetting.detect_simultaneous_clusters(events, time_tolerance_days=0.04, min_co_occurring_stars=3)
    
    assert len(clusters) == 1
    assert clusters[0]['star_count'] == 4
    assert abs(clusters[0]['mean_btjd'] - 100.01125) < 1e-4
    
    # Check flagged events
    assert ('star1', 100.0) in flagged
    assert ('star2', 100.01) in flagged
    assert ('star3', 100.02) in flagged
    assert ('star4', 100.015) in flagged
    
    # 105.0 (2 stars) and isolated events must NOT be flagged
    assert ('star1', 105.0) not in flagged
    assert ('star2', 105.01) not in flagged
    assert ('star3', 110.5) not in flagged
    assert ('star4', 115.2) not in flagged


def test_compute_star_vetting_metrics_clean_multi_dip():
    """Verify vetting metrics for clean multi-dip candidate."""
    cands = [
        {'id': 1, 'tic': 'clean_star', 'sector': 1, 'center_btjd': 10.0, 'score': 150.0, 'width_hr': 2.0, 'depth_ppm': 15000.0},
        {'id': 2, 'tic': 'clean_star', 'sector': 1, 'center_btjd': 15.0, 'score': 180.0, 'width_hr': 2.0, 'depth_ppm': 18000.0},
        {'id': 3, 'tic': 'clean_star', 'sector': 1, 'center_btjd': 20.0, 'score': 160.0, 'width_hr': 2.0, 'depth_ppm': 16000.0},
    ]
    flagged = set() # no systematics
    metrics = vetting.compute_star_vetting_metrics(cands, flagged)

    assert metrics['tic'] == 'clean_star'
    assert metrics['n_dips'] == 3
    assert metrics['n_clean_dips'] == 3
    assert metrics['systematic_fraction'] == 0.0
    assert metrics['max_clean_score'] == 180.0
    assert metrics['tier'] == 'Tier 2: Clean Multi-Dip Non-Periodic'
    assert metrics['priority_score'] > 0.0


def test_compute_star_vetting_metrics_systematic_dominated():
    """Verify stars dominated by spacecraft glitches are classified into Tier 5."""
    cands = [
        {'id': 1, 'tic': 'glitch_star', 'sector': 10, 'center_btjd': 100.0, 'score': 120.0, 'width_hr': 2.0, 'depth_ppm': 10000.0},
        {'id': 2, 'tic': 'glitch_star', 'sector': 10, 'center_btjd': 105.0, 'score': 110.0, 'width_hr': 2.0, 'depth_ppm': 9000.0},
        {'id': 3, 'tic': 'glitch_star', 'sector': 10, 'center_btjd': 110.0, 'score': 130.0, 'width_hr': 2.0, 'depth_ppm': 11000.0},
        {'id': 4, 'tic': 'glitch_star', 'sector': 10, 'center_btjd': 115.0, 'score': 105.0, 'width_hr': 2.0, 'depth_ppm': 8500.0},
    ]
    # All 4 dips match simultaneous spacecraft glitches
    flagged = {('glitch_star', 100.0), ('glitch_star', 105.0), ('glitch_star', 110.0), ('glitch_star', 115.0)}
    metrics = vetting.compute_star_vetting_metrics(cands, flagged)

    assert metrics['n_clean_dips'] == 0
    assert metrics['systematic_fraction'] == 1.0
    assert metrics['tier'] == 'Tier 5: Instrumental Systematics'
    assert metrics['priority_score'] == 0.0


def test_compute_star_vetting_metrics_unverified_periodic():
    """Verify unverified periodic systems receive Tier 4 and audit boost."""
    cands = [
        {'id': 1, 'tic': 'TIC_UNVERIFIED', 'sector': 4, 'center_btjd': 10.0, 'score': 200.0, 'width_hr': 2.0, 'depth_ppm': 40000.0},
        {'id': 2, 'tic': 'TIC_UNVERIFIED', 'sector': 4, 'center_btjd': 15.89, 'score': 210.0, 'width_hr': 2.0, 'depth_ppm': 41000.0},
        {'id': 3, 'tic': 'TIC_UNVERIFIED', 'sector': 4, 'center_btjd': 21.78, 'score': 195.0, 'width_hr': 2.0, 'depth_ppm': 39000.0},
    ]
    sys_info = {'is_periodic': True, 'period_days': 5.89, 'classification': 'periodic'}
    cat_match = {'catalog_identifier': 'Uncataloged', 'match_category': 'Unverified Identification'}
    
    metrics = vetting.compute_star_vetting_metrics(cands, set(), sys_info, cat_match)
    
    assert metrics['tier'] == 'Tier 4: Unverified Periodic System'
    assert metrics['is_periodic'] is True
    assert metrics['priority_score'] > 0.0


def test_rank_vetted_shortlist():
    """Verify shortlist produces balanced selection with target counts."""
    summaries = []
    # Create mock stars for different tiers
    for i in range(10):
        summaries.append({
            'tic': f'T1_{i}', 'tier': 'Tier 1: Extreme Non-Periodic Anomaly', 'priority_score': 15.0 - i
        })
    for i in range(10):
        summaries.append({
            'tic': f'T2_{i}', 'tier': 'Tier 2: Clean Multi-Dip Non-Periodic', 'priority_score': 10.0 - i
        })
    for i in range(8):
        summaries.append({
            'tic': f'T3_{i}', 'tier': 'Tier 3: Clean Single-Dip Candidate', 'priority_score': 5.0 - i
        })
    summaries.append({'tic': 'TIC_294261093', 'tier': 'Tier 4: Unverified Periodic System', 'priority_score': 8.5})
    summaries.append({'tic': 'TIC_4629065', 'tier': 'Tier 4: Unverified Periodic System', 'priority_score': 7.2})
    summaries.append({'tic': 'WASP-35', 'tier': 'Tier 5: Confirmed Catalog Object', 'priority_score': 0.0})

    shortlist = vetting.rank_vetted_shortlist(summaries, target_count=24)
    
    assert len(shortlist) == 24
    tics = [s['tic'] for s in shortlist]
    # Unverified periodic candidates must be included
    assert 'TIC_294261093' in tics
    assert 'TIC_4629065' in tics
    # Top Tier 1 and Tier 2 stars included
    assert 'T1_0' in tics
    assert 'T2_0' in tics
    # Confirmed catalog objects in Tier 5 excluded from discovery shortlist
    assert 'WASP-35' not in tics


def test_verify_event_photometry_genuine_transit():
    """Verify that a genuine negative transit dip passes verification."""
    import numpy as np
    t = np.linspace(100.0, 101.0, 720) # 2-minute cadence over 1 day
    g = np.ones_like(t, dtype=bool)
    f = np.zeros_like(t)
    
    # Inject 1.5% genuine dip at t=100.5 of width 2h
    half_w = (2.0 / 24.0) / 2.0
    m_dip = (t >= 100.5 - half_w) & (t <= 100.5 + half_w)
    f[m_dip] = -0.015
    
    res = vetting.verify_event_photometry(t, f, g, center_btjd=100.5, width_hr=2.0, mf_depth_ppm=15000.0)
    assert res['is_verified'] is True
    assert res['reason'] == 'genuine_dimming_confirmed'
    assert abs(res['true_depth_min_pct'] - 1.5) < 0.1
    assert res['has_flare'] is False


def test_verify_event_photometry_flare_shoulder_phantom():
    """Verify that a phantom dip created by an elevated flare shoulder is rejected."""
    import numpy as np
    t = np.linspace(100.0, 101.0, 720)
    g = np.ones_like(t, dtype=bool)
    f = np.zeros_like(t)
    
    # Candidate at t=100.5 has flat flux f ~ 0.0, but right shoulder (t=100.6) has a +10% flare
    m_flare = (t >= 100.58) & (t <= 100.65)
    f[m_flare] = 0.10
    
    # Matched filter would perceive an apparent dip because shoulder is elevated
    res = vetting.verify_event_photometry(t, f, g, center_btjd=100.5, width_hr=2.0, mf_depth_ppm=30000.0)
    assert res['is_verified'] is False
    assert res['reason'] in ('flare_shoulder_artifact', 'no_significant_negative_excursion', 'phantom_shoulder_elevation')
    assert res['has_flare'] is True
    assert res['true_depth_mean_pct'] <= 0.1


def test_verify_event_photometry_insufficient_data():
    """Verify graceful handling when data coverage is sparse."""
    import numpy as np
    t = np.array([100.0, 100.1, 100.2])
    f = np.zeros_like(t)
    g = np.ones_like(t, dtype=bool)
    
    res = vetting.verify_event_photometry(t, f, g, center_btjd=100.1)
    assert res['outcome'] == 'inconclusive'
    assert res['is_verified'] is False
    assert res['reason'] == 'insufficient_local_cadences'


def test_verify_event_photometry_single_cadence_outlier():
    """Verify that isolated single-cadence negative spikes are rejected as probable phantoms."""
    import numpy as np
    t = np.linspace(100.0, 101.0, 720) # 2-min cadence
    g = np.ones_like(t, dtype=bool)
    f = np.zeros_like(t)
    
    # Single isolated negative point at index 360
    f[360] = -0.05
    
    res = vetting.verify_event_photometry(t, f, g, center_btjd=t[360], width_hr=2.0, mf_depth_ppm=20000.0)
    assert res['outcome'] == 'probable_phantom'
    assert res['is_verified'] is False
    assert res['reason'] == 'single_cadence_negative_outlier'
    assert res['consecutive_depressed'] <= 1


def test_verify_event_photometry_broad_eclipse():
    """Verify that broad eclipses (tau > 2h) with depressed shoulders pass verification."""
    import numpy as np
    t = np.linspace(100.0, 101.0, 720)
    g = np.ones_like(t, dtype=bool)
    f = np.zeros_like(t)
    
    # 6-hour broad eclipse centered at 100.5, depth 20%
    half_dur = 3.0 / 24.0
    m_ecl = abs(t - 100.5) <= half_dur
    f[m_ecl] = -0.20
    
    res = vetting.verify_event_photometry(t, f, g, center_btjd=100.5, width_hr=2.0, mf_depth_ppm=40000.0)
    assert res['outcome'] == 'verified_dimming'
    assert res['is_verified'] is True
    assert res['true_depth_min_pct'] > 15.0
    assert res['consecutive_depressed'] >= 20


def test_verify_event_photometry_transit_near_flare():
    """Verify that a genuine transit occurring near a flare passes if multi-cadence coherence is present."""
    import numpy as np
    t = np.linspace(100.0, 101.0, 720)
    g = np.ones_like(t, dtype=bool)
    f = np.zeros_like(t)
    
    # Flare in left shoulder (t=100.4 to 100.44) peaking at +8%
    m_flare = (t >= 100.40) & (t <= 100.44)
    f[m_flare] = 0.08
    
    # Genuine 1.5% transit in core (t=100.48 to 100.52)
    m_dip = (t >= 100.48) & (t <= 100.52)
    f[m_dip] = -0.015
    
    res = vetting.verify_event_photometry(t, f, g, center_btjd=100.5, width_hr=2.0, mf_depth_ppm=15000.0)
    assert res['outcome'] == 'verified_dimming'
    assert res['is_verified'] is True
    assert res['consecutive_depressed'] >= 5


