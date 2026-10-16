import os
import csv
import glob
import pytest
import numpy as np

pd = pytest.importorskip("pandas")

CATALOG_A_PATH = 'research/phase6a/outputs/catalog_astrophysical_validation.csv' if os.path.isfile('research/phase6a/outputs/catalog_astrophysical_validation.csv') else 'scratch/catalog_astrophysical_validation.csv'
CATALOG_B_PATH = 'research/phase6a/outputs/catalog_unresolved_research.csv' if os.path.isfile('research/phase6a/outputs/catalog_unresolved_research.csv') else 'scratch/catalog_unresolved_research.csv'
PROVENANCE_PATH = 'research/phase6a/outputs/target_fits_metadata_provenance.csv' if os.path.isfile('research/phase6a/outputs/target_fits_metadata_provenance.csv') else 'scratch/target_fits_metadata_provenance.csv'

pytestmark = pytest.mark.skipif(
    not os.path.isfile(CATALOG_A_PATH) or not os.path.isfile(CATALOG_B_PATH),
    reason="Phase 6A research catalog artifacts not present in standalone software distribution"
)

def test_catalogs_exist_and_non_empty():
    assert os.path.isfile(CATALOG_A_PATH), f"Catalog A missing: {CATALOG_A_PATH}"
    assert os.path.isfile(CATALOG_B_PATH), f"Catalog B missing: {CATALOG_B_PATH}"
    
    df_a = pd.read_csv(CATALOG_A_PATH)
    df_b = pd.read_csv(CATALOG_B_PATH)
    
    assert len(df_a) == 8, f"Catalog A should have 8 rows, got {len(df_a)}"
    assert len(df_b) == 9, f"Catalog B should have 9 rows, got {len(df_b)}"

def test_catalog_identifiers_format():
    df_a = pd.read_csv(CATALOG_A_PATH)
    df_b = pd.read_csv(CATALOG_B_PATH)
    
    for df in [df_a, df_b]:
        for tic in df['tic']:
            assert tic.startswith('TIC:'), f"Invalid TIC format: {tic}"
            raw_id = tic.replace('TIC:', '')
            assert raw_id.isdigit(), f"TIC ID not numeric: {tic}"

def test_sectors_match_fits_provenance():
    df_a = pd.read_csv(CATALOG_A_PATH)
    df_b = pd.read_csv(CATALOG_B_PATH)
    df_prov = pd.read_csv(PROVENANCE_PATH).drop_duplicates(subset=['tic'])
    prov_map = {r['tic']: str(r['header_sector']) for r in df_prov.to_dict(orient='records')}
    
    # Specific known ground truth sectors:
    ground_truth_sectors = {
        'TIC:71298511': '9',
        'TIC:247585034': '6',
        'TIC:229810745': '2',
        'TIC:4629065': '4',
        'TIC:197908268': '2',
        'TIC:279955276': '2',
        'TIC:262958558': '1',
        'TIC:234497727': '1',
        'TIC:4207747': '6',
        'TIC:281636513': '8',
        'TIC:341419283': '8',
        'TIC:30315127': '8',
        'TIC:4847760': '8',
        'TIC:53748370': '7',
        'TIC:277696601': '8',
        'TIC:294261093': '10'
    }
    
    # Check Catalog A
    for idx, row in df_a.iterrows():
        tic = row['tic']
        sec = str(row['sectors'])
        expected_sec = ground_truth_sectors.get(tic)
        assert sec == expected_sec, f"Sector mismatch for {tic} in Catalog A: expected {expected_sec}, got {sec}"
        assert sec == prov_map.get(tic), f"Provenance sector mismatch for {tic}: {sec} vs {prov_map.get(tic)}"
        
    # Check Catalog B
    for idx, row in df_b.iterrows():
        tic = row['tic']
        sec = str(row['sectors'])
        if tic == 'TIC:141126559':
            assert sec in ['5;6', '5,6', '5'], f"Unexpected multi-sector format for {tic}: {sec}"
        else:
            expected_sec = ground_truth_sectors.get(tic)
            assert sec == expected_sec, f"Sector mismatch for {tic} in Catalog B: expected {expected_sec}, got {sec}"

def test_stellar_parameters_match_fits_headers():
    df_a = pd.read_csv(CATALOG_A_PATH)
    df_b = pd.read_csv(CATALOG_B_PATH)
    df_prov = pd.read_csv(PROVENANCE_PATH).drop_duplicates(subset=['tic'])
    prov_map = {r['tic']: r for r in df_prov.to_dict(orient='records')}
    
    for df in [df_a, df_b]:
        for idx, row in df.iterrows():
            tic = row['tic']
            p = prov_map[tic]
            
            # Tmag check (within 0.005)
            cat_tmag = float(row['tmag'])
            fits_tmag = float(p['tmag'])
            assert abs(cat_tmag - fits_tmag) < 0.005, f"Tmag mismatch for {tic}: {cat_tmag} vs {fits_tmag}"
            
            # Teff check
            if pd.notnull(p['teff']) and str(p['teff']).strip() != '' and str(p['teff']).lower() != 'nan':
                fits_teff = float(p['teff'])
                assert pd.notnull(row['teff']), f"Teff unexpectedly empty in catalog for {tic}"
                cat_teff = float(row['teff'])
                assert abs(cat_teff - fits_teff) < 1.0, f"Teff mismatch for {tic}: {cat_teff} vs {fits_teff}"
            else:
                assert pd.isna(row['teff']) or str(row['teff']).strip() == '', f"Teff should be empty for {tic}"

def test_literature_crossmatch_records():
    df_a = pd.read_csv(CATALOG_A_PATH)
    df_b = pd.read_csv(CATALOG_B_PATH)
    
    row_yhyi = df_a[df_a['tic'] == 'TIC:262958558'].iloc[0]
    assert 'Y Hyi' in row_yhyi['literature_crossmatch'], "TIC 262958558 should be matched to Y Hyi"
    assert 'Mkrtichian' in row_yhyi['literature_crossmatch'], "TIC 262958558 should cite Mkrtichian et al."
    assert abs(float(row_yhyi['period_days_lit']) - 3.53597) < 0.01, "Period should match published 3.53597d"
    
    row_sttau = df_a[df_a['tic'] == 'TIC:247585034'].iloc[0]
    assert 'ST Tau' in row_sttau['literature_crossmatch'], "TIC 247585034 should be matched to ST Tau"
    assert abs(float(row_sttau['period_days_lit']) - 4.0343) < 0.01, "ST Tau lit period should be ~4.034d"
    assert abs(float(row_sttau['period_days_survey']) - 4.06757) < 0.01, "ST Tau survey period should be ~4.068d"
    
    row_420 = df_b[df_b['tic'] == 'TIC:4207747'].iloc[0]
    assert 'HD 35317' in row_420['literature_crossmatch'], "TIC 4207747 should be matched to HD 35317"
    assert 'Howard' in row_420['literature_crossmatch'], "TIC 4207747 should cite Howard et al. (2022)"
    assert 'RNAAS 6 96' in row_420['literature_crossmatch'], "TIC 4207747 should cite RNAAS 6 96"

def test_tic_4207747_transit_duration_bounds():
    df_b = pd.read_csv(CATALOG_B_PATH)
    row_420 = df_b[df_b['tic'] == 'TIC:4207747'].iloc[0]
    
    t14 = float(row_420['transit_duration_t14_hr'])
    fwhm = float(row_420['fwhm_hr'])
    
    # Must be physical: T14 between 2.5 and 6 hours, FWHM between 1.5 and 3 hours
    assert 2.5 <= t14 <= 6.0, f"Unphysical T14 duration for TIC 4207747: {t14} hr"
    assert 1.5 <= fwhm <= 3.0, f"Unphysical FWHM duration for TIC 4207747: {fwhm} hr"
    assert fwhm < t14, f"FWHM should be smaller than T14: {fwhm} >= {t14}"
