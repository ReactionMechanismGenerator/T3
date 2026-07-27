#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_cache
"""

import os

import pytest

from arc.common import save_yaml_file

from t3.pdep.cache import (SA_CACHE_METADATA_FILE_NAME,
                          hash_file,
                          max_abs_ts_coefficient,
                          sa_cache_metadata_path,
                          validate_sa_cache,
                          write_sa_cache_metadata,
                          )
from t3.pdep.selector import CACHE_STATUS_CACHED_REJECTED, CACHE_STATUS_CACHED_VALID

# A TS coefficient safely above the default absolute floor (~1.2e-7 mol/J for the default
# perturbation and min_delta_ln_k), standing in for real Arkane SA output that carries signal.
REALISTIC_TS_COEFFICIENT = 0.05


def _write(path, content):
    """Write a small text file at ``path``."""
    with open(path, 'w') as f:
        f.write(content)


def _write_realistic_sa_yaml(path):
    """Write a minimal but realistic Arkane SA YAML, with a TS coefficient above the absolute floor."""
    save_yaml_file(path=path, content={
        'structures': dict(),
        'reactant1 + reactant2 <=> product': {
            (1000.0, 'K', 1.0, 'bar'): {
                '(TS) TS1': REALISTIC_TS_COEFFICIENT,
            },
        },
    })


# --- 14. sa_cache_metadata_path -----------------------------------------------------------------

def test_sa_cache_metadata_path_sits_next_to_sa_yaml(tmp_path):
    """Test that the sidecar path is in the same directory as the SA YAML."""
    sa_path = str(tmp_path / 'sensitivity' / 'sa_coefficients.yml')
    metadata_path = sa_cache_metadata_path(sa_path=sa_path)
    assert os.path.dirname(metadata_path) == os.path.dirname(sa_path)
    assert os.path.basename(metadata_path) == SA_CACHE_METADATA_FILE_NAME


# --- 15. Round trip -------------------------------------------------------------------------------

def test_round_trip_write_then_validate_is_cached_valid(tmp_path):
    """Test that a freshly written sidecar validates as ``'cached_valid'`` with no warnings."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path, method='MSC')
    assert status == CACHE_STATUS_CACHED_VALID
    assert warnings == []


# --- 16. Missing sidecar ---------------------------------------------------------------------------

def test_missing_sidecar_is_rejected(tmp_path):
    """Test that an SA YAML with no sidecar (P1 case) is never trusted."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write(sa_path, 'some: sa\n')

    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path)
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert warnings


# --- 17. Network file changed after sidecar written ------------------------------------------------

def test_network_file_changed_after_write_is_rejected(tmp_path):
    """Test that a network file modified after the sidecar was written invalidates the cache."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write(sa_path, 'some: sa\n')

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    _write(network_path, 'reaction(label="reaction1_modified")\n')

    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path)
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert warnings


# --- 18. Sidecar records a different selector_version -----------------------------------------------

def test_sidecar_with_mismatched_selector_version_is_rejected(tmp_path):
    """Test that a sidecar recording a different ``selector_version`` is rejected."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write(sa_path, 'some: sa\n')

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    metadata_path = sa_cache_metadata_path(sa_path=sa_path)
    from arc.common import read_yaml_file, save_yaml_file
    metadata = read_yaml_file(path=metadata_path)
    metadata['selector_version'] = metadata['selector_version'] + 1
    save_yaml_file(path=metadata_path, content=metadata)

    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path)
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert warnings


# --- 19. Sidecar records a mismatched perturbation ---------------------------------------------------

def test_sidecar_with_mismatched_perturbation_is_rejected(tmp_path):
    """Test that a sidecar recording a mismatched ``perturbation`` is rejected."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write(sa_path, 'some: sa\n')

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC',
                            perturbation=8368.0)

    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path, perturbation=99999.0)
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert warnings


# --- 20. Missing SA YAML entirely ---------------------------------------------------------------------

def test_missing_sa_yaml_is_rejected(tmp_path):
    """Test that a missing SA YAML altogether is rejected."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')  # never written

    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path)
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert warnings


# --- 21. Malformed sidecar --------------------------------------------------------------------------

def test_malformed_sidecar_is_rejected_without_raising(tmp_path):
    """Test that an unparseable sidecar is rejected, and does not raise."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write(sa_path, 'some: sa\n')
    metadata_path = sa_cache_metadata_path(sa_path=sa_path)
    _write(metadata_path, '{not: valid: yaml: [')

    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path)
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert warnings


# --- max_abs_ts_coefficient helper (FIX1) ------------------------------------------------------------

def test_max_abs_ts_coefficient_finds_largest_across_reactions_and_conditions():
    """Test that the scan finds the largest abs TS coefficient across all reactions/conditions."""
    sa_dict = {
        'structures': {'A': 'adjlist'},  # must be skipped, not treated as a reaction key
        'A + B <=> C': {
            (1000.0, 'K', 1.0, 'bar'): {'(TS) TS1': -0.01, 'well1': 0.9},
            (1200.0, 'K', 1.0, 'bar'): {'(TS) TS1': 0.2},
        },
        'A + B <=> D': {
            (1000.0, 'K', 1.0, 'bar'): {'(TS) TS2': 0.05},
        },
    }
    assert max_abs_ts_coefficient(sa_dict) == 0.2


def test_max_abs_ts_coefficient_ignores_non_finite_values():
    """Test that non-finite TS coefficients (inf/nan) are ignored, not treated as the max."""
    sa_dict = {
        'A + B <=> C': {
            (1000.0, 'K', 1.0, 'bar'): {'(TS) TS1': float('inf'), '(TS) TS2': float('nan'), '(TS) TS3': 0.01},
        },
    }
    assert max_abs_ts_coefficient(sa_dict) == 0.01


def test_max_abs_ts_coefficient_returns_none_when_no_ts_entries():
    """Test that the scan returns None when there is no TS-prefixed coefficient anywhere."""
    sa_dict = {'A + B <=> C': {(1000.0, 'K', 1.0, 'bar'): {'well1': 0.9}}}
    assert max_abs_ts_coefficient(sa_dict) is None


def test_max_abs_ts_coefficient_returns_none_for_non_dict_input():
    """Test that the scan tolerates a non-dict sa_dict, returning None rather than raising."""
    assert max_abs_ts_coefficient('not a dict') is None
    assert max_abs_ts_coefficient(None) is None


# --- FIX1: cache rejected when TS rows carry no signal ------------------------------------------------

def test_cache_with_only_structural_zero_ts_coefficients_is_rejected(tmp_path):
    """Test that a cache whose TS rows are all below the absolute floor is rejected (P1)."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    save_yaml_file(path=sa_path, content={
        'structures': dict(),
        'reactant1 + reactant2 <=> product': {
            (1000.0, 'K', 1.0, 'bar'): {'(TS) TS1': 1e-18},  # a structural-zero denormal, below the floor
        },
    })

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path)
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert any('transition state' in w.lower() for w in warnings)


def test_cache_with_no_ts_coefficient_recorded_is_rejected(tmp_path):
    """Test that a cache with no recorded max_abs_ts_coefficient at all is rejected (P1)."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    save_yaml_file(path=sa_path, content={'structures': dict(), 'reactant1 + reactant2 <=> product': {}})

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path)
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert any('transition state' in w.lower() for w in warnings)


# --- FIX2: cache bound to the specific SA YAML and method ---------------------------------------------

def test_cache_rejected_when_sa_yaml_content_changed_since_sidecar_written(tmp_path):
    """Test that editing the SA YAML after the sidecar was written invalidates the cache."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    # Rewrite the SA YAML with different (but still realistic) content, without updating the sidecar.
    save_yaml_file(path=sa_path, content={
        'structures': dict(),
        'reactant1 + reactant2 <=> product': {
            (1000.0, 'K', 1.0, 'bar'): {'(TS) TS1': 0.09},
        },
    })

    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path)
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert any('changed' in w.lower() for w in warnings)


def test_cache_rejected_when_requested_method_disagrees_with_recorded_method(tmp_path):
    """Test that a cache generated with a different ME method than requested is rejected."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path, method='CSE')
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert any('method' in w.lower() for w in warnings)


# --- 21c. Threshold validation happens before any early return ---------------------------------------

def test_validate_sa_cache_raises_for_bad_perturbation_even_when_sa_yaml_is_missing(tmp_path):
    """Test that ``validate_sa_cache`` raises ``ValueError`` for a non-finite/non-positive
    ``perturbation`` before it ever reaches the early 'no cached sensitivity output' return -- a bad
    threshold is a caller bug that should surface immediately, not be masked by whichever early-return
    branch a particular call happens to hit first (here, a simply-missing SA file)."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')  # never written

    with pytest.raises(ValueError):
        validate_sa_cache(sa_path=sa_path, network_path=network_path, perturbation=0.0)


def test_validate_sa_cache_raises_for_bad_min_delta_ln_k_even_when_sidecar_is_missing(tmp_path):
    """Test that ``validate_sa_cache`` raises ``ValueError`` for a non-finite/non-positive
    ``min_delta_ln_k`` before it ever reaches the early 'no T3 metadata sidecar' return."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write(sa_path, 'some: sa\n')  # sidecar deliberately not written

    with pytest.raises(ValueError):
        validate_sa_cache(sa_path=sa_path, network_path=network_path, min_delta_ln_k=-1.0)


# --- 22. hash_file stability --------------------------------------------------------------------------

def test_hash_file_stable_for_identical_content_and_differs_for_changed_content(tmp_path):
    """Test that ``hash_file`` is stable for identical content and differs when content changes."""
    path_a = str(tmp_path / 'a.py')
    path_b = str(tmp_path / 'b.py')
    _write(path_a, 'identical content\n')
    _write(path_b, 'identical content\n')
    assert hash_file(path=path_a) == hash_file(path=path_b)

    _write(path_b, 'different content\n')
    assert hash_file(path=path_a) != hash_file(path=path_b)
