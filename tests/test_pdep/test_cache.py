#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_cache
"""

import os

import pytest

from arc.common import read_yaml_file, save_yaml_file

from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.cache import (SA_CACHE_METADATA_FILE_NAME,
                          hash_file,
                          max_abs_ts_coefficient,
                          read_arkane_log_rmg_py_commit,
                          read_t_grid_clamp_record,
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
    # Pin the REASON, not just the verdict: 'rejected' is reachable by seven other branches.
    assert any('no T3 metadata sidecar' in warning for warning in warnings), warnings


# --- 17. Network file changed after sidecar written ------------------------------------------------

def test_network_file_changed_after_write_is_rejected(tmp_path):
    """Test that a network file modified after the sidecar was written invalidates the cache."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    # A REALISTIC SA YAML, not a stub: this keeps the test meaningful in spirit even though
    # ``validate_sa_cache`` no longer gates on TS signal at all (that judgment moved to
    # ``t3.pdep.selector.select_from_sa_dict`` -- see FIX1 in cache.py). Using a realistic YAML here
    # still isolates this test to the network_file_hash check alone.
    _write_realistic_sa_yaml(sa_path)

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    _write(network_path, 'reaction(label="reaction1_modified")\n')

    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path, method='MSC')
    assert status == CACHE_STATUS_CACHED_REJECTED
    # Pin the REASON, not just the verdict: 'rejected' is reachable by seven other branches. Match on
    # both words together since 'changed' alone also appears in the (different) sa_file_hash message.
    assert any('network file' in warning.lower() and 'changed' in warning.lower()
              for warning in warnings), warnings


# --- 17b. Network file missing entirely -------------------------------------------------------------

def test_network_file_missing_is_rejected(tmp_path):
    """Test that a network file removed after the sidecar was written invalidates the cache.

    This is distinct from test_network_file_changed_after_write_is_rejected above: here the file is
    gone entirely (hash_file() cannot even be computed), not merely modified.
    """
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    os.remove(network_path)

    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path, method='MSC')
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert any('network file' in warning.lower() and 'missing' in warning.lower()
              for warning in warnings), warnings


def test_an_unreadable_network_file_is_rejected_rather_than_raising(tmp_path, monkeypatch):
    """Test that an OSError from hashing the network file is a cache rejection, not an exception.

    ``os.path.isfile`` answers "does a file exist here", not "can this process read it", so the
    guard above it does not cover a present-but-unreadable file. Letting the ``OSError`` escape
    would break this function's documented contract of returning a status, and at the ``t3/main.py``
    call site it is not caught by the surrounding ``(OSError, ValueError)`` handler either -- it
    reaches the outer ``except Exception``, which records an ``internal_error`` and re-raises,
    ending the whole campaign over one unreadable file.
    """
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)
    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')

    def refuse_the_network_file(path):
        if os.path.realpath(path) == os.path.realpath(network_path):
            raise PermissionError(13, 'Permission denied', path)
        return hash_file(path)

    monkeypatch.setattr('t3.pdep.cache.hash_file', refuse_the_network_file)
    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path, method='MSC')
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert any('network file' in warning.lower() and 'could not be read' in warning.lower()
               for warning in warnings), warnings


def test_an_unreadable_sa_file_is_rejected_rather_than_raising(tmp_path, monkeypatch):
    """Test that an OSError from hashing the SA YAML is a cache rejection, not an exception.

    The same hole as the network file above, at the second of the two hash gates: the SA output is
    confirmed to exist at the top of the function but is not read until much later.
    """
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)
    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')

    def refuse_the_sa_file(path):
        if os.path.realpath(path) == os.path.realpath(sa_path):
            raise PermissionError(13, 'Permission denied', path)
        return hash_file(path)

    monkeypatch.setattr('t3.pdep.cache.hash_file', refuse_the_sa_file)
    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path, method='MSC')
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert any('sensitivity output' in warning.lower() and 'could not be read' in warning.lower()
               for warning in warnings), warnings


@pytest.mark.skipif(hasattr(os, 'geteuid') and os.geteuid() == 0,
                    reason='root bypasses file permissions, so the unreadable file would be readable')
def test_a_genuinely_chmodded_network_file_is_rejected_rather_than_raising(tmp_path):
    """Test the same fail-closed behaviour against a real unreadable file, not a patched hash.

    The two monkeypatched tests above pin the contract; this one exists so that contract is known
    to hold for the filesystem condition it was written for, rather than only for a stubbed
    ``hash_file`` that might not raise what ``open()`` actually raises.
    """
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)
    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    os.chmod(network_path, 0o000)
    try:
        assert os.path.isfile(network_path), 'the guard this test is about must still see the file'
        status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path, method='MSC')
    finally:
        os.chmod(network_path, 0o644)
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert any('could not be read' in warning.lower() for warning in warnings), warnings


# --- 18. Sidecar records a different sa_cache_contract_version ---------------------------------------

def test_sidecar_with_mismatched_sa_cache_contract_version_is_rejected(tmp_path):
    """Test that a sidecar recording a different ``sa_cache_contract_version`` is rejected."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    # A REALISTIC SA YAML, not a stub: ``validate_sa_cache`` no longer gates on TS signal at all
    # (that judgment moved to ``t3.pdep.selector.select_from_sa_dict`` -- see FIX1 in cache.py), so
    # a stub would work here too now; kept realistic to match the sidecar's real-world shape.
    _write_realistic_sa_yaml(sa_path)

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    metadata_path = sa_cache_metadata_path(sa_path=sa_path)
    metadata = read_yaml_file(path=metadata_path)
    metadata['sa_cache_contract_version'] = metadata['sa_cache_contract_version'] + 1
    save_yaml_file(path=metadata_path, content=metadata)

    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path, method='MSC')
    assert status == CACHE_STATUS_CACHED_REJECTED
    # Pin the REASON, not just the verdict: 'rejected' is reachable by nine other branches.
    assert any('cache contract version' in warning for warning in warnings), warnings


def test_sidecar_missing_sa_cache_contract_version_key_is_rejected(tmp_path):
    """Test that a sidecar present but entirely missing the ``sa_cache_contract_version`` key is
    rejected (fail-closed on an unversioned file), distinct from test_missing_sidecar_is_rejected
    above (which covers the sidecar FILE being absent, not just this one key within it)."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    metadata_path = sa_cache_metadata_path(sa_path=sa_path)
    metadata = read_yaml_file(path=metadata_path)
    del metadata['sa_cache_contract_version']
    save_yaml_file(path=metadata_path, content=metadata)

    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path, method='MSC')
    assert status == CACHE_STATUS_CACHED_REJECTED
    assert any('cache contract version' in warning for warning in warnings), warnings


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
    # Pin the REASON, not just the verdict: 'rejected' is reachable by seven other branches.
    assert any('No cached sensitivity output' in warning for warning in warnings), warnings


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


# --- 21b. Sidecar parses but is not a dict -----------------------------------------------------------

def test_sidecar_parsing_to_non_dict_is_rejected(tmp_path):
    """Test that a sidecar which is valid YAML but not a mapping (e.g. a list) is rejected.

    This is distinct from test_malformed_sidecar_is_rejected_without_raising above: here
    read_sa_yaml_file() succeeds (the YAML is syntactically valid), but what it returns is not a dict,
    so none of the metadata.get(...) lookups below it would mean anything.
    """
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write(sa_path, 'some: sa\n')
    metadata_path = sa_cache_metadata_path(sa_path=sa_path)
    _write(metadata_path, '- item1\n- item2\n')  # parses to a list, not a dict

    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path)
    assert status == CACHE_STATUS_CACHED_REJECTED
    # Pin the REASON, not just the verdict: 'rejected' is reachable by seven other branches.
    assert any('malformed' in warning.lower() for warning in warnings), warnings


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

def test_cache_with_only_structural_zero_ts_coefficients_is_valid(tmp_path):
    """Test that a cache whose TS rows are all below the absolute floor is still CACHED_VALID (FIX1).

    Whether the SA output is USEFUL for criterion (b) (does it carry live TS signal) is a different
    question from whether the CACHE is trustworthy (right hashes, right method, right contract
    version, parseable). This case used to conflate the two: validate_sa_cache() rejected a cache
    that was, by every cache-validity measure, perfectly good -- it simply held a below-floor TS
    reading. That floor judgment now lives in ``t3.pdep.selector.select_from_sa_dict`` (evaluated
    per-reaction-key, not by scanning the whole dict), which is where a real trial run
    (T3-pdep-qm-trial-001) showed the granularity actually matters: a genuine Arkane SA landed
    real, non-zero-but-below-floor TS coefficients that this rejection would have thrown away and
    silently regenerated forever, instead of letting the selector report a considered verdict."""
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
    assert status == CACHE_STATUS_CACHED_VALID, warnings


def test_cache_with_no_ts_coefficient_recorded_is_valid(tmp_path):
    """Test that a cache with no recorded max_abs_ts_coefficient at all is still CACHED_VALID (FIX1).

    ``max_abs_ts_coefficient`` is recorded in the sidecar as provenance only -- useful for a human
    to read -- not as a cache-validity gate; see the module docstring and test above for why this
    guarantee moved to ``t3.pdep.selector.select_from_sa_dict``."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    save_yaml_file(path=sa_path, content={'structures': dict(), 'reactant1 + reactant2 <=> product': {}})

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    status, warnings = validate_sa_cache(sa_path=sa_path, network_path=network_path)
    assert status == CACHE_STATUS_CACHED_VALID, warnings


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


# --- 23. t_grid_clamp provenance in the sidecar -------------------------------------------------

def test_write_sa_cache_metadata_persists_t_grid_clamp_when_supplied(tmp_path):
    """A clamp decision made at write time (see t3.utils.network_thermo.TGridClampRecord) must
    survive past the run: when a caller supplies ``t_grid_clamp``, write_sa_cache_metadata must
    persist it into the sidecar under the ``t_grid_clamp`` key so a saved PDepNetworkSelection can
    later be traced back to whether its SA evidence rested on a clamped T grid."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)
    clamp_record = {'clamped': True, 'requested_t_max': 3200.0, 'thermo_ceiling': 3000.0,
                    'written_t_max': 3000.0, 'tlist_dropped': False, 'tlist_original_highest': None,
                    'skipped_species': []}

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2',
                            method='MSC', t_grid_clamp=clamp_record)

    metadata = read_yaml_file(sa_cache_metadata_path(sa_path=sa_path))
    assert metadata['t_grid_clamp'] == clamp_record


def test_write_sa_cache_metadata_omits_t_grid_clamp_key_entirely_when_not_supplied(tmp_path):
    """When no ``t_grid_clamp`` is supplied (the default), the sidecar must NOT gain a
    ``t_grid_clamp`` key at all -- not a ``None``/null value. A missing key and an explicit null
    both currently read back as 'unknown' via read_t_grid_clamp_record, but the sidecar itself
    must reflect 'this caller never supplied one', matching how a pre-existing (old) sidecar
    written before this feature existed would look, so the two situations stay indistinguishable
    from a caller who never asked for the field to begin with."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')

    metadata = read_yaml_file(sa_cache_metadata_path(sa_path=sa_path))
    assert 't_grid_clamp' not in metadata


def test_read_t_grid_clamp_record_round_trips_a_written_record(tmp_path):
    """The read side of the sidecar provenance must recover exactly what the write side put
    there, for both the clamped and the explicit not-clamped cases -- otherwise a
    PDepNetworkSelection reconstructed from this record could misreport its own T-grid
    provenance."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)
    not_clamped_record = {'clamped': False, 'requested_t_max': 2500.0, 'thermo_ceiling': 3000.0,
                          'written_t_max': 2500.0, 'tlist_dropped': False,
                          'tlist_original_highest': None, 'skipped_species': []}

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2',
                            method='MSC', t_grid_clamp=not_clamped_record)

    assert read_t_grid_clamp_record(sa_path=sa_path) == not_clamped_record


def test_read_t_grid_clamp_record_returns_none_for_a_sidecar_missing_the_key(tmp_path):
    """An 'old' sidecar -- written before this feature existed, i.e. lacking the ``t_grid_clamp``
    key entirely -- must read back as ``None`` (unknown provenance), NEVER as a value that could
    be mistaken for 'no clamp happened'. Conflating 'the key is absent' with 'clamped=False' would
    silently claim a verified-unclamped T grid for data that was never checked at all."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')

    assert read_t_grid_clamp_record(sa_path=sa_path) is None


def test_read_t_grid_clamp_record_returns_none_when_no_sidecar_exists(tmp_path):
    """No sidecar at all (e.g. SA output produced entirely outside T3) must also read as
    ``None``, exactly like the missing-key case -- read_t_grid_clamp_record never raises for this,
    since it exists purely to disclose provenance, not to gate anything."""
    sa_path = str(tmp_path / 'sa_coefficients.yml')  # no sidecar ever written next to this
    assert read_t_grid_clamp_record(sa_path=sa_path) is None


def test_read_t_grid_clamp_record_returns_none_for_unparseable_sidecar(tmp_path):
    """A sidecar that exists but is not valid YAML must also collapse to ``None`` rather than
    raising -- read_t_grid_clamp_record is a best-effort disclosure helper, and a caller consulting
    provenance for a broken sidecar should see 'unknown', the same outcome as if the sidecar were
    simply absent, not an exception it would need to specially handle."""
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write(sa_path, 'irrelevant\n')
    metadata_path = sa_cache_metadata_path(sa_path=sa_path)
    _write(metadata_path, '{not: [valid, yaml')

    assert read_t_grid_clamp_record(sa_path=sa_path) is None


def test_read_t_grid_clamp_record_returns_none_when_metadata_is_not_a_dict(tmp_path):
    """A sidecar whose top-level YAML content parses but is not a mapping (e.g. a bare list or
    scalar -- corrupted or hand-edited) must also collapse to ``None``, not raise or attempt a
    ``.get()`` on a non-dict."""
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write(sa_path, 'irrelevant\n')
    metadata_path = sa_cache_metadata_path(sa_path=sa_path)
    save_yaml_file(path=metadata_path, content=['not', 'a', 'dict'])

    assert read_t_grid_clamp_record(sa_path=sa_path) is None


def test_read_t_grid_clamp_record_returns_none_when_t_grid_clamp_value_is_not_a_dict(tmp_path):
    """If the ``t_grid_clamp`` key is present but its value is not itself a dict (e.g. hand-edited
    to a string or a bare boolean), that is not usable provenance -- read_t_grid_clamp_record must
    treat it the same as absent (``None``) rather than returning a malformed value a caller would
    then try to index into."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)
    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    metadata_path = sa_cache_metadata_path(sa_path=sa_path)
    metadata = read_yaml_file(metadata_path)
    metadata['t_grid_clamp'] = 'not-a-dict'
    save_yaml_file(path=metadata_path, content=metadata)

    assert read_t_grid_clamp_record(sa_path=sa_path) is None


@pytest.mark.parametrize('not_provenance', [
    {'requested_t_max': 2500.0, 'written_t_max': 2000.0},  # never says whether a clamp happened
    {'clamped': 'yes'},                                    # a verdict nobody actually stated
    {'clamped': 1},                                        # ditto: a truthy int is not an explicit bool
    {'clamped': True, 'thermo_ceiling': '2000'},           # a temperature that compares against nothing
    {'clamped': True, 'skipped_species': 'CH4'},           # would read back as three species, C, H and 4
    {'clamped': True, 'skipped_species': [42]},
])
def test_read_t_grid_clamp_record_returns_none_when_the_dict_is_not_clamp_provenance(tmp_path,
                                                                                     not_provenance):
    """A sidecar whose ``t_grid_clamp`` IS a dict but is not a TGridClampRecord rendering is exactly
    as useless as one whose value is a bare string, and must collapse the same way: to ``None``,
    unknown provenance. This is the last hole in this function's contract -- until now any dict at
    all was handed straight back, and t3.pdep.api copies whatever it returns into up to four live
    PDepNetworkSelection records, which then persist it as their own provenance. A malformed dict
    reaching that far is not caught anywhere downstream either, because nothing reads a key back out
    of it; it would simply sit in a saved decision record, misinforming the next human to open it.

    Collapsing rather than raising is this function's own documented contract: it exists purely to
    disclose provenance, so unknown provenance must never, by itself, cause a refusal."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)
    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    metadata_path = sa_cache_metadata_path(sa_path=sa_path)
    metadata = read_yaml_file(metadata_path)
    metadata['t_grid_clamp'] = not_provenance
    save_yaml_file(path=metadata_path, content=metadata)

    assert read_t_grid_clamp_record(sa_path=sa_path) is None


def test_read_t_grid_clamp_record_still_discloses_a_sidecar_a_newer_t3_wrote(tmp_path):
    """The mirror of the test above, and the one that keeps the new check from becoming a version
    tripwire. Sidecars outlive the version that wrote them: if a future T3 adds an eighth field to
    TGridClampRecord, every older T3 reading that sidecar must still disclose the provenance it CAN
    read, not silently discard the whole record as malformed. Losing provenance is the failure this
    function was written to prevent, so the check must not become a new way to cause it."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)
    from_a_newer_t3 = {'clamped': True, 'requested_t_max': 2500.0, 'thermo_ceiling': 2000.0,
                       'written_t_max': 2000.0, 'tlist_dropped': False,
                       'tlist_original_highest': None, 'skipped_species': [],
                       'clamp_reason': 'a field a future version added'}

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2',
                            method='MSC', t_grid_clamp=from_a_newer_t3)

    assert read_t_grid_clamp_record(sa_path=sa_path) == from_a_newer_t3


@pytest.mark.parametrize('not_provenance', [
    {'requested_t_max': 2500.0},
    {'clamped': 'yes'},
    {'clamped': True, 'thermo_ceiling': '2000'},
    {'clamped': True, 'skipped_species': 'CH4'},
])
def test_write_sa_cache_metadata_refuses_provenance_its_own_reader_would_discard(tmp_path,
                                                                                 not_provenance):
    """The write side of the leniency, and the half that makes it safe. read_t_grid_clamp_record
    collapses anything it cannot trust to None, which is right for a sidecar that may be foreign or
    predate a field -- but aimed at T3's OWN writer it would mean provenance this process actually
    had disappearing between the write and the read, with nothing logged on either side. So the
    asymmetry is deliberate and enforced from both ends: a possibly-foreign file is read leniently, a
    file T3 writes is written strictly."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)

    with pytest.raises(ValueError, match='provenance'):
        write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2',
                                method='MSC', t_grid_clamp=not_provenance)


def test_write_sa_cache_metadata_leaves_no_sidecar_behind_when_it_refuses(tmp_path):
    """A refused write must not have replaced a sidecar that was already there. The refusal happens
    before anything is read or written, so the previously recorded provenance survives intact rather
    than being traded for a malformed one."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)
    good = {'clamped': False, 'requested_t_max': 2500.0, 'thermo_ceiling': 3000.0,
            'written_t_max': 2500.0, 'tlist_dropped': False, 'tlist_original_highest': None,
            'skipped_species': []}
    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2',
                            method='MSC', t_grid_clamp=good)

    with pytest.raises(ValueError, match='provenance'):
        write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2',
                                method='MSC', t_grid_clamp={'clamped': 'yes'})

    assert read_t_grid_clamp_record(sa_path=sa_path) == good

# --- read_arkane_log_rmg_py_commit (RMG-Py provenance) --------------------------------------------
#
# Arkane does NOT run in T3's process: ARC's run_arkane shells out to
# ``micromamba run -n rmg_env python -m arkane input.py``, in a different conda environment. So the
# RMG-Py that produced an SA cannot be identified by introspecting this interpreter -- the only
# witness is the log Arkane itself writes into the output directory T3 controls.

ARKANE_LOG_PATH = os.path.join(TEST_DATA_BASE_PATH, 'pdep_arkane_log', 'arkane.log')
REAL_RMG_PY_COMMIT = 'e720866ae94eca51652978c15a0fb33c6827be67'


def test_read_arkane_log_rmg_py_commit_parses_a_real_arkane_log():
    """Test that the real commit is recovered from a genuine Arkane log. The SHA is not on the
    label line but on the line AFTER it, and the line after THAT is a date -- so a parser that
    takes the label line, or greps a single line, or takes every following line, gets nothing or
    garbage."""
    assert read_arkane_log_rmg_py_commit(ARKANE_LOG_PATH) == REAL_RMG_PY_COMMIT


def test_read_arkane_log_rmg_py_commit_does_not_return_the_rmg_database_commit():
    """Test that the RMG-database stanza, which shares the 'current git HEAD for' prefix and
    appears a few lines later in the same file, is not mistaken for the RMG-Py commit."""
    result = read_arkane_log_rmg_py_commit(ARKANE_LOG_PATH)
    assert result != '608f412ed7c109ed155dd877e13cd8959324e424'


def test_read_arkane_log_rmg_py_commit_returns_none_when_the_log_is_missing(tmp_path):
    """Test that a missing log yields None rather than raising: an Arkane job that died before
    writing its log must not crash the T3 iteration that is merely recording provenance."""
    assert read_arkane_log_rmg_py_commit(str(tmp_path / 'nonexistent' / 'arkane.log')) is None


def test_read_arkane_log_rmg_py_commit_returns_none_when_the_value_is_not_a_sha(tmp_path):
    """Test that a label followed by something that is not a commit hash yields None. Recording an
    arbitrary log line as a commit would be worse than recording nothing, because it would look
    like real provenance to whoever audits the sidecar later."""
    log_path = str(tmp_path / 'arkane.log')
    _write(log_path, 'The current git HEAD for RMG-Py is:\n\tnot-a-commit-hash\n')
    assert read_arkane_log_rmg_py_commit(log_path) is None


def test_read_arkane_log_rmg_py_commit_returns_none_when_the_label_is_absent(tmp_path):
    """Test that a log with no RMG-Py git-HEAD stanza at all yields None."""
    log_path = str(tmp_path / 'arkane.log')
    _write(log_path, 'Arkane execution initiated at Sat Aug  1 17:28:25 2026\n\nSome other line\n')
    assert read_arkane_log_rmg_py_commit(log_path) is None


def test_write_sa_cache_metadata_records_an_explicit_rmg_py_commit_verbatim(tmp_path):
    """Test that a commit supplied by the caller is recorded as-is. The commit now arrives from
    Arkane's log rather than from a local checkout, so it must NOT be re-derived from a path."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2',
                            method='MSC', rmg_py_commit=REAL_RMG_PY_COMMIT)

    metadata = read_yaml_file(sa_cache_metadata_path(sa_path=sa_path))
    assert metadata['rmg_py_commit'] == REAL_RMG_PY_COMMIT


def test_write_sa_cache_metadata_still_falls_back_to_the_checkout_when_no_commit_is_given(tmp_path):
    """Test that the pre-existing rmg_py_path route still works when no explicit commit is passed,
    so a caller that has a real checkout on hand is not forced through the log."""
    network_path = str(tmp_path / 'network4_2.py')
    _write(network_path, 'reaction(label="reaction1")\n')
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    _write_realistic_sa_yaml(sa_path)

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2',
                            method='MSC', rmg_py_path=str(tmp_path / 'not-a-git-checkout'))

    metadata = read_yaml_file(sa_cache_metadata_path(sa_path=sa_path))
    assert metadata['rmg_py_commit'] is None
