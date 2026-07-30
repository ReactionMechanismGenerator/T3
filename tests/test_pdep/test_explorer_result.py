#!/usr/bin/env python3
# encoding: utf-8

"""
t3 test_pdep test_explorer_result module

Tests ``PDepExplorationResult``.
"""

import yaml

import pytest

from t3.pdep.explorer.result import (
    EXPLORATION_STATUS_FAILED,
    EXPLORATION_STATUS_SKIPPED,
    EXPLORATION_STATUS_SUCCEEDED,
    PDepExplorationResult,
)
from t3.pdep.parser import PDepArkaneReaction
from t3.pdep.selector import PDepNetworkSelection


def _make_reaction():
    return PDepArkaneReaction(
        reactants=('A',),
        products=('B', 'C'),
        kinetics_type='Chebyshev',
        kinetics_params={'coeffs': [[1.0, 2.0], [None, 4.0]], 'Tmin': (300, 'K')},
        numeric_values=(1.0, 2.0, None, 4.0, 300),
    )


def test_status_succeeded_is_accepted_with_a_network_path():
    """Test that status='succeeded' with a non-empty network_paths is accepted (the minimal
    coherent 'succeeded' shape)."""
    result = PDepExplorationResult(
        network_id='network4_2', status=EXPLORATION_STATUS_SUCCEEDED,
        network_paths=('pdep/final/network1_full.py',))
    assert result.status == EXPLORATION_STATUS_SUCCEEDED


def test_status_failed_is_accepted_with_a_reason():
    """Test that status='failed' with a non-empty reasons is accepted (the minimal coherent
    'failed' shape)."""
    result = PDepExplorationResult(
        network_id='network4_2', status=EXPLORATION_STATUS_FAILED, reasons=('ME solve diverged',))
    assert result.status == EXPLORATION_STATUS_FAILED


def test_status_skipped_is_accepted_with_a_reason_and_no_artifacts():
    """Test that status='skipped' with a non-empty reasons and no artifact fields is accepted
    (the minimal coherent 'skipped' shape)."""
    result = PDepExplorationResult(
        network_id='network4_2', status=EXPLORATION_STATUS_SKIPPED,
        reasons=('network did not qualify for QM refinement',))
    assert result.status == EXPLORATION_STATUS_SKIPPED


def test_unknown_status_is_refused():
    """Test that a status string outside the three EXPLORATION_STATUS_* constants is refused,
    fail-closed, rather than silently accepted as some fourth undocumented state."""
    with pytest.raises(ValueError, match='status'):
        PDepExplorationResult(network_id='network4_2', status='maybe')


def test_skipped_with_network_paths_is_refused():
    """Test that a 'skipped' result carrying network_paths is refused: nothing ran, so there can
    be no final network artifact.

    Regression this guards: collapsing 'skipped' into 'failed'-shaped bookkeeping could leave a
    stale network_paths value from a reused/copied result on a status that says nothing ran.
    """
    with pytest.raises(ValueError, match='skipped'):
        PDepExplorationResult(
            network_id='network4_2', status=EXPLORATION_STATUS_SKIPPED,
            reasons=('did not qualify',), network_paths=('pdep/final/network1_full.py',))


def test_skipped_with_output_paths_is_refused():
    """Test that a 'skipped' result carrying output_paths is refused: nothing ran, so there can
    be no raw Arkane output artifact."""
    with pytest.raises(ValueError, match='skipped'):
        PDepExplorationResult(
            network_id='network4_2', status=EXPLORATION_STATUS_SKIPPED,
            reasons=('did not qualify',), output_paths=('output.py',))


def test_skipped_with_k_tp_as_written_is_refused():
    """Test that a 'skipped' result carrying k_tp_as_written is refused: nothing ran, so there
    can be no parsed k(T,P) entries."""
    with pytest.raises(ValueError, match='skipped'):
        PDepExplorationResult(
            network_id='network4_2', status=EXPLORATION_STATUS_SKIPPED,
            reasons=('did not qualify',), k_tp_as_written=(_make_reaction(),))


def test_failed_with_no_reasons_is_refused():
    """Test that a 'failed' result with empty reasons is refused: a negative outcome with no
    stated reason is unusable to a caller deciding what to do next."""
    with pytest.raises(ValueError, match='reasons'):
        PDepExplorationResult(network_id='network4_2', status=EXPLORATION_STATUS_FAILED)


def test_skipped_with_no_reasons_is_refused():
    """Test that a 'skipped' result with empty reasons is refused, mirroring the 'failed' case:
    a skip with no stated reason is equally unusable to a caller."""
    with pytest.raises(ValueError, match='reasons'):
        PDepExplorationResult(network_id='network4_2', status=EXPLORATION_STATUS_SKIPPED)


def test_succeeded_with_no_network_paths_is_refused():
    """Test that a 'succeeded' result with empty network_paths is refused.

    Regression this guards: success with no artifact is the fail-open shape this branch has
    closed repeatedly elsewhere (e.g. ArkaneExplorerAdapter's own success gating) -- a caller
    reading status=='succeeded' must be able to trust that network_paths is non-empty without
    re-checking it itself.
    """
    with pytest.raises(ValueError, match='network_paths'):
        PDepExplorationResult(network_id='network4_2', status=EXPLORATION_STATUS_SUCCEEDED)


def test_as_dict_round_trips_through_yaml_safe_dump():
    """Test that a fully-populated PDepExplorationResult's as_dict() output is accepted by
    yaml.safe_dump, i.e. it contains only plain, safe-dumpable types -- proven rather than
    assumed, since a stray tuple or dataclass instance would make safe_dump raise."""
    selection = PDepNetworkSelection(network_id='network4_2', qualified=True)
    result = PDepExplorationResult(
        network_id='network4_2',
        status=EXPLORATION_STATUS_SUCCEEDED,
        reasons=('exploration ran',),
        network_paths=('pdep/final/network1_full.py',),
        output_paths=('output.py',),
        k_tp_as_written=(_make_reaction(),),
        manifest={'input_hash': 'abc123', 'tool_versions': {'arkane': '1.0'}, 'nested': [(1, 2)]},
        selection=selection,
    )
    dumped = yaml.safe_dump(result.as_dict())
    assert isinstance(dumped, str)
    assert 'network4_2' in dumped


def test_as_dict_contains_no_tuples_or_dataclass_instances():
    """Test that as_dict() output contains no tuples and no dataclass instances anywhere in the
    structure, walked recursively rather than spot-checked at the top level -- a nested tuple
    (e.g. inside manifest, or inside a k_tp_as_written entry's kinetics_params) would otherwise
    survive undetected by a shallow check.
    """
    selection = PDepNetworkSelection(network_id='network4_2', qualified=True)
    result = PDepExplorationResult(
        network_id='network4_2',
        status=EXPLORATION_STATUS_SUCCEEDED,
        reasons=('exploration ran',),
        network_paths=('pdep/final/network1_full.py',),
        output_paths=('output.py',),
        k_tp_as_written=(_make_reaction(),),
        manifest={'nested': [(1, 2), {'deep': (3, 4)}]},
        selection=selection,
    )
    as_dict = result.as_dict()

    def _walk(value):
        assert not isinstance(value, tuple), f'Found a tuple: {value!r}'
        assert not is_dataclass_instance(value), f'Found a dataclass instance: {value!r}'
        if isinstance(value, dict):
            for item in value.values():
                _walk(item)
        elif isinstance(value, list):
            for item in value:
                _walk(item)

    def is_dataclass_instance(value) -> bool:
        return hasattr(value, '__dataclass_fields__') and not isinstance(value, type)

    _walk(as_dict)


def test_as_dict_renders_none_k_tp_direction_leaf():
    """Test that a None leaf inside a k_tp_as_written entry's kinetics_params survives all the
    way through PDepExplorationResult.as_dict(), not just PDepArkaneReaction.as_dict() alone --
    guarding the composition of the two as_dict() calls, not just each in isolation."""
    result = PDepExplorationResult(
        network_id='network4_2',
        status=EXPLORATION_STATUS_SUCCEEDED,
        network_paths=('pdep/final/network1_full.py',),
        k_tp_as_written=(_make_reaction(),),
    )
    as_dict = result.as_dict()
    assert as_dict['k_tp_as_written'][0]['kinetics_params']['coeffs'][1][0] is None


def test_a_failed_result_may_not_claim_a_usable_network_or_rate():
    """
    status='failed' carrying network_paths or k(T,P) must be refused.

    Those fields assert a usable result exists, which is the claim the failure denies. Would catch a
    caller (or a future refactor of explore_pdep_network) that forwards the adapter's artifacts
    uniformly on both outcomes, producing a record that reads as a success to anyone checking the
    artifact fields rather than the status.
    """
    with pytest.raises(ValueError, match='assert a usable result exists'):
        PDepExplorationResult(network_id='n1', status=EXPLORATION_STATUS_FAILED,
                              reasons=('did not converge',), network_paths=('/x/net.py',))


def test_a_failed_result_may_still_carry_its_output_paths():
    """
    The refusal above must not over-refuse into discarding failure diagnostics.

    A failed run's own logs are what a human needs to diagnose it. Refusing them alongside
    network_paths would be the tidier-looking rule and the wrong one -- the same over-refusal shape
    (refusing something that prevents nothing) this package has now corrected three times.
    """
    result = PDepExplorationResult(network_id='n1', status=EXPLORATION_STATUS_FAILED,
                                   reasons=('did not converge',), output_paths=('/x/output.py',))
    assert result.output_paths == ('/x/output.py',)


def test_the_manifest_is_copied_so_a_reported_result_cannot_be_rewritten():
    """
    Mutating the dict passed in must not change the result's manifest.

    A frozen dataclass holding someone else's live dict is frozen in name only: whoever passed it in
    can keep editing the provenance of a run that has already been reported.
    """
    manifest = {'rmgpy_revision': 'abc123'}
    result = PDepExplorationResult(network_id='n1', status=EXPLORATION_STATUS_SUCCEEDED,
                                   network_paths=('/x/net.py',), manifest=manifest)
    manifest['rmgpy_revision'] = 'tampered'
    manifest['injected'] = True
    assert result.manifest == {'rmgpy_revision': 'abc123'}
