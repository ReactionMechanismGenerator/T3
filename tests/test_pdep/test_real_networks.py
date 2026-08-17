#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_real_networks module

Exercises the PDep qualification gate against two complete, unmodified real chains rather than
against data a test constructed for itself. The fixtures and the reason this specific pair was
chosen are documented in ``tests/data/pdep_real_networks/README.md``.

A synthetic fixture can only ever assert that the gate does what the fixture's author already
believed. These two do something a synthetic pair almost never does: they disagree. The network
with the STRONGER sensitivity response is the one that does NOT qualify, because its kinetics are
not uncertain. That is the whole criterion -- sensitivity times uncertainty -- and it is why the
assertions below are on the semantics of each decision (which transition states, how many, certain
or uncertain, at what thresholds, with which warnings) and not merely on ``qualified``.
"""

import collections
import os

import pytest

from arc.common import read_yaml_file

from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.api import select_pdep_network
from t3.pdep.cache import read_arkane_log_rmg_py_commit, sa_cache_metadata_path
from t3.pdep.selector import (CACHE_STATUS_CACHED_VALID,
                              EVALUATION_STATUS_EVALUATED,
                              STRUCTURES_KEY,
                              TS_ENTRY_PREFIX,
                              )
from t3.pdep.yaml_safe import read_sa_yaml_file

REAL_NETWORKS_PATH = os.path.join(TEST_DATA_BASE_PATH, 'pdep_real_networks')

# The RMG-Py commit both runs were produced by. It contains the two ILT sensitivity fixes without
# which every transition-state coefficient in these files would be meaningless noise -- so this is
# not decoration, it is the reason the numbers below are worth asserting at all.
FIXTURE_RMG_PY_COMMIT = 'e720866ae94eca51652978c15a0fb33c6827be67'

# The floor below which a coefficient is indistinguishable from solver noise, derived from the
# 8368.0 J/mol perturbation Arkane applied. Recomputed by the gate on every call; pinned here
# because the ratios asserted further down are only meaningful against a known floor.
EXPECTED_COEFFICIENT_FLOOR = 1.1950286806883366e-07

EXPECTED_THRESHOLDS = {'relative_threshold': 0.001,
                       'min_delta_ln_k': 0.001,
                       'perturbation': 8368.0,
                       'coefficient_floor': EXPECTED_COEFFICIENT_FLOOR,
                       }


def network_paths(network_id: str) -> tuple:
    """
    Locate a vendored real network chain.

    Args:
        network_id (str): The network file stem, e.g. ``'network21_1'``.

    Returns:
        tuple: (network path (str), sa path (str), arkane log path (str)).
    """
    base = os.path.join(REAL_NETWORKS_PATH, network_id)
    return (os.path.join(base, f'{network_id}.py'),
            os.path.join(base, 'sensitivity', 'sa_coefficients.yml'),
            os.path.join(base, 'arkane.log'),
            )


def selection_triples(selection) -> dict:
    """
    Summarize a decision as (transition state, path reaction, uncertainty verdict) -> count.

    Args:
        selection (PDepNetworkSelection): The decision to summarize.

    Returns:
        dict: A count per distinct (ts_label, path_reaction_str, uncertain) triple.
    """
    return dict(collections.Counter((entry.ts_label, entry.path_reaction_str, entry.uncertain)
                                    for entry in selection.selected_ts))


def count_ts_rows(network_id: str) -> tuple:
    """
    Count the transition-state rows the SA file offers, and how many fall below the floor.

    Reads the vendored ``sa_coefficients.yml`` directly rather than through the selector, so that
    "the floor excluded N rows" is measured against the raw input instead of being restated from
    the selector's own output.

    Args:
        network_id (str): The network file stem.

    Returns:
        tuple: (total TS rows (int), TS rows below the coefficient floor (int)).
    """
    _, sa_path, _ = network_paths(network_id)
    sa_dict = read_sa_yaml_file(path=sa_path)
    total, below_floor = 0, 0
    for reaction_key, conditions in sa_dict.items():
        if reaction_key == STRUCTURES_KEY or not isinstance(reaction_key, str):
            continue
        for rows in conditions.values():
            for entry, coefficient in rows.items():
                if isinstance(entry, str) and entry.startswith(TS_ENTRY_PREFIX):
                    total += 1
                    below_floor += int(abs(float(coefficient)) < EXPECTED_COEFFICIENT_FLOOR)
    return total, below_floor


def select(network_id: str, method: str):
    """
    Run the gate over a vendored real chain exactly as a caller would.

    Cache validation is left ON deliberately: the vendored sidecar is part of the fixture, so a
    ``cached_valid`` status is itself an assertion that the copied files are unmodified and that
    the cache contract still accepts real T3 output.

    Args:
        network_id (str): The network file stem.
        method (str): The master-equation method the SA was generated with.

    Returns:
        PDepNetworkSelection: The decision.
    """
    network_path, sa_path, _ = network_paths(network_id)
    return select_pdep_network(network=network_path, sa_path=sa_path, method=method)


@pytest.mark.parametrize('network_id, method, network_source_hash, reactions_examined', [
    ('network21_1', 'CSE', 'sha256:bcbca07ce19d6d37b890ac1c82c3bb211d2ab6f4fa7243354287c7d789f04ab6', 3),
    ('network799_1', 'MSC', 'sha256:72efb347a7fde63f1589b0ece9d5dc05e357e2c5ca057ffb40088920c87db998', 3),
])
def test_a_real_chain_is_evaluated_from_a_valid_cache_without_warnings(network_id,
                                                                      method,
                                                                      network_source_hash,
                                                                      reactions_examined):
    """Test that a real, unmodified chain is actually evaluated rather than fail-closed."""
    selection = select(network_id, method)
    # `qualified` carries no signal unless the decision was genuinely computed, so this is the
    # precondition for every other assertion in this module, not a separate nicety.
    assert selection.evaluation_status == EVALUATION_STATUS_EVALUATED
    assert selection.cache_status == CACHE_STATUS_CACHED_VALID
    assert selection.warnings == list()
    assert selection.network_id == network_id
    assert selection.network_source_hash == network_source_hash
    assert selection.method == method
    assert selection.thresholds == EXPECTED_THRESHOLDS
    assert selection.network_reactions_examined == reactions_examined


def test_a_sensitive_network_with_certain_kinetics_does_not_qualify():
    """Test the negative real chain: strong sensitivity, no uncertainty, therefore no QM."""
    selection = select('network21_1', 'CSE')
    assert selection.qualified is False

    # Refusing to qualify must NOT mean the gate found nothing -- it selected 29 transition states
    # across three path reactions and then declined anyway. Asserting the selection alongside the
    # boolean is what distinguishes this from an over-refusing gate that returns False by
    # collapsing `selected_ts` to empty.
    assert len(selection.selected_ts) == 29
    # Counted as (label, path reaction, verdict) triples rather than by label alone: a label count
    # would survive a join that attached the right transition state to the WRONG path reaction,
    # which is the failure that would silently corrupt the uncertainty verdict below.
    assert selection_triples(selection) == {('TS1', 'N2H4(357) <=> NH3NH(403)', False): 11,
                                            ('TS2', 'N2H4(357) <=> H2(16) + H2NN(S)(380)', False): 6,
                                            ('TS3', 'N2H4(357) <=> H2(16) + N2H2(362)', False): 12,
                                            }
    assert sorted(selection.direction_keys) == ['N2H4(357) <=> H2(16) + H2NN(S)(380)',
                                                'N2H4(357) <=> H2(16) + N2H2(362)',
                                                'N2H4(357) <=> NH3NH(403)',
                                                ]

    # The reason it does not qualify, stated as the data rather than as the verdict.
    assert selection.uncertain_path_reactions == list()
    assert selection.uncertain_ts_labels() == list()
    assert not any(ts.uncertain for ts in selection.selected_ts)


def test_an_uncertain_network_qualifies_with_every_selected_transition_state_uncertain():
    """Test the positive real chain: weaker sensitivity, wholly uncertain kinetics, so QM."""
    selection = select('network799_1', 'MSC')
    assert selection.qualified is True

    assert len(selection.selected_ts) == 20
    assert selection_triples(selection) == {('TS1', 'O-2(13598) + CO2(11) <=> O=C1OO1(21648)', True): 4,
                                            ('TS2', '[O]C([O])=O(8160) <=> O=C1OO1(21648)', True): 12,
                                            ('TS3', '[O]O[C]=O(8100) <=> O=C1OO1(21648)', True): 4,
                                            }
    assert sorted(selection.direction_keys) == ['O=C1OO1(21648) <=> O-2(13598) + CO2(11)',
                                                'O=C1OO1(21648) <=> [O]C([O])=O(8160)',
                                                'O=C1OO1(21648) <=> [O]O[C]=O(8100)',
                                                ]

    # Every single selected transition state is uncertain here, which is the exact mirror of
    # network21_1 and the reason the pair is worth keeping together.
    assert len(selection.uncertain_path_reactions) == 20
    assert selection.uncertain_ts_labels() == ['TS1', 'TS2', 'TS3']
    assert all(ts.uncertain for ts in selection.selected_ts)


def test_the_more_sensitive_network_is_the_one_that_does_not_qualify():
    """Test that qualification is sensitivity AND uncertainty, not sensitivity alone."""
    unqualified = select('network21_1', 'CSE')
    qualified = select('network799_1', 'MSC')

    strongest_unqualified = max(abs(ts.coefficient) for ts in unqualified.selected_ts)
    strongest_qualified = max(abs(ts.coefficient) for ts in qualified.selected_ts)
    assert strongest_unqualified == pytest.approx(0.0005007416716395361, rel=1e-9)
    assert strongest_qualified == pytest.approx(9.515582575266414e-05, rel=1e-9)

    # The load-bearing comparison: the network that does NOT qualify responds five times more
    # strongly than the one that does. Any gate that ranked on sensitivity alone would invert
    # these two verdicts, and no assertion on a single network can catch that.
    assert strongest_unqualified > 5 * strongest_qualified
    assert unqualified.qualified is False
    assert qualified.qualified is True

    # Both are far above the noise floor, so the verdicts turn on uncertainty rather than on one
    # of them being marginal.
    assert strongest_unqualified / EXPECTED_COEFFICIENT_FLOOR == pytest.approx(4190.2, rel=1e-4)
    assert strongest_qualified / EXPECTED_COEFFICIENT_FLOOR == pytest.approx(796.3, rel=1e-4)


@pytest.mark.parametrize('network_id, method, total_rows, below_floor', [
    ('network21_1', 'CSE', 36, 7),
    ('network799_1', 'MSC', 36, 16),
])
def test_the_coefficient_floor_excludes_real_rows_rather_than_admitting_everything(network_id,
                                                                                  method,
                                                                                  total_rows,
                                                                                  below_floor):
    """Test that the floor actually removes transition-state rows a real SA file offered."""
    total, measured_below_floor = count_ts_rows(network_id)
    assert (total, measured_below_floor) == (total_rows, below_floor)

    # The floor accounts for EVERY exclusion here: what the selector kept is exactly what the raw
    # SA file offered minus what falls below the floor. Asserting the arithmetic rather than just
    # "some rows were dropped" is what makes this a claim about the floor specifically -- a gate
    # that discarded rows for some other reason would break this equality even if the count of
    # survivors happened to match.
    selection = select(network_id, method)
    assert len(selection.selected_ts) == total - below_floor


def test_a_real_selection_straddles_the_floor_instead_of_clearing_it_comfortably():
    """Test that the surviving rows span the floor closely, not by orders of magnitude."""
    selection = select('network21_1', 'CSE')
    magnitudes = sorted(abs(entry.coefficient) for entry in selection.selected_ts)

    # The weakest selected transition state is only ~2.2x the floor while the strongest is 4190x.
    # The selection is therefore neither "everything survived" nor "only the obvious ones did",
    # which is precisely the regime a hand-built fixture tends to miss.
    assert magnitudes[0] == pytest.approx(2.5993580396977083e-07, rel=1e-9)
    assert magnitudes[0] / EXPECTED_COEFFICIENT_FLOOR == pytest.approx(2.175, rel=1e-3)
    assert magnitudes[-1] / magnitudes[0] > 1000


@pytest.mark.parametrize('network_id, requested_t_max, thermo_ceiling', [
    ('network21_1', 2500.0, 3000.0),
    ('network799_1', 3200.0, 5000.0),
])
def test_a_real_chain_records_an_unclamped_temperature_grid(network_id, requested_t_max, thermo_ceiling):
    """Test that the T-grid provenance of a real run is carried onto the decision."""
    method = 'CSE' if network_id == 'network21_1' else 'MSC'
    selection = select(network_id, method)
    # `None` here would mean unknown provenance, which is a different (and weaker) claim than
    # "explicitly not clamped" -- so the dict itself is asserted, not just its `clamped` flag.
    assert selection.t_grid_clamp == {'clamped': False,
                                      'requested_t_max': requested_t_max,
                                      'skipped_species': list(),
                                      'thermo_ceiling': thermo_ceiling,
                                      'tlist_dropped': False,
                                      'tlist_original_highest': None,
                                      'written_t_max': requested_t_max,
                                      }


@pytest.mark.parametrize('network_id', ['network21_1', 'network799_1'])
def test_the_vendored_arkane_log_names_the_rmg_py_that_carries_the_ilt_fixes(network_id):
    """Test that the log shipped beside each chain names the expected RMG-Py commit."""
    _, _, arkane_log_path = network_paths(network_id)
    # Scope, stated honestly: this asserts what the vendored LOG says, and nothing links that log
    # to the SA file beside it -- no hash, no shared identifier, and `validate_sa_cache()` never
    # looks at a log. So this catches a fixture swapped for one from a different run only insofar
    # as the log travelled with it. It is evidence about the fixture's origin, not a binding proof
    # that this log produced this sa_coefficients.yml.
    assert read_arkane_log_rmg_py_commit(arkane_log_path) == FIXTURE_RMG_PY_COMMIT


@pytest.mark.parametrize('network_id', ['network21_1', 'network799_1'])
def test_the_vendored_sidecars_still_record_an_unknown_rmg_py_commit(network_id):
    """Test that the sidecars' null provenance is the fixture's state, not an accident."""
    _, sa_path, _ = network_paths(network_id)
    metadata = read_yaml_file(sa_cache_metadata_path(sa_path))
    # These sidecars were written before T3 learned to read the commit out of Arkane's log, so the
    # field is null. That is asserted rather than merely explained in the README, so that anyone
    # "helpfully" filling it in has to confront this test and the reason -- the sidecar hashes are
    # what make the cache validate, and a hand-written commit here would be invented provenance.
    assert 'rmg_py_commit' in metadata
    assert metadata['rmg_py_commit'] is None
    assert metadata['rmg_py_path'] is None
