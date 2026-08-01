#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_api
"""

import datetime
import json
import math
import os
from pathlib import Path

import pytest
import yaml

from arc.common import read_yaml_file, save_yaml_file

import t3.main as t3_main
from t3.common import TEST_DATA_BASE_PATH
from t3.logger import Logger
import t3.pdep.api as t3_pdep_api
from t3.pdep.api import (explore_pdep_network,
                         load_pdep_budget_record,
                         load_pdep_exploration_results,
                         load_pdep_network_selections,
                         rank_pdep_networks,
                         save_pdep_budget_record,
                         save_pdep_exploration_results,
                         save_pdep_network_selections,
                         select_pdep_network,
                         )
from t3.pdep.budget import (BUDGET_ALGORITHM_VERSION,
                            BUDGET_OUTCOME_ADMITTED,
                            BUDGET_OUTCOME_REFUSED,
                            BUDGET_RECORD_SCHEMA_VERSION,
                            BUDGET_SKIP_EXCEEDS_BUDGET,
                            PDepBudgetNetworkOutcome,
                            PDepBudgetRecord,
                            )
from t3.pdep.cache import hash_file, sa_cache_metadata_path, write_sa_cache_metadata
from t3.pdep.explorer.config import PDepExplorerConfig
from t3.pdep.explorer.result import (ADMISSION_POLICY_CALLER_ADMITTED,
                                     ADMISSION_POLICY_QUALIFIED_SELECTION,
                                     ADMISSION_POLICY_UNGATED,
                                     EXPLORATION_RESULT_SCHEMA_VERSION,
                                     EXPLORATION_STATUS_FAILED,
                                     EXPLORATION_STATUS_SKIPPED,
                                     EXPLORATION_STATUS_SUCCEEDED,
                                     PDepExplorationResult,
                                     )
from t3.pdep.parser import PDepArkaneReaction, parse_pdep_network_file
from t3.pdep.selector import (CACHE_STATUS_CACHED_REJECTED,
                              CACHE_STATUS_CACHED_VALID,
                              CACHE_STATUS_GENERATED,
                              CACHE_STATUS_UNVALIDATED,
                              EVALUATION_STATUS_EVALUATED,
                              EVALUATION_STATUS_NOT_EVALUATED,
                              SELECTION_ALGORITHM_VERSION,
                              SELECTION_SCHEMA_VERSION,
                              PDepNetworkSelection,
                              SensitiveTransitionState,
                              select_from_sa_dict,
                              )
from t3.schema import T3Sensitivity

from tests.test_pdep._wiring_helpers import (NETWORK_NAME,
                                             NETWORK_REACTION_STR,
                                             build_pdep_rxns_to_explore as _build_pdep_rxns_to_explore,
                                             build_sa_dict as _build_sa_dict,
                                             build_t3 as _build_t3,
                                             candidate_sa_path as _candidate_sa_path,
                                             network_path as _network_path,
                                             )

PDEP_NETWORK_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'iteration_1', 'RMG', 'pdep')
NETWORK_PATH = os.path.join(PDEP_NETWORK_DIR, 'network4_2.py')
# Computed with the same primitive production code uses, rather than hardcoded: a literal here
# would testify that a selection matches this file while actually only matching a constant.
NETWORK_SOURCE_HASH = hash_file(path=NETWORK_PATH)

# Deliberately kept OUTSIDE pdep_network/: test_main.py::test_determine_species_from_pdep_network
# sets t3.paths['PDep SA'] to tests/data/pdep_network/iteration_1/PDep_SA and its teardown runs
# shutil.rmtree() on that path, and this SA fixture must not live there.
SA_PATH = os.path.join(TEST_DATA_BASE_PATH, 'pdep_sa', 'network4_2_MSC', 'sa_coefficients.yml')

TARGET_REACTION = 'C1rad(2) + C3ene(27) <=> C2ene(29) + C2rad(3)'

# A minimal synthetic network with one uncertain (estimated) path reaction, reused across the
# rank_pdep_networks tests below with different reactant/product labels and TS coefficients.
SYNTHETIC_NETWORK_TEXT_A = '''
reaction(
    label = 'reactionA',
    reactants = ['S1', 'S2'],
    products = ['S3'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated template [x]"""),
)
'''

SYNTHETIC_NETWORK_TEXT_B = '''
reaction(
    label = 'reactionB',
    reactants = ['S4', 'S5'],
    products = ['S6'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated template [x]"""),
)
'''


@pytest.fixture(scope='module')
def sa_dict():
    """The real Arkane sensitivity dictionary for network4_2."""
    return read_yaml_file(path=SA_PATH)


def _write(path, content):
    """Write a small text file at ``path``."""
    with open(path, 'w') as f:
        f.write(content)


# --- 1. FIX A: parity with a real T3.determine_species_from_pdep_network() run -----------------

def test_select_pdep_network_matches_t3_in_run_decision(tmp_path, monkeypatch):
    """Test that select_pdep_network() agrees with a real T3.determine_species_from_pdep_network()
    run on the same inputs -- not merely with select_from_sa_dict() called directly (that would
    only prove api.py calls the core function, not that it produces the same decision T3's actual
    in-run path produces for a live run, which is the guarantee this module exists to provide)."""
    t3 = _build_t3(tmp_path)
    pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
    sa_dict_value = _build_sa_dict(t3)
    sa_path = _candidate_sa_path(t3, method='CSE')
    os.makedirs(os.path.dirname(sa_path), exist_ok=True)
    save_yaml_file(sa_path, sa_dict_value)
    write_sa_cache_metadata(sa_path=sa_path, network_path=_network_path(t3), network_id=NETWORK_NAME, method='CSE')

    def _fail_if_called(*args, **kwargs):
        pytest.fail('run_arkane_job should not be invoked when a valid cache is present.')

    monkeypatch.setattr(t3_main, 'run_arkane_job', _fail_if_called)

    t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)
    t3_selection = t3.pdep_network_selections[0]

    api_selection = select_pdep_network(
        network=_network_path(t3),
        sa_path=sa_path,
        network_reaction=NETWORK_REACTION_STR,
        relative_threshold=t3.t3['sensitivity']['pdep_SA_threshold'],
        min_delta_ln_k=t3.t3['sensitivity']['pdep_min_delta_ln_k'],
        method='CSE',
    )

    assert api_selection.qualified == t3_selection.qualified
    assert api_selection.network_id == t3_selection.network_id
    assert api_selection.direction_key == t3_selection.direction_key
    assert api_selection.method == t3_selection.method
    assert {(entry.ts_label, entry.condition) for entry in api_selection.selected_ts} == \
           {(entry.ts_label, entry.condition) for entry in t3_selection.selected_ts}


# --- 1b. FIX B: network accepts a str path or an already-parsed PDepNetwork identically ---------

def test_select_pdep_network_accepts_str_path_or_parsed_network_identically(sa_dict):
    """Test that passing an already-parsed PDepNetwork produces the same decision as its path."""
    parsed = parse_pdep_network_file(path=NETWORK_PATH)
    via_path = select_pdep_network(network=NETWORK_PATH, sa_dict=sa_dict, network_reaction=TARGET_REACTION,
                                   relative_threshold=0.001)
    via_parsed = select_pdep_network(network=parsed, sa_dict=sa_dict, network_reaction=TARGET_REACTION,
                                     relative_threshold=0.001)
    assert via_path.as_dict() == via_parsed.as_dict()


# --- 1c. FIX E: module default constants are derived from the schema, not copied literals -------

def test_default_constants_match_schema_defaults():
    """Test that DEFAULT_RELATIVE_THRESHOLD/DEFAULT_MIN_DELTA_LN_K track the schema's own defaults."""
    from t3.pdep.api import DEFAULT_MIN_DELTA_LN_K, DEFAULT_RELATIVE_THRESHOLD
    assert DEFAULT_RELATIVE_THRESHOLD == T3Sensitivity.model_fields['pdep_SA_threshold'].default
    assert DEFAULT_MIN_DELTA_LN_K == T3Sensitivity.model_fields['pdep_min_delta_ln_k'].default


# --- 2. Aggregation over every reaction key when network_reaction is None ----------------------

def test_select_pdep_network_aggregates_all_reaction_keys(sa_dict):
    """Test that network_reaction=None combines a decision per reaction key (45 in this fixture)."""
    combined = select_pdep_network(network=NETWORK_PATH, sa_dict=sa_dict, network_reaction=None,
                                   relative_threshold=0.001)
    assert combined.network_reactions_examined == 45
    # 20 of 45 reaction keys qualify in this fixture (see test_selector.py's sweep test), so the
    # combined (any-of) decision must also qualify.
    assert combined.qualified is True
    assert combined.network_reaction is None


# --- 3. ValueError for neither or both of sa_path/sa_dict ---------------------------------------

def test_select_pdep_network_raises_value_error_when_neither_sa_source_given():
    """Test that omitting both sa_path and sa_dict raises ValueError."""
    with pytest.raises(ValueError):
        select_pdep_network(network=NETWORK_PATH)


def test_select_pdep_network_raises_value_error_when_both_sa_sources_given(sa_dict):
    """Test that giving both sa_path and sa_dict raises ValueError."""
    with pytest.raises(ValueError):
        select_pdep_network(network=NETWORK_PATH, sa_path=SA_PATH, sa_dict=sa_dict)


# --- 3b. FIX K: 'how' without explore=True raises ValueError ------------------------------------

def test_select_pdep_network_raises_value_error_for_how_without_explore(sa_dict):
    """Test that passing 'how' without explore=True raises ValueError (it has no meaning there)."""
    with pytest.raises(ValueError):
        select_pdep_network(network=NETWORK_PATH, sa_dict=sa_dict, how='some_value')


# --- 4. explore=True redirects to explore_pdep_network() ----------------------------------------

def test_select_pdep_network_raises_value_error_for_explore_naming_explore_pdep_network(sa_dict):
    """Test that explore=True raises ValueError naming explore_pdep_network() as the real entry
    point (select_pdep_network() never explores; explore/how are retired placeholders)."""
    with pytest.raises(ValueError, match='explore_pdep_network'):
        select_pdep_network(network=NETWORK_PATH, sa_dict=sa_dict, explore=True)


# --- 5. Cache validate/reject, then valid after writing the sidecar -----------------------------

def test_select_pdep_network_rejects_missing_cache_then_accepts_after_sidecar_written(tmp_path):
    """Test the cache-rejection safety behavior: no sidecar -> rejected & unqualified; after writing
    the sidecar, the same sa_path is trusted and produces a real decision."""
    network_path = str(tmp_path / 'network4_2.py')
    with open(NETWORK_PATH, 'r') as f:
        _write(network_path, f.read())
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    save_yaml_file(path=sa_path, content=read_yaml_file(path=SA_PATH))

    rejected = select_pdep_network(network=network_path, sa_path=sa_path, network_reaction=TARGET_REACTION,
                                   relative_threshold=0.001, method='MSC')
    assert rejected.cache_status == CACHE_STATUS_CACHED_REJECTED
    assert rejected.qualified is False
    assert rejected.warnings

    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')
    accepted = select_pdep_network(network=network_path, sa_path=sa_path, network_reaction=TARGET_REACTION,
                                   relative_threshold=0.001, method='MSC')
    assert accepted.cache_status == CACHE_STATUS_CACHED_VALID
    assert accepted.selected_ts  # a real decision was computed, not a cache-rejection placeholder


# --- 5b. FIX G: validate_cache=False records 'unvalidated', not a fabricated real cache status ---

def test_select_pdep_network_validate_cache_false_reports_unvalidated(tmp_path):
    """Test that validate_cache=False trusts sa_path without checking provenance, and records that
    honestly as CACHE_STATUS_UNVALIDATED rather than as though a real check had taken place."""
    network_path = str(tmp_path / 'network4_2.py')
    with open(NETWORK_PATH, 'r') as f:
        _write(network_path, f.read())
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    save_yaml_file(path=sa_path, content=read_yaml_file(path=SA_PATH))
    # Deliberately no sidecar written: with validate_cache=True this would be CACHE_STATUS_CACHED_REJECTED.

    selection = select_pdep_network(network=network_path, sa_path=sa_path, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001, method='MSC', validate_cache=False)
    assert selection.cache_status == CACHE_STATUS_UNVALIDATED
    assert selection.selected_ts  # a real decision was computed, not a cache-rejection placeholder


# --- 5d. t_grid_clamp provenance threading through select_pdep_network --------------------------

def test_select_pdep_network_records_t_grid_clamp_from_the_sa_sidecar(tmp_path):
    """Test that select_pdep_network(sa_path=...) reads t_grid_clamp out of the SA cache sidecar
    (via read_t_grid_clamp_record) and carries it on the returned decision -- the whole point of
    persisting the clamp record is for a caller acting on a saved selection to be able to recover
    it without re-reading the sidecar itself."""
    network_path = str(tmp_path / 'network4_2.py')
    with open(NETWORK_PATH, 'r') as f:
        _write(network_path, f.read())
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    save_yaml_file(path=sa_path, content=read_yaml_file(path=SA_PATH))
    clamp_record = {'clamped': True, 'requested_t_max': 3200.0, 'thermo_ceiling': 3000.0}
    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2',
                            method='MSC', t_grid_clamp=clamp_record)

    selection = select_pdep_network(network=network_path, sa_path=sa_path, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001, method='MSC')
    assert selection.t_grid_clamp == clamp_record


def test_select_pdep_network_records_no_t_grid_clamp_when_the_sidecar_has_none(tmp_path):
    """Test that a valid sidecar that simply predates t_grid_clamp (the common case for every
    sidecar written before this feature existed) yields t_grid_clamp=None on the decision, not a
    refusal and not a fabricated 'not clamped' -- unknown provenance must stay unknown."""
    network_path = str(tmp_path / 'network4_2.py')
    with open(NETWORK_PATH, 'r') as f:
        _write(network_path, f.read())
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    save_yaml_file(path=sa_path, content=read_yaml_file(path=SA_PATH))
    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2', method='MSC')

    selection = select_pdep_network(network=network_path, sa_path=sa_path, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001, method='MSC')
    assert selection.t_grid_clamp is None


def test_select_pdep_network_records_no_t_grid_clamp_for_a_direct_sa_dict_call(sa_dict):
    """Test that select_pdep_network(sa_dict=...) -- the no-sa_path, no-sidecar-at-all path used
    e.g. when a caller already has SA data in memory -- records t_grid_clamp=None, since there is
    no sidecar to read the provenance from in the first place."""
    selection = select_pdep_network(network=NETWORK_PATH, sa_dict=sa_dict, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001)
    assert selection.t_grid_clamp is None


def test_select_pdep_network_drops_t_grid_clamp_on_the_cache_rejected_placeholder(tmp_path):
    """Test that the CACHE_STATUS_CACHED_REJECTED early-return placeholder carries NO t_grid_clamp.

    This test previously pinned the opposite, on the reasoning that surfacing what the sidecar
    claimed beats silently dropping it. That was reversed deliberately: the field means "the T grid
    this decision rests on", and a ``not_evaluated`` decision rests on nothing -- the SA was never
    read. Worse, the value would come from the very sidecar ``validate_sa_cache`` just refused to
    trust, so it is provenance borrowed from a source declared untrustworthy one line earlier.
    ``None`` (unknown provenance) is already the honest, well-tested answer on the two neighbouring
    paths above, and this path is no better informed than they are.
    """
    network_path = str(tmp_path / 'network4_2.py')
    with open(NETWORK_PATH, 'r') as f:
        _write(network_path, f.read())
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    save_yaml_file(path=sa_path, content=read_yaml_file(path=SA_PATH))
    clamp_record = {'clamped': False, 'requested_t_max': 2500.0}
    write_sa_cache_metadata(sa_path=sa_path, network_path=network_path, network_id='network4_2',
                            method='MSC', t_grid_clamp=clamp_record)
    # Deliberately rewrite the sidecar with a mismatching hash so the cache is rejected: force the
    # CACHE_STATUS_CACHED_REJECTED early-return construction site rather than the normal path.
    metadata_path = sa_cache_metadata_path(sa_path)
    metadata = read_yaml_file(path=metadata_path)
    metadata['network_file_hash'] = 'sha256:' + '0' * 64
    save_yaml_file(path=metadata_path, content=metadata)

    rejected = select_pdep_network(network=network_path, sa_path=sa_path, network_reaction=TARGET_REACTION,
                                   relative_threshold=0.001, method='MSC')
    assert rejected.cache_status == CACHE_STATUS_CACHED_REJECTED
    assert rejected.t_grid_clamp is None
    # The sidecar really did hold a clamp record: this asserts the value was dropped on the way to
    # the decision, not that the fixture never wrote one.
    assert read_yaml_file(path=sa_cache_metadata_path(sa_path))['t_grid_clamp'] == clamp_record


# --- 5c. FIX 3: a bad threshold raises even on the cache-rejected early-return path --------------

def test_select_pdep_network_raises_for_bad_threshold_even_when_cache_would_be_rejected(tmp_path):
    """Test that a NaN relative_threshold raises ValueError from select_pdep_network(), even for a
    sa_path/network pair that would otherwise hit the CACHE_STATUS_CACHED_REJECTED early return (no
    sidecar written here) -- that path bypasses select_from_sa_dict()/coefficient_floor() entirely,
    so select_pdep_network() must validate the thresholds itself before reaching it."""
    network_path = str(tmp_path / 'network4_2.py')
    with open(NETWORK_PATH, 'r') as f:
        _write(network_path, f.read())
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    save_yaml_file(path=sa_path, content=read_yaml_file(path=SA_PATH))
    # Deliberately no sidecar written: this sa_path/network pair would be CACHE_STATUS_CACHED_REJECTED.

    with pytest.raises(ValueError, match='relative_threshold'):
        select_pdep_network(network=network_path, sa_path=sa_path, network_reaction=TARGET_REACTION,
                            relative_threshold=math.nan, method='MSC')


# --- 6. rank_pdep_networks ordering, and unparseable networks included as not_evaluated ----------

def test_rank_pdep_networks_orders_and_includes_unparseable_as_not_evaluated(tmp_path):
    """Test that qualifying networks are ranked by descending strongest evidence, and a network
    whose file cannot be found/parsed is included (not dropped) as a not_evaluated placeholder,
    rather than aborting or silently vanishing from the returned list (FIX C)."""
    network_a_path = str(tmp_path / 'networkA.py')
    _write(network_a_path, SYNTHETIC_NETWORK_TEXT_A)
    network_b_path = str(tmp_path / 'networkB.py')
    _write(network_b_path, SYNTHETIC_NETWORK_TEXT_B)
    missing_path = str(tmp_path / 'networkMissing.py')  # never written

    sa_dict_a = {'S1 + S2 <=> S3': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 0.5}}}
    sa_dict_b = {'S4 + S5 <=> S6': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 0.1}}}

    networks = [
        {'network_path': network_b_path, 'sa_dict': sa_dict_b},
        {'network_path': network_a_path, 'sa_dict': sa_dict_a},
        {'network_path': missing_path, 'sa_dict': {}},
    ]
    ranked = rank_pdep_networks(networks, relative_threshold=0.001)

    assert [selection.network_id for selection in ranked] == ['networkA', 'networkB', 'networkMissing']
    missing_selection = ranked[-1]
    assert missing_selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert missing_selection.qualified is False
    assert missing_selection.warnings


# --- 6b. FIX J: a malformed network entry does not abort rank_pdep_networks ----------------------

def test_rank_pdep_networks_survives_a_malformed_entry():
    """Test that a malformed entry (neither a mapping nor a (network_path, sa_path) tuple/list, and
    a mapping missing 'network_path') produces a not_evaluated placeholder instead of raising."""
    networks = [42, {'sa_path': SA_PATH}]
    ranked = rank_pdep_networks(networks, relative_threshold=0.001)

    assert len(ranked) == 2
    for selection in ranked:
        assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
        assert selection.qualified is False
        assert selection.warnings


# --- 6c. FIX 3: a bad threshold raises once, not as N per-network placeholders -------------------

def test_rank_pdep_networks_raises_for_bad_threshold_instead_of_swallowing_it():
    """Test that a NaN relative_threshold raises ValueError directly from rank_pdep_networks(),
    rather than being caught by the per-network ``except Exception`` and turned into misleading
    not_evaluated placeholders for every network in the batch -- a bad threshold is a single config
    error affecting the whole call, not a per-network failure."""
    networks = [{'network_path': NETWORK_PATH, 'sa_path': SA_PATH},
                {'network_path': NETWORK_PATH, 'sa_path': SA_PATH}]
    with pytest.raises(ValueError, match='relative_threshold'):
        rank_pdep_networks(networks, relative_threshold=math.nan)


# --- 6c. rank_pdep_networks orders qualified networks by descending delta_ln_k ------------------

def test_rank_pdep_networks_orders_qualified_by_delta_ln_k_descending(tmp_path):
    """Test that among qualified decisions, the one with the larger delta_ln_k (among its
    uncertain_path_reactions) sorts first, even when that contradicts alphabetical network_id
    order (Mutation C4: max_delta_ln_k -> 0.0 would make this fall through to the network_id
    tiebreak and report the alphabetical order instead)."""
    text_low = '''
reaction(
    label = 'reactionLow',
    reactants = ['S1', 'S2'],
    products = ['S3'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated template [x]"""),
)
'''
    text_high = '''
reaction(
    label = 'reactionHigh',
    reactants = ['S4', 'S5'],
    products = ['S6'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated template [x]"""),
)
'''
    # 'aLowDelta' sorts before 'zHighDelta' alphabetically, but carries the SMALLER delta_ln_k, so
    # the two orderings disagree -- this is what lets the test distinguish real ranking from the
    # network_id tiebreak alone.
    low_path = str(tmp_path / 'aLowDelta.py')
    _write(low_path, text_low)
    high_path = str(tmp_path / 'zHighDelta.py')
    _write(high_path, text_high)

    sa_dict_low = {'S1 + S2 <=> S3': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 0.1}}}
    sa_dict_high = {'S4 + S5 <=> S6': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS9': 0.9}}}

    networks = [
        {'network_path': low_path, 'sa_dict': sa_dict_low},
        {'network_path': high_path, 'sa_dict': sa_dict_high},
    ]
    ranked = rank_pdep_networks(networks, relative_threshold=0.001)

    assert all(selection.qualified for selection in ranked)
    assert [selection.network_id for selection in ranked] == ['zHighDelta', 'aLowDelta']


# --- 6d. rank_pdep_networks: not_evaluated ranks ahead of evaluated-but-unqualified --------------

def test_rank_pdep_networks_not_evaluated_ranks_before_unqualified(tmp_path):
    """Test that a not_evaluated decision (tier 1) sorts before an evaluated-but-unqualified
    decision (tier 2), even when that contradicts alphabetical network_id order (Mutation C3:
    forcing not_evaluated's tier to 2 would collapse the two tiers together and report the
    alphabetical order instead)."""
    text_unqualified = '''
reaction(
    label = 'reactionUnqualified',
    reactants = ['P', 'Q'],
    products = ['R'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated template [x]"""),
)
'''
    # 'aUnqualified' sorts before 'zMissing' alphabetically, but is evaluated-but-unqualified
    # (tier 2) while 'zMissing' is not_evaluated (tier 1) -- the two orderings disagree.
    unqualified_path = str(tmp_path / 'aUnqualified.py')
    _write(unqualified_path, text_unqualified)
    missing_path = str(tmp_path / 'zMissing.py')  # never written

    sa_dict_unqualified = {'P + Q <=> R': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS5': 0.05}}}
    networks = [
        {'network_path': unqualified_path, 'sa_dict': sa_dict_unqualified},
        {'network_path': missing_path, 'sa_dict': {}},
    ]
    # relative_threshold=2.0 guarantees the single TS row's coefficient never clears its own
    # (2x-scaled) cutoff, so the unqualified network is genuinely evaluated-but-unqualified.
    ranked = rank_pdep_networks(networks, relative_threshold=2.0)

    assert [selection.network_id for selection in ranked] == ['zMissing', 'aUnqualified']
    assert ranked[0].evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert ranked[1].evaluation_status == EVALUATION_STATUS_EVALUATED
    assert ranked[1].qualified is False


# --- 6e. rank_pdep_networks: tiebreak is network_id ascending, not input order -------------------

def test_rank_pdep_networks_ties_broken_by_network_id_ascending(tmp_path):
    """Test that two qualified decisions with equal delta_ln_k are ordered by network_id ascending,
    not by their order in the input list (Mutation C5: dropping network_id from the sort key would
    leave Python's stable sort to preserve input order instead, which here is the reverse)."""
    text = '''
reaction(
    label = 'reactionTie',
    reactants = ['{r1}', '{r2}'],
    products = ['{p1}'],
    transitionState = '{ts}',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated template [x]"""),
)
'''
    zzz_path = str(tmp_path / 'zzz.py')
    _write(zzz_path, text.format(r1='P1', r2='Q1', p1='R1', ts='TS5'))
    aaa_path = str(tmp_path / 'aaa.py')
    _write(aaa_path, text.format(r1='P2', r2='Q2', p1='R2', ts='TS6'))

    sa_dict_zzz = {'P1 + Q1 <=> R1': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS5': 0.5}}}
    sa_dict_aaa = {'P2 + Q2 <=> R2': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS6': 0.5}}}

    # 'zzz' appears FIRST in the input list, so if the tiebreak were dropped, stable sort would
    # report ['zzz', 'aaa'] instead of the network_id-ascending order asserted below.
    networks = [
        {'network_path': zzz_path, 'sa_dict': sa_dict_zzz},
        {'network_path': aaa_path, 'sa_dict': sa_dict_aaa},
    ]
    ranked = rank_pdep_networks(networks, relative_threshold=0.001)

    assert all(selection.qualified for selection in ranked)
    deltas = [selection.uncertain_path_reactions[0].delta_ln_k for selection in ranked]
    assert deltas[0] == pytest.approx(deltas[1])  # confirms this is genuinely a tie
    assert [selection.network_id for selection in ranked] == ['aaa', 'zzz']


# --- 6f. rank_pdep_networks: _unpack_network_entry's own diagnoses survive as warning text -------

def test_rank_pdep_networks_dict_entry_missing_network_path_names_the_key():
    """Test that a mapping-form entry missing 'network_path' is recorded as not_evaluated with a
    warning naming the missing key specifically (Mutation C6: disabling the ``network_path is
    None`` guard would let ``entry.get('network_path')`` (``None``) flow on into
    ``select_pdep_network``, which fails for a different reason and would no longer name
    'network_path' at all)."""
    networks = [{'sa_path': SA_PATH}]
    ranked = rank_pdep_networks(networks, relative_threshold=0.001)

    assert len(ranked) == 1
    selection = ranked[0]
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert len(selection.warnings) == 1
    assert "missing a 'network_path' key" in selection.warnings[0]


def test_rank_pdep_networks_short_tuple_entry_names_the_required_shape():
    """Test that a tuple/list-form entry shorter than (network_path, sa_path) is recorded as
    not_evaluated with a warning naming the required shape and the actual element count (Mutation
    C7: disabling the ``len(entry) < 2`` guard would let a 1-element entry unpack ``sa_path`` as
    ``None`` and fall through to a different, unrelated failure downstream)."""
    networks = [('only_network_path',)]
    ranked = rank_pdep_networks(networks, relative_threshold=0.001)

    assert len(ranked) == 1
    selection = ranked[0]
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert len(selection.warnings) == 1
    assert 'must have at least (network_path, sa_path)' in selection.warnings[0]
    assert 'got 1 element(s)' in selection.warnings[0]


# --- 7. save_pdep_network_selections round trip and JSON-serializability ------------------------

def test_save_pdep_network_selections_round_trips_and_is_json_serializable(tmp_path, sa_dict):
    """Test that saved selections round-trip through read_yaml_file and are JSON-serializable, and
    carry a selection_schema_version envelope marker (FIX M)."""
    network = parse_pdep_network_file(path=NETWORK_PATH)
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001)
    path = str(tmp_path / 'selections.yml')
    returned_path = save_pdep_network_selections(path=path, selections=[selection])

    assert returned_path == path
    loaded = read_yaml_file(path=path)
    assert loaded == {'selection_schema_version': SELECTION_SCHEMA_VERSION, 'selections': [selection.as_dict()]}
    serialized = json.dumps(loaded)
    assert isinstance(serialized, str)


# --- 7b. Every PDepNetworkSelection carries non-None schema/algorithm version markers ------------

def test_select_pdep_network_cache_rejected_carries_version_markers(tmp_path):
    """Test that a decision from the cache-rejected path (api.py's cache-rejection thresholds
    site) carries non-None selection_schema_version/selection_algorithm_version, even though that
    site never mentions either key explicitly -- they come from the dataclass field defaults, not
    from anything a caller must remember to pass."""
    network_path = str(tmp_path / 'network4_2.py')
    with open(NETWORK_PATH, 'r') as f:
        _write(network_path, f.read())
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    save_yaml_file(path=sa_path, content=read_yaml_file(path=SA_PATH))
    # Deliberately no sidecar written, so this takes the cache-rejected path.

    rejected = select_pdep_network(network=network_path, sa_path=sa_path, network_reaction=TARGET_REACTION,
                                   relative_threshold=0.001, method='MSC')
    assert rejected.cache_status == CACHE_STATUS_CACHED_REJECTED
    assert rejected.selection_schema_version == SELECTION_SCHEMA_VERSION
    assert rejected.selection_algorithm_version == SELECTION_ALGORITHM_VERSION
    # The version markers above come from field DEFAULTS, so they cannot catch a construction site
    # that forgets them. network_source_hash has no such default -- it defaults to None, and
    # explore_pdep_network refuses a None hash BEFORE it reaches the evaluation gate. So dropping it
    # here would turn a clean 'skipped' into a ValueError blaming the decision for the constructor's
    # omission. Mutation-checked: deleting the kwarg at api.py:165 left the whole suite green.
    assert rejected.network_source_hash == parse_pdep_network_file(path=network_path).source_hash


def test_select_pdep_network_no_reaction_keys_carries_version_markers():
    """Test that a decision from the no-SA-keys path (api.py's "no reaction keys found" thresholds
    site) carries non-None selection_schema_version/selection_algorithm_version, for the same
    reason as the cache-rejected case above."""
    selection = select_pdep_network(network=NETWORK_PATH, sa_dict={}, network_reaction=None,
                                    relative_threshold=0.001)
    assert selection.selection_schema_version == SELECTION_SCHEMA_VERSION
    assert selection.selection_algorithm_version == SELECTION_ALGORITHM_VERSION
    # See the sibling test above: this is the assertion the defaults cannot make.
    assert selection.network_source_hash == parse_pdep_network_file(path=NETWORK_PATH).source_hash


# --- 7c. selection_schema_version survives nesting inside PDepExplorationResult.as_dict() --------

def test_exploration_result_as_dict_nests_selection_schema_version(sa_dict):
    """Test that PDepExplorationResult.as_dict()['selection']['selection_schema_version'] is
    present: a version marker on the enclosing envelope alone cannot describe a nested record, so
    this must survive as a field on the nested selection itself."""
    network = parse_pdep_network_file(path=NETWORK_PATH)
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001)
    result = PDepExplorationResult(network_id=network.network_id, status=EXPLORATION_STATUS_SKIPPED,
                                   reasons=('not qualified',), selection=selection)
    rendered = result.as_dict()
    assert rendered['selection']['selection_schema_version'] == SELECTION_SCHEMA_VERSION
    assert rendered['selection']['selection_algorithm_version'] == SELECTION_ALGORITHM_VERSION


# --- 7d. save_pdep_exploration_results round trip, YAML-safety, and envelope shape ---------------

def _make_exploration_result(sa_dict, status=EXPLORATION_STATUS_SUCCEEDED, selection=None,
                             include_selection=True):
    """Build a realistic ``PDepExplorationResult`` for the save-round-trip tests below: a nested
    selection (real, from ``select_from_sa_dict``, unless the caller overrides it), a nested
    ``PDepArkaneReaction`` k(T,P) entry, and a free-form manifest containing tuples -- so saving it
    actually exercises the nested-selection, nested-reaction, and tuple-flattening paths at once."""
    network = parse_pdep_network_file(path=NETWORK_PATH)
    if include_selection and selection is None:
        selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=TARGET_REACTION,
                                        relative_threshold=0.001)
    reaction = PDepArkaneReaction(
        reactants=('A',),
        products=('B', 'C'),
        kinetics_type='Chebyshev',
        kinetics_params={'coeffs': [[1.0, 2.0], [None, 4.0]], 'Tmin': (300, 'K')},
        numeric_values=(1.0, 2.0, None, 4.0, 300),
    )
    if status == EXPLORATION_STATUS_SUCCEEDED:
        return PDepExplorationResult(network_id=network.network_id, status=status,
                                     network_paths=('pdep/final/network1_full.py',),
                                     k_tp_as_written=(reaction,),
                                     manifest={'nested': [(1, 2), {'deep': (3, 4)}]},
                                     selection=selection)
    return PDepExplorationResult(network_id=network.network_id, status=status,
                                 reasons=('did not qualify',), selection=selection)


def test_save_pdep_exploration_results_round_trips_and_is_json_serializable(tmp_path, sa_dict):
    """Test that saved exploration results round-trip through read_yaml_file, carry an
    exploration_result_schema_version envelope marker, and are JSON-serializable -- mirroring
    test_save_pdep_network_selections_round_trips_and_is_json_serializable above."""
    result = _make_exploration_result(sa_dict)
    path = str(tmp_path / 'exploration_results.yml')
    returned_path = save_pdep_exploration_results(path=path, results=[result])

    assert returned_path == path
    loaded = read_yaml_file(path=path)
    assert loaded == {'exploration_result_schema_version': EXPLORATION_RESULT_SCHEMA_VERSION,
                      'results': [result.as_dict()]}
    serialized = json.dumps(loaded)
    assert isinstance(serialized, str)


def test_save_pdep_exploration_results_produces_no_python_object_tags(tmp_path, sa_dict):
    """Test that the saved file is YAML-safe: since a result nests a selection, nested
    PDepArkaneReaction k(T,P) entries, and a free-form manifest, this is the assertion that
    actually protects the on-disk format from silently regressing to an unsafe dump."""
    result = _make_exploration_result(sa_dict)
    path = str(tmp_path / 'exploration_results.yml')
    save_pdep_exploration_results(path=path, results=[result])

    with open(path, 'r') as f:
        raw = f.read()
    assert '!!python/' not in raw


def test_save_pdep_exploration_results_handles_selection_none(tmp_path):
    """Test that a result with selection=None (the skipped-before-any-decision case) serializes
    without blowing up, and round-trips with a null selection."""
    result = PDepExplorationResult(network_id='network4_2', status=EXPLORATION_STATUS_SKIPPED,
                                   reasons=('no selection available',), selection=None)
    path = str(tmp_path / 'exploration_results.yml')
    save_pdep_exploration_results(path=path, results=[result])

    loaded = read_yaml_file(path=path)
    assert loaded['results'] == [result.as_dict()]
    assert loaded['results'][0]['selection'] is None


def test_saved_exploration_result_selection_carries_its_own_schema_version(tmp_path, sa_dict):
    """Test that the nested selection in the SAVED file carries its own selection_schema_version --
    the assertion that proves omitting a selection-version key from the envelope (a deliberate
    design choice; see save_pdep_exploration_results' docstring) is safe, because each nested
    selection record already self-describes."""
    result = _make_exploration_result(sa_dict)
    path = str(tmp_path / 'exploration_results.yml')
    save_pdep_exploration_results(path=path, results=[result])

    loaded = read_yaml_file(path=path)
    saved_selection = loaded['results'][0]['selection']
    assert saved_selection['selection_schema_version'] == SELECTION_SCHEMA_VERSION
    assert saved_selection['selection_algorithm_version'] == SELECTION_ALGORITHM_VERSION
    assert 'selection_schema_version' not in loaded  # not duplicated on the envelope


def test_save_pdep_exploration_results_writes_a_valid_envelope_for_an_empty_list(tmp_path):
    """Test that saving an empty list still writes a valid envelope, not an omitted/empty file."""
    path = str(tmp_path / 'exploration_results.yml')
    returned_path = save_pdep_exploration_results(path=path, results=[])

    assert returned_path == path
    loaded = read_yaml_file(path=path)
    assert loaded == {'exploration_result_schema_version': EXPLORATION_RESULT_SCHEMA_VERSION, 'results': []}


# --- 8. log_pdep_network_summary does not raise for populated or empty lists ---------------------

def test_log_pdep_network_summary_does_not_raise(tmp_path, sa_dict):
    """Test that Logger.log_pdep_network_summary handles both a populated and an empty list."""
    network = parse_pdep_network_file(path=NETWORK_PATH)
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001)
    qualifying_selection = select_pdep_network(network=NETWORK_PATH, sa_dict=sa_dict, network_reaction=None,
                                               relative_threshold=0.001)
    assert qualifying_selection.qualified  # sanity check: this is the populated-list case

    logger = Logger(project='test_pdep_api', project_directory=str(tmp_path), verbose=None,
                    t0=datetime.datetime.now())
    logger.log_pdep_network_summary(selections=[selection, qualifying_selection])
    logger.log_pdep_network_summary(selections=[])


# --- 9. explore_pdep_network() -------------------------------------------------------------------
#
# All of these use a fake explorer adapter/factory (monkeypatched onto t3.pdep.api.explorer_factory,
# the name actually looked up inside explore_pdep_network()) rather than running real Arkane: the
# point of these tests is to pin explore_pdep_network()'s OWN guards -- the budget gate, the
# filesystem containment re-checks, and the no-catch-all rule around adapter.explore() -- not to
# exercise a real PES exploration.

class _FakeAdapter:
    """
    A minimal stand-in for a PESExplorerAdapter, recording what it was constructed with.

    Deliberately models the ABC's CONTRACT (``get_networks()``, ``reasons``, ``output_paths``,
    ``manifest``) rather than ArkaneExplorerAdapter's internals. An earlier version of this double
    exposed ``final_network_paths``, which is Arkane's own attribute and is not part of the adapter
    interface at all -- so the suite's model of a "conforming adapter" was one no second explorer
    would have had to satisfy, and the API's coupling to a single concrete class was invisible here.
    A double that testifies to a contract the real interface does not define proves nothing.
    """

    def __init__(self, succeed=True, reasons=(), raise_error=None, on_construct=None,
                 network_paths=('/fake/final_network.py',), **kwargs):
        self._succeed = succeed
        self.reasons = reasons
        self._raise_error = raise_error
        self._network_paths = network_paths
        self.output_paths = ('/fake/output1.py',)
        self.manifest = {'fake': True}
        self.construct_kwargs = kwargs
        if on_construct is not None:
            on_construct(kwargs)

    def explore(self) -> bool:
        if self._raise_error is not None:
            raise self._raise_error
        return self._succeed

    def get_networks(self) -> tuple:
        if not self._succeed:
            raise RuntimeError('get_networks() was called after a failed explore().')
        return self._network_paths

    def get_k_tp(self):
        return tuple()


def _make_fake_factory(monkeypatch, **adapter_kwargs):
    """Monkeypatch t3.pdep.api.explorer_factory with a fake that records whether it was called and
    returns a _FakeAdapter built from adapter_kwargs. Returns the mutable call-record list."""
    calls = []

    def _fake_factory(**kwargs):
        calls.append(kwargs)
        return _FakeAdapter(**adapter_kwargs)

    monkeypatch.setattr(t3_pdep_api, 'explorer_factory', _fake_factory)
    return calls


def _make_config(trusted_output_root, output_directory, method='CSE'):
    """Build a valid PDepExplorerConfig pointed at the given (not-yet-existing) directories."""
    return PDepExplorerConfig(
        explorer='FakeExplorer',
        trusted_output_root=trusted_output_root,
        output_directory=output_directory,
        seed_species=('H',),
        method=method,
    )


def test_explore_pdep_network_refuses_parsed_network_object(tmp_path, monkeypatch):
    """Test that a parsed PDepNetwork object (rather than a path string) is refused, with the reason
    stated in the message."""
    _make_fake_factory(monkeypatch)
    parsed = parse_pdep_network_file(path=NETWORK_PATH)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'))
    with pytest.raises(ValueError, match='path string'):
        explore_pdep_network(network_path=parsed, config=config)


def test_explore_pdep_network_refuses_selection_method_mismatch(tmp_path, monkeypatch):
    """Test that a selection whose method disagrees with config.method is refused."""
    _make_fake_factory(monkeypatch)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=True, method='MSC')
    with pytest.raises(ValueError, match='method'):
        explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection)


def test_explore_pdep_network_refuses_selection_network_id_mismatch(tmp_path, monkeypatch):
    """Test that a selection whose network_id disagrees with the parsed network is refused."""
    _make_fake_factory(monkeypatch)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='some_other_network', qualified=True, method='CSE')
    with pytest.raises(ValueError, match='network_id'):
        explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection)


def test_explore_pdep_network_skips_non_qualifying_selection_without_constructing_adapter(tmp_path, monkeypatch):
    """Test that a non-qualifying selection produces status='skipped' carrying selection.reason(),
    and -- the whole point of the budget gate -- that the explorer factory is never called."""
    calls = _make_fake_factory(monkeypatch)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    output_directory = os.path.join(trusted_root, 'run1')
    config = _make_config(trusted_root, output_directory, method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=False, method='CSE',
                                     network_source_hash=NETWORK_SOURCE_HASH)

    result = explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection)

    assert result.status == EXPLORATION_STATUS_SKIPPED
    assert result.reasons == (selection.reason(),)
    # Identity is deliberately no longer guaranteed here: PDepExplorationResult.__post_init__ now
    # deep-copies ``selection`` (the same defense it already gave ``manifest``), so the result
    # snapshots the decision rather than aliasing the caller's mutable object. Equality still
    # tests what this assertion means: the result carries the decision it was given.
    assert result.selection == selection
    assert calls == []  # the factory (and therefore any adapter) was never invoked
    # A skipped exploration must not have touched the filesystem either.
    assert not os.path.exists(output_directory)


def test_explore_pdep_network_forwards_the_parsed_content_hash_to_the_factory(tmp_path, monkeypatch):
    """Test the FIRST link of the TOCTOU chain: this function's own hash check proves only what the
    bytes were when it parsed them, and the writer reopens the same pathname later. Forwarding
    parsed_network.source_hash is what binds the writer's read to the approved bytes -- and deleting
    that one keyword argument is invisible to every other test here, since they all assert on the
    result rather than on what the factory was handed (confirmed by mutation)."""
    calls = _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')

    explore_pdep_network(network_path=NETWORK_PATH, config=config)

    assert len(calls) == 1
    assert calls[0]['expected_source_hash'] == NETWORK_SOURCE_HASH


def test_explore_pdep_network_accepts_a_qualified_but_partially_evaluated_selection(tmp_path, monkeypatch):
    """Test that a selection which qualified but whose coverage was partial is ACCEPTED. combine()
    marks any partially covered aggregate 'not_evaluated' because that is the truth about coverage,
    but a positive verdict there is backed by whichever evaluated component qualified -- and
    combine() only counts evaluated components' votes. Refusing it would be an over-refusal, this
    branch's other failure mode. The asymmetry is deliberate: a partial 'no' is unsupported, a
    partial 'yes' is not."""
    calls = _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=True, method='CSE',
                                     network_source_hash=NETWORK_SOURCE_HASH)
    selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED

    result = explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection)

    assert len(calls) == 1
    assert result.status == EXPLORATION_STATUS_SUCCEEDED


def test_explore_pdep_network_refuses_a_selection_with_no_source_hash(tmp_path, monkeypatch):
    """Test that a selection carrying no content binding is refused rather than gated on its file
    stem. network_id matches every revision of network4_2.py, so accepting a hash-less selection
    would be the fail-open this check exists to close."""
    calls = _make_fake_factory(monkeypatch)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=True, method='CSE')
    # Matched on the wording specific to the no-hash branch, NOT merely on 'network_source_hash':
    # a hash-less selection also trips the mismatch check below it (None != any real hash), so a
    # loose match would pass with the no-hash branch deleted -- confirmed by mutation. The branch
    # earns its place by DIAGNOSIS: the mismatch message says the network file has changed since
    # the decision was made, which is a false statement about a decision that never recorded one.
    with pytest.raises(ValueError, match='carries no binding to the content'):
        explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection)
    assert calls == []


def test_explore_pdep_network_refuses_a_selection_made_against_different_content(tmp_path, monkeypatch):
    """Test that a QUALIFYING selection whose recorded hash does not match the file on disk is
    refused: the decision describes a revision of the network that is not the one about to be
    explored. This is the case network_id structurally cannot catch -- same stem, different bytes."""
    calls = _make_fake_factory(monkeypatch)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=True, method='CSE',
                                     network_source_hash='sha256:' + 'e' * 64)
    with pytest.raises(ValueError, match='has changed since'):
        explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection)
    assert calls == []


def test_explore_pdep_network_refuses_a_selection_after_the_network_file_is_edited(tmp_path, monkeypatch):
    """Test the end-to-end shape of the hazard: a selection made against a real network file, the
    file then edited on disk, and the exploration refused. Both the selection and the check compute
    their hashes from real files, so nothing here is asserted against a constant."""
    calls = _make_fake_factory(monkeypatch)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    network_copy = tmp_path / 'network4_2.py'
    network_copy.write_bytes(open(NETWORK_PATH, 'rb').read())
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(
        network_id='network4_2', qualified=True, method='CSE',
        network_source_hash=parse_pdep_network_file(path=str(network_copy)).source_hash)

    # The exact hazard: the network file is rewritten (as RMG does every iteration) after the
    # decision was made, but before it is acted on.
    network_copy.write_bytes(open(NETWORK_PATH, 'rb').read() + b'\n# RMG rewrote this network\n')

    with pytest.raises(ValueError, match='has changed since'):
        explore_pdep_network(network_path=str(network_copy), config=config, selection=selection)
    assert calls == []


def test_explore_pdep_network_explores_anyway_when_selection_is_none(tmp_path, monkeypatch):
    """Test that omitting selection (None) runs the exploration unconditionally."""
    calls = _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')

    result = explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=None)

    assert len(calls) == 1
    assert result.status == EXPLORATION_STATUS_SUCCEEDED


def test_explore_pdep_network_accepts_a_real_selection_end_to_end(tmp_path, monkeypatch):
    """Test that a selection produced by select_pdep_network() against the SAME file passes every gate
    and reaches the adapter. This is the over-refusal guard for the content-hash check: every other
    test of it asserts a refusal, so without this one the check could have made the gated path
    impossible to use and the suite would still be green."""
    calls = _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='MSC')

    # network_reaction=None on purpose: this is the combined path, so the test also covers combine()
    # carrying network_source_hash through to the aggregate decision the gate actually receives.
    selection = select_pdep_network(network=NETWORK_PATH, sa_path=SA_PATH, network_reaction=None,
                                    relative_threshold=0.001, method='MSC', validate_cache=False)
    assert selection.qualified, 'fixture no longer qualifies; this test can no longer exercise the gate'
    assert selection.network_source_hash == NETWORK_SOURCE_HASH

    result = explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection)

    assert len(calls) == 1
    assert result.status == EXPLORATION_STATUS_SUCCEEDED


def test_explore_pdep_network_refuses_symlinked_trusted_output_root(tmp_path, monkeypatch):
    """Test that a trusted_output_root that is itself a symlink is refused, even though it resolves
    to a real directory and containment checks would otherwise pass."""
    _make_fake_factory(monkeypatch)
    real_root = tmp_path / 'real_root'
    os.makedirs(real_root)
    linked_root = tmp_path / 'root_link'
    linked_root.symlink_to(real_root, target_is_directory=True)
    config = _make_config(str(linked_root), os.path.join(str(linked_root), 'run1'), method='CSE')
    with pytest.raises(ValueError, match='symlink'):
        explore_pdep_network(network_path=NETWORK_PATH, config=config)


def test_explore_pdep_network_refuses_nonexistent_trusted_output_root(tmp_path, monkeypatch):
    """Test that a trusted_output_root that does not exist on disk is refused (explore_pdep_network()
    never creates the root itself)."""
    _make_fake_factory(monkeypatch)
    trusted_root = str(tmp_path / 'does_not_exist')  # deliberately never created
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    with pytest.raises(ValueError, match='trusted_output_root'):
        explore_pdep_network(network_path=NETWORK_PATH, config=config)


def test_explore_pdep_network_refuses_output_directory_outside_root(tmp_path, monkeypatch):
    """Test that explore_pdep_network() independently re-verifies output_directory containment
    rather than only trusting PDepExplorerConfig's construction-time check. A well-behaved config
    cannot even be constructed with an out-of-root output_directory (PDepExplorerConfig.__post_init__
    already refuses that), so this simulates a config tampered with after construction -- e.g. by a
    stale reference or a bug elsewhere -- to prove explore_pdep_network()'s own re-check is real and
    not merely decorative."""
    _make_fake_factory(monkeypatch)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    other_root = str(tmp_path / 'other_root')
    os.makedirs(other_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    # PDepExplorerConfig is frozen; object.__setattr__ is the documented way to bypass that for a
    # deliberate tamper, mirroring PDepExplorerConfig.__post_init__'s own use of object.__setattr__
    # for coerced field values.
    object.__setattr__(config, 'output_directory', os.path.join(other_root, 'run1'))
    with pytest.raises(ValueError, match='trusted_output_root'):
        explore_pdep_network(network_path=NETWORK_PATH, config=config)


def test_explore_pdep_network_creates_intermediate_dirs_but_not_output_directory_itself(tmp_path, monkeypatch):
    """Test that explore_pdep_network() creates the intermediate directories between the trusted
    root and output_directory, but never pre-creates output_directory itself -- that leaf directory
    must stay the adapter's own rule-0 atomic os.mkdir claim (see ArkaneExplorerAdapter's
    _claim_run_directory docstring)."""
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    output_directory = os.path.join(trusted_root, 'nested', 'deeper', 'run1')
    observed = {}

    def _on_construct(kwargs):
        # Captured at the moment the fake adapter is constructed -- i.e. right where the real
        # adapter's rule-0 os.mkdir claim would happen next.
        observed['parent_exists'] = os.path.isdir(os.path.dirname(output_directory))
        observed['output_directory_exists'] = os.path.isdir(output_directory)

    _make_fake_factory(monkeypatch, succeed=True, on_construct=_on_construct)
    config = _make_config(trusted_root, output_directory, method='CSE')

    explore_pdep_network(network_path=NETWORK_PATH, config=config)

    assert observed['parent_exists'] is True
    assert observed['output_directory_exists'] is False


def test_explore_pdep_network_reports_failed_status_with_reasons_no_exception(tmp_path, monkeypatch):
    """Test that an ordinary Arkane run failure is a RECORDED result (status='failed' with reasons),
    not an exception."""
    _make_fake_factory(monkeypatch, succeed=False, reasons=('Arkane reported job failure.',))
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')

    result = explore_pdep_network(network_path=NETWORK_PATH, config=config)

    assert result.status == EXPLORATION_STATUS_FAILED
    assert result.reasons == ('Arkane reported job failure.',)


def test_explore_pdep_network_propagates_runtime_error_from_adapter(tmp_path, monkeypatch):
    """Test that a RuntimeError raised out of adapter.explore() (e.g. a rule-0 directory-claim
    collision surfaced as an exception) propagates verbatim -- explore_pdep_network() must not wrap
    adapter.explore() in a catch-all that would relabel it as a 'failed' result."""
    _make_fake_factory(monkeypatch, raise_error=RuntimeError('rule-0 collision'))
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')

    with pytest.raises(RuntimeError, match='rule-0 collision'):
        explore_pdep_network(network_path=NETWORK_PATH, config=config)


def test_explore_pdep_network_reports_an_adapter_that_fails_without_saying_why(monkeypatch, tmp_path):
    """
    An adapter returning False with no reasons must yield a FAILED result naming the contract breach.

    Would catch two opposite mistakes. Letting the empty tuple through reaches
    PDepExplorationResult's empty-reasons guard, which raises a ValueError that reads as though this
    function built a malformed result -- burying the real fact, which is that an adapter broke the
    contract documented on PESExplorerAdapter.reasons. Silently substituting a bland "exploration
    failed" would lose it just as thoroughly. The exploration genuinely failed and the result must
    say so, while naming the adapter that could not diagnose itself.
    """
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    network_path = NETWORK_PATH
    _make_fake_factory(monkeypatch, succeed=False, reasons=tuple())

    result = explore_pdep_network(network_path=network_path, config=config)

    assert result.status == EXPLORATION_STATUS_FAILED
    assert len(result.reasons) == 1
    assert '_FakeAdapter' in result.reasons[0]
    assert 'contract violation' in result.reasons[0]


def test_explore_pdep_network_reads_networks_through_the_abc_contract_not_arkanes_attribute(
        monkeypatch, tmp_path):
    """
    The success path must read artifacts via get_networks(), the ABC's method.

    Would catch a regression to ``adapter.final_network_paths``, which is ArkaneExplorerAdapter's
    own state and appears nowhere in PESExplorerAdapter. Reading it makes the factory's pluggability
    a fiction: a second explorer implementing every abstract method would satisfy the ABC and then
    raise AttributeError here. This double has no such attribute, so the coupling cannot come back
    unnoticed.
    """
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    network_path = NETWORK_PATH
    _make_fake_factory(monkeypatch, network_paths=('/fake/from_get_networks.py',))

    result = explore_pdep_network(network_path=network_path, config=config)

    assert result.status == EXPLORATION_STATUS_SUCCEEDED
    assert result.network_paths == ('/fake/from_get_networks.py',)


def test_explore_pdep_network_copies_the_adapters_manifest(monkeypatch, tmp_path):
    """
    Mutating the adapter's manifest afterwards must not rewrite the reported result.

    The result advertises itself as a frozen record of what happened. Aliasing the adapter's live
    dict would let anything still holding the adapter edit the provenance of a run already reported
    -- fake immutability, the same defect class as a frozen config holding a caller's mutable dict.
    """
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    network_path = NETWORK_PATH
    constructed = []

    def _capturing_factory(**kwargs):
        adapter = _FakeAdapter()
        constructed.append(adapter)
        return adapter

    monkeypatch.setattr(t3_pdep_api, 'explorer_factory', _capturing_factory)

    result = explore_pdep_network(network_path=network_path, config=config)
    assert result.manifest == {'fake': True}

    constructed[0].manifest['fake'] = 'tampered'
    constructed[0].manifest['injected'] = True
    assert result.manifest == {'fake': True}


def test_explore_pdep_network_refuses_a_selection_that_was_never_evaluated(monkeypatch, tmp_path):
    """
    A not-evaluated selection must be REFUSED as a gate, not read as "does not qualify".

    This is the defect the test exists for. When select_pdep_network() rejects a stale SA cache it
    returns qualified=False AND evaluation_status='not_evaluated', and reason() still renders the
    full "does not qualify: no transition state the network is sensitive to ..." sentence. Gating on
    `qualified` alone therefore turns a missing evaluation into a negative verdict: the exploration
    is silently skipped, and the caller is handed a decision that was never made, phrased as though
    it had been. Raising is deliberate -- exploring anyway and skipping are both guesses about what
    the caller meant.
    """
    _make_fake_factory(monkeypatch)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=False, method='CSE',
                                     network_source_hash=NETWORK_SOURCE_HASH)
    selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED

    with pytest.raises(ValueError) as exc_info:
        explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection)
    message = str(exc_info.value)
    assert 'not_evaluated' in message
    assert 'carries no verdict' in message


def test_explore_pdep_network_reports_every_reason_a_selection_cannot_gate(monkeypatch, tmp_path):
    """Test that a selection which is BOTH stale (its recorded hash no longer matches the network
    file) AND unevaluated (evaluation_status != 'evaluated', qualified=False) is refused with a
    SINGLE error that reports BOTH diagnoses, not just whichever guard happened to run first. Each
    of these two problems is independently sufficient to invalidate the gate, and each is diagnosed
    by different evidence (a content hash vs. a coverage flag) -- silently discarding one diagnosis
    because the other fired first would send the caller off to fix one problem, re-run, and hit the
    second unannounced problem on the very next attempt."""
    _make_fake_factory(monkeypatch)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=False, method='CSE',
                                     network_source_hash='sha256:' + 'e' * 64)
    selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED

    with pytest.raises(ValueError) as exc_info:
        explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection)
    message = str(exc_info.value)
    # Distinctive substrings from the hash-mismatch diagnosis:
    assert 'has changed since' in message
    # Distinctive substrings from the evaluation-status diagnosis:
    assert 'not_evaluated' in message
    assert 'carries no verdict' in message


def test_explore_pdep_network_reports_method_and_network_id_alongside_evaluation_status(monkeypatch, tmp_path):
    """Test that a selection failing method, network_id, AND evaluation-status simultaneously is
    refused with a SINGLE error naming all three -- proving method/network_id are folded into the
    same accumulator as the hash/evaluation-status checks (187b28f only accumulated the latter two;
    method/network_id used to raise immediately on whichever was checked first, silently discarding
    every other reason). network_id already proves this is a different network entirely, so the
    hash-mismatch diagnosis (which reads "the network file has changed since the decision was made")
    would be misleading here and must be short-circuited -- its wording must NOT appear."""
    _make_fake_factory(monkeypatch)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='some_other_network', qualified=False, method='MSC',
                                     network_source_hash='sha256:' + 'e' * 64)
    selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED

    with pytest.raises(ValueError) as exc_info:
        explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection)
    message = str(exc_info.value)
    assert 'for 3 independent reasons' in message
    assert 'method' in message
    assert 'network_id' in message
    assert 'not_evaluated' in message
    assert 'carries no verdict' in message
    # Short-circuited: network_id already proves it, so the (misleading) hash-mismatch wording must
    # not also appear.
    assert 'has changed since' not in message


def test_explore_pdep_network_carries_output_paths_on_failure(monkeypatch, tmp_path):
    """
    A failed run must carry output_paths, and must NOT carry network_paths or k(T,P).

    'reasons' says what T3 concluded; output_paths say where Arkane's own logs and partial output
    are. Dropping them leaves whoever is diagnosing the failure to rediscover the run directory by
    hand, and this result is the only record they get. network_paths/k_tp are withheld for the
    opposite reason: they would assert a usable result exists, which is exactly what the failure
    denies.
    """
    _make_fake_factory(monkeypatch, succeed=False, reasons=('Arkane exploration did not converge.',))
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')

    result = explore_pdep_network(network_path=NETWORK_PATH, config=config)

    assert result.status == EXPLORATION_STATUS_FAILED
    assert result.output_paths == ('/fake/output1.py',)
    assert result.network_paths == tuple()
    assert result.k_tp_as_written == tuple()


# --- 8b. explore_pdep_network()'s admission policy -----------------------------------------------
#
# Qualification was originally a GATE: a selection that did not qualify meant "do not explore". It is
# now a RANKING input, and the spend decision can legitimately be made by the caller instead (T3 makes
# it in t3.main via t3.pdep.budget.apply_pdep_qm_budget). ``admission_policy`` is how a caller states
# WHICH of those two happened, and the tests below pin the asymmetry it creates: the qualification
# checks are policy and step aside under 'caller_admitted', while the method/network_id/hash checks
# are provenance and never do.

def test_explore_pdep_network_runs_an_unqualified_selection_when_the_caller_admitted_it(tmp_path, monkeypatch):
    """Test that an EVALUATED, non-qualifying selection explores under 'caller_admitted' instead of
    being skipped. This is the one cell of the gate's truth table the ranking reframe changes: the
    caller has ranked this network against the whole field and decided to spend on it anyway, so
    `qualified` is a tier here, not a veto. The selection is still passed -- and still fully checked
    below -- because it is what binds the run to the content the decision was made about."""
    calls = _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=False, method='CSE',
                                     network_source_hash=NETWORK_SOURCE_HASH)
    selection.evaluation_status = EVALUATION_STATUS_EVALUATED

    result = explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection,
                                  admission_policy=ADMISSION_POLICY_CALLER_ADMITTED)

    assert len(calls) == 1  # the adapter WAS constructed; under the default policy it would not be
    assert result.status == EXPLORATION_STATUS_SUCCEEDED
    # The result must say what admitted it. Without this, a succeeded result carrying an unqualified
    # selection is unreadable after the fact: it looks exactly like the gate having been bypassed.
    assert result.admission_policy == ADMISSION_POLICY_CALLER_ADMITTED


def test_explore_pdep_network_runs_a_never_evaluated_selection_when_the_caller_admitted_it(tmp_path, monkeypatch):
    """Test that 'caller_admitted' also stands down the not_evaluated refusal, not only the skip.

    That refusal exists because `qualified` carries no verdict unless the decision was evaluated, and
    reading a missing evaluation as a negative one is silent corruption. But its whole subject is
    using `qualified` AS A GATE -- and a caller declaring it admitted the network elsewhere is not
    using it as one. Keeping the raise here would force that caller back to selection=None, which
    throws away the method/network_id/hash binding that is the entire reason to pass a selection.
    That is the over-refusal failure mode, and it is the one this branch keeps rediscovering.
    """
    calls = _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=False, method='CSE',
                                     network_source_hash=NETWORK_SOURCE_HASH)
    selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED

    result = explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection,
                                  admission_policy=ADMISSION_POLICY_CALLER_ADMITTED)

    assert len(calls) == 1
    assert result.status == EXPLORATION_STATUS_SUCCEEDED


def test_explore_pdep_network_still_refuses_a_stale_selection_when_the_caller_admitted_it(tmp_path, monkeypatch):
    """Test that 'caller_admitted' stands down the QUALIFICATION checks and nothing else. A caller
    can decide to spend on a network it ranked; it cannot decide that a decision made about
    different bytes describes this run. Admitting a network is a budget statement, and no budget
    statement makes a stale selection current -- so the content binding stays unconditional."""
    calls = _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=True, method='CSE',
                                     network_source_hash='sha256:' + 'e' * 64)

    with pytest.raises(ValueError) as exc_info:
        explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection,
                             admission_policy=ADMISSION_POLICY_CALLER_ADMITTED)
    assert 'has changed since' in str(exc_info.value)
    assert calls == []


def test_explore_pdep_network_refuses_an_unknown_admission_policy(tmp_path, monkeypatch):
    """Test that a misspelled policy is refused rather than silently treated as the default. Falling
    back to the default would silently REINSTATE the gate for a caller who explicitly asked for it to
    stand aside, so the network they ranked and paid for would be skipped and reported as declined."""
    _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')

    with pytest.raises(ValueError, match='admission_policy'):
        explore_pdep_network(network_path=NETWORK_PATH, config=config,
                             admission_policy='caller-admitted')


def test_explore_pdep_network_refuses_caller_admitted_without_a_selection(tmp_path, monkeypatch):
    """Test that 'caller_admitted' with no selection is refused rather than quietly inert.

    The policy's entire job is to keep the provenance binding while overriding the qualification
    verdict. With selection=None there is no binding to keep, so the parameter would do nothing at
    all -- and the most likely way to write this call is by forgetting the selection you meant to
    pass, which loses exactly the binding you were reaching for. Note the asymmetry with the
    over-refusal the previous tests guard against: refusing here rejects an incoherent CALL, not a
    run whose data supports it, and the message names the escape (drop the argument).
    """
    calls = _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')

    with pytest.raises(ValueError) as exc_info:
        explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=None,
                             admission_policy=ADMISSION_POLICY_CALLER_ADMITTED)
    message = str(exc_info.value)
    assert 'selection' in message
    assert calls == []


def test_explore_pdep_network_records_the_default_admission_policy_on_a_skipped_result(tmp_path, monkeypatch):
    """Test that the default policy still skips a non-qualifying selection AND says so on the result.
    'skipped' alone does not distinguish "the qualification gate declined this" from any other
    not-run outcome; the recorded policy is what makes the pair readable."""
    calls = _make_fake_factory(monkeypatch)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=False, method='CSE',
                                     network_source_hash=NETWORK_SOURCE_HASH)

    result = explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection)

    assert result.status == EXPLORATION_STATUS_SKIPPED
    assert result.admission_policy == ADMISSION_POLICY_QUALIFIED_SELECTION
    assert calls == []


def test_explore_pdep_network_records_ungated_when_no_selection_was_given(tmp_path, monkeypatch):
    """Test that a selection-less exploration records 'ungated', not the argument's default.

    admission_policy answers "what admitted this run". With selection=None nothing did -- so
    recording 'qualified_selection' would have the result assert that a qualified selection admitted
    an exploration for which no selection ever existed. That is a false provenance claim, and the
    most consequential kind: it is exactly the claim someone auditing an expensive QM run would lean
    on. The value is therefore DERIVED from whether a selection was given, never echoed from the
    argument.
    """
    calls = _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')

    result = explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=None)

    assert len(calls) == 1
    assert result.status == EXPLORATION_STATUS_SUCCEEDED
    assert result.admission_policy == ADMISSION_POLICY_UNGATED


def test_explore_pdep_network_refuses_ungated_as_a_requested_policy(tmp_path, monkeypatch):
    """Test that 'ungated' cannot be REQUESTED, only derived. It is not a policy a caller chooses --
    it is what passing no selection means -- so accepting it as an argument would let a caller claim
    'ungated' while passing a selection, i.e. state two contradictory things about one run."""
    _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')

    with pytest.raises(ValueError, match='admission_policy'):
        explore_pdep_network(network_path=NETWORK_PATH, config=config,
                             admission_policy=ADMISSION_POLICY_UNGATED)


def test_explore_pdep_network_refuses_a_provenance_mismatch_even_when_the_caller_admitted_it(tmp_path, monkeypatch):
    """Test that method and network_id -- not only the content hash -- stay enforced under
    'caller_admitted'. Each of the three is a different way for a decision to be about something
    other than this run, and admitting a network says nothing about any of them."""
    calls = _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='some_other_network', qualified=False, method='MSC',
                                     network_source_hash=NETWORK_SOURCE_HASH)

    with pytest.raises(ValueError) as exc_info:
        explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection,
                             admission_policy=ADMISSION_POLICY_CALLER_ADMITTED)
    message = str(exc_info.value)
    assert 'method' in message
    assert 'network_id' in message
    assert calls == []


def test_explore_pdep_network_refuses_an_unreadable_evaluation_status_under_every_policy(tmp_path, monkeypatch):
    """Test that an evaluation_status outside the known set is refused even under 'caller_admitted'.

    That is malformed DATA, not a verdict this function may decline to consult: the policy stands
    down the qualification checks, and an unreadable coverage flag is not one. Without this check the
    'caller_admitted' path would carry a typo'd or invented status straight through to the recorded
    provenance of an expensive run, unexamined.
    """
    calls = _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=True, method='CSE',
                                     network_source_hash=NETWORK_SOURCE_HASH)
    selection.evaluation_status = 'evaluted'  # a typo, not a status

    with pytest.raises(ValueError) as exc_info:
        explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection,
                             admission_policy=ADMISSION_POLICY_CALLER_ADMITTED)
    assert 'evaluation_status' in str(exc_info.value)
    assert calls == []


def test_explore_pdep_network_records_the_admission_policy_on_a_failed_run(tmp_path, monkeypatch):
    """Test that a run which was admitted and then FAILED still records what admitted it. A failure
    is exactly when someone goes back to ask why this network was explored at all, so the failed path
    is the one where dropping the field costs the most -- and it is a separate construction site from
    the succeeded one, so it can rot independently."""
    _make_fake_factory(monkeypatch, succeed=False, reasons=('Arkane exploration did not converge.',))
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='CSE')
    selection = PDepNetworkSelection(network_id='network4_2', qualified=False, method='CSE',
                                     network_source_hash=NETWORK_SOURCE_HASH)
    selection.evaluation_status = EVALUATION_STATUS_EVALUATED

    result = explore_pdep_network(network_path=NETWORK_PATH, config=config, selection=selection,
                                  admission_policy=ADMISSION_POLICY_CALLER_ADMITTED)

    assert result.status == EXPLORATION_STATUS_FAILED
    assert result.admission_policy == ADMISSION_POLICY_CALLER_ADMITTED


def test_load_pdep_exploration_results_round_trips_the_admission_policy(tmp_path):
    """Test that admission_policy survives save -> load.

    It is recorded precisely so someone reading the results file later can tell why an unqualified
    network was explored, and that reader is downstream of the SAVE. as_dict() writing the key while
    the loader ignores it would leave every reloaded record claiming the default policy -- the
    reconstruction silently rewriting the one fact the field exists to preserve. A round-trip over a
    NON-default value is the only assertion that catches it; the existing round-trip tests all use
    the default and pass either way.
    """
    result = PDepExplorationResult(
        network_id='network4_2', status=EXPLORATION_STATUS_FAILED,
        reasons=('Arkane exploration did not converge.',),
        output_paths=('/fake/output1.py',),
        manifest={}, selection=_build_full_selection(),
        admission_policy=ADMISSION_POLICY_CALLER_ADMITTED,
    )
    path = str(tmp_path / 'exploration_results.yml')
    save_pdep_exploration_results(path=path, results=[result])

    loaded = load_pdep_exploration_results(path=path)

    assert loaded == [result]
    assert loaded[0].admission_policy == ADMISSION_POLICY_CALLER_ADMITTED


@pytest.mark.parametrize('with_selection, expected_policy', [
    (True, ADMISSION_POLICY_QUALIFIED_SELECTION),
    (False, ADMISSION_POLICY_UNGATED),
])
def test_load_pdep_exploration_results_derives_the_policy_for_a_record_written_before_the_field(
        tmp_path, with_selection, expected_policy):
    """Test that a record predating admission_policy is DERIVED, not defaulted and not refused.

    Such a record predates 'caller_admitted' entirely, so its basis is not a guess: carrying a
    selection means the qualification gate admitted it, carrying none means it was ungated. Refusing
    instead would make those files unloadable while still calling them
    exploration_result_schema_version 1 -- two incompatible shapes under one version number, which is
    what that version exists to prevent. A blanket default would be worse still: it would relabel
    selection-less records as gate-admitted, the exact false claim the field was added to stop.
    """
    result = PDepExplorationResult(
        network_id='network4_2', status=EXPLORATION_STATUS_FAILED,
        reasons=('Arkane exploration did not converge.',),
        output_paths=('/fake/output1.py',),
        selection=_build_full_selection() if with_selection else None,
    )
    path = str(tmp_path / 'exploration_results.yml')
    save_pdep_exploration_results(path=path, results=[result])
    # Strip the key back out, reproducing a file written before the field existed.
    content = read_yaml_file(path)
    del content['results'][0]['admission_policy']
    save_yaml_file(path=path, content=content)

    loaded = load_pdep_exploration_results(path=path)

    assert loaded[0].admission_policy == expected_policy


# --- 9. load_pdep_network_selections / load_pdep_exploration_results ----------------------------
#
# The two loaders are the read side of save_pdep_network_selections()/save_pdep_exploration_results()
# above. Both are STRICT (no migration path for an older on-disk shape -- see the module comment in
# t3/pdep/api.py above _sensitive_transition_state_from_dict for why that is safe to assert). These
# tests cover: round-trip equality (built from live objects, not hand-written expected dicts), the
# documented manifest-tuple-normalization caveat, network_source_hash preservation specifically, one
# distinctive-wording test per refusal branch, and a reconstruction-mutation sanity target.

def _build_full_selection():
    """A PDepNetworkSelection with every field set to a non-default value, including nested
    SensitiveTransitionState entries in both selected_ts and uncertain_path_reactions with
    non-default (tuple) conditions -- for the round-trip equality test below."""
    ts_a = SensitiveTransitionState(
        ts_label='TS3', coefficient=0.5, condition=(1000.0, 'K', 1.0, 'bar'),
        path_reaction_label='rxn1', path_reaction_str='A + B <=> C', kinetics_comment='estimate',
        uncertain=True, delta_ln_k=0.05)
    ts_b = SensitiveTransitionState(
        ts_label='TS7', coefficient=-0.3, condition=(1200.0, 'K', 2.0, 'bar'),
        path_reaction_label='rxn2', path_reaction_str='D + E <=> F', kinetics_comment='library',
        uncertain=False, delta_ln_k=0.03)
    return PDepNetworkSelection(
        network_id='network4_2',
        network_source_hash='sha256:deadbeef',
        qualified=True,
        network_reaction='A + B <=> C',
        direction_key='A + B <=> C',
        direction_keys=['A + B <=> C', 'D + E <=> F'],
        direction_ambiguous=True,
        method='MSC',
        sa_path='/some/path/sa.yml',
        cache_status=CACHE_STATUS_CACHED_VALID,
        thresholds={'relative_threshold': 0.01, 'min_delta_ln_k': 0.001, 'perturbation': 1000.0},
        selected_ts=[ts_a, ts_b],
        uncertain_path_reactions=[ts_a],
        warnings=['a warning'],
        network_reactions_examined=2,
        evaluation_status=EVALUATION_STATUS_EVALUATED,
        selection_schema_version=SELECTION_SCHEMA_VERSION,
        selection_algorithm_version=SELECTION_ALGORITHM_VERSION,
        t_grid_clamp={'clamped': True, 'requested_t_max': 3200.0, 'thermo_ceiling': 3000.0,
                     'written_t_max': 3000.0, 'tlist_dropped': False, 'tlist_original_highest': None,
                     'skipped_species': []},
    )


def test_load_pdep_network_selections_round_trip_equality(tmp_path):
    """Test that load(save(x)) equals x for a fully-populated PDepNetworkSelection, derived from the
    live object rather than a hand-written expected dict -- so the test cannot silently pass by
    omitting a field neither side checks."""
    selection = _build_full_selection()
    path = str(tmp_path / 'selections.yml')
    save_pdep_network_selections(path=path, selections=[selection])

    loaded = load_pdep_network_selections(path=path)

    assert loaded == [selection]


def test_load_pdep_network_selections_preserves_network_source_hash(tmp_path):
    """Test that network_source_hash survives load(save(x)) specifically: this is the field
    explore_pdep_network() refuses a selection for lacking (see
    test_explore_pdep_network_refuses_a_selection_with_no_source_hash above), so a loader that
    silently dropped or mangled it would defeat that guard for any selection read back from disk."""
    selection = PDepNetworkSelection(network_id='network4_2', network_source_hash='sha256:abc123')
    path = str(tmp_path / 'selections.yml')
    save_pdep_network_selections(path=path, selections=[selection])

    loaded = load_pdep_network_selections(path=path)

    assert loaded[0].network_source_hash == 'sha256:abc123'


def _build_full_arkane_reaction():
    """A PDepArkaneReaction with every field set to a non-default value. kinetics_params
    deliberately contains no tuples: PDepArkaneReaction.as_dict() passes kinetics_params through
    to_json_safe just like PDepExplorationResult.as_dict() does for manifest, so a tuple nested in
    either can never round-trip to an equal object (see
    test_save_pdep_exploration_results_normalizes_manifest_tuples_to_lists below, which documents
    that caveat for manifest; the same reasoning applies here)."""
    return PDepArkaneReaction(
        reactants=('A', 'B'),
        products=('C',),
        kinetics_type='Chebyshev',
        kinetics_params={'coeffs': [[1.0, 2.0], [3.0, None]], 'Tmin': 300, 'Tmin_unit': 'K'},
        numeric_values=(1.0, 2.0, 3.0, None, 300),
        rate_payload_numeric_values=(1.0, 2.0, 3.0, None),
        missing_kinetics_keys=('Pmin',),
    )


def test_load_pdep_exploration_results_round_trip_equality_succeeded(tmp_path):
    """Test that load(save(x)) equals x for a fully-populated, status='succeeded'
    PDepExplorationResult: non-empty network_paths/output_paths/k_tp_as_written, a nested
    PDepArkaneReaction with every field set, a nested selection, and a manifest without tuples
    (tuples inside manifest are a documented, separate caveat -- see the dedicated test below)."""
    selection = _build_full_selection()
    reaction = _build_full_arkane_reaction()
    result = PDepExplorationResult(
        network_id='network4_2', status=EXPLORATION_STATUS_SUCCEEDED,
        reasons=(),
        network_paths=('pdep/final/network1_full.py', 'pdep/final/network2_full.py'),
        output_paths=('/fake/output1.py',),
        k_tp_as_written=(reaction,),
        manifest={'note': 'ok', 'seeds': [1, 2, 3], 'meta': {'k': 'v'}},
        selection=selection,
    )
    path = str(tmp_path / 'exploration_results.yml')
    save_pdep_exploration_results(path=path, results=[result])

    loaded = load_pdep_exploration_results(path=path)

    assert loaded == [result]


def test_load_pdep_exploration_results_round_trip_equality_failed_with_reasons(tmp_path):
    """Test that load(save(x)) equals x for a status='failed' PDepExplorationResult with a non-empty
    reasons tuple and a non-empty output_paths tuple -- the two fields the succeeded-status test
    above cannot exercise non-trivially, since 'succeeded' requires reasons to stay empty."""
    result = PDepExplorationResult(
        network_id='network4_2', status=EXPLORATION_STATUS_FAILED,
        reasons=('Arkane exploration did not converge.', 'timed out'),
        output_paths=('/fake/output1.py',),
        manifest={}, selection=None,
    )
    path = str(tmp_path / 'exploration_results.yml')
    save_pdep_exploration_results(path=path, results=[result])

    loaded = load_pdep_exploration_results(path=path)

    assert loaded == [result]


def test_save_pdep_exploration_results_normalizes_manifest_tuples_to_lists(tmp_path):
    """Document (not merely tolerate) that PDepExplorationResult.as_dict() passes manifest through
    to_json_safe, which converts tuples anywhere inside it to lists -- so a manifest containing a
    tuple can never round-trip to an equal object. This is a known and intended limitation of the
    on-disk format, not a loader bug: assert the normalized (list) shape explicitly, and assert the
    loaded object is NOT equal to the original, rather than choosing a tuple-free fixture everywhere
    else in this file that would hide the limitation entirely."""
    result = PDepExplorationResult(
        network_id='network4_2', status=EXPLORATION_STATUS_SKIPPED,
        reasons=('skipped',), manifest={'nested': [(1, 2), {'deep': (3, 4)}]}, selection=None,
    )
    path = str(tmp_path / 'exploration_results.yml')
    save_pdep_exploration_results(path=path, results=[result])

    loaded = load_pdep_exploration_results(path=path)

    assert loaded[0].manifest == {'nested': [[1, 2], {'deep': [3, 4]}]}
    assert loaded[0] != result  # tuple normalization means it cannot equal the original


# --- 9a. load_pdep_network_selections refusal branches, one test per branch ---------------------

def test_load_pdep_network_selections_refuses_non_mapping_top_level(tmp_path):
    """Test the top-level-not-a-mapping refusal, on its own distinctive wording."""
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content=[1, 2, 3])

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    assert 'does not contain a mapping' in str(exc_info.value)


def test_load_pdep_network_selections_refuses_missing_version_key(tmp_path):
    """Test the absent-envelope-version-key refusal, on its own distinctive wording -- distinct from
    the unknown-version wording checked below (an unversioned file is not "version 1", it is
    unknown)."""
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selections': []})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    assert "no 'selection_schema_version' key" in str(exc_info.value)
    assert 'unknown shape' in str(exc_info.value)


def test_load_pdep_network_selections_refuses_unknown_version(tmp_path):
    """Test the unknown-envelope-version refusal, on its own distinctive wording."""
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': 999, 'selections': []})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    assert 'only understands version' in str(exc_info.value)


def test_load_pdep_network_selections_refuses_missing_selections_key(tmp_path):
    """Test the missing/non-list 'selections' key refusal."""
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    assert "no 'selections' list" in str(exc_info.value)


def test_load_pdep_network_selections_refuses_non_list_selections_key(tmp_path):
    """Test the 'selections' key being present but not a list -- the other half of the same refusal
    branch as the missing-key test above."""
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': 'not-a-list'})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    assert "no 'selections' list" in str(exc_info.value)


def test_load_pdep_network_selections_refuses_non_mapping_record(tmp_path):
    """Test the per-record not-a-mapping refusal."""
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [42]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    assert 'selection record that is not a mapping' in str(exc_info.value)


def test_load_pdep_network_selections_refuses_unknown_record_schema_version(tmp_path):
    """Test the per-record unknown-selection_schema_version refusal."""
    rendered = _build_full_selection().as_dict()
    rendered['selection_schema_version'] = 999
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    assert 'which this code does not understand' in str(exc_info.value)


# --- 9b. load_pdep_exploration_results refusal branches, one test per branch --------------------

def test_load_pdep_exploration_results_refuses_non_mapping_top_level(tmp_path):
    """Test the top-level-not-a-mapping refusal, on its own distinctive wording."""
    path = str(tmp_path / 'exploration_results.yml')
    save_yaml_file(path=path, content=[1, 2])

    with pytest.raises(ValueError) as exc_info:
        load_pdep_exploration_results(path=path)
    assert 'does not contain a mapping' in str(exc_info.value)


def test_load_pdep_exploration_results_refuses_missing_version_key(tmp_path):
    """Test the absent-envelope-version-key refusal, distinct wording from the unknown-version case
    checked below."""
    path = str(tmp_path / 'exploration_results.yml')
    save_yaml_file(path=path, content={'results': []})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_exploration_results(path=path)
    assert "no 'exploration_result_schema_version' key" in str(exc_info.value)
    assert 'unknown shape' in str(exc_info.value)


def test_load_pdep_exploration_results_refuses_unknown_version(tmp_path):
    """Test the unknown-envelope-version refusal, on its own distinctive wording."""
    path = str(tmp_path / 'exploration_results.yml')
    save_yaml_file(path=path, content={'exploration_result_schema_version': 999, 'results': []})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_exploration_results(path=path)
    assert 'only understands version' in str(exc_info.value)


def test_load_pdep_exploration_results_refuses_missing_results_key(tmp_path):
    """Test the missing/non-list 'results' key refusal."""
    path = str(tmp_path / 'exploration_results.yml')
    save_yaml_file(path=path,
                   content={'exploration_result_schema_version': EXPLORATION_RESULT_SCHEMA_VERSION})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_exploration_results(path=path)
    assert "no 'results' list" in str(exc_info.value)


def test_load_pdep_exploration_results_refuses_non_mapping_record(tmp_path):
    """Test the per-record not-a-mapping refusal."""
    path = str(tmp_path / 'exploration_results.yml')
    save_yaml_file(path=path, content={'exploration_result_schema_version': EXPLORATION_RESULT_SCHEMA_VERSION,
                                       'results': [42]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_exploration_results(path=path)
    assert 'result record that is not a mapping' in str(exc_info.value)


def test_load_pdep_exploration_results_refuses_unknown_nested_selection_version(tmp_path):
    """Test that a result carrying a nested selection this code does not understand is refused via
    the same per-record version check load_pdep_network_selections uses, proving the two loaders
    share -- not duplicate -- that logic."""
    result = PDepExplorationResult(network_id='network4_2', status=EXPLORATION_STATUS_SKIPPED,
                                   reasons=('x',), selection=PDepNetworkSelection(network_id='network4_2'))
    rendered = result.as_dict()
    rendered['selection']['selection_schema_version'] = 999
    path = str(tmp_path / 'exploration_results.yml')
    save_yaml_file(path=path, content={'exploration_result_schema_version': EXPLORATION_RESULT_SCHEMA_VERSION,
                                       'results': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_exploration_results(path=path)
    assert 'which this code does not understand' in str(exc_info.value)


def test_load_pdep_network_selections_refuses_null_record(tmp_path):
    """Test that a null (None) entry in the selections list is refused, not silently returned as
    None in the loaded list for a caller to later dereference."""
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [None]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    assert 'selections[0]' in str(exc_info.value)
    assert 'null' in str(exc_info.value)


def test_load_pdep_network_selections_refuses_record_missing_required_field(tmp_path):
    """Test that a record missing a required field raises a diagnostic ValueError naming the file,
    the record's index, and the missing field -- instead of a raw KeyError('network_reaction') out
    of a public API."""
    rendered = _build_full_selection().as_dict()
    del rendered['network_reaction']
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "missing required field 'network_reaction'" in message


def test_load_pdep_network_selections_refuses_a_string_for_a_list_typed_field(tmp_path):
    """Test that a list-typed field given a bare string is refused rather than silently coerced
    character by character (e.g. warnings: 'AB' silently becoming ('A', 'B'))."""
    rendered = _build_full_selection().as_dict()
    rendered['warnings'] = 'AB'
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'warnings' that must be a list" in message


def test_load_pdep_network_selections_refuses_unknown_record_algorithm_version(tmp_path):
    """Test the per-record unknown-selection_algorithm_version refusal, on wording distinguishable
    from the unknown-selection_schema_version refusal above (both must be tellable apart)."""
    rendered = _build_full_selection().as_dict()
    rendered['selection_algorithm_version'] = 999
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert 'selection_algorithm_version=999' in message
    assert 'cannot interpret' in message
    assert 'which this code does not understand' not in message  # distinct from the schema-version wording


def test_load_pdep_network_selections_refuses_the_forged_budget_gate_bypass(tmp_path):
    """Named regression test for the exact forgery the budget gate in explore_pdep_network is
    vulnerable to without type validation: qualified: 'yes' is a non-empty STRING, therefore
    truthy, so a bare presence check lets it satisfy 'not selection.qualified'; combined with
    evaluation_status: 'bogus' (which is neither EVALUATION_STATUS_EVALUATED nor
    EVALUATION_STATUS_NOT_EVALUATED, so it does not equal EVALUATION_STATUS_EVALUATED and would
    also happen to satisfy that half of the guard), this hand-edited record would sail straight
    through the gate. Both must be refused at LOAD time so a forged record can never even become a
    PDepNetworkSelection object."""
    rendered = _build_full_selection().as_dict()
    rendered['qualified'] = 'yes'
    rendered['evaluation_status'] = 'bogus'
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'qualified'" in message


def test_load_pdep_network_selections_refuses_qualified_true_with_no_uncertain_evidence(tmp_path):
    """Named regression test for the cross-field forgery the type/enum checks in
    _selection_from_dict do not catch: a record with a well-typed bool qualified=True, a
    recognized evaluation_status, a correct method/network/hash, and empty
    uncertain_path_reactions. Neither select_from_sa_dict() (which sets
    qualified = bool(uncertain_path_reactions) for a single decision) nor combine() (whose
    qualified = any(...) over evaluated components is only ever True because the qualifying
    component's own non-empty uncertain_path_reactions is part of the unioned result) can ever
    produce qualified=True alongside an empty uncertain_path_reactions -- so this combination can
    only arise from a hand-edited file, and must be refused at LOAD time. Because the record is
    refused before load_pdep_network_selections returns, no PDepNetworkSelection object bearing
    this lie is ever constructed, so explore_pdep_network's budget gate can never be reached with
    it -- there is no selection object to hand it."""
    rendered = _build_full_selection().as_dict()
    rendered['qualified'] = True
    rendered['uncertain_path_reactions'] = []
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert 'network4_2' in message
    assert 'uncertain_path_reactions' in message


def test_load_pdep_network_selections_refuses_uncertain_path_reactions_entry_not_marked_uncertain(
        tmp_path):
    """Named regression test for a subtler forgery of the same hole: qualified=True with a
    NON-empty uncertain_path_reactions that nonetheless is not real evidence, because the entry it
    contains has uncertain=False (or was never marked at all). select_from_sa_dict() builds
    uncertain_path_reactions as exactly ``[entry for entry in selected_ts if entry.uncertain]``, so
    every real entry has uncertain is True by construction; a record whose uncertain_path_reactions
    entry fails that must be hand-edited, not decision output."""
    rendered = _build_full_selection().as_dict()
    rendered['qualified'] = True
    # ts_b (the second selected_ts entry in _build_full_selection) has uncertain=False -- swap it in
    # as the sole "uncertain" entry to forge evidence that was never actually uncertain.
    rendered['uncertain_path_reactions'] = [rendered['selected_ts'][1]]
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert 'network4_2' in message
    assert 'uncertain_path_reactions[0]' in message


def test_load_pdep_network_selections_refuses_uncertain_path_reactions_entry_absent_from_selected_ts(
        tmp_path):
    """Named regression test for a third variant of the same hole: qualified=True with a
    uncertain_path_reactions entry that is well-formed (uncertain=True) but does not appear in
    selected_ts at all -- evidence conjured from nothing rather than drawn from the network's own
    selected transition states. select_from_sa_dict() builds uncertain_path_reactions as a subset
    of selected_ts by construction (a list comprehension filtering selected_ts itself), and
    combine() only ever unions uncertain_path_reactions alongside the matching union of
    selected_ts, so a real decision can never produce an uncertain_path_reactions entry that is not
    also present in selected_ts. A record with such an entry must be hand-edited, not decision
    output, and is refused at LOAD time."""
    rendered = _build_full_selection().as_dict()
    fabricated_entry = dict(rendered['uncertain_path_reactions'][0])
    fabricated_entry['ts_label'] = 'TS_NOT_IN_SELECTED_TS'
    rendered['uncertain_path_reactions'] = [fabricated_entry]
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert 'network4_2' in message
    assert 'uncertain_path_reactions[0]' in message
    assert 'selected_ts' in message


def test_load_pdep_network_selections_accepts_a_qualified_but_partially_evaluated_selection(tmp_path):
    """Over-refusal guard: a record with qualified=True, evaluation_status='not_evaluated', and
    real, non-empty uncertain_path_reactions evidence must still LOAD. This is the documented
    "partial yes is supported" case (see PDepNetworkSelection.reason() and the gate in
    explore_pdep_network): select_from_sa_dict() and combine() both allow a positive verdict to
    stand on whatever evaluated evidence qualified it, even when other components/coverage never
    got evaluated. If the cross-field check added to _selection_from_dict() required
    evaluation_status == EVALUATION_STATUS_EVALUATED whenever qualified is True, this legitimate
    record would wrongly be refused -- exactly the over-refusal failure mode this test pins."""
    rendered = _build_full_selection().as_dict()
    assert rendered['qualified'] is True
    assert rendered['uncertain_path_reactions'], 'fixture must carry real evidence for this to test anything'
    rendered['evaluation_status'] = EVALUATION_STATUS_NOT_EVALUATED
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    loaded = load_pdep_network_selections(path=path)

    assert len(loaded) == 1
    assert loaded[0].qualified is True
    assert loaded[0].evaluation_status == EVALUATION_STATUS_NOT_EVALUATED


def test_load_pdep_network_selections_accepts_unqualified_with_evidence_from_a_combined_component(
        tmp_path):
    """Over-refusal guard for the asymmetry combine() bakes in: qualified is voted over EVALUATED
    components only, but selected_ts/uncertain_path_reactions are unioned over ALL components
    (t3/pdep/selector.py combine(), _union_preserving_order calls). So a not-evaluated component
    can contribute real uncertain evidence to the aggregate's uncertain_path_reactions while the
    aggregate's qualified stays False, because no EVALUATED component voted True. Built through the
    real PDepNetworkSelection.combine() classmethod, not assembled by hand, so this pins the actual
    production behaviour rather than an assumption about it. If the cross-field check required
    qualified == bool(uncertain_path_reactions) (the single-decision equality), this legitimate
    combined record would wrongly be refused."""
    ts_x = SensitiveTransitionState(
        ts_label='TSX', coefficient=0.4, condition=(900.0, 'K', 1.0, 'bar'), path_reaction_label='rxnX',
        path_reaction_str='X + Y <=> Z', kinetics_comment='estimate', uncertain=True, delta_ln_k=0.02)
    not_evaluated_component = PDepNetworkSelection(
        network_id='network4_2', network_source_hash=NETWORK_SOURCE_HASH, qualified=False,
        selected_ts=[ts_x], uncertain_path_reactions=[ts_x])
    not_evaluated_component.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
    evaluated_negative_component = PDepNetworkSelection(
        network_id='network4_2', network_source_hash=NETWORK_SOURCE_HASH, qualified=False)
    assert evaluated_negative_component.evaluation_status == EVALUATION_STATUS_EVALUATED

    combined = PDepNetworkSelection.combine([not_evaluated_component, evaluated_negative_component])
    assert combined.qualified is False, 'fixture setup no longer exercises the intended combine() branch'
    assert combined.uncertain_path_reactions == [ts_x], \
        'fixture setup no longer exercises the intended combine() branch'

    path = str(tmp_path / 'selections.yml')
    save_pdep_network_selections(path=path, selections=[combined])
    loaded = load_pdep_network_selections(path=path)

    assert len(loaded) == 1
    assert loaded[0].qualified is False
    assert loaded[0].uncertain_path_reactions == [ts_x]
    assert loaded[0].evaluation_status == EVALUATION_STATUS_NOT_EVALUATED


def test_load_pdep_network_selections_round_trips_a_real_select_from_sa_dict_decision(tmp_path, sa_dict):
    """Test that a PDepNetworkSelection produced by the real select_from_sa_dict() -- not a
    hand-built fixture that happens to agree with the cross-field validator -- survives
    save/load unchanged. Built through the real constructor so the validator is exercised against
    genuine decision output rather than against a fixture engineered to satisfy it."""
    network = parse_pdep_network_file(path=NETWORK_PATH)
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001)

    path = str(tmp_path / 'selections.yml')
    save_pdep_network_selections(path=path, selections=[selection])
    loaded = load_pdep_network_selections(path=path)

    assert len(loaded) == 1
    assert loaded[0] == selection


def test_load_pdep_network_selections_refuses_a_string_for_qualified(tmp_path):
    """Test that qualified: 'yes' -- a non-empty string, therefore truthy -- is refused rather than
    silently accepted as a bool, since a truthy non-bool here would forge its way past
    explore_pdep_network's 'not selection.qualified' budget gate."""
    rendered = _build_full_selection().as_dict()
    rendered['qualified'] = 'yes'
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'qualified' that must be a bool" in message


def test_load_pdep_network_selections_refuses_an_int_for_qualified(tmp_path):
    """Test that qualified: 1 is refused: bool is a subclass of int in Python, so a naive
    isinstance(value, int) check would wrongly accept this -- isinstance(value, bool) must be
    checked instead (and first)."""
    rendered = _build_full_selection().as_dict()
    rendered['qualified'] = 1
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'qualified' that must be a bool" in message


def test_load_pdep_network_selections_refuses_none_for_qualified(tmp_path):
    """Test that qualified: null is refused: PDepNetworkSelection.qualified is a plain bool field
    (never optional), so None must not be silently accepted."""
    rendered = _build_full_selection().as_dict()
    rendered['qualified'] = None
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'qualified' that must be a bool" in message


def test_load_pdep_network_selections_refuses_an_unrecognized_evaluation_status(tmp_path):
    """Test that evaluation_status: 'bogus' is refused rather than silently accepted as some
    unrecognized third state -- only EVALUATION_STATUS_EVALUATED/EVALUATION_STATUS_NOT_EVALUATED
    are understood."""
    rendered = _build_full_selection().as_dict()
    rendered['evaluation_status'] = 'bogus'
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'evaluation_status'" in message
    assert 'bogus' in message


def test_load_pdep_network_selections_refuses_an_unrecognized_cache_status(tmp_path):
    """Test that cache_status: 'bogus' is refused; cache_status is optional (None is a legitimate
    value when there is no cache to validate), but a non-null value must be one of the four
    recognized CACHE_STATUS_* constants."""
    rendered = _build_full_selection().as_dict()
    rendered['cache_status'] = 'bogus'
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'cache_status'" in message
    assert 'bogus' in message


def test_load_pdep_network_selections_accepts_a_null_cache_status(tmp_path):
    """Test that cache_status: null keeps loading rather than being refused: unlike
    evaluation_status, cache_status is a legitimately optional field."""
    rendered = _build_full_selection().as_dict()
    rendered['cache_status'] = None
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    loaded = load_pdep_network_selections(path=path)

    assert loaded[0].cache_status is None


def test_load_pdep_network_selections_reconstructs_a_present_t_grid_clamp(tmp_path):
    """Test that a record carrying a real t_grid_clamp mapping reconstructs it as a fresh dict equal
    to what was written -- the ordinary, fully-provenanced case."""
    rendered = _build_full_selection().as_dict()
    assert rendered['t_grid_clamp'] == {
        'clamped': True, 'requested_t_max': 3200.0, 'thermo_ceiling': 3000.0,
        'written_t_max': 3000.0, 'tlist_dropped': False, 'tlist_original_highest': None,
        'skipped_species': [],
    }
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    loaded = load_pdep_network_selections(path=path)

    assert loaded[0].t_grid_clamp == rendered['t_grid_clamp']


def test_load_pdep_network_selections_defaults_a_missing_t_grid_clamp_key_to_none(tmp_path):
    """Test that a record predating this field -- the key ABSENT entirely, not merely null --
    reconstructs t_grid_clamp as None (unknown provenance) rather than raising or defaulting to
    some other value. This is the exact scenario of an old on-disk selection written before
    t_grid_clamp existed; it must keep loading, not be refused."""
    rendered = _build_full_selection().as_dict()
    del rendered['t_grid_clamp']
    assert 't_grid_clamp' not in rendered
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    loaded = load_pdep_network_selections(path=path)

    assert loaded[0].t_grid_clamp is None
    assert loaded[0].evaluation_status == EVALUATION_STATUS_EVALUATED, (
        'a missing t_grid_clamp key must not gate evaluation_status')


def test_load_pdep_network_selections_accepts_an_explicit_null_t_grid_clamp(tmp_path):
    """Test that t_grid_clamp: null (explicit, as opposed to the key being absent) also
    reconstructs to None -- the same unknown-provenance outcome as the absent-key case, since a
    caller writing None explicitly and a caller from before the field existed must be
    indistinguishable to a reader."""
    rendered = _build_full_selection().as_dict()
    rendered['t_grid_clamp'] = None
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    loaded = load_pdep_network_selections(path=path)

    assert loaded[0].t_grid_clamp is None


def test_load_pdep_network_selections_refuses_a_non_mapping_t_grid_clamp(tmp_path):
    """Test that a t_grid_clamp value that is neither a mapping nor null is refused as a genuine
    malformation -- TGridClampRecord.as_dict() never renders anything else, so a string/list/number
    here means the file was hand-edited or corrupted, not merely old."""
    rendered = _build_full_selection().as_dict()
    rendered['t_grid_clamp'] = 'not-a-dict'
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError, match='t_grid_clamp'):
        load_pdep_network_selections(path=path)


def test_load_pdep_network_selections_refuses_a_non_string_network_id(tmp_path):
    """Test that network_id (a required, non-optional str) given a non-string value is refused."""
    rendered = _build_full_selection().as_dict()
    rendered['network_id'] = 123
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'network_id' that must be a string" in message


def test_load_pdep_network_selections_refuses_a_non_string_network_source_hash(tmp_path):
    """Test that network_source_hash (optional str) given a non-string, non-null value is refused
    -- this is the field explore_pdep_network binds its budget gate to, so a wrong-typed value here
    must never be silently accepted."""
    rendered = _build_full_selection().as_dict()
    rendered['network_source_hash'] = 123
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'network_source_hash' that must be a string or null" in message


def test_load_pdep_network_selections_refuses_a_non_string_method(tmp_path):
    """Test that method (optional str) given a non-string, non-null value is refused -- this is the
    field explore_pdep_network compares against config.method, so a wrong-typed value must never be
    silently accepted."""
    rendered = _build_full_selection().as_dict()
    rendered['method'] = 42
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'method' that must be a string or null" in message


def test_load_pdep_network_selections_refuses_a_bool_network_reactions_examined(tmp_path):
    """Test that network_reactions_examined: True is refused: bool is a subclass of int, so this
    field's isinstance(value, bool) check must run before (and instead of) a bare
    isinstance(value, int) check that would wrongly accept it."""
    rendered = _build_full_selection().as_dict()
    rendered['network_reactions_examined'] = True
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'network_reactions_examined' that must be an int" in message


def test_load_pdep_network_selections_refuses_a_string_network_reactions_examined(tmp_path):
    """Test that network_reactions_examined: '2' (a numeric-looking string) is refused rather than
    silently coerced."""
    rendered = _build_full_selection().as_dict()
    rendered['network_reactions_examined'] = '2'
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'network_reactions_examined' that must be an int" in message


def test_load_pdep_network_selections_refuses_a_non_bool_direction_ambiguous(tmp_path):
    """Test that direction_ambiguous given a non-bool value is refused."""
    rendered = _build_full_selection().as_dict()
    rendered['direction_ambiguous'] = 'yes'
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "field 'direction_ambiguous' that must be a bool" in message


def test_load_pdep_network_selections_refuses_a_non_numeric_threshold_value(tmp_path):
    """Test that a thresholds dict entry given a non-numeric (bool) value is refused, WITHOUT
    requiring a fixed key set -- select_from_sa_dict's malformed-sa_dict branch records a
    coefficient_floor key absent from every other thresholds dict, so pinning an exact key set
    would refuse a legitimately-produced record."""
    rendered = _build_full_selection().as_dict()
    rendered['thresholds']['relative_threshold'] = True
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0]' in message
    assert "'thresholds' entry 'relative_threshold' that must be numeric" in message


def test_load_pdep_network_selections_accepts_an_unlisted_threshold_key(tmp_path):
    """Test that a thresholds dict with an extra key not present on any other record (e.g.
    coefficient_floor) is accepted, confirming thresholds validation does not enforce a fixed key
    set."""
    rendered = _build_full_selection().as_dict()
    rendered['thresholds']['coefficient_floor'] = 1e-8
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    loaded = load_pdep_network_selections(path=path)

    assert loaded[0].thresholds['coefficient_floor'] == 1e-8


def test_load_pdep_network_selections_refuses_a_non_bool_uncertain_on_a_sensitive_ts(tmp_path):
    """Test that a selected_ts entry's 'uncertain' field given a non-bool, non-null value is
    refused -- sibling fix to _selection_from_dict's type validation, applied to
    _sensitive_transition_state_from_dict, naming both the record's index AND the nested list
    index/field so a malformed nested entry is diagnosable."""
    rendered = _build_full_selection().as_dict()
    rendered['selected_ts'][0]['uncertain'] = 'yes'
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0].selected_ts[0]' in message
    assert "field 'uncertain' that must be a bool or null" in message


def test_load_pdep_network_selections_refuses_a_non_numeric_coefficient_on_a_sensitive_ts(tmp_path):
    """Test that a selected_ts entry's 'coefficient' field given a non-numeric value is refused
    rather than raising an opaque KeyError-adjacent TypeError deeper in the dataclass."""
    rendered = _build_full_selection().as_dict()
    rendered['selected_ts'][0]['coefficient'] = 'not-a-number'
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0].selected_ts[0]' in message
    assert "'coefficient' that must be numeric" in message


def test_load_pdep_network_selections_refuses_a_missing_field_on_a_sensitive_ts(tmp_path):
    """Test that a selected_ts entry missing a required field (e.g. ts_label) raises a diagnostic
    ValueError naming the nested context and the field, instead of a raw KeyError('ts_label') --
    the same fix _require_record_field already gives the top-level selection fields, now applied
    consistently to the nested SensitiveTransitionState entries."""
    rendered = _build_full_selection().as_dict()
    del rendered['selected_ts'][0]['ts_label']
    path = str(tmp_path / 'selections.yml')
    save_yaml_file(path=path, content={'selection_schema_version': SELECTION_SCHEMA_VERSION,
                                       'selections': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_selections(path=path)
    message = str(exc_info.value)
    assert 'selections[0].selected_ts[0]' in message
    assert "missing required field 'ts_label'" in message


def test_load_pdep_exploration_results_refuses_a_non_string_kinetics_type(tmp_path):
    """Test that a k_tp_as_written entry's kinetics_type field given a non-string, non-null value
    is refused."""
    reaction = _build_full_arkane_reaction()
    result = PDepExplorationResult(network_id='network4_2', status=EXPLORATION_STATUS_SUCCEEDED,
                                   reasons=(), network_paths=('pdep/final/network4_2_full.py',),
                                   k_tp_as_written=(reaction,), selection=None)
    rendered = result.as_dict()
    rendered['k_tp_as_written'][0]['kinetics_type'] = 42
    path = str(tmp_path / 'exploration_results.yml')
    save_yaml_file(path=path, content={'exploration_result_schema_version': EXPLORATION_RESULT_SCHEMA_VERSION,
                                       'results': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_exploration_results(path=path)
    message = str(exc_info.value)
    assert 'results[0].k_tp_as_written[0]' in message
    assert "field 'kinetics_type' that must be a string or null" in message


def test_load_pdep_exploration_results_refuses_a_non_mapping_kinetics_params(tmp_path):
    """Test that a k_tp_as_written entry's kinetics_params field given a non-mapping value is
    refused."""
    reaction = _build_full_arkane_reaction()
    result = PDepExplorationResult(network_id='network4_2', status=EXPLORATION_STATUS_SUCCEEDED,
                                   reasons=(), network_paths=('pdep/final/network4_2_full.py',),
                                   k_tp_as_written=(reaction,), selection=None)
    rendered = result.as_dict()
    rendered['k_tp_as_written'][0]['kinetics_params'] = 'not-a-mapping'
    path = str(tmp_path / 'exploration_results.yml')
    save_yaml_file(path=path, content={'exploration_result_schema_version': EXPLORATION_RESULT_SCHEMA_VERSION,
                                       'results': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_exploration_results(path=path)
    message = str(exc_info.value)
    assert 'results[0].k_tp_as_written[0]' in message
    assert "field 'kinetics_params' that must be a mapping" in message


def test_load_pdep_exploration_results_refuses_record_missing_required_field(tmp_path):
    """Test that a result record missing a required field raises a diagnostic ValueError naming the
    file, the record's index, and the missing field."""
    result = PDepExplorationResult(network_id='network4_2', status=EXPLORATION_STATUS_SKIPPED,
                                   reasons=('x',), selection=None)
    rendered = result.as_dict()
    del rendered['network_id']
    path = str(tmp_path / 'exploration_results.yml')
    save_yaml_file(path=path, content={'exploration_result_schema_version': EXPLORATION_RESULT_SCHEMA_VERSION,
                                       'results': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_exploration_results(path=path)
    message = str(exc_info.value)
    assert 'results[0]' in message
    assert "missing required field 'network_id'" in message


def test_load_pdep_exploration_results_refuses_a_string_for_a_list_typed_field(tmp_path):
    """Test that a list-typed field on a result record given a bare string is refused rather than
    silently coerced character by character."""
    result = PDepExplorationResult(network_id='network4_2', status=EXPLORATION_STATUS_FAILED,
                                   reasons=('x',), output_paths=('/fake/output1.py',), selection=None)
    rendered = result.as_dict()
    rendered['output_paths'] = 'AB'
    path = str(tmp_path / 'exploration_results.yml')
    save_yaml_file(path=path, content={'exploration_result_schema_version': EXPLORATION_RESULT_SCHEMA_VERSION,
                                       'results': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_exploration_results(path=path)
    message = str(exc_info.value)
    assert 'results[0]' in message
    assert "field 'output_paths' that must be a list" in message


def test_load_pdep_exploration_results_refuses_unknown_nested_selection_algorithm_version(tmp_path):
    """Test that a result carrying a nested selection with an unknown selection_algorithm_version is
    refused, on wording distinguishable from the unknown-schema-version refusal."""
    result = PDepExplorationResult(network_id='network4_2', status=EXPLORATION_STATUS_SKIPPED,
                                   reasons=('x',), selection=PDepNetworkSelection(network_id='network4_2'))
    rendered = result.as_dict()
    rendered['selection']['selection_algorithm_version'] = 999
    path = str(tmp_path / 'exploration_results.yml')
    save_yaml_file(path=path, content={'exploration_result_schema_version': EXPLORATION_RESULT_SCHEMA_VERSION,
                                       'results': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_exploration_results(path=path)
    message = str(exc_info.value)
    assert 'results[0].selection' in message
    assert 'cannot interpret' in message


def test_load_pdep_exploration_results_refuses_a_result_missing_the_selection_key(tmp_path):
    """Test that a result record with NO 'selection' key at all is refused, distinct from a record
    carrying an explicit 'selection: null' (which is the documented "explored without a gating
    decision" state and must keep loading -- see
    test_save_pdep_exploration_results_handles_selection_none above)."""
    result = PDepExplorationResult(network_id='network4_2', status=EXPLORATION_STATUS_SKIPPED,
                                   reasons=('x',), selection=None)
    rendered = result.as_dict()
    del rendered['selection']
    path = str(tmp_path / 'exploration_results.yml')
    save_yaml_file(path=path, content={'exploration_result_schema_version': EXPLORATION_RESULT_SCHEMA_VERSION,
                                       'results': [rendered]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_exploration_results(path=path)
    message = str(exc_info.value)
    assert 'results[0]' in message
    assert "missing the required 'selection' key" in message


def test_a_cache_rejected_selection_is_refused_for_the_RIGHT_reason(tmp_path, monkeypatch):
    """
    Test that a selection from select_pdep_network's cache-rejected path is refused by
    explore_pdep_network for being UNEVALUATED, not for carrying no content hash.

    Both refusals are correct behaviour and both raise ValueError, so a test asserting merely that
    ValueError is raised would pass either way. The distinction is the whole point: the hash check
    runs BEFORE the evaluation check, so if the cache-rejected construction site ever stopped
    recording network_source_hash, this call would still raise -- but with a message asserting the
    network file changed under a decision that never recorded one, blaming the wrong thing and
    sending a reader to the wrong file. Mutation-checked: deleting the kwarg at api.py:165 left the
    entire suite green before this test existed.
    """
    _make_fake_factory(monkeypatch, succeed=True)
    trusted_root = str(tmp_path / 'root')
    os.makedirs(trusted_root)
    config = _make_config(trusted_root, os.path.join(trusted_root, 'run1'), method='MSC')

    network_path = str(tmp_path / 'network4_2.py')
    with open(NETWORK_PATH, 'r') as f:
        _write(network_path, f.read())
    sa_path = str(tmp_path / 'sa_coefficients.yml')
    save_yaml_file(path=sa_path, content=read_yaml_file(path=SA_PATH))

    # No sidecar is written, so the cache is rejected and the selection is never evaluated.
    rejected = select_pdep_network(network=network_path, sa_path=sa_path,
                                   network_reaction=TARGET_REACTION, relative_threshold=0.001,
                                   method='MSC')
    assert rejected.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED

    with pytest.raises(ValueError) as excinfo:
        explore_pdep_network(network_path=network_path, config=config, selection=rejected)

    message = str(excinfo.value)
    assert 'carries no verdict' in message, message
    assert 'carries no binding to the content' not in message, message


def test_the_loaders_refuse_a_python_object_tag(tmp_path):
    """
    Test that the public loaders cannot be made to construct a Python object from a YAML tag.

    `t3/pdep/yaml_safe.py`'s module docstring already states the rule these loaders were written
    against and broke: `t3.pdep.api` is a PUBLIC entrypoint reading a CALLER-SUPPLIED path, so a
    `yaml.FullLoader` read (which is what `arc.common.read_yaml_file` does) hands whoever influences
    that path control over what gets constructed. Nothing T3 writes needs any tag at all --
    `as_dict()` renders tuples as plain lists -- so a tag here means the file was not written by T3.

    Both loaders are checked: routing only one of them through the safe reader would leave the hole
    open, and no other test distinguishes them.
    """
    payload = ("selection_schema_version: 1\n"
               "selections:\n"
               "  - !!python/object/apply:os.system ['echo pwned']\n")
    path = str(tmp_path / 'evil_selections.yml')
    _write(path, payload)
    with pytest.raises(ValueError, match='plain YAML'):
        load_pdep_network_selections(path=path)

    result_payload = ("exploration_result_schema_version: 1\n"
                      "results:\n"
                      "  - !!python/object/apply:os.system ['echo pwned']\n")
    result_path = str(tmp_path / 'evil_results.yml')
    _write(result_path, result_payload)
    with pytest.raises(ValueError, match='plain YAML'):
        load_pdep_exploration_results(path=result_path)


def _build_full_budget_record():
    """A PDepBudgetRecord with every field set to a non-default value, and two network_outcomes
    entries -- one admitted (a named network) and one refused (an unnamed offer) -- so a loader that
    rebuilt the record field-by-field but forgot a field (on either the record or the nested outcome)
    would be caught by round-trip equality rather than silently dropping it."""
    admitted = PDepBudgetNetworkOutcome(
        network_id='network1_1', outcome=BUDGET_OUTCOME_ADMITTED, cost=5, rank=0,
        network_source_hash='sha256:full_budget_record', method='master_equation',
        remaining_transition_states=3,
    )
    refused = PDepBudgetNetworkOutcome(
        network_id=None, unnamed_offer_index=7, outcome=BUDGET_OUTCOME_REFUSED, cost=2, rank=1,
        reason_code=BUDGET_SKIP_EXCEEDS_BUDGET, reason='exceeds the remaining transition-state budget',
    )
    return PDepBudgetRecord(
        iteration=3, max_transition_states=10, max_networks=4, total_cost=5,
        network_outcomes=(admitted, refused),
    )


def test_save_pdep_budget_record_round_trips_and_is_json_serializable(tmp_path):
    """Test that save_pdep_budget_record() writes something json.dumps() can handle AND that
    load_pdep_budget_record() (the STRICT loader real callers use, not a permissive read_yaml_file()
    call) accepts it back as the identical record. A prior version of this test only checked
    json.dumps(read_yaml_file(path)) -- but json.dumps((1, 2)) succeeds, so that assertion could not
    have caught a YAML-only type (e.g. a tuple) surviving into the saved file; and reading back with
    read_yaml_file() rather than the strict loader proved nothing about loader compatibility. Assert
    both properties directly: no non-plain (list/dict/str/int/float/bool/None) values anywhere in the
    saved structure, and successful round-trip through the strict loader."""
    record = _build_full_budget_record()
    path = str(tmp_path / 'budget.yml')

    save_pdep_budget_record(path=path, record=record)

    content = read_yaml_file(path=path)
    json.dumps(content)  # must not raise

    def _assert_only_plain_types(value):
        if isinstance(value, dict):
            for key, item in value.items():
                assert isinstance(key, str), f'non-string mapping key {key!r} in saved record'
                _assert_only_plain_types(item)
        elif isinstance(value, list):
            for item in value:
                _assert_only_plain_types(item)
        else:
            assert value is None or isinstance(value, (str, int, float, bool)), (
                f'non-plain value {value!r} of type {type(value).__name__} in saved record; '
                f'a tuple, set, or other YAML-only type here would satisfy json.dumps((1, 2)) but '
                f'still be incompatible with what a strict JSON-based caller expects')

    _assert_only_plain_types(content)

    loaded = load_pdep_budget_record(path=path)
    assert loaded == record


def test_load_pdep_budget_record_round_trip_equality(tmp_path):
    """Test that load(save(x)) equals x for a fully-populated PDepBudgetRecord, derived from the live
    object rather than a hand-written expected dict -- so the test cannot silently pass by omitting a
    field neither side checks. This is the exact defect class this branch already shipped once: a
    field present in as_dict() but forgotten in the loader is silently dropped on round-trip, and a
    round-trip test using only default field values would not catch it."""
    record = _build_full_budget_record()
    path = str(tmp_path / 'budget.yml')
    save_pdep_budget_record(path=path, record=record)

    loaded = load_pdep_budget_record(path=path)

    assert loaded == record


def test_load_pdep_budget_record_round_trip_equality_empty_record(tmp_path):
    """Test that a budget record with zero network_outcomes (the "nothing was refused and nothing
    was admitted" case -- e.g. the budget ran over zero candidates) round-trips correctly. Absence of
    any outcomes must not be confused with absence of the file itself."""
    record = PDepBudgetRecord(iteration=0, max_transition_states=None, max_networks=None,
                              total_cost=0, network_outcomes=())
    path = str(tmp_path / 'budget.yml')
    save_pdep_budget_record(path=path, record=record)

    loaded = load_pdep_budget_record(path=path)

    assert loaded == record
    assert loaded.network_outcomes == ()


def test_load_pdep_budget_record_refuses_non_mapping_top_level(tmp_path):
    """Test the top-level-not-a-mapping refusal, on its own distinctive wording."""
    path = str(tmp_path / 'budget.yml')
    save_yaml_file(path=path, content=[1, 2, 3])

    with pytest.raises(ValueError) as exc_info:
        load_pdep_budget_record(path=path)
    assert 'does not contain a mapping' in str(exc_info.value)


def test_load_pdep_budget_record_refuses_missing_version_key(tmp_path):
    """Test the absent-version-key refusal, on its own distinctive wording -- distinct from the
    unknown-version wording checked below (an unversioned file is not "version 1", it is unknown)."""
    path = str(tmp_path / 'budget.yml')
    save_yaml_file(path=path, content={'iteration': 0, 'total_cost': 0, 'network_outcomes': []})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_budget_record(path=path)
    assert "no 'budget_record_schema_version' key" in str(exc_info.value)
    assert 'unknown shape' in str(exc_info.value)


def test_load_pdep_budget_record_refuses_unknown_schema_version(tmp_path):
    """Test the unknown-schema-version refusal, on its own distinctive wording."""
    record = _build_full_budget_record()
    rendered = record.as_dict()
    rendered['budget_record_schema_version'] = 999
    path = str(tmp_path / 'budget.yml')
    save_yaml_file(path=path, content=rendered)

    with pytest.raises(ValueError) as exc_info:
        load_pdep_budget_record(path=path)
    assert 'only understands version' in str(exc_info.value)


def test_load_pdep_budget_record_refuses_unknown_algorithm_version(tmp_path):
    """Test the unknown-algorithm-version refusal, on its own distinctive wording, distinct from the
    schema-version refusal above."""
    record = _build_full_budget_record()
    rendered = record.as_dict()
    rendered['budget_algorithm_version'] = 999
    path = str(tmp_path / 'budget.yml')
    save_yaml_file(path=path, content=rendered)

    with pytest.raises(ValueError) as exc_info:
        load_pdep_budget_record(path=path)
    assert 'budget_algorithm_version' in str(exc_info.value)


def test_load_pdep_budget_record_refuses_missing_network_outcomes_key(tmp_path):
    """Test the missing/non-list 'network_outcomes' key refusal."""
    path = str(tmp_path / 'budget.yml')
    save_yaml_file(path=path, content={'budget_record_schema_version': BUDGET_RECORD_SCHEMA_VERSION,
                                       'budget_algorithm_version': BUDGET_ALGORITHM_VERSION,
                                       'iteration': 0, 'max_transition_states': None,
                                       'max_networks': None, 'total_cost': 0})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_budget_record(path=path)
    assert "no 'network_outcomes' list" in str(exc_info.value)


def test_load_pdep_budget_record_refuses_non_mapping_outcome(tmp_path):
    """Test the per-outcome not-a-mapping refusal."""
    path = str(tmp_path / 'budget.yml')
    save_yaml_file(path=path, content={'budget_record_schema_version': BUDGET_RECORD_SCHEMA_VERSION,
                                       'budget_algorithm_version': BUDGET_ALGORITHM_VERSION,
                                       'iteration': 0, 'max_transition_states': None,
                                       'max_networks': None, 'total_cost': 0,
                                       'network_outcomes': [42]})

    with pytest.raises(ValueError) as exc_info:
        load_pdep_budget_record(path=path)
    assert 'not a mapping' in str(exc_info.value)


def test_load_pdep_budget_record_refuses_a_string_for_cost(tmp_path):
    """Test per-field type validation on a nested network outcome: a string where an int (cost) is
    required must be refused, not silently coerced."""
    record = _build_full_budget_record()
    rendered = record.as_dict()
    rendered['network_outcomes'][0]['cost'] = 'five'
    path = str(tmp_path / 'budget.yml')
    save_yaml_file(path=path, content=rendered)

    with pytest.raises(ValueError) as exc_info:
        load_pdep_budget_record(path=path)
    assert 'cost' in str(exc_info.value)


def test_load_pdep_budget_record_refuses_an_unrecognized_outcome(tmp_path):
    """Test that an 'outcome' value outside VALID_BUDGET_OUTCOMES is refused."""
    record = _build_full_budget_record()
    rendered = record.as_dict()
    rendered['network_outcomes'][0]['outcome'] = 'maybe'
    path = str(tmp_path / 'budget.yml')
    save_yaml_file(path=path, content=rendered)

    with pytest.raises(ValueError) as exc_info:
        load_pdep_budget_record(path=path)
    assert 'outcome' in str(exc_info.value)


def test_load_pdep_budget_record_refuses_a_python_object_tag(tmp_path):
    """Test that load_pdep_budget_record() cannot be made to construct a Python object from a YAML
    tag -- mirroring test_the_loaders_refuse_a_python_object_tag above for the other two loaders."""
    payload = ("budget_record_schema_version: 1\n"
               "budget_algorithm_version: 1\n"
               "iteration: 0\n"
               "max_transition_states: null\n"
               "max_networks: null\n"
               "total_cost: 0\n"
               "network_outcomes:\n"
               "  - !!python/object/apply:os.system ['echo pwned']\n")
    path = str(tmp_path / 'evil_budget.yml')
    _write(path, payload)
    with pytest.raises(ValueError, match='plain YAML'):
        load_pdep_budget_record(path=path)


def test_load_pdep_budget_record_refuses_a_bool_schema_version(tmp_path):
    """Test that ``budget_record_schema_version: true`` is refused rather than silently accepted as
    version 1. In Python ``bool`` is a subclass of ``int`` and ``True == 1``, so a bare ``!=`` version
    comparison (with no ``isinstance(..., bool)`` exclusion) would let this through as if it were the
    real integer version -- exactly the gap ``_require_int_field``/``_validate_budget`` already guard
    against for other fields in this module."""
    record = _build_full_budget_record()
    rendered = record.as_dict()
    rendered['budget_record_schema_version'] = True
    path = str(tmp_path / 'budget.yml')
    save_yaml_file(path=path, content=rendered)

    with pytest.raises(ValueError) as exc_info:
        load_pdep_budget_record(path=path)
    assert 'only understands version' in str(exc_info.value)


def test_load_pdep_budget_record_refuses_a_bool_algorithm_version(tmp_path):
    """Test that ``budget_algorithm_version: true`` is refused rather than silently accepted as
    version 1, mirroring the schema-version bool-bypass test above."""
    record = _build_full_budget_record()
    rendered = record.as_dict()
    rendered['budget_algorithm_version'] = True
    path = str(tmp_path / 'budget.yml')
    save_yaml_file(path=path, content=rendered)

    with pytest.raises(ValueError) as exc_info:
        load_pdep_budget_record(path=path)
    assert 'budget_algorithm_version' in str(exc_info.value)


def test_pdep_budget_record_constructor_refuses_a_bool_schema_version():
    """Test the same bool-bypass gap directly on PDepBudgetRecord.__post_init__ (t3/pdep/budget.py),
    independent of the YAML loader in t3/pdep/api.py -- both gates must reject a bool, since a caller
    can construct a PDepBudgetRecord directly without going through the loader at all."""
    admitted = PDepBudgetNetworkOutcome(
        network_id='network1_1', outcome=BUDGET_OUTCOME_ADMITTED, cost=5, rank=0,
        network_source_hash='sha256:x', method='master_equation',
    )
    with pytest.raises(ValueError, match='PDepBudgetRecord.schema_version'):
        PDepBudgetRecord(iteration=0, max_transition_states=None, max_networks=None, total_cost=5,
                         network_outcomes=(admitted,), schema_version=True)


def test_pdep_budget_record_constructor_refuses_a_bool_algorithm_version():
    """Test the same bool-bypass gap directly on PDepBudgetRecord.__post_init__ for
    algorithm_version, mirroring the schema-version test above."""
    admitted = PDepBudgetNetworkOutcome(
        network_id='network1_1', outcome=BUDGET_OUTCOME_ADMITTED, cost=5, rank=0,
        network_source_hash='sha256:x', method='master_equation',
    )
    with pytest.raises(ValueError, match='PDepBudgetRecord.algorithm_version'):
        PDepBudgetRecord(iteration=0, max_transition_states=None, max_networks=None, total_cost=5,
                         network_outcomes=(admitted,), algorithm_version=True)


def test_pdep_budget_network_outcome_refuses_a_negative_remaining_transition_states():
    """Test that PDepBudgetNetworkOutcome.__post_init__ rejects a negative
    remaining_transition_states, mirroring the non-negativity checks it already applies to cost and
    rank. apply_pdep_qm_budget() itself never produces a negative remaining_transition_states, but a
    hand-built or loaded-from-disk record could carry one, and a negative "transition states left"
    value is nonsensical and must not load as authoritative."""
    with pytest.raises(ValueError, match='remaining_transition_states'):
        PDepBudgetNetworkOutcome(
            network_id='network1_1', outcome=BUDGET_OUTCOME_ADMITTED, cost=5, rank=0,
            network_source_hash='sha256:x', method='master_equation',
            remaining_transition_states=-1,
        )


def test_load_pdep_budget_record_refuses_a_negative_remaining_transition_states(tmp_path):
    """Test the same negative-remaining_transition_states rejection from the loader side: a record
    file on disk with remaining_transition_states=-1 must be refused, not accepted as authoritative."""
    record = _build_full_budget_record()
    rendered = record.as_dict()
    rendered['network_outcomes'][0]['remaining_transition_states'] = -1
    path = str(tmp_path / 'budget.yml')
    save_yaml_file(path=path, content=rendered)

    with pytest.raises(ValueError, match='remaining_transition_states'):
        load_pdep_budget_record(path=path)


def test_save_pdep_budget_record_replaces_a_pre_existing_symlink_rather_than_writing_through_it(tmp_path):
    """Test the symlink write-through vulnerability directly: if `path` is a symlink pointing outside
    the destination directory, saving a budget record there must replace the symlink itself with a
    regular file, never write through the link and mutate whatever it points to. This is the
    concrete attack save_pdep_budget_record's atomic-write fix (mkstemp + os.replace) is meant to
    close -- os.replace() unlinks/replaces the symlink entry rather than following it, unlike a naive
    open(path, 'w')."""
    outside_dir = tmp_path / 'outside'
    outside_dir.mkdir()
    target = outside_dir / 'victim.yml'
    target.write_text('sentinel: do-not-touch\n')

    iteration_dir = tmp_path / 'iteration_1'
    iteration_dir.mkdir()
    link_path = iteration_dir / 'budget.yml'
    os.symlink(str(target), str(link_path))

    record = _build_full_budget_record()
    save_pdep_budget_record(path=str(link_path), record=record)

    # The outside target must be untouched -- the write must not have gone through the symlink.
    assert target.read_text() == 'sentinel: do-not-touch\n'
    # The path where the record was requested must no longer be a symlink: it was replaced by a
    # regular file carrying the saved record.
    assert not os.path.islink(str(link_path))
    loaded = load_pdep_budget_record(path=str(link_path))
    assert loaded == record


# --- The writers refuse to write what the loaders could never read -------------------------------
#
# `_read_persisted_yaml_file`'s docstring states the invariant these tests pin: "Nothing T3 writes
# here [needs tag support] ... any file containing one is not a file T3 wrote." That was an
# assertion about the writers, made on the read side, and nothing enforced it on the write side.
# `arc.common.save_yaml_file` dumps with a full representer, so ONE non-plain value anywhere in a
# record -- a `pathlib.Path` handed to a `str`-annotated field is the realistic way in -- is written
# out as a `!!python/object/apply:...` tag, the write reports success, and the file is then refused
# by every T3 loader forever. The corruption is silent, durable, and discovered only later, by which
# point the iteration the record described is gone.

def test_save_pdep_network_selections_refuses_a_record_the_loader_could_never_read(tmp_path):
    """Test that a non-plain field value is refused at write time rather than written as a Python
    object tag. A `Path` in `sa_path` is the realistic case: the field is annotated `str`, nothing
    enforces that, and `os.path.join`/`open` all accept a Path happily, so it reaches the dumper."""
    selection = PDepNetworkSelection(network_id='network4_2', sa_path=Path(tmp_path) / 'sa.yml')
    path = str(tmp_path / 'selections.yml')

    with pytest.raises(ValueError, match='plain YAML'):
        save_pdep_network_selections(path=path, selections=[selection])
    assert not os.path.exists(path), 'the refused write must not leave a file behind'


def test_a_refused_write_does_not_destroy_the_record_already_on_disk(tmp_path):
    """Test that refusing the write leaves the PREVIOUS record intact. This is the half that makes
    the refusal worth having: `save_pdep_network_selections` writes straight to `path`, so without
    the check the corrupt dump replaces a good file's contents, and the run loses both the new
    record and the one it already had."""
    path = str(tmp_path / 'selections.yml')
    good = PDepNetworkSelection(network_id='network4_2', sa_path=str(tmp_path / 'sa.yml'))
    save_pdep_network_selections(path=path, selections=[good])

    # Constructed OUTSIDE the `raises` block on purpose: when `PDepNetworkSelection` grows the
    # per-field validation its sibling record types already have, the refusal moves to this line and
    # the test fails loudly rather than passing while no longer exercising clobber resistance.
    bad = PDepNetworkSelection(network_id='network5_1', sa_path=Path(tmp_path) / 'sa.yml')
    with pytest.raises(ValueError, match='plain YAML'):
        save_pdep_network_selections(path=path, selections=[bad])

    assert load_pdep_network_selections(path=path) == [good]


def test_save_pdep_exploration_results_refuses_a_record_the_loader_could_never_read(tmp_path):
    """Test the same contract on the exploration-result writer. Each result nests a serialized
    selection, so the offending value here is one level down -- the check has to look at the whole
    rendered content, not just its top-level keys."""
    result = PDepExplorationResult(
        network_id='network4_2',
        status=EXPLORATION_STATUS_SKIPPED,
        reasons=('not qualified',),
        selection=PDepNetworkSelection(network_id='network4_2', sa_path=Path(tmp_path) / 'sa.yml'),
    )
    path = str(tmp_path / 'results.yml')

    with pytest.raises(ValueError, match='plain YAML'):
        save_pdep_exploration_results(path=path, results=[result])
    assert not os.path.exists(path)


def test_the_check_draws_its_line_at_parseability_not_at_schema_validity(tmp_path):
    """Test the BOUNDARY of the write-time check, so the next reader does not mistake it for
    validation. A non-numeric threshold is perfectly good YAML: it is written, and refused later per
    record with a message naming the field. A `Path` in the same dict is not YAML at all: it costs
    the whole file, and is refused now. The two are different failures, and only the second one is
    this check's business -- giving `PDepNetworkSelection` the per-field validation its sibling
    record types already have is the separate, named piece of work that closes the first."""
    path = str(tmp_path / 'selections.yml')
    schema_invalid = PDepNetworkSelection(network_id='network4_2',
                                          thresholds={'relative_threshold': 'not a number'})
    save_pdep_network_selections(path=path, selections=[schema_invalid])
    with pytest.raises(ValueError, match='relative_threshold'):
        load_pdep_network_selections(path=path)

    unparseable = PDepNetworkSelection(network_id='network4_2',
                                       thresholds={'relative_threshold': 0.001, 'source': Path(tmp_path)})
    with pytest.raises(ValueError, match='plain YAML'):
        save_pdep_network_selections(path=str(tmp_path / 'other.yml'), selections=[unparseable])


def test_the_check_runs_on_the_bytes_that_will_actually_be_written(tmp_path, monkeypatch):
    """Test the budget writer's call site, which no record it can legitimately hold can reach:
    `PDepBudgetNetworkOutcome.__post_init__` type-checks every `str` field, so the guard there is a
    contract for future fields rather than a live gate. Reaching it means forcing `as_dict()` to
    return something it never would -- which is also the only way to pin that the refusal happens
    BEFORE `mkstemp`, and so leaves no `.pdep-budget-*` dropping behind."""
    monkeypatch.setattr(PDepBudgetRecord, 'as_dict', lambda self: {'iteration': Path(tmp_path)})
    path = str(tmp_path / 'budget.yml')

    with pytest.raises(ValueError, match='plain YAML'):
        save_pdep_budget_record(path=path, record=_build_full_budget_record())
    assert os.listdir(str(tmp_path)) == [], 'no file and no staged temp file may survive'


def test_the_check_uses_the_same_renderer_the_write_uses(tmp_path):
    """Test that the check is not merely `yaml.safe_dump`. `arc.common.save_yaml_file` renders via
    `to_yaml`, which installs its own `str` representer on the global dumper -- so a check performed
    with a different dumper would be approving bytes other than the ones that reach the file, which
    is the exact shape of the bug being prevented. Pinning it here because the two happen to agree
    today, and a check that silently stops matching the writer is worse than no check."""
    content = {'a': 'plain', 'b': ['multi\nline\nstring', 0.5, None, True]}
    t3_pdep_api._refuse_content_that_would_not_parse_back(content=content, path=str(tmp_path / 'x.yml'))

    path = str(tmp_path / 'x.yml')
    save_yaml_file(path=path, content=content)
    with open(path, 'r') as f:
        assert yaml.safe_load(f) == content
