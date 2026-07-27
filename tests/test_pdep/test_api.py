#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_api
"""

import datetime
import json
import math
import os

import pytest
from arc.common import read_yaml_file, save_yaml_file

import t3.main as t3_main
from t3.common import TEST_DATA_BASE_PATH
from t3.logger import Logger
from t3.pdep.api import rank_pdep_networks, save_pdep_network_selections, select_pdep_network
from t3.pdep.cache import write_sa_cache_metadata
from t3.pdep.parser import parse_pdep_network_file
from t3.pdep.selector import (CACHE_STATUS_CACHED_REJECTED,
                              CACHE_STATUS_CACHED_VALID,
                              CACHE_STATUS_UNVALIDATED,
                              EVALUATION_STATUS_NOT_EVALUATED,
                              SELECTOR_VERSION,
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


# --- 4. NotImplementedError for explore=True -----------------------------------------------------

def test_select_pdep_network_raises_not_implemented_for_explore(sa_dict):
    """Test that explore=True raises NotImplementedError (Commit 4 is not implemented yet)."""
    with pytest.raises(NotImplementedError):
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


# --- 7. save_pdep_network_selections round trip and JSON-serializability ------------------------

def test_save_pdep_network_selections_round_trips_and_is_json_serializable(tmp_path, sa_dict):
    """Test that saved selections round-trip through read_yaml_file and are JSON-serializable, and
    carry a selector_version schema marker (FIX M)."""
    network = parse_pdep_network_file(path=NETWORK_PATH)
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001)
    path = str(tmp_path / 'selections.yml')
    returned_path = save_pdep_network_selections(path=path, selections=[selection])

    assert returned_path == path
    loaded = read_yaml_file(path=path)
    assert loaded == {'selector_version': SELECTOR_VERSION, 'selections': [selection.as_dict()]}
    serialized = json.dumps(loaded)
    assert isinstance(serialized, str)


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
