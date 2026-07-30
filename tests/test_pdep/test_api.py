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
import t3.pdep.api as t3_pdep_api
from t3.pdep.api import (explore_pdep_network,
                         rank_pdep_networks,
                         save_pdep_network_selections,
                         select_pdep_network,
                         )
from t3.pdep.cache import hash_file, write_sa_cache_metadata
from t3.pdep.explorer.config import PDepExplorerConfig
from t3.pdep.explorer.result import (EXPLORATION_STATUS_FAILED,
                                     EXPLORATION_STATUS_SKIPPED,
                                     EXPLORATION_STATUS_SUCCEEDED,
                                     )
from t3.pdep.parser import parse_pdep_network_file
from t3.pdep.selector import (CACHE_STATUS_CACHED_REJECTED,
                              CACHE_STATUS_CACHED_VALID,
                              CACHE_STATUS_UNVALIDATED,
                              EVALUATION_STATUS_NOT_EVALUATED,
                              SELECTOR_VERSION,
                              PDepNetworkSelection,
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
    assert result.selection is selection
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
