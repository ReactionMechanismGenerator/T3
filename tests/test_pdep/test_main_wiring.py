#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_main_wiring module

Acceptance test for T3.determine_species_from_pdep_network()'s in-run wiring: both cache paths
(reuse a valid cached SA, and generate + persist a new one), the cache-rejected path (an SA
output present without T3's own sidecar must not be trusted), and a cross-check that the wiring's
`qualified`/`network_id` verdict agrees with a direct, independent call into
`t3.pdep.selector.select_from_sa_dict()` for the same inputs.

Real Arkane is never invoked here: `t3.main.run_arkane_job` is monkeypatched throughout, and a
synthetic Arkane sensitivity dictionary stands in for real output. The real fixture network file
(`tests/data/pdep_network/iteration_1/RMG/pdep/network4_2.py`) is used, but the whole
`iteration_1` fixture tree is copied into `tmp_path` first, so nothing under `tests/data/` is
ever written to.
"""

import os
import shutil

import pytest

import t3.main as t3_main
from t3.chem import T3Reaction
from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.cache import sa_cache_metadata_path, write_sa_cache_metadata
from t3.pdep.parser import parse_pdep_network_file
from t3.pdep.selector import (CACHE_STATUS_CACHED_VALID,
                              CACHE_STATUS_GENERATED,
                              select_from_sa_dict,
                              )
from t3.utils.rmg_shim import PDepNetwork
from arc.common import save_yaml_file

from tests.common import run_minimal

FIXTURE_ITERATION_1 = os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'iteration_1')

# The real Arkane sensitivity fixture for network4_2/MSC lives OUTSIDE pdep_network/, at
# tests/data/pdep_sa/network4_2_MSC/sa_coefficients.yml, because test_main.py's
# test_determine_species_from_pdep_network sets t3.paths['PDep SA'] to
# tests/data/pdep_network/iteration_1/PDep_SA and its teardown does shutil.rmtree() on that path,
# which would delete it if it lived there. It is staged into the tmp_path copy below.
SA_FIXTURE_PATH = os.path.join(TEST_DATA_BASE_PATH, 'pdep_sa', 'network4_2_MSC', 'sa_coefficients.yml')

# The network reaction under test, in network species labels: reaction2 / TS2 in network4_2.py,
# i.e. H(34) + C4ene(26) <=> C4rad(5). Real RMG species indices established by the existing
# `tests/test_main.py::test_determine_species_from_pdep_network` test.
H_INDEX, C4ENE_INDEX, C4RAD_INDEX = 35, 27, 6
NETWORK_REACTION_STR = 'H(34) + C4ene(26) <=> C4rad(5)'
NETWORK_NAME = 'network4_2'
CONDITION = (1000.0, 'K', 1.0, 'bar')

# Sensitive to two transition states of the same network at this condition:
#  - TS2 (reaction2, the network reaction itself): kinetics is a training-reaction exact match,
#    i.e. certain.
#  - TS1 (reaction1, a different path reaction of the same network): kinetics is an estimate,
#    i.e. uncertain. This is what should make the network qualify for QM refinement.
TS2_COEFFICIENT = 0.05
TS1_COEFFICIENT = 0.04


def _build_t3(tmp_path):
    """Build a real T3 instance against a tmp_path copy of the iteration_1 fixture tree.

    Also stages the network4_2/MSC sensitivity sidecar fixture (which lives outside
    pdep_network/, see ``SA_FIXTURE_PATH`` above) into the tmp_path copy, at the exact
    location it would occupy in a real run.
    """
    project_directory = str(tmp_path)
    shutil.copytree(FIXTURE_ITERATION_1, os.path.join(project_directory, 'iteration_1'))
    msc_sensitivity_dir = os.path.join(project_directory, 'iteration_1', 'PDep_SA', 'network4_2',
                                       'MSC', 'sensitivity')
    os.makedirs(msc_sensitivity_dir, exist_ok=True)
    shutil.copy(SA_FIXTURE_PATH, os.path.join(msc_sensitivity_dir, 'sa_coefficients.yml'))
    t3 = run_minimal(project_directory=project_directory, iteration=1, set_paths=True)
    t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_yaml_file()
    return t3


def _build_pdep_rxns_to_explore(t3):
    """A T3Reaction (the production path) referencing real species, pressure-dependent, network 4."""
    reaction = T3Reaction(
        r_species=[t3.rmg_species[H_INDEX], t3.rmg_species[C4ENE_INDEX]],
        p_species=[t3.rmg_species[C4RAD_INDEX]],
        is_pressure_dependent=True,
        network=PDepNetwork(index=4),
    )
    return [(reaction, 2, t3.rmg_species[C4RAD_INDEX].label)]


def _build_sa_dict(t3):
    """A synthetic Arkane PDep sensitivity dictionary standing in for real Arkane SA output."""
    return {
        'structures': {
            'H(34)': t3.rmg_species[H_INDEX].mol.to_adjacency_list(),
            'C4ene(26)': t3.rmg_species[C4ENE_INDEX].mol.to_adjacency_list(),
            'C4rad(5)': t3.rmg_species[C4RAD_INDEX].mol.to_adjacency_list(),
        },
        NETWORK_REACTION_STR: {
            CONDITION: {
                '(TS) TS2': TS2_COEFFICIENT,
                '(TS) TS1': TS1_COEFFICIENT,
            },
        },
    }


def _network_path(t3):
    return os.path.join(t3.paths['RMG PDep'], f'{NETWORK_NAME}.py')


def _candidate_sa_path(t3, method='CSE'):
    return os.path.join(t3.paths['PDep SA'], NETWORK_NAME, method, 'sensitivity', 'sa_coefficients.yml')


class TestDetermineSpeciesFromPdepNetworkWiring(object):

    def test_cached_valid_path_never_invokes_arkane(self, tmp_path, monkeypatch):
        """A valid, T3-vouched-for cache must be reused; Arkane must never be invoked."""
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_path = _candidate_sa_path(t3, method='CSE')
        os.makedirs(os.path.dirname(sa_path), exist_ok=True)
        save_yaml_file(sa_path, _build_sa_dict(t3))
        write_sa_cache_metadata(sa_path=sa_path,
                                network_path=_network_path(t3),
                                network_id=NETWORK_NAME,
                                method='CSE',
                                )

        def _fail_if_called(*args, **kwargs):
            pytest.fail('run_arkane_job should not be invoked when a valid cache is present.')

        monkeypatch.setattr(t3_main, 'run_arkane_job', _fail_if_called)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert len(t3.pdep_network_selections) == 1
        selection = t3.pdep_network_selections[0]
        assert selection.network_id == NETWORK_NAME
        assert selection.cache_status == CACHE_STATUS_CACHED_VALID
        assert selection.qualified is True

    def test_generate_path_invokes_arkane_and_persists_sidecar(self, tmp_path, monkeypatch):
        """No cache present: Arkane must be invoked, and a sidecar must be written on success."""
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_dict = _build_sa_dict(t3)
        calls = list()

        def _fake_run_arkane_job(input_file, output_directory, plot=True, logger=None, **kwargs):
            calls.append(output_directory)
            sensitivity_dir = os.path.join(output_directory, 'sensitivity')
            os.makedirs(sensitivity_dir, exist_ok=True)
            save_yaml_file(os.path.join(sensitivity_dir, 'sa_coefficients.yml'), sa_dict)
            return True

        monkeypatch.setattr(t3_main, 'run_arkane_job', _fake_run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert len(calls) == 1, 'run_arkane_job should have been invoked exactly once (CSE, the first ME method).'
        assert len(t3.pdep_network_selections) == 1
        selection = t3.pdep_network_selections[0]
        assert selection.network_id == NETWORK_NAME
        assert selection.cache_status == CACHE_STATUS_GENERATED
        assert selection.qualified is True

        sa_path = _candidate_sa_path(t3, method='CSE')
        assert os.path.isfile(sa_path)
        assert os.path.isfile(sa_cache_metadata_path(sa_path)), \
            'A T3 sidecar should have been written alongside the freshly generated SA output.'

    def test_sa_output_without_sidecar_is_rejected_and_arkane_is_rerun(self, tmp_path, monkeypatch):
        """An SA output present without T3's own sidecar must not be trusted; Arkane re-runs."""
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_dict = _build_sa_dict(t3)

        # Pre-existing SA output from some prior, untracked run: no t3_sa_cache.yml sidecar.
        stale_sa_path = _candidate_sa_path(t3, method='CSE')
        os.makedirs(os.path.dirname(stale_sa_path), exist_ok=True)
        save_yaml_file(stale_sa_path, {'structures': {}, 'stale reaction': {}})
        assert not os.path.isfile(sa_cache_metadata_path(stale_sa_path))

        calls = list()

        def _fake_run_arkane_job(input_file, output_directory, plot=True, logger=None, **kwargs):
            calls.append(output_directory)
            sensitivity_dir = os.path.join(output_directory, 'sensitivity')
            os.makedirs(sensitivity_dir, exist_ok=True)
            save_yaml_file(os.path.join(sensitivity_dir, 'sa_coefficients.yml'), sa_dict)
            return True

        monkeypatch.setattr(t3_main, 'run_arkane_job', _fake_run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert len(calls) == 1, 'Arkane must be re-run when an SA output lacks a T3 sidecar.'
        selection = t3.pdep_network_selections[0]
        assert selection.cache_status == CACHE_STATUS_GENERATED
        assert selection.qualified is True

    def test_qualified_and_network_id_agree_with_direct_selector_call(self, tmp_path, monkeypatch):
        """The wiring's decision must agree with select_from_sa_dict() called directly."""
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_dict = _build_sa_dict(t3)
        sa_path = _candidate_sa_path(t3, method='CSE')
        os.makedirs(os.path.dirname(sa_path), exist_ok=True)
        save_yaml_file(sa_path, sa_dict)
        write_sa_cache_metadata(sa_path=sa_path,
                                network_path=_network_path(t3),
                                network_id=NETWORK_NAME,
                                method='CSE',
                                )

        def _fail_if_called(*args, **kwargs):
            pytest.fail('run_arkane_job should not be invoked when a valid cache is present.')

        monkeypatch.setattr(t3_main, 'run_arkane_job', _fail_if_called)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)
        selection = t3.pdep_network_selections[0]

        network = parse_pdep_network_file(_network_path(t3))
        direct_selection = select_from_sa_dict(
            sa_dict=sa_dict,
            network=network,
            network_reaction=NETWORK_REACTION_STR,
            relative_threshold=t3.t3['sensitivity']['pdep_SA_threshold'],
            min_delta_ln_k=t3.t3['sensitivity']['pdep_min_delta_ln_k'],
        )

        assert selection.network_id == direct_selection.network_id == NETWORK_NAME
        assert selection.qualified == direct_selection.qualified is True
