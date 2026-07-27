#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep _wiring_helpers module

Shared fixture-construction helpers for tests that need a real, working T3 instance wired up
against the ``network4_2`` PDep fixture tree (``tests/data/pdep_network/iteration_1``).

This is used by both ``test_main_wiring.py`` (the in-run wiring acceptance tests) and
``test_api.py`` (FIX A's parity test, which compares ``t3.pdep.api``'s decision against a real
``T3.determine_species_from_pdep_network()`` run) so the two test files exercise literally the
same T3 construction rather than two independently-maintained copies that could quietly drift
apart and stop actually testing the same thing.

Not itself a test module (no ``test_*`` functions), so pytest does not collect it directly.
"""

import os

from t3.chem import T3Reaction
from t3.common import TEST_DATA_BASE_PATH
from t3.utils.rmg_shim import PDepNetwork

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


def build_t3(tmp_path):
    """Build a real T3 instance against a tmp_path copy of the iteration_1 fixture tree.

    Also stages the network4_2/MSC sensitivity sidecar fixture (which lives outside
    pdep_network/, see ``SA_FIXTURE_PATH`` above) into the tmp_path copy, at the exact
    location it would occupy in a real run.
    """
    import shutil

    project_directory = str(tmp_path)
    shutil.copytree(FIXTURE_ITERATION_1, os.path.join(project_directory, 'iteration_1'))
    msc_sensitivity_dir = os.path.join(project_directory, 'iteration_1', 'PDep_SA', 'network4_2',
                                       'MSC', 'sensitivity')
    os.makedirs(msc_sensitivity_dir, exist_ok=True)
    shutil.copy(SA_FIXTURE_PATH, os.path.join(msc_sensitivity_dir, 'sa_coefficients.yml'))
    t3 = run_minimal(project_directory=project_directory, iteration=1, set_paths=True)
    t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_yaml_file()
    return t3


def build_pdep_rxns_to_explore(t3):
    """A T3Reaction (the production path) referencing real species, pressure-dependent, network 4."""
    reaction = T3Reaction(
        r_species=[t3.rmg_species[H_INDEX], t3.rmg_species[C4ENE_INDEX]],
        p_species=[t3.rmg_species[C4RAD_INDEX]],
        is_pressure_dependent=True,
        network=PDepNetwork(index=4),
    )
    return [(reaction, 2, t3.rmg_species[C4RAD_INDEX].label)]


def build_sa_dict(t3):
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


def network_path(t3):
    return os.path.join(t3.paths['RMG PDep'], f'{NETWORK_NAME}.py')


def candidate_sa_path(t3, method='CSE'):
    return os.path.join(t3.paths['PDep SA'], NETWORK_NAME, method, 'sensitivity', 'sa_coefficients.yml')
