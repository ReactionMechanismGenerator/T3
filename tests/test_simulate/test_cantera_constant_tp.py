#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_cantera_constantTP module
"""

import os
import shutil

import numpy as np
from arc.common import read_yaml_file

from t3.chem import T3Reaction
from t3.common import (SIMULATE_TEST_DATA_BASE_PATH, TEST_DATA_BASE_PATH,
                        sa_dict_to_yaml, sa_dict_from_yaml)
from tests.common import run_minimal
from t3.simulate.cantera_constant_tp import CanteraConstantTP
from t3.simulate.factory import simulate_factory


TEST_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'cantera_simulator_test')
MINIMAL_DATA_DIR = os.path.join(TEST_DATA_BASE_PATH, 'minimal_data')
NH3_MODEL_PATH = os.path.join(TEST_DATA_BASE_PATH, 'models', 'NH3.yaml')


def test_set_up_no_sa():
    """
    Simulate RMG's minimal example without SA. By setting observable_list = list(), no SA is run.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    cantera_simulator_adapter = CanteraConstantTP(t3=t3.t3,
                                                  rmg=t3.rmg,
                                                  paths=t3.paths,
                                                  logger=t3.logger,
                                                  atol=t3.rmg['model']['atol'],
                                                  rtol=t3.rmg['model']['rtol'],
                                                  observable_list=list(),
                                                  )
    cantera_simulator_adapter.simulate()
    # check that there are over 20 time steps (as an arbitrary number) to indicate that the solver completed
    time = cantera_simulator_adapter.all_data[0][0].data
    assert len(time) > 20


def test_get_sa_coefficients():
    """
    Simulate RMG's minimal example with SA.
    Then run the `get_sa_coefficients()` method to test that Cantera correctly performed SA.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    observable_list = ['OH', 'H']
    cantera_simulator_adapter = CanteraConstantTP(t3=t3.t3,
                                                  rmg=t3.rmg,
                                                  paths=t3.paths,
                                                  logger=t3.logger,
                                                  atol=t3.rmg['model']['atol'],
                                                  rtol=t3.rmg['model']['rtol'],
                                                  observable_list=observable_list,
                                                  sa_atol=t3.t3['sensitivity']['atol'],
                                                  sa_rtol=t3.t3['sensitivity']['rtol'],
                                                  )
    cantera_simulator_adapter.simulate()
    sa_dict = cantera_simulator_adapter.get_sa_coefficients()
    # check that there are over 20 time steps (as an arbitrary number) to indicate that the solver completed
    assert len(sa_dict['time'][0]) > 20
    # check that there are SA data for the 2 requested species
    assert len(sa_dict['kinetics'][0]) == 2
    assert len(sa_dict['thermo'][0]) == 2


def test_get_idt_by_T():
    """
    Calculate the ignition delay time for RMG's minimal example.
    Since this adapter simulates at constant T, this method returns a dictionary whose values are empty lists.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    cantera_simulator_adapter = CanteraConstantTP(t3=t3.t3,
                                                  rmg=t3.rmg,
                                                  paths=t3.paths,
                                                  logger=t3.logger,
                                                  atol=t3.rmg['model']['atol'],
                                                  rtol=t3.rmg['model']['rtol'],
                                                  observable_list=list(),
                                                  )
    cantera_simulator_adapter.simulate()
    idt_dict = cantera_simulator_adapter.get_idt_by_T()
    assert len(idt_dict['idt']) == 0
    assert len(idt_dict['idt_index']) == 0


def test_find_equilibrium():
    """
    Find equilibrium mole fractions of RMG's minimal example.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    cantera_simulator_adapter = CanteraConstantTP(t3=t3.t3,
                                                  rmg=t3.rmg,
                                                  paths=t3.paths,
                                                  logger=t3.logger,
                                                  atol=t3.rmg['model']['atol'],
                                                  rtol=t3.rmg['model']['rtol'],
                                                  observable_list=list(),
                                                  )
    cantera_simulator_adapter.conditions
    equilibrium_dictionaries = cantera_simulator_adapter.find_equilibrium('TP')
    assert len(equilibrium_dictionaries[0].keys()) == 8


def test_get_t50():
    """
    Find half life of give species from RMG's minimal example.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    cantera_simulator_adapter = CanteraConstantTP(t3=t3.t3,
                                                  rmg=t3.rmg,
                                                  paths=t3.paths,
                                                  logger=t3.logger,
                                                  atol=t3.rmg['model']['atol'],
                                                  rtol=t3.rmg['model']['rtol'],
                                                  observable_list=list(),
                                                  )
    cantera_simulator_adapter.simulate()
    t50_list = cantera_simulator_adapter.get_t50('H2(1)')
    assert len(t50_list) == 1


def test_determine_reactions_based_on_sa():
    """
    Run SA via CanteraConstantTP on the minimal_data model, then call
    determine_reactions_based_on_sa to verify that T3 correctly identifies
    reactions whose rate coefficients need refinement.
    """
    t3 = run_minimal(project_directory=MINIMAL_DATA_DIR, iteration=1, set_paths=True)
    try:
        t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_yaml_file()
        sa_observables = ['H2', 'OH']
        simulate_adapter = simulate_factory(simulate_method='CanteraConstantTP',
                                            t3=t3.t3,
                                            rmg=t3.rmg,
                                            paths=t3.paths,
                                            logger=t3.logger,
                                            atol=t3.rmg['model']['atol'],
                                            rtol=t3.rmg['model']['rtol'],
                                            observable_list=sa_observables,
                                            sa_atol=t3.t3['sensitivity']['atol'],
                                            sa_rtol=t3.t3['sensitivity']['rtol'],
                                            )
        simulate_adapter.simulate()
        t3.sa_dict = simulate_adapter.get_sa_coefficients()
        assert t3.sa_dict is not None
        assert len(t3.sa_dict['kinetics'][0]) == 2
        reaction_keys = t3.determine_reactions_based_on_sa()
        assert len(reaction_keys) > 0, "Expected at least one reaction to require refinement"
        for key in reaction_keys:
            assert key in t3.reactions
            assert isinstance(t3.reactions[key], T3Reaction)
            assert len(t3.reactions[key].reasons) > 0
            assert 'most sensitive reaction' in t3.reactions[key].reasons[0]
    finally:
        _cleanup_minimal_data(t3)


def test_determine_species_based_on_sa():
    """
    Run SA via CanteraConstantTP on the minimal_data model, then call
    determine_species_based_on_sa to verify that T3 correctly identifies
    species whose thermo needs refinement.
    """
    t3 = run_minimal(project_directory=MINIMAL_DATA_DIR, iteration=1, set_paths=True)
    try:
        t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_yaml_file()
        sa_observables = ['H2', 'OH']
        simulate_adapter = simulate_factory(simulate_method='CanteraConstantTP',
                                            t3=t3.t3,
                                            rmg=t3.rmg,
                                            paths=t3.paths,
                                            logger=t3.logger,
                                            atol=t3.rmg['model']['atol'],
                                            rtol=t3.rmg['model']['rtol'],
                                            observable_list=sa_observables,
                                            sa_atol=t3.t3['sensitivity']['atol'],
                                            sa_rtol=t3.t3['sensitivity']['rtol'],
                                            )
        simulate_adapter.simulate()
        t3.sa_dict = simulate_adapter.get_sa_coefficients()
        assert t3.sa_dict is not None
        assert len(t3.sa_dict['thermo'][0]) == 2
        species_keys = t3.determine_species_based_on_sa()
        assert isinstance(species_keys, list)
        assert len(species_keys) > 0, "Expected at least one species to require refinement"
        for key in species_keys:
            assert key in t3.species
    finally:
        _cleanup_minimal_data(t3)


def _cleanup_minimal_data(t3):
    """Remove SA / PDep_SA dirs and t3.log left under MINIMAL_DATA_DIR by a test."""
    for path in [t3.paths.get('SA'), t3.paths.get('PDep SA'), os.path.join(MINIMAL_DATA_DIR, 't3.log')]:
        if not path:
            continue
        if os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
        elif os.path.isfile(path):
            os.remove(path)


NH3_SA_TEST_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'nh3_sa_test')


def _run_nh3_sa():
    """
    Helper: run Cantera SA on the NH3 model with H and NH2 as observables.
    Returns (t3_object, sa_dict).
    """
    t3 = run_minimal(project_directory=NH3_SA_TEST_DIR, iteration=1, set_paths=True)
    # Override the cantera path to point at the NH3 model
    t3.paths['cantera annotated'] = NH3_MODEL_PATH
    # Build a minimal rmg species list matching the NH3 model.
    # The adapter matches by label (without RMG index) against Cantera species names.
    # Observable species (H, NH2) must be in the species list for the adapter to track them.
    t3.rmg['species'] = [
        {'label': 'NH3', 'smiles': 'N', 'concentration': 0.5,
         'observable': False, 'SA_observable': False, 'balance': False, 'reactive': True},
        {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0.25,
         'observable': False, 'SA_observable': False, 'balance': False, 'reactive': True},
        {'label': 'N2', 'smiles': 'N#N', 'concentration': 0.25,
         'observable': False, 'SA_observable': False, 'balance': False, 'reactive': True},
        {'label': 'H', 'smiles': '[H]', 'concentration': 0,
         'observable': False, 'SA_observable': True, 'balance': False, 'reactive': True},
        {'label': 'NH2', 'smiles': '[NH2]', 'concentration': 0,
         'observable': False, 'SA_observable': True, 'balance': False, 'reactive': True},
    ]
    observable_list = ['H', 'NH2']
    adapter = CanteraConstantTP(t3=t3.t3,
                                rmg=t3.rmg,
                                paths=t3.paths,
                                logger=t3.logger,
                                atol=1e-16,
                                rtol=1e-8,
                                observable_list=observable_list,
                                sa_atol=1e-6,
                                sa_rtol=1e-4,
                                )
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    return t3, sa_dict


def test_sa_coefficients_nh3():
    """
    Run Cantera SA on the NH3 model with H and NH2 observables and verify
    the sa_dict structure.
    """
    print(f'Running {NH3_SA_TEST_DIR}')
    t3, sa_dict = _run_nh3_sa()
    assert sa_dict is not None
    # Two observables
    assert 'H(24)' in sa_dict['kinetics'][0]
    assert 'NH2(6)' in sa_dict['kinetics'][0]
    # Time series should have entries
    assert len(sa_dict['time'][0]) > 10
    # Kinetics SA should have entries for reactions (int keys)
    for obs in ['H(24)', 'NH2(6)']:
        params = sa_dict['kinetics'][0][obs]
        assert len(params) > 0
        for param, values in params.items():
            assert isinstance(param, int)
            assert isinstance(values, np.ndarray)
            assert len(values) == len(sa_dict['time'][0])
    # shutil.rmtree(NH3_SA_TEST_DIR, ignore_errors=True)


def test_sa_yaml_round_trip_nh3():
    """
    Run Cantera SA on the NH3 model, serialize to YAML, load back,
    and verify the round-trip preserves all data.
    """
    from t3.common import save_yaml_file

    t3, sa_dict = _run_nh3_sa()
    # Serialize
    metadata = {'adapter': 'CanteraConstantTP', 'iteration': 1,
                'observables': list(sa_dict['kinetics'][0].keys())}
    serialized = sa_dict_to_yaml(sa_dict, metadata=metadata)

    # Verify serialized structure uses conditions list
    assert 'conditions' in serialized
    assert len(serialized['conditions']) == 1
    cond0 = serialized['conditions'][0]
    assert isinstance(cond0['time'], list)
    assert isinstance(cond0['time'][0], float)
    for obs, params in cond0['kinetics'].items():
        for param, values in params.items():
            assert isinstance(values, list)
    assert serialized['metadata'] == metadata

    # Write to disk and read back
    sa_path = os.path.join(NH3_SA_TEST_DIR, 'sa_coefficients.yml')
    os.makedirs(os.path.dirname(sa_path), exist_ok=True)
    save_yaml_file(path=sa_path, content=serialized, top_keys=['metadata'])
    assert os.path.isfile(sa_path)

    # Verify metadata appears first in the file
    with open(sa_path, 'r') as f:
        first_line = f.readline()
    assert first_line.startswith('metadata:'), f"Expected metadata at top, got: {first_line!r}"

    loaded_raw = read_yaml_file(path=sa_path)
    loaded_sa_dict = sa_dict_from_yaml(loaded_raw)

    # Verify round-trip: time (per-condition lists)
    assert len(loaded_sa_dict['time']) == 1
    np.testing.assert_array_almost_equal(loaded_sa_dict['time'][0], sa_dict['time'][0])

    # Verify round-trip: kinetics
    assert set(loaded_sa_dict['kinetics'][0].keys()) == set(sa_dict['kinetics'][0].keys())
    for obs in sa_dict['kinetics'][0]:
        assert set(loaded_sa_dict['kinetics'][0][obs].keys()) == set(sa_dict['kinetics'][0][obs].keys())
        for param in sa_dict['kinetics'][0][obs]:
            np.testing.assert_array_almost_equal(
                loaded_sa_dict['kinetics'][0][obs][param],
                sa_dict['kinetics'][0][obs][param],
            )

    # Thermo section should exist as a per-condition list
    assert 'thermo' in loaded_sa_dict
    assert isinstance(loaded_sa_dict['thermo'], list)
    assert len(loaded_sa_dict['thermo']) == 1
    # shutil.rmtree(NH3_SA_TEST_DIR, ignore_errors=True)


def test_dump_sa_coefficients_nh3():
    """
    Run Cantera SA on the NH3 model and test the T3.dump_sa_coefficients /
    load_sa_coefficients methods end-to-end.
    """
    t3, sa_dict = _run_nh3_sa()
    t3.sa_dict = sa_dict
    t3.dump_sa_coefficients()

    sa_path = t3.paths['SA coefficients']
    assert os.path.isfile(sa_path)

    # Read the file and check structure
    raw = read_yaml_file(path=sa_path)
    assert 'metadata' in raw
    assert raw['metadata']['adapter'] == 'CanteraConstantTP'
    assert raw['metadata']['iteration'] == 1
    assert set(raw['metadata']['observables']) == {'H(24)', 'NH2(6)'}
    assert 'conditions' in raw
    assert len(raw['conditions']) == 1

    # Test load round-trip
    t3.sa_dict = None
    t3.load_sa_coefficients()
    assert t3.sa_dict is not None
    assert len(t3.sa_dict['time'][0]) > 10
    assert 'H(24)' in t3.sa_dict['kinetics'][0]
    assert 'NH2(6)' in t3.sa_dict['kinetics'][0]
    for obs in ['H(24)', 'NH2(6)']:
        for param, values in t3.sa_dict['kinetics'][0][obs].items():
            assert isinstance(param, int)
            assert isinstance(values, np.ndarray)
            np.testing.assert_array_almost_equal(values, sa_dict['kinetics'][0][obs][param])
    shutil.rmtree(NH3_SA_TEST_DIR, ignore_errors=True)


def teardown_module():
    """
    Safety net that removes any artifacts the tests in this module may have
    left behind, even if a test failed mid-run before its inline cleanup.
    """
    for path in [os.path.join(TEST_DIR, 'log_archive'),
                 os.path.join(MINIMAL_DATA_DIR, 'iteration_1', 'SA'),
                 os.path.join(MINIMAL_DATA_DIR, 'iteration_1', 'PDep_SA'),
                 os.path.join(MINIMAL_DATA_DIR, 'log_archive'),
                 NH3_SA_TEST_DIR]:
        if os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
    for file in [os.path.join(TEST_DIR, 't3.log'),
                 os.path.join(MINIMAL_DATA_DIR, 't3.log')]:
        if os.path.isfile(file):
            os.remove(file)

