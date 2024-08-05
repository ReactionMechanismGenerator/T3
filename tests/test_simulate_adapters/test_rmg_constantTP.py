#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_rmg_constantTP module
"""

import os
import shutil

from arc.common import read_yaml_file

from t3.common import SIMULATE_DATA_BASE_PATH
from tests.common import run_minimal
from t3.main import T3
from t3.simulate.rmg_constantTP import RMGConstantTP


TEST_DIR_1 = os.path.join(SIMULATE_DATA_BASE_PATH, 'rmg_simulator_test')
TEST_DIR_2 = os.path.join(SIMULATE_DATA_BASE_PATH, 'rmg_simulator_test_ranges_3a')
TEST_DIR_3 = os.path.join(SIMULATE_DATA_BASE_PATH, 'rmg_simulator_test_ranges_3b')
TEST_DIR_4 = os.path.join(SIMULATE_DATA_BASE_PATH, 'rmg_simulator_test_ranges_3c')
RANGED_INPUT_DICT_1 = {'verbose': 10, 'project_directory': TEST_DIR_2,
                       'project': 'test_get_species_concentration_lists_from_ranged_params_case_3a',
                       't3': {'options': {'max_T3_iterations': 2, 'max_RMG_walltime': '00:00:05:00'},
                              'sensitivity': {'adapter': 'RMGConstantTP', 'top_SA_species': 10}},
                       'rmg': {'database': {'thermo_libraries': ['primaryThermoLibrary'], 'kinetics_libraries': []},
                               'reactors': [{'type': 'gas batch constant T P',
                                             'T': 1000, 'P': 1, 'termination_rate_ratio': 0.1}],
                               'model': {'core_tolerance': [0.01, 0.001]}},
                       'qm': {'adapter': 'ARC', 'level_of_theory': 'b3lyp/6-31g(d,p)'}}
RANGED_INPUT_DICT_2 = {'verbose': 10, 'project_directory': TEST_DIR_3,
                       'project': 'test_get_species_concentration_lists_from_ranged_params_case_3b',
                       't3': {'options': {'max_T3_iterations': 2, 'max_RMG_walltime': '00:00:05:00'},
                              'sensitivity': {'adapter': 'RMGConstantTP', 'top_SA_species': 10}},
                       'rmg': {'database': {'thermo_libraries': ['primaryThermoLibrary'], 'kinetics_libraries': []},
                               'reactors': [{'type': 'gas batch constant T P',
                                             'T': 1000, 'P': 1, 'termination_rate_ratio': 0.1}],
                               'model': {'core_tolerance': [0.01, 0.001]}},
                       'qm': {'adapter': 'ARC', 'level_of_theory': 'b3lyp/6-31g(d,p)'}}


def test_set_up_no_sa():
    """
    Run RMG's minimal example without SA by testing the `set_up` method within the RMGSimulator init method.
    By setting observable_list = list(), no SA is run. Instead, RMG is just used to simulate the mechanism.
    """
    t3 = run_minimal(project_directory=TEST_DIR_1)
    t3.set_paths()

    rmg_simulator_adapter = RMGConstantTP(t3=t3.t3,
                                          rmg=t3.rmg,
                                          paths=t3.paths,
                                          logger=t3.logger,
                                          atol=t3.rmg['model']['atol'],
                                          rtol=t3.rmg['model']['rtol'],
                                          observable_list=list(),
                                          )
    rmg_simulator_adapter.simulate()
    # check that RMG produced the concentration profiles
    concentration_profiles = os.path.isfile(os.path.join(t3.paths['RMG'], 'solver', 'simulation_1_12.csv'))
    assert concentration_profiles


def test_get_sa_coefficients():
    """
    Run RMG's minimal example with SA by testing the `set_up` method within the RMGConstantTP init method.
    Then run the `get_sa_coefficients()` method to test that RMGConstantTP correctly parses the SA csv files
    to obtain sa_dict.
    """
    t3 = run_minimal(project_directory=TEST_DIR_1)
    t3.set_paths()
    observable_list = ['OH', 'H']
    rmg_simulator_adapter = RMGConstantTP(t3=t3.t3,
                                          rmg=t3.rmg,
                                          paths=t3.paths,
                                          logger=t3.logger,
                                          atol=t3.rmg['model']['atol'],
                                          rtol=t3.rmg['model']['rtol'],
                                          observable_list=observable_list,
                                          sa_atol=t3.t3['sensitivity']['atol'],
                                          sa_rtol=t3.t3['sensitivity']['rtol'],
                                          )
    rmg_simulator_adapter.simulate()
    # check that RMG ran SA
    performed_sa = os.path.isdir(t3.paths['SA'])
    assert performed_sa
    # check that RMG produced the corresponding csv files
    files = ['sensitivity_1_SPC_3.csv', 'sensitivity_1_SPC_4.csv', 'simulation_1_12.csv']
    for file in files:
        exist = os.path.isfile(os.path.join(t3.paths['SA solver'], file))
        assert exist
    # check that RMG can correctly parse the csv files to product the SA dictionary
    sa_dict = rmg_simulator_adapter.get_sa_coefficients()
    # check that there are over 100 time steps (as an arbitrary number) to indicate that the solver completed
    assert len(sa_dict['time']) > 100
    # check that there are SA data for the 2 requested species
    assert len(sa_dict['kinetics']) == 2
    assert len(sa_dict['thermo']) == 2


def test_get_idt_by_T():
    """
    Calculate the ignition delay time for RMG's minimal example.
    Since this adapter simulates at constant T, this method returns a dictionary whose values are empty lists.
    """
    t3 = run_minimal(project_directory=TEST_DIR_1)
    t3.set_paths()
    rmg_simulator_adapter = RMGConstantTP(t3=t3.t3,
                                          rmg=t3.rmg,
                                          paths=t3.paths,
                                          logger=t3.logger,
                                          atol=t3.rmg['model']['atol'],
                                          rtol=t3.rmg['model']['rtol'],
                                          observable_list=list(),
                                          )
    rmg_simulator_adapter.simulate()
    idt_dict = rmg_simulator_adapter.get_idt_by_T()
    assert len(idt_dict['idt']) == 0
    assert len(idt_dict['idt_index']) == 0


def test_get_species_concentration_lists_from_ranged_params_case_3():
    """
    Test the get_species_concentration_lists_from_ranged_params() method
    For case 3 where ``modify_concentration_ranges_together`` is set to ``True``
    """
    input_dict = RANGED_INPUT_DICT_1.copy()
    input_dict['rmg']['species'] = [
        {'label': 'fuel', 'smiles': 'CC', 'concentration': 1, 'reactive': True, 'balance': False},
        {'label': 'O2', 'smiles': '[O][O]', 'concentration': [3.75, 11.25], 'reactive': True, 'balance': False},
        {'label': 'N2', 'smiles': 'N#N', 'concentration': [14.1, 42.3], 'reactive': False, 'balance': False}]
    input_dict['t3']['options']['modify_concentration_ranges_together'] = True
    input_dict['t3']['options']['modify_concentration_ranges_in_reverse'] = False
    t3 = T3(**input_dict)
    t3.iteration = 1
    t3.set_paths()
    rmg_simulator_adapter = RMGConstantTP(t3=t3.t3,
                                          rmg=t3.rmg,
                                          paths=t3.paths,
                                          logger=t3.logger,
                                          atol=t3.rmg['model']['atol'],
                                          rtol=t3.rmg['model']['rtol'],
                                          observable_list=list(),
                                          )
    species_lists = rmg_simulator_adapter.get_species_concentration_lists_from_ranged_params()
    assert species_lists == [[{'concentration': 14.1, 'label': 'N2'},
                              {'concentration': 3.75, 'label': 'O2'},
                              {'concentration': 1, 'label': 'fuel'}],
                             [{'concentration': 28.199999999999996, 'label': 'N2'},
                              {'concentration': 7.5, 'label': 'O2'},
                              {'concentration': 1, 'label': 'fuel'}],
                             [{'concentration': 42.3, 'label': 'N2'},
                              {'concentration': 11.25, 'label': 'O2'},
                              {'concentration': 1, 'label': 'fuel'}]]

    input_dict = RANGED_INPUT_DICT_2.copy()
    input_dict['rmg']['species'] = [
        {'label': 'FA', 'smiles': 'OC=O', 'concentration': [1, 4], 'reactive': True, 'balance': False},
        {'label': 'O2', 'smiles': '[O][O]', 'concentration': 1, 'reactive': True, 'balance': False},
        {'label': 'N2', 'smiles': 'N#N', 'concentration': 3.76, 'reactive': False, 'balance': False}]
    input_dict['t3']['options']['modify_concentration_ranges_together'] = True
    input_dict['t3']['options']['modify_concentration_ranges_in_reverse'] = False
    t3 = T3(**input_dict)
    t3.iteration = 1
    t3.set_paths()
    rmg_simulator_adapter = RMGConstantTP(t3=t3.t3,
                                          rmg=t3.rmg,
                                          paths=t3.paths,
                                          logger=t3.logger,
                                          atol=t3.rmg['model']['atol'],
                                          rtol=t3.rmg['model']['rtol'],
                                          observable_list=list(),
                                          )
    species_lists = rmg_simulator_adapter.get_species_concentration_lists_from_ranged_params()
    assert species_lists == [[{'label': 'N2', 'concentration': 3.76},
                              {'label': 'O2', 'concentration': 1.0},
                              {'label': 'FA', 'concentration': 1.0}],
                             [{'label': 'N2', 'concentration': 3.76},
                              {'label': 'FA', 'concentration': 2.5},
                              {'label': 'O2', 'concentration': 1.0}],
                             [{'label': 'FA', 'concentration': 4.0},
                              {'label': 'N2', 'concentration': 3.76},
                              {'label': 'O2', 'concentration': 1.0}]]

    input_dict = RANGED_INPUT_DICT_1.copy()
    input_dict['rmg']['species'] = [
        {'label': 'fuel', 'smiles': 'CC', 'concentration': 1, 'reactive': True, 'balance': False},
        {'label': 'O2', 'smiles': '[O][O]', 'concentration': [3.75, 11.25], 'reactive': True, 'balance': False},
        {'label': 'N2', 'smiles': 'N#N', 'concentration': [14.1, 42.3], 'reactive': False, 'balance': False}]
    input_dict['t3']['options']['modify_concentration_ranges_together'] = False
    input_dict['t3']['options']['modify_concentration_ranges_in_reverse'] = False
    t3 = T3(**input_dict)
    t3.iteration = 1
    t3.set_paths()
    rmg_simulator_adapter = RMGConstantTP(t3=t3.t3,
                                          rmg=t3.rmg,
                                          paths=t3.paths,
                                          logger=t3.logger,
                                          atol=t3.rmg['model']['atol'],
                                          rtol=t3.rmg['model']['rtol'],
                                          observable_list=list(),
                                          )
    species_lists = rmg_simulator_adapter.get_species_concentration_lists_from_ranged_params()
    assert len(species_lists) == 9
    assert len(species_lists[0]) == 19


def test_run_sa_via_rmg():
    """Test running sensitivity analysis via RMG"""
    input_dict = read_yaml_file(os.path.join(TEST_DIR_4, 'input.yml'))
    input_dict['project_directory'] = TEST_DIR_4
    t3_object = T3(**input_dict)
    t3_object.iteration = 1
    t3_object.set_paths()
    t3_object.sa_observables = ['FA']
    simulate_adapter = RMGConstantTP(t3=t3_object.t3,
                                     rmg=t3_object.rmg,
                                     paths=t3_object.paths,
                                     logger=t3_object.logger,
                                     atol=t3_object.rmg['model']['atol'],
                                     rtol=t3_object.rmg['model']['rtol'],
                                     observable_list=t3_object.sa_observables,
                                     sa_atol=t3_object.t3['sensitivity']['atol'],
                                     sa_rtol=t3_object.t3['sensitivity']['rtol'],
                                     global_observables=None,
                                     )
    simulate_adapter.simulate()
    t3_object.sa_dict = simulate_adapter.get_sa_coefficients()
    assert list(t3_object.sa_dict.keys()) == ['kinetics', 'thermo', 'time']
    assert list(t3_object.sa_dict['kinetics'].keys()) == ['FA(1)']
    assert list(t3_object.sa_dict['thermo'].keys()) == ['FA(1)']


def teardown_module():
    """
    A method that is run after all unit tests in this class.
    Delete all project directories created during these unit tests
    """
    for test_dir in [TEST_DIR_1, TEST_DIR_2, TEST_DIR_3, TEST_DIR_4]:
        for i in [0, 1]:
            solver_directory = os.path.join(test_dir, f'iteration_{i}', 'RMG', 'solver')
            species_directory = os.path.join(test_dir, f'iteration_{i}', 'RMG', 'species')
            sa_directory = os.path.join(test_dir, f'iteration_{i}', 'SA')
            log_archive = os.path.join(test_dir, 'log_archive')
            dirs = [solver_directory, species_directory, sa_directory, log_archive]
            for dir in dirs:
                if os.path.isdir(dir):
                    shutil.rmtree(dir, ignore_errors=True)
            files = [os.path.join(test_dir, 't3.log')]
            for file in files:
                if os.path.isfile(file):
                    os.remove(file)
