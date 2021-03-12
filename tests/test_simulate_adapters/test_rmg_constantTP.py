#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_rmg_constantTP module
"""

import os
import shutil

from t3.common import SIMULATE_DATA_BASE_PATH
from tests.common import run_minimal
from t3.simulate.rmg_constantTP import RMGConstantTP


TEST_DIR = os.path.join(SIMULATE_DATA_BASE_PATH, 'rmg_simulator_test')


def test_set_up_no_sa():
    """
    Run RMG's minimal example without SA by testing the `set_up` method within the RMGSimulator init method.
    By setting observable_list = list(), no SA is run. Instead, RMG is just used to simulate the mechanism.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
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
    t3 = run_minimal(project_directory=TEST_DIR)
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
    t3 = run_minimal(project_directory=TEST_DIR)
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


def teardown_module():
    """
    A method that is run after all unit tests in this class.
    Delete all project directories created during these unit tests
    """
    solver_directory = os.path.join(TEST_DIR, 'iteration_0', 'RMG', 'solver')
    species_directory = os.path.join(TEST_DIR, 'iteration_0', 'RMG', 'species')
    sa_directory = os.path.join(TEST_DIR, 'iteration_0', 'SA')
    log_archive = os.path.join(TEST_DIR, 'log_archive')
    dirs = [solver_directory, species_directory, sa_directory, log_archive]
    for dir in dirs:
        if os.path.isdir(dir):
            shutil.rmtree(dir, ignore_errors=True)
    files = [os.path.join(TEST_DIR, 't3.log')]
    for file in files:
        if os.path.isfile(file):
            os.remove(file)
