#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_rms_constantHP module
"""

import os
import shutil

from t3.common import SIMULATE_DATA_BASE_PATH
from tests.common import run_minimal
from t3.simulate.rms_constantHP import RMSConstantHP


TEST_DIR = os.path.join(SIMULATE_DATA_BASE_PATH, 'rms_simulator_test')


def test_set_up_no_sa():
    """
    Run RMG's minimal example without SA by testing the `set_up` method within the RMSSimulator init method.
    By setting observable_list = list(), no SA is run. Instead, RMS is used to just simulate the mechanism.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()

    rms_simulator_adapter = RMSConstantHP(t3=t3.t3,
                                          rmg=t3.rmg,
                                          paths=t3.paths,
                                          logger=t3.logger,
                                          atol=t3.rmg['model']['atol'],
                                          rtol=t3.rmg['model']['rtol'],
                                          observable_list=list(),
                                          )
    rms_simulator_adapter.simulate()
    # check that there are over 20 time steps (as an arbitrary number) to indicate that the solver completed
    assert len(rms_simulator_adapter.sol.t) > 20


def test_get_sa_coefficients():
    """
    Run RMG's minimal example with SA.
    Then run the `get_sa_coefficients()` method to test that the adapter correctly parses the bsol object
    to obtain sa_dict.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    observable_list = ['OH', 'H']
    rms_simulator_adapter = RMSConstantHP(t3=t3.t3,
                                          rmg=t3.rmg,
                                          paths=t3.paths,
                                          logger=t3.logger,
                                          atol=t3.rmg['model']['atol'],
                                          rtol=t3.rmg['model']['rtol'],
                                          observable_list=observable_list,
                                          sa_atol=t3.t3['sensitivity']['atol'],
                                          sa_rtol=t3.t3['sensitivity']['rtol'],
                                          )
    rms_simulator_adapter.simulate()
    sa_dict = rms_simulator_adapter.get_sa_coefficients()

    # check that there are over 20 time steps (as an arbitrary number) to indicate that the solver completed
    assert len(sa_dict['time']) > 20
    # check that there are SA data for the 2 requested species
    assert len(sa_dict['kinetics']) == 2
    assert len(sa_dict['thermo']) == 2


def test_get_idt_by_T():
    """
    Calculate the ignition delay time for RMG's minimal example.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    rms_simulator_adapter = RMSConstantHP(t3=t3.t3,
                                          rmg=t3.rmg,
                                          paths=t3.paths,
                                          logger=t3.logger,
                                          atol=t3.rmg['model']['atol'],
                                          rtol=t3.rmg['model']['rtol'],
                                          observable_list=list(),
                                          sa_atol=t3.t3['sensitivity']['atol'],
                                          sa_rtol=t3.t3['sensitivity']['rtol'],
                                          )
    rms_simulator_adapter.simulate()
    idt_dict = rms_simulator_adapter.get_idt_by_T()
    assert len(idt_dict['idt']) == 1
    assert len(idt_dict['idt_index']) == 1


def teardown_module():
    """
    A method that is run after all unit tests in this class.
    Delete all project directories created during these unit tests
    """
    log_archive = os.path.join(TEST_DIR, 'log_archive')
    dirs = [log_archive]
    for dir in dirs:
        if os.path.isdir(dir):
            shutil.rmtree(dir, ignore_errors=True)
    files = [os.path.join(TEST_DIR, 't3.log')]
    for file in files:
        if os.path.isfile(file):
            os.remove(file)
