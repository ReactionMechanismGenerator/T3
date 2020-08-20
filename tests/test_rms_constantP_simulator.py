#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_rms_constantP_simulator module
"""

import os
import shutil

from t3.common import DATA_BASE_PATH
from tests.common import run_minimal
from t3.simulate.rms_constantP_simulator import RMSConstantP


TEST_DIR = os.path.join(DATA_BASE_PATH, 'rms_constantP_simulator_test')


def test_set_up_no_sa():
    """
    Run RMG's minimal example without SA by testing the `set_up` method within the RMSConstantP init method.
    By setting observable_list = list(), no SA is run. Instead, RMS is used to just simulate the mechanism.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()

    RMSSimulator_adapter = RMSConstantP(t3=t3.t3,
                                        rmg=t3.rmg,
                                        paths=t3.paths,
                                        logger=t3.logger,
                                        atol=t3.rmg['model']['atol'],
                                        rtol=t3.rmg['model']['rtol'],
                                        observable_list=list(),
                                        )
    # check that there are over 20 time steps (as an arbitrary number) to indicate that the solver completed
    assert len(RMSSimulator_adapter.sol.t) > 20


def test_get_sa_coefficients():
    """
    Run RMG's minimal example with SA by testing the `set_up` method within the RMSConstantP init method.
    Then run the `get_sa_coefficients()` method to test that RMSConstantP correctly parses the bsol object
    to obtain sa_dict.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    observable_list = ['OH', 'H']
    RMSSimulator_adapter = RMSConstantP(t3=t3.t3,
                                        rmg=t3.rmg,
                                        paths=t3.paths,
                                        logger=t3.logger,
                                        atol=t3.rmg['model']['atol'],
                                        rtol=t3.rmg['model']['rtol'],
                                        observable_list=observable_list,
                                        sa_atol=t3.t3['sensitivity']['atol'],
                                        sa_rtol=t3.t3['sensitivity']['atol'],
                                        )
    sa_dict = RMSSimulator_adapter.get_sa_coefficients()

    # check that there are over 20 time steps (as an arbitrary number) to indicate that the solver completed
    assert len(sa_dict['time']) > 20


def teardown_module():
    """
    A method that is run after all unit tests in this class.
    Delete all project directories created during these unit tests
    """
    log_archive = os.path.join(TEST_DIR, 'log_archive')
    dirs = [log_archive]
    for dir in dirs:
        if os.path.isdir(dir):
            shutil.rmtree(dir)
    files = [os.path.join(TEST_DIR, 't3.log')]
    for file in files:
        if os.path.isfile(file):
            os.remove(file)
