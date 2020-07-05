#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_rms_constantTP_simulator module
"""

import os
import shutil

from t3.common import DATA_BASE_PATH
from tests.common import run_minimal
from t3.simulate.rms_constantTP_simulator import RMSConstantTP


TEST_DIR = os.path.join(DATA_BASE_PATH, 'rms_constantTP_simulator_test')


def test_set_up_no_sa():
    """
    Runs RMG's minimal example without SA by testing the `set_up` method in the RMSConstantTP class.
    By setting observable_list = list(), no SA is run. Instead, RMS is used to just simulate the mechanism.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()

    RMSSimulator_adapter = RMSConstantTP(t3=t3.t3,
                                         rmg=t3.rmg,
                                         paths=t3.paths,
                                         logger=t3.logger,
                                         atol=t3.rmg['model']['atol'],
                                         rtol=t3.rmg['model']['rtol'],
                                         observable_list=list(),
                                         )
    # check that there are over 100 time steps (as an arbitrary number) to indicate that the solver completed
    assert len(RMSSimulator_adapter.sol.t) > 100


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
