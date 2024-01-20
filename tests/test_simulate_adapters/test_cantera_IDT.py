#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_cantera_IDT module
"""

import os
import shutil

from arc.common import almost_equal_lists

from t3.common import SIMULATE_DATA_BASE_PATH
from tests.common import run_minimal
from t3.simulate.cantera_IDT import CanteraIDT, get_t_and_p_lists


TEST_DIR = os.path.join(SIMULATE_DATA_BASE_PATH, 'cantera_simulator_test')


def test_determine_radical_label():
    """
    Test the `determine_radical_label()` method.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    ct_adapter = CanteraIDT(t3=t3.t3,
                            rmg=t3.rmg,
                            paths=t3.paths,
                            logger=t3.logger,
                            atol=t3.rmg['model']['atol'],
                            rtol=t3.rmg['model']['rtol'],
                            )
    label = ct_adapter.determine_radical_label()
    assert label == 'OH(4)'


def test_get_t_and_p_lists():
    """
    Test the `get_t_and_p_lists()` method.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': 1000,
                           'P': 1e0,
                           'termination_conversion': {'ethane': 0.2},
                           'termination_time': [5, 's'],
                           'termination_rate_ratio': 0.01,
                           'conditions_per_iteration': 12},
                          {'type': 'gas batch constant T P',
                           'T': [800, 1750],
                           'P': [1, 100],
                           'termination_conversion': {'ethane': 0.2},
                           'termination_time': [5, 's'],
                           'termination_rate_ratio': 0.01,
                           'conditions_per_iteration': 12},
                          ]
    ct_adapter = CanteraIDT(t3=t3.t3,
                            rmg=t3.rmg,
                            paths=t3.paths,
                            logger=t3.logger,
                            atol=t3.rmg['model']['atol'],
                            rtol=t3.rmg['model']['rtol'],
                            )
    T_list, P_list = get_t_and_p_lists(ct_adapter.rmg['reactors'][0])
    print(T_list, P_list)
    assert T_list == [1000.0]
    assert P_list == [1.0]
    T_list, P_list = get_t_and_p_lists(ct_adapter.rmg['reactors'][1])
    print(T_list, P_list)
    assert almost_equal_lists(T_list, [800.0, 832.27, 867.26, 905.31, 946.86, 992.41, 1042.55, 1098.04, 1159.76,
                                       1228.84, 1306.67, 1395.02, 1496.18, 1613.17, 1750.00], rtol=0.001, atol=0.1)
    assert P_list == [1.0, 10.0, 100.0]


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
