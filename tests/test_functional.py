#!/usr/bin/env python3
# encoding: utf-8

"""
functional test that runs T3's minimal example
"""

import os
import shutil

from arc.common import read_yaml_file

from t3 import T3
from t3.common import DATA_BASE_PATH
from t3.utils.dependencies import check_dependencies


functional_minimal_directory = os.path.join(DATA_BASE_PATH, 'functional_minimal_example')


def test_minimal_example():
    """Tests that the minimal example is functional (ARC is not being called)"""
    input_file = os.path.join(functional_minimal_directory, 'input.yml')
    input_dict = read_yaml_file(path=input_file)
    verbose = 20
    input_dict['verbose'] = verbose
    input_dict['project_directory'] = functional_minimal_directory

    # check that RMG and ARC are available
    check_dependencies()

    # run the minimal example
    t3_object = T3(**input_dict)
    t3_object.execute()

    # check that the minimal example ran to completion by checking for existence of some files and directories
    assert os.path.isfile(os.path.join(functional_minimal_directory, 't3.log'))
    assert os.path.isdir(os.path.join(functional_minimal_directory, 'iteration_0'))
    assert os.path.isfile(os.path.join(functional_minimal_directory, 'iteration_0', 'RMG', 'chemkin', 'species_dictionary.txt'))
    assert os.path.isdir(os.path.join(functional_minimal_directory, 'iteration_1'))
    assert os.path.isfile(os.path.join(functional_minimal_directory, 'iteration_1', 'RMG', 'chemkin', 'species_dictionary.txt'))

    # remove directories created by this functional test
    log_file = os.path.join(functional_minimal_directory, 't3.log')
    if os.path.isfile(log_file):
        os.remove(log_file)
    iteration_directories = [os.path.join(functional_minimal_directory, f'iteration_{i}') for i in range(2)]
    for iteration_dir in iteration_directories:
        if os.path.isdir(iteration_dir):
            shutil.rmtree(iteration_dir)
