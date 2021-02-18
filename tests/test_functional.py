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


def test_no_t3_no_qm():
    """Test proper execution of T3 without specifying neither of the t3 nor qm args"""
    rmg_args = {'database': {'thermo_libraries': ['primaryThermoLibrary', 'BurkeH2O2'],
                             'kinetics_libraries': ['BurkeH2O2inN2']},
                'species': [{'label': 'H2',
                             'smiles': '[H][H]',
                             'concentration': 0.67},
                            {'label': 'O2',
                             'smiles': '[O][O]',
                             'concentration': 0.33}],
                'reactors': [{'type': 'gas batch constant T P',
                              'T': 1000,
                              'P': 1,
                              'termination_conversion': {'H2': 0.1},
                              'termination_time': [1, 'ms']}],
                'model': {'core_tolerance': 0.01}}

    t3_object = T3(project='T3_functional_test_1',
                   rmg=rmg_args,
                   clean_dir=True,
                   )
    t3_object.execute()
    with open(os.path.join(t3_object.project_directory, 't3.log'), 'r') as f:
        lines = f.readlines()
    for line in ['#                  The   Tandem   Tool   (T3)                  #\n',
                 'Starting project T3_functional_test_1\n',
                 'T3 iteration 1:\n',
                 ]:
        assert line in lines
    assert not any('T3 iteration 0:\n' in line for line in lines)
    assert any('T3 execution terminated' in line for line in lines)
    shutil.rmtree(t3_object.project_directory, ignore_errors=True)


def test_minimal_example():
    """Tests that the minimal example is functional (ARC is not being called)"""
    delete_minimal_example_dirs()
    functional_minimal_directory = os.path.join(DATA_BASE_PATH, 'functional_minimal_example')
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
    assert os.path.isdir(os.path.join(functional_minimal_directory, 'iteration_1'))
    assert os.path.isfile(os.path.join(functional_minimal_directory, 'iteration_1',
                                       'RMG', 'chemkin', 'species_dictionary.txt'))
    assert os.path.isdir(os.path.join(functional_minimal_directory, 'iteration_2'))
    assert os.path.isfile(os.path.join(functional_minimal_directory, 'iteration_2',
                                       'RMG', 'chemkin', 'species_dictionary.txt'))
    delete_minimal_example_dirs()


def delete_minimal_example_dirs():
    """remove directories created by the test_minimal_example functional test"""
    functional_minimal_directory = os.path.join(DATA_BASE_PATH, 'functional_minimal_example')
    log_file = os.path.join(functional_minimal_directory, 't3.log')
    if os.path.isfile(log_file):
        os.remove(log_file)
    iteration_directories = [os.path.join(functional_minimal_directory, f'iteration_{i + 1}') for i in range(2)]
    for iteration_dir in iteration_directories:
        if os.path.isdir(iteration_dir):
            shutil.rmtree(iteration_dir, ignore_errors=True)
