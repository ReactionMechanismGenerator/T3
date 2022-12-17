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
                             'kinetics_libraries': ['primaryH2O2']},
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
                'model': {'core_tolerance': 0.01},
                'rmg_execution_type': 'incore',
                }

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


def test_computing_thermo():
    """Tests computing thermo for two species and running RMG with the updated data"""
    functional_test_directory = os.path.join(DATA_BASE_PATH, 'functional_2_thermo')
    #delete_selective_content_from_test_dirs(test_dir=functional_test_directory)
    input_file = os.path.join(functional_test_directory, 'input.yml')
    input_dict = read_yaml_file(path=input_file)
    input_dict['verbose'] = 20
    input_dict['project_directory'] = functional_test_directory

    # check that RMG and ARC are available
    check_dependencies()

    # run the minimal example
    t3_object = T3(**input_dict)
    t3_object.execute()
    assert os.path.isfile(os.path.join(functional_test_directory, 't3.log'))
    assert os.path.isfile(os.path.join(functional_test_directory, 'species.yml'))
    assert os.path.isfile(os.path.join(functional_test_directory, 'reactions.yml'))
    assert os.path.isdir(os.path.join(functional_test_directory, 'iteration_1'))
    assert os.path.isfile(os.path.join(functional_test_directory, 'iteration_1',
                                       'RMG', 'chemkin', 'species_dictionary.txt'))
    assert os.path.isdir(os.path.join(functional_test_directory, 'iteration_2'))
    assert os.path.isfile(os.path.join(functional_test_directory, 'iteration_2', 'RMG', 'input.py'))
    assert os.path.isfile(os.path.join(functional_test_directory, 'iteration_2', 'RMG', 'RMG.log'))

    expected_lines = [{'line': 'T3 iteration 1:', 'exists': False},
                      {'line': 'Running RMG (tolerance = 0.1, iteration 1)...', 'exists': False},
                      {'line': 'Running RMG (tolerance = 0.1, iteration 2)...', 'exists': False},
                      {'line': 'Running a simulation with SA using RMGConstantTP for 1 conditions...', 'exists': False},
                      {'line': 'Additional calculations required: True', 'exists': False},
                      {'line': 's0_2-propyl  C[CH]C      SA observable', 'exists': False},
                      {'line': 'All species thermodynamic calculations in this iteration successfully converged.', 'exists': False},
                      {'line': 'T3 iteration 2 (just generating a model using RMG):', 'exists': False},
                      {'line': 'SPECIES SUMMARY', 'exists': False},
                      {'line': 'Species for which thermodynamic data was calculate:', 'exists': False},
                      {'line': 'All species calculated by ARC successfully converged', 'exists': False},
                      ]
    with open(os.path.join(functional_test_directory, 't3.log'), 'r') as f:
        lines = f.readlines()
    for line in lines:
        for expected_line_dict in expected_lines:
            if expected_line_dict['line'] in line:
                expected_line_dict['exists'] = True
    print('\n\n********************** test_computing_thermo:\n')
    for expected_line_dict in expected_lines:
        print(expected_line_dict['line'])  # assists in debugging this test, otherwise error messages aren't informative
        assert expected_line_dict['exists'] is True


def delete_selective_content_from_test_dirs(test_dir: str):
    """remove directories created by the test_minimal_example functional test"""
    for (_, dirs, files) in os.walk(test_dir):
        for file_ in files:
            if file_ != 'input.yml' and os.path.isfile(os.path.join(test_dir, file_)):
                os.remove(os.path.join(test_dir, file_))
        for dir_ in dirs:
            if os.path.isdir(os.path.join(test_dir, dir_)):
                shutil.rmtree(os.path.join(test_dir, dir_), ignore_errors=True)
        break


def teardown_module():
    """teardown any state that was previously setup with a setup_module method."""
    test_dirs_to_selectively_delete = [os.path.join(DATA_BASE_PATH, 'functional_2_thermo'),
                                       ]
    for test_dir in test_dirs_to_selectively_delete:
        delete_selective_content_from_test_dirs(test_dir)
    test_dirs = [os.path.join(DATA_BASE_PATH, 'T3_functional_test_1'),
                 ]
    for test_dir in test_dirs:
        shutil.rmtree(test_dir, ignore_errors=True)
