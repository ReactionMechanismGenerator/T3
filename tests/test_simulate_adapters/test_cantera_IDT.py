#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_cantera_IDT module
"""

import os
import shutil

from arc.common import almost_equal_lists

from t3.common import SIMULATE_DATA_BASE_PATH
from tests.common import almost_equal, run_minimal
from t3.simulate.cantera_IDT import CanteraIDT, get_t_and_p_lists
from t3.utils.fix_cantera import fix_cantera


TEST_DIR = os.path.join(SIMULATE_DATA_BASE_PATH, 'cantera_simulator_test')
TEST_DIR_IDT = os.path.join(SIMULATE_DATA_BASE_PATH, 'cantera_idt_test')


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


def test_get_cantera_species_label():
    """
    Test the `get_cantera_species_label()` method.
    """
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 2  # CH4 model
    t3.set_paths()
    ct_adapter = CanteraIDT(t3=t3.t3,
                            rmg=t3.rmg,
                            paths=t3.paths,
                            logger=t3.logger,
                            atol=t3.rmg['model']['atol'],
                            rtol=t3.rmg['model']['rtol'],
                            )
    assert ct_adapter.get_cantera_species_label('methane') == 'methane(1)'
    assert ct_adapter.get_cantera_species_label('O2') == 'O2(2)'
    assert ct_adapter.get_cantera_species_label('OH') == 'OH(6)'
    assert ct_adapter.get_cantera_species_label('CH3CH2OO') == 'CH3CH2OO(40)'


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


def test_simulate_seiser():
    """
    Test the ``simulate()`` method for computing IDT using the Seiser model.
    """
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 1  # Seiser model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': 1000, 'P': 1,
                           'termination_rate_ratio': 0.01},
                          ]
    t3.rmg['species'] = [{'label': 'nc7h16', 'smiles': 'CCCCCCC', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [1.0]},
                         {'label': 'o2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxygen'},
                         {'label': 'n2', 'smiles': 'N#N', 'concentration': 0, 'role': 'nitrogen'}]
    ct_adapter = CanteraIDT(t3=t3.t3,
                            rmg=t3.rmg,
                            paths=t3.paths,
                            logger=t3.logger,
                            atol=t3.rmg['model']['atol'],
                            rtol=t3.rmg['model']['rtol'],
                            )
    assert ct_adapter.idt_dict == dict()
    ct_adapter.simulate()
    print(ct_adapter.idt_dict)
    assert almost_equal(ct_adapter.idt_dict[(1.0, 1, 1000)], 0.0328465)


def test_simulate_rmg_ammonia():
    """
    Test the ``simulate()`` method for computing IDT for NH3.
    """
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 0  # NH3 model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': [1500, 1600], 'P': [1, 100],
                           'termination_rate_ratio': 0.01},
                          ]
    t3.rmg['species'] = [{'label': 'NH3', 'smiles': 'N', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [0.5, 1.0]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxygen'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'nitrogen'}]
    ct_adapter = CanteraIDT(t3=t3.t3,
                            rmg=t3.rmg,
                            paths=t3.paths,
                            logger=t3.logger,
                            atol=t3.rmg['model']['atol'],
                            rtol=t3.rmg['model']['rtol'],
                            )
    assert ct_adapter.idt_dict == dict()
    ct_adapter.simulate()
    assert len(ct_adapter.idt_dict.keys()) == 90
    for val in ct_adapter.idt_dict.values():
        assert val is not None
    assert almost_equal(ct_adapter.idt_dict[(1.0, 100.0, 1600.0)], 2.12081e-05)


def test_simulate_rmg_methane():
    """
    Test the ``simulate()`` method for computing IDT for CH4.
    """
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 2  # CH4 model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': [1000, 1200], 'P': 10,
                           'termination_rate_ratio': 0.01},
                          ]
    t3.rmg['species'] = [{'label': 'methane', 'smiles': 'C', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [1.0]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxygen'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'nitrogen'}]
    ct_adapter = CanteraIDT(t3=t3.t3,
                            rmg=t3.rmg,
                            paths=t3.paths,
                            logger=t3.logger,
                            atol=t3.rmg['model']['atol'],
                            rtol=t3.rmg['model']['rtol'],
                            )
    assert ct_adapter.idt_dict == dict()
    ct_adapter.simulate()
    assert len(ct_adapter.idt_dict.keys()) == 15
    for val in ct_adapter.idt_dict.values():
        assert val is not None
    assert almost_equal(ct_adapter.idt_dict[(1.0, 10, 1150.6849315068491)], 0.009160914)


def test_simulate_rmg_heptane():
    """
    Test the ``simulate()`` method for computing IDT for and RMG heptane model.
    """
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 3  # C7H16 model
    t3.set_paths()
    fix_cantera(model_path=t3.paths['cantera annotated'])
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': [1000, 2000], 'P': 10,
                           'termination_rate_ratio': 0.01},
                          ]
    t3.rmg['species'] = [{'label': 'n-heptane', 'smiles': 'CCCCCCC', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [1.0, 2.0]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxygen'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'nitrogen'}]
    ct_adapter = CanteraIDT(t3=t3.t3,
                            rmg=t3.rmg,
                            paths=t3.paths,
                            logger=t3.logger,
                            atol=t3.rmg['model']['atol'],
                            rtol=t3.rmg['model']['rtol'],
                            )
    assert ct_adapter.idt_dict == dict()
    ct_adapter.simulate()
    print(ct_adapter.idt_dict)
    assert len(ct_adapter.idt_dict.keys()) == 30
    for val in ct_adapter.idt_dict.values():
        assert val is not None
    assert almost_equal(ct_adapter.idt_dict[(1.0, 10, 2000.0)], 1.32002925e-08)


def test_get_concentration_combinations():
    """
    Test the ``get_concentration_combinations()`` method.
    """
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 1  # Seiser model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': 1000, 'P': 1,
                           'termination_rate_ratio': 0.01},
                          ]
    t3.rmg['species'] = [{'label': 'nc7h16', 'smiles': 'CCCCCCC', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [1.0]},
                         {'label': 'o2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxygen'},
                         {'label': 'n2', 'smiles': 'N#N', 'concentration': 0, 'role': 'nitrogen'}]
    ct_adapter = CanteraIDT(t3=t3.t3,
                            rmg=t3.rmg,
                            paths=t3.paths,
                            logger=t3.logger,
                            atol=t3.rmg['model']['atol'],
                            rtol=t3.rmg['model']['rtol'],
                            )
    equivalence_ratios, concentration_combinations = ct_adapter.get_concentration_combinations()
    assert equivalence_ratios == [1.0]
    assert concentration_combinations == [{'nc7h16': 1.0, 'o2': 11.0, 'n2': 41.36}]

    t3.iteration = 0  # NH3 model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': [800, 1750], 'P': [1, 100],
                           'termination_rate_ratio': 0.01},
                          ]
    t3.rmg['species'] = [{'label': 'NH3', 'smiles': 'N', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [0.5, 1.0, 1.5]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxygen'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'nitrogen'}]
    ct_adapter = CanteraIDT(t3=t3.t3,
                            rmg=t3.rmg,
                            paths=t3.paths,
                            logger=t3.logger,
                            atol=t3.rmg['model']['atol'],
                            rtol=t3.rmg['model']['rtol'],
                            )
    equivalence_ratios, concentration_combinations = ct_adapter.get_concentration_combinations()
    assert equivalence_ratios == [0.5, 1.0, 1.5]
    assert concentration_combinations == [{'NH3(1)': 1, 'O2(3)': 0.375, 'N2(2)': 1.41},
                                          {'NH3(1)': 1, 'O2(3)': 0.75, 'N2(2)': 2.82},
                                          {'NH3(1)': 1, 'O2(3)': 1.125, 'N2(2)': 4.2299999999999995}]


def teardown_module():
    """
    A method that is run after all unit tests in this class.
    Delete all project directories created during these unit tests
    """
    for test_dir_path in [TEST_DIR, TEST_DIR_IDT]:
        log_archive = os.path.join(test_dir_path, 'log_archive')
        dirs = [log_archive]
        for dir_ in dirs:
            if os.path.isdir(dir_):
                shutil.rmtree(dir_, ignore_errors=True)
        files = [os.path.join(test_dir_path, 't3.log')]
        for file in files:
            if os.path.isfile(file):
                os.remove(file)
        for iteration in [0, 1, 2, 3, 4, 5]:
            figs_path = os.path.join(test_dir_path, f'iteration_{iteration}', 'Figures')
            if os.path.isdir(figs_path):
                shutil.rmtree(figs_path, ignore_errors=True)
