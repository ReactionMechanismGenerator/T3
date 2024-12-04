#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_cantera_IDT module
"""

import os
import shutil

import cantera as ct

from arc.common import almost_equal_lists, read_yaml_file

from t3.common import SIMULATE_TEST_DATA_BASE_PATH, TEST_DATA_BASE_PATH
from tests.common import almost_equal, run_minimal
from t3.simulate.cantera_IDT import (CanteraIDT, DELTA_H, DELTA_K, calculate_arrhenius_rate_coefficient,
                                     calculate_troe_rate_coefficient, calculate_chebyshev_rate_coefficient,
                                     calculate_plog_rate_coefficient, get_Ea_units, get_h298, get_pressure_from_cantera,
                                     get_t_and_p_lists, get_top_sa_coefficients, perturb_enthalpy,
                                     perturb_reaction_rate_coefficient, plot_idt_vs_temperature)
from t3.utils.fix_cantera import fix_cantera


TEST_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'cantera_simulator_test')
TEST_DIR_IDT = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'cantera_idt_test')


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
    assert T_list == [1000.0]
    assert P_list == [1.0]
    T_list, P_list = get_t_and_p_lists(ct_adapter.rmg['reactors'][1])
    assert almost_equal_lists(T_list, [800.0, 818.51, 837.91, 858.24, 879.58, 902.01, 925.62, 950.495, 976.74, 1004.48,
                                       1033.85, 1064.98, 1098.04, 1133.22, 1170.73, 1210.81, 1253.73, 1299.81, 1349.40,
                                       1402.92, 1460.87, 1523.81, 1592.42, 1667.49, 1750.00], rtol=0.001, atol=0.1)
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
    assert ct_adapter.reactor_idt_dict is None
    ct_adapter.simulate()
    assert almost_equal(ct_adapter.reactor_idt_dict[0][1.0][1.0][1000.0], 0.0328465)
    os.remove(t3.paths['SA IDT dict']) if os.path.exists(t3.paths['SA IDT dict']) else None
    os.remove(t3.paths['SA IDT dict top X']) if os.path.exists(t3.paths['SA IDT dict top X']) else None


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
    assert ct_adapter.reactor_idt_dict is None
    ct_adapter.simulate()
    assert list(ct_adapter.reactor_idt_dict[0].keys()) == [0.5, 1.0]
    assert list(ct_adapter.reactor_idt_dict[0][0.5].keys()) == [1.0, 10.0, 100.0]
    assert almost_equal_lists(ct_adapter.reactor_idt_dict[0][0.5][10.0].keys(), [1500.0, 1510.49, 1521.13, 1531.91, 1542.86,
                                                                      1553.96, 1565.22, 1576.64, 1588.24, 1600.0],
                              rtol=0.001, atol=0.1)
    assert abs(ct_adapter.reactor_idt_dict[0][0.5][10.0][1600.0] - 0.000182883) < 1e-7
    for val in ct_adapter.reactor_idt_dict[0][0.5][10.0].values():
        assert val is not None


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
    assert ct_adapter.reactor_idt_dict is None
    ct_adapter.simulate()
    assert list(ct_adapter.reactor_idt_dict[0].keys()) == [1.0]
    assert list(ct_adapter.reactor_idt_dict[0][1.0].keys()) == [10.0]
    assert almost_equal_lists(ct_adapter.reactor_idt_dict[0][1.0][10.0].keys(),
                              [1000.0, 1008.85, 1017.86, 1027.03, 1036.36, 1045.87, 1055.56, 1065.42, 1075.47,
                               1085.71, 1096.15, 1106.80, 1117.65, 1128.71, 1140.0, 1151.52, 1163.27, 1175.26,
                               1187.5, 1200.0], rtol=0.001, atol=0.1)
    for val in ct_adapter.reactor_idt_dict[0].values():
        assert val is not None
    assert almost_equal(ct_adapter.reactor_idt_dict[0][1.0][10.0][1200.0], 0.00526738418)


def test_simulate_rmg_methane_reduced():
    """
    Test the ``simulate()`` method for computing IDT for the reduced CH4 model.
    """
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 4  # Reduced CH4 model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': [800, 2000], 'P': 10,
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
    assert ct_adapter.reactor_idt_dict is None
    ct_adapter.simulate()
    assert list(ct_adapter.reactor_idt_dict[0].keys()) == [1.0]
    assert list(ct_adapter.reactor_idt_dict[0][1.0].keys()) == [10.0]
    assert almost_equal_lists(ct_adapter.reactor_idt_dict[0][1.0][10.0].keys(),
                              [820.5128, 842.1052, 864.8648, 888.88888, 914.28571, 941.176470, 969.6969, 1000.0,
                               1032.2580, 1066.66666, 1103.4482, 1142.85714, 1185.1851, 1230.76923, 1280.0, 1333.33333,
                               1391.30434, 1454.54545, 1523.8095, 1600.0, 1684.2105, 1777.77777, 1882.3529, 2000.0], rtol=0.001, atol=0.1)
    for val in ct_adapter.reactor_idt_dict[0].values():
        assert val is not None
    assert almost_equal(ct_adapter.reactor_idt_dict[0][1.0][10.0][1600.0], 2.8396e-05)


def test_get_sa_coefficients():
    """
    Test the ``get_sa_coefficients()`` method.
    """
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 4  # A reduced CH4 model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': [800, 2000], 'P': 10,
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
    sa_dict = ct_adapter.get_sa_coefficients()
    assert len(sa_dict) == 2
    assert len(sa_dict['kinetics']) == 1
    assert almost_equal_lists(sa_dict['kinetics']['IDT'][0][1.0][10.0].keys(),
                              [820.5128, 842.1052, 864.8648, 888.88888, 914.28571, 941.176470, 969.6969, 1000.0,
                               1032.2580, 1066.66666, 1103.4482, 1142.85714, 1185.1851, 1230.76923, 1280.0, 1333.33333,
                               1391.30434, 1454.54545, 1523.8095, 1600.0, 1684.2105, 1777.77777, 1882.3529, 2000.0], rtol=0.001, atol=0.1)
    assert almost_equal(sa_dict['kinetics']['IDT'][0][1.0][10.0][1000.0][68], -0.31029061953466325)
    assert almost_equal(sa_dict['kinetics']['IDT'][0][1.0][10.0][1000.0][110], -0.2865901229666261)
    assert almost_equal(sa_dict['kinetics']['IDT'][0][1.0][10.0][1000.0][77], 0.15626917713911037)
    assert almost_equal_lists(sa_dict['thermo']['IDT'][0][1.0][10.0].keys(),
                              [820.5128, 842.1052, 864.8648, 888.88888, 914.28571, 941.176470, 969.6969, 1000.0,
                               1032.2580, 1066.66666, 1103.4482, 1142.85714, 1185.1851, 1230.76923, 1280.0, 1333.33333,
                               1391.30434, 1454.54545, 1523.8095, 1600.0, 1684.2105, 1777.77777, 1882.3529, 2000.0], rtol=0.001, atol=0.1)
    assert almost_equal(sa_dict['thermo']['IDT'][0][1.0][10.0][1000.0][4], -0.16204949148906997)
    assert almost_equal(sa_dict['thermo']['IDT'][0][1.0][10.0][1000.0][11], 0.14341172469660043)
    assert almost_equal(sa_dict['thermo']['IDT'][0][1.0][10.0][1000.0][12], 0.0012104558725444484)
    assert os.path.isfile(t3.paths['SA IDT dict'])
    assert os.path.isfile(t3.paths['SA IDT dict top X'])
    os.remove(t3.paths['SA IDT dict']) if os.path.exists(t3.paths['SA IDT dict']) else None
    os.remove(t3.paths['SA IDT dict top X']) if os.path.exists(t3.paths['SA IDT dict top X']) else None


def test_simulate_rmg_heptane():
    """
    Test the ``simulate()`` method for computing IDT for and RMG heptane model.
    """
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 3  # C7H16 model
    t3.set_paths()
    fix_cantera(model_path=t3.paths['cantera annotated'])
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': 1500, 'P': 10,
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
    assert ct_adapter.reactor_idt_dict is None
    ct_adapter.simulate()
    assert almost_equal(ct_adapter.reactor_idt_dict[0][1.0][10.0][1500.0], 8.8617387e-08)
    assert almost_equal(ct_adapter.reactor_idt_dict[0][2.0][10.0][1500.0], 1.6509341e-07)


def test_get_concentration_combinations():
    """
    Test the ``get_concentration_combinations()`` method.
    """
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 1  # Seiser test model
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


def test_perturb_enthalpy():
    """
    Test the ``perturb_enthalpy()`` method.
    """
    model_path = os.path.join(TEST_DATA_BASE_PATH, 'models', 'N2H4.yaml')
    perturbed_model_path = os.path.join(TEST_DATA_BASE_PATH, 'models', 'N2H4_perturbed.yaml')
    model = ct.Solution(model_path)
    assert almost_equal(model.species('N2(2)').thermo.h(298), -0.0405687)
    assert almost_equal(model.species('N2(2)').thermo.s(298), 191501.4639)
    assert almost_equal(model.species('N2(2)').thermo.cp(298), 29083.8382)
    success = perturb_enthalpy(original_path=model_path, perturbed_path=perturbed_model_path, species_index=3)
    assert success is True
    perturbed_model = ct.Solution(perturbed_model_path)
    assert almost_equal(perturbed_model.species('N2(2)').thermo.h(298), 99999.959431)
    assert almost_equal(perturbed_model.species('N2(2)').thermo.s(298), 191501.4639)
    assert almost_equal(perturbed_model.species('N2(2)').thermo.cp(298), 29083.8382)
    assert almost_equal(perturbed_model.species("N2(2)").thermo.h(298) - model.species("N2(2)").thermo.h(298), 1e5)
    assert almost_equal(perturbed_model.species("N2(2)").thermo.h(2000) - model.species("N2(2)").thermo.h(2000), 1e5)
    assert almost_equal(perturbed_model.species("N2(2)").thermo.s(298) - model.species("N2(2)").thermo.s(298), 0.0)
    assert almost_equal(perturbed_model.species("N2(2)").thermo.s(2000) - model.species("N2(2)").thermo.s(2000), 0.0)
    assert almost_equal(perturbed_model.species("N2(2)").thermo.cp(298) - model.species("N2(2)").thermo.cp(298), 0.0)
    assert almost_equal(perturbed_model.species("N2(2)").thermo.cp(2000) - model.species("N2(2)").thermo.cp(2000), 0.0)

    original, modified = 87426041.87140828, 87526041.8714083
    assert almost_equal(model.species('H4N2(1)').thermo.h(0.1), original)
    success = perturb_enthalpy(original_path=model_path, perturbed_path=perturbed_model_path, species_index=2)
    assert success is True
    perturbed_model = ct.Solution(perturbed_model_path)
    assert almost_equal(perturbed_model.species('H4N2(1)').thermo.h(0.1), modified)
    assert almost_equal(modified - original, DELTA_H * 1e6)  # this is 0.1 kJ/mol in J/kmol units
    os.remove(perturbed_model_path) if os.path.exists(perturbed_model_path) else None


def test_perturb_reaction_rate_coefficient():
    """
    Test the ``perturb_reaction_rate_coefficient()`` method.
    """
    model_path = os.path.join(TEST_DATA_BASE_PATH, 'models', 'N2H4.yaml')
    perturbed_model_path = os.path.join(TEST_DATA_BASE_PATH, 'models', 'N2H4_perturbed.yaml')
    original_model = ct.Solution(model_path)

    # Test Simple Arrhenius
    reaction_index = 1
    original_A = original_model.reaction(reaction_index).rate.pre_exponential_factor
    success = perturb_reaction_rate_coefficient(model_path, perturbed_model_path, reaction_index)
    assert success is True
    perturbed_model = ct.Solution(perturbed_model_path)
    perturbed_A = perturbed_model.reaction(reaction_index).rate.pre_exponential_factor
    assert almost_equal(perturbed_A, original_A * (1 + DELTA_K))

    # Test Falloff reaction
    reaction_index = 2
    original_low_A = original_model.reaction(reaction_index).low_rate.pre_exponential_factor
    original_high_A = original_model.reaction(reaction_index).high_rate.pre_exponential_factor
    success = perturb_reaction_rate_coefficient(model_path, perturbed_model_path, reaction_index)
    assert success is True
    perturbed_model = ct.Solution(perturbed_model_path)
    perturbed_low_A = perturbed_model.reaction(reaction_index).low_rate.pre_exponential_factor
    perturbed_high_A = perturbed_model.reaction(reaction_index).high_rate.pre_exponential_factor
    assert almost_equal(perturbed_low_A, original_low_A * (1 + DELTA_K))
    assert almost_equal(perturbed_high_A, original_high_A * (1 + DELTA_K))

    # Test Pressure-Dependent Arrhenius
    reaction_index = 5
    original_rates = [rc[1].pre_exponential_factor for rc in original_model.reaction(reaction_index).rates]
    success = perturb_reaction_rate_coefficient(model_path, perturbed_model_path, reaction_index)
    assert success is True
    perturbed_model = ct.Solution(perturbed_model_path)
    perturbed_rates = [rc[1].pre_exponential_factor for rc in perturbed_model.reaction(reaction_index).rates]
    for original_A, perturbed_A in zip(original_rates, perturbed_rates):
        assert almost_equal(perturbed_A, original_A * (1 + DELTA_K))

    os.remove(perturbed_model_path) if os.path.exists(perturbed_model_path) else None


def test_get_h298():
    """Test the ``get_h298()`` function."""
    model = ct.Solution(os.path.join(TEST_DATA_BASE_PATH, 'models', 'N2H4.yaml'))
    assert almost_equal(get_h298(model, 2), 97.5703)  # N2H4
    assert almost_equal(get_h298(model, 3), 0.0)  # N2
    assert almost_equal(get_h298(model, 4), 217.98615)  # H


def test_get_Ea_units():
    """
    Test the ``get_Ea_units()`` method.
    """
    assert get_Ea_units(os.path.join(TEST_DATA_BASE_PATH, 'models', 'N2H4.yaml')) == 'kcal/mol'
    assert get_Ea_units(os.path.join(TEST_DATA_BASE_PATH, 'models', 'eA_units.yaml')) == 'J/mol'


def test_get_pressure_from_cantera():
    """
    Test the ``get_pressure_from_cantera()`` method.
    """
    assert get_pressure_from_cantera('1 bar') == 1.0
    assert get_pressure_from_cantera('3 atm') == 3 * 1.01325
    assert get_pressure_from_cantera('100 Pa') == 1e-3


def test_calculate_arrhenius_rate_coefficient():
    """
    Test the ``calculate_arrhenius_rate()`` function.
    """
    assert almost_equal(calculate_arrhenius_rate_coefficient(
        A=1e12, n=0.5, Ea=30, T=1000, Ea_units='cal/mol') / 3.11e+13, 1.0, places=2)
    assert almost_equal(calculate_arrhenius_rate_coefficient(
        A=1e12, n=0.5, Ea=30, T=1000, Ea_units='kcal/mol') / 8.78e+06, 1.0, places=2)
    assert almost_equal(calculate_arrhenius_rate_coefficient(
        A=1e12, n=0.5, Ea=8.2, T=1000, Ea_units='kJ/mol') / 1.18e+13, 1.0, places=2)
    assert almost_equal(calculate_arrhenius_rate_coefficient(
        A=1e12, n=0.5, Ea=8200, T=1000, Ea_units='J/mol') / 1.18e+13, 1.0, places=2)


def test_calculate_troe_rate_coefficient():
    """
    Test the ``calculate_troe_rate()`` function.
    """
    reaction_data = {'equation': 'NH2(5) + NH2(5) (+M) <=> H4N2(1) (+M)',
                     'type': 'falloff',
                     'low-P-rate-constant': {'A': 1.6e+34, 'b': -5.49, 'Ea': 1.987},
                     'high-P-rate-constant': {'A': 5.6e+14, 'b': -0.414, 'Ea': 0.066},
                     'Troe': {'A': 0.31, 'T3': 1.0e-30, 'T1': 1.0e+30},
                     'efficiencies': {'N2(2)': 1.0, 'Ar': 0.5, 'ammonia(9)': 2.93}}
    assert almost_equal(calculate_troe_rate_coefficient(reaction_data=reaction_data, T=1000, P=0.1, Ea_units='cal/mol') / 3.69e+11, 1.0, places=2)
    assert almost_equal(calculate_troe_rate_coefficient(reaction_data=reaction_data, T=1000, P=100, Ea_units='cal/mol') / 1.72e+13, 1.0, places=2)
    reaction_data = {'equation': 'H4N2(1) (+M) <=> H(3) + H3N2(6) (+M)',
                     'type': 'falloff',
                     'low-P-rate-constant': {'A': 1.95e+47, 'b': -8.5, 'Ea': 82.384},
                     'high-P-rate-constant': {'A': 5.69e+14, 'b': -0.28, 'Ea': 81.034}}
    assert almost_equal(calculate_troe_rate_coefficient(reaction_data=reaction_data, T=1000, P=0.1, Ea_units='cal/mol') / 7.81e+13, 1.0, places=2)
    assert almost_equal(calculate_troe_rate_coefficient(reaction_data=reaction_data, T=1000, P=100, Ea_units='cal/mol') / 7.90e+13, 1.0, places=2)


def test_calculate_plog_rate_coefficient():
    """
    Test the ``calculate_plog_rate()`` function.
    """
    reaction_data = {'equation': 'NH2(5) + NH2(5) <=> H(3) + H3N2(6)',
                     'type': 'pressure-dependent-Arrhenius',
                     'rate-constants': [{'P': '0.1 atm', 'A': 9.2e+11, 'b': -0.01, 'Ea': 10.014},
                                        {'P': '1.0 atm', 'A': 1.2e+12, 'b': -0.03, 'Ea': 10.084},
                                        {'P': '10.0 atm', 'A': 4.7e+12, 'b': -0.2, 'Ea': 10.62}]}
    assert almost_equal(calculate_plog_rate_coefficient(reaction_data=reaction_data, T=1000, P=0.5, Ea_units='cal/mol') / 9.33e+11, 1.0, places=2)
    assert almost_equal(calculate_plog_rate_coefficient(reaction_data=reaction_data, T=1000, P=7.5, Ea_units='cal/mol') / 1.15e+12, 1.0, places=2)


def test_calculate_chebyshev_rate_coefficient():
    """
    Test the ``calculate_chebyshev_rate()`` function.
    """
    reaction_data = {'equation': 'N2H3 <=> N2H2 + H',
                     'type': 'Chebyshev',
                     'temperature-range': [300.0, 3000.0],
                     'pressure-range': ['0.01 bar', '100 bar'],
                     'data': [[-6.40114, 0.924384, -0.148344, 0.00833621],
                              [14.5099, 0.809594, 0.0374586, -0.0271928],
                              [-0.620903, 0.192692, 0.0680562, 0.00277832],
                              [-0.312241, 0.0135261, 0.0200691, 0.00804536],
                              [-0.133503, -0.0148586, -0.000903153, 0.00233355],
                              [-0.0440289, -0.012249, -0.00402369, -0.000514046]]}
    assert almost_equal(calculate_chebyshev_rate_coefficient(reaction_data=reaction_data, T=1000, P=100) / 2.48e+03, 1.0, places=2)
    assert almost_equal(calculate_chebyshev_rate_coefficient(reaction_data=reaction_data, T=1500, P=0.5) / 8.19e+04, 1.0, places=2)
    assert almost_equal(calculate_chebyshev_rate_coefficient(reaction_data=reaction_data, T=1000, P=20) / 1.49e+03, 1.0, places=2)
    assert almost_equal(calculate_chebyshev_rate_coefficient(reaction_data=reaction_data, T=2000, P=10) / 1.26e+07, 1.0, places=2)


def test_get_top_sa_coefficients():
    """Test the get_top_sa_coefficients() function"""
    data_path = os.path.join(TEST_DATA_BASE_PATH, 'sa_idt_1.yaml')
    top_sa_dict = get_top_sa_coefficients(idt_sa_dict=read_yaml_file(data_path),
                                          top_species=10,
                                          top_reactions=10)
    assert top_sa_dict['kinetics']['IDT'][0][1.0][10.0][1000.0] == {13: 0.3205566319720274,
                                                                    14: -0.22644432990213972,
                                                                    45: -0.2211180397468687,
                                                                    46: -0.24633871316234654,
                                                                    47: 0.2964208373351383,
                                                                    51: -0.11147805613602418,
                                                                    67: 0.22325779078467647,
                                                                    68: -0.31029061953466325,
                                                                    77: 0.15626917713911037,
                                                                    110: -0.2865901229666261}
    assert top_sa_dict['kinetics']['IDT'][0][1.0][10.0][1280.0] == {14: 0.15935729904061308,
                                                                    44: -0.12459281014675062,
                                                                    47: -0.6053434095960982,
                                                                    51: -0.44917550899211955,
                                                                    53: -0.24370751128422777,
                                                                    64: -0.21607974243243627,
                                                                    67: -0.20653334574996596,
                                                                    75: -0.21197372330142064,
                                                                    77: -0.20221059540821848,
                                                                    81: -0.21607829777848703}
    assert top_sa_dict['thermo']['IDT'][0][1.0][10.0][1600.0] == {3: 0.021847208262576294,
                                                                  4: -0.01674391971202145,
                                                                  5: 0.00706845583932834,
                                                                  6: -0.007669745192204748,
                                                                  9: 0.002504932522953768,
                                                                  11: 0.01786371289736624,
                                                                  18: -0.05249629811869966,
                                                                  19: 0.03296233019706528,
                                                                  20: 0.0022711498018319824,
                                                                  21: -0.005124399631177829}
    assert top_sa_dict['thermo']['IDT'][0][1.0][10.0][2000.0] == {3: 0.0056417741975125545,
                                                                  4: -0.006821261336671581,
                                                                  5: 0.01210001530866201,
                                                                  6: -0.003575609137733701,
                                                                  10: 0.0012149600916090855,
                                                                  11: 0.003058814430615631,
                                                                  18: -0.027582072400766624,
                                                                  19: 0.01896165943937473,
                                                                  24: -0.001310798627769729,
                                                                  27: -0.0008972909873407294}

    data_path = os.path.join(TEST_DATA_BASE_PATH, 'sa_idt_2.yaml')
    top_sa_dict = get_top_sa_coefficients(idt_sa_dict=read_yaml_file(data_path),
                                          top_species=10,
                                          top_reactions=10)
    assert top_sa_dict['kinetics']['IDT'][0][1.0][10.0][1200.0] == {0: -0.1753106434982762,
                                                                    15: 0.12127708624343236,
                                                                    62: -0.1253435501225756,
                                                                    65: -0.1845257082104551,
                                                                    66: -0.1334340281883038,
                                                                    85: -0.10664923956319464,
                                                                    90: 0.12718281741288295,
                                                                    108: -0.15507533285747638,
                                                                    165: -0.11471694565718685,
                                                                    245: -0.11446478011585266}
    assert top_sa_dict['thermo']['IDT'][0][1.0][10.0][1200.0] == {2: -0.005988036940960675,
                                                                  3: -0.09572951383236826,
                                                                  4: 0.0210129385535655,
                                                                  5: -0.04500208755278678,
                                                                  9: 0.08203778520740615,
                                                                  10: 0.005072337466243485,
                                                                  12: -0.01147549061238161,
                                                                  24: 0.011809516976394461,
                                                                  34: 0.005113408917627866,
                                                                  35: 0.005070908258511206}

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
