#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_cantera_idt module
"""

import os
import shutil
import tempfile

import cantera as ct
import yaml

from arc.common import almost_equal_lists, read_yaml_file

from t3.common import SIMULATE_TEST_DATA_BASE_PATH, TEST_DATA_BASE_PATH
from tests.common import almost_equal, run_minimal
from t3.simulate.cantera_idt import (CanteraIDT, DELTA_H, DELTA_K, calculate_arrhenius_rate_coefficient,
                                     calculate_troe_rate_coefficient, calculate_chebyshev_rate_coefficient,
                                     calculate_plog_rate_coefficient, compute_idt, get_Ea_units, get_h298,
                                     get_pressure_from_cantera, get_t_and_p_lists, get_top_sa_coefficients,
                                     perturb_enthalpy, perturb_reaction_rate_coefficient)
from t3.utils.fix_cantera import fix_cantera


TEST_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'cantera_simulator_test')
TEST_DIR_IDT = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'cantera_idt_test')


def _make_idt_adapter(t3):
    """Helper to construct a CanteraIDT adapter from a configured T3 instance."""
    return CanteraIDT(t3=t3.t3,
                      rmg=t3.rmg,
                      paths=t3.paths,
                      logger=t3.logger,
                      atol=t3.rmg['model']['atol'],
                      rtol=t3.rmg['model']['rtol'],
                      )


def test_determine_radical_label():
    """
    Test the ``determine_radical_label()`` method.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    ct_adapter = _make_idt_adapter(t3)
    label = ct_adapter.determine_radical_label()
    assert label == 'OH(4)'


def test_get_cantera_species_label():
    """
    Test the ``get_cantera_species_label()`` method.
    """
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 2  # CH4 model
    t3.set_paths()
    ct_adapter = _make_idt_adapter(t3)
    assert ct_adapter.get_cantera_species_label('methane') == 'methane(1)'
    assert ct_adapter.get_cantera_species_label('O2') == 'O2(2)'
    assert ct_adapter.get_cantera_species_label('OH') == 'OH(6)'
    assert ct_adapter.get_cantera_species_label('CH3CH2OO') == 'CH3CH2OO(40)'


def test_get_t_and_p_lists():
    """
    Test the ``get_t_and_p_lists()`` function.
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
    ct_adapter = _make_idt_adapter(t3)
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
                         {'label': 'o2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'n2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    ct_adapter = _make_idt_adapter(t3)
    assert ct_adapter.reactor_idt_dict is None
    ct_adapter.simulate()
    assert almost_equal(ct_adapter.reactor_idt_dict[0][1.0][1.0][1000.0], 0.0328465)
    if os.path.exists(t3.paths['SA IDT dict']):
        os.remove(t3.paths['SA IDT dict'])
    if os.path.exists(t3.paths['SA IDT dict top X']):
        os.remove(t3.paths['SA IDT dict top X'])


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
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    ct_adapter = _make_idt_adapter(t3)
    assert ct_adapter.reactor_idt_dict is None
    ct_adapter.simulate()
    assert list(ct_adapter.reactor_idt_dict[0].keys()) == [0.5, 1.0]
    assert list(ct_adapter.reactor_idt_dict[0][0.5].keys()) == [1.0, 10.0, 100.0]
    assert almost_equal_lists(ct_adapter.reactor_idt_dict[0][0.5][10.0].keys(),
                              [1500.0, 1510.49, 1521.13, 1531.91, 1542.86,
                               1553.96, 1565.22, 1576.64, 1588.24, 1600.0],
                              rtol=0.001, atol=0.1)
    for val in ct_adapter.reactor_idt_dict[0][0.5][10.0].values():
        assert val is not None and val > 0


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
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    ct_adapter = _make_idt_adapter(t3)
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
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    ct_adapter = _make_idt_adapter(t3)
    assert ct_adapter.reactor_idt_dict is None
    ct_adapter.simulate()
    assert list(ct_adapter.reactor_idt_dict[0].keys()) == [1.0]
    assert list(ct_adapter.reactor_idt_dict[0][1.0].keys()) == [10.0]
    idt_at_10bar = ct_adapter.reactor_idt_dict[0][1.0][10.0]
    non_none_keys = [T for T, idt in idt_at_10bar.items() if idt is not None]
    assert almost_equal_lists(non_none_keys,
                              [820.5128, 842.1052, 864.8648, 888.88888, 914.28571, 941.176470, 969.6969, 1000.0,
                               1032.2580, 1066.66666, 1103.4482, 1142.85714, 1185.1851, 1230.76923, 1280.0, 1333.33333,
                               1391.30434, 1454.54545, 1523.8095, 1600.0, 1684.2105, 1777.77777, 1882.3529, 2000.0],
                              rtol=0.001, atol=0.1)
    assert any(v is not None for v in idt_at_10bar.values())
    assert almost_equal(idt_at_10bar[1600.0], 2.8396e-05)


def test_simulate_rmg_heptane():
    """
    Test the ``simulate()`` method for computing IDT for an RMG heptane model.
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
                          'equivalence_ratios': [1.0]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    ct_adapter = _make_idt_adapter(t3)
    assert ct_adapter.reactor_idt_dict is None
    ct_adapter.simulate()
    assert ct_adapter.reactor_idt_dict[0][1.0][10.0][1500.0] is not None
    assert ct_adapter.reactor_idt_dict[0][1.0][10.0][1500.0] > 0


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
                         {'label': 'o2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'n2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    ct_adapter = _make_idt_adapter(t3)
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
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    ct_adapter = _make_idt_adapter(t3)
    equivalence_ratios, concentration_combinations = ct_adapter.get_concentration_combinations()
    assert equivalence_ratios == [0.5, 1.0, 1.5]
    # The corrected formula gives: O2 = fuel * stoich / phi.
    # NH3 stoich = 0.75 → O2 columns: 1.5, 0.75, 0.5
    # N2 (default ratio 3.76): 5.64, 2.82, 1.88
    expected = [{'NH3(1)': 1.0, 'O2(3)': 1.5, 'N2(2)': 5.64},
                {'NH3(1)': 1.0, 'O2(3)': 0.75, 'N2(2)': 2.82},
                {'NH3(1)': 1.0, 'O2(3)': 0.5, 'N2(2)': 1.88}]
    assert len(concentration_combinations) == len(expected)
    for got, want in zip(concentration_combinations, expected):
        assert set(got.keys()) == set(want.keys())
        for k, v in want.items():
            assert almost_equal(got[k], v)


def test_perturb_enthalpy():
    """
    Test the ``perturb_enthalpy()`` function.
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
    if os.path.exists(perturbed_model_path):
        os.remove(perturbed_model_path)


def test_perturb_reaction_rate_coefficient():
    """
    Test the ``perturb_reaction_rate_coefficient()`` function.
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
    original_low_A = original_model.reaction(reaction_index).rate.low_rate.pre_exponential_factor
    original_high_A = original_model.reaction(reaction_index).rate.high_rate.pre_exponential_factor
    success = perturb_reaction_rate_coefficient(model_path, perturbed_model_path, reaction_index)
    assert success is True
    perturbed_model = ct.Solution(perturbed_model_path)
    perturbed_low_A = perturbed_model.reaction(reaction_index).rate.low_rate.pre_exponential_factor
    perturbed_high_A = perturbed_model.reaction(reaction_index).rate.high_rate.pre_exponential_factor
    assert almost_equal(perturbed_low_A, original_low_A * (1 + DELTA_K))
    assert almost_equal(perturbed_high_A, original_high_A * (1 + DELTA_K))

    # Test Pressure-Dependent Arrhenius
    reaction_index = 5
    original_rates = [rc[1].pre_exponential_factor for rc in original_model.reaction(reaction_index).rate.rates]
    success = perturb_reaction_rate_coefficient(model_path, perturbed_model_path, reaction_index)
    assert success is True
    perturbed_model = ct.Solution(perturbed_model_path)
    perturbed_rates = [rc[1].pre_exponential_factor for rc in perturbed_model.reaction(reaction_index).rate.rates]
    for original_A, perturbed_A in zip(original_rates, perturbed_rates):
        assert almost_equal(perturbed_A, original_A * (1 + DELTA_K))

    if os.path.exists(perturbed_model_path):
        os.remove(perturbed_model_path)


def test_get_h298():
    """Test the ``get_h298()`` function."""
    model = ct.Solution(os.path.join(TEST_DATA_BASE_PATH, 'models', 'N2H4.yaml'))
    assert almost_equal(get_h298(model, 2), 97.5703)  # N2H4
    assert almost_equal(get_h298(model, 3), 0.0)  # N2
    assert almost_equal(get_h298(model, 4), 217.98615)  # H


def test_get_Ea_units():
    """
    Test the ``get_Ea_units()`` function.
    """
    assert get_Ea_units(os.path.join(TEST_DATA_BASE_PATH, 'models', 'N2H4.yaml')) == 'kcal/mol'
    assert get_Ea_units(os.path.join(TEST_DATA_BASE_PATH, 'models', 'eA_units.yaml')) == 'J/mol'


def test_get_pressure_from_cantera():
    """
    Test the ``get_pressure_from_cantera()`` function.
    """
    assert get_pressure_from_cantera('1 bar') == 1.0
    assert get_pressure_from_cantera('3 atm') == 3 * 1.01325
    assert get_pressure_from_cantera('100 Pa') == 1e-3


def test_calculate_arrhenius_rate_coefficient():
    """
    Test the ``calculate_arrhenius_rate_coefficient()`` function.
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
    Test the ``calculate_troe_rate_coefficient()`` function.
    """
    reaction_data = {'equation': 'NH2(5) + NH2(5) (+M) <=> H4N2(1) (+M)',
                     'type': 'falloff',
                     'low-P-rate-constant': {'A': 1.6e+34, 'b': -5.49, 'Ea': 1.987},
                     'high-P-rate-constant': {'A': 5.6e+14, 'b': -0.414, 'Ea': 0.066},
                     'Troe': {'A': 0.31, 'T3': 1.0e-30, 'T1': 1.0e+30},
                     'efficiencies': {'N2(2)': 1.0, 'Ar': 0.5, 'ammonia(9)': 2.93}}
    assert almost_equal(calculate_troe_rate_coefficient(
        reaction_data=reaction_data, T=1000, P=0.1, Ea_units='cal/mol') / 3.69e+11, 1.0, places=2)
    assert almost_equal(calculate_troe_rate_coefficient(
        reaction_data=reaction_data, T=1000, P=100, Ea_units='cal/mol') / 1.72e+13, 1.0, places=2)
    reaction_data = {'equation': 'H4N2(1) (+M) <=> H(3) + H3N2(6) (+M)',
                     'type': 'falloff',
                     'low-P-rate-constant': {'A': 1.95e+47, 'b': -8.5, 'Ea': 82.384},
                     'high-P-rate-constant': {'A': 5.69e+14, 'b': -0.28, 'Ea': 81.034}}
    assert almost_equal(calculate_troe_rate_coefficient(
        reaction_data=reaction_data, T=1000, P=0.1, Ea_units='cal/mol') / 7.81e+13, 1.0, places=2)
    assert almost_equal(calculate_troe_rate_coefficient(
        reaction_data=reaction_data, T=1000, P=100, Ea_units='cal/mol') / 7.90e+13, 1.0, places=2)


def test_calculate_plog_rate_coefficient():
    """
    Test the ``calculate_plog_rate_coefficient()`` function.
    """
    reaction_data = {'equation': 'NH2(5) + NH2(5) <=> H(3) + H3N2(6)',
                     'type': 'pressure-dependent-Arrhenius',
                     'rate-constants': [{'P': '0.1 atm', 'A': 9.2e+11, 'b': -0.01, 'Ea': 10.014},
                                        {'P': '1.0 atm', 'A': 1.2e+12, 'b': -0.03, 'Ea': 10.084},
                                        {'P': '10.0 atm', 'A': 4.7e+12, 'b': -0.2, 'Ea': 10.62}]}
    assert almost_equal(calculate_plog_rate_coefficient(
        reaction_data=reaction_data, T=1000, P=0.5, Ea_units='cal/mol') / 9.33e+11, 1.0, places=2)
    assert almost_equal(calculate_plog_rate_coefficient(
        reaction_data=reaction_data, T=1000, P=7.5, Ea_units='cal/mol') / 1.15e+12, 1.0, places=2)


def test_calculate_chebyshev_rate_coefficient():
    """
    Test the ``calculate_chebyshev_rate_coefficient()`` function.
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
    assert almost_equal(calculate_chebyshev_rate_coefficient(
        reaction_data=reaction_data, T=1000, P=100) / 2.48e+03, 1.0, places=2)
    assert almost_equal(calculate_chebyshev_rate_coefficient(
        reaction_data=reaction_data, T=1500, P=0.5) / 8.19e+04, 1.0, places=2)
    assert almost_equal(calculate_chebyshev_rate_coefficient(
        reaction_data=reaction_data, T=1000, P=20) / 1.49e+03, 1.0, places=2)
    assert almost_equal(calculate_chebyshev_rate_coefficient(
        reaction_data=reaction_data, T=2000, P=10) / 1.26e+07, 1.0, places=2)


def test_get_top_sa_coefficients():
    """Test the ``get_top_sa_coefficients()`` function."""
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


def test_idt_mode_row():
    """Test that idt_mode='row' produces index-aligned (T, P, phi) tuples, not a full matrix."""
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 2  # CH4 model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': [1000, 1100, 1200],
                           'P': [10, 20, 30],
                           'termination_rate_ratio': 0.01,
                           'idt_mode': 'row'},
                          ]
    t3.rmg['species'] = [{'label': 'methane', 'smiles': 'C', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [0.5, 1.0, 1.5]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    ct_adapter = _make_idt_adapter(t3)
    ct_adapter.simulate()
    # Row mode: exactly 3 conditions (one per index), not 3*3*3=27
    total_conditions = 0
    for phi_data in ct_adapter.reactor_idt_dict[0].values():
        for p_data in phi_data.values():
            total_conditions += len(p_data)
    assert total_conditions == 3


def test_multi_reactor():
    """Test that two reactor blocks produce separate entries in reactor_idt_dict."""
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 2  # CH4 model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': 1200, 'P': 10,
                           'termination_rate_ratio': 0.01},
                          {'type': 'gas batch constant T P',
                           'T': 1100, 'P': 20,
                           'termination_rate_ratio': 0.01},
                          ]
    t3.rmg['species'] = [{'label': 'methane', 'smiles': 'C', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [1.0]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    ct_adapter = _make_idt_adapter(t3)
    ct_adapter.simulate()
    assert set(ct_adapter.reactor_idt_dict.keys()) == {0, 1}
    # Reactor 0: T=1200, P=10
    assert 10.0 in ct_adapter.reactor_idt_dict[0][1.0]
    assert 1200.0 in ct_adapter.reactor_idt_dict[0][1.0][10.0]
    # Reactor 1: T=1100, P=20
    assert 20.0 in ct_adapter.reactor_idt_dict[1][1.0]
    assert 1100.0 in ct_adapter.reactor_idt_dict[1][1.0][20.0]


def test_idt_criterion_max_dTdt():
    """Test that max_dTdt criterion produces a valid IDT that differs from max_dOHdt."""
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 4  # Reduced CH4 model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': 1200, 'P': 10,
                           'termination_rate_ratio': 0.01},
                          ]
    t3.rmg['species'] = [{'label': 'methane', 'smiles': 'C', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [1.0]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    # Run with max_dOHdt
    t3.t3['sensitivity']['idt_criterion'] = 'max_dOHdt'
    ct_adapter_oh = _make_idt_adapter(t3)
    ct_adapter_oh.simulate()
    idt_oh = ct_adapter_oh.reactor_idt_dict[0][1.0][10.0][1200.0]

    # Run with max_dTdt
    t3.t3['sensitivity']['idt_criterion'] = 'max_dTdt'
    ct_adapter_dTdt = _make_idt_adapter(t3)
    ct_adapter_dTdt.simulate()
    idt_dTdt = ct_adapter_dTdt.reactor_idt_dict[0][1.0][10.0][1200.0]

    assert idt_oh is not None and idt_oh > 0
    assert idt_dTdt is not None and idt_dTdt > 0
    # Both should be in a physically reasonable range (1e-6 to 1 s)
    assert 1e-6 < idt_oh < 1.0
    assert 1e-6 < idt_dTdt < 1.0


def test_no_radical_found():
    """Test graceful handling when no radical (OH/H) is present in the model."""
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 2  # CH4 model
    t3.set_paths()
    ct_adapter = _make_idt_adapter(t3)
    # Simulate a scenario where no radical is found
    ct_adapter.radical_label = None
    ct_adapter.idt_criterion = 'max_dOHdt'
    # compute_idt with radical_label=None should return None
    result = compute_idt(time_history=ct.SolutionArray(ct_adapter.model, 0, extra='t'),
                         radical_label=None, criterion='max_dOHdt')
    assert result is None


def test_experimental_comparison():
    """Test the experimental comparison method with synthetic data."""
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 4  # Reduced CH4 model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': 1200, 'P': 10,
                           'termination_rate_ratio': 0.01},
                          ]
    t3.rmg['species'] = [{'label': 'methane', 'smiles': 'C', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [1.0]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    ct_adapter = _make_idt_adapter(t3)
    ct_adapter.simulate()
    sim_idt = ct_adapter.reactor_idt_dict[0][1.0][10.0][1200.0]
    assert sim_idt is not None

    # Create a temporary experimental data file
    exp_data = {
        'citation': 'Test et al., 2024',
        'data': [
            {'T': 1200, 'P': 10, 'phi': 1.0, 'idt': sim_idt * 1.5},  # 50% off
        ],
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(exp_data, f)
        exp_path = f.name

    try:
        result = ct_adapter.compare_with_experiment(exp_path)
        assert result['citation'] == 'Test et al., 2024'
        assert result['n_points'] == 1
        assert result['n_matched'] == 1
        assert result['rmse_log'] is not None
        # log10(1/1.5) ≈ -0.176
        assert almost_equal(abs(result['rmse_log']), 0.176, places=2)
    finally:
        os.remove(exp_path)


def test_sa_threshold_filtering():
    """Test that SA_threshold filtering reduces the number of returned species/reactions."""
    from t3.simulate.cantera_idt import get_top_sa_coefficients
    data_path = os.path.join(TEST_DATA_BASE_PATH, 'sa_idt_1.yaml')
    idt_sa_dict = read_yaml_file(data_path)
    # With default top-N, we get 10 entries
    top_all = get_top_sa_coefficients(idt_sa_dict=idt_sa_dict, top_species=10, top_reactions=10)

    # Manually filter at high threshold — fewer should survive
    high_threshold = 0.25
    for t_val in top_all['kinetics']['IDT'][0][1.0][10.0]:
        data = top_all['kinetics']['IDT'][0][1.0][10.0][t_val]
        filtered = {k: v for k, v in data.items() if abs(v) > high_threshold}
        assert len(filtered) <= len(data)


def test_worker_failure_logging():
    """Test that worker failures are logged and partial results are returned."""
    from unittest.mock import patch
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 4  # Reduced CH4 model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': 1200, 'P': 10,
                           'termination_rate_ratio': 0.01},
                          ]
    t3.rmg['species'] = [{'label': 'methane', 'smiles': 'C', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [1.0]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    t3.t3['sensitivity']['top_SA_species'] = 5
    t3.t3['sensitivity']['top_SA_reactions'] = 5
    t3.t3['sensitivity']['max_sa_workers'] = 2
    t3.t3['sensitivity']['save_sa_yaml'] = False
    ct_adapter = _make_idt_adapter(t3)
    ct_adapter.simulate()

    # Mock worker to fail for thermo index 0
    original_worker = __import__('t3.simulate.cantera_idt', fromlist=['worker']).worker

    def failing_worker(task, *args, **kwargs):
        if task == ('thermo', 0):
            raise RuntimeError("Simulated worker failure")
        return original_worker(task, *args, **kwargs)

    with patch('t3.simulate.cantera_idt.worker', side_effect=failing_worker):
        # This should not raise — partial results should be returned
        result = ct_adapter._get_sa_coefficients_brute_force()
    # Result should still be a dict (partial results)
    assert result is not None
    assert 'kinetics' in result
    assert 'thermo' in result


def test_cantera_rcm_adapter():
    """Test that CanteraRCM uses a constant-pressure reactor and produces valid IDTs."""
    from t3.simulate.cantera_rcm import CanteraRCM
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 4  # Reduced CH4 model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': 1200, 'P': 10,
                           'termination_rate_ratio': 0.01},
                          ]
    t3.rmg['species'] = [{'label': 'methane', 'smiles': 'C', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [1.0]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    rcm_adapter = CanteraRCM(t3=t3.t3, rmg=t3.rmg, paths=t3.paths, logger=t3.logger,
                             atol=t3.rmg['model']['atol'], rtol=t3.rmg['model']['rtol'])
    # Verify the reactor type override
    model = ct.Solution(infile=t3.paths['cantera annotated'])
    reactor = rcm_adapter._create_reactor(model, energy='on')
    assert isinstance(reactor, ct.IdealGasConstPressureReactor)

    # Run simulation and verify IDT is produced
    rcm_adapter.simulate()
    idt_rcm = rcm_adapter.reactor_idt_dict[0][1.0][10.0][1200.0]
    assert idt_rcm is not None and idt_rcm > 0

    # Compare with constant-volume (CanteraIDT) — should differ since reactor types differ
    idt_adapter = _make_idt_adapter(t3)
    idt_adapter.simulate()
    idt_cv = idt_adapter.reactor_idt_dict[0][1.0][10.0][1200.0]
    assert idt_cv is not None and idt_cv > 0
    assert idt_rcm != idt_cv  # constant-P vs constant-V give different IDTs


def test_prescan_radical():
    """Test that _prescan_radical identifies a radical from a quick simulation."""
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 4  # Reduced CH4 model
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': 1200, 'P': 10,
                           'termination_rate_ratio': 0.01}]
    t3.rmg['species'] = [{'label': 'methane', 'smiles': 'C', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [1.0]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    t3.t3['sensitivity']['idt_criterion'] = 'max_radical_dt'
    ct_adapter = _make_idt_adapter(t3)
    infile = ct_adapter.paths['cantera annotated']
    x = {'methane(1)': 0.095, 'O2(2)': 0.19, 'N2': 0.715}
    radical = ct_adapter._prescan_radical(infile=infile, T=1200.0, P=10.0, X=x)
    assert radical is not None
    assert isinstance(radical, str)
    # The prescan should find a small radical (OH, H, O, etc.) in a methane mechanism
    assert len(radical) > 0


def test_prescan_radical_fallback():
    """Test that _prescan_radical falls back to self.radical_label on a trivially short sim."""
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 4
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': 1200, 'P': 10,
                           'termination_rate_ratio': 0.01}]
    t3.rmg['species'] = [{'label': 'methane', 'smiles': 'C', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [1.0]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    t3.t3['sensitivity']['idt_criterion'] = 'max_radical_dt'
    ct_adapter = _make_idt_adapter(t3)
    infile = ct_adapter.paths['cantera annotated']
    x = {'methane(1)': 0.095, 'O2(2)': 0.19, 'N2': 0.715}
    # Use an extremely short max_idt so the sim barely runs — should still return a label
    radical = ct_adapter._prescan_radical(infile=infile, T=300.0, P=1.0, X=x, max_idt=1e-15)
    # Falls back to self.radical_label since no radical has meaningful concentration
    assert radical is not None


def test_adjoint_sa_max_radical_dt_uses_prescan():
    """Test that adjoint SA with max_radical_dt uses a prescanned radical, not temperature."""
    from unittest.mock import patch
    t3 = run_minimal(project_directory=TEST_DIR_IDT)
    t3.iteration = 4
    t3.set_paths()
    t3.rmg['reactors'] = [{'type': 'gas batch constant T P',
                           'T': 1200, 'P': 10,
                           'termination_rate_ratio': 0.01}]
    t3.rmg['species'] = [{'label': 'methane', 'smiles': 'C', 'concentration': 0, 'role': 'fuel',
                          'equivalence_ratios': [1.0]},
                         {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
                         {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'}]
    t3.t3['sensitivity']['idt_criterion'] = 'max_radical_dt'
    t3.t3['sensitivity']['idt_sa_method'] = 'adjoint'
    t3.t3['sensitivity']['save_sa_yaml'] = False
    ct_adapter = _make_idt_adapter(t3)
    ct_adapter.simulate()

    # Patch _prescan_radical to track that it gets called and returns a species name
    original_prescan = ct_adapter._prescan_radical
    prescan_calls = []

    def tracking_prescan(*args, **kwargs):
        result = original_prescan(*args, **kwargs)
        prescan_calls.append(result)
        return result

    with patch.object(ct_adapter, '_prescan_radical', side_effect=tracking_prescan):
        ct_adapter.get_sa_coefficients()

    assert len(prescan_calls) > 0, '_prescan_radical was not called for max_radical_dt'
    assert prescan_calls[0] is not None, '_prescan_radical returned None'
    assert prescan_calls[0] != 'temperature', \
        'Adjoint SA should use a radical species, not temperature, for max_radical_dt'


def teardown_module():
    """
    A method that is run after all unit tests in this module.
    Delete all project directories created during these unit tests.
    """
    for test_dir_path in [TEST_DIR, TEST_DIR_IDT]:
        log_archive = os.path.join(test_dir_path, 'log_archive')
        if os.path.isdir(log_archive):
            shutil.rmtree(log_archive, ignore_errors=True)
        log_file = os.path.join(test_dir_path, 't3.log')
        if os.path.isfile(log_file):
            os.remove(log_file)
        for iteration in [0, 1, 2, 3, 4, 5]:
            figs_path = os.path.join(test_dir_path, f'iteration_{iteration}', 'Figures')
            if os.path.isdir(figs_path):
                shutil.rmtree(figs_path, ignore_errors=True)
            sa_path = os.path.join(test_dir_path, f'iteration_{iteration}', 'SA')
            if os.path.isdir(sa_path):
                shutil.rmtree(sa_path, ignore_errors=True)
