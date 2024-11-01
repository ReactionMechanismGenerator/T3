#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_utils module
"""

import os
import pytest

from rmgpy.species import Species

import t3.common as common
from t3.common import TEST_DATA_BASE_PATH
from t3.schema import RMGSpecies
from tests.common import run_minimal


def test_dict_to_str():
    """Test prettifying a dictionary"""
    dictionary = {'label1': {'spc': 'Species1', 'reason': 'Reason1'},
                  'label2': {'spc': 'Species2', 'reason': 'Reason2'}}
    output = common.dict_to_str(dictionary)
    expected_output = """label1:
  spc: Species1
  reason: Reason1
label2:
  spc: Species2
  reason: Reason2
"""
    assert output == expected_output


def test_get_species_by_label():
    """Test getting species by label"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    rmg_species, rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    label = 'H2O'
    species = common.get_species_by_label(label, rmg_species)
    assert species.label == label
    assert species.index == 7


def test_get_rmg_species_from_a_species_dict():
    """Test getting an RMG species from a species dictionary"""
    species = common.get_rmg_species_from_a_species_dict(
        species_dict=RMGSpecies(**{'label': 'spc', 'smiles': 'C=O'}).dict())
    assert isinstance(species, Species)
    assert species.label == 'spc'
    assert species.molecule[0].to_smiles() == 'C=O'

    adj = """1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 O u0 p2 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}"""
    species = common.get_rmg_species_from_a_species_dict(
        species_dict=RMGSpecies(**{'label': 'spc', 'adjlist': adj}).dict())
    assert isinstance(species, Species)
    assert species.label == 'spc'
    assert species.molecule[0].to_smiles() == 'C=O'

    species = common.get_rmg_species_from_a_species_dict(
        species_dict=RMGSpecies(**{'label': 'spc', 'inchi': 'InChI=1S/CH2O/c1-2/h1H2'}).dict())
    assert isinstance(species, Species)
    assert species.label == 'spc'
    assert species.molecule[0].to_smiles() == 'C=O'

    xyz = """O  0.0000000  0.0000000  0.7047750
C  0.0000000  0.0000000 -0.5593030
H  0.0000000  0.9470590 -1.1411940
H  0.0000000 -0.9470590 -1.1411940"""
    species = common.get_rmg_species_from_a_species_dict(
        species_dict=RMGSpecies(**{'label': 'spc', 'xyz': [xyz]}).dict())
    assert isinstance(species, Species)
    assert species.label == 'spc'
    assert species.molecule[0].to_smiles() == 'C=O'

    species = common.get_rmg_species_from_a_species_dict(species_dict=RMGSpecies(**{'label': 'spc'}).dict(),
                                                         raise_error=False)
    assert species is None
    with pytest.raises(ValueError):
        common.get_rmg_species_from_a_species_dict(species_dict=RMGSpecies(**{'label': 'spc'}).dict(),
                                                   raise_error=True)


def test_convert_termination_time_to_seconds():
    """Test converting the termination_time from an RMG reactor to seconds"""
    t_final = (1, 'micro-s')
    assert common.convert_termination_time_to_seconds(t_final) == 1e-6

    t_final = (1, 'ms')
    assert common.convert_termination_time_to_seconds(t_final) == 1e-3

    t_final = (1, 's')
    assert common.convert_termination_time_to_seconds(t_final) == 1

    t_final = (1, 'hours')
    assert common.convert_termination_time_to_seconds(t_final) == 3600

    t_final = (1, 'hrs')
    assert common.convert_termination_time_to_seconds(t_final) == 3600

    t_final = (1, 'days')
    assert common.convert_termination_time_to_seconds(t_final) == 3600*24


def test_get_values_within_range():
    """Test the get_values_within_range() function"""
    assert common.get_values_within_range([0, 10], 1) == [5]
    assert common.get_values_within_range([1, 10], 2) == [4, 7]
    assert common.get_values_within_range([0, 10], 3) == [0, 5, 10]
    assert common.get_values_within_range([1, 10], 4) == [1, 4, 7, 10]
    assert common.get_values_within_range(7.2, 40) == [7.2]
    assert common.get_values_within_range([7.2], 40) == [7.2]
    assert common.get_values_within_range([1, 25], 3, use_log_scale=True) == [1, 10]
    assert common.get_values_within_range([1, 500], 3, use_log_scale=True) == [1, 10, 100]
    assert common.get_values_within_range([0.1, 1000], 5, use_log_scale=True) == [0.1, 1, 10, 100, 1000]


def test_get_interval():
    """Test thew get_interval() function"""
    assert common.get_interval([0, 10], 1) == 5
    assert common.get_interval([1, 10], 2) == 3
    assert common.get_interval([0, 10], 3) == 5
    assert common.get_interval([0, 9], 4) == 3
    assert common.get_interval([500, 1200], 3) == 350


def test_get_chem_to_rmg_rxn_index_map():
    """Test the get_chem_to_rmg_rxn_index_map() function"""
    chemkin_path = os.path.join(TEST_DATA_BASE_PATH, 'chem_annotated_1.inp')
    rxn_map = common.get_chem_to_rmg_rxn_index_map(chem_annotated_path=chemkin_path)
    assert rxn_map == {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9, 11: 9, 12: 10, 13: 11, 14: 12,
                       15: 13, 16: 14, 17: 15, 18: 16, 19: 17, 20: 18, 21: 19, 22: 20, 23: 21, 24: 22, 25: 23, 26: 24,
                       27: 25, 28: 26, 29: 27, 30: 28, 31: 29, 32: 30, 33: 30, 34: 31, 35: 32, 36: 33, 37: 34, 38: 35,
                       39: 36, 40: 37, 41: 38, 42: 39, 43: 40, 44: 41, 45: 41, 46: 42, 47: 43, 48: 44, 49: 45, 50: 46,
                       51: 47, 52: 48, 53: 49, 54: 50, 55: 51, 56: 52, 57: 53, 58: 54, 59: 55, 60: 56, 61: 57, 62: 58,
                       63: 59, 64: 60, 65: 61}


def test_get_observable_label_from_header():
    """
    Test that the `get_observable_label_from_header` function correctly parses the header of the RMG simulation csv file
    to obtain the observable labels.
    """
    label = common.get_observable_label_from_header('dln[H(15)]/dln[k2]: O2(14)+H(15)(+M)<=>HO2(17)(+M)')
    assert label == 'H(15)'

    label = common.get_observable_label_from_header('dln[C2H4(12)]/dln[k16]: C8H17(11)<=>C8H17(5)')
    assert label == 'C2H4(12)'

    label = common.get_observable_label_from_header('dln[H(15)]/dG[CC(C)CO[O](237)]')
    assert label == 'H(15)'

    label = common.get_observable_label_from_header('dln[ethane(1)]/dln[k8]: H(6)+ethane(1)=H2(12)+C2H5(5)')
    assert label == 'ethane(1)'


def test_get_parameter_from_header():
    """
    Test that the `get_parameter_from_header` function correctly parses the header of the RMG simulation csv file
    to obtain the parameter labels.
    """
    label = common.get_parameter_from_header('dln[H(15)]/dln[k2]: O2(14)+H(15)(+M)<=>HO2(17)(+M)')
    assert label == 'k2'

    label = common.get_parameter_from_header('dln[C2H4(12)]/dln[k16]: C8H17(11)<=>C8H17(5)')
    assert label == 'k16'

    label = common.get_parameter_from_header('dln[H(15)]/dG[t-C4H9(60)]')
    assert label == 't-C4H9(60)'

    label = common.get_parameter_from_header('dln[H(15)]/dG[CC(C)CO[O](237)]')
    assert label == 'CC(C)CO[O](237)'

    label = common.get_parameter_from_header('dln[ethane(1)]/dG[C2H4(8)]')
    assert label == 'C2H4(8)'

    label = common.get_parameter_from_header('dln[ethane(1)]/dln[k8]: H(6)+ethane(1)=H2(12)+C2H5(5)')
    assert label == 'k8'
