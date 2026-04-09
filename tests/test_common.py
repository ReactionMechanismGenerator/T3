#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_utils module
"""

import os
import shutil
import tempfile

import pytest

from arc.common import read_yaml_file

from t3.chem import T3Species
import t3.common as common
from t3.common import TEST_DATA_BASE_PATH
from t3.schema import RMGSpecies
from t3.utils.writer import get_species_obj_from_a_species_dict


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
    t3_species = [T3Species(label='H2O', key=7, smiles='O'),
                  T3Species(label='CH4', key=1, smiles='C')]
    label = 'H2O'
    species = common.get_species_by_label(label, t3_species)
    assert species.label == label
    assert species.key == 7

    species = common.get_species_by_label('CH4', t3_species)
    assert species.label == 'CH4'
    assert species.key == 1


def test_get_rmg_species_from_a_species_dict():
    """Test getting an RMG species from a species dictionary"""
    species = get_species_obj_from_a_species_dict(
        species_dict=RMGSpecies(**{'label': 'spc', 'smiles': 'C=O'}).dict())
    assert isinstance(species, T3Species)
    assert species.label == 'spc'
    assert species.mol.to_smiles() == 'C=O'

    adj = """1 C u0 p0 c0 {2,D} {3,S} {4,S}
2 O u0 p2 c0 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}"""
    species = get_species_obj_from_a_species_dict(
        species_dict=RMGSpecies(**{'label': 'spc', 'adjlist': adj}).dict())
    assert isinstance(species, T3Species)
    assert species.label == 'spc'
    assert species.mol.to_smiles() == 'C=O'

    species = get_species_obj_from_a_species_dict(
        species_dict=RMGSpecies(**{'label': 'spc', 'inchi': 'InChI=1S/CH2O/c1-2/h1H2'}).dict())
    assert isinstance(species, T3Species)
    assert species.label == 'spc'
    assert species.mol.to_smiles() == 'C=O'

    xyz = """O  0.0000000  0.0000000  0.7047750
C  0.0000000  0.0000000 -0.5593030
H  0.0000000  0.9470590 -1.1411940
H  0.0000000 -0.9470590 -1.1411940"""
    species = get_species_obj_from_a_species_dict(
        species_dict=RMGSpecies(**{'label': 'spc', 'xyz': [xyz]}).dict())
    assert isinstance(species, T3Species)
    assert species.label == 'spc'
    assert species.mol.to_smiles() == 'C=O'

    species = get_species_obj_from_a_species_dict(species_dict=RMGSpecies(**{'label': 'spc'}).dict(),
                                                  raise_error=False)
    assert species is None
    with pytest.raises(ValueError):
        get_species_obj_from_a_species_dict(species_dict=RMGSpecies(**{'label': 'spc'}).dict(),
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

    # Cantera thermo SA uses /dH[ (enthalpy perturbation)
    label = common.get_parameter_from_header('dln[H2(1)]/dH[OH(3)]')
    assert label == 'OH(3)'

    label = common.get_parameter_from_header('dln[NH2(6)]/dH[NH3(1)]')
    assert label == 'NH3(1)'


def test_save_yaml_file_top_keys():
    """Test that save_yaml_file writes top_keys first in the specified order."""
    tmp_dir = tempfile.mkdtemp()
    path = os.path.join(tmp_dir, 'test.yml')
    content = {
        'zebra': [1, 2, 3],
        'alpha': 'first',
        'metadata': {'adapter': 'CanteraConstantTP', 'iteration': 1},
        'beta': {'nested': True},
    }
    common.save_yaml_file(path=path, content=content, top_keys=['metadata', 'alpha'])
    with open(path, 'r') as f:
        lines = f.readlines()
    # Metadata block must start on the very first line
    assert lines[0].startswith('metadata:')
    # Find where each top-level key appears
    top_level_keys = [line.split(':')[0] for line in lines
                      if line and not line[0].isspace() and ':' in line]
    assert top_level_keys[0] == 'metadata'
    assert top_level_keys[1] == 'alpha'
    # The remaining keys (zebra, beta) follow in insertion order
    assert 'zebra' in top_level_keys
    assert 'beta' in top_level_keys

    # Verify the file is valid YAML that round-trips correctly
    loaded = read_yaml_file(path=path)
    assert loaded['metadata'] == content['metadata']
    assert loaded['alpha'] == content['alpha']
    assert loaded['zebra'] == content['zebra']
    assert loaded['beta'] == content['beta']
    shutil.rmtree(tmp_dir)


def test_save_yaml_file_no_top_keys():
    """Test that save_yaml_file without top_keys preserves insertion order."""
    tmp_dir = tempfile.mkdtemp()
    path = os.path.join(tmp_dir, 'test.yml')
    content = {'charlie': 3, 'alpha': 1, 'beta': 2}
    common.save_yaml_file(path=path, content=content)
    with open(path, 'r') as f:
        lines = f.readlines()
    top_level_keys = [line.split(':')[0] for line in lines
                      if line.strip() and not line[0].isspace() and ':' in line]
    # Without top_keys, keys should be in insertion order (not alphabetized)
    assert top_level_keys == ['charlie', 'alpha', 'beta']
    shutil.rmtree(tmp_dir)


def test_save_yaml_file_list_content():
    """save_yaml_file should accept a list as the top-level content."""
    tmp_dir = tempfile.mkdtemp()
    path = os.path.join(tmp_dir, 'list.yml')
    content = [{'label': 'OH', 'idx': 1}, {'label': 'H2O', 'idx': 2}]
    common.save_yaml_file(path=path, content=content)
    loaded = read_yaml_file(path=path)
    assert loaded == content
    shutil.rmtree(tmp_dir)


def test_save_yaml_file_list_content_with_top_keys_ignored():
    """When content is a list, top_keys must be ignored (no crash)."""
    tmp_dir = tempfile.mkdtemp()
    path = os.path.join(tmp_dir, 'list_top.yml')
    content = [1, 2, 3]
    common.save_yaml_file(path=path, content=content, top_keys=['anything'])
    loaded = read_yaml_file(path=path)
    assert loaded == content
    shutil.rmtree(tmp_dir)


def test_save_yaml_file_creates_missing_parent_dir():
    """save_yaml_file should create a missing parent directory."""
    tmp_dir = tempfile.mkdtemp()
    nested = os.path.join(tmp_dir, 'a', 'b', 'c')
    path = os.path.join(nested, 'out.yml')
    assert not os.path.exists(nested)
    common.save_yaml_file(path=path, content={'k': 'v'})
    assert os.path.isfile(path)
    assert read_yaml_file(path=path) == {'k': 'v'}
    shutil.rmtree(tmp_dir)


def test_save_yaml_file_no_parent_dir_in_path():
    """save_yaml_file should not crash when path has no directory component."""
    tmp_dir = tempfile.mkdtemp()
    cwd = os.getcwd()
    try:
        os.chdir(tmp_dir)
        common.save_yaml_file(path='bare.yml', content={'k': 'v'})
        assert os.path.isfile(os.path.join(tmp_dir, 'bare.yml'))
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp_dir)


def test_save_yaml_file_top_keys_covering_all_keys():
    """When top_keys covers every key, no remainder block should be written."""
    tmp_dir = tempfile.mkdtemp()
    path = os.path.join(tmp_dir, 'all_top.yml')
    content = {'b': 2, 'a': 1}
    common.save_yaml_file(path=path, content=content, top_keys=['a', 'b'])
    with open(path, 'r') as f:
        lines = f.readlines()
    top_level_keys = [line.split(':')[0] for line in lines
                      if line.strip() and not line[0].isspace() and ':' in line]
    assert top_level_keys == ['a', 'b']
    assert read_yaml_file(path=path) == content
    shutil.rmtree(tmp_dir)


def test_save_yaml_file_top_keys_missing_from_content():
    """top_keys may name keys that aren't in content; they should be skipped."""
    tmp_dir = tempfile.mkdtemp()
    path = os.path.join(tmp_dir, 'missing_top.yml')
    content = {'real': 1, 'also_real': 2}
    common.save_yaml_file(path=path, content=content,
                          top_keys=['ghost', 'real', 'phantom'])
    with open(path, 'r') as f:
        lines = f.readlines()
    top_level_keys = [line.split(':')[0] for line in lines
                      if line.strip() and not line[0].isspace() and ':' in line]
    assert top_level_keys[0] == 'real'
    assert 'ghost' not in top_level_keys
    assert 'phantom' not in top_level_keys
    assert read_yaml_file(path=path) == content
    shutil.rmtree(tmp_dir)


def test_save_yaml_file_empty_dict():
    """save_yaml_file should handle an empty dict without error."""
    tmp_dir = tempfile.mkdtemp()
    path = os.path.join(tmp_dir, 'empty.yml')
    common.save_yaml_file(path=path, content={})
    assert os.path.isfile(path)
    loaded = read_yaml_file(path=path)
    # Empty YAML may load as {} or None depending on the writer
    assert loaded in ({}, None)
    shutil.rmtree(tmp_dir)


def test_save_yaml_file_nested_round_trip():
    """save_yaml_file output should round-trip nested structures faithfully."""
    tmp_dir = tempfile.mkdtemp()
    path = os.path.join(tmp_dir, 'nested.yml')
    content = {
        'meta': {'iter': 3, 'tags': ['a', 'b']},
        'data': [{'x': 1.5, 'y': [10, 20]}, {'x': 2.5, 'y': [30, 40]}],
    }
    common.save_yaml_file(path=path, content=content, top_keys=['meta'])
    loaded = read_yaml_file(path=path)
    assert loaded == content
    shutil.rmtree(tmp_dir)


# ---------------------------------------------------------------------------
# Tests for the new helpers (remove_numeric_parentheses, numpy_to_list,
# get_atom_counts, get_oxidizer_stoichiometry, determine_concentrations_by_equivalence_ratios)
# ---------------------------------------------------------------------------

def test_remove_numeric_parentheses():
    assert common.remove_numeric_parentheses('OH(12)') == 'OH'
    assert common.remove_numeric_parentheses('CH3CH2(245)') == 'CH3CH2'
    assert common.remove_numeric_parentheses('plain') == 'plain'
    assert common.remove_numeric_parentheses('') == ''
    # Only trailing parens are stripped:
    assert common.remove_numeric_parentheses('A(1)B(2)') == 'A(1)B'
    assert common.remove_numeric_parentheses('S(40)') == 'S'


def test_numpy_to_list():
    import numpy as np
    assert common.numpy_to_list(np.array([1, 2, 3])) == [1, 2, 3]
    assert common.numpy_to_list(np.array([[1, 2], [3, 4]])) == [[1, 2], [3, 4]]
    assert common.numpy_to_list([1, 2, 3]) == [1, 2, 3]
    assert common.numpy_to_list(None) is None
    assert common.numpy_to_list(7.5) == 7.5


def test_get_atom_counts():
    """Atom counts should match formulas regardless of identifier kind."""
    assert common.get_atom_counts(smiles='C') == {'C': 1, 'H': 4, 'N': 0, 'O': 0, 'other': 0}
    assert common.get_atom_counts(smiles='CCO') == {'C': 2, 'H': 6, 'N': 0, 'O': 1, 'other': 0}
    assert common.get_atom_counts(smiles='N') == {'C': 0, 'H': 3, 'N': 1, 'O': 0, 'other': 0}
    assert common.get_atom_counts(smiles='[O][O]') == {'C': 0, 'H': 0, 'N': 0, 'O': 2, 'other': 0}
    assert common.get_atom_counts(smiles='OO') == {'C': 0, 'H': 2, 'N': 0, 'O': 2, 'other': 0}  # H2O2
    assert common.get_atom_counts(smiles='N#N') == {'C': 0, 'H': 0, 'N': 2, 'O': 0, 'other': 0}
    # Adjacency list path:
    assert common.get_atom_counts(adjlist="""1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}""") == {'C': 1, 'H': 4, 'N': 0, 'O': 0, 'other': 0}
    # Missing identifier:
    with pytest.raises(ValueError):
        common.get_atom_counts()


def test_get_oxidizer_stoichiometry_o2():
    """Hydrocarbon + O2 stoichiometry should match hand-balanced reactions."""
    # CH4 + 2 O2
    assert common.get_oxidizer_stoichiometry(fuel_smiles='C', oxidizer_smiles='[O][O]') == 2.0
    # C3H8 + 5 O2
    assert common.get_oxidizer_stoichiometry(fuel_smiles='CCC', oxidizer_smiles='[O][O]') == 5.0
    # C5H12 + 8 O2
    assert common.get_oxidizer_stoichiometry(fuel_smiles='CCCCC', oxidizer_smiles='[O][O]') == 8.0
    # C6H14 + 9.5 O2
    assert common.get_oxidizer_stoichiometry(fuel_smiles='CCCCCC', oxidizer_smiles='[O][O]') == 9.5
    # C2H5OH (CCO) + 3 O2 → 2 CO2 + 3 H2O
    assert common.get_oxidizer_stoichiometry(fuel_smiles='CCO', oxidizer_smiles='[O][O]') == 3.0
    # NH3 + 0.75 O2 → 0.5 N2 + 1.5 H2O   (4 NH3 + 3 O2 → 2 N2 + 6 H2O)
    assert common.get_oxidizer_stoichiometry(fuel_smiles='N', oxidizer_smiles='[O][O]') == 0.75
    # Ethylamine NCC: a=2, b=7, c=1, d=0  ⇒ (4 + 3.5)/2 = 3.75
    assert common.get_oxidizer_stoichiometry(fuel_smiles='NCC', oxidizer_smiles='[O][O]') == 3.75


def test_get_oxidizer_stoichiometry_h2o2():
    """Fuel + H2O2 stoichiometry: each H2O2 supplies 1 net active O."""
    # C3H8 + 10 H2O2 → 3 CO2 + 14 H2O
    assert common.get_oxidizer_stoichiometry(fuel_smiles='CCC', oxidizer_smiles='OO') == 10.0
    # CH4 + 4 H2O2 → CO2 + 6 H2O
    assert common.get_oxidizer_stoichiometry(fuel_smiles='C', oxidizer_smiles='OO') == 4.0
    # 2 NH3 + 3 H2O2 → N2 + 6 H2O   ⇒ 1.5 H2O2 per NH3
    assert common.get_oxidizer_stoichiometry(fuel_smiles='N', oxidizer_smiles='OO') == 1.5


def test_get_oxidizer_stoichiometry_n2o():
    """N2O has 1 O atom per molecule and no H — each supplies 1 net active O."""
    # CH4 + 4 N2O → CO2 + 2 H2O + 4 N2
    assert common.get_oxidizer_stoichiometry(fuel_smiles='C', oxidizer_smiles='[N-]=[N+]=O') == 4.0


def test_get_oxidizer_stoichiometry_rejects_bad_inputs():
    # Sulfur-containing fuel: not C/H/N/O
    with pytest.raises(ValueError, match='not C/H/N/O'):
        common.get_oxidizer_stoichiometry(fuel_smiles='CS', oxidizer_smiles='[O][O]')
    # Water as "oxidizer" supplies no net active O (k_O - k_H/2 = 1 - 1 = 0)
    with pytest.raises(ValueError, match='no net active oxygen'):
        common.get_oxidizer_stoichiometry(fuel_smiles='C', oxidizer_smiles='O')


def test_determine_concentrations_no_fuel():
    """No fuel role → returns None."""
    species = [{'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'}]
    assert common.determine_concentrations_by_equivalence_ratios(species) is None


def test_determine_concentrations_no_equivalence_ratios():
    """Fuel role but no equivalence_ratios → returns None."""
    species = [
        {'label': 'CH4', 'smiles': 'C', 'concentration': 1, 'role': 'fuel'},
        {'label': 'O2', 'smiles': '[O][O]', 'concentration': 2, 'role': 'oxidizer'},
    ]
    assert common.determine_concentrations_by_equivalence_ratios(species) is None


def test_determine_concentrations_propane_air():
    """Propane in air at φ ∈ {0.5, 1.0, 1.5} — checks the corrected (stoich/φ) formula."""
    species = [
        {'label': 'propane', 'smiles': 'CCC', 'role': 'fuel',
         'equivalence_ratios': [0.5, 1.0, 1.5], 'concentration': 0},
        {'label': 'O2', 'smiles': '[O][O]', 'role': 'oxidizer', 'concentration': 0},
        {'label': 'N2', 'smiles': 'N#N', 'role': 'diluent', 'concentration': 0},
    ]
    sweep = common.determine_concentrations_by_equivalence_ratios(species)
    assert sweep is not None
    assert sweep['equivalence_ratios'] == [0.5, 1.0, 1.5]
    assert sweep['concentrations']['propane'] == [1.0, 1.0, 1.0]
    # Stoich for C3H8 = 5; corrected formula: [O2] = 5 / φ
    assert sweep['concentrations']['O2'] == pytest.approx([10.0, 5.0, 10.0 / 3.0])
    # Diluent N2 at default 3.76 ratio
    expected_n2 = [10.0 * 3.76, 5.0 * 3.76, (10.0 / 3.0) * 3.76]
    assert sweep['concentrations']['N2'] == pytest.approx(expected_n2)


def test_determine_concentrations_methane_air_single_phi():
    """Single equivalence ratio still produces a 1-element column."""
    species = [
        {'label': 'methane', 'smiles': 'C', 'role': 'fuel',
         'equivalence_ratios': [1.0], 'concentration': 0},
        {'label': 'O2', 'smiles': '[O][O]', 'role': 'oxidizer', 'concentration': 0},
        {'label': 'N2', 'smiles': 'N#N', 'role': 'diluent', 'concentration': 0},
    ]
    sweep = common.determine_concentrations_by_equivalence_ratios(species)
    assert sweep['equivalence_ratios'] == [1.0]
    assert sweep['concentrations']['methane'] == [1.0]
    assert sweep['concentrations']['O2'] == [2.0]
    assert sweep['concentrations']['N2'] == pytest.approx([2.0 * 3.76])


def test_determine_concentrations_ammonia_air():
    """Ammonia / air — small stoich, verify lean/rich direction is correct."""
    species = [
        {'label': 'NH3', 'smiles': 'N', 'role': 'fuel',
         'equivalence_ratios': [0.5, 1.0, 1.5], 'concentration': 0},
        {'label': 'O2', 'smiles': '[O][O]', 'role': 'oxidizer', 'concentration': 0},
    ]
    sweep = common.determine_concentrations_by_equivalence_ratios(species)
    # Stoich for NH3 = 0.75; corrected: [O2] = 0.75 / φ
    assert sweep['concentrations']['O2'] == pytest.approx([1.5, 0.75, 0.5])


def test_determine_concentrations_explicit_fuel_concentration():
    """Non-default fuel concentration scales every column linearly."""
    species = [
        {'label': 'CH4', 'smiles': 'C', 'role': 'fuel',
         'equivalence_ratios': [1.0, 2.0], 'concentration': 3.0},
        {'label': 'O2', 'smiles': '[O][O]', 'role': 'oxidizer', 'concentration': 0},
    ]
    sweep = common.determine_concentrations_by_equivalence_ratios(species)
    assert sweep['concentrations']['CH4'] == [3.0, 3.0]
    # 3 * 2 / 1 = 6 ; 3 * 2 / 2 = 3
    assert sweep['concentrations']['O2'] == pytest.approx([6.0, 3.0])


def test_determine_concentrations_h2o2_oxidizer():
    """Propane + H2O2 (no O2) — alkane / non-O2 oxidizer mixture."""
    species = [
        {'label': 'C3H8', 'smiles': 'CCC', 'role': 'fuel',
         'equivalence_ratios': [0.5, 1.0, 2.0], 'concentration': 0},
        {'label': 'H2O2', 'smiles': 'OO', 'role': 'oxidizer', 'concentration': 0},
    ]
    sweep = common.determine_concentrations_by_equivalence_ratios(species)
    # Stoich C3H8 + 10 H2O2; corrected: [H2O2] = 10 / φ
    assert sweep['concentrations']['H2O2'] == pytest.approx([20.0, 10.0, 5.0])


def test_determine_concentrations_multi_oxidizer():
    """Multi-oxidizer mixture: 60% O2 + 40% H2O2 share the demand."""
    species = [
        {'label': 'CH4', 'smiles': 'C', 'role': 'fuel',
         'equivalence_ratios': [1.0], 'concentration': 0},
        {'label': 'O2', 'smiles': '[O][O]', 'role': 'oxidizer',
         'concentration': 0, 'oxidizer_fraction': 0.6},
        {'label': 'H2O2', 'smiles': 'OO', 'role': 'oxidizer',
         'concentration': 0, 'oxidizer_fraction': 0.4},
    ]
    sweep = common.determine_concentrations_by_equivalence_ratios(species)
    # CH4: stoich for O2 = 2, for H2O2 = 4. Each gets its own stoich * fraction / φ.
    assert sweep['concentrations']['O2'] == pytest.approx([2.0 * 0.6])
    assert sweep['concentrations']['H2O2'] == pytest.approx([4.0 * 0.4])


def test_determine_concentrations_multi_oxidizer_bad_fractions():
    """Oxidizer fractions must sum to 1.0."""
    species = [
        {'label': 'CH4', 'smiles': 'C', 'role': 'fuel',
         'equivalence_ratios': [1.0], 'concentration': 0},
        {'label': 'O2', 'smiles': '[O][O]', 'role': 'oxidizer',
         'concentration': 0, 'oxidizer_fraction': 0.7},
        {'label': 'H2O2', 'smiles': 'OO', 'role': 'oxidizer',
         'concentration': 0, 'oxidizer_fraction': 0.5},
    ]
    with pytest.raises(ValueError, match='must sum to 1.0'):
        common.determine_concentrations_by_equivalence_ratios(species)


def test_determine_concentrations_multi_diluent_custom_ratios():
    """Two diluents, each with its own diluent_to_oxidizer_ratio."""
    species = [
        {'label': 'CH4', 'smiles': 'C', 'role': 'fuel',
         'equivalence_ratios': [1.0], 'concentration': 0},
        {'label': 'O2', 'smiles': '[O][O]', 'role': 'oxidizer', 'concentration': 0},
        {'label': 'N2', 'smiles': 'N#N', 'role': 'diluent',
         'concentration': 0, 'diluent_to_oxidizer_ratio': 3.76},
        {'label': 'Ar', 'smiles': '[Ar]', 'role': 'diluent',
         'concentration': 0, 'diluent_to_oxidizer_ratio': 1.0},
    ]
    sweep = common.determine_concentrations_by_equivalence_ratios(species)
    assert sweep['concentrations']['O2'] == [2.0]
    assert sweep['concentrations']['N2'] == pytest.approx([2.0 * 3.76])
    assert sweep['concentrations']['Ar'] == pytest.approx([2.0 * 1.0])


def test_determine_concentrations_diluent_default_ratio():
    """A diluent without an explicit ratio defaults to 3.76 (air-like)."""
    species = [
        {'label': 'CH4', 'smiles': 'C', 'role': 'fuel',
         'equivalence_ratios': [1.0], 'concentration': 0},
        {'label': 'O2', 'smiles': '[O][O]', 'role': 'oxidizer', 'concentration': 0},
        {'label': 'N2', 'smiles': 'N#N', 'role': 'diluent', 'concentration': 0},
    ]
    sweep = common.determine_concentrations_by_equivalence_ratios(species)
    assert sweep['concentrations']['N2'] == pytest.approx([2.0 * 3.76])


def test_determine_concentrations_inert_passthrough():
    """A non-zero, no-role species is passed through as a constant column."""
    species = [
        {'label': 'CH4', 'smiles': 'C', 'role': 'fuel',
         'equivalence_ratios': [1.0, 0.5], 'concentration': 0},
        {'label': 'O2', 'smiles': '[O][O]', 'role': 'oxidizer', 'concentration': 0},
        {'label': 'He', 'smiles': '[He]', 'role': None, 'concentration': 0.5},
    ]
    sweep = common.determine_concentrations_by_equivalence_ratios(species)
    assert sweep['concentrations']['He'] == [0.5, 0.5]


def test_determine_concentrations_two_fuels_rejected():
    species = [
        {'label': 'CH4', 'smiles': 'C', 'role': 'fuel',
         'equivalence_ratios': [1.0], 'concentration': 0},
        {'label': 'C2H6', 'smiles': 'CC', 'role': 'fuel',
         'equivalence_ratios': [1.0], 'concentration': 0},
        {'label': 'O2', 'smiles': '[O][O]', 'role': 'oxidizer', 'concentration': 0},
    ]
    with pytest.raises(ValueError, match='Exactly one species'):
        common.determine_concentrations_by_equivalence_ratios(species)


def test_determine_concentrations_fuel_without_oxidizer_rejected():
    species = [
        {'label': 'CH4', 'smiles': 'C', 'role': 'fuel',
         'equivalence_ratios': [1.0], 'concentration': 0},
    ]
    with pytest.raises(ValueError, match='no.*role="oxidizer"'):
        common.determine_concentrations_by_equivalence_ratios(species)


def test_determine_concentrations_empty_eq_ratios_rejected():
    species = [
        {'label': 'CH4', 'smiles': 'C', 'role': 'fuel',
         'equivalence_ratios': [], 'concentration': 0},
        {'label': 'O2', 'smiles': '[O][O]', 'role': 'oxidizer', 'concentration': 0},
    ]
    with pytest.raises(ValueError, match='empty equivalence_ratios'):
        common.determine_concentrations_by_equivalence_ratios(species)
