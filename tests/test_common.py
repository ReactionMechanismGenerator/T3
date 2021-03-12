#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_utils module
"""

import os
import pytest

from rmgpy.species import Species

import t3.common as common
from t3.common import DATA_BASE_PATH
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
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'minimal_data'),
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
    t_final = [1, 'micro-s']
    assert common.convert_termination_time_to_seconds(t_final) == 1e-6

    t_final = [1, 'ms']
    assert common.convert_termination_time_to_seconds(t_final) == 1e-3

    t_final = [1, 's']
    assert common.convert_termination_time_to_seconds(t_final) == 1

    t_final = [1, 'hours']
    assert common.convert_termination_time_to_seconds(t_final) == 3600

    t_final = [1, 'hrs']
    assert common.convert_termination_time_to_seconds(t_final) == 3600

    t_final = [1, 'days']
    assert common.convert_termination_time_to_seconds(t_final) == 3600*24
