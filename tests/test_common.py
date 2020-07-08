#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_utils module
"""
import os

import t3.common as common
from t3.common import DATA_BASE_PATH
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
