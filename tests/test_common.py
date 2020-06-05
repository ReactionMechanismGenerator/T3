#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_utils module
"""

import t3.common as common


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
