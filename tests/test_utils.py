"""
t3 tests test_utils module
"""

import tandem.utils as utils


def test_dict_to_str():
    """Test prettifying a dictionary"""
    dictionary = {'label1': {'spc': 'Species1', 'reason': 'Reason1'},
                  'label2': {'spc': 'Species2', 'reason': 'Reason2'}}
    output = utils.dict_to_str(dictionary)
    expected_output = """label1:
  spc: Species1
  reason: Reason1
label2:
  spc: Species2
  reason: Reason2
"""
    assert output == expected_output


def test_combine_dicts():
    """Test combining two dictionaries"""
    dict1 = {1: 7, 3: 9}
    dict2 = {5: 2, 4: 8, 1: 7}

    combined = utils.combine_dicts(dict1=dict1, dict2=None)
    assert combined == dict1

    combined = utils.combine_dicts(dict1=dict1, dict2=dict())
    assert combined == dict1

    combined = utils.combine_dicts(dict1=None, dict2=dict1)
    assert combined == dict1

    combined = utils.combine_dicts(dict1=dict1, dict2=dict2)
    assert combined == {1: 7, 3: 9, 5: 2, 4: 8}
