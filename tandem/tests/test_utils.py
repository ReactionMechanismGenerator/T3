"""
ts tests test_utils module
"""

import unittest

import tandem.utils as utils


class TestUtils(unittest.TestCase):
    """
    Contains unit tests for the T3 utils module
    """

    def test_dict_to_str(self):
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
        self.assertEqual(output, expected_output)

    def test_combine_dicts(self):
        """Test combining two dictionaries"""
        dict1 = {1: 7, 3: 9}
        dict2 = {5: 2, 4: 8, 1: 7}
        combined = utils.combine_dicts(dict1=dict1, dict2=None)
        self.assertEqual(combined, dict1)
        combined = utils.combine_dicts(dict1=dict1, dict2=dict())
        self.assertEqual(combined, dict1)
        combined = utils.combine_dicts(dict1=None, dict2=dict1)
        self.assertEqual(combined, dict1)
        combined = utils.combine_dicts(dict1=dict1, dict2=dict2)
        self.assertEqual(combined, {1: 7, 3: 9, 5: 2, 4: 8})
