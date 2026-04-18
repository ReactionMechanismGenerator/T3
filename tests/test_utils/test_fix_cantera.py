#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_fix_cantera module
"""

import os
import shutil

from t3.common import TEST_DATA_BASE_PATH
from tests.common import copy_model
import t3.utils.fix_cantera as fix


def test_fix_cantera_1():
    """Test fixing a Cantera file, case 1."""
    model_path = os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', '1_test_fix_cantera.yaml')
    shutil.copyfile(os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', '1.yaml'), model_path)
    done = fix.fix_cantera(model_path)
    assert done is True
    os.remove(model_path)
    os.remove(model_path + '.bak')


def test_fix_cantera_2():
    """Test fixing a Cantera file, case 2."""
    model_path = os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', '2_test_fix_cantera.yaml')
    shutil.copyfile(os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', '2.yaml'), model_path)
    done = fix.fix_cantera(model_path)
    assert done is True
    os.remove(model_path)
    os.remove(model_path + '.bak')


def test_fix_cantera_NH3():
    """Test fixing a Cantera file, NH3 case."""
    model_path = os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', 'NH3_test_fix_cantera.yaml')
    shutil.copyfile(os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', 'NH3.yaml'), model_path)
    done = fix.fix_cantera(model_path)
    assert done is True
    os.remove(model_path)
    os.remove(model_path + '.bak')


def test_fix_undeclared_duplicate_reactions():
    """Test fixing undeclared duplicate reactions."""
    model_path_1 = os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', '1_test_fix_undeclared_duplicate_reactions.yaml')
    shutil.copyfile(os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', '1.yaml'), model_path_1)
    tb_1 = fix.get_traceback(model_path_1)
    model_path_2 = copy_model(os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', '1.yaml'))
    fix.fix_undeclared_duplicate_reactions(model_path_2, tb_1, marked_dups=[])
    tb = fix.get_traceback(model_path_2)
    assert tb is None
    os.remove(model_path_1)
    os.remove(model_path_2)


def test_get_traceback():
    """Test getting a traceback."""
    model_path = os.path.join(TEST_DATA_BASE_PATH, 'models', 'HOCHO.yaml')
    tb = fix.get_traceback(model_path)
    assert tb is None

    # Case 1: Undeclared duplicates
    model_path = copy_model(os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', '1.yaml'), name='C')
    shutil.copyfile(os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', '1.yaml'), model_path)
    tb = fix.get_traceback(model_path)

    # Robust assertions: Check for key error components instead of exact string match
    assert "CanteraError" in tb
    assert "Undeclared duplicate reactions detected" in tb
    # Check that the specific problematic equations are mentioned, regardless of order
    assert "HNO(63) + HNO(63) <=> HNO(T)(117) + HNO(T)(117)" in tb
    assert "HNO(T)(117) <=> HNO(63)" in tb

    os.remove(model_path)

    # Case 2: No duplicate found for declared duplicate
    model_path = copy_model(os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', '2.yaml'), name='B')
    tb = fix.get_traceback(model_path)

    # Robust assertions
    assert "CanteraError" in tb
    assert "No duplicate found for declared duplicate reaction" in tb
    assert "NH2O(93) + O2(3) <=> HNO(T)(117) + HO2(10)" in tb

    os.remove(model_path)


def test_get_dup_rxn_indices():
    """Test getting the reaction indices from the traceback."""
    model_path = copy_model(os.path.join(TEST_DATA_BASE_PATH, 'models', 'HOCHO.yaml'))
    tb = fix.get_traceback(model_path)
    rxns = fix.get_dup_rxn_indices(tb)
    assert rxns == []
    os.remove(model_path)

    model_path = copy_model(os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', '1.yaml'))
    tb = fix.get_traceback(model_path)
    rxns = fix.get_dup_rxn_indices(tb)
    assert rxns == [1, 2]
    os.remove(model_path)


def test_get_mistakenly_marked_dup_rxns():
    """Test getting the reaction indices from the traceback."""
    model_path = copy_model(os.path.join(TEST_DATA_BASE_PATH, 'models', 'HOCHO.yaml'))
    tb = fix.get_traceback(model_path)
    rxns = fix.get_mistakenly_marked_dup_rxns(tb)
    assert rxns == []
    os.remove(model_path)

    tb = 'No duplicate found for declared duplicate reaction number 51 (NH2O(93) + O2(3) <=> HNO(T)(117) + HO2(10))'
    rxns = fix.get_mistakenly_marked_dup_rxns(tb)
    assert rxns == [51]

    model_path = copy_model(os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups', '2.yaml'), name='A')
    tb = fix.get_traceback(model_path)
    rxns = fix.get_mistakenly_marked_dup_rxns(tb)
    assert rxns == [0]
    os.remove(model_path)


def teardown_module():
    """Safety net: remove any temp_*.yaml / *_test_fix_cantera*.yaml /
    *.bak files left under tests/data/models/ if a test failed mid-run."""
    models_dirs = [os.path.join(TEST_DATA_BASE_PATH, 'models'),
                   os.path.join(TEST_DATA_BASE_PATH, 'models', 'dups')]
    for d in models_dirs:
        if not os.path.isdir(d):
            continue
        for name in os.listdir(d):
            if name.startswith('temp_') or 'test_fix_cantera' in name or name.endswith('.bak'):
                try:
                    os.remove(os.path.join(d, name))
                except OSError:
                    pass

