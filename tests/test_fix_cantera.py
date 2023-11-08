#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_flux module
"""

import os
import shutil

from t3.common import t3_path, DATA_BASE_PATH
from tests.common import copy_model
import t3.utils.fix_cantera as fix


def test_fix_cantera_1():
    """Test fixing a Cantera file, case 1."""
    model_path = os.path.join(DATA_BASE_PATH, 'models', 'dups', '1_test_fix_cantera.yaml')
    shutil.copyfile(os.path.join(DATA_BASE_PATH, 'models', 'dups', '1.yaml'), model_path)
    done = fix.fix_cantera(model_path)
    assert done is True
    os.remove(model_path)
    os.remove(model_path + '.bak')


def test_fix_cantera_2():
    """Test fixing a Cantera file, case 2."""
    model_path = os.path.join(DATA_BASE_PATH, 'models', 'dups', '2_test_fix_cantera.yaml')
    shutil.copyfile(os.path.join(DATA_BASE_PATH, 'models', 'dups', '2.yaml'), model_path)
    done = fix.fix_cantera(model_path)
    assert done is True
    os.remove(model_path)
    os.remove(model_path + '.bak')


def test_fix_cantera_NH3():
    """Test fixing a Cantera file, NH3 case."""
    model_path = os.path.join(DATA_BASE_PATH, 'models', 'dups', 'NH3_test_fix_cantera.yaml')
    shutil.copyfile(os.path.join(DATA_BASE_PATH, 'models', 'dups', 'NH3.yaml'), model_path)
    done = fix.fix_cantera(model_path)
    assert done is True
    os.remove(model_path)
    os.remove(model_path + '.bak')


def test_fix_undeclared_duplicate_reactions():
    """Test fixing undeclared duplicate reactions."""
    model_path_1 = os.path.join(DATA_BASE_PATH, 'models', 'dups', '1_test_fix_undeclared_duplicate_reactions.yaml')
    shutil.copyfile(os.path.join(DATA_BASE_PATH, 'models', 'dups', '1.yaml'), model_path_1)
    tb_1 = fix.get_traceback(model_path_1)
    model_path_2 = copy_model(os.path.join(DATA_BASE_PATH, 'models', 'dups', '1.yaml'))
    fix.fix_undeclared_duplicate_reactions(model_path_2, tb_1)
    tb = fix.get_traceback(model_path_2)
    assert tb is None
    os.remove(model_path_1)
    os.remove(model_path_2)


def test_get_traceback():
    """Test getting a traceback."""
    model_path = os.path.join(DATA_BASE_PATH, 'models', 'HOCHO.yaml')
    tb = fix.get_traceback(model_path)
    assert tb is None

    model_path = copy_model(os.path.join(DATA_BASE_PATH, 'models', 'dups', '1.yaml'), name='C')
    shutil.copyfile(os.path.join(DATA_BASE_PATH, 'models', 'dups', '1.yaml'), model_path)
    tb = fix.get_traceback(model_path)
    expected_tb = f"""in get_traceback
    ct.Solution(model_path)
  File "build/python/cantera/base.pyx", line 71, in cantera._cantera._SolutionBase.__cinit__
  File "build/python/cantera/base.pyx", line 128, in cantera._cantera._SolutionBase._cinit
  File "build/python/cantera/base.pyx", line 215, in cantera._cantera._SolutionBase._init_yaml
cantera._cantera.CanteraError: 
*******************************************************************************
InputFileError thrown by Kinetics::checkDuplicates:
Error on lines 72 and 87 of {t3_path}/tests/data/models/dups/temp_C_1.yaml:
Undeclared duplicate reactions detected:
Reaction 2: 2 HNO(63) <=> 2 HNO(T)(117)
Reaction 1: HNO(T)(117) <=> HNO(63)

|  Line |
|    67 |     note: Epsilon & sigma estimated with Tc=544.31 K, Pc=43.8 bar (from
|    68 |       Joback method)
|    69 |   note: HNO(T)(117)
|    70 | 
|    71 | reactions:
>    72 > - equation: HNO(63) + HNO(63) <=> HNO(T)(117) + HNO(T)(117)  # Reaction 515
            ^
|    73 |   type: Chebyshev
|    74 |   temperature-range: [500.0, 1600.0]
|    75 |   pressure-range: [0.099 atm, 98.692 atm]
...
|    82 |   - [8.737e-04, -2.09e-10, -1.281e-10, -6.01e-11]
|    83 |   note: |-
|    84 |     Reaction index: Chemkin #515; RMG #1957
|    85 |     PDep reaction: PDepNetwork #87
|    86 |     Flux pairs: HNO(63), HNO(T)(117); HNO(63), HNO(T)(117);
>    87 > - equation: HNO(T)(117) <=> HNO(63)  # Reaction 561
            ^
|    88 |   type: Chebyshev
|    89 |   temperature-range: [500.0, 1600.0]
|    90 |   pressure-range: [0.099 atm, 98.692 atm]
*******************************************************************************

"""
    assert expected_tb in tb
    os.remove(model_path)

    model_path = copy_model(os.path.join(DATA_BASE_PATH, 'models', 'dups', '2.yaml'), name='B')
    tb = fix.get_traceback(model_path)
    expected_tb = f"""in get_traceback
    ct.Solution(model_path)
  File "build/python/cantera/base.pyx", line 71, in cantera._cantera._SolutionBase.__cinit__
  File "build/python/cantera/base.pyx", line 128, in cantera._cantera._SolutionBase._cinit
  File "build/python/cantera/base.pyx", line 215, in cantera._cantera._SolutionBase._init_yaml
cantera._cantera.CanteraError: 
*******************************************************************************
InputFileError thrown by Kinetics::checkDuplicates:
Error on line 149 of {t3_path}/tests/data/models/dups/temp_B_2.yaml:
No duplicate found for declared duplicate reaction number 0 (NH2O(93) + O2(3) <=> HNO(T)(117) + HO2(10))
|  Line |
|   144 |     note: Epsilon & sigma estimated with Tc=544.31 K, Pc=43.8 bar (from
|   145 |       Joback method)
|   146 |   note: HNO(T)(117)
|   147 | 
|   148 | reactions:
>   149 > - equation: O2(3) + NH2O(93) <=> HO2(10) + HNO(T)(117)  # Reaction 503
            ^
|   150 |   duplicate: true
|   151 |   rate-constant: {{A: 4429.0, b: 2.578, Ea: 29.877}}
|   152 |   note: |-
*******************************************************************************

"""
    assert expected_tb in tb
    os.remove(model_path)


def test_get_dup_rxn_indices():
    """Test getting the reaction indices from the traceback."""
    model_path = copy_model(os.path.join(DATA_BASE_PATH, 'models', 'HOCHO.yaml'))
    tb = fix.get_traceback(model_path)
    rxns = fix.get_dup_rxn_indices(tb)
    assert rxns == []
    os.remove(model_path)

    model_path = copy_model(os.path.join(DATA_BASE_PATH, 'models', 'dups', '1.yaml'))
    tb = fix.get_traceback(model_path)
    rxns = fix.get_dup_rxn_indices(tb)
    assert rxns == [2, 1]
    os.remove(model_path)


def test_get_mistakenly_marked_dup_rxns():
    """Test getting the reaction indices from the traceback."""
    model_path = copy_model(os.path.join(DATA_BASE_PATH, 'models', 'HOCHO.yaml'))
    tb = fix.get_traceback(model_path)
    rxns = fix.get_mistakenly_marked_dup_rxns(tb)
    assert rxns == []
    os.remove(model_path)

    tb = 'No duplicate found for declared duplicate reaction number 51 (NH2O(93) + O2(3) <=> HNO(T)(117) + HO2(10))'
    rxns = fix.get_mistakenly_marked_dup_rxns(tb)
    assert rxns == [51]

    model_path = copy_model(os.path.join(DATA_BASE_PATH, 'models', 'dups', '2.yaml'), name='A')
    tb = fix.get_traceback(model_path)
    rxns = fix.get_mistakenly_marked_dup_rxns(tb)
    assert rxns == [0]
    os.remove(model_path)
