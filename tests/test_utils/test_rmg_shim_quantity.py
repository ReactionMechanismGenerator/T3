#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_utils test_rmg_shim_quantity module
"""

import pytest

from t3.utils.rmg_shim import Quantity


@pytest.mark.parametrize('units', ['K', 'Pa', 's'])
def test_quantity_accepts_si_units(units):
    """A Quantity in an SI unit keeps its value as value_si."""
    quantity = Quantity(2.5, units)
    assert quantity.value_si == 2.5
    assert quantity.units == units


@pytest.mark.parametrize('units', ['bar', 'atm', 'ms', 'hours', None])
def test_quantity_rejects_non_si_units(units):
    """A Quantity refuses any unit it cannot hold as SI, naming the unit."""
    with pytest.raises(ValueError, match=repr(units)):
        Quantity(1.0, units)
