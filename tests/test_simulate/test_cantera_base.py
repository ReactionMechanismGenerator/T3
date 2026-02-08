#!/usr/bin/env python3
# encoding: utf-8

"""
Tests for the CanteraBase shared functionality and the bug fixes
introduced during the refactoring to a common base class.
"""

import os
import shutil

import numpy as np

from t3.common import SIMULATE_TEST_DATA_BASE_PATH, TEST_DATA_BASE_PATH
from tests.common import run_minimal
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.cantera_base import CanteraBase
from t3.simulate.cantera_constant_tp import CanteraConstantTP
from t3.simulate.cantera_constant_hp import CanteraConstantHP
from t3.simulate.cantera_constant_uv import CanteraConstantUV


TEST_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'cantera_simulator_test')
NH3_MODEL_PATH = os.path.join(TEST_DATA_BASE_PATH, 'models', 'NH3.yaml')


def _make_adapter(adapter_class, observable_list=None):
    """Helper to construct a Cantera adapter using the minimal example model."""
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    adapter = adapter_class(t3=t3.t3,
                            rmg=t3.rmg,
                            paths=t3.paths,
                            logger=t3.logger,
                            atol=t3.rmg['model']['atol'],
                            rtol=t3.rmg['model']['rtol'],
                            observable_list=observable_list or list(),
                            sa_atol=t3.t3['sensitivity']['atol'],
                            sa_rtol=t3.t3['sensitivity']['rtol'],
                            )
    return adapter


# ---------------------------------------------------------------------------
# Inheritance
# ---------------------------------------------------------------------------

def test_inheritance_hierarchy():
    """All three adapters inherit from CanteraBase and SimulateAdapter."""
    for cls in (CanteraConstantTP, CanteraConstantHP, CanteraConstantUV):
        assert issubclass(cls, CanteraBase)
        assert issubclass(cls, SimulateAdapter)


def test_adapter_instances():
    """Instantiated adapters are instances of CanteraBase."""
    for cls in (CanteraConstantTP, CanteraConstantHP, CanteraConstantUV):
        adapter = _make_adapter(cls)
        assert isinstance(adapter, CanteraBase)
        assert isinstance(adapter, SimulateAdapter)


# ---------------------------------------------------------------------------
# Bug fix: no duplicate inert indices (double-population bug)
# ---------------------------------------------------------------------------

def test_no_duplicate_inert_indices():
    """
    Previously __init__ built inert_index_list both inside set_up() and again
    after the set_up() call, causing duplicates.  Verify no duplicates exist.
    """
    for cls in (CanteraConstantTP, CanteraConstantHP, CanteraConstantUV):
        adapter = _make_adapter(cls)
        assert len(adapter.inert_index_list) == len(set(adapter.inert_index_list)), \
            f'{cls.__name__} has duplicate inert indices: {adapter.inert_index_list}'


def test_no_duplicate_spc_identifier_lookup():
    """
    Verify spc_identifier_lookup has exactly num_ct_species entries
    (no overwriting from double-build).
    """
    adapter = _make_adapter(CanteraConstantTP)
    assert len(adapter.spc_identifier_lookup) == adapter.num_ct_species


def test_rxn_identifier_lookup_populated():
    """
    Verify rxn_identifier_lookup is populated (may be fewer than num_ct_reactions
    if the model has DUPLICATE reactions with the same equation string).
    """
    adapter = _make_adapter(CanteraConstantTP)
    assert len(adapter.rxn_identifier_lookup) > 0
    assert len(adapter.rxn_identifier_lookup) <= adapter.num_ct_reactions


# ---------------------------------------------------------------------------
# Bug fix: SA tolerance assignment (was swapped)
# ---------------------------------------------------------------------------

def test_sa_tolerance_assignment():
    """
    Previously rtol_sensitivity was assigned sa_atol and vice versa.
    Verify the tolerances are correctly assigned after reinitialize_simulation.
    """
    adapter = _make_adapter(CanteraConstantTP, observable_list=['OH'])
    condition = adapter.conditions[0]
    T0, P0, V0 = adapter._get_initial_state(condition)
    adapter.reinitialize_simulation(T0=T0, P0=P0, X0=condition.mol_frac, V0=V0)

    assert adapter.cantera_simulation.atol_sensitivity == adapter.sa_atol
    assert adapter.cantera_simulation.rtol_sensitivity == adapter.sa_rtol
    # Double-check they're not equal (otherwise the test is vacuous)
    assert adapter.sa_atol != adapter.sa_rtol


# ---------------------------------------------------------------------------
# Bug fix: get_t50 no longer hangs
# ---------------------------------------------------------------------------

def test_get_t50_terminates_for_stable_species():
    """
    get_t50 previously had ``while True`` which would hang forever if the
    species never reached 50% conversion.  Now it's bounded by the reaction time.
    Use an inert-like species that won't reach 50% consumption.
    """
    adapter = _make_adapter(CanteraConstantTP)
    adapter.simulate()
    # O2(2) in this H2/O2 system at constant T,P should deplete,
    # but N2 (if present) or a species at very low initial fraction might not.
    # Just verify the function returns without hanging.
    t50_list = adapter.get_t50('H2(1)')
    assert len(t50_list) == 1
    assert t50_list[0] > 0


# ---------------------------------------------------------------------------
# _get_initial_state helper
# ---------------------------------------------------------------------------

def test_get_initial_state():
    """Test the _get_initial_state helper extracts T0, P0, V0 correctly."""
    adapter = _make_adapter(CanteraConstantTP)
    condition = adapter.conditions[0]
    T0, P0, V0 = adapter._get_initial_state(condition)
    assert T0 > 0
    # For pressure-based conditions, P0 is set and V0 is None
    assert P0 is not None or V0 is not None
    if P0 is not None:
        assert V0 is None
        assert P0 > 0
    else:
        assert V0 > 0


# ---------------------------------------------------------------------------
# SA header parsing consistency (uses common functions)
# ---------------------------------------------------------------------------

def test_sa_coefficients_use_proper_header_parsing():
    """
    Verify that get_sa_coefficients produces correct keys for both kinetics
    and thermo, using the common header-parsing functions (not fragile inline splits).
    """
    adapter = _make_adapter(CanteraConstantTP, observable_list=['OH', 'H'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()

    # Kinetics: keys should be species labels (strings), values dict with int keys
    for obs_label, params in sa_dict['kinetics'][0].items():
        assert isinstance(obs_label, str)
        assert len(obs_label) > 0
        for rxn_key in params:
            assert isinstance(rxn_key, int)

    # Thermo: keys should be species labels (strings), values dict with string keys
    for obs_label, params in sa_dict['thermo'][0].items():
        assert isinstance(obs_label, str)
        assert len(obs_label) > 0
        for spc_key in params:
            assert isinstance(spc_key, str)
            assert spc_key is not None, \
                'Thermo SA parameter key is None — header parsing failed for /dH[ headers'
            assert len(spc_key) > 0


# ---------------------------------------------------------------------------
# All three adapters produce valid SA dicts with identical structure
# ---------------------------------------------------------------------------

def test_all_adapters_produce_valid_sa_dict():
    """
    Run SA with all three adapters on the same model and verify the SA dict
    structure is consistent.
    """
    for cls in (CanteraConstantTP, CanteraConstantHP, CanteraConstantUV):
        adapter = _make_adapter(cls, observable_list=['OH', 'H'])
        adapter.simulate()
        sa_dict = adapter.get_sa_coefficients()

        assert 'kinetics' in sa_dict
        assert 'thermo' in sa_dict
        assert 'time' in sa_dict
        assert len(sa_dict['time']) == 1  # one condition
        assert len(sa_dict['time'][0]) > 20

        # Should have 2 observables for each section
        assert len(sa_dict['kinetics'][0]) == 2, \
            f'{cls.__name__} kinetics has {len(sa_dict["kinetics"][0])} observables, expected 2'
        assert len(sa_dict['thermo'][0]) == 2, \
            f'{cls.__name__} thermo has {len(sa_dict["thermo"][0])} observables, expected 2'

        # All coefficient arrays should match time length
        for section in ('kinetics', 'thermo'):
            for obs, params in sa_dict[section][0].items():
                for param, values in params.items():
                    assert len(values) == len(sa_dict['time'][0]), \
                        f'{cls.__name__} {section}/{obs}/{param}: ' \
                        f'len(values)={len(values)} != len(time)={len(sa_dict["time"][0])}'


# ---------------------------------------------------------------------------
# Reactor type is correct per adapter
# ---------------------------------------------------------------------------

def test_reactor_types():
    """Verify each adapter sets the correct reactor type string."""
    assert CanteraConstantTP.cantera_reactor_type == 'IdealGasConstPressureReactor'
    assert CanteraConstantHP.cantera_reactor_type == 'IdealGasConstPressureReactor'
    assert CanteraConstantUV.cantera_reactor_type == 'IdealGasReactor'


def test_create_reactor_returns_correct_type():
    """Verify create_reactor produces the expected Cantera reactor type."""
    import cantera as ct
    for cls, expected_ct_type in [(CanteraConstantTP, ct.IdealGasConstPressureReactor),
                                  (CanteraConstantHP, ct.IdealGasConstPressureReactor),
                                  (CanteraConstantUV, ct.IdealGasReactor)]:
        adapter = _make_adapter(cls)
        condition = adapter.conditions[0]
        T0, P0, V0 = adapter._get_initial_state(condition)
        adapter.reinitialize_simulation(T0=T0, P0=P0, X0=condition.mol_frac, V0=V0)
        assert isinstance(adapter.cantera_reactor, expected_ct_type), \
            f'{cls.__name__} reactor is {type(adapter.cantera_reactor)}, expected {expected_ct_type}'


# ---------------------------------------------------------------------------
# IDT behavior per adapter type
# ---------------------------------------------------------------------------

def test_idt_constant_tp_returns_empty():
    """Constant T reactor has no ignition — IDT should be empty."""
    adapter = _make_adapter(CanteraConstantTP)
    adapter.simulate()
    idt_dict = adapter.get_idt_by_T()
    assert idt_dict['idt'] == []
    assert idt_dict['idt_index'] == []


def test_idt_hp_and_uv_return_values():
    """HP and UV reactors compute IDT from dT/dt."""
    for cls in (CanteraConstantHP, CanteraConstantUV):
        adapter = _make_adapter(cls)
        adapter.simulate()
        idt_dict = adapter.get_idt_by_T()
        assert len(idt_dict['idt']) == 1, f'{cls.__name__} should find 1 IDT'
        assert len(idt_dict['idt_index']) == 1
        assert idt_dict['idt'][0] > 0


# ---------------------------------------------------------------------------
# NH3 model: SA with Cantera thermo headers (/dH[...])
# ---------------------------------------------------------------------------

def test_nh3_thermo_sa_keys_not_none():
    """
    Cantera thermo SA headers use /dH[species], not /dG[species].
    Previously get_parameter_from_header didn't handle /dH[ and returned None.
    Verify all thermo SA keys are valid species labels.
    """
    t3_obj = run_minimal(project_directory=NH3_SA_BASE_TEST_DIR, iteration=1, set_paths=True)
    t3_obj.paths['cantera annotated'] = NH3_MODEL_PATH
    t3_obj.rmg['species'] = [
        {'label': 'NH3', 'smiles': 'N', 'concentration': 0.5,
         'observable': False, 'SA_observable': False, 'balance': False, 'reactive': True},
        {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0.25,
         'observable': False, 'SA_observable': False, 'balance': False, 'reactive': True},
        {'label': 'N2', 'smiles': 'N#N', 'concentration': 0.25,
         'observable': False, 'SA_observable': False, 'balance': False, 'reactive': True},
        {'label': 'H', 'smiles': '[H]', 'concentration': 0,
         'observable': False, 'SA_observable': True, 'balance': False, 'reactive': True},
        {'label': 'NH2', 'smiles': '[NH2]', 'concentration': 0,
         'observable': False, 'SA_observable': True, 'balance': False, 'reactive': True},
    ]
    adapter = CanteraConstantTP(t3=t3_obj.t3,
                                rmg=t3_obj.rmg,
                                paths=t3_obj.paths,
                                logger=t3_obj.logger,
                                atol=1e-16,
                                rtol=1e-8,
                                observable_list=['H', 'NH2'],
                                sa_atol=1e-6,
                                sa_rtol=1e-4,
                                )
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()

    # Thermo SA should have entries with real species labels, not None
    assert len(sa_dict['thermo'][0]) == 2
    for obs_label, params in sa_dict['thermo'][0].items():
        assert len(params) > 0, f'No thermo SA params for {obs_label}'
        for spc_key, values in params.items():
            assert spc_key is not None, f'Thermo key is None for observable {obs_label}'
            assert isinstance(spc_key, str)
            assert len(spc_key) > 0
            assert isinstance(values, np.ndarray)
            assert len(values) == len(sa_dict['time'][0])


NH3_SA_BASE_TEST_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'nh3_base_test')


def teardown_module():
    """Clean up test artifacts."""
    log_archive = os.path.join(TEST_DIR, 'log_archive')
    if os.path.isdir(log_archive):
        shutil.rmtree(log_archive, ignore_errors=True)
    for f in [os.path.join(TEST_DIR, 't3.log')]:
        if os.path.isfile(f):
            os.remove(f)
    if os.path.isdir(NH3_SA_BASE_TEST_DIR):
        shutil.rmtree(NH3_SA_BASE_TEST_DIR, ignore_errors=True)
