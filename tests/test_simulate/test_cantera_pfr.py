#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_cantera_pfr module
"""

import os
import shutil
from unittest import mock

import cantera as ct
import numpy as np

from t3.common import SIMULATE_TEST_DATA_BASE_PATH, TEST_DATA_BASE_PATH
from tests.common import run_minimal
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.cantera_base import CanteraBase
from t3.simulate.cantera_constant_tp import CanteraConstantTP
from t3.simulate.cantera_pfr import CanteraPFR
import t3.simulate.cantera_pfr as pfr_module


TEST_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'cantera_simulator_test')
NH3_MODEL_PATH = os.path.join(TEST_DATA_BASE_PATH, 'models', 'NH3.yaml')
PFR_TEST_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'pfr_test')
PFR_NH3_CHAIN_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'pfr_nh3_chain_test')


def _make_pfr(observable_list=None):
    """Helper to construct a CanteraPFR adapter using the minimal example model."""
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    return CanteraPFR(t3=t3.t3,
                      rmg=t3.rmg,
                      paths=t3.paths,
                      logger=t3.logger,
                      atol=t3.rmg['model']['atol'],
                      rtol=t3.rmg['model']['rtol'],
                      observable_list=observable_list or list(),
                      sa_atol=t3.t3['sensitivity']['atol'],
                      sa_rtol=t3.t3['sensitivity']['rtol'],
                      )


# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

def test_module_constants_defaults():
    """Verify the module-level constants have expected default values."""
    assert pfr_module.METHOD == 'lagrangian'
    assert pfr_module.N_CELLS == 100
    assert pfr_module.AREA == 1e-4
    assert pfr_module.LENGTH == 1.0


# ---------------------------------------------------------------------------
# Inheritance, registration, class attributes
# ---------------------------------------------------------------------------

def test_inheritance():
    """CanteraPFR inherits from CanteraBase and SimulateAdapter."""
    assert issubclass(CanteraPFR, CanteraBase)
    assert issubclass(CanteraPFR, SimulateAdapter)


def test_cantera_reactor_type_attribute():
    """CanteraPFR has the correct cantera_reactor_type class attribute."""
    assert CanteraPFR.cantera_reactor_type == 'IdealGasConstPressureReactor'


def test_pfr_registered():
    """CanteraPFR is registered in the factory."""
    from t3.simulate.factory import _registered_simulate_adapters
    assert 'CanteraPFR' in _registered_simulate_adapters
    assert _registered_simulate_adapters['CanteraPFR'] is CanteraPFR


def test_pfr_via_factory():
    """CanteraPFR can be instantiated through the simulate_factory."""
    from t3.simulate.factory import simulate_factory
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    adapter = simulate_factory(simulate_method='CanteraPFR',
                               t3=t3.t3,
                               rmg=t3.rmg,
                               paths=t3.paths,
                               logger=t3.logger,
                               atol=t3.rmg['model']['atol'],
                               rtol=t3.rmg['model']['rtol'],
                               observable_list=['OH'],
                               )
    assert isinstance(adapter, CanteraPFR)


def test_create_reactor_returns_correct_type():
    """create_reactor returns a ct.IdealGasConstPressureReactor with energy='off'."""
    adapter = _make_pfr()
    condition = adapter.conditions[0]
    T0, P0, V0 = adapter._get_initial_state(condition)
    adapter.reinitialize_simulation(T0=T0, P0=P0, X0=condition.mol_frac, V0=V0)
    assert isinstance(adapter.cantera_reactor, ct.IdealGasConstPressureReactor)


def test_pfr_is_instance_of_base():
    """An instantiated PFR adapter is an instance of CanteraBase and SimulateAdapter."""
    adapter = _make_pfr()
    assert isinstance(adapter, CanteraBase)
    assert isinstance(adapter, SimulateAdapter)
    assert isinstance(adapter, CanteraPFR)


# ---------------------------------------------------------------------------
# Lagrangian method — basic simulation
# ---------------------------------------------------------------------------

def test_lagrangian_simulate_no_sa():
    """Run PFR Lagrangian simulation without SA."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    time = adapter.all_data[0][0].data
    assert len(time) > 20


def test_lagrangian_all_data_structure():
    """Verify all_data has the correct tuple structure after Lagrangian simulation."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    assert len(adapter.all_data) == 1  # one condition
    time_gd, condition_data, rxn_sa_data, thermo_sa_data = adapter.all_data[0]
    # time GenericData
    assert time_gd.label == 'Time'
    assert len(time_gd.data) > 20
    # condition_data: T, P, then one per species
    assert len(condition_data) == 2 + adapter.num_ct_species
    assert condition_data[0].label == 'Temperature'
    assert condition_data[1].label == 'Pressure'
    # No SA data when no observables
    assert rxn_sa_data == []
    assert thermo_sa_data == []


def test_lagrangian_time_monotonically_increasing():
    """Time data from Lagrangian simulation should be strictly increasing."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    times = np.array(adapter.all_data[0][0].data)
    assert np.all(np.diff(times) > 0)


def test_lagrangian_species_profiles():
    """Verify Lagrangian species profiles are physically reasonable."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    condition_data = adapter.all_data[0][1]
    # Temperature should be constant (isothermal)
    T_data = np.array(condition_data[0].data)
    np.testing.assert_allclose(T_data, T_data[0], atol=1e-6)
    # Pressure should be constant (isobaric)
    P_data = np.array(condition_data[1].data)
    np.testing.assert_allclose(P_data, P_data[0], rtol=1e-3)
    # Mole fractions should sum to ~1 at every step
    n_spc = adapter.num_ct_species
    for i in range(len(T_data)):
        x_sum = sum(condition_data[2 + s].data[i] for s in range(n_spc))
        np.testing.assert_allclose(x_sum, 1.0, atol=1e-10)
    # All mole fractions should be non-negative
    for s in range(n_spc):
        assert np.all(np.array(condition_data[2 + s].data) >= -1e-15), \
            f'Negative mole fraction for {condition_data[2 + s].label}'


def test_lagrangian_species_labels_match_model():
    """Species GenericData labels should match Cantera model species names."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    condition_data = adapter.all_data[0][1]
    model_species = [s.name for s in adapter.model.species()]
    for s in range(adapter.num_ct_species):
        assert condition_data[2 + s].label == model_species[s]


# ---------------------------------------------------------------------------
# Lagrangian method — distance computation
# ---------------------------------------------------------------------------

def test_lagrangian_distance_computed():
    """Verify that distance_data is populated after Lagrangian simulation."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    assert hasattr(adapter, 'distance_data')
    assert len(adapter.distance_data) == 1  # one condition
    distance = adapter.distance_data[0]
    assert len(distance) == len(adapter.all_data[0][0].data)
    # Distance should be monotonically increasing
    assert np.all(np.diff(distance) > 0)
    # Final distance should be approximately LENGTH
    np.testing.assert_allclose(distance[-1], pfr_module.LENGTH, rtol=0.05)


def test_lagrangian_distance_starts_near_zero():
    """The first distance point should be very close to zero."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    distance = adapter.distance_data[0]
    # First point = u * dt where dt is the first time step — should be small
    assert distance[0] < 0.01 * pfr_module.LENGTH


def test_lagrangian_distance_all_positive():
    """All distance values should be strictly positive."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    distance = adapter.distance_data[0]
    assert np.all(distance > 0)


def test_lagrangian_distance_with_custom_length():
    """Distance scales with the LENGTH constant."""
    adapter = _make_pfr()
    custom_length = 5.0
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'), \
         mock.patch.object(pfr_module, 'LENGTH', custom_length):
        adapter.simulate()
    distance = adapter.distance_data[0]
    np.testing.assert_allclose(distance[-1], custom_length, rtol=0.05)


# ---------------------------------------------------------------------------
# Lagrangian method — SA
# ---------------------------------------------------------------------------

def test_lagrangian_sa():
    """Run PFR Lagrangian simulation with SA and verify the sa_dict structure."""
    adapter = _make_pfr(observable_list=['OH', 'H'])
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    assert len(sa_dict['time'][0]) > 20
    assert len(sa_dict['kinetics'][0]) == 2
    assert len(sa_dict['thermo'][0]) == 2
    for obs_label, params in sa_dict['kinetics'][0].items():
        for rxn_key, values in params.items():
            assert isinstance(rxn_key, int)
            assert isinstance(values, np.ndarray)
    for obs_label, params in sa_dict['thermo'][0].items():
        for spc_key, values in params.items():
            assert isinstance(spc_key, str)
            assert spc_key is not None


def test_lagrangian_sa_single_observable():
    """SA with a single observable produces correct structure."""
    adapter = _make_pfr(observable_list=['OH'])
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    assert len(sa_dict['kinetics'][0]) == 1
    assert len(sa_dict['thermo'][0]) == 1
    obs = list(sa_dict['kinetics'][0].keys())[0]
    assert 'OH' in obs


def test_lagrangian_sa_coefficient_array_lengths():
    """All SA coefficient arrays must have the same length as the time array."""
    adapter = _make_pfr(observable_list=['OH', 'H'])
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    n_time = len(sa_dict['time'][0])
    for section in ('kinetics', 'thermo'):
        for obs, params in sa_dict[section][0].items():
            for param, values in params.items():
                assert len(values) == n_time, \
                    f'{section}/{obs}/{param}: len={len(values)} != n_time={n_time}'


def test_lagrangian_sa_has_nonzero_kinetics_coefficients():
    """At least some kinetics SA coefficients should be non-zero (reactions matter)."""
    adapter = _make_pfr(observable_list=['OH'])
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    max_abs = 0.0
    for obs, params in sa_dict['kinetics'][0].items():
        for param, values in params.items():
            max_abs = max(max_abs, np.max(np.abs(values)))
    assert max_abs > 1e-10, 'All kinetics SA coefficients are zero'


def test_lagrangian_sa_has_nonzero_thermo_coefficients():
    """At least some thermo SA coefficients should be non-zero (thermodynamics matter)."""
    adapter = _make_pfr(observable_list=['OH'])
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    max_abs = 0.0
    for obs, params in sa_dict['thermo'][0].items():
        for param, values in params.items():
            max_abs = max(max_abs, np.max(np.abs(values)))
    assert max_abs > 1e-10, 'All thermo SA coefficients are zero'


def test_lagrangian_sa_thermo_keys_are_valid_species():
    """All thermo SA keys should be valid species labels (not None, not empty)."""
    adapter = _make_pfr(observable_list=['OH', 'H'])
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    model_species = {s.name for s in adapter.model.species()}
    for obs_label, params in sa_dict['thermo'][0].items():
        assert len(params) > 0, f'No thermo SA params for {obs_label}'
        for spc_key in params:
            assert spc_key is not None, f'Thermo key is None for observable {obs_label}'
            assert isinstance(spc_key, str)
            assert len(spc_key) > 0
            assert spc_key in model_species, \
                f'Thermo SA key {spc_key!r} not in model species'


def test_lagrangian_sa_kinetics_keys_are_valid_reactions():
    """All kinetics SA reaction indices should be valid (1-based, within num_reactions)."""
    adapter = _make_pfr(observable_list=['OH'])
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    for obs_label, params in sa_dict['kinetics'][0].items():
        for rxn_key in params:
            assert isinstance(rxn_key, int)
            assert 1 <= rxn_key <= adapter.num_ct_reactions, \
                f'Reaction index {rxn_key} out of range [1, {adapter.num_ct_reactions}]'


def test_lagrangian_sa_all_data_raw_structure():
    """Verify all_data raw SA data lists are populated when SA is active."""
    adapter = _make_pfr(observable_list=['OH', 'H'])
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    _, _, rxn_sa_data, thermo_sa_data = adapter.all_data[0]
    # Should have entries: n_observables * n_reactions for kinetics
    assert len(rxn_sa_data) == 2 * adapter.num_ct_reactions
    # Should have entries: n_observables * n_species for thermo
    assert len(thermo_sa_data) == 2 * adapter.num_ct_species
    # Each entry should be a GenericData with correct data length
    n_time = len(adapter.all_data[0][0].data)
    for gd in rxn_sa_data:
        assert len(gd.data) == n_time
    for gd in thermo_sa_data:
        assert len(gd.data) == n_time


def test_lagrangian_matches_constant_tp():
    """
    An isothermal isobaric PFR is mathematically equivalent to a constant-TP
    batch reactor.  Verify that both produce identical SA coefficients.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    common_kwargs = dict(t3=t3.t3, rmg=t3.rmg, paths=t3.paths, logger=t3.logger,
                         atol=t3.rmg['model']['atol'], rtol=t3.rmg['model']['rtol'],
                         observable_list=['OH', 'H'],
                         sa_atol=t3.t3['sensitivity']['atol'],
                         sa_rtol=t3.t3['sensitivity']['rtol'])

    tp = CanteraConstantTP(**common_kwargs)
    tp.simulate()
    sa_tp = tp.get_sa_coefficients()

    pfr = CanteraPFR(**common_kwargs)
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        pfr.simulate()
    sa_pfr = pfr.get_sa_coefficients()

    np.testing.assert_array_equal(sa_tp['time'][0], sa_pfr['time'][0])
    for obs in sa_tp['kinetics'][0]:
        for param in sa_tp['kinetics'][0][obs]:
            np.testing.assert_array_almost_equal(
                sa_tp['kinetics'][0][obs][param], sa_pfr['kinetics'][0][obs][param])
    for obs in sa_tp['thermo'][0]:
        for param in sa_tp['thermo'][0][obs]:
            np.testing.assert_array_almost_equal(
                sa_tp['thermo'][0][obs][param], sa_pfr['thermo'][0][obs][param])


def test_lagrangian_matches_constant_tp_species_profiles():
    """
    Lagrangian PFR and constant-TP batch reactor should produce identical
    species profiles (not just SA).
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    common_kwargs = dict(t3=t3.t3, rmg=t3.rmg, paths=t3.paths, logger=t3.logger,
                         atol=t3.rmg['model']['atol'], rtol=t3.rmg['model']['rtol'],
                         observable_list=[],
                         sa_atol=t3.t3['sensitivity']['atol'],
                         sa_rtol=t3.t3['sensitivity']['rtol'])

    tp = CanteraConstantTP(**common_kwargs)
    tp.simulate()
    tp_data = tp.all_data[0][1]

    pfr = CanteraPFR(**common_kwargs)
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        pfr.simulate()
    pfr_data = pfr.all_data[0][1]

    # Same number of data points
    assert len(tp_data[0].data) == len(pfr_data[0].data)
    # Same species mole fractions
    for s in range(tp.num_ct_species):
        np.testing.assert_array_almost_equal(
            tp_data[2 + s].data, pfr_data[2 + s].data,
            err_msg=f'Mismatch for {tp_data[2 + s].label}')


# ---------------------------------------------------------------------------
# Chain-of-reactors method — basic simulation
# ---------------------------------------------------------------------------

def test_chain_simulate_no_sa():
    """Run PFR chain-of-reactors simulation without SA."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'chain'), \
         mock.patch.object(pfr_module, 'N_CELLS', 50):
        adapter.simulate()
    time = adapter.all_data[0][0].data
    assert len(time) == 50  # one data point per cell


def test_chain_all_data_structure():
    """Verify all_data has the correct tuple structure after chain simulation."""
    n_cells = 30
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'chain'), \
         mock.patch.object(pfr_module, 'N_CELLS', n_cells):
        adapter.simulate()
    assert len(adapter.all_data) == 1
    time_gd, condition_data, rxn_sa_data, thermo_sa_data = adapter.all_data[0]
    assert time_gd.label == 'Time'
    assert len(time_gd.data) == n_cells
    assert len(condition_data) == 2 + adapter.num_ct_species
    assert condition_data[0].label == 'Temperature'
    assert condition_data[1].label == 'Pressure'
    # Chain produces no SA data
    assert rxn_sa_data == []
    assert thermo_sa_data == []


def test_chain_time_monotonically_increasing():
    """Cumulative residence time from chain should be strictly increasing."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'chain'), \
         mock.patch.object(pfr_module, 'N_CELLS', 50):
        adapter.simulate()
    times = np.array(adapter.all_data[0][0].data)
    assert np.all(np.diff(times) > 0)


def test_chain_distance_computed():
    """Verify that distance_data is populated after chain simulation."""
    n_cells = 50
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'chain'), \
         mock.patch.object(pfr_module, 'N_CELLS', n_cells):
        adapter.simulate()
    assert len(adapter.distance_data) == 1
    distance = adapter.distance_data[0]
    assert len(distance) == n_cells
    # Distance should be evenly spaced
    dz = pfr_module.LENGTH / n_cells
    expected = np.arange(1, n_cells + 1) * dz
    np.testing.assert_array_almost_equal(distance, expected)


def test_chain_distance_final_equals_length():
    """The last distance in the chain should exactly equal LENGTH."""
    n_cells = 40
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'chain'), \
         mock.patch.object(pfr_module, 'N_CELLS', n_cells):
        adapter.simulate()
    distance = adapter.distance_data[0]
    np.testing.assert_allclose(distance[-1], pfr_module.LENGTH, rtol=1e-10)


def test_chain_species_profiles():
    """Verify species profiles from the chain method are physically reasonable."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'chain'), \
         mock.patch.object(pfr_module, 'N_CELLS', 50):
        adapter.simulate()
    # condition_data: [T, P, species...]
    condition_data = adapter.all_data[0][1]
    # Temperature should be constant (isothermal)
    T_data = condition_data[0].data
    assert all(abs(T - T_data[0]) < 1e-6 for T in T_data)
    # Pressure should be constant (isobaric)
    P_data = condition_data[1].data
    assert all(abs(P - P_data[0]) / P_data[0] < 1e-3 for P in P_data)
    # Mole fractions should sum to ~1 at every cell
    n_spc = adapter.num_ct_species
    for i in range(len(T_data)):
        x_sum = sum(condition_data[2 + s].data[i] for s in range(n_spc))
        np.testing.assert_allclose(x_sum, 1.0, atol=1e-10)


def test_chain_species_labels_match_model():
    """Species GenericData labels from chain should match Cantera model species names."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'chain'), \
         mock.patch.object(pfr_module, 'N_CELLS', 20):
        adapter.simulate()
    condition_data = adapter.all_data[0][1]
    model_species = [s.name for s in adapter.model.species()]
    for s in range(adapter.num_ct_species):
        assert condition_data[2 + s].label == model_species[s]


def test_chain_mole_fractions_nonnegative():
    """All mole fractions from chain should be non-negative."""
    adapter = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'chain'), \
         mock.patch.object(pfr_module, 'N_CELLS', 50):
        adapter.simulate()
    condition_data = adapter.all_data[0][1]
    for s in range(adapter.num_ct_species):
        data = np.array(condition_data[2 + s].data)
        assert np.all(data >= -1e-15), \
            f'Negative mole fraction for {condition_data[2 + s].label}'


# ---------------------------------------------------------------------------
# Chain-of-reactors method — SA fallback
# ---------------------------------------------------------------------------

def test_chain_sa_fallback_to_lagrangian():
    """When SA is requested with chain method, adapter falls back to Lagrangian."""
    adapter = _make_pfr(observable_list=['OH', 'H'])
    with mock.patch.object(pfr_module, 'METHOD', 'chain'):
        adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    # SA should have been produced (via Lagrangian fallback)
    assert len(sa_dict['kinetics'][0]) == 2
    assert len(sa_dict['thermo'][0]) == 2
    # Distance should also be computed (Lagrangian computes it)
    assert hasattr(adapter, 'distance_data')
    assert len(adapter.distance_data) == 1


def test_chain_sa_fallback_matches_direct_lagrangian():
    """
    SA produced via chain fallback should be identical to direct Lagrangian SA,
    since the fallback literally runs the Lagrangian method.
    """
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    common_kwargs = dict(t3=t3.t3, rmg=t3.rmg, paths=t3.paths, logger=t3.logger,
                         atol=t3.rmg['model']['atol'], rtol=t3.rmg['model']['rtol'],
                         observable_list=['OH'],
                         sa_atol=t3.t3['sensitivity']['atol'],
                         sa_rtol=t3.t3['sensitivity']['rtol'])

    direct = CanteraPFR(**common_kwargs)
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        direct.simulate()
    sa_direct = direct.get_sa_coefficients()

    fallback = CanteraPFR(**common_kwargs)
    with mock.patch.object(pfr_module, 'METHOD', 'chain'):
        fallback.simulate()
    sa_fallback = fallback.get_sa_coefficients()

    np.testing.assert_array_equal(sa_direct['time'][0], sa_fallback['time'][0])
    for obs in sa_direct['kinetics'][0]:
        for param in sa_direct['kinetics'][0][obs]:
            np.testing.assert_array_equal(
                sa_direct['kinetics'][0][obs][param], sa_fallback['kinetics'][0][obs][param])
    for obs in sa_direct['thermo'][0]:
        for param in sa_direct['thermo'][0][obs]:
            np.testing.assert_array_equal(
                sa_direct['thermo'][0][obs][param], sa_fallback['thermo'][0][obs][param])


def test_chain_sa_fallback_also_computes_distance():
    """When chain falls back to Lagrangian for SA, distance should still be computed."""
    adapter = _make_pfr(observable_list=['OH'])
    with mock.patch.object(pfr_module, 'METHOD', 'chain'):
        adapter.simulate()
    assert hasattr(adapter, 'distance_data')
    assert len(adapter.distance_data) == 1
    distance = adapter.distance_data[0]
    assert len(distance) == len(adapter.all_data[0][0].data)
    assert np.all(np.diff(distance) > 0)


# ---------------------------------------------------------------------------
# Chain-of-reactors method — convergence & residence time
# ---------------------------------------------------------------------------

def test_chain_cumulative_residence_time():
    """
    Verify that the cumulative residence time at the last cell approximately
    equals the total residence time specified in the input.
    """
    adapter = _make_pfr()
    total_time = adapter.conditions[0].reaction_time.value_si
    with mock.patch.object(pfr_module, 'METHOD', 'chain'), \
         mock.patch.object(pfr_module, 'N_CELLS', 100):
        adapter.simulate()
    times = adapter.all_data[0][0].data
    # Final cumulative residence time should be close to total_time
    np.testing.assert_allclose(times[-1], total_time, rtol=0.1)


def test_chain_converges_to_lagrangian():
    """
    With enough cells, the chain method should converge to the Lagrangian solution.
    Compare final species mole fractions.
    """
    adapter_lag = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter_lag.simulate()
    # Get final mole fractions from Lagrangian
    lag_data = adapter_lag.all_data[0][1]
    lag_final = {lag_data[2 + s].label: lag_data[2 + s].data[-1]
                 for s in range(adapter_lag.num_ct_species)}

    adapter_chain = _make_pfr()
    with mock.patch.object(pfr_module, 'METHOD', 'chain'), \
         mock.patch.object(pfr_module, 'N_CELLS', 200):
        adapter_chain.simulate()
    chain_data = adapter_chain.all_data[0][1]
    chain_final = {chain_data[2 + s].label: chain_data[2 + s].data[-1]
                   for s in range(adapter_chain.num_ct_species)}

    for spc, x_lag in lag_final.items():
        x_chain = chain_final[spc]
        if x_lag > 1e-6:  # only compare non-trace species
            np.testing.assert_allclose(x_chain, x_lag, rtol=0.05,
                                       err_msg=f'Chain/Lagrangian mismatch for {spc}')


def test_chain_with_different_n_cells():
    """Chain simulation works with various N_CELLS values and always produces correct output length."""
    for n_cells in [10, 50]:
        adapter = _make_pfr()
        with mock.patch.object(pfr_module, 'METHOD', 'chain'), \
             mock.patch.object(pfr_module, 'N_CELLS', n_cells):
            adapter.simulate()
        assert len(adapter.all_data[0][0].data) == n_cells
        assert len(adapter.distance_data[0]) == n_cells
        # Mole fracs still sum to 1
        condition_data = adapter.all_data[0][1]
        for i in range(n_cells):
            x_sum = sum(condition_data[2 + s].data[i] for s in range(adapter.num_ct_species))
            np.testing.assert_allclose(x_sum, 1.0, atol=1e-10)


# ---------------------------------------------------------------------------
# IDT and other inherited methods
# ---------------------------------------------------------------------------

def test_get_idt_returns_empty():
    """Isothermal PFR has no ignition."""
    adapter = _make_pfr()
    adapter.simulate()
    idt_dict = adapter.get_idt_by_T()
    assert idt_dict == {'idt': [], 'idt_index': []}


def test_get_idt_returns_empty_before_simulate():
    """get_idt_by_T should return empty even without running simulate first."""
    adapter = _make_pfr()
    idt_dict = adapter.get_idt_by_T()
    assert idt_dict == {'idt': [], 'idt_index': []}


def test_find_equilibrium():
    """PFR adapter can find equilibrium (inherited from base)."""
    adapter = _make_pfr()
    equilibrium_dicts = adapter.find_equilibrium('TP')
    assert len(equilibrium_dicts) == 1
    assert len(equilibrium_dicts[0].keys()) == 8


def test_get_t50():
    """PFR adapter can compute t50 (residence time to 50% conversion)."""
    adapter = _make_pfr()
    adapter.simulate()
    t50_list = adapter.get_t50('H2(1)')
    assert len(t50_list) == 1
    assert t50_list[0] > 0


def test_get_t50_within_residence_time():
    """t50 should be approximately within the total residence time (may overshoot by one step)."""
    adapter = _make_pfr()
    adapter.simulate()
    total_time = adapter.conditions[0].reaction_time.value_si
    t50_list = adapter.get_t50('H2(1)')
    # The integrator may overshoot by one step, so allow a small margin
    assert t50_list[0] < total_time * 1.1


# ---------------------------------------------------------------------------
# NH3 model — Lagrangian
# ---------------------------------------------------------------------------

def test_pfr_nh3_lagrangian():
    """Run PFR Lagrangian SA on the NH3 model."""
    t3_obj = run_minimal(project_directory=PFR_TEST_DIR, iteration=1, set_paths=True)
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
    ]
    adapter = CanteraPFR(t3=t3_obj.t3,
                         rmg=t3_obj.rmg,
                         paths=t3_obj.paths,
                         logger=t3_obj.logger,
                         atol=1e-16,
                         rtol=1e-8,
                         observable_list=['H'],
                         sa_atol=1e-6,
                         sa_rtol=1e-4,
                         )
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    assert len(sa_dict['kinetics'][0]) == 1
    assert len(sa_dict['thermo'][0]) == 1
    # Distance should be computed
    assert len(adapter.distance_data) == 1
    assert len(adapter.distance_data[0]) == len(sa_dict['time'][0])


def test_pfr_nh3_lagrangian_thermo_keys_not_none():
    """
    NH3 model uses Cantera /dH[] headers.
    Verify all thermo SA keys are valid species labels (not None).
    """
    t3_obj = run_minimal(project_directory=PFR_TEST_DIR, iteration=1, set_paths=True)
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
    ]
    adapter = CanteraPFR(t3=t3_obj.t3,
                         rmg=t3_obj.rmg,
                         paths=t3_obj.paths,
                         logger=t3_obj.logger,
                         atol=1e-16,
                         rtol=1e-8,
                         observable_list=['H'],
                         sa_atol=1e-6,
                         sa_rtol=1e-4,
                         )
    with mock.patch.object(pfr_module, 'METHOD', 'lagrangian'):
        adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    for obs_label, params in sa_dict['thermo'][0].items():
        assert len(params) > 0
        for spc_key, values in params.items():
            assert spc_key is not None, f'Thermo key is None for observable {obs_label}'
            assert isinstance(spc_key, str)
            assert len(spc_key) > 0
            assert isinstance(values, np.ndarray)


# ---------------------------------------------------------------------------
# NH3 model — Chain
# ---------------------------------------------------------------------------

def test_pfr_nh3_chain():
    """Run PFR chain-of-reactors on the NH3 model (no SA)."""
    t3_obj = run_minimal(project_directory=PFR_NH3_CHAIN_DIR, iteration=1, set_paths=True)
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
    ]
    adapter = CanteraPFR(t3=t3_obj.t3,
                         rmg=t3_obj.rmg,
                         paths=t3_obj.paths,
                         logger=t3_obj.logger,
                         atol=1e-16,
                         rtol=1e-8,
                         observable_list=[],
                         sa_atol=1e-6,
                         sa_rtol=1e-4,
                         )
    n_cells = 50
    with mock.patch.object(pfr_module, 'METHOD', 'chain'), \
         mock.patch.object(pfr_module, 'N_CELLS', n_cells):
        adapter.simulate()
    assert len(adapter.all_data) == 1
    assert len(adapter.all_data[0][0].data) == n_cells
    assert len(adapter.distance_data) == 1
    assert len(adapter.distance_data[0]) == n_cells
    # Mole fractions should be physical
    condition_data = adapter.all_data[0][1]
    for i in range(n_cells):
        x_sum = sum(condition_data[2 + s].data[i] for s in range(adapter.num_ct_species))
        np.testing.assert_allclose(x_sum, 1.0, atol=1e-10)


def teardown_module():
    """Clean up test artifacts."""
    log_archive = os.path.join(TEST_DIR, 'log_archive')
    if os.path.isdir(log_archive):
        shutil.rmtree(log_archive, ignore_errors=True)
    for f in [os.path.join(TEST_DIR, 't3.log')]:
        if os.path.isfile(f):
            os.remove(f)
    for d in [PFR_TEST_DIR, PFR_NH3_CHAIN_DIR]:
        if os.path.isdir(d):
            shutil.rmtree(d, ignore_errors=True)
