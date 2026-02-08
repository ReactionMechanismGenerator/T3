#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_cantera_jsr module
"""

import os
import shutil

import cantera as ct
import numpy as np
import pytest

from t3.common import SIMULATE_TEST_DATA_BASE_PATH, TEST_DATA_BASE_PATH
from tests.common import run_minimal
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.cantera_base import CanteraBase
from t3.simulate.cantera_constant_tp import CanteraConstantTP
from t3.simulate.cantera_jsr import CanteraJSR
import t3.simulate.cantera_jsr as jsr_module


TEST_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'cantera_simulator_test')
JSR_TEST_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'jsr_test')
JSR_NH3_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'jsr_nh3_test')
NH3_MODEL_PATH = os.path.join(TEST_DATA_BASE_PATH, 'models', 'NH3.yaml')


def _make_jsr(observable_list=None, project_directory=None):
    """Helper to construct a CanteraJSR adapter using the minimal example model."""
    t3 = run_minimal(project_directory=project_directory or TEST_DIR)
    t3.set_paths()
    return CanteraJSR(t3=t3.t3,
                      rmg=t3.rmg,
                      paths=t3.paths,
                      logger=t3.logger,
                      atol=t3.rmg['model']['atol'],
                      rtol=t3.rmg['model']['rtol'],
                      observable_list=observable_list or list(),
                      sa_atol=t3.t3['sensitivity']['atol'],
                      sa_rtol=t3.t3['sensitivity']['rtol'],
                      )


def _make_nh3_jsr(observable_list=None):
    """Helper to construct a CanteraJSR adapter using the NH3 model."""
    t3 = run_minimal(project_directory=JSR_NH3_DIR, iteration=1, set_paths=True)
    t3.paths['cantera annotated'] = NH3_MODEL_PATH
    t3.rmg['species'] = [
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
    return CanteraJSR(t3=t3.t3,
                      rmg=t3.rmg,
                      paths=t3.paths,
                      logger=t3.logger,
                      atol=1e-16,
                      rtol=1e-8,
                      observable_list=observable_list or list(),
                      sa_atol=1e-6,
                      sa_rtol=1e-4,
                      )


# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

def test_module_constants_defaults():
    """Verify the module-level constants have expected default values."""
    assert jsr_module.VOLUME == 1e-4
    assert jsr_module.PRESSURE_COEFF == 0.01


# ---------------------------------------------------------------------------
# Inheritance, registration, class attributes
# ---------------------------------------------------------------------------

def test_inheritance():
    """CanteraJSR inherits from CanteraBase and SimulateAdapter."""
    assert issubclass(CanteraJSR, CanteraBase)
    assert issubclass(CanteraJSR, SimulateAdapter)


def test_cantera_reactor_type_attribute():
    """CanteraJSR has the correct cantera_reactor_type class attribute."""
    assert CanteraJSR.cantera_reactor_type == 'IdealGasReactor'


def test_jsr_registered():
    """CanteraJSR is registered in the factory."""
    from t3.simulate.factory import _registered_simulate_adapters
    assert 'CanteraJSR' in _registered_simulate_adapters
    assert _registered_simulate_adapters['CanteraJSR'] is CanteraJSR


def test_jsr_via_factory():
    """CanteraJSR can be instantiated through the simulate_factory."""
    from t3.simulate.factory import simulate_factory
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    adapter = simulate_factory(simulate_method='CanteraJSR',
                               t3=t3.t3,
                               rmg=t3.rmg,
                               paths=t3.paths,
                               logger=t3.logger,
                               atol=t3.rmg['model']['atol'],
                               rtol=t3.rmg['model']['rtol'],
                               observable_list=['OH'],
                               )
    assert isinstance(adapter, CanteraJSR)


def test_create_reactor_returns_ideal_gas_reactor():
    """create_reactor returns a bare ct.IdealGasReactor (used by inherited helpers)."""
    adapter = _make_jsr()
    reactor = adapter.create_reactor()
    assert isinstance(reactor, ct.IdealGasReactor)


def test_jsr_is_instance_of_base():
    """An instantiated JSR adapter is an instance of CanteraBase and SimulateAdapter."""
    adapter = _make_jsr()
    assert isinstance(adapter, CanteraBase)
    assert isinstance(adapter, SimulateAdapter)
    assert isinstance(adapter, CanteraJSR)


# ---------------------------------------------------------------------------
# Basic simulation
# ---------------------------------------------------------------------------

def test_simulate_no_sa():
    """Run JSR simulation without SA — populates all_data."""
    adapter = _make_jsr()
    adapter.simulate()
    assert len(adapter.all_data) == 1  # one condition
    time = adapter.all_data[0][0].data
    assert len(time) > 5  # JSR converges fast; modest threshold


def test_all_data_structure():
    """Verify all_data has the correct tuple structure after JSR simulation."""
    adapter = _make_jsr()
    adapter.simulate()
    assert len(adapter.all_data) == 1
    time_gd, condition_data, rxn_sa_data, thermo_sa_data = adapter.all_data[0]
    # time GenericData
    assert time_gd.label == 'Time'
    assert time_gd.units == 's'
    assert len(time_gd.data) > 5
    # condition_data: T, P, then one per species
    assert len(condition_data) == 2 + adapter.num_ct_species
    assert condition_data[0].label == 'Temperature'
    assert condition_data[0].units == 'K'
    assert condition_data[1].label == 'Pressure'
    assert condition_data[1].units == 'Pa'
    # No SA data when no observables
    assert rxn_sa_data == []
    assert thermo_sa_data == []


def test_time_monotonically_increasing():
    """Time data from JSR simulation should be strictly increasing."""
    adapter = _make_jsr()
    adapter.simulate()
    times = np.array(adapter.all_data[0][0].data)
    assert np.all(np.diff(times) > 0)


def test_time_starts_above_zero():
    """First recorded time should be > 0 (after the first solver step)."""
    adapter = _make_jsr()
    adapter.simulate()
    times = adapter.all_data[0][0].data
    assert times[0] > 0


def test_time_reaches_residence_time():
    """The integration loop runs until the network time reaches the residence time."""
    adapter = _make_jsr()
    residence_time = adapter.conditions[0].reaction_time.value_si
    adapter.simulate()
    times = adapter.all_data[0][0].data
    # The base-class loop terminates when sim.time >= reaction_time, so the
    # final recorded time must be at or beyond the residence time.
    assert times[-1] >= residence_time


def test_temperature_constant():
    """JSR is isothermal — temperature should not vary."""
    adapter = _make_jsr()
    adapter.simulate()
    T_data = np.array(adapter.all_data[0][1][0].data)
    np.testing.assert_allclose(T_data, T_data[0], atol=1e-6)


def test_pressure_approximately_constant():
    """JSR is isobaric — pressure should be approximately constant."""
    adapter = _make_jsr()
    adapter.simulate()
    P_data = np.array(adapter.all_data[0][1][1].data)
    np.testing.assert_allclose(P_data, P_data[0], rtol=1e-3)


def test_mole_fractions_sum_to_one():
    """Mole fractions should sum to ~1 at every time step."""
    adapter = _make_jsr()
    adapter.simulate()
    condition_data = adapter.all_data[0][1]
    n_spc = adapter.num_ct_species
    n_steps = len(condition_data[0].data)
    for i in range(n_steps):
        x_sum = sum(condition_data[2 + s].data[i] for s in range(n_spc))
        np.testing.assert_allclose(x_sum, 1.0, atol=1e-10)


def test_mole_fractions_nonnegative():
    """All mole fractions should be non-negative (within numerical tolerance)."""
    adapter = _make_jsr()
    adapter.simulate()
    condition_data = adapter.all_data[0][1]
    for s in range(adapter.num_ct_species):
        data = np.array(condition_data[2 + s].data)
        assert np.all(data >= -1e-15), \
            f'Negative mole fraction for {condition_data[2 + s].label}'


def test_species_labels_match_model():
    """Species GenericData labels should match Cantera model species names."""
    adapter = _make_jsr()
    adapter.simulate()
    condition_data = adapter.all_data[0][1]
    model_species = [s.name for s in adapter.model.species()]
    for s in range(adapter.num_ct_species):
        assert condition_data[2 + s].label == model_species[s]


# ---------------------------------------------------------------------------
# Flow network construction
# ---------------------------------------------------------------------------

def test_flow_network_objects_attached_after_simulate():
    """After simulate(), the JSR flow-network objects are attached to the adapter."""
    adapter = _make_jsr()
    adapter.simulate()
    assert isinstance(adapter._jsr_inlet, ct.Reservoir)
    assert isinstance(adapter._jsr_exhaust, ct.Reservoir)
    assert isinstance(adapter._jsr_mfc, ct.MassFlowController)
    assert isinstance(adapter._jsr_pc, ct.PressureController)
    assert isinstance(adapter.cantera_reactor, ct.IdealGasReactor)


def test_mass_flow_rate_matches_residence_time():
    """mdot should equal m_reactor(t=0) / tau (within roundoff).

    The MassFlowController is built with ``mdot = cantera_reactor.mass / tau``
    at construction time.  Cantera holds that value constant during integration
    even though the reactor mass evolves slightly.
    """
    adapter = _make_jsr()
    adapter.simulate()
    residence_time = adapter.conditions[0].reaction_time.value_si
    mdot = adapter._jsr_mfc.mass_flow_rate
    # The reactor mass equals (density * VOLUME) and the inlet reservoir is
    # constructed from the same Solution object, so its density is the inlet
    # state density.  At t=0 (and for the isothermal/isobaric JSR network,
    # at all times) the controller mdot = inlet_density * VOLUME / tau.
    expected_mdot = adapter._jsr_inlet.thermo.density * jsr_module.VOLUME / residence_time
    np.testing.assert_allclose(mdot, expected_mdot, rtol=1e-3)


def test_pressure_controller_primary_is_mfc():
    """The PressureController's primary should be the MassFlowController."""
    adapter = _make_jsr()
    adapter.simulate()
    # Primary device of the PC is the MFC; can check via mass_flow_rate equivalence
    assert adapter._jsr_pc.mass_flow_rate == pytest.approx(adapter._jsr_mfc.mass_flow_rate,
                                                            rel=1e-3, abs=1e-12)


# ---------------------------------------------------------------------------
# SA
# ---------------------------------------------------------------------------

def test_sa_two_observables():
    """JSR with two SA observables produces a correctly shaped sa_dict."""
    adapter = _make_jsr(observable_list=['OH', 'H'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    assert len(sa_dict['time']) == 1
    assert len(sa_dict['kinetics']) == 1
    assert len(sa_dict['thermo']) == 1
    assert len(sa_dict['kinetics'][0]) == 2
    assert len(sa_dict['thermo'][0]) == 2


def test_sa_single_observable():
    """SA with a single observable produces correct structure."""
    adapter = _make_jsr(observable_list=['OH'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    assert len(sa_dict['kinetics'][0]) == 1
    assert len(sa_dict['thermo'][0]) == 1
    obs = list(sa_dict['kinetics'][0].keys())[0]
    assert 'OH' in obs


def test_sa_coefficient_array_lengths():
    """All SA coefficient arrays must have the same length as the time array."""
    adapter = _make_jsr(observable_list=['OH', 'H'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    n_time = len(sa_dict['time'][0])
    for section in ('kinetics', 'thermo'):
        for obs, params in sa_dict[section][0].items():
            for param, values in params.items():
                assert len(values) == n_time, \
                    f'{section}/{obs}/{param}: len={len(values)} != n_time={n_time}'


def test_sa_kinetics_keys_are_valid_reaction_indices():
    """All kinetics SA reaction indices should be 1-based and within bounds."""
    adapter = _make_jsr(observable_list=['OH'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    for obs_label, params in sa_dict['kinetics'][0].items():
        for rxn_key in params:
            assert isinstance(rxn_key, int)
            assert 1 <= rxn_key <= adapter.num_ct_reactions, \
                f'Reaction index {rxn_key} out of range [1, {adapter.num_ct_reactions}]'


def test_sa_thermo_keys_are_valid_species():
    """All thermo SA keys should be valid species labels in the model."""
    adapter = _make_jsr(observable_list=['OH', 'H'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    model_species = {s.name for s in adapter.model.species()}
    for obs_label, params in sa_dict['thermo'][0].items():
        assert len(params) > 0
        for spc_key in params:
            assert spc_key is not None
            assert isinstance(spc_key, str)
            assert len(spc_key) > 0
            assert spc_key in model_species


def test_sa_has_nonzero_kinetics_coefficients():
    """At least some kinetics SA coefficients should be non-zero."""
    adapter = _make_jsr(observable_list=['OH', 'H'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    max_abs = 0.0
    for obs, params in sa_dict['kinetics'][0].items():
        for param, values in params.items():
            max_abs = max(max_abs, float(np.max(np.abs(values))))
    assert max_abs > 1e-12, 'All kinetics SA coefficients are zero'


def test_sa_all_data_raw_structure():
    """Verify all_data raw SA data lists are populated when SA is active."""
    adapter = _make_jsr(observable_list=['OH', 'H'])
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


def test_sa_kinetics_values_are_ndarray():
    """SA kinetics values should be numpy arrays."""
    adapter = _make_jsr(observable_list=['OH'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    for obs, params in sa_dict['kinetics'][0].items():
        for param, values in params.items():
            assert isinstance(values, np.ndarray)


# ---------------------------------------------------------------------------
# IDT and t50
# ---------------------------------------------------------------------------

def test_get_idt_returns_empty():
    """Steady-state isothermal JSR has no ignition delay."""
    adapter = _make_jsr()
    adapter.simulate()
    idt_dict = adapter.get_idt_by_T()
    assert idt_dict == {'idt': [], 'idt_index': []}


def test_get_idt_returns_empty_before_simulate():
    """get_idt_by_T should return empty even without running simulate first."""
    adapter = _make_jsr()
    idt_dict = adapter.get_idt_by_T()
    assert idt_dict == {'idt': [], 'idt_index': []}


def test_get_t50_raises_not_implemented():
    """get_t50 is not meaningful for a JSR and must raise NotImplementedError."""
    adapter = _make_jsr()
    adapter.simulate()
    with pytest.raises(NotImplementedError, match='JSR'):
        adapter.get_t50('H2(1)')


# ---------------------------------------------------------------------------
# find_equilibrium (inherited)
# ---------------------------------------------------------------------------

def test_find_equilibrium():
    """JSR adapter can find equilibrium (inherited from base, reactor wiring is unused)."""
    adapter = _make_jsr()
    equilibrium_dicts = adapter.find_equilibrium('TP')
    assert len(equilibrium_dicts) == 1
    assert len(equilibrium_dicts[0].keys()) == 8


def test_find_equilibrium_then_simulate():
    """find_equilibrium and simulate can be called in any order on the same adapter."""
    adapter = _make_jsr()
    adapter.find_equilibrium('TP')
    adapter.simulate()
    assert len(adapter.all_data) == 1
    assert len(adapter.all_data[0][0].data) > 5


def test_simulate_then_find_equilibrium():
    """simulate followed by find_equilibrium works."""
    adapter = _make_jsr()
    adapter.simulate()
    eq = adapter.find_equilibrium('TP')
    assert len(eq) == 1


def test_simulate_twice_in_a_row():
    """simulate can be called multiple times in a row (residence-time iterator resets)."""
    adapter = _make_jsr()
    adapter.simulate()
    n1 = len(adapter.all_data[0][0].data)
    adapter.simulate()
    n2 = len(adapter.all_data[0][0].data)
    # Same input should give the same number of steps
    assert n1 == n2


# ---------------------------------------------------------------------------
# Module constant override
# ---------------------------------------------------------------------------

def test_volume_override_runs_without_error():
    """Overriding the VOLUME constant should produce a valid simulation."""
    original_volume = jsr_module.VOLUME
    try:
        jsr_module.VOLUME = 1e-3   # 1 L (10x default)
        adapter = _make_jsr()
        adapter.simulate()
        assert len(adapter.all_data[0][0].data) > 5
    finally:
        jsr_module.VOLUME = original_volume


def test_volume_does_not_change_steady_state_composition():
    """
    For a fixed inlet state and residence time, the JSR steady-state composition
    is independent of reactor volume (volume only scales the absolute mass flow
    rate, not mdot/m_reactor = 1/tau).  Verify this physical invariant.
    """
    original_volume = jsr_module.VOLUME
    try:
        # Run with default volume
        adapter1 = _make_jsr()
        adapter1.simulate()
        final1 = {adapter1.all_data[0][1][2 + s].label: adapter1.all_data[0][1][2 + s].data[-1]
                  for s in range(adapter1.num_ct_species)}

        # Run with 10x volume
        jsr_module.VOLUME = original_volume * 10
        adapter2 = _make_jsr()
        adapter2.simulate()
        final2 = {adapter2.all_data[0][1][2 + s].label: adapter2.all_data[0][1][2 + s].data[-1]
                  for s in range(adapter2.num_ct_species)}

        # Two integrator runs at different volumes pick slightly different
        # internal step sizes, so use a generous absolute tolerance.  The
        # underlying ODEs (mdot/m_reactor = 1/tau, chemistry term independent
        # of volume) are exactly volume-invariant — any drift here is solver
        # noise, not physics.
        for spc in final1:
            np.testing.assert_allclose(final2[spc], final1[spc], rtol=1e-2, atol=1e-9,
                                       err_msg=f'Volume changed steady state for {spc}')
    finally:
        jsr_module.VOLUME = original_volume


# ---------------------------------------------------------------------------
# Physical sanity: short tau ≈ inlet, long tau ≈ batch result
# ---------------------------------------------------------------------------

def test_jsr_long_residence_time_approaches_batch_final_state():
    """
    With a long residence time, the JSR steady-state composition should
    approach the final composition of an equivalent batch reactor (because
    the gas spends enough time in the JSR for the chemistry to fully proceed).
    """
    # Batch reactor at the same T, P, X_in for the same physical time
    batch = CanteraConstantTP(
        t3=run_minimal(project_directory=TEST_DIR).t3,
        rmg=run_minimal(project_directory=TEST_DIR).rmg,
        paths={'cantera annotated': _make_jsr().paths['cantera annotated']},
        logger=_make_jsr().logger,
        atol=1e-16, rtol=1e-8,
        observable_list=list(),
    )
    batch.simulate()
    batch_final = {batch.all_data[0][1][2 + s].label: batch.all_data[0][1][2 + s].data[-1]
                   for s in range(batch.num_ct_species)}

    # JSR at the same conditions, with the residence time = batch reaction time
    jsr = _make_jsr()
    jsr.simulate()
    jsr_final = {jsr.all_data[0][1][2 + s].label: jsr.all_data[0][1][2 + s].data[-1]
                 for s in range(jsr.num_ct_species)}

    # For non-trace species, JSR steady state should be close to batch final
    # (qualitative sanity — both reach the same near-equilibrium state)
    for spc in batch_final:
        if batch_final[spc] > 1e-3:
            np.testing.assert_allclose(jsr_final[spc], batch_final[spc], rtol=0.20,
                                       err_msg=f'JSR ≠ batch for {spc}')


# ---------------------------------------------------------------------------
# NH3 model — end-to-end with SA
# ---------------------------------------------------------------------------

def test_jsr_nh3_no_sa():
    """Run JSR on the NH3 model without SA."""
    adapter = _make_nh3_jsr()
    adapter.simulate()
    assert len(adapter.all_data) == 1
    assert len(adapter.all_data[0][0].data) > 5
    # Mole fractions sum to 1
    n_spc = adapter.num_ct_species
    condition_data = adapter.all_data[0][1]
    for i in range(len(condition_data[0].data)):
        x_sum = sum(condition_data[2 + s].data[i] for s in range(n_spc))
        np.testing.assert_allclose(x_sum, 1.0, atol=1e-10)


def test_jsr_nh3_with_sa():
    """Run JSR on the NH3 model with H/NH2 as SA observables."""
    adapter = _make_nh3_jsr(observable_list=['H', 'NH2'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    assert sa_dict is not None
    assert 'H(24)' in sa_dict['kinetics'][0]
    assert 'NH2(6)' in sa_dict['kinetics'][0]
    assert len(sa_dict['time'][0]) > 5
    # Each kinetics observable has at least one reaction with non-trivial coefficient
    for obs in ['H(24)', 'NH2(6)']:
        params = sa_dict['kinetics'][0][obs]
        assert len(params) > 0
        for param, values in params.items():
            assert isinstance(param, int)
            assert isinstance(values, np.ndarray)
            assert len(values) == len(sa_dict['time'][0])


def test_jsr_nh3_sa_thermo_keys_match_model():
    """NH3 JSR SA thermo keys should be a subset of NH3 model species."""
    adapter = _make_nh3_jsr(observable_list=['H'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    model_species = {s.name for s in adapter.model.species()}
    for obs, params in sa_dict['thermo'][0].items():
        for spc_key in params:
            assert spc_key in model_species


# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------

def teardown_module():
    """Remove project directories created during these unit tests."""
    for d in (JSR_TEST_DIR, JSR_NH3_DIR):
        if os.path.isdir(d):
            shutil.rmtree(d, ignore_errors=True)
    log_archive = os.path.join(TEST_DIR, 'log_archive')
    if os.path.isdir(log_archive):
        shutil.rmtree(log_archive, ignore_errors=True)
    for f in (os.path.join(TEST_DIR, 't3.log'),):
        if os.path.isfile(f):
            os.remove(f)
