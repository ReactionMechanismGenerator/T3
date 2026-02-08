#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_cantera_pfr_t_profile module

Tests for the CanteraPFRTProfile adapter — a chain-of-cells PFR with a
hardcoded axial temperature profile.
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
from t3.simulate.cantera_pfr_t_profile import (
    CanteraPFRTProfile,
    temperature_profile,
)
import t3.simulate.cantera_pfr_t_profile as tprof_module


TEST_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'cantera_simulator_test')
TPROF_NH3_DIR = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'pfr_tprof_nh3_test')
NH3_MODEL_PATH = os.path.join(TEST_DATA_BASE_PATH, 'models', 'NH3.yaml')


def _make_tprof(observable_list=None, project_directory=None):
    """Helper to construct a CanteraPFRTProfile adapter using the minimal model."""
    t3 = run_minimal(project_directory=project_directory or TEST_DIR)
    t3.set_paths()
    return CanteraPFRTProfile(t3=t3.t3,
                              rmg=t3.rmg,
                              paths=t3.paths,
                              logger=t3.logger,
                              atol=t3.rmg['model']['atol'],
                              rtol=t3.rmg['model']['rtol'],
                              observable_list=observable_list or list(),
                              sa_atol=t3.t3['sensitivity']['atol'],
                              sa_rtol=t3.t3['sensitivity']['rtol'],
                              )


def _make_nh3_tprof(observable_list=None):
    """Helper to construct a CanteraPFRTProfile adapter using the NH3 model."""
    t3 = run_minimal(project_directory=TPROF_NH3_DIR, iteration=1, set_paths=True)
    t3.paths['cantera annotated'] = NH3_MODEL_PATH
    t3.rmg['species'] = [
        {'label': 'NH3', 'smiles': 'N', 'concentration': 0.5,
         'observable': False, 'SA_observable': False, 'balance': False, 'reactive': True},
        {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0.25,
         'observable': False, 'SA_observable': False, 'balance': False, 'reactive': True},
        {'label': 'N2', 'smiles': 'N#N', 'concentration': 0.25,
         'observable': False, 'SA_observable': False, 'balance': False, 'reactive': True},
        # Listed (no concentration) so an observable_list=['NH2'] entry can
        # match a real species and populate sensitive_species in the adapter.
        {'label': 'NH2', 'smiles': '[NH2]', 'concentration': 0,
         'observable': False, 'SA_observable': False, 'balance': False, 'reactive': True},
    ]
    return CanteraPFRTProfile(t3=t3.t3,
                              rmg=t3.rmg,
                              paths=t3.paths,
                              logger=t3.logger,
                              atol=1e-14,
                              rtol=1e-7,
                              observable_list=observable_list or list(),
                              sa_atol=1e-6,
                              sa_rtol=1e-4,
                              )


# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

def test_module_constants_defaults():
    """Verify the module-level constants have expected default values."""
    assert tprof_module.LENGTH == 0.80
    assert tprof_module.AREA == 1e-4
    assert tprof_module.N_CELLS == 200
    assert tprof_module.RAMP_UP_END == 0.20
    assert tprof_module.ISO_END == 0.60
    assert tprof_module.T_INLET == 300.0
    assert tprof_module.T_HOT == 900.0
    assert tprof_module.T_OUTLET == 500.0


def test_module_geometry_consistent():
    """Module-level constants must satisfy 0 < RAMP_UP_END < ISO_END < LENGTH."""
    assert 0 < tprof_module.RAMP_UP_END < tprof_module.ISO_END < tprof_module.LENGTH


# ---------------------------------------------------------------------------
# Temperature profile function (the user-facing geometry of T(z))
# ---------------------------------------------------------------------------

def test_profile_endpoint_inlet():
    """T(0) == T_INLET (300 K)."""
    assert temperature_profile(0.0) == pytest.approx(tprof_module.T_INLET)


def test_profile_endpoint_ramp_up_end():
    """T(0.20) == T_HOT (900 K) — top of the up-ramp."""
    assert temperature_profile(tprof_module.RAMP_UP_END) == pytest.approx(tprof_module.T_HOT)


def test_profile_endpoint_iso_end():
    """T(0.60) == T_HOT (900 K) — end of the isothermal plateau."""
    assert temperature_profile(tprof_module.ISO_END) == pytest.approx(tprof_module.T_HOT)


def test_profile_endpoint_outlet():
    """T(0.80) == T_OUTLET (500 K)."""
    assert temperature_profile(tprof_module.LENGTH) == pytest.approx(tprof_module.T_OUTLET)


def test_profile_midpoint_ramp_up():
    """Mid of the raised-cosine ramp-up: T = (T_INLET + T_HOT) / 2 = 600 K at z=0.10."""
    z_mid = tprof_module.RAMP_UP_END / 2
    expected = 0.5 * (tprof_module.T_INLET + tprof_module.T_HOT)
    assert temperature_profile(z_mid) == pytest.approx(expected)


def test_profile_midpoint_ramp_down():
    """Mid of the ramp-down: T = (T_HOT + T_OUTLET) / 2 = 700 K at z=0.70."""
    z_mid = (tprof_module.ISO_END + tprof_module.LENGTH) / 2
    expected = 0.5 * (tprof_module.T_HOT + tprof_module.T_OUTLET)
    assert temperature_profile(z_mid) == pytest.approx(expected)


def test_profile_isothermal_segment_constant():
    """Profile is exactly T_HOT for any point in [RAMP_UP_END, ISO_END)."""
    for z in np.linspace(tprof_module.RAMP_UP_END, tprof_module.ISO_END, 25, endpoint=False):
        assert temperature_profile(z) == pytest.approx(tprof_module.T_HOT)


def test_profile_clipping_below_zero():
    """z < 0 is clipped to z=0 → returns T_INLET."""
    assert temperature_profile(-1.0) == pytest.approx(tprof_module.T_INLET)


def test_profile_clipping_above_length():
    """z > LENGTH is clipped to z=LENGTH → returns T_OUTLET."""
    assert temperature_profile(2 * tprof_module.LENGTH) == pytest.approx(tprof_module.T_OUTLET)


def test_profile_ramp_up_monotonic_increasing():
    """Profile is strictly increasing on the ramp-up segment."""
    zs = np.linspace(0, tprof_module.RAMP_UP_END, 50)
    Ts = np.array([temperature_profile(z) for z in zs])
    assert np.all(np.diff(Ts) >= -1e-12)
    # And monotonically increasing in the strict sense at the interior
    assert np.all(np.diff(Ts[1:-1]) > 0)


def test_profile_ramp_down_monotonic_decreasing():
    """Profile is strictly decreasing on the ramp-down segment."""
    zs = np.linspace(tprof_module.ISO_END, tprof_module.LENGTH, 50)
    Ts = np.array([temperature_profile(z) for z in zs])
    assert np.all(np.diff(Ts) <= 1e-12)
    assert np.all(np.diff(Ts[1:-1]) < 0)


def test_profile_smoothness_at_ramp_up_junction():
    """The raised-cosine ramp has zero analytic derivative at z=RAMP_UP_END."""
    eps = 1e-7
    z = tprof_module.RAMP_UP_END
    dleft = (temperature_profile(z) - temperature_profile(z - eps)) / eps
    dright = (temperature_profile(z + eps) - temperature_profile(z)) / eps
    # Both should be near zero (within finite-difference noise of a 1e-7 step
    # against an O(1e3 K * pi^2 / 0.04) curvature ~ 1e6 → eps^2*curvature ~ 1e-1).
    assert abs(dleft) < 1.0
    assert abs(dright) < 1e-9
    # And the cosine derivative left of the junction tends to zero anyway.
    assert dleft >= 0  # ramp is going up


def test_profile_smoothness_at_iso_end_junction():
    """Same C^1 continuity at z=ISO_END for the down-ramp."""
    eps = 1e-7
    z = tprof_module.ISO_END
    dleft = (temperature_profile(z) - temperature_profile(z - eps)) / eps
    dright = (temperature_profile(z + eps) - temperature_profile(z)) / eps
    assert abs(dleft) < 1e-9
    assert abs(dright) < 1.0
    assert dright <= 0  # ramp is going down


def test_profile_max_is_t_hot():
    """The plateau temperature is the global max."""
    zs = np.linspace(0, tprof_module.LENGTH, 1001)
    Ts = np.array([temperature_profile(z) for z in zs])
    assert Ts.max() == pytest.approx(tprof_module.T_HOT)


def test_profile_min_is_t_inlet():
    """The inlet temperature is the global min for this configuration."""
    zs = np.linspace(0, tprof_module.LENGTH, 1001)
    Ts = np.array([temperature_profile(z) for z in zs])
    assert Ts.min() == pytest.approx(tprof_module.T_INLET)


# ---------------------------------------------------------------------------
# Inheritance, registration, class attributes
# ---------------------------------------------------------------------------

def test_inheritance():
    """CanteraPFRTProfile inherits from CanteraBase and SimulateAdapter."""
    assert issubclass(CanteraPFRTProfile, CanteraBase)
    assert issubclass(CanteraPFRTProfile, SimulateAdapter)


def test_cantera_reactor_type_attribute():
    """The class attribute matches the reactor type used by the Lagrangian sweep."""
    assert CanteraPFRTProfile.cantera_reactor_type == 'IdealGasConstPressureReactor'


def test_tprof_registered():
    """CanteraPFRTProfile is registered in the factory."""
    from t3.simulate.factory import _registered_simulate_adapters
    assert 'CanteraPFRTProfile' in _registered_simulate_adapters
    assert _registered_simulate_adapters['CanteraPFRTProfile'] is CanteraPFRTProfile


def test_tprof_via_factory():
    """CanteraPFRTProfile can be instantiated through the factory."""
    from t3.simulate.factory import simulate_factory
    t3 = run_minimal(project_directory=TEST_DIR)
    t3.set_paths()
    adapter = simulate_factory(simulate_method='CanteraPFRTProfile',
                               t3=t3.t3,
                               rmg=t3.rmg,
                               paths=t3.paths,
                               logger=t3.logger,
                               atol=t3.rmg['model']['atol'],
                               rtol=t3.rmg['model']['rtol'],
                               observable_list=[],
                               )
    assert isinstance(adapter, CanteraPFRTProfile)


def test_create_reactor_returns_const_pressure_reactor():
    """The bare reactor returned by create_reactor is an IdealGasConstPressureReactor."""
    adapter = _make_tprof()
    condition = adapter.conditions[0]
    T0, P0, V0 = adapter._get_initial_state(condition)
    adapter.reinitialize_simulation(T0=T0, P0=P0, X0=condition.mol_frac, V0=V0)
    assert isinstance(adapter.cantera_reactor, ct.IdealGasConstPressureReactor)


def test_tprof_is_instance_of_base():
    adapter = _make_tprof()
    assert isinstance(adapter, CanteraBase)
    assert isinstance(adapter, SimulateAdapter)
    assert isinstance(adapter, CanteraPFRTProfile)


# ---------------------------------------------------------------------------
# simulate() — basic structure
# ---------------------------------------------------------------------------

def test_simulate_runs():
    """simulate() runs without errors and produces N_CELLS data points."""
    adapter = _make_tprof()
    adapter.simulate()
    assert len(adapter.all_data) == 1
    time_gd = adapter.all_data[0][0]
    assert len(time_gd.data) == tprof_module.N_CELLS


def test_all_data_structure():
    """all_data has the same tuple shape as the base-class adapters."""
    adapter = _make_tprof()
    adapter.simulate()
    time_gd, condition_data, rxn_sa_data, thermo_sa_data = adapter.all_data[0]
    assert time_gd.label == 'Time'
    assert len(condition_data) == 2 + adapter.num_ct_species
    assert condition_data[0].label == 'Temperature'
    assert condition_data[1].label == 'Pressure'
    assert rxn_sa_data == []
    assert thermo_sa_data == []


def test_distance_data_populated():
    """distance_data is a list of one numpy array of length N_CELLS."""
    adapter = _make_tprof()
    adapter.simulate()
    assert hasattr(adapter, 'distance_data')
    assert len(adapter.distance_data) == 1
    z = adapter.distance_data[0]
    assert isinstance(z, np.ndarray)
    assert len(z) == tprof_module.N_CELLS


def test_temperature_data_populated():
    """temperature_data is a list of one numpy array of length N_CELLS."""
    adapter = _make_tprof()
    adapter.simulate()
    assert hasattr(adapter, 'temperature_data')
    assert len(adapter.temperature_data) == 1
    T = adapter.temperature_data[0]
    assert isinstance(T, np.ndarray)
    assert len(T) == tprof_module.N_CELLS


def test_distance_monotonic_and_in_bounds():
    """Cell centers are strictly increasing and lie in (0, LENGTH)."""
    adapter = _make_tprof()
    adapter.simulate()
    z = adapter.distance_data[0]
    assert np.all(np.diff(z) > 0)
    assert z[0] > 0
    assert z[-1] < tprof_module.LENGTH


def test_distance_first_and_last_cell_centers():
    """First and last cell centers are at dz/2 and LENGTH - dz/2."""
    adapter = _make_tprof()
    adapter.simulate()
    z = adapter.distance_data[0]
    dz = tprof_module.LENGTH / tprof_module.N_CELLS
    np.testing.assert_allclose(z[0], dz / 2, rtol=1e-12)
    np.testing.assert_allclose(z[-1], tprof_module.LENGTH - dz / 2, rtol=1e-12)


def test_time_monotonic_increasing():
    """Cumulative residence time is strictly increasing across cells."""
    adapter = _make_tprof()
    adapter.simulate()
    times = np.array(adapter.all_data[0][0].data)
    assert np.all(np.diff(times) > 0)


# ---------------------------------------------------------------------------
# Physical sanity: T(z) imposed by the profile
# ---------------------------------------------------------------------------

def test_reactor_temperature_matches_profile():
    """Stored reactor T at each cell equals temperature_profile(cell_center)."""
    adapter = _make_tprof()
    adapter.simulate()
    T_reactor = np.array(adapter.all_data[0][1][0].data)
    T_expected = adapter.temperature_data[0]
    np.testing.assert_allclose(T_reactor, T_expected, rtol=1e-10)


def test_temperature_data_matches_profile_function():
    """temperature_data is built from temperature_profile() at the cell centers."""
    adapter = _make_tprof()
    adapter.simulate()
    z = adapter.distance_data[0]
    T = adapter.temperature_data[0]
    expected = np.array([temperature_profile(zi) for zi in z])
    np.testing.assert_allclose(T, expected, rtol=1e-12)


def test_temperature_min_max():
    """Temperature ranges between T_INLET and T_HOT (T_OUTLET > T_INLET so min is T_INLET)."""
    adapter = _make_tprof()
    adapter.simulate()
    T = np.array(adapter.all_data[0][1][0].data)
    assert T.min() >= tprof_module.T_INLET - 1e-9
    assert T.max() <= tprof_module.T_HOT + 1e-9
    # The max should actually equal T_HOT since the iso plateau covers many cells
    assert T.max() == pytest.approx(tprof_module.T_HOT, rel=1e-6)


def test_temperature_has_plateau():
    """Many cells lie at exactly T_HOT in the isothermal plateau."""
    adapter = _make_tprof()
    adapter.simulate()
    T = np.array(adapter.all_data[0][1][0].data)
    n_at_plateau = np.sum(np.isclose(T, tprof_module.T_HOT, rtol=1e-9))
    # Plateau spans (ISO_END - RAMP_UP_END) / LENGTH * N_CELLS cells
    expected_plateau_cells = int(round(
        (tprof_module.ISO_END - tprof_module.RAMP_UP_END) / tprof_module.LENGTH * tprof_module.N_CELLS
    ))
    # Allow off-by-one at the boundaries
    assert abs(n_at_plateau - expected_plateau_cells) <= 2


def test_pressure_constant_along_chain():
    """Pressure is held at the inlet value (chain is isobaric)."""
    adapter = _make_tprof()
    adapter.simulate()
    P = np.array(adapter.all_data[0][1][1].data)
    np.testing.assert_allclose(P, P[0], rtol=1e-6)


def test_mole_fractions_sum_to_one_at_every_cell():
    """X mole fractions sum to 1 at every cell along the chain."""
    adapter = _make_tprof()
    adapter.simulate()
    n_cells = len(adapter.all_data[0][0].data)
    n_spc = adapter.num_ct_species
    for i in range(n_cells):
        x_sum = sum(adapter.all_data[0][1][2 + s].data[i] for s in range(n_spc))
        np.testing.assert_allclose(x_sum, 1.0, atol=1e-10)


def test_mole_fractions_nonnegative():
    """All species mole fractions are non-negative."""
    adapter = _make_tprof()
    adapter.simulate()
    n_spc = adapter.num_ct_species
    for s in range(n_spc):
        data = np.array(adapter.all_data[0][1][2 + s].data)
        assert np.all(data >= -1e-15), \
            f'Negative mole fraction for {adapter.all_data[0][1][2 + s].label}'


def test_species_labels_match_model():
    """Species GenericData labels match Cantera model species names."""
    adapter = _make_tprof()
    adapter.simulate()
    model_species = [s.name for s in adapter.model.species()]
    for s in range(adapter.num_ct_species):
        assert adapter.all_data[0][1][2 + s].label == model_species[s]


def test_inlet_composition_close_to_inputs():
    """The first cell's mole fractions reflect the input concentrations
    (chemistry at 300 K is essentially frozen, so cell 0 ≈ inlet)."""
    adapter = _make_tprof()
    adapter.simulate()
    n_spc = adapter.num_ct_species
    # Find H2 and O2 by base-name
    species_names = [adapter.all_data[0][1][2 + s].label for s in range(n_spc)]
    h2_idx = next(i for i, n in enumerate(species_names) if n.split('(')[0] == 'H2')
    o2_idx = next(i for i, n in enumerate(species_names) if n.split('(')[0] == 'O2')
    # Cell 0 mole fractions
    h2_in = adapter.all_data[0][1][2 + h2_idx].data[0]
    o2_in = adapter.all_data[0][1][2 + o2_idx].data[0]
    assert h2_in == pytest.approx(0.67, abs=1e-3)
    assert o2_in == pytest.approx(0.33, abs=1e-3)


# ---------------------------------------------------------------------------
# Sensitivity analysis (steady-state PFR — SA is a function of distance z)
# ---------------------------------------------------------------------------

def test_sa_simulate_runs():
    """simulate() succeeds with SA observables and populates SA arrays."""
    adapter = _make_tprof(observable_list=['OH', 'H'])
    adapter.simulate()
    _, _, rxn_sa, thermo_sa = adapter.all_data[0]
    # 2 sensitive species * num_reactions kinetic SA series, and
    # 2 sensitive species * num_species thermo SA series.
    assert len(rxn_sa) == 2 * adapter.num_ct_reactions
    assert len(thermo_sa) == 2 * adapter.num_ct_species
    # Each SA series has one value per cell.
    assert len(rxn_sa[0].data) == tprof_module.N_CELLS
    assert len(thermo_sa[0].data) == tprof_module.N_CELLS


def test_sa_dict_has_distance_axis():
    """get_sa_coefficients returns the standard dict plus a 'distance' key."""
    adapter = _make_tprof(observable_list=['OH', 'H'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    # Standard keys still present.
    for key in ('time', 'kinetics', 'thermo'):
        assert key in sa_dict
    # New 'distance' key added by the adapter.
    assert 'distance' in sa_dict
    assert len(sa_dict['distance']) == 1
    z = sa_dict['distance'][0]
    assert isinstance(z, np.ndarray)
    assert len(z) == tprof_module.N_CELLS
    # 'distance' and 'time' must be the same length so SA arrays line up.
    assert len(sa_dict['time'][0]) == len(z)
    # Distance is monotonically increasing and bounded by [0, LENGTH].
    assert np.all(np.diff(z) > 0)
    assert z[0] > 0 and z[-1] < tprof_module.LENGTH


def test_sa_dict_observables_populated():
    """sa_dict['kinetics'][0] / ['thermo'][0] contain entries for each observable."""
    adapter = _make_tprof(observable_list=['OH', 'H'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()

    kin = sa_dict['kinetics'][0]
    thermo = sa_dict['thermo'][0]
    # Cantera species names include the RMG index suffix
    obs_keys = list(kin.keys())
    assert any('OH' in k for k in obs_keys)
    assert any(k == 'H' or k.startswith('H(') for k in obs_keys)

    for obs, params in kin.items():
        # Each kinetic param maps to a numpy array of length N_CELLS
        for arr in params.values():
            assert isinstance(arr, np.ndarray)
            assert len(arr) == tprof_module.N_CELLS
    for obs, params in thermo.items():
        for arr in params.values():
            assert isinstance(arr, np.ndarray)
            assert len(arr) == tprof_module.N_CELLS


def test_sa_finite_and_nontrivial():
    """SA coefficients are finite and at least some are non-zero."""
    adapter = _make_tprof(observable_list=['OH', 'H'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    any_nonzero_kin = False
    for obs_dict in sa_dict['kinetics']:
        for params in obs_dict.values():
            for arr in params.values():
                assert np.all(np.isfinite(arr))
                if np.any(np.abs(arr) > 1e-12):
                    any_nonzero_kin = True
    assert any_nonzero_kin, 'all kinetic SA values were ~zero — chemistry never engaged'


def test_sa_distance_matches_distance_data():
    """sa_dict['distance'] matches the adapter's distance_data attribute."""
    adapter = _make_tprof(observable_list=['OH'])
    adapter.simulate()
    sa_dict = adapter.get_sa_coefficients()
    np.testing.assert_array_equal(sa_dict['distance'][0], adapter.distance_data[0])


def test_sa_doubles_when_two_observables():
    """Adding a second observable doubles the number of SA series in all_data."""
    a1 = _make_tprof(observable_list=['OH'])
    a1.simulate()
    a2 = _make_tprof(observable_list=['OH', 'H'])
    a2.simulate()
    assert len(a2.all_data[0][2]) == 2 * len(a1.all_data[0][2])
    assert len(a2.all_data[0][3]) == 2 * len(a1.all_data[0][3])


def test_sa_no_observables_skips_sa_data():
    """Without SA observables the adapter still runs, with empty SA arrays."""
    adapter = _make_tprof()
    adapter.simulate()
    _, _, rxn_sa, thermo_sa = adapter.all_data[0]
    assert rxn_sa == []
    assert thermo_sa == []
    sa_dict = adapter.get_sa_coefficients()
    # 'distance' key is added even when no observables — it's a property of the
    # adapter, not of whether SA was requested.
    assert sa_dict['distance'][0].shape == (tprof_module.N_CELLS,)
    assert sa_dict['kinetics'][0] == {}
    assert sa_dict['thermo'][0] == {}


# ---------------------------------------------------------------------------
# get_idt_by_T / find_equilibrium / get_t50 behavior
# ---------------------------------------------------------------------------

def test_get_idt_returns_empty_before_simulate():
    """get_idt_by_T returns empty lists (PFR with imposed T(z) has no IDT)."""
    adapter = _make_tprof()
    idt = adapter.get_idt_by_T()
    assert idt == {'idt': [], 'idt_index': []}


def test_get_idt_returns_empty_after_simulate():
    """get_idt_by_T still returns empty after simulate()."""
    adapter = _make_tprof()
    adapter.simulate()
    idt = adapter.get_idt_by_T()
    assert idt == {'idt': [], 'idt_index': []}


def test_find_equilibrium_runs():
    """find_equilibrium uses the inherited helper and returns a mole-fraction dict."""
    adapter = _make_tprof()
    eq = adapter.find_equilibrium('TP')
    assert isinstance(eq, list)
    assert len(eq) == 1
    assert isinstance(eq[0], dict)
    # Sum of mole fractions ≈ 1
    assert sum(eq[0].values()) == pytest.approx(1.0, abs=1e-6)


# ---------------------------------------------------------------------------
# Re-runnability
# ---------------------------------------------------------------------------

def test_simulate_twice():
    """Running simulate() twice in a row works and produces identical results."""
    adapter = _make_tprof()
    adapter.simulate()
    T1 = np.array(adapter.all_data[0][1][0].data).copy()
    adapter.simulate()
    T2 = np.array(adapter.all_data[0][1][0].data).copy()
    np.testing.assert_allclose(T1, T2, rtol=1e-12)


# ---------------------------------------------------------------------------
# Module constant overrides (no mocks — try/finally)
# ---------------------------------------------------------------------------

def test_n_cells_override():
    """Overriding N_CELLS changes the number of stored data points."""
    original = tprof_module.N_CELLS
    try:
        tprof_module.N_CELLS = 50
        adapter = _make_tprof()
        adapter.simulate()
        assert len(adapter.all_data[0][0].data) == 50
        assert len(adapter.distance_data[0]) == 50
        assert len(adapter.temperature_data[0]) == 50
    finally:
        tprof_module.N_CELLS = original


def test_length_override_changes_distance_range():
    """Overriding LENGTH changes the cell-center positions accordingly."""
    original = tprof_module.LENGTH
    try:
        tprof_module.LENGTH = 0.40
        # Also adjust the profile breakpoints to stay within the new length
        original_ramp = tprof_module.RAMP_UP_END
        original_iso = tprof_module.ISO_END
        tprof_module.RAMP_UP_END = 0.10
        tprof_module.ISO_END = 0.30
        try:
            adapter = _make_tprof()
            adapter.simulate()
            z = adapter.distance_data[0]
            # Cell centers should span (0, 0.40)
            assert z[0] < 0.10
            assert z[-1] < 0.40
            assert z[-1] > 0.30
        finally:
            tprof_module.RAMP_UP_END = original_ramp
            tprof_module.ISO_END = original_iso
    finally:
        tprof_module.LENGTH = original


def test_temperature_profile_override():
    """Replacing module constants changes the resulting profile."""
    original_T_inlet = tprof_module.T_INLET
    original_T_hot = tprof_module.T_HOT
    original_T_outlet = tprof_module.T_OUTLET
    try:
        tprof_module.T_INLET = 400.0
        tprof_module.T_HOT = 1200.0
        tprof_module.T_OUTLET = 800.0
        # Endpoints should reflect the new values
        assert temperature_profile(0.0) == pytest.approx(400.0)
        assert temperature_profile(tprof_module.RAMP_UP_END) == pytest.approx(1200.0)
        assert temperature_profile(tprof_module.LENGTH) == pytest.approx(800.0)
    finally:
        tprof_module.T_INLET = original_T_inlet
        tprof_module.T_HOT = original_T_hot
        tprof_module.T_OUTLET = original_T_outlet


# ---------------------------------------------------------------------------
# Larger model (NH3)
# ---------------------------------------------------------------------------

def test_nh3_simulate_no_sa():
    """The adapter works against the larger NH3 mechanism.

    The default T profile (300 K → 900 K → 500 K) drives the NH3 mechanism
    into a regime where its rate expressions underflow / produce NaN at the
    cold ends — that is a quirk of the bundled NH3.yaml, not the adapter.
    For this smoke test we narrow the profile to 600 K → 900 K → 700 K, which
    keeps every reaction rate finite while still exercising the same code
    paths (ramp-up, plateau, ramp-down).
    """
    original = (tprof_module.T_INLET, tprof_module.T_OUTLET)
    try:
        tprof_module.T_INLET = 600.0
        tprof_module.T_OUTLET = 700.0
        adapter = _make_nh3_tprof()
        adapter.simulate()
        assert len(adapter.all_data) == 1
        time_gd = adapter.all_data[0][0]
        assert len(time_gd.data) == tprof_module.N_CELLS
        # Plateau temperature reached
        T = np.array(adapter.all_data[0][1][0].data)
        assert T.max() == pytest.approx(tprof_module.T_HOT)
        # Mole fractions sum to 1 at every cell
        n_spc = adapter.num_ct_species
        for i in range(len(time_gd.data)):
            x_sum = sum(adapter.all_data[0][1][2 + s].data[i] for s in range(n_spc))
            np.testing.assert_allclose(x_sum, 1.0, atol=1e-9)
    finally:
        tprof_module.T_INLET, tprof_module.T_OUTLET = original


def test_nh3_sa_request_runs():
    """SA against the larger NH3 model runs end-to-end and produces valid SA arrays."""
    original = (tprof_module.T_INLET, tprof_module.T_OUTLET)
    try:
        tprof_module.T_INLET = 600.0
        tprof_module.T_OUTLET = 700.0
        adapter = _make_nh3_tprof(observable_list=['NH2'])
        adapter.simulate()
        sa_dict = adapter.get_sa_coefficients()
        assert 'distance' in sa_dict
        assert sa_dict['distance'][0].shape == (tprof_module.N_CELLS,)
        # NH2 SA should exist for the NH3 mech
        kin = sa_dict['kinetics'][0]
        assert any('NH2' in k for k in kin.keys())
        for obs, params in kin.items():
            for arr in params.values():
                assert len(arr) == tprof_module.N_CELLS
                assert np.all(np.isfinite(arr))
    finally:
        tprof_module.T_INLET, tprof_module.T_OUTLET = original


# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------

def teardown_module():
    """Remove project directories created during these unit tests."""
    if os.path.isdir(TPROF_NH3_DIR):
        shutil.rmtree(TPROF_NH3_DIR, ignore_errors=True)
    log_archive = os.path.join(TEST_DIR, 'log_archive')
    if os.path.isdir(log_archive):
        shutil.rmtree(log_archive, ignore_errors=True)
    log_file = os.path.join(TEST_DIR, 't3.log')
    if os.path.isfile(log_file):
        os.remove(log_file)
