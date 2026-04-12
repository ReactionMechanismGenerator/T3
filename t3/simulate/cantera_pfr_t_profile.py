"""
Cantera Simulator Adapter module
================================

Plug flow reactor (PFR) with a hardcoded **axial temperature profile**.

This adapter is the *prescribed-temperature* (Scenario A) version of a
non-isothermal PFR -- the wall temperature ``T(z)`` is known a priori
(e.g., a measured profile from a flow-tube experiment) and the gas
adopts that temperature at every axial position.  The chemistry then
evolves under that imposed schedule.

Method
------

A single Lagrangian particle is advanced through the reactor:

* The gas is wrapped in a ``ct.IdealGasConstPressureReactor`` with
  ``energy='off'`` -- with the energy equation disabled, the temperature
  is whatever we set externally and the integrator only advances the
  composition.
* The reactor length is discretized into ``N_CELLS`` equal segments.
  At each segment center ``z_n = (n + 1/2) * dz`` we look up
  ``T_n = temperature_profile(z_n)``, override ``self.model.TP``, call
  ``reactor.syncState()``, and advance the integrator by
  ``dt_n = dz / u_n``, where the local velocity ``u_n`` is recomputed at
  the new density: ``u_n = mdot / (rho_n * AREA)``.
* The mass flow rate is constant along the chain
  ``mdot = rho_inlet * u_inlet * AREA`` with
  ``u_inlet = LENGTH / total_residence_time`` taken from the reactor
  entry's ``termination_time``.

Because the gas density rises and falls with ``T``, the velocity (and
therefore ``dt`` per segment) varies dynamically along the tube -- this
is exactly the recommendation in
https://cantera.org/examples/python/reactors/pfr.py.html and the
write-up in Cantera/cantera#523 for non-isothermal flow-tube modelling.

Defaults (the configuration this file ships with):

    Total length     : 80 cm
    z in [0,    0.20]: smooth (raised-cosine) ramp 300 K -> 900 K
    z in [0.20, 0.60]: isothermal at 900 K
    z in [0.60, 0.80]: smooth (raised-cosine) ramp 900 K -> 500 K

The raised-cosine ramp ``(1 - cos(pi t)) / 2`` has zero derivative at both
ends, so the profile is C^1 continuous at the junctions between the ramp
and isothermal segments -- no kink in ``T(z)`` for the integrator to
stumble over.

Sensitivity analysis
--------------------

SA **is** supported.  Because the reactor is steady-state in the lab frame,
the natural independent variable for the SA coefficients is the axial
position ``z`` (not the residence time ``t`` of the Lagrangian particle).
Each cell sets up the gas at its imposed ``T_n`` and the reactor network
advances by ``dt_n``; Cantera's per-step sensitivity matrix is preserved
across the ``gas.TP`` override (verified empirically) so the sensitivity
ODEs continue evolving consistently from cell to cell.

After ``simulate()`` runs, ``get_sa_coefficients()`` returns the usual
``{'time', 'kinetics', 'thermo'}`` dictionary (so all of T3's existing
SA-consuming code works unchanged), with one extra key:

* ``sa_dict['distance']`` -- a list of one numpy array per condition giving
  the segment-center axial position [m] for every SA sample.  ``time`` and
  ``distance`` have the same length.  Code that wants to plot SA against
  ``z`` should index by ``sa_dict['distance']``; the ``time`` axis is the
  cumulative residence time of the Lagrangian particle and is provided for
  interface compatibility with the other adapters.

Customising the profile
-----------------------

Advanced users have two equally good options:

1.  Tune the module-level constants below (geometry and segment
    temperatures) and/or replace the body of :func:`temperature_profile`
    with their own ``T(z)``.  This file ships with a smooth raised-cosine
    profile, but a piecewise-linear ramp, a measured CSV interpolated via
    ``np.interp``, or any other function works the same way as long as
    it returns a temperature [K] for an axial position [m].
2.  Copy this file as the starting point for a brand-new adapter (and
    rename the registered ``CanteraPFRTProfile`` accordingly).
"""

import cantera as ct
import numpy as np

from t3.simulate.cantera_base import CanteraBase
from t3.simulate.factory import register_simulate_adapter
from t3.utils.rmg_shim import GenericData


# ---------------------------------------------------------------------------
# Module-level constants -- users modify these for their application
# ---------------------------------------------------------------------------
LENGTH = 0.80          # m, total reactor length (= 80 cm)
AREA = 1e-4            # m^2, cross-sectional area
N_CELLS = 200          # number of CSTR cells along the chain

# Temperature-profile geometry (must satisfy 0 < RAMP_UP_END < ISO_END < LENGTH)
RAMP_UP_END = 0.20     # m, end of the ramp-up region
ISO_END = 0.60         # m, end of the isothermal region
T_INLET = 300.0        # K, inlet temperature
T_HOT = 900.0          # K, isothermal plateau temperature
T_OUTLET = 500.0       # K, outlet temperature

if not (0 < RAMP_UP_END < ISO_END < LENGTH):
    raise ValueError(
        f'Temperature profile constants must satisfy 0 < RAMP_UP_END ({RAMP_UP_END}) < ISO_END ({ISO_END}) < LENGTH ({LENGTH})'
    )


def temperature_profile(z):
    """
    Return the wall/gas temperature [K] at axial position ``z`` [m].

    The default profile uses raised-cosine ramps (C^1 continuous at the
    junctions) for the up-ramp and down-ramp segments::

         T(z) = T_INLET + (T_HOT  - T_INLET) * 0.5 * (1 - cos(pi * z / RAMP_UP_END))
                                                     for 0   <= z < RAMP_UP_END
         T(z) = T_HOT                                for RAMP_UP_END <= z < ISO_END
         T(z) = T_HOT + (T_OUTLET - T_HOT) * 0.5
                * (1 - cos(pi * (z - ISO_END) / (LENGTH - ISO_END)))
                                                     for ISO_END <= z <= LENGTH

    Args:
        z (float): Axial position [m].  Clipped to [0, LENGTH].

    Returns:
        float: Temperature [K] at ``z``.
    """
    z = float(np.clip(z, 0.0, LENGTH))
    if z < RAMP_UP_END:
        t = z / RAMP_UP_END
        return T_INLET + (T_HOT - T_INLET) * 0.5 * (1 - np.cos(np.pi * t))
    if z < ISO_END:
        return T_HOT
    t = (z - ISO_END) / (LENGTH - ISO_END)
    return T_HOT + (T_OUTLET - T_HOT) * 0.5 * (1 - np.cos(np.pi * t))


class CanteraPFRTProfile(CanteraBase):
    """
    Plug flow reactor with an axial temperature profile (Lagrangian).

    A single ``ct.IdealGasConstPressureReactor`` with ``energy='off'`` is
    advected through ``N_CELLS`` axial segments.  At each segment we sample
    ``T = temperature_profile(z_center)``, override ``self.model.TP``, call
    ``reactor.syncState()``, and advance the network by ``dt = dz / u``.
    The local velocity ``u = mdot / (rho * AREA)`` is recomputed at every
    segment because the density of an isobaric ideal gas changes with T.

    See the module docstring for the modelling rationale and references.

    Notes
    -----
    *   The energy equation is **off** -- the temperature is imposed
        externally by ``temperature_profile`` rather than evolved by the
        integrator.  This is the same idealization as a measured-T(z) flow
        reactor where the wall temperature is the controlling variable and
        the gas equilibrates to the wall thermally.
    *   Pressure is held constant at the inlet value (``rmg.reactors[0]['P']``).
        Density evolves with T as the parcel marches through the reactor.
    *   ``self.distance_data`` (one numpy array per condition) holds the
        segment-center positions [m].
    *   ``self.temperature_data`` (one numpy array per condition) holds the
        imposed segment temperatures [K] -- handy for plotting alongside
        the species profiles.
    *   ``self.all_data`` matches the shape produced by the other Cantera
        adapters: a list of
        ``(time_gd, [T_gd, P_gd, *species_gd], rxn_sa, thermo_sa)`` tuples.
    *   Sensitivity analysis **is** supported.  The reactor is steady-state
        in the lab frame, so the natural independent variable for the SA
        coefficients is the axial position ``z``.
        :meth:`get_sa_coefficients` adds a ``'distance'`` key to the usual
        ``sa_dict`` (alongside ``'time'``) carrying the cell-center
        positions [m].

    Args:
        t3 (dict): The T3.t3 attribute.
        rmg (dict): The T3.rmg attribute.
        paths (dict): The T3.paths attribute.
        logger (Logger): Instance of T3's Logger class.
        atol (float): Absolute tolerance for the time integrator.
        rtol (float): Relative tolerance for the time integrator.
        observable_list (Optional[list]): Species used for SA. Not supported.
        sa_atol (float): Absolute tolerance for SA. Not supported.
        sa_rtol (float): Relative tolerance for SA. Not supported.
        global_observables (Optional[List[str]]): Global observables. Not supported.
    """

    cantera_reactor_type = 'IdealGasConstPressureReactor'

    def create_reactor(self):
        """
        Bare reactor used by inherited helpers (e.g. :meth:`find_equilibrium`).

        The actual Lagrangian sweep is built inside :meth:`simulate`.  This
        method just returns a const-pressure reactor with the energy
        equation disabled, matching the type used during the sweep.
        """
        return ct.IdealGasConstPressureReactor(self.model, energy='off')

    def get_idt_by_T(self):
        """A non-isothermal PFR with imposed wall T(z) has no chemistry-driven IDT."""
        return {'idt': list(), 'idt_index': list()}

    def simulate(self):
        """
        Sweep the Lagrangian particle through the reactor.

        For each axial segment ``n in [0, N_CELLS)``:

        1. Compute the segment center ``z_n = (n + 1/2) * dz`` and the
           imposed segment temperature ``T_n = temperature_profile(z_n)``.
        2. Override ``self.model.TP = (T_n, P_in)`` (the composition is
           inherited from the previous segment) and call
           ``reactor.syncState()`` so the integrator picks up the new state.
        3. Recompute the local velocity ``u_n = mdot / (rho_n * AREA)``
           at the new density and step forward in time by
           ``dt_n = dz / u_n``.
        4. If SA observables are configured, read the per-cell sensitivity
           matrix from Cantera (it is preserved across the ``gas.TP``
           override) and convert it from the Cantera mass-fraction basis to
           the mole-fraction basis used by the rest of T3 -- exactly the
           same correction the base class applies in :meth:`CanteraBase.simulate`.

        After the sweep, ``self.all_data`` follows the same shape as the
        other Cantera adapters and ``self.distance_data`` /
        ``self.temperature_data`` carry the matching axial axis.
        """
        sa_enabled = bool(self.sensitive_species)
        self.logger.info(
            f'Running a PFR-with-T-profile Lagrangian-particle simulation '
            f'{"(with SA) " if sa_enabled else ""}'
            f'using {self.__class__.__name__} (LENGTH={LENGTH} m, N_CELLS={N_CELLS})...')

        species_names_list = [species.name for species in self.model.species()]
        self.all_data = list()
        self.distance_data = list()
        self.temperature_data = list()

        dz = LENGTH / N_CELLS
        # Segment centers (where T is sampled)
        cell_centers = np.array([(n + 0.5) * dz for n in range(N_CELLS)])
        cell_temperatures = np.array([temperature_profile(z) for z in cell_centers])

        # Index mask for the non-inert species (used by the SA correction).
        non_inert_mask = np.ones(self.num_ct_species, dtype=bool)
        non_inert_mask[self.inert_index_list] = False

        for condition in self.conditions:
            T0, P0, V0 = self._get_initial_state(condition)
            if P0 is None:
                raise ValueError(
                    'CanteraPFRTProfile requires a pressure-based initial state '
                    '(T, P, X).  V0 was specified instead.')
            P_in = P0

            # Inlet density and constant mass flow rate.  The inlet velocity
            # is set so that an isothermal-at-T0 trip down the tube would take
            # exactly the configured termination_time -- this fixes mdot.
            self.model.TPX = T0, P_in, condition.mol_frac
            total_residence_time = condition.reaction_time.value_si
            u_inlet = LENGTH / total_residence_time
            mass_flow_rate = self.model.density * u_inlet * AREA

            # Build the Lagrangian particle reactor (a constant-pressure
            # ideal-gas batch with energy='off') seeded at the first segment's
            # imposed temperature.
            self.model.TPX = cell_temperatures[0], P_in, condition.mol_frac
            reactor = ct.IdealGasConstPressureReactor(self.model, energy='off')
            sim = ct.ReactorNet([reactor])
            sim.atol = self.atol
            sim.rtol = self.rtol

            if sa_enabled:
                # Register every reaction and every species enthalpy for SA
                # before the first integration step (matches the base class).
                for i in range(self.num_ct_reactions):
                    reactor.add_sensitivity_reaction(i)
                for i in range(self.num_ct_species):
                    reactor.add_sensitivity_species_enthalpy(i)
                sim.atol_sensitivity = self.sa_atol
                sim.rtol_sensitivity = self.sa_rtol

            # Keep base-class attributes in sync (used by inherited helpers).
            self.cantera_reactor = reactor
            self.cantera_simulation = sim

            times = list()
            temperatures = list()
            pressures = list()
            species_data = list()
            kinetic_sensitivity_data = list()
            thermo_sensitivity_data = list()

            for n in range(N_CELLS):
                T_cell = cell_temperatures[n]

                # Override the temperature at the start of this segment.  The
                # composition is whatever the previous step left in the
                # reactor (or the inlet for n == 0).  syncState pushes the
                # new T into the reactor without touching the SA matrix.
                self.model.TP = T_cell, P_in
                reactor.syncState()

                # Local velocity at the new (T, rho) and the corresponding
                # time step to traverse one segment.
                rho = self.model.density
                if rho * AREA <= 0:
                    raise ValueError(f'Invalid density ({rho}) or area ({AREA}) in PFR T-profile cell {n}')
                u = mass_flow_rate / (rho * AREA)
                if u <= 0:
                    raise ValueError(f'Invalid velocity ({u}) in PFR T-profile cell {n}')
                dt = dz / u

                sim.advance(sim.time + dt)

                times.append(sim.time)
                temperatures.append(reactor.T)
                pressures.append(reactor.thermo.P)
                species_data.append(reactor.thermo[species_names_list].X)

                if sa_enabled:
                    # Cantera returns mass-based sensitivities; convert to
                    # mole-fraction basis using the same identity the base
                    # class uses:
                    #   d ln x_k = d ln w_k - sum_i (d ln w_i) * x_i
                    # over non-inert i.
                    mass_frac_sa = sim.sensitivities()
                    # IdealGasConstPressureReactor: rows are
                    # [mass, enthalpy/temperature, then species].
                    mass_frac_sa = mass_frac_sa[2:, :]
                    x = np.array(species_data[-1])
                    mass_frac_sa *= x[:, np.newaxis]

                    kin_correction = mass_frac_sa[non_inert_mask, :self.num_ct_reactions].sum(axis=0)
                    thermo_correction = mass_frac_sa[non_inert_mask, self.num_ct_reactions:].sum(axis=0)

                    kin_sa = np.zeros(len(self.sensitive_species) * self.num_ct_reactions)
                    thermo_sa = np.zeros(len(self.sensitive_species) * self.num_ct_species)

                    for index, species in enumerate(self.sensitive_species):
                        raw_kin = np.array([sim.sensitivity(species, j)
                                            for j in range(self.num_ct_reactions)])
                        start = self.num_ct_reactions * index
                        kin_sa[start:start + self.num_ct_reactions] = raw_kin - kin_correction

                        raw_thermo = np.array([
                            sim.sensitivity(species, j + self.num_ct_reactions)
                            for j in range(self.num_ct_species)])
                        start = self.num_ct_species * index
                        thermo_sa[start:start + self.num_ct_species] = raw_thermo - thermo_correction

                    kinetic_sensitivity_data.append(kin_sa)
                    thermo_sensitivity_data.append(thermo_sa)

            species_data = np.array(species_data)
            kinetic_sensitivity_data = np.array(kinetic_sensitivity_data)
            thermo_sensitivity_data = np.array(thermo_sensitivity_data)

            time_gd = GenericData(label='Time', data=times, units='s')
            temp_gd = GenericData(label='Temperature', data=temperatures, units='K')
            pres_gd = GenericData(label='Pressure', data=pressures, units='Pa')
            condition_data = [temp_gd, pres_gd]

            for index, species in enumerate(self.model.species()):
                condition_data.append(GenericData(label=species.name,
                                                  species=species,
                                                  data=species_data[:, index],
                                                  index=index))

            reaction_sensitivity_data = []
            for index, species in enumerate(self.sensitive_species):
                for j in range(self.num_ct_reactions):
                    reaction_sensitivity_data.append(GenericData(
                        label='dln[{0}]/dln[k{1}]: {2}'.format(species, j + 1, self.model.reactions()[j]),
                        species=species,
                        reaction=self.model.reactions()[j],
                        data=kinetic_sensitivity_data[:, self.num_ct_reactions * index + j],
                        index=j + 1,
                    ))

            thermodynamic_sensitivity_data = []
            for index, species in enumerate(self.sensitive_species):
                for j in range(self.num_ct_species):
                    thermodynamic_sensitivity_data.append(GenericData(
                        label='dln[{0}]/dH[{1}]'.format(species, self.model.species()[j].name),
                        species=species,
                        data=thermo_sensitivity_data[:, self.num_ct_species * index + j],
                        index=j + 1,
                    ))

            self.all_data.append((time_gd, condition_data,
                                  reaction_sensitivity_data,
                                  thermodynamic_sensitivity_data))
            self.distance_data.append(cell_centers.copy())
            self.temperature_data.append(cell_temperatures.copy())

    def get_sa_coefficients(self):
        """
        Return the standard ``sa_dict`` plus an extra ``'distance'`` key.

        For this adapter the natural independent variable is the axial
        position ``z`` (the reactor is steady-state in the lab frame); we
        therefore add a ``'distance'`` key alongside the inherited
        ``'time'`` key, listing one numpy array of cell-center positions
        [m] per simulation condition.  ``time`` and ``distance`` have the
        same length, so existing T3 code that consumes ``sa_dict['time']``
        keeps working.
        """
        sa_dict = super().get_sa_coefficients()
        sa_dict['distance'] = [np.asarray(d) for d in self.distance_data]
        return sa_dict


register_simulate_adapter("CanteraPFRTProfile", CanteraPFRTProfile)
