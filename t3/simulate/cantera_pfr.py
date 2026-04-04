"""
Cantera Simulator Adapter module
Used to run mechanism analysis with Cantera as an ideal gas in an isothermal isobaric
plug flow reactor (PFR).

Two solution methods are available, controlled by the ``METHOD`` constant:

**Lagrangian** (default):
    A fluid element is tracked in time using a constant-T constant-P batch reactor.
    Distance along the reactor is reconstructed from the instantaneous velocity at
    each time step:  z += u * dt,  where u = mass_flow_rate / (A * rho).
    Sensitivity analysis is fully supported.

**Chain of reactors**:
    The PFR is discretized into ``N_CELLS`` well-mixed cells connected by
    ``MassFlowController`` / ``PressureController`` pairs.  Each cell is solved
    to steady state before marching to the next.
    If SA is requested with this method, the adapter automatically falls back to
    the Lagrangian method (results are equivalent for isothermal isobaric PFR).

Users can modify the constants below to configure the PFR for their application.
"""

import cantera as ct
import numpy as np

from t3.simulate.cantera_base import CanteraBase
from t3.simulate.factory import register_simulate_adapter
from t3.utils.rmg_shim import GenericData


# ---------------------------------------------------------------------------
# Module-level constants — users modify these for their application
# ---------------------------------------------------------------------------
METHOD = 'lagrangian'   # 'lagrangian' or 'chain'
N_CELLS = 100           # number of cells for the chain method
AREA = 1e-4             # m^2, cross-sectional area
LENGTH = 1.0            # m, total reactor length


class CanteraPFR(CanteraBase):
    """
    Simulates an ideal gas in a steady-state, isothermal, isobaric plug flow reactor.

    The independent variable is residence time (equivalent to axial distance / velocity).
    Uses ``ct.IdealGasConstPressureReactor`` with ``energy='off'`` for the Lagrangian
    method, and ``ct.IdealGasReactor`` cells with flow controllers for the chain method.

    After simulation, ``self.distance_data`` contains a list of numpy arrays (one per
    condition) giving the axial distance [m] at each data point.

    Args:
        t3 (dict): The T3.t3 attribute.
        rmg (dict): The T3.rmg attribute.
        paths (dict): The T3.paths attribute.
        logger (Logger): Instance of T3's Logger class.
        atol (float): The absolute tolerance used when integrating.
        rtol (float): The relative tolerance used when integrating.
        observable_list (Optional[list]): Species used for SA.
        sa_atol (float): The absolute tolerance for sensitivity analysis.
        sa_rtol (float): The relative tolerance for sensitivity analysis.
        global_observables (Optional[List[str]]): Global observables.
    """

    cantera_reactor_type = 'IdealGasConstPressureReactor'

    def create_reactor(self):
        """Create a constant-pressure reactor with the energy equation disabled (isothermal PFR)."""
        return ct.IdealGasConstPressureReactor(self.model, energy='off')

    def get_idt_by_T(self):
        """Isothermal PFR has no ignition."""
        return {'idt': list(), 'idt_index': list()}

    def simulate(self):
        """
        Simulate the PFR using the method specified by the ``METHOD`` module constant.

        If ``METHOD == 'chain'`` and SA is requested, the Lagrangian method is used
        instead (the governing equations are identical for isothermal isobaric PFR).
        """
        if METHOD == 'chain' and self.sensitive_species:
            self.logger.info(
                'SA is not natively supported with the chain-of-reactors method. '
                'Using the Lagrangian method for this run '
                '(results are equivalent for isothermal isobaric PFR).')
            self._simulate_lagrangian()
        elif METHOD == 'chain':
            self._simulate_chain()
        else:
            self._simulate_lagrangian()

    def _simulate_lagrangian(self):
        """
        Lagrangian particle simulation.

        Delegates to the base-class ``simulate()`` for time integration and SA,
        then post-processes to compute axial distance from the velocity field.
        """
        super().simulate()
        self._compute_distance()

    def _compute_distance(self):
        """
        Compute axial distance from time-series data stored by the base-class
        ``simulate()``.

        For an isothermal isobaric ideal gas the velocity is inversely
        proportional to mean molecular weight (T, P, and area cancel):

            u(t) = u_0 * mean_MW(0) / mean_MW(t)
            z    = cumulative integral of u * dt

        Results are stored in ``self.distance_data``.
        """
        self.distance_data = []
        molecular_weights = self.model.molecular_weights  # kg/kmol per species

        for idx, condition_data in enumerate(self.all_data):
            time_gd, data_list, _, _ = condition_data
            times = np.array(time_gd.data)

            # Build mole-fraction matrix (n_steps x n_species)
            X = np.column_stack([data_list[2 + s].data for s in range(self.num_ct_species)])

            # Mean molecular weight at each step [kg/kmol]
            mean_MW = X @ molecular_weights

            # Velocity: u = u_0 * mean_MW_0 / mean_MW  (P, T, A all cancel)
            total_residence_time = self.conditions[idx].reaction_time.value_si
            u_0 = LENGTH / total_residence_time
            u = u_0 * mean_MW[0] / mean_MW

            # Axial distance [m] via cumulative integration
            dt = np.diff(times, prepend=0.0)
            distance = np.cumsum(u * dt)

            self.distance_data.append(distance)

    def _simulate_chain(self):
        """
        Chain-of-reactors simulation.

        The PFR is divided into ``N_CELLS`` equal-length cells.  Each cell is a
        well-mixed reactor (``ct.IdealGasReactor`` with ``energy='off'``) connected
        to upstream and downstream reservoirs via a ``MassFlowController`` and a
        ``PressureController``.  Cells are solved to steady state sequentially,
        with the outlet of cell *n* feeding the inlet of cell *n+1*.

        The results are stored in ``self.all_data`` in the same format as the
        base-class ``simulate()`` so that ``get_sa_coefficients()`` works
        (SA data will be empty).
        """
        self.logger.info(f'Running a PFR chain-of-reactors simulation using {self.__class__.__name__}...')

        species_names_list = [species.name for species in self.model.species()]
        self.all_data = list()
        self.distance_data = []

        dz = LENGTH / N_CELLS

        for condition in self.conditions:
            T0, P0, V0 = self._get_initial_state(condition)
            if V0 is None:
                self.model.TPX = T0, P0, condition.mol_frac
            elif P0 is None:
                self.model.TDX = T0, 1 / V0, condition.mol_frac

            total_residence_time = condition.reaction_time.value_si
            u_0 = LENGTH / total_residence_time
            mass_flow_rate = self.model.density * u_0 * AREA
            cell_volume = AREA * dz

            # Create the reactor network (reused across cells)
            reactor = ct.IdealGasReactor(self.model, energy='off', volume=cell_volume)
            upstream = ct.Reservoir(self.model, name='upstream')
            downstream = ct.Reservoir(self.model, name='downstream')
            # Cantera requires these objects to exist (they register with the reactor network)
            mfc = ct.MassFlowController(upstream, reactor, mdot=mass_flow_rate)  # noqa: F841
            ct.PressureController(reactor, downstream, primary=mfc, K=1e-12)
            sim = ct.ReactorNet([reactor])
            sim.atol = self.atol
            sim.rtol = self.rtol

            # Keep base-class attributes in sync so downstream code can access them
            self.cantera_reactor = reactor
            self.cantera_simulation = sim

            # Storage
            times = []
            temperature = []
            pressure = []
            species_data = []
            cumulative_time = 0.0
            distances = []

            for n in range(N_CELLS):
                if n > 0:
                    # Feed outlet of previous cell as inlet to next cell
                    upstream.phase.TPX = reactor.T, reactor.thermo.P, reactor.thermo.X
                    sim.reinitialize()

                sim.advance_to_steady_state()

                # Residence time in this cell
                cell_residence_time = reactor.mass / mass_flow_rate
                cumulative_time += cell_residence_time

                times.append(cumulative_time)
                temperature.append(reactor.T)
                pressure.append(reactor.thermo.P)
                species_data.append(reactor.thermo[species_names_list].X)
                distances.append((n + 1) * dz)

            # Convert to numpy arrays
            species_data = np.array(species_data)

            # Pack into GenericData objects (same format as base class)
            time_gd = GenericData(label='Time', data=times, units='s')
            temp_gd = GenericData(label='Temperature', data=temperature, units='K')
            pres_gd = GenericData(label='Pressure', data=pressure, units='Pa')
            condition_data = [temp_gd, pres_gd]

            for index, species in enumerate(self.model.species()):
                species_gd = GenericData(label=species.name,
                                         species=species,
                                         data=species_data[:, index],
                                         index=index)
                condition_data.append(species_gd)

            # No SA data for chain method (simulate() handles fallback)
            self.all_data.append((time_gd, condition_data, [], []))
            self.distance_data.append(np.array(distances))


register_simulate_adapter("CanteraPFR", CanteraPFR)
