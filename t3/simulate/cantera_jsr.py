"""
Cantera Simulator Adapter module
Used to run mechanism analysis with Cantera as an ideal gas in a steady-state
jet stirred reactor (JSR), also known as a perfectly stirred reactor (PSR) or
continuously stirred tank reactor (CSTR).

A JSR is an open, well-mixed reactor with a continuous feed and outflow.  T3
models the JSR as isothermal and isobaric: feed gas enters at constant T and P,
the reactor contents are perfectly mixed at the same T and P, and the residence
time tau (set from the reactor's ``termination_time``) determines the mass flow
rate via ``mdot = m_reactor / tau``.

The reactor is approached transiently from the inlet state and integrated
forward in time over tau, recording the trajectory at every solver step.
Sensitivity analysis is supported and uses the standard ``CanteraBase``
machinery (the JSR is an ``IdealGasReactor`` with the energy equation off,
so the SA mass-fraction-to-mole-fraction correction applies unchanged).

Inspiration for the flow-network construction comes from
``set_jsr`` / ``run_jsr`` in ``t3/utils/flux.py`` (used for flux-based species
selection).  Users can modify the module-level ``VOLUME`` constant below to
change the JSR volume for their application.
"""

import cantera as ct

from t3.simulate.cantera_base import CanteraBase
from t3.simulate.factory import register_simulate_adapter


# ---------------------------------------------------------------------------
# Module-level constants — users modify these for their application
# ---------------------------------------------------------------------------
VOLUME = 1e-4              # m^3, JSR volume (= 100 cm^3, matching set_jsr default)
PRESSURE_COEFF = 0.01      # PressureController pressure_coeff: proportionality
                           # constant for pressure-difference flow correction
                           # (kg/s/Pa).  Larger -> faster pressure equilibration.


class CanteraJSR(CanteraBase):
    """
    Simulates an ideal gas in a steady-state jet stirred reactor (JSR / PSR / CSTR).

    The JSR is isothermal and isobaric: ``ct.IdealGasReactor`` with
    ``energy='off'``, fed by an upstream Reservoir through a
    ``MassFlowController`` and exhausted through a ``PressureController`` to a
    downstream Reservoir.  The mass flow rate is set so that the residence time
    of the reactor equals the input ``termination_time`` for each condition:

        mdot = m_reactor(0) / tau

    The base-class :meth:`simulate` time loop integrates from t = 0 to tau,
    capturing the transient approach to steady state and (optionally) the
    sensitivity coefficients at every time point.

    The reactor volume is controlled by the module-level :data:`VOLUME` constant.
    Because ``mdot`` is computed from the reactor mass, the steady-state
    composition is independent of ``VOLUME`` for a fixed inlet state and tau —
    only the absolute mass flow rate scales with volume.

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

    cantera_reactor_type = 'IdealGasReactor'

    def create_reactor(self):
        """
        Create a bare ``ct.IdealGasReactor`` matching the JSR reactor type.

        This is used by inherited helpers such as :meth:`find_equilibrium` that
        only need a reactor object (not the full flow network).  The actual JSR
        flow network is built in :meth:`reinitialize_simulation`.
        """
        return ct.IdealGasReactor(self.model, energy='off', volume=VOLUME)

    def reinitialize_simulation(self, T0=None, P0=None, X0=None, V0=None):
        """
        Build the JSR flow network for the next condition:

            inlet Reservoir -> MassFlowController -> JSR -> PressureController -> outlet Reservoir

        The mass flow rate is set so the residence time of the reactor equals
        the input ``termination_time`` for the current condition.  Sensitivity
        reactions and species enthalpies are registered when SA is requested.

        The residence time for the current condition is consumed from
        ``self._jsr_residence_time_iter``, which is initialized by
        :meth:`simulate` (or any other entry point that loops over
        ``self.conditions``).
        """
        # Set inlet conditions on self.model so that the inlet Reservoir
        # captures the feed state.
        if V0 is None:
            self.model.TPX = T0, P0, X0
        elif P0 is None:
            self.model.TDX = T0, 1 / V0, X0

        residence_time = next(self._jsr_residence_time_iter)

        # Build flow devices.  Reservoirs and devices are stashed on self
        # to keep them alive for the duration of the simulation (Cantera does
        # not retain Python references to flow-network objects).
        self._jsr_inlet = ct.Reservoir(self.model)
        self._jsr_exhaust = ct.Reservoir(self.model)
        self.cantera_reactor = ct.IdealGasReactor(self.model, energy='off', volume=VOLUME)

        if residence_time <= 0:
            raise ValueError(f'Invalid residence time: {residence_time}')
        self._jsr_mfc = ct.MassFlowController(
            upstream=self._jsr_inlet,
            downstream=self.cantera_reactor,
            mdot=self.cantera_reactor.mass / residence_time,
        )

        self._jsr_pc = ct.PressureController(
            upstream=self.cantera_reactor,
            downstream=self._jsr_exhaust,
            primary=self._jsr_mfc,
        )
        self._jsr_pc.pressure_coeff = PRESSURE_COEFF

        self.cantera_simulation = ct.ReactorNet([self.cantera_reactor])
        self.cantera_simulation.atol = self.atol
        self.cantera_simulation.rtol = self.rtol

        if self.sensitive_species:
            for i in range(self.num_ct_reactions):
                self.cantera_reactor.add_sensitivity_reaction(i)
            for i in range(self.num_ct_species):
                self.cantera_reactor.add_sensitivity_species_enthalpy(i)
            self.cantera_simulation.atol_sensitivity = self.sa_atol
            self.cantera_simulation.rtol_sensitivity = self.sa_rtol

    def _reset_residence_time_iterator(self):
        """Reset the residence-time iterator to the start of ``self.conditions``."""
        self._jsr_residence_time_iter = iter(
            [c.reaction_time.value_si for c in self.conditions]
        )

    def simulate(self):
        """
        Run the JSR transient simulation.

        Resets the residence-time iterator and delegates to the base-class
        :meth:`CanteraBase.simulate` for the per-condition time loop and SA
        collection.  The base-class loop calls :meth:`reinitialize_simulation`
        once per condition, which builds the JSR flow network described above.
        """
        self._reset_residence_time_iterator()
        super().simulate()

    def find_equilibrium(self, constrained_state_vars):
        """
        Compute equilibrium mole fractions for each condition (inherited
        behavior).  Equilibrium is a property of the gas state, not the reactor
        wiring, so the JSR flow network is built but unused.

        The residence-time iterator is reset before delegating to the base-class
        implementation to keep ``find_equilibrium`` callable independently of
        :meth:`simulate`.
        """
        self._reset_residence_time_iterator()
        return super().find_equilibrium(constrained_state_vars)

    def get_idt_by_T(self):
        """A steady-state isothermal JSR has no ignition delay."""
        return {'idt': list(), 'idt_index': list()}

    def get_t50(self, species, criteria='mass_frac'):
        """
        Half-life is not meaningful for a steady-state JSR with continuous feed:
        the species concentration approaches a steady-state value rather than
        decaying monotonically from its initial amount.

        Raises:
            NotImplementedError: Always.
        """
        raise NotImplementedError(
            'get_t50 is not meaningful for a JSR — the reactor is fed continuously '
            'so species concentrations approach a steady-state value rather than '
            'decaying monotonically. Use a closed batch reactor adapter '
            '(e.g. CanteraConstantTP) for half-life calculations.'
        )


register_simulate_adapter("CanteraJSR", CanteraJSR)
