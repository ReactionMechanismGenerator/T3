"""
Cantera Simulator Adapter module for rapid compression machine (RCM) IDT calculations.

Models the post-compression state of an RCM as a closed, adiabatic,
constant-pressure batch reactor (``IdealGasConstPressureReactor``).
Heat loss to walls is not modeled. Inherits all IDT simulation and
sensitivity analysis capabilities from :class:`CanteraIDT`, overriding
only the reactor type.
"""

import cantera as ct

from t3.simulate.cantera_idt import CanteraIDT
from t3.simulate.factory import register_simulate_adapter


class CanteraRCM(CanteraIDT):
    """
    CanteraRCM is a SimulateAdapter that runs ignition-delay-time (IDT) simulations
    using a constant-pressure reactor (``IdealGasConstPressureReactor``), modelling
    rapid compression machine experiments.

    Inherits all functionality from :class:`CanteraIDT` — only the reactor type
    differs.
    """

    def _create_reactor(self, model: ct.Solution, energy: str = 'on'):
        """Create a constant-pressure reactor for RCM simulation."""
        return ct.IdealGasConstPressureReactor(model, energy=energy)


register_simulate_adapter("CanteraRCM", CanteraRCM)
