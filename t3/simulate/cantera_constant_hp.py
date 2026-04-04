"""
Cantera Simulator Adapter module
Used to run mechanism analysis with Cantera as an ideal gas in a batch reactor at constant H-P
"""

import cantera as ct

from t3.simulate.cantera_base import CanteraBase
from t3.simulate.factory import register_simulate_adapter


class CanteraConstantHP(CanteraBase):
    """
    Simulates ideal gases in a batch reactor at constant pressure and enthalpy (adiabatic, constant P).
    """

    cantera_reactor_type = 'IdealGasConstPressureReactor'

    def create_reactor(self):
        """Create a constant-pressure reactor with the energy equation enabled (adiabatic)."""
        return ct.IdealGasConstPressureReactor(self.model)


register_simulate_adapter("CanteraConstantHP", CanteraConstantHP)
