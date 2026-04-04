"""
Cantera Simulator Adapter module
Used to run mechanism analysis with Cantera as an ideal gas in a batch reactor at constant V-U (adiabatic)
"""

import cantera as ct

from t3.simulate.cantera_base import CanteraBase
from t3.simulate.factory import register_simulate_adapter


class CanteraConstantUV(CanteraBase):
    """
    Simulates ideal gases in a batch reactor adiabatically at constant volume.
    """

    cantera_reactor_type = 'IdealGasReactor'

    def create_reactor(self):
        """Create a constant-volume ideal gas reactor (adiabatic)."""
        return ct.IdealGasReactor(self.model)


register_simulate_adapter("CanteraConstantUV", CanteraConstantUV)
