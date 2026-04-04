"""
Cantera Simulator Adapter module
Used to run mechanism analysis with Cantera as ideal gases in a batch reactor at constant T and P
"""

import cantera as ct

from t3.simulate.cantera_base import CanteraBase
from t3.simulate.factory import register_simulate_adapter


class CanteraConstantTP(CanteraBase):
    """
    Simulates ideal gases in a batch reactor at constant temperature and pressure.
    Uses ``ct.IdealGasConstPressureReactor`` with ``energy='off'`` to maintain constant T.
    """

    cantera_reactor_type = 'IdealGasConstPressureReactor'

    def create_reactor(self):
        """Create a constant-pressure reactor with the energy equation disabled (constant T)."""
        return ct.IdealGasConstPressureReactor(self.model, energy='off')

    def get_idt_by_T(self):
        """
        Since this adapter simulates at constant T, there is no ignition.

        Returns:
            idt_dict (dict): Dictionary whose values are empty lists.
        """
        return {'idt': list(), 'idt_index': list()}


register_simulate_adapter("CanteraConstantTP", CanteraConstantTP)
