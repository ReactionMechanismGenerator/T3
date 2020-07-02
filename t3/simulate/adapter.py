"""
A module for the abstract SimulateAdapter class.
This class allows users to easily simulate a mechanism and perform analysis with RMG, RMS, or Cantera.
"""

from abc import ABC, abstractmethod


class SimulateAdapter(ABC):
    """
    An abstract Simulate Adapter class
    """

    @abstractmethod
    def set_up(self):
        """
        Simulate a job to obtain species profiles. Run sensitivity analysis if requested.
        """
        pass

    @abstractmethod
    def get_sa_coefficients(self):
        """
        Return sensitivity analysis coefficients in a standard dictionary format.
        """
        pass
