"""
A module for the abstract SimulateAdapter class.
This class allows users to easily simulate a mechanism and perform analysis with external packages such as
RMG, RMS, or Cantera.
"""

from abc import ABC, abstractmethod


class SimulateAdapter(ABC):
    """
    An abstract Simulate Adapter class
    """

    @abstractmethod
    def set_up(self):
        """
        Set up the class, possibly by reading an input file, initializing values, and/or setting up the domain.
        """
        pass

    @abstractmethod
    def simulate(self):
        """
        Simulate a job to obtain species profiles. Run sensitivity analysis if requested.
        """
        pass

    @abstractmethod
    def get_sa_coefficients(self):
        """
        Obtain the sensitivity analysis coefficients.

        Returns:
             sa_dict (dict): a SA dictionary, whose structure is given in the docstring for T3/t3/main.py
        """
        pass
