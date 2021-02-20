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

    @abstractmethod
    def get_idt_by_T(self):
        """
        Finds the ignition point by approximating dT/dt as a first order forward difference
        and then finds the point of maximum slope.

        Returns:
            idt_dict (dict): Dictionary whose keys include 'idt' and 'idt_index' and whose values are lists of
                             the ignition delay time in seconds and index at which idt occurs respectively.
        """
        pass
