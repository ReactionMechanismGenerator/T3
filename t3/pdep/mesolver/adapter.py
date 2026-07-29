"""
t3 pdep mesolver adapter module.

This module defines the abstract MESolverAdapter class. This class allows users to easily solve a
given pressure-dependent reaction network with a master-equation solver (e.g., Arkane, MESS,
Mesmer) and obtain k(T,P) results.
"""

from abc import ABC, abstractmethod


class MESolverAdapter(ABC):
    """
    The abstract MESolverAdapter class.

    Attributes:
        supports_ilt_complement (bool): Whether this solver can solve a network in which only
                                        a subset of wells and transition states carry QM data
                                        while the remaining path reactions stay as RMG/ILT
                                        estimates. Arkane can (it does not require every node
                                        in the network to be QM-computed); MESS and Mesmer
                                        cannot, they require the entire PES to be QM-calculated.
    """

    supports_ilt_complement: bool = False

    @abstractmethod
    def set_up(self):
        """
        Set up the job directory and write the solver's input file for the given network.
        """
        pass

    @abstractmethod
    def solve(self):
        """
        Execute the solver and determine whether the solve genuinely succeeded.

        Each adapter is responsible for its own success criterion and must not assume that
        any other family's success criterion applies.

        Returns:
            bool: Whether the solve succeeded.
        """
        pass

    @abstractmethod
    def get_k_tp(self):
        """
        Obtain the solved k(T,P) results.

        Implementations return entries directed as the underlying tool wrote them, which is NOT
        guaranteed to match the direction the caller asked about. A consumer needing a specific
        direction must resolve it and reverse the entry itself. See the full reasoning on
        ``t3.pdep.explorer.arkane.ArkaneExplorerAdapter.get_k_tp``.

        Returns:
            k_tp: The solved k(T,P) results.
        """
        pass
