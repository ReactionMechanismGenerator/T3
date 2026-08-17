"""
A module for generating mesolver adapters.
"""

from typing import TYPE_CHECKING

from t3.pdep.mesolver.adapter import MESolverAdapter

if TYPE_CHECKING:
    from t3.logger import Logger

_registered_mesolver_adapters = {}


def register_mesolver_adapter(me_solver: str,
                              mesolver_class: type[MESolverAdapter],
                              ) -> None:
    """
    A register for the mesolver adapters.

    Args:
        me_solver (str): A string representation for a mesolver adapter.
        mesolver_class (Type[MESolverAdapter]): The mesolver adapter class.

    Raises:
        TypeError: If ``mesolver_class`` is not a ``MESolverAdapter`` instance.
    """
    if not issubclass(mesolver_class, MESolverAdapter):
        raise TypeError(f'mesolver_class is not a MESolverAdapter, got {mesolver_class} which is a {type(mesolver_class)}')
    # Registering the same key twice silently overwrites the earlier entry. This deliberately
    # matches the behavior of ``t3/simulate/factory.py::register_simulate_adapter`` (a plain
    # dict assignment with no duplicate check) rather than inventing a stricter policy here.
    _registered_mesolver_adapters[me_solver] = mesolver_class


def mesolver_factory(me_solver: str,
                     network_path: str,
                     output_directory: str,
                     method: str,
                     logger: 'Logger | None' = None,
                     allow_ilt_complement: bool = False,
                     expected_reactions: int | set | None = None,
                     ) -> MESolverAdapter:
    """
    A factory generating the mesolver adapter corresponding to ``me_solver``.

    Args:
        me_solver (str): The mesolver adapter name. Examples: 'Arkane', 'MESS', or 'Mesmer'.
        network_path (str): The path to the pressure-dependent reaction network to solve.
        output_directory (str): The path to the directory in which to write the solver's
                                input/output files.
        method (str): The approximation method to use for solving the master equation
                      (e.g., 'CSE', 'MSC', 'RS').
        logger (Logger, optional): The current T3 Logger instance.
        allow_ilt_complement (bool, optional): Whether the network being solved is only
                                              partially QM-computed, with the remaining path
                                              reactions estimated via RMG/ILT.
        expected_reactions (int | set, optional): The expected coverage of solved-reaction
                                                 entries in the solver's output, passed through
                                                 to the adapter (see
                                                 ``ArkaneMESolverAdapter.__init__``). When
                                                 ``None`` (the default), output completeness is
                                                 not verified.

    Raises:
        ValueError: If the provided me_solver is not in the keys for the
                    _registered_mesolver_adapters dictionary, or if allow_ilt_complement is
                    True while the resolved adapter does not support solving a partially
                    QM-computed network.

    Returns:
        MESolverAdapter: The requested MESolverAdapter child, initialized with the respective
                         arguments.
    """

    if me_solver not in _registered_mesolver_adapters.keys():
        raise ValueError(f'The "me_solver" argument of {me_solver} was not present in the keys for the '
                         f'_registered_mesolver_adapters dictionary: {list(_registered_mesolver_adapters.keys())}'
                         f'\nPlease check that the mesolver adapter was registered properly.')

    mesolver_class = _registered_mesolver_adapters[me_solver]

    # NOTE: this allow_ilt_complement rule is currently enforced in the factory ONLY, so direct
    # construction of an adapter bypasses it. That is harmless today because the sole registered
    # adapter (Arkane) supports the ILT complement anyway; the check becomes load-bearing once a
    # non-Arkane adapter exists (MESS/Mesmer, a later increment), at which point it should move
    # into the adapters themselves (e.g., their __init__), not live here alone.
    if allow_ilt_complement and not mesolver_class.supports_ilt_complement:
        raise ValueError(f'The "{me_solver}" mesolver adapter does not support solving a network in which only a '
                         f'subset of wells and transition states carry QM data (allow_ilt_complement=True was '
                         f'requested). This solver requires the entire PES to be QM-calculated and cannot consume '
                         f'a partially-computed network.')

    adapter = mesolver_class(network_path=network_path,
                             output_directory=output_directory,
                             method=method,
                             logger=logger,
                             allow_ilt_complement=allow_ilt_complement,
                             expected_reactions=expected_reactions,
                             )
    return adapter
