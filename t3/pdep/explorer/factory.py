"""
A module for generating explorer adapters.
"""

from typing import TYPE_CHECKING

from t3.pdep.explorer.adapter import PESExplorerAdapter

if TYPE_CHECKING:
    from t3.logger import Logger

_registered_explorer_adapters = {}


def register_explorer_adapter(explorer: str,
                              explorer_class: type[PESExplorerAdapter],
                              ) -> None:
    """
    A register for the explorer adapters.

    Args:
        explorer (str): A string representation for an explorer adapter.
        explorer_class (Type[PESExplorerAdapter]): The explorer adapter class.

    Raises:
        TypeError: If ``explorer_class`` is not a ``PESExplorerAdapter`` instance.
    """
    if not issubclass(explorer_class, PESExplorerAdapter):
        raise TypeError(f'explorer_class is not a PESExplorerAdapter, got {explorer_class} which is a '
                        f'{type(explorer_class)}')
    # Registering the same key twice silently overwrites the earlier entry. This deliberately
    # matches the behavior of ``t3/pdep/mesolver/factory.py::register_mesolver_adapter`` (a plain
    # dict assignment with no duplicate check) rather than inventing a stricter policy here.
    _registered_explorer_adapters[explorer] = explorer_class


def explorer_factory(explorer: str,
                     seed_species,
                     output_directory: str,
                     bath_gas: dict = None,
                     explore_tol: float = None,
                     energy_tol: float = None,
                     flux_tol: float = None,
                     maximum_radical_electrons: int = None,
                     logger: 'Logger | None' = None,
                     transition_state_seeds: tuple = None,
                     database_kwargs: dict = None,
                     ) -> PESExplorerAdapter:
    """
    A factory generating the explorer adapter corresponding to ``explorer``.

    Unlike ``t3.pdep.mesolver.factory.mesolver_factory``, this factory ROUTES ONLY: it performs
    no seed/capability validation of its own. The seed/capability rules (``max_source_species``,
    ``supports_transition_state_seeds``) are enforced inside each concrete adapter's own
    ``__init__``, so no rule here is reachable-by-bypass via direct construction. This is a
    deliberate divergence from ``mesolver_factory``, whose ``allow_ilt_complement`` check lives
    in the factory only (see the note at ``t3/pdep/mesolver/factory.py:83-87``).

    Args:
        explorer (str): The explorer adapter name. Example: 'Arkane'.
        seed_species (list | tuple): The source (seed) species labels to explore from.
        output_directory (str): The path to the directory in which to write the explorer's
                                input/output files.
        bath_gas (dict, optional): The bath gas composition, mapping species labels to mole
                                  fractions.
        explore_tol (float, optional): The energy tolerance for exploring new isomers/reactions.
        energy_tol (float, optional): The energy tolerance for including a well/transition state
                                      in the output network.
        flux_tol (float, optional): The flux tolerance for including a well/transition state in
                                    the output network.
        maximum_radical_electrons (int, optional): The maximum number of radical electrons
                                                   allowed in an explored species.
        logger (Logger, optional): The current T3 Logger instance.
        transition_state_seeds (tuple, optional): Transition state label(s) to seed the
                                                  exploration from, in addition to (or instead
                                                  of) ``seed_species``. Only meaningful for an
                                                  adapter with ``supports_transition_state_seeds``
                                                  True.
        database_kwargs (dict, optional): Keyword arguments describing the RMG database
                                          settings to use for the exploration.

    Raises:
        ValueError: If the provided explorer is not in the keys for the
                    _registered_explorer_adapters dictionary.

    Returns:
        PESExplorerAdapter: The requested PESExplorerAdapter child, initialized with the
                            respective arguments.
    """

    if explorer not in _registered_explorer_adapters.keys():
        raise ValueError(f'The "explorer" argument of {explorer} was not present in the keys for the '
                         f'_registered_explorer_adapters dictionary: {list(_registered_explorer_adapters.keys())}'
                         f'\nPlease check that the explorer adapter was registered properly.')

    explorer_class = _registered_explorer_adapters[explorer]

    adapter = explorer_class(seed_species=seed_species,
                             output_directory=output_directory,
                             bath_gas=bath_gas,
                             explore_tol=explore_tol,
                             energy_tol=energy_tol,
                             flux_tol=flux_tol,
                             maximum_radical_electrons=maximum_radical_electrons,
                             logger=logger,
                             transition_state_seeds=transition_state_seeds,
                             database_kwargs=database_kwargs,
                             )
    return adapter
