"""
A module for generating explorer adapters.
"""

from typing import TYPE_CHECKING

from t3.pdep.explorer.adapter import PESExplorerAdapter, validate_explorer_seed

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
                     network_path: str,
                     method: str,
                     bath_gas: dict = None,
                     explore_tol: float = None,
                     energy_tol: float = None,
                     flux_tol: float = None,
                     maximum_radical_electrons: int = None,
                     logger: 'Logger | None' = None,
                     transition_state_seeds: tuple = None,
                     database_kwargs: dict = None,
                     expected_source_hash: str = None,
                     timeout: float = None,
                     ) -> PESExplorerAdapter:
    """
    A factory generating the explorer adapter corresponding to ``explorer``.

    Unlike ``t3.pdep.mesolver.factory.mesolver_factory``, this factory ROUTES ONLY: it performs
    no seed/capability validation of its own. The seed/capability rules (``max_source_species``,
    ``supports_transition_state_seeds``) are enforced inside each concrete adapter's own
    ``__init__``, so no rule here is reachable-by-bypass via direct construction. That is still
    true and still the reason the rules live in the adapter, but it is not sufficient on its own:
    enforcement inside ``__init__`` is bypassed by any subclass that overrides ``__init__`` without
    calling ``super()``, and ``register_explorer_adapter`` cannot detect that (``issubclass`` remains
    True). So the seed rules are re-asserted HERE as well, against the arguments this call was given.
    The two checks guard two different entry paths -- direct construction and the factory -- and
    neither makes the other redundant. Compare ``mesolver_factory``, whose ``allow_ilt_complement``
    check lives in the factory only (see the note at ``t3/pdep/mesolver/factory.py:83-87``).

    Args:
        explorer (str): The explorer adapter name. Example: 'Arkane'.
        seed_species (list | tuple): The source (seed) species labels to explore from.
        output_directory (str): The path to the directory in which to write the explorer's
                                input/output files.
        network_path (str): The path to the RMG P-dep network file (or Arkane network input file)
                            to explore from.
        method (str): The master-equation method, e.g. 'CSE', 'MSC' or 'RS'.
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
        expected_source_hash (str, optional): The content hash (``t3.pdep.hashing`` format)
                                              ``network_path``'s bytes must match when the adapter's
                                              ``set_up()`` reads them, forwarded verbatim to the
                                              adapter's constructor.
        timeout (float, optional): The wall-clock deadline, in seconds, for the explorer's
                                   underlying process, forwarded verbatim to the adapter's
                                   constructor (``None`` means no deadline; validated and enforced
                                   at the runner layer, see
                                   ``t3.runners.rmg_runner.run_arkane_job``).

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

    # Re-asserted here, against the arguments THIS call was given, rather than trusted to the
    # adapter. Registration checks only ``issubclass``, which is a claim about ancestry and not about
    # whether the class's ``__init__`` ever reached ``PESExplorerAdapter.__init__``: a subclass that
    # overrides ``__init__`` and forgets ``super().__init__(...)`` inherits every seed rule and
    # enforces none of them, and ``issubclass`` is still True, so registration blesses it. The shared
    # module-level rule function is used deliberately -- it cannot be overridden by the class being
    # checked, which a method or a classmethod could be.
    # The return value is kept, not discarded: it is the NORMALIZED pair (both tuples), which is
    # exactly what a well-behaved adapter stores, and so it is what the post-condition below compares
    # against. Today that is not a behavioural difference -- normalization is ``tuple(x or tuple())``
    # (adapter.py:50-51), so re-deriving it here from the raw arguments would give the same answer, and
    # a mutation that does exactly that is indistinguishable by test. It is kept anyway to hold ONE
    # source of truth: the day normalization gains a step (dedup, case-folding, ordering), a re-derived
    # copy here would silently disagree with what the adapter actually stored, and this check would
    # start refusing well-behaved adapters or, worse, accepting forgetful ones.
    expected_seed, expected_ts_seeds = validate_explorer_seed(
        seed_species=seed_species,
        transition_state_seeds=transition_state_seeds,
        max_source_species=explorer_class.max_source_species,
        supports_transition_state_seeds=explorer_class.supports_transition_state_seeds,
    )

    # Constructed with the NORMALIZED tuples, not the caller's original objects. The seed is handed
    # onward twice -- once to the validation above, once to the constructor, whose own
    # super().__init__ validates it again -- and a one-shot iterable does not survive that: the
    # generator arrives exhausted and the adapter reports 'requires at least one source species'
    # about a seed that was demonstrably supplied. Beyond the message, the two calls have to agree on
    # what was validated, and they cannot if the first one consumes it.
    adapter = explorer_class(seed_species=expected_seed,
                             output_directory=output_directory,
                             network_path=network_path,
                             method=method,
                             bath_gas=bath_gas,
                             explore_tol=explore_tol,
                             energy_tol=energy_tol,
                             flux_tol=flux_tol,
                             maximum_radical_electrons=maximum_radical_electrons,
                             logger=logger,
                             transition_state_seeds=expected_ts_seeds,
                             database_kwargs=database_kwargs,
                             expected_source_hash=expected_source_hash,
                             timeout=timeout,
                             )

    # The quiet half of the same defect: re-running the rules above only catches a forgetful adapter
    # when the seed happens to VIOLATE one. A valid seed that the adapter simply never stored passes
    # every rule and yields an object with no ``seed_species`` at all, which surfaces much later as an
    # AttributeError inside explore() -- far from the cause, and looking like a bug in the exploration
    # rather than in the adapter's constructor. So the outcome is checked, not just the inputs.
    #
    # BOTH fields are checked, against everything ``validate_explorer_seed`` returned, rather than
    # ``seed_species`` alone. An adapter that records the seed and drops the transition-state seeds
    # fails no rule and passes a seed-only check, and then runs an ordinary unseeded exploration --
    # a different and cheaper calculation than the one requested, completing successfully, so nothing
    # ever raises and the caller believes the result was TS-seeded. Enumerating one field made this a
    # statement about that field; checking the returned pair makes it a statement about the object.
    for field_name, expected in (('seed_species', expected_seed), ('transition_state_seeds', expected_ts_seeds)):
        found = getattr(adapter, field_name, None)
        if found != expected:
            raise TypeError(f'{explorer_class.__name__} did not record the {field_name} it was constructed '
                            f'with: expected {field_name} {expected!r}, found '
                            f'{getattr(adapter, field_name, "<unset>")!r}. A concrete PESExplorerAdapter must '
                            f'call super().__init__(...), which is what stores and validates the seed.')
    return adapter
