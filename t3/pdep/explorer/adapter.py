"""
t3 pdep explorer adapter module.

This module defines the abstract PESExplorerAdapter class. This class allows users to explore a
pressure-dependent potential energy surface (PES) starting from a seed of one or more source
species (e.g., via Arkane's ``explorer`` block) and obtain the resulting reaction network(s) and
their k(T,P) kinetics.
"""

from abc import ABC, abstractmethod


def refuse_bare_string_seed(field_name: str, value) -> None:
    """
    Refuse a bare ``str``/``bytes`` offered where a sequence of labels is expected.

    A module-level function rather than an inline check, because this rule MUST run before anything
    calls ``tuple()`` on the value, and there are now two places that do: ``validate_explorer_seed``
    below, and ``t3.pdep.explorer.config.PDepExplorerConfig``, which normalizes its seed fields to
    tuples at construction. Whichever of them runs first destroys the evidence for the other -- once
    ``'OH'`` has become ``('O', 'H')`` it is indistinguishable from a caller who deliberately asked
    for a bimolecular O + H source channel, so a later check cannot recover the distinction. Sharing
    one function keeps the refusal at every point of entry instead of only the first one that
    happened to be written.

    Args:
        field_name (str): The name of the field being checked, for the error message.
        value: The offered value.

    Raises:
        ValueError: If ``value`` is a ``str`` or ``bytes``.
    """
    if isinstance(value, (str, bytes)):
        raise ValueError(f"'{field_name}' must be a sequence of label strings, not a bare "
                         f"{type(value).__name__} ({value!r}). A string is itself a sequence, so this "
                         f"would silently be read character by character -- {tuple(value)!r} -- rather "
                         f"than as the single label it was meant to be.")


def validate_explorer_seed(seed_species, transition_state_seeds, max_source_species: int,
                           supports_transition_state_seeds: bool) -> tuple:
    """
    Apply the seed and capability rules, returning the normalized seed.

    A module-level function rather than a method, because it has two callers on two different entry
    paths and must be the same rule on both: ``PESExplorerAdapter.__init__`` (so direct construction
    is safe) and ``t3.pdep.explorer.factory.explorer_factory`` (so the blessed entry point is safe
    even for an adapter whose ``__init__`` never reached the base class -- see the Note on
    ``PESExplorerAdapter``). Being module-level is the point: a subclass cannot override it, so the
    factory's enforcement cannot be weakened by the very class it is checking. The capability values
    are passed in rather than read off a class, so the caller decides whose limits apply.

    Args:
        seed_species (list | tuple): The labels of the species forming the source configuration.
        transition_state_seeds (list | tuple): The labels of any transition states offered as seeds.
        max_source_species (int): The maximum number of source species this explorer allows.
        supports_transition_state_seeds (bool): Whether this explorer can seed from transition states.

    Raises:
        ValueError: If the seed is empty, has more than ``max_source_species`` entries, or carries
                   transition states when the explorer does not support them.

    Returns:
        tuple: ``(seed_species, transition_state_seeds)``, each normalized to a tuple.
    """
    # Refused BEFORE tuple(), because a str is a sequence and ``tuple('OH')`` is ``('O', 'H')``:
    # ``seed_species='OH'`` -- a plausible thing to write for a single seed -- silently became a
    # two-species bimolecular source channel, which passes the count rule below (2 <= 2) and, if the
    # network happens to declare species named 'O' and 'H', every later check as well. Arkane then
    # explores a different reaction than the one asked for. bytes is refused for the same reason.
    for name, value in (('seed_species', seed_species), ('transition_state_seeds', transition_state_seeds)):
        refuse_bare_string_seed(field_name=name, value=value)
    seed_species = tuple(seed_species or tuple())
    transition_state_seeds = tuple(transition_state_seeds or tuple())
    for name, values in (('seed_species', seed_species), ('transition_state_seeds', transition_state_seeds)):
        for value in values:
            if not isinstance(value, str) or not value:
                raise ValueError(f"Every entry of '{name}' must be a non-empty label string, got "
                                 f"{value!r} of type {type(value).__name__} in {values!r}.")

    # An EMPTY seed is refused separately from a too-large one. Arkane reaches the same else
    # branch for both and reports each as "Reactant channels with greater than 2 reactants not
    # supported" (arkane/explorer.py:159-169) -- a message that is simply false for an empty
    # seed, and would send a caller looking for a seed they never supplied.
    if not seed_species:
        raise ValueError(f'A PES exploration requires at least one source species, got '
                         f'{seed_species}.')
    if len(seed_species) > max_source_species:
        raise ValueError(f'A PES exploration source is a single reactant channel, so it may '
                         f'name at most {max_source_species} source species; got '
                         f'{len(seed_species)}: {seed_species}.')
    if transition_state_seeds and not supports_transition_state_seeds:
        raise ValueError(f'This explorer does not support transition-state-seeded exploration, '
                         f'yet {len(transition_state_seeds)} transition state(s) were given as '
                         f'seeds: {transition_state_seeds}.')

    return seed_species, transition_state_seeds


class PESExplorerAdapter(ABC):
    """
    The abstract PESExplorerAdapter class.

    Note:
        Unlike ``t3.pdep.mesolver.factory.mesolver_factory``'s ``allow_ilt_complement`` rule
        (which is enforced in the factory ONLY, so direct construction of a mesolver adapter
        bypasses it), the seed/capability rules below (``max_source_species`` and
        ``supports_transition_state_seeds``) are enforced HERE, by this base class's own
        ``__init__``, and not by ``t3.pdep.explorer.factory.explorer_factory``. The factory routes
        a name to a registered class only; it performs no capability checks of its own. This makes
        direct construction of an adapter just as safe as going through the factory.

        The enforcement deliberately lives in concrete base-class code rather than in a
        "every concrete adapter MUST validate this" instruction. Moving the check out of the
        factory closes a bypass, but leaving it unimplemented in the base would replace that
        bypass with something worse: a rule with NO enforcement site at all, honoured only by
        whichever adapter author remembers to re-implement it. Concrete adapters therefore call
        ``super().__init__(...)`` and inherit the rules instead of restating them.

        ``get_networks()`` and ``get_k_tp()`` MUST raise ``RuntimeError`` when called before a
        successful ``explore()`` (mirroring the contract already enforced by
        ``ArkaneMESolverAdapter.get_k_tp`` at ``t3/pdep/mesolver/arkane.py:179-183``). A caller
        must never silently receive empty or placeholder results as though they were the outcome
        of a genuine exploration.

    Attributes:
        max_source_species (int): The maximum number of source (seed) species an exploration may
                                  start from. Defaults to 2, Arkane's hard limit: Arkane's
                                  ``explorer`` block accepts a ``source`` of one or two species
                                  and raises for anything else (it has no concept of a
                                  disconnected seed set beyond a single bimolecular reactant
                                  pair).
        supports_transition_state_seeds (bool): Whether this explorer can seed an exploration
                                                from one or more transition states in addition
                                                to (or instead of) stable species. Defaults to
                                                False: Arkane's ``source`` is a list of
                                                ``species_dict`` labels only; a
                                                ``transitionState()`` label in ``source`` raises
                                                a ``KeyError``, so Arkane never supports this.
    """

    max_source_species: int = 2
    supports_transition_state_seeds: bool = False

    def __init__(self, seed_species=None, transition_state_seeds=None, **kwargs):
        """
        Validate and store the seed. Concrete adapters must call this via ``super().__init__()``.

        Args:
            seed_species (list | tuple, optional): The labels of the species forming the source
                                                   configuration to explore from.
            transition_state_seeds (list | tuple, optional): The labels of any transition states
                                                             offered as seeds. Refused unless the
                                                             adapter sets
                                                             ``supports_transition_state_seeds``.
            **kwargs: Ignored here; accepted so a concrete adapter can pass its own arguments
                      through a single ``super().__init__(...)`` call.

        Raises:
            ValueError: If the seed is empty, has more than ``max_source_species`` entries, or
                        carries transition states when the adapter does not support them.
        """
        # The rules live in ``validate_explorer_seed``, not inline here, because ``explorer_factory``
        # has to apply the identical rules on its own path -- a subclass that overrides ``__init__``
        # and forgets ``super().__init__(...)`` never reaches this line at all. Two call sites, one
        # rule; duplicating the checks would let the two entry points drift apart silently.
        self.seed_species, self.transition_state_seeds = validate_explorer_seed(
            seed_species=seed_species,
            transition_state_seeds=transition_state_seeds,
            max_source_species=self.max_source_species,
            supports_transition_state_seeds=self.supports_transition_state_seeds,
        )

    @abstractmethod
    def set_up(self):
        """
        Set up the job directory and write the explorer's input file for the given seed.

        MUST be called by ``explore()``, never from ``__init__``. Constructing an adapter is
        required to be side-effect-free: an explorer that judges its run directory before
        launching (to prove no stale artifact can be mistaken for a fresh result) cannot make that
        judgement honestly if construction has already written into the directory. Several
        adapters elsewhere in T3 do call ``set_up()`` from ``__init__``; that idiom is a known
        sharp edge there and must not be inherited here, where the ordering is load-bearing.
        """
        pass

    @abstractmethod
    def explore(self) -> bool:
        """
        Execute the exploration and determine whether it genuinely succeeded.

        Each adapter is responsible for its own success criterion and must not assume that any
        other family's success criterion applies.

        Returns:
            bool: Whether the exploration succeeded.
        """
        pass

    @abstractmethod
    def get_networks(self) -> tuple:
        """
        Obtain the explored network(s).

        Raises:
            RuntimeError: If ``explore()`` has not been called yet, or was called and did not
                         succeed. A caller must never silently receive empty results as if they
                         were the outcome of a genuine exploration.

        Returns:
            tuple: The explored network(s).
        """
        pass

    @abstractmethod
    def get_k_tp(self):
        """
        Obtain the parsed k(T,P) results for the explored network(s).

        Implementations return entries directed as the underlying tool wrote them, which is NOT
        guaranteed to match the direction the caller asked about. A consumer needing a specific
        direction must resolve it and reverse the entry itself. See the full reasoning on
        ``t3.pdep.explorer.arkane.ArkaneExplorerAdapter.get_k_tp``.

        Raises:
            RuntimeError: If ``explore()`` has not been called yet, or was called and did not
                         succeed. A caller must never silently receive empty results as if they
                         were the outcome of a genuine exploration.

        Returns:
            k_tp: The parsed k(T,P) results.
        """
        pass
