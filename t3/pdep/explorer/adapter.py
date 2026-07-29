"""
t3 pdep explorer adapter module.

This module defines the abstract PESExplorerAdapter class. This class allows users to explore a
pressure-dependent potential energy surface (PES) starting from a seed of one or more source
species (e.g., via Arkane's ``explorer`` block) and obtain the resulting reaction network(s) and
their k(T,P) kinetics.
"""

from abc import ABC, abstractmethod


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
        seed_species = tuple(seed_species or tuple())
        transition_state_seeds = tuple(transition_state_seeds or tuple())

        # An EMPTY seed is refused separately from a too-large one. Arkane reaches the same else
        # branch for both and reports each as "Reactant channels with greater than 2 reactants not
        # supported" (arkane/explorer.py:159-169) -- a message that is simply false for an empty
        # seed, and would send a caller looking for a seed they never supplied.
        if not seed_species:
            raise ValueError(f'A PES exploration requires at least one source species, got '
                             f'{seed_species}.')
        if len(seed_species) > self.max_source_species:
            raise ValueError(f'A PES exploration source is a single reactant channel, so it may '
                             f'name at most {self.max_source_species} source species; got '
                             f'{len(seed_species)}: {seed_species}.')
        if transition_state_seeds and not self.supports_transition_state_seeds:
            raise ValueError(f'This explorer does not support transition-state-seeded exploration, '
                             f'yet {len(transition_state_seeds)} transition state(s) were given as '
                             f'seeds: {transition_state_seeds}.')

        self.seed_species = seed_species
        self.transition_state_seeds = transition_state_seeds

    @abstractmethod
    def set_up(self):
        """
        Set up the job directory and write the explorer's input file for the given seed.
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

        Raises:
            RuntimeError: If ``explore()`` has not been called yet, or was called and did not
                         succeed. A caller must never silently receive empty results as if they
                         were the outcome of a genuine exploration.

        Returns:
            k_tp: The parsed k(T,P) results.
        """
        pass
