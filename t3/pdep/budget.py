#!/usr/bin/env python3
# encoding: utf-8

"""
t3 pdep budget module

Choose which of the P-dep networks that qualified for QM refinement actually get it, when there is
a budget.

Qualification and affordability are different questions, and this module answers only the second.
:mod:`t3.pdep.selector` decides whether a network's rate coefficient is sensitive to kinetics that
are merely an estimate -- a per-network judgement about the physics. This module takes the whole
field of such judgements and, given a bound on how much quantum chemistry the user is willing to
buy this iteration, decides which of them fit. With no bound set, every qualified network fits and
this module only imposes an order.

Three properties are load-bearing:

* **It is pure.** It reads decisions and returns indices. Projecting a network's cost by actually
  queueing it would mutate T3's species, reactions and QM state, which is exactly what a budget
  must be able to decide against BEFORE committing to it.
* **It never slices a network.** A network that does not fit is skipped whole. Its transition
  states are ranked against each other nowhere -- ``PDepNetworkSelection.uncertain_ts_labels()``
  returns them SORTED BY LABEL, not by sensitivity -- so taking "the first k" of them would mean
  choosing quantum chemistry alphabetically. A half-computed network is also not half a result: the
  master-equation solve that consumes it is hybrid QM+ILT, and which transition states are QM is a
  physics decision owned by the selection, not by whatever happened to fit.
* **Its cost is an upper bound, never an under-estimate.** The cost charged for a network is its
  number of uncertain transition states. Some of those will not become ARC jobs at queue time (an
  unsafe label, missing structures, a transition state shared by several path reactions, a reaction
  T3 already knows). Charging for them anyway means the budget can only ever over-estimate the spend
  it authorizes -- the safe direction for a bound on an expensive resource.

Networks offered more than once in one iteration -- the normal case, since several sensitive
reactions can belong to a single network -- are charged once, for the UNION of what their offers
name. Two offers of one network are two answers to two different questions ("which transition states
is THIS observable reaction's rate sensitive to?"), so they can name different transition states; the
cost of admitting both is what queueing both actually produces, since
``T3.queue_pdep_transition_states`` de-duplicates against records it already wrote. For the same
reason such a network is ranked by the BEST of its offers: if any one of them found a strong rate
response, the network has that response.

A network that could not be evaluated is refused outright rather than ranked. It is not a cheap
candidate, it is an absent measurement -- there is nothing to rank it by and nothing to justify
spending on it. Being refused for that reason is recorded like any other refusal, never dropped.
"""

from dataclasses import dataclass, field

from t3.pdep.selector import EVALUATION_STATUS_NOT_EVALUATED, PDepNetworkSelection, selection_rank_key

__all__ = ['PDepBudgetDecision',
           'PDepBudgetSkip',
           'apply_pdep_qm_budget',
           'projected_ts_cost',
           ]


@dataclass(frozen=True)
class PDepBudgetSkip:
    """
    One network the budget refused, and why.

    A refusal is recorded rather than merely dropped: a network that qualified for QM refinement and
    then did not get it is a decision the run took, and it has to be reportable. Silently shrinking
    the work would read, downstream and in the logs, exactly like "there was nothing more to do".

    Attributes:
        network_id (str): The refused network's identifier.
        cost (int): The transition-state cost that was charged for it.
        remaining_transition_states (int | None): The transition-state budget still unspent when it
            was considered, or ``None`` if no transition-state budget was set.
        reason (str): A human-readable explanation naming the bound that refused it.
    """
    network_id: str
    cost: int
    remaining_transition_states: int | None
    reason: str


@dataclass(frozen=True)
class PDepBudgetDecision:
    """
    The outcome of applying a QM budget to a field of network decisions.

    Attributes:
        admitted_indices (tuple): Indices into the input sequence, in the order the networks should
            be queued: most deserving first. Every offer of an admitted network is present, so a
            network offered twice contributes two indices.
        skipped (tuple): One ``PDepBudgetSkip`` per distinct refused network, in the order the
            networks were considered.
        total_cost (int): The transition-state cost charged for everything admitted.
    """
    admitted_indices: tuple = field(default_factory=tuple)
    skipped: tuple = field(default_factory=tuple)
    total_cost: int = 0


def projected_ts_cost(selection: PDepNetworkSelection) -> int:
    """
    The transition-state cost to charge for refining one network.

    This is the number of distinct transition states the selection names as uncertain, which is an
    upper bound on the ARC jobs the network will produce: see the module docstring on why the bound
    is deliberately loose in that direction.

    Args:
        selection (PDepNetworkSelection): The network's decision.

    Returns:
        int: The number of uncertain transition states.
    """
    return len(selection.uncertain_ts_labels())


def _validate_budget(value, name: str) -> None:
    """
    Refuse a budget that is not ``None`` or a positive whole number.

    ``bool`` is rejected explicitly: it is a subclass of ``int``, so ``max_networks=True`` would
    otherwise silently mean "at most one network".

    Args:
        value: The budget to validate.
        name (str): The parameter name, for the error message.

    Raises:
        ValueError: If the value is neither ``None`` nor an ``int`` greater than zero.
    """
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f'{name} must be None (no limit) or a positive integer, got {value!r}.')


def apply_pdep_qm_budget(selections,
                         max_transition_states: int | None = None,
                         max_networks: int | None = None,
                         ) -> PDepBudgetDecision:
    """
    Rank qualified network decisions and admit as many as the budget allows.

    Networks are considered most-deserving-first, by :func:`t3.pdep.selector.selection_rank_key`,
    and admitted while they fit. A network that does not fit the remaining transition-state budget
    is skipped WHOLE and the walk continues, so a later, cheaper network can still be admitted --
    the alternative, stopping at the first network that does not fit, would leave budget unspent for
    no benefit. It never takes part of a network; see the module docstring.

    With both budgets ``None`` every network is admitted and only the order changes, which is what
    makes the budget opt-in: setting nothing refines exactly the networks T3 refined before.

    A network that no offer could evaluate is refused before it is ever ranked, and a network whose
    cost exceeds the ENTIRE transition-state budget is refused with a distinct reason, since no
    remaining budget can ever be larger than the whole of it and no future iteration will change
    that on its own.

    Args:
        selections: A sequence of ``PDepNetworkSelection`` objects to consider, in the order they
            were found. Repeats of one ``network_id`` are charged once, for the union of what they
            name, and ranked by the best of them (see the module docstring). A selection with no
            ``network_id`` is treated as its own network, since nothing identifies it as a repeat.
        max_transition_states (int, optional): The most transition states to admit this iteration,
            or ``None`` for no limit.
        max_networks (int, optional): The most networks to admit this iteration, or ``None`` for no
            limit.

    Returns:
        PDepBudgetDecision: The admitted indices, in ranked order, and the refusals.

    Raises:
        ValueError: If either budget is neither ``None`` nor a positive integer.
    """
    _validate_budget(max_transition_states, 'max_transition_states')
    _validate_budget(max_networks, 'max_networks')

    indices_by_network = dict()  # Keys are network identities, values are every index offering them.
    for index, selection in enumerate(selections):
        # An unnamed network cannot be recognized as the same network twice, so it is grouped alone.
        # Collapsing every unnamed decision into one pseudo-network would charge several distinct
        # networks as one and admit or refuse them together.
        identity = selection.network_id if selection.network_id else (None, index)
        indices_by_network.setdefault(identity, list()).append(index)

    # Rank each network by the BEST of its offers and charge it for the UNION of what they name; see
    # the module docstring on why two offers of one network are not the same decision twice.
    ranked = sorted(indices_by_network.items(),
                    key=lambda item: min(selection_rank_key(selections[index]) for index in item[1]))

    admitted_indices, skipped = list(), list()
    remaining = max_transition_states
    admitted_networks, total_cost = 0, 0
    for identity, indices in ranked:
        offers = [selections[index] for index in indices]
        network_id = offers[0].network_id
        cost = len({ts_label for selection in offers for ts_label in selection.uncertain_ts_labels()})
        if all(selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED for selection in offers):
            skipped.append(PDepBudgetSkip(network_id=network_id,
                                          cost=cost,
                                          remaining_transition_states=remaining,
                                          reason='it could not be evaluated, so there is no measured '
                                                 'sensitivity to justify spending the QM budget on it'))
            continue
        if max_transition_states is not None and cost > max_transition_states:
            # Distinguished from "does not fit right now" on purpose: no budget will ever remain
            # larger than the whole budget, so this network is refused every iteration, forever,
            # and the only thing that changes that is the user raising the limit.
            skipped.append(PDepBudgetSkip(network_id=network_id,
                                          cost=cost,
                                          remaining_transition_states=remaining,
                                          reason=f'it needs {cost} transition state(s), which exceeds the '
                                                 f'entire budget of {max_transition_states}, so it can '
                                                 f'never be refined until that limit is raised'))
            continue
        if max_networks is not None and admitted_networks >= max_networks:
            skipped.append(PDepBudgetSkip(network_id=network_id,
                                          cost=cost,
                                          remaining_transition_states=remaining,
                                          reason=f'the limit of {max_networks} network(s) per iteration '
                                                 f'was already reached'))
            continue
        if remaining is not None and cost > remaining:
            skipped.append(PDepBudgetSkip(network_id=network_id,
                                          cost=cost,
                                          remaining_transition_states=remaining,
                                          reason=f'it needs {cost} transition state(s) and only {remaining} '
                                                 f'remain(s) of the budget of {max_transition_states}'))
            continue
        admitted_indices.extend(indices)
        admitted_networks += 1
        total_cost += cost
        if remaining is not None:
            remaining -= cost
    return PDepBudgetDecision(admitted_indices=tuple(admitted_indices),
                              skipped=tuple(skipped),
                              total_cost=total_cost,
                              )
