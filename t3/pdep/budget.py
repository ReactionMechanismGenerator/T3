#!/usr/bin/env python3
# encoding: utf-8

"""
t3.pdep.budget module

Once P-dep networks have been qualified for QM refinement, this module decides which of them the
QM budget actually pays for.

Qualification and budgeting are deliberately separate questions. :mod:`t3.pdep.selector` decides
whether a network's rate coefficient sensitivity estimate is trustworthy -- that is per-network
physics judgement. This module decides how many of the qualified networks the caller is willing
to spend quantum chemistry on this iteration, and which ones fit -- that is a resource-allocation
decision, and it should not know or care why any one network was qualified, only how much it
costs and how it ranks against the others.

Three things are load-bearing here:

**This module stays pure.** It only makes decisions; it never touches disk, mutates T3's species
or reactions state, or queues QM work. Everything downstream -- whether to actually persist a
decision, or act on it -- is decided by the caller, never here.

**A cost is charged once per network, not per transition state.** The cost charged to a network
is the number of its uncertain transition states (see ``PDepNetworkSelection.uncertain_ts_labels()``,
which already sorts and de-duplicates by label), not a sensitivity-weighted count and not the
number of QM jobs it will eventually produce -- half-computing a network's cost from only its
most sensitive path would let a network with one very sensitive transition state and nine
merely-uncertain ones masquerade as ten times cheaper than it is. Whatever hybrid QM+ILT physics
decision consumes this cost later owns the question of which transition states within an admitted
network actually get queued; that decision has not happened here and this module does not
anticipate it.

**Cost is an upper bound, never an under-estimate.** The cost charged to a network is a ceiling on
the QM work it could produce, not a promise of what ``T3.queue_pdep_transition_states`` will
actually queue -- that step can queue fewer jobs (unsafe labels, missing structures, a transition
state already shared with an admitted reaction) but never more. Charging the budget on the
over-estimate means the budget can only under-spend against its stated limit, never overspend it.

Networks can be offered more than once in one iteration -- e.g. because more than one sensitive
reaction belongs to the same underlying network -- and a repeat offer of the same ``network_id``
is charged once, for the UNION of transition-state labels across all its offers, not once per
offer. Two offers of one network answer the same question ("which transition states is THIS
observable reaction's rate sensitive to?"), not two independent questions, so charging their union
once is the honest cost of the queueing ``T3.queue_pdep_transition_states`` actually produces --
it already de-duplicates a transition state shared by more than one offer of a network. A network
offered more than once is ranked by the BEST of its offers, since one strong sensitivity response
justifies the spend regardless of how many weaker responses also named the network.

A network that could not be evaluated is still ranked -- ``selection_rank_key`` gives it its own
tier, neither trusted nor dismissed, between qualified and evaluated-and-unqualified -- but it is
refused inside this module's walk, immediately once reached, rather than admitted: it is not a
cheap candidate that lost a fair fight for space against the budget, it is an absent measurement,
so there is nothing to justify spending QM time on it. Its persisted rank records only where it
sat in the walk, not that it was preferred or disfavored. Being refused for that reason is
recorded like any other refusal, never silently dropped.
"""

import os
from dataclasses import dataclass, field

from t3.pdep.selector import EVALUATION_STATUS_NOT_EVALUATED, PDepNetworkSelection, selection_rank_key


__all__ = ['BUDGET_ALGORITHM_VERSION',
           'BUDGET_OUTCOME_ADMITTED',
           'BUDGET_OUTCOME_REFUSED',
           'BUDGET_RECORD_FILE_NAME',
           'BUDGET_RECORD_SCHEMA_VERSION',
           'BUDGET_SKIP_DOES_NOT_FIT_REMAINING',
           'BUDGET_SKIP_EXCEEDS_BUDGET',
           'BUDGET_SKIP_MAX_NETWORKS_REACHED',
           'BUDGET_SKIP_NOT_EVALUATED',
           'PDepBudgetConsideration',
           'PDepBudgetDecision',
           'PDepBudgetNetworkOutcome',
           'PDepBudgetRecord',
           'PDepBudgetSkip',
           'VALID_BUDGET_OUTCOMES',
           'VALID_BUDGET_SKIP_REASON_CODES',
           'apply_pdep_qm_budget',
           'budget_record_path',
           'build_pdep_budget_record',
           'projected_ts_cost',
           ]


# A machine-readable code for why a network was skipped, so PDepBudgetSkip.reason's prose can be
# reworded freely without breaking anything that switches on the CODE. These are exactly
# apply_pdep_qm_budget's four refusal sites, in the order the walk can reach them.
BUDGET_SKIP_NOT_EVALUATED = 'not_evaluated'
BUDGET_SKIP_EXCEEDS_BUDGET = 'exceeds_budget'
BUDGET_SKIP_MAX_NETWORKS_REACHED = 'max_networks_reached'
BUDGET_SKIP_DOES_NOT_FIT_REMAINING = 'does_not_fit_remaining'
VALID_BUDGET_SKIP_REASON_CODES = (BUDGET_SKIP_NOT_EVALUATED,
                                  BUDGET_SKIP_EXCEEDS_BUDGET,
                                  BUDGET_SKIP_MAX_NETWORKS_REACHED,
                                  BUDGET_SKIP_DOES_NOT_FIT_REMAINING,
                                  )

# The two things a network's budget outcome can be. A network is either admitted (the budget
# authorizes it to be queued this iteration -- queueing itself can still yield zero ARC jobs, for
# the same over-estimate reasons the module docstring gives for cost) or refused (it will not be
# queued this iteration) -- there is no third state, so PDepBudgetNetworkOutcome.__post_init__
# enforces that admitted and refused are mutually exclusive with carrying refusal detail.
BUDGET_OUTCOME_ADMITTED = 'admitted'
BUDGET_OUTCOME_REFUSED = 'refused'
VALID_BUDGET_OUTCOMES = (BUDGET_OUTCOME_ADMITTED, BUDGET_OUTCOME_REFUSED)

# ONE version for the SHAPE PDepBudgetRecord.as_dict() serializes -- not what the record means, only
# what it looks like on disk. Bump this when a field is added, removed, renamed, or changes type;
# do not bump it for a change to the admission/refusal LOGIC, which is what BUDGET_ALGORITHM_VERSION
# below is for. Sibling on-disk-shape versions elsewhere: the selection this record's networks got
# in (SELECTION_SCHEMA_VERSION, in t3.pdep.selector), the SA cache this record's costs were computed
# from (SA_CACHE_CONTRACT_VERSION, in t3.pdep.cache), and the exploration result this budget's
# networks may eventually produce (EXPLORATION_RESULT_SCHEMA_VERSION, in t3.pdep.explorer.result).
# Nothing here migrates an old record forward; a version mismatch is a signal to re-derive, not
# to guess -- so pre-ship, these constants should still be bumped freely rather than treated as frozen.
BUDGET_RECORD_SCHEMA_VERSION = 1

# The SEMANTICS of the budgeting decision: which refusal/admission rules apply and in what order
# (not-evaluated refused before ranking is trusted, then exceeds-the-whole-budget, then
# max-networks-reached, then does-not-fit-remaining; union-charging and best-of-offers ranking for
# a network offered more than once). Bump this when that LOGIC changes in a way that could flip a
# past decision, so a saved budget record can be told apart from one a newer apply_pdep_qm_budget
# would produce given the same inputs. Do NOT bump this for a change to the on-disk SHAPE of a
# record (see BUDGET_RECORD_SCHEMA_VERSION above) -- that does not change what the decision means.
BUDGET_ALGORITHM_VERSION = 1

# The budget record is written under the ITERATION directory (not the ARC project directory, unlike
# t3.pdep.join's TS_JOIN_SIDECAR_FILE_NAME): it describes a decision this iteration took over EVERY
# qualified network, admitted or refused, so it belongs beside the iteration's other top-level
# artifacts rather than nested inside the ARC run that only some admitted networks ever reach.
BUDGET_RECORD_FILE_NAME = 't3_pdep_qm_budget.yml'


def budget_record_path(iteration_directory: str) -> str:
    """
    Return the path a ``PDepBudgetRecord`` for one iteration is (or would be) written to.

    Args:
        iteration_directory (str): The T3 iteration directory (``self.paths['iteration']``).

    Returns:
        str: ``iteration_directory``, joined with ``BUDGET_RECORD_FILE_NAME``.
    """
    return os.path.join(iteration_directory, BUDGET_RECORD_FILE_NAME)


@dataclass(frozen=True)
class PDepBudgetSkip:
    """
    One network the QM budget refused to admit this iteration.

    Recorded so a refusal is legible downstream -- in logs, in the durable record -- rather than
    the network simply vanishing with "the budget did exactly nothing to it".

    Attributes:
        network_id (str, optional): The network's identifier, or a falsy value for an unnamed network.
        cost (int): The transition-state cost charged to it.
        remaining_transition_states (int, optional): The transition-state budget still available
            immediately BEFORE this network was considered, or ``None`` if no transition-state
            budget was configured.
        reason_code (str): One of ``VALID_BUDGET_SKIP_REASON_CODES``, machine-readable.
        reason (str): A human-readable explanation, not bound to stay stable across versions.
    """
    network_id: str | None
    cost: int
    remaining_transition_states: int | None
    reason_code: str
    reason: str

    def __post_init__(self):
        if self.reason_code not in VALID_BUDGET_SKIP_REASON_CODES:
            raise ValueError(f'PDepBudgetSkip.reason_code must be one of {VALID_BUDGET_SKIP_REASON_CODES}, '
                             f'got {self.reason_code!r}.')


@dataclass(frozen=True)
class PDepBudgetConsideration:
    """
    One network's full passage through apply_pdep_qm_budget's single walk, computed there once
    and consumed by build_pdep_budget_record -- never recomputed from ``selections`` a second
    time, so a durable record can never drift from, or be handed a mismatched pair with, the
    decision it is meant to describe.

    Attributes:
        identity: The grouping key this network's offers were grouped under during the walk --
            either the shared, truthy ``network_id``, or ``(None, first_offer_index)`` for an
            unnamed network. Not ``network_id`` alone, because an unnamed network's ``network_id``
            (``None`` or ``''``) is not unique across distinct unnamed networks.
        network_id (str, optional): The network's identifier, or a falsy value if unnamed.
        offer_indices (tuple): Every index into the walk's ``selections`` input that offered this
            network, in the order offered.
        cost (int): The transition-state cost charged to it (the union of uncertain transition
            state labels across all its offers).
        rank (int): Position in the most-deserving-first walk, 0 first.
        remaining_before (int, optional): Transition-state budget remaining immediately BEFORE this
            network was considered, or ``None`` if no transition-state budget was configured.
        outcome (str): One of ``VALID_BUDGET_OUTCOMES``.
        reason_code (str, optional): One of ``VALID_BUDGET_SKIP_REASON_CODES`` if refused, else
            ``None``.
        reason (str, optional): Human-readable, if refused, else ``None``.
        network_source_hash (str, optional): This network's ``network_source_hash``, resolved
            across its offers the same way ``PDepNetworkSelection.combine()`` resolves it: ``None``
            if its offers disagree or only some recorded one, rather than silently binding the
            aggregate to one offer's bytes.
        method (str, optional): This network's ``method``, resolved across its offers the same way
            ``combine()`` resolves ``method``: ``None`` on disagreement, since it is a weaker signal
            than identity.
    """
    identity: object
    network_id: str | None
    offer_indices: tuple
    cost: int
    rank: int
    remaining_before: int | None
    outcome: str
    reason_code: str | None = None
    reason: str | None = None
    network_source_hash: str | None = None
    method: str | None = None


@dataclass(frozen=True)
class PDepBudgetDecision:
    """
    The outcome of applying the QM budget to one iteration's field of qualified networks.

    Attributes:
        admitted_indices (tuple): Indices into the input ``selections`` list that the budget
            admits: an admitted network contributes every one of its offering indices, since a
            repeat offer of one network is charged and ranked once but still identifies each
            reaction that named it.
        skipped (tuple): One ``PDepBudgetSkip`` per refused network, in the order the walk refused
            them -- not necessarily rank order, since refusal happens inline as the walk reaches
            each network.
        total_cost (int): The summed transition-state cost of every admitted network.
        considered (tuple): One ``PDepBudgetConsideration`` per distinct network considered, in
            most-deserving-first order -- the single authoritative walk's own record of identity,
            offers, cost, rank, remaining-before-consideration, outcome, and (for a refusal)
            reason, plus the network_source_hash/method this same walk resolved.
            ``build_pdep_budget_record`` consumes ONLY this field (and the two below) rather than
            re-deriving grouping, ranking, or remaining-budget bookkeeping a second time from
            ``selections``, so the two can never drift apart or be handed a mismatched pair.
        max_transition_states (int, optional): The transition-state budget this decision was
            actually computed under, carried on the decision itself so a consumer never has to
            pass the same budget configuration to two different functions and trust they agree.
        max_networks (int, optional): The network-count budget this decision was actually computed
            under, for the same reason.
    """
    admitted_indices: tuple = field(default_factory=tuple)
    skipped: tuple = field(default_factory=tuple)
    total_cost: int = 0
    considered: tuple = field(default_factory=tuple)
    max_transition_states: int | None = None
    max_networks: int | None = None


def projected_ts_cost(selection: PDepNetworkSelection) -> int:
    """
    The transition-state cost a network's own selection would charge the budget, standing alone.

    ARC job counts are deliberately NOT part of this: cost is a bound on the QM work a network
    could produce, not a promise of what will actually be queued (see the module docstring).

    Args:
        selection (PDepNetworkSelection): The network's qualification decision.

    Returns:
        int: The number of distinct uncertain transition-state labels this selection names.
    """
    return len(selection.uncertain_ts_labels())


def _validate_budget(value, name: str) -> None:
    """
    Validate one of the two optional budget limits (``max_transition_states``/``max_networks``).

    ``None`` means no limit -- that must stay distinguishable from ``0``, which would mean "admit
    nothing", so this rejects a bool (a subclass of ``int`` in Python -- ``max_networks=True`` must
    not silently become "admit exactly one network").

    Args:
        value: The value to validate.
        name (str): The parameter name, for the error message.

    Raises:
        ValueError: If ``value`` is neither ``None`` nor a positive, non-bool ``int``.
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
    Decide which qualified networks the QM budget admits this iteration.

    Networks are considered most-deserving-first, by :func:`t3.pdep.selector.selection_rank_key`,
    and an admitted network never yields its transition-state budget back even if a later,
    cheaper network would not fit in what it left behind -- once a network is admitted the WHOLE
    of its cost is committed. A network that does not fit is refused, not queued partially: there
    is no notion of admitting half a network's transition states here.

    Both budgets default to ``None``, meaning no limit -- opt-in budgeting, not opt-out.

    Refusal reasons, in the order the walk can reach them: not evaluated (refused before ranking
    is even trusted); cost exceeds the ENTIRE transition-state budget (a distinct reason from
    merely not fitting in what remains, since no future remaining-budget could ever admit it,
    this iteration or otherwise); the network-count limit already reached; or the network's cost
    does not fit in what is left of the transition-state budget.

    Args:
        selections: An iterable of ``PDepNetworkSelection`` objects, in the order they were found.
            Repeats of one ``network_id`` are charged once, for the union of offers by that name,
            and ranked by the best of those offers. A falsy ``network_id`` is treated as its own,
            unmergeable network, since nothing identifies it as the same network as another.
        max_transition_states (int, optional): Transition-state budget to admit this iteration, or
            ``None`` for no limit.
        max_networks (int, optional): Number of networks to admit this iteration, or ``None`` for
            no limit.

    Returns:
        PDepBudgetDecision: Admitted indices, in ranked order, skips, and this walk's own detail.

    Raises:
        ValueError: If either budget is not ``None`` or a positive integer, or if a repeated
            network's own offers disagree on a real (non-``None``) ``network_source_hash``.
    """
    _validate_budget(max_transition_states, 'max_transition_states')
    _validate_budget(max_networks, 'max_networks')

    # Group offers by network identity, index by index, so a repeated network's several offers are
    # considered exactly once, together.
    indices_by_network = dict()
    for index, selection in enumerate(selections):
        # A falsy network_id cannot be recognized as the same network twice, so each unnamed
        # offer is its own, singleton pseudo-network -- never merged with another unnamed network.
        identity = selection.network_id if selection.network_id else (None, index)
        indices_by_network.setdefault(identity, list()).append(index)

    # Each network is ranked by the BEST of its offers, and charged the UNION of their costs.
    ranked = sorted(indices_by_network.items(),
                    key=lambda item: min(selection_rank_key(selections[index]) for index in item[1]))

    admitted_indices, skipped, considered = list(), list(), list()
    remaining = max_transition_states
    admitted_networks, total_cost = 0, 0
    for rank, (identity, indices) in enumerate(ranked):
        offers = [selections[index] for index in indices]
        network_id = offers[0].network_id
        cost = len({ts_label for selection in offers for ts_label in selection.uncertain_ts_labels()})
        remaining_before = remaining

        # network_source_hash/method mirror PDepNetworkSelection.combine()'s identity/weaker-signal
        # rules exactly, since these several offers of one network are precisely what combine()
        # exists to fold into one: a real hash disagreement means these offers do not actually
        # describe the same network at all (fail closed, raise), while a hash partially missing, or
        # a method disagreement, is a weaker signal recorded as None rather than silently bound to
        # whichever offer happened to be listed first.
        source_hashes = {selection.network_source_hash for selection in offers}
        known_source_hashes = {value for value in source_hashes if value is not None}
        if len(known_source_hashes) > 1:
            raise ValueError(f'Cannot budget offers of network {network_id!r} computed from different '
                             f'revisions of the network: {sorted(known_source_hashes)}.')
        network_source_hash = offers[0].network_source_hash
        if len(source_hashes) > 1:
            network_source_hash = None
        methods = {selection.method for selection in offers}
        method = offers[0].method
        if len(methods) > 1:
            method = None

        if all(selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED for selection in offers):
            skip = PDepBudgetSkip(network_id=network_id,
                                  cost=cost,
                                  remaining_transition_states=remaining,
                                  reason_code=BUDGET_SKIP_NOT_EVALUATED,
                                  reason='it could not be evaluated, so there is no measured '
                                         'sensitivity to justify spending the QM budget on it')
            skipped.append(skip)
            considered.append(PDepBudgetConsideration(identity=identity, network_id=network_id,
                                                      offer_indices=tuple(indices), cost=cost, rank=rank,
                                                      remaining_before=remaining_before,
                                                      outcome=BUDGET_OUTCOME_REFUSED,
                                                      reason_code=skip.reason_code, reason=skip.reason,
                                                      network_source_hash=network_source_hash, method=method))
            continue
        if max_transition_states is not None and cost > max_transition_states:
            skip = PDepBudgetSkip(network_id=network_id,
                                  cost=cost,
                                  remaining_transition_states=remaining,
                                  reason_code=BUDGET_SKIP_EXCEEDS_BUDGET,
                                  reason=f'its cost ({cost}) exceeds the entire transition-state '
                                         f'budget ({max_transition_states}), not just what remains, so it can '
                                         f'never be admitted until the limit itself is raised')
            skipped.append(skip)
            considered.append(PDepBudgetConsideration(identity=identity, network_id=network_id,
                                                      offer_indices=tuple(indices), cost=cost, rank=rank,
                                                      remaining_before=remaining_before,
                                                      outcome=BUDGET_OUTCOME_REFUSED,
                                                      reason_code=skip.reason_code, reason=skip.reason,
                                                      network_source_hash=network_source_hash, method=method))
            continue
        if max_networks is not None and admitted_networks >= max_networks:
            skip = PDepBudgetSkip(network_id=network_id,
                                  cost=cost,
                                  remaining_transition_states=remaining,
                                  reason_code=BUDGET_SKIP_MAX_NETWORKS_REACHED,
                                  reason=f'the network-count budget ({max_networks}) was already reached')
            skipped.append(skip)
            considered.append(PDepBudgetConsideration(identity=identity, network_id=network_id,
                                                      offer_indices=tuple(indices), cost=cost, rank=rank,
                                                      remaining_before=remaining_before,
                                                      outcome=BUDGET_OUTCOME_REFUSED,
                                                      reason_code=skip.reason_code, reason=skip.reason,
                                                      network_source_hash=network_source_hash, method=method))
            continue
        if remaining is not None and cost > remaining:
            skip = PDepBudgetSkip(network_id=network_id,
                                  cost=cost,
                                  remaining_transition_states=remaining,
                                  reason_code=BUDGET_SKIP_DOES_NOT_FIT_REMAINING,
                                  reason=f'its cost ({cost}) does not fit in the transition-state budget '
                                         f'remaining ({remaining})')
            skipped.append(skip)
            considered.append(PDepBudgetConsideration(identity=identity, network_id=network_id,
                                                      offer_indices=tuple(indices), cost=cost, rank=rank,
                                                      remaining_before=remaining_before,
                                                      outcome=BUDGET_OUTCOME_REFUSED,
                                                      reason_code=skip.reason_code, reason=skip.reason,
                                                      network_source_hash=network_source_hash, method=method))
            continue

        admitted_indices.extend(indices)
        admitted_networks += 1
        total_cost += cost
        considered.append(PDepBudgetConsideration(identity=identity, network_id=network_id,
                                                  offer_indices=tuple(indices), cost=cost, rank=rank,
                                                  remaining_before=remaining_before,
                                                  outcome=BUDGET_OUTCOME_ADMITTED,
                                                  network_source_hash=network_source_hash, method=method))
        if remaining is not None:
            remaining -= cost

    return PDepBudgetDecision(admitted_indices=tuple(admitted_indices),
                              skipped=tuple(skipped),
                              total_cost=total_cost,
                              considered=tuple(considered),
                              max_transition_states=max_transition_states,
                              max_networks=max_networks,
                              )


@dataclass(frozen=True)
class PDepBudgetNetworkOutcome:
    """
    What the QM budget decided for one network, admitted or refused.

    A durable record without an entry for every considered network -- admitted or not -- could not
    be told apart from "the budget never ran"; recording only refusals would make "no entry" mean
    both "admitted" and "not considered", which is exactly the kind of stale/incomplete record
    this module otherwise refuses to produce silently.

    Attributes:
        network_id (str, optional): The network's identifier, or a falsy value if unnamed.
        outcome (str): One of ``VALID_BUDGET_OUTCOMES``.
        cost (int): The transition-state cost charged to this network.
        rank (int): Position in the most-deserving-first walk (0 first), from
            :func:`t3.pdep.selector.selection_rank_key`.
        network_source_hash (str, optional): The admitted/refused network's ``network_source_hash``,
            resolved across its offers, so the record can be traced back to the network file it
            was computed from.
        method (str, optional): The admitted/refused network's master-equation method.
        reason_code (str, optional): One of ``VALID_BUDGET_SKIP_REASON_CODES`` if refused, else
            ``None``. Mirrors ``PDepBudgetSkip.reason_code``.
        reason (str, optional): Human-readable if refused, else ``None``. Mirrors
            ``PDepBudgetSkip.reason``.
        remaining_transition_states (int, optional): Transition-state budget still available
            immediately BEFORE this network was considered, or ``None`` if no transition-state
            budget was configured.
        unnamed_offer_index (int, optional): The first index that offered this network, but ONLY
            if ``network_id`` is falsy -- an unnamed network's identity, since a bare, repeated
            falsy ``network_id`` cannot tell two distinct unnamed networks apart on its own.
            ``None`` for a named network, where ``network_id`` alone is already identity.
    """
    network_id: str | None
    outcome: str
    cost: int
    rank: int
    network_source_hash: str | None = None
    method: str | None = None
    reason_code: str | None = None
    reason: str | None = None
    remaining_transition_states: int | None = None
    unnamed_offer_index: int | None = None

    def __post_init__(self):
        if self.outcome not in VALID_BUDGET_OUTCOMES:
            raise ValueError(f'PDepBudgetNetworkOutcome.outcome must be one of {VALID_BUDGET_OUTCOMES}, '
                             f'got {self.outcome!r}.')
        if self.outcome == BUDGET_OUTCOME_ADMITTED and (self.reason_code is not None or self.reason is not None):
            raise ValueError('An admitted PDepBudgetNetworkOutcome must not carry refusal detail, got '
                             f'reason_code={self.reason_code!r}, reason={self.reason!r}.')
        if self.outcome == BUDGET_OUTCOME_REFUSED:
            if self.reason_code is None or self.reason is None:
                raise ValueError('A refused PDepBudgetNetworkOutcome must carry a reason_code and reason, got '
                                 f'reason_code={self.reason_code!r}, reason={self.reason!r}.')
            if self.reason_code not in VALID_BUDGET_SKIP_REASON_CODES:
                raise ValueError(f'PDepBudgetNetworkOutcome.reason_code must be one of '
                                 f'{VALID_BUDGET_SKIP_REASON_CODES} when outcome is refused, got '
                                 f'{self.reason_code!r}.')
            if not isinstance(self.reason, str) or not self.reason:
                raise ValueError('PDepBudgetNetworkOutcome.reason must be a non-empty string when outcome is '
                                 f'refused, got {self.reason!r}.')
        if self.network_id is not None and not isinstance(self.network_id, str):
            raise ValueError(f'PDepBudgetNetworkOutcome.network_id must be a str or None, got {self.network_id!r}.')
        if self.method is not None and not isinstance(self.method, str):
            raise ValueError(f'PDepBudgetNetworkOutcome.method must be a str or None, got {self.method!r}.')
        if self.network_source_hash is not None and not isinstance(self.network_source_hash, str):
            raise ValueError('PDepBudgetNetworkOutcome.network_source_hash must be a str or None, got '
                             f'{self.network_source_hash!r}.')
        if isinstance(self.cost, bool) or not isinstance(self.cost, int) or self.cost < 0:
            raise ValueError(f'PDepBudgetNetworkOutcome.cost must be a non-negative integer, got {self.cost!r}.')
        if isinstance(self.rank, bool) or not isinstance(self.rank, int) or self.rank < 0:
            raise ValueError(f'PDepBudgetNetworkOutcome.rank must be a non-negative integer, got {self.rank!r}.')
        if self.remaining_transition_states is not None and (
                isinstance(self.remaining_transition_states, bool)
                or not isinstance(self.remaining_transition_states, int)
                or self.remaining_transition_states < 0):
            raise ValueError('PDepBudgetNetworkOutcome.remaining_transition_states must be a non-negative '
                             f'int or None, got {self.remaining_transition_states!r}.')
        if self.unnamed_offer_index is not None and (
                isinstance(self.unnamed_offer_index, bool) or not isinstance(self.unnamed_offer_index, int)
                or self.unnamed_offer_index < 0):
            raise ValueError('PDepBudgetNetworkOutcome.unnamed_offer_index must be a non-negative integer or '
                             f'None, got {self.unnamed_offer_index!r}.')
        if self.network_id and self.unnamed_offer_index is not None:
            raise ValueError('PDepBudgetNetworkOutcome.unnamed_offer_index is only for a falsy network_id (an '
                             f'unnamed network), got network_id={self.network_id!r} alongside '
                             f'unnamed_offer_index={self.unnamed_offer_index!r}.')

    def as_dict(self) -> dict:
        """Return a YAML/JSON-safe dict -- every field is already a JSON-safe primitive."""
        return {'network_id': self.network_id,
               'outcome': self.outcome,
               'cost': self.cost,
               'network_source_hash': self.network_source_hash,
               'method': self.method,
               'reason_code': self.reason_code,
               'reason': self.reason,
               'remaining_transition_states': self.remaining_transition_states,
               'rank': self.rank,
               'unnamed_offer_index': self.unnamed_offer_index,
               }


@dataclass(frozen=True)
class PDepBudgetRecord:
    """
    The durable, on-disk record of one iteration's QM budget decision.

    ``BUDGET_RECORD_SCHEMA_VERSION`` records what SHAPE this record is; ``algorithm_version``
    (``BUDGET_ALGORITHM_VERSION`` by default) records what admission/refusal LOGIC produced it --
    two different questions, so a shape-identical record from an updated algorithm can still be
    told apart from one an older ``apply_pdep_qm_budget`` would have produced given the same inputs.

    Attributes:
        iteration (int): The T3 iteration this budget was applied in.
        max_transition_states (int, optional): The transition-state budget configured, or ``None``
            if unbounded.
        max_networks (int, optional): The network-count budget configured, or ``None`` if unbounded.
        total_cost (int): The summed transition-state cost of every admitted network.
        network_outcomes (tuple): One ``PDepBudgetNetworkOutcome`` per considered network, in
            most-deserving-first order. No two outcomes may share the same network identity.
        schema_version (int): ``BUDGET_RECORD_SCHEMA_VERSION`` this record's shape matches.
        algorithm_version (int): ``BUDGET_ALGORITHM_VERSION`` this record's decision was computed
            under.
    """
    iteration: int
    max_transition_states: int | None
    max_networks: int | None
    total_cost: int
    network_outcomes: tuple
    schema_version: int = BUDGET_RECORD_SCHEMA_VERSION
    algorithm_version: int = BUDGET_ALGORITHM_VERSION

    def __post_init__(self):
        _validate_budget(self.max_transition_states, 'max_transition_states')
        _validate_budget(self.max_networks, 'max_networks')
        if isinstance(self.iteration, bool) or not isinstance(self.iteration, int) or self.iteration < 0:
            raise ValueError(f'PDepBudgetRecord.iteration must be a non-negative integer, got '
                             f'{self.iteration!r}.')
        if isinstance(self.total_cost, bool) or not isinstance(self.total_cost, int) or self.total_cost < 0:
            raise ValueError(f'PDepBudgetRecord.total_cost must be a non-negative integer, got '
                             f'{self.total_cost!r}.')
        if isinstance(self.schema_version, bool) or self.schema_version != BUDGET_RECORD_SCHEMA_VERSION:
            raise ValueError(f'PDepBudgetRecord.schema_version must be {BUDGET_RECORD_SCHEMA_VERSION}, got '
                             f'{self.schema_version!r}.')
        if isinstance(self.algorithm_version, bool) or self.algorithm_version != BUDGET_ALGORITHM_VERSION:
            raise ValueError(f'PDepBudgetRecord.algorithm_version must be {BUDGET_ALGORITHM_VERSION}, got '
                             f'{self.algorithm_version!r}.')
        for outcome in self.network_outcomes:
            if not isinstance(outcome, PDepBudgetNetworkOutcome):
                raise ValueError('PDepBudgetRecord.network_outcomes must contain only PDepBudgetNetworkOutcome '
                                 f'instances, got {outcome!r}.')
        # The identity used for the duplicate check is NOT bare network_id: two distinct unnamed
        # networks legitimately share the same falsy network_id (both None, say), and only
        # unnamed_offer_index tells them apart -- see PDepBudgetNetworkOutcome's own docstring. A
        # real, non-empty network_id repeated with unnamed_offer_index unset both times is still
        # invalid: budget semantics intentionally merge one network_id into a single outcome, so a
        # duplicate real id means identity was already corrupted before this record was built.
        identities = [(outcome.network_id, outcome.unnamed_offer_index) for outcome in self.network_outcomes]
        if len(identities) != len(set(identities)):
            raise ValueError('PDepBudgetRecord.network_outcomes must not repeat a network identity '
                             f'(network_id, unnamed_offer_index), got {identities!r}.')
        ranks = [outcome.rank for outcome in self.network_outcomes]
        if sorted(ranks) != list(range(len(ranks))):
            raise ValueError(f'PDepBudgetRecord.network_outcomes ranks must be exactly 0..{len(ranks) - 1} '
                             f'with no gaps or repeats, got {ranks!r}.')
        admitted_cost = sum(outcome.cost for outcome in self.network_outcomes
                           if outcome.outcome == BUDGET_OUTCOME_ADMITTED)
        if self.total_cost != admitted_cost:
            raise ValueError(f'PDepBudgetRecord.total_cost ({self.total_cost!r}) must equal the sum of admitted '
                             f'outcome costs ({admitted_cost!r}).')
        object.__setattr__(self, 'network_outcomes', tuple(self.network_outcomes))

    def as_dict(self) -> dict:
        """Return a YAML/JSON-safe dict -- ``network_outcomes`` recurses into each outcome's own
        ``as_dict()``, since a list of dataclasses is not YAML-safe on its own."""
        return {'budget_record_schema_version': self.schema_version,
               'budget_algorithm_version': self.algorithm_version,
               'iteration': self.iteration,
               'max_transition_states': self.max_transition_states,
               'max_networks': self.max_networks,
               'total_cost': self.total_cost,
               'network_outcomes': [outcome.as_dict() for outcome in self.network_outcomes],
               }


def build_pdep_budget_record(decision: PDepBudgetDecision, iteration: int) -> PDepBudgetRecord:
    """
    Turn one budget decision into its durable record, network outcome by network outcome.

    Consumes ONLY ``decision.considered`` (plus ``decision.max_transition_states``,
    ``decision.max_networks``, and ``decision.total_cost``) -- everything already computed once,
    in ``apply_pdep_qm_budget``'s single authoritative walk. There is no ``selections`` parameter
    here and no grouping/ranking/remaining-budget bookkeeping to keep in lockstep with anything: a
    second, independently-maintained copy of that walk was exactly the drift risk -- and, worse,
    the risk of being handed a ``decision``/``selections`` pair that never actually came from the
    same call -- that an adversarial review of an earlier version of this function found.

    Args:
        decision (PDepBudgetDecision): The decision ``apply_pdep_qm_budget`` returned.
        iteration (int): The T3 iteration this budget was applied in.

    Returns:
        PDepBudgetRecord: One outcome per considered network, most-deserving-first.

    Raises:
        ValueError: If ``decision.considered`` does not account for every admission and refusal the
            decision reports.
    """
    # ``considered`` is defaulted, so a decision can be constructed without it -- and a record built
    # from such a decision would say, with the full authority of a durable artifact, that nothing was
    # considered, while the same decision reports networks admitted and refused. Silence is the one
    # thing this record exists to prevent, so an incoherent decision is refused here, at the point
    # where authority is claimed, rather than rendered.
    admitted_offer_indices = {index for consideration in decision.considered
                              if consideration.outcome == BUDGET_OUTCOME_ADMITTED
                              for index in consideration.offer_indices}
    refused_count = sum(1 for consideration in decision.considered
                        if consideration.outcome == BUDGET_OUTCOME_REFUSED)
    if admitted_offer_indices != set(decision.admitted_indices) or refused_count != len(decision.skipped):
        raise ValueError(f'Refusing to build a budget record from a decision whose considered networks do '
                         f'not account for what it reports. The decision admits offer indices '
                         f'{sorted(set(decision.admitted_indices))} and refuses {len(decision.skipped)} '
                         f'network(s), but its considered networks admit {sorted(admitted_offer_indices)} '
                         f'and refuse {refused_count}.')
    network_outcomes = tuple(
        PDepBudgetNetworkOutcome(network_id=consideration.network_id,
                                 outcome=consideration.outcome,
                                 cost=consideration.cost,
                                 rank=consideration.rank,
                                 network_source_hash=consideration.network_source_hash,
                                 method=consideration.method,
                                 reason_code=consideration.reason_code,
                                 reason=consideration.reason,
                                 remaining_transition_states=consideration.remaining_before,
                                 unnamed_offer_index=(consideration.offer_indices[0]
                                                      if not consideration.network_id else None),
                                 )
        for consideration in decision.considered
        )
    return PDepBudgetRecord(iteration=iteration,
                            max_transition_states=decision.max_transition_states,
                            max_networks=decision.max_networks,
                            total_cost=decision.total_cost,
                            network_outcomes=network_outcomes,
                            )
