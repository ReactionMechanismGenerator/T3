#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_budget
"""

from dataclasses import replace

import pytest

from t3.pdep.budget import (BUDGET_ALGORITHM_VERSION,
                            BUDGET_OUTCOME_ADMITTED,
                            BUDGET_OUTCOME_REFUSED,
                            BUDGET_RECORD_SCHEMA_VERSION,
                            BUDGET_SKIP_DOES_NOT_FIT_REMAINING,
                            BUDGET_SKIP_EXCEEDS_BUDGET,
                            BUDGET_SKIP_MAX_NETWORKS_REACHED,
                            BUDGET_SKIP_NOT_EVALUATED,
                            PDepBudgetDecision,
                            PDepBudgetNetworkOutcome,
                            PDepBudgetRecord,
                            PDepBudgetSkip,
                            VALID_BUDGET_SKIP_REASON_CODES,
                            apply_pdep_qm_budget,
                            build_pdep_budget_record,
                            projected_ts_cost,
                            )
from t3.pdep.selector import (EVALUATION_STATUS_EVALUATED,
                              EVALUATION_STATUS_NOT_EVALUATED,
                              PDepNetworkSelection,
                              SensitiveTransitionState,
                              )


def _make_ts(ts_label, delta_ln_k=1.0, uncertain=True):
    """Build a minimal SensitiveTransitionState with a controllable ts_label and delta_ln_k."""
    return SensitiveTransitionState(ts_label=ts_label, coefficient=0.1, condition=(300.0, 'K', 1.0, 'bar'),
                                    path_reaction_label=f'reaction_{ts_label}', path_reaction_str='A + B <=> C',
                                    kinetics_comment='Estimated using template [x]', uncertain=uncertain,
                                    delta_ln_k=delta_ln_k)


def _make_selection(network_id, n_ts=1, delta_ln_k=1.0, qualified=True,
                    evaluation_status=EVALUATION_STATUS_EVALUATED):
    """Build a qualified PDepNetworkSelection through its real constructor, with ``n_ts`` distinct,
    uncertain transition state labels (each providing the given ``delta_ln_k``) -- this both
    controls the network's transition-state cost (``n_ts``) and, via ``delta_ln_k``, its rank."""
    entries = [_make_ts(f'{network_id}_TS{i}', delta_ln_k=delta_ln_k) for i in range(n_ts)]
    return PDepNetworkSelection(
        network_id=network_id,
        qualified=qualified,
        uncertain_path_reactions=list(entries),
        selected_ts=list(entries),
        evaluation_status=evaluation_status,
    )


# --- 1. No budgets: every index admitted, nothing skipped, total_cost sums all costs -----------

def test_no_budgets_admits_everything():
    """Test that with both budgets None, every index is admitted, skipped is empty, and total_cost
    is the sum of all network costs -- if this regresses, opting out of budgeting (the documented
    ``None``-means-no-limit default) would silently start dropping networks."""
    selections = [_make_selection('netA', n_ts=2, delta_ln_k=1.0),
                 _make_selection('netB', n_ts=3, delta_ln_k=2.0),
                 _make_selection('netC', n_ts=1, delta_ln_k=3.0)]
    decision = apply_pdep_qm_budget(selections)
    assert set(decision.admitted_indices) == {0, 1, 2}
    assert decision.skipped == ()
    assert decision.total_cost == 2 + 3 + 1


# --- 2. No budgets: admitted order is RANKED order, not input order -----------------------------

def test_no_budgets_admitted_order_is_ranked_not_input_order():
    """Test that with no budgets, the admitted order follows the ranked (most-deserving-first)
    order, not the input order -- if the module stopped sorting before walking, an unbudgeted run
    would still silently reorder which networks are queued first downstream."""
    # netA is offered first but is LESS deserving (lower delta_ln_k); netB is offered second but is
    # MORE deserving. Ranked order must be [netB, netA], the reverse of input order.
    selections = [_make_selection('netA', n_ts=1, delta_ln_k=1.0),
                 _make_selection('netB', n_ts=1, delta_ln_k=5.0)]
    decision = apply_pdep_qm_budget(selections)
    assert decision.admitted_indices == (1, 0)


# --- 3. max_transition_states admits in ranked order while it fits, cost never exceeds budget ---

def test_max_transition_states_admits_in_ranked_order_within_budget():
    """Test that max_transition_states admits ranked networks while they fit and that total
    admitted cost never exceeds the budget -- if the walk stopped honoring rank or overspent, a
    less-deserving network could be queued ahead of (or alongside costs beyond) a more-deserving one."""
    selections = [_make_selection('netA', n_ts=2, delta_ln_k=1.0),
                 _make_selection('netB', n_ts=2, delta_ln_k=5.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=2)
    # netB (rank first, cost 2) fits exactly; netA (cost 2) would push total to 4 and is refused.
    assert decision.admitted_indices == (1,)
    assert decision.total_cost == 2
    assert decision.total_cost <= 2


# --- 4. A network that does not fit is skipped WHOLE, never partially taken --------------------

def test_network_that_does_not_fit_is_skipped_whole():
    """Test that a network whose cost exceeds the remaining transition-state budget is skipped
    WHOLE -- none of its indices are admitted -- rather than partially taken. A half-admitted
    network would leave some of its transition states queued for QM without the master-equation
    update that consumes them all together, per the module's 'never slices a network' guarantee."""
    selections = [_make_selection('netA', n_ts=5, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=3)
    assert decision.admitted_indices == ()
    assert len(decision.skipped) == 1
    assert decision.skipped[0].network_id == 'netA'


# --- 5. The walk continues past a network that does not fit ------------------------------------

def test_walk_continues_past_a_skipped_network():
    """Test that skipping a network that does not fit does not stop the walk: a later, cheaper,
    less-deserving network must still be admitted if it fits the remaining budget -- stopping at
    the first refusal would leave budget unspent for no benefit, which the module docstring
    explicitly rejects as the alternative."""
    # netA (rank first, cost 5) does not fit budget 3; netB (rank second, cost 2) does.
    selections = [_make_selection('netA', n_ts=5, delta_ln_k=10.0),
                 _make_selection('netB', n_ts=2, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=3)
    assert 1 in decision.admitted_indices
    assert 0 not in decision.admitted_indices
    assert decision.total_cost == 2


# --- 6. max_networks caps the number of distinct networks, independent of TS cost ---------------

def test_max_networks_caps_distinct_network_count():
    """Test that max_networks caps the number of DISTINCT networks admitted regardless of
    transition-state cost -- a network-count cap must bind even when the transition-state budget
    (here, None/unlimited) would happily admit more."""
    selections = [_make_selection('netA', n_ts=1, delta_ln_k=3.0),
                 _make_selection('netB', n_ts=1, delta_ln_k=2.0),
                 _make_selection('netC', n_ts=1, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections, max_networks=2)
    admitted_network_ids = {selections[index].network_id for index in decision.admitted_indices}
    assert len(admitted_network_ids) == 2
    assert admitted_network_ids == {'netA', 'netB'}


# --- 7. Both budgets together: whichever binds first refuses -----------------------------------

def test_both_budgets_transition_state_budget_binds_first():
    """Test that when both budgets are set and the transition-state budget is the one that would
    refuse a network first (max_networks alone would still allow it), that network is skipped and
    the reason names the transition-state bound."""
    selections = [_make_selection('netA', n_ts=5, delta_ln_k=2.0),
                 _make_selection('netB', n_ts=1, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=3, max_networks=5)
    assert decision.admitted_indices == (1,)
    assert decision.skipped[0].network_id == 'netA'
    assert 'state' in decision.skipped[0].reason


def test_both_budgets_max_networks_binds_first():
    """Test that when both budgets are set and max_networks is the one that would refuse a network
    first (the transition-state budget alone would still have room), that network is skipped and
    the reason names the network-count bound."""
    selections = [_make_selection('netA', n_ts=1, delta_ln_k=3.0),
                 _make_selection('netB', n_ts=1, delta_ln_k=2.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=100, max_networks=1)
    assert decision.admitted_indices == (0,)
    assert decision.skipped[0].network_id == 'netB'
    assert 'network' in decision.skipped[0].reason


# --- 8. A network offered twice is charged once, and both indices are admitted ------------------

def test_repeated_network_id_charged_once_both_indices_admitted():
    """Test that a network offered twice (two entries sharing one network_id, the normal case of
    several sensitive reactions belonging to one network) is charged ONCE against the budget, and,
    if admitted, BOTH of its indices appear in admitted_indices -- double-charging would starve the
    budget for networks that legitimately appear more than once."""
    shared = _make_selection('netA', n_ts=2, delta_ln_k=1.0)
    other = _make_selection('netB', n_ts=2, delta_ln_k=1.0)
    selections = [shared, other, shared]
    decision = apply_pdep_qm_budget(selections, max_transition_states=2)
    # Only one network's worth of cost (2) can be admitted; netA's two indices (0 and 2) must both
    # be admitted together, charged once, or neither -- assert whichever network was admitted has
    # exactly the number of indices its offers contributed.
    assert decision.total_cost == 2
    admitted_network_ids = {selections[index].network_id for index in decision.admitted_indices}
    assert len(admitted_network_ids) == 1
    if 'netA' in admitted_network_ids:
        assert set(decision.admitted_indices) == {0, 2}


def test_repeated_network_id_both_indices_admitted_when_it_fits():
    """Test directly (unambiguous case) that when the repeated network fits the budget alongside
    another network, both of its indices are admitted and it is charged only once."""
    shared = _make_selection('netA', n_ts=1, delta_ln_k=5.0)
    other = _make_selection('netB', n_ts=1, delta_ln_k=1.0)
    selections = [shared, other, shared]
    decision = apply_pdep_qm_budget(selections, max_transition_states=2)
    assert set(decision.admitted_indices) == {0, 1, 2}
    assert decision.total_cost == 2


# --- 9. Every refusal is recorded in skipped with network_id, cost, and a reason ----------------

def test_skipped_entries_record_network_id_cost_and_reason():
    """Test that every refusal is recorded in skipped with the refused network's network_id, the
    cost that was charged against it, and a non-empty reason -- a refusal that dropped any of these
    would be a silently-shrunk decision the module docstring says must never happen."""
    selections = [_make_selection('netA', n_ts=5, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=3)
    assert len(decision.skipped) == 1
    skip = decision.skipped[0]
    assert isinstance(skip, PDepBudgetSkip)
    assert skip.network_id == 'netA'
    assert skip.cost == 5
    assert skip.reason and isinstance(skip.reason, str)


def test_nothing_silently_dropped_admitted_plus_skipped_covers_all_networks():
    """Test that every distinct network offered is accounted for as either admitted or skipped --
    nothing vanishes without a trace."""
    selections = [_make_selection('netA', n_ts=5, delta_ln_k=3.0),
                 _make_selection('netB', n_ts=1, delta_ln_k=2.0),
                 _make_selection('netC', n_ts=1, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=2)
    admitted_network_ids = {selections[index].network_id for index in decision.admitted_indices}
    skipped_network_ids = {skip.network_id for skip in decision.skipped}
    assert admitted_network_ids | skipped_network_ids == {'netA', 'netB', 'netC'}
    assert admitted_network_ids.isdisjoint(skipped_network_ids)


# --- 10. apply_pdep_qm_budget rejects a budget that is 0, negative, a float, or a bool ----------

@pytest.mark.parametrize('bad_value', [0, -1, 1.5, True, False])
def test_apply_pdep_qm_budget_rejects_bad_max_transition_states(bad_value):
    """Test that a zero, negative, float, or bool max_transition_states raises ValueError -- a bool
    is an int subclass, so without an explicit isinstance(value, bool) rejection,
    max_transition_states=True would silently mean 'admit at most one transition state'."""
    selections = [_make_selection('netA', n_ts=1)]
    with pytest.raises(ValueError):
        apply_pdep_qm_budget(selections, max_transition_states=bad_value)


@pytest.mark.parametrize('bad_value', [0, -1, 1.5, True, False])
def test_apply_pdep_qm_budget_rejects_bad_max_networks(bad_value):
    """Test that a zero, negative, float, or bool max_networks raises ValueError -- a bool is an
    int subclass, so without an explicit isinstance(value, bool) rejection, max_networks=True
    would silently mean 'admit at most one network'."""
    selections = [_make_selection('netA', n_ts=1)]
    with pytest.raises(ValueError):
        apply_pdep_qm_budget(selections, max_networks=bad_value)


def test_apply_pdep_qm_budget_accepts_none_and_positive_int():
    """Test that None and a positive int are accepted for both budgets (the valid cases bracketing
    the rejections above), so the rejection tests are not vacuously passing on everything."""
    selections = [_make_selection('netA', n_ts=1)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=None, max_networks=None)
    assert decision.admitted_indices == (0,)
    decision = apply_pdep_qm_budget(selections, max_transition_states=1, max_networks=1)
    assert decision.admitted_indices == (0,)


# --- 11. projected_ts_cost equals the number of distinct uncertain TS labels --------------------

def test_projected_ts_cost_equals_distinct_uncertain_ts_label_count():
    """Test that projected_ts_cost equals len(selection.uncertain_ts_labels()) -- if it instead
    counted evidence entries, a network sensitive to the same transition state via several path
    reactions or conditions would be charged as if it required more independent QM jobs than it
    actually does."""
    selection = _make_selection('netA', n_ts=3, delta_ln_k=1.0)
    assert projected_ts_cost(selection) == 3
    assert projected_ts_cost(selection) == len(selection.uncertain_ts_labels())


def test_projected_ts_cost_dedupes_shared_ts_label_across_entries():
    """Test that a selection whose several uncertain evidence entries share ONE ts_label costs 1,
    not the number of entries -- a transition state shared by two path reactions (the module
    docstring's own example) becomes one ARC job, not two."""
    shared_label_entries = [_make_ts('TS1', uncertain=True), _make_ts('TS1', uncertain=True)]
    selection = PDepNetworkSelection(network_id='netA', qualified=True,
                                     uncertain_path_reactions=shared_label_entries)
    assert projected_ts_cost(selection) == 1
    assert len(selection.uncertain_path_reactions) == 2


# --- 12. Deterministic tie-break: identical rank evidence orders by network_id ------------------

def test_ties_broken_by_network_id():
    """Test that two networks with identical rank evidence (same qualified tier, same
    max delta_ln_k) come out in network_id order -- without a deterministic tie-break, admission
    order (and, under a tight budget, which of two equally-deserving networks gets admitted) would
    depend on incidental input order or dict/set iteration, not a documented rule."""
    selections = [_make_selection('netZ', n_ts=1, delta_ln_k=1.0),
                 _make_selection('netA', n_ts=1, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections)
    # netA and netZ tie on tier and delta_ln_k; network_id order puts netA first regardless of
    # input order (netZ was offered first here).
    assert decision.admitted_indices == (1, 0)


def test_ties_broken_by_network_id_under_a_binding_budget():
    """Test the tie-break under a budget tight enough that only one of two tied networks can be
    admitted: the one admitted must be the one with the earlier network_id, not the one offered
    first."""
    selections = [_make_selection('netZ', n_ts=1, delta_ln_k=1.0),
                 _make_selection('netA', n_ts=1, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections, max_networks=1)
    admitted_network_ids = {selections[index].network_id for index in decision.admitted_indices}
    assert admitted_network_ids == {'netA'}


# --- 13. Repeat offers naming different ts_labels are charged for the UNION, not the first -----

def test_repeated_network_id_charged_for_union_of_ts_labels_not_first_offer():
    """Test that two entries sharing one network_id but naming DIFFERENT ts_labels (offer A names
    {TS1, TS2}, offer B names {TS2, TS3}) are charged for the UNION of those labels (3), not for
    the first offer alone (2) -- if only the first offer's labels were counted, a budget that
    cannot afford the union would wrongly admit the network as if it only needed 2 transition
    states, under-spending the very budget the module promises never to under-estimate."""
    offer_a = PDepNetworkSelection(
        network_id='netA', qualified=True,
        uncertain_path_reactions=[_make_ts('TS1'), _make_ts('TS2')],
    )
    offer_b = PDepNetworkSelection(
        network_id='netA', qualified=True,
        uncertain_path_reactions=[_make_ts('TS2'), _make_ts('TS3')],
    )
    selections = [offer_a, offer_b]
    # Budget fits 2 transition states but not the true union of 3.
    decision = apply_pdep_qm_budget(selections, max_transition_states=2)
    assert decision.admitted_indices == ()
    assert len(decision.skipped) == 1
    assert decision.skipped[0].network_id == 'netA'
    assert decision.skipped[0].cost == 3


# --- 14. Repeat offers are ranked by the BEST of them, not the first ----------------------------

def test_repeated_network_id_ranked_by_best_offer_not_first():
    """Test that a network offered twice is ranked by the BEST (most deserving) of its offers, not
    the first one seen -- a weak first offer (small delta_ln_k) followed by a strong second offer
    (large delta_ln_k) must still let the network out-rank a single-offer competitor whose strength
    sits strictly between the two; ranking by the first offer alone would wrongly bury this network
    behind the competitor."""
    weak_first = PDepNetworkSelection(
        network_id='netA', qualified=True,
        uncertain_path_reactions=[_make_ts('netA_TS1', delta_ln_k=1.0)],
    )
    strong_second = PDepNetworkSelection(
        network_id='netA', qualified=True,
        uncertain_path_reactions=[_make_ts('netA_TS2', delta_ln_k=10.0)],
    )
    competitor = _make_selection('netB', n_ts=1, delta_ln_k=5.0)
    selections = [weak_first, competitor, strong_second]
    decision = apply_pdep_qm_budget(selections, max_networks=1)
    admitted_network_ids = {selections[index].network_id for index in decision.admitted_indices}
    assert admitted_network_ids == {'netA'}


# --- 15. A selection with no network_id is its own network, distinct from other unnamed ones ---

@pytest.mark.parametrize('empty_id', [None, ''])
def test_selection_without_network_id_is_its_own_network(empty_id):
    """Test that two selections both carrying no network_id (None, and separately '') are treated
    as TWO distinct networks, not collapsed into one -- with max_networks=1, exactly one is
    admitted and the other is skipped. Collapsing unnamed decisions by their (shared, empty)
    network_id would charge several genuinely distinct networks as if they were a single one,
    silently starving the budget for the rest."""
    first = PDepNetworkSelection(
        network_id=empty_id, qualified=True,
        uncertain_path_reactions=[_make_ts('unnamed_TS1', delta_ln_k=5.0)],
    )
    second = PDepNetworkSelection(
        network_id=empty_id, qualified=True,
        uncertain_path_reactions=[_make_ts('unnamed_TS2', delta_ln_k=1.0)],
    )
    selections = [first, second]
    decision = apply_pdep_qm_budget(selections, max_networks=1)
    assert len(decision.admitted_indices) == 1
    assert len(decision.skipped) == 1


# --- 16. A network no offer could evaluate is refused before it is ranked ----------------------

def test_network_that_could_not_be_evaluated_is_refused_with_a_measurement_reason():
    """Test that a network none of whose offers could be evaluated (evaluation_status ==
    EVALUATION_STATUS_NOT_EVALUATED) is refused outright, with a reason that says it could not be
    evaluated -- per the module docstring, such a network carries no signal and must never be
    ranked as if it had one."""
    selections = [_make_selection('netA', n_ts=1, delta_ln_k=1.0,
                                  qualified=False,
                                  evaluation_status=EVALUATION_STATUS_NOT_EVALUATED)]
    decision = apply_pdep_qm_budget(selections)
    assert decision.admitted_indices == ()
    assert len(decision.skipped) == 1
    assert decision.skipped[0].network_id == 'netA'
    assert 'could not be evaluated' in decision.skipped[0].reason


def test_network_with_one_not_evaluated_and_one_evaluated_offer_is_not_refused_for_that_reason():
    """Test the converse of the above: a network with ONE not-evaluated offer and ONE evaluated
    offer has a measurement (from the evaluated offer) and must NOT be refused with the
    'could not be evaluated' reason -- only a network for which every single offer failed to
    evaluate carries no signal at all."""
    not_evaluated_offer = PDepNetworkSelection(
        network_id='netA', qualified=False,
        uncertain_path_reactions=[_make_ts('netA_TS1', delta_ln_k=1.0)],
        evaluation_status=EVALUATION_STATUS_NOT_EVALUATED,
    )
    evaluated_offer = PDepNetworkSelection(
        network_id='netA', qualified=True,
        uncertain_path_reactions=[_make_ts('netA_TS2', delta_ln_k=1.0)],
        evaluation_status=EVALUATION_STATUS_EVALUATED,
    )
    selections = [not_evaluated_offer, evaluated_offer]
    decision = apply_pdep_qm_budget(selections)
    assert decision.admitted_indices == (0, 1)
    assert decision.skipped == ()


# --- 17. A network whose cost exceeds the WHOLE transition-state budget is refused permanently --

def test_network_exceeding_entire_budget_gets_a_distinct_permanent_refusal():
    """Test that a network whose cost exceeds the WHOLE transition-state budget (not just what
    currently remains) is skipped with a reason distinct from, and never confusable with, the
    ordinary 'does not fit right now' refusal -- and that this reason says the network can never be
    refined until the limit itself is raised, since no future iteration's remaining budget can ever
    exceed the whole budget."""
    exceeds_whole_budget = _make_selection('netA', n_ts=5, delta_ln_k=1.0)
    fits_but_loses_this_round = _make_selection('netB', n_ts=3, delta_ln_k=10.0)
    selections = [exceeds_whole_budget, fits_but_loses_this_round]
    decision = apply_pdep_qm_budget(selections, max_transition_states=3)
    skips_by_id = {skip.network_id: skip for skip in decision.skipped}
    assert 'netA' in skips_by_id
    permanent_reason = skips_by_id['netA'].reason
    assert 'never' in permanent_reason or 'raise' in permanent_reason or 'limit' in permanent_reason

    # Contrast with an ordinary same-iteration "only N remain" refusal: a network that fits within
    # the whole budget but not within what remains after a more-deserving network was admitted.
    admits_one = _make_selection('netX', n_ts=2, delta_ln_k=10.0)
    loses_to_remaining_budget = _make_selection('netY', n_ts=2, delta_ln_k=1.0)
    ordinary_decision = apply_pdep_qm_budget([admits_one, loses_to_remaining_budget],
                                             max_transition_states=2)
    ordinary_reason = ordinary_decision.skipped[0].reason
    assert ordinary_reason != permanent_reason


# --- 18. Each of the four refusal reasons carries the matching machine-stable reason_code --------

def test_not_evaluated_skip_carries_not_evaluated_reason_code():
    """Test that a network refused because no offer could be evaluated carries
    reason_code == BUDGET_SKIP_NOT_EVALUATED -- a machine reader must be able to tell this refusal
    apart from the other three without parsing prose."""
    selections = [_make_selection('netA', n_ts=1, delta_ln_k=1.0, qualified=False,
                                  evaluation_status=EVALUATION_STATUS_NOT_EVALUATED)]
    decision = apply_pdep_qm_budget(selections)
    assert decision.skipped[0].reason_code == BUDGET_SKIP_NOT_EVALUATED


def test_exceeds_whole_budget_skip_carries_exceeds_budget_reason_code():
    """Test that a network whose cost exceeds the entire transition-state budget carries
    reason_code == BUDGET_SKIP_EXCEEDS_BUDGET."""
    selections = [_make_selection('netA', n_ts=5, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=3)
    assert decision.skipped[0].reason_code == BUDGET_SKIP_EXCEEDS_BUDGET


def test_max_networks_reached_skip_carries_max_networks_reached_reason_code():
    """Test that a network refused because the per-iteration network-count limit was already
    reached carries reason_code == BUDGET_SKIP_MAX_NETWORKS_REACHED."""
    selections = [_make_selection('netA', n_ts=1, delta_ln_k=3.0),
                 _make_selection('netB', n_ts=1, delta_ln_k=2.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=100, max_networks=1)
    assert decision.skipped[0].reason_code == BUDGET_SKIP_MAX_NETWORKS_REACHED


def test_does_not_fit_remaining_skip_carries_does_not_fit_remaining_reason_code():
    """Test that a network that fits the whole transition-state budget but not what currently
    remains carries reason_code == BUDGET_SKIP_DOES_NOT_FIT_REMAINING, distinct from the
    exceeds-entire-budget reason code even though both are transition-state refusals."""
    selections = [_make_selection('netA', n_ts=2, delta_ln_k=10.0),
                 _make_selection('netB', n_ts=3, delta_ln_k=5.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=3)
    skips_by_id = {skip.network_id: skip for skip in decision.skipped}
    assert skips_by_id['netB'].reason_code == BUDGET_SKIP_DOES_NOT_FIT_REMAINING


def test_every_reason_code_is_one_of_the_valid_tuple():
    """Test that VALID_BUDGET_SKIP_REASON_CODES contains exactly the four codes above -- catches a
    typo'd or forgotten code silently falling outside the validated set."""
    assert set(VALID_BUDGET_SKIP_REASON_CODES) == {BUDGET_SKIP_NOT_EVALUATED, BUDGET_SKIP_EXCEEDS_BUDGET,
                                                    BUDGET_SKIP_MAX_NETWORKS_REACHED,
                                                    BUDGET_SKIP_DOES_NOT_FIT_REMAINING}


def test_pdep_budget_skip_rejects_unknown_reason_code():
    """Test that constructing a PDepBudgetSkip with a reason_code outside the valid tuple raises
    ValueError -- an unrecognized code is worse than no code at all, since a downstream reader would
    silently treat it as a fifth, undocumented reason."""
    with pytest.raises(ValueError):
        PDepBudgetSkip(network_id='netA', cost=1, remaining_transition_states=None,
                      reason_code='some_made_up_reason', reason='irrelevant')


def test_pdep_budget_skip_still_carries_prose_reason():
    """Test that PDepBudgetSkip.reason (the prose t3/main.py logs) still works exactly as before --
    reason_code is additive, not a replacement, so the existing warning log must not degrade."""
    selections = [_make_selection('netA', n_ts=5, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=3)
    assert isinstance(decision.skipped[0].reason, str) and decision.skipped[0].reason


# --- 19. PDepBudgetNetworkOutcome enforces fail-closed coherence between outcome and reason -------

def test_network_outcome_admitted_forbids_reason_and_reason_code():
    """Test that an 'admitted' PDepBudgetNetworkOutcome cannot carry a reason or reason_code -- an
    admitted network was never refused, so carrying refusal detail would be a false record."""
    with pytest.raises(ValueError):
        PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_ADMITTED, cost=1,
                                 reason_code=BUDGET_SKIP_NOT_EVALUATED, reason='should not be here',
                                 rank=0)


def test_network_outcome_refused_requires_reason_and_reason_code():
    """Test that a 'refused' PDepBudgetNetworkOutcome must carry both a reason and a reason_code --
    a refusal with no stated reason is unusable to a reader trying to understand what happened."""
    with pytest.raises(ValueError):
        PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_REFUSED, cost=1,
                                 reason_code=None, reason=None, rank=0)


def test_network_outcome_rejects_unknown_outcome_value():
    """Test that an outcome value outside {'admitted', 'refused'} raises ValueError."""
    with pytest.raises(ValueError):
        PDepBudgetNetworkOutcome(network_id='netA', outcome='maybe', cost=1, rank=0)


def test_network_outcome_as_dict_is_json_safe():
    """Test that PDepBudgetNetworkOutcome.as_dict() returns a plain dict of JSON/YAML-safe types."""
    outcome = PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_REFUSED, cost=3,
                                       network_source_hash='sha256:abc', method='MSC',
                                       reason_code=BUDGET_SKIP_NOT_EVALUATED, reason='could not be evaluated',
                                       remaining_transition_states=5, rank=0)
    as_dict = outcome.as_dict()
    assert as_dict == {
        'network_id': 'netA',
        'outcome': BUDGET_OUTCOME_REFUSED,
        'cost': 3,
        'network_source_hash': 'sha256:abc',
        'method': 'MSC',
        'reason_code': BUDGET_SKIP_NOT_EVALUATED,
        'reason': 'could not be evaluated',
        'remaining_transition_states': 5,
        'rank': 0,
        'unnamed_offer_index': None,
    }


# --- 20. build_pdep_budget_record resolves indices to identities, admitted AND refused -----------

def test_build_pdep_budget_record_covers_admitted_and_refused_networks():
    """Test that build_pdep_budget_record produces exactly one outcome per distinct network, with
    the correct outcome/cost for both the admitted and the refused network -- this is the whole
    point of the record: without the admitted side, a reader could not tell 'nothing refused' from
    'nothing recorded'."""
    selections = [_make_selection('netA', n_ts=5, delta_ln_k=10.0),
                 _make_selection('netB', n_ts=2, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=3)
    record = build_pdep_budget_record(decision, iteration=2)
    outcomes_by_id = {outcome.network_id: outcome for outcome in record.network_outcomes}
    assert set(outcomes_by_id) == {'netA', 'netB'}
    assert outcomes_by_id['netA'].outcome == BUDGET_OUTCOME_REFUSED
    assert outcomes_by_id['netA'].cost == 5
    assert outcomes_by_id['netA'].reason_code == BUDGET_SKIP_EXCEEDS_BUDGET
    assert outcomes_by_id['netB'].outcome == BUDGET_OUTCOME_ADMITTED
    assert outcomes_by_id['netB'].cost == 2
    assert outcomes_by_id['netB'].reason_code is None
    assert outcomes_by_id['netB'].reason is None
    assert record.iteration == 2
    assert record.max_transition_states == 3
    assert record.max_networks is None
    assert record.total_cost == decision.total_cost


def test_build_pdep_budget_record_collapses_repeated_network_id_to_one_outcome():
    """Test that a network offered twice (sharing one network_id) produces exactly ONE outcome
    entry in the record, charged for the union cost apply_pdep_qm_budget actually charged -- per
    the module's own union-charging semantics (budget.py's ``cost = len({ts_label for selection in
    offers for ts_label in selection.uncertain_ts_labels()})``), a duplicate entry per offer would
    misrepresent one network as two."""
    shared = _make_selection('netA', n_ts=2, delta_ln_k=1.0)
    other = _make_selection('netB', n_ts=2, delta_ln_k=1.0)
    selections = [shared, other, shared]
    decision = apply_pdep_qm_budget(selections, max_transition_states=2)
    record = build_pdep_budget_record(decision, iteration=1)
    assert len(record.network_outcomes) == 2
    netA_outcomes = [outcome for outcome in record.network_outcomes if outcome.network_id == 'netA']
    assert len(netA_outcomes) == 1
    assert netA_outcomes[0].cost == 2


def test_build_pdep_budget_record_pulls_network_source_hash_and_method():
    """Test that build_pdep_budget_record carries network_source_hash and method through from
    decision.considered -- these do not exist anywhere on PDepBudgetSkip, only on the
    PDepBudgetConsideration apply_pdep_qm_budget's own walk already resolved them onto."""
    selection = PDepNetworkSelection(network_id='netA', qualified=True,
                                     network_source_hash='sha256:deadbeef', method='MSC',
                                     uncertain_path_reactions=[_make_ts('netA_TS1')],
                                     selected_ts=[_make_ts('netA_TS1')])
    decision = apply_pdep_qm_budget([selection])
    record = build_pdep_budget_record(decision, iteration=0)
    assert record.network_outcomes[0].network_source_hash == 'sha256:deadbeef'
    assert record.network_outcomes[0].method == 'MSC'


def test_build_pdep_budget_record_as_dict_is_json_safe_and_deep_copied():
    """Test that PDepBudgetRecord.as_dict() renders a plain dict with a list of plain dicts for
    network_outcomes, and that mutating the returned dict does not affect the record."""
    selections = [_make_selection('netA', n_ts=1, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections)
    record = build_pdep_budget_record(decision, iteration=3)
    as_dict = record.as_dict()
    assert as_dict['budget_record_schema_version'] == BUDGET_RECORD_SCHEMA_VERSION
    assert as_dict['iteration'] == 3
    assert isinstance(as_dict['network_outcomes'], list)
    assert isinstance(as_dict['network_outcomes'][0], dict)
    as_dict['network_outcomes'].append('mutated')
    assert len(record.network_outcomes) == 1


def test_build_pdep_budget_record_outcomes_are_in_ranked_order():
    """Test that record.network_outcomes is in the same most-deserving-first rank order
    apply_pdep_qm_budget itself used -- both the admitted network (rank 0) and the refused one
    (rank 1) must appear in that order, and each outcome's rank field must say so."""
    selections = [_make_selection('netB', n_ts=1, delta_ln_k=10.0),
                 _make_selection('netA', n_ts=5, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections, max_transition_states=3)
    record = build_pdep_budget_record(decision, iteration=0)
    assert [outcome.network_id for outcome in record.network_outcomes] == ['netB', 'netA']
    assert [outcome.rank for outcome in record.network_outcomes] == [0, 1]


# --- 21. Two unnamed networks sharing the same falsy network_id survive the WHOLE builder --------

def test_two_unnamed_networks_with_same_falsy_network_id_survive_the_whole_builder():
    """Test that two DISTINCT unnamed networks that both carry network_id=None do not collide
    anywhere in the pipeline, from apply_pdep_qm_budget through build_pdep_budget_record -- if the
    duplicate check ever went back to bare network_id, this would crash on PDepBudgetRecord
    construction (two entries with network_id=None), rather than telling the two apart by
    unnamed_offer_index the way the whole builder is documented to."""
    selections = [_make_selection(None, n_ts=1, delta_ln_k=5.0),
                 _make_selection(None, n_ts=2, delta_ln_k=1.0)]
    decision = apply_pdep_qm_budget(selections)
    record = build_pdep_budget_record(decision, iteration=0)
    assert len(record.network_outcomes) == 2
    unnamed_offer_indices = {outcome.unnamed_offer_index for outcome in record.network_outcomes}
    assert unnamed_offer_indices == {0, 1}
    assert all(outcome.network_id is None for outcome in record.network_outcomes)


# --- 22. Repeated network with disagreeing network_source_hash raises, never silently picks one --

def test_apply_pdep_qm_budget_raises_on_conflicting_network_source_hash():
    """Test that two offers of the same network_id carrying two DIFFERENT, both-real
    network_source_hash values raise -- mirroring PDepNetworkSelection.combine(), for which a real
    hash disagreement means the offers do not actually describe the same network revision at all,
    so silently picking one (e.g. offers[0]'s) would misattribute the record to the wrong file."""
    offers = [_make_selection('netA', n_ts=1, delta_ln_k=1.0),
             _make_selection('netA', n_ts=1, delta_ln_k=2.0)]
    offers[0] = replace(offers[0], network_source_hash='sha256:aaaa')
    offers[1] = replace(offers[1], network_source_hash='sha256:bbbb')
    with pytest.raises(ValueError, match='different'):
        apply_pdep_qm_budget(offers)


# --- 23. Repeated network with hash present on only some offers resolves to None, never adopted --

def test_apply_pdep_qm_budget_partially_missing_hash_resolves_to_none():
    """Test that when only SOME offers of a repeated network record a network_source_hash (the
    others leaving it None), the considered network_source_hash resolves to None rather than
    silently adopting the one real hash that happened to be recorded -- mirroring combine()'s rule
    that a partially-missing hash is not evidence, just an absent one."""
    offers = [_make_selection('netA', n_ts=1, delta_ln_k=1.0),
             _make_selection('netA', n_ts=1, delta_ln_k=2.0)]
    offers[0] = replace(offers[0], network_source_hash='sha256:aaaa')
    decision = apply_pdep_qm_budget(offers)
    record = build_pdep_budget_record(decision, iteration=0)
    assert record.network_outcomes[0].network_source_hash is None


# --- 24. Repeated network with disagreeing method resolves to None, no raise (weaker signal) -----

def test_apply_pdep_qm_budget_conflicting_method_resolves_to_none_without_raising():
    """Test that two offers of one network disagreeing on method (a weaker signal than identity,
    per combine()) resolve to None rather than raising -- method disagreement alone must never be
    treated as the fatal 'these are different networks' condition network_source_hash disagreement is."""
    offers = [_make_selection('netA', n_ts=1, delta_ln_k=1.0),
             _make_selection('netA', n_ts=1, delta_ln_k=2.0)]
    offers[0] = replace(offers[0], method='MSC')
    offers[1] = replace(offers[1], method='RRKM')
    decision = apply_pdep_qm_budget(offers)  # must not raise
    record = build_pdep_budget_record(decision, iteration=0)
    assert record.network_outcomes[0].method is None


# --- 25. A not-evaluated network's recorded rank/outcome survives to the final record ------------

def test_not_evaluated_network_survives_to_the_record_with_its_walk_rank():
    """Test that a not-evaluated network is not silently dropped between apply_pdep_qm_budget and
    build_pdep_budget_record -- it must appear in the final record, refused, with the rank it held
    in the walk (not omitted the way 'nothing to rank it by' would wrongly imply, per the review's
    finding that selection_rank_key does give it a tier)."""
    selections = [_make_selection('netA', n_ts=1, delta_ln_k=5.0),
                 _make_selection('netB', n_ts=1, evaluation_status=EVALUATION_STATUS_NOT_EVALUATED)]
    decision = apply_pdep_qm_budget(selections)
    record = build_pdep_budget_record(decision, iteration=0)
    outcomes_by_id = {outcome.network_id: outcome for outcome in record.network_outcomes}
    assert outcomes_by_id['netB'].outcome == BUDGET_OUTCOME_REFUSED
    assert outcomes_by_id['netB'].reason_code == BUDGET_SKIP_NOT_EVALUATED
    assert outcomes_by_id['netB'].rank == 1


# --- 26. PDepBudgetNetworkOutcome direct-construction validation, one invalid field at a time -----

def test_network_outcome_rejects_negative_cost():
    """Test that a negative cost is rejected -- a negative transition-state count can never be a
    real charge."""
    with pytest.raises(ValueError, match='cost'):
        PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_ADMITTED, cost=-1, rank=0)


def test_network_outcome_rejects_boolean_cost():
    """Test that a bool cost is rejected -- bool is an int subclass in Python, and
    outcome=True/False must never silently pass as a valid cost."""
    with pytest.raises(ValueError, match='cost'):
        PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_ADMITTED, cost=True, rank=0)


def test_network_outcome_rejects_non_int_cost():
    """Test that a non-int cost (e.g. a float) is rejected."""
    with pytest.raises(ValueError, match='cost'):
        PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_ADMITTED, cost=1.5, rank=0)


def test_network_outcome_rejects_negative_rank():
    """Test that a negative rank is rejected -- rank is a position in a walk and can never be
    negative."""
    with pytest.raises(ValueError, match='rank'):
        PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_ADMITTED, cost=1, rank=-1)


def test_network_outcome_rejects_boolean_rank():
    """Test that a bool rank is rejected, for the same reason bool cost is."""
    with pytest.raises(ValueError, match='rank'):
        PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_ADMITTED, cost=1, rank=False)


def test_network_outcome_rejects_invalid_reason_code():
    """Test that a refused outcome with a reason_code outside VALID_BUDGET_SKIP_REASON_CODES is
    rejected -- an unrecognized code would silently break anything downstream that switches on it."""
    with pytest.raises(ValueError, match='reason_code'):
        PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_REFUSED, cost=1, rank=0,
                                 reason_code='not_a_real_code', reason='some reason')


def test_network_outcome_rejects_empty_reason():
    """Test that a refused outcome with an empty-string reason is rejected -- an empty reason is
    indistinguishable from a reason that was silently dropped."""
    with pytest.raises(ValueError, match='reason'):
        PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_REFUSED, cost=1, rank=0,
                                 reason_code=BUDGET_SKIP_NOT_EVALUATED, reason='')


def test_network_outcome_rejects_non_string_network_id():
    """Test that a non-str, truthy network_id (e.g. an int) is rejected -- network_id must be a str
    or a falsy value (None/''), never some other truthy non-str type masquerading as an identifier."""
    with pytest.raises(ValueError, match='network_id'):
        PDepBudgetNetworkOutcome(network_id=12345, outcome=BUDGET_OUTCOME_ADMITTED, cost=1, rank=0)


def test_network_outcome_rejects_non_string_method():
    """Test that a non-str, non-None method is rejected."""
    with pytest.raises(ValueError, match='method'):
        PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_ADMITTED, cost=1, rank=0, method=123)


def test_network_outcome_rejects_malformed_network_source_hash():
    """Test that a non-str, non-None network_source_hash is rejected."""
    with pytest.raises(ValueError, match='network_source_hash'):
        PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_ADMITTED, cost=1, rank=0,
                                 network_source_hash=b'sha256:deadbeef')


def test_network_outcome_rejects_unnamed_offer_index_alongside_real_network_id():
    """Test that unnamed_offer_index is rejected alongside a real, truthy network_id -- it is only
    meaningful as the identity of an UNNAMED network, and setting both would leave two conflicting
    notions of identity on one outcome."""
    with pytest.raises(ValueError, match='unnamed_offer_index'):
        PDepBudgetNetworkOutcome(network_id='netA', outcome=BUDGET_OUTCOME_ADMITTED, cost=1, rank=0,
                                 unnamed_offer_index=0)


# --- 27. PDepBudgetRecord direct-construction validation, one invalid field at a time -------------

def _make_outcome(network_id, outcome=BUDGET_OUTCOME_ADMITTED, cost=1, rank=0, **kwargs):
    """Build a minimal, valid PDepBudgetNetworkOutcome, defaulting to admitted."""
    return PDepBudgetNetworkOutcome(network_id=network_id, outcome=outcome, cost=cost, rank=rank, **kwargs)


def test_record_rejects_wrong_schema_version():
    """Test that a schema_version other than BUDGET_RECORD_SCHEMA_VERSION is rejected -- this field
    exists precisely so an old, differently-shaped record can never be silently read as current."""
    with pytest.raises(ValueError, match='schema_version'):
        PDepBudgetRecord(iteration=0, max_transition_states=None, max_networks=None, total_cost=1,
                        network_outcomes=(_make_outcome('netA'),), schema_version=BUDGET_RECORD_SCHEMA_VERSION + 1)


def test_record_rejects_wrong_algorithm_version():
    """Test that an algorithm_version other than BUDGET_ALGORITHM_VERSION is rejected -- for the
    same reason as schema_version, but for the admission/refusal LOGIC rather than the shape."""
    with pytest.raises(ValueError, match='algorithm_version'):
        PDepBudgetRecord(iteration=0, max_transition_states=None, max_networks=None, total_cost=1,
                        network_outcomes=(_make_outcome('netA'),), algorithm_version=BUDGET_ALGORITHM_VERSION + 1)


def test_record_rejects_non_outcome_entries():
    """Test that network_outcomes containing something other than a PDepBudgetNetworkOutcome (e.g.
    a plain dict) is rejected, rather than silently accepted and only failing later at serialization."""
    with pytest.raises(ValueError, match='PDepBudgetNetworkOutcome'):
        PDepBudgetRecord(iteration=0, max_transition_states=None, max_networks=None, total_cost=0,
                        network_outcomes=({'network_id': 'netA'},))


def test_record_rejects_duplicate_rank():
    """Test that two network_outcomes sharing the same rank are rejected -- rank is meant to be
    each network's unique position in one walk."""
    with pytest.raises(ValueError, match='rank'):
        PDepBudgetRecord(iteration=0, max_transition_states=None, max_networks=None, total_cost=2,
                        network_outcomes=(_make_outcome('netA', rank=0), _make_outcome('netB', rank=0)))


def test_record_rejects_gapped_rank():
    """Test that ranks with a gap (0, 2 but no 1) are rejected -- a gap would mean some
    considered network's outcome went missing between the walk and the record."""
    with pytest.raises(ValueError, match='rank'):
        PDepBudgetRecord(iteration=0, max_transition_states=None, max_networks=None, total_cost=2,
                        network_outcomes=(_make_outcome('netA', rank=0), _make_outcome('netB', rank=2)))


def test_record_rejects_total_cost_not_matching_sum_of_admitted_costs():
    """Test that a total_cost disagreeing with the sum of admitted outcomes' costs is rejected --
    this cross-check is what would catch a caller building a record by hand from a stale or
    mismatched decision."""
    with pytest.raises(ValueError, match='total_cost'):
        PDepBudgetRecord(iteration=0, max_transition_states=None, max_networks=None, total_cost=999,
                        network_outcomes=(_make_outcome('netA', cost=1),))


def test_record_rejects_duplicate_unnamed_network_identity():
    """Test that two unnamed outcomes (falsy network_id) sharing the same unnamed_offer_index are
    still rejected as a duplicate identity -- the composite key must actually distinguish, not just
    exist."""
    with pytest.raises(ValueError, match='identity'):
        PDepBudgetRecord(iteration=0, max_transition_states=None, max_networks=None, total_cost=0,
                        network_outcomes=(_make_outcome(None, rank=0, outcome=BUDGET_OUTCOME_REFUSED,
                                                        reason_code=BUDGET_SKIP_NOT_EVALUATED,
                                                        reason='x', unnamed_offer_index=0),
                                         _make_outcome(None, rank=1, outcome=BUDGET_OUTCOME_REFUSED,
                                                      reason_code=BUDGET_SKIP_NOT_EVALUATED,
                                                      reason='x', unnamed_offer_index=0)))


def test_build_refuses_a_decision_whose_considered_does_not_account_for_its_admissions():
    """Test that a decision claiming admissions it carries no consideration detail for is refused
    rather than rendered as a record saying nothing was considered at all.

    ``considered`` is defaulted, so a decision can legitimately be constructed without it -- but a
    record built from such a decision would be an authoritative-looking lie, silently reporting an
    empty field where two networks were in fact admitted. The record is where authority is claimed,
    so that is where the incoherence has to be caught."""
    decision = PDepBudgetDecision(admitted_indices=(0, 1), skipped=tuple(), total_cost=0)
    with pytest.raises(ValueError, match='considered'):
        build_pdep_budget_record(decision, iteration=3)


def test_build_refuses_a_decision_whose_considered_does_not_account_for_its_refusals():
    """Test that the same guard fires on the refusal side, where a dropped consideration would erase
    the very refusal the record exists to preserve."""
    skip = PDepBudgetSkip(network_id='netA', cost=1, remaining_transition_states=None,
                          reason_code=BUDGET_SKIP_NOT_EVALUATED, reason='it could not be evaluated')
    decision = PDepBudgetDecision(admitted_indices=tuple(), skipped=(skip,), total_cost=0)
    with pytest.raises(ValueError, match='considered'):
        build_pdep_budget_record(decision, iteration=3)
