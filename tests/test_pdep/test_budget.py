#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_budget
"""

import pytest

from t3.pdep.budget import (PDepBudgetDecision,
                            PDepBudgetSkip,
                            apply_pdep_qm_budget,
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
