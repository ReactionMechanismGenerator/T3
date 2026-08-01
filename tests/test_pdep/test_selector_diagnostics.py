#!/usr/bin/env python3
# encoding: utf-8

"""
Tests for t3.pdep.selector.select_from_sa_dict_with_diagnostics.

The selector distinguishes its nine refusal outcomes only by the prose it appends to
``selection.warnings``. The diagnostics wrapper reports a machine-readable code alongside the
decision so the durable assessment record can COUNT them; these tests pin each code to the
condition that actually produces it, because a code that silently reports the wrong outcome is
worse than no code at all -- it would make a miscount look authoritative.
"""

import pytest

from t3.pdep.parser import parse_pdep_network_text
from t3.pdep.reason_codes import (PRE_SELECTOR_REASON_CODES,
                                  REASON_EVALUATED_NO_UNCERTAIN_TS,
                                  REASON_QUALIFIED_UNCERTAIN_TS,
                                  REASON_SELECTOR_DIRECTION_ENTRY_MALFORMED,
                                  REASON_SELECTOR_DIRECTION_UNRESOLVED,
                                  REASON_SELECTOR_MALFORMED_CONDITIONS_NO_TS_ROWS,
                                  REASON_SELECTOR_NEGATIVE_INCOMPLETE_DATA,
                                  REASON_SELECTOR_NO_TS_ROWS,
                                  REASON_SELECTOR_NO_USABLE_TS_ROWS,
                                  REASON_SELECTOR_SA_PAYLOAD_MALFORMED,
                                  REASON_SELECTOR_TS_PROVENANCE_UNASSESSED,
                                  REASON_SELECTOR_TS_RESPONSE_BELOW_FLOOR,
                                  SELECTION_BEARING_REASON_CODES,
                                  )
from t3.pdep.selector import (EVALUATION_STATUS_EVALUATED,
                              EVALUATION_STATUS_NOT_EVALUATED,
                              select_from_sa_dict,
                              select_from_sa_dict_with_diagnostics,
                              )


def network_text(comment):
    """A one-reaction network whose single path reaction carries ``comment`` as its kinetics."""
    return f'''
reaction(
    label = 'reactionA',
    reactants = ['S1', 'S2'],
    products = ['S3'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""{comment}"""),
)
'''


UNCERTAIN_COMMENT = 'Estimated using template [x]'
CERTAIN_COMMENT = "Reaction library: primaryNitrogenLibrary"
CONDITION = (300.0, 'K', 1.0, 'bar')
REACTION = 'A + B <=> C'
# Comfortably above the ~1.195e-7 mol/J coefficient floor for these thresholds.
ABOVE_FLOOR = 1.0e-3


def diagnose(sa_dict, comment=UNCERTAIN_COMMENT, network_reaction=REACTION, network_id='synthetic'):
    """Run the diagnostics wrapper over a one-reaction synthetic network."""
    network = parse_pdep_network_text(text=network_text(comment), network_id=network_id)
    return select_from_sa_dict_with_diagnostics(sa_dict=sa_dict,
                                                network=network,
                                                network_reaction=network_reaction,
                                                relative_threshold=0.001,
                                                min_delta_ln_k=1e-3,
                                                perturbation=8368.0,
                                                )


# --- The two verdicts ----------------------------------------------------------------------------

def test_a_sensitive_uncertain_transition_state_reports_the_qualified_code():
    selection, reason_code, secondary = diagnose({REACTION: {CONDITION: {'(TS) TS1': ABOVE_FLOOR}}})
    assert selection.qualified is True
    assert selection.evaluation_status == EVALUATION_STATUS_EVALUATED
    assert [entry.ts_label for entry in selection.uncertain_path_reactions] == ['TS1']
    assert reason_code == REASON_QUALIFIED_UNCERTAIN_TS
    assert secondary == ()


def test_a_sensitive_but_certain_transition_state_reports_the_evaluated_negative_code():
    selection, reason_code, secondary = diagnose({REACTION: {CONDITION: {'(TS) TS1': ABOVE_FLOOR}}},
                                                 comment=CERTAIN_COMMENT)
    assert selection.qualified is False
    assert selection.evaluation_status == EVALUATION_STATUS_EVALUATED
    assert [entry.ts_label for entry in selection.selected_ts] == ['TS1']
    assert selection.uncertain_path_reactions == []
    assert reason_code == REASON_EVALUATED_NO_UNCERTAIN_TS
    assert secondary == ()


# --- The nine refusals ---------------------------------------------------------------------------
#
# Each asserts the decision state alongside the code: a code that agreed with nothing in the
# selection it accompanies would be a label rather than a diagnosis.

def test_a_non_mapping_sensitivity_payload_reports_the_malformed_payload_code():
    # Reachable only through the selector's PUBLIC entry points. T3's own funnel inspects
    # `structures` first and so can never get a non-mapping payload this far -- which is why this
    # code is named for the selector and the funnel has its own `sa_output_malformed`.
    selection, reason_code, _ = diagnose(['not', 'a', 'mapping'])
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert selection.selected_ts == []
    assert reason_code == REASON_SELECTOR_SA_PAYLOAD_MALFORMED


def test_a_reaction_absent_from_the_sensitivity_output_reports_the_unresolved_direction_code():
    selection, reason_code, _ = diagnose({'X <=> Y': {CONDITION: {'(TS) TS1': ABOVE_FLOOR}}})
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert selection.direction_key is None
    assert selection.network_reactions_examined == 0
    assert reason_code == REASON_SELECTOR_DIRECTION_UNRESOLVED


def test_a_non_mapping_condition_block_reports_the_malformed_direction_entry_code():
    selection, reason_code, _ = diagnose({REACTION: 'not a mapping of conditions'})
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert selection.direction_key == REACTION
    assert reason_code == REASON_SELECTOR_DIRECTION_ENTRY_MALFORMED


def test_a_condition_with_only_well_rows_reports_the_no_transition_state_rows_code():
    selection, reason_code, _ = diagnose({REACTION: {CONDITION: {'well1': 1.0}}})
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert selection.selected_ts == []
    assert any('No transition state sensitivity rows were found' in warning
               for warning in selection.warnings)
    assert reason_code == REASON_SELECTOR_NO_TS_ROWS


def test_malformed_conditions_and_no_transition_state_rows_report_their_own_code():
    # Distinct from the case above: rows could not be read at all, rather than read and found to
    # contain no transition state. Conflating them would hide a parsing failure as an absence.
    selection, reason_code, _ = diagnose({REACTION: {CONDITION: 'not a mapping of entries'}})
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert reason_code == REASON_SELECTOR_MALFORMED_CONDITIONS_NO_TS_ROWS
    assert reason_code != REASON_SELECTOR_NO_TS_ROWS


def test_transition_state_rows_that_are_all_non_finite_report_the_no_usable_rows_code():
    selection, reason_code, _ = diagnose({REACTION: {CONDITION: {'(TS) TS1': float('nan'),
                                                                 '(TS) TS2': float('inf')}}})
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert selection.selected_ts == []
    assert reason_code == REASON_SELECTOR_NO_USABLE_TS_ROWS


def test_transition_state_rows_below_the_absolute_floor_report_the_below_floor_code():
    selection, reason_code, _ = diagnose({REACTION: {CONDITION: {'(TS) TS1': 1.0e-18,
                                                                 '(TS) TS2': 1.0e-15}}})
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert selection.selected_ts == []
    assert any('did not reach the absolute floor' in warning for warning in selection.warnings)
    assert reason_code == REASON_SELECTOR_TS_RESPONSE_BELOW_FLOOR


def test_a_negative_verdict_resting_on_discarded_rows_reports_the_incomplete_data_code():
    # One usable, certain, above-floor row would give a real negative -- but a non-finite row was
    # discarded alongside it, so the negative rests on an incomplete picture.
    selection, reason_code, _ = diagnose({REACTION: {CONDITION: {'(TS) TS1': ABOVE_FLOOR,
                                                                 '(TS) TS2': float('nan')}}},
                                         comment=CERTAIN_COMMENT)
    assert selection.qualified is False
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert [entry.ts_label for entry in selection.selected_ts] == ['TS1']
    assert reason_code == REASON_SELECTOR_NEGATIVE_INCOMPLETE_DATA


def test_a_selected_transition_state_with_no_joined_path_reaction_reports_the_unassessed_code():
    # '(TS) TS_ORPHAN' joins to no path reaction, so its provenance is unknown rather than certain;
    # a negative built on it would be a fail-OPEN.
    selection, reason_code, _ = diagnose({REACTION: {CONDITION: {'(TS) TS_ORPHAN': ABOVE_FLOOR}}})
    assert selection.qualified is False
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert [entry.uncertain for entry in selection.selected_ts] == [None]
    assert reason_code == REASON_SELECTOR_TS_PROVENANCE_UNASSESSED


def test_when_two_refusals_both_apply_neither_is_dropped():
    # TS_ORPHAN is selected but unjoined (provenance unassessed) AND a non-finite row was discarded
    # alongside it. These are INDEPENDENT defects in the evidence -- a transition state can fail to
    # join whether or not other rows were unreadable -- so reporting only the first would undercount
    # the second wherever these are tallied.
    selection, reason_code, secondary = diagnose(
        {REACTION: {CONDITION: {'(TS) TS_ORPHAN': ABOVE_FLOOR, '(TS) TS2': float('nan')}}})
    assert selection.qualified is False
    assert any(entry.uncertain is None for entry in selection.selected_ts)
    assert reason_code == REASON_SELECTOR_NEGATIVE_INCOMPLETE_DATA
    assert secondary == (REASON_SELECTOR_TS_PROVENANCE_UNASSESSED,)


# --- Contract-level properties -------------------------------------------------------------------

ALL_SCENARIOS = [
    ({REACTION: {CONDITION: {'(TS) TS1': ABOVE_FLOOR}}}, UNCERTAIN_COMMENT),
    ({REACTION: {CONDITION: {'(TS) TS1': ABOVE_FLOOR}}}, CERTAIN_COMMENT),
    (['not', 'a', 'mapping'], UNCERTAIN_COMMENT),
    ({'X <=> Y': {}}, UNCERTAIN_COMMENT),
    ({REACTION: 'not a mapping'}, UNCERTAIN_COMMENT),
    ({REACTION: {CONDITION: {'well1': 1.0}}}, UNCERTAIN_COMMENT),
    ({REACTION: {CONDITION: 'not a mapping'}}, UNCERTAIN_COMMENT),
    ({REACTION: {CONDITION: {'(TS) TS1': float('nan')}}}, UNCERTAIN_COMMENT),
    ({REACTION: {CONDITION: {'(TS) TS1': 1.0e-18}}}, UNCERTAIN_COMMENT),
    ({REACTION: {CONDITION: {'(TS) TS_ORPHAN': ABOVE_FLOOR}}}, UNCERTAIN_COMMENT),
    ({REACTION: {CONDITION: {'(TS) TS_ORPHAN': ABOVE_FLOOR, '(TS) TS2': float('nan')}}},
     UNCERTAIN_COMMENT),
]


@pytest.mark.parametrize('sa_dict, comment', ALL_SCENARIOS)
def test_every_reported_code_is_selection_bearing_and_never_none(sa_dict, comment):
    # The assessment record refuses a selector-produced code without a nested selection and refuses
    # a pre-selector code WITH one, so a wrapper that returned a pre-selector code here would make
    # every assessment built from it unconstructable.
    selection, reason_code, secondary = diagnose(sa_dict, comment=comment)
    assert reason_code is not None
    assert reason_code in SELECTION_BEARING_REASON_CODES
    assert reason_code not in PRE_SELECTOR_REASON_CODES
    for code in secondary:
        assert code in SELECTION_BEARING_REASON_CODES
        assert code != reason_code


@pytest.mark.parametrize('sa_dict, comment', ALL_SCENARIOS)
def test_the_code_agrees_with_the_evaluation_status_it_accompanies(sa_dict, comment):
    selection, reason_code, _ = diagnose(sa_dict, comment=comment)
    is_verdict = reason_code in (REASON_QUALIFIED_UNCERTAIN_TS, REASON_EVALUATED_NO_UNCERTAIN_TS)
    assert is_verdict == (selection.evaluation_status == EVALUATION_STATUS_EVALUATED)
    assert (reason_code == REASON_QUALIFIED_UNCERTAIN_TS) == selection.qualified


@pytest.mark.parametrize('sa_dict, comment', ALL_SCENARIOS)
def test_the_plain_entry_point_returns_exactly_what_the_diagnostics_one_decides(sa_dict, comment):
    # One implementation behind both, so the decision cannot drift from its own diagnosis.
    plain = select_from_sa_dict(
        sa_dict=sa_dict,
        network=parse_pdep_network_text(text=network_text(comment), network_id='synthetic'),
        network_reaction=REACTION, relative_threshold=0.001, min_delta_ln_k=1e-3, perturbation=8368.0)
    diagnosed, _, _ = diagnose(sa_dict, comment=comment)
    assert plain.as_dict() == diagnosed.as_dict()
