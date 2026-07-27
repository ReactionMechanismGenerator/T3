#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_selector
"""

import json
import math
import os

import pytest
from arc.common import read_yaml_file

from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.parser import parse_pdep_network_file, parse_pdep_network_text
from t3.pdep.selector import (SensitiveTransitionState,
                              coefficient_floor,
                              resolve_direction_key,
                              select_from_sa_dict,
                              select_sensitive_wells,
                              )

PDEP_NETWORK_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'iteration_1', 'RMG', 'pdep')
# Deliberately kept OUTSIDE pdep_network/: test_main.py::test_determine_species_from_pdep_network
# sets t3.paths['PDep SA'] to tests/data/pdep_network/iteration_1/PDep_SA and its teardown does
# shutil.rmtree() on that path, which would delete this fixture if it lived there.
SA_PATH = os.path.join(TEST_DATA_BASE_PATH, 'pdep_sa', 'network4_2_MSC', 'sa_coefficients.yml')

TARGET_REACTION = 'C1rad(2) + C3ene(27) <=> C2ene(29) + C2rad(3)'
CERTAIN_TS_LABELS = {'TS2', 'TS3', 'TS5', 'TS8', 'TS9', 'TS10'}
UNCERTAIN_TS_LABELS = {'TS1', 'TS4', 'TS6', 'TS7'}


@pytest.fixture(scope='module')
def sa_dict():
    """The real Arkane sensitivity dictionary for network4_2."""
    return read_yaml_file(path=SA_PATH)


@pytest.fixture(scope='module')
def network():
    """The real parsed network4_2 network."""
    return parse_pdep_network_file(path=os.path.join(PDEP_NETWORK_DIR, 'network4_2.py'))


# --- 1. Integration: the non-qualifying target reaction --------------------------------------

def test_target_reaction_does_not_qualify(sa_dict, network):
    """Test that the target reaction does not qualify and its selected TS set is exactly the certain ones."""
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001)
    assert selection.qualified is False
    labels = {entry.ts_label for entry in selection.selected_ts}
    assert labels == CERTAIN_TS_LABELS


# --- 2. Integration: sweep of all 45 keys -----------------------------------------------------

def test_sweep_all_keys_yields_20_qualified(sa_dict, network):
    """Test that sweeping all 45 network reaction keys yields exactly 20 qualified networks."""
    keys = [key for key in sa_dict.keys() if key != 'structures']
    assert len(keys) == 45
    qualified_count = 0
    for key in keys:
        selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=key,
                                        relative_threshold=0.001)
        if selection.qualified:
            qualified_count += 1
    assert qualified_count == 20


# --- 3. Integration: a qualifying key's evidence ----------------------------------------------

def test_qualifying_key_evidence(sa_dict, network):
    """Test the evidence recorded for a qualifying key: non-empty, all uncertain, labels a subset
    of the known uncertain TS set, and a reason mentioning the network id."""
    keys = [key for key in sa_dict.keys() if key != 'structures']
    qualifying_key = None
    for key in keys:
        selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=key,
                                        relative_threshold=0.001)
        if selection.qualified:
            qualifying_key = key
            break
    assert qualifying_key is not None
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=qualifying_key,
                                    relative_threshold=0.001)
    assert selection.uncertain_path_reactions
    assert all(entry.uncertain is True for entry in selection.uncertain_path_reactions)
    assert set(selection.uncertain_ts_labels()) <= UNCERTAIN_TS_LABELS
    assert network.network_id in selection.reason()


# --- 4. thresholds round-trip ------------------------------------------------------------------

def test_thresholds_round_trip(sa_dict, network):
    """Test that ``thresholds`` records the values passed in plus the derived ``coefficient_floor``."""
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001, min_delta_ln_k=1e-3, perturbation=8368.0)
    assert selection.thresholds['relative_threshold'] == 0.001
    assert selection.thresholds['min_delta_ln_k'] == 1e-3
    assert selection.thresholds['perturbation'] == 8368.0
    assert selection.thresholds['coefficient_floor'] == pytest.approx(1e-3 / 8368.0)


# --- 5. Path reaction labels are not unique -----------------------------------------------------

def test_path_reaction_labels_not_unique(network):
    """Test directly that path reaction labels are not unique within network4_2 (e.g. 'reaction1' appears
    for both TS1 and TS10)."""
    labels = [path_reaction.label for path_reaction in network.path_reactions]
    assert len(set(labels)) < len(labels)
    assert labels.count('reaction1') == 2


# --- 6. coefficient_floor ------------------------------------------------------------------------

def test_coefficient_floor_value():
    """Test that ``coefficient_floor`` returns ``min_delta_ln_k / perturbation``."""
    assert coefficient_floor(min_delta_ln_k=1e-3, perturbation=8368.0) == pytest.approx(1e-3 / 8368.0)


def test_coefficient_floor_raises_for_non_positive_perturbation():
    """Test that ``coefficient_floor`` raises ``ValueError`` for a non-positive perturbation."""
    with pytest.raises(ValueError):
        coefficient_floor(min_delta_ln_k=1e-3, perturbation=0)
    with pytest.raises(ValueError):
        coefficient_floor(min_delta_ln_k=1e-3, perturbation=-1.0)


# --- 7. The absolute floor bites -------------------------------------------------------------------

SYNTHETIC_NETWORK_TEXT = '''
reaction(
    label = 'reactionA',
    reactants = ['S1', 'S2'],
    products = ['S3'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [x]"""),
)
'''


def test_absolute_floor_rejects_denormal_rows_despite_relative_gate():
    """Test that the absolute floor rejects a condition where every TS coefficient is a real,
    finite, measured response (~1e-18/1e-15, i.e. a denormal, NOT a structural zero) that is
    smaller than the floor, even though one coefficient is 1000x the others and would pass the
    relative gate alone. This USED TO be treated as an answer to criterion (b) -- a real response
    too small to act on, evaluation_status staying 'evaluated' with qualified=False (the old case
    (d) of the below-floor split). The first PDep-QM trial run showed that reasoning does not
    survive contact with real data: a live network's TS coefficients (~1e-16-1e-14 mol/J) sat many
    orders of magnitude below the floor (~1.2e-7 mol/J) while its well coefficients sat at or above
    it, even though Arkane's log confirmed it had perturbed the TS -- a below-floor TS reading is a
    dead instrument, not a small answer. Cases (c) and (d) are now merged: evaluation_status must be
    'not_evaluated', and the warning must still say the largest coefficient did not reach the floor,
    but must now also say criterion (b) cannot be evaluated."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_floor')
    # TS2 is 1000x TS1, so the relative gate alone (cutoff = 0.001 * 1e-15 = 1e-18) would select
    # both; the absolute floor (~1.195e-7 mol/J) must reject both, since neither reaches it.
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 1.0e-18, '(TS) TS2': 1.0e-15}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001, min_delta_ln_k=1e-3, perturbation=8368.0)
    assert selection.selected_ts == []
    assert selection.qualified is False
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert any('did not reach the absolute floor' in warning for warning in selection.warnings)
    assert any('cannot be evaluated' in warning for warning in selection.warnings)


def test_below_floor_no_ts_rows_at_all_is_not_evaluated():
    """Test case (a) of the below-floor split: a condition with no TS-prefixed rows at all (only
    well entries) is NOT a below-floor measurement -- there was nothing to measure -- so
    evaluation_status must be 'not_evaluated', with wording distinct from the 'seen but unusable'
    and 'below floor' cases."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_no_ts_rows')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {'well1': 1.0}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert any('No transition state sensitivity rows were found' in warning
               for warning in selection.warnings)
    assert not any('did not reach the absolute floor' in warning for warning in selection.warnings)
    assert not any('exactly zero' in warning for warning in selection.warnings)


def test_below_floor_ts_rows_seen_but_none_usable_is_not_evaluated():
    """Test case (b) of the below-floor split: TS-prefixed rows were seen but every one of them was
    non-finite, so none were usable -- the negative verdict this would otherwise produce is
    unsupported, so evaluation_status must be 'not_evaluated' with wording naming that rows were
    seen but unusable, distinct from case (a)'s 'none found' wording."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_none_usable')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': float('nan'),
                                                          '(TS) TS2': float('inf'),
                                                          }}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert any('none were usable' in warning for warning in selection.warnings)
    assert not any('No transition state sensitivity rows were found' in warning
                   for warning in selection.warnings)


def test_below_floor_exact_zero_is_not_evaluated():
    """Test case (c) of the below-floor split: usable TS rows exist and every one is EXACTLY 0.0,
    the structural-zero signature of an Arkane that perturbs a synthesized TS E0 without reaching
    the ILT rate expression. This is a dead instrument, not a measurement, so evaluation_status
    must be 'not_evaluated'. Cases (c) (exact zero) and the former case (d) (finite-but-below-floor)
    are now merged into one 'below the absolute floor' diagnosis, so the warning must report BOTH
    that the coefficient did not reach the floor AND that it was exactly zero."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_exact_zero')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 0.0, '(TS) TS2': 0.0}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001, min_delta_ln_k=1e-3, perturbation=8368.0)
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert any('exactly zero' in warning for warning in selection.warnings)
    assert any('did not reach the absolute floor' in warning for warning in selection.warnings)
    # The single most common cause of this diagnosis is an Arkane predating RMG-Py PR #2990, which is
    # what makes a TS perturbation reach an ILT-derived rate at all. Without naming it, the reader is
    # told their instrument is dead but not where the fix is -- and the likeliest reading of "this
    # network carries no TS signal" is "this network is uninteresting", which is the wrong action.
    assert any('#2990' in warning for warning in selection.warnings), \
        f'The diagnosis did not name the RMG-Py fix that resolves it: {selection.warnings}'


def test_below_floor_exact_zero_reports_well_row_response():
    """Test case (c)'s enriched message: when the TS channel is dead (every usable TS coefficient
    exactly 0.0) but a finite non-TS (well) row for the SAME direction key DID respond, the warning
    must name that the wells responded and report their largest absolute coefficient -- this is what
    makes the 'dead instrument, not a dead network' hypothesis credible, distinct from case (c)'s
    base wording alone."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_exact_zero_wells')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 0.0, '(TS) TS2': 0.0, 'well1': 5.0}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001, min_delta_ln_k=1e-3, perturbation=8368.0)
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert any('wells did respond' in warning and '5.000e+00' in warning for warning in selection.warnings)


def test_below_floor_exact_zero_notes_absence_of_well_rows_too():
    """Test case (c)'s enriched message when there are NO finite well rows either: the message must
    say so explicitly rather than silently omitting the well context."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_exact_zero_no_wells')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 0.0, '(TS) TS2': 0.0}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001, min_delta_ln_k=1e-3, perturbation=8368.0)
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert any('No finite non-transition-state (well) rows were found' in warning
               for warning in selection.warnings)


def test_below_floor_negative_unsupported_by_discarded_malformed_and_nonfinite_rows():
    """Test that a below-floor TS reading (case (c), merged) computed alongside a malformed
    condition block still reaches 'not_evaluated', and via the merged case's OWN diagnosis rather
    than the separate downstream 'rests on incomplete data' discard-check. Before the case (c)/(d)
    merge, this exact input reached NOT_EVALUATED only through that downstream check, since the old
    case (d) left evaluation_status 'evaluated' at dispatch time; now the merged case already sets
    NOT_EVALUATED before the discard-check's guard (``evaluation_status == EVALUATION_STATUS_EVALUATED``)
    can ever be true, so the 'rests on incomplete data' wording is not the reason for this verdict
    any more -- the malformed block is still surfaced, but only via the generic 'Ignored ... malformed
    condition' warning. The verdict itself (not_evaluated) is unchanged; only which code path -- and
    which message -- produces it has moved earlier."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_discarded_negative')
    sa_dict = {'A + B <=> C': {
        (300.0, 'K', 1.0, 'bar'): 'not a dict',
        (400.0, 'K', 1.0, 'bar'): {'(TS) TS1': 1.0e-18},
    }}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001, min_delta_ln_k=1e-3, perturbation=8368.0)
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert any('did not reach the absolute floor' in warning for warning in selection.warnings)
    assert any('malformed condition' in warning for warning in selection.warnings)
    assert not any('rests on incomplete data' in warning for warning in selection.warnings)


def test_below_floor_all_malformed_condition_blocks_distinct_from_no_rows_found():
    """Test case (1b): when EVERY condition block for a direction key is malformed (so
    ts_rows_seen == 0 for a reason distinct from case (a)), the diagnosis must say whether TS rows
    existed cannot be determined -- NOT case (a)'s 'no transition state sensitivity rows were found'
    wording, which asserts a settled fact this data cannot support."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_all_malformed')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): 'not a dict'}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert any('cannot be determined' in warning for warning in selection.warnings)
    assert not any('No transition state sensitivity rows were found' in warning
                   for warning in selection.warnings)


# --- 8. Non-finite handling -----------------------------------------------------------------------

def test_non_finite_ts_rows_ignored_and_warned():
    """Test that ``nan`` and ``inf`` TS rows are ignored and counted in a warning."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_nonfinite')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': float('nan'),
                                                          '(TS) TS2': float('inf'),
                                                          }}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.selected_ts == []
    assert any('non-finite' in warning and '2' in warning for warning in selection.warnings)


# --- 9. resolve_direction_key -------------------------------------------------------------------

def test_resolve_direction_key_exact_match():
    """Test that an exact key match returns the key unchanged, unambiguous, with no warnings."""
    sa_dict = {'A + B <=> C': dict(), 'structures': dict()}
    key, ambiguous, warnings = resolve_direction_key(sa_dict=sa_dict, network_reaction='A + B <=> C')
    assert key == 'A + B <=> C'
    assert ambiguous is False
    assert warnings == []


def test_resolve_direction_key_legalized_label_match():
    """Test that an ARC-legalized label (bracket vs parenthesis style) is matched with a warning,
    and is not reported as ambiguous."""
    sa_dict = {'A[1] + B <=> C': dict(), 'structures': dict()}
    key, ambiguous, warnings = resolve_direction_key(sa_dict=sa_dict, network_reaction='A(1) + B <=> C')
    assert key == 'A[1] + B <=> C'
    assert ambiguous is False
    assert len(warnings) == 1


def test_resolve_direction_key_only_reverse_present():
    """Test that when only the reverse direction is present, it is returned with ``direction_ambiguous``
    True and a warning."""
    sa_dict = {'C <=> A + B': dict(), 'structures': dict()}
    key, ambiguous, warnings = resolve_direction_key(sa_dict=sa_dict, network_reaction='A + B <=> C')
    assert key == 'C <=> A + B'
    assert ambiguous is True
    assert len(warnings) == 1


def test_resolve_direction_key_absent():
    """Test that an absent reaction returns ``(None, False, [warning])``."""
    sa_dict = {'X + Y <=> Z': dict(), 'structures': dict()}
    key, ambiguous, warnings = resolve_direction_key(sa_dict=sa_dict, network_reaction='A + B <=> C')
    assert key is None
    assert ambiguous is False
    assert len(warnings) == 1


# --- 10. Reaction absent from SA dict -------------------------------------------------------------

def test_select_from_sa_dict_reaction_absent_does_not_raise():
    """Test that a network reaction absent from the SA dict returns ``qualified is False`` with a
    warning and an empty ``selected_ts``, without raising."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_absent')
    sa_dict = {'X + Y <=> Z': dict(), 'structures': dict()}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.qualified is False
    assert selection.selected_ts == []
    assert selection.warnings


# --- 11. TS join failure --------------------------------------------------------------------------

def test_ts_join_failure_recorded_as_evidence_with_none_uncertain():
    """Test that a TS entry with no matching ``transitionState`` in the network produces an evidence
    record with ``uncertain is None``, a warning, and does NOT cause qualification -- but since the
    only evidence is an unassessed TS (no ``True`` entry anywhere), the negative verdict is
    unsupported: ``evaluation_status`` must be ``EVALUATION_STATUS_NOT_EVALUATED``, not a confident
    ``qualified=False``."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_unjoined')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS_UNKNOWN': 1.0}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert len(selection.selected_ts) == 1
    entry = selection.selected_ts[0]
    assert isinstance(entry, SensitiveTransitionState)
    assert entry.ts_label == 'TS_UNKNOWN'
    assert entry.uncertain is None
    assert entry.path_reaction_label is None
    assert any('Could not join' in warning for warning in selection.warnings)
    assert selection.qualified is False
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert any('provenance could not be assessed' in warning and 'criterion (b)' in warning
               for warning in selection.warnings)


def test_ts_join_failure_mixed_with_certain_entry_still_not_evaluated():
    """Test that when the selected TS set has one certain (``uncertain=False``) entry and one
    unjoinable (``uncertain=None``) entry -- and no ``True`` entry anywhere -- the aggregate verdict
    is still unsupported: a certain-but-negative entry does not paper over an unassessed one."""
    network = parse_pdep_network_text(text=SHARED_TS_CERTAIN_ONLY_NETWORK_TEXT,
                                      network_id='synthetic_mixed_none_false')
    sa_dict = {'S1 + S2 <=> S3': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 1.0, '(TS) TS_UNKNOWN': 1.0}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='S1 + S2 <=> S3',
                                    relative_threshold=0.001)
    uncertain_flags = {entry.ts_label: entry.uncertain for entry in selection.selected_ts}
    assert uncertain_flags == {'TS1': False, 'TS_UNKNOWN': None}
    assert selection.qualified is False
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED


# --- 12. Shared TS --------------------------------------------------------------------------------

SHARED_TS_NETWORK_TEXT = '''
reaction(
    label = 'reactionA',
    reactants = ['S1', 'S2'],
    products = ['S3'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Exact match found for rate rule [x]"""),
)

reaction(
    label = 'reactionB',
    reactants = ['S3'],
    products = ['S4', 'S5'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [x]"""),
)
'''


def test_shared_ts_both_assessed_and_network_qualifies():
    """Test a TS shared by two path reactions (one certain, one uncertain): both are assessed, a
    warning about the shared label is recorded, and the network qualifies."""
    network = parse_pdep_network_text(text=SHARED_TS_NETWORK_TEXT, network_id='synthetic_shared')
    sa_dict = {'S1 + S2 <=> S3': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 1.0}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='S1 + S2 <=> S3',
                                    relative_threshold=0.001)
    assert len(selection.selected_ts) == 2
    labels = {entry.path_reaction_label for entry in selection.selected_ts}
    assert labels == {'reactionA', 'reactionB'}
    uncertain_flags = {entry.path_reaction_label: entry.uncertain for entry in selection.selected_ts}
    assert uncertain_flags == {'reactionA': False, 'reactionB': True}
    assert any('shared by' in warning for warning in selection.warnings)
    assert selection.qualified is True
    assert selection.evaluation_status == EVALUATION_STATUS_EVALUATED


SHARED_TS_CERTAIN_ONLY_NETWORK_TEXT = '''
reaction(
    label = 'reactionA',
    reactants = ['S1', 'S2'],
    products = ['S3'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Exact match found for rate rule [x]"""),
)
'''


def test_all_certain_no_true_no_none_stays_evaluated_and_unqualified():
    """Test that when every selected TS is assessed and certain (``uncertain=False``, no ``True`` and
    no ``None`` anywhere), the negative verdict IS supported: ``evaluation_status`` stays
    ``EVALUATION_STATUS_EVALUATED`` and ``qualified`` stays ``False`` -- unchanged from before FIX 2."""
    network = parse_pdep_network_text(text=SHARED_TS_CERTAIN_ONLY_NETWORK_TEXT,
                                      network_id='synthetic_certain_only')
    sa_dict = {'S1 + S2 <=> S3': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 1.0}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='S1 + S2 <=> S3',
                                    relative_threshold=0.001)
    assert len(selection.selected_ts) == 1
    assert selection.selected_ts[0].uncertain is False
    assert selection.qualified is False
    assert selection.evaluation_status == EVALUATION_STATUS_EVALUATED


# --- 13. select_sensitive_wells ranks by absolute value ---------------------------------------------

def test_select_sensitive_wells_ranks_by_absolute_value():
    """Test that ``select_sensitive_wells`` selects a strongly negative coefficient (ranking is by
    absolute value, not signed value), excludes ``(TS) `` rows, and applies the floor."""
    entries_by_condition = {
        (300.0, 'K', 1.0, 'bar'): {
            'well_negative': -1.0,
            'well_small': 1e-4,
            '(TS) TS1': -2.0,
            'well_below_floor': 1e-18,
        },
    }
    sensitive_wells = select_sensitive_wells(entries_by_condition=entries_by_condition,
                                            relative_threshold=0.5,
                                            min_delta_ln_k=1e-3,
                                            perturbation=8368.0)
    assert 'well_negative' in sensitive_wells
    assert '(TS) TS1' not in sensitive_wells
    assert 'well_below_floor' not in sensitive_wells
    assert 'well_small' not in sensitive_wells
    assert sensitive_wells['well_negative'] == [(300.0, 'K', 1.0, 'bar')]


# --- FIX4: malformed/corrupted SA data is rejected, not crashed on ------------------------------------

def test_select_from_sa_dict_tolerates_non_dict_sa_dict():
    """Test that a non-dict sa_dict is rejected with a warning, not raised."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_malformed_top')
    selection = select_from_sa_dict(sa_dict='not a dict', network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.qualified is False
    assert selection.selected_ts == []
    assert any('malformed' in warning.lower() for warning in selection.warnings)


def test_select_from_sa_dict_tolerates_non_dict_resolved_entry():
    """Test that a resolved SA entry that is not a dict-of-conditions is rejected with a warning."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_malformed_entry')
    sa_dict = {'A + B <=> C': 'not a dict', 'structures': dict()}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.qualified is False
    assert selection.selected_ts == []
    assert any('malformed' in warning.lower() for warning in selection.warnings)


def test_select_from_sa_dict_tolerates_non_dict_condition_entries():
    """Test that a condition mapping to something other than a dict is skipped with a warning,
    while other, well-formed conditions are still evaluated."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_malformed_condition')
    sa_dict = {'A + B <=> C': {
        (300.0, 'K', 1.0, 'bar'): 'not a dict',
        (400.0, 'K', 1.0, 'bar'): {'(TS) TS1': 1.0},
    }}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert any('malformed' in warning.lower() for warning in selection.warnings)
    assert len(selection.selected_ts) == 1


def test_select_from_sa_dict_tolerates_non_string_entry_keys():
    """Test that entry keys that are not strings (e.g. a stray non-string key) are skipped, not
    raised on, when checking the TS-entry prefix."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_malformed_keys')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {42: 1.0, '(TS) TS1': 1.0}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert len(selection.selected_ts) == 1
    assert selection.selected_ts[0].ts_label == 'TS1'


def test_resolve_direction_key_tolerates_non_dict_sa_dict():
    """Test that resolve_direction_key tolerates a non-dict sa_dict, returning (None, False, [warning])."""
    key, ambiguous, warnings = resolve_direction_key(sa_dict='not a dict', network_reaction='A + B <=> C')
    assert key is None
    assert ambiguous is False
    assert warnings


def test_resolve_direction_key_tolerates_non_string_keys():
    """Test that resolve_direction_key ignores non-string keys rather than raising on them."""
    sa_dict = {42: dict(), 'A + B <=> C': dict(), 'structures': dict()}
    key, ambiguous, warnings = resolve_direction_key(sa_dict=sa_dict, network_reaction='A + B <=> C')
    assert key == 'A + B <=> C'
    assert ambiguous is False


# --- FIX6: delta_ln_k field ----------------------------------------------------------------------

def test_delta_ln_k_equals_abs_coefficient_times_perturbation():
    """Test that every selected TS evidence record carries delta_ln_k == abs(coefficient) * perturbation."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_delta_ln_k')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': -0.03}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001, perturbation=8368.0)
    assert len(selection.selected_ts) == 1
    entry = selection.selected_ts[0]
    assert entry.delta_ln_k == pytest.approx(abs(-0.03) * 8368.0)


# --- FIX7: as_dict() JSON/YAML-safe rendering -----------------------------------------------------

def test_selection_as_dict_is_json_serializable(sa_dict, network):
    """Test that select_from_sa_dict()'s result renders to a JSON-serializable plain dict, with
    conditions rendered as T/T_unit/P/P_unit dicts rather than tuples."""
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001)
    rendered = selection.as_dict()
    serialized = json.dumps(rendered)
    assert isinstance(serialized, str)
    assert rendered['selected_ts']
    for entry in rendered['selected_ts']:
        assert isinstance(entry, dict)
        assert set(entry['condition'].keys()) == {'T', 'T_unit', 'P', 'P_unit'}
        assert entry['delta_ln_k'] == pytest.approx(abs(entry['coefficient']) * selection.thresholds['perturbation'])


def test_sensitive_transition_state_as_dict_shape():
    """Test SensitiveTransitionState.as_dict()'s exact shape directly."""
    entry = SensitiveTransitionState(ts_label='TS1', coefficient=0.1, condition=(300.0, 'K', 1.0, 'bar'),
                                     path_reaction_label='reactionA', path_reaction_str='S1 + S2 <=> S3',
                                     kinetics_comment='comment', uncertain=True, delta_ln_k=0.1 * 8368.0)
    rendered = entry.as_dict()
    assert rendered == {
        'ts_label': 'TS1',
        'coefficient': 0.1,
        'condition': {'T': 300.0, 'T_unit': 'K', 'P': 1.0, 'P_unit': 'bar'},
        'path_reaction_label': 'reactionA',
        'path_reaction_str': 'S1 + S2 <=> S3',
        'kinetics_comment': 'comment',
        'uncertain': True,
        'delta_ln_k': 0.1 * 8368.0,
    }
    json.dumps(rendered)


# --- FIX8: reworded reverse-direction warning, pinned against the real fixture ------------------------

def test_resolve_direction_key_reverse_warning_wording():
    """Test the reworded reverse-direction warning: it must state the requested direction is
    absent, name the opposite-direction entry being used instead, and assert coefficient equality
    to solver tolerance between directions -- not imply a different derivative question is asked."""
    sa_dict = {'C <=> A + B': dict(), 'structures': dict()}
    key, ambiguous, warnings = resolve_direction_key(sa_dict=sa_dict, network_reaction='A + B <=> C')
    assert key == 'C <=> A + B'
    assert ambiguous is True
    assert len(warnings) == 1
    warning = warnings[0].lower()
    assert 'no entry for the requested direction' in warning
    assert 'opposite-direction entry' in warning
    assert 'solver tolerance' in warning


def test_forward_reverse_coefficients_agree_on_real_fixture(sa_dict):
    """Premise guard (FIX8): pin the empirical claim the reworded warning now makes -- that
    Arkane's sensitivity coefficients for the two directions of a reaction agree to within solver
    tolerance -- against the real network4_2 fixture. For every reaction present in the SA output
    under both directions, every TS/well entry at or above the absolute floor on either side must
    agree between directions to within 1e-6 relative. If this premise were ever false, FIX8's
    reworded warning would be asserting something untrue about the underlying Arkane data.
    """
    keys = [key for key in sa_dict.keys() if key != 'structures']

    def reverse_of(reaction_str):
        lhs, rhs = reaction_str.split('<=>')
        return f'{rhs.strip()} <=> {lhs.strip()}'

    floor = coefficient_floor(min_delta_ln_k=1e-3, perturbation=8368.0)
    checked_pairs = 0
    for key in keys:
        reverse_key = reverse_of(key)
        if reverse_key not in keys:
            continue
        forward_entries_by_condition = sa_dict[key]
        reverse_entries_by_condition = sa_dict[reverse_key]
        for condition, forward_entries in forward_entries_by_condition.items():
            reverse_entries = reverse_entries_by_condition.get(condition)
            if not reverse_entries:
                continue
            for entry, forward_value in forward_entries.items():
                if entry not in reverse_entries:
                    continue
                reverse_value = reverse_entries[entry]
                if abs(forward_value) < floor and abs(reverse_value) < floor:
                    continue
                checked_pairs += 1
                relative_difference = abs(forward_value - reverse_value) / max(abs(forward_value),
                                                                                abs(reverse_value))
                assert relative_difference < 1e-6, (
                    f'{key} vs {reverse_key}, entry {entry}, condition {condition}: '
                    f'{forward_value} vs {reverse_value} (relative difference {relative_difference})')
    assert checked_pairs > 0, 'No forward/reverse pair was found in the fixture to check the premise against.'


def test_reason_does_not_report_a_verdict_on_a_decision_that_was_never_made():
    """`reason()` used to branch on `qualified` alone. `qualified` defaults to False, so a
    `not_evaluated` decision -- one where criterion (b) could not be computed at all -- rendered as
    a confident "does not qualify", and that string is what lands in the provenance record a human
    later reads. A decision that was never made must not read as a negative one."""
    selection = PDepNetworkSelection(network_id='network1',
                                     qualified=False,
                                     evaluation_status=EVALUATION_STATUS_NOT_EVALUATED,
                                     warnings=['The SA output could not be read.'])
    reason = selection.reason()
    assert 'does not qualify' not in reason, f'A never-made decision was reported as a negative: {reason}'
    assert 'NOT evaluated' in reason, f'Reason did not say the criterion was not evaluated: {reason}'
    assert 'The SA output could not be read.' in reason, \
        f'Reason did not surface the recorded diagnosis: {reason}'


def test_reason_still_reports_a_real_negative_as_a_negative():
    """The over-refusal counterpart: an EVALUATED negative is a real, computed finding and must keep
    reading as one. Fixing the not_evaluated case must not make every negative unspeakable."""
    selection = PDepNetworkSelection(network_id='network1',
                                     qualified=False,
                                     evaluation_status=EVALUATION_STATUS_EVALUATED)
    reason = selection.reason()
    assert 'does not qualify' in reason, f'A real evaluated negative stopped reading as a negative: {reason}'
    assert 'NOT evaluated' not in reason, f'A real evaluated negative was reported as unevaluated: {reason}'


def test_reason_reports_a_partial_yes_as_a_yes():
    """`combine()` unions `qualified` over EVALUATED components only, but sets `evaluation_status` to
    not_evaluated if ANY component was not evaluated -- so `qualified=True, not_evaluated` is a
    reachable state, not a defensive one. The asymmetry says a partial yes is supported: it must
    still read as a qualification, with the incompleteness noted rather than the verdict withdrawn."""
    selection = PDepNetworkSelection(
        network_id='network1',
        qualified=True,
        evaluation_status=EVALUATION_STATUS_NOT_EVALUATED,
        uncertain_path_reactions=[SensitiveTransitionState(ts_label='TS1',
                                                          coefficient=1.0,
                                                          condition=(300.0, 'K', 1.0, 'bar'),
                                                          path_reaction_label='r1',
                                                          path_reaction_str='A <=> B',
                                                          kinetics_comment='Estimated',
                                                          uncertain=True,
                                                          delta_ln_k=1.0)],
    )
    reason = selection.reason()
    assert 'qualifies for QM refinement' in reason, f'A partial yes stopped reading as a yes: {reason}'
    assert 'TS1' in reason, f'The evidence for the yes was not named: {reason}'
    assert 'could not be evaluated' in reason, f'The incompleteness was not disclosed: {reason}'


def test_reason_on_a_not_evaluated_decision_with_no_recorded_diagnosis():
    """Every not_evaluated path in `select_from_sa_dict` appends a warning, but `evaluation_status`
    is a plain mutable field with no invariant tying it to `warnings`, and a persisted record is
    validated for type, not for that pairing. Silence must read as "no diagnosis", not as a verdict."""
    selection = PDepNetworkSelection(network_id='network1',
                                     qualified=False,
                                     evaluation_status=EVALUATION_STATUS_NOT_EVALUATED)
    reason = selection.reason()
    assert 'does not qualify' not in reason, f'A never-made decision was reported as a negative: {reason}'
    assert 'No diagnosis was recorded' in reason, f'Missing diagnosis was not disclosed: {reason}'


def test_reason_on_a_qualified_record_that_names_no_evidence():
    """A record loaded off disk is validated for field TYPE, not for the cross-field invariant that
    `qualified=True` implies a non-empty `uncertain_path_reactions`. Without a guard this rendered
    "sensitive to , whose kinetics are estimates" -- a sentence with a hole where the evidence goes."""
    selection = PDepNetworkSelection(network_id='network1',
                                     qualified=True,
                                     evaluation_status=EVALUATION_STATUS_EVALUATED)
    reason = selection.reason()
    assert 'sensitive to ,' not in reason, f'Rendered an empty evidence list as prose: {reason}'
    assert 'names no uncertain transition state' in reason, \
        f'Did not disclose that the evidence is missing: {reason}'
