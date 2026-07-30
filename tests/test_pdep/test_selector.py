#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_selector
"""

import dataclasses
import json
import math
import os

import pytest
from arc.common import read_yaml_file

from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.parser import parse_pdep_network_file, parse_pdep_network_text
from t3.pdep.selector import (EVALUATION_STATUS_EVALUATED,
                              EVALUATION_STATUS_NOT_EVALUATED,
                              PDepNetworkSelection,
                              SensitiveTransitionState,
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


def test_coefficient_floor_raises_for_non_positive_min_delta_ln_k():
    """Test that ``coefficient_floor`` raises ``ValueError`` for a non-positive ``min_delta_ln_k``,
    naming that parameter (not ``perturbation``) in the message -- a zero or negative
    ``min_delta_ln_k`` would silently produce a floor of 0.0 (or a negative floor), which is a
    fail-open, not a deliberate 'accept anything' choice."""
    with pytest.raises(ValueError, match='min_delta_ln_k'):
        coefficient_floor(min_delta_ln_k=0, perturbation=8368.0)
    with pytest.raises(ValueError, match='min_delta_ln_k'):
        coefficient_floor(min_delta_ln_k=-1e-3, perturbation=8368.0)


def test_coefficient_floor_raises_for_non_finite_values():
    """Test that ``coefficient_floor`` raises ``ValueError`` for a NaN or infinite
    ``min_delta_ln_k`` or ``perturbation`` -- a non-finite input would silently produce a NaN or
    inf floor (nothing, or everything, qualifies) rather than a real gate."""
    with pytest.raises(ValueError, match='min_delta_ln_k'):
        coefficient_floor(min_delta_ln_k=math.nan, perturbation=8368.0)
    with pytest.raises(ValueError, match='min_delta_ln_k'):
        coefficient_floor(min_delta_ln_k=math.inf, perturbation=8368.0)
    with pytest.raises(ValueError, match='perturbation'):
        coefficient_floor(min_delta_ln_k=1e-3, perturbation=math.nan)
    with pytest.raises(ValueError, match='perturbation'):
        coefficient_floor(min_delta_ln_k=1e-3, perturbation=math.inf)


def test_coefficient_floor_raises_when_derived_floor_overflows():
    """Test that ``coefficient_floor`` raises ``ValueError`` when ``min_delta_ln_k`` and
    ``perturbation`` are each individually finite and positive, but the DERIVED floor
    ``min_delta_ln_k / perturbation`` overflows to infinity -- validating only the two inputs is
    not enough, since finite inputs can still divide to a non-finite result that would silently
    make the absolute gate reject nothing (an inf floor) or accept everything (a 0.0 floor from
    the symmetric underflow case)."""
    with pytest.raises(ValueError):
        coefficient_floor(min_delta_ln_k=1e300, perturbation=1e-300)


def test_coefficient_floor_raises_when_derived_floor_underflows_to_zero():
    """Test that ``coefficient_floor`` raises ``ValueError`` when the derived floor underflows to
    exactly 0.0 despite finite, positive inputs -- a 0.0 floor silently disables the absolute gate
    (every finite coefficient, including a denormal structural zero, would clear it)."""
    with pytest.raises(ValueError):
        coefficient_floor(min_delta_ln_k=1e-300, perturbation=1e300)


def test_select_from_sa_dict_raises_when_relative_cutoff_overflows():
    """Test that ``select_from_sa_dict`` raises ``ValueError`` rather than silently computing a
    non-finite cutoff when ``relative_threshold * max_abs_ts`` overflows -- an inf cutoff would
    fail every row closed (a wrong, unreported negative) instead of surfacing the degenerate
    threshold that produced it."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_cutoff_overflow')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 1.0e250}}}
    with pytest.raises(ValueError):
        select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                            relative_threshold=1.0e250, min_delta_ln_k=1e-3, perturbation=8368.0)


def test_select_from_sa_dict_raises_for_bad_relative_threshold():
    """Test that ``select_from_sa_dict`` raises ``ValueError`` for a NaN or negative
    ``relative_threshold`` before doing any work -- a NaN threshold makes every relative-gate
    comparison silently False, which would fall through to whatever the absolute floor decides
    rather than raising on the bad input."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_bad_rel')
    sa_dict = {'S1 + S2 <=> S3': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 1.0}}}
    with pytest.raises(ValueError, match='relative_threshold'):
        select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='S1 + S2 <=> S3',
                            relative_threshold=math.nan)
    with pytest.raises(ValueError, match='relative_threshold'):
        select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='S1 + S2 <=> S3',
                            relative_threshold=-0.001)


def test_select_sensitive_wells_raises_for_bad_relative_threshold():
    """Test that ``select_sensitive_wells`` raises ``ValueError`` for a NaN or negative
    ``relative_threshold``, for the same reason as ``select_from_sa_dict`` above."""
    entries_by_condition = {(300.0, 'K', 1.0, 'bar'): {'well1': 1.0}}
    with pytest.raises(ValueError, match='relative_threshold'):
        select_sensitive_wells(entries_by_condition=entries_by_condition, relative_threshold=math.nan)
    with pytest.raises(ValueError, match='relative_threshold'):
        select_sensitive_wells(entries_by_condition=entries_by_condition, relative_threshold=-0.001)


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


def test_select_sensitive_wells_tolerates_non_string_entry_keys():
    """Test that ``select_sensitive_wells`` tolerates a non-string entry key (e.g. a corrupted SA
    dict with a numeric or tuple key) by treating it as a non-TS row that is simply excluded, the
    same tolerance :func:`select_from_sa_dict` already applies on its TS-selection path -- without
    this guard, ``entry.startswith(TS_ENTRY_PREFIX)`` raises ``AttributeError`` on a non-string key
    instead of producing a decision."""
    entries_by_condition = {
        (300.0, 'K', 1.0, 'bar'): {
            'well_negative': -1.0,
            42: -5.0,
        },
    }
    sensitive_wells = select_sensitive_wells(entries_by_condition=entries_by_condition,
                                            relative_threshold=0.5,
                                            min_delta_ln_k=1e-3,
                                            perturbation=8368.0)
    assert 'well_negative' in sensitive_wells
    assert 42 not in sensitive_wells


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


# --- FIX L: SensitiveTransitionState.as_dict() falls back for a non-4-tuple condition -----------

def test_sensitive_transition_state_as_dict_falls_back_for_non_4_tuple_condition():
    """Test that as_dict() renders a malformed (non-4-element) condition as a plain list rather
    than crashing on the T/T_unit/P/P_unit unpack."""
    entry = SensitiveTransitionState(ts_label='TS1', coefficient=0.1, condition=(300.0, 'K'),
                                     path_reaction_label='reactionA', path_reaction_str='S1 + S2 <=> S3',
                                     kinetics_comment='comment', uncertain=True, delta_ln_k=0.1 * 8368.0)
    rendered = entry.as_dict()
    assert rendered['condition'] == [300.0, 'K']
    json.dumps(rendered)


# --- FIX C: evaluation_status --------------------------------------------------------------------

def test_evaluation_status_evaluated_for_a_normal_decision():
    """Test that a normal, fully-computed decision carries evaluation_status == 'evaluated'."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_eval_status_ok')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 1.0}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.evaluation_status == EVALUATION_STATUS_EVALUATED


def test_evaluation_status_not_evaluated_for_malformed_sa_dict():
    """Test that a non-dict sa_dict sets evaluation_status to 'not_evaluated'."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_eval_status_top')
    selection = select_from_sa_dict(sa_dict='not a dict', network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED


def test_evaluation_status_not_evaluated_when_direction_key_unresolved():
    """Test that a network reaction absent from the SA dict sets evaluation_status to 'not_evaluated',
    since no decision could actually be computed."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_eval_status_absent')
    sa_dict = {'X + Y <=> Z': dict(), 'structures': dict()}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED


def test_evaluation_status_not_evaluated_for_malformed_resolved_entry():
    """Test that a resolved SA entry that is not a dict-of-conditions sets evaluation_status to
    'not_evaluated'."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_eval_status_entry')
    sa_dict = {'A + B <=> C': 'not a dict', 'structures': dict()}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED


# --- FIX I: network_reactions_examined for single-reaction decisions -----------------------------

def test_network_reactions_examined_is_1_when_direction_key_resolves():
    """Test that network_reactions_examined is 1 once a direction key was resolved."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_examined_ok')
    sa_dict = {'A + B <=> C': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 1.0}}}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.network_reactions_examined == 1


def test_network_reactions_examined_is_0_when_direction_key_unresolved():
    """Test that network_reactions_examined is 0 when the requested reaction could not be located
    in the SA output at all."""
    network = parse_pdep_network_text(text=SYNTHETIC_NETWORK_TEXT, network_id='synthetic_examined_absent')
    sa_dict = {'X + Y <=> Z': dict(), 'structures': dict()}
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction='A + B <=> C',
                                    relative_threshold=0.001)
    assert selection.network_reactions_examined == 0


# --- FIX H: combine() -----------------------------------------------------------------------------

def _make_decision(network_id='net1', qualified=False, method='MSC', cache_status='generated',
                   direction_key='A + B <=> C', selected_ts=None, uncertain_path_reactions=None,
                   warnings=None):
    """Build a minimal PDepNetworkSelection for combine() tests."""
    return PDepNetworkSelection(
        network_id=network_id,
        qualified=qualified,
        network_reaction='A + B <=> C',
        direction_key=direction_key,
        method=method,
        sa_path='sa.yml',
        cache_status=cache_status,
        thresholds={'relative_threshold': 0.001},
        selected_ts=list(selected_ts) if selected_ts is not None else list(),
        uncertain_path_reactions=list(uncertain_path_reactions) if uncertain_path_reactions is not None else list(),
        warnings=list(warnings) if warnings is not None else list(),
    )


def test_combine_raises_value_error_for_empty_list():
    """Test that combine() raises ValueError on an empty list of decisions."""
    with pytest.raises(ValueError):
        PDepNetworkSelection.combine([])


def test_combine_raises_value_error_for_network_id_mismatch():
    """Test that combine() raises ValueError when decisions disagree on network_id, since silently
    combining decisions for different networks would fabricate confidence."""
    decisions = [_make_decision(network_id='net1'), _make_decision(network_id='net2')]
    with pytest.raises(ValueError):
        PDepNetworkSelection.combine(decisions)


def test_combine_single_component_is_equivalent():
    """Test that combining a single decision is equivalent to that decision, modulo network_reaction
    becoming None and direction_keys being populated."""
    decision = _make_decision(qualified=True, direction_key='A + B <=> C')
    combined = PDepNetworkSelection.combine([decision])
    assert combined.network_id == decision.network_id
    assert combined.qualified == decision.qualified
    assert combined.method == decision.method
    assert combined.cache_status == decision.cache_status
    assert combined.network_reaction is None
    assert combined.direction_keys == ['A + B <=> C']
    assert combined.network_reactions_examined == 1


def test_combine_method_and_cache_status_disagreement_sets_none_and_warns():
    """Test that combine() sets method/cache_status to None and records a warning, rather than
    silently adopting the first value, when components disagree."""
    decisions = [_make_decision(method='MSC', cache_status='generated'),
                 _make_decision(method='RS', cache_status='cached_valid')]
    combined = PDepNetworkSelection.combine(decisions)
    assert combined.method is None
    assert combined.cache_status is None
    assert any('disagree on method' in warning for warning in combined.warnings)
    assert any('disagree on cache_status' in warning for warning in combined.warnings)


def test_combine_records_direction_keys_and_unions_qualified():
    """Test that combine() records per-component direction_keys in input order, and that qualified
    is True iff any input qualified."""
    decisions = [_make_decision(qualified=False, direction_key='A + B <=> C'),
                 _make_decision(qualified=True, direction_key='D + E <=> F')]
    combined = PDepNetworkSelection.combine(decisions)
    assert combined.direction_key is None
    assert combined.direction_keys == ['A + B <=> C', 'D + E <=> F']
    assert combined.qualified is True


# --- ROUND-34 P2: as_dict() and combine() enumerate the fields BY HAND ------------------------------
#
# Both ``PDepNetworkSelection.as_dict()`` and ``PDepNetworkSelection.combine()`` list the dataclass's
# fields literally, so a field added to the dataclass is silently absent from saved YAML (as_dict)
# and silently reset to its default on an aggregate decision (combine). Neither failure is visible
# from any test that asserts on the fields it already knows about -- the whole point is the field
# nobody thought to assert on. These two tests close that by deriving the expected field set from
# ``dataclasses.fields()``, so adding a field to the dataclass FAILS them until it is handled.

FULL_DECISION_SOURCE_HASH = 'sha256:' + '0' * 64


def _make_full_decision(direction_key='A + B <=> C', direction_ambiguous=False, evaluation_status=None,
                        selected_ts=None, uncertain_path_reactions=None, warnings=None,
                        network_source_hash=FULL_DECISION_SOURCE_HASH):
    """Build a PDepNetworkSelection with every field set to a non-default value, for the field-coverage
    tests below."""
    decision = _make_decision(
        network_id='net1',
        qualified=False,
        method='MSC',
        cache_status='generated',
        direction_key=direction_key,
        selected_ts=selected_ts,
        uncertain_path_reactions=uncertain_path_reactions,
        warnings=warnings,
    )
    decision.direction_ambiguous = direction_ambiguous
    decision.network_reactions_examined = 1
    decision.network_source_hash = network_source_hash
    if evaluation_status is not None:
        decision.evaluation_status = evaluation_status
    return decision


def _make_ts_entry(ts_label):
    """Build a SensitiveTransitionState for the field-coverage tests below."""
    return SensitiveTransitionState(
        ts_label=ts_label,
        coefficient=1e-5,
        condition=(1000.0, 'K', 1.0, 'bar'),
        path_reaction_label='reaction1',
        path_reaction_str=f'A + B <=> C via {ts_label}',
        kinetics_comment='Estimated using template',
        uncertain=True,
        delta_ln_k=0.5,
    )


def test_as_dict_renders_every_dataclass_field():
    """Test that as_dict() renders EVERY PDepNetworkSelection field, so a newly added field cannot
    silently vanish from the saved YAML. as_dict() enumerates its keys by hand, so nothing but this
    equality stops a field from being dropped."""
    declared = {dataclass_field.name for dataclass_field in dataclasses.fields(PDepNetworkSelection)}
    rendered = set(PDepNetworkSelection(network_id='net1').as_dict().keys())
    assert rendered == declared, (
        f'as_dict() does not render every field. Missing: {sorted(declared - rendered)}; '
        f'unexpected: {sorted(rendered - declared)}.')


def test_as_dict_containers_are_isolated_from_the_live_selection():
    """Test that mutating any container in as_dict()'s returned dict cannot rewrite the
    PDepNetworkSelection's own state. as_dict() is the on-disk/reported shape, so a caller
    mutating what it returns (directly, or by mutating something nested inside it) must never be
    able to corrupt the live decision object."""
    ts_entry = SensitiveTransitionState(
        ts_label='TS1', coefficient=0.1, condition=(300.0, 'K', 1.0, 'bar'),
        path_reaction_label='R1', path_reaction_str='A + B <=> C',
        kinetics_comment='estimate', uncertain=True, delta_ln_k=0.1 * 8368.0)
    decision = PDepNetworkSelection(
        network_id='net1',
        qualified=True,
        thresholds={'relative_threshold': 0.001, 'nested': [1, 2, 3]},
        direction_keys=['A + B <=> C'],
        selected_ts=[ts_entry],
        uncertain_path_reactions=[ts_entry],
        warnings=['a warning'],
    )
    rendered = decision.as_dict()

    # Mutate every container in the rendered dict, including nested containers.
    rendered['thresholds']['nested'].append(999)
    rendered['thresholds']['relative_threshold'] = -1.0
    rendered['thresholds']['injected'] = True
    rendered['direction_keys'].append('injected')
    rendered['warnings'].append('injected')
    rendered['selected_ts'].append({'injected': True})
    rendered['selected_ts'][0]['ts_label'] = 'tampered'
    rendered['uncertain_path_reactions'].append({'injected': True})
    rendered['uncertain_path_reactions'][0]['ts_label'] = 'tampered'

    assert decision.thresholds == {'relative_threshold': 0.001, 'nested': [1, 2, 3]}, (
        f'thresholds mutated via as_dict() output: {decision.thresholds}')
    assert decision.direction_keys == ['A + B <=> C'], (
        f'direction_keys mutated via as_dict() output: {decision.direction_keys}')
    assert decision.warnings == ['a warning'], (
        f'warnings mutated via as_dict() output: {decision.warnings}')
    assert len(decision.selected_ts) == 1 and decision.selected_ts[0].ts_label == 'TS1', (
        f'selected_ts mutated via as_dict() output: {decision.selected_ts}')
    assert len(decision.uncertain_path_reactions) == 1 and decision.uncertain_path_reactions[0].ts_label == 'TS1', (
        f'uncertain_path_reactions mutated via as_dict() output: {decision.uncertain_path_reactions}')


# How ``combine()`` must treat each field, as a predicate on the COMBINED decision produced from the
# two components built in ``test_combine_handles_every_dataclass_field`` below. A field added to
# ``PDepNetworkSelection`` without an entry here fails ``test_combine_states_a_policy_for_every_field``,
# which is deliberate: the author has to state what an aggregate decision means for it, because the
# default -- silently taking the dataclass default -- is exactly the failure this guards.
_COMBINE_EXPECTED = {
    # Identity: taken from the first decision; all components must agree or combine() raises.
    'network_id': lambda value: value == 'net1',
    # Identity too: propagated when the components agree, refused when two real hashes differ,
    # None + warning when only some components recorded one.
    'network_source_hash': lambda value: value == FULL_DECISION_SOURCE_HASH,
    # Unioned across components: True iff any component qualified.
    'qualified': lambda value: value is False,
    # Deliberately RESET: the aggregate no longer refers to one network reaction.
    'network_reaction': lambda value: value is None,
    'direction_key': lambda value: value is None,
    # Recorded per component, in input order.
    'direction_keys': lambda value: value == ['A + B <=> C', 'D + E <=> F'],
    # Unioned: True iff any component had it True.
    'direction_ambiguous': lambda value: value is True,
    # Taken from the first decision (None + warning if components disagree).
    'method': lambda value: value == 'MSC',
    'sa_path': lambda value: value == 'sa.yml',
    'cache_status': lambda value: value == 'generated',
    'thresholds': lambda value: value == {'relative_threshold': 0.001},
    # Unioned, de-duplicated, first-seen order.
    'selected_ts': lambda value: [entry.ts_label for entry in value] == ['TS1', 'TS2'],
    'uncertain_path_reactions': lambda value: [entry.ts_label for entry in value] == ['TS1'],
    'warnings': lambda value: 'component warning' in value,
    # Set to the number of decisions combined.
    'network_reactions_examined': lambda value: value == 2,
    # Fail-closed: an aggregate 'does not qualify' is only backed if EVERY component was actually
    # evaluated. One never-evaluated component makes the aggregate not_evaluated.
    'evaluation_status': lambda value: value == EVALUATION_STATUS_NOT_EVALUATED,
}


def test_combine_states_a_policy_for_every_field():
    """Test that _COMBINE_EXPECTED covers every PDepNetworkSelection field. This is a conscience
    test: it does not check combine()'s behaviour, it forces whoever adds a field to state what an
    aggregate decision means for it, instead of inheriting the dataclass default by omission."""
    declared = {dataclass_field.name for dataclass_field in dataclasses.fields(PDepNetworkSelection)}
    assert set(_COMBINE_EXPECTED) == declared, (
        f'_COMBINE_EXPECTED does not cover every field. Missing: {sorted(declared - set(_COMBINE_EXPECTED))}; '
        f'stale: {sorted(set(_COMBINE_EXPECTED) - declared)}.')


def test_combine_handles_every_dataclass_field():
    """Test that combine() produces the documented value for EVERY field, from two components that
    set every field to a non-default value. A field combine() forgets falls back to the dataclass
    default and fails its predicate here."""
    first = _make_full_decision(
        direction_key='A + B <=> C',
        direction_ambiguous=False,
        selected_ts=[_make_ts_entry('TS1')],
        uncertain_path_reactions=[_make_ts_entry('TS1')],
        warnings=['component warning'],
    )
    second = _make_full_decision(
        direction_key='D + E <=> F',
        direction_ambiguous=True,
        evaluation_status=EVALUATION_STATUS_NOT_EVALUATED,
        selected_ts=[_make_ts_entry('TS2')],
    )
    combined = PDepNetworkSelection.combine([first, second])
    failures = [name for name, predicate in _COMBINE_EXPECTED.items()
                if not predicate(getattr(combined, name))]
    assert not failures, (
        f'combine() produced an unexpected value for: '
        f'{ {name: getattr(combined, name) for name in failures} }.')


def test_combine_reports_partial_coverage_even_when_the_aggregate_qualified():
    """Test that evaluation_status states a FACT about coverage, not a judgement about usability: an
    aggregate that qualified is still 'not_evaluated' when a component was not evaluated, because it
    was not. Whether such an aggregate may be acted on is the consumer's call -- see
    test_api.py::test_explore_pdep_network_accepts_a_qualified_but_partially_evaluated_selection."""
    qualifying = _make_full_decision(uncertain_path_reactions=[_make_ts_entry('TS1')])
    qualifying.qualified = True
    unevaluated = _make_full_decision(direction_key=None, evaluation_status=EVALUATION_STATUS_NOT_EVALUATED)
    combined = PDepNetworkSelection.combine([qualifying, unevaluated])
    assert combined.qualified is True
    assert combined.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED
    assert any('were not evaluated' in warning for warning in combined.warnings)


def test_combine_ignores_the_qualified_flag_of_an_unevaluated_component():
    """Test that a component carrying qualified=True alongside evaluation_status='not_evaluated' does
    not carry the aggregate. PDepNetworkSelection is mutable and enforces no invariant, so that state
    is constructible; counting its vote would let a flag its own status declares meaningless decide
    whether the network gets expensive QM."""
    unevaluated = _make_full_decision(evaluation_status=EVALUATION_STATUS_NOT_EVALUATED)
    unevaluated.qualified = True
    evaluated = _make_full_decision(direction_key=None)
    combined = PDepNetworkSelection.combine([unevaluated, evaluated])
    assert combined.qualified is False


def test_combine_warns_about_unevaluated_components_and_counts_them():
    """Test that combine() records HOW MANY components were not evaluated, so a partial evaluation is
    never invisible, whichever way evaluation_status lands."""
    decisions = [_make_full_decision(),
                 _make_full_decision(direction_key=None, evaluation_status=EVALUATION_STATUS_NOT_EVALUATED),
                 _make_full_decision(direction_key=None, evaluation_status=EVALUATION_STATUS_NOT_EVALUATED)]
    combined = PDepNetworkSelection.combine(decisions)
    assert any('2 of 3 combined decisions were not evaluated' in warning for warning in combined.warnings)


def test_combine_of_all_evaluated_components_records_no_unevaluated_warning():
    """Test that the not-evaluated warning is NOT emitted when every component was evaluated, so the
    warning stays a signal rather than noise on every aggregate."""
    combined = PDepNetworkSelection.combine([_make_full_decision(), _make_full_decision()])
    assert combined.evaluation_status == EVALUATION_STATUS_EVALUATED
    assert not any('were not evaluated' in warning for warning in combined.warnings)


# --- ROUND-34 P1: network identity bound to CONTENT, not the file stem ------------------------------

def test_combine_raises_for_two_different_source_hashes():
    """Test that combine() refuses two components computed from different revisions of the network.
    network_id alone cannot catch this: both components name the same file stem, and an aggregate
    over two revisions describes no network that exists."""
    decisions = [_make_full_decision(network_source_hash='sha256:' + 'a' * 64),
                 _make_full_decision(network_source_hash='sha256:' + 'b' * 64)]
    with pytest.raises(ValueError):
        PDepNetworkSelection.combine(decisions)


def test_combine_records_none_and_warns_when_only_some_components_have_a_source_hash():
    """Test that a hash missing on some components yields None plus a warning, rather than adopting
    the one real hash: missing means 'not recorded', not 'different', but binding the aggregate to
    bytes part of it was never computed from would be a fail-open."""
    decisions = [_make_full_decision(), _make_full_decision(network_source_hash=None)]
    combined = PDepNetworkSelection.combine(decisions)
    assert combined.network_source_hash is None
    assert any('network_source_hash' in warning for warning in combined.warnings)


def test_combine_of_components_without_a_source_hash_does_not_warn():
    """Test that components that all lack a source hash combine to None WITHOUT a warning, so the
    warning stays a signal about disagreement rather than noise about absence."""
    combined = PDepNetworkSelection.combine([_make_full_decision(network_source_hash=None),
                                             _make_full_decision(network_source_hash=None)])
    assert combined.network_source_hash is None
    assert not any('network_source_hash' in warning for warning in combined.warnings)


def test_select_from_sa_dict_records_the_parsed_network_source_hash(sa_dict, network):
    """Test that a decision carries the content hash of the network file it was computed from, so a
    consumer acting on it later can tell whether the file has changed since."""
    selection = select_from_sa_dict(sa_dict=sa_dict, network=network, network_reaction=TARGET_REACTION,
                                    relative_threshold=0.001)
    assert selection.network_source_hash == network.source_hash
    assert selection.network_source_hash.startswith('sha256:')
def test_as_dict_does_not_hand_out_the_conditions_nested_objects():
    """Test that mutating a nested object inside the ``condition`` rendered by
    ``SensitiveTransitionState.as_dict()`` cannot reach the live record.

    Same defect class as the manifest/selection/thresholds leaks already closed elsewhere on this
    branch: a frozen dataclass that hands out its own nested containers is frozen in name only.
    A condition reaching this record from a loaded selection is not guaranteed to be a flat tuple
    of scalars (``_sensitive_transition_state_from_dict`` coerces whatever it finds), so both the
    4-element branch and the fallback branch must copy rather than alias.
    """
    entry = SensitiveTransitionState(ts_label='TS1', coefficient=1.0, condition=([300.0], 'K', 1.0, 'bar'),
                                     path_reaction_label='r1', path_reaction_str='A <=> B',
                                     kinetics_comment='', uncertain=True, delta_ln_k=1.0)
    rendered = entry.as_dict()
    rendered['condition']['T'].append('tampered')
    assert entry.condition[0] == [300.0], f'The live condition was rewritten: {entry.condition!r}.'


def test_as_dict_does_not_hand_out_a_malformed_condition_itself():
    """Test that the malformed-condition fallback in ``as_dict()`` copies too.

    A condition that is neither a tuple nor a list was previously returned by reference, so the
    caller received the record's own object rather than a rendering of it.
    """
    entry = SensitiveTransitionState(ts_label='TS1', coefficient=1.0, condition={'T': [300.0]},
                                     path_reaction_label='r1', path_reaction_str='A <=> B',
                                     kinetics_comment='', uncertain=True, delta_ln_k=1.0)
    rendered = entry.as_dict()
    assert rendered['condition'] is not entry.condition, 'as_dict() handed out the live condition object.'
    rendered['condition']['T'].append('tampered')
    assert entry.condition == {'T': [300.0]}, f'The live condition was rewritten: {entry.condition!r}.'


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
