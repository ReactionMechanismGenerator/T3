#!/usr/bin/env python3
# encoding: utf-8

"""
Tests for the t3.pdep.assessment module.

The assessment record is the durable answer to "what did T3 decide about every PDep network it
looked at this iteration?", including the ones it could not decide about. Its whole value is that
it cannot quietly disagree with itself, so most of what is tested here is the cross-field
invariants rather than the field values.
"""

import dataclasses
import os

import pytest

from arc.common import read_yaml_file, save_yaml_file

from t3.pdep.assessment import (ASSESSMENT_RECORD_FILE_NAME,
                                ASSESSMENT_RECORD_SCHEMA_VERSION,
                                PDepNetworkAssessment,
                                assessments_record_path,
                                )
from t3.pdep.reason_codes import (INTERNAL_ERROR_REASON_CODES,
                                  PRE_SELECTOR_REASON_CODES,
                                  REASON_CODE_STATUS,
                                  REASON_EVALUATED_NO_UNCERTAIN_TS,
                                  REASON_INTERNAL_ERROR,
                                  REASON_NETWORK_PARSE_FAILED,
                                  REASON_QUALIFIED_UNCERTAIN_TS,
                                  REASON_SA_ALL_METHODS_FAILED,
                                  REASON_SA_OUTPUT_MALFORMED,
                                  REASON_SELECTOR_NO_TS_ROWS,
                                  REASON_SELECTOR_SA_PAYLOAD_MALFORMED,
                                  REASON_SELECTOR_TS_PROVENANCE_UNASSESSED,
                                  SELECTION_BEARING_REASON_CODES,
                                  STATUS_EVALUATED_NEGATIVE,
                                  STATUS_INTERNAL_ERROR,
                                  STATUS_NOT_EVALUATED,
                                  STATUS_QUALIFIED,
                                  VALID_ASSESSMENT_REASON_CODES,
                                  VALID_ASSESSMENT_STATUSES,
                                  )
from t3.pdep.selector import (CACHE_STATUS_CACHED_VALID,
                              PDepNetworkSelection,
                              SensitiveTransitionState,
                              )


def uncertain_ts(ts_label='TS1'):
    """One piece of the evidence a qualified selection rests on."""
    return SensitiveTransitionState(ts_label=ts_label, coefficient=1.0e-3, condition=(1000.0, 1.0),
                                    path_reaction_label='rxn1', path_reaction_str='A <=> B',
                                    kinetics_comment='Estimated using template', uncertain=True,
                                    delta_ln_k=0.5)


def qualified_selection(network_id='network4_2'):
    """A selection that qualifies, with the evidence a qualified one must carry."""
    selection = PDepNetworkSelection(network_id=network_id, qualified=True)
    selection.selected_ts = [uncertain_ts()]
    selection.uncertain_path_reactions = [uncertain_ts()]
    return selection


def negative_selection(network_id='network4_2'):
    """A selection that was evaluated and did not qualify."""
    return PDepNetworkSelection(network_id=network_id, qualified=False)


def selection_for(reason_code):
    """The nested selection a given reason code requires, or None if it forbids one."""
    if reason_code == REASON_QUALIFIED_UNCERTAIN_TS:
        return qualified_selection()
    if reason_code in SELECTION_BEARING_REASON_CODES:
        return negative_selection()
    return None


def assessment(**overrides):
    """A minimal valid assessment, overridable field by field."""
    kwargs = {'network_id': 'network4_2',
              'iteration': 1,
              'status': STATUS_NOT_EVALUATED,
              'reason_code': REASON_SA_ALL_METHODS_FAILED,
              }
    kwargs.update(overrides)
    return PDepNetworkAssessment(**kwargs)


def test_the_record_file_name_is_iteration_scoped_and_joined_onto_the_iteration_directory():
    assert assessments_record_path('/tmp/iteration_3') == os.path.join('/tmp/iteration_3',
                                                                       ASSESSMENT_RECORD_FILE_NAME)


def test_the_closed_reason_code_set_is_partitioned_into_three_disjoint_categories():
    # The categories are what decide whether a nested selection is required, forbidden, or optional,
    # so a code belonging to none of them -- or to two -- would escape that check entirely.
    assert len(VALID_ASSESSMENT_REASON_CODES) == len(set(VALID_ASSESSMENT_REASON_CODES))
    categories = [set(SELECTION_BEARING_REASON_CODES), set(PRE_SELECTOR_REASON_CODES),
                  set(INTERNAL_ERROR_REASON_CODES)]
    assert set.union(*categories) == set(VALID_ASSESSMENT_REASON_CODES)
    for i, first in enumerate(categories):
        for second in categories[i + 1:]:
            assert not first & second


def test_the_two_malformed_payload_codes_sit_in_the_categories_their_names_promise():
    # These two exist because ONE code could not serve both callers. T3's funnel inspects
    # `structures` before it can build a network reaction string, so a non-mapping SA payload is
    # caught there and never reaches the selector -- while the selector's public entry points can
    # be handed one directly. Putting either in the other's category would make the funnel unable
    # to record its own failure without either fabricating a selection or violating the type.
    assert REASON_SA_OUTPUT_MALFORMED in PRE_SELECTOR_REASON_CODES
    assert REASON_SA_OUTPUT_MALFORMED not in SELECTION_BEARING_REASON_CODES
    assert REASON_SELECTOR_SA_PAYLOAD_MALFORMED in SELECTION_BEARING_REASON_CODES
    assert REASON_SELECTOR_SA_PAYLOAD_MALFORMED not in PRE_SELECTOR_REASON_CODES


@pytest.mark.parametrize('reason_code', ['sa_all_methods_failed', 'sa_output_missing',
                                         'sa_output_unreadable', 'sa_structures_missing',
                                         'species_label_mapping_failed', 'network_parse_failed',
                                         'network_input_write_failed', 'network_discovery_failed'])
def test_the_funnels_own_failures_all_forbid_a_nested_selection(reason_code):
    # Every one of these names a point BEFORE the selector ran. If any drifted into the
    # selection-bearing half, the funnel would be forced to attach a selection it does not have.
    assert reason_code in PRE_SELECTOR_REASON_CODES
    with pytest.raises(ValueError, match='selection'):
        assessment(status=STATUS_NOT_EVALUATED, reason_code=reason_code,
                   selection=negative_selection())


@pytest.mark.parametrize('reason_code', [code for code in SELECTION_BEARING_REASON_CODES])
def test_every_selection_bearing_code_actually_requires_a_selection(reason_code):
    with pytest.raises(ValueError, match='selection'):
        assessment(status=REASON_CODE_STATUS[reason_code], reason_code=reason_code, selection=None)


def test_every_reason_code_has_a_status_and_every_status_is_reachable():
    assert set(REASON_CODE_STATUS) == set(VALID_ASSESSMENT_REASON_CODES)
    assert set(REASON_CODE_STATUS.values()) == set(VALID_ASSESSMENT_STATUSES)


def test_a_crash_is_not_counted_as_a_network_that_could_not_be_evaluated():
    # An internal error is a BUG, not an outcome. Bucketing it with the legitimately unevaluable
    # networks would make a T3 defect read as the ordinary cost of doing business.
    assert REASON_CODE_STATUS[REASON_INTERNAL_ERROR] == STATUS_INTERNAL_ERROR
    assert STATUS_INTERNAL_ERROR != STATUS_NOT_EVALUATED


@pytest.mark.parametrize('status', ['Qualified', 'unknown', '', None, 0])
def test_an_unrecognized_status_is_refused(status):
    with pytest.raises(ValueError, match='status'):
        assessment(status=status)


@pytest.mark.parametrize('reason_code', ['selector_made_up', '', None, 0])
def test_an_unrecognized_reason_code_is_refused(reason_code):
    with pytest.raises(ValueError, match='reason_code'):
        assessment(reason_code=reason_code)


# --- The cross-field invariant: a status may not disagree with its reason code -------------------

@pytest.mark.parametrize('reason_code', VALID_ASSESSMENT_REASON_CODES)
def test_every_reason_code_accepts_only_the_status_it_implies(reason_code):
    # Parametrized over the WHOLE closed set rather than a sample: a single mis-mapped code is
    # exactly the defect this mapping exists to prevent, and sampling would miss most of them.
    selection = selection_for(reason_code)
    implied = REASON_CODE_STATUS[reason_code]
    assessment(status=implied, reason_code=reason_code, selection=selection)
    for status in VALID_ASSESSMENT_STATUSES:
        if status == implied:
            continue
        with pytest.raises(ValueError, match='status'):
            assessment(status=status, reason_code=reason_code, selection=selection)


# --- The selection-presence invariant ------------------------------------------------------------

def test_a_selector_produced_reason_code_requires_the_selection_it_came_from():
    with pytest.raises(ValueError, match='selection'):
        assessment(status=STATUS_NOT_EVALUATED, reason_code=REASON_SELECTOR_NO_TS_ROWS,
                   selection=None)


def test_a_reason_code_from_before_the_selector_ran_must_not_carry_a_selection():
    # A network that never reached the selector cannot have produced one; a selection here would be
    # borrowed from somewhere else and would misattribute evidence.
    with pytest.raises(ValueError, match='selection'):
        assessment(status=STATUS_NOT_EVALUATED, reason_code=REASON_SA_ALL_METHODS_FAILED,
                   selection=negative_selection())


def test_an_internal_error_may_carry_a_selection_or_not():
    # A crash can happen before the selector ran or after it returned. Forbidding the selection
    # would discard the evidence already obtained -- the most useful part of the breadcrumb.
    assessment(status=STATUS_INTERNAL_ERROR, reason_code=REASON_INTERNAL_ERROR, selection=None)
    assessment(status=STATUS_INTERNAL_ERROR, reason_code=REASON_INTERNAL_ERROR,
               selection=qualified_selection())


def test_a_qualified_assessment_must_carry_a_selection_that_actually_qualified():
    with pytest.raises(ValueError, match='selection'):
        assessment(status=STATUS_QUALIFIED, reason_code=REASON_QUALIFIED_UNCERTAIN_TS,
                   selection=negative_selection())


def test_an_evaluated_negative_assessment_must_not_carry_a_qualified_selection():
    with pytest.raises(ValueError, match='selection'):
        assessment(status=STATUS_EVALUATED_NEGATIVE, reason_code=REASON_EVALUATED_NO_UNCERTAIN_TS,
                   selection=qualified_selection())


@pytest.mark.parametrize('impostor', [object(), 'a selection', 42, {'qualified': True},
                                      dict(), ['not', 'a', 'selection']])
def test_only_a_real_selection_can_be_nested(impostor):
    # Duck-typing on `.qualified` alone would let any object through construction and then fail at
    # serialization time -- after the record had already been treated as valid.
    with pytest.raises(ValueError, match='selection'):
        assessment(status=STATUS_NOT_EVALUATED, reason_code=REASON_SELECTOR_NO_TS_ROWS,
                   selection=impostor)


def test_the_nested_selection_is_snapshotted_so_a_later_mutation_cannot_rewrite_history():
    # Freezing the outer record does not freeze the selection handed to it; a caller still holding
    # the original could otherwise flip `qualified` and leave this validated record self-contradictory.
    selection = qualified_selection()
    record = assessment(status=STATUS_QUALIFIED, reason_code=REASON_QUALIFIED_UNCERTAIN_TS,
                        selection=selection)
    selection.qualified = False
    selection.network_id = 'mutated'
    assert record.selection.qualified is True
    assert record.selection.network_id == 'network4_2'
    assert record.as_dict()['selection']['qualified'] is True


def test_the_record_is_frozen_so_a_validated_field_cannot_be_reassigned():
    record = assessment()
    for field_name, value in [('status', STATUS_QUALIFIED), ('reason_code', REASON_INTERNAL_ERROR),
                              ('network_id', 'other'), ('selection', qualified_selection())]:
        with pytest.raises(dataclasses.FrozenInstanceError):
            setattr(record, field_name, value)


# --- Field validation ----------------------------------------------------------------------------

@pytest.mark.parametrize('iteration', [-1, 1.5, True, None, '1'])
def test_a_non_integral_or_negative_iteration_is_refused(iteration):
    with pytest.raises(ValueError, match='iteration'):
        assessment(iteration=iteration)


@pytest.mark.parametrize('schema_version', [ASSESSMENT_RECORD_SCHEMA_VERSION + 1, 0, True, None])
def test_a_schema_version_this_code_did_not_write_is_refused(schema_version):
    # `True == 1` in Python, so a bool would otherwise slip past an `!=` comparison against 1 and
    # be rendered as `true` in the YAML.
    with pytest.raises(ValueError, match='schema_version'):
        assessment(schema_version=schema_version)


@pytest.mark.parametrize('network_id', ['', None, 42])
def test_an_empty_or_non_string_network_id_is_refused(network_id):
    with pytest.raises(ValueError, match='network_id'):
        assessment(network_id=network_id)


@pytest.mark.parametrize('field_name', ['requested_me_methods', 'warnings'])
def test_a_bare_string_sequence_is_refused_rather_than_split_into_characters(field_name):
    # tuple('MSC') is ('M', 'S', 'C') -- a silent corruption that would record three nonsense ME
    # methods and read back as though T3 had genuinely tried them.
    with pytest.raises(ValueError, match=field_name):
        assessment(**{field_name: 'MSC'})


@pytest.mark.parametrize('field_name', ['requested_me_methods', 'warnings'])
def test_a_non_string_entry_in_a_string_sequence_is_refused(field_name):
    with pytest.raises(ValueError, match=field_name):
        assessment(**{field_name: ['CSE', 7]})


@pytest.mark.parametrize('sa_rank_index', [-1, 1.5, True, 'first'])
def test_a_non_integral_or_negative_sa_rank_index_is_refused(sa_rank_index):
    with pytest.raises(ValueError, match='sa_rank_index'):
        assessment(sa_rank_index=sa_rank_index)


@pytest.mark.parametrize('field_name', ['network_path', 'network_source_hash', 'observable_label',
                                        'chemkin_reaction', 'network_reaction', 'final_method',
                                        'sa_path', 'cache_status'])
def test_a_non_string_optional_field_is_refused(field_name):
    with pytest.raises(ValueError, match=field_name):
        assessment(**{field_name: 42})


# --- Secondary reason codes ----------------------------------------------------------------------

def test_secondary_reason_codes_record_further_defects_rather_than_dropping_them():
    record = assessment(status=STATUS_NOT_EVALUATED, reason_code=REASON_SELECTOR_NO_TS_ROWS,
                        selection=negative_selection(),
                        secondary_reason_codes=(REASON_SELECTOR_TS_PROVENANCE_UNASSESSED,))
    assert record.secondary_reason_codes == (REASON_SELECTOR_TS_PROVENANCE_UNASSESSED,)
    assert record.as_dict()['secondary_reason_codes'] == [REASON_SELECTOR_TS_PROVENANCE_UNASSESSED]


def test_a_secondary_code_that_repeats_the_primary_is_refused():
    with pytest.raises(ValueError, match='reason_code'):
        assessment(secondary_reason_codes=(REASON_SA_ALL_METHODS_FAILED,))


def test_an_unrecognized_or_repeated_secondary_code_is_refused():
    with pytest.raises(ValueError, match='reason_code'):
        assessment(secondary_reason_codes=('made_up',))
    with pytest.raises(ValueError, match='reason_code'):
        assessment(secondary_reason_codes=(REASON_NETWORK_PARSE_FAILED, REASON_NETWORK_PARSE_FAILED))


# --- Serialization -------------------------------------------------------------------------------

def test_as_dict_renders_every_field_and_nests_the_selection_rather_than_flattening_it():
    record = assessment(status=STATUS_QUALIFIED, reason_code=REASON_QUALIFIED_UNCERTAIN_TS,
                        selection=qualified_selection(), network_path='/runs/network4_2.py',
                        network_source_hash='sha256:abc', observable_label='OH',
                        sa_rank_index=2, chemkin_reaction='A(1) <=> B(2)',
                        network_reaction='A <=> B', requested_me_methods=('CSE', 'MSC'),
                        final_method='MSC', sa_path='/runs/sa.yml', cache_status='cached_valid',
                        warnings=('a warning',))
    content = record.as_dict()
    assert content == {'assessment_record_schema_version': ASSESSMENT_RECORD_SCHEMA_VERSION,
                       'network_id': 'network4_2',
                       'iteration': 1,
                       'status': STATUS_QUALIFIED,
                       'reason_code': REASON_QUALIFIED_UNCERTAIN_TS,
                       'secondary_reason_codes': [],
                       'network_path': '/runs/network4_2.py',
                       'network_source_hash': 'sha256:abc',
                       'observable_label': 'OH',
                       'sa_rank_index': 2,
                       'chemkin_reaction': 'A(1) <=> B(2)',
                       'network_reaction': 'A <=> B',
                       'requested_me_methods': ['CSE', 'MSC'],
                       'final_method': 'MSC',
                       'sa_path': '/runs/sa.yml',
                       'cache_status': 'cached_valid',
                       'warnings': ['a warning'],
                       'selection': record.selection.as_dict(),
                       }
    # Nested, not flattened: the selection self-describes via its own schema/algorithm version
    # fields, which only survive if it is rendered as a nested mapping.
    assert content['selection']['selection_schema_version'] is not None


def test_as_dict_renders_a_missing_selection_as_none_rather_than_omitting_the_key():
    # An absent key and a null value read the same to a careless consumer but not to a strict one;
    # the key is always present so "no selection" is stated rather than implied.
    content = assessment().as_dict()
    assert 'selection' in content
    assert content['selection'] is None


def test_every_as_dict_value_survives_a_yaml_round_trip(tmp_path):
    record = assessment(status=STATUS_EVALUATED_NEGATIVE,
                        reason_code=REASON_EVALUATED_NO_UNCERTAIN_TS,
                        selection=negative_selection(), requested_me_methods=('CSE',),
                        warnings=('w1', 'w2'))
    path = os.path.join(str(tmp_path), 'record.yml')
    save_yaml_file(path=path, content=record.as_dict())
    assert read_yaml_file(path) == record.as_dict()


# --- The nested selection cannot drift away from the record after construction -------------------
# Rounds 55-56 established that a record must not be able to contradict itself. Round 57 found the
# hole left in that: freezing the record froze the REFERENCE, not the mutable object behind it, so
# every invariant checked at construction could be undone immediately afterwards and the corrupted
# state was what got serialized.

def test_mutating_the_nested_selection_cannot_change_what_gets_persisted():
    record = assessment(status=STATUS_QUALIFIED, reason_code=REASON_QUALIFIED_UNCERTAIN_TS,
                        selection=qualified_selection())
    persisted_before = record.as_dict()['selection']
    record.selection.selected_ts.append(uncertain_ts(ts_label='TS_SMUGGLED'))
    with pytest.raises(ValueError, match='no longer matches'):
        record.as_dict()
    # And what was captured at construction is untouched by the mutation.
    assert [entry['ts_label'] for entry in persisted_before['selected_ts']] == ['TS1']


def test_flipping_qualified_on_the_nested_selection_is_caught_at_serialization():
    # The exact corruption the deep copy was believed to prevent: the copy stopped the CALLER's
    # object from reaching in, but the record handed out its own copy just as freely.
    record = assessment(status=STATUS_QUALIFIED, reason_code=REASON_QUALIFIED_UNCERTAIN_TS,
                        selection=qualified_selection())
    record.selection.qualified = False
    with pytest.raises(ValueError, match='no longer matches'):
        record.as_dict()


def test_repointing_the_nested_selection_at_another_network_is_caught_at_serialization():
    record = assessment(status=STATUS_QUALIFIED, reason_code=REASON_QUALIFIED_UNCERTAIN_TS,
                        selection=qualified_selection())
    record.selection.network_id = 'network12_1'
    with pytest.raises(ValueError, match='no longer matches'):
        record.as_dict()


def test_an_unmutated_record_serializes_the_same_selection_every_time():
    # The drift check must not fire on the ordinary path, or it converts a safety net into an
    # outage: as_dict() is called once per network per iteration by the funnel.
    record = assessment(status=STATUS_QUALIFIED, reason_code=REASON_QUALIFIED_UNCERTAIN_TS,
                        selection=qualified_selection())
    assert record.as_dict() == record.as_dict()


def test_the_persisted_selection_cannot_be_edited_through_the_rendering_it_returns():
    record = assessment(status=STATUS_QUALIFIED, reason_code=REASON_QUALIFIED_UNCERTAIN_TS,
                        selection=qualified_selection())
    record.as_dict()['selection']['qualified'] = 'tampered'
    assert record.as_dict()['selection']['qualified'] is True


# --- The nested selection must be about the SAME network as the record ---------------------------

def test_a_selection_for_another_network_is_refused():
    # t3/pdep/reason_codes.py's own module docstring says the category rules exist so that evidence
    # from another network or another pass is never attached to this one -- but nothing enforced the
    # network part of that claim until now.
    with pytest.raises(ValueError, match='is about network'):
        assessment(status=STATUS_QUALIFIED, reason_code=REASON_QUALIFIED_UNCERTAIN_TS,
                   network_id='network4_2', selection=qualified_selection(network_id='network12_1'))


@pytest.mark.parametrize('record_field, selection_field', [('network_source_hash',
                                                            'network_source_hash'),
                                                           ('network_reaction', 'network_reaction'),
                                                           ('sa_path', 'sa_path'),
                                                           ('cache_status', 'cache_status'),
                                                           ('final_method', 'method')])
def test_provenance_the_record_and_its_selection_both_carry_may_not_disagree(record_field,
                                                                             selection_field):
    # Every field duplicated across the two is a chance for them to contradict, and a record that
    # names one SA run while its evidence names another is worse than one that names neither.
    selection = negative_selection()
    setattr(selection, selection_field, 'from_the_selection')
    with pytest.raises(ValueError, match='disagree'):
        assessment(status=STATUS_EVALUATED_NEGATIVE, reason_code=REASON_EVALUATED_NO_UNCERTAIN_TS,
                   selection=selection, **{record_field: 'from_the_record'})


# The value is carried per field rather than reusing one literal: `cache_status` is an enum on the
# selection (`PDepNetworkSelection.validate` checks it against `VALID_CACHE_STATUSES`), so a
# placeholder string is no longer constructible there. What is under test is whether ONE side
# carrying a value is accepted, not what the value is.
@pytest.mark.parametrize('record_field, selection_field, value',
                         [('network_source_hash', 'network_source_hash', 'sha256:only_here'),
                          ('network_reaction', 'network_reaction', 'A <=> B'),
                          ('sa_path', 'sa_path', '/runs/only_here.yml'),
                          ('cache_status', 'cache_status', CACHE_STATUS_CACHED_VALID),
                          ('final_method', 'method', 'MSC')])
def test_provenance_only_one_of_the_two_carries_is_accepted(record_field, selection_field, value):
    # Absence is not disagreement. The funnel legitimately knows things the selector never recorded
    # and vice versa, so requiring both to be set would refuse ordinary records.
    selection = negative_selection()
    setattr(selection, selection_field, value)
    assessment(status=STATUS_EVALUATED_NEGATIVE, reason_code=REASON_EVALUATED_NO_UNCERTAIN_TS,
               selection=selection)
    other = negative_selection()
    assessment(status=STATUS_EVALUATED_NEGATIVE, reason_code=REASON_EVALUATED_NO_UNCERTAIN_TS,
               selection=other, **{record_field: value})


def test_provenance_that_agrees_is_accepted():
    selection = negative_selection()
    selection.sa_path = '/runs/sa.yml'
    selection.method = 'MSC'
    assessment(status=STATUS_EVALUATED_NEGATIVE, reason_code=REASON_EVALUATED_NO_UNCERTAIN_TS,
               selection=selection, sa_path='/runs/sa.yml', final_method='MSC')


def test_an_internal_error_may_not_nest_another_networks_evidence():
    # An internal error is the one category where the two halves of a record are allowed to disagree
    # about the WORK -- a crash part-way through populating one legitimately leaves them mid-update.
    # It is not licence to disagree about the SUBJECT: "these fields are mid-update" is a coherent
    # thing for a crash record to say, "this is evidence about a different network" is not.
    with pytest.raises(ValueError, match='does not make another network'):
        PDepNetworkAssessment(network_id='network4_2',
                              status=STATUS_INTERNAL_ERROR,
                              reason_code=REASON_INTERNAL_ERROR,
                              selection=PDepNetworkSelection(network_id='network9_1'))


def test_an_internal_error_may_still_disagree_about_the_work_it_was_interrupted_during():
    # The other side of the same line: this must NOT be refused, or the breadcrumb would be discarded
    # over exactly the inconsistency it exists to report.
    record = PDepNetworkAssessment(network_id='network4_2',
                                   status=STATUS_INTERNAL_ERROR,
                                   reason_code=REASON_INTERNAL_ERROR,
                                   sa_path='/runs/half/written.yml',
                                   final_method='CSE',
                                   selection=PDepNetworkSelection(network_id='network4_2',
                                                                  sa_path='/runs/a/different.yml',
                                                                  method='MSC'))
    assert record.as_dict()['sa_path'] == '/runs/half/written.yml'
