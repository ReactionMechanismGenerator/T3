#!/usr/bin/env python3
# encoding: utf-8

"""
Tests for saving and loading the PDep network assessment record
(``t3.pdep.api.save_pdep_network_assessments`` / ``load_pdep_network_assessments``).

The record exists to be believed later: it is the only durable answer to "what did T3 decide about
every network it looked at, including the ones it could not decide about?". So the tests here are
mostly about the two ways a durable file betrays that trust -- being written half-way, and being
read back as something other than what was written.
"""

import os
from pathlib import Path

import pytest

from arc.common import save_yaml_file

from t3.pdep.api import (load_pdep_network_assessments,
                         save_pdep_network_assessments,
                         save_pdep_network_selections,
                         )
from t3.pdep.assessment import (ASSESSMENT_ENVELOPE_SCHEMA_VERSION,
                                ASSESSMENT_RECORD_SCHEMA_VERSION,
                                PDepNetworkAssessment,
                                )
from t3.pdep.reason_codes import (REASON_EVALUATED_NO_UNCERTAIN_TS,
                                  REASON_INTERNAL_ERROR,
                                  REASON_NETWORK_PARSE_FAILED,
                                  REASON_QUALIFIED_UNCERTAIN_TS,
                                  REASON_SA_ALL_METHODS_FAILED,
                                  REASON_SELECTOR_NO_TS_ROWS,
                                  STATUS_EVALUATED_NEGATIVE,
                                  STATUS_INTERNAL_ERROR,
                                  STATUS_NOT_EVALUATED,
                                  STATUS_QUALIFIED,
                                  )
from t3.pdep.selector import PDepNetworkSelection, SensitiveTransitionState


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


def unevaluable_assessment(**overrides):
    """A network that was never assessed -- the case that had no durable record before increment 38."""
    kwargs = {'network_id': 'network4_2',
              'iteration': 1,
              'status': STATUS_NOT_EVALUATED,
              'reason_code': REASON_SA_ALL_METHODS_FAILED,
              }
    kwargs.update(overrides)
    return PDepNetworkAssessment(**kwargs)


def qualified_assessment(**overrides):
    """A network that qualified, carrying the selection that is the evidence for the verdict."""
    kwargs = {'network_id': 'network4_2',
              'iteration': 1,
              'status': STATUS_QUALIFIED,
              'reason_code': REASON_QUALIFIED_UNCERTAIN_TS,
              }
    kwargs.update(overrides)
    # The evidence has to be about the network the record is about, so the selection follows any
    # network_id the caller overrode rather than being pinned to the default.
    kwargs.setdefault('selection', qualified_selection(network_id=kwargs['network_id']))
    return PDepNetworkAssessment(**kwargs)


def fully_populated_assessment():
    """An assessment with every optional field set, so a round trip cannot pass by dropping one."""
    return PDepNetworkAssessment(network_id='network12_1',
                                 iteration=3,
                                 status=STATUS_EVALUATED_NEGATIVE,
                                 reason_code=REASON_EVALUATED_NO_UNCERTAIN_TS,
                                 secondary_reason_codes=[REASON_SELECTOR_NO_TS_ROWS],
                                 network_path='/runs/t3/network12_1.py',
                                 network_source_hash='sha256:abc123',
                                 observable_label='OH',
                                 sa_rank_index=2,
                                 chemkin_reaction='C2H5(4)+O2(2)=C2H4(6)+HO2(7)',
                                 network_reaction='C2H5 + O2 <=> C2H4 + HO2',
                                 requested_me_methods=['CSE', 'MSC'],
                                 final_method='MSC',
                                 sa_path='/runs/t3/sa_coefficients.yml',
                                 cache_status='cached_valid',
                                 warnings=['a warning', 'another warning'],
                                 selection=negative_selection(network_id='network12_1'),
                                 )


def saved(tmp_path, assessments, file_name='t3_pdep_network_assessments.yml', complete=True):
    """Save ``assessments`` under ``tmp_path`` and return the path written."""
    return save_pdep_network_assessments(os.path.join(str(tmp_path), file_name), assessments,
                                         complete=complete)


def envelope(records, version=ASSESSMENT_ENVELOPE_SCHEMA_VERSION, complete=True):
    """A hand-built on-disk envelope, for the malformed-file cases a saver would never produce."""
    return {'assessment_envelope_schema_version': version,
            'complete': complete,
            'assessments': records}


def write_envelope(tmp_path, content, file_name='hand_written.yml'):
    """Write an arbitrary payload where an assessment record file is expected."""
    path = os.path.join(str(tmp_path), file_name)
    save_yaml_file(path=path, content=content)
    return path


# --- Round trips ---------------------------------------------------------------------------------

def test_a_fully_populated_assessment_survives_a_round_trip_field_for_field(tmp_path):
    original = fully_populated_assessment()
    loaded = load_pdep_network_assessments(saved(tmp_path, [original]))
    assert loaded == [original]


def test_an_unevaluable_network_round_trips_with_no_nested_selection(tmp_path):
    original = unevaluable_assessment()
    loaded = load_pdep_network_assessments(saved(tmp_path, [original]))
    assert loaded == [original]
    assert loaded[0].selection is None


def test_a_qualified_assessment_round_trips_with_its_nested_selection_evidence(tmp_path):
    original = qualified_assessment()
    loaded = load_pdep_network_assessments(saved(tmp_path, [original]))
    assert loaded == [original]
    # The nested selection is the evidence for the verdict; a round trip that kept the verdict and
    # dropped the evidence would leave a record that asserts a conclusion it can no longer support.
    assert loaded[0].selection is not None
    assert loaded[0].selection.qualified is True
    assert [ts.ts_label for ts in loaded[0].selection.selected_ts] == ['TS1']
    assert loaded[0].selection.selected_ts[0].condition == (1000.0, 1.0)


def test_an_internal_error_round_trips_whether_or_not_it_captured_a_selection(tmp_path):
    # A crash can land on either side of the selector, and the record has to survive both.
    with_evidence = PDepNetworkAssessment(network_id='network1_1', status=STATUS_INTERNAL_ERROR,
                                          reason_code=REASON_INTERNAL_ERROR,
                                          selection=negative_selection(network_id='network1_1'))
    without_evidence = PDepNetworkAssessment(network_id='network2_1', status=STATUS_INTERNAL_ERROR,
                                             reason_code=REASON_INTERNAL_ERROR)
    loaded = load_pdep_network_assessments(saved(tmp_path, [with_evidence, without_evidence]))
    assert loaded == [with_evidence, without_evidence]


def test_an_empty_iteration_round_trips_as_an_empty_list(tmp_path):
    # This is the whole reason the file carries an envelope rather than being a bare list: an empty
    # list has no record of its own to carry a version, so the envelope has to.
    assert load_pdep_network_assessments(saved(tmp_path, [])) == []


def test_load_preserves_file_order(tmp_path):
    originals = [unevaluable_assessment(network_id=f'network{index}_1') for index in range(5)]
    loaded = load_pdep_network_assessments(saved(tmp_path, originals))
    assert [record.network_id for record in loaded] == [record.network_id for record in originals]


def test_every_status_and_its_reason_survive_a_round_trip_together(tmp_path):
    originals = [unevaluable_assessment(network_id='n1'),
                 qualified_assessment(network_id='n2'),
                 PDepNetworkAssessment(network_id='n3', status=STATUS_EVALUATED_NEGATIVE,
                                       reason_code=REASON_EVALUATED_NO_UNCERTAIN_TS,
                                       selection=negative_selection(network_id='n3')),
                 PDepNetworkAssessment(network_id='n4', status=STATUS_INTERNAL_ERROR,
                                       reason_code=REASON_INTERNAL_ERROR)]
    loaded = load_pdep_network_assessments(saved(tmp_path, originals))
    assert [(record.status, record.reason_code) for record in loaded] == \
           [(record.status, record.reason_code) for record in originals]


def test_the_loader_is_not_stricter_than_the_constructor_about_cache_status(tmp_path):
    # The constructor takes cache_status as a free optional string, so the loader must too. A loader
    # that enum-restricted it would make a record T3 legitimately wrote unreadable by T3 -- the one
    # failure a durable record must never have.
    original = unevaluable_assessment(cache_status='some_future_cache_state')
    assert load_pdep_network_assessments(saved(tmp_path, [original])) == [original]


# --- The on-disk shape ---------------------------------------------------------------------------

def test_the_saver_returns_the_path_it_wrote_so_callers_can_chain_it(tmp_path):
    path = os.path.join(str(tmp_path), 'assessments.yml')
    assert save_pdep_network_assessments(path, [unevaluable_assessment()], complete=True) == path


def test_the_file_is_a_versioned_envelope_around_an_assessments_list(tmp_path):
    from t3.pdep.api import _read_persisted_yaml_file
    content = _read_persisted_yaml_file(path=saved(tmp_path, [unevaluable_assessment()]))
    assert content['assessment_envelope_schema_version'] == ASSESSMENT_ENVELOPE_SCHEMA_VERSION
    assert isinstance(content['assessments'], list) and len(content['assessments']) == 1
    assert content['assessments'][0]['network_id'] == 'network4_2'


def test_the_file_contains_only_plain_yaml_types(tmp_path):
    # A record written with Python object tags could only be read back by executing them.
    with open(saved(tmp_path, [fully_populated_assessment()]), 'r') as f:
        text = f.read()
    assert '!!python' not in text


def test_saving_replaces_a_previous_record_rather_than_appending_to_it(tmp_path):
    path = os.path.join(str(tmp_path), 'assessments.yml')
    save_pdep_network_assessments(path, [unevaluable_assessment(network_id='old')], complete=True)
    save_pdep_network_assessments(path, [unevaluable_assessment(network_id='new')], complete=True)
    assert [record.network_id for record in load_pdep_network_assessments(path)] == ['new']


def test_the_write_is_atomic_and_leaves_no_staging_file_behind(tmp_path):
    # The funnel rewrites this file once per network, so a crash mid-write is a real scenario, not a
    # theoretical one: a truncated but still parseable record would be trusted as authoritative.
    save_pdep_network_assessments(os.path.join(str(tmp_path), 'assessments.yml'),
                                  [fully_populated_assessment()], complete=True)
    assert sorted(os.listdir(str(tmp_path))) == ['assessments.yml']


def test_a_failed_write_leaves_the_previous_record_intact_rather_than_truncated(tmp_path):
    # The discriminating half of atomicity: a direct write onto `path` would have destroyed the
    # previous record before failing, leaving a truncated-but-parseable file that the next reader
    # would trust. Staging elsewhere and renaming means a failure is a no-op instead.
    import t3.pdep.api as api_module
    path = os.path.join(str(tmp_path), 'assessments.yml')
    save_pdep_network_assessments(path, [unevaluable_assessment(network_id='survivor')], complete=True)
    real_save_yaml_file = api_module.save_yaml_file

    def save_then_die(path, content):
        real_save_yaml_file(path=path, content=content)
        raise OSError('the disk filled up mid-write')

    api_module.save_yaml_file = save_then_die
    try:
        with pytest.raises(OSError):
            save_pdep_network_assessments(path, [unevaluable_assessment(network_id='casualty')],
                                          complete=True)
    finally:
        api_module.save_yaml_file = real_save_yaml_file
    assert [record.network_id for record in load_pdep_network_assessments(path)] == ['survivor']
    assert sorted(os.listdir(str(tmp_path))) == ['assessments.yml']


def test_saving_onto_a_symlink_replaces_the_link_rather_than_writing_through_it(tmp_path):
    outside = os.path.join(str(tmp_path), 'outside.yml')
    with open(outside, 'w') as f:
        f.write('do not clobber me\n')
    link = os.path.join(str(tmp_path), 'assessments.yml')
    os.symlink(outside, link)
    save_pdep_network_assessments(link, [unevaluable_assessment()], complete=True)
    assert not os.path.islink(link)
    with open(outside, 'r') as f:
        assert f.read() == 'do not clobber me\n'


# --- Refusing a file that cannot be trusted ------------------------------------------------------

def test_a_non_mapping_top_level_is_refused(tmp_path):
    path = write_envelope(tmp_path, [{'network_id': 'network4_2'}])
    with pytest.raises(ValueError, match='mapping at its top level'):
        load_pdep_network_assessments(path)


def test_an_unversioned_file_is_refused_rather_than_assumed_to_be_version_one(tmp_path):
    path = write_envelope(tmp_path, {'assessments': []})
    with pytest.raises(ValueError, match='assessment_envelope_schema_version'):
        load_pdep_network_assessments(path)


def test_a_file_of_an_unrecognized_version_is_refused(tmp_path):
    path = write_envelope(tmp_path, envelope([], version=ASSESSMENT_ENVELOPE_SCHEMA_VERSION + 1))
    with pytest.raises(ValueError, match='only understands version'):
        load_pdep_network_assessments(path)


def test_a_boolean_version_is_refused_even_though_true_equals_one(tmp_path):
    # `True == 1` in Python, so a `version: true` file would otherwise pass a plain equality check
    # against schema version 1 and be read as though it were correctly versioned.
    path = write_envelope(tmp_path, envelope([], version=True))
    with pytest.raises(ValueError, match='only understands version'):
        load_pdep_network_assessments(path)


def test_a_missing_assessments_key_is_refused(tmp_path):
    path = write_envelope(tmp_path,
                          {'assessment_envelope_schema_version': ASSESSMENT_ENVELOPE_SCHEMA_VERSION,
                           'complete': True})
    with pytest.raises(ValueError, match="no 'assessments' list"):
        load_pdep_network_assessments(path)


def test_an_assessments_key_that_is_not_a_list_is_refused(tmp_path):
    path = write_envelope(tmp_path, envelope({'network4_2': {}}))
    with pytest.raises(ValueError, match="no 'assessments' list"):
        load_pdep_network_assessments(path)


def test_an_entry_that_is_not_a_mapping_is_refused_and_named_by_its_index(tmp_path):
    path = write_envelope(tmp_path, envelope([['not', 'a', 'mapping']]))
    with pytest.raises(ValueError, match=r'assessments\[0\]'):
        load_pdep_network_assessments(path)


def test_a_null_entry_is_refused(tmp_path):
    path = write_envelope(tmp_path, envelope([None]))
    with pytest.raises(ValueError, match=r'assessments\[0\]'):
        load_pdep_network_assessments(path)


def test_a_record_missing_a_required_field_is_refused_naming_the_field(tmp_path):
    record = fully_populated_assessment().as_dict()
    del record['observable_label']
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match='observable_label'):
        load_pdep_network_assessments(path)


def test_a_record_missing_the_selection_key_is_refused(tmp_path):
    # as_dict() ALWAYS writes this key, rendering it null when there is no selection. An absent key
    # is therefore not the same as an explicit null -- it means the record was not written by T3.
    record = unevaluable_assessment().as_dict()
    del record['selection']
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match='selection'):
        load_pdep_network_assessments(path)


def test_a_record_with_an_unrecognized_reason_code_is_refused(tmp_path):
    record = unevaluable_assessment().as_dict()
    record['reason_code'] = 'the_dog_ate_it'
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match='reason_code'):
        load_pdep_network_assessments(path)


def test_a_record_with_an_unrecognized_status_is_refused(tmp_path):
    record = unevaluable_assessment().as_dict()
    record['status'] = 'probably_fine'
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match='status'):
        load_pdep_network_assessments(path)


def test_a_record_with_an_unrecognized_secondary_reason_code_is_refused(tmp_path):
    record = unevaluable_assessment().as_dict()
    record['secondary_reason_codes'] = ['the_dog_ate_it']
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match='secondary_reason_codes'):
        load_pdep_network_assessments(path)


def test_a_record_whose_status_contradicts_its_reason_code_is_refused_at_load(tmp_path):
    # The invariant the record type exists to enforce has to hold on the way IN as well as on the
    # way out -- a hand-edited file claiming a verdict for a network whose SA never parsed would
    # otherwise mislead with the full authority of provenance.
    record = unevaluable_assessment().as_dict()
    record['status'] = STATUS_QUALIFIED
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_assessments(path)
    assert 'implies' in str(exc_info.value)


def test_a_pre_selector_record_carrying_a_selection_is_refused_at_load(tmp_path):
    record = unevaluable_assessment(reason_code=REASON_NETWORK_PARSE_FAILED).as_dict()
    record['selection'] = negative_selection().as_dict()
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError):
        load_pdep_network_assessments(path)


def test_a_selection_bearing_record_with_a_null_selection_is_refused_at_load(tmp_path):
    record = qualified_assessment().as_dict()
    record['selection'] = None
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError):
        load_pdep_network_assessments(path)


def test_a_refusal_from_the_record_type_still_names_the_file_and_the_entry(tmp_path):
    # The constructor's own message names the network but not the file it came from; on a run with
    # several iterations on disk, that is the difference between a fixable report and a hunt.
    record = unevaluable_assessment().as_dict()
    record['status'] = STATUS_QUALIFIED
    path = write_envelope(tmp_path, envelope([record]), file_name='which_one_was_it.yml')
    with pytest.raises(ValueError) as exc_info:
        load_pdep_network_assessments(path)
    assert 'which_one_was_it.yml' in str(exc_info.value)
    assert 'assessments[0]' in str(exc_info.value)


@pytest.mark.parametrize('key', ['warnings', 'requested_me_methods', 'secondary_reason_codes'])
def test_a_bare_string_where_a_list_belongs_is_refused_rather_than_shredded(tmp_path, key):
    # tuple('MSC') is ('M', 'S', 'C'): coercing a bare string would render three nonsense entries
    # into the record and read back as though T3 had genuinely produced them.
    record = fully_populated_assessment().as_dict()
    record[key] = 'MSC'
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match='string is not accepted'):
        load_pdep_network_assessments(path)


@pytest.mark.parametrize('key', ['network_path', 'observable_label', 'chemkin_reaction',
                                 'network_reaction', 'final_method', 'sa_path', 'cache_status',
                                 'network_source_hash'])
def test_a_non_string_optional_field_is_refused(tmp_path, key):
    record = fully_populated_assessment().as_dict()
    record[key] = 17
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match=key):
        load_pdep_network_assessments(path)


@pytest.mark.parametrize('key', ['iteration', 'sa_rank_index'])
def test_a_non_integer_where_an_integer_belongs_is_refused(tmp_path, key):
    record = fully_populated_assessment().as_dict()
    record[key] = 'two'
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match=key):
        load_pdep_network_assessments(path)


@pytest.mark.parametrize('key', ['iteration', 'sa_rank_index'])
def test_a_boolean_where_an_integer_belongs_is_refused(tmp_path, key):
    record = fully_populated_assessment().as_dict()
    record[key] = True
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match=key):
        load_pdep_network_assessments(path)


def test_a_record_level_version_the_code_does_not_understand_is_refused(tmp_path):
    # The envelope version and the record version are separate keys and are checked separately, so
    # a file whose envelope says 1 cannot smuggle in a record of another shape.
    record = unevaluable_assessment().as_dict()
    record['assessment_record_schema_version'] = ASSESSMENT_RECORD_SCHEMA_VERSION + 1
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match='assessment_record_schema_version'):
        load_pdep_network_assessments(path)


def test_a_nested_selection_of_an_unrecognized_schema_version_is_refused(tmp_path):
    record = qualified_assessment().as_dict()
    record['selection']['selection_schema_version'] = 99
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match='selection_schema_version'):
        load_pdep_network_assessments(path)


def test_a_python_tagged_yaml_file_is_refused_rather_than_executed(tmp_path):
    path = os.path.join(str(tmp_path), 'tagged.yml')
    with open(path, 'w') as f:
        f.write('assessment_envelope_schema_version: 1\n'
                'assessments: !!python/object/apply:os.system ["echo pwned"]\n')
    with pytest.raises(ValueError):
        load_pdep_network_assessments(path)


# --- The envelope version and the record version are independent ---------------------------------

def test_the_envelope_and_the_record_carry_separate_version_keys(tmp_path):
    # Sharing one key would make each version a hostage of the other: adding a field to a record
    # would force the envelope to claim a change it never underwent, and vice versa.
    from t3.pdep.api import _read_persisted_yaml_file
    content = _read_persisted_yaml_file(path=saved(tmp_path, [unevaluable_assessment()]))
    assert content['assessment_envelope_schema_version'] == ASSESSMENT_ENVELOPE_SCHEMA_VERSION
    assert content['assessments'][0]['assessment_record_schema_version'] == \
           ASSESSMENT_RECORD_SCHEMA_VERSION
    assert 'assessment_record_schema_version' not in content
    assert 'assessment_envelope_schema_version' not in content['assessments'][0]


def test_an_envelope_versioned_with_the_record_key_alone_is_refused(tmp_path):
    # The pre-round-57 shape, which reused one key for both roles. It is not readable now, which is
    # correct: nothing on this branch has shipped, so there is no such file to migrate.
    path = write_envelope(tmp_path, {'assessment_record_schema_version': 1, 'assessments': []})
    with pytest.raises(ValueError, match='assessment_envelope_schema_version'):
        load_pdep_network_assessments(path)


# --- Durability ----------------------------------------------------------------------------------

def test_the_staged_file_is_flushed_before_and_the_directory_after_the_rename(tmp_path):
    # Atomic is not durable: os.replace orders the rename, but says nothing about whether the bytes
    # behind it ever reached the disk. Without both flushes a power loss can leave the rename applied
    # and the contents lost -- the file present, and empty or torn.
    import t3.pdep.api as api_module
    fsynced = list()
    real_fsync_path = api_module._fsync_path
    api_module._fsync_path = lambda path: fsynced.append(path)
    try:
        path = save_pdep_network_assessments(os.path.join(str(tmp_path), 'assessments.yml'),
                                             [unevaluable_assessment()], complete=True)
    finally:
        api_module._fsync_path = real_fsync_path
    assert len(fsynced) == 2
    assert fsynced[0].endswith('assessments.yml') and fsynced[0] != path
    assert fsynced[1] == str(tmp_path)


def test_the_staging_directory_is_private_to_this_process(tmp_path):
    # Staging inside a fresh 0700 directory is what closes the window between creating the staged
    # file and reopening it by path to write: nothing else can substitute a file at a path only this
    # process can reach.
    import t3.pdep.api as api_module
    observed = dict()
    real_save_yaml_file = api_module.save_yaml_file

    def record_mode(path, content):
        observed['dir'] = os.path.dirname(path)
        observed['mode'] = os.stat(os.path.dirname(path)).st_mode & 0o777
        return real_save_yaml_file(path=path, content=content)

    api_module.save_yaml_file = record_mode
    try:
        save_pdep_network_assessments(os.path.join(str(tmp_path), 'assessments.yml'),
                                      [unevaluable_assessment()], complete=True)
    finally:
        api_module.save_yaml_file = real_save_yaml_file
    assert observed['mode'] == 0o700
    assert observed['dir'] != str(tmp_path)
    assert not os.path.exists(observed['dir'])


# --- Round-trip totality for the NESTED selection ------------------------------------------------

def test_a_nested_selection_carrying_a_boolean_schema_version_is_refused(tmp_path):
    # True == 1, so a `selection_schema_version: true` nested record would otherwise be read as
    # correctly versioned -- the same hole that was already closed for the assessment's own version.
    record = qualified_assessment().as_dict()
    record['selection']['selection_schema_version'] = True
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match='selection_schema_version'):
        load_pdep_network_assessments(path)


def test_a_nested_selection_carrying_a_boolean_algorithm_version_is_refused(tmp_path):
    record = qualified_assessment().as_dict()
    record['selection']['selection_algorithm_version'] = True
    path = write_envelope(tmp_path, envelope([record]))
    with pytest.raises(ValueError, match='selection_algorithm_version'):
        load_pdep_network_assessments(path)


@pytest.mark.parametrize('field_name, value', [('cache_status', 'some_future_cache_state'),
                                               ('warnings', 'a bare string'),
                                               ('network_reactions_examined', True)])
def test_a_nested_selection_mutated_into_an_unloadable_state_is_refused_at_the_write(tmp_path,
                                                                                     field_name,
                                                                                     value):
    """Test the gap this test used to merely PIN, now that it is closed.

    `PDepNetworkSelection` is mutable by design -- the selector builds it across ~20 assignments --
    so for a long time it validated nothing while `_selection_from_dict` was strict, and a selection
    mutated into one of these states saved cleanly and could then never be loaded back. The old
    comment here said validating at construction would not help, because the mutation happens after
    it, and that the real fix was one contract shared by both halves. That fix is
    `PDepNetworkSelection.validate()`, called from `__post_init__` AND from `as_dict()` -- the one
    funnel every persisted copy passes through, whatever was done to the record in between.

    The refusal now names the field, and lands even earlier than the write: `PDepNetworkAssessment`
    renders the selection when it captures its snapshot, so a mutated selection is refused the
    moment an assessment tries to carry it. A record that cannot be read back can no longer be
    built, never mind saved.
    """
    selection = negative_selection()
    setattr(selection, field_name, value)

    with pytest.raises(ValueError, match=field_name):
        PDepNetworkAssessment(network_id='network4_2', status=STATUS_EVALUATED_NEGATIVE,
                              reason_code=REASON_EVALUATED_NO_UNCERTAIN_TS, selection=selection)

    # And directly at the write, for a selection that is not nested inside anything.
    with pytest.raises(ValueError, match=field_name):
        save_pdep_network_selections(path=os.path.join(str(tmp_path), 'selections.yml'),
                                     selections=[selection])
    assert os.listdir(str(tmp_path)) == [], 'no partial file and no staging directory may survive'


def test_a_nested_selection_holding_a_path_is_refused_by_name_not_by_yaml_tag(tmp_path):
    """Test the value that used to be the WORST case, and no longer is.

    A `Path` in a nested `sa_path` was the one class the loader could not diagnose: `yaml.safe_load`
    refuses the `!!python/object/apply:` tag before any per-record check runs, so a single bad
    nested field cost the whole iteration's assessments and reported a YAML tag rather than a field.
    The selection's own contract now catches it at render time and names `sa_path`.

    The write-time plain-YAML check in `t3.pdep.api` is what remains underneath, for the dict fields
    no contract describes -- see `test_the_check_covers_exactly_what_no_field_validator_describes`.
    """
    selection = negative_selection()
    selection.sa_path = Path('/runs/t3/sa_coefficients.yml')

    with pytest.raises(ValueError, match='sa_path'):
        PDepNetworkAssessment(network_id='network4_2', status=STATUS_EVALUATED_NEGATIVE,
                              reason_code=REASON_EVALUATED_NO_UNCERTAIN_TS, selection=selection)

    with pytest.raises(ValueError, match='sa_path'):
        save_pdep_network_selections(path=os.path.join(str(tmp_path), 'selections.yml'),
                                     selections=[selection])
    assert os.listdir(str(tmp_path)) == [], 'no partial file and no staging directory may survive'


# --- A record cannot be corrupted between construction and the file ------------------------------

def test_mutating_a_records_selection_before_saving_refuses_the_save(tmp_path):
    record = qualified_assessment()
    record.selection.qualified = False
    with pytest.raises(ValueError, match='no longer matches'):
        saved(tmp_path, [record])


def test_a_refused_save_does_not_leave_a_partial_file_or_staging_directory(tmp_path):
    good = unevaluable_assessment(network_id='good')
    corrupted = qualified_assessment()
    corrupted.selection.qualified = False
    with pytest.raises(ValueError):
        saved(tmp_path, [good, corrupted])
    assert os.listdir(str(tmp_path)) == []

# --- Telling a finished record from a crash scene ------------------------------------------------
#
# T3 rewrites this file after every network, so a well-formed file holding four of an iteration's
# twelve networks is an ordinary thing to find on disk. It is also the one corruption that cannot be
# detected by looking at the records, since every one of them is perfectly valid.

def test_the_envelope_states_whether_the_record_is_finished(tmp_path):
    from t3.pdep.api import _read_persisted_yaml_file
    assert _read_persisted_yaml_file(path=saved(tmp_path, [unevaluable_assessment()]))['complete'] is True
    assert _read_persisted_yaml_file(
        path=saved(tmp_path, [unevaluable_assessment()], complete=False))['complete'] is False


def test_the_completeness_flag_has_no_default_so_it_cannot_be_forgotten(tmp_path):
    # Defaulting it either way would be wrong: a forgotten argument would either promote a crash
    # scene to a finished record, or demote every finished record to an unreadable one.
    with pytest.raises(TypeError):
        save_pdep_network_assessments(os.path.join(str(tmp_path), 'assessments.yml'),
                                      [unevaluable_assessment()])


@pytest.mark.parametrize('complete', [1, 0, 'yes', None, [], object()])
def test_a_non_boolean_completeness_flag_is_refused_at_the_write(tmp_path, complete):
    # A truthy stand-in would be written out as itself and refused on the way back in, which is a
    # confusing place to discover it -- and `complete: 1` reads as a value, not as a claim.
    with pytest.raises(ValueError, match='must be a bool'):
        save_pdep_network_assessments(os.path.join(str(tmp_path), 'assessments.yml'),
                                      [unevaluable_assessment()], complete=complete)


def test_an_unfinished_record_is_refused_by_default(tmp_path):
    path = saved(tmp_path, [unevaluable_assessment(network_id='as_far_as_it_got')], complete=False)
    with pytest.raises(ValueError, match='complete=False'):
        load_pdep_network_assessments(path)


def test_an_unfinished_record_can_be_read_deliberately(tmp_path):
    # Which is exactly what an operator investigating the crash wants: the records ARE the evidence.
    path = saved(tmp_path, [unevaluable_assessment(network_id='as_far_as_it_got')], complete=False)
    records = load_pdep_network_assessments(path, allow_incomplete=True)
    assert [record.network_id for record in records] == ['as_far_as_it_got']


def test_allow_incomplete_does_not_weaken_any_other_check(tmp_path):
    # The escape hatch is about the flag alone; a malformed record inside a partial file is still
    # malformed, and reading the partial file must not become a way to smuggle one past the loader.
    path = write_envelope(tmp_path, envelope([{'network_id': 'network4_2'}], complete=False))
    with pytest.raises(ValueError, match='assessment_record_schema_version'):
        load_pdep_network_assessments(path, allow_incomplete=True)


def test_a_file_with_no_completeness_flag_at_all_is_refused(tmp_path):
    # Absence is not "complete": a file written by something that did not know about this flag makes
    # no claim either way, and guessing the reassuring one is how a partial record gets believed.
    path = write_envelope(tmp_path, {'assessment_envelope_schema_version': ASSESSMENT_ENVELOPE_SCHEMA_VERSION,
                                     'assessments': []})
    with pytest.raises(ValueError, match="no 'complete' key"):
        load_pdep_network_assessments(path)


@pytest.mark.parametrize('complete', [1, 0, 'true', None])
def test_a_non_boolean_completeness_flag_is_refused_at_the_read(tmp_path, complete):
    # `complete: 1` would pass a truthiness test and quietly promote an unfinished record.
    path = write_envelope(tmp_path, envelope([], complete=complete))
    with pytest.raises(ValueError, match='must be a bool'):
        load_pdep_network_assessments(path)


def test_a_partial_file_whose_assessments_key_is_not_a_list_is_refused_cleanly(tmp_path):
    # The refusal message counts the records the partial file got to, so the list had to be validated
    # BEFORE it -- otherwise a hand-written `assessments: 1` raised a bare TypeError from inside the
    # message of the error it was about to raise.
    path = write_envelope(tmp_path, envelope(1, complete=False))
    with pytest.raises(ValueError, match="no 'assessments' list"):
        load_pdep_network_assessments(path)
    with pytest.raises(ValueError, match="no 'assessments' list"):
        load_pdep_network_assessments(path, allow_incomplete=True)


@pytest.mark.parametrize('allow_incomplete', [1, 'false', 'yes', None, []])
def test_a_non_boolean_allow_incomplete_is_refused(tmp_path, allow_incomplete):
    # Being strict about `complete` and lax about the flag that overrides it would put the whole
    # guarantee behind a truthiness test: allow_incomplete='false' would read the partial file.
    path = saved(tmp_path, [unevaluable_assessment()], complete=False)
    with pytest.raises(ValueError, match='must be a bool'):
        load_pdep_network_assessments(path, allow_incomplete=allow_incomplete)
