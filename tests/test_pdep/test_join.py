"""
Tests for the t3.pdep.join module, the network-TS <-> ARC-TS join.
"""

import os

import pytest

from t3.pdep.join import (ARC_TS_LABEL_PREFIX,
                          JOIN_STATUS_ALREADY_PRESENT,
                          JOIN_STATUS_NOT_QUEUED,
                          JOIN_STATUS_QUEUED,
                          TSJoinRecord,
                          TS_JOIN_SIDECAR_FILE_NAME,
                          arc_ts_label,
                          expected_ts_artifact_path,
                          merge_ts_join_records,
                          read_ts_join_networks,
                          read_ts_join_sidecar,
                          ts_join_sidecar_path,
                          validate_ts_join_networks,
                          validate_ts_join_records,
                          write_ts_join_sidecar,
                          )

# A syntactically valid source_sha256 in the one representation this module stores: the full
# 'sha256:<hex>' string exactly as t3.pdep.cache.hash_file returns it. join.py validates the
# format only (file existence and re-hashing are the capture step's job), so a literal is enough.
VALID_SHA = f'sha256:{"0123456789abcdef" * 4}'


def _record(network_id='network4_1',
            network_ts_label='TS3',
            status=JOIN_STATUS_QUEUED,
            **kwargs) -> TSJoinRecord:
    """Build a record with the label derived exactly as production code derives it."""
    return TSJoinRecord(network_id=network_id,
                        network_ts_label=network_ts_label,
                        arc_ts_label=kwargs.pop('arc_ts_label', arc_ts_label(network_id, network_ts_label)),
                        status=status,
                        **kwargs)


def _networks_for(records: list, method: str = 'MSC') -> dict:
    """Build a well-formed networks block covering exactly the networks the records reference."""
    return {network_id: {'source_path': f'/rmg/pdep/{network_id}.py',
                         'source_sha256': VALID_SHA,
                         'method': method}
            for network_id in {record.network_id for record in records}}


def test_arc_ts_label_is_deterministic_and_namespaced():
    """Test that the ARC label is a pure function of the network and TS, and is namespaced."""
    label = arc_ts_label(network_id='network4_1', network_ts_label='TS3')
    assert label == 'T3PDep_network4_1_TS3'
    assert label == arc_ts_label(network_id='network4_1', network_ts_label='TS3')
    assert label.startswith(ARC_TS_LABEL_PREFIX)


def test_arc_ts_label_cannot_collide_with_arcs_own_namespace():
    """Test that a T3 label is never of the ``TS<index>`` form ARC assigns to unlabeled reactions.

    If it were, ARC could hand the same label to a reaction T3 did not queue, and the artifact read
    back for a network transition state would be some other reaction's result.
    """
    label = arc_ts_label(network_id='network1_1', network_ts_label='TS1')
    assert not (label.startswith('TS') and label[2:].isdigit())


def test_arc_ts_label_distinguishes_networks_and_transition_states():
    """Test that the same TS label in different networks maps to different ARC labels, and vice versa."""
    assert arc_ts_label('network4_1', 'TS1') != arc_ts_label('network4_2', 'TS1')
    assert arc_ts_label('network4_1', 'TS1') != arc_ts_label('network4_1', 'TS2')


@pytest.mark.parametrize('network_id, network_ts_label', [
    ('net work', 'TS1'),
    ('net-1', 'TS1'),
    ('network4_1', 'TS 1'),
    ('network4_1', 'TS(1)'),
    ('net/../etc', 'TS1'),
])
def test_arc_ts_label_refuses_unsafe_components(network_id, network_ts_label):
    """Test that an unsafe component is refused rather than sanitized.

    Sanitizing is the dangerous option: ``net-1`` and ``net_1`` would both become ``net_1`` and two
    distinct networks would then share one ARC label. The path-traversal case matters separately,
    since the label becomes a file name under the ARC project directory.
    """
    with pytest.raises(ValueError):
        arc_ts_label(network_id=network_id, network_ts_label=network_ts_label)


@pytest.mark.parametrize('network_id, network_ts_label', [('', 'TS1'), ('network4_1', '')])
def test_arc_ts_label_refuses_empty_components(network_id, network_ts_label):
    """Test that an empty component is refused rather than producing a degenerate label."""
    with pytest.raises(ValueError):
        arc_ts_label(network_id=network_id, network_ts_label=network_ts_label)


def test_expected_ts_artifact_path_matches_arcs_layout():
    """Test the artifact path ARC writes per transition state."""
    path = expected_ts_artifact_path(arc_project_directory='/proj', label='T3PDep_network4_1_TS3')
    assert path == os.path.join('/proj', 'calcs', 'statmech', 'kinetics', 'TSs',
                                'T3PDep_network4_1_TS3.py')


def test_record_rejects_an_unrecognized_status():
    """Test that a record cannot be built with a status the discovery step would not understand."""
    with pytest.raises(ValueError, match='status'):
        _record(status='probably_fine')


def test_record_round_trips_through_as_dict():
    """Test that a record survives the YAML rendering it is stored as."""
    record = _record(path_reaction_labels=('reaction1', 'reaction4'),
                     path_reaction_strs=('A <=> B', 'A <=> C + D'),
                     t3_reaction_key=7,
                     expected_artifact_path='/proj/calcs/statmech/kinetics/TSs/T3PDep_network4_1_TS3.py',
                     reason='selected')
    assert TSJoinRecord.from_dict(record.as_dict()) == record


def test_record_as_dict_contains_no_tuples():
    """Test that the rendering is plain YAML-safe types, as the sidecar is written verbatim."""
    rendered = _record(path_reaction_labels=('reaction1',), path_reaction_strs=('A <=> B',)).as_dict()
    assert not any(isinstance(value, tuple) for value in rendered.values())


@pytest.mark.parametrize('missing', ['network_id', 'network_ts_label', 'status'])
def test_record_from_dict_refuses_a_missing_required_field(missing):
    """Test that a truncated sidecar entry is refused rather than silently read as a partial join.

    ``arc_ts_label`` is deliberately not parametrized here: it is now legitimately ``None`` for a
    transition state ``arc_ts_label()`` refused to label, so a missing/``None`` value is not itself
    an error (see ``test_record_round_trips_with_no_arc_ts_label``).
    """
    rendered = _record().as_dict()
    rendered[missing] = None
    with pytest.raises(ValueError, match=missing):
        TSJoinRecord.from_dict(rendered)


def test_record_round_trips_with_no_arc_ts_label():
    """Test that a record whose ``arc_ts_label`` is `None` survives the YAML rendering it is stored as.

    This is the ``arc_ts_label()`` refusal case: the transition state still gets a record, just with
    no ARC label, and that record must serialize and deserialize without being mistaken for a
    truncated one.
    """
    record = TSJoinRecord(network_id='network4_1',
                          network_ts_label='TS3',
                          status=JOIN_STATUS_NOT_QUEUED,
                          arc_ts_label=None,
                          reason='network_id contains unsafe characters')
    assert TSJoinRecord.from_dict(record.as_dict()) == record
    assert TSJoinRecord.from_dict(record.as_dict()).arc_ts_label is None


def test_validate_refuses_one_network_ts_mapped_twice():
    """Test that a duplicated network TS is refused: T3 would not know which result to read back."""
    records = [_record(), _record(arc_ts_label='T3PDep_network4_1_TS3_other')]
    with pytest.raises(ValueError, match='Ambiguous'):
        validate_ts_join_records(records)


def test_validate_refuses_one_arc_label_claimed_twice():
    """Test that a shared ARC label is refused.

    This is the direction that silently corrupts: both network transition states would read back
    the SAME artifact, attributing one network's quantum chemistry to another.
    """
    shared = arc_ts_label('network4_1', 'TS3')
    records = [_record(network_ts_label='TS3'),
               _record(network_ts_label='TS4', arc_ts_label=shared)]
    with pytest.raises(ValueError, match='Ambiguous'):
        validate_ts_join_records(records)


def test_validate_accepts_two_records_with_no_arc_ts_label():
    """Test that two `arc_ts_label=None` records do not collide with each other.

    Both `arc_ts_label()` refusals get `arc_ts_label=None`, and `None` is not a real ARC label two
    transition states could actually share, so it must be exempt from the ARC-label uniqueness check.
    """
    records = [_record(network_ts_label='TS3', arc_ts_label=None, status=JOIN_STATUS_NOT_QUEUED),
              _record(network_ts_label='TS4', arc_ts_label=None, status=JOIN_STATUS_NOT_QUEUED)]
    validate_ts_join_records(records)


def test_validate_still_refuses_a_duplicated_network_ts_with_no_arc_ts_label():
    """Test that network-key uniqueness is still enforced even when `arc_ts_label` is `None`."""
    records = [_record(network_ts_label='TS3', arc_ts_label=None, status=JOIN_STATUS_NOT_QUEUED),
              _record(network_ts_label='TS3', arc_ts_label=None, status=JOIN_STATUS_NOT_QUEUED)]
    with pytest.raises(ValueError, match='Ambiguous'):
        validate_ts_join_records(records)


def test_validate_accepts_a_one_to_one_mapping():
    """Test that a well-formed set of records passes."""
    validate_ts_join_records([_record(network_ts_label='TS3'),
                              _record(network_ts_label='TS4'),
                              _record(network_id='network4_2', network_ts_label='TS3')])


def test_merge_absorbs_an_identical_repeat():
    """Test that offering the same TS twice is absorbed, not duplicated.

    Several network reactions can point at one network, so the same transition state is legitimately
    offered more than once within a single iteration.
    """
    merged = merge_ts_join_records(existing=[_record()], new=[_record()])
    assert len(merged) == 1


def test_merge_refuses_a_conflicting_repeat():
    """Test that a repeat which disagrees is refused rather than last-write-wins."""
    with pytest.raises(ValueError, match='Conflicting'):
        merge_ts_join_records(existing=[_record(t3_reaction_key=1)],
                              new=[_record(t3_reaction_key=2)])


def test_merge_preserves_first_seen_order():
    """Test that the sidecar reads in the order the work was decided."""
    merged = merge_ts_join_records(existing=[_record(network_ts_label='TS3')],
                                   new=[_record(network_ts_label='TS1'),
                                        _record(network_ts_label='TS2')])
    assert [record.network_ts_label for record in merged] == ['TS3', 'TS1', 'TS2']


def test_sidecar_round_trips(tmp_path):
    """Test that records survive a write/read cycle through the sidecar."""
    project = str(tmp_path / 'ARC')
    records = [_record(network_ts_label='TS3', t3_reaction_key=4,
                       path_reaction_labels=('reaction1',), path_reaction_strs=('A <=> B',)),
               _record(network_ts_label='TS4', status=JOIN_STATUS_ALREADY_PRESENT, reason='already known'),
               _record(network_id='network1_1', network_ts_label='TS1',
                       status=JOIN_STATUS_NOT_QUEUED, reason='no structure')]
    networks = _networks_for(records)
    path = write_ts_join_sidecar(arc_project_directory=project, records=records, networks=networks)
    assert path == ts_join_sidecar_path(project)
    assert os.path.basename(path) == TS_JOIN_SIDECAR_FILE_NAME
    assert read_ts_join_sidecar(project) == records


def test_reading_an_absent_sidecar_is_empty_not_an_error(tmp_path):
    """Test that a project that queued no PDep transition states reads as empty.

    An absent sidecar means "nothing was queued", which is a normal outcome and must not be
    confused with a lost join.
    """
    assert read_ts_join_sidecar(str(tmp_path)) == list()


def test_reading_a_malformed_sidecar_raises(tmp_path):
    """Test that a sidecar which is not the expected mapping is refused rather than read as empty."""
    project = tmp_path / 'ARC'
    project.mkdir()
    (project / TS_JOIN_SIDECAR_FILE_NAME).write_text('- just\n- a list\n')
    with pytest.raises(ValueError, match='malformed'):
        read_ts_join_sidecar(str(project))


def test_reading_a_sidecar_with_an_ambiguous_mapping_raises(tmp_path):
    """Test that ambiguity is caught on the way IN as well as on the way out.

    The sidecar is a file on disk between two ARC-length steps; it can be edited or merged by hand
    in between, so the reader cannot assume the writer's guarantees still hold.
    """
    project = tmp_path / 'ARC'
    project.mkdir()
    shared = arc_ts_label('network4_1', 'TS3')
    (project / TS_JOIN_SIDECAR_FILE_NAME).write_text(
        'transition_states:\n'
        f'- {{network_id: network4_1, network_ts_label: TS3, arc_ts_label: {shared}, status: queued}}\n'
        f'- {{network_id: network4_1, network_ts_label: TS4, arc_ts_label: {shared}, status: queued}}\n'
        'networks:\n'
        f'  network4_1: {{source_path: /rmg/pdep/network4_1.py, source_sha256: {VALID_SHA}, method: MSC}}\n')
    with pytest.raises(ValueError, match='Ambiguous'):
        read_ts_join_sidecar(str(project))


def test_write_refuses_ambiguous_records_without_leaving_a_partial_sidecar(tmp_path):
    """Test that a refused write leaves nothing behind for a later step to misread as authoritative."""
    project = str(tmp_path / 'ARC')
    shared = arc_ts_label('network4_1', 'TS3')
    records = [_record(network_ts_label='TS3'),
               _record(network_ts_label='TS4', arc_ts_label=shared)]
    with pytest.raises(ValueError, match='Ambiguous'):
        write_ts_join_sidecar(arc_project_directory=project,
                              records=records,
                              networks=_networks_for(records))
    assert not os.path.isfile(ts_join_sidecar_path(project))


def test_networks_block_round_trips_losslessly(tmp_path):
    """Test that the networks block survives a write/read cycle byte-equal in content."""
    project = str(tmp_path / 'ARC')
    records = [_record(network_ts_label='TS3'),
               _record(network_id='network1_1', network_ts_label='TS1')]
    networks = {'network4_1': {'source_path': '/rmg/pdep/network4_1.py',
                               'source_sha256': VALID_SHA,
                               'method': 'MSC'},
                'network1_1': {'source_path': '/rmg/pdep/network1_1.py',
                               'source_sha256': f'sha256:{"fedcba9876543210" * 4}',
                               'method': 'CSE'}}
    write_ts_join_sidecar(arc_project_directory=project, records=records, networks=networks)
    assert read_ts_join_networks(project) == networks


def test_read_networks_of_an_absent_sidecar_is_empty_not_an_error(tmp_path):
    """Test that a project that queued nothing has an empty networks mapping, mirroring records."""
    assert read_ts_join_networks(str(tmp_path)) == dict()


def test_reading_a_sidecar_without_a_networks_key_raises(tmp_path):
    """Test fail-closed on a sidecar predating (or stripped of) the networks block.

    No backward compatibility is offered on this branch: a sidecar that names transition states but
    cannot say which network file (at which hash, refined with which ME method) they came from is
    not a partial answer, it is no answer -- the capture step could vendor the wrong network file.
    """
    project = tmp_path / 'ARC'
    project.mkdir()
    label = arc_ts_label('network4_1', 'TS3')
    (project / TS_JOIN_SIDECAR_FILE_NAME).write_text(
        'transition_states:\n'
        f'- {{network_id: network4_1, network_ts_label: TS3, arc_ts_label: {label}, status: queued}}\n')
    with pytest.raises(ValueError, match="'networks' key"):
        read_ts_join_sidecar(str(project))
    with pytest.raises(ValueError, match="'networks' key"):
        read_ts_join_networks(str(project))


def test_a_sidecar_with_empty_records_and_no_networks_key_still_raises(tmp_path):
    """Test that the missing-key refusal is not shadowed by the cross-reference check.

    With zero records the cross-reference between records and networks is vacuously satisfied, so
    only the explicit key-presence check stands between this file and being silently read as a
    valid empty sidecar -- and a sidecar stripped of its ``networks`` key is malformed regardless
    of how many records it holds.
    """
    project = tmp_path / 'ARC'
    project.mkdir()
    (project / TS_JOIN_SIDECAR_FILE_NAME).write_text('transition_states: []\n')
    with pytest.raises(ValueError, match="'networks' key"):
        read_ts_join_sidecar(str(project))


def test_a_record_referencing_an_unknown_network_raises():
    """Test the record -> networks direction: every referenced network must have an entry."""
    records = [_record(), _record(network_id='network1_1', network_ts_label='TS1')]
    networks = _networks_for([records[0]])  # network1_1 deliberately absent
    with pytest.raises(ValueError, match='network1_1'):
        validate_ts_join_networks(records=records, networks=networks)


def test_an_unreferenced_networks_entry_raises():
    """Test the networks -> record direction: an entry no record references is refused.

    An unreferenced entry means the sidecar asserts an identity for work it never actually
    queued -- either the records were dropped or the entry belongs to another iteration; both
    are corruption, not noise.
    """
    records = [_record()]
    networks = _networks_for(records)
    networks['network9_9'] = {'source_path': '/rmg/pdep/network9_9.py',
                              'source_sha256': VALID_SHA,
                              'method': 'MSC'}
    with pytest.raises(ValueError, match='network9_9'):
        validate_ts_join_networks(records=records, networks=networks)


def test_networks_entries_required_even_when_all_records_were_not_queued():
    """Test that not_queued records still pin their network's identity: they reference it too."""
    records = [_record(status=JOIN_STATUS_NOT_QUEUED, arc_ts_label=None, reason='refused')]
    with pytest.raises(ValueError, match='network4_1'):
        validate_ts_join_networks(records=records, networks=dict())


@pytest.mark.parametrize('method', ['msc', 'modified strong collision', 'default', '', None])
def test_an_unrecognized_method_raises(method):
    """Test that only a method key t3.utils.writer.METHOD_MAP accepts is allowed.

    The method is what the downstream hybrid Arkane input will be written with; an unknown value
    here would only fail later, deep inside Arkane input generation, far from its cause.
    """
    records = [_record()]
    networks = _networks_for(records)
    networks['network4_1']['method'] = method
    with pytest.raises(ValueError, match='method'):
        validate_ts_join_networks(records=records, networks=networks)


@pytest.mark.parametrize('bad_sha', [
    '0123456789abcdef' * 4,             # bare hex digest, missing the 'sha256:' prefix
    'sha256:' + '0123456789ABCDEF' * 4,  # uppercase hex
    'sha256:' + '0123456789abcdef' * 3,  # too short
    'sha256:' + 'ghijklmnopqrstuv' * 4,  # not hex at all
    'md5:' + '0123456789abcdef' * 4,     # wrong algorithm prefix
    '',
    None,
])
def test_a_malformed_source_sha256_raises(bad_sha):
    """Test strict validation of the stored hash representation: 'sha256:' + 64 lowercase hex."""
    records = [_record()]
    networks = _networks_for(records)
    networks['network4_1']['source_sha256'] = bad_sha
    with pytest.raises(ValueError, match='source_sha256'):
        validate_ts_join_networks(records=records, networks=networks)


@pytest.mark.parametrize('bad_path', ['', '   ', None, 42])
def test_a_blank_source_path_raises(bad_path):
    """Test that source_path must be a non-blank string."""
    records = [_record()]
    networks = _networks_for(records)
    networks['network4_1']['source_path'] = bad_path
    with pytest.raises(ValueError, match='source_path'):
        validate_ts_join_networks(records=records, networks=networks)


def test_write_refuses_a_bad_networks_block_without_leaving_a_partial_sidecar(tmp_path):
    """Test that networks validation happens before anything is written, like record validation."""
    project = str(tmp_path / 'ARC')
    records = [_record()]
    with pytest.raises(ValueError, match='network4_1'):
        write_ts_join_sidecar(arc_project_directory=project, records=records, networks=dict())
    assert not os.path.isfile(ts_join_sidecar_path(project))


def test_reading_a_sidecar_with_an_invalid_networks_entry_raises(tmp_path):
    """Test that networks validation is enforced on the way IN as well.

    The sidecar sits on disk between two ARC-length steps and can be edited in between, so the
    reader re-validates rather than assuming the writer's guarantees still hold.
    """
    project = tmp_path / 'ARC'
    project.mkdir()
    label = arc_ts_label('network4_1', 'TS3')
    (project / TS_JOIN_SIDECAR_FILE_NAME).write_text(
        'transition_states:\n'
        f'- {{network_id: network4_1, network_ts_label: TS3, arc_ts_label: {label}, status: queued}}\n'
        'networks:\n'
        f'  network4_1: {{source_path: /rmg/pdep/network4_1.py, source_sha256: {VALID_SHA}, '
        f'method: not_a_method}}\n')
    with pytest.raises(ValueError, match='method'):
        read_ts_join_sidecar(str(project))
