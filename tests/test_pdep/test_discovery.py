"""
Tests for t3.pdep.discovery.

Every fixture is built under ``tmp_path``; none reads from ``tests/data/``, since these tests
exercise how ``discovery`` reacts to the *shape* of an ARC project directory, not to any one
recorded run of it.
"""

import os

import pytest

from t3.pdep.discovery import (
    ARTIFACT_STATUS_MISSING,
    ARTIFACT_STATUS_NOT_QUEUED,
    ARTIFACT_STATUS_UNUSABLE,
    ARTIFACT_STATUS_UNVERIFIED,
    ARTIFACT_STATUS_USABLE,
    HYBRID_STATUS_COMPLETE,
    HYBRID_STATUS_NOT_EVALUATED,
    HYBRID_STATUS_PARTIAL_SELECTED_QM,
    TSArtifactRecord,
    discover_ts_artifacts,
    evaluate_pdep_hybrid,
    read_arc_convergence,
)
from t3.pdep.join import JOIN_STATUS_ALREADY_PRESENT, JOIN_STATUS_NOT_QUEUED, JOIN_STATUS_QUEUED, TSJoinRecord, \
    arc_ts_label, expected_ts_artifact_path
from t3.pdep.selector import SensitiveTransitionState


def _write_artifact(path: str, log_path: str | None = None) -> None:
    """Write a minimal artifact at ``path`` in the ARC statmech species-file shape
    (``arc/statmech/arkane.py::species_input_template``): module-level ``linear``,
    ``spinMultiplicity``, ``energy``, ``geometry`` and ``frequencies`` assignments, the
    ``Log(...)`` ones pointing at ``log_path`` (relative to the artifact's own directory).

    When ``log_path`` is ``None``, a stub log file is created next to the artifact so its
    ``Log(...)`` references still resolve.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if log_path is None:
        log_path = os.path.join(os.path.dirname(path), 'output.out')
        with open(log_path, 'w') as f:
            f.write('stub quantum chemistry log\n')
    rel = os.path.relpath(log_path, os.path.dirname(path))
    content = f"""linear = False

spinMultiplicity = 2

energy = Log('{rel}')

geometry = Log('{rel}')

frequencies = Log('{rel}')
"""
    with open(path, 'w') as f:
        f.write(content)


def _queued_record(network_id='network4_1', network_ts_label='TS3', arc_project_directory=None):
    label = arc_ts_label(network_id, network_ts_label)
    expected_path = expected_ts_artifact_path(arc_project_directory, label) if arc_project_directory else None
    return TSJoinRecord(
        network_id=network_id,
        network_ts_label=network_ts_label,
        status=JOIN_STATUS_QUEUED,
        arc_ts_label=label,
        expected_artifact_path=expected_path,
        reason='Queued to ARC.',
    )


class TestReadArcConvergence:

    def test_no_status_source_returns_none(self, tmp_path):
        assert read_arc_convergence(str(tmp_path)) is None

    def test_prefers_status_yml_over_output_yml(self, tmp_path):
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        # status.yml says explicitly NOT converged; output.yml says converged -- status.yml is
        # the richer, tri-state-preserving source and must win. If output.yml were preferred
        # instead, this test would see its verdict.
        (output_dir / 'status.yml').write_text(
            "T3PDep_network4_1_TS3:\n  convergence: false\n  job_types: {}\n  paths: {}\n  info: ''\n  errors: ''\n"
        )
        (output_dir / 'output.yml').write_text(
            "species: []\n"
            "transition_states:\n"
            "  - label: T3PDep_network4_1_TS3\n"
            "    converged: true\n"
        )
        convergence = read_arc_convergence(str(tmp_path))
        assert convergence is not None
        assert convergence['T3PDep_network4_1_TS3']['converged'] is False

    def test_stale_status_yml_converged_vs_output_yml_not_converged_fails_closed(self, tmp_path):
        # The one genuinely dangerous disagreement direction: status.yml (possibly stale, from an
        # interrupted earlier run) says converged while output.yml (ARC's final consolidated
        # summary) says NOT converged. Fail closed: treat the transition state as not converged,
        # with a reason naming both sources and the disagreement.
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        (output_dir / 'status.yml').write_text(
            "T3PDep_network4_1_TS3:\n  convergence: true\n  job_types: {}\n  paths: {}\n  info: ''\n  errors: ''\n"
        )
        (output_dir / 'output.yml').write_text(
            "species: []\n"
            "transition_states:\n"
            "  - label: T3PDep_network4_1_TS3\n"
            "    converged: false\n"
        )
        convergence = read_arc_convergence(str(tmp_path))
        assert convergence is not None
        entry = convergence['T3PDep_network4_1_TS3']
        assert entry['converged'] is False
        assert 'status.yml' in entry['reason']
        assert 'output.yml' in entry['reason']
        assert 'disagree' in entry['reason'].lower()

    def test_status_yml_unknown_vs_output_yml_not_converged_is_not_a_conflict(self, tmp_path):
        # status.yml: null vs output.yml: false is the EXPECTED lossy-collapse case (output.yml
        # writes `converged: entry.get('convergence') is True`, folding None into False), not a
        # real conflict. The tri-state 'unknown' must survive, and no disagreement is reported.
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        (output_dir / 'status.yml').write_text(
            "T3PDep_network4_1_TS3:\n  convergence: null\n  job_types: {}\n  paths: {}\n  info: ''\n  errors: ''\n"
        )
        (output_dir / 'output.yml').write_text(
            "species: []\n"
            "transition_states:\n"
            "  - label: T3PDep_network4_1_TS3\n"
            "    converged: false\n"
        )
        convergence = read_arc_convergence(str(tmp_path))
        assert convergence is not None
        entry = convergence['T3PDep_network4_1_TS3']
        assert entry['converged'] is None
        assert 'disagree' not in entry['reason'].lower()

    def test_fresh_status_yml_converged_beats_stale_output_yml_not_converged(self, tmp_path):
        # Finding 19: the reverse of the dangerous direction above. output.yml is written LAST in
        # a single ARC execute() call (status.yml unconditionally first, output.yml later via
        # write_output_yml) -- so when the two disagree, a status.yml that is strictly NEWER than
        # output.yml cannot be the stale one; it is evidence of a later, restarted run that got
        # further than whatever left output.yml behind (e.g. the restart crashed again before
        # write_output_yml was reached). Treating output.yml as automatically authoritative here
        # would reject a genuinely valid final result. status.yml's converged=True must win when
        # it is demonstrably fresher.
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        output_path = output_dir / 'output.yml'
        status_path = output_dir / 'status.yml'
        # Write output.yml first (the stale one, from an earlier interrupted run).
        output_path.write_text(
            "species: []\n"
            "transition_states:\n"
            "  - label: T3PDep_network4_1_TS3\n"
            "    converged: false\n"
        )
        # status.yml is written afterward, by a later restarted run, and reports convergence.
        status_path.write_text(
            "T3PDep_network4_1_TS3:\n  convergence: true\n  job_types: {}\n  paths: {}\n  info: ''\n  errors: ''\n"
        )
        # Make the freshness relationship explicit and unambiguous regardless of filesystem
        # mtime resolution.
        now = os.path.getmtime(status_path)
        os.utime(output_path, (now - 10, now - 10))
        os.utime(status_path, (now, now))
        convergence = read_arc_convergence(str(tmp_path))
        assert convergence is not None
        entry = convergence['T3PDep_network4_1_TS3']
        assert entry['converged'] is True

    def test_output_yml_non_convergence_fails_closed_with_caveat(self, tmp_path):
        # Only output.yml exists (no status.yml), and it reports not-converged. output.yml cannot
        # distinguish "explicitly failed" from "never assessed", so this must be treated
        # fail-closed as not converged -- but the reason text must say so honestly rather than
        # claiming a definite failure the source cannot actually support.
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        (output_dir / 'output.yml').write_text(
            "species: []\n"
            "transition_states:\n"
            "  - label: T3PDep_network4_1_TS3\n"
            "    converged: false\n"
        )
        convergence = read_arc_convergence(str(tmp_path))
        assert convergence is not None
        entry = convergence['T3PDep_network4_1_TS3']
        assert entry['converged'] is False
        assert 'could not distinguish' in entry['reason']

    def test_reads_status_yml_when_it_is_the_only_source(self, tmp_path):
        # Only status.yml exists (no output.yml): its verdict is read directly. (This is the
        # status-only case -- nothing here "falls back" to anything.)
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        (output_dir / 'status.yml').write_text(
            "T3PDep_network4_1_TS3:\n"
            "  convergence: true\n"
            "  job_types: {}\n"
            "  paths: {}\n"
            "  info: ''\n"
            "  errors: ''\n"
        )
        convergence = read_arc_convergence(str(tmp_path))
        assert convergence is not None
        assert convergence['T3PDep_network4_1_TS3']['converged'] is True

    @pytest.mark.parametrize('yaml_value, type_name',
                             [('"false"', 'str'), ('"true"', 'str'), ('1', 'int'), ('[true]', 'list')])
    def test_status_yml_non_boolean_convergence_is_unknown_never_coerced(self, tmp_path, yaml_value, type_name):
        # A YAML STRING like "false" (quoted, hand-edited, or written by a different ARC version)
        # is truthy, so a bool(...) coercion would silently promote a failed transition state to
        # converged. Anything that is not literally True/False/None must map to UNKNOWN (None),
        # with a reason naming the unexpected value and its type.
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        (output_dir / 'status.yml').write_text(
            f"T3PDep_network4_1_TS3:\n"
            f"  convergence: {yaml_value}\n"
            f"  job_types: {{}}\n"
            f"  paths: {{}}\n"
            f"  info: ''\n"
            f"  errors: ''\n"
        )
        convergence = read_arc_convergence(str(tmp_path))
        assert convergence is not None
        entry = convergence['T3PDep_network4_1_TS3']
        assert entry['converged'] is None
        assert 'unexpected' in entry['reason'].lower()
        assert type_name in entry['reason']

    @pytest.mark.parametrize('yaml_value, type_name',
                             [('"false"', 'str'), ('"true"', 'str'), ('1', 'int')])
    def test_output_yml_non_boolean_converged_is_unknown_never_coerced(self, tmp_path, yaml_value, type_name):
        # Same guarantee for the output.yml fallback: a non-boolean 'converged' value must never
        # be coerced with bool(...); it maps to UNKNOWN (None) with a reason naming value and type.
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        (output_dir / 'output.yml').write_text(
            f"species: []\n"
            f"transition_states:\n"
            f"  - label: T3PDep_network4_1_TS3\n"
            f"    converged: {yaml_value}\n"
        )
        convergence = read_arc_convergence(str(tmp_path))
        assert convergence is not None
        entry = convergence['T3PDep_network4_1_TS3']
        assert entry['converged'] is None
        assert 'unexpected' in entry['reason'].lower()
        assert type_name in entry['reason']

    def test_status_yml_unknown_convergence_is_none_not_false(self, tmp_path):
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        (output_dir / 'status.yml').write_text(
            "T3PDep_network4_1_TS3:\n"
            "  convergence: null\n"
            "  job_types: {}\n"
            "  paths: {}\n"
            "  info: ''\n"
            "  errors: ''\n"
        )
        convergence = read_arc_convergence(str(tmp_path))
        assert convergence is not None
        assert convergence['T3PDep_network4_1_TS3']['converged'] is None


class TestDiscoverTsArtifacts:

    def test_usable_artifact(self, tmp_path):
        # Finding 18: USABLE requires convergence to be explicitly True, not merely unknown -- so
        # this (the only test exercising a genuine positive-convergence -> USABLE outcome through
        # discover_ts_artifacts, as opposed to constructing TSArtifactRecord directly) must supply
        # a real status.yml entry reporting convergence.
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        _write_artifact(record.expected_artifact_path)
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        (output_dir / 'status.yml').write_text(
            f"{record.arc_ts_label}:\n"
            "  convergence: true\n"
            "  job_types: {}\n"
            "  paths: {}\n"
            "  info: ''\n"
            "  errors: ''\n"
        )
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert len(artifacts) == 1
        assert artifacts[0].status == ARTIFACT_STATUS_USABLE
        assert artifacts[0].converged is True
        assert artifacts[0].artifact_path is not None
        assert os.path.isfile(artifacts[0].artifact_path)

    def test_missing_artifact(self, tmp_path):
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert artifacts[0].status == ARTIFACT_STATUS_MISSING
        assert artifacts[0].artifact_path is None

    def test_artifact_present_but_referenced_log_deleted_is_unusable(self, tmp_path):
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        log_path = str(tmp_path / 'calcs' / 'statmech' / 'kinetics' / 'TSs' / 'some.log')
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        with open(log_path, 'w') as f:
            f.write('fake log')
        _write_artifact(record.expected_artifact_path, log_path=log_path)
        # Delete the log AFTER writing the artifact that references it.
        os.remove(log_path)
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert artifacts[0].status == ARTIFACT_STATUS_UNUSABLE
        assert 'Log' in artifacts[0].reason or 'log' in artifacts[0].reason

    def test_explicitly_unconverged_is_unusable(self, tmp_path):
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        _write_artifact(record.expected_artifact_path)
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        (output_dir / 'status.yml').write_text(
            f"{record.arc_ts_label}:\n"
            "  convergence: false\n"
            "  job_types: {}\n"
            "  paths: {}\n"
            "  info: 'optimization failed'\n"
            "  errors: 'SCF did not converge'\n"
        )
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert artifacts[0].status == ARTIFACT_STATUS_UNUSABLE
        assert artifacts[0].converged is False

    def test_no_status_source_at_all_is_not_promoted_to_usable(self, tmp_path):
        # Finding 18 (TDD red phase): no output/status.yml AND no output/output.yml exist at all
        # (read_arc_convergence returns None -- convergence entirely UNKNOWN, e.g. immediately
        # after arc.execute() left no status record whatsoever). Even though the artifact parses
        # cleanly, an unknown convergence must NOT be silently promoted to 'usable' -- that would
        # let a QM result T3 never actually confirmed converged flow straight into a hybrid model.
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        _write_artifact(record.expected_artifact_path)
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert artifacts[0].status == ARTIFACT_STATUS_UNVERIFIED
        assert artifacts[0].converged is None

    def test_unknown_convergence_is_unverified_not_usable_if_artifact_parses(self, tmp_path):
        # Finding 18: convergence unknown (label absent from the only status source that exists)
        # must not be silently promoted to USABLE -- it is UNVERIFIED: the artifact is fine, but
        # nothing ever confirmed it converged.
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        _write_artifact(record.expected_artifact_path)
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        # No entry at all for this label -- convergence must be treated as unknown, not unconverged.
        (output_dir / 'status.yml').write_text('SomeOtherLabel:\n  convergence: true\n  job_types: {}\n'
                                                "  paths: {}\n  info: ''\n  errors: ''\n")
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert artifacts[0].status == ARTIFACT_STATUS_UNVERIFIED
        assert artifacts[0].converged is None

    def test_already_present_record_preserved(self, tmp_path):
        arc_dir = str(tmp_path)
        record = TSJoinRecord(
            network_id='network4_1',
            network_ts_label='TS3',
            status=JOIN_STATUS_ALREADY_PRESENT,
            reason='This transition state was already known to T3.',
        )
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert len(artifacts) == 1
        assert artifacts[0].status == ARTIFACT_STATUS_MISSING
        assert artifacts[0].artifact_path is None
        assert artifacts[0].reason == record.reason

    def test_not_queued_record_preserved(self, tmp_path):
        arc_dir = str(tmp_path)
        record = TSJoinRecord(
            network_id='network4_1',
            network_ts_label='TS3',
            status=JOIN_STATUS_NOT_QUEUED,
            reason='Could not build the species for this transition state.',
        )
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert len(artifacts) == 1
        assert artifacts[0].status == ARTIFACT_STATUS_NOT_QUEUED
        assert artifacts[0].artifact_path is None
        assert artifacts[0].reason == record.reason

    def test_artifact_referencing_log_outside_project_is_unusable(self, tmp_path):
        # An artifact whose Log(...) resolves outside the ARC project directory is refused by
        # the reader (t3.pdep.hybrid._read_qm_artifact's confinement) and classified unusable
        # here -- never usable, and never a candidate for vendoring that outside file.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        outside_log = tmp_path / 'outside.out'
        outside_log.write_text('stub quantum chemistry log\n')
        record = _queued_record(arc_project_directory=arc_dir)
        _write_artifact(record.expected_artifact_path, log_path=str(outside_log))
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert artifacts[0].status == ARTIFACT_STATUS_UNUSABLE
        assert 'outside' in artifacts[0].reason.lower()

    def test_wrong_shape_artifact_is_unusable(self, tmp_path):
        # A file that merely parses as Python (and resolves its Log(...) calls, of which it has
        # none) is NOT enough: an Arkane-DSL-shaped `transitionState(label=..., E0=...)` file is
        # not the ARC statmech species-file shape and must be classified unusable, with a reason
        # naming the missing module-level assignments.
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        os.makedirs(os.path.dirname(record.expected_artifact_path), exist_ok=True)
        with open(record.expected_artifact_path, 'w') as f:
            f.write("transitionState(label='TS0', E0=(100.0, 'kJ/mol'))\n")
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert artifacts[0].status == ARTIFACT_STATUS_UNUSABLE
        assert artifacts[0].artifact_path is None
        for name in ('energy', 'geometry', 'frequencies'):
            assert name in artifacts[0].reason

    def test_partially_shaped_artifact_is_unusable(self, tmp_path):
        # Every one of the load-bearing assignments (energy, geometry, frequencies) is mandatory;
        # an artifact carrying only some of them is unusable, naming exactly the missing ones.
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        os.makedirs(os.path.dirname(record.expected_artifact_path), exist_ok=True)
        log_path = os.path.join(os.path.dirname(record.expected_artifact_path), 'output.out')
        with open(log_path, 'w') as f:
            f.write('stub quantum chemistry log\n')
        with open(record.expected_artifact_path, 'w') as f:
            f.write("linear = False\n\nspinMultiplicity = 2\n\nenergy = Log('output.out')\n")
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert artifacts[0].status == ARTIFACT_STATUS_UNUSABLE
        assert 'geometry' in artifacts[0].reason
        assert 'frequencies' in artifacts[0].reason

    def test_unrecognized_join_status_raises(self, tmp_path):
        # An unrecognized join status must raise loudly, never fall through into the queued
        # branch and get a recomputed artifact path (an `assert` would vanish under `python -O`).
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        record.status = 'weird_status'  # bypass TSJoinRecord.__init__'s own validation
        with pytest.raises(ValueError, match='weird_status'):
            discover_ts_artifacts([record], arc_dir)

    def test_stored_path_disagreement_recomputed_wins(self, tmp_path):
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        # Corrupt the stored path (as if a sidecar had been hand-edited); the artifact is only at
        # the recomputed location.
        record.expected_artifact_path = os.path.join(arc_dir, 'somewhere', 'else.py')
        _write_artifact(expected_ts_artifact_path(arc_dir, record.arc_ts_label))
        # Finding 18 requires convergence to be explicitly True for USABLE; supply that here so
        # this test still exercises USABLE while keeping its original path-disagreement intent.
        output_dir = tmp_path / 'output'
        output_dir.mkdir()
        (output_dir / 'status.yml').write_text(
            f"{record.arc_ts_label}:\n"
            "  convergence: true\n"
            "  job_types: {}\n"
            "  paths: {}\n"
            "  info: ''\n"
            "  errors: ''\n"
        )
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert artifacts[0].status == ARTIFACT_STATUS_USABLE
        assert artifacts[0].artifact_path == expected_ts_artifact_path(arc_dir, record.arc_ts_label)

    def test_stored_label_disagreement_is_refused(self, tmp_path):
        # Finding 20: unlike the stored-path disagreement (merely advisory -- the recomputed path
        # is used regardless), a stored `arc_ts_label` that disagrees with the recomputed one
        # means the join record's identity itself cannot be trusted (a hand-edited or stale
        # sidecar), so it must be refused (ARTIFACT_STATUS_UNUSABLE), not silently normalized to
        # the recomputed label as if nothing were wrong -- even though the artifact for the
        # RECOMPUTED label is present and would otherwise be usable.
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        recomputed_label = arc_ts_label(record.network_id, record.network_ts_label)
        _write_artifact(expected_ts_artifact_path(arc_dir, recomputed_label))
        record.arc_ts_label = 'SomeOtherLabelEntirely'
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert len(artifacts) == 1
        assert artifacts[0].status == ARTIFACT_STATUS_UNUSABLE
        assert 'SomeOtherLabelEntirely' in artifacts[0].reason
        assert recomputed_label in artifacts[0].reason

    def test_path_escaping_project_directory_is_refused(self, tmp_path, monkeypatch):
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        outside = tmp_path.parent / 'outside_the_project.py'
        with open(outside, 'w') as f:
            f.write("transitionState(label='TS0', E0=(100.0, 'kJ/mol'))\n")

        import t3.pdep.discovery as discovery_module

        def _escaping_path(_arc_project_directory, _label):
            return str(outside)

        monkeypatch.setattr(discovery_module, 'expected_ts_artifact_path', _escaping_path)
        artifacts = discover_ts_artifacts([record], arc_dir)
        assert artifacts[0].status == ARTIFACT_STATUS_UNUSABLE
        assert artifacts[0].artifact_path is None
        assert 'outside' in artifacts[0].reason.lower() or 'refus' in artifacts[0].reason.lower()

    def test_sensitivity_evidence_attached(self, tmp_path):
        arc_dir = str(tmp_path)
        record = _queued_record(arc_project_directory=arc_dir)
        sensitivity = SensitiveTransitionState(
            ts_label='TS3', coefficient=0.5, condition=(1000.0, 'K', 1.0, 'bar'),
            path_reaction_label='reaction1', path_reaction_str='A + B <=> C',
            kinetics_comment='', uncertain=True, delta_ln_k=0.75,
        )
        artifacts = discover_ts_artifacts([record], arc_dir,
                                          sensitivity_by_ts={('network4_1', 'TS3'): sensitivity})
        assert artifacts[0].coefficient == 0.5
        assert artifacts[0].delta_ln_k == 0.75

    def test_sensitivity_keyed_by_network_never_cross_contaminates(self, tmp_path):
        # Two different networks can both have a 'TS1'. The sensitivity mapping is keyed by
        # (network_id, network_ts_label) -- the same key TSJoinRecord.key uses -- so
        # network4_1/TS1 must never inherit network4_2/TS1's sensitivity, and vice versa.
        arc_dir = str(tmp_path)
        record_1 = _queued_record(network_id='network4_1', network_ts_label='TS1', arc_project_directory=arc_dir)
        record_2 = _queued_record(network_id='network4_2', network_ts_label='TS1', arc_project_directory=arc_dir)
        sensitivity = {
            ('network4_1', 'TS1'): SensitiveTransitionState(
                ts_label='TS1', coefficient=0.5, condition=(1000.0, 'K', 1.0, 'bar'),
                path_reaction_label='reaction1', path_reaction_str='A + B <=> C',
                kinetics_comment='', uncertain=True, delta_ln_k=0.75,
            ),
            ('network4_2', 'TS1'): SensitiveTransitionState(
                ts_label='TS1', coefficient=-0.2, condition=(1000.0, 'K', 1.0, 'bar'),
                path_reaction_label='reaction9', path_reaction_str='D + E <=> F',
                kinetics_comment='', uncertain=True, delta_ln_k=0.1,
            ),
        }
        artifacts = discover_ts_artifacts([record_1, record_2], arc_dir, sensitivity_by_ts=sensitivity)
        assert artifacts[0].coefficient == 0.5
        assert artifacts[0].delta_ln_k == 0.75
        assert artifacts[1].coefficient == -0.2
        assert artifacts[1].delta_ln_k == 0.1


class TestEvaluatePdepHybrid:

    def _artifact(self, status, network_ts_label='TS3', delta_ln_k=None, network_id='network4_1'):
        label = arc_ts_label(network_id, network_ts_label)
        return TSArtifactRecord(
            network_id=network_id,
            network_ts_label=network_ts_label,
            arc_ts_label=label,
            status=status,
            # A usable record always carries the path of its usable artifact; anything else never does.
            artifact_path=f'/arc/calcs/statmech/kinetics/TSs/{label}.py' if status == ARTIFACT_STATUS_USABLE else None,
            converged=None,
            reason='',
            delta_ln_k=delta_ln_k,
        )

    def test_dominance_flag_none_when_any_selected_ts_is_unscored(self):
        # If even one selected transition state carries no sensitivity value, the unscored one
        # could itself be the dominant one -- so the flag must be UNKNOWN (None), never a quiet
        # False computed from the scored subset alone.
        artifacts = [
            self._artifact(ARTIFACT_STATUS_USABLE, 'TS3', delta_ln_k=0.2),
            self._artifact(ARTIFACT_STATUS_MISSING, 'TS4', delta_ln_k=None),
        ]
        evaluation = evaluate_pdep_hybrid(artifacts, allow_partial_selected_qm=True)
        assert evaluation.dominant_sensitivity_ts_missing is None

    def test_dominance_flag_none_for_singleton_selected_set(self):
        # With a single selected transition state, "the most sensitive one" is meaningless:
        # the flag must be None, not a vacuous True/False.
        artifacts = [self._artifact(ARTIFACT_STATUS_MISSING, 'TS3', delta_ln_k=0.9)]
        evaluation = evaluate_pdep_hybrid(artifacts)
        assert evaluation.dominant_sensitivity_ts_missing is None

    def test_refuses_records_spanning_more_than_one_network(self):
        # An evaluation describes ONE network; mixing records from two networks must be refused
        # loudly rather than silently folded into one wrong verdict.
        artifacts = [self._artifact(ARTIFACT_STATUS_USABLE, network_id='network4_1'),
                     self._artifact(ARTIFACT_STATUS_USABLE, network_id='network4_2')]
        with pytest.raises(ValueError, match='network4_1.*network4_2|network4_2.*network4_1'):
            evaluate_pdep_hybrid(artifacts)

    def test_refuses_duplicate_network_ts_records(self):
        # The same (network_id, network_ts_label) appearing twice would double-count that
        # transition state in the coverage counts; refuse rather than mis-count.
        artifacts = [self._artifact(ARTIFACT_STATUS_USABLE), self._artifact(ARTIFACT_STATUS_USABLE)]
        with pytest.raises(ValueError, match='TS3'):
            evaluate_pdep_hybrid(artifacts)

    def test_refuses_usable_record_without_artifact_path(self):
        # A record claiming 'usable' while carrying no artifact path is internally inconsistent:
        # nothing downstream could actually consume it. Refuse rather than count it as usable.
        artifact = self._artifact(ARTIFACT_STATUS_USABLE)
        artifact.artifact_path = None
        with pytest.raises(ValueError, match='artifact_path|artifact path'):
            evaluate_pdep_hybrid([artifact])

    def test_all_usable_is_complete_and_accepted(self):
        artifacts = [self._artifact(ARTIFACT_STATUS_USABLE), self._artifact(ARTIFACT_STATUS_USABLE, 'TS4')]
        evaluation = evaluate_pdep_hybrid(artifacts)
        assert evaluation.hybrid_status == HYBRID_STATUS_COMPLETE
        assert evaluation.accepted is True

    def test_partial_rejected_by_default(self):
        artifacts = [self._artifact(ARTIFACT_STATUS_USABLE), self._artifact(ARTIFACT_STATUS_MISSING, 'TS4')]
        evaluation = evaluate_pdep_hybrid(artifacts)
        assert evaluation.hybrid_status == HYBRID_STATUS_PARTIAL_SELECTED_QM
        assert evaluation.accepted is False

    def test_partial_accepted_when_allowed(self):
        artifacts = [self._artifact(ARTIFACT_STATUS_USABLE), self._artifact(ARTIFACT_STATUS_MISSING, 'TS4')]
        evaluation = evaluate_pdep_hybrid(artifacts, allow_partial_selected_qm=True)
        assert evaluation.hybrid_status == HYBRID_STATUS_PARTIAL_SELECTED_QM
        assert evaluation.accepted is True

    def test_zero_usable_never_accepted_even_when_partial_allowed(self):
        artifacts = [self._artifact(ARTIFACT_STATUS_MISSING), self._artifact(ARTIFACT_STATUS_UNUSABLE, 'TS4')]
        evaluation = evaluate_pdep_hybrid(artifacts, allow_partial_selected_qm=True)
        assert evaluation.hybrid_status == HYBRID_STATUS_NOT_EVALUATED
        assert evaluation.accepted is False

    def test_dominant_sensitivity_missing_flag_set(self):
        artifacts = [
            self._artifact(ARTIFACT_STATUS_USABLE, 'TS3', delta_ln_k=0.2),
            self._artifact(ARTIFACT_STATUS_MISSING, 'TS4', delta_ln_k=0.9),
        ]
        evaluation = evaluate_pdep_hybrid(artifacts, allow_partial_selected_qm=True)
        assert evaluation.dominant_sensitivity_ts_missing is True

    def test_dominant_sensitivity_missing_flag_false_when_dominant_is_usable(self):
        artifacts = [
            self._artifact(ARTIFACT_STATUS_USABLE, 'TS3', delta_ln_k=0.9),
            self._artifact(ARTIFACT_STATUS_MISSING, 'TS4', delta_ln_k=0.2),
        ]
        evaluation = evaluate_pdep_hybrid(artifacts, allow_partial_selected_qm=True)
        assert evaluation.dominant_sensitivity_ts_missing is False

    def test_dominant_sensitivity_flag_none_without_sensitivity_data(self):
        artifacts = [self._artifact(ARTIFACT_STATUS_USABLE), self._artifact(ARTIFACT_STATUS_MISSING, 'TS4')]
        evaluation = evaluate_pdep_hybrid(artifacts, allow_partial_selected_qm=True)
        assert evaluation.dominant_sensitivity_ts_missing is None
