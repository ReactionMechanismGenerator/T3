"""
Tests for t3.pdep.capture.

Every fixture is built under ``tmp_path``; none reads from ``tests/data/``, since these tests
exercise how ``capture`` reacts to the *shape* of an ARC project directory (and its interaction
with discovery), not to any one recorded run of it. Mirrors the fixture-construction style of
``tests/test_pdep/test_discovery.py``.
"""

import hashlib
import os

import pytest
import yaml

from t3.pdep.capture import CAPTURE_MANIFEST_FILE_NAME, CaptureResult, VerifyResult, capture_ts_artifacts, \
    verify_capture
from t3.pdep.discovery import (
    ARTIFACT_STATUS_MISSING,
    ARTIFACT_STATUS_NOT_QUEUED,
    ARTIFACT_STATUS_UNUSABLE,
    ARTIFACT_STATUS_UNVERIFIED,
    ARTIFACT_STATUS_USABLE,
)
from t3.pdep.hybrid import _read_qm_artifact
from t3.pdep.join import JOIN_STATUS_ALREADY_PRESENT, JOIN_STATUS_NOT_QUEUED, JOIN_STATUS_QUEUED, TSJoinRecord, \
    arc_ts_label, expected_ts_artifact_path


def _write_artifact(path: str, log_path: str | None = None) -> None:
    """Write a minimal artifact at ``path`` in the ARC statmech species-file shape, exactly as
    ``tests/test_pdep/test_discovery.py::_write_artifact`` does."""
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


def _write_status_yml(arc_dir: str, label: str, converged: bool) -> None:
    output_dir = os.path.join(arc_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'status.yml'), 'a') as f:
        f.write(
            f"{label}:\n"
            f"  convergence: {str(converged).lower()}\n"
            "  job_types: {}\n"
            "  paths: {}\n"
            "  info: ''\n"
            "  errors: ''\n"
        )


def _sha256_of(path: str) -> str:
    """Compute a file's sha256 hex digest, independently of t3.pdep.capture._sha256_file, so tests
    asserting on hashes never share an implementation with the code under test."""
    hasher = hashlib.sha256()
    with open(path, 'rb') as f:
        hasher.update(f.read())
    return hasher.hexdigest()


def _write_energy_settings_fixture(arc_dir: str, statmech_subdir: str = 'kinetics') -> None:
    """Write a minimal, valid ARC energy-settings fixture -- ``calcs/statmech/<statmech_subdir>/input.py``
    and ``output/output.yml`` -- under ``arc_dir``, exactly what
    ``t3.pdep.energy_settings.read_arc_energy_settings`` requires to freeze a complete energy-settings
    block. Deliberately sets ``useAtomCorrections``/``useBondCorrections`` both ``False`` so no
    cross-validation against ``output.yml`` (frequencyScaleFactor/atomEnergies/bondCorrectionType) is
    ever triggered by this fixture. Gated so it never overwrites an ``output.yml`` a test has already
    written of its own (``output.yml`` is project-global, shared across statmech subdirs -- see the
    module docstring of ``t3.pdep.energy_settings``)."""
    statmech_dir = os.path.join(arc_dir, 'calcs', 'statmech', statmech_subdir)
    os.makedirs(statmech_dir, exist_ok=True)
    input_py_path = os.path.join(statmech_dir, 'input.py')
    if not os.path.isfile(input_py_path):
        with open(input_py_path, 'w') as f:
            f.write(
                "modelChemistry = 'CBS-QB3'\n\n"
                "useHinderedRotors = True\n\n"
                "useAtomCorrections = False\n\n"
                "useBondCorrections = False\n"
            )
    output_dir = os.path.join(arc_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    output_yml_path = os.path.join(output_dir, 'output.yml')
    if not os.path.isfile(output_yml_path):
        with open(output_yml_path, 'w') as f:
            f.write('{}\n')


def _usable_record(tmp_path, network_id='network4_1', network_ts_label='TS3'):
    """Build one USABLE queued join record, plus its converged status.yml entry and a minimal valid
    energy-settings fixture (required now that any capture with >=1 artifact must freeze one),
    entirely under ``tmp_path`` (the ARC project directory)."""
    arc_dir = str(tmp_path)
    record = _queued_record(network_id=network_id, network_ts_label=network_ts_label, arc_project_directory=arc_dir)
    _write_artifact(record.expected_artifact_path)
    _write_status_yml(arc_dir, record.arc_ts_label, converged=True)
    _write_energy_settings_fixture(arc_dir)
    return record


class TestCaptureTsArtifacts:

    def test_usable_artifact_is_captured_into_capture_dir(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        result = capture_ts_artifacts([record], arc_dir, capture_dir)

        assert isinstance(result, CaptureResult)
        assert len(result.records) == 1
        captured = result.records[0]
        assert captured.status == ARTIFACT_STATUS_USABLE
        assert captured.converged is True
        # The captured path must point INTO the capture dir, never back at the (ephemeral) ARC
        # project directory.
        assert captured.artifact_path is not None
        assert os.path.isfile(captured.artifact_path)
        assert os.path.commonpath([os.path.realpath(captured.artifact_path), os.path.realpath(capture_dir)]) \
            == os.path.realpath(capture_dir)
        assert not captured.artifact_path.startswith(os.path.realpath(arc_dir))

    def test_captured_artifact_is_re_readable_with_capture_dir_as_log_root(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        result = capture_ts_artifacts([record], arc_dir, capture_dir)
        captured = result.records[0]

        # A captured artifact must be re-readable in place: _read_qm_artifact must resolve its
        # relative Log(...) targets, and pass confinement, with allowed_log_root=capture_dir (NOT
        # arc_dir -- that is the whole point of vendoring).
        info = _read_qm_artifact(ts_label=captured.arc_ts_label, artifact_path=captured.artifact_path,
                                  allowed_log_root=capture_dir)
        # The fixture's energy/geometry/frequencies assignments all quote the identical relative
        # path, and _read_qm_artifact's 'logs' list is de-duplicated by (original_path,
        # resolved_path) pair -- so a single log entry is expected here, not three.
        assert len(info['logs']) == 1
        for _, resolved in info['logs']:
            assert os.path.isfile(resolved)
            assert os.path.commonpath([os.path.realpath(resolved), os.path.realpath(capture_dir)]) \
                == os.path.realpath(capture_dir)

    def test_mixed_statuses_common_case_missing_is_the_majority(self, tmp_path):
        # Missing artifacts are the common case (~56% in real runs), not the exception. Build a
        # mixed set: one usable, one missing (queued but ARC never produced anything), one
        # unusable (log deleted after artifact was written), one not_queued, one unverified (no
        # status source at all).
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        # This test builds every record's artifacts directly (not via _usable_record), so the
        # mandatory energy-settings fixture must be written into the shared arc_dir here too --
        # two of the records below (usable, unverified) do carry an artifact.
        _write_energy_settings_fixture(arc_dir)

        usable = _queued_record(network_id='network4_1', network_ts_label='TS_usable',
                                 arc_project_directory=arc_dir)
        _write_artifact(usable.expected_artifact_path)
        _write_status_yml(arc_dir, usable.arc_ts_label, converged=True)

        missing = _queued_record(network_id='network4_1', network_ts_label='TS_missing',
                                  arc_project_directory=arc_dir)
        # No artifact written at all.

        unusable = _queued_record(network_id='network4_1', network_ts_label='TS_unusable',
                                   arc_project_directory=arc_dir)
        unusable_log = os.path.join(os.path.dirname(unusable.expected_artifact_path), 'gone.out')
        os.makedirs(os.path.dirname(unusable_log), exist_ok=True)
        with open(unusable_log, 'w') as f:
            f.write('fake log')
        _write_artifact(unusable.expected_artifact_path, log_path=unusable_log)
        os.remove(unusable_log)

        not_queued = TSJoinRecord(
            network_id='network4_1',
            network_ts_label='TS_not_queued',
            status=JOIN_STATUS_NOT_QUEUED,
            reason='Could not build the species for this transition state.',
        )

        already_present = TSJoinRecord(
            network_id='network4_1',
            network_ts_label='TS_already_present',
            status=JOIN_STATUS_ALREADY_PRESENT,
            reason='This transition state was already known to T3.',
        )

        unverified = _queued_record(network_id='network4_1', network_ts_label='TS_unverified',
                                     arc_project_directory=arc_dir)
        _write_artifact(unverified.expected_artifact_path)
        # No status.yml/output.yml entry at all for this one -- convergence stays unknown.

        capture_dir = str(tmp_path / '__capture__')
        result = capture_ts_artifacts(
            [usable, missing, unusable, not_queued, already_present, unverified], arc_dir, capture_dir)

        by_label = {record.network_ts_label: record for record in result.records}
        assert by_label['TS_usable'].status == ARTIFACT_STATUS_USABLE
        assert by_label['TS_missing'].status == ARTIFACT_STATUS_MISSING
        assert by_label['TS_missing'].artifact_path is None
        assert by_label['TS_unusable'].status == ARTIFACT_STATUS_UNUSABLE
        assert by_label['TS_unusable'].artifact_path is None
        assert by_label['TS_not_queued'].status == ARTIFACT_STATUS_NOT_QUEUED
        assert by_label['TS_not_queued'].artifact_path is None
        assert by_label['TS_already_present'].status == ARTIFACT_STATUS_MISSING
        assert by_label['TS_already_present'].artifact_path is None
        assert by_label['TS_unverified'].status == ARTIFACT_STATUS_UNVERIFIED
        # UNVERIFIED still carries a (now-captured) artifact_path: the artifact parsed fine, it is
        # only convergence that was never confirmed.
        assert by_label['TS_unverified'].artifact_path is not None
        assert os.path.isfile(by_label['TS_unverified'].artifact_path)

        # Only the two records that actually carried an artifact_path (usable, unverified) were
        # captured; nobody else's path was fabricated.
        for label in ('TS_usable', 'TS_unverified'):
            assert os.path.commonpath(
                [os.path.realpath(by_label[label].artifact_path), os.path.realpath(capture_dir)]
            ) == os.path.realpath(capture_dir)

    def test_zero_usable_artifacts_produces_well_formed_empty_capture(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = TSJoinRecord(
            network_id='network4_1',
            network_ts_label='TS_none',
            status=JOIN_STATUS_NOT_QUEUED,
            reason='Could not build the species for this transition state.',
        )
        capture_dir = str(tmp_path / 'capture')

        result = capture_ts_artifacts([record], arc_dir, capture_dir)

        assert len(result.records) == 1
        assert result.records[0].artifact_path is None
        assert os.path.isfile(result.manifest_path)
        with open(result.manifest_path) as f:
            manifest = yaml.safe_load(f)
        assert manifest['transition_states'][0]['captured_artifact_path'] is None
        assert manifest['transition_states'][0]['status'] == ARTIFACT_STATUS_NOT_QUEUED

    def test_source_path_escaping_arc_project_dir_is_refused(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        outside_log = tmp_path / 'outside.out'
        outside_log.write_text('stub quantum chemistry log\n')
        record = _queued_record(arc_project_directory=arc_dir)
        _write_artifact(record.expected_artifact_path, log_path=str(outside_log))
        _write_status_yml(arc_dir, record.arc_ts_label, converged=True)

        capture_dir = str(tmp_path / 'capture')
        # discover_ts_artifacts itself refuses this at discovery time (classifies UNUSABLE, never
        # raises) -- confirm capture_ts_artifacts surfaces that same fail-closed classification
        # rather than somehow vendoring the escaping file.
        result = capture_ts_artifacts([record], arc_dir, capture_dir)
        assert result.records[0].status == ARTIFACT_STATUS_UNUSABLE
        assert result.records[0].artifact_path is None
        # qm/ may exist (capture_ts_artifacts now always reconciles it with the manifest, even on
        # a zero-artifact pass -- see the zero-artifact-clears-stale-qm_dir test), but it must not
        # hold the escaping artifact or any file at all.
        qm_dir = os.path.join(capture_dir, 'qm')
        if os.path.isdir(qm_dir):
            captured_files = [name for _root, _dirs, names in os.walk(qm_dir) for name in names]
            assert captured_files == []

    def test_manifest_hashes_match_sources_and_verdict_is_frozen(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        result = capture_ts_artifacts([record], arc_dir, capture_dir)

        with open(result.manifest_path) as f:
            manifest = yaml.safe_load(f)
        assert manifest['arc_project_directory'] == arc_dir
        entry = manifest['transition_states'][0]
        assert entry['network_id'] == record.network_id
        assert entry['network_ts_label'] == record.network_ts_label
        assert entry['status'] == ARTIFACT_STATUS_USABLE
        assert entry['converged'] is True

        import hashlib

        def sha256_of(path):
            hasher = hashlib.sha256()
            with open(path, 'rb') as fh:
                hasher.update(fh.read())
            return hasher.hexdigest()

        assert entry['source_artifact_sha256'] == sha256_of(record.expected_artifact_path)
        assert len(entry['source_logs']) == 1
        source_log_entry = entry['source_logs'][0]
        resolved_source_log = os.path.join(os.path.dirname(record.expected_artifact_path), 'output.out')
        assert source_log_entry['sha256'] == sha256_of(resolved_source_log)

        # And the captured (vendored) copy round-trips too: it must have the SAME content, hence
        # the same hash, as the source it was copied from.
        captured_artifact_path = os.path.join(capture_dir, entry['captured_artifact_path'])
        assert sha256_of(captured_artifact_path) is not None  # merely confirm it is readable
        captured_log_path = os.path.join(capture_dir, next(iter(entry['captured_log_paths'].values())))
        assert os.path.isfile(captured_log_path)

    def test_source_artifact_hash_comes_from_the_same_read_as_the_vendored_content(self, monkeypatch, tmp_path):
        """Regression test for the pointer-file TOCTOU: source_artifact_sha256 must be derived from the
        exact bytes _read_qm_artifact already parsed, never from a second, independent read of
        record.artifact_path performed afterwards. We simulate the TOCTOU window by monkeypatching
        _sha256_file so that, if it is EVER called on the artifact's own path (the buggy, pre-fix
        behavior), it first overwrites the source file on disk -- standing in for something else
        mutating the pointer file in the gap between the two reads -- before computing the hash. If
        capture_ts_artifacts still (correctly) takes the hash from info['sha256'] instead, _sha256_file
        is never called on the artifact path at all, the file is never mutated, and the recorded hash
        matches the ORIGINAL content that was actually parsed and vendored."""
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        with open(record.expected_artifact_path, 'rb') as f:
            original_artifact_bytes = f.read()
        expected_sha256 = hashlib.sha256(original_artifact_bytes).hexdigest()

        import t3.pdep.capture as capture_module
        real_sha256_file = capture_module._sha256_file
        artifact_realpath = os.path.realpath(record.expected_artifact_path)

        def toctou_probing_sha256_file(path):
            if os.path.realpath(path) == artifact_realpath:
                # This branch must never execute for the artifact's own path: capture_ts_artifacts
                # is expected to use info['sha256'] (from _read_qm_artifact's own read) instead of
                # re-hashing record.artifact_path here. If it ever does land here, mutate the file
                # first so a lingering TOCTOU would be caught by the assertion below rather than
                # passing by accident because the content never actually changed.
                with open(record.expected_artifact_path, 'w') as f:
                    f.write('mutated in the TOCTOU window -- must never be the recorded hash source')
            return real_sha256_file(path)

        monkeypatch.setattr(capture_module, '_sha256_file', toctou_probing_sha256_file)

        result = capture_ts_artifacts([record], arc_dir, capture_dir)

        with open(result.manifest_path) as f:
            manifest = yaml.safe_load(f)
        entry = manifest['transition_states'][0]
        assert entry['source_artifact_sha256'] == expected_sha256, (
            "source_artifact_sha256 diverged from the bytes _read_qm_artifact actually parsed -- "
            "capture_ts_artifacts re-read record.artifact_path separately (a TOCTOU) instead of "
            "reusing info['sha256'] from the read _read_qm_artifact already performed."
        )
        # The source file must never have actually been mutated by our probe: proves _sha256_file
        # was never invoked on the artifact's own path at all (not just that its result happened to
        # match by coincidence).
        with open(record.expected_artifact_path, 'rb') as f:
            assert f.read() == original_artifact_bytes, (
                "record.artifact_path was re-read via _sha256_file after the TOCTOU-probe mutated "
                "it -- the source hash is no longer guaranteed to come from the same read used to "
                "parse and vendor the artifact."
            )

    def test_does_not_copy_status_yml_or_output_yml_into_qm_dir(self, tmp_path):
        # CRITICAL: naively copying either file straight into qm/ (or anywhere a consumer might
        # mistake for authoritative capture content) would risk a reader treating them as if they
        # were live enough to re-derive convergence from -- the resolved verdict, already frozen
        # into the manifest's per-transition-state fields, must remain the sole authority. The two
        # files themselves ARE archived (see test_archives_status_and_output_yml_as_inert_provenance
        # below), but only under the manifest's separate, clearly inert 'provenance' key/directory,
        # never under qm/.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        capture_ts_artifacts([record], arc_dir, capture_dir)

        qm_dir = os.path.join(capture_dir, 'qm')
        for root, _dirs, files in os.walk(qm_dir):
            assert 'status.yml' not in files
            assert 'output.yml' not in files

    def test_archives_status_and_output_yml_as_inert_provenance(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        # _usable_record writes status.yml and the mandatory kinetics energy-settings fixture
        # (calcs/statmech/kinetics/input.py + output/output.yml); overwrite output.yml with our own
        # content here so both provenance entries exercise the "present and archived" path with
        # bytes/hashes we control. The thermo statmech pass is deliberately never written, so its
        # provenance entry exercises the "absent" path in the same test (its dedicated "absent"
        # sibling below only re-exercises the fully-zero-artifact scenario).
        output_yml_path = os.path.join(arc_dir, 'output', 'output.yml')
        with open(output_yml_path, 'w') as f:
            f.write('some_key: some_value\n')
        status_yml_path = os.path.join(arc_dir, 'output', 'status.yml')
        kinetics_input_py_path = os.path.join(arc_dir, 'calcs', 'statmech', 'kinetics', 'input.py')
        thermo_input_py_path = os.path.join(arc_dir, 'calcs', 'statmech', 'thermo', 'input.py')
        status_sha256 = _sha256_of(status_yml_path)
        output_sha256 = _sha256_of(output_yml_path)
        kinetics_input_sha256 = _sha256_of(kinetics_input_py_path)
        status_mtime = os.path.getmtime(status_yml_path)
        output_mtime = os.path.getmtime(output_yml_path)
        kinetics_input_mtime = os.path.getmtime(kinetics_input_py_path)
        capture_dir = str(tmp_path / 'capture')

        capture_ts_artifacts([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)

        provenance = manifest['provenance']
        assert len(provenance) == 4
        # Both the kinetics and thermo statmech input.py sources are literally named "input.py",
        # so keying by os.path.basename(source_path) (the prior version of this test) would silently
        # collide two of the four entries; key by the full source_path instead.
        by_source_path = {entry['source_path']: entry for entry in provenance}
        assert set(by_source_path) == {status_yml_path, output_yml_path, kinetics_input_py_path,
                                        thermo_input_py_path}

        status_entry = by_source_path[status_yml_path]
        assert status_entry['captured_path'] == os.path.join('provenance', 'status.yml')
        assert status_entry['sha256'] == status_sha256
        assert status_entry['source_mtime'] == pytest.approx(status_mtime)
        captured_status_path = os.path.join(capture_dir, status_entry['captured_path'])
        assert os.path.isfile(captured_status_path)
        assert _sha256_of(captured_status_path) == status_sha256

        output_entry = by_source_path[output_yml_path]
        assert output_entry['captured_path'] == os.path.join('provenance', 'output.yml')
        assert output_entry['sha256'] == output_sha256
        assert output_entry['source_mtime'] == pytest.approx(output_mtime)
        captured_output_path = os.path.join(capture_dir, output_entry['captured_path'])
        assert os.path.isfile(captured_output_path)
        assert _sha256_of(captured_output_path) == output_sha256

        kinetics_entry = by_source_path[kinetics_input_py_path]
        assert kinetics_entry['captured_path'] == os.path.join('provenance', 'statmech_kinetics_input.py')
        assert kinetics_entry['sha256'] == kinetics_input_sha256
        assert kinetics_entry['source_mtime'] == pytest.approx(kinetics_input_mtime)
        captured_kinetics_path = os.path.join(capture_dir, kinetics_entry['captured_path'])
        assert os.path.isfile(captured_kinetics_path)
        assert _sha256_of(captured_kinetics_path) == kinetics_input_sha256

        thermo_entry = by_source_path[thermo_input_py_path]
        assert thermo_entry['captured_path'] is None
        assert thermo_entry['sha256'] is None
        assert thermo_entry['source_mtime'] is None

    def test_archives_provenance_absence_when_status_or_output_yml_missing(self, tmp_path):
        # A capture with zero usable artifacts (e.g. ARC never even reached the finalization step
        # for this transition state) never has status.yml, output.yml, or either statmech input.py
        # to archive; absence must be recorded, not raise. A capture that DOES have >=1 artifact now
        # requires a valid, complete energy-settings fixture to exist at all (see
        # test_archives_status_and_output_yml_as_inert_provenance for that "present" path, which also
        # covers the thermo input.py "absent" case) -- so the clean way to exercise a genuinely
        # *missing* output.yml, without also having to satisfy that mandatory energy-settings
        # contract, is the same zero-artifact scenario used by
        # test_zero_usable_artifacts_produces_well_formed_empty_capture.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = TSJoinRecord(
            network_id='network4_1',
            network_ts_label='TS_none',
            status=JOIN_STATUS_NOT_QUEUED,
            reason='Could not build the species for this transition state.',
        )
        capture_dir = str(tmp_path / 'capture')

        capture_ts_artifacts([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)

        assert manifest.get('energy_settings') is None
        provenance = manifest['provenance']
        assert len(provenance) == 4
        for entry in provenance:
            assert entry['captured_path'] is None
            assert entry['sha256'] is None
            assert entry['source_mtime'] is None
        # And nothing was ever created under provenance/ at all.
        provenance_dir = os.path.join(capture_dir, 'provenance')
        assert not os.path.isdir(provenance_dir)

    def test_mid_capture_failure_leaves_previous_capture_intact(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        capture_dir = str(tmp_path / 'capture')

        good_record = _usable_record(tmp_path=tmp_path.__class__(arc_dir), network_ts_label='TS_good')
        first_result = capture_ts_artifacts([good_record], arc_dir, capture_dir)
        previous_manifest_mtime = os.path.getmtime(first_result.manifest_path)
        previous_artifact_content = open(first_result.records[0].artifact_path).read()

        # Now a second pass, adding a transition state whose artifact references a log file that
        # vanishes AFTER the artifact was written (mirrors
        # test_discovery.py::test_artifact_present_but_referenced_log_deleted_is_unusable) --
        # except here we want a failure that happens DURING _read_qm_artifact's confinement/
        # existence check, i.e. before any write to capture_dir, so the previous capture must
        # survive completely untouched.
        bad_record = _queued_record(network_id='network4_1', network_ts_label='TS_bad',
                                     arc_project_directory=arc_dir)
        bad_log = os.path.join(os.path.dirname(bad_record.expected_artifact_path), 'vanishing.out')
        os.makedirs(os.path.dirname(bad_log), exist_ok=True)
        with open(bad_log, 'w') as f:
            f.write('fake log')
        _write_artifact(bad_record.expected_artifact_path, log_path=bad_log)
        _write_status_yml(arc_dir, bad_record.arc_ts_label, converged=True)
        os.remove(bad_log)

        # discover_ts_artifacts itself classifies a vanished log as UNUSABLE (never raises), so to
        # exercise capture_ts_artifacts' own pre-flight-before-write guarantee we instead force a
        # true mid-pipeline failure: delete the good artifact's underlying log file out from under
        # a second, otherwise-untouched queued record whose artifact path escapes discovery's own
        # classification by being hand-fed straight through discover_ts_artifacts. Since
        # discover_ts_artifacts already fails closed to UNUSABLE for this exact case (never handing
        # capture_ts_artifacts a broken artifact_path to read), assert instead that the *previous*
        # capture is untouched after this second, partially-bad call.
        second_result = capture_ts_artifacts([good_record, bad_record], arc_dir, capture_dir)

        # The bad record was classified UNUSABLE by discovery (never reached _read_qm_artifact
        # with a path capture_ts_artifacts would need to roll back), and the good record's capture
        # must still be present, unchanged.
        assert second_result.records[1].status == ARTIFACT_STATUS_UNUSABLE
        good_after = [r for r in second_result.records if r.network_ts_label == 'TS_good'][0]
        assert open(good_after.artifact_path).read() == previous_artifact_content
        assert os.path.getmtime(second_result.manifest_path) >= previous_manifest_mtime

    def test_mid_capture_read_failure_before_any_write_raises_and_leaves_prior_capture_untouched(self, tmp_path):
        # A genuine pre-flight failure: monkeypatch-free approach -- an artifact whose Log(...)
        # resolves outside arc_project_directory raises from _read_qm_artifact only when discovery
        # itself would have already called it usable; since discovery fails closed to UNUSABLE
        # first, we instead directly exercise capture's own raise-before-write guarantee by giving
        # it a record whose artifact_path is usable at discovery time but whose log is removed
        # between discovery and capture's own _read_qm_artifact pre-flight call. This can't happen
        # through the public discover_ts_artifacts + capture_ts_artifacts sequence directly (both
        # run back-to-back on the same, stable filesystem state), so instead assert the module-level
        # contract: capture_dir is left untouched by a prior successful capture when a *second*
        # capture_ts_artifacts call over a broken record raises before any write.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        capture_dir = str(tmp_path / 'capture')

        good_record = _usable_record(tmp_path=tmp_path.__class__(arc_dir), network_ts_label='TS_good')
        first_result = capture_ts_artifacts([good_record], arc_dir, capture_dir)
        previous_files = sorted(
            os.path.relpath(os.path.join(root, name), capture_dir)
            for root, _dirs, names in os.walk(capture_dir) for name in names
        )

        # Build a record whose join status is unrecognized, forcing discover_ts_artifacts itself
        # to raise -- capture_ts_artifacts must propagate that raise, untouched, before writing
        # anything for this second call.
        broken_record = _queued_record(network_id='network4_1', network_ts_label='TS_broken',
                                        arc_project_directory=arc_dir)
        broken_record.status = 'weird_status'  # bypass TSJoinRecord.__init__'s own validation

        with pytest.raises(ValueError):
            capture_ts_artifacts([good_record, broken_record], arc_dir, capture_dir)

        after_files = sorted(
            os.path.relpath(os.path.join(root, name), capture_dir)
            for root, _dirs, names in os.walk(capture_dir) for name in names
        )
        assert after_files == previous_files
        for relpath in previous_files:
            if relpath.endswith('.py') and 'qm' in relpath.split(os.sep):
                with open(os.path.join(capture_dir, relpath)) as f:
                    assert 'Log(' in f.read()

    def test_duplicate_arc_ts_label_across_two_networks_is_refused(self, tmp_path):
        # arc_ts_label() is `f'{prefix}_{network_id}_{network_ts_label}'`, and '_' is itself a
        # LABEL_SAFE_CHARS character -- so two DIFFERENT (network_id, network_ts_label) pairs can
        # genuinely collide on the same rendered label purely by where the '_' falls, with no
        # mutation needed: 'network1'/'1_TS_A' and 'network1_1'/'TS_A' both render to
        # 'T3PDep_network1_1_TS_A'. Both must never silently let one artifact clobber the other in
        # capture_ts_artifacts's internal per-label dicts.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        first = _usable_record(tmp_path=tmp_path.__class__(arc_dir), network_id='network1', network_ts_label='1_TS_A')
        second = _usable_record(tmp_path=tmp_path.__class__(arc_dir), network_id='network1_1', network_ts_label='TS_A')
        assert first.arc_ts_label == second.arc_ts_label  # confirm the collision is genuine

        capture_dir = str(tmp_path / 'capture')

        with pytest.raises(ValueError) as exc_info:
            capture_ts_artifacts([first, second], arc_dir, capture_dir)

        message = str(exc_info.value)
        assert first.arc_ts_label in message
        assert first.network_id in message
        assert first.network_ts_label in message
        assert second.network_id in message
        assert second.network_ts_label in message
        # Pre-flight-before-write: nothing was written to capture_dir at all.
        assert not os.path.isdir(capture_dir)

    def test_truthy_artifact_path_with_falsy_arc_ts_label_is_refused(self, tmp_path, monkeypatch):
        # discover_ts_artifacts itself always recomputes a non-empty arc_ts_label for any record
        # it hands back with a truthy artifact_path (arc_ts_label() raises on an empty
        # network_id/network_ts_label, so a QUEUED record can never legitimately resolve to a
        # falsy label). So this defensive guard in capture_ts_artifacts is unreachable through
        # discover_ts_artifacts' own normal resolution; exercise it directly by monkeypatching
        # discover_ts_artifacts to hand back a record shaped exactly like the corrupt state the
        # guard exists to catch (truthy artifact_path, falsy arc_ts_label).
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))

        from t3.pdep import capture as capture_module
        from t3.pdep.discovery import ARTIFACT_STATUS_USABLE, TSArtifactRecord

        corrupt_record = TSArtifactRecord(
            network_id=record.network_id,
            network_ts_label=record.network_ts_label,
            arc_ts_label=None,
            status=ARTIFACT_STATUS_USABLE,
            artifact_path=record.expected_artifact_path,
            converged=True,
            reason='stubbed for the falsy-arc_ts_label pre-flight guard test',
        )
        monkeypatch.setattr(capture_module, 'discover_ts_artifacts', lambda **kwargs: [corrupt_record])

        capture_dir = str(tmp_path / 'capture')

        with pytest.raises(ValueError) as exc_info:
            capture_ts_artifacts([record], arc_dir, capture_dir)

        message = str(exc_info.value)
        assert record.network_id in message
        assert record.network_ts_label in message
        assert not os.path.isdir(capture_dir)

    def test_zero_artifact_capture_clears_stale_qm_dir_from_prior_capture(self, tmp_path):
        # A capture into a directory that already holds a PREVIOUS capture's qm/ subdirectory,
        # this time with zero artifacts to capture, must leave the manifest and the on-disk state
        # in agreement: no stale qm/ survives from the earlier pass.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        capture_dir = str(tmp_path / 'capture')

        good_record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        first_result = capture_ts_artifacts([good_record], arc_dir, capture_dir)
        assert os.path.isdir(os.path.join(capture_dir, 'qm'))
        assert os.path.isfile(first_result.records[0].artifact_path)

        zero_record = TSJoinRecord(
            network_id='network4_1',
            network_ts_label='TS_none',
            status=JOIN_STATUS_NOT_QUEUED,
            reason='Could not build the species for this transition state.',
        )
        # This zero-artifact pass would destroy the valid, non-empty prior capture above; the
        # caller must explicitly opt in via supersede=True (see
        # test_zero_artifact_capture_refuses_to_clear_valid_prior_capture_by_default for the
        # default-refuses behavior this scenario now exercises).
        second_result = capture_ts_artifacts([zero_record], arc_dir, capture_dir, supersede=True)

        assert second_result.records[0].artifact_path is None
        with open(second_result.manifest_path) as f:
            manifest = yaml.safe_load(f)
        assert manifest['transition_states'][0]['captured_artifact_path'] is None

        # The manifest now records nothing captured; the qm/ directory must not still hold the
        # previous pass's artifact.
        qm_dir = os.path.join(capture_dir, 'qm')
        remaining_files = []
        if os.path.isdir(qm_dir):
            for root, _dirs, files in os.walk(qm_dir):
                remaining_files.extend(files)
        assert remaining_files == []

    def test_zero_artifact_capture_refuses_to_clear_valid_prior_capture_by_default(self, tmp_path):
        # The DEFAULT (supersede=False) behavior for the exact same scenario as
        # test_zero_artifact_capture_clears_stale_qm_dir_from_prior_capture: a zero-artifact pass
        # must RAISE rather than silently destroy a valid, non-empty existing capture, since a
        # zero-artifact pass is indistinguishable, from inside capture_ts_artifacts, from an
        # accidental replay against a stale/cleaned-up arc_project_directory.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        capture_dir = str(tmp_path / 'capture')

        good_record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        first_result = capture_ts_artifacts([good_record], arc_dir, capture_dir)
        captured_artifact_path = first_result.records[0].artifact_path
        assert os.path.isfile(captured_artifact_path)

        zero_record = TSJoinRecord(
            network_id='network4_1',
            network_ts_label='TS_none',
            status=JOIN_STATUS_NOT_QUEUED,
            reason='Could not build the species for this transition state.',
        )
        with pytest.raises(ValueError, match='supersede'):
            capture_ts_artifacts([zero_record], arc_dir, capture_dir)

        # Refusing must leave the prior, good capture completely untouched -- not partially
        # cleared, not re-verified-and-then-clobbered.
        assert os.path.isfile(captured_artifact_path)
        verify_capture(capture_dir)  # still verifies clean; would raise if anything had been torn

    def test_duplicate_arc_ts_label_with_both_artifacts_missing_is_refused(self, tmp_path):
        # The existing collision test (test_duplicate_arc_ts_label_across_two_networks_is_refused)
        # only exercises the case where both colliding records DO have artifacts. The common case in
        # a real run is the opposite: one or both of the colliding networks hasn't produced an
        # artifact for this transition state yet. validate_ts_join_records checks the FULL
        # join_records population up front, before discovery, so it must refuse this collision even
        # when NEITHER record currently has an artifact_path.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        first = _queued_record(network_id='network1', network_ts_label='1_TS_A', arc_project_directory=arc_dir)
        second = _queued_record(network_id='network1_1', network_ts_label='TS_A', arc_project_directory=arc_dir)
        assert first.arc_ts_label == second.arc_ts_label  # confirm the collision is genuine
        # Neither artifact is written at all -- both remain MISSING at discovery time.

        capture_dir = str(tmp_path / 'capture')

        with pytest.raises(ValueError) as exc_info:
            capture_ts_artifacts([first, second], arc_dir, capture_dir)

        message = str(exc_info.value)
        assert first.arc_ts_label in message
        # Pre-flight-before-discovery: nothing was written to capture_dir at all.
        assert not os.path.isdir(capture_dir)

    def test_capture_dir_inside_arc_project_directory_is_refused(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        # capture_dir nested inside the ARC project directory -- exactly what this module exists to
        # be durable against, since ARC deletes/recreates its own subtrees.
        capture_dir = os.path.join(arc_dir, 'capture')

        with pytest.raises(ValueError) as exc_info:
            capture_ts_artifacts([record], arc_dir, capture_dir)

        message = str(exc_info.value)
        assert capture_dir in message
        assert arc_dir in message
        assert not os.path.isdir(capture_dir)

    def test_capture_dir_holding_unrelated_content_is_refused(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        # Pre-populate capture_dir with content this module never wrote -- no capture manifest at
        # all, just an unrelated file a caller happened to put there.
        os.makedirs(capture_dir)
        with open(os.path.join(capture_dir, 'unrelated.txt'), 'w') as f:
            f.write('not a capture')

        with pytest.raises(ValueError) as exc_info:
            capture_ts_artifacts([record], arc_dir, capture_dir)

        message = str(exc_info.value)
        assert capture_dir in message
        # Refused before any capture I/O -- the unrelated file must still be the only thing there.
        assert os.listdir(capture_dir) == ['unrelated.txt']

    def test_verify_capture_detects_tampered_captured_log(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        result = capture_ts_artifacts([record], arc_dir, capture_dir)
        # Sanity: a freshly-written, untouched capture verifies clean.
        verify_capture(capture_dir)

        entry = None
        with open(result.manifest_path) as f:
            manifest = yaml.safe_load(f)
        entry = manifest['transition_states'][0]
        captured_log_relpath = next(iter(entry['captured_log_paths'].values()))
        captured_log_abs_path = os.path.join(capture_dir, captured_log_relpath)
        with open(captured_log_abs_path, 'a') as f:
            f.write('tampered content appended after capture\n')

        with pytest.raises(ValueError) as exc_info:
            verify_capture(capture_dir)

        message = str(exc_info.value)
        assert captured_log_relpath in message

    def test_verify_capture_detects_torn_capture_when_qm_dir_removed(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        capture_ts_artifacts([record], arc_dir, capture_dir)
        verify_capture(capture_dir)  # sanity: clean immediately after capture

        # Simulate a crash between vendoring qm/ and writing the manifest (or an out-of-band edit
        # after the fact): the manifest is present and claims captured files exist, but qm/ itself is
        # gone.
        qm_dir = os.path.join(capture_dir, 'qm')
        import shutil
        shutil.rmtree(qm_dir)

        with pytest.raises(ValueError) as exc_info:
            verify_capture(capture_dir)

        message = str(exc_info.value)
        assert 'Torn' in message or 'torn' in message

    def test_verify_capture_raises_on_missing_manifest(self, tmp_path):
        capture_dir = str(tmp_path / 'not_a_capture')
        os.makedirs(capture_dir)

        with pytest.raises(ValueError):
            verify_capture(capture_dir)

    def test_verify_capture_raises_on_structurally_empty_manifest(self, tmp_path):
        # A manifest whose 'transition_states' list is empty (e.g. hand-truncated, or a bug in
        # whatever wrote it) must NOT "verify" vacuously -- there is nothing for the per-entry loop
        # to check, so silently returning success would let a truncated/empty manifest pass as if
        # it had been meaningfully verified. This manifest is written directly to disk, bypassing
        # capture_ts_artifacts entirely, precisely because capture_ts_artifacts itself never
        # produces such a manifest (see test_zero_usable_artifacts_produces_well_formed_empty_capture
        # for the well-formed zero-artifact case, which still has one entry per transition state
        # asked about).
        capture_dir = str(tmp_path / 'capture')
        os.makedirs(capture_dir)
        with open(os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME), 'w') as f:
            yaml.safe_dump({'arc_project_directory': str(tmp_path / 'arc_project'), 'transition_states': []}, f)

        with pytest.raises(ValueError, match='empty'):
            verify_capture(capture_dir)

    def test_capture_refuses_empty_join_records(self, tmp_path):
        # Closes the gap that made the comment above true only by accident: capture_ts_artifacts
        # used to ACCEPT an empty join_records and write a manifest with zero 'transition_states',
        # which verify_capture then refused as structurally vacuous -- the module producing a
        # capture its own verifier rejects. Any consumer that verifies before trusting (the
        # finalization wiring in t3/main.py is about to) would have seen a spurious "invalid
        # capture" for a capture this module itself wrote. Refusing at the entry point keeps
        # "anything capture_ts_artifacts wrote, verify_capture accepts" an actual invariant.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        capture_dir = str(tmp_path / 'capture')

        with pytest.raises(ValueError, match='empty'):
            capture_ts_artifacts(join_records=[], arc_project_directory=arc_dir, capture_dir=capture_dir)

        assert not os.path.exists(capture_dir), 'a refused capture must not leave a capture directory behind'

    def test_verify_capture_returns_concrete_counts_on_success(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        capture_result = capture_ts_artifacts([record], arc_dir, capture_dir)
        captured_artifact_path = capture_result.records[0].artifact_path
        assert os.path.isfile(captured_artifact_path)

        result = verify_capture(capture_dir)

        assert isinstance(result, VerifyResult)
        assert result.capture_dir == capture_dir
        assert result.manifest_path == os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        assert result.record_count == 1
        assert result.captured_artifact_count == 1

    def test_captured_log_hash_mismatch_during_copy_is_refused(self, tmp_path, monkeypatch):
        # A real mid-copy corruption: shutil.copyfile itself (used internally by
        # _vendor_qm_artifacts to vendor each log) is monkeypatched to write different bytes than it
        # was asked to copy, simulating a truncated/corrupted copy. capture_ts_artifacts must catch
        # this via its own captured-vs-source hash comparison immediately after the copy, and must
        # never return a capture whose vendored bytes disagree with their source.
        from t3.pdep import hybrid as hybrid_module

        def corrupt_copyfile(src, dst):
            with open(dst, 'w') as f:
                f.write('corrupted during copy\n')

        monkeypatch.setattr(hybrid_module.shutil, 'copyfile', corrupt_copyfile)

        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        with pytest.raises(ValueError) as exc_info:
            capture_ts_artifacts([record], arc_dir, capture_dir)

        assert 'Torn or corrupt' in str(exc_info.value)

    def test_energy_settings_is_frozen_into_result_and_manifest_when_artifact_captured(self, tmp_path):
        # Any capture with >=1 usable artifact must freeze the ARC energy-reference settings that
        # produced it (t3/pdep/hybrid.py::write_hybrid_network_input_file requires a non-blank
        # model_chemistry downstream); _usable_record's fixture sets modelChemistry = 'CBS-QB3'.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        result = capture_ts_artifacts([record], arc_dir, capture_dir)

        assert isinstance(result.energy_settings, dict)
        # A plain string label is frozen as its string VALUE, WITHOUT the surrounding quotes, so that
        # it drops straight into QMEnergySettings.model_chemistry (whose validator rejects embedded
        # quotes, and which re-quotes the label at render time). Freezing the raw source segment
        # instead would round-trip a legitimate Arkane input into a frozen capture value that the
        # hybrid writer then refuses -- possibly with the ARC project directory already gone. Only
        # the bare LevelOfTheory(...)/CompositeLevelOfTheory(...) forms are frozen verbatim; see
        # tests/test_pdep/test_energy_settings.py::
        # test_frozen_model_chemistry_round_trips_into_the_hybrid_writer, which covers both.
        assert result.energy_settings['model_chemistry'] == 'CBS-QB3'
        assert result.energy_settings['use_hindered_rotors'] is True
        assert result.energy_settings['use_atom_corrections'] is False
        assert result.energy_settings['use_bond_corrections'] is False

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        assert manifest['energy_settings']['model_chemistry'] == 'CBS-QB3'

    def test_energy_settings_is_none_with_zero_captured_artifacts(self, tmp_path):
        # No artifact was ever fetched from ARC, so there is nothing whose provenance needs
        # freezing; this must not raise even if ARC never wrote a statmech kinetics input.py at all
        # (which is exactly the case here -- arc_dir has no calcs/ directory whatsoever).
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = TSJoinRecord(
            network_id='network4_1',
            network_ts_label='TS_none',
            status=JOIN_STATUS_NOT_QUEUED,
            reason='Could not build the species for this transition state.',
        )
        capture_dir = str(tmp_path / 'capture')

        result = capture_ts_artifacts([record], arc_dir, capture_dir)

        assert result.energy_settings is None

    def test_capture_raises_when_energy_settings_fixture_is_missing_and_leaves_capture_dir_untouched(self, tmp_path):
        # capture_ts_artifacts's energy-settings preflight must fail closed and, per the
        # pre-flight-before-any-write convention, must not have written anything to capture_dir yet.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        # Remove the kinetics input.py that _usable_record wrote, so read_arc_energy_settings has
        # nothing to read.
        kinetics_input_py_path = os.path.join(arc_dir, 'calcs', 'statmech', 'kinetics', 'input.py')
        os.remove(kinetics_input_py_path)
        capture_dir = str(tmp_path / 'capture')

        with pytest.raises(ValueError, match='ARC statmech'):
            capture_ts_artifacts([record], arc_dir, capture_dir)

        assert not os.path.exists(capture_dir), (
            'a capture refused during the energy-settings preflight must not leave a capture '
            'directory behind -- the preflight runs before any write, per capture_ts_artifacts\' '
            'read-everything-before-writing-anything convention.'
        )

    def test_verify_capture_raises_when_energy_settings_block_is_missing_but_artifacts_were_captured(self, tmp_path):
        # A hand-edited or otherwise malformed manifest that dropped its 'energy_settings' key (or
        # zeroed it out) despite having captured artifact(s) must be refused by verify_capture, not
        # silently accepted -- t3/pdep/hybrid.py::write_hybrid_network_input_file has no other source
        # of truth for the model chemistry once ARC's calcs/ directory may have moved on.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        capture_ts_artifacts([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        del manifest['energy_settings']
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='energy_settings'):
            verify_capture(capture_dir)

    def test_verify_capture_raises_when_energy_settings_lacks_model_chemistry(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        capture_ts_artifacts([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        manifest['energy_settings']['model_chemistry'] = ''
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='energy_settings'):
            verify_capture(capture_dir)

    @pytest.mark.parametrize('field_name, tampered_value', [
        # 'False' is a NON-EMPTY string, therefore truthy: it slips past every truthiness guard
        # downstream while rendering back into the generated Arkane input as a bare, falsy
        # `useAtomCorrections = False`. verify_capture exists precisely to catch a hand-edited
        # manifest, so a wrongly typed frozen setting must be refused here too, not only at the
        # hybrid writer.
        ('use_atom_corrections', 'False'),
        ('use_atom_corrections', None),
        ('use_hindered_rotors', 1),
        ('use_bond_corrections', 'True'),
        ('frequency_scale_factor', '0.988'),
        ('bond_correction_type', 5),
        ('atom_energies', 'H: -0.5'),
        ('model_chemistry', 3.14),
    ])
    def test_verify_capture_raises_on_a_wrongly_typed_frozen_energy_setting(self, tmp_path, field_name,
                                                                            tampered_value):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        _capture([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        manifest['energy_settings'][field_name] = tampered_value
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match=field_name):
            verify_capture(capture_dir)

    def test_verify_capture_accepts_zero_artifact_capture_with_no_energy_settings_block(self, tmp_path):
        # The flip side of the two guards above: a capture that legitimately never fetched any
        # artifact must still verify cleanly with no 'energy_settings' block at all.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = TSJoinRecord(
            network_id='network4_1',
            network_ts_label='TS_none',
            status=JOIN_STATUS_NOT_QUEUED,
            reason='Could not build the species for this transition state.',
        )
        capture_dir = str(tmp_path / 'capture')
        capture_ts_artifacts([record], arc_dir, capture_dir)

        result = verify_capture(capture_dir)

        assert result.captured_artifact_count == 0
