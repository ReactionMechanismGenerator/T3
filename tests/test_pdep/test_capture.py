"""
Tests for t3.pdep.capture.

Every fixture is built under ``tmp_path``; none reads from ``tests/data/``, since these tests
exercise how ``capture`` reacts to the *shape* of an ARC project directory (and its interaction
with discovery), not to any one recorded run of it. Mirrors the fixture-construction style of
``tests/test_pdep/test_discovery.py``.
"""

import hashlib
import multiprocessing
import os
import shutil

import pytest
import yaml

from t3.pdep import capture as capture_module
from t3.pdep.cache import hash_file
from t3.pdep.capture import CAPTURE_MANIFEST_FILE_NAME, CaptureResult, VerifyResult, capture_ts_artifacts, \
    verify_capture
from t3.pdep.discovery import (
    ARTIFACT_STATUS_MISSING,
    ARTIFACT_STATUS_NOT_QUEUED,
    ARTIFACT_STATUS_UNUSABLE,
    ARTIFACT_STATUS_UNVERIFIED,
    ARTIFACT_STATUS_USABLE,
)
from arc.common import read_yaml_file
from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.hybrid import QMEnergySettings, _read_qm_artifact, write_hybrid_network_input_file
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


def _queued_record(network_id='network4_1', network_ts_label='TS3', arc_project_directory=None,
                   coefficient=0.05, delta_ln_k=0.02):
    label = arc_ts_label(network_id, network_ts_label)
    expected_path = expected_ts_artifact_path(arc_project_directory, label) if arc_project_directory else None
    return TSJoinRecord(
        network_id=network_id,
        network_ts_label=network_ts_label,
        status=JOIN_STATUS_QUEUED,
        arc_ts_label=label,
        expected_artifact_path=expected_path,
        reason='Queued to ARC.',
        coefficient=coefficient,
        delta_ln_k=delta_ln_k,
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
    block. ``useAtomCorrections`` is ``True`` to stay faithful to ARC, which computes it as
    ``bool(model_chemistry or atom_energies)`` and so can never emit ``False`` alongside a
    ``modelChemistry`` line; a ``False`` here would be an ARC-impossible state that the hybrid
    writer rightly refuses, making the end-to-end path untestable. ``atomEnergies`` is populated
    with a realistic ARC-shaped dict (not omitted) so ``use_atom_corrections=True`` resolves to a
    real ``atom_energies`` mapping rather than ``None`` -- the fail-open state
    ``write_hybrid_network_input_file`` now refuses, which would otherwise make the end-to-end
    capture-to-hybrid-write path untestable too. ``useBondCorrections`` stays ``False`` so no
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
                "useAtomCorrections = True\n\n"
                "atomEnergies = {'C': -37.844411, 'H': -0.499818, 'N': -54.581501, 'O': -75.062219}\n\n"
                "useBondCorrections = False\n"
            )
    output_dir = os.path.join(arc_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    output_yml_path = os.path.join(output_dir, 'output.yml')
    if not os.path.isfile(output_yml_path):
        with open(output_yml_path, 'w') as f:
            f.write('{}\n')


def _networks_for(join_records: list, networks_dir: str, method: str = 'MSC') -> dict:
    """Build the sidecar-shaped networks block for the given records, backed by real stub RMG
    network files under ``networks_dir`` (created if absent) so that capture's live-hash check has
    a real file to verify. Hashes are computed with ``t3.pdep.cache.hash_file``, the same primitive
    production code records them with."""
    os.makedirs(networks_dir, exist_ok=True)
    networks = dict()
    for network_id in sorted({record.network_id for record in join_records}):
        source_path = os.path.join(networks_dir, f'{network_id}.py')
        if not os.path.isfile(source_path):
            with open(source_path, 'w') as f:
                f.write(f"# stub RMG network file\nnetwork(label='{network_id}')\n")
        networks[network_id] = {'source_path': source_path,
                                'source_sha256': hash_file(source_path),
                                'method': method}
    return networks


def _capture(join_records, arc_project_directory, capture_dir, networks=None, **kwargs):
    """Call ``capture_ts_artifacts`` with a well-formed networks block derived from the records
    (stub network files live in a sibling directory of the ARC project, standing in for RMG's
    ``pdep/``), unless the test supplies its own."""
    if networks is None:
        networks = _networks_for(join_records, networks_dir=f'{arc_project_directory}_rmg_pdep')
    if 'sensitivity_by_ts' not in kwargs:
        # Mirror `main.py::_capture_pdep_ts_artifacts`'s own derivation: the sensitivity evidence a
        # captured transition state requires comes off the (durable) join record's own frozen
        # `coefficient`/`delta_ln_k`, not from some separate in-memory selection.
        kwargs['sensitivity_by_ts'] = {record.key: (record.coefficient, record.delta_ln_k) for record in join_records}
    return capture_ts_artifacts(join_records, arc_project_directory, capture_dir,
                                networks=networks, **kwargs)


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

        result = _capture([record], arc_dir, capture_dir)

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

        result = _capture([record], arc_dir, capture_dir)
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
        result = _capture(
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

        result = _capture([record], arc_dir, capture_dir)

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
        result = _capture([record], arc_dir, capture_dir)
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

        result = _capture([record], arc_dir, capture_dir)

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

        result = _capture([record], arc_dir, capture_dir)

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

        _capture([record], arc_dir, capture_dir)

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

        _capture([record], arc_dir, capture_dir)

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

        _capture([record], arc_dir, capture_dir)

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
        first_result = _capture([good_record], arc_dir, capture_dir)
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
        second_result = _capture([good_record, bad_record], arc_dir, capture_dir)

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
        _capture([good_record], arc_dir, capture_dir)
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
            _capture([good_record, broken_record], arc_dir, capture_dir)

        after_files = sorted(
            os.path.relpath(os.path.join(root, name), capture_dir)
            for root, _dirs, names in os.walk(capture_dir) for name in names
        )
        assert after_files == previous_files
        for relpath in previous_files:
            if relpath.endswith('.py') and 'qm' in relpath.split(os.sep):
                with open(os.path.join(capture_dir, relpath)) as f:
                    assert 'Log(' in f.read()

    def test_a_crash_during_the_atomic_swap_leaves_the_prior_good_capture_intact(self, tmp_path, monkeypatch):
        # The swap step (rename capture_dir aside, os.replace staging_dir into place, remove the
        # renamed-aside tree) is the one place a fully-staged, already-self-verified new capture
        # touches the previous good one at all. Simulate a crash exactly there -- os.replace itself
        # raising, e.g. an OSError from a full disk or a cross-device link -- and confirm
        # capture_ts_artifacts renames the old tree straight back rather than leaving capture_dir
        # empty or half-populated. A naive destroy-then-copy swap (no rename-aside, no restore-on-
        # failure) would pass every OTHER test in this file (none of them inject a failure inside
        # the swap itself) while silently destroying the prior good capture the instant os.replace
        # failed -- this test exists specifically to kill that mutation.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        capture_dir = str(tmp_path / 'capture')

        good_record = _usable_record(tmp_path=tmp_path.__class__(arc_dir), network_ts_label='TS_good')
        first_result = _capture([good_record], arc_dir, capture_dir)
        previous_manifest_mtime = os.path.getmtime(first_result.manifest_path)
        previous_files = sorted(
            os.path.relpath(os.path.join(root, name), capture_dir)
            for root, _dirs, names in os.walk(capture_dir) for name in names
        )
        previous_artifact_content = open(first_result.records[0].artifact_path).read()

        real_os_replace = capture_module.os.replace

        def failing_replace(src, dst):
            # Only fail the staging_dir -> capture_dir swap itself; anything else (there is
            # nothing else in this codepath, but fail closed rather than assume) uses the real
            # implementation.
            if os.path.basename(dst) == os.path.basename(capture_dir):
                raise OSError('simulated crash during atomic swap')
            return real_os_replace(src, dst)

        monkeypatch.setattr(capture_module.os, 'replace', failing_replace)

        with pytest.raises(OSError, match='simulated crash during atomic swap'):
            _capture([good_record], arc_dir, capture_dir)

        # The prior good capture must be exactly as it was: same files, same manifest mtime, same
        # artifact content -- never observed empty (old tree destroyed before the swap) or
        # half-populated (new tree partially applied).
        assert os.path.isdir(capture_dir)
        after_files = sorted(
            os.path.relpath(os.path.join(root, name), capture_dir)
            for root, _dirs, names in os.walk(capture_dir) for name in names
        )
        assert after_files == previous_files
        assert os.path.getmtime(first_result.manifest_path) == previous_manifest_mtime
        assert open(first_result.records[0].artifact_path).read() == previous_artifact_content
        # No leftover staging/old-capture scratch directories beside capture_dir.
        parent_dir = os.path.dirname(os.path.abspath(capture_dir))
        leftover = [name for name in os.listdir(parent_dir)
                    if name.startswith('.capture-staging-') or name.startswith('.capture-old-')]
        assert leftover == []

    def test_a_real_process_death_between_rename_aside_and_replace_is_recoverable(self, tmp_path):
        # Defect: the swap is rename-aside (capture_dir -> '.capture-old-*') THEN a separate
        # os.replace (staging_dir -> capture_dir) -- two syscalls, not one. The test above only
        # covers a caught in-process exception (os.replace raising); it proves nothing about a
        # process that dies for REAL (kill -9, OOM-kill, power loss) between the two syscalls,
        # since that skips every except/finally handler in _capture_ts_artifacts_locked entirely.
        # Simulate that for real: fork a child, let it perform the actual rename-aside, then call
        # os._exit() immediately after -- before it ever reaches os.replace. Confirm the defect
        # window is real (capture_dir observably absent, prior good capture sitting under
        # '.capture-old-*'), then confirm _recover_capture_swap_state -- run by the next call that
        # holds the lock -- restores the exact pre-crash capture.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        capture_dir = str(tmp_path / 'capture')

        good_record = _usable_record(tmp_path=tmp_path.__class__(arc_dir), network_ts_label='TS_good')
        first_result = _capture([good_record], arc_dir, capture_dir)
        previous_manifest_mtime = os.path.getmtime(first_result.manifest_path)
        previous_artifact_content = open(first_result.records[0].artifact_path).read()
        previous_files = sorted(
            os.path.relpath(os.path.join(root, name), capture_dir)
            for root, _dirs, names in os.walk(capture_dir) for name in names
        )

        second_arc_dir = str(tmp_path / 'arc_project_2')
        os.makedirs(second_arc_dir)
        second_record = _usable_record(tmp_path=tmp_path.__class__(second_arc_dir), network_ts_label='TS_second')

        def _child(arc_dir_, capture_dir_, record_):
            # Runs in a forked child (fork start method: no pickling, the child is a full memory
            # copy of this process). Patch os.rename so the SWAP's rename-aside call -- identified
            # by its destination being a '.capture-old-*' sibling, exactly as production code names
            # it -- performs the real rename and then calls os._exit immediately: a genuine process
            # death that skips every Python-level except/finally in capture.py, landing exactly
            # between the rename-aside and the subsequent os.replace.
            real_rename = os.rename

            def crash_after_rename_aside(src, dst):
                if os.path.basename(dst).startswith(capture_module._CAPTURE_OLD_DIR_PREFIX):
                    real_rename(src, dst)
                    os._exit(1)
                return real_rename(src, dst)

            capture_module.os.rename = crash_after_rename_aside
            try:
                _capture([record_], arc_dir_, capture_dir_)
            finally:
                # Only reached if the crash point above was never hit (a bug in this test's own
                # setup, not in capture.py) -- exit nonzero regardless so the parent never mistakes
                # this for a clean, non-crashing run.
                os._exit(2)

        ctx = multiprocessing.get_context('fork')
        child = ctx.Process(target=_child, args=(second_arc_dir, capture_dir, second_record))
        child.start()
        child.join(timeout=30)
        assert child.exitcode == 1, (
            f"child did not die at the intended crash point (exitcode {child.exitcode!r}); this "
            "test's own crash injection did not fire as expected, not a capture.py failure."
        )

        # Empirically confirm the defect window is real: capture_dir must be absent (rename-aside
        # completed and removed it; os.replace never ran to repopulate it), with exactly one valid
        # '.capture-old-*' sibling holding the untouched prior good capture.
        assert not os.path.isdir(capture_dir)
        parent_dir = os.path.dirname(os.path.abspath(capture_dir))
        old_dirs = [name for name in os.listdir(parent_dir)
                    if name.startswith(capture_module._CAPTURE_OLD_DIR_PREFIX)]
        assert len(old_dirs) == 1
        recovered_manifest_path = os.path.join(parent_dir, old_dirs[0], CAPTURE_MANIFEST_FILE_NAME)
        assert os.path.getmtime(recovered_manifest_path) == previous_manifest_mtime

        # The crashed child died while holding capture_dir's lock -- clean it up first so the
        # recovery call below isn't blocked by an unrelated stale-lock concern (that is covered by
        # its own dedicated test).
        stale_lock_path = capture_module._capture_lock_path(capture_dir)
        if os.path.isfile(stale_lock_path):
            os.remove(stale_lock_path)

        # The next lock-holding call recovers it: _recover_capture_swap_state restores capture_dir
        # to EXACTLY the pre-crash content -- not the second (never-swapped-in) attempt's content.
        capture_module._recover_capture_swap_state(capture_dir=capture_dir)
        assert os.path.isdir(capture_dir)
        after_files = sorted(
            os.path.relpath(os.path.join(root, name), capture_dir)
            for root, _dirs, names in os.walk(capture_dir) for name in names
        )
        assert after_files == previous_files
        assert os.path.getmtime(first_result.manifest_path) == previous_manifest_mtime
        assert open(first_result.records[0].artifact_path).read() == previous_artifact_content
        verify_result = verify_capture(capture_dir)
        assert verify_result.record_count == 1
        assert verify_result.captured_artifact_count == 1

        leftover = [name for name in os.listdir(parent_dir)
                    if name.startswith(capture_module._CAPTURE_OLD_DIR_PREFIX)
                    or name.startswith(capture_module._CAPTURE_STAGING_DIR_PREFIX)]
        assert leftover == []

    def test_a_second_capture_call_against_a_locked_capture_dir_is_refused(self, tmp_path):
        # Defect: without interprocess exclusion, two processes could both pass
        # _refuse_unowned_capture_dir, build two competing staging directories, and race the
        # rename-aside/os.replace swap against each other. Simulate a concurrent process already
        # capturing into capture_dir by taking its lock directly (with OUR pid recorded, so it is
        # never mistaken for stale), then assert a second, genuinely concurrent capture_ts_artifacts
        # call over the SAME capture_dir raises rather than proceeding.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        capture_dir = str(tmp_path / 'capture')
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir), network_ts_label='TS_locked')

        lock_path = capture_module._acquire_capture_lock(capture_dir=capture_dir)
        try:
            with pytest.raises(RuntimeError, match='already held'):
                _capture([record], arc_dir, capture_dir)
            # Refusing to take the lock must mean NOTHING was written: capture_dir must not even
            # be created, let alone partially populated, by the call that failed to get the lock.
            assert not os.path.exists(capture_dir)
        finally:
            capture_module._release_capture_lock(lock_path)

        # Once the lock is released, an otherwise-identical call succeeds normally.
        _capture([record], arc_dir, capture_dir)
        assert os.path.isdir(capture_dir)

    def test_a_lock_held_by_a_dead_process_is_reclaimed_rather_than_blocking_forever(self, tmp_path):
        # A process that dies while holding the lock (kill -9, OOM-kill -- anything that skips the
        # 'finally: _release_capture_lock(...)' in capture_ts_artifacts) must not permanently wedge
        # this capture_dir: the next call must recognize the recorded PID no longer corresponds to a
        # live process and reclaim the lock, rather than failing closed forever.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        capture_dir = str(tmp_path / 'capture')
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir), network_ts_label='TS_stale')

        lock_path = capture_module._capture_lock_path(capture_dir)
        os.makedirs(os.path.dirname(lock_path), exist_ok=True)
        # Fabricate a lock file recording a PID guaranteed to be dead: fork a short-lived child,
        # reap it, and use its now-freed PID -- a real dead PID, not a guessed-large/unlikely one.
        child_pid = os.fork()
        if child_pid == 0:
            os._exit(0)
        os.waitpid(child_pid, 0)
        with open(lock_path, 'w') as f:
            f.write(f'{child_pid}\n')

        _capture([record], arc_dir, capture_dir)
        assert os.path.isdir(capture_dir)

    def test_multiple_ambiguous_old_capture_siblings_refuses_to_guess_which_is_real(self, tmp_path):
        # If _recover_capture_swap_state ever found MORE than one sibling '.capture-old-*'
        # directory that each hold a valid manifest while capture_dir itself is missing/empty,
        # silently picking one would risk resurrecting the wrong (possibly stale) prior capture.
        # This should never arise from a single crashed swap -- but the guard must refuse outright
        # rather than guess, if it somehow does.
        # Fabricate two independently-valid '.capture-old-*' siblings directly -- a valid manifest
        # is just a dict with 'transition_states' and 'arc_project_directory' keys per
        # _has_valid_capture_manifest, so this does not need two full real captures (which, routed
        # through the SAME capture_dir, would trigger capture_dir's own crash-recovery in between
        # and silently consume the first one before the second existed).
        capture_dir = str(tmp_path / 'capture')
        manifest = {'transition_states': [], 'arc_project_directory': '/nonexistent'}
        old_a = os.path.join(str(tmp_path), capture_module._CAPTURE_OLD_DIR_PREFIX + 'a')
        old_b = os.path.join(str(tmp_path), capture_module._CAPTURE_OLD_DIR_PREFIX + 'b')
        os.makedirs(old_a)
        os.makedirs(old_b)
        with open(os.path.join(old_a, CAPTURE_MANIFEST_FILE_NAME), 'w') as f:
            yaml.dump(manifest, f)
        with open(os.path.join(old_b, CAPTURE_MANIFEST_FILE_NAME), 'w') as f:
            yaml.dump(manifest, f)

        # capture_dir is absent, with two independently-valid '.capture-old-*' siblings.
        assert not os.path.exists(capture_dir)
        with pytest.raises(ValueError, match='Refusing to recover'):
            capture_module._recover_capture_swap_state(capture_dir=capture_dir)

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
            _capture([first, second], arc_dir, capture_dir)

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
            _capture([record], arc_dir, capture_dir)

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
        first_result = _capture([good_record], arc_dir, capture_dir)
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
        second_result = _capture([zero_record], arc_dir, capture_dir, supersede=True)

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
        first_result = _capture([good_record], arc_dir, capture_dir)
        captured_artifact_path = first_result.records[0].artifact_path
        assert os.path.isfile(captured_artifact_path)

        zero_record = TSJoinRecord(
            network_id='network4_1',
            network_ts_label='TS_none',
            status=JOIN_STATUS_NOT_QUEUED,
            reason='Could not build the species for this transition state.',
        )
        with pytest.raises(ValueError, match='supersede'):
            _capture([zero_record], arc_dir, capture_dir)

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
            _capture([first, second], arc_dir, capture_dir)

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
            _capture([record], arc_dir, capture_dir)

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
            _capture([record], arc_dir, capture_dir)

        message = str(exc_info.value)
        assert capture_dir in message
        # Refused before any capture I/O -- the unrelated file must still be the only thing there.
        assert os.listdir(capture_dir) == ['unrelated.txt']

    def test_verify_capture_detects_tampered_captured_log(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        result = _capture([record], arc_dir, capture_dir)
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

        _capture([record], arc_dir, capture_dir)
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
            _capture(join_records=[], arc_project_directory=arc_dir, capture_dir=capture_dir)

        assert not os.path.exists(capture_dir), 'a refused capture must not leave a capture directory behind'

    def test_verify_capture_returns_concrete_counts_on_success(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        capture_result = _capture([record], arc_dir, capture_dir)
        captured_artifact_path = capture_result.records[0].artifact_path
        assert os.path.isfile(captured_artifact_path)

        result = verify_capture(capture_dir)

        assert isinstance(result, VerifyResult)
        assert result.capture_dir == capture_dir
        assert result.manifest_path == os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        assert result.record_count == 1
        assert result.captured_artifact_count == 1

    def test_verify_capture_returns_the_verified_ts_records_for_the_writer_to_consume(self, tmp_path):
        # P1 #1: the writer (t3.main._write_pdep_hybrid_network_inputs) used to re-read the
        # manifest itself after verify_capture already read and validated it, and rebuild
        # TSArtifactRecords from that second, unverified read. verify_capture must hand back the
        # records it already validated (with resolved, confined artifact paths) so the writer
        # never needs a second read of untrusted disk state.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        capture_result = _capture([record], arc_dir, capture_dir)
        captured_artifact_path = capture_result.records[0].artifact_path

        result = verify_capture(capture_dir)

        assert len(result.ts_records) == 1
        ts_record = result.ts_records[0]
        assert ts_record.network_id == record.network_id
        assert ts_record.network_ts_label == record.network_ts_label
        assert ts_record.status == ARTIFACT_STATUS_USABLE
        # The path handed back must already be resolved against capture_dir (an absolute,
        # existing path), not the raw manifest-relative string -- a consumer that only ever sees
        # this resolved value cannot reconstruct an unconfined path out of it.
        assert ts_record.artifact_path == captured_artifact_path
        assert os.path.isfile(ts_record.artifact_path)

    def test_verify_capture_refuses_a_usable_status_with_no_captured_artifact_path(self, tmp_path):
        # P1 #2: verify_capture never used to look at 'status' at all -- it counted a captured
        # artifact purely by 'captured_artifact_path is not None'. A hand-edited (or corrupted)
        # manifest claiming status: usable but with captured_artifact_path: null would verify
        # clean and silently be treated as zero captured artifacts by any consumer that (rightly)
        # trusts 'usable' to mean "there is a vendored artifact for this". That is exactly the
        # kind of torn/tampered state verify_capture exists to catch.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        _capture([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        entry = manifest['transition_states'][0]
        assert entry['status'] == ARTIFACT_STATUS_USABLE
        entry['captured_artifact_path'] = None
        entry['captured_artifact_sha256'] = None
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='usable'):
            verify_capture(capture_dir)

    def test_verify_capture_refuses_a_non_null_captured_artifact_path_on_a_status_forbidding_one(self, tmp_path):
        # The converse of the above: a status that discovery would never pair with an artifact
        # path (e.g. 'missing') showing up with a non-null captured_artifact_path is just as
        # untrustworthy a manifest as the reverse pairing -- both mean the status/path pairing no
        # longer matches what discovery would have produced.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        _capture([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        entry = manifest['transition_states'][0]
        entry['status'] = ARTIFACT_STATUS_MISSING
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='missing'):
            verify_capture(capture_dir)

    def test_verify_capture_refuses_an_unrecognized_status(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        _capture([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        manifest['transition_states'][0]['status'] = 'not_a_real_status'
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        # 'Refusing to verify' pins this to verify_capture's OWN status guard rather than the
        # downstream, differently-worded 'Unrecognized transition state artifact status' raised by
        # TSArtifactRecord.__init__ -- a loose 'unrecognized' match would pass even if
        # verify_capture's own check were deleted, since that constructor enforces the same
        # invariant redundantly a few lines later.
        with pytest.raises(ValueError, match='Refusing to verify.*unrecognized status'):
            verify_capture(capture_dir)

    def test_verify_capture_refuses_a_missing_status_field(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        _capture([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        del manifest['transition_states'][0]['status']
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        # Same reasoning as the sibling unrecognized-status test above: a missing status field
        # reads as `entry.get('status')` -> None, which is also not in TS_ARTIFACT_STATUSES, so the
        # same 'Refusing to verify ... unrecognized status' guard fires here. Pinning to that exact
        # phrasing (rather than a bare 'status' match) is what makes this test actually exercise
        # verify_capture's own check instead of passing on TSArtifactRecord's redundant one.
        with pytest.raises(ValueError, match='Refusing to verify.*unrecognized status'):
            verify_capture(capture_dir)

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
            _capture([record], arc_dir, capture_dir)

        assert 'Torn or corrupt' in str(exc_info.value)

    def test_energy_settings_is_frozen_into_result_and_manifest_when_artifact_captured(self, tmp_path):
        # Any capture with >=1 usable artifact must freeze the ARC energy-reference settings that
        # produced it (t3/pdep/hybrid.py::write_hybrid_network_input_file requires a non-blank
        # model_chemistry downstream); _usable_record's fixture sets modelChemistry = 'CBS-QB3'.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        result = _capture([record], arc_dir, capture_dir)

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
        # True, matching the fixture, which matches ARC: use_aec is computed as
        # bool(model_chemistry or atom_energies), so ARC cannot emit False alongside a modelChemistry
        # line. A False here would be an ARC-impossible state that the hybrid writer rightly refuses.
        assert result.energy_settings['use_atom_corrections'] is True
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

        result = _capture([record], arc_dir, capture_dir)

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
            _capture([record], arc_dir, capture_dir)

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
        _capture([record], arc_dir, capture_dir)

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
        _capture([record], arc_dir, capture_dir)

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

    def test_verify_capture_raises_when_a_single_frozen_energy_setting_key_is_missing(self, tmp_path):
        # Distinct from test_verify_capture_raises_when_energy_settings_block_is_missing_but_
        # artifacts_were_captured (which drops the WHOLE 'energy_settings' key) and from
        # test_verify_capture_raises_on_a_wrongly_typed_frozen_energy_setting (which always
        # re-assigns a field to a tampered VALUE, never removes it): this drops exactly one field
        # from an otherwise-intact 'energy_settings' dict. QMEnergySettings.from_frozen must fail
        # closed on a missing key rather than defaulting it to None, since a manifest hand-edited to
        # drop e.g. 'atom_energies' would otherwise silently reach write_hybrid_network_input_file
        # as if atom_energies had legitimately never been set.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        _capture([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        del manifest['energy_settings']['atom_energies']
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='atom_energies'):
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
        _capture([record], arc_dir, capture_dir)

        result = verify_capture(capture_dir)

        assert result.captured_artifact_count == 0

    def test_sensitivity_evidence_is_frozen_into_manifest_when_artifact_captured(self, tmp_path):
        # The evidence that justified selecting this transition state (frozen onto the join record
        # at queue time, per _queued_record's defaults) must be carried into the manifest entry, so
        # a partial-capture decision can be made from the capture alone -- never re-derived from an
        # in-memory selection, which is empty on a restart.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        _capture([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        entry = manifest['transition_states'][0]
        assert entry['coefficient'] == 0.05
        assert entry['delta_ln_k'] == 0.02

    def test_sensitivity_evidence_is_frozen_for_a_not_queued_entry_too(self, tmp_path):
        # A transition state that never reached ARC (NOT_QUEUED, no captured artifact) still carries
        # the evidence that made it a candidate: a partial capture needs to be able to tell whether a
        # transition state that never got queued was in fact the dominant one.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = TSJoinRecord(
            network_id='network4_1',
            network_ts_label='TS_none',
            status=JOIN_STATUS_NOT_QUEUED,
            reason='Could not build the species for this transition state.',
            coefficient=0.07,
            delta_ln_k=0.03,
        )
        capture_dir = str(tmp_path / 'capture')

        _capture([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        assert manifest['transition_states'], 'the not-queued entry must still be recorded in the manifest'
        entry = manifest['transition_states'][0]
        assert entry['coefficient'] == 0.07
        assert entry['delta_ln_k'] == 0.03
        assert entry['captured_artifact_path'] is None

    def test_verify_capture_raises_when_a_captured_entrys_sensitivity_evidence_is_missing(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        _capture([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        del manifest['transition_states'][0]['coefficient']
        del manifest['transition_states'][0]['delta_ln_k']
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='coefficient'):
            verify_capture(capture_dir)

    @pytest.mark.parametrize('field_name, tampered_value', [
        ('coefficient', None),
        ('delta_ln_k', None),
        ('coefficient', float('nan')),
        ('delta_ln_k', float('inf')),
        ('coefficient', 'not_a_number'),
        ('delta_ln_k', True),
    ])
    def test_verify_capture_raises_on_a_non_finite_or_missing_captured_sensitivity_value(
            self, tmp_path, field_name, tampered_value):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')
        _capture([record], arc_dir, capture_dir)

        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        manifest['transition_states'][0][field_name] = tampered_value
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='coefficient|delta_ln_k'):
            verify_capture(capture_dir)


class TestNetworkVendoring:
    """The vendored-network-source layer: the capture must carry its own copy of each RMG
    ``network.py`` (plus the frozen ME method), verified against the identity the join sidecar
    recorded, so downstream hybrid-network writing never needs the RMG ``pdep/`` directory or the
    ARC project directory to still exist."""

    def _captured_setup(self, tmp_path):
        """One usable record, its networks block (with a real stub network file), and a completed
        capture. Returns ``(record, networks, capture_dir)``."""
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        networks = _networks_for([record], networks_dir=str(tmp_path / 'rmg_pdep'))
        capture_dir = str(tmp_path / 'capture')
        result = _capture([record], arc_dir, capture_dir, networks=networks)
        return record, networks, capture_dir, result

    def test_network_file_is_vendored_byte_identical(self, tmp_path):
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        source_path = networks[record.network_id]['source_path']
        vendored_path = os.path.join(capture_dir, 'networks', f'{record.network_id}.py')

        assert os.path.isfile(vendored_path)
        with open(source_path, 'rb') as f:
            source_bytes = f.read()
        with open(vendored_path, 'rb') as f:
            vendored_bytes = f.read()
        assert len(source_bytes) > 0
        assert vendored_bytes == source_bytes

    def test_manifest_records_an_authoritative_networks_block(self, tmp_path):
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        with open(os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)) as f:
            manifest = yaml.safe_load(f)

        assert record.network_id in manifest['networks']
        entry = manifest['networks'][record.network_id]
        assert entry['source_path'] == networks[record.network_id]['source_path']
        # The recorded hash is the full 'sha256:<hex>' string, verified against an independent
        # hash implementation so the test never shares one with the code under test.
        assert entry['source_sha256'] == f"sha256:{_sha256_of(entry['source_path'])}"
        assert entry['method'] == 'MSC'
        assert entry['captured_path'] == os.path.join('networks', f'{record.network_id}.py')

    def test_capture_result_surfaces_the_networks_mapping(self, tmp_path):
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        assert isinstance(result.networks, dict)
        assert set(result.networks) == {record.network_id}
        assert result.networks[record.network_id]['method'] == 'MSC'
        assert result.networks[record.network_id]['captured_path'] == \
            os.path.join('networks', f'{record.network_id}.py')

    def test_capture_survives_deletion_of_rmg_pdep_and_arc_project(self, tmp_path):
        # THE acceptance property of the whole vendoring layer: after a successful capture, both
        # source trees can vanish entirely and the capture alone still verifies and yields a usable
        # network source (the vendored copy), the frozen ME method, and the frozen energy settings.
        record, networks, capture_dir, result = self._captured_setup(tmp_path)

        shutil.rmtree(str(tmp_path / 'rmg_pdep'))
        shutil.rmtree(str(tmp_path / 'arc_project'))

        verify_result = verify_capture(capture_dir)
        assert verify_result.captured_artifact_count == 1
        assert isinstance(verify_result.networks, dict)
        assert len(verify_result.networks) == 1
        entry = verify_result.networks[record.network_id]
        vendored_abs_path = os.path.join(capture_dir, entry['captured_path'])
        assert os.path.isfile(vendored_abs_path)
        with open(vendored_abs_path) as f:
            assert record.network_id in f.read()
        assert entry['method'] == 'MSC'
        assert isinstance(verify_result.energy_settings, dict)
        assert verify_result.energy_settings['model_chemistry'] == 'CBS-QB3'

    def test_capture_raises_when_live_network_hash_mismatches_the_recorded_one(self, tmp_path):
        # "The world changed under us": the RMG network file was rewritten between the SA/selection
        # pass that recorded its hash and this capture pass. The capture must refuse -- vendoring
        # the new bytes would freeze a network the selection never actually examined.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        networks = _networks_for([record], networks_dir=str(tmp_path / 'rmg_pdep'))
        with open(networks[record.network_id]['source_path'], 'a') as f:
            f.write('# the network file changed after selection ran\n')
        capture_dir = str(tmp_path / 'capture')

        with pytest.raises(ValueError, match='source_sha256|changed'):
            _capture([record], arc_dir, capture_dir, networks=networks)

        assert not os.path.exists(capture_dir), (
            'a capture refused during the network pre-flight must not leave a capture directory '
            'behind -- the pre-flight runs before any write.'
        )

    def test_capture_raises_specifically_when_a_network_file_is_missing(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        networks = _networks_for([record], networks_dir=str(tmp_path / 'rmg_pdep'))
        os.remove(networks[record.network_id]['source_path'])
        capture_dir = str(tmp_path / 'capture')

        with pytest.raises(ValueError, match=record.network_id):
            _capture([record], arc_dir, capture_dir, networks=networks)

        assert not os.path.exists(capture_dir)

    def test_capture_refuses_a_networks_block_missing_a_referenced_network(self, tmp_path):
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        capture_dir = str(tmp_path / 'capture')

        with pytest.raises(ValueError, match=record.network_id):
            _capture([record], arc_dir, capture_dir, networks=dict())

        assert not os.path.exists(capture_dir)

    def test_a_corrupted_network_copy_is_refused_at_capture_time(self, tmp_path, monkeypatch):
        # Simulates mid-copy corruption of the network source vendoring itself (the same failure
        # mode test_captured_log_hash_mismatch_during_copy_is_refused simulates for logs): the
        # vendored copy's hash must be re-checked against the sidecar's recorded source_sha256
        # immediately after the copy, and any disagreement refused on the spot.
        real_copyfile = shutil.copyfile

        def corrupt_network_copyfile(src, dst):
            if os.path.basename(src).startswith('network'):
                with open(dst, 'w') as f:
                    f.write('corrupted network copy\n')
                return dst
            return real_copyfile(src, dst)

        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        networks = _networks_for([record], networks_dir=str(tmp_path / 'rmg_pdep'))
        capture_dir = str(tmp_path / 'capture')
        monkeypatch.setattr(capture_module.shutil, 'copyfile', corrupt_network_copyfile)

        with pytest.raises(ValueError, match='Torn or corrupt'):
            _capture([record], arc_dir, capture_dir, networks=networks)

    def test_a_corrupted_network_copy_during_recapture_does_not_destroy_the_prior_good_capture(
            self, tmp_path, monkeypatch):
        # Defect A (round-22 P1 #2): capture supersession is not transactional. A RE-capture over
        # an already-good capture_dir that dies partway through _vendor_network_sources must leave
        # the ORIGINAL good networks/ file exactly as it was -- not a destroy-then-half-rebuild
        # torn state with no way back. This is the probe that distinguishes the finding from
        # test_a_corrupted_network_copy_is_refused_at_capture_time above, which only ever exercises
        # a FIRST capture into an empty capture_dir and never asserts anything about prior data.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        networks = _networks_for([record], networks_dir=str(tmp_path / 'rmg_pdep'))
        capture_dir = str(tmp_path / 'capture')

        # First capture: succeeds, leaves a good vendored network file behind.
        _capture([record], arc_dir, capture_dir, networks=networks)
        vendored_path = os.path.join(capture_dir, 'networks', f'{record.network_id}.py')
        with open(vendored_path, 'rb') as f:
            good_bytes = f.read()
        assert len(good_bytes) > 0

        # Second capture (re-capture over the same, already-good capture_dir): corrupt the
        # network copy mid-vendoring, simulating a crash/tamper partway through supersession.
        real_copyfile = shutil.copyfile

        def corrupt_network_copyfile(src, dst):
            if os.path.basename(src).startswith('network'):
                with open(dst, 'w') as f:
                    f.write('corrupted network copy\n')
                return dst
            return real_copyfile(src, dst)

        monkeypatch.setattr(capture_module.shutil, 'copyfile', corrupt_network_copyfile)

        with pytest.raises(ValueError, match='Torn or corrupt'):
            _capture([record], arc_dir, capture_dir, networks=networks)

        # The prior good capture must survive completely untouched -- transactional supersession
        # means this failed re-capture can never leave the old good state gone or half-overwritten.
        assert os.path.isfile(vendored_path), (
            'a failed re-capture destroyed the previously good vendored network file -- '
            'capture supersession is not transactional'
        )
        with open(vendored_path, 'rb') as f:
            assert f.read() == good_bytes, (
                'a failed re-capture left a corrupted/partial network file in place of the '
                'good one -- capture supersession is not transactional'
            )

    def test_pre_swap_self_check_refuses_a_staged_capture_whose_manifest_networks_block_is_incomplete(
            self, tmp_path, monkeypatch):
        # capture_ts_artifacts calls verify_capture(staging_dir) as a self-check BEFORE swapping the
        # staged tree into capture_dir's place -- this is distinct from the per-copy hash checks
        # _vendor_network_sources already performs during vendoring (covered by
        # test_a_corrupted_network_copy_is_refused_at_capture_time above): those catch a corrupted
        # COPY, not an internal bug that silently drops a network from the manifest's 'networks'
        # block after every copy already vendored and hashed correctly. Simulate exactly that: wrap
        # _vendor_network_sources so the real vendoring/hashing succeeds (every file on disk is
        # fine), but the dict handed back for the manifest is missing one entry -- the only way this
        # is ever caught is verify_capture's own 'networks' completeness check, run here as a
        # pre-swap self-check. Without that self-check call, this internally-inconsistent capture
        # would be swapped into capture_dir's place and returned as if it were good.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        record = _usable_record(tmp_path=tmp_path.__class__(arc_dir))
        networks = _networks_for([record], networks_dir=str(tmp_path / 'rmg_pdep'))
        capture_dir = str(tmp_path / 'capture')

        real_vendor_network_sources = capture_module._vendor_network_sources

        def dropping_vendor_network_sources(networks, capture_dir):
            result = real_vendor_network_sources(networks=networks, capture_dir=capture_dir)
            # The real vendoring/hashing already happened correctly; only the returned manifest
            # dict is truncated, simulating a bug that loses track of an already-vendored network.
            result.pop(record.network_id, None)
            return result

        monkeypatch.setattr(capture_module, '_vendor_network_sources', dropping_vendor_network_sources)

        with pytest.raises(ValueError, match=record.network_id):
            _capture([record], arc_dir, capture_dir, networks=networks)

        # First-time capture: a staged tree that fails its own pre-swap self-check must never be
        # swapped into capture_dir's place at all.
        assert not os.path.exists(capture_dir), (
            'a staged capture that fails its own pre-swap self-check (verify_capture(staging_dir)) '
            'was swapped into capture_dir anyway -- the pre-swap self-check is not load-bearing'
        )

    def test_a_network_id_that_escapes_the_networks_subdir_is_refused(self, tmp_path):
        # A network_id is only ever label-validated for QUEUED transition states (arc_ts_label
        # refuses unsafe components); a not_queued record can legitimately carry any string, so the
        # vendoring step itself must confine the vendored file name to the networks/ subdirectory
        # rather than let 'net/../x'-style ids write outside it.
        arc_dir = str(tmp_path / 'arc_project')
        os.makedirs(arc_dir)
        evil_id = os.path.join('..', 'escaped')
        record = TSJoinRecord(network_id=evil_id,
                              network_ts_label='TS1',
                              status=JOIN_STATUS_NOT_QUEUED,
                              reason='network_id contains unsafe characters')
        networks_dir = str(tmp_path / 'rmg_pdep')
        os.makedirs(networks_dir)
        source_path = os.path.join(networks_dir, 'evil.py')
        with open(source_path, 'w') as f:
            f.write('# stub\n')
        networks = {evil_id: {'source_path': source_path,
                              'source_sha256': hash_file(source_path),
                              'method': 'MSC'}}
        capture_dir = str(tmp_path / 'capture')

        with pytest.raises(ValueError, match='outside the capture'):
            _capture([record], arc_dir, capture_dir, networks=networks)

        assert not os.path.exists(os.path.join(capture_dir, '..', 'escaped.py'))
        assert not os.path.exists(str(tmp_path / 'escaped.py'))

    def test_verify_catches_a_deleted_vendored_network_file(self, tmp_path):
        # Mutation 2 of the brief, as a live test: the manifest names a vendored network that no
        # longer exists on disk. This must be a torn-capture failure, never a silent pass.
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        os.remove(os.path.join(capture_dir, 'networks', f'{record.network_id}.py'))

        with pytest.raises(ValueError, match='network'):
            verify_capture(capture_dir)

    def test_verify_catches_a_tampered_vendored_network_file(self, tmp_path):
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        with open(os.path.join(capture_dir, 'networks', f'{record.network_id}.py'), 'a') as f:
            f.write('# tampered after capture\n')

        with pytest.raises(ValueError, match='sha256'):
            verify_capture(capture_dir)

    def test_verify_refuses_an_empty_networks_block_when_transition_states_reference_networks(self, tmp_path):
        # The anti-vacuity case called out by the brief: transition_states is non-empty and every
        # entry names a network_id, so an empty (or absent) networks block cannot "pass because
        # there was nothing to iterate over".
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        manifest['networks'] = dict()
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='networks'):
            verify_capture(capture_dir)

    def test_verify_refuses_a_manifest_with_no_networks_key(self, tmp_path):
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        del manifest['networks']
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='networks'):
            verify_capture(capture_dir)

    def test_verify_refuses_an_unreferenced_networks_entry(self, tmp_path):
        # The other direction: a networks entry no transition state references asserts an identity
        # for work this capture never contained.
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        manifest['networks']['network9_9'] = dict(manifest['networks'][record.network_id])
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='network9_9'):
            verify_capture(capture_dir)

    def test_verify_refuses_an_invalid_method_in_the_networks_block(self, tmp_path):
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        manifest['networks'][record.network_id]['method'] = 'definitely_not_an_me_method'
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='method'):
            verify_capture(capture_dir)

    def test_a_stale_vendored_network_from_a_prior_capture_is_replaced(self, tmp_path):
        # A re-capture into an owned capture_dir must leave the vendored networks/ directory in
        # exact agreement with the new manifest: a corrupted or stale file with the same name is
        # overwritten (never silently kept), and a leftover file from a network the new capture
        # does not contain is removed.
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        vendored_path = os.path.join(capture_dir, 'networks', f'{record.network_id}.py')
        with open(vendored_path, 'w') as f:
            f.write('# stale bytes from an interrupted prior run\n')
        stale_extra_path = os.path.join(capture_dir, 'networks', 'network9_9.py')
        with open(stale_extra_path, 'w') as f:
            f.write('# leftover from a prior capture of a different network set\n')

        arc_dir = str(tmp_path / 'arc_project')
        _capture([record], arc_dir, capture_dir, networks=networks)

        with open(networks[record.network_id]['source_path'], 'rb') as f:
            source_bytes = f.read()
        with open(vendored_path, 'rb') as f:
            vendored_bytes = f.read()
        assert vendored_bytes == source_bytes
        assert not os.path.exists(stale_extra_path)
        verify_capture(capture_dir)

    def test_verify_refuses_an_absolute_captured_artifact_path_that_escapes_capture_dir(self, tmp_path):
        # Defect B (round-22 P1 #3): verify_capture must confine every manifest-supplied path to
        # capture_dir. os.path.join(a, '/abs/path') silently returns '/abs/path', so an absolute
        # captured_artifact_path in the manifest would otherwise make verification read (and
        # "pass") an arbitrary external file instead of the capture's own vendored copy.
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        outside_file = tmp_path / 'outside_capture_dir.py'
        outside_file.write_bytes(b'an external file the tampered manifest points at\n')
        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        manifest['transition_states'][0]['captured_artifact_path'] = str(outside_file)
        manifest['transition_states'][0]['captured_artifact_sha256'] = _sha256_of(str(outside_file))
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='outside'):
            verify_capture(capture_dir)

    def test_verify_refuses_a_captured_log_sha256_key_that_escapes_capture_dir_via_dotdot(self, tmp_path):
        # Same vector, but through a captured_log_sha256 KEY (a relpath, not the artifact path
        # itself) using '..' to walk back out of capture_dir, rather than an absolute path.
        #
        # The escaping relpath must resolve to a file that genuinely EXISTS and genuinely hashes
        # to the recorded value -- one directory level up from capture_dir, exactly where
        # outside_file is written below -- so that without confinement this would actually
        # "verify" (wrong file, right-looking hash), not merely raise some unrelated error. An
        # earlier version of this test used two '..' components, overshooting past tmp_path to a
        # location with no file at all; that made the test vacuous, since the resulting "no longer
        # exists on disk" ValueError coincidentally matched the pytest.raises 'outside' pattern
        # only because the filename itself (outside_log.out) appears in that unrelated message --
        # the test would have "passed" even with the confinement check deleted entirely.
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        outside_file_path = os.path.join(os.path.dirname(capture_dir), 'outside_log.out')
        with open(outside_file_path, 'wb') as f:
            f.write(b'an external log the tampered manifest points at\n')
        escaping_relpath = os.path.join('..', os.path.basename(outside_file_path))
        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        entry = manifest['transition_states'][0]
        entry['captured_log_sha256'] = {escaping_relpath: _sha256_of(outside_file_path)}
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='outside'):
            verify_capture(capture_dir)

    def test_verify_refuses_a_networks_block_captured_path_that_escapes_capture_dir(self, tmp_path):
        # The third vector called out by the finding: the networks block's own captured_path.
        record, networks, capture_dir, result = self._captured_setup(tmp_path)
        outside_file = tmp_path / 'outside_network.py'
        outside_file.write_bytes(b'an external network file the tampered manifest points at\n')
        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        with open(manifest_path) as f:
            manifest = yaml.safe_load(f)
        manifest['networks'][record.network_id]['captured_path'] = str(outside_file)
        manifest['networks'][record.network_id]['source_sha256'] = f'sha256:{_sha256_of(str(outside_file))}'
        with open(manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)

        with pytest.raises(ValueError, match='outside'):
            verify_capture(capture_dir)


REAL_NETWORK_FIXTURE = os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'iteration_1', 'RMG', 'pdep',
                                    'network4_1.py')


def test_a_capture_alone_can_drive_a_hybrid_network_write(tmp_path):
    """
    The end-to-end acceptance property this whole capture layer exists for: with BOTH the RMG pdep
    source directory and the ARC project directory deleted, the capture alone must still supply
    everything ``write_hybrid_network_input_file`` needs -- the network source, the ME method, the
    energy-reference settings, and the QM artifact (with its ``Log(...)`` files resolving inside the
    capture).

    Every other capture test asserts pieces of this in isolation: that the manifest verifies, that
    the vendored network is byte-identical, that the frozen settings round-trip. None of them proves
    the pieces actually COMPOSE into a usable hybrid write, which is the only thing the capture is
    for. This test deliberately uses the REAL network fixture rather than the stub written by
    ``_networks_for``: a stub has no ``reaction(...)`` entries, so it hashes and vendors perfectly
    while being unparseable by the writer that has to consume it -- byte-identity is not usability.
    """
    arc_dir = str(tmp_path / 'arc_project')
    os.makedirs(arc_dir)
    record = _usable_record(tmp_path=tmp_path.__class__(arc_dir), network_ts_label='TS1')

    rmg_pdep_dir = str(tmp_path / 'rmg_pdep')
    os.makedirs(rmg_pdep_dir)
    live_source_path = os.path.join(rmg_pdep_dir, f'{record.network_id}.py')
    shutil.copyfile(REAL_NETWORK_FIXTURE, live_source_path)
    networks = {record.network_id: {'source_path': live_source_path,
                                    'source_sha256': hash_file(live_source_path),
                                    'method': 'MSC'}}
    capture_dir = str(tmp_path / 'capture')
    _capture([record], arc_dir, capture_dir, networks=networks)

    shutil.rmtree(rmg_pdep_dir)
    shutil.rmtree(arc_dir)
    assert not os.path.exists(rmg_pdep_dir)
    assert not os.path.exists(arc_dir)

    manifest = read_yaml_file(path=os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME))
    energy_settings = manifest['energy_settings']
    network_entry = manifest['networks'][record.network_id]
    captured_entries = [entry for entry in manifest['transition_states']
                        if entry.get('captured_artifact_path')]
    assert len(captured_entries) == 1, captured_entries
    entry = captured_entries[0]

    # The writer consumes DOF-normalized conformer data (E0/frequencies/imaginary extracted from the
    # captured artifact through Arkane, in the real loop); this test uses a stub artifact, so it
    # hands the writer synthetic vibration-only data directly, exercising the capture->write handoff
    # without running Arkane.
    ts_conformer = {'label': entry['network_ts_label'], 'is_ts': True, 'E0_kJ_mol': -38.0,
                    'frequencies_cm_1': [500.0, 800.0, 1200.0], 'imaginary_frequency_cm_1': -1800.0,
                    'spin_multiplicity': 1, 'optical_isomers': 1, 'hindered_rotors': []}
    result = write_hybrid_network_input_file(
        source_path=os.path.join(capture_dir, network_entry['captured_path']),
        dest_path=str(tmp_path / 'hybrid' / 'input.py'),
        method=network_entry['method'],
        qm_transition_states={entry['network_ts_label']: ts_conformer},
        energy_settings=QMEnergySettings.from_frozen(energy_settings),
    )

    assert os.path.isfile(result.dest_path)
    assert result.qm_ts_labels == (entry['network_ts_label'],)
    # The point of a HYBRID network: the transition states that never got QM must survive as ILT
    # rather than being dropped. Asserting non-emptiness first keeps this from passing vacuously.
    assert len(result.ilt_ts_labels) > 0
    assert entry['network_ts_label'] not in result.ilt_ts_labels
    with open(result.dest_path, 'r') as f:
        written = f.read()
    assert 'useAtomCorrections = True' in written


class TestSensitivityAggregationMarker:
    """The 'sensitivity_aggregation' manifest marker: T3's in-run captures record nothing (their
    evidence is observable-selected and threshold-gated), the standalone PES loop marks its
    manifests 'all_directions_max_abs' -- the two conventions' admissible coefficient ranges are
    incompatible, and the marker is what prevents a silent comparison on the shared field names."""

    def test_default_writes_no_marker_and_verify_surfaces_none(self, tmp_path):
        record = _usable_record(tmp_path)
        result = _capture([record], str(tmp_path), f'{tmp_path}_capture')
        manifest = read_yaml_file(result.manifest_path)
        assert 'sensitivity_aggregation' not in manifest
        assert verify_capture(result.capture_dir).sensitivity_aggregation is None

    def test_marker_is_written_verbatim_and_surfaced_by_verify(self, tmp_path):
        record = _usable_record(tmp_path)
        result = _capture([record], str(tmp_path), f'{tmp_path}_capture',
                          sensitivity_aggregation=capture_module.SENSITIVITY_AGGREGATION_ALL_DIRECTIONS_MAX_ABS)
        manifest = read_yaml_file(result.manifest_path)
        assert manifest['sensitivity_aggregation'] == 'all_directions_max_abs'
        assert verify_capture(result.capture_dir).sensitivity_aggregation == 'all_directions_max_abs'

    def test_an_unrecognized_marker_is_refused_on_write(self, tmp_path):
        record = _usable_record(tmp_path)
        with pytest.raises(ValueError, match="Unrecognized 'sensitivity_aggregation'"):
            _capture([record], str(tmp_path), f'{tmp_path}_capture',
                     sensitivity_aggregation='made_up_convention')

    def test_a_hand_edited_marker_is_refused_by_verify(self, tmp_path):
        record = _usable_record(tmp_path)
        result = _capture([record], str(tmp_path), f'{tmp_path}_capture',
                          sensitivity_aggregation=capture_module.SENSITIVITY_AGGREGATION_ALL_DIRECTIONS_MAX_ABS)
        manifest = read_yaml_file(result.manifest_path)
        manifest['sensitivity_aggregation'] = 'tampered'
        with open(result.manifest_path, 'w') as f:
            yaml.safe_dump(manifest, f)
        with pytest.raises(ValueError, match="unrecognized 'sensitivity_aggregation'"):
            verify_capture(result.capture_dir)
