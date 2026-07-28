"""
t3 pdep capture module

Freeze, into a durable and self-contained directory, everything downstream P-dep processing
(hybrid network write, ME solve) needs about the QM transition state (TS) artifacts ARC produced
-- so that nothing downstream ever has to reach back into the live ARC project directory again.

Why this module exists: ARC deletes and recreates ``<project>/calcs/statmech/kinetics/`` on every
rate pass (``arc/scheduler.py``, ``delete_existing_subdir=True``), so the pointer ``.py`` files
under it are EPHEMERAL. They carry ARC-rendered metadata (``linear``/``spinMultiplicity``/rotor
selections, sp/freq path bindings -- see ``arc/statmech/arkane.py::species_input_template``) that
is not recoverable from the logs alone; once ARC's next pass wipes that directory, a T3 run that
deferred vendoring until later would find nothing left to vendor. So capture happens as a single,
mandatory step, run immediately after discovery resolves convergence against the still-live ARC
project directory.

There is deliberately only ONE public entry point that discovers and vendors, ``capture_ts_artifacts``,
which does discover -> capture -> freeze as one indivisible operation. There is no separate "call
capture later" helper: a caller that can forget to invoke it just relocates the very bug this
module exists to prevent. The other public function, ``verify_capture``, does not discover or
vendor anything -- it only re-checks a capture that already exists, so a consumer that did not
create the capture can still validate it before trusting it.

This module never duplicates the copying/atomic-replace machinery already hardened in
``t3.pdep.hybrid``: artifact reading (with ``Log(...)`` confinement) is ``_read_qm_artifact``, and
staged, all-or-nothing vendoring is ``_vendor_qm_artifacts``. This module only adds what those two
do not already provide: discovery-then-capture sequencing, a durable manifest, and sha256
provenance hashes of both the sources and the vendored copies.

The manifest's hashes are load-bearing, not decorative: ``capture_ts_artifacts`` hashes each
captured (vendored) file immediately after it is copied and compares it against the source hash
computed before the copy, so a copy that silently corrupted or truncated a file is caught and
raised on the spot -- never returned as a capture that disagrees with its own manifest. Those same
captured-file hashes are also recorded in the manifest so that ``verify_capture`` can later re-check
a capture using only the manifest and the capture directory, with no dependency on the (likely
long-gone) ARC project directory. This is also how a TORN capture is detected: the vendored ``qm/``
directory and the manifest are two separate objects written at different times (the manifest last),
so a crash or external edit between the two can leave them disagreeing; ``verify_capture`` refuses
rather than silently consuming whichever one it finds.

The recorded ``source_artifact_sha256`` for the pointer ``.py`` file is taken from the SAME read
``t3.pdep.hybrid._read_qm_artifact`` performs to parse it (its returned ``'sha256'`` field), never
from a second, later read of the file. Hashing via a separate read would be a TOCTOU: if the
pointer file changed on disk between the parse-read and a later hash-only read, the recorded hash
would describe bytes that were never the ones actually parsed or vendored, silently making the
manifest lie about provenance. The equivalent window for each referenced LOG file (a separate read
for its source hash, versus the read ``_vendor_qm_artifacts`` performs moments later to copy it) is
NOT closed the same way -- doing so would require reading every log file fully into memory here
just to hash it, which statmech logs are not guaranteed to be small enough to make cheap. That
narrower, log-only window is a documented residual, not an oversight.

The vendored pointer ``.py`` file is a special case for verification: it is intentionally REWRITTEN
during vendoring (``_vendor_qm_artifacts`` rewrites its ``Log(...)`` arguments to be relative to
``capture_dir``), so its captured bytes never equal its source bytes -- asserting byte-equality for
it would be asserting something false by construction. What IS verified for the pointer is a
structural property that survives the rewrite: read back the captured pointer with
``_read_qm_artifact`` (confined to ``capture_dir``) and confirm it resolves to exactly the same
number of logs as the source did, and that the content hash of each resolved captured log is one of
the already-verified, already-hashed vendored logs. In other words: not "same bytes," but "still
correctly and exclusively points at its own already-verified vendored logs."
"""

import hashlib
import os
import shutil
import tempfile
from dataclasses import dataclass

from arc.common import read_yaml_file, save_yaml_file

from t3.pdep.discovery import TSArtifactRecord, discover_ts_artifacts
from t3.pdep.hybrid import _read_qm_artifact, _vendor_qm_artifacts, _vendored_log_names
from t3.pdep.join import validate_ts_join_records

# The manifest file name, written at the top level of every capture directory.
CAPTURE_MANIFEST_FILE_NAME = 'capture_manifest.yml'

# The subdirectory (mirroring t3.pdep.hybrid's own 'qm/' convention) that vendored artifacts and
# logs are written under, inside the capture directory.
_CAPTURE_QM_SUBDIR = 'qm'

# The subdirectory that inert, non-authoritative provenance copies (of ARC's output/status.yml and
# output/output.yml) are archived under. Deliberately separate from _CAPTURE_QM_SUBDIR: nothing in
# the capture-consumption path may ever treat this directory's contents as authoritative -- see
# _capture_provenance_files and the module docstring.
_CAPTURE_PROVENANCE_SUBDIR = 'provenance'


@dataclass(frozen=True)
class CaptureResult:
    """
    The outcome of ``capture_ts_artifacts``.

    Args:
        capture_dir (str): The capture directory populated (its ``qm/`` subdirectory holds the
                           vendored artifacts/logs, its top level holds the manifest).
        manifest_path (str): The path to the written manifest file (see
                             ``CAPTURE_MANIFEST_FILE_NAME``).
        records (tuple): One ``TSArtifactRecord`` per transition state discovery inspected. For
                         every record that had an artifact (``ARTIFACT_STATUS_USABLE`` or
                         ``ARTIFACT_STATUS_UNVERIFIED``), ``artifact_path`` now points INTO
                         ``capture_dir``, never at the (ephemeral) ARC project directory. Every
                         other field is exactly what discovery resolved, frozen at capture time.
    """
    capture_dir: str
    manifest_path: str
    records: tuple


@dataclass(frozen=True)
class VerifyResult:
    """
    The outcome of ``verify_capture``.

    ``verify_capture`` raises on any invalid capture (see its docstring), so simply returning from
    it already means "this capture verified." What it CANNOT communicate through a bare, successful
    return is whether the capture is non-trivially populated or is a well-formed but entirely empty
    one (e.g. a run whose transition states were all ``not_queued``/``already_present``, so nothing
    was ever queued to ARC in the first place) -- both verify cleanly, but a caller (in particular
    the wiring guard this module exists to support) needs to be able to tell them apart WITHOUT
    relying on exception-vs-return-value ambiguity to do it. This result exposes that distinction as
    plain counts instead.

    Args:
        capture_dir (str): The capture directory that was verified.
        manifest_path (str): The manifest path that was read.
        record_count (int): The number of transition-state entries the manifest records. Always
                            >= 1 for anything ``verify_capture`` returns from -- a manifest with
                            zero entries is refused (see ``verify_capture``), never silently
                            accepted as vacuously valid.
        captured_artifact_count (int): The number of those entries that actually carry a captured
                                       artifact (``captured_artifact_path is not None``). May
                                       legitimately be zero -- that is the "verified but empty"
                                       case this result exists to make distinguishable from
                                       "verified and non-empty" for the caller.
    """
    capture_dir: str
    manifest_path: str
    record_count: int
    captured_artifact_count: int


def capture_ts_artifacts(join_records: list,
                         arc_project_directory: str,
                         capture_dir: str,
                         sensitivity_by_ts: dict | None = None,
                         supersede: bool = False,
                         ) -> CaptureResult:
    """
    Discover, capture and freeze every network transition state's QM statmech artifact, in one
    indivisible step.

    This is the ONLY entry point into this module. It runs ``discover_ts_artifacts`` against the
    LIVE ``arc_project_directory`` (resolving convergence now, while ARC's own state is still
    fresh), vendors every artifact discovery found (and the logs it references) into
    ``capture_dir``, and writes a manifest recording, per transition state, enough provenance to
    audit the capture after the fact -- all before returning. There is no way to discover without
    also capturing: a caller that could do one without the other could silently skip capture and
    only find out once ARC has already wiped the source files out from under it.

    Args:
        join_records (list): The ``t3.pdep.join.TSJoinRecord`` entries to resolve, exactly as
                             ``discover_ts_artifacts`` expects.
        arc_project_directory (str): The ARC project directory ARC just ran in. Read from, never
                                    written to.
        capture_dir (str): The directory to populate. Must be nonexistent, empty, or a directory
                           this module previously captured into (see ``_refuse_unowned_capture_dir``);
                           created if it does not already exist. Must NOT resolve inside
                           ``arc_project_directory`` (see ``_refuse_capture_dir_inside_arc_project``).
                           Repeated calls (e.g. a re-run) atomically replace its ``qm/``
                           subdirectory and its manifest; see ``_vendor_qm_artifacts`` for the
                           staging + atomic-replace guarantee that protects a prior capture from a
                           mid-vendoring failure.
        sensitivity_by_ts (dict, optional): Forwarded verbatim to ``discover_ts_artifacts``.
        supersede (bool, optional): Defaults to ``False``. If this discovery pass finds ZERO usable
                                   artifacts (``records_with_artifact`` empty) AND ``capture_dir``
                                   already holds a valid, non-empty prior capture, the default
                                   behavior is to RAISE rather than silently clear that prior
                                   capture's ``qm/`` -- a zero-artifact pass is the common, innocent
                                   case for a genuinely empty capture, but it is ALSO exactly what
                                   an accidental replay produces (e.g. against a stale or
                                   already-cleaned-up ``arc_project_directory``), and the two are
                                   indistinguishable from inside this function. Pass ``True`` only
                                   when the caller has independently confirmed this replacement is
                                   intended (e.g. a deliberate, deliberate-only re-run); it is never
                                   inferred from context. Has no effect when the existing
                                   ``capture_dir`` is empty, unowned-but-refused earlier, or itself
                                   fails ``verify_capture`` (nothing valid to protect in that case).

    Raises:
        ValueError: Anything ``discover_ts_artifacts``, ``_read_qm_artifact`` or
                   ``_vendor_qm_artifacts`` themselves raise (e.g. an artifact whose ``Log(...)``
                   path escapes ``arc_project_directory``, or references a log file that no longer
                   exists); a duplicate ``arc_ts_label`` or join key anywhere in ``join_records``
                   (via ``validate_ts_join_records``, checked BEFORE discovery runs, over ALL
                   records regardless of whether they currently have an artifact); a ``capture_dir``
                   that resolves inside ``arc_project_directory``; a ``capture_dir`` that holds
                   content this module did not put there; a record with a truthy ``artifact_path``
                   but a falsy ``arc_ts_label``; a zero-artifact pass that would destroy a valid,
                   non-empty existing capture without ``supersede=True``; or a captured file whose
                   hash does not match its source (or, for the rewritten pointer, whose resolved
                   logs' hashes do not match the vendored logs captured alongside it) -- a torn or
                   corrupted copy. All of this is raised before ``capture_dir`` is touched at all
                   where possible, and the post-copy hash checks fail closed the instant a mismatch
                   is found, so a single bad transition state cannot leave a partially-populated or
                   silently-corrupt capture directory behind, or destroy a prior, good capture.

    Returns:
        CaptureResult: The populated capture directory, its manifest path, and the frozen records.
    """
    # Pre-flight, before discovery even runs: refuse ANY duplicate arc_ts_label (or duplicate
    # (network_id, network_ts_label) key) across the FULL join_records population, regardless of
    # whether a given record currently has an artifact. Checking only records_with_artifact (as an
    # earlier version of this function did) misses the common case where one of the two colliding
    # records has no artifact YET -- e.g. one network's ARC run hasn't reached that transition state
    # -- and would let the surviving record's artifact silently claim the shared arc_ts_label with no
    # error at all today, only a confusing collision later once the second artifact appears.
    # t3.pdep.join.validate_ts_join_records already performs exactly this all-records check (over
    # both record.key and record.arc_ts_label), so it is called here rather than re-implemented.
    validate_ts_join_records(join_records)

    # Pre-flight: refuse an empty join_records outright. Capturing "nothing" is never a meaningful
    # request -- an iteration that queued no transition states must not reach this function at all
    # (its caller returns early) -- and allowing it breaks the invariant verify_capture depends on:
    # a manifest with zero 'transition_states' entries is exactly what verify_capture refuses as
    # structurally vacuous, so without this guard capture_ts_artifacts could write a capture that
    # its own verifier then rejects. That inconsistency would surface as a spurious "invalid
    # capture" in any consumer that verifies before trusting, for a capture this module itself
    # produced. Refusing here keeps "anything capture_ts_artifacts wrote, verify_capture accepts"
    # a real invariant rather than an aspiration.
    if not join_records:
        raise ValueError(f"Refusing to capture into '{capture_dir}': join_records is empty. An iteration that "
                         f"queued no P-dep transition states must not call capture_ts_artifacts at all; capturing "
                         f"zero transition states would write a manifest that verify_capture refuses as "
                         f"structurally vacuous.")

    # Pre-flight: refuse a capture_dir that resolves inside arc_project_directory. ARC deletes and
    # recreates its own subtrees (and the whole project directory can later be cleaned up) -- the
    # entire point of this module is durability against that, so a capture_dir nested inside the ARC
    # project directory would be destroyed right along with the thing it is meant to survive.
    _refuse_capture_dir_inside_arc_project(arc_project_directory=arc_project_directory, capture_dir=capture_dir)

    # Pre-flight: refuse a capture_dir that holds unrelated content this module did not put there.
    # A capture_dir is only ever allowed to be (a) nonexistent, (b) empty, or (c) a directory this
    # module previously captured into (identified by a valid capture manifest at its top level) --
    # anything else risks silently treating a caller's unrelated directory as "ours" and replacing
    # its qm/ subdirectory.
    _refuse_unowned_capture_dir(capture_dir=capture_dir)

    records = discover_ts_artifacts(
        join_records=join_records,
        arc_project_directory=arc_project_directory,
        sensitivity_by_ts=sensitivity_by_ts,
    )

    records_with_artifact = [record for record in records if record.artifact_path]

    # Pre-flight: refuse, before touching capture_dir at all, any record whose artifact_path is
    # truthy but whose arc_ts_label is falsy -- such a record would key artifact_infos/source_hashes
    # on None and hand _read_qm_artifact a None ts_label, silently corrupting the keying below.
    for record in records_with_artifact:
        if not record.arc_ts_label:
            raise ValueError(
                f"Refusing to capture transition state network_id={record.network_id!r} "
                f"network_ts_label={record.network_ts_label!r}: it has an artifact_path but no "
                f"arc_ts_label, so it cannot be captured without silently corrupting the "
                f"per-label keying capture_ts_artifacts relies on."
            )

    # Pre-flight: parse every QM artifact and resolve (and hash) every Log(...) file it references
    # BEFORE anything is written to capture_dir. Mirrors write_hybrid_network_input_file's own
    # pre-flight-before-write ordering, for the same reason: a missing/escaping log must raise
    # before any vendoring starts, never partway through it.
    artifact_infos = dict()
    source_hashes = dict()  # arc_ts_label -> {'artifact': sha256, 'logs': {original_path: sha256}}
    for record in records_with_artifact:
        info = _read_qm_artifact(
            ts_label=record.arc_ts_label,
            artifact_path=record.artifact_path,
            allowed_log_root=arc_project_directory,
        )
        artifact_infos[record.arc_ts_label] = info
        source_hashes[record.arc_ts_label] = {
            # 'artifact' is taken from info['sha256'] -- the hash _read_qm_artifact computed from
            # the SAME read it used to parse the pointer file -- rather than from a second, separate
            # read of record.artifact_path here. A separate re-read would be a TOCTOU: if the
            # pointer file changed on disk between _read_qm_artifact's read and a later hash-only
            # read, the recorded source_artifact_sha256 would describe bytes that were never the
            # ones actually parsed or vendored, silently making the manifest lie about provenance.
            'artifact': info['sha256'],
            # RESIDUAL WINDOW (documented, not closed): unlike the pointer file above, each log
            # file's source hash below IS a separate read from the copy _vendor_qm_artifacts performs
            # moments later (via shutil.copyfile), so a log file that changes on disk in that gap
            # would still go undetected. Closing this the same way as the pointer file would require
            # reading each full log file into memory here (statmech logs can be large, unlike the
            # small pointer .py file) just to hash it, then handing those same bytes to the vendoring
            # step instead of letting it stream-copy from disk -- a real memory/design cost this fix
            # does not take on. This window is intentionally left open.
            'logs': {original_path: _sha256_file(resolved_path) for original_path, resolved_path in info['logs']},
        }

    # Pre-flight: refuse to destroy a valid, non-empty EXISTING capture when this pass found ZERO
    # usable artifacts, unless the caller explicitly opts in via supersede=True. Zero artifacts is
    # the common, legitimate case for a first capture (nothing has reached ARC yet) or a genuine
    # re-run where every transition state is deliberately being cleared -- but it is ALSO exactly
    # what an accidental replay produces after ARC's tree was cleaned up between discovery passes
    # (e.g. a caller re-running capture against a stale arc_project_directory, or a bug that drops
    # every record's artifact_path). Without this guard, capture_ts_artifacts would silently accept
    # either case identically and unconditionally clear a prior capture's qm/ -- turning a real,
    # already-verified capture into an empty one with no error at all. supersede=True is how a
    # caller PROVES the replacement is intended; it is never inferred from context.
    if not records_with_artifact and not supersede and _has_valid_capture_manifest(capture_dir):
        try:
            existing = verify_capture(capture_dir)
        except ValueError:
            existing = None
        if existing is not None and existing.captured_artifact_count > 0:
            raise ValueError(
                f"Refusing to capture zero artifacts into '{capture_dir}': it already holds a valid, "
                f"non-empty capture ({existing.captured_artifact_count} captured artifact(s) across "
                f"{existing.record_count} transition state(s)), and this pass discovered no usable "
                f"artifact for any transition state. Overwriting it now would silently destroy a good "
                f"capture, indistinguishable from an accidental replay against a stale or cleaned-up "
                f"arc_project_directory. Pass supersede=True to confirm this replacement is intended."
            )

    if not os.path.isdir(capture_dir):
        os.makedirs(capture_dir)

    qm_dir = os.path.join(capture_dir, _CAPTURE_QM_SUBDIR)
    # Zero usable artifacts is the common case (missing/unusable/not_queued/unverified-with-no-path
    # records carry no artifact_path at all), not an exception to special-case: _vendor_qm_artifacts
    # is always called, even with an empty artifact_infos, so that a zero-artifact pass atomically
    # clears any stale qm/ a PREVIOUS capture into this same capture_dir left behind, keeping the
    # directory in agreement with the (now-empty) manifest. _vendor_qm_artifacts's own managed-qm_dir
    # guard (basename exactly 'qm', parent exactly dest_dir) still applies, so a foreign qm_dir is
    # never touched. (The guard above already refused this outright when it would destroy a valid,
    # non-empty prior capture without supersede=True; reaching this point means either there was
    # nothing worth protecting, or the caller explicitly opted in.)
    _vendor_qm_artifacts(artifact_infos=artifact_infos, qm_dir=qm_dir, dest_dir=capture_dir)

    captured_records = list()
    manifest_entries = list()
    for record in records:
        info = artifact_infos.get(record.arc_ts_label)
        if info is None:
            captured_records.append(record)
            manifest_entries.append(_manifest_entry(record=record, source_hashes=None, captured_artifact_path=None,
                                                     captured_log_paths=None))
            continue

        log_names = _vendored_log_names(resolved_paths=[resolved for _, resolved in info['logs']])
        captured_artifact_relpath = os.path.join(_CAPTURE_QM_SUBDIR, f'{record.arc_ts_label}.py')
        captured_log_relpaths = {
            original_path: os.path.join(_CAPTURE_QM_SUBDIR, 'logs', record.arc_ts_label, log_names[resolved_path])
            for original_path, resolved_path in info['logs']
        }

        # Verify captured-vs-source immediately after the copy, fail closed on any disagreement.
        # Logs are byte-for-byte copies (_vendor_qm_artifacts never rewrites them), so their captured
        # hash must equal the source hash computed before the copy; any mismatch means the copy
        # itself corrupted/truncated a file and this capture must never be returned as if it were
        # trustworthy.
        captured_log_sha256 = dict()
        for original_path, resolved_path in info['logs']:
            captured_relpath = captured_log_relpaths[original_path]
            captured_abs_path = os.path.join(capture_dir, captured_relpath)
            captured_hash = _sha256_file(captured_abs_path)
            source_hash = source_hashes[record.arc_ts_label]['logs'][original_path]
            if captured_hash != source_hash:
                raise ValueError(
                    f"Torn or corrupt capture detected for transition state {record.arc_ts_label!r}: "
                    f"the vendored log at '{captured_abs_path}' has sha256 {captured_hash}, but its "
                    f"source '{resolved_path}' hashed to {source_hash} before the copy. Refusing to "
                    f"return a capture whose vendored bytes disagree with their own source."
                )
            captured_log_sha256[captured_relpath] = captured_hash

        # The vendored pointer .py is intentionally REWRITTEN by _vendor_qm_artifacts (its Log(...)
        # arguments become relative to capture_dir), so its captured bytes can never equal the
        # source's -- asserting byte-equality here would be asserting something false by
        # construction. What IS verified instead is a structural property that survives the
        # rewrite: read the captured pointer back with _read_qm_artifact (confined to capture_dir)
        # and confirm it resolves to exactly as many logs as the source did, and that every one of
        # those resolved logs' content hash is one of the logs just hashed and verified above. That
        # is: not "same bytes," but "still correctly and exclusively points at its own
        # already-verified vendored logs."
        captured_artifact_abs_path = os.path.join(capture_dir, captured_artifact_relpath)
        captured_artifact_sha256 = _sha256_file(captured_artifact_abs_path)
        captured_pointer_info = _read_qm_artifact(
            ts_label=record.arc_ts_label,
            artifact_path=captured_artifact_abs_path,
            allowed_log_root=capture_dir,
        )
        if len(captured_pointer_info['logs']) != len(info['logs']):
            raise ValueError(
                f"Torn or corrupt capture detected for transition state {record.arc_ts_label!r}: the "
                f"captured pointer '{captured_artifact_abs_path}' resolves to "
                f"{len(captured_pointer_info['logs'])} log(s), but the source pointer resolved to "
                f"{len(info['logs'])}. Refusing to return a capture whose pointer no longer matches "
                f"the artifact it was vendored from."
            )
        expected_hashes = sorted(captured_log_sha256.values())
        actual_hashes = sorted(_sha256_file(resolved) for _, resolved in captured_pointer_info['logs'])
        if actual_hashes != expected_hashes:
            raise ValueError(
                f"Torn or corrupt capture detected for transition state {record.arc_ts_label!r}: the "
                f"captured pointer '{captured_artifact_abs_path}' resolves to logs whose content "
                f"hashes do not match the vendored logs captured alongside it. Refusing to return a "
                f"capture whose pointer no longer correctly and exclusively references its own "
                f"already-verified logs."
            )

        record_dict = record.as_dict()
        # CRITICAL: the raw ARC pointer path must not cross this boundary -- only the vendored,
        # capture_dir-relative path is handed back to callers downstream of capture.
        record_dict['artifact_path'] = captured_artifact_abs_path
        captured_records.append(TSArtifactRecord.from_dict(record_dict))

        manifest_entries.append(_manifest_entry(
            record=record,
            source_hashes=source_hashes[record.arc_ts_label],
            captured_artifact_path=captured_artifact_relpath,
            captured_log_paths=captured_log_relpaths,
            captured_artifact_sha256=captured_artifact_sha256,
            captured_log_sha256=captured_log_sha256,
        ))

    # Archive status.yml/output.yml as inert provenance (see _capture_provenance_files) AFTER
    # vendoring/verification above has succeeded, so a provenance-archiving failure never leaves a
    # capture_dir in a state inconsistent with what was just verified, and BEFORE the manifest is
    # written, so the manifest's 'provenance' entries always describe files that are already on
    # disk under capture_dir by the time a reader could see the manifest.
    provenance_entries = _capture_provenance_files(arc_project_directory=arc_project_directory,
                                                   capture_dir=capture_dir)

    manifest_path = _write_manifest(capture_dir=capture_dir, arc_project_directory=arc_project_directory,
                                    manifest_entries=manifest_entries, provenance_entries=provenance_entries)

    return CaptureResult(capture_dir=capture_dir, manifest_path=manifest_path, records=tuple(captured_records))


def _manifest_entry(record: TSArtifactRecord, source_hashes: dict | None, captured_artifact_path: str | None,
                    captured_log_paths: dict | None, captured_artifact_sha256: str | None = None,
                    captured_log_sha256: dict | None = None) -> dict:
    """
    Build one transition state's manifest entry.

    Deliberately does NOT record ``output/status.yml`` or ``output/output.yml`` contents inline in
    this per-transition-state entry: the frozen ``status``/``converged``/``reason`` fields below
    already carry the resolved VERDICT, which is what everything downstream actually needs, and
    remain the sole authority. ``status.yml``/``output.yml`` themselves ARE archived, but only as
    inert, non-authoritative provenance under the manifest's separate top-level ``'provenance'``
    key and ``capture_dir``'s ``provenance/`` subdirectory -- never under ``qm/``, and never
    consulted by anything that consumes a capture. See ``_capture_provenance_files`` and the
    module-level docstring for why they are archived at all (auditing the mtime-based restart
    tiebreak in ``t3.pdep.discovery.read_arc_convergence``) and why they must never become
    authoritative here.

    The ``captured_*_sha256`` fields are what makes the manifest's hashes load-bearing rather than
    decorative: they let ``verify_capture`` re-check a capture using only the manifest and the
    capture directory (no dependency on the source/ARC project directory, which may be long gone),
    and they are how a TORN capture -- ``qm/`` and the manifest disagreeing because a crash or edit
    happened between the two being written -- is detected rather than silently consumed.

    Args:
        record (TSArtifactRecord): The discovery record for this transition state.
        source_hashes (dict, optional): ``{'artifact': sha256, 'logs': {original_path: sha256}}``
                                        for the SOURCE (pre-vendoring) files, or ``None`` if this
                                        record carried no artifact to capture.
        captured_artifact_path (str, optional): The vendored artifact's path, relative to
                                                ``capture_dir``, or ``None``.
        captured_log_paths (dict, optional): Original log path -> vendored log path (relative to
                                             ``capture_dir``), or ``None``.
        captured_artifact_sha256 (str, optional): sha256 of the vendored pointer ``.py``, as it
                                                  actually reads on disk after vendoring. This is a
                                                  fingerprint only -- it is never compared to
                                                  ``source_artifact_sha256``, since the pointer is
                                                  intentionally rewritten during vendoring and will
                                                  never match its source bytes.
        captured_log_sha256 (dict, optional): Vendored log path (relative to ``capture_dir``) ->
                                              sha256, already verified equal to the corresponding
                                              source log's hash at capture time.

    Returns:
        dict: The manifest entry.
    """
    source_logs = list()
    if source_hashes is not None:
        for original_path, sha256 in source_hashes['logs'].items():
            source_logs.append({'original_path': original_path, 'sha256': sha256})
    return {
        'network_id': record.network_id,
        'network_ts_label': record.network_ts_label,
        'arc_ts_label': record.arc_ts_label,
        'status': record.status,
        'converged': record.converged,
        'reason': record.reason,
        'source_artifact_path': record.artifact_path,
        'source_artifact_sha256': source_hashes['artifact'] if source_hashes is not None else None,
        'source_logs': source_logs,
        'captured_artifact_path': captured_artifact_path,
        'captured_log_paths': captured_log_paths or dict(),
        'captured_artifact_sha256': captured_artifact_sha256,
        'captured_log_sha256': captured_log_sha256 or dict(),
    }


def _write_manifest(capture_dir: str, arc_project_directory: str, manifest_entries: list,
                    provenance_entries: list | None = None) -> str:
    """
    Write the capture manifest atomically: stage it under a temp name inside ``capture_dir``, then
    ``os.replace`` it onto ``CAPTURE_MANIFEST_FILE_NAME``, so a write failure partway through (e.g.
    a full disk) cannot leave a truncated, half-written manifest in a prior capture's place.

    Args:
        capture_dir (str): The capture directory to write the manifest into.
        arc_project_directory (str): Recorded purely as audit provenance (where this capture came
                                    from); never read back by anything that consumes the capture.
        manifest_entries (list): One dict (see ``_manifest_entry``) per transition state.
        provenance_entries (list, optional): One dict (see ``_capture_provenance_files``) per
                                            archived ``output/status.yml``/``output/output.yml``
                                            file, or ``None``. Recorded under the manifest's
                                            ``'provenance'`` key. INERT: these entries describe
                                            files archived purely for post-hoc human audit of the
                                            status/output mtime tiebreak (see
                                            ``_capture_provenance_files`` and the module docstring);
                                            nothing in the capture-consumption path may ever read
                                            this key back as authority. The frozen ``status``/
                                            ``converged``/``reason`` fields already written per
                                            transition state above remain the sole authority.

    Returns:
        str: The path the manifest was written to.
    """
    manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
    content = {
        'arc_project_directory': arc_project_directory,
        'transition_states': manifest_entries,
        # INERT provenance only -- never read back as authority by anything in the
        # capture-consumption path. See _capture_provenance_files and the module docstring.
        'provenance': provenance_entries or list(),
    }

    fd, staged_path = tempfile.mkstemp(prefix='.capture-manifest-', dir=capture_dir)
    os.close(fd)
    try:
        save_yaml_file(path=staged_path, content=content)
        os.replace(staged_path, manifest_path)
    finally:
        if os.path.isfile(staged_path):
            os.remove(staged_path)
    return manifest_path


def _capture_provenance_files(arc_project_directory: str, capture_dir: str) -> list:
    """
    Archive ARC's ``output/status.yml`` and ``output/output.yml`` as INERT provenance.

    ``t3.pdep.discovery.read_arc_convergence`` uses ``getmtime(status.yml) > getmtime(output.yml)``
    as a restart tiebreak when the two files disagree about convergence -- a decision this module's
    caller (``capture_ts_artifacts``) has already resolved into the frozen ``status``/``converged``/
    ``reason`` fields recorded per transition state in the manifest (see ``_manifest_entry``). Those
    frozen fields remain the ONLY authority anything downstream of capture may consult; nothing in
    the capture-consumption path (``verify_capture``, or any future wiring guard) may ever read the
    files this function archives back as authority for anything. They are archived here purely so a
    human (or a future investigation) can audit, after the fact, what ``status.yml``/``output.yml``
    actually said -- and, critically, WHEN each was last written -- at the moment of capture, since
    without that a post-hoc audit of a restart-tiebreak decision would be impossible once ARC's own
    project directory is gone.

    Deliberately archived under ``_CAPTURE_PROVENANCE_SUBDIR`` ('provenance/'), a directory
    entirely separate from ``_CAPTURE_QM_SUBDIR`` ('qm/'): nothing that reads ``qm/`` as "the
    capture" should ever notice these files exist, and nothing named 'provenance' should ever look
    authoritative.

    ``shutil.copyfile`` does not preserve source mtimes, and copy-file mtimes would be worthless for
    auditing the tiebreak anyway (they would all read as "now", collapsing the very distinction the
    tiebreak depends on) -- so each source file's mtime is recorded explicitly, as a manifest VALUE,
    alongside its sha256, rather than relied upon to survive the copy.

    Either or both files may legitimately not exist (e.g. a transition state resolved from ARC
    without ever going through the ``output/`` finalization step); absence is recorded, not raised.

    Args:
        arc_project_directory (str): The ARC project directory to read ``output/status.yml`` and
                                     ``output/output.yml`` from.
        capture_dir (str): The capture directory being populated.

    Returns:
        list: One dict per file (``status.yml`` then ``output.yml``), each with keys
             ``source_path``, ``captured_path`` (relative to ``capture_dir``, or ``None`` if the
             source file did not exist), ``sha256`` (or ``None``), and ``source_mtime`` (the
             source file's ``os.path.getmtime`` result at archive time, or ``None``).
    """
    provenance_dir = os.path.join(capture_dir, _CAPTURE_PROVENANCE_SUBDIR)
    entries = list()
    for file_name in ('status.yml', 'output.yml'):
        source_path = os.path.join(arc_project_directory, 'output', file_name)
        if not os.path.isfile(source_path):
            entries.append({'source_path': source_path, 'captured_path': None, 'sha256': None,
                            'source_mtime': None})
            continue

        # Record the source mtime BEFORE the copy, for the same reason the pointer-file hash is
        # taken from the read that parsed it (see the module docstring): the value that matters is
        # the one describing the bytes actually archived, not a later re-stat that could observe a
        # file ARC has since rewritten.
        source_mtime = os.path.getmtime(source_path)
        source_sha256 = _sha256_file(source_path)

        if not os.path.isdir(provenance_dir):
            os.makedirs(provenance_dir)
        captured_relpath = os.path.join(_CAPTURE_PROVENANCE_SUBDIR, file_name)
        captured_abs_path = os.path.join(capture_dir, captured_relpath)
        fd, staged_path = tempfile.mkstemp(prefix=f'.provenance-{file_name}-', dir=provenance_dir)
        os.close(fd)
        try:
            shutil.copyfile(source_path, staged_path)
            os.replace(staged_path, captured_abs_path)
        finally:
            if os.path.isfile(staged_path):
                os.remove(staged_path)

        entries.append({'source_path': source_path, 'captured_path': captured_relpath,
                        'sha256': source_sha256, 'source_mtime': source_mtime})
    return entries


def _sha256_file(path: str) -> str:
    """
    Compute the sha256 hex digest of a file's contents.

    Args:
        path (str): The file path to hash.

    Returns:
        str: The hex digest.
    """
    hasher = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            hasher.update(chunk)
    return hasher.hexdigest()


def _refuse_capture_dir_inside_arc_project(arc_project_directory: str, capture_dir: str) -> None:
    """
    Refuse a ``capture_dir`` that resolves inside ``arc_project_directory``.

    Mirrors the ``os.path.realpath`` + ``os.path.commonpath`` confinement idiom used elsewhere in
    this module family (``t3.pdep.discovery._confine_to_project``, ``t3.pdep.hybrid._read_qm_artifact``),
    just inverted: those refuse a path that ESCAPES a root; this refuses a path that is NESTED
    inside one. Using ``realpath`` on both sides means a symlink cannot be used to make a
    capture_dir that is really inside the ARC project directory look like it is outside.

    Args:
        arc_project_directory (str): The ARC project directory, which ARC deletes/recreates
                                     subtrees of (and which may itself be cleaned up later).
        capture_dir (str): The directory this capture is about to populate.

    Raises:
        ValueError: If ``capture_dir`` resolves to ``arc_project_directory`` itself or to anything
                   inside it.
    """
    resolved_arc_project = os.path.realpath(arc_project_directory)
    resolved_capture_dir = os.path.realpath(capture_dir)
    if resolved_capture_dir == resolved_arc_project \
            or os.path.commonpath([resolved_arc_project, resolved_capture_dir]) == resolved_arc_project:
        raise ValueError(
            f"Refusing to capture into '{capture_dir}' (resolved to '{resolved_capture_dir}'): it is "
            f"inside the ARC project directory '{arc_project_directory}' (resolved to "
            f"'{resolved_arc_project}'). ARC deletes and recreates its own subtrees on every rate "
            f"pass, and the whole project directory can later be cleaned up -- the entire point of "
            f"this module is durability against exactly that, so a capture_dir nested inside "
            f"arc_project_directory would be destroyed right along with the thing it exists to "
            f"survive. Choose a capture_dir outside the ARC project directory."
        )


def _has_valid_capture_manifest(capture_dir: str) -> bool:
    """
    Whether ``capture_dir`` already contains a manifest this module could have written.

    Used only to distinguish "a directory this module previously captured into" (safe to treat as
    managed/replaceable) from "a directory holding unrelated content" (refused). Deliberately does
    NOT hash-verify the manifest's contents here -- that is ``verify_capture``'s job, and requiring
    it here would make re-capturing into a torn capture impossible to ever repair.

    Args:
        capture_dir (str): The directory to inspect.

    Returns:
        bool: True if a readable, well-formed capture manifest is present at the top level of
             ``capture_dir``.
    """
    manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
    if not os.path.isfile(manifest_path):
        return False
    try:
        content = read_yaml_file(manifest_path)
    except Exception:
        return False
    return isinstance(content, dict) and 'transition_states' in content and 'arc_project_directory' in content


def _refuse_unowned_capture_dir(capture_dir: str) -> None:
    """
    Refuse a ``capture_dir`` that holds content this module did not put there.

    A ``capture_dir`` is only ever allowed to be (a) nonexistent, (b) empty, or (c) a directory this
    module previously captured into. Anything else is refused up front, before any I/O, rather than
    unconditionally treating ``<capture_dir>/qm`` as managed/replaceable -- silently doing so would
    let a caller-supplied directory holding unrelated content have part of it deleted and replaced.

    Args:
        capture_dir (str): The directory this capture is about to populate.

    Raises:
        ValueError: If ``capture_dir`` exists, is non-empty, and does not contain a valid capture
                   manifest.
    """
    if not os.path.isdir(capture_dir):
        return
    if not os.listdir(capture_dir):
        return
    if _has_valid_capture_manifest(capture_dir):
        return
    raise ValueError(
        f"Refusing to capture into '{capture_dir}': it already exists, is non-empty, and does not "
        f"contain a capture manifest this module could have written. A capture_dir must be a "
        f"dedicated, per-capture directory -- either nonexistent, empty, or one this module "
        f"previously captured into -- never a directory holding unrelated content, since this "
        f"function treats its 'qm' subdirectory as freely replaceable."
    )


def verify_capture(capture_dir: str) -> VerifyResult:
    """
    Re-read a capture's manifest and confirm every captured file still matches its recorded hash.

    This is the counterpart a CONSUMER (not the code that created the capture) calls before
    trusting it: it re-derives every check from the manifest and the files actually on disk in
    ``capture_dir``, with no dependency on ``arc_project_directory`` (which may no longer exist by
    the time this runs). It also detects a TORN capture: the vendored ``qm/`` directory and the
    manifest are two separate objects, written at different times (the manifest last, via
    ``_write_manifest``'s atomic replace) -- a crash or an out-of-band edit between the two, or
    after either one, can leave them disagreeing, and this function refuses that rather than
    silently consuming whichever one it finds.

    A manifest that is structurally absent, empty, or unparseable is refused (raised on) rather
    than treated as vacuously "verified": a ``transition_states`` list with zero entries has
    nothing for the loop below to check, so silently returning success in that case would let an
    empty or truncated manifest "pass" verification with no work done at all -- exactly the
    fail-open failure mode this function exists to prevent. Note this is stricter than what
    ``capture_ts_artifacts`` itself allows: a WELL-FORMED zero-artifact capture (every transition
    state ``not_queued``/``already_present``, nothing ever queued into ARC) still writes at least
    one manifest entry per transition state it was asked about (with ``captured_artifact_path``
    ``None``), so ``transition_states`` being non-empty is a real invariant of anything
    ``capture_ts_artifacts`` ever wrote, not an accident of the happy path.

    Args:
        capture_dir (str): The capture directory to verify (as previously produced by
                           ``capture_ts_artifacts``).

    Returns:
        VerifyResult: Counts describing what was verified, so a caller can distinguish "verified
                     and non-empty" from anything else without relying on exception-vs-return-value
                     ambiguity. ``record_count`` is always >= 1 (a manifest with zero entries is
                     refused, see above).

    Raises:
        ValueError: If the manifest is missing, malformed, or structurally empty (no transition
                   state entries at all); a captured file the manifest references is missing from
                   disk; or a captured file's actual sha256 does not match the hash recorded for it
                   in the manifest.
    """
    manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
    if not os.path.isfile(manifest_path):
        raise ValueError(f"Refusing to verify '{capture_dir}': no capture manifest found at '{manifest_path}'.")
    try:
        content = read_yaml_file(manifest_path)
    except Exception as e:
        raise ValueError(f"Refusing to verify '{capture_dir}': the manifest at '{manifest_path}' could not be "
                         f"read ({e}).")
    if not isinstance(content, dict) or 'transition_states' not in content:
        raise ValueError(f"Refusing to verify '{capture_dir}': the manifest at '{manifest_path}' is malformed "
                         f"(expected a mapping with a 'transition_states' key).")

    transition_states = content['transition_states']
    if not transition_states:
        raise ValueError(
            f"Refusing to verify '{capture_dir}': the manifest at '{manifest_path}' has a "
            f"'transition_states' list that is empty (or missing/None). A manifest ever written by "
            f"capture_ts_artifacts records at least one entry per transition state it was asked "
            f"about, even when none carried an artifact -- an empty list here means the manifest is "
            f"structurally vacuous (truncated, hand-edited, or never really produced by a capture),"
            f" and treating it as 'verified' would silently check nothing at all."
        )

    record_count = 0
    captured_artifact_count = 0
    for entry in transition_states:
        record_count += 1
        ts_label = entry.get('arc_ts_label')
        captured_artifact_path = entry.get('captured_artifact_path')
        if captured_artifact_path is not None:
            captured_artifact_count += 1
            _verify_one_captured_file(
                capture_dir=capture_dir,
                relpath=captured_artifact_path,
                expected_sha256=entry.get('captured_artifact_sha256'),
                ts_label=ts_label,
            )
        for log_relpath, expected_log_sha256 in (entry.get('captured_log_sha256') or dict()).items():
            _verify_one_captured_file(
                capture_dir=capture_dir,
                relpath=log_relpath,
                expected_sha256=expected_log_sha256,
                ts_label=ts_label,
            )

    return VerifyResult(capture_dir=capture_dir, manifest_path=manifest_path, record_count=record_count,
                        captured_artifact_count=captured_artifact_count)


def _verify_one_captured_file(capture_dir: str, relpath: str, expected_sha256: str | None, ts_label) -> None:
    """
    Verify one captured file named in a manifest still exists and matches its recorded hash.

    Args:
        capture_dir (str): The capture directory this file is relative to.
        relpath (str): The captured file's path, relative to ``capture_dir``.
        expected_sha256 (str, optional): The hash recorded for it in the manifest.
        ts_label: The transition state this file belongs to, for error messages only.

    Raises:
        ValueError: If the file is missing, the manifest recorded no hash to check against, or the
                   file's actual hash does not match the recorded one.
    """
    abs_path = os.path.join(capture_dir, relpath)
    if not os.path.isfile(abs_path):
        raise ValueError(
            f"Torn capture detected in '{capture_dir}': the manifest for transition state "
            f"{ts_label!r} records a captured file at '{relpath}', but it no longer exists on disk. "
            f"The vendored qm/ directory and the manifest are written as two separate objects; this "
            f"disagreement means the capture was interrupted or altered after the fact. Refusing to "
            f"use it."
        )
    if not expected_sha256:
        raise ValueError(
            f"Torn capture detected in '{capture_dir}': the manifest for transition state "
            f"{ts_label!r} records captured file '{relpath}' with no hash to verify it against. "
            f"Refusing to trust an unverifiable capture."
        )
    actual_sha256 = _sha256_file(abs_path)
    if actual_sha256 != expected_sha256:
        raise ValueError(
            f"Torn or tampered capture detected in '{capture_dir}': transition state {ts_label!r}'s "
            f"captured file '{relpath}' has sha256 {actual_sha256}, but the manifest records "
            f"{expected_sha256}. Refusing to use it."
        )
