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
import math
import os
import shutil
import tempfile
from dataclasses import dataclass

from arc.common import read_yaml_file, save_yaml_file

from t3.pdep.cache import hash_file
from t3.pdep.discovery import (TS_ARTIFACT_STATUSES,
                               TS_ARTIFACT_STATUSES_REQUIRING_ARTIFACT_PATH,
                               TS_ARTIFACT_STATUSES_REQUIRING_NO_ARTIFACT_PATH,
                               TSArtifactRecord,
                               discover_ts_artifacts)
from t3.pdep.energy_settings import read_arc_energy_settings, validate_frozen_energy_settings
from t3.pdep.hybrid import _read_qm_artifact, _vendor_qm_artifacts, _vendored_log_names
from t3.pdep.join import validate_networks_block, validate_ts_join_records

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

# The subdirectory the vendored RMG network source files (one ``<network_id>.py`` per network the
# capture's transition states belong to) are written under, inside the capture directory. Unlike
# 'provenance/', this directory IS authoritative: it holds the exact network file the SA/selection
# pass examined, verified by hash against the join sidecar's frozen identity, and it is what
# downstream hybrid-network writing must read -- never the (long-gone) RMG ``pdep/`` directory.
_CAPTURE_NETWORKS_SUBDIR = 'networks'


def captured_qm_artifact_path(capture_dir: str, arc_ts_label: str) -> str:
    """
    Where ``capture_ts_artifacts`` vendors the QM artifact for one transition state.

    This layout -- ``<capture_dir>/qm/<arc_ts_label>.py`` -- is a stable invariant of this module
    (``_CAPTURE_QM_SUBDIR``; ``verify_capture`` re-verifies it, and ``t3.pdep.pes_qm`` already
    relies on it to recover an adopted artifact's own capture directory). Naming it as a function
    lets a caller that only knows a converged transition state's label (e.g.
    ``t3.pdep.pes_loop``, carrying QM artifacts across round boundaries) derive the vendored
    path without re-encoding the layout.

    Args:
        capture_dir (str): The capture directory the artifact was vendored into.
        arc_ts_label (str): The ARC transition state label (``t3.pdep.join.arc_ts_label``).

    Returns:
        str: The vendored artifact's path.
    """
    return os.path.join(capture_dir, _CAPTURE_QM_SUBDIR, f'{arc_ts_label}.py')


# The value the standalone PES loop records under a manifest's 'sensitivity_aggregation' key:
# its per-TS evidence is the maximum-|coefficient| row across EVERY network reaction direction
# key and every condition of the Arkane ME SA output, ungated (see t3.pdep.pes_sa). T3's in-run
# path records nothing under this key: its evidence is taken from ONE observable-selected
# direction key and has already passed the selector's relative/absolute gates. The two are NOT
# comparable numbers, and this marker is what keeps a loop-written manifest from ever being
# silently compared against a T3-written one on the same field name.
SENSITIVITY_AGGREGATION_ALL_DIRECTIONS_MAX_ABS = 'all_directions_max_abs'

# Every marker value this module will write or verify. An advisory field is still not a free-text
# field: a value outside this set reaching a consumer means a hand-edited or newer-format
# manifest, and both the write path and verify_capture refuse it rather than pass it through --
# the same fail-closed stance verify_capture already takes on an unrecognized artifact status.
SENSITIVITY_AGGREGATIONS = (SENSITIVITY_AGGREGATION_ALL_DIRECTIONS_MAX_ABS,)

# Sibling-of-capture_dir naming convention for the atomic-swap machinery in
# _capture_ts_artifacts_locked and its crash-recovery counterpart, _recover_capture_swap_state.
# Both live in the same parent_dir as capture_dir itself (never inside it), so a fresh capture's
# own qm/networks/provenance layout is never confused with this module's own scratch directories.
_CAPTURE_STAGING_DIR_PREFIX = '.capture-staging-'
_CAPTURE_OLD_DIR_PREFIX = '.capture-old-'


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
        energy_settings (dict, optional): The frozen, AUTHORITATIVE energy-reference settings block
                                          (see ``t3.pdep.energy_settings.read_arc_energy_settings``),
                                          or ``None`` if this capture found zero usable artifacts (no
                                          transition state had anything to freeze a model chemistry
                                          for). Also recorded under the manifest's ``'energy_settings'``
                                          key -- unlike ``'provenance'``, this key IS authoritative.
        networks (dict, optional): The frozen, AUTHORITATIVE network identity block:
                                   ``network_id -> {source_path, source_sha256, captured_path,
                                   method}``, one entry per network the captured transition states
                                   belong to, with ``captured_path`` relative to ``capture_dir``.
                                   Also recorded under the manifest's ``'networks'`` key. Like
                                   ``energy_settings`` (and unlike ``'provenance'``), this IS
                                   authoritative: downstream consumers must read the capture's own
                                   vendored network copy and recorded ME method, never re-derive
                                   them from the ARC project or the RMG ``pdep/`` directory.
    """
    capture_dir: str
    manifest_path: str
    records: tuple
    energy_settings: dict | None = None
    networks: dict | None = None


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
        networks (dict, optional): The manifest's AUTHORITATIVE ``'networks'`` block, re-verified:
                                   ``network_id -> {source_path, source_sha256, captured_path,
                                   method}`` with ``captured_path`` relative to ``capture_dir``.
                                   Surfaced here so a consumer can reach the vendored network
                                   source and its frozen ME method from the verified capture alone,
                                   with no dependency on the RMG ``pdep/`` or ARC project
                                   directories still existing.
        energy_settings (dict, optional): The manifest's AUTHORITATIVE ``'energy_settings'`` block
                                          (``None`` only for a verified zero-artifact capture),
                                          surfaced for the same reason as ``networks``.
        sensitivity_aggregation (str, optional): The manifest's ``'sensitivity_aggregation'``
                                                 marker, or ``None`` for a manifest written under
                                                 T3's in-run convention (see
                                                 ``capture_ts_artifacts``). Surfaced so a consumer
                                                 comparing sensitivity evidence across captures can
                                                 refuse to compare incompatible aggregations.
        ts_records (tuple): One ``TSArtifactRecord`` per transition-state manifest entry, built
                            from EXACTLY the fields this function itself just validated (including
                            the status/captured_artifact_path consistency check below) -- never a
                            second, independent read of the manifest. ``artifact_path`` on each
                            record is already resolved to an absolute path under ``capture_dir``
                            (confined via the same per-file check used to verify the artifact
                            exists), so a consumer building on ``ts_records`` cannot reconstruct an
                            unconfined path itself. A caller (in particular the hybrid-input writer)
                            must consume this instead of re-reading the manifest: a second,
                            untrusted read could observe a manifest mutated out from under the first
                            read (e.g. a concurrent re-capture), reintroducing exactly the
                            check-then-use race this function's own re-verification exists to close.
    """
    capture_dir: str
    manifest_path: str
    record_count: int
    captured_artifact_count: int
    networks: dict | None = None
    energy_settings: dict | None = None
    ts_records: tuple = ()
    sensitivity_aggregation: str | None = None


def capture_ts_artifacts(join_records: list,
                         arc_project_directory: str,
                         capture_dir: str,
                         networks: dict,
                         sensitivity_by_ts: dict | None = None,
                         supersede: bool = False,
                         sensitivity_aggregation: str | None = None,
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
                           The ENTIRE new capture (``qm/``, ``networks/``, ``provenance/`` and the
                           manifest) is built in a sibling staging directory and self-checked with
                           ``verify_capture`` before ``capture_dir`` is touched at all; only then is
                           it atomically swapped into ``capture_dir``'s place (a same-filesystem
                           rename-aside-then-``os.replace``, mirroring ``_vendor_qm_artifacts``'s own
                           staging + atomic-replace idiom, extended to the whole directory). A
                           repeated call (e.g. a re-run) therefore either fully replaces a prior
                           capture with a complete, self-verified new one, or -- on any failure
                           anywhere in discovery, vendoring, or verification -- leaves the prior
                           capture completely untouched; it can never observe or leave behind a torn
                           mix of old and new files. This whole call holds an interprocess lock on
                           ``capture_dir`` (see ``_acquire_capture_lock``), so two concurrent
                           processes capturing the same ``capture_dir`` can never interleave their
                           staging and swap. Honest limit: POSIX has no single syscall that replaces
                           a non-empty directory's contents atomically, so the rename-aside and the
                           ``os.replace`` above are still two separate steps; a process that dies
                           between them (``kill -9``, OOM-kill, power loss -- anything that skips
                           this function's own ``except``/``finally`` handlers) leaves ``capture_dir``
                           transiently ABSENT with the last good capture sitting under a sibling
                           ``.capture-old-*`` name. What IS guaranteed is that the very next call
                           against the same ``capture_dir`` (which must take the same lock first)
                           finds and restores it before doing anything else, via
                           ``_recover_capture_swap_state`` -- so no capture is ever lost, but a
                           concurrent reader that bypasses this module's lock could transiently
                           observe ``capture_dir`` missing (never corrupt) until that recovery runs.
        networks (dict): The join sidecar's frozen network identity block:
                        ``network_id -> {source_path, source_sha256, method}``, covering exactly
                        the networks ``join_records`` reference (see
                        ``t3.pdep.join.validate_networks_block``). Every named network's LIVE file
                        at ``source_path`` is re-hashed with ``t3.pdep.cache.hash_file`` and must
                        still match the recorded ``source_sha256`` -- a mismatch means the RMG
                        network file changed between the SA/selection pass and this capture, and
                        is refused rather than vendored. Each verified file is then copied into
                        ``capture_dir``'s ``networks/`` subdirectory and recorded under the
                        manifest's AUTHORITATIVE ``'networks'`` key.
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
        sensitivity_aggregation (str, optional): How the per-transition-state
                                                ``coefficient``/``delta_ln_k`` evidence in
                                                ``sensitivity_by_ts`` was aggregated, recorded
                                                verbatim under the manifest's
                                                ``'sensitivity_aggregation'`` key when given.
                                                ``None`` (the default, and T3's in-run
                                                convention) records nothing: that evidence comes
                                                from one observable-selected direction key and has
                                                passed the selector's gates. The standalone PES
                                                loop passes
                                                ``SENSITIVITY_AGGREGATION_ALL_DIRECTIONS_MAX_ABS``
                                                -- an ungated max over ALL direction keys -- whose
                                                values are NOT comparable to the in-run ones; the
                                                marker is what prevents a silent comparison.

    Raises:
        RuntimeError: If ``capture_dir``'s interprocess lock is already held by another (live, or
                     unidentifiable) process -- this function fails CLOSED rather than silently
                     proceeding without exclusion; see ``_acquire_capture_lock``.
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
                   silently-corrupt capture directory behind, or destroy a prior, good capture. All
                   of this happens against a staging directory rather than ``capture_dir`` itself
                   (see the ``capture_dir`` Args entry above), and a staged capture that fails its
                   own ``verify_capture`` self-check is refused the same way -- ``capture_dir`` is
                   never replaced with (or left as) anything less than a complete, self-verified
                   capture.

    Returns:
        CaptureResult: The populated capture directory, its manifest path, and the frozen records.
    """
    if sensitivity_aggregation is not None and sensitivity_aggregation not in SENSITIVITY_AGGREGATIONS:
        raise ValueError(f"Unrecognized 'sensitivity_aggregation' value {sensitivity_aggregation!r}; "
                         f"expected None (T3's unmarked in-run convention) or one of "
                         f"{sorted(SENSITIVITY_AGGREGATIONS)}. A free-text marker would defeat the "
                         f"comparison guard this field exists to provide.")
    lock_path = _acquire_capture_lock(capture_dir=capture_dir)
    try:
        # Recover from a crash inside a PRIOR call's atomic swap (see _recover_capture_swap_state)
        # before doing anything else. Safe to run unconditionally, now that this process holds the
        # lock above: no other process can be mid-build against capture_dir while this runs, so
        # there is no risk of "recovering" out from under a concurrent, still-in-flight capture.
        _recover_capture_swap_state(capture_dir=capture_dir)
        return _capture_ts_artifacts_locked(
            join_records=join_records,
            arc_project_directory=arc_project_directory,
            capture_dir=capture_dir,
            networks=networks,
            sensitivity_by_ts=sensitivity_by_ts,
            supersede=supersede,
            sensitivity_aggregation=sensitivity_aggregation,
        )
    finally:
        _release_capture_lock(lock_path)


def _capture_ts_artifacts_locked(join_records: list,
                                 arc_project_directory: str,
                                 capture_dir: str,
                                 networks: dict,
                                 sensitivity_by_ts: dict | None,
                                 supersede: bool,
                                 sensitivity_aggregation: str | None,
                                 ) -> CaptureResult:
    """
    The body of ``capture_ts_artifacts``, run only once the caller already holds ``capture_dir``'s
    interprocess lock and crash recovery has already run. Split out purely so those two concerns
    (see ``_acquire_capture_lock`` and ``_recover_capture_swap_state``) can wrap this body in a
    ``try``/``finally`` without reindenting it; see ``capture_ts_artifacts`` for the full contract,
    Args, Returns and Raises.
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

    # Pre-flight: validate the networks identity block against the join records (both directions,
    # plus per-entry field validation), then verify every named network's LIVE file still hashes to
    # the sidecar's recorded source_sha256 -- all BEFORE discovery runs and before anything is
    # written to capture_dir. A network file that changed (or vanished) between the SA/selection
    # pass and this capture is a real "the world changed under us" defect: vendoring the new bytes
    # would freeze a network the selection never actually examined.
    validate_networks_block(referenced_network_ids={record.network_id for record in join_records},
                            networks=networks,
                            context=f"the capture request for '{capture_dir}'")
    _preflight_network_sources(networks=networks)

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

    # Pre-flight: freeze the AUTHORITATIVE energy-reference settings block (model chemistry, atom/
    # bond corrections, frequency scale factor) from ARC's OWN statmech input directives, before any
    # vendoring starts -- mirrors the artifact/log pre-flight below for the same reason: a capture
    # that cannot freeze a complete energy-reference is worthless to downstream hybrid-network
    # writing regardless of whether its QM artifacts parsed cleanly, so this must fail before
    # touching capture_dir at all, never partway through. Only attempted when at least one record
    # actually carries an artifact: a zero-artifact capture (nothing ever reached ARC) has no energy
    # reference to freeze, and read_arc_energy_settings is fail-closed on ARC files that a
    # legitimately artifact-free capture has no reason to expect exist.
    energy_settings = None
    if records_with_artifact:
        energy_settings = read_arc_energy_settings(arc_project_directory, statmech_subdir='kinetics')

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

    # Build the ENTIRE new capture in a sibling staging directory, never touching capture_dir
    # itself until the staged tree is complete and has passed verify_capture's own self-check.
    # This is the directory-level analogue of _vendor_qm_artifacts's staging + atomic-replace
    # idiom (see t3.pdep.hybrid), extended to cover the whole capture rather than just its qm/
    # subdirectory: qm/, networks/, provenance/ and the manifest are all staged together, so a
    # crash or exception at ANY point below leaves capture_dir either completely absent (first
    # capture) or exactly as it was before this call (a re-capture) -- never a torn mix of some
    # new files and some stale ones. Without this, a partial write (e.g. a corrupted mid-copy that
    # _vendor_network_sources's rmtree-then-rebuild already began) could destroy a previously good,
    # already-verified capture and leave an incomplete one in its place, with no way back.
    parent_dir = os.path.dirname(os.path.abspath(capture_dir))
    if not os.path.isdir(parent_dir):
        os.makedirs(parent_dir)
    staging_dir = tempfile.mkdtemp(prefix=_CAPTURE_STAGING_DIR_PREFIX, dir=parent_dir)
    try:
        qm_dir = os.path.join(staging_dir, _CAPTURE_QM_SUBDIR)
        # Zero usable artifacts is the common case (missing/unusable/not_queued/unverified-with-no-
        # path records carry no artifact_path at all), not an exception to special-case:
        # _vendor_qm_artifacts is always called, even with an empty artifact_infos, so that a
        # zero-artifact pass still produces a (trivially empty) qm/ in the staged capture. Staging
        # is always fresh here, so _vendor_qm_artifacts's own managed-qm_dir/foreign-qm_dir branch
        # never triggers -- qm_dir cannot pre-exist inside a directory this function just created.
        _vendor_qm_artifacts(artifact_infos=artifact_infos, qm_dir=qm_dir, dest_dir=staging_dir)

        captured_records = list()
        manifest_entries = list()
        # (record, captured_artifact_relpath) pairs needing their final, capture_dir-rooted
        # artifact_path filled in AFTER the atomic swap below -- staging_dir will no longer exist
        # by then (it becomes capture_dir), so an absolute path built against it now would dangle.
        pending_records = list()
        for record in records:
            info = artifact_infos.get(record.arc_ts_label)
            if info is None:
                pending_records.append((record, None))
                manifest_entries.append(_manifest_entry(record=record, source_hashes=None,
                                                         captured_artifact_path=None, captured_log_paths=None))
                continue

            log_names = _vendored_log_names(resolved_paths=[resolved for _, resolved in info['logs']])
            captured_artifact_relpath = os.path.join(_CAPTURE_QM_SUBDIR, f'{record.arc_ts_label}.py')
            captured_log_relpaths = {
                original_path: os.path.join(_CAPTURE_QM_SUBDIR, 'logs', record.arc_ts_label,
                                            log_names[resolved_path])
                for original_path, resolved_path in info['logs']
            }

            # Verify captured-vs-source immediately after the copy, fail closed on any disagreement.
            # Logs are byte-for-byte copies (_vendor_qm_artifacts never rewrites them), so their
            # captured hash must equal the source hash computed before the copy; any mismatch means
            # the copy itself corrupted/truncated a file and this capture must never be returned as
            # if it were trustworthy.
            captured_log_sha256 = dict()
            for original_path, resolved_path in info['logs']:
                captured_relpath = captured_log_relpaths[original_path]
                captured_abs_path = os.path.join(staging_dir, captured_relpath)
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

            # The vendored pointer .py is intentionally REWRITTEN by _vendor_qm_artifacts (its
            # Log(...) arguments become relative to its own final directory), so its captured bytes
            # can never equal the source's -- asserting byte-equality here would be asserting
            # something false by construction. What IS verified instead is a structural property
            # that survives the rewrite: read the captured pointer back with _read_qm_artifact
            # (confined to staging_dir, standing in for the eventual capture_dir) and confirm it
            # resolves to exactly as many logs as the source did, and that every one of those
            # resolved logs' content hash is one of the logs just hashed and verified above. That
            # is: not "same bytes," but "still correctly and exclusively points at its own
            # already-verified vendored logs."
            captured_artifact_abs_path = os.path.join(staging_dir, captured_artifact_relpath)
            captured_artifact_sha256 = _sha256_file(captured_artifact_abs_path)
            captured_pointer_info = _read_qm_artifact(
                ts_label=record.arc_ts_label,
                artifact_path=captured_artifact_abs_path,
                allowed_log_root=staging_dir,
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

            pending_records.append((record, captured_artifact_relpath))
            manifest_entries.append(_manifest_entry(
                record=record,
                source_hashes=source_hashes[record.arc_ts_label],
                captured_artifact_path=captured_artifact_relpath,
                captured_log_paths=captured_log_relpaths,
                captured_artifact_sha256=captured_artifact_sha256,
                captured_log_sha256=captured_log_sha256,
            ))

        # Vendor the (already pre-flight-verified) network source files into networks/, after the
        # TS artifact vendoring/verification above, for the same all-verified-content-before-
        # provenance sequencing: each copy is re-hashed against the sidecar's recorded
        # source_sha256 immediately, failing closed on any disagreement, and the manifest below is
        # only ever written once every vendored network already sits verified on disk. This runs
        # against staging_dir, which _vendor_network_sources always finds networks/-less (a fresh
        # staging directory), so its rebuild-from-scratch rmtree never has anything stale to clear.
        manifest_networks = _vendor_network_sources(networks=networks, capture_dir=staging_dir)

        # Archive status.yml/output.yml as inert provenance (see _capture_provenance_files) AFTER
        # vendoring/verification above has succeeded, so a provenance-archiving failure never
        # leaves a capture_dir in a state inconsistent with what was just verified, and BEFORE the
        # manifest is written, so the manifest's 'provenance' entries always describe files that
        # are already on disk under the (eventual) capture_dir by the time a reader could see the
        # manifest.
        provenance_entries = _capture_provenance_files(arc_project_directory=arc_project_directory,
                                                       capture_dir=staging_dir)

        _write_manifest(capture_dir=staging_dir, arc_project_directory=arc_project_directory,
                        manifest_entries=manifest_entries, provenance_entries=provenance_entries,
                        energy_settings=energy_settings, networks=manifest_networks,
                        sensitivity_aggregation=sensitivity_aggregation)

        # Self-check the fully staged capture with the same verifier a downstream consumer would
        # use, BEFORE it is ever swapped into capture_dir's place. This is the natural post-write
        # check verify_capture already exists to perform; running it here means a staged capture
        # that would fail its own verification is refused (and staging_dir discarded, capture_dir
        # untouched) rather than ever becoming the thing callers see as "the capture".
        verify_capture(staging_dir)

        # Every file has been staged and self-verified; nothing below can fail on CONTENT grounds
        # any more. Swap staging_dir into capture_dir's place atomically: rename any existing
        # capture_dir aside first (a same-filesystem, near-instant rename, not a copy -- os.replace
        # cannot atomically replace a non-empty directory directly), then os.replace the staged
        # directory into the now-vacant capture_dir path (atomic, since the target is absent), and
        # only remove the renamed-aside old tree once the new one is confirmed in place. This
        # closes the gap a naive rmtree-then-replace would leave open: a crash between the two
        # steps would destroy the prior good capture and leave nothing in its place, even though
        # os.replace itself is atomic. If os.replace itself fails, the old tree is renamed straight
        # back so capture_dir is left exactly as it was, never observed empty or half-populated.
        old_capture_dir = None
        if os.path.isdir(capture_dir):
            old_capture_dir = tempfile.mkdtemp(prefix=_CAPTURE_OLD_DIR_PREFIX, dir=parent_dir)
            os.rmdir(old_capture_dir)
            os.rename(capture_dir, old_capture_dir)
        try:
            os.replace(staging_dir, capture_dir)
        except BaseException:
            if old_capture_dir is not None:
                os.rename(old_capture_dir, capture_dir)
            raise
        if old_capture_dir is not None:
            shutil.rmtree(old_capture_dir)
    finally:
        if os.path.isdir(staging_dir):
            shutil.rmtree(staging_dir, ignore_errors=True)

    # Fill in each captured record's final, capture_dir-rooted artifact_path now that capture_dir
    # holds the (formerly staged) files -- deferred from the loop above because staging_dir no
    # longer exists by this point (os.replace renamed it into capture_dir), so an absolute path
    # built against it earlier would have dangled.
    for record, captured_artifact_relpath in pending_records:
        if captured_artifact_relpath is None:
            captured_records.append(record)
            continue
        record_dict = record.as_dict()
        # CRITICAL: the raw ARC pointer path must not cross this boundary -- only the vendored,
        # capture_dir-relative path is handed back to callers downstream of capture.
        record_dict['artifact_path'] = os.path.join(capture_dir, captured_artifact_relpath)
        captured_records.append(TSArtifactRecord.from_dict(record_dict))

    return CaptureResult(capture_dir=capture_dir, manifest_path=os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME),
                         records=tuple(captured_records), energy_settings=energy_settings, networks=manifest_networks)


def _preflight_network_sources(networks: dict) -> None:
    """
    Verify every network's LIVE source file exists and still hashes to its recorded identity.

    Runs before anything is written to the capture directory. A missing file and a hash mismatch
    are refused separately, each with a message naming the network and what to do about it: both
    mean the RMG ``pdep/`` tree this capture was told about is not the one the SA/selection pass
    examined, and vendoring whatever is (or is not) there now would silently freeze the wrong
    network.

    Args:
        networks (dict): The sidecar's ``network_id -> {source_path, source_sha256, method}``
                        block, already validated by ``validate_networks_block``.

    Raises:
        ValueError: If a network's ``source_path`` does not exist on disk, or the file's current
                   ``t3.pdep.cache.hash_file`` result differs from the recorded ``source_sha256``.
    """
    for network_id, entry in networks.items():
        source_path = entry['source_path']
        if not os.path.isfile(source_path):
            raise ValueError(
                f"Cannot capture the network source for '{network_id}': the RMG network file "
                f"'{source_path}' recorded in the join sidecar no longer exists on disk. The file the "
                f"SA/selection pass examined is gone before it could be vendored -- capture must run "
                f"before the RMG pdep/ directory of this iteration is cleaned up or overwritten."
            )
        live_sha256 = hash_file(source_path)
        if live_sha256 != entry['source_sha256']:
            raise ValueError(
                f"Refusing to capture the network source for '{network_id}': the RMG network file "
                f"'{source_path}' now hashes to {live_sha256}, but the join sidecar recorded "
                f"{entry['source_sha256']} when the SA/selection pass examined it. The file changed "
                f"between selection and capture, so vendoring it would freeze a network the selection "
                f"never actually saw. This is a real defect in the run's sequencing, not a warning."
            )


def _vendor_network_sources(networks: dict, capture_dir: str) -> dict:
    """
    Copy every (already pre-flight-verified) network source file into ``capture_dir``'s
    ``networks/`` subdirectory, verifying each copy by hash, and build the manifest's
    AUTHORITATIVE ``'networks'`` block.

    The subdirectory is rebuilt from scratch on every capture: a stale vendored copy left by a
    prior capture into this same (owned) ``capture_dir`` -- including one for a network this
    capture no longer contains -- would otherwise survive in silent disagreement with the new
    manifest. Each file is staged via ``tempfile.mkstemp`` + ``os.replace`` (the same atomic-copy
    idiom as ``_capture_provenance_files``), never copied directly onto its final name, and then
    re-hashed with ``t3.pdep.cache.hash_file``: since network files are vendored byte-for-byte
    (unlike the rewritten QM pointer files), the vendored hash must equal the sidecar's recorded
    ``source_sha256`` exactly, and any disagreement is a torn/corrupt copy refused on the spot.

    Args:
        networks (dict): The sidecar's ``network_id -> {source_path, source_sha256, method}``
                        block, already validated and pre-flight-verified against the live files.
        capture_dir (str): The capture directory being populated.

    Raises:
        ValueError: If a ``network_id`` would escape the ``networks/`` subdirectory as a file name,
                   or a vendored copy's hash does not match the recorded ``source_sha256``.

    Returns:
        dict: ``network_id -> {source_path, source_sha256, captured_path, method}``, with
             ``captured_path`` relative to ``capture_dir``.
    """
    networks_dir = os.path.join(capture_dir, _CAPTURE_NETWORKS_SUBDIR)
    if os.path.isdir(networks_dir):
        shutil.rmtree(networks_dir)
    os.makedirs(networks_dir)
    manifest_networks = dict()
    for network_id, entry in networks.items():
        captured_relpath = os.path.join(_CAPTURE_NETWORKS_SUBDIR, f'{network_id}.py')
        captured_abs_path = os.path.join(capture_dir, captured_relpath)
        # Confinement, mirroring the realpath idiom used across this module family: a network_id
        # containing a path separator or '..' (possible for a never-queued network, whose id was
        # never label-validated) must not be able to place its vendored copy outside networks/.
        resolved_networks_dir = os.path.realpath(networks_dir)
        resolved_captured = os.path.realpath(captured_abs_path)
        if os.path.dirname(resolved_captured) != resolved_networks_dir:
            raise ValueError(
                f"Refusing to vendor the network source for network_id={network_id!r}: its vendored "
                f"file name would resolve to '{resolved_captured}', outside the capture's networks/ "
                f"subdirectory '{resolved_networks_dir}'. A network id must be a plain file name "
                f"component; refusing rather than writing outside the capture directory."
            )
        fd, staged_path = tempfile.mkstemp(prefix='.network-', dir=networks_dir)
        os.close(fd)
        try:
            shutil.copyfile(entry['source_path'], staged_path)
            os.replace(staged_path, captured_abs_path)
        finally:
            if os.path.isfile(staged_path):
                os.remove(staged_path)
        captured_sha256 = hash_file(captured_abs_path)
        if captured_sha256 != entry['source_sha256']:
            raise ValueError(
                f"Torn or corrupt capture detected for network '{network_id}': the vendored network "
                f"file '{captured_abs_path}' hashes to {captured_sha256}, but its source "
                f"'{entry['source_path']}' was recorded (and pre-flight-verified) as "
                f"{entry['source_sha256']}. Refusing to return a capture whose vendored bytes disagree "
                f"with their own source."
            )
        manifest_networks[network_id] = {
            'source_path': entry['source_path'],
            'source_sha256': entry['source_sha256'],
            'captured_path': captured_relpath,
            'method': entry['method'],
        }
    return manifest_networks


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
        'coefficient': record.coefficient,
        'delta_ln_k': record.delta_ln_k,
        'path_reaction_labels': list(record.path_reaction_labels),
        'captured_artifact_path': captured_artifact_path,
        'captured_log_paths': captured_log_paths or dict(),
        'captured_artifact_sha256': captured_artifact_sha256,
        'captured_log_sha256': captured_log_sha256 or dict(),
    }


def _write_manifest(capture_dir: str, arc_project_directory: str, manifest_entries: list,
                    provenance_entries: list | None = None, energy_settings: dict | None = None,
                    networks: dict | None = None,
                    sensitivity_aggregation: str | None = None) -> str:
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
        energy_settings (dict, optional): The frozen energy-reference settings block (see
                                          ``t3.pdep.energy_settings.read_arc_energy_settings``), or
                                          ``None`` if this capture found zero usable artifacts.
                                          Recorded under the manifest's ``'energy_settings'`` key.
                                          Unlike ``'provenance'``, this key IS authoritative: it is
                                          the sole source downstream hybrid-network writing may rely
                                          on for a captured transition state's model chemistry and
                                          atom/bond corrections.
        networks (dict, optional): The vendored-network identity block (see
                                  ``_vendor_network_sources``): ``network_id -> {source_path,
                                  source_sha256, captured_path, method}``. Recorded under the
                                  manifest's ``'networks'`` key. Like ``'energy_settings'`` (and
                                  unlike ``'provenance'``), this key IS authoritative: downstream
                                  consumers must read the capture's own vendored network copy at
                                  ``captured_path`` and its recorded ``method``, never re-derive
                                  either from the ARC project or the RMG ``pdep/`` directory --
                                  both of which are expected to be gone by consumption time.

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
        # AUTHORITATIVE, unlike 'provenance' above: the sole source downstream hybrid-network
        # writing may rely on for this capture's model chemistry and atom/bond corrections.
        'energy_settings': energy_settings,
        # AUTHORITATIVE, like 'energy_settings' above: the sole source downstream consumers may
        # rely on for each network's vendored source file and frozen ME method.
        'networks': networks or dict(),
    }
    if sensitivity_aggregation is not None:
        # How the per-TS coefficient/delta_ln_k evidence in 'transition_states' was aggregated,
        # when it was NOT T3's in-run convention (one observable-selected, threshold-gated
        # direction key -- the unmarked default). A manifest carrying
        # SENSITIVITY_AGGREGATION_ALL_DIRECTIONS_MAX_ABS holds ungated max-|coefficient| values
        # over ALL direction keys (the standalone PES loop's convention); the two conventions'
        # admissible ranges are incompatible and must never be compared on the shared field name.
        content['sensitivity_aggregation'] = sensitivity_aggregation

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
    Archive ARC's ``output/status.yml``, ``output/output.yml``, and the ``calcs/statmech/kinetics``/
    ``calcs/statmech/thermo`` ``input.py`` files as INERT provenance.

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

    Any of these files may legitimately not exist (e.g. a transition state resolved from ARC
    without ever going through the ``output/`` finalization step, or a project that never ran a
    ``thermo`` statmech pass at all); absence is recorded, not raised. The kinetics/thermo
    ``input.py`` files are archived here as inert provenance ONLY -- the AUTHORITATIVE energy
    settings frozen from ``calcs/statmech/kinetics/input.py`` live in the manifest's separate,
    authoritative ``'energy_settings'`` key (see ``read_arc_energy_settings`` and
    ``capture_ts_artifacts``); nothing may read these archived copies back as authority for
    anything, same as ``status.yml``/``output.yml`` above.

    Args:
        arc_project_directory (str): The ARC project directory to read ``output/status.yml``,
                                     ``output/output.yml``, ``calcs/statmech/kinetics/input.py``, and
                                     ``calcs/statmech/thermo/input.py`` from.
        capture_dir (str): The capture directory being populated.

    Returns:
        list: One dict per file (``status.yml``, ``output.yml``, ``statmech_kinetics_input.py``,
             then ``statmech_thermo_input.py``), each with keys ``source_path``, ``captured_path``
             (relative to ``capture_dir``, or ``None`` if the source file did not exist), ``sha256``
             (or ``None``), and ``source_mtime`` (the source file's ``os.path.getmtime`` result at
             archive time, or ``None``).
    """
    provenance_dir = os.path.join(capture_dir, _CAPTURE_PROVENANCE_SUBDIR)
    # (source_path, captured_file_name): captured_file_name is distinct for the two input.py files
    # (kinetics/thermo) so they never clobber each other under provenance/, since both are literally
    # named 'input.py' at their respective sources.
    files_to_archive = (
        (os.path.join(arc_project_directory, 'output', 'status.yml'), 'status.yml'),
        (os.path.join(arc_project_directory, 'output', 'output.yml'), 'output.yml'),
        (os.path.join(arc_project_directory, 'calcs', 'statmech', 'kinetics', 'input.py'),
         'statmech_kinetics_input.py'),
        (os.path.join(arc_project_directory, 'calcs', 'statmech', 'thermo', 'input.py'),
         'statmech_thermo_input.py'),
    )
    entries = list()
    for source_path, captured_file_name in files_to_archive:
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
        captured_relpath = os.path.join(_CAPTURE_PROVENANCE_SUBDIR, captured_file_name)
        captured_abs_path = os.path.join(capture_dir, captured_relpath)
        fd, staged_path = tempfile.mkstemp(prefix=f'.provenance-{captured_file_name}-',
                                           dir=provenance_dir)
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


def _confine_to_capture_dir(capture_dir: str, path: str) -> None:
    """
    Refuse a manifest-supplied path that resolves outside ``capture_dir``.

    Mirrors ``t3.pdep.discovery._confine_to_project``'s exact standalone-function shape and
    escape-polarity check (``os.path.realpath`` of both sides, then ``os.path.commonpath``): a
    capture is supposed to be self-contained, so nothing a manifest names (``captured_artifact_path``,
    a ``captured_log_sha256`` key, or the ``'networks'`` block's ``captured_path``) may ever be
    allowed to resolve outside the directory that manifest lives in. Without this check,
    ``os.path.join`` silently lets an absolute value in the manifest override ``capture_dir``
    entirely (``os.path.join(a, '/abs/path') == '/abs/path'``), and a relative value containing
    ``'..'`` components can walk back out of it -- either way, a tampered or corrupted manifest
    could make verification read (and "pass") an external file that happens to match, rather than
    the capture's own vendored copy.

    Args:
        capture_dir (str): The capture directory a manifest-supplied path must stay under.
        path (str): The candidate path (already joined onto ``capture_dir``) to check.

    Raises:
        ValueError: If ``path`` resolves to a location outside ``capture_dir``.
    """
    resolved_capture_dir = os.path.realpath(capture_dir)
    resolved_path = os.path.realpath(path)
    if resolved_path == resolved_capture_dir \
            or os.path.commonpath([resolved_capture_dir, resolved_path]) != resolved_capture_dir:
        raise ValueError(f"Refusing to verify '{path}': it resolves to '{resolved_path}', which is outside the "
                         f"capture directory '{capture_dir}' (resolved to '{resolved_capture_dir}'). A capture "
                         f"must be self-contained; a manifest-supplied path may never escape it.")


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


def _capture_lock_path(capture_dir: str) -> str:
    """
    The interprocess lock file path for a given ``capture_dir``.

    Always a sibling of ``capture_dir`` (never inside it), so the lock file itself can never be
    mistaken for capture content by ``_refuse_unowned_capture_dir``, ``_has_valid_capture_manifest``
    or ``verify_capture``, and so it survives across the rename-aside/``os.replace`` swap (which
    only ever touches ``capture_dir`` and its ``.capture-old-*``/``.capture-staging-*`` siblings).

    Args:
        capture_dir (str): The capture directory this lock guards.

    Returns:
        str: The lock file's path.
    """
    return os.path.abspath(capture_dir).rstrip(os.sep) + '.capture-lock'


def _is_stale_capture_lock(lock_path: str) -> bool:
    """
    Whether ``lock_path`` was written by a process that is now confirmed dead.

    Best-effort, and deliberately biased toward "not stale" whenever it cannot be sure: unreadable
    or unparsable lock content, or a PID belonging to a live process (whether or not this process
    can signal it), are all treated as "still held". Only ``os.kill(pid, 0)`` raising
    ``ProcessLookupError`` -- meaning the PID no longer identifies any process at all -- counts as
    proof of death. A ``PermissionError`` means the PID exists but is owned by a different user, so
    it is NOT stale from this function's point of view; treating "cannot confirm" as "confirmed
    dead" would turn an ordinary permissions quirk into a false stale-lock removal, which is exactly
    the failure mode this function must avoid given it directly precedes deleting the lock file.

    Args:
        lock_path (str): The lock file to inspect.

    Returns:
        bool: True only if the lock file's recorded PID is confirmed to no longer exist.
    """
    try:
        with open(lock_path) as f:
            content = f.read().strip()
        pid = int(content.splitlines()[0]) if content else None
    except (OSError, ValueError, IndexError):
        return False  # unreadable/unparsable: cannot confirm anything, so do not treat as stale
    if pid is None:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return True
    except PermissionError:
        return False
    else:
        return False


def _read_capture_lock_holder(lock_path: str) -> str:
    """Best-effort, human-readable description of a capture lock's holder, for error messages only."""
    try:
        with open(lock_path) as f:
            content = f.read().strip()
        return f'pid {content}' if content else '(lock file is empty)'
    except OSError as e:
        return f'(could not read lock file: {e})'


def _acquire_capture_lock(capture_dir: str) -> str:
    """
    Take an exclusive, interprocess lock on ``capture_dir`` before any capture work begins.

    Without this, two concurrent T3 processes capturing the same ``capture_dir`` could interleave
    their staging and swap (e.g. both build a staging directory, then race each other's
    rename-aside/``os.replace`` swap), producing a capture neither of them individually verified.

    Implemented with ``os.open(..., os.O_CREAT | os.O_EXCL)``: POSIX guarantees the
    create-if-absent check and the file's creation happen as a single atomic operation, so two
    processes racing to create the same lock path can never both "win" -- exactly one ``os.open``
    call succeeds, no matter how closely timed, and every other one raises ``FileExistsError``. The
    winning process's PID is written into the file afterward purely for diagnosability (so a human
    staring at a stuck lock can identify who holds it); the lock's exclusivity comes entirely from
    the atomic file creation, never from the PID written afterward.

    Stale-lock handling: if a lock file already exists but ``_is_stale_capture_lock`` confirms its
    recorded PID no longer corresponds to any live process (e.g. that process was killed without
    reaching its ``finally``), the stale lock is removed and acquisition is retried exactly once.
    This is inherently best-effort (there is an unavoidable race between the liveness check and the
    removal, and PIDs can in principle be recycled by the OS in that window), but the retried
    ``os.open`` is itself the atomic check, so a genuine concurrent live holder is still caught even
    if the stale-cleanup race is lost.

    Args:
        capture_dir (str): The capture directory about to be captured into.

    Returns:
        str: The lock file path, to be passed to ``_release_capture_lock`` once the caller is done
            (regardless of success or failure).

    Raises:
        RuntimeError: If the lock is already held by a live (or unidentifiable) process, or a
                     retried acquisition after stale-lock cleanup still fails. Fails CLOSED: this
                     function never returns normally without truly holding the lock.
    """
    lock_path = _capture_lock_path(capture_dir)
    parent_dir = os.path.dirname(lock_path)
    if parent_dir and not os.path.isdir(parent_dir):
        os.makedirs(parent_dir)

    for attempt in range(2):
        try:
            # 0o600, not 0o644: this file is written by _acquire_capture_lock and read only by
            # _read_capture_lock_holder, both inside one user's own capture directory, to recover
            # the PID of a possibly-dead holder. There is no cross-user reader to serve, so the
            # narrower mode is simply the right one (CodeQL alert 112, py/overly-permissive-file).
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        except FileExistsError:
            if attempt == 0 and _is_stale_capture_lock(lock_path):
                try:
                    os.remove(lock_path)
                except FileNotFoundError:
                    pass  # another process already reclaimed it; the retry below is still safe
                continue
            holder = _read_capture_lock_holder(lock_path)
            raise RuntimeError(
                f"Refusing to capture into '{capture_dir}': its lock file '{lock_path}' is already "
                f"held ({holder}). Another process appears to be capturing this same directory "
                f"concurrently. Failing closed rather than racing an in-progress capture; if the "
                f"holding process is confirmed dead, remove the lock file manually and retry."
            )
        else:
            try:
                os.write(fd, f'{os.getpid()}\n'.encode('utf-8'))
            finally:
                os.close(fd)
            return lock_path
    # Unreachable: every iteration above either returns or raises. Kept only so this function can
    # never silently fall through to "no lock, no error" if that analysis is ever wrong.
    raise RuntimeError(f"Refusing to capture into '{capture_dir}': could not acquire its lock "
                      f"'{lock_path}' after stale-lock recovery.")


def _release_capture_lock(lock_path: str) -> None:
    """
    Release a capture lock previously taken by ``_acquire_capture_lock``.

    Always called from a ``finally`` block, so the lock is released whether the capture it guarded
    succeeded or raised. Tolerates the file already being gone rather than raising during cleanup,
    since raising here would mask whatever real exception is already in flight in the caller.

    Args:
        lock_path (str): The lock file path returned by ``_acquire_capture_lock``.
    """
    try:
        os.remove(lock_path)
    except FileNotFoundError:
        # Already gone -- a stale-lock reclaim by another process removed it, or this is a second
        # release on the same path. The post-condition this function promises ("the lock is not
        # held") is satisfied either way, and raising here would mask whatever real exception is
        # already unwinding through the caller's finally.
        pass


def _recover_capture_swap_state(capture_dir: str) -> None:
    """
    Self-heal ``capture_dir`` from a crash inside a PRIOR call's atomic swap, and clear leftover
    staging scratch. Called at the start of every ``capture_ts_artifacts`` call, while this process
    already holds ``capture_dir``'s lock (see ``_acquire_capture_lock``), so nothing else can be
    concurrently mid-build or mid-swap while this runs.

    Why this exists: the swap that promotes a fully-staged, self-verified new capture into
    ``capture_dir``'s place (see ``_capture_ts_artifacts_locked``) works by renaming any existing
    ``capture_dir`` aside to a sibling ``.capture-old-*`` path, then ``os.replace``-ing the staging
    directory into the now-vacant ``capture_dir`` path. ``os.replace`` itself is atomic, but POSIX
    has no primitive that atomically replaces a non-empty directory's contents in one step, so the
    rename-aside and the ``os.replace`` are necessarily two separate syscalls. A process that dies
    (``kill -9``, OOM-kill, power loss -- anything that skips the in-process ``except``/``finally``
    handlers already in ``_capture_ts_artifacts_locked``) BETWEEN those two syscalls leaves
    ``capture_dir`` completely ABSENT, with the last known-good capture sitting under the
    ``.capture-old-*`` name, which nothing else recognizes or recovers on its own.

    What this function actually guarantees, stated honestly: ``capture_dir`` is never left
    permanently missing after such a crash, because the very next call that takes the lock finds
    and restores it here, before doing anything else. What it does NOT provide: a single atomic
    syscall that makes the swap itself crash-proof in real time -- some other process reading
    ``capture_dir`` directly, bypassing this module's lock, between the crash and this recovery
    would transiently see it absent (never corrupt, just momentarily missing) until this function
    next runs. True real-time atomicity would require ``capture_dir`` to be something other than an
    ordinary directory (e.g. a symlink swapped with ``os.replace``) -- a larger change to this
    module's on-disk contract, which every consumer of ``capture_dir`` would need to tolerate, that
    this fix does not take on.

    Args:
        capture_dir (str): The capture directory about to be populated.

    Raises:
        ValueError: If ``capture_dir`` itself is missing or empty AND more than one sibling
                   ``.capture-old-*`` directory holds a valid capture manifest -- an ambiguous state
                   that should never arise from a single crashed swap. Refuses to guess which one is
                   the real prior capture rather than silently picking one.
    """
    abs_capture_dir = os.path.abspath(capture_dir)
    parent_dir = os.path.dirname(abs_capture_dir)
    if not os.path.isdir(parent_dir):
        return  # capture_dir's parent does not even exist yet; nothing to recover

    def _sibling_dirs(prefix):
        return sorted(
            os.path.join(parent_dir, name) for name in os.listdir(parent_dir)
            if name.startswith(prefix) and os.path.isdir(os.path.join(parent_dir, name))
        )

    if _has_valid_capture_manifest(capture_dir):
        # The swap already completed (or never started). Any '.capture-old-*' siblings are
        # leftover garbage from a crash AFTER os.replace succeeded but before the old tree's
        # cleanup rmtree ran; any '.capture-staging-*' siblings are leftover garbage from a crash
        # before some capture's swap ever ran. Both are safe to discard now: this process holds
        # capture_dir's lock, so nothing else can be concurrently mid-build against it.
        for stale_dir in _sibling_dirs(_CAPTURE_OLD_DIR_PREFIX) + _sibling_dirs(_CAPTURE_STAGING_DIR_PREFIX):
            shutil.rmtree(stale_dir, ignore_errors=True)
        return

    # capture_dir itself is missing or invalid. Only ever treat it as recoverable-from when it is
    # actually absent or empty -- if it exists and holds unrelated, non-empty content, leave it
    # alone entirely and let _refuse_unowned_capture_dir raise its own clear error, rather than ever
    # deleting content this module cannot prove it owns.
    capture_dir_absent_or_empty = not os.path.isdir(capture_dir) or not os.listdir(capture_dir)
    if not capture_dir_absent_or_empty:
        return

    old_dirs = _sibling_dirs(_CAPTURE_OLD_DIR_PREFIX)
    valid_old_dirs = [d for d in old_dirs if _has_valid_capture_manifest(d)]
    if not valid_old_dirs:
        return  # nothing recoverable; leave any (manifest-less) old/staging siblings alone

    if len(valid_old_dirs) > 1:
        raise ValueError(
            f"Refusing to recover '{capture_dir}': found {len(valid_old_dirs)} sibling "
            f"'{_CAPTURE_OLD_DIR_PREFIX}*' directories that each hold a valid capture manifest "
            f"({valid_old_dirs!r}), while capture_dir itself is missing or empty. This should never "
            f"happen from a single crashed swap; refusing to guess which one is the real prior "
            f"capture rather than silently picking one -- resolve this manually."
        )

    recovered_dir = valid_old_dirs[0]
    if os.path.isdir(capture_dir):
        os.rmdir(capture_dir)  # confirmed empty above
    os.rename(recovered_dir, capture_dir)

    # Now that capture_dir holds the recovered, valid capture, clear every other leftover
    # '.capture-old-*'/'.capture-staging-*' sibling the same way the already-valid branch above
    # does (recovered_dir itself no longer appears in the listing, having just been renamed away).
    for stale_dir in _sibling_dirs(_CAPTURE_OLD_DIR_PREFIX) + _sibling_dirs(_CAPTURE_STAGING_DIR_PREFIX):
        shutil.rmtree(stale_dir, ignore_errors=True)


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
                   disk; a captured file's actual sha256 does not match the hash recorded for it
                   in the manifest; or the ``'networks'`` block is absent, does not exactly cover
                   the networks the transition state entries reference (in either direction), or
                   names a vendored network file that is missing or no longer hashes to its
                   recorded ``source_sha256``. The networks check is structurally non-vacuous:
                   ``transition_states`` is non-empty (enforced above) and every entry names a
                   ``network_id``, so the referenced-network set an empty ``'networks'`` block is
                   compared against is never empty -- an absent or empty block always fails here,
                   it can never pass because there was nothing to iterate over.
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
    referenced_network_ids = set()
    ts_records = list()
    for entry in transition_states:
        record_count += 1
        network_id = entry.get('network_id')
        if not network_id:
            raise ValueError(
                f"Refusing to verify '{capture_dir}': a transition state entry in the manifest at "
                f"'{manifest_path}' has no network_id ({entry!r}). Every entry capture_ts_artifacts "
                f"ever wrote names the network its transition state belongs to; a missing one means "
                f"the manifest is malformed or was hand-edited."
            )
        referenced_network_ids.add(network_id)
        ts_label = entry.get('arc_ts_label')
        status = entry.get('status')
        # Fail closed on the status itself before trusting any path derived from it: a manifest is
        # untrusted data re-read off disk, not something this function produced, so an unrecognized
        # status (hand-edited, corrupted, or from a newer/older capture format) must be refused
        # rather than silently falling through whichever branch below happens to match its path.
        if status not in TS_ARTIFACT_STATUSES:
            raise ValueError(
                f"Refusing to verify '{capture_dir}': the manifest entry for transition state "
                f"{ts_label!r} at '{manifest_path}' has an unrecognized status {status!r}; expected "
                f"one of {sorted(TS_ARTIFACT_STATUSES)}. A missing or unrecognized status means the "
                f"manifest is malformed, hand-edited, or predates this capture format."
            )
        captured_artifact_path = entry.get('captured_artifact_path')
        # Enforce the status/path pairing discovery itself would have produced (see the two
        # frozensets' definitions in t3.pdep.discovery): counting captured_artifact_count purely off
        # "is captured_artifact_path non-None" -- without ever looking at status -- let a manifest
        # hand-edited to 'status: usable' with a null captured_artifact_path verify cleanly and then
        # silently vanish from every downstream count, even though nothing was ever actually
        # captured for it.
        if status in TS_ARTIFACT_STATUSES_REQUIRING_ARTIFACT_PATH and captured_artifact_path is None:
            raise ValueError(
                f"Refusing to verify '{capture_dir}': the manifest entry for transition state "
                f"{ts_label!r} at '{manifest_path}' has status {status!r}, which must always carry a "
                f"captured_artifact_path, but captured_artifact_path is None. This status/path "
                f"combination is never produced by capture_ts_artifacts; the manifest is malformed or "
                f"was hand-edited."
            )
        if status in TS_ARTIFACT_STATUSES_REQUIRING_NO_ARTIFACT_PATH and captured_artifact_path is not None:
            raise ValueError(
                f"Refusing to verify '{capture_dir}': the manifest entry for transition state "
                f"{ts_label!r} at '{manifest_path}' has status {status!r}, which must never carry a "
                f"captured_artifact_path, but captured_artifact_path is {captured_artifact_path!r}. This "
                f"status/path combination is never produced by capture_ts_artifacts; the manifest is "
                f"malformed or was hand-edited."
            )
        resolved_artifact_path = None
        if captured_artifact_path is not None:
            captured_artifact_count += 1
            resolved_artifact_path = os.path.join(capture_dir, captured_artifact_path)
            _confine_to_capture_dir(capture_dir=capture_dir, path=resolved_artifact_path)
            _verify_one_captured_file(
                capture_dir=capture_dir,
                relpath=captured_artifact_path,
                expected_sha256=entry.get('captured_artifact_sha256'),
                ts_label=ts_label,
            )
            # A captured artifact must carry the sensitivity evidence that justified selecting this
            # transition state in the first place: without it, a partial-capture / "is the dominant
            # transition state missing?" decision cannot be made from the capture alone, including on
            # a restart where the in-memory selection that originally produced this evidence is gone.
            coefficient, delta_ln_k = entry.get('coefficient'), entry.get('delta_ln_k')
            if not _is_finite_manifest_number(coefficient) or not _is_finite_manifest_number(delta_ln_k):
                raise ValueError(
                    f"Refusing to verify '{capture_dir}': the manifest entry for transition state "
                    f"{ts_label!r} at '{manifest_path}' has a captured artifact but no valid "
                    f"'coefficient'/'delta_ln_k' sensitivity evidence (got coefficient={coefficient!r}, "
                    f"delta_ln_k={delta_ln_k!r}). Every captured transition state must have been queued "
                    f"with the sensitivity evidence that justified its selection frozen onto its join "
                    f"record -- a missing or malformed value here means the manifest is malformed or "
                    f"predates this capture format."
                )
        for log_relpath, expected_log_sha256 in (entry.get('captured_log_sha256') or dict()).items():
            _verify_one_captured_file(
                capture_dir=capture_dir,
                relpath=log_relpath,
                expected_sha256=expected_log_sha256,
                ts_label=ts_label,
            )
        ts_records.append(TSArtifactRecord(
            network_id=network_id,
            network_ts_label=entry.get('network_ts_label'),
            arc_ts_label=ts_label,
            status=status,
            artifact_path=resolved_artifact_path,
            converged=entry.get('converged'),
            reason=entry.get('reason') or '',
            coefficient=entry.get('coefficient'),
            delta_ln_k=entry.get('delta_ln_k'),
            path_reaction_labels=tuple(entry.get('path_reaction_labels') or ()),
        ))

    if captured_artifact_count > 0:
        energy_settings = content.get('energy_settings')
        if not isinstance(energy_settings, dict) or not energy_settings.get('model_chemistry'):
            raise ValueError(
                f"Refusing to verify '{capture_dir}': the manifest at '{manifest_path}' captured "
                f"{captured_artifact_count} artifact(s) but has no valid 'energy_settings' block "
                f"with a non-blank 'model_chemistry'. A capture with at least one captured artifact "
                f"must always have frozen its authoritative energy settings at capture time (see "
                f"capture_ts_artifacts/read_arc_energy_settings) -- a missing or incomplete block "
                f"here means the manifest is malformed or was hand-edited."
            )
        # Presence is not enough: the block came back off disk as YAML, so every value's TYPE has
        # to be re-checked against the one authoritative definition of the frozen block's shape.
        validate_frozen_energy_settings(energy_settings, context=f"the capture manifest at '{manifest_path}'")

    # Verify the AUTHORITATIVE networks block. Structurally non-vacuous by construction: the
    # transition_states loop above enforced record_count >= 1 and a network_id on every entry, so
    # referenced_network_ids is never empty here -- an absent or empty 'networks' block therefore
    # always fails the coverage check inside validate_networks_block; it can never slip through
    # as "nothing to iterate over".
    if 'networks' not in content:
        raise ValueError(
            f"Refusing to verify '{capture_dir}': the manifest at '{manifest_path}' has no 'networks' "
            f"key. Every capture freezes, per network its transition states belong to, the vendored "
            f"network source file, its hash, and the ME method (the authoritative identity downstream "
            f"consumers must read instead of the long-gone RMG pdep/ directory) -- a manifest without "
            f"it is malformed or predates this capture format, and there is no safe fallback."
        )
    if len(referenced_network_ids) == 0:
        raise ValueError(
            f"Refusing to verify '{capture_dir}': no transition state entry in the manifest at "
            f"'{manifest_path}' references a network, so the networks block has nothing to be checked "
            f"against and this verification would be vacuous."
        )
    networks = content['networks'] or dict()
    validate_networks_block(referenced_network_ids=referenced_network_ids,
                            networks=networks,
                            context=f"the capture manifest at '{manifest_path}'")
    for network_id, entry in networks.items():
        captured_path = entry.get('captured_path')
        if not isinstance(captured_path, str) or not captured_path.strip():
            raise ValueError(
                f"Refusing to verify '{capture_dir}': the networks entry for '{network_id}' in the "
                f"manifest at '{manifest_path}' has an invalid captured_path ({captured_path!r}); "
                f"expected the vendored network file's path relative to the capture directory."
            )
        captured_abs_path = os.path.join(capture_dir, captured_path)
        _confine_to_capture_dir(capture_dir=capture_dir, path=captured_abs_path)
        if not os.path.isfile(captured_abs_path):
            raise ValueError(
                f"Torn capture detected in '{capture_dir}': the manifest records a vendored network "
                f"source for '{network_id}' at '{captured_path}', but it no longer exists on disk. "
                f"Refusing to use a capture whose vendored networks disagree with its own manifest."
            )
        actual_sha256 = hash_file(captured_abs_path)
        if actual_sha256 != entry['source_sha256']:
            raise ValueError(
                f"Torn or tampered capture detected in '{capture_dir}': the vendored network source "
                f"for '{network_id}' at '{captured_path}' hashes to {actual_sha256}, but the manifest "
                f"records {entry['source_sha256']}. Refusing to use it."
            )

    sensitivity_aggregation = content.get('sensitivity_aggregation')
    if sensitivity_aggregation is not None and sensitivity_aggregation not in SENSITIVITY_AGGREGATIONS:
        raise ValueError(
            f"Refusing to verify '{capture_dir}': the manifest at '{manifest_path}' carries an "
            f"unrecognized 'sensitivity_aggregation' value {sensitivity_aggregation!r}; expected "
            f"the key to be absent (T3's unmarked in-run convention) or one of "
            f"{sorted(SENSITIVITY_AGGREGATIONS)}. An unrecognized marker means the manifest is "
            f"hand-edited or from a newer capture format."
        )
    return VerifyResult(capture_dir=capture_dir, manifest_path=manifest_path, record_count=record_count,
                        captured_artifact_count=captured_artifact_count,
                        networks=networks, energy_settings=content.get('energy_settings'),
                        ts_records=tuple(ts_records),
                        sensitivity_aggregation=sensitivity_aggregation)


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
    _confine_to_capture_dir(capture_dir=capture_dir, path=abs_path)
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


def _is_finite_manifest_number(value) -> bool:
    """
    Whether a value read back from a manifest is a usable finite real number.

    Mirrors ``t3.pdep.join._is_finite_number``/``t3.pdep.selector._is_finite``: ``bool`` is
    deliberately excluded even though it is an ``int`` subtype, and NaN/inf are refused since YAML
    can round-trip them but they carry no usable sensitivity information.

    Args:
        value: The value to check, as read back from the manifest (already YAML-deserialized).

    Returns:
        bool: ``True`` if ``value`` is a finite, non-boolean ``int``/``float``.
    """
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
