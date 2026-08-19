"""
t3 pdep discovery module

The post-ARC counterpart to ``t3.pdep.join``: given the join records written before ARC ran,
decide, per network transition state, whether ARC actually produced a usable quantum-chemistry
artifact -- and give an honest reason when it did not.

A queued transition state is NOT guaranteed to come back with an artifact. In one real ARC
project inspected while designing this module, only 2 of 5 queued transition states had an
Arkane statmech ``.py`` waiting for them: a queued-but-missing result is the common case here,
not the exception, and this module's every decision has to treat it that way rather than as a
surprise to special-case.

Two things are never trusted from the join step and are recomputed here instead:

* The ARC transition state label. ``t3.pdep.join.arc_ts_label`` is a pure function of
  ``(network_id, network_ts_label)``, so it is recomputed rather than read from the stored join
  record -- a hand-edited or stale sidecar entry must not silently redirect which artifact gets
  read back for a given network transition state.
* The expected artifact path, for the same reason, via ``t3.pdep.join.expected_ts_artifact_path``.
  Any disagreement between the stored path and the recomputed one is recorded in the reason text,
  but the recomputed path is always the one actually read.

ARC exposes two different sources for whether a transition state's job converged, and they
disagree in a way that matters:

* ``output/status.yml`` is the raw ``self.output`` mapping ARC keeps per label
  (``arc/main.py`` -> ``self.scheduler.output``), and it preserves the true tri-state:
  ``convergence`` is ``True`` (converged), ``False`` (explicitly failed), or ``None``/absent
  (never assessed -- e.g. the job never ran, or ARC's project was interrupted first).
* ``output/output.yml`` (``arc/output.py::write_output_yml``) is a later, consolidated summary,
  and it collapses that tri-state: it writes ``'converged': entry.get('convergence') is True``
  (``arc/output.py:436``), so ``False`` and ``None`` become indistinguishable. A transition state
  ARC never even attempted therefore reads identically, from ``output.yml`` alone, to one ARC
  tried and explicitly failed.

This module prefers ``output/status.yml`` when it exists, and falls back to ``output/output.yml``
only when ``status.yml`` is absent. This is the OPPOSITE of ARC's own "latest wins" convention,
and deliberately so: this module's every decision turns on distinguishing "explicitly failed"
from "not yet assessed", which is exactly the distinction ``output.yml`` cannot make. Preferring
the lossier source for the one question that depends on the lost bit would be backwards. When
only ``output.yml`` is available and it reports a transition state as not converged, that verdict
is still trusted fail-closed (treated as ``False``, not ``None``) per this module's standing
refuse-rather-than-degrade rule -- but the reason text says plainly that the source could not
distinguish "failed" from "not yet determined", so nothing downstream overstates what ARC
actually said. Every reason string also records which file the verdict came from.

One exception to the ``status.yml``-first preference, when BOTH sources exist: ``status.yml``
saying converged while ``output.yml`` (the later, consolidated summary) says NOT converged is
treated fail-closed as not converged, since a stale ``status.yml`` from an interrupted earlier
run could otherwise promote a transition state ARC's final output says failed -- UNLESS
``status.yml`` is strictly newer (by file mtime) than ``output.yml``. ARC writes ``status.yml``
unconditionally, before ``output.yml``, within a single ``execute()`` call, so ``output.yml`` is
never older than ``status.yml`` in any run that completes normally; a ``status.yml`` that is
nonetheless demonstrably newer can only be the product of a later, restarted run that got
further than whatever left ``output.yml`` behind, so its ``True`` is trusted instead. Only that
one direction is a conflict: ``status.yml: null`` vs ``output.yml: false`` is the expected lossy
collapse described above, and the reverse disagreement simply keeps ``status.yml``'s richer
verdict. See ``read_arc_convergence``.

Regardless of source, unknown convergence (``None``) is never treated as unconverged -- but nor
is it treated as converged. A transition state whose convergence status could not be determined
(no status source at all right after ``arc.execute()``, or the label simply absent from whatever
source does exist) is genuinely suspicious immediately after an ARC run, not neutral: silently
promoting it to ``usable`` would let a QM result T3 never actually confirmed as converged flow
straight into a hybrid kinetic model. So an artifact with unknown convergence that otherwise
parses cleanly is ``ARTIFACT_STATUS_UNVERIFIED``, a distinct outcome from both ``usable`` and
``unusable``: it is not silently discarded (the artifact and its path are still worth recording,
and a human -- or a re-run that produces a real status source -- can resolve it), but
``evaluate_pdep_hybrid`` never counts it toward ``usable_count``, so it can never on its own let a
hybrid network be accepted. Failing it hard as ``unusable`` would be just as wrong as accepting it
uncritically, since a not-yet-assessed job is not evidence of anything -- but it must not be
free-riding on ``usable`` either.

USABLE requires all of: convergence is explicitly ``True`` (not merely non-``False``), the
expected artifact file exists, it is actually shaped like an ARC statmech species-file
(``_check_artifact_shape`` -- merely parsing as Python is not enough), and it can actually be read
back --
``t3.pdep.hybrid._read_qm_artifact`` parses it and
resolves every ``Log(...)`` path it references. That last check is a deliberate preflight: an
artifact that exists on disk but whose referenced log files have gone missing (a very real
failure mode, since ARC wipes and rewrites ``calcs/statmech/kinetics`` on every restart) would
otherwise raise LATE, inside whatever ME-solving step tries to actually consume it, far from
where the fix -- re-running or re-vendoring that one transition state -- would be made.
"""

import ast
import os

from arc.common import read_yaml_file

from t3.pdep.hybrid import _parse_as_ast, _read_qm_artifact
from t3.pdep.join import (JOIN_STATUS_ALREADY_PRESENT,
                          JOIN_STATUS_NOT_QUEUED,
                          JOIN_STATUS_QUEUED,
                          arc_ts_label,
                          expected_ts_artifact_path,
                          )

# A transition state's expected artifact exists, is confined to the ARC project directory,
# parses, has the ARC statmech species-file shape (``REQUIRED_ARTIFACT_ASSIGNMENTS``), resolves
# every Log(...) it references, and is explicitly converged (``converged is True``).
ARTIFACT_STATUS_USABLE = 'usable'
# The transition state was queued (or already known) but no artifact was found for it.
ARTIFACT_STATUS_MISSING = 'missing'
# An artifact was found but is not safe to use: explicitly unconverged, unparseable, escapes the
# ARC project directory, or references a log file that no longer exists.
ARTIFACT_STATUS_UNUSABLE = 'unusable'
# The transition state could not even be queued to ARC; carried over unchanged from the join step.
ARTIFACT_STATUS_NOT_QUEUED = 'not_queued'
# An artifact was found, is confined to the ARC project directory, parses, has the ARC
# statmech species-file shape, and resolves every Log(...) it references -- but convergence is
# unknown (``converged is None``): no status source at all was found, or the label is absent from
# whatever status source does exist. This is deliberately distinct from both USABLE (never
# promoted from unknown convergence -- fail-closed) and UNUSABLE (an unassessed job is not
# evidence of failure either): the record is kept for visibility, but ``evaluate_pdep_hybrid``
# never counts it toward ``usable_count``.
ARTIFACT_STATUS_UNVERIFIED = 'unverified'

TS_ARTIFACT_STATUSES = frozenset({
    ARTIFACT_STATUS_USABLE, ARTIFACT_STATUS_MISSING, ARTIFACT_STATUS_UNUSABLE, ARTIFACT_STATUS_NOT_QUEUED,
    ARTIFACT_STATUS_UNVERIFIED,
})

# Which statuses this module ever pairs with a non-None artifact_path, and which it never does --
# read straight off this function's own branches below (the USABLE/UNVERIFIED tail sets
# artifact_path=recomputed_path; MISSING/UNUSABLE/NOT_QUEUED all set it None). A capture manifest
# is untrusted data re-read off disk by a consumer (``t3.pdep.capture.verify_capture``) that never
# ran this function itself, so it cannot assume a status/path pairing this module would never
# produce; these two sets let that consumer enforce the pairing is exactly what discovery would
# have written, rather than only checking the path when one happens to be present.
TS_ARTIFACT_STATUSES_REQUIRING_ARTIFACT_PATH = frozenset({ARTIFACT_STATUS_USABLE, ARTIFACT_STATUS_UNVERIFIED})
TS_ARTIFACT_STATUSES_REQUIRING_NO_ARTIFACT_PATH = frozenset({
    ARTIFACT_STATUS_MISSING, ARTIFACT_STATUS_UNUSABLE, ARTIFACT_STATUS_NOT_QUEUED,
})

# Every selected transition state came back usable: the hybrid network is fully QM-refined.
HYBRID_STATUS_COMPLETE = 'complete'
# At least one, but not all, selected transition states came back usable.
HYBRID_STATUS_PARTIAL_SELECTED_QM = 'partial_selected_qm'
# No selected transition state came back usable (including the case of no selected transition
# states at all): there is nothing here to build a hybrid network from.
HYBRID_STATUS_NOT_EVALUATED = 'not_evaluated'

HYBRID_STATUSES = frozenset({HYBRID_STATUS_COMPLETE, HYBRID_STATUS_PARTIAL_SELECTED_QM, HYBRID_STATUS_NOT_EVALUATED})

# The module-level assignments a usable artifact MUST carry to count as an ARC statmech
# species-file (``arc/statmech/arkane.py::species_input_template``). ``energy``, ``geometry`` and
# ``frequencies`` are the load-bearing trio Arkane's statmech actually consumes, so they are
# mandatory. A real artifact also carries ``linear`` and ``spinMultiplicity``, but those are not
# required here, and ``rotors`` never is: no real artifact inspected while designing this module
# carried rotors at all. A file that merely parses as Python -- e.g. one shaped like the Arkane
# network DSL's ``transitionState(label=..., E0=...)`` -- is NOT an artifact and is classified
# unusable rather than usable.
REQUIRED_ARTIFACT_ASSIGNMENTS = ('energy', 'geometry', 'frequencies')


def _status_yml_entry_to_convergence(label: str, entry: dict) -> tuple:
    """
    Build the ``(converged, reason)`` pair for one ``status.yml`` entry.

    ``status.yml`` preserves the tri-state faithfully, so its own ``convergence`` field is
    trusted -- but only under STRICT identity: literally ``True`` means converged, literally
    ``False`` means explicitly failed, and ANYTHING else (``None``, a string like ``'false'``
    from a hand-edited or differently-serialized file, an int, a list, ...) means UNKNOWN.
    ``bool(...)`` coercion is never applied: the string ``'false'`` is truthy, so coercing it
    would silently promote a failed transition state to converged.

    Args:
        label (str): The ARC transition state label this entry describes.
        entry (dict): The raw per-label mapping from ``status.yml``.

    Returns:
        tuple: ``(converged, reason)``, where ``converged`` is ``True``/``False``/``None``.
    """
    raw = entry.get('convergence')
    errors = (entry.get('errors') or '').strip()
    info = (entry.get('info') or '').strip()
    detail = '; '.join(part for part in (info, errors) if part)
    if raw is True:
        converged = True
        reason = f"ARC reports '{label}' as converged (status.yml)."
    elif raw is False:
        converged = False
        reason = f"ARC reports '{label}' as explicitly not converged (status.yml)."
    elif raw is None:
        converged = None
        reason = f"ARC has no convergence verdict for '{label}' yet (status.yml); treating it as unknown."
    else:
        converged = None
        reason = (f"ARC's status.yml carries an unexpected 'convergence' value {raw!r} "
                  f"(type {type(raw).__name__}) for '{label}'; refusing to coerce it, so convergence "
                  f"is treated as unknown.")
    if detail:
        reason = f'{reason} Detail: {detail}.'
    return converged, reason


def _output_yml_entry_to_convergence(label: str, entry: dict) -> tuple:
    """
    Build the ``(converged, reason)`` pair for one ``output.yml`` entry.

    ``output.yml``'s ``'converged'`` field is ``entry.get('convergence') is True``
    (``arc/output.py::_spc_to_dict``), which collapses ``status.yml``'s tri-state: an explicitly
    failed job and a job that was never assessed both read back as ``converged: False`` here.
    That loss of information is called out in the returned reason text rather than silently
    reported as if it were a genuine "not converged" verdict.

    The ``'converged'`` value is read under STRICT identity, never ``bool(...)`` coercion:
    literally ``True`` means converged, literally ``False`` means not converged (with the
    lossy-collapse caveat above), and ANYTHING else (a string like ``'false'``, an int, ``None``,
    ...) means UNKNOWN -- the string ``'false'`` is truthy, so coercing it would silently promote
    a failed transition state to converged.

    Args:
        label (str): The ARC transition state label this entry describes.
        entry (dict): The raw per-label mapping from ``output.yml``.

    Returns:
        tuple: ``(converged, reason)``. ``converged`` is ``True``, ``False``, or ``None`` (the
            latter only for an unexpected, non-boolean ``'converged'`` value); the reason text for
            a ``False`` explains that it may mean "unknown" rather than "failed".
    """
    raw = entry.get('converged')
    if raw is True:
        converged = True
        reason = f"ARC reports '{label}' as converged (output.yml)."
    elif raw is False:
        converged = False
        reason = (f"ARC reports '{label}' as not converged (output.yml); output.yml could not "
                  f"distinguish 'explicitly failed' from 'not yet determined' for this label, so "
                  f"this verdict is treated fail-closed as not converged, but it may only mean "
                  f"the job was never assessed.")
    else:
        converged = None
        reason = (f"ARC's output.yml carries an unexpected 'converged' value {raw!r} "
                  f"(type {type(raw).__name__}) for '{label}'; refusing to coerce it, so convergence "
                  f"is treated as unknown.")
    return converged, reason


def read_arc_convergence(arc_project_directory: str) -> dict | None:
    """
    Read whichever ARC status source exists and normalize it into one per-label mapping.

    ``output/status.yml`` is preferred when present: it is the richer source, preserving ARC's
    true tri-state convergence verdict (``True``/``False``/``None``). ``output/output.yml`` is
    read as a fallback when ``status.yml`` is absent -- it collapses "explicitly failed" and
    "never assessed" into the same ``False``, which is exactly the distinction this module's
    callers need. See the module docstring for the full rationale.

    One exception to that preference, when BOTH sources exist: if ``status.yml`` says a label
    converged while ``output.yml`` -- ARC's later, consolidated summary -- says it did NOT, a
    stale ``status.yml`` left behind by an interrupted earlier run could be promoting a
    transition state ARC's final output says failed. That one direction fails closed: the label
    is treated as not converged, with a reason naming both sources and the disagreement -- UNLESS
    ``status.yml`` is strictly newer (by file mtime) than ``output.yml``, in which case
    ``status.yml``'s converged verdict is trusted instead, since ARC's unconditional,
    always-first ``status.yml`` write means a newer ``status.yml`` can only come from a later
    restarted run, not a stale one. The reverse direction (``status.yml`` not converged vs
    ``output.yml`` converged) simply keeps ``status.yml``'s richer verdict, and
    ``status.yml: null`` vs ``output.yml: false`` is the EXPECTED lossy collapse (``output.yml``
    folds ``None`` into ``False``), not a conflict -- neither of those overrides anything.

    Args:
        arc_project_directory (str): The ARC project directory to read from.

    Returns:
        dict | None: A mapping of ARC transition state label to
            ``{'converged': True | False | None, 'reason': str}``, or ``None`` if neither status
            source exists -- meaning convergence is entirely UNKNOWN for every label, not that
            every label is unconverged.
    """
    output_yml_path = os.path.join(arc_project_directory, 'output', 'output.yml')
    status_yml_path = os.path.join(arc_project_directory, 'output', 'status.yml')

    if os.path.isfile(status_yml_path):
        content = read_yaml_file(status_yml_path) or dict()
        result = dict()
        for label, entry in content.items():
            converged, reason = _status_yml_entry_to_convergence(label, entry or dict())
            result[label] = {'converged': converged, 'reason': reason}
        if os.path.isfile(output_yml_path):
            output_content = read_yaml_file(output_yml_path) or dict()
            for entry in (output_content.get('transition_states') or list()):
                label = entry.get('label')
                if not label or label not in result:
                    continue
                output_converged, _ = _output_yml_entry_to_convergence(label, entry)
                if result[label]['converged'] is True and output_converged is False:
                    # ARC writes status.yml unconditionally, BEFORE output.yml, within a single
                    # execute() call -- so in any run that completes normally, output.yml is never
                    # older than status.yml. If status.yml is nonetheless demonstrably NEWER, it
                    # cannot be the stale side of this disagreement: it is evidence of a later,
                    # restarted run that got further than whatever left output.yml behind (e.g. a
                    # restart that reached the status.yml write but crashed again before
                    # write_output_yml). Trust status.yml's converged=True only in that case;
                    # otherwise fail closed as before, since a stale status.yml from an
                    # interrupted earlier run must not promote a transition state ARC's final
                    # consolidated output says failed.
                    try:
                        status_is_newer = os.path.getmtime(status_yml_path) > os.path.getmtime(output_yml_path)
                    except OSError:
                        status_is_newer = False
                    if status_is_newer:
                        result[label] = {
                            'converged': True,
                            'reason': (f"status.yml and output.yml disagree about '{label}': status.yml reports "
                                       f"it as converged, and output.yml (ARC's later, consolidated summary) "
                                       f"reports it as NOT converged. status.yml is strictly newer than "
                                       f"output.yml, so it is trusted as the result of a later restarted run; "
                                       f"the disagreement is not treated as a stale-status.yml conflict."),
                        }
                    else:
                        result[label] = {
                            'converged': False,
                            'reason': (f"status.yml and output.yml disagree about '{label}': status.yml reports it "
                                       f"as converged, but output.yml (ARC's later, consolidated summary) reports "
                                       f"it as NOT converged. A stale status.yml from an interrupted earlier run "
                                       f"could otherwise promote a failed transition state, so this is treated "
                                       f"fail-closed as not converged."),
                        }
        return result

    if os.path.isfile(output_yml_path):
        content = read_yaml_file(output_yml_path) or dict()
        result = dict()
        for entry in (content.get('transition_states') or list()):
            label = entry.get('label')
            if not label:
                continue
            converged, reason = _output_yml_entry_to_convergence(label, entry)
            result[label] = {'converged': converged, 'reason': reason}
        return result

    return None


def _confine_to_project(arc_project_directory: str, path: str) -> None:
    """
    Refuse a path that resolves outside the ARC project directory.

    Mirrors the confinement check in ``t3.runners.rmg_runner.run_arkane_job``: the check is done
    against ``os.path.realpath`` of both sides (so a traversal component or a symlink cannot
    escape it), but the caller keeps using the original, unresolved path afterward.

    Args:
        arc_project_directory (str): The ARC project directory a resolved artifact must stay under.
        path (str): The candidate path to check.

    Raises:
        ValueError: If ``path`` resolves to a location outside ``arc_project_directory``.
    """
    resolved_project = os.path.realpath(arc_project_directory)
    resolved_path = os.path.realpath(path)
    if resolved_path == resolved_project \
            or os.path.commonpath([resolved_project, resolved_path]) != resolved_project:
        raise ValueError(f"Refusing to read '{path}': it resolves to '{resolved_path}', which is outside the "
                         f"ARC project directory '{arc_project_directory}' (resolved to '{resolved_project}').")


def _check_artifact_shape(ts_label: str, artifact_path: str) -> None:
    """
    Refuse an artifact that is not shaped like an ARC statmech species-file.

    A file can parse cleanly as Python and resolve every ``Log(...)`` it references while still
    being the wrong KIND of file entirely (e.g. an Arkane network-DSL ``transitionState(...)``
    call). The ARC shape is a set of module-level assignments; the mandatory, load-bearing ones
    are ``REQUIRED_ARTIFACT_ASSIGNMENTS`` (see the comment there for what is deliberately NOT
    required).

    Args:
        ts_label (str): The transition state label the artifact belongs to (used in error text).
        artifact_path (str): The path to the candidate artifact file.

    Raises:
        ValueError: If any of ``REQUIRED_ARTIFACT_ASSIGNMENTS`` is not assigned at the artifact's
            module level, naming exactly the missing ones.
    """
    with open(artifact_path, 'r') as f:
        content = f.read()
    tree = _parse_as_ast(text=content, path=artifact_path)
    assigned = set()
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    assigned.add(target.id)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            assigned.add(node.target.id)
    missing = [name for name in REQUIRED_ARTIFACT_ASSIGNMENTS if name not in assigned]
    if missing:
        raise ValueError(f"Transition state '{ts_label}''s file at '{artifact_path}' is not shaped like an ARC "
                         f"statmech species-file: it is missing the mandatory module-level assignment(s) "
                         f"{missing} (a real artifact assigns "
                         f"{', '.join(repr(name) for name in REQUIRED_ARTIFACT_ASSIGNMENTS)}; see "
                         f"arc/statmech/arkane.py::species_input_template).")


def _extract_sensitivity(evidence) -> tuple:
    """
    Pull ``(coefficient, delta_ln_k)`` out of a sensitivity evidence object, if any.

    Args:
        evidence: Either ``None``, a ``t3.pdep.selector.SensitiveTransitionState`` (or any object
            exposing ``.coefficient``/``.delta_ln_k`` attributes), or a ``(coefficient,
            delta_ln_k)`` pair.

    Returns:
        tuple: ``(coefficient, delta_ln_k)``, either of which may be ``None``.
    """
    if evidence is None:
        return None, None
    if hasattr(evidence, 'coefficient') and hasattr(evidence, 'delta_ln_k'):
        return evidence.coefficient, evidence.delta_ln_k
    coefficient, delta_ln_k = evidence
    return coefficient, delta_ln_k


class TSArtifactRecord:
    """
    Whether one network transition state's ARC quantum chemistry is usable, and why.

    Args:
        network_id (str): The network file stem, e.g. ``'network4_1'``.
        network_ts_label (str): The transition state label in the network file's namespace.
        arc_ts_label (str, optional): The ARC transition state label this record was resolved
            against (always the recomputed one, never the one stored in the join record).
        status (str): One of ``TS_ARTIFACT_STATUSES``.
        artifact_path (str, optional): The path of the artifact, when one was resolved. Carried for
            BOTH ``ARTIFACT_STATUS_USABLE`` (converged, ready to use) and
            ``ARTIFACT_STATUS_UNVERIFIED`` (parses/resolves fine, but convergence could not be
            determined) records; ``None`` for every other status. WARNING: a truthy
            ``artifact_path`` is NOT a usability test. Code that gates on
            ``if record.artifact_path:`` instead of ``record.status == ARTIFACT_STATUS_USABLE``
            will silently treat an UNVERIFIED record -- one whose convergence is unknown -- as if
            it were USABLE. Always branch on ``status``.
        converged (bool, optional): ARC's convergence verdict for this transition state, or
            ``None`` when it is unknown (never assessed, or no status source exists at all).
        reason (str): Why this record has the status it has.
        coefficient (float, optional): The signed sensitivity coefficient, when sensitivity
            evidence was supplied for this transition state.
        delta_ln_k (float, optional): The corresponding dimensionless rate response, when
            sensitivity evidence was supplied for this transition state.
        path_reaction_labels (tuple): The labels of the path reactions this transition state owns,
            carried through from the ``t3.pdep.join.TSJoinRecord`` this record was resolved from.
            This is the structural identity of the transition state -- unlike ``network_ts_label``
            (positional: Arkane numbers ``TS1``, ``TS2``, ... by path-reaction index, so pruning a
            reaction shifts every later label) or ``network_id`` (Arkane's final network artifact is
            named by explorer run index, e.g. ``network0_full``, not by any stem a caller controls),
            this tuple identifies the actual saddle point and is what any cross-capture matching
            (e.g. reusing a prior capture's QM for a differently-named network) must key on.
    """

    def __init__(self,
                 network_id: str,
                 network_ts_label: str,
                 arc_ts_label: str | None,
                 status: str,
                 artifact_path: str | None = None,
                 converged: bool | None = None,
                 reason: str = '',
                 coefficient: float | None = None,
                 delta_ln_k: float | None = None,
                 path_reaction_labels: tuple = tuple(),
                 ):
        if status not in TS_ARTIFACT_STATUSES:
            raise ValueError(f"Unrecognized transition state artifact status '{status}'; "
                             f"expected one of {list(TS_ARTIFACT_STATUSES)}.")
        self.network_id = network_id
        self.network_ts_label = network_ts_label
        self.arc_ts_label = arc_ts_label
        self.status = status
        self.artifact_path = artifact_path
        self.converged = converged
        self.reason = reason
        self.coefficient = coefficient
        self.delta_ln_k = delta_ln_k
        self.path_reaction_labels = tuple(path_reaction_labels)

    @property
    def key(self) -> tuple:
        """
        The identity of this record: the network transition state it describes.

        Returns:
            tuple: ``(network_id, network_ts_label)``.
        """
        return self.network_id, self.network_ts_label

    def as_dict(self) -> dict:
        """
        Render this record as plain YAML-safe types.

        Returns:
            dict: A plain dict.
        """
        return {
            'network_id': self.network_id,
            'network_ts_label': self.network_ts_label,
            'arc_ts_label': self.arc_ts_label,
            'status': self.status,
            'artifact_path': self.artifact_path,
            'converged': self.converged,
            'reason': self.reason,
            'coefficient': self.coefficient,
            'delta_ln_k': self.delta_ln_k,
            'path_reaction_labels': list(self.path_reaction_labels),
        }

    @classmethod
    def from_dict(cls, record_dict: dict) -> 'TSArtifactRecord':
        """
        Reconstruct a record from its ``as_dict()`` rendering.

        Args:
            record_dict (dict): The rendered record.

        Raises:
            ValueError: If a required field is missing, or the status is unrecognized.

        Returns:
            TSArtifactRecord: The reconstructed record.
        """
        for required in ('network_id', 'network_ts_label', 'status'):
            if not record_dict.get(required):
                raise ValueError(f"A transition state artifact record is missing the required field "
                                 f"'{required}': {record_dict}")
        return cls(
            network_id=record_dict['network_id'],
            network_ts_label=record_dict['network_ts_label'],
            arc_ts_label=record_dict.get('arc_ts_label'),
            status=record_dict['status'],
            artifact_path=record_dict.get('artifact_path'),
            converged=record_dict.get('converged'),
            reason=record_dict.get('reason') or '',
            coefficient=record_dict.get('coefficient'),
            delta_ln_k=record_dict.get('delta_ln_k'),
            path_reaction_labels=tuple(record_dict.get('path_reaction_labels') or ()),
        )

    def __repr__(self) -> str:
        return f'TSArtifactRecord({self.network_id}/{self.network_ts_label} -> {self.arc_ts_label}, {self.status})'

    def __eq__(self, other) -> bool:
        if not isinstance(other, TSArtifactRecord):
            return NotImplemented
        return self.as_dict() == other.as_dict()


def discover_ts_artifacts(join_records: list, arc_project_directory: str, sensitivity_by_ts: dict | None = None,
                          ) -> list:
    """
    Decide, per network transition state, whether ARC produced a usable quantum-chemistry artifact.

    The ARC label and expected artifact path are always recomputed from ``(network_id,
    network_ts_label)`` -- never read from the join record's stored values -- so a stale or
    hand-edited sidecar cannot redirect which file gets read back for a given transition state.
    Any disagreement between the stored and recomputed *path* is merely noted in the reason text
    (the recomputed path is used regardless). A disagreement between the stored and recomputed
    *label*, when the stored label is present, is treated much more seriously: it means the join
    record's identity itself cannot be trusted, so the record is refused outright
    (``ARTIFACT_STATUS_UNUSABLE``) rather than silently normalized to the recomputed label. A
    stored label of ``None`` is not an error -- it is the legitimate, deliberate value for
    ``not_queued`` records -- only a present-and-different value is refused.

    ``already_present`` and ``not_queued`` join records are passed through with their original
    status and reason preserved, without inventing an artifact path for them: neither one was
    queued in this iteration, so there is nothing new here to discover about them.

    Args:
        join_records (list): The ``t3.pdep.join.TSJoinRecord`` entries describing what was queued.
        arc_project_directory (str): The ARC project directory the join records were queued into.
        sensitivity_by_ts (dict, optional): A mapping of ``(network_id, network_ts_label)`` (the
            same key ``TSJoinRecord.key`` uses -- never the bare transition state label, since two
            different networks can both have a ``'TS1'``) to sensitivity evidence for that
            transition state -- either a ``t3.pdep.selector.SensitiveTransitionState`` (or any
            object exposing ``.coefficient``/``.delta_ln_k``) or a ``(coefficient, delta_ln_k)``
            pair. Used only to carry that evidence through onto the returned records; it plays no
            part in the usability decision itself.

    Returns:
        list: One ``TSArtifactRecord`` per join record, in the same order.
    """
    sensitivity_by_ts = sensitivity_by_ts or dict()
    convergence = read_arc_convergence(arc_project_directory)
    records = list()

    for join_record in join_records:
        coefficient, delta_ln_k = _extract_sensitivity(sensitivity_by_ts.get(join_record.key))

        if join_record.status == JOIN_STATUS_NOT_QUEUED:
            records.append(TSArtifactRecord(
                network_id=join_record.network_id,
                network_ts_label=join_record.network_ts_label,
                path_reaction_labels=join_record.path_reaction_labels,
                arc_ts_label=join_record.arc_ts_label,
                status=ARTIFACT_STATUS_NOT_QUEUED,
                artifact_path=None,
                converged=None,
                reason=join_record.reason,
                coefficient=coefficient,
                delta_ln_k=delta_ln_k,
            ))
            continue

        if join_record.status == JOIN_STATUS_ALREADY_PRESENT:
            records.append(TSArtifactRecord(
                network_id=join_record.network_id,
                network_ts_label=join_record.network_ts_label,
                path_reaction_labels=join_record.path_reaction_labels,
                arc_ts_label=join_record.arc_ts_label,
                status=ARTIFACT_STATUS_MISSING,
                artifact_path=None,
                converged=None,
                reason=join_record.reason,
                coefficient=coefficient,
                delta_ln_k=delta_ln_k,
            ))
            continue

        # JOIN_STATUS_QUEUED: recompute the label and the expected path; never trust the stored ones.
        # An explicit raise, never an `assert` (which vanishes under `python -O`, after which an
        # unrecognized status would fall through into this queued branch and be handed a
        # recomputed artifact path it must never have).
        if join_record.status != JOIN_STATUS_QUEUED:
            raise ValueError(f"Unrecognized join record status '{join_record.status}' for "
                             f"'{join_record.network_id}/{join_record.network_ts_label}'; expected "
                             f"'{JOIN_STATUS_QUEUED}', '{JOIN_STATUS_ALREADY_PRESENT}' or "
                             f"'{JOIN_STATUS_NOT_QUEUED}'.")
        recomputed_label = arc_ts_label(join_record.network_id, join_record.network_ts_label)
        recomputed_path = expected_ts_artifact_path(arc_project_directory, recomputed_label)

        reason_parts = list()
        if join_record.expected_artifact_path and join_record.expected_artifact_path != recomputed_path:
            reason_parts.append(f"The stored expected artifact path '{join_record.expected_artifact_path}' "
                                f"disagrees with the recomputed path '{recomputed_path}'; the stored path is "
                                f"advisory only (a hand-edited sidecar could point anywhere), so the recomputed "
                                f"path is the one actually read.")

        # Fail-closed identity check: a stored `arc_ts_label` that is PRESENT but disagrees with
        # the recomputed one means this join record no longer names the transition state it
        # claims to (a hand-edited or stale sidecar). Unlike the path disagreement above (merely
        # advisory -- the recomputed path is used regardless), a label disagreement is refused
        # outright rather than silently normalized away: `arc_ts_label` is legitimately `None` for
        # `not_queued` records, so only a truthy-and-different value is an error here.
        if join_record.arc_ts_label and join_record.arc_ts_label != recomputed_label:
            reason_parts.append(f"The stored ARC transition state label '{join_record.arc_ts_label}' disagrees "
                                f"with the recomputed label '{recomputed_label}'; this join record's identity "
                                f"cannot be trusted, so it is refused rather than silently normalized to the "
                                f"recomputed label.")
            records.append(TSArtifactRecord(
                network_id=join_record.network_id,
                network_ts_label=join_record.network_ts_label,
                path_reaction_labels=join_record.path_reaction_labels,
                arc_ts_label=recomputed_label,
                status=ARTIFACT_STATUS_UNUSABLE,
                artifact_path=None,
                converged=None,
                reason=' '.join(reason_parts),
                coefficient=coefficient,
                delta_ln_k=delta_ln_k,
            ))
            continue

        try:
            _confine_to_project(arc_project_directory, recomputed_path)
        except ValueError as e:
            reason_parts.append(str(e))
            records.append(TSArtifactRecord(
                network_id=join_record.network_id,
                network_ts_label=join_record.network_ts_label,
                path_reaction_labels=join_record.path_reaction_labels,
                arc_ts_label=recomputed_label,
                status=ARTIFACT_STATUS_UNUSABLE,
                artifact_path=None,
                converged=None,
                reason=' '.join(reason_parts),
                coefficient=coefficient,
                delta_ln_k=delta_ln_k,
            ))
            continue

        converged_entry = (convergence or dict()).get(recomputed_label)
        if converged_entry is not None:
            converged_flag = converged_entry['converged']
            reason_parts.append(converged_entry['reason'])
        else:
            converged_flag = None
            if convergence is None:
                reason_parts.append(f"No ARC status source (output.yml or status.yml) was found under "
                                    f"'{arc_project_directory}'; convergence for '{recomputed_label}' is unknown.")
            else:
                reason_parts.append(f"'{recomputed_label}' does not appear in the ARC status source found under "
                                    f"'{arc_project_directory}'; convergence is unknown.")

        if converged_flag is False:
            records.append(TSArtifactRecord(
                network_id=join_record.network_id,
                network_ts_label=join_record.network_ts_label,
                path_reaction_labels=join_record.path_reaction_labels,
                arc_ts_label=recomputed_label,
                status=ARTIFACT_STATUS_UNUSABLE,
                artifact_path=None,
                converged=converged_flag,
                reason=' '.join(reason_parts),
                coefficient=coefficient,
                delta_ln_k=delta_ln_k,
            ))
            continue

        if not os.path.isfile(recomputed_path):
            reason_parts.append(f"No Arkane statmech artifact was found at the expected path "
                                f"'{recomputed_path}'.")
            records.append(TSArtifactRecord(
                network_id=join_record.network_id,
                network_ts_label=join_record.network_ts_label,
                path_reaction_labels=join_record.path_reaction_labels,
                arc_ts_label=recomputed_label,
                status=ARTIFACT_STATUS_MISSING,
                artifact_path=None,
                converged=converged_flag,
                reason=' '.join(reason_parts),
                coefficient=coefficient,
                delta_ln_k=delta_ln_k,
            ))
            continue

        try:
            _read_qm_artifact(ts_label=recomputed_label, artifact_path=recomputed_path,
                              allowed_log_root=arc_project_directory)
            _check_artifact_shape(ts_label=recomputed_label, artifact_path=recomputed_path)
        except ValueError as e:
            reason_parts.append(f'The artifact exists but could not be used: {e}')
            records.append(TSArtifactRecord(
                network_id=join_record.network_id,
                network_ts_label=join_record.network_ts_label,
                path_reaction_labels=join_record.path_reaction_labels,
                arc_ts_label=recomputed_label,
                status=ARTIFACT_STATUS_UNUSABLE,
                artifact_path=None,
                converged=converged_flag,
                reason=' '.join(reason_parts),
                coefficient=coefficient,
                delta_ln_k=delta_ln_k,
            ))
            continue

        reason_parts.append(f"'{recomputed_path}' parses, has the ARC statmech species-file shape, and every "
                            f"Log(...) it references resolves.")
        if converged_flag is True:
            artifact_status = ARTIFACT_STATUS_USABLE
        else:
            # converged_flag is None here (the ``is False`` case already returned above): the
            # artifact itself is fine, but no status source ever confirmed convergence. Fail
            # closed on USABLE without failing hard on UNUSABLE -- see ARTIFACT_STATUS_UNVERIFIED.
            reason_parts.append("Convergence was never confirmed by an ARC status source, so this "
                                "artifact is not promoted to 'usable' even though it parses cleanly.")
            artifact_status = ARTIFACT_STATUS_UNVERIFIED
        records.append(TSArtifactRecord(
            network_id=join_record.network_id,
            network_ts_label=join_record.network_ts_label,
            path_reaction_labels=join_record.path_reaction_labels,
            arc_ts_label=recomputed_label,
            status=artifact_status,
            artifact_path=recomputed_path,
            converged=converged_flag,
            reason=' '.join(reason_parts),
            coefficient=coefficient,
            delta_ln_k=delta_ln_k,
        ))

    return records


class PDepHybridEvaluation:
    """
    Whether a P-dep network's selected transition states, taken together, support a hybrid QM+ILT
    network.

    Args:
        hybrid_status (str): One of ``HYBRID_STATUSES``.
        accepted (bool): Whether the network should proceed to a hybrid build.
        usable_count (int): The number of selected transition states with a usable artifact.
        selected_count (int): The total number of selected transition states evaluated.
        artifacts (list): The ``TSArtifactRecord`` entries this evaluation was computed from.
        dominant_sensitivity_ts_missing (bool, optional): ``True`` if the selected transition
            state with the single largest ``delta_ln_k`` is NOT usable; ``False`` if it is;
            ``None`` when the question cannot honestly be answered -- i.e. unless EVERY selected
            transition state carries a sensitivity value (an unscored one could itself be the
            dominant one) and the selected set has more than one member (dominance is meaningless
            for a singleton).
        reason (str): A human-readable summary of how this evaluation was reached.
    """

    def __init__(self,
                 hybrid_status: str,
                 accepted: bool,
                 usable_count: int,
                 selected_count: int,
                 artifacts: list,
                 dominant_sensitivity_ts_missing: bool | None = None,
                 reason: str = '',
                 ):
        if hybrid_status not in HYBRID_STATUSES:
            raise ValueError(f"Unrecognized hybrid status '{hybrid_status}'; expected one of "
                             f"{list(HYBRID_STATUSES)}.")
        self.hybrid_status = hybrid_status
        self.accepted = accepted
        self.usable_count = usable_count
        self.selected_count = selected_count
        self.artifacts = list(artifacts)
        self.dominant_sensitivity_ts_missing = dominant_sensitivity_ts_missing
        self.reason = reason

    def as_dict(self) -> dict:
        """
        Render this evaluation as plain YAML-safe types.

        Returns:
            dict: A plain dict.
        """
        return {
            'hybrid_status': self.hybrid_status,
            'accepted': self.accepted,
            'usable_count': self.usable_count,
            'selected_count': self.selected_count,
            'artifacts': [artifact.as_dict() for artifact in self.artifacts],
            'dominant_sensitivity_ts_missing': self.dominant_sensitivity_ts_missing,
            'reason': self.reason,
        }

    @classmethod
    def from_dict(cls, evaluation_dict: dict) -> 'PDepHybridEvaluation':
        """
        Reconstruct an evaluation from its ``as_dict()`` rendering.

        Args:
            evaluation_dict (dict): The rendered evaluation.

        Returns:
            PDepHybridEvaluation: The reconstructed evaluation.
        """
        return cls(
            hybrid_status=evaluation_dict['hybrid_status'],
            accepted=evaluation_dict['accepted'],
            usable_count=evaluation_dict['usable_count'],
            selected_count=evaluation_dict['selected_count'],
            artifacts=[TSArtifactRecord.from_dict(entry) for entry in (evaluation_dict.get('artifacts') or list())],
            dominant_sensitivity_ts_missing=evaluation_dict.get('dominant_sensitivity_ts_missing'),
            reason=evaluation_dict.get('reason') or '',
        )

    def __repr__(self) -> str:
        return (f'PDepHybridEvaluation({self.hybrid_status}, accepted={self.accepted}, '
                f'{self.usable_count}/{self.selected_count} usable)')

    def __eq__(self, other) -> bool:
        if not isinstance(other, PDepHybridEvaluation):
            return NotImplemented
        return self.as_dict() == other.as_dict()


def evaluate_pdep_hybrid(artifacts: list, allow_partial_selected_qm: bool = False) -> PDepHybridEvaluation:
    """
    Decide whether a network's selected transition states support a hybrid QM+ILT network.

    Strict by default: every selected transition state must be usable
    (``HYBRID_STATUS_COMPLETE``) for ``accepted`` to be ``True``. Setting
    ``allow_partial_selected_qm=True`` additionally accepts ``HYBRID_STATUS_PARTIAL_SELECTED_QM``
    (at least one, but not all, usable) -- this is an opt-in relaxation, never the default,
    because a hybrid network missing QM for some of its most uncertain transition states is a
    real degradation of what was asked for, not a free win. A network with zero usable transition
    states is never accepted, regardless of ``allow_partial_selected_qm``: there is nothing to
    build a hybrid network from.

    This function does not read or write ``t3.pdep.selector.PDepNetworkSelection.evaluation_status``
    at all -- that field answers a PRE-ARC question (was this network's sensitivity evaluated?)
    while this one answers a POST-ARC question (did the QM refinement it asked for come back?),
    and overloading one status field with both would erase that distinction.

    Args:
        artifacts (list): The ``TSArtifactRecord`` entries for one network's selected transition
            states (as returned by ``discover_ts_artifacts``).
        allow_partial_selected_qm (bool, optional): Whether a partially-usable set of selected
            transition states should still be accepted.

    Raises:
        ValueError: If the records span more than one ``network_id`` (an evaluation describes ONE
            network); if the same ``(network_id, network_ts_label)`` appears more than once (the
            coverage counts would double-count it); or if a record claims to be usable but carries
            no ``artifact_path`` (nothing downstream could actually consume it). Refusing loudly
            beats silently producing a wrong verdict from malformed input.

    Returns:
        PDepHybridEvaluation: The coverage, status, acceptance decision, and (when sensitivity
            evidence is available) whether the most sensitivity-dominant selected transition
            state is among the ones missing.
    """
    network_ids = sorted({artifact.network_id for artifact in artifacts})
    if len(network_ids) > 1:
        raise ValueError(f'An evaluation describes ONE network, but the supplied records span '
                         f'{len(network_ids)} networks: {network_ids}. Refusing to evaluate them together.')
    seen_keys = set()
    for artifact in artifacts:
        if artifact.key in seen_keys:
            raise ValueError(f"The transition state '{artifact.network_id}/{artifact.network_ts_label}' appears "
                             f"more than once in the supplied records; refusing to evaluate a set that would "
                             f"double-count it.")
        seen_keys.add(artifact.key)
        if artifact.status == ARTIFACT_STATUS_USABLE and not artifact.artifact_path:
            raise ValueError(f"The record for '{artifact.network_id}/{artifact.network_ts_label}' claims status "
                             f"'{ARTIFACT_STATUS_USABLE}' but carries no artifact_path; nothing downstream could "
                             f"actually consume it, so this record is refused rather than counted as usable.")

    selected_count = len(artifacts)
    usable = [artifact for artifact in artifacts if artifact.status == ARTIFACT_STATUS_USABLE]
    usable_count = len(usable)

    if selected_count == 0 or usable_count == 0:
        hybrid_status = HYBRID_STATUS_NOT_EVALUATED
    elif usable_count == selected_count:
        hybrid_status = HYBRID_STATUS_COMPLETE
    else:
        hybrid_status = HYBRID_STATUS_PARTIAL_SELECTED_QM

    if hybrid_status == HYBRID_STATUS_COMPLETE:
        accepted = True
    elif hybrid_status == HYBRID_STATUS_PARTIAL_SELECTED_QM:
        accepted = bool(allow_partial_selected_qm)
    else:
        accepted = False

    # The transition state with the single largest sensitivity response (delta_ln_k) is the one
    # whose missing QM would matter most; flag it explicitly rather than leaving it buried in a
    # per-TS artifact list. The flag is only ever asserted (True/False) on COMPLETE evidence:
    # if even one selected transition state carries no sensitivity value, the unscored one could
    # itself be the dominant one, so the honest answer is UNKNOWN (None), never a quiet False
    # computed from the scored subset alone. It is likewise None for a singleton selected set,
    # where "the most sensitive one" is meaningless.
    dominant_sensitivity_ts_missing = None
    if len(artifacts) > 1 and all(artifact.delta_ln_k is not None for artifact in artifacts):
        dominant = max(artifacts, key=lambda artifact: artifact.delta_ln_k)
        dominant_sensitivity_ts_missing = dominant.status != ARTIFACT_STATUS_USABLE

    reason = (f'{usable_count}/{selected_count} selected transition state(s) have a usable QM artifact '
             f'({hybrid_status}); accepted={accepted}.')
    if dominant_sensitivity_ts_missing:
        reason += (' The selected transition state with the largest sensitivity response '
                  '(delta_ln_k) is among those without a usable artifact.')

    return PDepHybridEvaluation(
        hybrid_status=hybrid_status,
        accepted=accepted,
        usable_count=usable_count,
        selected_count=selected_count,
        artifacts=artifacts,
        dominant_sensitivity_ts_missing=dominant_sensitivity_ts_missing,
        reason=reason,
    )
