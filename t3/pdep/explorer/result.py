"""
t3 pdep explorer result module.

Defines ``PDepExplorationResult``, the plain-data record a PES exploration decision (skip / run-
and-fail / run-and-succeed) is packaged into for callers outside the adapter layer. It mirrors the
house style set by ``t3.pdep.selector.PDepNetworkSelection``: a frozen dataclass with a status
string, evidence fields, and an ``as_dict()`` that renders everything as plain JSON/YAML-safe
types.
"""

import copy

from dataclasses import dataclass, field

from t3.pdep.parser import to_json_safe
from t3.pdep.selector import PDepNetworkSelection

# Exploration status values for PDepExplorationResult.status. THREE statuses, not a bool, mirroring
# the never-ran-vs-failed distinction ArkaneExplorerAdapter already preserves with its
# ``succeeded is None`` sentinel (t3/pdep/explorer/arkane.py, ~line 147):
#
# - 'skipped'   the exploration was never attempted, because the qualification gate (a
#               non-qualifying PDepNetworkSelection) declined it. No Arkane job ran at all. See
#               ``admission_policy`` below: that gate is the ONLY thing that can decline a run.
# - 'failed'    an Arkane exploration actually ran and did not succeed.
# - 'succeeded' it ran and succeeded.
#
# Collapsing 'skipped' into 'failed' would report "we tried and Arkane failed" about a run that
# never happened, which is a wrong statement of fact a caller could act on (e.g. counting it
# against an Arkane failure-rate budget it never touched).
EXPLORATION_STATUS_SUCCEEDED = 'succeeded'
EXPLORATION_STATUS_FAILED = 'failed'
EXPLORATION_STATUS_SKIPPED = 'skipped'

_VALID_EXPLORATION_STATUSES = (EXPLORATION_STATUS_SUCCEEDED, EXPLORATION_STATUS_FAILED,
                               EXPLORATION_STATUS_SKIPPED)

# WHAT ADMITTED this exploration -- recorded because ``status`` alone stopped being able to say it.
# Qualification used to be a gate: an unqualified selection meant "do not explore", so a 'succeeded'
# result carrying a selection implied that selection had qualified. Under the ranking reframe the
# spend decision can instead be made by the caller (T3 makes it in ``t3.main`` via
# ``t3.pdep.budget.apply_pdep_qm_budget``, which ranks the whole field of qualified networks against
# a budget). A succeeded result carrying an UNQUALIFIED selection is therefore now a legitimate
# outcome -- and indistinguishable, after the fact, from the gate having been bypassed, unless the
# basis is recorded alongside it.
#
# - 'qualified_selection'  the exploration ran (or was skipped) under the qualification gate: a
#                          selection was the admission authority, and an unqualified one declined it.
# - 'caller_admitted'      the caller admitted this network itself and passed the selection only to
#                          bind the run to the content the decision was made about. The
#                          qualification checks stand aside; the provenance checks do not.
# - 'ungated'              no selection was given at all (``explore_pdep_network(selection=None)``),
#                          so nothing gated the run and there is no recorded decision behind it.
#
# 'ungated' exists because the first two do not cover the selection-less call, and defaulting it to
# 'qualified_selection' would have the record assert that a qualified selection admitted a run for
# which no selection existed -- a false provenance claim, and the most consequential kind, since it
# is exactly the claim someone auditing an expensive QM run would rely on.
ADMISSION_POLICY_QUALIFIED_SELECTION = 'qualified_selection'
ADMISSION_POLICY_CALLER_ADMITTED = 'caller_admitted'
ADMISSION_POLICY_UNGATED = 'ungated'

# Public, unlike the status tuple above: t3.pdep.api validates this same field before it ever builds
# a result, and one shared tuple is what keeps the two checks from drifting apart. Note that only the
# first two are REQUESTABLE as an ``explore_pdep_network()`` argument -- 'ungated' is not a policy a
# caller chooses, it is what passing no selection MEANS, so that function derives it rather than
# accepting it (see its own REQUESTABLE_ADMISSION_POLICIES).
VALID_ADMISSION_POLICIES = (ADMISSION_POLICY_QUALIFIED_SELECTION, ADMISSION_POLICY_CALLER_ADMITTED,
                            ADMISSION_POLICY_UNGATED)

# The SHAPE of one PDepExplorationResult.as_dict() record: its set of keys and their types. This is
# its ONE job. It is deliberately distinct from three other version constants that already exist
# for other jobs, none of which this one duplicates: ``t3.pdep.cache.SA_CACHE_CONTRACT_VERSION``
# (whether a cached SA YAML on disk may still be trusted), ``t3.pdep.selector.SELECTION_SCHEMA_VERSION``
# (the SHAPE of one saved *selection* record), and ``t3.pdep.selector.SELECTION_ALGORITHM_VERSION``
# (the selection *decision logic*). A change to any of those three does not change the shape of an
# exploration-result record, so it must not force a bump here -- and a change to this record's shape
# (a field added, removed, renamed, or rendered differently) does not touch any of those three.
# Version 1 means "the shape as of first ship"; pre-ship development churn on this shape is
# deliberately NOT versioned, because t3.pdep has never shipped a release -- no saved record has
# ever left this repo under a version number that would need to keep meaning the same thing.
EXPLORATION_RESULT_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class PDepExplorationResult:
    """
    The outcome of one PES exploration decision for one network: skipped, failed, or succeeded.

    ``__post_init__`` enforces fail-closed coherence between ``status`` and the other fields (see
    below), so a caller reading only ``status`` cannot be misled by leftover artifacts from a
    different code path, and a caller reading the artifact fields cannot be misled into thinking
    something ran when it did not.

    Args:
        network_id (str, optional): The network file stem this result is for, e.g. ``'network4_2'``,
            or ``None`` if the exploration never got far enough to resolve one.
        status (str): One of ``'succeeded'``, ``'failed'``, ``'skipped'`` (see the
            ``EXPLORATION_STATUS_*`` constants). ``'skipped'`` means the exploration was never
            attempted because the qualification gate declined it; ``'failed'`` means it ran and did
            not succeed; ``'succeeded'`` means it ran and succeeded.
        reasons (tuple): Human-readable reasons for a ``'failed'`` or ``'skipped'`` outcome. Required
            (non-empty) for both: a negative outcome with no stated reason is unusable to a caller
            trying to decide what to do next.
        network_paths (tuple): The resolved final network file path(s) an exploration wrote, e.g.
            from ``ArkaneExplorerAdapter.get_networks()``. Must be empty unless ``status`` is
            ``'succeeded'``, and non-empty when it is (success with no artifact is a fail-open shape).
        output_paths (tuple): The raw Arkane ``output*.py`` path(s) an exploration wrote. Must be
            empty when ``status`` is ``'skipped'`` (nothing ran), but is CARRIED on ``'failed'`` as
            well as ``'succeeded'``: a failed run's own logs and partial output are exactly what
            someone diagnosing the failure needs, and this result is the only record they get of
            where to find them.
        k_tp_as_written (tuple): The parsed ``PDepArkaneReaction`` entries for the explored network(s),
            e.g. from ``ArkaneExplorerAdapter.get_k_tp()``. DIRECTION IS ARKANE'S, NOT THE CALLER'S:
            each entry carries ``reactants``/``products`` exactly as Arkane wrote them, which is NOT
            guaranteed to match the direction of anything the caller asked about (see
            ``ArkaneExplorerAdapter.get_k_tp``'s docstring for the full reasoning). A consumer that
            needs a rate in a particular direction must resolve the direction itself and reverse the
            entry when required; these entries are deliberately not canonicalized or reversed here.
            Must be empty unless ``status`` is ``'succeeded'``.
        manifest (dict): The exploration run's manifest (provenance: input hashes, tool versions,
            etc.), if any.
        selection (PDepNetworkSelection, optional): The network decision this exploration result is
            downstream of, if any (e.g. the non-qualifying decision that produced a ``'skipped'``
            result). Whether that decision was the thing that ADMITTED the run is a separate
            question, answered by ``admission_policy``.
        admission_policy (str): What admitted this exploration -- one of the ``ADMISSION_POLICY_*``
            constants; see their definition above for the full reasoning. Under the default,
            ``'qualified_selection'``, the selection was the admission authority, so a ``'skipped'``
            result means it declined. Under ``'caller_admitted'`` the caller made the spend decision
            itself, which is what makes a ``'succeeded'`` result carrying an UNQUALIFIED selection a
            coherent record rather than evidence of a bypassed gate. ``'ungated'`` means no
            ``selection`` was given, so nothing gated the run. A ``'skipped'`` status is cross-checked
            against this field (see Raises): only the qualification gate can decline an exploration,
            so it is the only basis a skipped result can name.

    Raises:
        ValueError: If ``status`` is not one of the ``EXPLORATION_STATUS_*`` constants; if
            ``admission_policy`` is not one of the ``ADMISSION_POLICY_*`` constants; if ``status``
            is ``'skipped'`` and ``admission_policy`` is not ``'qualified_selection'``; if
            ``status`` is ``'skipped'`` and any of ``network_paths``/``output_paths``/
            ``k_tp_as_written`` is non-empty (nothing ran, so there can be no artifacts); if
            ``status`` is ``'failed'`` or ``'skipped'`` and ``reasons`` is empty; if ``status`` is
            ``'succeeded'`` and ``network_paths`` is empty; or if ``status`` is ``'failed'`` and
            ``network_paths``/``k_tp_as_written`` is non-empty (those assert a usable result the
            failure denies).
    """
    network_id: str | None
    status: str
    reasons: tuple = tuple()
    network_paths: tuple = tuple()
    output_paths: tuple = tuple()
    k_tp_as_written: tuple = tuple()
    manifest: dict = field(default_factory=dict)
    selection: 'PDepNetworkSelection | None' = None
    admission_policy: str = ADMISSION_POLICY_QUALIFIED_SELECTION

    def __post_init__(self):
        if self.status not in _VALID_EXPLORATION_STATUSES:
            raise ValueError(
                f'PDepExplorationResult.status must be one of {_VALID_EXPLORATION_STATUSES}, '
                f'got {self.status!r}.')
        if self.admission_policy not in VALID_ADMISSION_POLICIES:
            raise ValueError(
                f'PDepExplorationResult.admission_policy must be one of '
                f'{VALID_ADMISSION_POLICIES}, got {self.admission_policy!r}.')
        if self.status == EXPLORATION_STATUS_SKIPPED \
                and self.admission_policy != ADMISSION_POLICY_QUALIFIED_SELECTION:
            # ``admission_policy`` answers "what admitted this exploration". A skipped exploration was
            # not admitted by anything, and only the qualification gate can decline one -- a caller
            # that admitted a network did not skip it, and a run with no selection had nothing to
            # decline it. So any other value here is a false statement about why nothing ran, not
            # merely an unusual combination.
            raise ValueError(
                f"PDepExplorationResult with status='skipped' must carry "
                f"admission_policy={ADMISSION_POLICY_QUALIFIED_SELECTION!r}, got "
                f"{self.admission_policy!r}: the qualification gate is the only thing that can "
                f"decline an exploration, so it is the only thing a skipped result can name.")
        if self.status == EXPLORATION_STATUS_SKIPPED:
            if self.network_paths or self.output_paths or self.k_tp_as_written:
                raise ValueError(
                    "PDepExplorationResult with status='skipped' cannot carry network_paths / "
                    "output_paths / k_tp_as_written: a skipped exploration never ran, so it cannot "
                    "have produced artifacts.")
        if self.status in (EXPLORATION_STATUS_FAILED, EXPLORATION_STATUS_SKIPPED):
            if not self.reasons:
                raise ValueError(
                    f"PDepExplorationResult with status={self.status!r} must carry at least one "
                    f"entry in reasons: a negative outcome with no stated reason is unusable to a "
                    f"caller deciding what to do next.")
        if self.status == EXPLORATION_STATUS_SUCCEEDED and not self.network_paths:
            raise ValueError(
                "PDepExplorationResult with status='succeeded' must carry at least one entry in "
                "network_paths: success with no artifact is a fail-open shape.")
        if self.status == EXPLORATION_STATUS_FAILED and (self.network_paths or self.k_tp_as_written):
            raise ValueError(
                "PDepExplorationResult with status='failed' cannot carry network_paths or "
                "k_tp_as_written: those assert a usable result exists, which is the claim the "
                "failure denies. A failed run MAY carry output_paths -- Arkane's own logs and "
                "partial output are exactly what a human diagnosing the failure needs.")
        # Copied, not aliased. A frozen dataclass holding someone else's live dict/object is frozen
        # in name only: whoever passed it in can keep editing the provenance of a run that has
        # already been reported. ``selection`` is a mutable ``PDepNetworkSelection`` (callers
        # routinely do ``selection.warnings.append(...)`` / ``selection.evaluation_status = ...``),
        # so it needs the same defense as ``manifest`` -- a caller still holding it could otherwise
        # rewrite the recorded decision after the fact. Callers that want to correlate this result
        # with the selection object they passed in must therefore do so by equality, not identity
        # (``PDepNetworkSelection`` is a plain dataclass, so ``==`` compares field values).
        object.__setattr__(self, 'manifest', copy.deepcopy(self.manifest))
        object.__setattr__(self, 'selection', copy.deepcopy(self.selection))

    def as_dict(self) -> dict:
        """
        Render this result as plain JSON/YAML-safe types.

        Returns:
            dict: A plain dict, containing no dataclass instances or tuples anywhere. ``selection``
                is rendered via ``PDepNetworkSelection.as_dict()`` (or ``None``), ``k_tp_as_written``
                via ``PDepArkaneReaction.as_dict()``, and ``manifest`` is deep-copied with any nested
                tuples converted to lists (a manifest is free-form provenance data, so it is not
                assumed to already be shallow/flat).
        """
        return {
            'network_id': self.network_id,
            'status': self.status,
            'reasons': list(self.reasons),
            'network_paths': list(self.network_paths),
            'output_paths': list(self.output_paths),
            'k_tp_as_written': [entry.as_dict() for entry in self.k_tp_as_written],
            'manifest': to_json_safe(self.manifest),
            'selection': self.selection.as_dict() if self.selection is not None else None,
            'admission_policy': self.admission_policy,
        }
