"""
t3 pdep explorer result module.

Defines ``PDepExplorationResult``, the plain-data record a PES exploration decision (skip / run-
and-fail / run-and-succeed) is packaged into for callers outside the adapter layer. It mirrors the
house style set by ``t3.pdep.selector.PDepNetworkSelection``: a frozen dataclass with a status
string, evidence fields, and an ``as_dict()`` that renders everything as plain JSON/YAML-safe
types.
"""

from dataclasses import dataclass, field

from t3.pdep.parser import PDepArkaneReaction, to_json_safe
from t3.pdep.selector import PDepNetworkSelection

# Exploration status values for PDepExplorationResult.status. THREE statuses, not a bool, mirroring
# the never-ran-vs-failed distinction ArkaneExplorerAdapter already preserves with its
# ``succeeded is None`` sentinel (t3/pdep/explorer/arkane.py, ~line 147):
#
# - 'skipped'   the exploration was never attempted, because the caller's budget gate (a
#               non-qualifying PDepNetworkSelection) declined it. No Arkane job ran at all.
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
            attempted because the caller's budget gate declined it; ``'failed'`` means it ran and did
            not succeed; ``'succeeded'`` means it ran and succeeded.
        reasons (tuple): Human-readable reasons for a ``'failed'`` or ``'skipped'`` outcome. Required
            (non-empty) for both: a negative outcome with no stated reason is unusable to a caller
            trying to decide what to do next.
        network_paths (tuple): The resolved final network file path(s) an exploration wrote, e.g.
            from ``ArkaneExplorerAdapter.get_networks()``. Must be empty unless ``status`` is
            ``'succeeded'``, and non-empty when it is (success with no artifact is a fail-open shape).
        output_paths (tuple): The raw Arkane ``output*.py`` path(s) an exploration wrote. Must be
            empty unless ``status`` is ``'succeeded'``.
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
        selection (PDepNetworkSelection, optional): The budget-gate decision this exploration result
            is downstream of, if any (e.g. the non-qualifying decision that produced a ``'skipped'``
            result).

    Raises:
        ValueError: If ``status`` is not one of the ``EXPLORATION_STATUS_*`` constants; if
            ``status`` is ``'skipped'`` and any of ``network_paths``/``output_paths``/
            ``k_tp_as_written`` is non-empty (nothing ran, so there can be no artifacts); if
            ``status`` is ``'failed'`` or ``'skipped'`` and ``reasons`` is empty; if ``status`` is
            ``'succeeded'`` and ``network_paths`` is empty.
    """
    network_id: str | None
    status: str
    reasons: tuple = tuple()
    network_paths: tuple = tuple()
    output_paths: tuple = tuple()
    k_tp_as_written: tuple = tuple()
    manifest: dict = field(default_factory=dict)
    selection: 'PDepNetworkSelection | None' = None

    def __post_init__(self):
        if self.status not in _VALID_EXPLORATION_STATUSES:
            raise ValueError(
                f'PDepExplorationResult.status must be one of {_VALID_EXPLORATION_STATUSES}, '
                f'got {self.status!r}.')
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
        }
