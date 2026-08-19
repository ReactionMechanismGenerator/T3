"""
t3 pdep pes_rounds module

The per-round bookkeeping of the PES exploration loop: which path reactions in an explored network
still need quantum chemistry, and which must be left alone.

Three separate reasons disqualify a path reaction from QM, and they are kept distinct because they
mean different things to someone reading the round record:

* it is **barrierless** -- there is no saddle point to find, so QM cannot succeed (see
  ``t3.pdep.barrierless``);
* it **already has QM** from an earlier round or an adopted prior project -- re-queueing it would
  spend the budget twice for the same number;
* it **declares no transition state** in the network file -- there is nothing to point a job at.

Every disqualification is recorded with its reason rather than filtered out silently. A loop that
quietly drops channels produces a PES that looks complete and is not, which is precisely the
failure mode the whole exercise exists to remove.

Candidate order is the network file's own order, deterministically. The budget admits a prefix of
this list, so a non-deterministic order would admit a different subset on each run and make two
runs of the same input incomparable.
"""

from dataclasses import dataclass

from t3.pdep.barrierless import classify_barrierless
from t3.pdep.parser import PDepNetwork, PDepPathReaction


@dataclass(frozen=True)
class QMCandidate:
    """
    A path reaction whose transition state should be sent to ARC.

    Attributes:
        path_reaction (PDepPathReaction): The parsed path reaction.
        ts_label (str): The network-local transition state label (e.g. ``'TS1'``).
        family (str | None): The RMG family, when the kinetics comment named one.
    """
    path_reaction: PDepPathReaction
    ts_label: str
    family: str | None


@dataclass(frozen=True)
class SkippedChannel:
    """
    A path reaction deliberately not sent to ARC, and why.

    Attributes:
        label (str): The path reaction's label.
        reason (str): Why it was skipped, for the log and the round record.
    """
    label: str
    reason: str


@dataclass(frozen=True)
class CandidateSplit:
    """
    The outcome of splitting a network's path reactions into QM candidates and skips.

    Attributes:
        candidates (tuple): The ``QMCandidate`` objects, in network-file order.
        skipped (tuple): The ``SkippedChannel`` objects, in network-file order.
    """
    candidates: tuple = ()
    skipped: tuple = ()


def split_qm_candidates(network: PDepNetwork, computed_ts_labels: frozenset) -> CandidateSplit:
    """
    Split a network's path reactions into those needing QM and those that must not get it.

    Args:
        network (PDepNetwork): The parsed network, freshly explored or freshly generated.
        computed_ts_labels (frozenset): Network-local TS labels that already have QM, from this
                                        loop's earlier rounds or from an adopted prior project.

    Returns:
        CandidateSplit: Candidates and skips, each in the network file's own order.
    """
    candidates, skipped = [], []
    for path_reaction in network.path_reactions:
        ts_label = path_reaction.transition_state
        if ts_label is None:
            skipped.append(SkippedChannel(
                label=path_reaction.label,
                reason=f"'{path_reaction.label}': declares no transition state in the network "
                       f'file, so there is nothing to compute.'))
            continue
        if ts_label in computed_ts_labels:
            skipped.append(SkippedChannel(
                label=path_reaction.label,
                reason=f"'{path_reaction.label}': transition state {ts_label} already has QM; not "
                       f'queueing it again.'))
            continue
        verdict = classify_barrierless(path_reaction)
        if verdict.is_barrierless:
            skipped.append(SkippedChannel(label=path_reaction.label, reason=verdict.reason))
            continue
        candidates.append(QMCandidate(path_reaction=path_reaction, ts_label=ts_label,
                                      family=verdict.family))
    return CandidateSplit(candidates=tuple(candidates), skipped=tuple(skipped))
