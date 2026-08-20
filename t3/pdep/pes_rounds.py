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

import os
from dataclasses import dataclass, replace

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
        coefficient (float | None): The signed dln(k)/dE0 sensitivity coefficient (mol/J) this
            round's master-equation SA measured for this transition state
            (``t3.pdep.pes_sa.run_round_me_sensitivity``), stamped by
            ``attach_sensitivity_evidence``. ``None`` only before stamping;
            ``t3.pdep.pes_qm.arc_qm_runner`` refuses to queue a candidate left at ``None``,
            because ``t3.pdep.capture`` would (rightly) refuse its artifact after the QM spend.
        delta_ln_k (float | None): The corresponding dimensionless rate response,
            ``abs(coefficient) * perturbation`` -- same convention as
            ``t3.pdep.selector.SensitiveTransitionState``.
    """
    path_reaction: PDepPathReaction
    ts_label: str
    family: str | None
    coefficient: float | None = None
    delta_ln_k: float | None = None


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


def attach_sensitivity_evidence(split: CandidateSplit,
                                evidence_by_ts_label: dict,
                                min_delta_ln_k: float = 0.0) -> CandidateSplit:
    """
    Stamp each candidate with the real sensitivity evidence this round's ME SA measured for it.

    A candidate whose transition state has NO entry in ``evidence_by_ts_label`` is moved to
    ``skipped`` with its reason recorded, never queued: ``t3.pdep.capture`` refuses (fail-closed,
    by design) any captured artifact that carries no ``coefficient``/``delta_ln_k`` evidence, so
    queueing such a candidate would spend real QM and then lose the result at capture time -- and
    inventing a number instead would violate the standing constraint that no rate-determining
    parameter may be fabricated.

    A candidate whose measured ``delta_ln_k`` is BELOW ``min_delta_ln_k`` is also skipped (with
    the measured value in its reason): a response below the floor -- including a measured exact
    zero -- does not justify the QM spend, and folding its coefficient into a capture manifest as
    "the sensitivity evidence that justified selecting this transition state" would make that
    sentence false. This mirrors T3's in-run selection, where ``t3.pdep.selector``'s
    ``_bounded_cutoff`` floors at ``min_delta_ln_k / perturbation > 0`` and such a record is
    definitionally unreachable. The floor applies under BOTH ``qm.scope`` values -- 'sensitive'
    ranks and 'all' does not, but neither may queue below it.

    Args:
        split (CandidateSplit): The split to stamp, from ``split_qm_candidates``.
        evidence_by_ts_label (dict): Network-local TS label -> ``(coefficient, delta_ln_k)``,
            from ``t3.pdep.pes_sa.run_round_me_sensitivity``.
        min_delta_ln_k (float, optional): The smallest measured ln(k) response that justifies
            queueing (``config.qm.min_delta_ln_k``). ``0.0`` disables the floor.

    Returns:
        CandidateSplit: The same candidates in the same order, each carrying its evidence, minus
        the ones with no evidence or with a response below the floor (appended to ``skipped``,
        each with its reason).
    """
    candidates, skipped = [], list(split.skipped)
    for candidate in split.candidates:
        pair = evidence_by_ts_label.get(candidate.ts_label)
        if pair is None:
            skipped.append(SkippedChannel(
                label=candidate.path_reaction.label,
                reason=f"'{candidate.path_reaction.label}': transition state "
                       f'{candidate.ts_label} has no finite sensitivity evidence in this '
                       f"round's master-equation sensitivity analysis, and a captured artifact "
                       f'must carry the evidence that justified selecting it '
                       f'(t3.pdep.capture); not queueing it rather than inventing a number.'))
            continue
        coefficient, delta_ln_k = pair
        if delta_ln_k < min_delta_ln_k:
            skipped.append(SkippedChannel(
                label=candidate.path_reaction.label,
                reason=f"'{candidate.path_reaction.label}': transition state "
                       f'{candidate.ts_label} measured a ln(k) response of {delta_ln_k:.3e}, '
                       f'below the min_delta_ln_k floor ({min_delta_ln_k:.3e}); its leverage on '
                       f'the network is too small to justify the QM spend.'))
            continue
        candidates.append(replace(candidate, coefficient=coefficient, delta_ln_k=delta_ln_k))
    return CandidateSplit(candidates=tuple(candidates), skipped=tuple(skipped))


PES_LOOP_DIAGRAM_FILENAME = 'pes_diagram.png'


@dataclass(frozen=True)
class RoundPaths:
    """
    Where one round of the loop puts its artifacts.

    Attributes:
        root (str): The round's own directory.
        arc_project (str): The ARC project directory for this round's QM.
        explorer_output (str): Where the Arkane explorer writes.
        sa (str): Where this round's master-equation sensitivity analysis runs
            (``t3.pdep.pes_sa.run_round_me_sensitivity``).
        capture (str): Where this round's QM artifacts are frozen.
        hybrid (str): Where this round's hybrid network input file is written.
        diagram (str): The PES diagram for this round.
    """
    root: str
    arc_project: str
    explorer_output: str
    sa: str
    capture: str
    hybrid: str
    diagram: str


def round_paths(project_directory: str, round_index: int) -> RoundPaths:
    """
    Resolve the artifact layout for one round.

    Every round is self-contained, and in particular gets its OWN ARC project directory. That is
    what lets the loop run ARC more than once without fighting ``t3.pdep.capture``'s single-shot
    window: ARC recreates ``calcs/statmech/`` with ``delete_existing_subdir=True`` on every rate
    pass, so a second ARC run sharing one project directory would destroy the first round's
    un-captured artifacts. Separate projects make that structurally impossible rather than
    merely discouraged.

    The capture directory is deliberately a sibling of the ARC project, never nested inside it:
    ``capture_ts_artifacts`` refuses a ``capture_dir`` that resolves inside the ARC project
    directory, for the same reason.

    Args:
        project_directory (str): The loop's project directory. Must be absolute.
        round_index (int): The zero-based round number.

    Returns:
        RoundPaths: The resolved layout.

    Raises:
        ValueError: If ``project_directory`` is not absolute, or ``round_index`` is negative.
    """
    if not os.path.isabs(project_directory):
        raise ValueError(f"'project_directory' must be an absolute path, got "
                         f"{project_directory!r}.")
    if round_index < 0:
        raise ValueError(f"'round_index' must be non-negative, got {round_index}.")
    root = os.path.join(project_directory, f'round_{round_index}')
    return RoundPaths(root=root,
                      arc_project=os.path.join(root, 'ARC'),
                      explorer_output=os.path.join(root, 'explorer'),
                      sa=os.path.join(root, 'SA'),
                      capture=os.path.join(root, 'capture'),
                      hybrid=os.path.join(root, 'hybrid'),
                      diagram=os.path.join(root, PES_LOOP_DIAGRAM_FILENAME))


def hybrid_network_path(paths: RoundPaths, network_id: str) -> str:
    """
    Where a round's ``qm_runner`` must write its hybrid network input file.

    ``RoundPaths.hybrid`` is a directory, not a file, and the PES loop needs a file path to hand
    the next round's explorer. It also needs that file's stem to carry ``network_id`` rather than
    the literal ``'hybrid'``, because ``parse_pdep_network_file`` derives ``network_id =
    Path(path).stem`` (``t3/pdep/parser.py:729``) -- every round writing to a ``hybrid.py`` stem
    would collapse distinct networks onto one ``network_id``, and with it one ARC job-label
    namespace (the exact failure ruling C4 exists to prevent).

    This lives in ``t3.pdep.pes_rounds`` (not ``t3.pdep.pes_loop``, which re-exports it for
    backward compatibility) so that both ``t3.pdep.pes_loop`` and ``t3.pdep.pes_qm`` can import it
    without creating an import cycle between those two modules.

    Args:
        paths (RoundPaths): This round's paths.
        network_id (str): The network id to preserve in the file's stem (normally the just-explored
            network's own ``network_id``).

    Returns:
        str: The hybrid network file path this round's ``qm_runner`` must write to, and the path
        the next round explores from.
    """
    return os.path.join(paths.hybrid, f'{network_id}.py')
