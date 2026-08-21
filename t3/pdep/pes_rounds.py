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

import logging
import os
from dataclasses import dataclass, replace

from arc.molecule.molecule import Molecule

from t3.pdep.barrierless import classify_barrierless, rmg_family
from t3.pdep.parser import PDepNetwork, PDepPathReaction, canonical_channel_pair

_logger = logging.getLogger(__name__)

# adjacency-list text -> canonical SMILES (or None for unparseable text), memoized because one
# loop converts the same handful of structures once per round per consumer.
_CANONICAL_STRUCTURE_CACHE: dict = dict()


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


def _canonical_structure(adjlist: str) -> str | None:
    """
    Reduce one RMG adjacency list to a canonical, label- and atom-order-independent identity.

    ARC's own ``Molecule`` (never ``rmgpy``) parses the adjacency list and renders a canonical
    SMILES, so two projects that wrote the same molecule with different species labels or a
    different atom order still produce the same identity string.

    What the identity is NOT independent of -- deliberately: the spin multiplicity is appended
    explicitly (``...|m<multiplicity>``), because SMILES alone does not encode spin state and two
    genuinely distinct RMG species collapse onto one string without it -- singlet and triplet CH2
    are both ``'[CH2]'``. Adopting a singlet channel's barrier onto the triplet channel (or vice
    versa) is exactly the wrong-saddle-point misattribution structural keying exists to prevent,
    so the multiplicity is part of the identity. Charge and lone-pair differences already separate
    through the SMILES itself (``[OH]`` vs ``[OH-]``; O(3P) ``[O]`` vs O(1D) ``O``).

    Args:
        adjlist (str): The RMG adjacency-list text.

    Returns:
        str | None: The canonical identity (``<canonical SMILES>|m<multiplicity>``), or ``None``
        if the text could not be parsed into a molecule -- the caller must then refuse to key the
        channel rather than guess.
    """
    if adjlist in _CANONICAL_STRUCTURE_CACHE:
        return _CANONICAL_STRUCTURE_CACHE[adjlist]
    try:
        molecule = Molecule().from_adjacency_list(adjlist)
        identity = f'{molecule.to_smiles()}|m{molecule.multiplicity}'
    except Exception as e:
        _logger.warning(f'Could not canonicalize an adjacency list into a molecule ({e}); the '
                        f'channel using it cannot be structurally keyed.')
        identity = None
    _CANONICAL_STRUCTURE_CACHE[adjlist] = identity
    return identity


def structural_channel_key(path_reaction: PDepPathReaction,
                           species_structures: dict) -> tuple | None:
    """
    The direction-insensitive STRUCTURAL identity of the channel a path reaction connects.

    Why labels are not enough: every TS label Arkane writes is purely positional --
    ``rmgpy/rmg/pdep.py:854`` replaces every path reaction's transition state with a fresh,
    label-less object before every network file write, so ``TS3`` means nothing more than
    "whatever sat at ``path_reactions[2]`` when that file was written", and pruning or discovery
    between explorations shifts indices. Reaction labels (``reaction<i>``) are positional too,
    and can even collide within one file after a remove-then-append. Carrying anything across
    network-file writes by either label can attach a computed barrier to the WRONG saddle point.
    This key is what may be carried instead: the two channels the reaction connects, each species
    reduced to its canonical structure (``_canonical_structure``), canonicalized
    direction-insensitively through ``t3.pdep.parser.canonical_channel_pair`` -- the identity the
    network's own topology assigns the reaction, immune to renumbering AND to per-project species
    label indices.

    Args:
        path_reaction (PDepPathReaction): The path reaction to key.
        species_structures (dict): Species label -> RMG adjacency-list text, i.e. the owning
            network's ``PDepNetwork.species_structures``.

    Returns:
        tuple | None: The structural key, or ``None`` when any participating species has no
        parseable structure in ``species_structures`` -- the channel then has no structural
        identity and must not be carried (fail closed, never keyed by label instead).
    """
    sides = list()
    for labels in (path_reaction.reactants, path_reaction.products):
        side = list()
        for label in labels:
            adjlist = species_structures.get(label)
            if not adjlist:
                return None
            smiles = _canonical_structure(adjlist)
            if smiles is None:
                return None
            side.append(smiles)
        sides.append(side)
    return canonical_channel_pair(sides[0], sides[1])


def channel_keys_by_ts_label(network: PDepNetwork) -> dict:
    """
    Map each of a network's transition state labels to its channel's structural key.

    Fail-closed by construction, in three ways, because this map is what QM artifacts are carried
    across network-file writes on (see ``structural_channel_key`` for why labels cannot be):

    - a transition state shared by SEVERAL path reactions is omitted (which of its channels would
      be "the" identity is ambiguous -- the same reason ``t3.main`` refuses to queue such a TS);
    - a transition state whose channel cannot be structurally keyed (a species with no parseable
      structure) is omitted;
    - two transition states whose channels key IDENTICALLY are BOTH omitted -- a duplicate channel
      pair means the structural key cannot distinguish them, and attaching an artifact to
      "whichever matched first" is exactly the wrong-saddle-point failure this keying exists to
      prevent.

    An omitted transition state is simply re-decided (and, if selected, re-queued) each round
    rather than carried -- duplicated QM spend, never a misattributed barrier.

    Args:
        network (PDepNetwork): The parsed network.

    Returns:
        dict: Network-local TS label -> structural key, for every transition state that has
        exactly one, unambiguous structural identity.
    """
    keys = dict()
    for ts_label, path_reactions in network.path_reactions_by_ts().items():
        if len(path_reactions) != 1:
            _logger.warning(f'Transition state {ts_label!r} of network {network.network_id!r} is '
                            f'shared by {len(path_reactions)} path reactions; it has no single '
                            f'structural channel identity and will not be carried across rounds.')
            continue
        key = structural_channel_key(path_reactions[0], network.species_structures)
        if key is None:
            _logger.warning(f'Transition state {ts_label!r} of network {network.network_id!r} '
                            f'cannot be structurally keyed (a participating species has no '
                            f'parseable structure); it will not be carried across rounds.')
            continue
        keys[ts_label] = key
    by_key = dict()
    for ts_label, key in keys.items():
        by_key.setdefault(key, list()).append(ts_label)
    for key, ts_labels in by_key.items():
        if len(ts_labels) > 1:
            _logger.warning(f'Transition states {sorted(ts_labels)} of network '
                            f'{network.network_id!r} connect structurally identical channels; '
                            f'none of them can be carried across rounds unambiguously.')
            for ts_label in ts_labels:
                del keys[ts_label]
    return keys


def adoption_channel_keys_by_ts_label(network: PDepNetwork) -> dict:
    """
    Map each transition state to the identity CROSS-CAPTURE adoption may be matched on: its
    structural channel key PLUS a discriminator for the PATHWAY that produced it -- the RMG family
    when the kinetics comment names one, the comment itself otherwise.

    Why the endpoints alone are not enough here. ``structural_channel_key`` identifies the
    reactants and products a channel connects -- not a saddle point. ``channel_keys_by_ts_label``
    covers the case where two pathways between the same endpoints sit in ONE file (both are
    dropped as structurally duplicated), but that guard cannot reach across files: if a prior
    capture's network holds one A<=>B pathway and this run's network holds a DIFFERENT A<=>B
    pathway, each file holds a single, identical key, no duplicate is present, and the prior
    artifact is attached to a saddle point it was never computed for -- silently, as a computed
    barrier. That is the precise misattribution this whole keying exists to prevent.

    What actually discriminates a pathway, established empirically against the vendored real
    networks rather than assumed:

    - The RMG family (``t3.pdep.barrierless.rmg_family``) is the obvious candidate, and on its own
      it is not enough: it is present in only 1 of 3 path reactions in
      ``tests/data/pdep_real_networks/network799_1`` and 0 of 3 in ``network21_1``, whose kinetics
      come from a reaction library and name no family at all.
    - The kinetics comment the family is parsed OUT of is present in all six, and discriminates
      just as well -- ``'Estimated from node Root_N-1R!H->C'`` and
      ``"Reaction library: 'primaryNitrogenLibrary'"`` each name the specific rate rule or library
      a pathway came from. It is the fallback here.
    - Nothing else exists. The kinetics comment is the ONLY carrier of pathway provenance in a
      network file, and ``t3.pdep.hybrid`` DELIBERATELY drops the whole ``kinetics = ...`` entry
      of every QM'd path reaction (its invariant 2, ``t3/pdep/hybrid.py:18``), so a comment-derived
      discriminator is null on exactly the channels a WITHIN-RUN carry would need it for. The
      transition state's ``E0`` goes the same way, replaced by a statmech ``Log(...)`` reference.

    So the discriminator exists for cross-capture adoption -- where both sides are RMG-estimated
    networks that still carry their kinetics comments -- and a channel that carries neither is
    REFUSED (omitted from this map, logged) rather than matched on its endpoints. A refused match
    costs quantum chemistry that was already paid for once; a false match costs correctness of the
    barrier the whole loop exists to compute, presented as if it had been computed. That asymmetry
    decides it.

    This is NOT what the within-run, round-to-round carry uses -- see ``channel_keys_by_ts_label``,
    which stays endpoints-only precisely because no discriminator survives the hybrid write this
    loop itself performs between one round and the next.

    Args:
        network (PDepNetwork): The parsed network.

    Returns:
        dict: Network-local TS label -> ``(structural key, pathway discriminator)``, for every
        transition state that has both an unambiguous structural identity and a discriminator.
    """
    path_reactions_by_ts = network.path_reactions_by_ts()
    keys = dict()
    for ts_label, key in channel_keys_by_ts_label(network).items():
        path_reaction = path_reactions_by_ts[ts_label][0]
        # Family first, comment second. The family is the coarser and more stable of the two: it
        # survives rate-rule retraining and degeneracy changes that reword a comment, so preferring
        # it keeps a family-bearing channel adoptable across RMG database versions. The comment is
        # the fallback, not a second-class one -- 'Estimated from node Root_N-1R!H->C' and
        # "Reaction library: 'primaryNitrogenLibrary'" each name the specific rate rule or library
        # a pathway came from, which is exactly the pathway identity the endpoints cannot supply.
        discriminator = rmg_family(path_reaction) or (path_reaction.kinetics_comment or '').strip()
        if not discriminator:
            _logger.info(
                f'Transition state {ts_label!r} of network {network.network_id!r} names neither an '
                f'RMG family nor any kinetics comment, so its channel is identified by its '
                f'endpoints alone. A different pathway between the same endpoints would key '
                f'identically across captures, so prior QM will not be adopted for it -- it will '
                f'be recomputed instead.')
            continue
        keys[ts_label] = (key, discriminator)
    return keys


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
