"""
t3 pdep pes_loop module

The standalone PES exploration loop: explore a network with Arkane, draw the PES diagram Arkane's
own output supports, queue the barriered channels that still lack a transition-state QM result,
run QM, fold the results back into the network, and re-explore -- round after round until the
network converges, the round budget runs out, or a round cannot explore at all.

Why this lives here, standalone, rather than inside ``t3.main``'s own iteration loop: the whole
point is to let this loop be run and corroborated on its own, against a single network, before it
is trusted to influence what T3's kinetic-model iteration spends its QM budget on. ``t3.main``
folding directly into an unreviewed loop that sequences real ARC submissions would make every bug
here a live-campaign bug; a standalone entry point lets Task 6's real ARC-backed ``qm_runner`` and
this driver be exercised and debugged in isolation first.

Why the PES diagram is redrawn from EACH round's freshly explored network, never from the input
network some earlier round started with: T3's hybrid network file (``t3.pdep.hybrid``) represents
a QM'd transition state as a statmech artifact reference -- a ``Log(...)`` chain into a Gaussian
output file -- never as an inline ``E0`` (see ``t3/pdep/hybrid.py:46``). The computed barrier does
not exist as a NUMBER until Arkane actually runs and applies its atom-energy corrections to that
statmech data. Drawing from any input file, hybrid or otherwise, can therefore never show a
computed barrier; only Arkane's own output network can. This is the bug this whole loop exists to
fix: T3 has been observed to draw its PES diagram once, from the pre-QM network, and never redraw
it, producing byte-identical diagrams across two runs even though one of them converged a genuine
saddle point.

Why ``qm_runner`` is injected rather than called directly: this loop's own job is sequencing --
explore, draw, split candidates from skips, decide whether to stop -- and none of that requires an
ARC cluster to test. Injecting the QM step as a callable with the signature
``qm_runner(candidates, paths, config, network_id) -> tuple[frozenset[str], frozenset[str]]``
(``converged_ts_labels, queued_ts_labels``) lets every round-sequencing decision (no-candidates,
converged, stalled, diagram-only, max-rounds, failed-explore) be covered
by tests that never touch ARC. Task 6 supplies the real ARC-backed implementation; ``qm_runner=None``
is not a placeholder for "not implemented yet", it is the honest, permanent behaviour for a caller
who wants explore-and-draw only, with no QM spend at all -- so it must never raise, and it is
recorded as ``PES_LOOP_DIAGRAM_ONLY``, never ``PES_LOOP_STALLED``: the two look similar (a round
that computed no new TS) but differ in kind -- ``stalled`` is an operational fault worth
investigating (ARC ran and made no progress), ``diagram_only`` is complete success of exactly what
was asked (no QM was ever requested).

Why every round gets its own ARC project directory (``t3.pdep.pes_rounds.round_paths``): ARC
recreates its ``calcs/statmech/`` directory with ``delete_existing_subdir=True`` on every rate
pass, so a second ARC run sharing one project directory would destroy the first round's
un-captured artifacts. Giving each round a fresh, self-contained directory tree makes that
structurally impossible instead of merely a discipline to remember.

This module never imports ``rmgpy`` or ``arkane``, and never ``exec``s or ``import``s a network
file: everything at that boundary flows through ``t3.pdep.parser`` (``ast.parse`` only) and
``t3.pdep.diagram``.
"""

import logging
import os
from dataclasses import dataclass
from pathlib import Path

from t3.pdep.api import explore_pdep_network
from t3.pdep.capture import captured_qm_artifact_path
from t3.pdep.diagram import draw_pes_diagram
from t3.pdep.join import arc_ts_label
from t3.pdep.explorer.config import PDepExplorerConfig
from t3.pdep.explorer.result import EXPLORATION_STATUS_SUCCEEDED
from t3.pdep.parser import parse_pdep_network_file
from t3.pdep.pes_qm import adopt_prior_qm
from t3.pdep.pes_rounds import (RoundPaths, adoption_channel_keys_by_ts_label,
                                attach_sensitivity_evidence,
                                channel_keys_by_ts_label, hybrid_network_path, round_paths,
                                split_qm_candidates)
from t3.pdep.pes_sa import run_round_me_sensitivity
from t3.schema import PESLoopConfig

_logger = logging.getLogger(__name__)


# Status values for PESLoopResult.status / RoundRecord.status.
#
# - 'converged'      every barriered channel the network declares now has QM, or every channel
#                     still lacking it measured a ln(k) response below qm.min_delta_ln_k -- a
#                     measured "not worth the spend", recorded per channel in the round's skipped
#                     entries (round > 0).
# - 'no_candidates'  the very first round already had nothing to compute (round == 0): a network
#                     that is barrierless end to end, or whose every barriered channel measured a
#                     response below the floor, is not a loop that "converged" after doing
#                     work, it never had any to do.
# - 'stalled'        a round ran the QM runner and it returned no newly-converged TS labels, and
#                     ``config.termination.stop_when_no_new_ts`` says that is a reason to stop now
#                     rather than spend the rest of the round budget on a runner that is not making
#                     progress. This is an OPERATIONAL FAULT: ARC ran, spent time, and produced
#                     nothing usable, which is worth investigating.
# - 'diagram_only'    ``qm_runner=None``: no QM was ever going to be attempted, so the loop explored
#                     and drew one diagram and stopped, exactly as asked. This is the COMPLETE
#                     SUCCESS of a diagram-only request, not a fault, and must be distinguishable
#                     from 'stalled' so that a monitor alerting on 'stalled' does not page on every
#                     diagram-only run.
# - 'max_rounds'      the round budget was exhausted without converging, without stalling, and
#                     without a failed explore.
# - 'failed'          a round could not explore its network at all -- including a round whose
#                     expected hybrid network input file was never written by ``qm_runner`` -- never
#                     papered over as convergence.
PES_LOOP_CONVERGED = 'converged'
PES_LOOP_MAX_ROUNDS = 'max_rounds'
PES_LOOP_NO_CANDIDATES = 'no_candidates'
PES_LOOP_STALLED = 'stalled'
PES_LOOP_DIAGRAM_ONLY = 'diagram_only'
PES_LOOP_FAILED = 'failed'


@dataclass(frozen=True)
class RoundRecord:
    """
    What one round of the PES exploration loop did.

    Attributes:
        index (int): The zero-based round number.
        network_path (str): The network file this round explored from (before this round's own
            re-exploration), or the freshly explored network if exploration succeeded.
        diagram_path (str | None): Where this round's PES diagram was written, or ``None`` if
            drawing it was refused (a species with no ``E0``) or exploration itself failed.
        queued_ts_labels (tuple): The network-local TS labels ``qm_runner`` actually queued this
            round. For a real ``qm_runner`` this is its own reported ``queued_ts_labels`` (N3: it
            can silently drop a candidate the loop offered, e.g. for missing structure, so this is
            not simply every candidate the loop handed it); for ``qm_runner=None`` or a round that
            never reaches a runner, it is every candidate the loop would have offered.
        skipped (tuple): The ``SkippedChannel`` entries for this round, each carrying its reason.
        status (str): This round's outcome. ``'continuing'`` for a round that queued QM and is not
            the loop's last one; otherwise one of the ``PES_LOOP_*`` status constants, matching the
            ``PESLoopResult.status`` this round ended the loop with.
        reason (str): Why, for a non-obvious status (chiefly ``'failed'``). ``''`` otherwise.
    """
    index: int
    network_path: str
    diagram_path: str | None
    queued_ts_labels: tuple
    skipped: tuple
    status: str
    reason: str = ''


@dataclass(frozen=True)
class PESLoopResult:
    """
    The outcome of running the PES exploration loop to a stop.

    Attributes:
        rounds (tuple): The ``RoundRecord`` for every round that ran, in order.
        status (str): One of the ``PES_LOOP_*`` status constants.
        reason (str): Why the loop stopped, when ``status`` alone does not say (chiefly
            ``'failed'``). ``''`` otherwise.
        final_network_path (str | None): The last successfully explored network file, or ``None``
            if the very first round failed to explore. On ``'max_rounds'`` with the last round
            having converged (or adopted) anything, this is the re-exploration of THAT round's
            hybrid network -- the surface carrying the barriers that round computed -- not the
            pre-QM network its own ``RoundRecord`` carries.
        final_diagram_path (str | None): The PES diagram drawn from ``final_network_path``, or
            ``None`` if no diagram could be drawn. Same rule as above: on ``'max_rounds'`` after a
            productive final round this is the QM-informed diagram, drawn after that round's
            quantum chemistry was folded in, never the last ``RoundRecord``'s pre-QM one.
    """
    rounds: tuple
    status: str
    reason: str
    final_network_path: str | None
    final_diagram_path: str | None


def _build_explorer_config(config: PESLoopConfig, project_directory: str,
                           paths: RoundPaths) -> PDepExplorerConfig:
    """Build this round's ``PDepExplorerConfig`` from ``config.pes``."""
    return PDepExplorerConfig(
        explorer='Arkane',
        trusted_output_root=project_directory,
        output_directory=paths.explorer_output,
        seed_species=tuple(config.pes.source),
        method=config.pes.method,
        bath_gas=config.pes.bath_gas,
        explore_tol=config.pes.explore_tol,
        energy_tol=config.pes.energy_tol,
        flux_tol=config.pes.flux_tol,
        maximum_radical_electrons=config.pes.maximum_radical_electrons,
        timeout=config.pes.timeout,
    )


def _trim_candidates(candidates: tuple, config: PESLoopConfig, logger) -> tuple:
    """
    Apply ``config.qm.scope`` and ``config.qm.max_transition_states_per_round`` to ``candidates``.

    Every candidate reaching this function already carries the REAL sensitivity evidence this
    round's own master-equation SA measured for it (``t3.pdep.pes_sa.run_round_me_sensitivity``,
    stamped by ``t3.pdep.pes_rounds.attach_sensitivity_evidence``), so ``scope == 'sensitive'``
    ranks by that evidence -- descending ``abs(coefficient)``, ties keeping network-file order --
    rather than degrading to file order the way it had to before the loop had an SA of its own.
    ``scope == 'all'`` keeps the network file's own order. Both take the first ``limit``.

    Args:
        candidates (tuple): The ``QMCandidate`` entries to trim, in network-file order, each
                            already stamped with its sensitivity evidence.
        config (PESLoopConfig): The loop's configuration.
        logger: A T3 ``Logger``, or ``None``.

    Returns:
        tuple: The trimmed candidates -- evidence-ranked for ``scope == 'sensitive'``,
        network-file order for ``scope == 'all'``.
    """
    limit = config.qm.max_transition_states_per_round
    if config.qm.scope == 'sensitive':
        # sorted() is stable, so equal |coefficient| keeps network-file order -- the determinism
        # split_qm_candidates' own docstring promises.
        candidates = tuple(sorted(candidates, key=lambda candidate: -abs(candidate.coefficient)))
    if len(candidates) > limit:
        message = (f"PES loop: taking the {'most sensitive' if config.qm.scope == 'sensitive' else 'first'} "
                   f'{limit} of {len(candidates)} candidate(s) '
                   f'(qm.max_transition_states_per_round).')
        _logger.info(message)
        if logger is not None:
            logger.info(message)
    return candidates[:limit]


def _draw_round_diagram(explored_network_path: str, diagram_path: str, logger) -> str | None:
    """
    Draw one round's PES diagram -- strictly best-effort, exactly as
    ``t3.pdep.explorer.arkane.ArkaneExplorer._draw_pes_diagram`` already treats it.

    Every exception is swallowed and logged, not only the ``ValueError`` a parser/data refusal
    raises: a matplotlib backend problem, a full disk, or a read-only output directory has nothing
    to do with whether the exploration this diagram describes succeeded, and must never flip a
    round that genuinely explored into an aborted loop. The diagram describes the result; it is
    not part of it.

    Args:
        explored_network_path (str): The freshly explored network to draw from.
        diagram_path (str): Where to write the diagram.
        logger: A T3 ``Logger``, or ``None``.

    Returns:
        str | None: ``diagram_path`` if the diagram was drawn, ``None`` if drawing failed.
    """
    try:
        draw_pes_diagram(network_path=explored_network_path, output_path=diagram_path)
    except Exception as e:
        message = (f'PES loop: could not draw the diagram for {explored_network_path} '
                   f'({type(e).__name__}): {e}')
        _logger.warning(message)
        if logger is not None:
            logger.warning(message)
        return None
    return diagram_path


def run_pes_loop(config: PESLoopConfig, project_directory: str, qm_runner=None,
                 adopted_ts_labels: frozenset | None = None, logger=None) -> PESLoopResult:
    """
    Run the PES exploration loop to a stop.

    Args:
        config (PESLoopConfig): The loop's configuration.
        project_directory (str): The loop's project directory. Must be absolute (see
            ``t3.pdep.pes_rounds.round_paths``).
        qm_runner: An injected callable, ``qm_runner(candidates, paths, config, network_id,
            adopted=...) -> tuple[frozenset[str], frozenset[str]]`` -- ``(converged_ts_labels,
            queued_ts_labels)``, the network-local TS labels that converged this round and the
            ones the runner actually queued (not necessarily every candidate it was handed -- N3).
            ``adopted`` is passed EVERY round: a dict of network-local TS label -> artifact path
            for every QM artifact already in hand before this round ran (reused from prior
            projects, or converged by an earlier round of this loop), which the runner must fold
            into this round's hybrid network alongside whatever converges now -- every round's
            hybrid carries the CUMULATIVE QM, never just that round's (a TS that ever converged
            must never revert to an RMG/ILT rate line in a later hybrid). ``None`` (the default)
            means "no QM this run" -- explore and draw the first round's
            diagram only, which is the honest behaviour when no QM runner is configured, not a
            crash.
        adopted_ts_labels (dict | frozenset, optional): TS labels of the SEED network
            (``config.pes.network``) already computed before this loop started (e.g. reused from
            an earlier T3 project), preferably as a MAPPING of TS label -> QM artifact path --
            the same shape ``config.reuse.from_t3_projects`` resolves to. Each label is
            translated to its structural channel key through the seed network (labels are
            positional and do not survive re-exploration -- see the keying comment at the top of
            this function's body) and seeded as computed before round 0; a label with no
            unambiguous structural identity is warned about and re-decided each round instead.
            A label whose artifact is missing (a bare set of labels, or a path that does not
            exist) is NOT marked computed: without an artifact there is nothing to fold into any
            round's hybrid network, so the channel would keep its RMG/ILT rate line and the final
            diagram would show an ESTIMATED barrier while the run reported quantum chemistry for
            it -- and if such a label were the only candidate, round 0 would return
            ``'no_candidates'`` before ``qm_runner`` was ever called. It is warned about and
            queued normally instead. This is unioned with -- not replaced by -- whatever
            ``config.reuse.from_t3_projects`` itself resolves to (see below); passing this
            explicitly is for callers (e.g. tests) that already have the adopted set in hand and
            want to skip re-discovering it.
        logger: A T3 ``Logger``, or ``None``.

    When ``config.reuse.from_t3_projects`` is non-empty, this function itself calls
    ``t3.pdep.pes_qm.adopt_prior_qm`` before round 0 (``network_id`` derived the same way
    ``t3.pdep.parser.parse_pdep_network_file`` would derive it from ``config.pes.network``, i.e.
    ``Path(config.pes.network).stem``), handing it the seed network's structural channel keys
    (``t3.pdep.pes_rounds.channel_keys_by_ts_label``) to match prior captures against, and seeds
    every adopted channel as computed. A channel adopted this way is therefore never offered to
    ``qm_runner`` (its per-round label is translated into ``split_qm_candidates``'s computed set)
    -- it is not merely available for the caller to use, it actually prevents the redundant ARC
    submission.

    Returns:
        PESLoopResult: The rounds that ran and why the loop stopped.

    Raises:
        ValueError: If ``config.pes`` cannot build a valid ``PDepExplorerConfig`` (e.g. no
            ``bath_gas``), propagated from ``PDepExplorerConfig.__post_init__``.
        FileNotFoundError: If a successfully explored network's own output file cannot be read,
            propagated from ``parse_pdep_network_file``. A round-N ``qm_runner`` that fails to
            write its hybrid network file is NOT this case -- that is caught explicitly and
            returned as ``PES_LOOP_FAILED`` rather than left to raise (see the file-existence check
            at the top of every round after the first).
    """
    rounds = []
    # Cross-round memory is keyed STRUCTURALLY (t3.pdep.pes_rounds.structural_channel_key), never
    # by TS or reaction label: every label Arkane writes into a network file is purely positional
    # (rmgpy/rmg/pdep.py:854 replaces every path reaction's transition state with a fresh,
    # label-less object before each write, and path_reactions is append-with-remove), so a
    # label-keyed carry can attach a computed barrier to the WRONG saddle point after a pruning or
    # discovery renumbers the file. computed_channels is every channel whose QM is already in
    # hand; qm_artifacts_by_channel maps those that have an artifact to its path in the capture
    # where it was FIRST captured (or the prior project's own capture, for adopted ones) -- the
    # ORIGINATING capture, whose verified manifest still exists, never a later round's re-vendored
    # copy. Every round translates both through the freshly explored network's own
    # channel_keys_by_ts_label map, so the runner always speaks THAT round's labels (defect-3
    # fix: without the per-round fold, a TS converged in round N silently reverts to an RMG/ILT
    # Arrhenius line in round N+1's hybrid -- the fail-open shape t3.pdep.hybrid's invariant 2
    # exists to prevent).
    computed_channels = set()
    qm_artifacts_by_channel = dict()
    # Initialised unconditionally, not only inside the branch that needs it: the two consumers
    # below are guarded by exactly the same condition, so an uninitialised read is unreachable --
    # but a reader (and a static analyser: CodeQL flagged this) cannot see that without checking
    # three conditions against each other. The seed network is still parsed only when something
    # actually needs its keys; an empty map is the truthful value when nothing does.
    seed_channel_keys = dict()
    seed_adoption_keys = dict()
    if adopted_ts_labels or config.reuse.from_t3_projects:
        # require_reactions=False: config.pes.network is the seed, which is legitimately
        # source-only on round 0 (no reaction(...) yet -- the explorer creates them). This parse
        # only needs the seed's species/channel keys for adoption matching, never its reactions.
        seed_network = parse_pdep_network_file(path=config.pes.network, require_reactions=False)
        seed_channel_keys = channel_keys_by_ts_label(seed_network)
        seed_adoption_keys = adoption_channel_keys_by_ts_label(seed_network)
    if adopted_ts_labels:
        # A mapping supplies the artifact for each label; a bare iterable supplies none, and every
        # label in it therefore resolves to None below and is refused (see the fail-closed check).
        adopted_artifacts = (dict(adopted_ts_labels) if isinstance(adopted_ts_labels, dict)
                             else {ts_label: None for ts_label in adopted_ts_labels})
        for ts_label, artifact_path in sorted(adopted_artifacts.items()):
            key = seed_channel_keys.get(ts_label)
            if key is None:
                message = (f'PES loop: adopted_ts_labels names {ts_label!r}, but that transition '
                           f'state has no unambiguous structural channel identity in '
                           f'{config.pes.network!r}; it cannot be carried and will be re-decided '
                           f'each round.')
                _logger.warning(message)
                if logger is not None:
                    logger.warning(message)
                continue
            # Fail closed. Marking a channel computed WITHOUT an artifact to fold is a claim the
            # loop cannot honour: qm_runner is never handed anything for it, so the hybrid network
            # keeps that channel's RMG/ILT rate line and the final diagram shows an ESTIMATED
            # barrier -- while the run reports the channel as having quantum chemistry. Worse, if
            # such a label is the only candidate, the round returns 'no_candidates' before the
            # runner is even called. Refusing to mark it costs a re-computation; accepting it
            # costs the correctness of the number this whole loop exists to produce.
            if not artifact_path:
                message = (f'PES loop: adopted_ts_labels names {ts_label!r} with no artifact path, '
                           f'so there is nothing to fold into this run\'s hybrid network. Marking '
                           f'it computed would leave that channel on its RMG estimate while '
                           f'reporting it as computed, so it is NOT marked computed and will be '
                           f'queued normally. Pass a mapping of TS label -> artifact path (the '
                           f'way config.reuse.from_t3_projects resolves one) to actually reuse it.')
                _logger.warning(message)
                if logger is not None:
                    logger.warning(message)
                continue
            if not os.path.isfile(artifact_path):
                message = (f'PES loop: adopted_ts_labels gives {ts_label!r} the artifact '
                           f'{artifact_path!r}, which does not exist; it is NOT marked computed '
                           f'and will be queued normally.')
                _logger.warning(message)
                if logger is not None:
                    logger.warning(message)
                continue
            computed_channels.add(key)
            qm_artifacts_by_channel[key] = artifact_path
    if config.reuse.from_t3_projects:
        network_id = Path(config.pes.network).stem
        # Adoption is matched on the FAMILY-QUALIFIED key, never the endpoints-only one the
        # within-run carry uses: the two files compared here are unrelated, so a different pathway
        # between the same endpoints would key identically and the prior artifact would land on a
        # saddle point it was never computed for. See
        # t3.pdep.pes_rounds.adoption_channel_keys_by_ts_label. freq_level goes with sp_level
        # because ARC's own model_chemistry string depends on both.
        adopted = adopt_prior_qm(config.reuse.from_t3_projects, network_id, config.qm.sp_level,
                                 seed_adoption_keys,
                                 freq_level=config.qm.freq_level)
        if adopted:
            message = (f"PES loop: reusing {len(adopted)} prior QM result(s) for network "
                      f"{network_id!r}: {sorted(adopted)}.")
            _logger.info(message)
            if logger is not None:
                logger.info(message)
        for ts_label, artifact_path in adopted.items():
            # adopt_prior_qm only ever returns labels taken from the keys of seed_channel_keys.
            key = seed_channel_keys[ts_label]
            computed_channels.add(key)
            qm_artifacts_by_channel[key] = artifact_path
    # Whether round 0 starts with adopted ARTIFACTS to fold (reuse), as opposed to labels merely
    # marked computed (adopted_ts_labels carries no artifact paths) -- the round-0 fold decision.
    prior_adopted_artifacts = bool(qm_artifacts_by_channel)
    current_network_path = config.pes.network
    max_rounds = config.termination.max_rounds
    # Whether the round that ran last folded QM into a hybrid network -- read only by the
    # final-draw pass below, which is reachable only after at least one round ran.
    last_round_made_progress = False

    for round_index in range(max_rounds):
        paths = round_paths(project_directory, round_index)
        os.makedirs(paths.root, exist_ok=True)

        # Rounds after the first explore a hybrid network that a previous round's qm_runner is
        # contractually obliged to have written (see hybrid_network_path). If it did not, letting
        # that reach parse_pdep_network_file would raise a bare FileNotFoundError here -- after the
        # PREVIOUS round already spent real ARC time -- with no round record explaining why. Fail
        # the round explicitly instead, naming the contract that was broken.
        if round_index > 0 and not os.path.isfile(current_network_path):
            reason = (f"round {round_index}: the expected hybrid network file "
                     f"{current_network_path!r} does not exist. The previous round's qm_runner "
                     "must write it (see t3.pdep.pes_rounds.hybrid_network_path) before returning.")
            prior = rounds[-1] if rounds else None
            rounds.append(RoundRecord(index=round_index, network_path=current_network_path,
                                      diagram_path=None, queued_ts_labels=(), skipped=(),
                                      status=PES_LOOP_FAILED, reason=reason))
            return PESLoopResult(rounds=tuple(rounds), status=PES_LOOP_FAILED, reason=reason,
                                 final_network_path=prior.network_path if prior else None,
                                 final_diagram_path=prior.diagram_path if prior else None)

        explorer_config = _build_explorer_config(config, project_directory, paths)
        # No selection is available in standalone mode: there is no PDepNetworkSelection to bind
        # this run to, only a network path and a config. The default admission policy
        # ('qualified_selection' with selection=None) is what explore_pdep_network's own contract
        # calls the 'ungated' outcome -- nothing gated the run because nothing was there to gate it
        # -- which is the truthful record here, not 'caller_admitted': that policy REQUIRES a
        # selection to bind the run to (explore_pdep_network raises ValueError without one), because
        # its whole point is recording that the caller looked at a selection and admitted it anyway.
        result = explore_pdep_network(network_path=current_network_path, config=explorer_config,
                                      logger=logger)

        if result.status != EXPLORATION_STATUS_SUCCEEDED:
            reason = '; '.join(result.reasons) if result.reasons else \
                f"exploration ended with status {result.status!r} and no stated reason."
            prior = rounds[-1] if rounds else None
            rounds.append(RoundRecord(index=round_index, network_path=current_network_path,
                                      diagram_path=None, queued_ts_labels=(), skipped=(),
                                      status=PES_LOOP_FAILED, reason=reason))
            # The same rule as the missing-hybrid branch above: this round failed, but rounds
            # 0..N-1 explored real networks and drew real diagrams, and those remain the best
            # result this run has. Reporting None here throws them away and makes the caller (e.g.
            # PES.py, which logs both paths) report a run that produced nothing at all.
            return PESLoopResult(rounds=tuple(rounds), status=PES_LOOP_FAILED, reason=reason,
                                 final_network_path=prior.network_path if prior else None,
                                 final_diagram_path=prior.diagram_path if prior else None)

        explored_network_path = result.network_paths[0]
        network = parse_pdep_network_file(explored_network_path)

        diagram_path = _draw_round_diagram(explored_network_path, paths.diagram, logger)

        # Translate the structural cross-round memory into THIS round's own labels through the
        # freshly explored network -- the only labels this round's file, runner, and hybrid speak.
        keys_by_ts = channel_keys_by_ts_label(network)
        declared_ts_labels = tuple(sorted(network.path_reactions_by_ts()))
        if declared_ts_labels and not keys_by_ts:
            # The network declares transition states but NOT ONE of them could be structurally
            # keyed (no parseable structures, or every channel duplicated). The refusal semantics
            # alone are fail-closed but operationally silent: nothing is ever "already computed",
            # so every round would re-submit the SAME transition states to ARC until max_rounds --
            # up to a full budget of duplicate QM spend, with only a log line to explain it. This
            # is a diagnosable fault of the network/exploration, so the ROUND fails, carrying the
            # unkeyable labels in its reason.
            reason = (f'round {round_index}: none of the transition states '
                      f'{list(declared_ts_labels)} declared by network {network.network_id!r} has '
                      f'an unambiguous structural channel identity (missing/unparseable species '
                      f'structures, or structurally duplicated channels). QM computed for them '
                      f'could not be carried across rounds, so every round would re-submit the '
                      f'same transition states to ARC; refusing to spend the budget on that.')
            rounds.append(RoundRecord(index=round_index, network_path=explored_network_path,
                                      diagram_path=diagram_path, queued_ts_labels=(),
                                      skipped=split_qm_candidates(network, frozenset()).skipped,
                                      status=PES_LOOP_FAILED, reason=reason))
            return PESLoopResult(rounds=tuple(rounds), status=PES_LOOP_FAILED, reason=reason,
                                 final_network_path=explored_network_path,
                                 final_diagram_path=diagram_path)
        computed_ts_labels = frozenset(
            ts_label for ts_label, key in keys_by_ts.items() if key in computed_channels)
        split = split_qm_candidates(network, computed_ts_labels)

        # C2: a round that adopted prior QM at round 0 must not early-return here just because it
        # has no NEW candidates to queue -- every candidate having been satisfied by reuse is the
        # headline success case for reuse, not "nothing to do". qm_runner is still called below
        # (with an empty candidates tuple) so it can fold the adopted artifacts into a hybrid
        # network; only when there is also no qm_runner to do that (no way to write a hybrid at
        # all) does this fall through to the diagram-only return below.
        round0_has_adopted = round_index == 0 and prior_adopted_artifacts
        if not split.candidates and not round0_has_adopted:
            status = PES_LOOP_NO_CANDIDATES if round_index == 0 else PES_LOOP_CONVERGED
            rounds.append(RoundRecord(index=round_index, network_path=explored_network_path,
                                      diagram_path=diagram_path, queued_ts_labels=(),
                                      skipped=split.skipped, status=status))
            return PESLoopResult(rounds=tuple(rounds), status=status, reason='',
                                 final_network_path=explored_network_path,
                                 final_diagram_path=diagram_path)

        if qm_runner is None:
            # Offered, not queued -- no runner ever runs, so nothing can be dropped, and this is
            # the honest record. No ME SA runs on this path either (there is no QM spend to
            # justify and no capture to evidence), so scope='sensitive' cannot rank here: the
            # record simply reports the first `limit` candidates in network-file order.
            offered = split.candidates[:config.qm.max_transition_states_per_round]
            rounds.append(RoundRecord(index=round_index, network_path=explored_network_path,
                                      diagram_path=diagram_path,
                                      queued_ts_labels=tuple(candidate.ts_label
                                                             for candidate in offered),
                                      skipped=split.skipped, status=PES_LOOP_DIAGRAM_ONLY,
                                      reason='no qm_runner configured: explored and drew the '
                                            'diagram only, nothing was computed.'))
            return PESLoopResult(rounds=tuple(rounds), status=PES_LOOP_DIAGRAM_ONLY,
                                 reason=rounds[-1].reason, final_network_path=explored_network_path,
                                 final_diagram_path=diagram_path)

        if split.candidates:
            # Defect-1 fix: every queued transition state must carry REAL sensitivity evidence --
            # capture_ts_artifacts refuses (fail-closed, by design) any captured artifact without
            # it, so a loop with no evidence source could never complete a round that converged
            # anything. The evidence is measured, never invented: this round's own Arkane ME SA
            # on the freshly explored network, through the same production seams T3's in-iteration
            # pipeline uses (see t3.pdep.pes_sa). An SA that fails, or leaves every candidate
            # without a finite row, fails the ROUND loudly -- a stall reported honestly, never a
            # fabricated coefficient shimmed past the capture guard.
            try:
                evidence = run_round_me_sensitivity(network_path=explored_network_path,
                                                    sa_dir=paths.sa,
                                                    method=config.pes.method,
                                                    timeout=config.pes.timeout,
                                                    logger=logger)
            except (OSError, ValueError) as e:
                reason = (f'round {round_index}: could not obtain master-equation sensitivity '
                          f'evidence for network {network.network_id!r}: {e} Every queued '
                          f'transition state must carry real sensitivity evidence '
                          f'(t3.pdep.capture refuses a captured artifact without it), so this '
                          f'round cannot queue QM.')
                rounds.append(RoundRecord(index=round_index, network_path=explored_network_path,
                                          diagram_path=diagram_path, queued_ts_labels=(),
                                          skipped=split.skipped, status=PES_LOOP_FAILED,
                                          reason=reason))
                return PESLoopResult(rounds=tuple(rounds), status=PES_LOOP_FAILED, reason=reason,
                                     final_network_path=explored_network_path,
                                     final_diagram_path=diagram_path)
            any_finite_evidence = any(candidate.ts_label in evidence
                                      for candidate in split.candidates)
            split = attach_sensitivity_evidence(split, evidence,
                                                min_delta_ln_k=config.qm.min_delta_ln_k)
            if not split.candidates and not round0_has_adopted:
                if not any_finite_evidence:
                    # The SA ran but measured nothing finite for ANY candidate -- an operational
                    # fault in the evidence pipeline, not a legitimate "nothing worth computing".
                    reason = (f'round {round_index}: the master-equation sensitivity analysis for '
                              f'network {network.network_id!r} measured no finite sensitivity '
                              f'evidence for any QM candidate, so nothing can be queued (see the '
                              f"round record's skipped entries for the per-candidate reasons).")
                    rounds.append(RoundRecord(index=round_index, network_path=explored_network_path,
                                              diagram_path=diagram_path, queued_ts_labels=(),
                                              skipped=split.skipped, status=PES_LOOP_FAILED,
                                              reason=reason))
                    return PESLoopResult(rounds=tuple(rounds), status=PES_LOOP_FAILED,
                                         reason=reason,
                                         final_network_path=explored_network_path,
                                         final_diagram_path=diagram_path)
                # Every remaining candidate measured a real response BELOW the floor: a measured
                # "nothing here is worth the QM spend", which is the same honest terminal outcome
                # as having no candidates at all -- each skip carries its measured value in the
                # round record. Never 'failed' (nothing malfunctioned) and never a queue of
                # below-floor transition states (their manifests would record a coefficient that
                # justified nothing). The reason clause is what separates this, at the summary
                # level, from "everything was already computed".
                status = PES_LOOP_NO_CANDIDATES if round_index == 0 else PES_LOOP_CONVERGED
                reason = (f'round {round_index}: every remaining QM candidate of network '
                          f'{network.network_id!r} measured a ln(k) response below the '
                          f'min_delta_ln_k floor ({config.qm.min_delta_ln_k:.3e}); see the round '
                          f"record's skipped entries for the measured values.")
                rounds.append(RoundRecord(index=round_index, network_path=explored_network_path,
                                          diagram_path=diagram_path, queued_ts_labels=(),
                                          skipped=split.skipped, status=status, reason=reason))
                return PESLoopResult(rounds=tuple(rounds), status=status, reason=reason,
                                     final_network_path=explored_network_path,
                                     final_diagram_path=diagram_path)

        candidates = _trim_candidates(split.candidates, config, logger)
        # The candidates are OFFERED, not yet queued: the round record's queued_ts_labels comes
        # from qm_runner's own reported queued half (below), never from this list, because
        # build_arc_input can silently drop a candidate the loop offered (N3).

        # Defect-3 fix: EVERY round's qm_runner call carries every previously computed artifact
        # (adopted from prior projects at round 0, converged by any earlier round of this loop),
        # uniformly, translated into THIS round's own labels through keys_by_ts -- the previous
        # round0-only special case is exactly what let QM die at the round boundary. A channel
        # carried in qm_artifacts_by_channel but absent from this round's network (pruned, or no
        # longer unambiguously keyable) is simply not folded this round; it stays in the memory.
        adopted_for_round = {
            ts_label: qm_artifacts_by_channel[key]
            for ts_label, key in keys_by_ts.items() if key in qm_artifacts_by_channel
        }
        newly_computed, actually_queued = qm_runner(candidates, paths, config, network.network_id,
                                                    adopted=adopted_for_round)
        for ts_label in newly_computed:
            key = keys_by_ts.get(ts_label)
            if key is None:
                # Converged but structurally unkeyable (see channel_keys_by_ts_label): it cannot
                # be carried without risking a wrong-saddle-point attribution, so it will be
                # re-decided (and possibly re-queued) next round -- duplicated spend, loudly, in
                # preference to a misattributed barrier.
                message = (f'PES loop: transition state {ts_label!r} of network '
                           f'{network.network_id!r} converged but has no unambiguous structural '
                           f'channel identity; its artifact cannot be carried across rounds.')
                _logger.warning(message)
                if logger is not None:
                    logger.warning(message)
                continue
            computed_channels.add(key)
            # Where this round's capture vendored the converged artifact -- a stable layout
            # invariant of t3.pdep.capture (see captured_qm_artifact_path). Recorded off the
            # ORIGINATING round's capture so later rounds always re-vendor from a directory whose
            # own verified manifest still exists.
            qm_artifacts_by_channel[key] = captured_qm_artifact_path(
                paths.capture, arc_ts_label(network.network_id, ts_label))
        queued_ts_labels = tuple(actually_queued)

        # Important: qm_runner's returned newly_computed excludes adopted TS labels by contract
        # (t3.pdep.pes_qm.arc_qm_runner's converged_ts_labels is what ARC itself converged THIS
        # round; adopted artifacts were already usable before this round even ran). Folding
        # `adopted` in was still real progress -- a hybrid network carrying them was written this
        # round (see t3.pdep.pes_qm.arc_qm_runner's C2) -- so "no progress" and "advance to the
        # hybrid file" must both look at newly_computed OR round0_has_adopted, not newly_computed
        # alone.
        made_progress = bool(newly_computed) or round0_has_adopted

        if config.termination.stop_when_no_new_ts and not made_progress:
            reason = ('qm_runner returned no newly-converged transition states and '
                     'termination.stop_when_no_new_ts is set: continuing would not make progress.')
            rounds.append(RoundRecord(index=round_index, network_path=explored_network_path,
                                      diagram_path=diagram_path, queued_ts_labels=queued_ts_labels,
                                      skipped=split.skipped, status=PES_LOOP_STALLED, reason=reason))
            return PESLoopResult(rounds=tuple(rounds), status=PES_LOOP_STALLED, reason=reason,
                                 final_network_path=explored_network_path,
                                 final_diagram_path=diagram_path)

        # A round that converges nothing AND adopted nothing correctly does NOT write a hybrid file
        # (qm_runner's empty-convergence contract -- see t3.pdep.pes_qm.arc_qm_runner). Advancing
        # current_network_path to hybrid_network_path anyway would make the next round's
        # file-existence check (above) blame qm_runner for "failing to write" a file it was never
        # supposed to write. So only advance to the hybrid file when qm_runner actually converged or
        # folded in something this round; otherwise stay on the network just explored and let
        # max_rounds bound the (otherwise honest) repetition.
        no_progress_reason = ('' if made_progress else
                              'qm_runner converged no new transition states this round; '
                              're-exploring the same network next round rather than advancing to '
                              'a hybrid file that was never written.')
        rounds.append(RoundRecord(index=round_index, network_path=explored_network_path,
                                  diagram_path=diagram_path, queued_ts_labels=queued_ts_labels,
                                  skipped=split.skipped, status='continuing',
                                  reason=no_progress_reason))

        # The next round explores the QM-informed surface: a real qm_runner (Task 6) writes this
        # round's hybrid network input file (t3.pdep.hybrid) to hybrid_network_path(paths,
        # network.network_id) from its captured statmech artifacts, folding the newly converged
        # (and any adopted) transition states in as Log(...) references. This driver trusts that
        # contract rather than writing the hybrid file itself -- that write belongs with the capture
        # step that knows what it captured -- and enforces it explicitly at the top of the next
        # round (above). When nothing converged or was adopted this round, there is no hybrid file
        # to trust -- see no_progress_reason.
        current_network_path = (hybrid_network_path(paths, network.network_id) if made_progress
                                else explored_network_path)
        last_round_made_progress = made_progress

    reason = f'the round budget ({max_rounds}) was exhausted without converging.'
    final_network_path = rounds[-1].network_path
    final_diagram_path = rounds[-1].diagram_path

    # THE defect this whole loop exists to eliminate, re-entering by the back door. The last
    # permitted round's hybrid network -- the one carrying every barrier ARC just converged -- is
    # written AFTER that round's own diagram was already drawn from its pre-QM explored network.
    # Returning rounds[-1] here therefore hands the caller (PES.py logs both paths, and a campaign
    # reads the diagram) the RMG-estimate surface, silently, in exactly the configuration a user is
    # most likely to run first (max_rounds: 1). So the loop completes one final
    # exploration-and-draw from that hybrid before returning. This charges NO further QM: no
    # candidate split, no ME SA, no qm_runner call, no ARC project directory -- only the explore
    # and the draw. It gets round_paths(max_rounds) as its own directory so its explorer output and
    # diagram never overwrite the last real round's, and it appends no RoundRecord: no round of QM
    # happened, and claiming one would misreport the budget that was actually spent.
    if last_round_made_progress and os.path.isfile(current_network_path):
        final_paths = round_paths(project_directory, max_rounds)
        os.makedirs(final_paths.root, exist_ok=True)
        final_result = explore_pdep_network(
            network_path=current_network_path,
            config=_build_explorer_config(config, project_directory, final_paths), logger=logger)
        if final_result.status == EXPLORATION_STATUS_SUCCEEDED and final_result.network_paths:
            final_network_path = final_result.network_paths[0]
            final_diagram_path = _draw_round_diagram(final_network_path, final_paths.diagram,
                                                     logger)
            reason = (f'{reason} The final network and diagram are the QM-informed re-exploration '
                      f'of round {max_rounds - 1}\'s hybrid network {current_network_path!r}, '
                      f'drawn after that round\'s quantum chemistry was folded in.')
        else:
            # Fail OPEN, loudly: the last round's own explored network and diagram remain the best
            # result this run has, and reporting them is truthful -- they simply predate this
            # round's QM. Turning a bonus re-exploration into PES_LOOP_FAILED would throw away
            # every round that did work.
            failure = ('; '.join(final_result.reasons) if final_result.reasons
                       else f'status {final_result.status!r}')
            message = (f'PES loop: could not re-explore the final round\'s hybrid network '
                       f'{current_network_path!r} to draw a QM-informed diagram ({failure}); '
                       f'reporting round {max_rounds - 1}\'s own explored network and diagram, '
                       f'which predate that round\'s quantum chemistry.')
            _logger.warning(message)
            if logger is not None:
                logger.warning(message)
            reason = f'{reason} {message}'

    return PESLoopResult(rounds=tuple(rounds), status=PES_LOOP_MAX_ROUNDS, reason=reason,
                         final_network_path=final_network_path,
                         final_diagram_path=final_diagram_path)
