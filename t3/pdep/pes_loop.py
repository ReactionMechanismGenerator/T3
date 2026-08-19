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
``qm_runner(candidates, paths, config, network_id) -> frozenset[str]`` lets every round-sequencing
decision (no-candidates, converged, stalled, diagram-only, max-rounds, failed-explore) be covered
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

from t3.pdep.api import explore_pdep_network
from t3.pdep.diagram import draw_pes_diagram
from t3.pdep.explorer.config import PDepExplorerConfig
from t3.pdep.explorer.result import EXPLORATION_STATUS_SUCCEEDED
from t3.pdep.parser import parse_pdep_network_file
from t3.pdep.pes_rounds import RoundPaths, round_paths, split_qm_candidates
from t3.schema import PESLoopConfig

_logger = logging.getLogger(__name__)


def hybrid_network_path(paths: RoundPaths, network_id: str) -> str:
    """
    Where a round's ``qm_runner`` must write its hybrid network input file.

    ``RoundPaths.hybrid`` (``t3.pdep.pes_rounds``) is a directory, not a file, and this loop needs
    a file path to hand the next round's explorer. It also needs that file's stem to carry
    ``network_id`` rather than the literal ``'hybrid'``, because ``parse_pdep_network_file`` derives
    ``network_id = Path(path).stem`` (``t3/pdep/parser.py:729``) -- every round writing to a
    ``hybrid.py`` stem would collapse distinct networks onto one ``network_id``, and with it one ARC
    job-label namespace (the exact failure ruling C4 exists to prevent).

    This is a module-level function, not an inline expression, so Task 6's real ``qm_runner`` can
    import and call it directly rather than re-deriving the convention -- the contract lives in
    code, not in a comment.

    Args:
        paths (RoundPaths): This round's paths.
        network_id (str): The network id to preserve in the file's stem (normally the just-explored
            network's own ``network_id``).

    Returns:
        str: The hybrid network file path this round's ``qm_runner`` must write to, and the path
        the next round explores from.
    """
    return os.path.join(paths.hybrid, f'{network_id}.py')

# Status values for PESLoopResult.status / RoundRecord.status.
#
# - 'converged'      every barriered channel the network declares now has QM (round > 0).
# - 'no_candidates'  the very first round already had nothing to compute (round == 0): a network
#                     that is barrierless end to end is not a loop that "converged" after doing
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
        queued_ts_labels (tuple): The network-local TS labels sent to ``qm_runner`` this round.
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
            if the very first round failed to explore.
        final_diagram_path (str | None): The last round's PES diagram, or ``None`` if no round ever
            managed to draw one.
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
        explorer='arkane',
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

    Standalone mode has no sensitivity-analysis dict to rank candidates by -- that ranking only
    exists inside a T3 iteration, where ``t3.pdep.budget`` has a whole field of networks and an SA
    result to weigh them against. So ``scope == 'sensitive'`` degrades to the network file's own
    order here, and that degradation is logged at WARNING level rather than applied silently: a
    caller who asked for sensitivity-ranked spend and silently got file order instead has no way to
    notice unless the degradation is loud.

    Args:
        candidates (tuple): The ``QMCandidate`` entries to trim, in network-file order.
        config (PESLoopConfig): The loop's configuration.
        logger: A T3 ``Logger``, or ``None``.

    Returns:
        tuple: The trimmed candidates, in network-file order.
    """
    limit = config.qm.max_transition_states_per_round
    if config.qm.scope == 'sensitive' and len(candidates) > limit:
        message = ("PES loop: qm.scope='sensitive' was requested, but standalone mode has no "
                   "sensitivity-analysis dict to rank candidates by. Degrading to the network "
                   f"file's own order and taking the first {limit} of {len(candidates)} "
                   f"candidate(s).")
        _logger.warning(message)
        if logger is not None:
            logger.warning(message)
    return candidates[:limit]


def run_pes_loop(config: PESLoopConfig, project_directory: str, qm_runner=None,
                 adopted_ts_labels: frozenset | None = None, logger=None) -> PESLoopResult:
    """
    Run the PES exploration loop to a stop.

    Args:
        config (PESLoopConfig): The loop's configuration.
        project_directory (str): The loop's project directory. Must be absolute (see
            ``t3.pdep.pes_rounds.round_paths``).
        qm_runner: An injected callable, ``qm_runner(candidates, paths, config, network_id) ->
            frozenset[str]``, returning the network-local TS labels that converged this round.
            ``None`` (the default) means "no QM this run" -- explore and draw the first round's
            diagram only, which is the honest behaviour when no QM runner is configured, not a
            crash.
        adopted_ts_labels (frozenset, optional): Network-local TS labels already computed before
            this loop started (e.g. reused from an earlier T3 project), seeded into
            ``computed_ts_labels`` before round 0.
        logger: A T3 ``Logger``, or ``None``.

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
    computed_ts_labels = frozenset(adopted_ts_labels) if adopted_ts_labels else frozenset()
    current_network_path = config.pes.network
    max_rounds = config.termination.max_rounds

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
                     "must write it (see t3.pdep.pes_loop.hybrid_network_path) before returning.")
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
            rounds.append(RoundRecord(index=round_index, network_path=current_network_path,
                                      diagram_path=None, queued_ts_labels=(), skipped=(),
                                      status=PES_LOOP_FAILED, reason=reason))
            return PESLoopResult(rounds=tuple(rounds), status=PES_LOOP_FAILED, reason=reason,
                                 final_network_path=None, final_diagram_path=None)

        explored_network_path = result.network_paths[0]
        network = parse_pdep_network_file(explored_network_path)

        diagram_path = None
        try:
            draw_pes_diagram(network_path=explored_network_path, output_path=paths.diagram)
            diagram_path = paths.diagram
        except ValueError as e:
            # A legitimate refusal (a species with no E0), not a bug -- record it and keep going
            # rather than let a diagram problem kill an otherwise-working exploration round.
            message = f'PES loop: could not draw the diagram for {explored_network_path}: {e}'
            _logger.warning(message)
            if logger is not None:
                logger.warning(message)

        split = split_qm_candidates(network, computed_ts_labels)

        if not split.candidates:
            status = PES_LOOP_NO_CANDIDATES if round_index == 0 else PES_LOOP_CONVERGED
            rounds.append(RoundRecord(index=round_index, network_path=explored_network_path,
                                      diagram_path=diagram_path, queued_ts_labels=(),
                                      skipped=split.skipped, status=status))
            return PESLoopResult(rounds=tuple(rounds), status=status, reason='',
                                 final_network_path=explored_network_path,
                                 final_diagram_path=diagram_path)

        candidates = _trim_candidates(split.candidates, config, logger)
        queued_ts_labels = tuple(candidate.ts_label for candidate in candidates)

        if qm_runner is None:
            rounds.append(RoundRecord(index=round_index, network_path=explored_network_path,
                                      diagram_path=diagram_path, queued_ts_labels=queued_ts_labels,
                                      skipped=split.skipped, status=PES_LOOP_DIAGRAM_ONLY,
                                      reason='no qm_runner configured: explored and drew the '
                                            'diagram only, nothing was computed.'))
            return PESLoopResult(rounds=tuple(rounds), status=PES_LOOP_DIAGRAM_ONLY,
                                 reason=rounds[-1].reason, final_network_path=explored_network_path,
                                 final_diagram_path=diagram_path)

        newly_computed = qm_runner(candidates, paths, config, network.network_id)
        computed_ts_labels = computed_ts_labels | newly_computed

        if config.termination.stop_when_no_new_ts and not newly_computed:
            reason = ('qm_runner returned no newly-converged transition states and '
                     'termination.stop_when_no_new_ts is set: continuing would not make progress.')
            rounds.append(RoundRecord(index=round_index, network_path=explored_network_path,
                                      diagram_path=diagram_path, queued_ts_labels=queued_ts_labels,
                                      skipped=split.skipped, status=PES_LOOP_STALLED, reason=reason))
            return PESLoopResult(rounds=tuple(rounds), status=PES_LOOP_STALLED, reason=reason,
                                 final_network_path=explored_network_path,
                                 final_diagram_path=diagram_path)

        # A round that converges nothing correctly does NOT write a hybrid file (qm_runner's
        # empty-convergence contract -- see t3.pdep.pes_qm.arc_qm_runner). Advancing
        # current_network_path to hybrid_network_path anyway would make the next round's
        # file-existence check (above) blame qm_runner for "failing to write" a file it was never
        # supposed to write. So only advance to the hybrid file when qm_runner actually converged
        # something this round; otherwise stay on the network just explored and let max_rounds
        # bound the (otherwise honest) repetition.
        no_progress_reason = ('' if newly_computed else
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
        # transition states in as Log(...) references. This driver trusts that contract rather than
        # writing the hybrid file itself -- that write belongs with the capture step that knows what
        # it captured -- and enforces it explicitly at the top of the next round (above). When
        # nothing converged this round, there is no hybrid file to trust -- see no_progress_reason.
        current_network_path = (hybrid_network_path(paths, network.network_id) if newly_computed
                                else explored_network_path)

    reason = f'the round budget ({max_rounds}) was exhausted without converging.'
    return PESLoopResult(rounds=tuple(rounds), status=PES_LOOP_MAX_ROUNDS, reason=reason,
                         final_network_path=rounds[-1].network_path,
                         final_diagram_path=rounds[-1].diagram_path)
