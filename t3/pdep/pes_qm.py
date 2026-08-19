"""
t3 pdep pes_qm module

The real ARC-backed ``qm_runner`` that ``t3.pdep.pes_loop.run_pes_loop`` injects and calls as
``qm_runner(candidates, paths, config, network_id) -> tuple[frozenset[str], frozenset[str]]`` --
``(converged_ts_labels, queued_ts_labels)``. Everything else in ``t3.pdep`` is exercised without
ever touching ARC (see ``pes_loop.py``'s own module docstring for why); this module is the one
piece the loop cannot run without an ARC cluster behind it.

``build_arc_input`` is kept a pure function -- no I/O, no ARC import at call time beyond the plain
dict it returns -- specifically so the levels of theory, job types and TS label namespacing it
decides are unit-testable without running anything. ``arc_qm_runner`` is the thin, impure shell
around it: it writes that input, runs ARC via the exact mechanism ``t3.main.T3.run_arc`` already
uses, captures whatever transition-state QM artifacts converged, folds them into a hybrid Arkane
network input file, and reports back which network-local transition states are now usable AND
which were actually queued -- ``build_arc_input`` can silently drop a candidate for missing
structure, so the two are not always the caller's own candidate list (N3).

Every ``QMCandidate`` this module receives has already been screened by
``t3.pdep.pes_rounds.split_qm_candidates``: barrierless channels, channels that already have QM,
and channels with no declared transition state are filtered out before candidates ever reach
``build_arc_input``. This module does not re-check any of that -- it trusts the split it was
handed, the same way ``t3.main`` trusts ``add_reaction`` upstream of ``queue_pdep_transition_states``.

TS labels are never handed to ARC in the network's own namespace: ``t3.pdep.join.arc_ts_label``
namespaces every label by ``network_id`` before it reaches ``build_arc_input``'s ``reactions``
list, because a network-local ``TS1`` is almost never ARC's ``TS1`` and guessing from list order
has no basis (see ``t3/pdep/join.py``'s module docstring).

This module never imports ``rmgpy`` or ``arkane``. It imports ``arc`` directly, the same as
``t3.main`` does -- ARC is an allowed dependency; only RMG-Py and Arkane are not.
"""

import logging
import os

from arc.common import save_yaml_file
from arc.main import ARC

from t3.pdep.capture import capture_ts_artifacts
from t3.pdep.discovery import ARTIFACT_STATUS_USABLE
from t3.pdep.hashing import hash_file
from t3.pdep.hybrid import QMEnergySettings, write_hybrid_network_input_file
from t3.pdep.join import JOIN_STATUS_QUEUED, TSJoinRecord, arc_ts_label
from t3.pdep.parser import parse_pdep_network_file
from t3.pdep.pes_loop import hybrid_network_path
from t3.pdep.pes_rounds import RoundPaths
from t3.schema import PESLoopConfig

_logger = logging.getLogger(__name__)

# The name ARC's own input YAML is written under, inside the round's ARC project directory --
# mirrors t3.main.T3.run_arc's 'ARC input' path, kept local here since this loop's ARC project
# directory is per-round rather than per-T3-project.
ARC_INPUT_FILE_NAME = 'arc_input.yml'


def build_arc_input(candidates: tuple, paths: RoundPaths, config: PESLoopConfig, network_id: str,
                    species_structures: dict) -> dict:
    """
    Build the ARC input dict for one round's quantum chemistry, with no file or network I/O.

    Kept pure (beyond logging a warning for a candidate that cannot be built) so the levels of
    theory, job types, and TS label namespacing this function decides are unit-testable without
    running ARC, without writing a file, and without an ARC project directory existing on disk.

    Every reaction dict this function emits mirrors ``t3.main.T3.build_pdep_path_reaction`` --
    ARC's own ``ARCReaction.from_dict`` (``arc/reaction/reaction.py``), when handed a top-level
    ``species_list``, resolves a reaction's ``r_species``/``p_species`` entries by matching their
    ``'label'`` against that list rather than building fresh ``ARCSpecies`` from the reaction
    dict itself -- so a reaction dict alone, without a top-level ``'species'`` list, can never
    construct a valid ``ARCReaction``. ``species_structures`` (``t3.pdep.parser.PDepNetwork``'s
    field of the same name) is what supplies the adjacency list every ``'species'`` entry needs.

    A candidate whose reactant or product labels have no entry in ``species_structures`` is
    dropped from the returned ``reactions`` list (with a warning logged) rather than raising --
    the same fail-open choice ``build_pdep_path_reaction`` makes, so one network species with an
    unreadable structure loses only the reactions it participates in, not every candidate in the
    round.

    ``compute_thermo`` is explicitly set to ``False`` on every emitted species dict. ARC defaults
    a bare species dict's ``compute_thermo`` to ``True``, which would queue a full thermodynamics
    job for every reactant and product species reaching ARC through this path -- unbounded by
    ``config.qm.max_transition_states_per_round`` and outside what this round's TS-only QM budget
    is meant to cover. This loop only ever needs the transition state itself; reactant/product
    thermo, if wanted, is a distinct, separately-budgeted concern.

    Args:
        candidates (tuple): The ``t3.pdep.pes_rounds.QMCandidate`` entries to queue, already
                            screened by ``split_qm_candidates`` (no barrierless channels, none
                            that already have QM, none with no declared transition state).
        paths (RoundPaths): This round's artifact layout (``t3.pdep.pes_rounds.round_paths``).
        config (PESLoopConfig): The loop's configuration; ``config.qm`` supplies the levels of
                                theory, job types and TS adapters, ``config.pes`` the ME method.
        network_id (str): The network file stem this round explored, e.g. ``'network1_1'``. Used
                          to namespace every candidate's TS label via ``arc_ts_label`` so it can
                          never collide with another network's transition state of the same
                          network-local name.
        species_structures (dict): Network species label -> RMG adjacency-list text, i.e. the
                                   explored network's own ``PDepNetwork.species_structures``.

    Returns:
        dict: The ARC input, suitable for ``arc.main.ARC(**arc_input)``.
    """
    qm = config.qm
    species_by_label = dict()
    reactions = []
    for candidate in candidates:
        path_reaction = candidate.path_reaction
        labels = tuple(path_reaction.reactants) + tuple(path_reaction.products)
        missing = [label for label in labels if not species_structures.get(label)]
        if missing:
            _logger.warning(f"Network {network_id!r} carries no adjacency list for species "
                            f"{missing} of transition state '{candidate.ts_label}', so that "
                            f"reaction cannot be sent to ARC.")
            continue
        for label in labels:
            species_by_label.setdefault(label, species_structures[label])
        reactions.append({
            # N1: no 'label' key here, deliberately. ARCReaction.from_dict only synthesizes
            # self.label from reactants/products (arc/reaction/reaction.py's
            # set_label_reactants_products) when self.label starts out falsy -- if 'reactants'
            # and 'products' are already populated (as they always are here) AND 'label' carries
            # the network's raw path-reaction label (e.g. 'reaction1', no '<=>'), that synthesis
            # is skipped and check_atom_balance later crashes doing self.label.split('<=>')[well].
            # t3.main.T3.build_pdep_path_reaction never sets a label either -- this matches that
            # precedent.
            'ts_label': arc_ts_label(network_id, candidate.ts_label),
            'family': candidate.family,
            'reactants': list(path_reaction.reactants),
            'products': list(path_reaction.products),
            'r_species': [{'label': label} for label in path_reaction.reactants],
            'p_species': [{'label': label} for label in path_reaction.products],
        })
    species = [
        {'label': label, 'adjlist': adjlist, 'compute_thermo': False}
        for label, adjlist in species_by_label.items()
    ]
    return {
        'project_directory': paths.arc_project,
        'opt_level': qm.opt_level,
        'freq_level': qm.freq_level,
        'sp_level': qm.sp_level,
        'scan_level': qm.scan_level,
        'irc_level': qm.irc_level,
        'ts_adapters': qm.ts_adapters,
        'job_types': {
            'rotors': qm.rotors,
            'irc': qm.irc,
        },
        'species': species,
        'reactions': reactions,
    }


def _explored_network_path(paths: RoundPaths, network_id: str) -> str:
    """
    Reconstruct this round's explored network file path from its ``network_id`` alone.

    Arkane's own explorer (``t3.pdep.explorer.arkane``) always writes the network it explored to
    ``paths.explorer_output/pdep/final/<name>.py``, and ``network_id`` is, by construction
    (``t3.pdep.parser.parse_pdep_network_file``), exactly that file's stem -- so the path is
    deterministic from ``paths`` and ``network_id`` alone, with no need to thread the network
    object itself through ``qm_runner``'s signature.

    Args:
        paths (RoundPaths): This round's artifact layout.
        network_id (str): The network file stem this round explored.

    Returns:
        str: The path to this round's explored network file.
    """
    return os.path.join(paths.explorer_output, 'pdep', 'final', f'{network_id}.py')


def arc_qm_runner(candidates: tuple, paths: RoundPaths, config: PESLoopConfig, network_id: str) -> frozenset:
    """
    Run ARC on one round's QM candidates and fold whatever converges into a hybrid network file.

    This is the real ``qm_runner`` ``t3.pdep.pes_loop.run_pes_loop`` injects. It: parses this
    round's explored network to get species structures, builds the ARC input (``build_arc_input``),
    runs ARC via the same construction ``t3.main.T3.run_arc`` uses (``ARC(**arc_kwargs)`` then
    ``arc.execute()``), captures whatever transition-state QM artifacts converged
    (``t3.pdep.capture.capture_ts_artifacts``) against this round's own capture directory, and
    writes this round's hybrid Arkane network input file
    (``t3.pdep.hybrid.write_hybrid_network_input_file``) to exactly
    ``t3.pdep.pes_loop.hybrid_network_path(paths, network_id)`` -- the loop checks for that exact
    file at the top of every round after the first and fails the round if it is absent.

    No ``_refuse_if_queued_ts_lack_artifacts``-equivalent guard exists here, unlike
    ``t3.main.T3._capture_pdep_ts_artifacts``. That guard exists because ``t3.main`` can capture
    the SAME ARC project directory twice across separate iterations (a fresh capture, or a
    skip-re-capture replay over an existing manifest), and ARC deletes and recreates
    ``calcs/statmech/kinetics/`` on its next rate pass -- so a second capture over a project ARC
    has since wiped can silently see fewer artifacts than the first one did. This loop never
    reuses a project directory: every round gets its own freshly created ``paths.arc_project``
    (``t3.pdep.pes_rounds.round_paths`` namespaces it by round index), and capture happens exactly
    once, synchronously, immediately after ``arc.execute()`` returns within this same call -- there
    is no second call, no later iteration, and no window in which ARC could have cleaned up this
    round's own artifacts before they are captured. The artifact-wipe race the guard defends
    against is structurally excluded here, not merely unobserved.

    Args:
        candidates (tuple): The ``t3.pdep.pes_rounds.QMCandidate`` entries to queue, already
                            screened by ``split_qm_candidates``.
        paths (RoundPaths): This round's artifact layout.
        config (PESLoopConfig): The loop's configuration.
        network_id (str): The network file stem this round explored.

    Returns:
        tuple: ``(converged_ts_labels, queued_ts_labels)``, both ``frozenset``. ``converged_ts_labels``
               is the network-local TS labels ARC converged this round (``status ==
               t3.pdep.discovery.ARTIFACT_STATUS_USABLE`` only -- an ``UNVERIFIED`` artifact's
               convergence is, by design, unknown, and is never treated as usable).
               ``queued_ts_labels`` is the network-local TS labels that actually reached ARC --
               N3: ``build_arc_input`` silently drops a candidate whose reactant/product has no
               entry in ``species_structures`` (logging a warning and continuing rather than
               crashing), so the candidates this function was HANDED are not always the ones that
               were actually queued. Reporting this back, rather than letting the caller assume
               every offered candidate was queued, is what lets
               ``t3.pdep.pes_loop.RoundRecord.queued_ts_labels`` stay honest.
    """
    source_path = _explored_network_path(paths, network_id)
    network = parse_pdep_network_file(path=source_path)

    arc_input = build_arc_input(candidates, paths, config, network_id, network.species_structures)
    # N3: build_arc_input may have silently dropped some candidates (missing structure). Only the
    # candidates whose namespaced TS label actually made it into arc_input['reactions'] were
    # queued -- recompute the set from what build_arc_input actually returned rather than
    # trusting that every candidate handed in survived.
    queued_arc_ts_labels = frozenset(reaction['ts_label'] for reaction in arc_input['reactions'])
    queued_candidates = tuple(
        candidate for candidate in candidates
        if arc_ts_label(network_id, candidate.ts_label) in queued_arc_ts_labels
    )
    queued_ts_labels = frozenset(candidate.ts_label for candidate in queued_candidates)

    if not os.path.isdir(paths.arc_project):
        os.makedirs(paths.arc_project)
    arc_kwargs = dict(arc_input)
    arc_kwargs.setdefault('project', network_id)

    arc = ARC(**arc_kwargs)
    arc_input_path = os.path.join(paths.arc_project, ARC_INPUT_FILE_NAME)
    if not os.path.isfile(arc_input_path):
        save_yaml_file(path=arc_input_path, content=arc.as_dict())
    try:
        arc.execute()
    except Exception as e:
        _logger.error(f'ARC crashed with {e.__class__}. Got the following error message:\n{e}')
        raise

    networks = {
        network_id: {
            'source_path': source_path,
            'source_sha256': hash_file(source_path),
            'method': config.pes.method,
        },
    }
    join_records = [
        TSJoinRecord(
            network_id=network_id,
            network_ts_label=candidate.ts_label,
            status=JOIN_STATUS_QUEUED,
            arc_ts_label=arc_ts_label(network_id, candidate.ts_label),
            path_reaction_labels=(candidate.path_reaction.label,) if candidate.path_reaction is not None else (),
        )
        # N3: only candidates that actually reached ARC (queued_candidates), not every candidate
        # this function was handed -- a candidate build_arc_input dropped for missing structure
        # was never queued, so recording a JOIN_STATUS_QUEUED join record for it would be false.
        for candidate in queued_candidates
    ]

    capture_result = capture_ts_artifacts(
        join_records=join_records,
        arc_project_directory=paths.arc_project,
        capture_dir=paths.capture,
        networks=networks,
    )

    converged_labels = frozenset(
        record.network_ts_label for record in capture_result.records
        if record.status == ARTIFACT_STATUS_USABLE
    )
    if not converged_labels:
        _logger.info(f'ARC converged no transition states for network {network_id!r} this round.')
        return frozenset(), queued_ts_labels

    qm_transition_states = {
        record.network_ts_label: record.artifact_path
        for record in capture_result.records
        if record.status == ARTIFACT_STATUS_USABLE
    }
    energy_settings = QMEnergySettings.from_frozen(capture_result.energy_settings)
    # I1: the hybrid network is built from the CAPTURE's own vendored network copy and recorded
    # method (CaptureResult.networks is the authoritative source -- see t3.pdep.capture's own
    # docstring), never from the live/original explored network -- capture_result.networks may
    # differ from `networks` above if a concurrent/prior capture already vendored this network
    # under a different method or path.
    captured_network = capture_result.networks.get(network_id) if capture_result.networks else None
    if captured_network is None:
        raise KeyError(
            f"capture_ts_artifacts did not vendor a copy of network {network_id!r} even though "
            f"ARC converged {sorted(converged_labels)} for it -- CaptureResult.networks is "
            f"missing an entry this call needs to write the round's hybrid network file.")
    captured_source_path = os.path.join(capture_result.capture_dir, captured_network['captured_path'])
    result = write_hybrid_network_input_file(
        source_path=captured_source_path,
        dest_path=hybrid_network_path(paths, network_id),
        method=captured_network['method'],
        qm_transition_states=qm_transition_states,
        energy_settings=energy_settings,
        qm_artifacts_root=capture_result.capture_dir,
    )
    _logger.info(f"Wrote a hybrid P-dep network input for {network_id!r} to {result.dest_path}: "
                f"QM/RRKM for {list(result.qm_ts_labels)}, RMG/ILT for {list(result.ilt_ts_labels)}.")
    for warning in result.warnings:
        _logger.warning(warning)
    return converged_labels, queued_ts_labels
