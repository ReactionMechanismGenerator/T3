"""
t3 pdep pes_qm module

The real ARC-backed ``qm_runner`` that ``t3.pdep.pes_loop.run_pes_loop`` injects and calls as
``qm_runner(candidates, paths, config, network_id) -> frozenset[str]``. Everything else in
``t3.pdep`` is exercised without ever touching ARC (see ``pes_loop.py``'s own module docstring for
why); this module is the one piece the loop cannot run without an ARC cluster behind it.

``build_arc_input`` is kept a pure function -- no I/O, no ARC import at call time beyond the plain
dict it returns -- specifically so the levels of theory, job types and TS label namespacing it
decides are unit-testable without running anything. ``arc_qm_runner`` is the thin, impure shell
around it: it writes that input, runs ARC via the exact mechanism ``t3.main.T3.run_arc`` already
uses, captures whatever transition-state QM artifacts converged, folds them into a hybrid Arkane
network input file, and reports back which network-local transition states are now usable.

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
from t3.pdep.pes_loop import hybrid_network_path
from t3.pdep.pes_rounds import RoundPaths
from t3.schema import PESLoopConfig

_logger = logging.getLogger(__name__)

# The name ARC's own input YAML is written under, inside the round's ARC project directory --
# mirrors t3.main.T3.run_arc's 'ARC input' path, kept local here since this loop's ARC project
# directory is per-round rather than per-T3-project.
ARC_INPUT_FILE_NAME = 'arc_input.yml'


def build_arc_input(candidates: tuple, paths: RoundPaths, config: PESLoopConfig, network_id: str) -> dict:
    """
    Build the ARC input dict for one round's quantum chemistry, with no I/O.

    Kept pure so the levels of theory, job types, and TS label namespacing this function decides
    are unit-testable without running ARC, without writing a file, and without an ARC project
    directory existing on disk.

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

    Returns:
        dict: The ARC input, suitable for ``arc.main.ARC(**arc_input)``.
    """
    qm = config.qm
    reactions = []
    for candidate in candidates:
        reactions.append({
            'ts_label': arc_ts_label(network_id, candidate.ts_label),
            'network_ts_label': candidate.ts_label,
            'family': candidate.family,
        })
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

    This is the real ``qm_runner`` ``t3.pdep.pes_loop.run_pes_loop`` injects. It: builds the ARC
    input (``build_arc_input``), runs ARC via the same construction ``t3.main.T3.run_arc`` uses
    (``ARC(**arc_kwargs)`` then ``arc.execute()``), captures whatever transition-state QM
    artifacts converged (``t3.pdep.capture.capture_ts_artifacts``) against this round's own
    capture directory, and writes this round's hybrid Arkane network input file
    (``t3.pdep.hybrid.write_hybrid_network_input_file``) to exactly
    ``t3.pdep.pes_loop.hybrid_network_path(paths, network_id)`` -- the loop checks for that exact
    file at the top of every round after the first and fails the round if it is absent.

    Args:
        candidates (tuple): The ``t3.pdep.pes_rounds.QMCandidate`` entries to queue, already
                            screened by ``split_qm_candidates``.
        paths (RoundPaths): This round's artifact layout.
        config (PESLoopConfig): The loop's configuration.
        network_id (str): The network file stem this round explored.

    Returns:
        frozenset: The network-local TS labels ARC converged this round (``status ==
                   t3.pdep.discovery.ARTIFACT_STATUS_USABLE`` only -- an ``UNVERIFIED`` artifact's
                   convergence is, by design, unknown, and is never treated as usable).
    """
    arc_input = dict(build_arc_input(candidates, paths, config, network_id))
    if not os.path.isdir(paths.arc_project):
        os.makedirs(paths.arc_project)
    arc_input.setdefault('project', network_id)
    reactions = arc_input.pop('reactions')

    arc_kwargs = {key: value for key, value in arc_input.items()}
    arc_input_path = os.path.join(paths.arc_project, ARC_INPUT_FILE_NAME)
    if not os.path.isfile(arc_input_path):
        save_yaml_file(path=arc_input_path, content={**arc_kwargs, 'reactions': reactions})

    arc = ARC(**arc_kwargs)
    try:
        arc.execute()
    except Exception as e:
        _logger.error(f'ARC crashed with {e.__class__}. Got the following error message:\n{e}')
        raise

    source_path = _explored_network_path(paths, network_id)
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
        for candidate in candidates
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
        return frozenset()

    qm_transition_states = {
        record.network_ts_label: record.artifact_path
        for record in capture_result.records
        if record.status == ARTIFACT_STATUS_USABLE
    }
    energy_settings = QMEnergySettings.from_frozen(capture_result.energy_settings)
    write_hybrid_network_input_file(
        source_path=source_path,
        dest_path=hybrid_network_path(paths, network_id),
        method=config.pes.method,
        qm_transition_states=qm_transition_states,
        energy_settings=energy_settings,
    )
    return converged_labels
