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
from arc.level import Level, assign_frequency_scale_factor
from arc.main import ARC
from arc.statmech.arkane import get_arkane_model_chemistry

from t3.pdep.capture import CAPTURE_MANIFEST_FILE_NAME, capture_ts_artifacts, verify_capture
from t3.pdep.discovery import ARTIFACT_STATUS_USABLE
from t3.pdep.hashing import hash_file
from t3.pdep.hybrid import (QMEnergySettings, _read_qm_artifact, _vendor_qm_artifacts,
                            write_hybrid_network_input_file)
from t3.pdep.join import JOIN_STATUS_QUEUED, TSJoinRecord, arc_ts_label
from t3.pdep.parser import parse_pdep_network_file
from t3.pdep.pes_rounds import RoundPaths, channel_keys_by_ts_label, hybrid_network_path
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


def _vendor_adopted_artifacts(adopted: dict, capture_dir: str) -> dict[str, str]:
    """
    Copy adopted prior-QM artifacts, AND every log file they reference, into this round's own
    capture directory.

    ``write_hybrid_network_input_file``'s ``qm_artifacts_root`` guard only trusts artifact paths
    that live under the capture directory it is given -- a deliberate fail-closed security guard
    (see ``t3.pdep.hybrid``'s own docstring), which this function does not widen. An adopted
    artifact lives under a DIFFERENT run's capture directory, so it is vendored (copied) into an
    ``adopted/`` subdirectory of THIS round's own capture directory before being folded into
    ``qm_transition_states`` -- the existing guard then covers it for free, since ``adopted/``
    lives under the same ``capture_dir`` passed as ``qm_artifacts_root``.

    A captured artifact's own ``Log(...)`` arguments are relative to the artifact's OWN directory
    (see ``t3.pdep.hybrid.write_hybrid_network_input_file``'s docstring), so copying only the
    ``.py`` file -- without also copying the log tree it points at -- leaves those references
    dangling: ``write_hybrid_network_input_file``'s own pre-flight, run later on the vendored copy,
    would then raise a ``ValueError`` for a ``Log(...)`` file that "does not exist". This function
    therefore reuses ``t3.pdep.hybrid``'s own artifact-reading and vendoring machinery
    (``_read_qm_artifact``, ``_vendor_qm_artifacts``) -- the same functions ``t3.pdep.capture``
    already imports for the same reason -- instead of hand-copying just the ``.py`` file.

    Every adopted artifact was, by construction, itself produced by a PRIOR call to
    ``t3.pdep.capture.capture_ts_artifacts``, which always vendors an artifact to exactly
    ``<that capture's own capture_dir>/qm/<arc_ts_label>.py`` (see
    ``t3.pdep.capture``'s ``_CAPTURE_QM_SUBDIR`` and ``verify_capture``). So the artifact's own
    ``allowed_log_root`` for confinement -- the prior run's own capture directory -- is always
    recoverable as two directories up from its ``artifact_path``, with no need for
    ``adopt_prior_qm`` to separately track or pass it.

    Args:
        adopted (dict): Network-local TS label -> source artifact path, as returned by
                        ``adopt_prior_qm``.
        capture_dir (str): This round's own capture directory (``capture_result.capture_dir``).

    Returns:
        dict[str, str]: Network-local TS label -> the vendored (copied) artifact path, now inside
        ``capture_dir``.

    Raises:
        FileNotFoundError: An adopted artifact's source path no longer exists -- fail closed
                           rather than silently dropping it from the hybrid network.
        ValueError: An adopted artifact references a ``Log(...)`` file that no longer exists, or
                   that resolves outside its own prior capture directory -- see
                   ``t3.pdep.hybrid._read_qm_artifact``.
    """
    if not adopted:
        return dict()
    for ts_label, source_path in adopted.items():
        if not os.path.isfile(source_path):
            raise FileNotFoundError(
                f"PES loop: adopted artifact for {ts_label!r} no longer exists at "
                f"{source_path!r}; refusing to fold a missing artifact into this round's hybrid "
                f"network.")
    artifact_infos = {
        ts_label: _read_qm_artifact(
            ts_label=ts_label,
            artifact_path=source_path,
            allowed_log_root=os.path.dirname(os.path.dirname(os.path.abspath(source_path))),
        )
        for ts_label, source_path in adopted.items()
    }
    adopted_dir = os.path.join(capture_dir, 'adopted')
    _vendor_qm_artifacts(artifact_infos=artifact_infos, qm_dir=adopted_dir, dest_dir=capture_dir)
    return {ts_label: os.path.join(adopted_dir, f'{ts_label}.py') for ts_label in adopted}


def _adopted_energy_settings(adopted: dict) -> dict | None:
    """
    Recover the frozen energy-reference settings adopted artifacts were computed under, from their
    OWN prior capture manifest(s).

    ``t3.pdep.capture.capture_ts_artifacts`` only reads ``energy_settings`` from the live ARC
    project when THIS round captured at least one new artifact (``records_with_artifact``) -- a
    round that adopted a prior result but converged nothing new gets ``CaptureResult.energy_settings
    = None`` from that call, even though the adopted artifacts plainly WERE computed under some
    settings. Those settings are recorded in each adopted artifact's own (prior) capture manifest --
    the same manifest ``adopt_prior_qm`` already re-validated with ``verify_capture`` before
    adopting from it -- so this function re-derives the artifact's prior capture directory the same
    way ``_vendor_adopted_artifacts`` does (two directories up from ``artifact_path``, since a
    captured artifact always lives at ``<capture_dir>/qm/<label>.py``) and reads that manifest's
    settings directly, rather than treating "nothing new this round" as "no settings available".

    Args:
        adopted (dict): Network-local TS label -> ORIGINAL (pre-vendoring) adopted artifact path,
                        as returned by ``adopt_prior_qm``.

    Returns:
        dict, optional: The frozen ``energy_settings`` block, or ``None`` if ``adopted`` is empty.

    Raises:
        ValueError: Two adopted artifacts come from prior captures whose manifests disagree on
                   ``energy_settings`` -- refuse to guess which one is authoritative for this
                   round's hybrid network rather than silently picking one.
    """
    if not adopted:
        return None
    settings_by_capture_dir = dict()
    for source_path in adopted.values():
        capture_dir = os.path.dirname(os.path.dirname(os.path.abspath(source_path)))
        if capture_dir in settings_by_capture_dir:
            continue
        settings_by_capture_dir[capture_dir] = verify_capture(capture_dir).energy_settings
    settings_values = list(settings_by_capture_dir.values())
    first = settings_values[0]
    for capture_dir, other in list(settings_by_capture_dir.items())[1:]:
        if other != first:
            raise ValueError(
                f"PES loop: adopted artifacts come from prior captures with conflicting "
                f"energy_settings ({first!r} vs {other!r} at {capture_dir!r}); refusing to guess "
                f"which one is authoritative for this round's hybrid network.")
    return first


def arc_qm_runner(candidates: tuple, paths: RoundPaths, config: PESLoopConfig,
                  network_id: str, adopted: dict | None = None) -> tuple[frozenset[str], frozenset[str]]:
    """
    Run ARC on one round's QM candidates and fold whatever converges into a hybrid network file.

    This is the real ``qm_runner`` ``t3.pdep.pes_loop.run_pes_loop`` injects. It: parses this
    round's explored network to get species structures, builds the ARC input (``build_arc_input``),
    runs ARC via the same construction ``t3.main.T3.run_arc`` uses (``ARC(**arc_kwargs)`` then
    ``arc.execute()``), captures whatever transition-state QM artifacts converged
    (``t3.pdep.capture.capture_ts_artifacts``) against this round's own capture directory, and
    writes this round's hybrid Arkane network input file
    (``t3.pdep.hybrid.write_hybrid_network_input_file``) to exactly
    ``t3.pdep.pes_rounds.hybrid_network_path(paths, network_id)`` -- the loop checks for that exact
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
        adopted (dict, optional): Network-local TS label -> artifact path for every QM artifact
                                  already in hand before this round ran -- reused from a prior
                                  project (``t3.pdep.pes_qm.adopt_prior_qm``) or converged by an
                                  earlier round of this loop (``t3.pdep.pes_loop`` passes its
                                  cumulative map every round). Vendored into this round's own
                                  capture directory and folded into ``qm_transition_states``
                                  alongside whatever ARC converged this round, so every round's
                                  hybrid carries the CUMULATIVE QM. With no candidates to queue at
                                  all, ARC and capture are skipped entirely and this round's whole
                                  job is vendoring these and writing the hybrid.

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

    Raises:
        ValueError: If any candidate that would be queued carries no sensitivity evidence
                   (``QMCandidate.coefficient``/``delta_ln_k`` is ``None``) -- refused BEFORE ARC
                   runs, because ``capture_ts_artifacts`` would refuse the artifact AFTER the QM
                   spend; anything ``capture_ts_artifacts`` itself raises; anything the hybrid
                   write raises.
        KeyError: If the capture did not vendor a copy of ``network_id`` even though something
                 converged for it.
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

    # Defect-1 fix (fail-closed, BEFORE any QM is spent): every queued transition state must carry
    # the real sensitivity evidence that justified selecting it -- capture_ts_artifacts (below)
    # refuses, by deliberate design, any captured artifact without a finite coefficient/delta_ln_k,
    # and that refusal would land AFTER ARC already ran, losing the QM result with no second
    # chance. run_pes_loop stamps this evidence from the round's own master-equation SA
    # (t3.pdep.pes_sa via t3.pdep.pes_rounds.attach_sensitivity_evidence); a caller driving this
    # function directly must supply QMCandidate.coefficient/delta_ln_k itself. Refused here rather
    # than defaulted: no rate-determining parameter may ever be invented.
    missing_evidence = sorted(candidate.ts_label for candidate in queued_candidates
                              if candidate.coefficient is None or candidate.delta_ln_k is None)
    if missing_evidence:
        raise ValueError(
            f'Refusing to queue transition state(s) {missing_evidence} of network {network_id!r}: '
            f'they carry no sensitivity evidence (QMCandidate.coefficient/delta_ln_k is None). '
            f'capture_ts_artifacts would refuse their artifacts AFTER the QM spend, so this is '
            f'refused BEFORE ARC runs. Evidence comes from the round\'s own master-equation '
            f'sensitivity analysis (t3.pdep.pes_sa.run_round_me_sensitivity, stamped by '
            f't3.pdep.pes_rounds.attach_sensitivity_evidence); it is never defaulted or invented.')

    join_records = [
        TSJoinRecord(
            network_id=network_id,
            network_ts_label=candidate.ts_label,
            status=JOIN_STATUS_QUEUED,
            arc_ts_label=arc_ts_label(network_id, candidate.ts_label),
            path_reaction_labels=(candidate.path_reaction.label,) if candidate.path_reaction is not None else (),
            # The sensitivity evidence that justified selecting this transition state, frozen onto
            # the join record at queue time exactly as t3.main.T3.queue_pdep_transition_states
            # freezes it -- the guaranteed-present-by-now evidence the refusal above enforced.
            coefficient=candidate.coefficient,
            delta_ln_k=candidate.delta_ln_k,
        )
        # N3: only candidates that actually reached ARC (queued_candidates), not every candidate
        # this function was handed -- a candidate build_arc_input dropped for missing structure
        # was never queued, so recording a JOIN_STATUS_QUEUED join record for it would be false.
        for candidate in queued_candidates
    ]

    if join_records:
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
        capture_result = capture_ts_artifacts(
            join_records=join_records,
            arc_project_directory=paths.arc_project,
            capture_dir=paths.capture,
            networks=networks,
            # Derived off the join records' own frozen evidence, mirroring
            # t3.main.T3._capture_pdep_ts_artifacts byte for byte -- the argument whose omission
            # was defect 1: without it every captured artifact carries
            # coefficient=None/delta_ln_k=None and capture's own verify_capture self-check
            # (rightly) refuses the staged capture.
            sensitivity_by_ts={record.key: (record.coefficient, record.delta_ln_k)
                               for record in join_records},
        )
    else:
        # Defect-2 fix: a round with nothing to queue (every candidate satisfied by adoption, or
        # every candidate dropped by build_arc_input) must not run ARC on an empty reaction list,
        # and MUST NOT call capture_ts_artifacts at all -- capture refuses an empty join_records
        # outright ("an iteration that queued no P-dep transition states must not call
        # capture_ts_artifacts at all"), and that refusal is correct: there is nothing to capture.
        # With adopted artifacts in hand, this round's whole job is vendoring them and writing the
        # hybrid; with none, there is nothing to do at all.
        capture_result = None
        if not adopted:
            _logger.info(f'Nothing was queued and nothing was adopted for network {network_id!r} '
                         f'this round; no hybrid network to write.')
            return frozenset(), queued_ts_labels

    converged_labels = frozenset(
        record.network_ts_label for record in capture_result.records
        if record.status == ARTIFACT_STATUS_USABLE
    ) if capture_result is not None else frozenset()
    # C2: vendor adopted artifacts into THIS round's own capture directory before deciding
    # whether there is anything to write -- a round that adopted a prior result but converged
    # nothing new still has a hybrid network worth writing (an adopted TS must not silently
    # revert to an RMG/ILT rate estimate; see this module's caller, t3.pdep.pes_loop). When this
    # round captured nothing at all (capture_result is None), the round's capture directory is
    # still the vendoring root -- created here, holding only the adopted/ subtree.
    vendor_root = capture_result.capture_dir if capture_result is not None else paths.capture
    if capture_result is None and not os.path.isdir(vendor_root):
        os.makedirs(vendor_root)
    vendored_adopted = _vendor_adopted_artifacts(adopted or dict(), vendor_root)
    if not converged_labels and not vendored_adopted:
        _logger.info(f'ARC converged no transition states for network {network_id!r} this round.')
        return frozenset(), queued_ts_labels

    qm_transition_states = {
        record.network_ts_label: record.artifact_path
        for record in capture_result.records
        if record.status == ARTIFACT_STATUS_USABLE
    } if capture_result is not None else dict()
    qm_transition_states.update(vendored_adopted)
    # C3: capture_ts_artifacts only reads energy_settings when THIS round captured at least one
    # new artifact (records_with_artifact) -- a round that adopted a prior result but converged
    # nothing new gets capture_result.energy_settings=None (or no capture_result at all), which
    # would otherwise crash QMEnergySettings.from_frozen. The adopted artifacts were computed
    # under settings recorded in their OWN prior manifest -- and those, not a guess, are the
    # correct provenance for them.
    energy_settings_frozen = capture_result.energy_settings if capture_result is not None else None
    if energy_settings_frozen is None:
        energy_settings_frozen = _adopted_energy_settings(adopted or dict())
    energy_settings = QMEnergySettings.from_frozen(energy_settings_frozen)
    if capture_result is not None:
        # I1: the hybrid network is built from the CAPTURE's own vendored network copy and
        # recorded method (CaptureResult.networks is the authoritative source -- see
        # t3.pdep.capture's own docstring), never from the live/original explored network --
        # capture_result.networks may differ from `networks` above if a concurrent/prior capture
        # already vendored this network under a different method or path.
        captured_network = capture_result.networks.get(network_id) if capture_result.networks else None
        if captured_network is None:
            raise KeyError(
                f"capture_ts_artifacts did not vendor a copy of network {network_id!r} even though "
                f"ARC converged {sorted(converged_labels)} for it -- CaptureResult.networks is "
                f"missing an entry this call needs to write the round's hybrid network file.")
        hybrid_source_path = os.path.join(capture_result.capture_dir, captured_network['captured_path'])
        hybrid_method = captured_network['method']
    else:
        # No capture happened, so no vendored network copy exists; the round's own explored
        # network (the file this function parsed above) is the only -- and the truthful -- source.
        hybrid_source_path = source_path
        hybrid_method = config.pes.method
    result = write_hybrid_network_input_file(
        source_path=hybrid_source_path,
        dest_path=hybrid_network_path(paths, network_id),
        method=hybrid_method,
        qm_transition_states=qm_transition_states,
        energy_settings=energy_settings,
        qm_artifacts_root=vendor_root,
    )
    _logger.info(f"Wrote a hybrid P-dep network input for {network_id!r} to {result.dest_path}: "
                f"QM/RRKM for {list(result.qm_ts_labels)}, RMG/ILT for {list(result.ilt_ts_labels)}.")
    for warning in result.warnings:
        _logger.warning(warning)
    return converged_labels, queued_ts_labels


def _normalized_model_chemistry(sp_level: str) -> str | None:
    """
    Normalize a T3 ``sp_level`` string into ARC's own ``model_chemistry`` text.

    A capture manifest's ``model_chemistry`` is ARC's byte-for-byte ``LevelOfTheory(...)`` repr,
    written via ``arc.statmech.arkane.get_arkane_model_chemistry`` -- including year-suffix
    resolution (e.g. ``'wb97xd'`` -> ``'wb97xd2023'``) no T3-side string normalizer could
    reproduce. Comparing a manifest's ``model_chemistry`` against a raw ``config.qm.sp_level``
    string is therefore never true on real ARC data; this function runs the same ARC
    normalization chain ARC itself used to write the manifest, so the comparison is apples to
    apples.

    ``None`` -- never an exception -- is the answer whenever ARC cannot resolve one. In
    particular, ``get_arkane_model_chemistry`` itself RAISES (``ValueError: freq_level required
    when freq_scale_factor isn't provided``) for a non-composite level with no tabulated
    frequency scale factor -- e.g. ``'wb97xd2023/def2tzvp'``, ARC's own normalized form, exactly
    what a user copying a level out of a prior capture manifest would write. This function only
    ever has the single ``sp_level`` string, no frequency level to hand ARC in that case, so the
    raise is converted to the documented "could not resolve" answer: the caller
    (``adopt_prior_qm``) then warns that nothing will match, rather than the whole run dying
    before round 0.

    Args:
        sp_level (str): A T3 level-of-theory string, e.g. ``config.qm.sp_level``, undashed.

    Returns:
        str | None: The normalized ``model_chemistry`` text, or ``None`` if ARC could not
        resolve one.
    """
    level = Level(repr=sp_level)
    scale = assign_frequency_scale_factor(level)
    try:
        return get_arkane_model_chemistry(level, freq_scale_factor=scale)
    except ValueError as e:
        _logger.warning(f'PES loop: ARC could not normalize the level of theory {sp_level!r} '
                        f'into a model_chemistry ({e}); treating it as unresolvable.')
        return None


def adopt_prior_qm(from_t3_projects: list, network_id: str, level_of_theory: str,
                   channel_key_by_ts_label: dict) -> dict[str, str]:
    """
    Reuse transition states a previous T3 (or standalone PES loop) run already computed.

    Walks each configured project directory looking for capture manifests
    (``t3.pdep.capture.CAPTURE_MANIFEST_FILE_NAME``) written anywhere under it -- a capture
    directory is nested under a run's own iteration subdirectory (e.g. ``t3.main.T3``'s
    ``self.paths['PDep capture']``) whose exact name and depth this function cannot assume, so it
    is discovered by walking rather than by a fixed relative path. A directory that yields a
    manifest is not descended into further -- a capture directory is never itself nested inside
    another capture.

    Each manifest found is re-validated with ``verify_capture`` before anything in it is trusted
    (the same re-hashing, tamper/torn-capture detection a live consumer gets -- see
    ``t3.main``'s own use of it). An entry is adopted only when ALL of the following hold:

    - the capture's recorded ``energy_settings['model_chemistry']``, once normalized through
      ``_normalized_model_chemistry``, matches the requested ``level_of_theory``: mixing levels
      of theory inside one barrier's rate makes it inconsistent, so a prior result computed at a
      different level is refused even though it exists and converged (never adopted "close
      enough"). A refusal is logged with the manifest path, its raw ``model_chemistry``, and the
      requested level, so a silent "nothing was reused" is always explainable.
    - ``record.network_id`` is NEVER a gate -- it is logged only. Arkane names TS artifacts
      positionally (``'TS{i+1}'``, see ``arkane/pdep.py``), so a pruned or reordered exploration
      can re-issue the same label for a different transition state, and ``network_id`` itself
      never matches across independent PES-loop runs (Arkane names its own output by index, e.g.
      ``network0_full.py``, disjoint from T3's ``network<digits>_<digits>`` convention) -- gating
      on it would refuse the exact case this function exists to serve. The real match is
      STRUCTURAL, and deliberately not by ``path_reaction_labels`` either: reaction labels are
      just as positional as TS labels (``arkane/pdep.py:665``), and a remove-then-append can even
      leave two ``reaction(label='reaction3')`` blocks in one file (which ``arkane/input.py``
      silently renames). Instead, the record's transition state is resolved to its channel in the
      capture's own VENDORED network copy (the authoritative ``networks`` block every verified
      capture carries), reduced to a structural channel key
      (``t3.pdep.pes_rounds.channel_keys_by_ts_label`` -- canonical species structures, direction
      insensitive, immune to renumbering and to per-project species label indices), and matched
      against THIS run's own keys, supplied via ``channel_key_by_ts_label``. A record whose
      vendored network cannot be parsed, or whose transition state has no unambiguous structural
      identity there, or whose key matches no local transition state, is refused (logged).

    If two prior captures offer conflicting artifacts for what structurally matches the same
    local TS label, adoption for that label is refused (not a last-write-wins guess) and a
    warning is logged. Likewise, if two of THIS run's own transition states share one structural
    key, neither can be matched unambiguously and both are excluded up front.

    A project directory that does not exist, or whose manifest fails ``verify_capture`` (corrupt,
    torn, tampered, or simply absent), is skipped with a logged warning rather than raised -- a
    stale path left in a config must not kill a run that would otherwise work.

    Args:
        from_t3_projects (list): T3 (or standalone PES loop) project directories to search for
                                 reusable prior QM, e.g. ``config.reuse.from_t3_projects``.
        network_id (str): This run's network id. Logged for traceability only -- never a gate
                          (see above).
        level_of_theory (str): The level of theory ``model_chemistry`` must match to be adopted
                               (e.g. ``config.qm.sp_level``), undashed.
        channel_key_by_ts_label (dict): THIS run's own network-local TS label -> structural
                                        channel key, from
                                        ``t3.pdep.pes_rounds.channel_keys_by_ts_label`` over the
                                        run's own (seed) network. The identity adoption is
                                        matched against.

    Returns:
        dict[str, str]: Network-local TS label -> the adopted artifact path (already resolved,
        by ``verify_capture``, to an absolute path inside its capture directory).
    """
    normalized_requested = _normalized_model_chemistry(level_of_theory)
    if normalized_requested is None:
        _logger.warning(
            f"PES loop: ARC could not resolve a model_chemistry for the requested level of "
            f"theory {level_of_theory!r}; every prior capture's model_chemistry will therefore "
            f"fail to match and nothing will be adopted from {from_t3_projects!r}.")
    # Invert THIS run's map up front, refusing local transition states whose keys collide: a key
    # shared by two local labels cannot be matched unambiguously, and picking either would be the
    # wrong-saddle-point misattribution this structural matching exists to prevent.
    local_label_by_key = dict()
    ambiguous_local_keys = set()
    for local_label, key in channel_key_by_ts_label.items():
        if key in local_label_by_key:
            ambiguous_local_keys.add(key)
        local_label_by_key[key] = local_label
    for key in ambiguous_local_keys:
        _logger.warning(f'PES loop: several of this run\'s own transition states share the '
                        f'structural channel key {key!r}; none of them can adopt prior QM '
                        f'unambiguously.')
        del local_label_by_key[key]
    adopted = dict()
    conflicted = set()
    for project_directory in from_t3_projects:
        if not os.path.isdir(project_directory):
            _logger.warning(f"PES loop: reuse project directory {project_directory!r} does not "
                            f"exist; skipping it rather than failing this run.")
            continue
        for root, dirs, files in os.walk(project_directory):
            if CAPTURE_MANIFEST_FILE_NAME not in files:
                continue
            dirs[:] = []  # a capture directory is never nested inside another capture
            manifest_path = os.path.join(root, CAPTURE_MANIFEST_FILE_NAME)
            try:
                verified = verify_capture(root)
            except ValueError as e:
                _logger.warning(f"PES loop: skipping capture at {root!r} while looking for prior "
                                f"QM to reuse -- it failed verification: {e}")
                continue
            energy_settings = verified.energy_settings or {}
            manifest_model_chemistry = energy_settings.get('model_chemistry')
            if manifest_model_chemistry != normalized_requested:
                _logger.info(
                    f"PES loop: refusing prior QM at {manifest_path!r} -- its model_chemistry "
                    f"{manifest_model_chemistry!r} does not match the requested level of theory "
                    f"{level_of_theory!r} (normalized: {normalized_requested!r}).")
                continue
            # The capture's own VENDORED network copy is the authoritative structural record of
            # what each of its transition states connects (verified.networks -- see
            # t3.pdep.capture); parse each referenced network once and reduce its transition
            # states to structural channel keys. A network copy that cannot be parsed yields no
            # keys, so every record referencing it is refused below (logged per record).
            capture_keys_by_network_id = dict()
            for capture_network_id, network_entry in (verified.networks or dict()).items():
                captured_path = (network_entry or dict()).get('captured_path')
                if not captured_path:
                    _logger.warning(
                        f"PES loop: the capture at {root!r} records no vendored copy "
                        f"('captured_path') for network {capture_network_id!r}; no prior QM can "
                        f"be adopted from that network's records.")
                    capture_keys_by_network_id[capture_network_id] = dict()
                    continue
                vendored_network_path = os.path.join(root, captured_path)
                try:
                    capture_keys_by_network_id[capture_network_id] = channel_keys_by_ts_label(
                        parse_pdep_network_file(path=vendored_network_path))
                except (OSError, ValueError) as e:
                    _logger.warning(
                        f"PES loop: could not parse the vendored network copy "
                        f"{vendored_network_path!r} of the capture at {root!r} ({e}); no prior QM "
                        f"can be adopted from that network's records.")
                    capture_keys_by_network_id[capture_network_id] = dict()
            for record in verified.ts_records:
                if record.status != ARTIFACT_STATUS_USABLE:
                    continue
                if record.network_id != network_id:
                    # Never a gate -- see this function's own docstring: Arkane's network_id is
                    # per-run and never matches across independent PES-loop runs by construction.
                    # Logged only, so a structural match is always explainable even when the
                    # network ids visibly differ.
                    _logger.debug(
                        f"PES loop: matching prior QM at {manifest_path!r} for a different "
                        f"network_id ({record.network_id!r} vs this run's {network_id!r}) on "
                        f"the structural channel key alone -- network_id is a log-only field, "
                        f"never a gate (see this function's docstring).")
                record_key = capture_keys_by_network_id.get(record.network_id, dict()).get(
                    record.network_ts_label)
                if record_key is None:
                    _logger.info(
                        f"PES loop: refusing prior QM for {record.network_ts_label!r} at "
                        f"{manifest_path!r} -- its transition state has no unambiguous structural "
                        f"channel identity in the capture's own vendored network copy, so it "
                        f"cannot be matched to any of this run's own transition states.")
                    continue
                local_label = local_label_by_key.get(record_key)
                if local_label is None:
                    _logger.info(
                        f"PES loop: refusing prior QM for {record.network_ts_label!r} at "
                        f"{manifest_path!r} -- its structural channel key {record_key!r} matches "
                        f"none of this run's own transition states.")
                    continue
                if local_label in conflicted:
                    _logger.info(
                        f"PES loop: refusing prior QM for {local_label!r} at {manifest_path!r} -- "
                        f"already refused earlier this call due to a conflicting prior capture.")
                    continue
                if local_label in adopted and adopted[local_label] != record.artifact_path:
                    _logger.warning(
                        f"PES loop: refusing to adopt prior QM for {local_label!r} -- multiple "
                        f"prior captures offer conflicting artifacts ({adopted[local_label]!r} "
                        f"vs {record.artifact_path!r}); reuse for this transition state is "
                        f"refused.")
                    del adopted[local_label]
                    conflicted.add(local_label)
                    continue
                adopted[local_label] = record.artifact_path
    return adopted
