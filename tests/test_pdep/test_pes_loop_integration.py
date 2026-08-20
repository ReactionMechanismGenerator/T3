"""
t3 tests test_pdep test_pes_loop_integration module

N6 (round 3 of the Task 6 review): nothing anywhere paired the REAL ``t3.pdep.pes_loop.run_pes_loop``
with the REAL ``t3.pdep.pes_qm.arc_qm_runner``. Every ``test_pes_loop.py`` test injects a
hand-written ``qm_runner`` double whose two-tuple return can never disagree with what the loop
offered (that gap is exactly what round 3's N3 fix addressed on the loop side). Every
``test_pes_qm.py`` test drives ``arc_qm_runner`` directly, with no ``run_pes_loop`` around it, so it
can never observe whether the hybrid network file ``arc_qm_runner`` writes for round N is actually
the file ``run_pes_loop`` hands its round-(N+1) explorer -- each half's tests only ever check that
convention against itself.

This module is that missing pairing. It fakes only the three boundaries ``arc_qm_runner`` itself
cannot cross in a test process -- ``ARC`` (no cluster here), ``capture_ts_artifacts`` (no real QM
job to capture from), and ``write_hybrid_network_input_file`` (no real QM artifacts to fold into a
hybrid network) -- with argument-recording doubles. ``explore_pdep_network`` is also faked, the same
way every other ``pes_loop`` test fakes it (Arkane itself is out of scope for this loop's own test
suite, and always has been); everything downstream of that -- the real ``parse_pdep_network_file``,
the real ``split_qm_candidates``, the real ``draw_pes_diagram``, the real ``build_arc_input``, and
every join/capture computation inside ``arc_qm_runner`` -- runs unfaked, against the real
``network799_1`` fixture also used by ``test_pes_qm.py``.

Wiring ``run_pes_loop``'s ``qm_runner`` argument into ``t3.main`` stays deliberately out of scope --
this loop is standalone until corroborated (see ``pes_loop.py``'s own module docstring) -- this
module only corroborates that the loop and the real runner agree with each other.
"""

import os
import shutil

from t3.pdep.capture import (CaptureResult, capture_ts_artifacts, captured_qm_artifact_path,
                             verify_capture)
from t3.pdep.discovery import ARTIFACT_STATUS_USABLE, TSArtifactRecord
from t3.pdep.explorer.result import EXPLORATION_STATUS_SUCCEEDED, PDepExplorationResult
from t3.pdep.hashing import hash_file
from t3.pdep.hybrid import HybridNetworkResult
from t3.pdep.join import (JOIN_STATUS_QUEUED, TSJoinRecord, arc_ts_label,
                          expected_ts_artifact_path)
from t3.pdep.parser import parse_pdep_network_file
import t3.pdep.pes_qm as pes_qm
from t3.pdep.pes_loop import run_pes_loop
from t3.pdep.pes_qm import _explored_network_path, arc_qm_runner
from t3.pdep.pes_rounds import hybrid_network_path, round_paths, split_qm_candidates
from t3.schema import PESLoopConfig

_FIXTURE_NETWORK_PATH = os.path.join(
    os.path.dirname(__file__), '..', 'data', 'pdep_real_networks', 'network799_1', 'network799_1.py')

# The same frozen energy-settings shape t3.pdep.pes_qm.QMEnergySettings.from_frozen expects,
# mirroring test_pes_qm.py's own _FROZEN_ENERGY_SETTINGS.
_ENERGY_SETTINGS = {
    'model_chemistry': 'wb97xd/def2tzvp',
    'frequency_scale_factor': None,
    'use_hindered_rotors': True,
    'use_bond_corrections': False,
    'bond_correction_type': None,
    'atom_energies': None,
    'use_atom_corrections': True,
    'bond_additivity_corrections': None,
}


class _FakeARC(object):
    """Stands in for arc.main.ARC: no cluster, no execute -- just records what it was built with."""

    constructions = []
    execute_calls = 0

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        _FakeARC.constructions.append(kwargs)

    def execute(self):
        _FakeARC.execute_calls += 1

    def as_dict(self):
        return dict(self.kwargs)


def _config(network_path, max_rounds=3):
    return PESLoopConfig(pes={'network': network_path, 'source': ['A'],
                              'bath_gas': {'He': 1.0}},
                         termination={'max_rounds': max_rounds, 'stop_when_no_new_ts': False})


def test_real_run_pes_loop_wires_the_real_arc_qm_runner_across_rounds(tmp_path, monkeypatch):
    """Drive the real run_pes_loop with the real arc_qm_runner injected as qm_runner, faking only
    ARC, capture_ts_artifacts, and write_hybrid_network_input_file. Assert the loop completes at
    least two rounds, and -- the actual N6 contract -- that the exact hybrid file arc_qm_runner
    writes for round 0 is the exact file run_pes_loop hands its round-1 explorer, observed end to
    end through the real loop rather than inferred from two separately-passing unit tests."""
    # The fixture file's own stem ('network799_1') follows T3's own hybrid-file convention
    # (network<digits>_<digits>), not the real Arkane explorer's
    # ('^network(\\d+)_(full|reduced)\\.py$', t3.pdep.explorer.arkane._FINAL_NETWORK_FILENAME_RE).
    # This fake explorer stands in for that real explorer, so it must hand downstream a network_id
    # shaped the way Arkane's really is -- parse from a renamed copy rather than the fixture path
    # directly, so real_network.network_id and every path derived from it looks like what Arkane
    # would actually produce.
    fixture_network_path = os.path.join(str(tmp_path), 'network0_full.py')
    shutil.copyfile(_FIXTURE_NETWORK_PATH, fixture_network_path)
    real_network = parse_pdep_network_file(path=fixture_network_path)
    real_split = split_qm_candidates(real_network, frozenset())
    assert real_split.candidates, 'fixture must offer at least one real QM candidate'
    target = real_split.candidates[0]

    project_directory = str(tmp_path)
    # What run_pes_loop actually handed this fake as `network_path` each round -- the N6 contract
    # (round 1 gets round 0's hybrid file) is about THIS list, not about wherever the fake itself
    # chooses to write.
    received_network_paths = []

    def _fake_explore(*, network_path, config, logger=None):
        received_network_paths.append(network_path)
        # Not one of the three edges N6 names -- Arkane itself stays faked here exactly as it is
        # in every other pes_loop test -- but the file this hands downstream is a real, on-disk
        # copy of network799_1, written to the SAME canonical location the real Arkane explorer
        # would use (paths.explorer_output/pdep/final/<network_id>.py, round_paths(...)'s own
        # convention), because arc_qm_runner's real, unfaked _explored_network_path(paths,
        # network_id) recomputes that exact path independently rather than trusting whatever this
        # fake returns in network_paths -- so writing anywhere else would desync the two halves.
        round_index = len(received_network_paths) - 1
        paths = round_paths(project_directory, round_index)
        dest_path = _explored_network_path(paths, real_network.network_id)
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        shutil.copyfile(fixture_network_path, dest_path)
        return PDepExplorationResult(network_id=real_network.network_id,
                                     status=EXPLORATION_STATUS_SUCCEEDED,
                                     network_paths=(dest_path,))

    monkeypatch.setattr('t3.pdep.pes_loop.explore_pdep_network', _fake_explore)

    _FakeARC.constructions = []
    _FakeARC.execute_calls = 0
    monkeypatch.setattr(pes_qm, 'ARC', _FakeARC)

    write_hybrid_calls = []

    def _fake_capture_ts_artifacts(*, join_records, arc_project_directory, capture_dir, networks,
                                   sensitivity_by_ts):
        # Defect 1: the real capture_ts_artifacts refuses (via its verify_capture self-check) any
        # captured artifact whose join record carries no finite sensitivity evidence -- so even
        # this double must insist the runner actually passed it, keyed and finite, or the pairing
        # this module corroborates would still be one that cannot run against the real capture.
        assert sensitivity_by_ts, 'arc_qm_runner must pass the sensitivity evidence to capture'
        for key, (coefficient, delta_ln_k) in sensitivity_by_ts.items():
            assert coefficient is not None and delta_ln_k is not None, key
        os.makedirs(capture_dir, exist_ok=True)
        # Defect 3: the loop carries a converged TS across the round boundary as an adopted
        # artifact at the capture's own vendored path (captured_qm_artifact_path), and the REAL
        # _vendor_adopted_artifacts in round 1 fails closed on a missing file -- so this double
        # must actually write the artifact it claims was captured, not merely name it in a result
        # object (the same rule the hybrid double below already follows).
        vendored_path = captured_qm_artifact_path(
            capture_dir, arc_ts_label(real_network.network_id, target.ts_label))
        os.makedirs(os.path.dirname(vendored_path), exist_ok=True)
        with open(os.path.join(os.path.dirname(vendored_path), 'output.out'), 'w') as f:
            f.write('# stub ARC log output\n')
        with open(vendored_path, 'w') as f:
            f.write("geometry = Log('output.out')\n")
        return CaptureResult(
            capture_dir=capture_dir,
            manifest_path=os.path.join(capture_dir, 'manifest.yml'),
            records=(
                TSArtifactRecord(network_id=real_network.network_id,
                                 network_ts_label=target.ts_label,
                                 arc_ts_label=f'{real_network.network_id}_{target.ts_label}',
                                 status=ARTIFACT_STATUS_USABLE,
                                 artifact_path=os.path.join(capture_dir, f'{target.ts_label}.py'),
                                 converged=True,
                                 path_reaction_labels=(target.path_reaction.label,)),
            ),
            energy_settings=_ENERGY_SETTINGS,
            networks={
                real_network.network_id: {
                    'source_path': _FIXTURE_NETWORK_PATH,
                    'captured_path': f'{target.ts_label}.py',
                    'method': 'MSC',
                },
            },
        )

    def _fake_write_hybrid_network_input_file(*, source_path, dest_path, method,
                                              qm_transition_states, energy_settings,
                                              qm_artifacts_root):
        write_hybrid_calls.append({'source_path': source_path, 'dest_path': dest_path})
        # run_pes_loop's own round>0 guard checks for this file's existence on disk for real, so
        # the double must actually write it, not merely return a result object claiming to.
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        open(dest_path, 'w').close()
        return HybridNetworkResult(dest_path=dest_path, qm_ts_labels=tuple(qm_transition_states),
                                   ilt_ts_labels=(), vendored_files=(), warnings=())

    monkeypatch.setattr(pes_qm, 'capture_ts_artifacts', _fake_capture_ts_artifacts)
    monkeypatch.setattr(pes_qm, 'write_hybrid_network_input_file',
                        _fake_write_hybrid_network_input_file)

    config = _config(fixture_network_path, max_rounds=3)
    result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=arc_qm_runner)

    assert len(result.rounds) >= 2, (
        f'expected at least two rounds, got {len(result.rounds)}: '
        f'{[r.status for r in result.rounds]} ({result.status}: {result.reason})')
    assert len(write_hybrid_calls) >= 1
    round0_hybrid_path = write_hybrid_calls[0]['dest_path']
    assert os.path.isfile(round0_hybrid_path)
    assert len(received_network_paths) >= 2
    # The actual N6 contract: round 1's explorer must be handed EXACTLY the file arc_qm_runner
    # wrote for round 0.
    assert received_network_paths[1] == round0_hybrid_path


# ---------------------------------------------------------------------------------------------
# Real-capture pairing (this fix round). Everything below drives the REAL run_pes_loop with the
# REAL arc_qm_runner, the REAL capture_ts_artifacts (no capture double -- the gap every prior
# test on this branch shared, which is how defects 1-3 stayed green), the REAL Arkane
# master-equation SA (t3.pdep.pes_sa -- Arkane actually runs, once per round with candidates),
# the REAL adopt_prior_qm against a REAL prior capture, the REAL vendoring, and the REAL hybrid
# writer. Faked: ONLY the explorer (as everywhere in this suite) and ARC (execute() writes
# converged statmech artifacts into the round's own ARC project instead of submitting cluster
# jobs).
# ---------------------------------------------------------------------------------------------

_MODEL_CHEMISTRY = "LevelOfTheory(method='wb97xd2023',basis='def2tzvp',software='gaussian')"


def _write_artifact(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(os.path.join(os.path.dirname(path), 'output.out'), 'w') as f:
        f.write('stub quantum chemistry log\n')
    with open(path, 'w') as f:
        f.write("linear = False\n\nspinMultiplicity = 2\n\nenergy = Log('output.out')\n\n"
                "geometry = Log('output.out')\n\nfrequencies = Log('output.out')\n")


def _write_status(arc_dir, label, converged):
    output_dir = os.path.join(arc_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'status.yml'), 'a') as f:
        f.write(f"{label}:\n  convergence: {str(converged).lower()}\n  job_types: {{}}\n"
                f"  paths: {{}}\n  info: ''\n  errors: ''\n")


def _write_energy_settings(arc_dir):
    statmech_dir = os.path.join(arc_dir, 'calcs', 'statmech', 'kinetics')
    os.makedirs(statmech_dir, exist_ok=True)
    with open(os.path.join(statmech_dir, 'input.py'), 'w') as f:
        f.write(f"modelChemistry = {_MODEL_CHEMISTRY}\n\nuseHinderedRotors = True\n\n"
                "useAtomCorrections = True\n\n"
                "atomEnergies = {'C': -37.844411, 'H': -0.499818, 'O': -75.062219}\n\n"
                "useBondCorrections = False\n")
    output_dir = os.path.join(arc_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    output_yml = os.path.join(output_dir, 'output.yml')
    if not os.path.isfile(output_yml):
        with open(output_yml, 'w') as f:
            f.write('{}\n')


def _build_prior_capture(prior_project_dir, ts_labels):
    """A REAL prior capture for ``ts_labels`` (each joined to the fixture's own
    ``reaction<i>``), built by capture_ts_artifacts itself so adopt_prior_qm, the vendoring, and
    _adopted_energy_settings all run against the exact formats production writes."""
    arc_dir = os.path.join(prior_project_dir, 'arc_project')
    network_id = 'prior_run_network'
    records = []
    for ts_label in ts_labels:
        label = arc_ts_label(network_id, ts_label)
        record = TSJoinRecord(network_id=network_id, network_ts_label=ts_label,
                              status=JOIN_STATUS_QUEUED, arc_ts_label=label,
                              expected_artifact_path=expected_ts_artifact_path(arc_dir, label),
                              reason='Queued to ARC.', coefficient=-9.5e-05, delta_ln_k=0.795,
                              path_reaction_labels=(f'reaction{ts_label[2:]}',))
        _write_artifact(record.expected_artifact_path)
        _write_status(arc_dir, label, converged=True)
        records.append(record)
    _write_energy_settings(arc_dir)
    networks_dir = os.path.join(prior_project_dir, 'networks')
    os.makedirs(networks_dir, exist_ok=True)
    source_path = os.path.join(networks_dir, f'{network_id}.py')
    with open(source_path, 'w') as f:
        f.write(f"# stub RMG network file\nnetwork(label='{network_id}')\n")
    capture_ts_artifacts(
        records, arc_dir, os.path.join(prior_project_dir, 'capture'),
        networks={network_id: {'source_path': source_path,
                               'source_sha256': hash_file(source_path), 'method': 'MSC'}},
        sensitivity_by_ts={record.key: (record.coefficient, record.delta_ln_k)
                           for record in records})


class _ArtifactWritingFakeARC(object):
    """Stands in for arc.main.ARC at exactly the cluster boundary: execute() writes the statmech
    artifacts, convergence statuses, and energy settings a real ARC run would leave in the round's
    own project directory, per the class-level ``converge_plan`` (one set of network-local TS
    labels per execute() call)."""

    converge_plan = ()
    constructions = []
    executions = 0

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        _ArtifactWritingFakeARC.constructions.append(kwargs)

    def as_dict(self):
        return dict(self.kwargs)

    def execute(self):
        arc_dir = self.kwargs['project_directory']
        to_converge = _ArtifactWritingFakeARC.converge_plan[_ArtifactWritingFakeARC.executions]
        _ArtifactWritingFakeARC.executions += 1
        _write_energy_settings(arc_dir)
        for reaction in self.kwargs['reactions']:
            label = reaction['ts_label']
            if label.rsplit('_', 1)[-1] in to_converge:
                _write_artifact(expected_ts_artifact_path(arc_dir, label))
                _write_status(arc_dir, label, converged=True)
            else:
                _write_status(arc_dir, label, converged=False)


def _hybrid_qm_ilt(hybrid_path):
    """Which of the fixture's three TSs a hybrid file carries as QM/RRKM vs RMG/ILT: a QM'd TS's
    reaction block drops its ``kinetics = Arrhenius(...)`` line, an ILT one keeps it. Scanned
    textually because the hybrid's consumer is Arkane (which accepts its positional
    ``transitionState(...)`` calls), not ``t3.pdep.parser``."""
    with open(hybrid_path) as f:
        content = f.read()
    qm, ilt = [], []
    for reaction_label, ts_label in (('reaction1', 'TS1'), ('reaction2', 'TS2'),
                                     ('reaction3', 'TS3')):
        start = content.index(f"label = '{reaction_label}'")
        end = content.find('reaction(', start)
        block = content[start:end if end != -1 else len(content)]
        (ilt if 'kinetics = Arrhenius(' in block else qm).append(ts_label)
    return sorted(qm), sorted(ilt)


def _fixture_explorer(monkeypatch, project_directory, network_id):
    """The canonical fake explorer of this module: each round's 'exploration' deposits a copy of
    the real fixture at the exact path arc_qm_runner's real _explored_network_path recomputes."""
    explore_calls = []

    def _fake_explore(*, network_path, config, logger=None):
        explore_calls.append(network_path)
        paths = round_paths(project_directory, len(explore_calls) - 1)
        dest_path = _explored_network_path(paths, network_id)
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        shutil.copyfile(_FIXTURE_NETWORK_PATH, dest_path)
        return PDepExplorationResult(network_id=network_id,
                                     status=EXPLORATION_STATUS_SUCCEEDED,
                                     network_paths=(dest_path,))

    monkeypatch.setattr('t3.pdep.pes_loop.explore_pdep_network', _fake_explore)
    return explore_calls


def test_real_loop_real_capture_carries_cumulative_qm_across_rounds(tmp_path, monkeypatch):
    """The reviewer's exact scenario, against the REAL capture_ts_artifacts: round 0 adopts TS1
    from a real prior capture and converges TS2; round 1 converges TS3. Before this fix round the
    loop could not complete a single such round (defect 1), and when it did complete, round N+1's
    hybrid dropped round N's QM back to Arrhenius lines (defect 3)."""
    prior_project_dir = os.path.join(str(tmp_path), 'prior_project')
    _build_prior_capture(prior_project_dir, ts_labels=('TS1',))
    project_directory = os.path.join(str(tmp_path), 'loop_project')
    os.makedirs(project_directory)
    seed_path = os.path.join(str(tmp_path), 'network0_full.py')
    shutil.copyfile(_FIXTURE_NETWORK_PATH, seed_path)

    _fixture_explorer(monkeypatch, project_directory, 'network0_full')
    _ArtifactWritingFakeARC.converge_plan = ({'TS2'}, {'TS3'})
    _ArtifactWritingFakeARC.constructions = []
    _ArtifactWritingFakeARC.executions = 0
    monkeypatch.setattr(pes_qm, 'ARC', _ArtifactWritingFakeARC)

    config = PESLoopConfig(
        pes={'network': seed_path, 'source': ['O-2(13598)', 'CO2(11)'],
             'bath_gas': {'Ar': 1.0}, 'method': 'MSC'},
        termination={'max_rounds': 4},
        reuse={'from_t3_projects': [prior_project_dir]},
    )
    result = run_pes_loop(config, project_directory=project_directory, qm_runner=arc_qm_runner)

    assert result.status == 'converged', f'{result.status}: {result.reason}'
    assert [record.status for record in result.rounds] == ['continuing', 'continuing', 'converged']
    assert sorted(result.rounds[0].queued_ts_labels) == ['TS2', 'TS3']
    assert sorted(result.rounds[1].queued_ts_labels) == ['TS3']

    round0_hybrid = hybrid_network_path(round_paths(project_directory, 0), 'network0_full')
    round1_hybrid = hybrid_network_path(round_paths(project_directory, 1), 'network0_full')
    # Round 0: the adopted TS1 and the freshly converged TS2 are QM; TS3 is still ILT.
    assert _hybrid_qm_ilt(round0_hybrid) == (['TS1', 'TS2'], ['TS3'])
    # Round 1 -- THE defect-3 pin: the hybrid carries the CUMULATIVE QM, not just this round's.
    # The defective loop wrote QM=['TS3'], ILT=['TS1', 'TS2'] here.
    assert _hybrid_qm_ilt(round1_hybrid) == (['TS1', 'TS2', 'TS3'], [])

    # Both rounds' captures are REAL and re-verify cleanly, sensitivity evidence included --
    # exactly what defect 1 made impossible.
    for round_index in (0, 1):
        verified = verify_capture(round_paths(project_directory, round_index).capture)
        for record in verified.ts_records:
            if record.artifact_path is not None:
                assert record.coefficient is not None and record.delta_ln_k is not None


def test_real_loop_round_0_full_adoption_completes_with_real_vendoring(tmp_path, monkeypatch):
    """Defect 2, end to end against the real vendoring and the real _adopted_energy_settings:
    when EVERY candidate is adopted from a prior run, round 0 queues nothing, never touches ARC or
    capture, and still writes a hybrid carrying all three TSs as QM; round 1 then converges."""
    prior_project_dir = os.path.join(str(tmp_path), 'prior_project')
    _build_prior_capture(prior_project_dir, ts_labels=('TS1', 'TS2', 'TS3'))
    project_directory = os.path.join(str(tmp_path), 'loop_project')
    os.makedirs(project_directory)
    seed_path = os.path.join(str(tmp_path), 'network0_full.py')
    shutil.copyfile(_FIXTURE_NETWORK_PATH, seed_path)

    _fixture_explorer(monkeypatch, project_directory, 'network0_full')
    _ArtifactWritingFakeARC.converge_plan = ()
    _ArtifactWritingFakeARC.constructions = []
    _ArtifactWritingFakeARC.executions = 0
    monkeypatch.setattr(pes_qm, 'ARC', _ArtifactWritingFakeARC)

    config = PESLoopConfig(
        pes={'network': seed_path, 'source': ['O-2(13598)', 'CO2(11)'],
             'bath_gas': {'Ar': 1.0}, 'method': 'MSC'},
        termination={'max_rounds': 2},
        reuse={'from_t3_projects': [prior_project_dir]},
    )
    result = run_pes_loop(config, project_directory=project_directory, qm_runner=arc_qm_runner)

    assert result.status == 'converged', f'{result.status}: {result.reason}'
    assert [record.status for record in result.rounds] == ['continuing', 'converged']
    # Nothing was ever queued, so ARC must never even have been constructed.
    assert _ArtifactWritingFakeARC.constructions == []
    round0_hybrid = hybrid_network_path(round_paths(project_directory, 0), 'network0_full')
    assert _hybrid_qm_ilt(round0_hybrid) == (['TS1', 'TS2', 'TS3'], [])
    # The hybrid's energy settings came from the adopted artifacts' own prior manifest -- the
    # _adopted_energy_settings path a capture-less round is forced through.
    with open(round0_hybrid) as f:
        assert 'wb97xd2023' in f.read()
