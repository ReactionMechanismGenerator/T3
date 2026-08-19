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

from t3.pdep.capture import CaptureResult
from t3.pdep.discovery import ARTIFACT_STATUS_USABLE, TSArtifactRecord
from t3.pdep.explorer.result import EXPLORATION_STATUS_SUCCEEDED, PDepExplorationResult
from t3.pdep.hybrid import HybridNetworkResult
from t3.pdep.parser import parse_pdep_network_file
import t3.pdep.pes_qm as pes_qm
from t3.pdep.pes_loop import run_pes_loop
from t3.pdep.pes_qm import _explored_network_path, arc_qm_runner
from t3.pdep.pes_rounds import round_paths, split_qm_candidates
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


def _config(max_rounds=3):
    return PESLoopConfig(pes={'network': _FIXTURE_NETWORK_PATH, 'source': ['A'],
                              'bath_gas': {'He': 1.0}},
                         termination={'max_rounds': max_rounds, 'stop_when_no_new_ts': False})


def test_real_run_pes_loop_wires_the_real_arc_qm_runner_across_rounds(tmp_path, monkeypatch):
    """Drive the real run_pes_loop with the real arc_qm_runner injected as qm_runner, faking only
    ARC, capture_ts_artifacts, and write_hybrid_network_input_file. Assert the loop completes at
    least two rounds, and -- the actual N6 contract -- that the exact hybrid file arc_qm_runner
    writes for round 0 is the exact file run_pes_loop hands its round-1 explorer, observed end to
    end through the real loop rather than inferred from two separately-passing unit tests."""
    real_network = parse_pdep_network_file(path=_FIXTURE_NETWORK_PATH)
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
        shutil.copyfile(_FIXTURE_NETWORK_PATH, dest_path)
        return PDepExplorationResult(network_id=real_network.network_id,
                                     status=EXPLORATION_STATUS_SUCCEEDED,
                                     network_paths=(dest_path,))

    monkeypatch.setattr('t3.pdep.pes_loop.explore_pdep_network', _fake_explore)

    _FakeARC.constructions = []
    _FakeARC.execute_calls = 0
    monkeypatch.setattr(pes_qm, 'ARC', _FakeARC)

    write_hybrid_calls = []

    def _fake_capture_ts_artifacts(*, join_records, arc_project_directory, capture_dir, networks):
        os.makedirs(capture_dir, exist_ok=True)
        return CaptureResult(
            capture_dir=capture_dir,
            manifest_path=os.path.join(capture_dir, 'manifest.yml'),
            records=(
                TSArtifactRecord(network_id=real_network.network_id,
                                 network_ts_label=target.ts_label,
                                 arc_ts_label=f'{real_network.network_id}_{target.ts_label}',
                                 status=ARTIFACT_STATUS_USABLE,
                                 artifact_path=os.path.join(capture_dir, f'{target.ts_label}.py'),
                                 converged=True),
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

    config = _config(max_rounds=3)
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
