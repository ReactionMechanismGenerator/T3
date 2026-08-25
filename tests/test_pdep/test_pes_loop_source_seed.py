"""
End-to-end PES-loop test from a SOURCE-ONLY seed, driving the REAL Arkane explorer.

This is the test the I-009 Verifier asks for: no double at the explorer seam. It drives the real
``t3.pdep.pes_loop.run_pes_loop`` with ``qm_runner=None`` (the ``PES.py --diagram-only`` path) on a
source-only CHO2 seed -- one well, a bath gas, no hand-written ``reaction()`` -- and asserts the
loop runs Arkane's own explorer, accepts the genuinely-explored network, and writes its own PES
diagram. It exercises, through the real callee, four of the defects this ticket cleared: the seed
parser accepting a reaction-less seed, the ``explorer='Arkane'`` case fix, the ``_check_final_network_payload``
label normalization, and the scipy/RMG ``LinAlgWarning`` stderr false-positive. A fully-mocked seam
hid every one of them, which is exactly why this test refuses a double.

It runs Arkane under ``rmg_env`` via ``run_arkane_job`` (same as ``test_pes_sa.py``'s real-Arkane
test) and costs ~30 s; it is the single end-to-end explorer test in the suite.
"""

import json
import os

from t3.pdep.hybrid import QMEnergySettings, write_hybrid_network_input_file
from t3.pdep.pes_loop import PES_LOOP_DIAGRAM_ONLY, PES_LOOP_MAX_ROUNDS, run_pes_loop
from t3.pdep.pes_qm import _explored_network_path
from t3.pdep.pes_rounds import hybrid_network_path
from t3.schema import PESLoopConfig

_SOURCE_SEED_PATH = os.path.join(
    os.path.dirname(__file__), '..', 'data', 'pdep_source_seeds', 'cho2_source',
    'network_cho2_source.py')

# The REAL committed CHO2 round-0 conformer data (Arkane-extracted E0/frequencies/imaginary for
# both wells and the QM'd transition state, on one self-consistent energy reference), used to
# replay the QM the pilot run computed on the cluster without submitting anything.
_CHO2_CONFORMERS_PATH = os.path.join(
    os.path.dirname(__file__), '..', 'data', 'pdep_hybrid', 'cho2_round0', 'conformers.json')
_CHO2_ENERGY_SETTINGS = QMEnergySettings(
    model_chemistry="LevelOfTheory(method='wb97xd2023',basis='def2tzvp',software='gaussian')",
    atom_energies={
        'Br': -2574.174533595486, 'C': -37.84706210301937, 'Cl': -460.1467876783656,
        'F': -99.73955550924293, 'H': -0.5006557872395249, 'N': -54.584995947182875,
        'O': -75.07252406126821, 'S': -398.1105530401693,
    },
    tunneling='Eckart',
)
_CHO2_QM_TS = 'TS2'
_CHO2_QM_WELLS = ('O=[C]O(8)', '[O]C=O(1)')


def _source_only_config(network_path):
    """A PES-loop config seeded from the source-only CHO2 well, mirroring examples/pes_loop_cho2."""
    return PESLoopConfig(
        pes={'network': network_path,
             'source': ['[O]C=O'],
             'method': 'MSC',
             'bath_gas': {'Ar': 1.0},
             'explore_tol': 0.01,
             'energy_tol': 1.0e1,
             'flux_tol': 1.0e-6,
             'maximum_radical_electrons': 2,
             'timeout': 1200.0},
        termination={'max_rounds': 1, 'stop_when_no_new_ts': True},
    )


def test_diagram_only_runs_end_to_end_from_a_source_only_seed(tmp_path):
    """The real loop, real explorer, no double: a source-only seed explores to a genuine network,
    the payload check accepts it, and the loop writes its own diagram."""
    result = run_pes_loop(_source_only_config(os.path.abspath(_SOURCE_SEED_PATH)),
                          str(tmp_path), qm_runner=None)

    assert result.status == PES_LOOP_DIAGRAM_ONLY, (result.status, result.reason)
    assert len(result.rounds) == 1

    # The loop drew its OWN diagram (t3.pdep.pes_loop._draw_round_diagram), not merely the adapter's.
    assert result.final_diagram_path is not None, result.reason
    assert os.path.isfile(result.final_diagram_path)
    assert os.path.getsize(result.final_diagram_path) > 0

    # The accepted final network is the real explored artifact Arkane wrote.
    assert result.final_network_path is not None
    assert os.path.isfile(result.final_network_path)
    assert os.path.getsize(result.final_network_path) > 0


def _replay_qm_runner(candidates, paths, config, network_id, adopted=None):
    """A ``qm_runner`` that replays the pilot's already-computed QM instead of submitting to the
    cluster: it folds the REAL committed CHO2 conformers into a hybrid network through the REAL
    ``write_hybrid_network_input_file`` -- exactly what ``t3.pdep.pes_qm.arc_qm_runner`` does after
    ARC converges, minus the (gated) cluster submission and capture. The hybrid's source is this
    round's own explored network, so the labels it carries are the ones the loop itself produced.
    """
    with open(_CHO2_CONFORMERS_PATH) as f:
        conf = json.load(f)
    write_hybrid_network_input_file(
        source_path=_explored_network_path(paths, network_id),
        dest_path=hybrid_network_path(paths, network_id),
        method='MSC',
        qm_transition_states={_CHO2_QM_TS: conf[_CHO2_QM_TS]},
        qm_wells={well: conf[well] for well in _CHO2_QM_WELLS},
        energy_settings=_CHO2_ENERGY_SETTINGS,
    )
    return frozenset({_CHO2_QM_TS}), (_CHO2_QM_TS,)


def test_round_one_reexploration_of_a_real_hybrid_completes(tmp_path):
    """The loop's OWN round-0 -> round-1 path, through the REAL explorer and the REAL committed
    hybrid data, must complete -- not trip the output-vs-network label check (I-025).

    Round 0 explores the source-only CHO2 seed (real Arkane); a replay ``qm_runner`` folds the real
    committed QM conformers into a hybrid through the real ``write_hybrid_network_input_file``; the
    loop then re-explores that hybrid (real Arkane) and runs the real
    ``_check_final_network_payload``. The hybrid's species labels already carry one RMG index
    (``O=C=O(5)``) and Arkane appends a second when it writes the round's network file
    (``O=C=O(5)(3)``) while ``output.py`` keeps the single-indexed label; a normalization that
    strips only one index level refuses this genuinely single-run pairing. This test is red before
    the fix (final network stays round 0's pre-QM surface and the reason carries the mismatch) and
    green after it. Only the cluster submission is replayed -- every explore, the ME sensitivity
    analysis, the hybrid write and the check run unfaked, against real data.
    """
    result = run_pes_loop(_source_only_config(os.path.abspath(_SOURCE_SEED_PATH)),
                          str(tmp_path), qm_runner=_replay_qm_runner)

    # The guard must NOT refuse this pairing: the two artifacts describe the same network.
    assert 'do not describe the same network' not in (result.reason or ''), result.reason
    assert result.status == PES_LOOP_MAX_ROUNDS, (result.status, result.reason)

    # The reported final network is the QM-informed round-1 re-exploration of the hybrid (under
    # round_1/), not round 0's pre-QM surface -- i.e. the bonus re-exploration was ACCEPTED.
    assert result.final_network_path is not None, result.reason
    assert os.path.isfile(result.final_network_path)
    assert os.path.join('round_1', 'explorer') in result.final_network_path, result.final_network_path
    assert 'QM-informed re-exploration' in (result.reason or ''), result.reason
