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

import os

from t3.pdep.pes_loop import PES_LOOP_DIAGRAM_ONLY, run_pes_loop
from t3.schema import PESLoopConfig

_SOURCE_SEED_PATH = os.path.join(
    os.path.dirname(__file__), '..', 'data', 'pdep_source_seeds', 'cho2_source',
    'network_cho2_source.py')


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
