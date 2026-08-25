"""
t3 tests test_pdep test_round1_regression module

End-to-end regression test that PINS the round-0 -> round-1 handoff of the PES exploration loop.

Why this test exists, and why it is shaped the way it is
--------------------------------------------------------
The handoff -- T3 splices QM'd conformer data onto RMG estimates to write a *hybrid* Arkane
network (``t3.pdep.hybrid.write_hybrid_network_input_file``), then round 1 re-explores that hybrid
through the real Arkane explorer (``t3.pdep.api.explore_pdep_network``) -- was driven from EIGHT
consecutive contact-defects to zero. Every one of those eight was live while the whole test suite
was green, because each lived in the *contact* between a producer and a consumer that no test ever
put together on a real artifact; three were literally "the loop wrote a file its own reader
rejected." A fast fixture test that only inspects the written file cannot catch a consumer-side
defect, and consumer-side is exactly where these bugs live. So this test drives the REAL consumer.

Level, and the trade-off chosen
-------------------------------
This is a single slow (~30 s) end-to-end test. It shells out to Arkane in the ``rmg_env``
environment (Python 3.9) via the same ``run_arkane_job`` path production uses. It sits at the
*loop* level on purpose: the whole point of the ticket is catching what fast, in-process fixture
tests structurally cannot -- a producer/consumer contact defect. A file-only assertion would be
fast but blind to precisely the failure class being protected against.

The expensive regenerate->explore is run ONCE, in a module-scoped fixture, and four separate tests
assert against its result -- so each of the four assertions the ticket asks for shows up as its own
line in the pytest output (completes / rate coefficient / DOF consistency / parser-readable)
without paying the ~30 s cost four times over.

It REGENERATES the hybrid through the production writer rather than replaying a stored one: a test
seeded from a stored artifact can only exercise code *downstream* of where that artifact was
produced, and the stored ``round_0/hybrid/`` in the run directory was written by an older writer
and still carries a since-fixed defect. Regeneration is the only way to exercise the writer/reader
contact these defects lived in.

Gating
------
The test is marked ``slow`` and is SKIPPED (never silently: with an explicit reason) only where
``rmg_env`` is genuinely unavailable. Where ``rmg_env`` exists -- the developer boxes and CI that
run the master equation -- it RUNS as part of the normal ``tests/test_pdep`` suite. CI must keep
``rmg_env`` on the runner (it already needs it for the Arkane/ME tests) so this test executes
rather than skips; ``-m "not slow"`` deselects it for a fast local loop.

Real committed artifacts driven (NOT synthetic fixtures)
--------------------------------------------------------
- ``tests/data/pdep_hybrid/cho2_round0/source_network0_reduced.py`` -- the round-0 RMG network the
  hybrid is built from.
- ``tests/data/pdep_hybrid/cho2_round0/conformers.json`` -- the vibration-only QM conformer data
  (TS2 and two wells) spliced in, as ``t3.runners.statmech_conformer_extract`` produces it.

The energy-reference settings and the round-1 explorer configuration are the exact values the real
``r001_m1-cho2-pilot`` pilot used (its stored hybrid header and its ``input.yml``); they are inputs
to the writer/explorer, not the artifact under test.
"""

import json
import os
import re
import shutil
import subprocess

import pytest

from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.api import explore_pdep_network
from t3.pdep.explorer.config import PDepExplorerConfig
from t3.pdep.explorer.result import EXPLORATION_STATUS_SUCCEEDED
from t3.pdep.hybrid import QMEnergySettings, write_hybrid_network_input_file
# The single authoritative list of statmech modes that carry translation / overall rotation. A
# DOF-consistent p-dep network species must carry NONE of these (see the constant's own comment in
# t3.pdep.hybrid). Imported rather than re-listed so this test tracks the production definition.
from t3.pdep.hybrid import _EXTERNAL_DOF_MODE_NAMES
from t3.pdep.parser import parse_pdep_network_file

CHO2_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_hybrid', 'cho2_round0')
SOURCE_NETWORK = os.path.join(CHO2_DIR, 'source_network0_reduced.py')
CONFORMERS_JSON = os.path.join(CHO2_DIR, 'conformers.json')

# The exact energy-reference settings the r001_m1-cho2-pilot run used (its stored hybrid header:
# ~/runs/t3-pes/r001_m1-cho2-pilot/run/round_0/hybrid/network0_reduced.py). These are INPUTS to the
# writer, not the thing under test.
PILOT_ENERGY_SETTINGS = QMEnergySettings(
    model_chemistry="LevelOfTheory(method='wb97xd2023',basis='def2tzvp',software='gaussian')",
    frequency_scale_factor=0.988,
    atom_energies={
        'Br': -2574.174533595486, 'C': -37.84706210301937, 'Cl': -460.1467876783656,
        'F': -99.73955550924293, 'H': -0.5006557872395249, 'N': -54.584995947182875,
        'O': -75.07252406126821, 'S': -398.1105530401693,
    },
    tunneling='Eckart',
)

# The QM'd transition state and the two wells adopted onto the same energy reference, keyed exactly
# as they appear in the source network / conformers.json.
QM_TS_LABEL = 'TS2'
QM_WELL_LABELS = ('[O]C=O(1)', 'O=[C]O(8)')

# The measured, independently-verified reverse rate coefficient the working handoff produces, as an
# Arrhenius fit Arkane writes into the explorer's output.py for the QM'd reaction. See the tolerance
# discussion on test_rate_coefficient_matches below.
EXPECTED_K_REV_TST = {'A': 8.73653e+11, 'n': 0.308264, 'Ea': 106.537}
EXPECTED_K_REV_TST_TUNNELING = {'A': 0.201422, 'n': 3.90893, 'Ea': 71.7449}

# Tolerance: 1% relative on every fitted parameter. See test_rate_coefficient_matches for why.
K_REV_REL_TOL = 1.0e-2


def _rmg_env_available() -> bool:
    """Whether an ``rmg_env`` conda/micromamba environment exists on this box.

    The explorer shells Arkane out into ``rmg_env`` (see ``t3.runners.rmg_runner.run_arkane_job``);
    without it, round 1 cannot run. Checked by listing envs -- fast, and it does NOT import Arkane
    (which lives in that other environment). A False here produces an explicit skip with a reason,
    never a silent pass.
    """
    for mgr in ('micromamba', 'mamba', 'conda'):
        exe = shutil.which(mgr)
        if not exe:
            continue
        try:
            out = subprocess.run([exe, 'env', 'list'], capture_output=True, text=True, timeout=60)
        except (OSError, subprocess.SubprocessError):
            continue
        for line in out.stdout.splitlines():
            first = line.split()[0] if line.split() else ''
            if first == 'rmg_env' or line.rstrip().endswith('/rmg_env'):
                return True
    return False


def _parse_k_rev(output_py_text: str, which: str) -> dict:
    """Parse the ``# k_rev (<which>) = Arrhenius(A=(..,), n=.., Ea=(..,'kJ/mol'), ...)`` header line
    Arkane writes above the QM'd reaction block, returning ``{'A':, 'n':, 'Ea':}`` (Ea in kJ/mol).
    """
    num = r"([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
    pattern = (rf"# k_rev \({re.escape(which)}\) = Arrhenius\(A=\({num},'[^']*'\), "
               rf"n={num}, Ea=\({num},'kJ/mol'\)")
    m = re.search(pattern, output_py_text)
    assert m is not None, f"no 'k_rev ({which})' Arrhenius line found in output.py"
    return {'A': float(m.group(1)), 'n': float(m.group(2)), 'Ea': float(m.group(3))}


@pytest.fixture(scope='module')
def round1_handoff(tmp_path_factory):
    """Regenerate the real committed hybrid and drive it through the real round-1 explorer, ONCE.

    Returns a dict with the regenerated hybrid path, the explorer's output.py path (carrying the
    rate coefficient), the resolved final network path, and the exploration result.
    """
    if not _rmg_env_available():
        pytest.skip("rmg_env (conda/micromamba) is required to run the real Arkane round-1 "
                    "explorer; this end-to-end regression test cannot run without it. CI that runs "
                    "the master equation must keep rmg_env on the runner so this test executes.")

    root = str(tmp_path_factory.mktemp('cho2_round1'))
    hybrid_path = os.path.join(root, 'hybrid', 'network0_reduced.py')
    with open(CONFORMERS_JSON) as f:
        conformers = json.load(f)

    write_result = write_hybrid_network_input_file(
        source_path=SOURCE_NETWORK,
        dest_path=hybrid_path,
        method='MSC',
        qm_transition_states={QM_TS_LABEL: conformers[QM_TS_LABEL]},
        energy_settings=PILOT_ENERGY_SETTINGS,
        qm_wells={label: conformers[label] for label in QM_WELL_LABELS},
    )

    output_directory = os.path.join(root, 'round_1', 'explorer')
    config = PDepExplorerConfig(
        explorer='Arkane',
        trusted_output_root=root,
        output_directory=output_directory,
        seed_species=('[O]C=O',),
        method='MSC',
        bath_gas={'Ar': 1.0},
        explore_tol=0.01,
        energy_tol=1.0e+01,
        flux_tol=1.0e-06,
        maximum_radical_electrons=2,
        timeout=7200.0,
    )
    result = explore_pdep_network(network_path=hybrid_path, config=config)

    return {
        'write_result': write_result,
        'hybrid_path': hybrid_path,
        'result': result,
        'output_py_path': os.path.join(output_directory, 'output.py'),
    }


@pytest.mark.slow
def test_round1_completes(round1_handoff):
    """Assertion 1 (highest value): round 1 completes -- a final network is produced and the loop
    does not report failure. This alone catches every historical break that failed the explore
    outright (a singular master equation, a file the reader rejected, a bad seed resolution)."""
    result = round1_handoff['result']
    assert result.status == EXPLORATION_STATUS_SUCCEEDED, (
        f"round 1 did not complete: status={result.status!r}, reasons={result.reasons}")
    assert result.network_paths, "a succeeded exploration must carry at least one final network path"
    assert os.path.isfile(result.network_paths[0]), \
        f"final network path does not exist: {result.network_paths[0]!r}"


@pytest.mark.slow
def test_rate_coefficient_matches(round1_handoff):
    """Assertion 2: the reverse rate coefficient is what it should be, within tolerance.

    TOLERANCE CHOICE -- 1% relative on every fitted Arrhenius parameter (A, n, Ea), for both the
    bare-TST and the TST+tunneling fits. Rationale a future reader hitting this failure needs:

    * These six numbers are a deterministic TST + three-parameter Arrhenius least-squares fit; on a
      FIXED toolchain (the same rmg_env Arkane/numpy/scipy) they reproduce bit-for-bit, so 1% is
      enormous headroom for the re-fit itself and the test does not flap.
    * Every one of the eight historical breaks this test protects against moved the answer by ORDERS
      of magnitude (the DOF mismatch made k(E) explode by ~1e14) or produced no answer at all
      (singular matrix / failed explore). A 1% band catches all of those with room to spare and also
      flags an energy-reference or tunneling regression that shifts Ea by more than ~1 kJ/mol.
    * Therefore: a mismatch WITHIN 1% means "the numbers legitimately moved" -- most likely an
      Arkane/numpy/scipy version bump changed the fit; re-measure and update the constants. A
      mismatch BEYOND 1%, or a failed explore (assertion 1), means "you broke the handoff."
    """
    with open(round1_handoff['output_py_path']) as f:
        text = f.read()

    tst = _parse_k_rev(text, 'TST')
    for key, expected in EXPECTED_K_REV_TST.items():
        assert tst[key] == pytest.approx(expected, rel=K_REV_REL_TOL), \
            f"k_rev (TST) {key}={tst[key]!r} differs from expected {expected!r} by more than {K_REV_REL_TOL:.0%}"

    tst_t = _parse_k_rev(text, 'TST+T')
    for key, expected in EXPECTED_K_REV_TST_TUNNELING.items():
        assert tst_t[key] == pytest.approx(expected, rel=K_REV_REL_TOL), \
            f"k_rev (TST+T) {key}={tst_t[key]!r} differs from expected {expected!r} by more than {K_REV_REL_TOL:.0%}"


@pytest.mark.slow
def test_generated_hybrid_is_dof_consistent(round1_handoff):
    """Assertion 3: no isomer or transition state in the generated hybrid carries a translational or
    overall-rotational mode. This is the invariant whose violation -- a QM conformer spliced in with
    its full IdealGasTranslation/NonlinearRotor modes onto vibration-only estimated wells -- made
    Q_TS/Q_reactant fail to cancel and the modified-strong-collision matrix go singular.

    Checked here by an INDEPENDENT AST scan of the produced file (reusing only the production list of
    external-DOF mode names), not by trusting the writer's own outbound self-check -- so a mutation
    that both restores the modes AND disables that self-check is still caught."""
    import ast

    with open(round1_handoff['hybrid_path']) as f:
        text = f.read()
    # Scope the scan the way the production invariant scopes itself: isomers and transition states
    # only. A bimolecular reactant/product channel species or the bath gas is declared with
    # species(...) too, but is NOT a well -- it never enters the isomer<->isomer RRKM density of
    # states and so legitimately carries a full conformer. Checking every species() block would
    # make a perfectly valid network fail this test, which is a false alarm on the one assertion
    # that exists to catch a real and subtle defect.
    isomer_labels = set(parse_pdep_network_file(round1_handoff['hybrid_path']).isomers)
    offenders = []
    for node in ast.parse(text).body:
        if not (isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)):
            continue
        call = node.value
        name = getattr(call.func, 'id', None)
        if name not in ('species', 'transitionState'):
            continue
        keywords = {kw.arg: kw.value for kw in call.keywords if kw.arg is not None}
        modes_node = keywords.get('modes')
        if not isinstance(modes_node, (ast.List, ast.Tuple)):
            continue
        mode_names = {getattr(el.func, 'id', None) for el in modes_node.elts if isinstance(el, ast.Call)}
        label_node = keywords.get('label')
        label = label_node.value if isinstance(label_node, ast.Constant) else '<?>'
        if name == 'species' and label not in isomer_labels:
            continue  # a channel species or the bath gas: not a well, not part of the invariant.
        bad = sorted(mode_names & set(_EXTERNAL_DOF_MODE_NAMES))
        if bad:
            offenders.append(f"{name} {label!r}: {bad}")
    assert not offenders, (
        "the generated hybrid carries translational/rotational modes on isomers/transition states "
        f"(the singular-master-equation invariant is violated): {offenders}")


@pytest.mark.slow
def test_final_network_is_parseable_by_the_loop(round1_handoff):
    """Assertion 4: the generated final network is readable by the loop's OWN parser. Three of the
    eight historical defects were exactly this -- the explorer wrote a network file that
    ``parse_pdep_network_file`` then rejected."""
    result = round1_handoff['result']
    # Same guards as test_round1_completes, repeated rather than assumed: when the explore fails,
    # network_paths is empty, and indexing it here would raise IndexError instead of reporting the
    # actual round-1 failure -- the more confusing the more narrowly the test is selected.
    assert result.status == EXPLORATION_STATUS_SUCCEEDED, (
        f"round 1 did not complete, so there is no final network to parse: status={result.status!r}, "
        f"reasons={result.reasons}")
    assert result.network_paths, "a succeeded exploration must carry at least one final network path"
    final_network = result.network_paths[0]
    network = parse_pdep_network_file(final_network)
    # A parsed round-1 network must carry real content: species and at least one path reaction (the
    # QM'd channel). An empty parse would be a reader that "succeeds" on nothing.
    assert network.species_labels, f"parsed final network {final_network!r} carries no species"
    assert network.path_reactions, f"parsed final network {final_network!r} carries no reactions"
