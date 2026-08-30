"""
Tests for ``t3.pdep.explorer.seed`` -- generating a source-only PES seed from a species identifier.

The fast tests here drive the REAL entry point for the negative control (a bimolecular identifier
is refused before any subprocess spawns -- I-038 verifier 3) and the REAL assembly/templating for
the ``network()``/``pressureDependence()`` boilerplate. The heavy round-trip and non-CHO2 explorer
verifiers (1 and 2) drive the real RMG estimator subprocess and a real Arkane explorer round; they
are gated on ``rmg_env`` being importable and marked slow, mirroring ``test_pes_loop_source_seed``.
"""
import json
import math
import os
import re
import shutil
import subprocess

import pytest

from t3.pdep.explorer.input_file import write_arkane_explorer_input_file
from t3.pdep.explorer import seed as seed_module
from t3.pdep.explorer.seed import (
    BimolecularSourceError,
    _assemble_seed_file,
    _render_pdep_block,
    _run_driver,
    _DEFAULT_PDEP_SETTINGS,
    generate_source_seed,
)


def _rmg_env_available(env_name=None) -> bool:
    """True only if the ``rmg_env`` conda environment the estimator subprocess needs actually EXISTS.

    Checking merely for a conda/mamba launcher on ``PATH`` (or ``MAMBA_EXE``) is a false positive:
    the manager can be present while the env is not, and the slow tests then *fail* (the subprocess
    exits non-zero) instead of *skip*. Resolve the env name the way
    ``arc.job.env_run.rmg_env_command`` does, then ask a launcher to list its environments by name
    and look for it. Any failure -> treat the env as absent, so a machine that simply lacks it
    fail-SKIPS rather than fail-errors.
    """
    if env_name is None:
        try:
            from arc.imports import settings
            env_name = settings.get('RMG_ENV_NAME', 'rmg_env')
        except Exception:
            env_name = 'rmg_env'
    for launcher in _rmg_launcher_candidates():
        try:
            proc = subprocess.run([launcher, 'env', 'list', '--json'],
                                  capture_output=True, text=True, timeout=60)
        except Exception:
            continue
        if proc.returncode != 0:
            continue
        try:
            envs = json.loads(proc.stdout).get('envs', [])
        except (ValueError, AttributeError):
            continue
        if any(os.path.basename(str(p).rstrip('/')) == env_name for p in envs):
            return True
    return False


def _rmg_launcher_candidates():
    """The conda-family launchers to try, mirroring ``env_run._detect_launcher``'s preference order:
    the active one (``CONDA_EXE`` / ``MAMBA_EXE``) first, then conda -> mamba -> micromamba on PATH."""
    candidates = []
    for env_var in ('CONDA_EXE', 'MAMBA_EXE'):
        path = os.environ.get(env_var)
        if path and os.path.isfile(path) and path not in candidates:
            candidates.append(path)
    for name in ('conda', 'mamba', 'micromamba'):
        found = shutil.which(name)
        if found and found not in candidates:
            candidates.append(found)
    return candidates


_HAS_RMG_ENV = _rmg_env_available()

_CHO2_HAND_BUILT = os.path.join(
    os.path.dirname(__file__), '..', 'data', 'pdep_source_seeds', 'cho2_source',
    'network_cho2_source.py')


def _code_lines(text: str) -> list:
    """The scientific/structural lines of a seed file: every non-blank line that is not a ``#``
    comment. Comments are excluded so the generated file's header (which differs by design) does not
    count as a difference; every ``species()``/``network()``/``pressureDependence()`` field line and
    every adjacency-list line is kept and compared verbatim, indentation included."""
    return [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith('#')]


# --- Split invariant for the fitted thermo line (I-038) ---
#
# Seed generation is NOT bit-reproducible across environments: the Wilhoit->NASA least-squares
# fit converges slightly differently under a different BLAS/LAPACK/scipy build (~1e-5 relative on
# the coefficients; the fit breakpoint moves ~0.01 K out of ~994 K). Byte-identity was therefore
# never the right invariant for the *fitted* quantity. Everything else -- structure, E0, all
# frequencies, TransportData, SingleExponentialDown, and the group-additivity comment -- IS a real
# input that agrees exactly across environments, so those stay byte-exact and a change is a
# regression. Only the NASA polynomial coefficients and the fit breakpoint are compared numerically.

_NASA_POLY_RE = re.compile(
    r"NASAPolynomial\(coeffs=\[([^\]]*)\], Tmin=\(([^,]+),'K'\), Tmax=\(([^,]+),'K'\)\)")
_OUTER_BOUNDS_RE = re.compile(r"\], Tmin=\(([^,]+),'K'\), Tmax=\(([^,]+),'K'\), E0=")
_E0_RE = re.compile(r"E0=\(([^,]+),'kJ/mol'\)")
_CP0_RE = re.compile(r"Cp0=\(([^,]+),'J/\(mol\*K\)'\)")
_CPINF_RE = re.compile(r"CpInf=\(([^,]+),'J/\(mol\*K\)'\)")
_COMMENT_RE = re.compile(r'comment="""(.*)"""')

# Tolerances, justified from the observed cross-environment spread (largest coefficient relative
# difference ~5.8e-5 on coeffs[1]; breakpoint moved 0.011 K). A modest ~9x margin: tight enough that
# a real thermo change -- different groups, a shifted E0, a swapped frequency -- moves coefficients
# far more than this, loose enough to absorb pure BLAS-level fit noise.
_COEFF_RTOL = 5.0e-4
_BREAK_ATOL_K = 0.1


def _parse_source_thermo(line: str) -> dict:
    """Parse the source species' ``thermo = NASA(...)`` line into fitted (numeric) and non-fitted
    (exact) parts. The two ``NASAPolynomial`` coefficient sets and the fit breakpoint (poly-1 Tmax,
    which equals poly-2 Tmin) are the fit outputs; everything else is an exact input."""
    polys = _NASA_POLY_RE.findall(line)
    assert len(polys) == 2, f"expected two NASA polynomials, found {len(polys)} in {line!r}"
    coeffs = [[float(x) for x in cstr.split(',')] for cstr, _tmin, _tmax in polys]
    break_temps = [float(polys[0][2]), float(polys[1][1])]  # poly-1 Tmax and poly-2 Tmin: the fit knot
    fixed_poly_bounds = (polys[0][1], polys[1][2])           # poly-1 Tmin (100) & poly-2 Tmax (5000): fixed
    outer = _OUTER_BOUNDS_RE.search(line).groups()           # NASA outer Tmin/Tmax (100/5000): fixed
    return dict(
        coeffs=coeffs, break_temps=break_temps, fixed_poly_bounds=fixed_poly_bounds, outer=outer,
        e0=_E0_RE.search(line).group(1), cp0=_CP0_RE.search(line).group(1),
        cpinf=_CPINF_RE.search(line).group(1), comment=_COMMENT_RE.search(line).group(1))


def _assert_source_thermo(hand_line: str, gen_line: str) -> None:
    """The fitted thermo line: non-fit fields byte-exact, coefficients within ``_COEFF_RTOL``, the
    fit breakpoint within ``_BREAK_ATOL_K``."""
    h, g = _parse_source_thermo(hand_line), _parse_source_thermo(gen_line)
    # Non-fitted -> a change here is a real regression, never absorbed into the tolerance.
    assert h['comment'] == g['comment'], f"group-additivity comment changed: {h['comment']!r} vs {g['comment']!r}"
    assert h['e0'] == g['e0'], f"E0 inside NASA changed: {h['e0']!r} vs {g['e0']!r}"
    assert h['cp0'] == g['cp0'], f"Cp0 changed: {h['cp0']!r} vs {g['cp0']!r}"
    assert h['cpinf'] == g['cpinf'], f"CpInf changed: {h['cpinf']!r} vs {g['cpinf']!r}"
    assert h['outer'] == g['outer'], f"outer NASA T bounds changed: {h['outer']} vs {g['outer']}"
    assert h['fixed_poly_bounds'] == g['fixed_poly_bounds'], (
        f"fixed polynomial T bounds changed: {h['fixed_poly_bounds']} vs {g['fixed_poly_bounds']}")
    # Fitted -> numeric.
    for p, (hc, gc) in enumerate(zip(h['coeffs'], g['coeffs'])):
        assert len(hc) == len(gc) == 7, f"NASA poly{p} does not have 7 coeffs: {hc} vs {gc}"
        for k, (a, b) in enumerate(zip(hc, gc)):
            assert math.isclose(b, a, rel_tol=_COEFF_RTOL, abs_tol=0.0), (
                f"NASA poly{p} coeff[{k}] differs beyond rtol={_COEFF_RTOL:g}: "
                f"hand={a!r} gen={b!r} (rel={abs(b - a) / abs(a):.2e})")
    for a, b in zip(h['break_temps'], g['break_temps']):
        assert abs(b - a) <= _BREAK_ATOL_K, (
            f"fit breakpoint differs beyond {_BREAK_ATOL_K} K: hand={a} gen={b} (|d|={abs(b - a):.3g} K)")


def _assert_seed_matches(hand_text: str, gen_text: str) -> None:
    """Compare a generated seed against the hand-built reference under the split invariant: every
    non-fitted code line byte-exact (barring the two cosmetic ``label`` lines), the fitted source
    thermo line numeric. Raises ``AssertionError`` on any real difference, which is what the
    field-for-field test and the mutation demonstration both drive."""
    hand, gen = _code_lines(hand_text), _code_lines(gen_text)
    assert len(hand) == 57, f"hand-built seed changed: {len(hand)} code lines, expected 57"
    assert len(gen) == 57, f"generated seed has {len(gen)} code lines, expected 57"

    # The one fitted line is the source well's group-additivity thermo (NOT the Ar library thermo).
    src_idx = [i for i, ln in enumerate(hand)
               if ln.lstrip().startswith('thermo = NASA') and 'group additivity' in ln.lower()]
    assert len(src_idx) == 1, f"expected exactly one group-additivity thermo line, found {len(src_idx)}"
    (si,) = src_idx
    assert gen[si].lstrip().startswith('thermo = NASA') and 'group additivity' in gen[si].lower(), (
        f"generated line {si} is not the source thermo line: {gen[si]!r}")

    # Every other code line is byte-exact except the two label lines (network + pressureDependence).
    exact_diffs = [(i, h, g) for i, (h, g) in enumerate(zip(hand, gen)) if i != si and h != g]
    assert len(exact_diffs) == 2, (
        "expected exactly two differing non-thermo code lines (the network and pressureDependence "
        f"labels); got {len(exact_diffs)}:\n"
        + '\n'.join(f'  line {i}: hand={h!r} gen={g!r}' for i, h, g in exact_diffs))
    for _i, hand_line, gen_line in exact_diffs:
        assert hand_line == "    label = 'PDepNetwork CHO2 source',", \
            f"unexpected differing hand-built line: {hand_line!r}"
        assert gen_line == "    label = 'PDepNetwork [O]C=O source',", \
            f"unexpected differing generated line: {gen_line!r}"

    # The fitted thermo line: exact non-fit fields, numeric coeffs + breakpoint.
    _assert_source_thermo(hand[si], gen[si])


# --- Negative control (verifier 3): a bimolecular source is refused at the entry point ---

def test_bimolecular_dotted_smiles_is_refused(tmp_path):
    """A dotted (disconnected) SMILES is two molecules -- refused, naming the disjoint-network reason."""
    dest = str(tmp_path / 'never-written.py')
    with pytest.raises(BimolecularSourceError) as exc:
        generate_source_seed('[H].O=C=O', dest)
    msg = str(exc.value)
    assert 'disjoint' in msg
    assert 'two' in msg  # names the two-network consequence
    assert not os.path.isfile(dest)


def test_bimolecular_refusal_happens_before_any_subprocess(tmp_path):
    """The refusal must not depend on the estimator subprocess: a fake driver path would still never
    be reached, because the check precedes it."""
    # A three-fragment identifier is still refused (the rule is "exactly one connected molecule").
    dest = str(tmp_path / 'never-written-2.py')
    with pytest.raises(BimolecularSourceError):
        generate_source_seed('[H].[H].O', dest)
    assert not os.path.isfile(dest)


def test_unparseable_identifier_raises_valueerror(tmp_path):
    dest = str(tmp_path / 'never-written-3.py')
    with pytest.raises(ValueError):
        generate_source_seed('not a smiles !!!', dest)
    assert not os.path.isfile(dest)


def test_unimolecular_source_passes_the_gate():
    """A single connected molecule passes the check (the gate does not over-refuse)."""
    from t3.pdep.explorer.seed import _assert_unimolecular
    mol = _assert_unimolecular('[O]C=O')
    assert len(mol.split()) == 1


# --- Assembly / templating (fast, no subprocess) ---

_FAKE_SOURCE_BLOCK = "species(\n    label = '[O]C=O',\n    thermo = 'stub',\n)"
_FAKE_BATH_BLOCK = "species(\n    label = 'Ar',\n    thermo = 'stub',\n)"


def test_assembled_file_has_no_reaction_or_transition_state():
    """One source well is the contract: the assembled seed must carry NO reaction()/transitionState()
    call blocks (the header comment naming them is fine -- check the AST, not the raw text)."""
    import ast
    text = _assemble_seed_file(_FAKE_SOURCE_BLOCK, _FAKE_BATH_BLOCK, '[O]C=O', 'Ar', 'CHO2')
    called = {node.func.id for node in ast.walk(ast.parse(text))
              if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)}
    assert 'reaction' not in called
    assert 'transitionState' not in called
    assert 'species' in called and 'network' in called and 'pressureDependence' in called


def test_assembled_network_names_the_source_as_sole_isomer():
    text = _assemble_seed_file(_FAKE_SOURCE_BLOCK, _FAKE_BATH_BLOCK, '[O]C=O', 'Ar', 'CHO2')
    assert "label = 'PDepNetwork CHO2 source'," in text
    assert "isomers = [\n        '[O]C=O',\n    ]," in text
    assert "bathGas = {\n        'Ar': 1.0," in text
    # Both species blocks are spliced in verbatim.
    assert _FAKE_SOURCE_BLOCK in text
    assert _FAKE_BATH_BLOCK in text


def test_pdep_block_reproduces_rmg_comma_style():
    """RMG's Quantity repr has no space after the comma inside a quantity tuple; the template must
    match so a generated file diffs cleanly against an RMG-written one."""
    block = _render_pdep_block('PDepNetwork CHO2 source', dict(_DEFAULT_PDEP_SETTINGS))
    assert "Tmin = (700,'K')," in block
    assert "Pmax = (100,'bar')," in block
    assert "method = 'modified strong collision'," in block
    assert "interpolationModel = ('pdeparrhenius',)," in block
    assert "rmgmode = True," in block


def test_generated_seed_parses_as_python():
    """The assembled text must be valid Python source (Arkane execs it)."""
    import ast
    text = _assemble_seed_file(_FAKE_SOURCE_BLOCK, _FAKE_BATH_BLOCK, '[O]C=O', 'Ar', 'CHO2')
    ast.parse(text)


# --- RuntimeError contract (verifier 2): every driver failure surfaces diagnosably, not raw ---

def _fake_completed(stdout='', stderr='', returncode=0):
    return subprocess.CompletedProcess(args=['bash'], returncode=returncode, stdout=stdout, stderr=stderr)


def test_driver_timeout_surfaces_as_runtimeerror(monkeypatch):
    """A subprocess timeout must not escape as ``subprocess.TimeoutExpired``: it is re-raised as a
    ``RuntimeError`` naming the timeout and carrying whatever stdout/stderr was captured."""
    monkeypatch.setattr(seed_module, 'rmg_env_command', lambda **kwargs: 'true')

    def _raise_timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(cmd='bash', timeout=kwargs.get('timeout', 0),
                                        output='partial stdout', stderr='partial stderr from rmg')

    monkeypatch.setattr(seed_module.subprocess, 'run', _raise_timeout)
    with pytest.raises(RuntimeError) as exc:
        _run_driver({'source': {}, 'bath_gas': {}}, timeout=0.5)
    msg = str(exc.value)
    assert 'timed out' in msg
    assert 'partial stderr from rmg' in msg  # actionable context is carried through


def test_driver_malformed_json_surfaces_as_runtimeerror(monkeypatch):
    """A framed payload that is not valid JSON must not escape as ``json.JSONDecodeError``."""
    monkeypatch.setattr(seed_module, 'rmg_env_command', lambda **kwargs: 'true')
    framed = ('rmg log noise\n' + seed_module._BEGIN_MARKER
              + '\n{ this is not valid json ,,, \n' + seed_module._END_MARKER + '\n')
    monkeypatch.setattr(seed_module.subprocess, 'run',
                        lambda *a, **k: _fake_completed(stdout=framed, stderr='a warning line'))
    with pytest.raises(RuntimeError) as exc:
        _run_driver({'source': {}, 'bath_gas': {}}, timeout=1.0)
    msg = str(exc.value)
    assert 'not valid JSON' in msg
    assert 'a warning line' in msg


def test_driver_missing_species_blocks_key_surfaces_as_runtimeerror(monkeypatch):
    """A payload without the ``species_blocks`` key must not escape as ``KeyError``."""
    monkeypatch.setattr(seed_module, 'rmg_env_command', lambda **kwargs: 'true')
    framed = (seed_module._BEGIN_MARKER + '\n' + json.dumps({'wrong_key': []}) + '\n'
              + seed_module._END_MARKER + '\n')
    monkeypatch.setattr(seed_module.subprocess, 'run',
                        lambda *a, **k: _fake_completed(stdout=framed, stderr=''))
    with pytest.raises(RuntimeError) as exc:
        _run_driver({'source': {}, 'bath_gas': {}}, timeout=1.0)
    assert 'species_blocks' in str(exc.value)


def test_driver_wrong_length_species_blocks_surfaces_as_runtimeerror(monkeypatch):
    """A payload whose ``species_blocks`` is not exactly two entries must be validated in the driver
    and re-raised as ``RuntimeError`` -- NOT reach the caller as a confusing tuple-unpack error."""
    monkeypatch.setattr(seed_module, 'rmg_env_command', lambda **kwargs: 'true')
    framed = (seed_module._BEGIN_MARKER + '\n' + json.dumps({'species_blocks': ['only one block']})
              + '\n' + seed_module._END_MARKER + '\n')
    monkeypatch.setattr(seed_module.subprocess, 'run',
                        lambda *a, **k: _fake_completed(stdout=framed, stderr=''))
    with pytest.raises(RuntimeError) as exc:
        _run_driver({'source': {}, 'bath_gas': {}}, timeout=1.0)
    assert 'exactly two' in str(exc.value)


# --- Round-trip (verifier 1): a generated seed is accepted unmodified by the explorer writer ---

@pytest.mark.skipif(not _HAS_RMG_ENV, reason="needs rmg_env for the RMG estimator subprocess")
@pytest.mark.slow
def test_generated_cho2_seed_matches_the_hand_built_seed_field_for_field(tmp_path):
    """Regenerate the CHO2 seed from the bare identifier ``[O]C=O`` and compare it field-for-field
    against ``_CHO2_HAND_BUILT`` -- the hand-assembled seed that produced a successful end-to-end run
    (real surface, real QM, converged). This is the strongest available evidence the generator is
    correct, and it had only ever been done by hand and recorded in prose (I-038).

    Split invariant (I-038): the two files must match on EVERY scientific/structural code line --
    ``structure``, ``E0`` (standalone and inside ``NASA``), all six frequencies, ``spinMultiplicity``,
    ``opticalIsomers``, ``molecularWeight``, ``TransportData``, ``SingleExponentialDown``, the NASA
    ``Cp0``/``CpInf`` and the group-additivity ``comment`` -- BYTE-FOR-BYTE, differing ONLY on the two
    ``label`` lines (the ``network()``/``pressureDependence()`` labels: the hand-built file names the
    network 'CHO2', a bare-identifier regeneration names it after the identifier). A change in any of
    those exact fields (a shifted frequency, a moved ``E0``, different groups) is a real regression and
    fails this test.

    The ONLY line compared numerically is the source well's fitted ``thermo = NASA(...)`` polynomial:
    seed generation is not bit-reproducible across BLAS/LAPACK builds, so the least-squares
    coefficients drift ~1e-5 relative and the fit breakpoint ~0.01 K. Coefficients are compared at
    ``rtol=_COEFF_RTOL`` (5e-4) and the breakpoint at ``atol=_BREAK_ATOL_K`` (0.1 K) -- ~9x the
    observed cross-environment spread, tight enough that a real thermo change fails. This is a split
    of the invariant into exact-vs-numerical parts, NOT a blanket relaxation."""
    seed_path = str(tmp_path / 'network_source.py')
    # No network_name is passed: the label is derived from the identifier, which is the ONLY
    # tolerated non-fit difference against the hand-built seed (whose network is named 'CHO2').
    generate_source_seed('[O]C=O', seed_path, timeout=1800.0)

    with open(seed_path) as f:
        generated_text = f.read()
    with open(_CHO2_HAND_BUILT) as f:
        hand_built_text = f.read()

    _assert_seed_matches(hand_built_text, generated_text)


@pytest.mark.skipif(not _HAS_RMG_ENV, reason="needs rmg_env for the RMG estimator subprocess")
@pytest.mark.slow
def test_generated_cho2_seed_is_accepted_by_explorer_writer(tmp_path):
    """Generate the CHO2 seed from its identifier alone, then feed it UNMODIFIED to
    write_arkane_explorer_input_file. Acceptance = it parses, validates and writes a valid Arkane
    explorer input (verifier 1)."""
    seed_path = str(tmp_path / 'network_source.py')
    generate_source_seed('[O]C=O', seed_path, network_name='CHO2', timeout=1800.0)
    assert os.path.isfile(seed_path)

    dest = str(tmp_path / 'explorer_input.py')
    summary = write_arkane_explorer_input_file(
        source_path=seed_path, dest_path=dest, seed_species=['[O]C=O'],
        method='MSC', bath_gas={'Ar': 1.0}, explore_tol=0.01, energy_tol=1.0e1, flux_tol=1.0e-6,
        maximum_radical_electrons=2)
    assert os.path.isfile(dest)
    assert summary is not None


@pytest.mark.skipif(not _HAS_RMG_ENV, reason="needs rmg_env for the RMG estimator + Arkane explorer")
@pytest.mark.slow
def test_generated_noncho2_seed_explores_a_real_network(tmp_path):
    """Generalization (verifier 2): a generated seed for a NON-CHO2 unimolecular species must drive
    a real Arkane explorer round to a genuine network. CHO2 passing proves nothing here -- every
    prior run used the same hand-built CHO2 seed. n-propyl radical is a classic unimolecular pdep
    system (isomerization + beta-scission), pure hydrocarbon, so it exercises C/H group coverage the
    CHO2 seed never touched, and its source is a single connected molecule (so it must NOT trip the
    disjoint-network refusal)."""
    from t3.pdep.pes_loop import PES_LOOP_DIAGRAM_ONLY, run_pes_loop
    from t3.schema import PESLoopConfig

    seed = str(tmp_path / 'network_nc3h7_source.py')
    generate_source_seed('[CH2]CC', seed, network_name='nC3H7', timeout=1800.0)

    cfg = PESLoopConfig(
        pes={'network': os.path.abspath(seed), 'source': ['[CH2]CC'], 'method': 'MSC',
             'bath_gas': {'Ar': 1.0}, 'explore_tol': 0.01, 'energy_tol': 1.0e1,
             'flux_tol': 1.0e-6, 'maximum_radical_electrons': 2, 'timeout': 1800.0},
        termination={'max_rounds': 1, 'stop_when_no_new_ts': True})

    result = run_pes_loop(cfg, str(tmp_path / 'run'), qm_runner=None)

    assert result.status == PES_LOOP_DIAGRAM_ONLY, (result.status, result.reason)
    assert result.final_network_path is not None and os.path.isfile(result.final_network_path)
    # The explorer grew a genuine multi-well/multi-channel network from the generated seed.
    with open(result.final_network_path) as f:
        explored = f.read()
    assert len(re.findall(r'^reaction\(', explored, re.M)) >= 1, "explorer discovered no reactions"
    assert len(re.findall(r'^transitionState\(', explored, re.M)) >= 1
    assert result.final_diagram_path is not None and os.path.isfile(result.final_diagram_path)
