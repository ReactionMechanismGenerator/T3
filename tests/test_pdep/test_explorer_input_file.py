#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_explorer_input_file module

Tests for ``t3.pdep.explorer.input_file.write_arkane_explorer_input_file``: transforming an RMG
P-dep network file (or a T3 hybrid network input, both plain Arkane DSL) into a valid Arkane
PES-explorer input file.

Most fixtures are hand-built strings so individual paths (kinetics emission, each refusal) can be
exercised in isolation. The real RMG-written fixture at ``tests/data/pdep_network/.../network4_2.py``
is exercised end-to-end as well (``test_real_rmg_network_file_writes_end_to_end``): its bath-gas
species carry no ``reactive`` keyword at all -- the exact shape ``arkane/pdep.py:654``
(``save_input_file``) emits for every species, and the shape issue #183 showed this writer refused
wholesale. All test outputs are written to pytest's ``tmp_path``, never into ``tests/data/``.
"""

import ast
import os
import pathlib

import pytest

from t3.pdep.explorer.input_file import (
    _ARKANE_SIDE_EFFECT_NAMES,
    _build_database_block,
    _validate_generated_statements,
    _ARKANE_VALUE_NAMES,
    _get_call_name,
    _render_literal,
    _TOP_LEVEL_CALL_NAMES,
    _validate_source_statements,
    ExplorerInputSummary,
    write_arkane_explorer_input_file,
)
from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.hashing import hash_file

# A minimal but structurally realistic RMG pdep network source: two species (one of them a bath
# gas with reactive=False), one transition state, one reaction with no explicit kinetics=, a
# network(...) block (to be dropped), and a pressureDependence(...) block (to be kept, with its
# method line rewritten).
SOURCE_NO_KINETICS = """species(
    label='methoxy',
    structure=SMILES('[CH2]O'),
    E0=(0, 'kJ/mol'),
)
species(
    label='CH2O',
    structure=SMILES('C=O'),
    E0=(-50, 'kJ/mol'),
)
species(
    label='H',
    structure=SMILES('[H]'),
    E0=(0, 'kJ/mol'),
)
species(
    label='He',
    structure=SMILES('[He]'),
    E0=(0, 'kJ/mol'),
    reactive=False,
)
transitionState(
    label='TS1',
    E0=(50, 'kJ/mol'),
)
reaction(
    label='CH2O+H=methoxy',
    reactants=['CH2O', 'H'],
    products=['methoxy'],
    transitionState='TS1',
)
network(
    label='PDepNetwork #1',
    isomers=['methoxy'],
    reactants=[('CH2O', 'H')],
    bathGas={'He': 1.0},
)
pressureDependence(
    label='PDepNetwork #1',
    Tmin=(300, 'K'),
    Tmax=(2000, 'K'),
    Tcount=8,
    Pmin=(0.01, 'bar'),
    Pmax=(100, 'bar'),
    Pcount=5,
    maximumGrainSize=(0.5, 'kcal/mol'),
    maximumGrainCount=250,
    method = 'modified strong collision',
    activeKRotor=True,
    activeJRotor=True,
    rmgmode=False,
)
"""

# Same shape, but the reaction() block carries an explicit kinetics= entry (as a hybrid network
# input would for an ILT-complement reaction), so no kinetics(...) job directive should be emitted.
SOURCE_WITH_KINETICS = SOURCE_NO_KINETICS.replace(
    "    transitionState='TS1',\n)",
    "    transitionState='TS1',\n"
    "    kinetics=Arrhenius(A=(1e13, 's^-1'), n=0, Ea=(50, 'kJ/mol')),\n)",
)

# A variant whose bath gas species carries no 'reactive' keyword at all -- the shape every
# RMG-written network file has (arkane/pdep.py:654 save_input_file never emits 'reactive' for any
# species; issue #183) -- to exercise the reactive=False injection path.
SOURCE_BATH_GAS_NO_REACTIVE_KEYWORD = SOURCE_NO_KINETICS.replace(
    "    reactive=False,\n)",
    ")",
)

# A variant whose bath gas species carries an EXPLICIT reactive=True: a self-contradictory
# declaration (named as bath gas while asserted reactive) that must be refused, never silently
# rewritten.
SOURCE_BATH_GAS_EXPLICITLY_REACTIVE = SOURCE_NO_KINETICS.replace(
    "    reactive=False,\n",
    "    reactive=True,\n",
)

# A variant whose bath gas species carries a NON-LITERAL reactive value, which cannot be verified
# either way and must be refused.
SOURCE_BATH_GAS_NON_LITERAL_REACTIVE = SOURCE_NO_KINETICS.replace(
    "    reactive=False,\n",
    "    reactive=SMILES('[He]'),\n",
)

# A variant where the reaction() block has no 'label' keyword at all, so a kinetics(...) directive
# cannot be targeted at it.
SOURCE_REACTION_NO_LABEL = SOURCE_NO_KINETICS.replace(
    "    label='CH2O+H=methoxy',\n", "",
)

# A variant carrying two pressureDependence(...) blocks (ambiguous method selection).
SOURCE_TWO_PDEP_BLOCKS = SOURCE_NO_KINETICS + SOURCE_NO_KINETICS.split("pressureDependence(")[1].join(
    ["pressureDependence(", ""]
)

# A variant with exactly one pressureDependence(...) block (so the pdep-node-count guard is
# satisfied), but a second, stray top-level 'method = ...' assignment outside that block, so the
# file has two lines matching METHOD_LINE_CANDIDATE_RE. This isolates the method-rewrite-count
# guard from the pdep-node-count guard (the inverse of
# test_refuses_second_pressure_dependence_block_lacking_its_own_method_line).
SOURCE_STRAY_METHOD_LINE = SOURCE_NO_KINETICS + "method = 'stray'\n"

DEFAULT_KWARGS = dict(
    seed_species=('methoxy',),
    bath_gas={'He': 1.0},
    method='MSC',
    explore_tol=0.01,
    energy_tol=80.0,
    flux_tol=1e-6,
    maximum_radical_electrons=1,
)


def _write_source(tmp_path, text, name='network.py'):
    source_path = str(tmp_path / name)
    with open(source_path, 'w') as f:
        f.write(text)
    return source_path


def _species_reactive_map(path: str) -> dict:
    """
    Parse a generated input file and map every species() label to its literal ``reactive`` keyword
    value, or ``None`` where the block carries no ``reactive`` keyword at all.

    Args:
        path (str): The generated input file path.

    Returns:
        dict: label -> literal reactive value (or None when absent).
    """
    with open(path) as f:
        tree = ast.parse(f.read())
    reactive_map = dict()
    for node in tree.body:
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call) \
                or not isinstance(node.value.func, ast.Name) or node.value.func.id != 'species':
            continue
        kwargs = {kw.arg: kw.value for kw in node.value.keywords if kw.arg is not None}
        label = ast.literal_eval(kwargs['label'])
        assert [kw.arg for kw in node.value.keywords].count('reactive') <= 1, \
            f'species {label!r} carries a duplicated reactive keyword'
        reactive_map[label] = ast.literal_eval(kwargs['reactive']) if 'reactive' in kwargs else None
    return reactive_map


def test_writes_valid_python_with_all_required_blocks(tmp_path):
    """Happy path: database, no network(...), pressureDependence kept+rewritten, kinetics(...)
    emitted for the kinetics-less reaction, explorer(...) last."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    result = write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)

    assert isinstance(result, ExplorerInputSummary)
    assert result.dest_path == dest_path
    assert os.path.isfile(dest_path)

    with open(dest_path) as f:
        text = f.read()

    # Self-parseable.
    tree = ast.parse(text)
    top_level_calls = [node.value.func.id for node in tree.body
                       if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                       and isinstance(node.value.func, ast.Name)]

    assert 'database' in top_level_calls
    assert 'network' not in top_level_calls
    assert 'pressureDependence' in top_level_calls
    assert 'explorer' in top_level_calls
    assert top_level_calls.index('explorer') == len(top_level_calls) - 1
    assert top_level_calls.index('explorer') > top_level_calls.index('pressureDependence')
    assert top_level_calls.count('kinetics') == 1
    assert result.kinetics_labels_emitted == ('CH2O+H=methoxy',)
    # Sanity: the method line was rewritten via rewrite_arkane_method_line, not merely copied
    # verbatim -- confirmed separately (with a distinguishable method) in
    # test_rewrites_method_line_for_msc.


def test_writes_when_expected_source_hash_matches(tmp_path):
    """
    A CORRECT ``expected_source_hash`` must not refuse the write.

    This is the over-refusal guard for the TOCTOU fix: without it, a bug that always treats the
    hash as mismatched (e.g. comparing against the wrong bytes, or a stale digest format) would
    make ``expected_source_hash`` unusable for every real caller while the suite stayed green,
    because every OTHER test in this module omits the parameter entirely.
    """
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')
    expected_hash = hash_file(source_path)

    result = write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                              expected_source_hash=expected_hash, **DEFAULT_KWARGS)

    assert isinstance(result, ExplorerInputSummary)
    assert os.path.isfile(dest_path)


def test_refuses_when_source_file_changed_after_hash_was_computed(tmp_path):
    """
    A MISMATCHED ``expected_source_hash`` must refuse the write -- this closes the TOCTOU gap.

    A caller (``t3.pdep.api.explore_pdep_network``) content-hashes the network file and then hands
    only the PATHNAME down to this writer, which reopens it. If the file is rewritten in between,
    the writer would otherwise silently explore bytes nobody approved. The expected hash is
    computed from the REAL, still-original file via ``t3.pdep.hashing.hash_file`` -- not a
    hardcoded digest -- so this test fails loudly if the hash format or algorithm ever drifts,
    rather than passing against a stale magic string.
    """
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')
    expected_hash = hash_file(source_path)

    # Rewrite the file after the hash was computed, simulating the TOCTOU window.
    with open(source_path, 'w') as f:
        f.write(SOURCE_WITH_KINETICS)

    with pytest.raises(ValueError, match='changed after it was validated'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                         expected_source_hash=expected_hash, **DEFAULT_KWARGS)

    assert not os.path.isfile(dest_path), \
        'A hash mismatch must be caught before anything is written to dest_path.'


def test_expected_source_hash_none_is_backwards_compatible(tmp_path):
    """
    Omitting ``expected_source_hash`` (the default, ``None``) must still write, unconditionally.

    Every caller that predates this parameter -- and any caller with no hash to bind to, e.g. one
    invoking this writer directly with a freshly-authored file -- must keep working exactly as
    before.
    """
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    result = write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                              expected_source_hash=None, **DEFAULT_KWARGS)

    assert isinstance(result, ExplorerInputSummary)
    assert os.path.isfile(dest_path)


def test_rewrites_method_line_for_msc(tmp_path):
    """The method line is rewritten via rewrite_arkane_method_line, not left as the source's."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                     **{**DEFAULT_KWARGS, 'method': 'CSE'})

    with open(dest_path) as f:
        text = f.read()
    assert 'chemically-significant eigenvalues' in text
    assert 'modified strong collision' not in text


# A source whose atomic-O species is written the way RMG's explorer writes triplet O(3P):
# spinMultiplicity=1 contradicting its own 'multiplicity 3' adjacency header. [OH] beside it is
# self-consistent (multiplicity 2 / spinMultiplicity 2) and must be left untouched. Outer ''' with
# inner """ so the adjacency-list triple-quoted strings nest cleanly.
SOURCE_TRIPLET_O_INCONSISTENT = '''species(
    label='[O]',
    structure=adjacencyList("""multiplicity 3
1 O u2 p2 c0
"""),
    E0=(243.034, 'kJ/mol'),
    spinMultiplicity=1,
)
species(
    label='[OH]',
    structure=adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0=(28.0, 'kJ/mol'),
    spinMultiplicity=2,
)
species(
    label='He',
    structure=SMILES('[He]'),
    E0=(0, 'kJ/mol'),
    reactive=False,
)
pressureDependence(
    label='PDepNetwork #1',
    Tmin=(300, 'K'),
    Tmax=(2000, 'K'),
    Tcount=8,
    Pmin=(0.01, 'bar'),
    Pmax=(100, 'bar'),
    Pcount=5,
    maximumGrainSize=(0.5, 'kcal/mol'),
    maximumGrainCount=250,
    method = 'modified strong collision',
    activeKRotor=True,
    activeJRotor=True,
    rmgmode=False,
)
'''


def _species_spin_multiplicity(text, label):
    """The spinMultiplicity literal written for ``label`` in a generated explorer input file."""
    tree = ast.parse(text)
    for node in tree.body:
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call) \
                or not isinstance(node.value.func, ast.Name) or node.value.func.id != 'species':
            continue
        kwargs = {kw.arg: kw.value for kw in node.value.keywords if kw.arg is not None}
        if 'label' in kwargs and ast.literal_eval(kwargs['label']) == label:
            return ast.literal_eval(kwargs['spinMultiplicity'])
    raise AssertionError(f'no species {label!r} with a spinMultiplicity in the generated file')


def test_corrects_triplet_o_spin_multiplicity_that_contradicts_its_adjacency_header(tmp_path):
    """Defect 3: RMG's explorer writes triplet atomic O with spinMultiplicity=1 while its adjacency
    header says 'multiplicity 3'. Splicing that verbatim makes the next round's Arkane crash, so the
    writer rewrites spinMultiplicity to match the (authoritative) adjacency list. The self-consistent
    [OH] beside it is left untouched -- this corrects a contradiction, it does not rewrite spins."""
    source_path = _write_source(tmp_path, SOURCE_TRIPLET_O_INCONSISTENT)
    dest_path = str(tmp_path / 'input.py')

    result = write_arkane_explorer_input_file(
        source_path=source_path, dest_path=dest_path,
        **{**DEFAULT_KWARGS, 'seed_species': ('[O]',), 'maximum_radical_electrons': 2})

    assert result.spin_multiplicity_corrected == (('[O]', 1, 3),)
    with open(dest_path) as f:
        text = f.read()
    assert _species_spin_multiplicity(text, '[O]') == 3      # corrected
    assert _species_spin_multiplicity(text, '[OH]') == 2     # consistent -> untouched
    # The correction is surfaced, never silent.
    assert any('[O]' in w and 'spinMultiplicity' in w for w in result.warnings)


def test_does_not_touch_a_consistent_spin_multiplicity(tmp_path):
    """A source whose spinMultiplicity already matches its adjacency header triggers no correction --
    the over-correction guard, so the fix cannot start rewriting good files."""
    consistent = SOURCE_TRIPLET_O_INCONSISTENT.replace('    spinMultiplicity=1,\n',
                                                       '    spinMultiplicity=3,\n')
    source_path = _write_source(tmp_path, consistent)
    dest_path = str(tmp_path / 'input.py')
    result = write_arkane_explorer_input_file(
        source_path=source_path, dest_path=dest_path,
        **{**DEFAULT_KWARGS, 'seed_species': ('[O]',), 'maximum_radical_electrons': 2})
    assert result.spin_multiplicity_corrected == ()


def test_does_not_emit_kinetics_directive_when_reaction_has_explicit_kinetics(tmp_path):
    """Behavior: a reaction() with an explicit kinetics= keyword gets no kinetics(...) directive."""
    source_path = _write_source(tmp_path, SOURCE_WITH_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    result = write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)

    with open(dest_path) as f:
        text = f.read()
    tree = ast.parse(text)
    kinetics_calls = [node for node in tree.body if isinstance(node, ast.Expr)
                     and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name)
                     and node.value.func.id == 'kinetics']
    assert kinetics_calls == []
    assert result.kinetics_labels_emitted == ()


@pytest.mark.parametrize('bad_label', [
    # A single quote plus a semicolon: naive interpolation ("kinetics('{label}')") yields
    # "kinetics('x'); import os; z=('a')" -- three semicolon-separated, individually valid
    # top-level statements, all on ONE line. This is the case the self-check's ast.parse alone
    # CANNOT catch: the injected text is still syntactically valid Python, so only explicit
    # validation of the label (not merely rendering it with !r) closes this hole.
    "x'); import os; z=('a",
    # A newline plus a quote: naive interpolation yields "kinetics('x')\nimport os\nkinetics('a')"
    # -- three separate top-level statements across three lines, again all valid syntax.
    "x')\nimport os\nkinetics('a",
    # A trailing backslash: naive interpolation yields "kinetics('x\\')" -- the backslash escapes
    # the closing quote, leaving the string literal unterminated. Unlike the two cases above this
    # one happens to also break the self-check's ast.parse (a SyntaxError), but validation must
    # refuse it up front regardless, exactly like the other injection characters.
    'x\\',
])
def test_refuses_reaction_label_with_injection_characters(tmp_path, bad_label):
    """Fail-closed (defect A): a reaction label read out of the (untrusted) source network file is
    interpolated into a kinetics('<label>') job directive in a NEW, executable Arkane input file.
    A label containing a quote, newline or backslash could inject or corrupt directives in that
    file; this must be refused rather than silently written. The first two parametrized cases are
    deliberately crafted so the injected text remains valid Python SYNTAX (still passes
    ast.parse), which is exactly why explicit validation -- not merely the self-check -- is
    required."""
    text = SOURCE_NO_KINETICS.replace("label='CH2O+H=methoxy',", f"label={bad_label!r},")
    source_path = _write_source(tmp_path, text)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError, match='newline, quote character or backslash'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
    assert not os.path.isfile(dest_path)


def test_benign_reaction_label_round_trips_with_correct_quoting(tmp_path):
    """A benign label (no injection characters) still produces exactly one kinetics(...) directive,
    rendered with !r so it is a well-formed Python string literal."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    result = write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)

    assert result.kinetics_labels_emitted == ('CH2O+H=methoxy',)
    with open(dest_path) as f:
        text = f.read()
    assert "kinetics('CH2O+H=methoxy')\n" in text


def test_refuses_reaction_with_no_kinetics_and_no_label(tmp_path):
    """Fail-closed: a reaction() with neither explicit kinetics= nor a usable label to target a
    kinetics(...) directive at must be refused, not silently written with null kinetics."""
    source_path = _write_source(tmp_path, SOURCE_REACTION_NO_LABEL)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
    assert not os.path.isfile(dest_path)


def test_refuses_more_than_one_pressure_dependence_block(tmp_path):
    """Fail-closed: Arkane binds to the FIRST pressureDependence(...) block via break, so a second
    one would be silently ignored; refuse rather than let the caller believe both apply."""
    source_path = _write_source(tmp_path, SOURCE_TWO_PDEP_BLOCKS)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
    assert not os.path.isfile(dest_path)


def test_refuses_second_pressure_dependence_block_lacking_its_own_method_line(tmp_path):
    """Isolates the pdep-node-count guard from the method-line-count guard: a second
    pressureDependence(...) block with no method line of its own still leaves exactly one
    'method = ...' line in the file, so only the dedicated pdep-count check (not the method-line
    count check) can catch this case."""
    stub_second_block = (
        "pressureDependence(\n"
        "    label='PDepNetwork #1',\n"
        "    Tmin=(300, 'K'),\n"
        "    Tmax=(2000, 'K'),\n"
        "    Tcount=8,\n"
        ")\n"
    )
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS + stub_second_block)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
    assert not os.path.isfile(dest_path)


def test_a_stray_method_line_outside_the_pressure_dependence_block_is_ignored(tmp_path):
    """
    DELIBERATE INVERSION: this used to be refused, and the refusal was over-refusal.

    The rewrite scanned every line in the file, so a stray top-level ``method = ...`` made the match
    count 2 and the "exactly one" guard fired. But the file is not actually ambiguous -- only the
    ``method`` inside the ``pressureDependence(...)`` block is the one Arkane reads, and a top-level
    ``method`` binding does nothing once the block has been called with its own keyword. Now that the
    rewrite is anchored to the block's AST node, the stray line is simply not a candidate, and the
    real method line is rewritten as intended.
    """
    source_path = _write_source(tmp_path, SOURCE_STRAY_METHOD_LINE)
    dest_path = str(tmp_path / 'input.py')

    write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
    generated = open(dest_path).read()
    assert "method = 'modified strong collision'" in generated, generated
    assert "method = 'stray'" in generated, 'The stray assignment is inert and is left alone.'


def test_refuses_a_method_line_that_is_only_outside_the_pressure_dependence_block(tmp_path):
    """
    The case the line-scan got silently WRONG, which the "exactly one match" guard could not catch.

    Here the ``pressureDependence(...)`` block carries NO method line while a docstring elsewhere
    contains a line that matches the pattern. The scan found exactly one candidate, so the count
    guard was satisfied -- and the rewrite landed inside the docstring, leaving Arkane to run the
    block with no resolved method at all. Anchoring to the node turns this into the same refusal as
    any other missing method line.
    """
    without_method = SOURCE_NO_KINETICS.replace("    method = 'modified strong collision',\n", '')
    assert "method = " not in without_method
    source = '"""Doc.\nmethod = \'CSE\'\n"""\n' + without_method
    source_path = _write_source(tmp_path, source)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError, match='found 0'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
    assert not os.path.isfile(dest_path)


@pytest.mark.parametrize('seed', [(), ('methoxy', 'CH2O', 'H'), ('nonexistent',)])
def test_refuses_invalid_seed_species_count_or_membership(tmp_path, seed):
    """Arkane's source must be exactly 1 or 2 species that exist as species() blocks (P2, P3)."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                         **{**DEFAULT_KWARGS, 'seed_species': seed})
    assert not os.path.isfile(dest_path)


def test_refuses_transition_state_label_as_seed(tmp_path):
    """P2: source names are resolved from species_dict only; a transitionState() label raises
    KeyError in Arkane, so it must be refused here first."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError, match='transitionState\\(\\) block, not a'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                         **{**DEFAULT_KWARGS, 'seed_species': ('TS1',)})
    assert not os.path.isfile(dest_path)


def test_refuses_label_defined_as_both_species_and_transition_state(tmp_path):
    """A label that is BOTH a species() and a transitionState() is ambiguous; refuse it."""
    text = SOURCE_NO_KINETICS.replace("label='TS1',", "label='methoxy',")
    source_path = _write_source(tmp_path, text)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                         **{**DEFAULT_KWARGS, 'seed_species': ('methoxy',)})
    assert not os.path.isfile(dest_path)


def test_accepts_two_species_seed(tmp_path):
    """A bimolecular (2-species) seed configuration is valid per P3."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    result = write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                              **{**DEFAULT_KWARGS, 'seed_species': ('CH2O', 'H')})
    assert result.seed_species == ('CH2O', 'H')
    with open(dest_path) as f:
        assert "source=['CH2O', 'H']" in f.read()


@pytest.mark.parametrize('bad_value', [float('inf'), float('-inf'), float('nan')])
def test_refuses_a_non_finite_tolerance(tmp_path, bad_value):
    """
    A non-finite tolerance is refused, and the message tells the caller to omit the argument.

    This is not a hypothetical input: ``float('inf')`` is Arkane's OWN default for ``energy_tol``
    (``arkane/input.py:496``) and the natural way to spell "no filter". Rendered with ``repr`` it
    becomes a bare ``inf`` in the generated input file -- valid Python SYNTAX, so this module's
    ast.parse self-check passes it, but an undefined NAME that raises when Arkane loads the file.
    Refusing here is what keeps "the self-check passed" from meaning less than it appears to.
    """
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError, match='finite'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                         **{**DEFAULT_KWARGS, 'energy_tol': bad_value})
    assert not os.path.isfile(dest_path)


def test_refuses_bath_gas_label_not_a_species_block(tmp_path):
    """P6: bath-gas species must exist as species() blocks, or Arkane KeyErrors mid-run."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError, match='not defined as a species\\(\\) block'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                         **{**DEFAULT_KWARGS, 'bath_gas': {'Ar': 1.0}})
    assert not os.path.isfile(dest_path)


def test_injects_reactive_false_for_bath_gas_species_without_reactive_keyword(tmp_path):
    """
    Issue #183, the core defect: RMG's ``save_input_file`` (arkane/pdep.py:654) never writes a
    ``reactive`` keyword for ANY species, so requiring a literal ``reactive=False`` in the SOURCE
    refused every real RMG-written network file. The fact "this species is the bath gas" is
    declared by the caller (``bath_gas``), and this writer GENERATES the destination file -- so it
    emits ``reactive = False`` on the declared bath-gas species itself: a faithful translation
    between two encodings of one fact, not a rewrite of user input.
    """
    source_path = _write_source(tmp_path, SOURCE_BATH_GAS_NO_REACTIVE_KEYWORD)
    dest_path = str(tmp_path / 'input.py')

    result = write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)

    assert result.reactive_false_injected == ('He',)
    assert _species_reactive_map(dest_path) == {'methoxy': None, 'CH2O': None, 'H': None, 'He': False}


def test_does_not_reinject_reactive_false_when_source_already_carries_it(tmp_path):
    """A source bath-gas species already marked reactive=False is left verbatim: no duplicate
    keyword (which would be a SyntaxError Arkane could never load), and no injection recorded."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    result = write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)

    assert result.reactive_false_injected == tuple()
    assert _species_reactive_map(dest_path)['He'] is False


def test_refuses_bath_gas_species_explicitly_marked_reactive(tmp_path):
    """A species named as bath gas while carrying an explicit reactive=True is a contradiction:
    honouring the name would rewrite the source's own claim, honouring the claim would run Arkane
    statmech on the collider. Refuse loudly rather than pick a side."""
    source_path = _write_source(tmp_path, SOURCE_BATH_GAS_EXPLICITLY_REACTIVE)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError, match='explicit.*reactive=True|reactive=True.*explicit'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
    assert not os.path.isfile(dest_path)


def test_refuses_bath_gas_species_with_non_literal_reactive_value(tmp_path):
    """A bath-gas species whose reactive value cannot be literally evaluated can be verified
    neither as a conflict nor as already-unreactive; injecting next to it could contradict what the
    expression evaluates to at Arkane load time."""
    source_path = _write_source(tmp_path, SOURCE_BATH_GAS_NON_LITERAL_REACTIVE)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError, match='literal'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
    assert not os.path.isfile(dest_path)


def test_real_rmg_network_file_writes_end_to_end(tmp_path):
    """
    THE issue #183 regression test: a real, RMG-written network file (species blocks carrying no
    ``reactive`` keyword anywhere, bath gas declared only in the ``network(...)`` block RMG writes)
    must be writable into a loadable explorer input. Before the injection fix, this writer refused
    every such file, i.e. the entire feature path was unusable on real RMG output.
    """
    source_path = os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'iteration_1', 'RMG', 'pdep',
                               'network4_2.py')
    with open(source_path) as f:
        assert 'reactive' not in f.read(), \
            'fixture precondition: a real RMG-written network file carries no reactive keyword'
    dest_path = str(tmp_path / 'input.py')

    result = write_arkane_explorer_input_file(
        source_path=source_path, dest_path=dest_path,
        seed_species=('C4rad(5)',), method='MSC',
        bath_gas={'He': 0.5, 'Ne': 0.5},
    )

    assert sorted(result.reactive_false_injected) == ['He', 'Ne']
    reactive_map = _species_reactive_map(dest_path)
    assert reactive_map['He'] is False
    assert reactive_map['Ne'] is False
    assert all(value is None for label, value in reactive_map.items() if label not in ('He', 'Ne'))
    # The generated text must still be a loadable input: valid Python with exactly one
    # pressureDependence block, no network block, and the explorer block appended last.
    with open(dest_path) as f:
        tree = ast.parse(f.read())
    top_level_calls = [node.value.func.id for node in tree.body
                       if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                       and isinstance(node.value.func, ast.Name)]
    assert 'network' not in top_level_calls
    assert top_level_calls.count('pressureDependence') == 1
    assert top_level_calls[-1] == 'explorer'


def test_multiple_bath_gases_warns_and_records_fractions_without_honoring_them(tmp_path):
    """P16: requested fractions do not survive an RMG network update (overwritten with equal
    weights); T3 must not promise fraction fidelity, only record it and warn."""
    text = SOURCE_NO_KINETICS.replace(
        "species(\n    label='He',\n    structure=SMILES('[He]'),\n    E0=(0, 'kJ/mol'),\n    reactive=False,\n)",
        "species(\n    label='He',\n    structure=SMILES('[He]'),\n    E0=(0, 'kJ/mol'),\n    reactive=False,\n)\n"
        "species(\n    label='Ne',\n    structure=SMILES('[Ne]'),\n    E0=(0, 'kJ/mol'),\n    reactive=False,\n)",
    )
    source_path = _write_source(tmp_path, text)
    dest_path = str(tmp_path / 'input.py')

    result = write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                              **{**DEFAULT_KWARGS, 'bath_gas': {'He': 0.5, 'Ne': 0.5}})
    assert any('equal weights' in warning or 'overwrit' in warning for warning in result.warnings)
    assert result.bath_gas == {'He': 0.5, 'Ne': 0.5}
    with open(dest_path) as f:
        text = f.read()
    assert "'He': 0.5" in text and "'Ne': 0.5" in text


def test_single_bath_gas_no_warning(tmp_path):
    """A single bath gas does not trigger the multi-bath-gas fraction warning."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    result = write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
    assert not any('equal weights' in warning or 'overwrit' in warning for warning in result.warnings)


def test_database_directive_uses_caller_supplied_kwargs(tmp_path):
    """database_kwargs overrides are rendered into the prepended database(...) directive."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                     database_kwargs={'thermoLibraries': ['BurkeH2O2']},
                                     **DEFAULT_KWARGS)
    with open(dest_path) as f:
        text = f.read()
    assert "['BurkeH2O2']" in text
    assert text.lstrip().startswith('database(')


def test_refuses_database_kwargs_value_with_unparseable_repr(tmp_path):
    """A database_kwargs value whose repr() is not valid Python syntax at all is refused by
    _render_literal (it is not one of the accepted literal types), before the generated text is
    ever built -- so this never reaches the self-check at all (see
    test_refuses_to_write_when_generated_block_value_fails_the_load_self_check for that path)."""
    class _UnparseableRepr:
        def __repr__(self):
            return '<broken>'

    # Asserted at BOTH layers, because they answer different questions and the field contract now
    # answers first: 'thermoLibraries' must be a list of str, which this object is not, so the
    # end-to-end refusal quotes the field's contract rather than the rendering one. _render_literal's
    # own rule is still the backstop for any field that has no contract of its own (and for the items
    # inside a list), so it is pinned directly rather than left to be inferred from the write.
    with pytest.raises(ValueError, match='cannot be rendered as a Python literal'):
        _render_literal('thermoLibraries', _UnparseableRepr())

    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError, match='thermoLibraries'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                         database_kwargs={'thermoLibraries': _UnparseableRepr()},
                                         **DEFAULT_KWARGS)
    assert not os.path.isfile(dest_path)
    assert list(tmp_path.glob('.explorer-input-*')) == []


def test_refuses_database_kwargs_value_that_is_a_path_object(tmp_path):
    """A pathlib.Path (or similar constructor-repr object) renders as valid Python SYNTAX
    (``PosixPath('/x')``) via plain repr(), so it would pass the old ast.parse-only self-check --
    but it is a call to a name this generated file never imports, so it cannot LOAD. _render_literal
    refuses it up front rather than relying on the self-check to catch it.

    Pinned at both layers for the reason given in
    ``test_refuses_database_kwargs_value_with_unparseable_repr``: a Path is also not a list of str, so
    'thermoLibraries' own contract refuses the write first. The Path is nested INSIDE the list in the
    end-to-end case, which is the shape where _render_literal is the check that actually fires."""
    with pytest.raises(ValueError, match='cannot be rendered as a Python literal'):
        _render_literal('thermoLibraries', pathlib.Path('/x'))

    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError, match='thermoLibraries'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                         database_kwargs={'thermoLibraries': [pathlib.Path('/x')]},
                                         **DEFAULT_KWARGS)
    assert not os.path.isfile(dest_path)


def test_database_kwargs_list_of_strings_still_renders_correctly(tmp_path):
    """A valid database_kwargs value (a list of strings) still round-trips through _render_literal
    and renders as the expected Python literal."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                     database_kwargs={'reactionLibraries': ['lib1', 'lib2']},
                                     **DEFAULT_KWARGS)
    with open(dest_path) as f:
        text = f.read()
    assert "reactionLibraries=['lib1', 'lib2']" in text


def test_refuses_maximum_radical_electrons_non_finite(tmp_path):
    """maximum_radical_electrons=float('inf') renders as the bare name 'inf' via plain repr() --
    valid syntax, undefined name at load time -- so it must be refused, matching the coverage
    _validate_tolerances already gives explore_tol/energy_tol/flux_tol."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError, match='finite'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                         **{**DEFAULT_KWARGS, 'maximum_radical_electrons': float('inf')})
    assert not os.path.isfile(dest_path)


@pytest.mark.parametrize('bad_fraction', [float('nan'), float('inf')])
def test_refuses_non_finite_bath_gas_fraction(tmp_path, bad_fraction):
    """A bath-gas mole fraction of nan/inf renders as an undefined bare name via plain repr();
    _render_literal must refuse it rather than let it reach Arkane as a broken bathGas dict."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError, match='finite'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                         **{**DEFAULT_KWARGS, 'bath_gas': {'He': bad_fraction}})
    assert not os.path.isfile(dest_path)


def test_refuses_to_write_when_generated_block_value_fails_the_load_self_check(tmp_path, monkeypatch):
    """Defense in depth: even if a bug elsewhere let a non-literal value slip past _render_literal,
    the strengthened self-check (which walks the GENERATED database(...)/explorer(...) blocks and
    asserts every keyword value is ast.literal_eval-able) must still catch it before anything is
    written to dest_path. Simulated here by monkeypatching _render_literal itself to reintroduce a
    plain, unguarded repr().

    Two upstream guards are disabled here, not one. The per-field contracts now refuse
    ``maximum_radical_electrons=float('inf')`` before any rendering happens, so monkeypatching
    ``_render_literal`` alone no longer reaches the self-check -- the value never gets that far. That
    is the correct behaviour and it is pinned elsewhere; what this test is for is the LAST line of
    defense, so it has to simulate a bug in both layers above it to get there. Stacking the two
    patches is what keeps this a test of the self-check rather than a second test of the field
    contracts."""
    import t3.pdep.explorer.input_file as input_file_module

    monkeypatch.setattr(input_file_module, '_render_literal', lambda field_name, value: repr(value))
    # The stub returns its keywords unchanged rather than None: the real function's contract is to
    # hand back the validated (possibly coerced) values for the caller to render, so a None-returning
    # stub would simulate a DIFFERENT bug -- a crash in the caller -- and never reach the self-check.
    monkeypatch.setattr(input_file_module, 'validate_explorer_field_values', lambda **kwargs: kwargs)

    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(RuntimeError, match='self-check'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path,
                                         **{**DEFAULT_KWARGS, 'maximum_radical_electrons': float('inf')})
    assert not os.path.isfile(dest_path)
    assert list(tmp_path.glob('.explorer-input-*')) == []


def test_a_refused_source_is_rejected_before_anything_is_written(tmp_path):
    """
    Validation runs before the write, so a refused source never reaches the filesystem.

    This test used to claim it covered the atomic write's temp-file cleanup. It did not, and could
    not: the source is refused during validation, long before a temp file is ever created, so it
    would have passed just as happily against a naive ``open(dest_path, 'w')``. What it actually
    proves -- that refusal is ordered before any write -- is worth keeping, under its real name.
    The cleanup path is covered by the test below.
    """
    source_path = _write_source(tmp_path, SOURCE_REACTION_NO_LABEL)
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)

    assert not os.path.isfile(dest_path)
    assert list(tmp_path.glob('.explorer-input-*')) == []


def test_atomic_write_cleans_up_its_temp_file_when_the_final_replace_fails(tmp_path, monkeypatch):
    """
    A failure at the atomic ``os.replace`` must leave neither a partial dest_path nor a temp file.

    This is the case that actually exercises the staging/cleanup path: the source is VALID, so the
    writer proceeds all the way to a fully staged temp file, and only the final replace fails. A
    naive ``open(dest_path, 'w')`` implementation fails this test too -- it never calls
    ``os.replace``, so the injected failure never fires -- which is what makes the test load-bearing
    rather than decorative.
    """
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')

    def _failing_replace(src, dst):
        raise OSError('simulated failure at the atomic replace')

    monkeypatch.setattr('t3.pdep.explorer.input_file.os.replace', _failing_replace)

    with pytest.raises(OSError, match='simulated failure'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)

    assert not os.path.isfile(dest_path)
    assert list(tmp_path.glob('.explorer-input-*')) == []


def test_refuses_source_that_is_not_valid_python(tmp_path):
    """The source file is parsed via ast.parse only (never exec/eval); a source that is not valid
    Python at all must be refused with a ValueError naming the parse failure, not raise a raw
    SyntaxError out of this module."""
    source_path = _write_source(tmp_path, "def broken(:\n    pass\n")
    dest_path = str(tmp_path / 'input.py')

    with pytest.raises(ValueError, match='Could not parse'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
    assert not os.path.isfile(dest_path)


def test_creates_destination_directory_if_missing(tmp_path):
    """Mirrors write_hybrid_network_input_file: the destination directory is created if absent."""
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'nested' / 'dir' / 'input.py')

    write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
    assert os.path.isfile(dest_path)


def test_refuses_an_empty_bath_gas(tmp_path):
    """
    An empty bath gas is refused, not passed through as ``bathGas={}``.

    Every bath-gas check is a loop over the dict, so an empty one satisfies all of them vacuously
    -- the guard reported "valid" precisely when there was nothing to validate. It is also
    unusable: this writer removes the source's ``network(...)`` block, which is Arkane's only other
    source of a bath gas, so ``bathGas={}`` fails deep inside the Arkane run itself with
    ``InputError('bathGas not specified in explorer block')``.
    """
    source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
    dest_path = str(tmp_path / 'input.py')
    kwargs = dict(DEFAULT_KWARGS)
    kwargs['bath_gas'] = {}

    with pytest.raises(ValueError, match='bath gas is required'):
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **kwargs)
    assert not os.path.isfile(dest_path)


class TestGetCallNameRefusesAttributeCallees:
    """
    ``_get_call_name`` used to return ``call.func.attr`` for an ``ast.Attribute`` callee, so
    ``foo.species(...)`` read back as ``'species'`` -- a spoof that would have let a non-DSL call
    impersonate a DSL one at both of the helper's call sites (the source collection loop and the
    self-check that re-walks the generated ``database(...)``/``explorer(...)`` blocks).

    These are DIRECT unit tests on the helper, deliberately not routed through
    ``write_arkane_explorer_input_file``. Mutation testing showed why: reverting the helper to trust
    ``.attr`` breaks NO end-to-end test, because ``_validate_source_statements`` now refuses an
    attribute call at the top level before the collection loop ever runs. The helper's contract is
    therefore only pinned here -- without these, the docstring's promise would be untested.
    """

    def test_a_plain_name_callee_returns_its_name(self):
        """The ordinary case: ``species(...)`` reads back as ``'species'``."""
        call = ast.parse('species(label="x")').body[0].value
        assert _get_call_name(call) == 'species'

    def test_an_attribute_callee_returns_none_rather_than_the_attribute_name(self):
        """``foo.species(...)`` must NOT be mistaken for a top-level ``species(...)`` DSL call."""
        call = ast.parse('foo.species(label="x")').body[0].value
        assert _get_call_name(call) is None

    def test_a_deeply_nested_attribute_callee_returns_none(self):
        """A longer dotted path must not surface its final attribute either."""
        call = ast.parse('a.b.c.species(label="x")').body[0].value
        assert _get_call_name(call) is None

    def test_a_non_name_non_attribute_callee_returns_none(self):
        """A computed callee (subscript, lambda result) has no plain name to report."""
        for source in ('handlers["species"](label="x")', '(lambda **kw: None)(label="x")'):
            call = ast.parse(source).body[0].value
            assert _get_call_name(call) is None, source


class TestSourceExpressionWhitelistRefusesBareGadgetNodes:
    """
    The recursive expression whitelist must refuse ``ast.Attribute``, ``ast.Subscript``, comparisons
    and comprehensions wherever they appear -- not merely when they happen to sit under a call whose
    callee is itself already refused.

    Mutation testing motivated this class. Allowing those four node kinds through the whitelist did
    NOT break the end-to-end escape-gadget test, because that gadget's outermost node is a ``Call``
    with an ``Attribute`` callee, which the separate call rule rejects anyway. The node kinds were
    consequently unpinned in every shape that does NOT end in such a call -- a bare attribute load,
    an indexed lookup, a comprehension -- each of which is a working first hop of a sandbox escape.
    """

    @staticmethod
    def _refuses(source: str) -> str:
        """Assert ``_validate_source_statements`` refuses ``source``, returning the message."""
        with pytest.raises(ValueError) as exc_info:
            _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')
        return str(exc_info.value)

    @pytest.mark.parametrize('value', [
        '().__class__',
        '().__class__.__base__.__subclasses__',
        'species.__globals__',
    ])
    def test_a_bare_attribute_load_as_an_assigned_value_is_refused(self, value):
        """``title = ().__class__`` is the first hop of the verified escape, so it must be refused."""
        assert 'source.py' in self._refuses(f'title = {value}\n')

    @pytest.mark.parametrize('value', [
        'data[0]',
        '{"a": 1}["a"]',
        '[1, 2, 3][0]',
    ])
    def test_a_subscript_as_an_assigned_value_is_refused(self, value):
        """Indexing is how the escape reaches ``__builtins__``; a literal container is not enough."""
        assert 'source.py' in self._refuses(f'title = {value}\n')

    @pytest.mark.parametrize('value', ['[c for c in (1, 2)]', '{c for c in (1, 2)}', '(c for c in (1, 2))'])
    def test_a_comprehension_as_an_assigned_value_is_refused(self, value):
        """Comprehensions are how the escape searches ``__subclasses__()`` for its target class."""
        assert 'source.py' in self._refuses(f'title = {value}\n')

    @pytest.mark.parametrize('value', ['1 == 1', '1 if True else 2', 'True and False', 'lambda: 0'])
    def test_other_non_literal_expression_kinds_are_refused(self, value):
        """Comparisons, conditionals, boolean ops and lambdas have no place in a literal assignment."""
        assert 'source.py' in self._refuses(f'title = {value}\n')

    @pytest.mark.parametrize('argument', [
        'structure=().__class__',
        'structure=data[0]',
        'structure=[c for c in (1, 2)]',
        'thermo=NASA(polynomials=[p.__class__ for p in (1, 2)])',
    ])
    def test_a_gadget_node_hidden_in_a_dsl_call_argument_is_refused(self, argument):
        """The whitelist recurses into arguments, including inside a legitimate nested constructor."""
        assert 'source.py' in self._refuses(f'species(label="x", {argument})\n')

    @pytest.mark.parametrize('source', [
        'title = "a name"',
        'useHinderedRotors = True',
        'modelChemistry = "CBS-QB3"',
        'description = """multi\nline"""',
        'species(label="x", structure=SMILES("CC"))',
        'species(label="x", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[1.0, -2.0])]))',
        'pressureDependence(Tmin=(300, "K"), grainSize=(0.5, "kcal/mol"))',
        'network(isomers=["a"], reactants=[("b", "c")], bathGas={"He": 1.0})',
    ])
    def test_the_shapes_real_inputs_actually_use_are_still_accepted(self, source):
        """
        Over-refusal guard: every shape here appears in a real RMG-Py or Arkane input file. The
        arithmetic shapes are in a separate test below, because measurement showed which arithmetic
        real files actually contain and it is narrower than one would guess.
        """
        _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')

    @pytest.mark.parametrize('source', [
        # Measured across all 100 example input files in RMG-Py: '-' (45 occurrences) and '*' (2) are
        # the ONLY arithmetic operators that appear, always over numeric operands.
        'species(label="x", energy=(204.06 - -213.357, "kJ/mol"))',
        'species(label="x", energy=(447.5 * 0.011962, "kJ/mol"))',
        # '/' and '+' are NOT in that corpus. They are permitted anyway, deliberately: over numeric
        # operands the worst either can do is produce a number or raise ZeroDivisionError while Arkane
        # loads the file, which fails closed. Unit arithmetic like this is plausible in a
        # hand-written network file, and refusing a safe shape on the strength of an examples
        # directory being exhaustive is a guess, not a guard. '**' is the one refused, because it
        # alone turns short source into an arbitrarily large value.
        'pressureDependence(maximumGrainSize=(1.0e5 / 2, "J/mol"))',
        'pressureDependence(Tmin=(300 + 15, "K"))',
    ])
    def test_numeric_arithmetic_is_accepted(self, source):
        """Arithmetic over numeric literals stays legal -- the DoS rule must not cost real inputs."""
        _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')


class TestModelChemistryCallFormIsAccepted:
    """
    A real ARC/hybrid network source expresses ``modelChemistry`` as a bare
    ``LevelOfTheory(...)``/``CompositeLevelOfTheory(...)`` call, which Arkane exec's at load time.
    That is not ``ast.literal_eval``-able, so the literal-assignment rule refused it. This branch is
    now bridged to ``t3.pdep.hybrid._validate_model_chemistry_expression`` -- the exact checker T3
    uses when it EMITS this directive -- so the two paths agree and no second allowlist is grown.
    The exception is keyed strictly to the ``modelChemistry`` target and to a genuine call node; any
    other target, and any non-call value, still falls through to the existing refusal.
    """

    @staticmethod
    def _refuses(source: str) -> str:
        with pytest.raises(ValueError) as exc_info:
            _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')
        return str(exc_info.value)

    @pytest.mark.parametrize('source', [
        "modelChemistry = LevelOfTheory(method='wb97xd', basis='def2tzvp')",
        "modelChemistry = LevelOfTheory(method='wb97xd', basis='def2tzvp', software='qchem')",
        "modelChemistry = CompositeLevelOfTheory("
        "freq=LevelOfTheory(method='wb97xd', basis='def2tzvp'), "
        "energy=LevelOfTheory(method='dlpno-ccsd(t)', basis='cc-pvtz'))",
    ])
    def test_the_call_form_for_model_chemistry_is_accepted(self, source):
        """The bare-call form real ARC/hybrid sources use must load, not be refused as non-literal."""
        _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')

    @pytest.mark.parametrize('source', [
        # A positional arg -- the checker refuses these; the refusal must surface as a source refusal.
        "modelChemistry = LevelOfTheory('wb97xd')",
        # An unknown keyword.
        "modelChemistry = LevelOfTheory(method='wb97xd', bogus='x')",
        # A CompositeLevelOfTheory missing a required field.
        "modelChemistry = CompositeLevelOfTheory(freq=LevelOfTheory(method='wb97xd'))",
        # A non-literal keyword value smuggled into the call.
        "modelChemistry = LevelOfTheory(method=().__class__)",
    ])
    def test_a_malformed_call_form_for_model_chemistry_is_still_refused(self, source):
        assert 'source.py' in self._refuses(source)

    @pytest.mark.parametrize('target', ['title', 'basis', 'level_of_theory'])
    def test_the_call_form_under_any_other_target_is_still_refused(self, target):
        """The exception is keyed to ``modelChemistry`` alone; the call form elsewhere is non-literal."""
        assert 'source.py' in self._refuses(f"{target} = LevelOfTheory(method='wb97xd')")

    @pytest.mark.parametrize('value', ['1 + 2', '().__class__', 'foo("x")'])
    def test_a_computed_non_call_model_chemistry_value_is_still_refused(self, value):
        """Guards the trap: the checker treats any non-call string as a plain label, so the branch
        must gate on a real call node -- a computed ``modelChemistry`` value stays refused."""
        assert 'source.py' in self._refuses(f'modelChemistry = {value}')


class TestSourceIsNarrowedToNetworkSourceSyntax:
    """
    Codex's round-29 P1 B and P1/P2 D: being a name Arkane defines is not a reason to splice a call
    to it. Every source here is syntactically valid Arkane input that the previous, whole-namespace
    whitelist accepted, and each one either changes the generated run or is a load-time bomb.
    """

    @pytest.mark.parametrize('source, expected', [
        # Codex's named attack: the generated database(...) block is PREPENDED, so a database(...)
        # spliced out of the source lands after it and wins.
        ('database(thermoLibraries=["attacker"], kineticsFamilies="none")\n', 'database'),
        # Each of these appends an extra Arkane job to the file the explorer job is generated into.
        ('explorer(source=["x"], explore_tol=0.1)\n', 'explorer'),
        ('kinetics("some_reaction")\n', 'kinetics'),
        ('statmech("x")\n', 'statmech'),
        ('thermo("x", "NASA")\n', 'thermo'),
        ('bac(level_of_theory="x")\n', 'bac'),
        ('ae(level_of_theory="x")\n', 'ae'),
    ])
    def test_a_job_or_state_directive_is_refused_as_a_top_level_statement(self, source, expected):
        """A directive Arkane defines but this module does not understand must not be carried along."""
        with pytest.raises(ValueError) as excinfo:
            _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')
        message = str(excinfo.value)
        assert expected in message, message
        # The refusal has to say WHY this specific name is refused, not merely that it is unlisted:
        # 'Arkane has never heard of this' and 'Arkane defines this and that is the problem' send a
        # reader to completely different places.
        assert 'that is exactly the problem' in message, message

    @pytest.mark.parametrize('source', [
        # A job directive reached from inside a value: the top-level statement is one this module DOES
        # understand, so only the recursive expression check can catch these.
        'species(label=kinetics("x"))\n',
        'species(label="x", thermo=thermo("x", "NASA"))\n',
        'reaction(label="r", kinetics=Arrhenius(A=(1.0, "s^-1"), Ea=(database(thermoLibraries=[]), "kJ/mol")))\n',
        # ... including the state-mutating '# Functions' names, which are legal only as whole
        # top-level statements.
        'species(label=species(label="nested"))\n',
        'network(isomers=[network(isomers=[])])\n',
    ])
    def test_a_job_or_state_directive_nested_inside_a_value_is_refused(self, source):
        """The five top-level directives are legal AS statements, never as parts of an expression."""
        with pytest.raises(ValueError, match='from inside an expression'):
            _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')

    @pytest.mark.parametrize('source, match', [
        # Arkane's directives take positional arguments too (its own examples call thermo('x', 'NASA')
        # that way), so a payload does not have to arrive as a keyword. Every other test in this class
        # passes its payload as a keyword, which left the positional path of the top-level call
        # unexercised -- a mutation that skipped validating args entirely passed the whole suite.
        ('species("x" * 1000000000)\n', 'non-numeric operand'),
        ('species(kinetics("x"))\n', 'from inside an expression'),
        ('network([0.0] * 100000000)\n', 'non-numeric operand'),
        ('reaction(range(1000000000))\n', 'not a recognized Arkane value constructor'),
        ('transitionState(().__class__)\n', 'Attribute'),
    ])
    def test_a_payload_in_a_positional_argument_is_refused(self, source, match):
        """A top-level directive's positional arguments are validated exactly like its keywords."""
        with pytest.raises(ValueError, match=match):
            _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')

    def test_range_is_refused(self):
        """
        ``range`` is in Arkane's namespace but is not an input-file value constructor, and it is the
        cheapest bomb in it: no real input file uses it, and ``range(10 ** 1000000)`` is evaluated
        when Arkane loads the generated file.
        """
        source = 'species(label="x", energies=range(1000000000))\n'
        with pytest.raises(ValueError, match='not a recognized Arkane value constructor'):
            _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')

    @pytest.mark.parametrize('source, match', [
        # A bare '**' over numeric operands: the ONLY one of these the operator rule itself catches.
        # The operand rule cannot -- the operands really are numbers.
        ('species(label="x", energy=(10 ** 9, "kJ/mol"))\n', 'Pow'),
        # The payloads Codex named. Note these are caught by the OPERAND rule, not the operator rule,
        # even though they contain '**': the outer operator is '*' and its left operand is a string,
        # so the refusal happens before the recursion ever reaches the '**'. Each rule is therefore
        # load-bearing on its own; neither subsumes the other.
        ('species(label="x" * 10 ** 9)\n', 'non-numeric operand'),
        ('species(label="x", thermo=NASA(polynomials=[0] * 10 ** 8))\n', 'non-numeric operand'),
        # ... and the same payloads written without '**' at all, which an operator-only rule would
        # miss completely.
        ('species(label="x" * 1000000000)\n', 'non-numeric operand'),
        ('species(label="x", thermo=NASA(polynomials=[0.0] * 100000000))\n', 'non-numeric operand'),
        ('species(label="x", comment=("a", "b") * 100000000)\n', 'non-numeric operand'),
    ])
    def test_a_literal_blowup_is_refused(self, source, match):
        """Repetition shares its syntax with multiplication; only the operands tell them apart."""
        with pytest.raises(ValueError, match=match):
            _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')

    def test_arithmetic_on_a_boolean_is_refused(self):
        """
        ``bool`` is an ``int`` subclass, so a type check written as ``isinstance(value, int)`` would
        silently accept ``True * 8``. Refusing it costs nothing (no input file means arithmetic on a
        flag) and keeps the numeric rule stateable in one sentence.
        """
        source = 'species(label="x", energy=(True * 8, "kJ/mol"))\n'
        with pytest.raises(ValueError, match='non-numeric operand'):
            _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')

    @pytest.mark.parametrize('source, match', [
        # The docstring/code mismatch Codex found: the comments always said 'literal assignment', but
        # the recursive expression whitelist happily allowed a call on the right-hand side.
        ('x = kinetics("some_reaction")\n', 'not a literal'),
        ('x = explorer(source=["y"])\n', 'not a literal'),
        ('x = SMILES("CC")\n', 'not a literal'),
        ('x = 10 ** 9\n', 'not a literal'),
    ])
    def test_a_computed_top_level_assignment_is_refused(self, source, match):
        """A top-level assignment must be a plain literal, as the surrounding comments always said."""
        with pytest.raises(ValueError, match=match):
            _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')

    @pytest.mark.parametrize('source', [
        "title = 'my network'\n",
        "modelChemistry = 'CBS-QB3'\n",
        'useHinderedRotors = True\n',
        'frequencies = [1.0, 2.0, 3.0]\n',
        'bathGas = {"He": 1.0}\n',
    ])
    def test_a_literal_top_level_assignment_is_still_accepted(self, source):
        """Over-refusal guard: the assignment forms real hand-written Arkane inputs use (84 of them
        across RMG-Py's examples) are all plain literals, and must keep working."""
        _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')

    @pytest.mark.parametrize('source', [
        # The verified sandbox escape, rewritten so its OUTERMOST node is a comprehension rather than
        # a call. The call-shaped spelling of the same gadget is refused by the call rule, which is
        # why this shape has to be pinned separately: a statement rule that permitted any bare
        # expression -- rather than only a bare CONSTANT -- would let this exact payload through
        # while every other test still passed.
        ("[c.__init__.__globals__['__builtins__']['__import__']('os').system('touch /tmp/PWNED')\n"
         " for c in ().__class__.__base__.__subclasses__()]\n"),
        # Other bare non-constant expression statements, none of which a network source ever contains.
        "().__class__.__base__.__subclasses__()[0]\n",
        "some_name\n",
        "().__class__\n",
        'f"{open}"\n',
    ])
    def test_a_bare_non_constant_expression_statement_is_refused(self, source):
        """
        Only a bare CONSTANT is inert. Any other bare expression can evaluate attribute chains,
        subscripts and comprehensions -- the whole vocabulary of the escape -- so 'it is just an
        expression statement, it cannot assign anything' is not a safety argument.
        """
        with pytest.raises(ValueError):
            _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')

    def test_a_module_docstring_is_accepted(self):
        """
        A bare literal statement is permitted: evaluating a constant and discarding it cannot do
        anything, so refusing a docstring would buy nothing and would refuse a plausible source.
        """
        source = '"""A hand-annotated network file."""\nspecies(label="x")\n'
        _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')


class TestGeneratedFileCannotBeInjectedThrough:
    """
    Codex's round-30 P0 and P1s. The source file is not the only untrusted text reaching the
    generated input: this module interpolates its OWN arguments into it too, and a keyword NAME is
    interpolated outside any quoting, where the careful literal-rendering of values does not reach.
    """

    def test_a_database_kwarg_key_cannot_inject_a_top_level_statement(self, tmp_path):
        """
        The verified P0. This exact key closed the generated ``database(...)``, inserted a top-level
        statement, and reopened a ``database(`` that the rendered ``=[...],\\n)`` tail completed -- so
        the result was valid Python and passed the generated-file self-check, which only inspected the
        keyword VALUES of ``database``/``explorer`` calls and so never looked at the injected
        statement. It executed when Arkane loaded the file.
        """
        evil_key = ("thermoLibraries=['x'])\n"
                    "__PWNED__ = ().__class__.__base__.__subclasses__()\n"
                    "database(zzz")
        with pytest.raises(ValueError, match='OUTSIDE any quoting'):
            _build_database_block({evil_key: ['primaryThermoLibrary']})

    @pytest.mark.parametrize('key', [
        "x)\nimport os\ndatabase(y",       # statement injection
        'has spaces',                      # not a keyword name at all
        'class',                           # a Python keyword: renders as a SyntaxError
        '__PWNED__',                       # a dunder
        'trailing\n',                      # a bare newline is enough to leave the keyword position
    ])
    def test_a_database_kwarg_key_must_be_a_plain_identifier(self, key):
        """Not a blocklist of dangerous characters: an identifier simply cannot leave its position."""
        with pytest.raises(ValueError):
            _build_database_block({key: ['x']})

    def test_the_real_database_kwargs_still_render(self):
        """
        Over-refusal guard: every default key, and a plausible override, must still work.

        The override used to be ``frequencyScaleFactor=0.99``, which was a fabricated fixture: that is
        not a keyword Arkane's ``database()`` accepts (``arkane/input.py:83-84``) -- it is a top-level
        Arkane setting elsewhere. So this "over-refusal guard" was pinning the writer's freedom to emit
        a keyword that would have made Arkane raise ``TypeError: database() got an unexpected keyword
        argument`` at load time. Replaced with a real one, measured off Arkane's signature.
        """
        block = _build_database_block({'thermoLibraries': ['primaryThermoLibrary'],
                                       'kineticsFamilies': 'default',
                                       'transportLibraries': ['primaryTransportLibrary']})
        assert "transportLibraries=['primaryTransportLibrary']" in block
        ast.parse(block)

    def test_the_generated_file_self_check_refuses_an_unexpected_top_level_statement(self, tmp_path):
        """
        The self-check's blind spot, closed generally rather than per-payload.

        The key guard above closes the one known route into the generated file. This asks the question
        the self-check should always have asked -- is every statement in the file one this module meant
        to emit? -- so a future splice or rendering defect that manufactures a statement is caught even
        though nobody predicted it. Simulated here by injecting into the finished text, because the
        argument-level guards now (correctly) prevent constructing such a file through the API.
        """
        source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
        dest_path = str(tmp_path / 'input.py')
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
        good = open(dest_path).read()
        for injected in ('__PWNED__ = ().__class__.__base__.__subclasses__()\n',
                         'import os\n',
                         'foo.bar()\n'):
            tree = ast.parse(good + injected)
            with pytest.raises(RuntimeError, match='failed its own self-check'):
                _validate_generated_statements(tree=tree, text=good + injected, dest_path=dest_path)

    def test_the_writer_actually_runs_the_generated_statement_self_check(self, tmp_path, monkeypatch):
        """
        Pins the WIRING, not just the function. Calling ``_validate_generated_statements`` directly
        (as the test above does) proves it works; it does not prove ``write_arkane_explorer_input_file``
        calls it -- deleting the call passed the entire suite.

        Since the argument-level guards now prevent building a poisoned file through the API, the
        defect the self-check exists to catch is simulated at its source: a block builder that emits
        an extra top-level statement, which is exactly what a future splice or rendering bug looks
        like. Nothing may be written to disk.
        """
        source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
        dest_path = str(tmp_path / 'input.py')
        monkeypatch.setattr(
            't3.pdep.explorer.input_file._build_database_block',
            lambda database_kwargs: "database(thermoLibraries=[])\n__PWNED__ = ().__class__\n\n")
        with pytest.raises(RuntimeError, match='failed its own self-check'):
            write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
        assert not os.path.isfile(dest_path), 'A file that failed its own self-check must never be written.'

    def test_the_self_check_accepts_the_file_this_module_actually_generates(self):
        """Over-refusal guard for the check above -- it must pass on real generated output."""
        # Covered end-to-end by every writer test (they all run the self-check), pinned here directly
        # so a future tightening that breaks the normal case fails with a pointed message.
        text = ('"""Doc."""\ntitle = \'x\'\ndatabase(thermoLibraries=[])\nspecies(label="a")\n'
                'transitionState(label="ts")\nreaction(label="r")\nkinetics(\'r\')\n'
                'network(isomers=["a"])\npressureDependence(method="MSC")\nexplorer(source=["a"])\n')
        _validate_generated_statements(tree=ast.parse(text), text=text, dest_path='/nonexistent/input.py')

    @pytest.mark.parametrize('target', ['explorer', 'database', 'species', 'kinetics', 'network',
                                        'Arrhenius', 'NASA', 'range'])
    def test_a_source_may_not_shadow_a_name_arkane_defines(self, target):
        """
        Codex's round-30 P1: a literal VALUE is inert, but the NAME it binds need not be.
        ``explorer = None`` is spliced above the generated ``explorer(...)`` directive and turns it
        into a TypeError at load time; ``species = None`` does the same to every species block after
        it. Allowing literal assignments (which real hand-written inputs need) is only safe if the
        target name is not one Arkane put in the namespace.
        """
        source = f'{target} = None\n'
        with pytest.raises(ValueError, match='a name Arkane itself defines'):
            _validate_source_statements(ast.parse(source), source, '/nonexistent/source.py')


class TestArkaneDSLNameTranscription:
    """
    The two name sets are hand transcriptions of the ``local_context`` dict Arkane injects when it
    ``exec``s an input file (RMG-Py ``arkane/input.py:607``). T3 must never import ``arkane``, so the
    sets cannot be compared against Arkane at run time -- these canaries pin the transcription.

    What they pin is a PARTITION, not a copy. An earlier version of this module allowed every name
    Arkane defines, on the reasoning that refusing a name Arkane does define would reject a
    legitimate input. That reasoning was wrong, and it is worth recording why, because it reads as
    conservative: Arkane's namespace is not one kind of thing. The names it groups under ``# Jobs``
    and ``# Functions`` exist FOR their side effects -- calling ``kinetics(...)`` appends a job,
    ``database(...)`` reconfigures the run -- so splicing one out of an untrusted source hands that
    source control over the explorer job this module generates. Being a real Arkane name is what
    makes those dangerous, not what makes them safe. Only the value constructors are inert.

    So the sets matter in both directions. A value constructor missing from ``_ARKANE_VALUE_NAMES``
    makes the writer refuse a legitimate network source that nests it (an over-refusal). A
    side-effect name wrongly listed there re-opens the hole. And a name in neither set is refused
    with the generic message, which is correct but less actionable.
    """

    def test_the_two_sets_partition_arkanes_callable_namespace(self):
        """
        Together the sets must be exactly Arkane's ``local_context`` callables, and must not overlap.

        Overlap would be a live contradiction rather than a style problem: ``_ARKANE_SIDE_EFFECT_NAMES``
        is tested first in ``_validate_source_expression``, so a name in both sets would be refused
        inside expressions while the reader of ``_ARKANE_VALUE_NAMES`` would reasonably expect it to
        pass.
        """
        # Transcribed from RMG-Py arkane/input.py:607, using Arkane's own group comments there.
        # Excluded: '__builtins__' (set to None precisely so the exec has no builtins), the legacy
        # 'True'/'False' string keys (aliases of the bare constants, not callable), and 'range',
        # which is a builtin Arkane re-exposes rather than an input-file value constructor -- and
        # which is how 'range(10 ** 1000000)' would otherwise pass as a "recognized DSL call".
        value_constructors = {
            'SMILES', 'adjacencyList', 'InChI', 'LevelOfTheory', 'CompositeLevelOfTheory',
            'TransportData', 'SingleExponentialDown', 'Arrhenius', 'RateUncertainty',
            'IdealGasTranslation', 'LinearRotor', 'NonlinearRotor', 'KRotor', 'SphericalTopRotor',
            'HarmonicOscillator', 'HinderedRotor', 'FreeRotor',
            'ThermoData', 'Wilhoit', 'NASA', 'NASAPolynomial',
        }
        side_effecting = {
            # Arkane's own '# Jobs' group: each call appends a job to the run.
            'kinetics', 'statmech', 'thermo', 'pressureDependence', 'explorer', 'bac', 'ae',
            # Arkane's own '# Functions' group: each call mutates the loader's species/reaction/
            # database state.
            'species', 'transitionState', 'reaction', 'network', 'database',
        }
        assert set(_ARKANE_VALUE_NAMES) == value_constructors
        assert set(_ARKANE_SIDE_EFFECT_NAMES) == side_effecting
        assert not (set(_ARKANE_VALUE_NAMES) & set(_ARKANE_SIDE_EFFECT_NAMES))
        assert 'range' not in set(_ARKANE_VALUE_NAMES) | set(_ARKANE_SIDE_EFFECT_NAMES)

    def test_the_names_a_real_rmg_written_network_uses_are_allowed_where_it_uses_them(self):
        """
        The five directives RMG writes as top-level statements must be permitted THERE, and the
        constructors it nests must be permitted where they are nested -- the distinction the two
        sets exist to draw.
        """
        # Verified against 33 real RMG-written examples/rmg/polystyrene/pdep/network*.py files: every
        # top-level statement is an expression-statement calling one of these five, and every nested
        # call is one of the constructors below (measured; zero job directives, zero 'range').
        for name in ('species', 'transitionState', 'reaction', 'network', 'pressureDependence'):
            assert name in _TOP_LEVEL_CALL_NAMES
        for name in ('Arrhenius', 'HarmonicOscillator', 'HinderedRotor', 'NASA', 'NASAPolynomial',
                     'SingleExponentialDown', 'TransportData', 'SMILES', 'IdealGasTranslation',
                     'LinearRotor', 'NonlinearRotor'):
            assert name in _ARKANE_VALUE_NAMES

    def test_a_job_directive_is_never_a_permitted_top_level_statement(self):
        """
        The narrowing itself: no name whose purpose is a side effect may be spliced as a top-level
        statement. ``_TOP_LEVEL_CALL_NAMES`` does contain the five ``# Functions`` names -- a network
        source IS a sequence of those -- but must contain none of the ``# Jobs`` names.
        """
        for name in ('kinetics', 'statmech', 'thermo', 'explorer', 'bac', 'ae', 'database'):
            assert name not in _TOP_LEVEL_CALL_NAMES

    def test_the_escape_gadgets_building_blocks_are_not_allowed(self):
        """No dunder or builtin used by the verified sandbox escape may be an allowed callee."""
        for name in ('__import__', 'eval', 'exec', 'open', 'getattr', 'compile', 'input', '__build_class__'):
            assert name not in _ARKANE_VALUE_NAMES
            assert name not in _TOP_LEVEL_CALL_NAMES


class TestSourceStatementValidationClosesTheExecHole:
    """
    Regression tests for the remote-code-execution hole this module used to have:
    ``write_arkane_explorer_input_file`` only ``ast.parse``s the source (never ``exec``s it), but
    then spliced the source's text VERBATIM into a generated ``input.py`` that Arkane later
    ``exec``s. The collection loop only recognized ``Expr(Call(Name))`` nodes whose callee was one
    of ``_TOP_LEVEL_CALL_NAMES``; every other top-level statement -- an ``import``, a ``def``, an
    assignment computed from an arbitrary call -- was silently skipped by that loop and therefore
    survived into the written file untouched.

    These tests exercise the fix, ``_validate_source_statements``/``_validate_source_expression``,
    which now runs immediately after ``ast.parse`` and before anything else (including any
    filesystem write).
    """

    def test_refuses_top_level_import_statement(self, tmp_path):
        source_path = _write_source(tmp_path, "import os\n")
        dest_path = str(tmp_path / 'input.py')

        with pytest.raises(ValueError, match='top-level Import statement'):
            write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
        assert not os.path.isfile(dest_path)

    def test_refuses_top_level_from_import_statement(self, tmp_path):
        source_path = _write_source(tmp_path, "from os import system\n")
        dest_path = str(tmp_path / 'input.py')

        with pytest.raises(ValueError, match='top-level ImportFrom statement'):
            write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
        assert not os.path.isfile(dest_path)

    def test_refuses_the_real_builtins_sandbox_escape_gadget(self, tmp_path):
        """Arkane's ``global_context = {'__builtins__': None}`` (RMG-Py arkane/input.py:607) is not
        a sandbox: a pure attribute/subscript/comprehension chain escapes it and reaches
        '__import__'. This is that exact gadget, as a top-level assignment."""
        gadget = (
            "pwn = [c for c in ().__class__.__base__.__subclasses__() if c.__name__ == "
            "'catch_warnings'][0].__init__.__globals__['__builtins__']['__import__']"
            "('os').system('touch /tmp/PWNED')\n"
        )
        source_path = _write_source(tmp_path, gadget)
        dest_path = str(tmp_path / 'input.py')

        # Matched against this module's own refusal wording rather than a bare ValueError, so this
        # fails for the RIGHT reason pre-fix: the old code did not validate a top-level assignment's
        # value at all, so it would instead (wrongly) fail later, deep in unrelated collection logic
        # (a missing pressureDependence(...) block), not here.
        with pytest.raises(ValueError, match='spliced verbatim into a NEW file Arkane will exec'):
            write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
        assert not os.path.isfile(dest_path)

    def test_refuses_attribute_call_and_does_not_mistake_it_for_a_dsl_call(self, tmp_path):
        """foo.species(...) must not be treated as a recognized 'species' DSL call -- this is the
        defect in the OLD _get_call_name, which returned call.func.attr for an ast.Attribute."""
        source_path = _write_source(tmp_path, "foo.species(label='x')\n")
        dest_path = str(tmp_path / 'input.py')

        with pytest.raises(ValueError, match='non-plain-name callee') as excinfo:
            write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
        assert "recognized 'species'" not in str(excinfo.value)
        assert not os.path.isfile(dest_path)

    @pytest.mark.parametrize('bad_statement,expected_node_type', [
        ("def helper():\n    pass\n", 'FunctionDef'),
        ("class Helper:\n    pass\n", 'ClassDef'),
        ("if True:\n    species(label='x')\n", 'If'),
        ("for i in range(1):\n    pass\n", 'For'),
    ])
    def test_refuses_def_class_if_for_top_level_statements(self, tmp_path, bad_statement, expected_node_type):
        source_path = _write_source(tmp_path, bad_statement)
        dest_path = str(tmp_path / 'input.py')

        # Matched against the specific offending node type, so this fails for the RIGHT reason
        # pre-fix: the old collection loop simply 'continue'd past any node that was not
        # Expr(Call(...)), so these would otherwise have raised the unrelated 'pressureDependence'
        # or 'seed_species' ValueErrors instead, not one naming this statement.
        with pytest.raises(ValueError, match=f'top-level {expected_node_type} statement'):
            write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
        assert not os.path.isfile(dest_path)

    def test_accepts_name_equals_literal_assignments_and_they_survive_verbatim(self, tmp_path):
        """Over-refusal regression guard: RMG/Arkane inputs legitimately carry top-level
        'Name = <literal>' assignments (title, description, modelChemistry, useHinderedRotors,
        useBondCorrections -- see the hand-written RMG-Py/examples/arkane fixtures); these must be
        ACCEPTED and must survive into the generated file untouched."""
        preamble = (
            "title = 'A test title'\n"
            "description = 'A test description'\n"
            "modelChemistry = 'CBS-QB3'\n"
            "useHinderedRotors = True\n"
            "useBondCorrections = False\n"
        )
        source_path = _write_source(tmp_path, preamble + SOURCE_NO_KINETICS)
        dest_path = str(tmp_path / 'input.py')

        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)

        assert os.path.isfile(dest_path)
        with open(dest_path) as f:
            text = f.read()
        assert "title = 'A test title'" in text
        assert "description = 'A test description'" in text
        assert "modelChemistry = 'CBS-QB3'" in text
        assert "useHinderedRotors = True" in text
        assert "useBondCorrections = False" in text

    def test_accepts_nested_legitimate_constructor_calls(self, tmp_path):
        """A species() block nesting SMILES(...)/NASA(...)/NASAPolynomial(...) -- all real Arkane
        DSL functions -- must be accepted, not merely a bare top-level DSL call."""
        extra_species = (
            "species(\n"
            "    label='extra',\n"
            "    structure=SMILES('CC'),\n"
            "    thermo=NASA(polynomials=[NASAPolynomial(coeffs=[1.0, -2.0])]),\n"
            ")\n"
        )
        source_path = _write_source(tmp_path, SOURCE_NO_KINETICS + extra_species)
        dest_path = str(tmp_path / 'input.py')

        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)

        assert os.path.isfile(dest_path)
        with open(dest_path) as f:
            text = f.read()
        assert 'NASAPolynomial' in text

    @pytest.mark.parametrize('bad_call', [
        "species(label=eval('1'))\n",
        "species(label=__import__('os'))\n",
    ])
    def test_refuses_nested_non_dsl_calls(self, tmp_path, bad_call):
        source_path = _write_source(tmp_path, bad_call)
        dest_path = str(tmp_path / 'input.py')

        with pytest.raises(ValueError):
            write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
        assert not os.path.isfile(dest_path)

    def test_refuses_kwargs_unpacking_in_a_dsl_call(self, tmp_path):
        source_path = _write_source(tmp_path, "species(**{'label': 'x'})\n")
        dest_path = str(tmp_path / 'input.py')

        with pytest.raises(ValueError, match='keyword-argument unpacking'):
            write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)
        assert not os.path.isfile(dest_path)

    def test_accepts_negative_number_and_arithmetic_inside_a_dsl_call(self, tmp_path):
        """Regression guard for the UnaryOp/BinOp allowances: a negative number (already present in
        SOURCE_NO_KINETICS's CH2O species E0) and an arithmetic expression inside a DSL call
        argument must both be accepted, not refused as 'not a literal'."""
        source = SOURCE_NO_KINETICS.replace(
            "    maximumGrainSize=(0.5, 'kcal/mol'),\n",
            "    maximumGrainSize=(0.5 + 0.1, 'kcal/mol'),\n",
        )
        source_path = _write_source(tmp_path, source)
        dest_path = str(tmp_path / 'input.py')

        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)

        assert os.path.isfile(dest_path)

    def test_refusal_happens_before_anything_is_written_to_disk(self, tmp_path):
        source_path = _write_source(tmp_path, "import os\n")
        dest_path = str(tmp_path / 'input.py')

        with pytest.raises(ValueError):
            write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **DEFAULT_KWARGS)

        assert not os.path.isfile(dest_path)
        assert list(tmp_path.glob('.explorer-input-*')) == []


class TestFieldValuesAreValidatedPerFieldNotJustPerLiteralType:
    """
    Every rendered value is checked against ITS OWN field's contract, not just "is a literal".

    ``_render_literal`` answers "can this be written as a Python literal that will LOAD?" -- and it
    answers it correctly. The bug was treating that as the whole question. It is not the caller's
    question, which is "is this a plausible value for THIS field?", and the two have different
    answers per field: ``'2'`` and ``True`` are both perfectly renderable literals, and both are
    nonsense as ``maximumRadicalElectrons``. The consequence is not a crash in T3 -- it is Arkane
    accepting the file, running, and either dying deep inside the job or (worse) silently computing
    the wrong thing:

    * ``maximumRadicalElectrons=True`` renders as ``True`` and silently means 1
      (``rmgpy/constraints.py:144-147`` compares with ``>``, and ``True == 1``).
    * ``maximumRadicalElectrons='2'`` renders as ``'2'`` and reaches that same ``>`` as a str,
      raising a TypeError after the exploration has already been set up.
    * ``explore_tol=True`` passes ``math.isfinite`` (``math.isfinite(True)`` is True) and silently
      means an explore tolerance of 1.0 (``arkane/explorer.py:235``).
    * ``explore_tol=0`` or a negative value makes ``leak > explore_tol * kchar`` true for every
      isomer with any leak at all -- an unbounded exploration, not a filter.
    * ``bathGas={'He': 0.5}`` sums to 0.5, and the fractions are used UNNORMALIZED as a linear
      weight for sigma/mw and as an EXPONENT for epsilon (``rmgpy/pdep/configuration.pyx:190-197``),
      so the collision parameters come out physically wrong with nothing raising anywhere.

    Values are measured against real Arkane usage, not invented: see ``arkane/input.py:496`` for the
    defaults and ``examples/arkane/explorer/`` for the observed ranges.
    """

    @staticmethod
    def _write(tmp_path, **overrides):
        """Attempt a write with DEFAULT_KWARGS plus overrides, returning the generated text."""
        source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
        dest_path = str(tmp_path / 'input.py')
        kwargs = dict(DEFAULT_KWARGS)
        kwargs.update(overrides)
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **kwargs)
        with open(dest_path, 'r') as f:
            return f.read()

    @pytest.mark.parametrize('field,value,reason', [
        ('explore_tol', True, 'a bool passes math.isfinite and silently means a tolerance of 1.0'),
        ('explore_tol', False, 'False silently means a tolerance of 0.0'),
        ('explore_tol', '0.01', 'a str reaches Arkane and dies on explore_tol * kchar'),
        ('explore_tol', -0.01, 'a negative threshold is exceeded by a leak of zero: not a tolerance'),
        ('energy_tol', True, 'a bool silently means an energy tolerance of 1.0 RT'),
        ('energy_tol', '80', 'a str dies inside get_energy_filtered_reactions on tol * R * T'),
        ('energy_tol', -80.0, 'a negative dE removes channels BELOW the source: the opposite filter'),
        ('flux_tol', True, 'a bool silently means a flux tolerance of 1.0'),
        ('flux_tol', '1e-6', 'a str dies on the flux comparison'),
        ('flux_tol', -1e-6, 'a negative flux tolerance cannot filter anything'),
        ('maximum_radical_electrons', True, 'renders as True and silently means 1'),
        ('maximum_radical_electrons', '2', "reaches rmgpy/constraints.py's '>' as a str"),
        ('maximum_radical_electrons', 2.5, 'a fractional count of electrons is not a count'),
        ('maximum_radical_electrons', -1, 'a negative maximum forbids every species'),
    ])
    def test_a_renderable_literal_that_is_wrong_for_its_field_is_refused(self, tmp_path, field, value, reason):
        """Each of these is a valid Python literal that renders fine and is wrong for its field."""
        with pytest.raises(ValueError, match=field):
            self._write(tmp_path, **{field: value})

    @pytest.mark.parametrize('field,value', [
        ('explore_tol', 0.01),          # examples/arkane/explorer/methoxy/input.py:209
        ('explore_tol', 1.0e-1),        # examples/arkane/explorer/methyl+formaldehyde/input.py:48
        ('energy_tol', 8e1),            # examples/arkane/explorer/methoxy/input.py:210
        ('energy_tol', 45),             # an int is a legitimate spelling of a dimensionless multiple of RT
        ('flux_tol', 1e-15),            # arkane/data/methoxy_explore.py:243
        ('flux_tol', 0.0),              # Arkane's OWN default (arkane/input.py:496): "no flux filter"
        ('maximum_radical_electrons', 1),   # examples/arkane/explorer/methyl+formaldehyde/input.py:50
        ('maximum_radical_electrons', 2),   # examples/arkane/explorer/methoxy/input.py:213
        ('maximum_radical_electrons', 0),   # "no radicals at all" is a coherent request
    ])
    def test_every_value_real_arkane_inputs_actually_use_is_still_accepted(self, tmp_path, field, value):
        """
        The contracts must not reject the corpus they were derived from.

        ``flux_tol=0.0`` is the trap in this list: it is Arkane's own default and means "apply no
        flux filter", so a blanket "tolerances must be positive" rule -- which is right for the
        other two -- would refuse a legitimate, documented value. ``energy_tol=45`` is the second:
        the tolerances are dimensionless multiples of RT, not quantities, so an int spelling is
        exactly as valid as a float and rejecting non-floats would be arbitrary.
        """
        text = self._write(tmp_path, **{field: value})
        assert 'explorer(' in text

    @pytest.mark.parametrize('field', ['energy_tol', 'explore_tol'])
    @pytest.mark.parametrize('value', [0, 0.0])
    def test_a_degenerate_but_well_defined_zero_tolerance_is_accepted(self, tmp_path, field, value):
        """
        ``energy_tol=0`` is an exact cutoff, not a broken one, and must not be refused.

        The mechanics, verified against RMG-Py rather than assumed: Arkane gates the filter on
        ``self.energy_tol != np.inf`` (arkane/explorer.py:303,309), so zero is passed straight through;
        ``get_energy_filtered_reactions`` forms ``dE = tol * R * T`` and collects every reaction with
        ``E0 - E0source > dE`` (rmgpy/rmg/pdep.py:331,353); and that collected set is what
        ``remove_reactions`` DELETES (arkane/explorer.py:336). So at zero the filter removes every
        channel strictly above the source -- destructive, but exactly what was asked for.

        ``explore_tol=0`` is the same case at the other end: Arkane explores whenever
        ``get_leak_coefficient(...) > explore_tol * kchar`` (arkane/explorer.py:235), so zero means
        "explore on any leak at all" -- an effectively unbounded exploration, but a defined one.

        The reason refusing either was wrong is not that zero is harmless. It is that the bound bought
        nothing: ``energy_tol=1e-300`` is just as destructive, ``explore_tol=1e-300`` is just as
        unbounded, and both were always accepted -- so the rule refused one particular spelling of a
        value Arkane itself admits while providing no protection against the failure mode at all.
        Bounding the runaway cost of a near-zero ``explore_tol`` is a budget question for the caller,
        not something a lower bound on this field can express.

        A negative value is still refused for both, and that is not the same judgement call. A
        negative ``dE`` removes channels BELOW the source energy, inverting the documented meaning of
        the parameter rather than tightening it; a negative ``explore_tol`` makes the threshold
        negative, so a leak coefficient of exactly zero "exceeds" it and a network with no leak at all
        is explored -- neither is a looser tolerance, both are a different comparison.
        """
        text = self._write(tmp_path, **{field: value})
        assert f'{field}=' in text

    @pytest.mark.parametrize('value', [2.0, 0.0, 1.0e1])
    def test_an_integral_float_count_is_accepted_and_rendered_as_an_int(self, tmp_path, value):
        """
        ``maximum_radical_electrons`` is integer-only, but a float that IS a whole number is a
        legitimate spelling of one and must not be refused.

        This is not a hypothetical spelling: the value reaches this writer from T3's YAML config, and
        YAML parses ``2.0`` as a float, so refusing it would refuse a config file that asks for
        exactly what an accepted ``2`` asks for. The float is coerced rather than merely tolerated,
        because ``maximumRadicalElectrons=2.0`` in the generated file is a count written as a
        non-count -- it works, but it is not what a human writes and it is not what Arkane's own
        fixtures contain.
        """
        text = self._write(tmp_path, maximum_radical_electrons=value)
        assert f'maximumRadicalElectrons={int(value)}' in text
        assert f'maximumRadicalElectrons={value!r}' not in text

    @pytest.mark.parametrize('value', [2.5, -0.5, True, False, '2', float('inf'), float('nan')])
    def test_a_non_integral_count_is_still_refused_after_integral_floats_are_admitted(self, tmp_path, value):
        """
        Admitting an integral float must not widen the field to floats in general.

        ``True`` and ``False`` are the ones that would slip through a naive ``float(value).is_integer()``
        coercion: a bool is an int subclass, ``float(True).is_integer()`` is True, and ``int(True)`` is
        1 -- so the coercion would render a silent ``maximumRadicalElectrons=1`` from a caller who
        passed a flag. ``-0.5`` is integral in neither direction and also out of range, so it must be
        refused whichever check reaches it first.
        """
        with pytest.raises(ValueError, match='maximum_radical_electrons'):
            self._write(tmp_path, maximum_radical_electrons=value)

    @pytest.mark.parametrize('bath_gas,should_pass', [
        ({'He': 1.0}, True),            # examples/arkane/explorer/methoxy/input.py:212
        ({'He': 1}, True),              # the int spelling, as arkane/data/methoxy_explore.py:216 uses
        ({'He': 0.5, 'Ar': 0.5}, True),
        ({'He': 0.5}, False),           # sums to 0.5: silently wrong collision parameters
        ({'He': 1.5}, False),           # refused by the sum rule, which is the only rule that can see it
        ({'He': 0.5, 'Ar': 0.6}, False),
        ({'He': True}, False),          # renders as True, silently means 1.0
        ({'He': '1.0'}, False),         # a str fraction dies inside the collision-frequency maths
        ({'He': 0.0}, False),           # sum 0: caught by the sum rule
        ({'He': -1.0}, False),          # sum -1: caught by the sum rule
        # These two are the cases that isolate the PER-FRACTION rule: each composition sums to
        # exactly 1.0, so the sum rule is satisfied and only the per-species bound can refuse them.
        # Without them the per-fraction lower bound is unpinned -- neutralising it to -1e99 leaves the
        # suite green, because every other bad-fraction case here happens to break the sum as well.
        # They are not hypothetical: a negative fraction reaches ``epsilon ** frac``
        # (rmgpy/pdep/configuration.pyx:193) as a negative exponent, and a zero fraction makes that
        # term exactly 1.0 -- both silently wrong, neither raising.
        ({'He': -0.5, 'Ar': 1.5}, False),
        ({'He': 0.0, 'Ar': 1.0}, False),
    ])
    def test_bath_gas_fractions_are_validated_as_mole_fractions(self, tmp_path, bath_gas, should_pass):
        """
        The fractions must be mole fractions, and they must sum to 1.

        This is deliberately STRICTER than Arkane, which never checks either (``arkane/input.py:508``
        is a pure passthrough). That is the point: because nothing downstream normalizes them
        (``rmgpy/pdep/configuration.pyx:190-197`` uses them as raw weights and as an exponent), a set
        of fractions that does not sum to 1 does not fail -- it produces a physically wrong sigma,
        epsilon and molecular weight for the bath gas and reports success. A validator that only
        refuses what Arkane would refuse cannot catch that class of error at all.
        """
        # Every species named here must exist in the source and be marked reactive=False (P16), or
        # the bath-gas identity checks refuse the write before any fraction is ever examined -- which
        # would make these cases pass for entirely the wrong reason.
        source = SOURCE_NO_KINETICS.replace(
            "    reactive=False,\n)",
            "    reactive=False,\n)\n"
            "species(\n"
            "    label='Ar',\n"
            "    structure=SMILES('[Ar]'),\n"
            "    E0=(0, 'kJ/mol'),\n"
            "    reactive=False,\n)",
        )
        assert "label='Ar'" in source, 'fixture precondition: the multi-species cases need a second bath gas'
        source_path = _write_source(tmp_path, source)
        dest_path = str(tmp_path / 'input.py')
        kwargs = dict(DEFAULT_KWARGS)
        kwargs['bath_gas'] = bath_gas
        if should_pass:
            write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **kwargs)
            assert os.path.isfile(dest_path)
        else:
            with pytest.raises(ValueError, match='bath'):
                write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **kwargs)
            assert not os.path.isfile(dest_path)


class TestDatabaseKwargsAreValidatedAgainstArkanesRealSignature:
    """
    An unknown or wrongly-typed ``database(...)`` keyword is refused at write time, not at load time.

    ``_build_database_block`` closed the injection hole in the KEY position (a key must be a plain
    identifier), which is a syntactic guarantee, and stopped there. A key can be a perfectly
    well-formed identifier and still be one Arkane's ``database()`` does not accept -- in which case
    the generated file is valid Python that raises ``TypeError: database() got an unexpected keyword
    argument`` the moment Arkane loads it, after the whole run has been set up. Arkane's real
    signature is ``arkane/input.py:83-84``; the two magic-string keywords have their own
    ``isinstance(..., list)`` checks at ``:92-109``, so a wrong type there is an InputError rather
    than a silent misconfiguration -- but either way it is knowable here.
    """

    @staticmethod
    def _write(tmp_path, database_kwargs):
        source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
        dest_path = str(tmp_path / 'input.py')
        kwargs = dict(DEFAULT_KWARGS)
        kwargs['database_kwargs'] = database_kwargs
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **kwargs)
        with open(dest_path, 'r') as f:
            return f.read()

    @pytest.mark.parametrize('database_kwargs', [
        {'thermoLibrary': ['primaryThermoLibrary']},     # singular: a plausible typo
        {'thermolibraries': ['primaryThermoLibrary']},   # wrong case
        {'seedMechanisms': ['x']},                       # a real RMG keyword, not an Arkane one
        {'kineticsFamiles': 'default'},                  # transposed letters
    ])
    def test_a_keyword_arkanes_database_does_not_accept_is_refused(self, tmp_path, database_kwargs):
        """Each of these is a valid identifier, and each makes Arkane raise TypeError at load time."""
        with pytest.raises(ValueError, match='database'):
            self._write(tmp_path, database_kwargs)

    @pytest.mark.parametrize('database_kwargs', [
        {'kineticsFamilies': 'default'},
        {'kineticsFamilies': 'all'},
        {'kineticsFamilies': 'none'},
        {'kineticsFamilies': ['H_Abstraction', 'R_Recombination']},
        {'kineticsFamilies': ['!Intra_Disproportionation']},
        {'kineticsDepositories': 'default'},
        {'kineticsDepositories': 'all'},
        {'kineticsDepositories': ['training']},
        {'thermoLibraries': []},
        {'transportLibraries': None},
        {'frequenciesLibraries': ['x']},
        {'kineticsEstimator': 'rate rules'},
    ])
    def test_every_form_arkane_documents_is_accepted(self, tmp_path, database_kwargs):
        """
        The accepted set is read off Arkane's own signature and its own validation, not guessed.

        ``kineticsFamilies`` and ``kineticsDepositories`` each accept EITHER a magic string or a
        list (``arkane/input.py:92-109``), which is why they cannot share the plain list contract of
        the library keywords. ``frequenciesLibraries`` and ``kineticsEstimator`` are accepted and
        then never used by Arkane -- refusing them would be stricter than Arkane for no benefit,
        since a caller passing them is asking for a documented no-op, not a broken run.
        """
        text = self._write(tmp_path, database_kwargs)
        assert 'database(' in text

    @pytest.mark.parametrize('database_kwargs', [
        {'kineticsFamilies': 'defualt'},                 # a typo'd magic string
        {'kineticsFamilies': 3},
        {'kineticsDepositories': 'training'},            # a list element passed bare
        {'thermoLibraries': 'primaryThermoLibrary'},     # a bare str where a list is meant
        {'thermoLibraries': [1, 2]},
        {'reactionLibraries': True},
    ])
    def test_a_known_keyword_with_a_wrong_type_is_refused(self, tmp_path, database_kwargs):
        """A known key does not vouch for its value; each of these breaks or misconfigures Arkane."""
        with pytest.raises(ValueError, match='database'):
            self._write(tmp_path, database_kwargs)


class TestTheMethodFieldHasAContractToo:
    """
    ``method`` is the fifth field with a small closed set of legal values, and it had no contract.

    An unrecognized method previously reached ``METHOD_MAP[method]`` inside
    ``t3.utils.writer.rewrite_arkane_method_line`` and raised a bare ``KeyError`` -- no field name, no
    valid set, no indication that the caller had chosen the value -- and it raised only after the
    source had been read, validated and edited. ``None`` is the case that matters in practice, because
    it is what any caller that forgets the argument passes, and ``explorer_factory`` now forwards it
    straight through.
    """

    @pytest.mark.parametrize('method', [None, 'cse', 'CSE ', 'chemically-significant eigenvalues',
                                        'MSCC', '', 3, True])
    def test_an_unrecognized_method_is_refused_by_name(self, tmp_path, method):
        """Including the lowercase and already-expanded spellings, which look plausible."""
        source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
        dest_path = str(tmp_path / 'input.py')
        kwargs = dict(DEFAULT_KWARGS)
        kwargs['method'] = method
        with pytest.raises(ValueError, match='method'):
            write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **kwargs)
        assert not os.path.isfile(dest_path)

    @pytest.mark.parametrize('method', ['CSE', 'MSC', 'RS'])
    def test_every_method_the_map_defines_is_accepted(self, tmp_path, method):
        """Over-refusal guard, read off t3.utils.writer.METHOD_MAP rather than hardcoded here."""
        from t3.utils.writer import METHOD_MAP
        assert method in METHOD_MAP, 'fixture precondition: this list must track METHOD_MAP'
        source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
        dest_path = str(tmp_path / 'input.py')
        kwargs = dict(DEFAULT_KWARGS)
        kwargs['method'] = method
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **kwargs)
        assert os.path.isfile(dest_path)


class TestAMultiSpeciesBathGasIsRefusedUnlessTheRunCanHonourIt:
    """
    Enforcing sum-to-1 on a multi-species bath gas validated a quantity that is then DISCARDED.

    Verified in RMG-Py rather than assumed, because two sources disagreed and one was wrong:
    ``arkane/explorer.py:187`` assigns our requested ``bathGas`` onto the network, and then
    ``:274``/``:302``/``:346`` call ``network.update(...)`` -- which is
    ``rmgpy/rmg/pdep.py:740``, whose body at ``:857-863`` unconditionally does
    ``self.bath_gas = {}`` followed by ``self.bath_gas[spec] = 1.0 / len(bath_gas)`` over every
    UNREACTIVE CORE SPECIES. There is no rmgmode guard on that block (an earlier probe claimed there
    was; there is not).

    So ``bathGas={'He': 0.9, 'Ar': 0.1}`` passes a sum-to-1 check, is written into the input file,
    is recorded in the manifest as what ran -- and the run uses 0.5/0.5. A warning was already
    emitted, but a warning is the wrong instrument for a request that CANNOT be satisfied: the file
    still says 0.9/0.1 and the provenance still records it.

    Equal fractions are accepted because they match what the run will actually do. The warning is
    kept for that case, because the SET can still differ: equal weights are spread over the core's
    unreactive species, which may not be the species we named.
    """

    @staticmethod
    def _write(tmp_path, bath_gas):
        source = SOURCE_NO_KINETICS.replace(
            "    reactive=False,\n)",
            "    reactive=False,\n)\n"
            "species(\n"
            "    label='Ar',\n"
            "    structure=SMILES('[Ar]'),\n"
            "    E0=(0, 'kJ/mol'),\n"
            "    reactive=False,\n)",
        )
        source_path = _write_source(tmp_path, source)
        dest_path = str(tmp_path / 'input.py')
        kwargs = dict(DEFAULT_KWARGS)
        kwargs['bath_gas'] = bath_gas
        return write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **kwargs)

    @pytest.mark.parametrize('bath_gas', [
        {'He': 0.9, 'Ar': 0.1},
        {'He': 0.75, 'Ar': 0.25},
        {'He': 0.6, 'Ar': 0.4},
    ])
    def test_unequal_multi_species_fractions_are_refused(self, tmp_path, bath_gas):
        """Each of these sums to 1 and is still guaranteed to be replaced by 0.5/0.5."""
        with pytest.raises(ValueError, match='equal'):
            self._write(tmp_path, bath_gas)

    def test_equal_multi_species_fractions_are_accepted_with_a_warning(self, tmp_path):
        """Equal weights are what the run will apply, so the request is honourable."""
        summary = self._write(tmp_path, {'He': 0.5, 'Ar': 0.5})
        assert any('unreactive' in warning or 'equal' in warning for warning in summary.warnings), \
            f'the set can still be overridden, so the warning must survive: {summary.warnings}'

    def test_a_single_species_bath_gas_is_still_accepted(self, tmp_path):
        """Over-refusal guard: the shape every real Arkane example uses."""
        summary = self._write(tmp_path, {'He': 1.0})
        assert summary is not None

    def test_three_equal_species_are_accepted_despite_binary_rounding(self, tmp_path):
        """
        1/3 + 1/3 + 1/3 is not exactly 1.0 in binary floating point, and is exactly equal-weighted.

        This is the case where an equality-based rule and a sum rule could each misfire: the sum needs
        a tolerance, and the equality check needs to compare with a tolerance too rather than by
        ``==``, or a correctly-written three-species composition is refused.
        """
        source = SOURCE_NO_KINETICS.replace(
            "    reactive=False,\n)",
            "    reactive=False,\n)\n"
            "species(\n    label='Ar',\n    structure=SMILES('[Ar]'),\n    E0=(0, 'kJ/mol'),\n    reactive=False,\n)\n"
            "species(\n    label='Ne',\n    structure=SMILES('[Ne]'),\n    E0=(0, 'kJ/mol'),\n    reactive=False,\n)",
        )
        source_path = _write_source(tmp_path, source)
        dest_path = str(tmp_path / 'input.py')
        kwargs = dict(DEFAULT_KWARGS)
        # Hand-rounded decimals, NOT ``1 / 3``. ``1 / 3`` is bit-identical to the ``1.0 / len(bath_gas)``
        # the check computes, so an ``==`` comparison passes it and the tolerance goes untested -- that
        # mutation survived until this fixture replaced it. These are the spelling a person actually
        # writes, each within 1e-6 of an equal share and summing to 1.0.
        kwargs['bath_gas'] = {'He': 0.333333, 'Ar': 0.333333, 'Ne': 0.333334}
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **kwargs)
        assert os.path.isfile(dest_path)


class TestASeedGivenAsABareStringIsRefused:
    """
    ``tuple('OH')`` is ``('O', 'H')``, and a str is a sequence, so nothing complained.

    ``seed_species='OH'`` -- a plausible thing for a caller with a single seed to write -- silently
    became a TWO-species bimolecular source channel of 'O' and 'H'. It passes the 1-or-2 length rule
    (2), and if the source network happens to declare species named 'O' and 'H' it passes the
    label-existence checks too, and Arkane then explores an entirely different reaction than the one
    requested. Where the letters do NOT happen to be declared labels, the error message names 'O' and
    'H' and sends the caller looking for species they never asked about.

    A str and a bytes object are therefore refused as the CONTAINER, and every element is required to
    be a non-empty str, since a nested list or a None element is the same class of mistake.
    """

    @staticmethod
    def _write(tmp_path, seed_species):
        source_path = _write_source(tmp_path, SOURCE_NO_KINETICS)
        dest_path = str(tmp_path / 'input.py')
        kwargs = dict(DEFAULT_KWARGS)
        kwargs['seed_species'] = seed_species
        write_arkane_explorer_input_file(source_path=source_path, dest_path=dest_path, **kwargs)

    @pytest.mark.parametrize('seed_species', [
        'methoxy',        # the realistic one: 7 characters, so it fails on LENGTH, confusingly
        'OH',             # the dangerous one: exactly 2 characters, so it passes the length rule
        b'OH',
        ['methoxy', None],
        ['methoxy', ['nested']],
        ['methoxy', ''],  # an empty label is not a label
        [1],
    ])
    def test_a_string_or_a_non_string_element_is_refused_as_a_seed(self, tmp_path, seed_species):
        """Every one of these is currently either silently wrong or wrong with a misleading message."""
        with pytest.raises(ValueError, match='seed_species'):
            self._write(tmp_path, seed_species)

    def test_a_proper_one_element_sequence_is_still_accepted(self, tmp_path):
        """Over-refusal guard: the shape every caller is supposed to use."""
        self._write(tmp_path, ('methoxy',))
        assert os.path.isfile(str(tmp_path / 'input.py'))
