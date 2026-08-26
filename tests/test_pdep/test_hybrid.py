"""
t3 tests test_pdep test_hybrid module

Tests for t3.pdep.hybrid.write_hybrid_network_input_file: rewriting an RMG-generated P-dep
network file into a hybrid Arkane input file where a caller-selected subset of transition
states are switched to by-reference QM/RRKM statmech while every other transition state and
reaction is left untouched.

Fixtures live under ``tests/data/pdep_hybrid/`` (never under a ``t3.paths[...]`` target
directory, which other tests' teardown ``shutil.rmtree``s) and reuse the real network fixture
at ``tests/data/pdep_network/iteration_1/RMG/pdep/network4_1.py`` (11 species, 5 transition
states TS1-TS5, 5 path reactions) as the source network, per the sibling
``tests/test_pdep/test_parser.py``/``tests/test_pdep/test_arkane_mesolver.py`` convention.

The QM artifact fixtures under ``tests/data/pdep_hybrid/arc_ts/`` mirror the shape ARC's
``arc/statmech/arkane.py::species_input_template`` writes for a transition state, but the
``Log(...)`` targets they reference (``tests/data/pdep_hybrid/arc_ts/logs/*.out``) are
placeholder stub text files, NOT real, Arkane-parseable quantum chemistry output: this module
only ever inspects these files with ``ast.parse``/copies their bytes, it never runs Arkane on
them.

All test outputs are written to pytest's ``tmp_path``, never into ``tests/data/``.
"""

import ast
import dataclasses
import os

import pytest

import t3.pdep.hybrid as hybrid_module
from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.hybrid import (
    HybridNetworkResult,
    QMEnergySettings,
    write_hybrid_network_input_file,
)

NETWORK_FIXTURE = os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'iteration_1', 'RMG', 'pdep', 'network4_1.py')
ARC_TS_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_hybrid', 'arc_ts')
TS1_ARTIFACT = os.path.join(ARC_TS_DIR, 'TS1.py')
TS2_ARTIFACT = os.path.join(ARC_TS_DIR, 'TS2.py')
TS_MISSING_LOG_ARTIFACT = os.path.join(ARC_TS_DIR, 'TS_missing_log.py')
TS_COLLIDING_ARTIFACT = os.path.join(ARC_TS_DIR, 'TS_colliding.py')
TS_DUP_PATH_TEXT_ARTIFACT = os.path.join(ARC_TS_DIR, 'TS_dup_path_text.py')

# A realistic ARC atom-energies dict (see tests/test_pdep/test_energy_settings.py's xl1001
# fixture), not an empty/absent one: QMEnergySettings(use_atom_corrections=True) with
# atom_energies=None is exactly the fail-open state write_hybrid_network_input_file now refuses,
# and that state is reachable from a real ARC project (atomEnergies is only ever emitted by ARC
# for gfn2/torchani model chemistries), so this fixture must carry values, not omit them.
DEFAULT_ENERGY_SETTINGS = QMEnergySettings(
    model_chemistry='wb97xd/def2tzvp',
    atom_energies={
        'Br': -2574.174533595486,
        'C': -37.84706210301937,
        'Cl': -460.1467876783656,
        'F': -99.73955550924293,
        'H': -0.5006557872395249,
        'N': -54.584995947182875,
        'O': -75.07252406126821,
        'S': -398.1105530401693,
    },
)


def ts_conformer(e0=-38.0, frequencies=(500.0, 800.0, 1200.0, 1900.0, 2200.0),
                 imaginary=-1800.0, spin_multiplicity=1, optical_isomers=1, hindered_rotors=()):
    """A DOF-normalized (vibration-only) transition-state conformer-data dict, the shape
    ``t3.runners.statmech_conformer_extract`` produces and ``write_hybrid_network_input_file``
    now consumes (in place of a by-reference artifact path)."""
    return {
        'label': 'ts', 'is_ts': True, 'E0_kJ_mol': e0,
        'frequencies_cm_1': list(frequencies), 'imaginary_frequency_cm_1': imaginary,
        'spin_multiplicity': spin_multiplicity, 'optical_isomers': optical_isomers,
        'hindered_rotors': list(hindered_rotors),
    }


def well_conformer(e0=-170.0, frequencies=(500.0, 800.0, 1200.0, 1500.0, 1900.0, 2400.0),
                   spin_multiplicity=2, optical_isomers=1, hindered_rotors=()):
    """A DOF-normalized (vibration-only) well conformer-data dict (the ``is_ts=False`` shape)."""
    return {
        'label': 'well', 'is_ts': False, 'E0_kJ_mol': e0,
        'frequencies_cm_1': list(frequencies),
        'spin_multiplicity': spin_multiplicity, 'optical_isomers': optical_isomers,
        'hindered_rotors': list(hindered_rotors),
    }


def test_refuses_unknown_ts_label(tmp_path):
    """Behavior 1: an unknown transition state label raises ValueError."""
    with pytest.raises(ValueError, match='TS_NOT_IN_NETWORK'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=str(tmp_path / 'input.py'),
            method='CSE',
            qm_transition_states={'TS_NOT_IN_NETWORK': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
        )


def test_refuses_empty_qm_transition_states(tmp_path):
    """Behavior 2: an empty qm_transition_states dict raises ValueError."""
    with pytest.raises(ValueError, match='write_arkane_network_input_file'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=str(tmp_path / 'input.py'),
            method='CSE',
            qm_transition_states={},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
        )


@pytest.mark.parametrize('model_chemistry', [None, '', '   '])
def test_refuses_missing_or_blank_model_chemistry(tmp_path, model_chemistry):
    """Behavior 4: a missing/blank model_chemistry raises ValueError."""
    energy_settings = QMEnergySettings(model_chemistry=model_chemistry)
    with pytest.raises(ValueError, match='model_chemistry'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=str(tmp_path / 'input.py'),
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=energy_settings,
        )


def test_refuses_use_bond_corrections_true_unconditionally(tmp_path):
    """Item 1: use_bond_corrections=True is refused unconditionally, regardless of
    bond_correction_type. Arkane's input DSL has no directive to pin BAC VALUES -- only
    useBondCorrections (bool) and bondCorrectionType ('p'/'m', which SCHEME to apply); the actual
    numeric corrections are looked up from Arkane's own database at solve time (see
    arkane/input.py's local_context.get('atomEnergies', ...) being the only corrections-VALUES
    directive that exists). So turning this on can silently apply a BAC that was fit against a
    different reference than the RMG wells around a QM'd TS, with no way for this file to pin what
    was actually applied."""
    energy_settings = QMEnergySettings(
        model_chemistry='wb97xd/def2tzvp', use_bond_corrections=True, bond_correction_type='p',
        atom_energies={'C': -37.79, 'H': -0.5})
    with pytest.raises(ValueError, match='use_bond_corrections'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=str(tmp_path / 'input.py'),
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=energy_settings,
        )


def test_refuses_bond_correction_type_set_when_bond_corrections_disabled(tmp_path):
    """Item 1: bond_correction_type set while use_bond_corrections is False raises ValueError,
    since it would silently have no effect."""
    energy_settings = QMEnergySettings(
        model_chemistry='wb97xd/def2tzvp', use_bond_corrections=False, bond_correction_type='p',
        atom_energies={'C': -37.79, 'H': -0.5})
    with pytest.raises(ValueError, match='bond_correction_type'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=str(tmp_path / 'input.py'),
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=energy_settings,
        )


def test_refuses_use_atom_corrections_true_with_atom_energies_none(tmp_path):
    """Item 5: use_atom_corrections=True (the default) with atom_energies left None raises
    ValueError. This is the currently-live fail-open: _build_energy_header only emits
    atomEnergies = ... when atom_energies is not None, so a caller could otherwise turn atom
    corrections ON in the generated file without pinning the VALUES Arkane will actually look up
    -- exactly the gap that leaves a QM'd TS's E0 off the RMG wells' enthalpy-of-formation scale."""
    energy_settings = QMEnergySettings(model_chemistry='wb97xd/def2tzvp', atom_energies=None)
    with pytest.raises(ValueError, match='atom_energies'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=str(tmp_path / 'input.py'),
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=energy_settings,
        )


@pytest.mark.parametrize('bad_model_chemistry', [
    'wb97xd/def2tzvp\nuseAtomCorrections = False',
    "wb97xd\\'/def2tzvp",
    'wb97xd"/def2tzvp',
    'wb97xd\\/def2tzvp',
])
def test_rejects_model_chemistry_with_injection_characters(tmp_path, bad_model_chemistry):
    """Item 2: model_chemistry containing a newline, quote or backslash is rejected up front,
    rather than being raw-interpolated into the generated (executable) Arkane input."""
    energy_settings = QMEnergySettings(model_chemistry=bad_model_chemistry)
    with pytest.raises(ValueError, match='model_chemistry'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=str(tmp_path / 'input.py'),
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=energy_settings,
        )


@pytest.mark.parametrize('bad_tunneling', [
    "Eckart\nkinetics = Arrhenius(A=(1,'s^-1'))",
    "Eckart'",
    'Eckart"',
    'Eckart\\',
])
def test_rejects_tunneling_with_injection_characters(tmp_path, bad_tunneling):
    """Item 2: tunneling containing a newline, quote or backslash is rejected up front, rather
    than being raw-interpolated into the generated (executable) Arkane input."""
    energy_settings = QMEnergySettings(
        model_chemistry='wb97xd/def2tzvp', tunneling=bad_tunneling, atom_energies={'C': -37.79, 'H': -0.5})
    with pytest.raises(ValueError, match='tunneling'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=str(tmp_path / 'input.py'),
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=energy_settings,
        )


def test_refuses_use_atom_corrections_false(tmp_path):
    """Part 2: use_atom_corrections=False raises ValueError, since that directive silently
    disables Arkane's atom energy corrections."""
    energy_settings = QMEnergySettings(model_chemistry='wb97xd/def2tzvp', use_atom_corrections=False)
    with pytest.raises(ValueError, match='use_atom_corrections'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=str(tmp_path / 'input.py'),
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=energy_settings,
        )


def test_refuses_use_atom_corrections_false_even_with_atom_energies_present(tmp_path):
    """Part 2, isolated: use_atom_corrections=False must be refused on its own, not merely as a
    side effect of atom_energies also being unset. Supplying a non-None atom_energies here means
    the ONLY guard that can fire is the use_atom_corrections=False check itself; without this test,
    that guard could be deleted entirely and test_refuses_use_atom_corrections_false above would
    still pass (via the atom_energies-is-None guard raising instead), silently dissolving this
    refusal."""
    energy_settings = QMEnergySettings(model_chemistry='wb97xd/def2tzvp', use_atom_corrections=False,
                                       atom_energies={'H': -0.5})
    with pytest.raises(ValueError, match='use_atom_corrections'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=str(tmp_path / 'input.py'),
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=energy_settings,
        )


def test_emits_use_atom_corrections_true_by_default(tmp_path):
    """Part 2: the header unconditionally emits useAtomCorrections = True (the only value that
    ever reaches this point, since False is rejected up front)."""
    dest_path = str(tmp_path / 'input.py')
    write_hybrid_network_input_file(
        source_path=NETWORK_FIXTURE,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    assert 'useAtomCorrections = True' in text


@pytest.mark.parametrize('model_chemistry_expr', [
    "LevelOfTheory(method='wb97xd2023', basis='def2tzvp', software='gaussian')",
    "CompositeLevelOfTheory(freq=LevelOfTheory(method='wb97xd', basis='def2tzvp', software='gaussian'), "
    "energy=LevelOfTheory(method='dlpnoccsd(t)f122023', basis='ccpvtzf12', software='orca'))",
    # software_version is the one genuinely non-str LevelOfTheory field (Union[int, float, str]
    # in arkane/encorr/modelchem.py); it must be allowed, not blanket-rejected as non-str.
    "LevelOfTheory(method='wb97xd', basis='def2tzvp', software='gaussian', software_version=16)",
    "LevelOfTheory(method='wb97xd', basis='def2tzvp', software='gaussian', software_version=16.0)",
])
def test_accepts_and_renders_bare_level_of_theory_model_chemistry(tmp_path, model_chemistry_expr):
    """Part 2: a real ARC-shaped LevelOfTheory(...)/CompositeLevelOfTheory(...) object expression
    is accepted (it legitimately contains quote characters that _validate_no_injection_chars alone
    would reject) and rendered bare (unquoted), not as a quoted string containing that syntax."""
    dest_path = str(tmp_path / 'input.py')
    energy_settings = QMEnergySettings(model_chemistry=model_chemistry_expr, atom_energies={'C': -37.79, 'H': -0.5})
    write_hybrid_network_input_file(
        source_path=NETWORK_FIXTURE,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=energy_settings,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    assert f'modelChemistry = {model_chemistry_expr}' in text
    assert f'modelChemistry = "{model_chemistry_expr}"' not in text
    # The generated file must itself remain valid, executable-shaped Arkane syntax.
    ast.parse(text)


@pytest.mark.parametrize('bad_model_chemistry_expr', [
    "LevelOfTheory('wb97xd2023', basis='def2tzvp')",  # positional arg
    "LevelOfTheory(**{'method': 'wb97xd2023'})",  # **kwargs, dict-valued
    "LevelOfTheory(**'not_a_mapping')",  # **kwargs, constant-valued: isolates the `kw.arg is None`
                                         # guard, since kw.value here IS an ast.Constant and would
                                         # otherwise slip past the isinstance(Constant) check
    "LevelOfTheory(method=some_name)",  # non-literal keyword value
    "CompositeLevelOfTheory(LevelOfTheory(method='a'), energy=LevelOfTheory(method='b'))",  # positional arg
    "CompositeLevelOfTheory(freq=LevelOfTheory(method='a'), sp=LevelOfTheory(method='b'))",  # bad keyword
    "CompositeLevelOfTheory(freq='not_a_level_of_theory_call', energy=LevelOfTheory(method='b'))",
    # round 22 P1 #4: AST-structural allowlist accepted any keyword NAME and any ast.Constant
    # VALUE type, and a vacuous CompositeLevelOfTheory() with no freq/energy at all.
    "LevelOfTheory(__class__='x')",  # not a real LevelOfTheory field name
    "LevelOfTheory(method='a', extra_bogus_field='x')",  # not a real LevelOfTheory field name
    "LevelOfTheory(method=True)",  # method must be str, not bool
    "LevelOfTheory(method=...)",  # method must be str, not Ellipsis
    "LevelOfTheory(method=b'x')",  # method must be str, not bytes
    "LevelOfTheory(method=3)",  # method must be str, not int
    "LevelOfTheory(basis='def2tzvp')",  # method is a required LevelOfTheory field
    "CompositeLevelOfTheory()",  # vacuous: no freq/energy at all
    "CompositeLevelOfTheory(freq=LevelOfTheory(method='a'))",  # missing required energy keyword
    "CompositeLevelOfTheory(energy=LevelOfTheory(method='b'))",  # missing required freq keyword
    "CompositeLevelOfTheory(freq=LevelOfTheory(method='a'), energy=LevelOfTheory(method='b'), "
    "extra=LevelOfTheory(method='c'))",  # not a real CompositeLevelOfTheory field name
    "LevelOfTheory(method='a', method='b')",  # duplicate keyword: ast.parse alone permits this
                                               # even though CPython's compile() step rejects it
    "CompositeLevelOfTheory(freq=LevelOfTheory(method='a'), freq=LevelOfTheory(method='b'), "
    "energy=LevelOfTheory(method='c'))",  # duplicate keyword
])
def test_rejects_structurally_invalid_model_chemistry_object_expression(tmp_path, bad_model_chemistry_expr):
    """Part 2: a LevelOfTheory(...)/CompositeLevelOfTheory(...)-shaped model_chemistry that fails
    structural validation (positional args, **kwargs, unexpected keyword, non-literal/non-nested-
    call keyword value) is rejected rather than silently interpolated bare into the generated
    (executable) Arkane input."""
    energy_settings = QMEnergySettings(model_chemistry=bad_model_chemistry_expr)
    with pytest.raises(ValueError, match='model_chemistry'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=str(tmp_path / 'input.py'),
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=energy_settings,
        )


def test_replaces_qm_ts_block_and_leaves_others_byte_identical(tmp_path):
    """Behaviors 5 and 6: the QM'd TS's block is replaced by-reference; all other TS blocks
    (and the rest of the file) are left byte-identical."""
    dest_path = str(tmp_path / 'input.py')
    result = write_hybrid_network_input_file(
        source_path=NETWORK_FIXTURE,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
    )
    assert isinstance(result, HybridNetworkResult)
    assert result.qm_ts_labels == ('TS1',)
    assert result.ilt_ts_labels == ('TS2', 'TS3', 'TS4', 'TS5')

    with open(dest_path, 'r') as f:
        text = f.read()

    # New inline form for TS1: the whole conformer is spliced in directly, no by-reference
    # pointer into a vendored qm/ directory.
    assert "label = 'TS1'" in text
    assert "E0 = (-38.0,'kJ/mol')" in text
    assert 'HarmonicOscillator(' in text
    assert "frequency = (-1800.0,'cm^-1')" in text
    assert 'qm/' not in text
    assert 'E0 = (505.143' not in text

    # TS2-TS5's inline blocks are untouched, byte-identical to the source.
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    for label, e0 in (('TS2', '(205.409'), ('TS3', '(170.919'), ('TS4', '(469.777'), ('TS5', '(210.753')):
        block_start = source_text.index(f"label = '{label}'")
        block_end = source_text.index(')', block_start) + 1
        original_block = source_text[block_start:block_end]
        assert original_block in text
        assert e0 in original_block


def test_qm_reaction_drops_kinetics_and_adds_tunneling(tmp_path):
    """Behavior 7: dropping kinetics= only from the QM'd reaction, with tunneling added; other
    reactions are untouched."""
    dest_path = str(tmp_path / 'input.py')
    write_hybrid_network_input_file(
        source_path=NETWORK_FIXTURE,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
    )
    with open(dest_path, 'r') as f:
        text = f.read()

    reaction1_start = text.index("label = 'reaction1'")
    reaction1_end = text.index('\n)\n', reaction1_start)
    reaction1_block = text[reaction1_start:reaction1_end]
    assert 'kinetics' not in reaction1_block
    assert "tunneling = 'Eckart'" in reaction1_block

    # reaction2 (TS2, not QM'd) keeps its kinetics untouched.
    reaction2_start = text.index("label = 'reaction2'")
    reaction2_end = text.index('\n)\n', reaction2_start)
    reaction2_block = text[reaction2_start:reaction2_end]
    assert 'kinetics = Arrhenius' in reaction2_block
    assert 'tunneling' not in reaction2_block


def test_qm_reaction_omits_tunneling_when_not_configured(tmp_path):
    """tunneling=None on QMEnergySettings means no tunneling = ... is written."""
    dest_path = str(tmp_path / 'input.py')
    energy_settings = QMEnergySettings(
        model_chemistry='wb97xd/def2tzvp', tunneling=None, atom_energies={'C': -37.79, 'H': -0.5})
    write_hybrid_network_input_file(
        source_path=NETWORK_FIXTURE,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=energy_settings,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    reaction1_start = text.index("label = 'reaction1'")
    reaction1_end = text.index('\n)\n', reaction1_start)
    reaction1_block = text[reaction1_start:reaction1_end]
    assert 'kinetics' not in reaction1_block
    assert 'tunneling' not in reaction1_block


# Item 3: ``_kinetics_removal_edit`` must locate the splice span structurally (via AST keyword/
# call node positions), not by regex-hunting for a trailing comma, so it stays correct however
# ``kinetics`` is positioned among a reaction(...) call's other keywords and regardless of a
# trailing inline comment. The original ``reaction1`` block (kinetics last, its own line, no
# inline comment) is reused as the substring to replace, so the rest of the fixture network
# (other species/transition states/reactions) is left completely intact.
_ORIGINAL_REACTION1_BLOCK = (
    "reaction(\n"
    "    label = 'reaction1',\n"
    "    reactants = ['CH2(S)(53)', 'C3rad(4)'],\n"
    "    products = ['C4rad(5)'],\n"
    "    transitionState = 'TS1',\n"
    "    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), "
    "comment=\"\"\"Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]\n"
    "Euclidian distance = 1.0\n"
    "Multiplied by reaction path degeneracy 3.0\n"
    "family: 1,2_Insertion_carbene\n"
    "Ea raised from -1.5 to 0 kJ/mol.\"\"\"),\n"
    ")\n"
)


def _write_network_with_reaction1_variant(tmp_path, replacement_block):
    """Write a copy of NETWORK_FIXTURE with its reaction1 block swapped for ``replacement_block``."""
    with open(NETWORK_FIXTURE, 'r') as f:
        text = f.read()
    assert _ORIGINAL_REACTION1_BLOCK in text, 'fixture reaction1 block text has drifted'
    text = text.replace(_ORIGINAL_REACTION1_BLOCK, replacement_block)
    source_path = str(tmp_path / 'network_variant.py')
    with open(source_path, 'w') as f:
        f.write(text)
    return source_path


def test_kinetics_removal_handles_kinetics_as_first_keyword(tmp_path):
    """Item 3: kinetics as the FIRST keyword (followed by other keywords on later lines) must not
    leave a dangling/double comma behind -- the generated file must parse as valid Python."""
    variant = (
        "reaction(\n"
        "    kinetics = Arrhenius(A=(1,'s^-1'), n=0, Ea=(0,'kJ/mol')),\n"
        "    label = 'reaction1',\n"
        "    reactants = ['CH2(S)(53)', 'C3rad(4)'],\n"
        "    products = ['C4rad(5)'],\n"
        "    transitionState = 'TS1',\n"
        ")\n"
    )
    source_path = _write_network_with_reaction1_variant(tmp_path, variant)
    dest_path = str(tmp_path / 'input.py')
    write_hybrid_network_input_file(
        source_path=source_path,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    ast.parse(text)  # must not raise -- this is the crux of item 3
    reaction1_start = text.index("label = 'reaction1'")
    reaction1_end = text.index('\n)\n', reaction1_start)
    reaction1_block = text[reaction1_start:reaction1_end]
    assert 'kinetics' not in reaction1_block


def test_kinetics_removal_handles_trailing_inline_comment(tmp_path):
    """Item 3: a ``# ...`` comment trailing the comma after kinetics' value (on the same line) must
    not defeat the splice -- the old trailing-comma regex required an immediate newline after the
    comma and left a dangling comma (invalid Python) when a comment came between them."""
    variant = (
        "reaction(\n"
        "    label = 'reaction1',\n"
        "    reactants = ['CH2(S)(53)', 'C3rad(4)'],\n"
        "    products = ['C4rad(5)'],\n"
        "    transitionState = 'TS1',\n"
        "    kinetics = Arrhenius(A=(1,'s^-1'), n=0, Ea=(0,'kJ/mol')),  # RMG estimate\n"
        ")\n"
    )
    source_path = _write_network_with_reaction1_variant(tmp_path, variant)
    dest_path = str(tmp_path / 'input.py')
    write_hybrid_network_input_file(
        source_path=source_path,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    ast.parse(text)  # must not raise -- this is the crux of item 3
    reaction1_start = text.index("label = 'reaction1'")
    reaction1_end = text.index('\n)\n', reaction1_start)
    reaction1_block = text[reaction1_start:reaction1_end]
    assert 'kinetics' not in reaction1_block


def test_kinetics_removal_handles_single_line_reaction_call(tmp_path):
    """Item 3: a single-line reaction(...) call, with kinetics in the middle, must not leave a
    double comma (invalid Python) behind."""
    variant = (
        "reaction(label = 'reaction1', reactants = ['CH2(S)(53)', 'C3rad(4)'], "
        "products = ['C4rad(5)'], kinetics = Arrhenius(A=(1,'s^-1'), n=0, Ea=(0,'kJ/mol')), "
        "transitionState = 'TS1')\n"
    )
    source_path = _write_network_with_reaction1_variant(tmp_path, variant)
    dest_path = str(tmp_path / 'input.py')
    write_hybrid_network_input_file(
        source_path=source_path,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    ast.parse(text)  # must not raise -- this is the crux of item 3
    reaction1_start = text.index("label = 'reaction1'")
    reaction1_end = text.index(")\n", reaction1_start)
    reaction1_block = text[reaction1_start:reaction1_end]
    assert 'kinetics' not in reaction1_block
    assert ', ,' not in reaction1_block
    assert ',,' not in reaction1_block


def test_kinetics_removal_handles_non_ascii_text_earlier_on_the_same_line(tmp_path):
    """Item 4: ``ast`` node ``col_offset`` values are UTF-8 BYTE offsets, not character indices.
    A non-ASCII character earlier on the same line as an edited keyword (here, a label containing
    'e' with an acute accent, encoded as 2 UTF-8 bytes but 1 character) must not throw off the
    splice position -- if the byte offset is applied directly as a character index, the computed
    position lands one character too far into the line, corrupting the edit."""
    variant = (
        "reaction(label = 'é', reactants = ['CH2(S)(53)', 'C3rad(4)'], "
        "products = ['C4rad(5)'], kinetics = Arrhenius(A=(1,'s^-1'), n=0, Ea=(0,'kJ/mol')), "
        "transitionState = 'TS1')\n"
    )
    source_path = _write_network_with_reaction1_variant(tmp_path, variant)
    dest_path = str(tmp_path / 'input.py')
    write_hybrid_network_input_file(
        source_path=source_path,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    ast.parse(text)  # must not raise -- this is the crux of item 4
    reaction1_start = text.index("label = 'é'")
    reaction1_end = text.index(")\n", reaction1_start)
    reaction1_block = text[reaction1_start:reaction1_end]
    assert 'kinetics' not in reaction1_block
    assert 'transitionState' in reaction1_block


def test_injects_energy_header_before_first_species_and_never_disables_atom_corrections(tmp_path):
    """Behavior 8: the energy header is injected before the first species(...) block, and
    useAtomCorrections = False is never emitted."""
    dest_path = str(tmp_path / 'input.py')
    energy_settings = QMEnergySettings(
        model_chemistry='wb97xd/def2tzvp',
        frequency_scale_factor=0.988,
        use_hindered_rotors=True,
        use_bond_corrections=False,
        atom_energies={'C': -37.79, 'H': -0.5}
    )
    write_hybrid_network_input_file(
        source_path=NETWORK_FIXTURE,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=energy_settings,
    )
    with open(dest_path, 'r') as f:
        text = f.read()

    assert 'useAtomCorrections = False' not in text
    header_idx = text.index('modelChemistry = "wb97xd/def2tzvp"')
    first_species_idx = text.index("species(\n    label = 'C4rad(5)'")
    assert header_idx < first_species_idx
    assert 'frequencyScaleFactor = 0.988' in text
    assert 'useHinderedRotors = True' in text
    assert 'useBondCorrections = False' in text
    assert "atomEnergies" in text


def test_dest_path_write_is_atomic_no_torn_file_on_crash(tmp_path, monkeypatch):
    """The 'ARC finalization marker' (t3.main.T3._mark_arc_finalization_complete) is written via
    tempfile.mkstemp + fsync + os.replace precisely so a crash mid-write can never leave behind a
    marker that certifies unfinished work. Before this fix, write_hybrid_network_input_file wrote
    dest_path with a plain open(dest_path, 'w'); f.write(text) -- a crash partway through that
    write (e.g. disk full, power loss, OOM-kill) can leave a TRUNCATED input.py at dest_path. On
    the very next pass, _mark_arc_finalization_complete then runs and marks the iteration fully
    finalized, so check_arc_finalization_complete() reports success while a downstream consumer
    reads a torn Arkane input as if T3 stood behind it -- the exact fail-open shape the marker
    itself was built to refuse.

    This test simulates that crash by making the write's flush() raise partway through, and
    asserts:
      1. No partial dest_path is left behind when dest_path did not exist before the call.
      2. A pre-existing GOOD dest_path survives a failed rewrite completely intact (the atomic
         write must stage-then-replace, never truncate-in-place)."""
    dest_path = str(tmp_path / 'network' / 'input.py')
    os.makedirs(os.path.dirname(dest_path))
    pre_existing_content = "# a previously-written, good hybrid network input\n"
    with open(dest_path, 'w') as f:
        f.write(pre_existing_content)

    real_flush = os.fdopen

    class _CrashingFile:
        """Wraps a real file object but raises on flush(), simulating a crash mid-write after
        some bytes may already be buffered/written but before the write is complete."""

        def __init__(self, real_file):
            self._real_file = real_file

        def write(self, data):
            return self._real_file.write(data)

        def flush(self):
            raise OSError('simulated crash mid-write of dest_path')

        def __enter__(self):
            return self

        def __exit__(self, *exc_info):
            self._real_file.close()
            return False

    def _crashing_fdopen(fd, *args, **kwargs):
        # Only intercept writes destined for a temp file staged next to dest_path (this module's
        # own atomic-write staging file); anything else (e.g. vendoring's own temp files) must
        # behave normally so only the dest_path write itself is made to fail.
        real_file = real_flush(fd, *args, **kwargs)
        return _CrashingFile(real_file)

    monkeypatch.setattr(hybrid_module.os, 'fdopen', _crashing_fdopen)

    with pytest.raises(OSError, match='simulated crash mid-write of dest_path'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
        )

    # No leftover staging debris next to dest_path (would trip
    # T3._prune_stale_pdep_hybrid_outputs, which raises on any non-plain-directory entry it finds
    # directly under the 'PDep hybrid' root -- but a staging file lives ONE level deeper, inside a
    # network's own directory, so it can never be such a top-level entry; still, it must not be
    # left behind at all).
    # No '.hybrid-input-*' staging file may survive a failed write (the network is now spliced
    # inline with no vendoring step, so nothing else is created alongside dest_path either).
    dest_dir = os.path.dirname(dest_path)
    leftover_staging_files = [name for name in os.listdir(dest_dir) if name.startswith('.hybrid-input-')]
    assert leftover_staging_files == [], (
        f"Staging debris left behind after a failed write: {leftover_staging_files}")

    # The pre-existing good file must survive a failed rewrite completely intact -- a torn write
    # must never truncate or corrupt what was already there.
    with open(dest_path, 'r') as f:
        assert f.read() == pre_existing_content


def test_multiple_qm_transition_states(tmp_path):
    """Two QM'd transition states are both switched, and the third is left as ILT."""
    dest_path = str(tmp_path / 'input.py')
    result = write_hybrid_network_input_file(
        source_path=NETWORK_FIXTURE,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer(e0=-38.0), 'TS2': ts_conformer(e0=-50.0)},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
    )
    assert result.qm_ts_labels == ('TS1', 'TS2')
    assert result.ilt_ts_labels == ('TS3', 'TS4', 'TS5')
    with open(dest_path, 'r') as f:
        text = f.read()
    assert "label = 'TS1'" in text
    assert "E0 = (-38.0,'kJ/mol')" in text
    assert "label = 'TS2'" in text
    assert "E0 = (-50.0,'kJ/mol')" in text
    assert 'qm/' not in text


def test_rewrites_method_line_same_as_write_arkane_network_input_file(tmp_path):
    """Behavior 11: the method = ... line is rewritten exactly as
    write_arkane_network_input_file already does (reusing its rewrite_arkane_method_line
    helper)."""
    dest_path = str(tmp_path / 'input.py')
    write_hybrid_network_input_file(
        source_path=NETWORK_FIXTURE,
        dest_path=dest_path,
        method='RS',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    assert "method = 'reservoir state'" in text
    assert "method = 'modified strong collision'" not in text


def test_rewrites_unspaced_method_line(tmp_path):
    """B4(b): an unspaced method='...' line (no space around '=') must still be rewritten by
    _rewrite_method_and_sensitivity's own call site, which historically guarded on the literal
    substring 'method = ' (with a space) and so never fired for this shape."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    assert "    method = 'modified strong collision',\n" in source_text
    unspaced_text = source_text.replace(
        "    method = 'modified strong collision',\n", "    method='modified strong collision',\n")
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(unspaced_text)

    dest_path = str(tmp_path / 'input.py')
    write_hybrid_network_input_file(
        source_path=source_path,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    assert 'chemically-significant eigenvalues' in text
    assert 'modified strong collision' not in text


def test_raises_if_no_method_line(tmp_path):
    """B4(c): a source network with no 'method = ...' line at all must RAISE, never silently
    write a hybrid network file whose method was never actually rewritten."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    assert "    method = 'modified strong collision',\n" in source_text
    no_method_text = source_text.replace("    method = 'modified strong collision',\n", '')
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(no_method_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='method'):
        write_hybrid_network_input_file(
            source_path=source_path,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
        )


def test_clamps_pressure_dependence_and_sensitivity_conditions_tmax_to_the_species_thermo_ceiling(tmp_path):
    """_rewrite_method_and_sensitivity gets the same species-thermo-ceiling clamp as
    write_arkane_network_input_file (t3/utils/writer.py): a hybrid network is still text that
    contains the source network's intact species blocks, so it is just as exposed to RMG writing a
    pressureDependence(... Tmax=(X,'K') ...) block that asks for a temperature no species' NASA
    thermo actually supports -- standalone Arkane refuses that with 'No valid NASA polynomial at
    temperature X K.' even though RMG itself tolerated the extrapolation at generation time. Here
    NETWORK_FIXTURE's own 'C4rad(5)' isomer species is narrowed from a 5000 K outer NASA ceiling to
    3000 K, and the pdep block's requested Tmax is raised from 2100 K to 3200 K, so the network-wide
    ceiling (the MIN over species, 3000 K here since every other species is still valid to 5000 K
    or 6000 K) is below the requested Tmax and a clamp must fire in the WRITTEN hybrid file, in both
    the pressureDependence(...) block's own Tmax line and the injected sensitivity_conditions'
    high-T entries."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    assert "Tmax=(5000,'K'), E0=(63.0573,'kJ/mol')" in source_text
    assert "    Tmax = (2100,'K'),\n" in source_text
    narrowed_text = source_text.replace(
        "Tmax=(5000,'K'), E0=(63.0573,'kJ/mol')", "Tmax=(3000,'K'), E0=(63.0573,'kJ/mol')")
    narrowed_text = narrowed_text.replace("    Tmax = (2100,'K'),\n", "    Tmax = (3200,'K'),\n")
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(narrowed_text)

    dest_path = str(tmp_path / 'input.py')
    write_hybrid_network_input_file(
        source_path=source_path,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
        sensitivity=True,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    assert "Tmax = (3200,'K')" not in text
    assert "Tmax = (3000,'K')" in text
    assert "sensitivity_conditions" in text
    assert "(3200, 'K')" not in text
    assert "(3000, 'K')" in text


def test_raises_when_clamping_the_hybrid_networks_tmax_would_leave_it_at_or_below_tmin(tmp_path):
    """Same degenerate-range guard as write_arkane_network_input_file: if the narrowest species
    thermo ceiling (3000 K here) is at or below the pdep block's own Tmin (also raised to 3000 K in
    this fixture), clamping Tmax down to that ceiling would leave no valid temperature point to
    solve at all, so _rewrite_method_and_sensitivity must raise a ValueError naming both numbers
    rather than silently writing a degenerate T range."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    assert "Tmax=(5000,'K'), E0=(63.0573,'kJ/mol')" in source_text
    assert "    Tmin = (300,'K'),\n" in source_text
    assert "    Tmax = (2100,'K'),\n" in source_text
    degenerate_text = source_text.replace(
        "Tmax=(5000,'K'), E0=(63.0573,'kJ/mol')", "Tmax=(3000,'K'), E0=(63.0573,'kJ/mol')")
    degenerate_text = degenerate_text.replace("    Tmin = (300,'K'),\n", "    Tmin = (3000,'K'),\n")
    degenerate_text = degenerate_text.replace("    Tmax = (2100,'K'),\n", "    Tmax = (3200,'K'),\n")
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(degenerate_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='3000'):
        write_hybrid_network_input_file(
            source_path=source_path,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
        )


def test_drops_an_out_of_range_tlist_line_from_the_hybrid_network(tmp_path):
    """Arkane's pdep.py only calls generate_T_list() (building the grid from Tmin/Tmax/Tcount)
    when self.Tlist is None; if an explicit Tlist is present it wins outright and Tmin/Tmax are
    ignored. NETWORK_FIXTURE's pdep block already carries a real RMG-style Tlist
    (Tlist = ([302.491,...,1985.54],'K')); clamping the Tmax line alone (as the earlier fix did)
    would be a no-op once that Tlist's own high entry exceeds the narrowed thermo ceiling, since
    Arkane would still solve at the stale Tlist value and still fail. Here the ceiling is narrowed
    to 1500 K (below the Tlist's own 1985.54 K entry, and below the also-raised pdep Tmax of
    3200 K) so both the Tmax clamp AND the Tlist drop must fire in the written hybrid file, while
    Tmin/Tmax/Tcount remain."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    assert "Tmax=(5000,'K'), E0=(63.0573,'kJ/mol')" in source_text
    assert "    Tmax = (2100,'K'),\n" in source_text
    assert "Tlist = ([302.491,323.355,370.585,457.988,614.983,900.017,1394.8,1985.54],'K')" in source_text
    narrowed_text = source_text.replace(
        "Tmax=(5000,'K'), E0=(63.0573,'kJ/mol')", "Tmax=(1500,'K'), E0=(63.0573,'kJ/mol')")
    narrowed_text = narrowed_text.replace("    Tmax = (2100,'K'),\n", "    Tmax = (3200,'K'),\n")
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(narrowed_text)

    dest_path = str(tmp_path / 'input.py')
    write_hybrid_network_input_file(
        source_path=source_path,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
        sensitivity=True,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    assert 'Tlist' not in text
    assert "Tmin = (300,'K')" in text
    assert "Tmax = (1500,'K')" in text
    assert 'Tcount = 8,' in text


def test_leaves_an_in_range_tlist_line_untouched_in_the_hybrid_network(tmp_path):
    """When every Tlist entry is already within the (possibly clamped) species thermo ceiling, the
    line must be left byte-for-byte as written by RMG -- there is nothing to fix, and dropping it
    anyway would gratuitously throw away RMG's own precomputed grid."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    assert "Tlist = ([302.491,323.355,370.585,457.988,614.983,900.017,1394.8,1985.54],'K')" in source_text

    dest_path = str(tmp_path / 'input.py')
    write_hybrid_network_input_file(
        source_path=NETWORK_FIXTURE,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
        sensitivity=True,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    assert "Tlist = ([302.491,323.355,370.585,457.988,614.983,900.017,1394.8,1985.54],'K')" in text


def test_drops_tlist_even_when_hybrid_tmax_already_within_ceiling(tmp_path):
    """The drop decision must be made on the Tlist entries themselves, NOT on whether the Tmax
    clamp fired: a file can carry a Tmax already at/below the ceiling while Tlist still contains a
    stale higher entry -- exactly the shape a partial or earlier clamp would leave behind -- and
    that file must still be corrected."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    already_clamped_text = source_text.replace(
        "Tmax=(5000,'K'), E0=(63.0573,'kJ/mol')", "Tmax=(1500,'K'), E0=(63.0573,'kJ/mol')")
    already_clamped_text = already_clamped_text.replace(
        "    Tmax = (2100,'K'),\n", "    Tmax = (1500,'K'),\n")
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(already_clamped_text)

    dest_path = str(tmp_path / 'input.py')
    write_hybrid_network_input_file(
        source_path=source_path,
        dest_path=dest_path,
        method='CSE',
        qm_transition_states={'TS1': ts_conformer()},
        energy_settings=DEFAULT_ENERGY_SETTINGS,
        sensitivity=True,
    )
    with open(dest_path, 'r') as f:
        text = f.read()
    assert 'Tlist' not in text
    assert "Tmax = (1500,'K')" in text


def test_raises_on_out_of_range_tlist_without_tcount_in_the_hybrid_network(tmp_path):
    """Dropping an out-of-range Tlist is only safe if Tmin, Tmax and Tcount were all found in the
    block, since Arkane needs all three to regenerate the grid. If Tcount is missing, silently
    leaving the stale out-of-range Tlist would reproduce the bug being fixed, so this must raise a
    ValueError naming what's missing rather than writing a file Arkane will still fail on."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    narrowed_text = source_text.replace(
        "Tmax=(5000,'K'), E0=(63.0573,'kJ/mol')", "Tmax=(1500,'K'), E0=(63.0573,'kJ/mol')")
    narrowed_text = narrowed_text.replace("    Tmax = (2100,'K'),\n", "    Tmax = (3200,'K'),\n")
    assert "    Tcount = 8,\n" in narrowed_text
    no_tcount_text = narrowed_text.replace('    Tcount = 8,\n', '')
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(no_tcount_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='Tcount'):
        write_hybrid_network_input_file(
            source_path=source_path,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
            sensitivity=True,
        )


def test_propagates_a_non_unparseable_valueerror_from_network_thermo_t_max(tmp_path):
    """Finding C: _rewrite_method_and_sensitivity's try/except around network_thermo_t_max(text)
    exists ONLY to defer the 'text does not parse as Python at all' case to
    write_hybrid_network_input_file's own ast.parse(text) self-check further down the pipeline
    (see the comment above the try/except in hybrid.py). It must NOT also swallow a ValueError
    network_thermo_t_max raises for an entirely different reason -- e.g. a nested NASA(...) call
    (inside a species(...)'s thermo= keyword) that uses a positional argument, which
    network_thermo.py's _call_keywords refuses since a positional argument's name cannot be
    recovered from the AST. Before this fix, a blanket 'except ValueError: thermo_t_max = None'
    silently disabled BOTH the Tmax clamp AND the Tlist drop for such a file instead of surfacing
    the real defect.

    Deliberately targets the NASA(...) call rather than the outer species(...) call: hybrid.py's
    own pre-flight walk calls parse_pdep_network_file(source_path) unconditionally BEFORE
    _rewrite_method_and_sensitivity runs, and that parser already refuses a positional argument on
    any RECOGNIZED_TOP_LEVEL_CALLS member (species/transitionState/reaction/network/
    pressureDependence) -- so a positional arg on the top-level species(...) call itself would be
    caught earlier, by a different mechanism entirely, and this test would pass even with Finding
    C unfixed (as a first draft of this test discovered via mutation testing). The parser does NOT
    recurse into a NASA(...) call nested inside thermo=..., so a positional argument there reaches
    network_thermo_t_max's _call_keywords -- and only there -- uninterrupted.

    Mutates NETWORK_FIXTURE's 'H(34)' species(...) call so its thermo=NASA(...) passes
    `polynomials` positionally instead of as `polynomials=...`, which still parses as valid Python
    (ast.parse succeeds, and parse_pdep_network_file's top-level-only check does not object) but
    must make network_thermo_t_max -- and therefore write_hybrid_network_input_file -- raise, not
    silently degrade to 'no ceiling'."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    needle = (
        "thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], "
        "Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], "
        "Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), "
        "E0=(211.8,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), "
        "label=\"\"\"H\"\"\", comment=\"\"\"Thermo library: Klippenstein_Glarborg2016\"\"\"),"
    )
    assert needle in source_text
    positional_thermo = needle.replace('NASA(polynomials=[', 'NASA([', 1)
    assert positional_thermo != needle
    positional_text = source_text.replace(needle, positional_thermo)
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(positional_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='positional'):
        write_hybrid_network_input_file(
            source_path=source_path,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
        )


def test_warns_when_the_hybrid_networks_ceiling_skips_a_species(tmp_path, caplog):
    """Finding B: same visibility requirement as write_arkane_network_input_file -- when
    network_thermo_t_max's ceiling silently excludes a species whose thermo it could not read,
    the ceiling used to clamp the hybrid network may be too high, so this must not be silent.
    Replaces 'H(34)' species' NASA thermo with ThermoData (which carries no NASA Tmax to read),
    so H(34) is skipped while 'C4rad(5)' (5000 K outer NASA ceiling) still drives the clamp;
    asserts a WARNING is logged naming H(34) and that the ceiling could not account for it."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    h34_thermo = (
        "thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], "
        "Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], "
        "Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), "
        "E0=(211.8,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), "
        "label=\"\"\"H\"\"\", comment=\"\"\"Thermo library: Klippenstein_Glarborg2016\"\"\"),"
    )
    assert h34_thermo in source_text
    mixed_text = source_text.replace(
        h34_thermo,
        "thermo = ThermoData(Tdata=([300,400,500],'K'), Cpdata=([20.8,20.8,20.8],'J/(mol*K)'), "
        "H298=(0,'kJ/mol'), S298=(114.7,'J/(mol*K)')),",
    )
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(mixed_text)

    dest_path = str(tmp_path / 'input.py')
    import logging
    with caplog.at_level(logging.WARNING):
        write_hybrid_network_input_file(
            source_path=source_path,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
        )
    assert any("'H(34)'" in record.message and 'could not account for' in record.message
               for record in caplog.records)


def test_raises_on_hybrid_tlist_with_non_kelvin_unit(tmp_path):
    """Finding D: same shape validation as write_arkane_network_input_file -- a Tlist unit other
    than 'K' must not be silently compared as if it were Kelvin. A clamp is in play here (ceiling
    narrowed to 1500 K, pdep Tmax raised to 3200 K), so this must raise a ValueError naming the
    file and what wasn't understood."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    narrowed_text = source_text.replace(
        "Tmax=(5000,'K'), E0=(63.0573,'kJ/mol')", "Tmax=(1500,'K'), E0=(63.0573,'kJ/mol')")
    narrowed_text = narrowed_text.replace("    Tmax = (2100,'K'),\n", "    Tmax = (3200,'K'),\n")
    assert "],'K')" in narrowed_text
    bad_unit_text = narrowed_text.replace(
        "Tlist = ([302.491,323.355,370.585,457.988,614.983,900.017,1394.8,1985.54],'K'),\n",
        "Tlist = ([302.491,323.355,370.585,457.988,614.983,900.017,1394.8,1985.54],'degC'),\n",
    )
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(bad_unit_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='Tlist'):
        write_hybrid_network_input_file(
            source_path=source_path,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
            sensitivity=True,
        )


def test_raises_on_hybrid_tlist_with_non_2_tuple_shape(tmp_path):
    """Finding D: a Tlist RHS that literal-evaluates to a 2-tuple whose first element is not a
    sequence of numbers (a single scalar here) must not be indexed/iterated blindly. A clamp is
    in play (ceiling narrowed to 1500 K, pdep Tmax raised to 3200 K), so this must raise a
    ValueError naming the file and what wasn't understood."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    narrowed_text = source_text.replace(
        "Tmax=(5000,'K'), E0=(63.0573,'kJ/mol')", "Tmax=(1500,'K'), E0=(63.0573,'kJ/mol')")
    narrowed_text = narrowed_text.replace("    Tmax = (2100,'K'),\n", "    Tmax = (3200,'K'),\n")
    bad_shape_text = narrowed_text.replace(
        "Tlist = ([302.491,323.355,370.585,457.988,614.983,900.017,1394.8,1985.54],'K'),\n",
        "Tlist = (300.0,'K'),\n",
    )
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(bad_shape_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='Tlist'):
        write_hybrid_network_input_file(
            source_path=source_path,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
            sensitivity=True,
        )


def test_raises_on_hybrid_tlist_with_non_numeric_entries(tmp_path):
    """Finding D: a Tlist sequence containing a non-numeric entry must not be compared against
    the numeric thermo ceiling with a bare '>'. A clamp is in play (ceiling narrowed to 1500 K,
    pdep Tmax raised to 3200 K), so this must raise a ValueError naming the file and what wasn't
    understood, rather than an opaque TypeError."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    narrowed_text = source_text.replace(
        "Tmax=(5000,'K'), E0=(63.0573,'kJ/mol')", "Tmax=(1500,'K'), E0=(63.0573,'kJ/mol')")
    narrowed_text = narrowed_text.replace("    Tmax = (2100,'K'),\n", "    Tmax = (3200,'K'),\n")
    bad_entry_text = narrowed_text.replace(
        "Tlist = ([302.491,323.355,370.585,457.988,614.983,900.017,1394.8,1985.54],'K'),\n",
        "Tlist = (['not-a-number',323.355,370.585,457.988,614.983,900.017,1394.8,1985.54],'K'),\n",
    )
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(bad_entry_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='Tlist'):
        write_hybrid_network_input_file(
            source_path=source_path,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
            sensitivity=True,
        )


def test_raises_on_multiline_hybrid_tlist(tmp_path):
    """Finding D: the Tlist-drop logic's line-scan parsing assumes a single-line Tlist RHS. A
    multiline Tlist must be refused clearly with a ValueError naming the file, not parsed
    truncated/incorrectly by the line-scan or crashed on opaquely elsewhere. A clamp is in play
    (ceiling narrowed to 1500 K, pdep Tmax raised to 3200 K)."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    narrowed_text = source_text.replace(
        "Tmax=(5000,'K'), E0=(63.0573,'kJ/mol')", "Tmax=(1500,'K'), E0=(63.0573,'kJ/mol')")
    narrowed_text = narrowed_text.replace("    Tmax = (2100,'K'),\n", "    Tmax = (3200,'K'),\n")
    multiline_text = narrowed_text.replace(
        "    Tlist = ([302.491,323.355,370.585,457.988,614.983,900.017,1394.8,1985.54],'K'),\n",
        "    Tlist = ([302.491,323.355,370.585,457.988,614.983,900.017,1394.8,\n"
        "        1985.54],'K'),\n",
    )
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(multiline_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='Tlist'):
        write_hybrid_network_input_file(
            source_path=source_path,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
            sensitivity=True,
        )


def test_raises_if_transition_state_kwarg_is_present_but_unevaluable(tmp_path):
    """B3: a path reaction's transitionState=... kwarg that IS present but is not a literal value
    (e.g. transitionState=some_variable) must RAISE, not silently fall through
    _literal_or_none(...) to None. Falling through to None means `ts_label in qm_transition_states`
    is False, so the reaction's 'kinetics = ...' entry is neither dropped nor warned about, even
    when the reaction's real (unevaluable) transition state is actually one of the QM'd ones --
    breaking this module's invariant that a QM'd path reaction must never keep a 'kinetics = ...'
    entry (Arkane would otherwise silently fall back to ILT while HybridNetworkResult.qm_ts_labels
    falsely claims QM was used).

    Note: as of this writing, t3.pdep.parser.parse_pdep_network_file's own preflight parse (called
    by write_hybrid_network_input_file before this module's own AST loop ever runs) already raises
    for this same input via its own '_literal_or_raise' fail-closed helper -- so this end-to-end
    test alone cannot prove write_hybrid_network_input_file's OWN reaction-branch logic (the code
    this test is nominally about) independently raises. See
    test_raises_if_transition_state_kwarg_is_present_but_unevaluable_in_isolation below, which
    monkeypatches that preflight call out of the way to exercise this module's own fix directly."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    assert "    transitionState = 'TS1',\n" in source_text
    unevaluable_text = source_text.replace(
        "    transitionState = 'TS1',\n", "    transitionState = some_undefined_variable,\n")
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(unevaluable_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='transitionState'):
        write_hybrid_network_input_file(
            source_path=source_path,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
        )


def test_raises_if_transition_state_kwarg_is_present_but_unevaluable_in_isolation(tmp_path, monkeypatch):
    """B3, isolated from t3.pdep.parser's own (separately fixed) preflight check: monkeypatch
    parse_pdep_network_file out of the way so this test exercises ONLY
    write_hybrid_network_input_file's own reaction-branch handling of a present-but-unevaluable
    transitionState=... value, proving this module's own fix (not parser.py's) is what raises."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    assert "    transitionState = 'TS1',\n" in source_text
    unevaluable_text = source_text.replace(
        "    transitionState = 'TS1',\n", "    transitionState = some_undefined_variable,\n")
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(unevaluable_text)

    class _FakeNetwork:
        network_id = 'network4_1'
        transition_state_labels = ('TS1', 'TS2', 'TS3', 'TS4', 'TS5')
        species_labels = ()
        isomers = ()  # wells are validated against the isomers, not every species() block

    monkeypatch.setattr(hybrid_module, 'parse_pdep_network_file', lambda path: _FakeNetwork())

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='transitionState'):
        hybrid_module.write_hybrid_network_input_file(
            source_path=source_path,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
        )


def test_raises_if_generated_network_file_would_fail_to_parse(tmp_path, monkeypatch):
    """Item 8: before returning success, the generated network file's own text must be
    self-checked with ast.parse(...); if it would not parse (e.g. a splice went wrong), the write
    must raise a clear RuntimeError rather than silently handing back a broken (or worse, partially
    written) hybrid network."""
    real_edit = hybrid_module._kinetics_removal_edit

    def _corrupting_edit(**kwargs):
        start, end, replacement = real_edit(**kwargs)
        return start, end, replacement + 'this is not valid python !!! (\n'

    monkeypatch.setattr(hybrid_module, '_kinetics_removal_edit', _corrupting_edit)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(RuntimeError, match='input.py'):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
        )

    assert not os.path.exists(dest_path), 'a network file that fails its own self-check must never be written'


@pytest.mark.parametrize('field_name, value', [
    # The quoted 'False' is the dangerous shape: non-empty, hence TRUTHY, so
    # `if not energy_settings.use_atom_corrections` never fires -- while the f-string in
    # _build_energy_header renders it back out as a bare `useAtomCorrections = False` that Arkane
    # evaluates as the boolean. The guard is bypassed and the directive inverted in one step.
    ('use_atom_corrections', 'False'),
    ('use_atom_corrections', 0),
    ('use_atom_corrections', None),
    ('use_hindered_rotors', 'False'),
    ('use_hindered_rotors', 1),
    ('use_bond_corrections', 'True'),
    ('frequency_scale_factor', '0.988'),
    ('frequency_scale_factor', True),
    ('bond_correction_type', 5),
    ('atom_energies', 'H: -0.5'),
    ('model_chemistry', 3.14),
    ('tunneling', 7),
])
def test_refuses_wrongly_typed_energy_settings_fields(tmp_path, field_name, value):
    """QMEnergySettings is a plain dataclass, so it enforces nothing: the values reaching it come
    from a YAML-loaded capture manifest, which can carry any type at all. The write must therefore
    type-check every field itself, at the boundary, rather than trusting the annotation."""
    energy_settings = dataclasses.replace(DEFAULT_ENERGY_SETTINGS, **{field_name: value})
    with pytest.raises(ValueError, match=field_name):
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=str(tmp_path / 'input.py'),
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=energy_settings,
        )


class TestKwargsUnpackingIsRefusedByTheRewriterToo:
    """
    A ``**kwargs`` unpacking must never reach the writer -- and it is UPSTREAM that stops it.

    ``hybrid.py`` has two sites that read a recognized top-level call's keywords with
    ``{kw.arg: ... for kw in call.keywords if kw.arg is not None}``, discarding any ``**kwargs``
    unpacking. They were reported as carrying the same live defect as
    ``t3.pdep.parser._call_keywords``. They do not: ``write_hybrid_network_input_file`` calls
    ``parse_pdep_network_file`` unconditionally before that walk, and the parser's
    RECOGNIZED_TOP_LEVEL_CALLS is a strict superset of this module's _TOP_LEVEL_CALL_NAMES, so every
    call reaching those comprehensions has already been refused if it carried an unpacking. Adding a
    second guard in hybrid.py would be dead code, so none was added.

    That makes these tests tests of a COMPOSITION, and composition is what quietly breaks: the
    correctness of those two lines lives in another module. So the premise is asserted directly
    (``test_the_refusal_really_does_come_from_the_upstream_parse``) rather than left implicit -- if the
    parse is ever moved, made conditional, or its recognized set narrowed, that test fails and names
    the reason instead of leaving a silently reopened hole.

    What would happen if the unpacking DID get through, which is why this matters at all:

    * ``transitionState(**payload)`` yields ``label is None``, which is not in
      ``qm_transition_states``, so the QM artifact is never vendored in -- the written network keeps
      the RMG estimate while T3's records say the transition state was computed by QM.
    * ``reaction(**payload)`` yields no ``transitionState`` keyword at all, so the reaction looks like
      one with no transition state and its ``kinetics = ...`` entry is never dropped -- the ILT
      estimate survives into a file T3 believes carries QM kinetics.

    Neither raises, and the written file is valid Python either way. The refusal immediately below the
    second site already states the principle exactly ("Refuse rather than degrade ... the same
    fail-open bug already fixed for 'products='/'isomers='/'reactants='").
    """

    def test_the_refusal_really_does_come_from_the_upstream_parse(self, tmp_path, monkeypatch):
        """
        The premise the other two tests rest on, asserted instead of assumed.

        With the parser's refusal neutralised, the write must go THROUGH -- which is what proves
        hybrid.py's own comprehensions are not doing this work, and that the two tests below are
        pinning a composition rather than a local guard. If someone later adds a guard inside
        hybrid.py, this test fails and says so, which is the correct outcome: two guards for one rule
        is a decision to make deliberately, not to discover.
        """
        import t3.pdep.parser as parser_module

        original = parser_module._call_keywords
        monkeypatch.setattr(parser_module, '_call_keywords',
                            lambda call, path='', call_name=None: {kw.arg: kw.value
                                                                   for kw in call.keywords
                                                                   if kw.arg is not None})
        assert parser_module._call_keywords is not original, 'the patch must actually be in place'

        variant = ("reaction(\n"
                   "    label = 'reaction1',\n"
                   "    **payload\n"
                   ")\n")
        source_path = _write_network_with_reaction1_variant(tmp_path, variant)
        dest_path = str(tmp_path / 'input.py')
        write_hybrid_network_input_file(
            source_path=source_path,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
        )
        assert os.path.isfile(dest_path), \
            'with the upstream refusal removed the write succeeds -- hybrid.py does not refuse on its own'

    def test_a_kwargs_unpacking_on_a_reaction_is_refused(self, tmp_path):
        """The reaction site: an unpacking must not read as "this reaction has no transitionState"."""
        variant = (
            "reaction(\n"
            "    label = 'reaction1',\n"
            "    **payload\n"
            ")\n"
        )
        source_path = _write_network_with_reaction1_variant(tmp_path, variant)
        dest_path = str(tmp_path / 'input.py')
        with pytest.raises(ValueError, match=r'\*\*kwargs'):
            write_hybrid_network_input_file(
                source_path=source_path,
                dest_path=dest_path,
                method='CSE',
                qm_transition_states={'TS1': ts_conformer()},
                energy_settings=DEFAULT_ENERGY_SETTINGS,
            )
        assert not os.path.isfile(dest_path), 'nothing may be written once the source is refused'

    def test_a_kwargs_unpacking_on_a_transition_state_is_refused(self, tmp_path):
        """
        The transitionState site, which the reaction test does NOT cover.

        The two comprehensions are separate code, in separate branches, so a fix applied to one leaves
        the other open. Pinning both is what makes that impossible to half-do.
        """
        with open(NETWORK_FIXTURE, 'r') as f:
            text = f.read()
        assert "transitionState(" in text, 'fixture precondition: the source must declare a transitionState'
        text = text + "\ntransitionState(**payload)\n"
        source_path = str(tmp_path / 'network_variant.py')
        with open(source_path, 'w') as f:
            f.write(text)
        dest_path = str(tmp_path / 'input.py')
        with pytest.raises(ValueError, match=r'\*\*kwargs'):
            write_hybrid_network_input_file(
                source_path=source_path,
                dest_path=dest_path,
                method='CSE',
                qm_transition_states={'TS1': ts_conformer()},
                energy_settings=DEFAULT_ENERGY_SETTINGS,
            )
        assert not os.path.isfile(dest_path)

    def test_an_ordinary_network_without_any_unpacking_still_writes(self, tmp_path):
        """Over-refusal guard: the real fixture, which uses explicit keywords throughout."""
        dest_path = str(tmp_path / 'input.py')
        write_hybrid_network_input_file(
            source_path=NETWORK_FIXTURE,
            dest_path=dest_path,
            method='CSE',
            qm_transition_states={'TS1': ts_conformer()},
            energy_settings=DEFAULT_ENERGY_SETTINGS,
        )
        assert os.path.isfile(dest_path)


# ---------------------------------------------------------------------------------------------------
# I-023: adopt QM wells, and make the degrees-of-freedom treatment a network-wide invariant.
#
# These tests drive the REAL committed CHO2 round-0 data, not a synthetic stub: the source network
# is the real RMG-explored network the pilot run produced, and the conformer data is the real
# Arkane-extracted (E0, frequencies, imaginary) for both wells and the transition state, on the one
# self-consistent energy reference. Every prior defect on this producer/consumer handoff was
# invisible to synthetic fixtures, so the standing rule for this ticket is that at least one test
# drives the real artifact end to end.
# ---------------------------------------------------------------------------------------------------

CHO2_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_hybrid', 'cho2_round0')
CHO2_SOURCE = os.path.join(CHO2_DIR, 'source_network0_reduced.py')
CHO2_CONFORMERS = os.path.join(CHO2_DIR, 'conformers.json')
CHO2_ENERGY_SETTINGS = QMEnergySettings(
    model_chemistry="LevelOfTheory(method='wb97xd2023',basis='def2tzvp',software='gaussian')",
    atom_energies={
        'Br': -2574.174533595486, 'C': -37.84706210301937, 'Cl': -460.1467876783656,
        'F': -99.73955550924293, 'H': -0.5006557872395249, 'N': -54.584995947182875,
        'O': -75.07252406126821, 'S': -398.1105530401693,
    },
    tunneling='Eckart',
)
# The two isomers of the real CHO2 network and its QM'd transition state.
CHO2_WELL_QM = 'O=[C]O(8)'
CHO2_WELL_ESTIMATED = '[O]C=O(1)'
CHO2_TS = 'TS2'


def _load_cho2_conformers():
    import json
    with open(CHO2_CONFORMERS) as f:
        return json.load(f)


def _external_dof_mode_names(text):
    """Every statmech mode class *constructed* in the network text that carries translation or
    overall rotation. Matches the ``Name(`` constructor form so it never trips on the unrelated
    ``activeKRotor``/``activeJRotor`` pressureDependence flags."""
    return [name for name in ('IdealGasTranslation', 'LinearRotor', 'NonlinearRotor', 'KRotor',
                              'SphericalTopRotor') if f'{name}(' in text]


def _species_block(text, label):
    start = text.index(f"label = '{label}'")
    return text[start:text.index('\n)\n', start)]


@pytest.mark.parametrize('non_isomer_label', ['Ar', 'O=C=O(5)', '[H](6)'])
def test_only_an_isomer_can_be_adopted_as_a_well(tmp_path, non_isomer_label):
    """`qm_wells` must be validated against the network's ISOMERS, not against every species block.

    The CHO2 source declares five species(...) blocks but only two isomers: the other three are the
    bath gas `Ar` and the bimolecular channel species `O=C=O(5)` and `[H](6)`. Validating against
    `species_labels` accepted all five, so `qm_wells={'Ar': ...}` would have rewritten the BATH
    GAS's E0 and modes -- and nothing downstream would have caught it, because the write-time
    degrees-of-freedom invariant deliberately skips non-isomer species (channel species legitimately
    carry full conformers).
    """
    conf = _load_cho2_conformers()
    with pytest.raises(ValueError) as exc_info:
        write_hybrid_network_input_file(
            source_path=CHO2_SOURCE,
            dest_path=str(tmp_path / 'input.py'),
            method='MSC',
            qm_transition_states={CHO2_TS: conf[CHO2_TS]},
            qm_wells={non_isomer_label: conf[CHO2_WELL_QM]},
            energy_settings=CHO2_ENERGY_SETTINGS,
        )
    message = str(exc_info.value)
    assert non_isomer_label in message
    assert 'not isomers' in message


def test_i023_mixed_network_adopts_one_qm_well_and_stays_dof_consistent(tmp_path):
    """The mixed case, on the real CHO2 network: one well is adopted from QM, the other keeps its
    RMG estimate, and the transition state is QM'd. The produced network must be
    degrees-of-freedom-consistent -- every well and the TS vibration-only, no translation/rotation
    anywhere -- which is the invariant this ticket exists to guarantee even for a partially-explored
    network."""
    conf = _load_cho2_conformers()
    dest_path = str(tmp_path / 'input.py')
    result = write_hybrid_network_input_file(
        source_path=CHO2_SOURCE,
        dest_path=dest_path,
        method='MSC',
        qm_transition_states={CHO2_TS: conf[CHO2_TS]},
        qm_wells={CHO2_WELL_QM: conf[CHO2_WELL_QM]},  # only ONE well has QM data
        energy_settings=CHO2_ENERGY_SETTINGS,
    )
    assert result.qm_ts_labels == (CHO2_TS,)
    assert result.qm_well_labels == (CHO2_WELL_QM,)

    with open(dest_path) as f:
        text = f.read()

    # The invariant: no isomer or TS carries a translational/rotational mode.
    assert _external_dof_mode_names(text) == [], \
        f"DOF-inconsistent network leaked external modes: {_external_dof_mode_names(text)}"
    # The network is self-contained: it never falls back to a by-reference qm/ statmech file.
    assert 'qm/' not in text

    # The adopted well carries the QM E0 and QM frequencies (vibration-only, HarmonicOscillator).
    qm_well_block = _species_block(text, CHO2_WELL_QM)
    assert f"({conf[CHO2_WELL_QM]['E0_kJ_mol']!r},'kJ/mol')" in qm_well_block
    assert 'HarmonicOscillator' in qm_well_block

    # The un-adopted well keeps its RMG-estimated E0 -- it was NOT rewritten -- yet is still
    # vibration-only, so the mixed-provenance network is DOF-consistent.
    estimated_block = _species_block(text, CHO2_WELL_ESTIMATED)
    assert "E0 = (-169.935,'kJ/mol')" in estimated_block
    assert 'HarmonicOscillator' in estimated_block

    # The TS is inline (not by-reference) with its imaginary frequency and QM E0.
    assert f"({conf[CHO2_TS]['E0_kJ_mol']!r},'kJ/mol')" in text
    assert f"frequency = ({conf[CHO2_TS]['imaginary_frequency_cm_1']!r},'cm^-1')" in text


def test_i023_full_qm_network_adopts_both_wells_and_is_dof_consistent(tmp_path):
    """Both wells adopted from QM plus the QM'd TS: the whole isomer set is on one QM energy
    reference and every conformer is vibration-only. This is the network whose round-1 master
    equation solves (verified out of band by driving Arkane's explorer on it)."""
    conf = _load_cho2_conformers()
    dest_path = str(tmp_path / 'input.py')
    result = write_hybrid_network_input_file(
        source_path=CHO2_SOURCE,
        dest_path=dest_path,
        method='MSC',
        qm_transition_states={CHO2_TS: conf[CHO2_TS]},
        qm_wells={CHO2_WELL_QM: conf[CHO2_WELL_QM], CHO2_WELL_ESTIMATED: conf[CHO2_WELL_ESTIMATED]},
        energy_settings=CHO2_ENERGY_SETTINGS,
    )
    assert result.qm_well_labels == tuple(sorted((CHO2_WELL_QM, CHO2_WELL_ESTIMATED)))

    with open(dest_path) as f:
        text = f.read()
    assert _external_dof_mode_names(text) == []
    # Both wells now carry their QM E0 (on the same reference as the TS).
    assert f"({conf[CHO2_WELL_QM]['E0_kJ_mol']!r},'kJ/mol')" in text
    assert f"({conf[CHO2_WELL_ESTIMATED]['E0_kJ_mol']!r},'kJ/mol')" in text
    # The reaction whose TS was QM'd drops its ILT kinetics and takes tunneling.
    assert "tunneling = 'Eckart'" in text


def test_i023_refuses_a_conformer_with_a_non_hindered_rotor_mode(tmp_path):
    """The DOF invariant is enforced up front: a conformer-data dict whose internal-rotor entry is
    not a HinderedRotor (e.g. the extractor regressing and leaking an external rotor) must be
    refused before anything reaches disk, rather than emitting a DOF-inconsistent network."""
    conf = _load_cho2_conformers()
    bad = dict(conf[CHO2_WELL_QM])
    bad['hindered_rotors'] = [{'type': 'NonlinearRotor'}]
    with pytest.raises(ValueError, match='HinderedRotor'):
        write_hybrid_network_input_file(
            source_path=CHO2_SOURCE,
            dest_path=str(tmp_path / 'input.py'),
            method='MSC',
            qm_transition_states={CHO2_TS: conf[CHO2_TS]},
            qm_wells={CHO2_WELL_QM: bad},
            energy_settings=CHO2_ENERGY_SETTINGS,
        )
