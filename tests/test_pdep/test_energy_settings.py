"""
Tests for t3.pdep.energy_settings.

Fixture projects (all under tests/data/pdep_energy_settings/) are copied verbatim (input.py) or
hand-written but data-faithful (output.yml) from real ARC-generated output:

* xl1001_project -- a real multi-TS rate project. Its 'kinetics' subdir sets
  useBondCorrections = False while its 'thermo' subdir sets useBondCorrections = True,
  bondCorrectionType = 'p'; both subdirs share one project-level output.yml whose
  bond_additivity_corrections dict is non-null. This is the concrete regression fixture for the
  "output.yml is project-global, not per-subdir" finding (see energy_settings.py's module
  docstring): reading the 'kinetics' subdir must NOT raise and must NOT adopt that dict.
* composite_level_project -- a real benzyl thermo project using a bare, multi-line
  CompositeLevelOfTheory(freq=..., energy=...) modelChemistry expression, with sp_level != freq_level
  in output.yml (compositeness signaled that way, not via composite_method) and no
  frequencyScaleFactor in either file.
* missing_model_chemistry_project -- a real e2e project whose kinetics/input.py has no
  modelChemistry directive at all, and whose corrections are off (useAtomCorrections =
  useBondCorrections = False) even though output.yml still carries a non-null bac_type -- the
  concrete regression fixture for "null means unknown, not off" (and its mirror: a non-null
  project-level field is not evidence of a contradiction either).
"""

import copy
import os

import pytest
import yaml

from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.energy_settings import read_arc_energy_settings
from t3.pdep.hybrid import QMEnergySettings, _build_energy_header, _validate_model_chemistry_expression

DATA_BASE_PATH = os.path.join(TEST_DATA_BASE_PATH, 'pdep_energy_settings')


def test_xl1001_kinetics_subdir_does_not_adopt_project_global_bond_corrections():
    """The 'kinetics' subdir sets useBondCorrections = False; output.yml's non-null
    bond_additivity_corrections dict belongs to the sibling 'thermo' subdir's pass and must not be
    adopted, nor treated as a contradiction."""
    result = read_arc_energy_settings(os.path.join(DATA_BASE_PATH, 'xl1001_project'), statmech_subdir='kinetics')
    assert result['model_chemistry'] == "LevelOfTheory(method='wb97xd2023',basis='def2tzvp',software='gaussian')"
    assert result['frequency_scale_factor'] == pytest.approx(0.988)
    assert result['use_hindered_rotors'] is True
    assert result['use_atom_corrections'] is True
    assert result['use_bond_corrections'] is False
    assert result['bond_correction_type'] is None
    assert result['bond_additivity_corrections'] is None
    # Atom corrections ARE on for this subdir, so they should be populated from output.yml.
    assert result['atom_energies'] == {
        'Br': -2574.174533595486,
        'C': -37.84706210301937,
        'Cl': -460.1467876783656,
        'F': -99.73955550924293,
        'H': -0.5006557872395249,
        'N': -54.584995947182875,
        'O': -75.07252406126821,
        'S': -398.1105530401693,
    }


def test_xl1001_thermo_subdir_adopts_matching_project_global_bond_corrections():
    """The 'thermo' subdir sets useBondCorrections = True, bondCorrectionType = 'p', matching
    output.yml's bac_type = 'p'; the corrections dict should be adopted and no mismatch raised."""
    result = read_arc_energy_settings(os.path.join(DATA_BASE_PATH, 'xl1001_project'), statmech_subdir='thermo')
    assert result['use_bond_corrections'] is True
    assert result['bond_correction_type'] == 'p'
    assert result['bond_additivity_corrections'] == {
        'C-C': 2.9585184472615715,
        'C-H': 0.7891610772258629,
        'C-N': 1.7702965975878437,
        'C-O': 2.5823551057599605,
        'C-S': 2.403228751492449,
        'H-H': -3.6243516474330164,
    }


def test_composite_level_of_theory_expression_preserved_verbatim():
    """A bare, multi-line CompositeLevelOfTheory(...) expression must be recovered exactly as
    written (not reformatted/reconstructed), and the absence of frequencyScaleFactor in both files
    must resolve to None without raising."""
    result = read_arc_energy_settings(os.path.join(DATA_BASE_PATH, 'composite_level_project'),
                                       statmech_subdir='thermo')
    assert result['model_chemistry'] == (
        "CompositeLevelOfTheory(\n"
        "    freq=LevelOfTheory(method='wb97xd',basis='def2tzvp',software='gaussian'),\n"
        "    energy=LevelOfTheory(method='dlpnoccsd(t)f122023',basis='ccpvtzf12',software='orca')\n"
        ")"
    )
    assert result['frequency_scale_factor'] is None
    assert result['use_bond_corrections'] is True
    assert result['bond_correction_type'] == 'p'
    assert result['bond_additivity_corrections'] == {
        'C-C': -459.69651810328156,
        'C-H': -99.65901652526594,
        'C#N': -0.3283125828210814,
        'C-N': 0.17286574008491584,
    }


def test_missing_model_chemistry_raises():
    """input.py has no modelChemistry assignment at all: required, must raise."""
    with pytest.raises(ValueError, match='modelChemistry'):
        read_arc_energy_settings(os.path.join(DATA_BASE_PATH, 'missing_model_chemistry_project'),
                                  statmech_subdir='kinetics')


def test_corrections_off_but_output_yml_bac_type_present_is_not_an_error(tmp_path):
    """Regression fixture for 'null means unknown, not off' and its mirror: this subdir has both
    corrections flags False, yet output.yml still records bac_type = 'p' (a project-level leftover
    unrelated to this pass, since bondCorrectionType is entirely absent from this subdir's
    input.py so there is nothing of this pass's to compare it against). This must not raise, and
    frozen corrections fields must be None (not populated from the unrelated project-level data)."""
    # This fixture's input.py has no modelChemistry, so add one via a temp copy to isolate the
    # "corrections off" behavior from the "modelChemistry required" behavior tested above.
    src_project = os.path.join(DATA_BASE_PATH, 'missing_model_chemistry_project')
    project_dir = tmp_path / 'project'
    for root, _dirs, files in os.walk(src_project):
        rel = os.path.relpath(root, src_project)
        dest_dir = project_dir / rel if rel != '.' else project_dir
        dest_dir.mkdir(parents=True, exist_ok=True)
        for name in files:
            with open(os.path.join(root, name), 'r') as f:
                content = f.read()
            if name == 'input.py':
                content = content.replace(
                    "frequencyScaleFactor = 1.0",
                    "modelChemistry = 'mock/mockbasis'\nfrequencyScaleFactor = 1.0")
            with open(dest_dir / name, 'w') as f:
                f.write(content)
    result = read_arc_energy_settings(str(project_dir), statmech_subdir='kinetics')
    assert result['use_atom_corrections'] is False
    assert result['use_bond_corrections'] is False
    assert result['atom_energies'] is None
    assert result['bond_correction_type'] is None
    assert result['bond_additivity_corrections'] is None


def test_bond_correction_type_mismatch_raises(tmp_path):
    """A genuine same-pass disagreement between bondCorrectionType and output.yml's bac_type must
    raise (unlike the project-global non-null-vs-off case, this is an actual same-file
    contradiction check on values that both concern this pass)."""
    src_project = os.path.join(DATA_BASE_PATH, 'xl1001_project')
    project_dir = tmp_path / 'project'
    for root, _dirs, files in os.walk(src_project):
        rel = os.path.relpath(root, src_project)
        dest_dir = project_dir / rel if rel != '.' else project_dir
        dest_dir.mkdir(parents=True, exist_ok=True)
        for name in files:
            with open(os.path.join(root, name), 'r') as f:
                content = f.read()
            if name == 'input.py' and rel.endswith('thermo'):
                content = content.replace("bondCorrectionType = 'p'", "bondCorrectionType = 'm'")
            with open(dest_dir / name, 'w') as f:
                f.write(content)
    with pytest.raises(ValueError, match='[Bb]ond correction type mismatch'):
        read_arc_energy_settings(str(project_dir), statmech_subdir='thermo')


def test_frequency_scale_factor_mismatch_raises(tmp_path):
    """A genuine numeric disagreement between frequencyScaleFactor and output.yml's
    freq_scale_factor must raise."""
    src_project = os.path.join(DATA_BASE_PATH, 'xl1001_project')
    project_dir = tmp_path / 'project'
    for root, _dirs, files in os.walk(src_project):
        rel = os.path.relpath(root, src_project)
        dest_dir = project_dir / rel if rel != '.' else project_dir
        dest_dir.mkdir(parents=True, exist_ok=True)
        for name in files:
            with open(os.path.join(root, name), 'r') as f:
                content = f.read()
            if name == 'input.py' and rel.endswith('kinetics'):
                content = content.replace('frequencyScaleFactor = 0.988', 'frequencyScaleFactor = 0.95')
            with open(dest_dir / name, 'w') as f:
                f.write(content)
    with pytest.raises(ValueError, match='[Ff]requency scale factor mismatch'):
        read_arc_energy_settings(str(project_dir), statmech_subdir='kinetics')


def test_use_bond_corrections_true_without_bond_correction_type_raises(tmp_path):
    """Same-file contradiction: useBondCorrections = True but bondCorrectionType is never assigned
    in that same input.py. This must raise even before output.yml is consulted."""
    src_project = os.path.join(DATA_BASE_PATH, 'xl1001_project')
    project_dir = tmp_path / 'project'
    for root, _dirs, files in os.walk(src_project):
        rel = os.path.relpath(root, src_project)
        dest_dir = project_dir / rel if rel != '.' else project_dir
        dest_dir.mkdir(parents=True, exist_ok=True)
        for name in files:
            with open(os.path.join(root, name), 'r') as f:
                content = f.read()
            if name == 'input.py' and rel.endswith('thermo'):
                content = content.replace("bondCorrectionType = 'p'\n", "")
            with open(dest_dir / name, 'w') as f:
                f.write(content)
    with pytest.raises(ValueError, match='does not assign'):
        read_arc_energy_settings(str(project_dir), statmech_subdir='thermo')


def test_use_bond_corrections_false_with_bond_correction_type_raises(tmp_path):
    """Same-file contradiction: useBondCorrections = False but bondCorrectionType is nonetheless
    assigned in that same input.py."""
    src_project = os.path.join(DATA_BASE_PATH, 'xl1001_project')
    project_dir = tmp_path / 'project'
    for root, _dirs, files in os.walk(src_project):
        rel = os.path.relpath(root, src_project)
        dest_dir = project_dir / rel if rel != '.' else project_dir
        dest_dir.mkdir(parents=True, exist_ok=True)
        for name in files:
            with open(os.path.join(root, name), 'r') as f:
                content = f.read()
            if name == 'input.py' and rel.endswith('kinetics'):
                content = content.replace(
                    'useBondCorrections = False', "useBondCorrections = False\nbondCorrectionType = 'p'")
            with open(dest_dir / name, 'w') as f:
                f.write(content)
    with pytest.raises(ValueError, match='contradictory'):
        read_arc_energy_settings(str(project_dir), statmech_subdir='kinetics')


def test_missing_required_literal_field_raises(tmp_path):
    """ARC always writes useHinderedRotors/useAtomCorrections/useBondCorrections; a subdir missing
    one of them entirely is ambiguous and must raise rather than assume a default."""
    src_project = os.path.join(DATA_BASE_PATH, 'xl1001_project')
    project_dir = tmp_path / 'project'
    for root, _dirs, files in os.walk(src_project):
        rel = os.path.relpath(root, src_project)
        dest_dir = project_dir / rel if rel != '.' else project_dir
        dest_dir.mkdir(parents=True, exist_ok=True)
        for name in files:
            with open(os.path.join(root, name), 'r') as f:
                content = f.read()
            if name == 'input.py' and rel.endswith('kinetics'):
                content = content.replace('useHinderedRotors = True\n', '')
            with open(dest_dir / name, 'w') as f:
                f.write(content)
    with pytest.raises(ValueError, match='does not assign .useHinderedRotors.'):
        read_arc_energy_settings(str(project_dir), statmech_subdir='kinetics')


def test_atom_corrections_off_does_not_adopt_project_global_atom_energies(tmp_path):
    """Regression fixture for the 'null means unknown, not off' guard's OTHER direction: unlike
    missing_model_chemistry_project (where output.yml's atom_energy_corrections happens to also be
    null), this uses xl1001_project's output.yml, which has a genuinely non-null
    atom_energy_corrections dict. With useAtomCorrections flipped to False for this subdir, that
    project-global dict must NOT be adopted -- atom_energies must resolve to None. Without this
    fixture, a mutant that deletes the 'if not use_atom_corrections: return None' early return
    would go undetected, since every other fixture with useAtomCorrections = False also happens to
    have a null output.yml dict."""
    src_project = os.path.join(DATA_BASE_PATH, 'xl1001_project')
    project_dir = tmp_path / 'project'
    for root, _dirs, files in os.walk(src_project):
        rel = os.path.relpath(root, src_project)
        dest_dir = project_dir / rel if rel != '.' else project_dir
        dest_dir.mkdir(parents=True, exist_ok=True)
        for name in files:
            with open(os.path.join(root, name), 'r') as f:
                content = f.read()
            if name == 'input.py' and rel.endswith('kinetics'):
                content = content.replace('useAtomCorrections = True', 'useAtomCorrections = False')
            with open(dest_dir / name, 'w') as f:
                f.write(content)
    result = read_arc_energy_settings(str(project_dir), statmech_subdir='kinetics')
    assert result['use_atom_corrections'] is False
    assert result['atom_energies'] is None


def test_result_is_yaml_round_trippable():
    """The returned dict must survive a yaml.safe_dump / yaml.safe_load round trip unchanged --
    this is the property capture.py's manifest writer/reader will rely on."""
    result = read_arc_energy_settings(os.path.join(DATA_BASE_PATH, 'xl1001_project'), statmech_subdir='thermo')
    dumped = yaml.safe_dump(result)
    reloaded = yaml.safe_load(dumped)
    assert reloaded == result


def test_missing_input_py_raises():
    with pytest.raises(ValueError, match='statmech input file not found'):
        read_arc_energy_settings(os.path.join(DATA_BASE_PATH, 'xl1001_project'), statmech_subdir='nonexistent')


def test_duplicate_assignment_raises(tmp_path):
    src_project = os.path.join(DATA_BASE_PATH, 'xl1001_project')
    project_dir = tmp_path / 'project'
    for root, _dirs, files in os.walk(src_project):
        rel = os.path.relpath(root, src_project)
        dest_dir = project_dir / rel if rel != '.' else project_dir
        dest_dir.mkdir(parents=True, exist_ok=True)
        for name in files:
            with open(os.path.join(root, name), 'r') as f:
                content = f.read()
            if name == 'input.py' and rel.endswith('kinetics'):
                content = content.replace(
                    "frequencyScaleFactor = 0.988",
                    "frequencyScaleFactor = 0.988\nfrequencyScaleFactor = 0.99")
            with open(dest_dir / name, 'w') as f:
                f.write(content)
    with pytest.raises(ValueError, match='more than once'):
        read_arc_energy_settings(str(project_dir), statmech_subdir='kinetics')


def _write_minimal_arc_project(project_dir, model_chemistry_line: str, directives: dict | None = None) -> str:
    """
    Write a minimal but structurally valid ARC project directory whose kinetics/input.py carries
    ``model_chemistry_line`` verbatim.

    Args:
        project_dir: The directory to build the project in (a pathlib path or str).
        model_chemistry_line (str): The full ``modelChemistry = ...`` line to write.
        directives (dict, optional): Overrides for the non-modelChemistry directive lines, mapping
                                     a directive name to the verbatim source text of its value.
                                     Used to write directives ARC itself would never emit (e.g. a
                                     quoted ``'False'``), which is exactly what the type guards
                                     exist to refuse.

    Returns:
        str: The ARC project directory path.
    """
    directive_values = {'frequencyScaleFactor': '0.988',
                        'useHinderedRotors': 'True',
                        'useAtomCorrections': 'True',
                        'useBondCorrections': 'False'}
    directive_values.update(directives or dict())
    statmech_dir = os.path.join(str(project_dir), 'calcs', 'statmech', 'kinetics')
    os.makedirs(statmech_dir, exist_ok=True)
    os.makedirs(os.path.join(str(project_dir), 'output'), exist_ok=True)
    with open(os.path.join(statmech_dir, 'input.py'), 'w') as f:
        f.write("#!/usr/bin/env python\n"
                "title = 'Arkane kinetics calculation'\n"
                f'{model_chemistry_line}\n'
                + ''.join(f'{name} = {value}\n' for name, value in directive_values.items()))
    with open(os.path.join(str(project_dir), 'output', 'output.yml'), 'w') as f:
        f.write('sp_level:\n  method: wb97xd\n  basis: def2tzvp\n'
                'freq_scale_factor: 0.988\nfreq_scale_factor_source: null\n'
                'bac_type: null\natom_energy_corrections: null\n'
                'bond_additivity_corrections: null\n')
    return str(project_dir)


@pytest.mark.parametrize('model_chemistry_line, expected_frozen, expected_rendered', [
    # The bare object-expression form must be frozen VERBATIM and rendered BARE, so that Arkane
    # evaluates it as executable DSL rather than reading it as an inert string.
    ("modelChemistry = LevelOfTheory(method='wb97xd2023',basis='def2tzvp',software='gaussian')",
     "LevelOfTheory(method='wb97xd2023',basis='def2tzvp',software='gaussian')",
     'modelChemistry = LevelOfTheory(method=\'wb97xd2023\',basis=\'def2tzvp\',software=\'gaussian\')'),
    # The plain string-label form must be frozen as its string VALUE (no surrounding quotes) and
    # re-quoted at render time.
    ("modelChemistry = 'CBS-QB3'", 'CBS-QB3', 'modelChemistry = "CBS-QB3"'),
])
def test_frozen_model_chemistry_round_trips_into_the_hybrid_writer(tmp_path, model_chemistry_line,
                                                                   expected_frozen, expected_rendered):
    """
    A frozen ``model_chemistry`` must be directly consumable by ``QMEnergySettings``.

    This closes a seam defect in which two individually-correct behaviours composed into a hole:
    ``read_arc_energy_settings`` froze the ``modelChemistry`` SOURCE SEGMENT (which, for the plain
    string-label form, carries the surrounding quote characters), while
    ``QMEnergySettings.model_chemistry`` holds an UNQUOTED label and its validator rejects embedded
    quotes. The result was that a legitimate Arkane input round-tripped into a frozen capture value
    that the hybrid writer then refused outright -- with the ARC project directory possibly already
    gone, making it unrecoverable. Asserting the RENDERED line (not merely that validation passes)
    is what makes this test also cover the bare-vs-quoted rendering branch.
    """
    project_dir = _write_minimal_arc_project(tmp_path / 'arc_project', model_chemistry_line)
    frozen = read_arc_energy_settings(arc_project_directory=project_dir, statmech_subdir='kinetics')
    assert frozen['model_chemistry'] == expected_frozen
    # Must not raise -- this is the exact call the hybrid writer makes on the frozen value.
    _validate_model_chemistry_expression('model_chemistry', frozen['model_chemistry'])
    header = _build_energy_header(QMEnergySettings(model_chemistry=frozen['model_chemistry'],
                                                   frequency_scale_factor=frozen['frequency_scale_factor']))
    assert header.splitlines()[0] == expected_rendered


def test_frozen_model_chemistry_rejects_a_non_string_literal(tmp_path):
    """A numeric modelChemistry is neither a label nor an object expression; refuse it rather than
    freezing a value nothing downstream could use."""
    project_dir = _write_minimal_arc_project(tmp_path / 'arc_project', 'modelChemistry = 3.14')
    with pytest.raises(ValueError, match='non-string literal'):
        read_arc_energy_settings(arc_project_directory=project_dir, statmech_subdir='kinetics')


@pytest.mark.parametrize('directive_name, directive_source', [
    # A quoted 'False' is the dangerous one: it is a NON-EMPTY string, so every truthiness test
    # downstream ('if not use_atom_corrections', 'if use_bond_corrections') reads it as True,
    # while an f-string renders it back into the generated Arkane input as a bare False that
    # Arkane then evaluates as the boolean. The guard is bypassed and the directive is inverted.
    ('useAtomCorrections', "'False'"),
    ('useHinderedRotors', "'False'"),
    ('useBondCorrections', "'True'"),
    # A truthy int is the same defect wearing a different type.
    ('useAtomCorrections', '1'),
    ('useBondCorrections', '0'),
])
def test_frozen_boolean_directives_reject_non_boolean_literals(tmp_path, directive_name, directive_source):
    """A boolean Arkane directive must freeze as an actual bool. ARC never writes anything else,
    so a non-bool here means the file is not what it claims to be -- and silently accepting it
    hands every downstream truthiness guard a value whose truthiness is the OPPOSITE of what the
    directive says."""
    project_dir = _write_minimal_arc_project(tmp_path / 'arc_project',
                                             "modelChemistry = 'CBS-QB3'",
                                             directives={directive_name: directive_source})
    with pytest.raises(ValueError, match=f"'{directive_name}'"):
        read_arc_energy_settings(arc_project_directory=project_dir, statmech_subdir='kinetics')


@pytest.mark.parametrize('directive_name, directive_source', [
    ('frequencyScaleFactor', "'0.988'"),
    ('frequencyScaleFactor', 'True'),
    ('bondCorrectionType', '5'),
    ('atomEnergies', "'H: -0.5'"),
])
def test_frozen_non_boolean_directives_reject_wrongly_typed_literals(tmp_path, directive_name, directive_source):
    """Every frozen field is type-checked, not just the booleans: a string frequency scale factor
    would render into the generated input as a bare, unquoted token, and ``True`` is an int in
    Python's numeric tower so it would pass a naive isinstance(value, (int, float)) check while
    scaling every frequency by 1.0."""
    directives = {directive_name: directive_source}
    if directive_name == 'bondCorrectionType':
        directives['useBondCorrections'] = 'True'
    project_dir = _write_minimal_arc_project(tmp_path / 'arc_project',
                                             "modelChemistry = 'CBS-QB3'",
                                             directives=directives)
    with pytest.raises(ValueError, match=f"'{directive_name}'"):
        read_arc_energy_settings(arc_project_directory=project_dir, statmech_subdir='kinetics')
