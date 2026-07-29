"""
t3 pdep energy_settings module

Freezes a complete, authoritative energy-reference settings block from an ARC project directory
at TS-artifact **capture time**, so that a later hybrid network write (see
``t3.pdep.hybrid.write_hybrid_network_input_file``, which requires a
``t3.pdep.hybrid.QMEnergySettings`` with a non-blank ``model_chemistry``) needs nothing further
from the ARC project directory.

Why this module exists: ARC recreates ``<project>/calcs/statmech/<kinetics|thermo>/`` from scratch
on every rate pass (``delete_existing_subdir=True`` in ``arc/scheduler.py``), so that directory's
``input.py`` -- the only place the ``modelChemistry`` directive Arkane actually used for that pass
lives -- is ephemeral. ARC's per-TS statmech artifact (``arc/statmech/arkane.py``'s
``species_input_template``) does not carry the level of theory at all. If the TS-artifact capture
step does not freeze this information immediately, it is gone by the time a hybrid network needs
to be written.

This module reads exactly two files per call:

* ``<arc_project_directory>/calcs/statmech/<statmech_subdir>/input.py`` -- an executable Arkane DSL
  file. **NEVER exec'd, never eval'd, never imported.** It is read as text and inspected with
  ``ast.parse`` only, mirroring the house pattern in ``t3/pdep/parser.py`` and ``t3/pdep/hybrid.py``.
  ``ast.literal_eval`` is applied to individual AST nodes (never to the whole file, never via
  ``exec``) to pull out plain-literal values, and ``ast.get_source_segment`` is used to recover the
  *exact* original text of the ``modelChemistry`` expression -- which may itself be a bare,
  executable ``LevelOfTheory(...)``/``CompositeLevelOfTheory(...)`` object-expression, not just a
  string literal -- so it is preserved byte-for-byte rather than reconstructed.
* ``<arc_project_directory>/output/output.yml`` -- ARC's flattened, non-executable output summary.
  Every real ``output.yml`` inspected while building this module (three real ARC-generated
  projects) contained zero ``!!python``-tagged nodes, so it is read with a plain ``yaml.safe_load``
  (unlike ``t3/pdep/yaml_safe.py``, which exists specifically for Arkane SA YAML files that DO carry
  a ``!!python/tuple`` tag).

Design note -- ``output.yml`` is project-global, not per-statmech-pass: empirically (see the real
fixtures under ``tests/data/pdep_energy_settings/xl1001_project/``), a single project's
``output.yml`` can report a non-null ``bond_additivity_corrections`` dict and a ``bac_type`` even
for a statmech subdir whose *own* ``input.py`` sets ``useBondCorrections = False`` -- because that
project's *other* statmech subdir (e.g. ``thermo``) enabled bond corrections and both subdirs share
the one project-level ``output.yml``. Consequently this module does **not** treat "flag is False in
this subdir's ``input.py`` but the dict is non-null in ``output.yml``" as a contradiction for either
atom or bond corrections: that pattern is expected, not an error. The only cross-checks performed
are ones where a mismatch cannot be explained by a sibling subdir's differing settings: numeric
``frequencyScaleFactor`` vs. ``freq_scale_factor`` equality (a property of the frequency level of
theory, consistent across every real fixture that sets it), and ``bondCorrectionType`` vs.
``bac_type`` equality -- but *only* when this subdir's own ``input.py`` explicitly enables bond
corrections with a specific type, since only then is the comparison actually about this pass.

Design note -- "null means unknown, not off": when ``output.yml``'s corrections dict is
``null``/absent for a subdir whose ``input.py`` flag is ``True``, that means the numeric
corrections were not recorded for this particular run, NOT that corrections are disabled. This
module never conflates a missing/null dict with a ``False`` flag.

Design note -- model chemistry NAME is deliberately never cross-checked between ``input.py`` and
``output.yml``: ``output.yml``'s ``sp_level``/``freq_level``/``arkane_level_of_theory`` blocks use a
different schema (parsed sub-fields) than ``input.py``'s verbatim ``modelChemistry`` expression
text, and reconstructing one from the other to compare would risk exactly the kind of subtle
reformatting this module exists to avoid. The frozen ``model_chemistry`` value is always sourced
from ``input.py`` alone -- the one file that is actually Arkane's *input* for this pass.
"""

import ast
import math
import os

import yaml

# The two Arkane DSL call forms a modelChemistry directive may legitimately use in place of a
# plain string label; kept in sync with t3.pdep.hybrid._MODEL_CHEMISTRY_CALL_NAMES.
_MODEL_CHEMISTRY_CALL_NAMES = ('LevelOfTheory', 'CompositeLevelOfTheory')

# Fields read directly (as ast.literal_eval-able literals) from calcs/statmech/<subdir>/input.py,
# mapped to the Python type each MUST have. The type is enforced, not merely documented: a
# directive that parses to the wrong type is the single most dangerous shape a frozen setting can
# take, because a non-empty string like 'False' is TRUTHY. Every downstream truthiness guard
# ('if not use_atom_corrections', 'if use_bond_corrections') would then read it as the OPPOSITE of
# what the directive says, while an f-string renders it back into the generated Arkane input as a
# bare False that Arkane evaluates as the boolean -- a guard that is bypassed and a directive that
# is inverted, at once. bool is checked identically-typed rather than via isinstance because
# isinstance(True, int) is True in Python's numeric tower, so a bool would otherwise sail through
# the frequencyScaleFactor check and silently scale every frequency by 1.0.
_LITERAL_FIELD_TYPES = {
    'frequencyScaleFactor': (int, float),
    'useHinderedRotors': (bool,),
    'useAtomCorrections': (bool,),
    'useBondCorrections': (bool,),
    'bondCorrectionType': (str,),
    'atomEnergies': (dict,),
}
_LITERAL_FIELD_NAMES = tuple(_LITERAL_FIELD_TYPES)


def read_arc_energy_settings(arc_project_directory: str, statmech_subdir: str = 'kinetics') -> dict:
    """
    Read and freeze a complete, authoritative energy-reference settings block from an ARC project
    directory, at TS-artifact capture time.

    Args:
        arc_project_directory (str): The ARC project directory (the directory that directly
                                     contains ``calcs/`` and ``output/``).
        statmech_subdir (str, optional): Which ephemeral statmech pass to read
                                         ``calcs/statmech/<statmech_subdir>/input.py`` from.
                                         Default: ``'kinetics'``.

    Raises:
        ValueError: If ``calcs/statmech/<statmech_subdir>/input.py`` or ``output/output.yml`` does
                   not exist; if ``input.py`` cannot be parsed with ``ast.parse``; if it assigns any
                   of ``modelChemistry``/``frequencyScaleFactor``/``useHinderedRotors``/
                   ``useAtomCorrections``/``useBondCorrections``/``bondCorrectionType``/
                   ``atomEnergies`` more than once at module level (ambiguous); if it is missing
                   ``modelChemistry``, ``useHinderedRotors``, ``useAtomCorrections`` or
                   ``useBondCorrections`` (all required, since ARC always writes them); if a literal
                   field's assigned value is not a plain Python literal; if ``useBondCorrections`` is
                   ``True`` but ``bondCorrectionType`` is unset, or ``False`` but
                   ``bondCorrectionType`` is set (contradictory within the same file); if
                   ``output.yml`` cannot be parsed with ``yaml.safe_load``, or is not a mapping; if
                   ``output.yml``'s ``freq_scale_factor``, ``atom_energy_corrections`` or
                   ``bond_additivity_corrections`` fields are present but not of the expected type
                   (float/dict); or if ``input.py``'s ``frequencyScaleFactor``/``bondCorrectionType``
                   numerically/textually disagrees with ``output.yml``'s ``freq_scale_factor``/
                   ``bac_type``.

    Returns:
        dict: A YAML-round-trippable, frozen energy-settings block with keys ``model_chemistry``
             (str, the verbatim ``modelChemistry`` expression text), ``frequency_scale_factor``
             (float | None), ``use_hindered_rotors`` (bool), ``use_atom_corrections`` (bool),
             ``atom_energies`` (dict | None), ``use_bond_corrections`` (bool),
             ``bond_correction_type`` (str | None), ``bond_additivity_corrections`` (dict | None),
             and ``source_paths`` (dict with ``input_py``/``output_yml`` keys, recorded for
             traceability only -- never re-read from later).
    """
    input_py_path = os.path.join(arc_project_directory, 'calcs', 'statmech', statmech_subdir, 'input.py')
    output_yml_path = os.path.join(arc_project_directory, 'output', 'output.yml')

    if not os.path.isfile(input_py_path):
        raise ValueError(f"ARC statmech input file not found at '{input_py_path}'.")
    if not os.path.isfile(output_yml_path):
        raise ValueError(f"ARC output file not found at '{output_yml_path}'.")

    with open(input_py_path, 'r') as f:
        input_py_text = f.read()

    try:
        tree = ast.parse(input_py_text, filename=input_py_path)
    except SyntaxError as e:
        raise ValueError(f"Could not parse '{input_py_path}' as Python: {e}.") from e

    assigns_by_name = _collect_module_level_assigns(tree)

    for name, nodes in assigns_by_name.items():
        if len(nodes) > 1:
            raise ValueError(f"'{input_py_path}' assigns '{name}' more than once at module level; this is "
                             f"ambiguous and refused rather than silently taking the last value.")

    model_chemistry_nodes = assigns_by_name.get('modelChemistry')
    if not model_chemistry_nodes:
        raise ValueError(f"'{input_py_path}' does not assign 'modelChemistry'; without it, Arkane cannot apply "
                         f"atom energy corrections, so a QM'd transition state's E0 would not be on the same "
                         f"energy reference scale as the RMG wells around it.")
    # A plain string label (modelChemistry = 'CBS-QB3') must be frozen as its string VALUE, not as
    # its source segment: the segment would carry the surrounding quote characters, and
    # t3.pdep.hybrid.QMEnergySettings.model_chemistry holds an UNQUOTED label in that case (its
    # validator rejects embedded quotes, and _build_energy_header re-quotes it when rendering).
    # Freezing the segment here would round-trip a legitimate Arkane input into a frozen value that
    # the hybrid writer then refuses -- two individually-correct behaviours composing into a hole.
    # Only the bare object-expression forms (LevelOfTheory(...)/CompositeLevelOfTheory(...)), which
    # must stay executable Arkane syntax, are frozen verbatim as source text.
    model_chemistry_node = model_chemistry_nodes[0]
    if isinstance(model_chemistry_node, ast.Constant):
        if not isinstance(model_chemistry_node.value, str):
            raise ValueError(f"'{input_py_path}' assigns a non-string literal to 'modelChemistry' "
                             f"({model_chemistry_node.value!r}); expected either a string label or a "
                             f"LevelOfTheory(...)/CompositeLevelOfTheory(...) expression.")
        model_chemistry = model_chemistry_node.value
    else:
        model_chemistry = ast.get_source_segment(input_py_text, model_chemistry_node)
    if model_chemistry is None or not model_chemistry.strip():
        raise ValueError(f"Could not recover the source text of 'modelChemistry' in '{input_py_path}'.")

    literals = dict()
    for field_name in _LITERAL_FIELD_NAMES:
        nodes = assigns_by_name.get(field_name)
        if not nodes:
            literals[field_name] = None
            continue
        try:
            value = ast.literal_eval(nodes[0])
        except (ValueError, TypeError, SyntaxError) as e:
            raise ValueError(f"'{input_py_path}' assigns '{field_name}' to a non-literal value that could not be "
                             f"safely evaluated with ast.literal_eval: {e}.") from e
        _validate_literal_type(field_name=field_name, value=value, input_py_path=input_py_path)
        literals[field_name] = value

    for required_field_name in ('useHinderedRotors', 'useAtomCorrections', 'useBondCorrections'):
        if literals[required_field_name] is None:
            raise ValueError(f"'{input_py_path}' does not assign '{required_field_name}'; ARC always writes this "
                             f"directive, so its absence is ambiguous rather than a meaningful default.")

    use_bond_corrections = literals['useBondCorrections']
    bond_correction_type = literals['bondCorrectionType']
    if use_bond_corrections and bond_correction_type is None:
        raise ValueError(f"'{input_py_path}' sets useBondCorrections = True but does not assign "
                         f"'bondCorrectionType'; Arkane would otherwise silently default it to 'p' "
                         f"(Petersson-type) rather than making that choice explicit.")
    if not use_bond_corrections and bond_correction_type is not None:
        raise ValueError(f"'{input_py_path}' sets useBondCorrections = False but assigns bondCorrectionType = "
                         f"{bond_correction_type!r}; this is contradictory within a single file.")

    with open(output_yml_path, 'r') as f:
        output_yml_text = f.read()
    try:
        output_yml = yaml.safe_load(output_yml_text)
    except yaml.YAMLError as e:
        raise ValueError(f"Could not parse '{output_yml_path}' as YAML: {e}.") from e
    if not isinstance(output_yml, dict):
        raise ValueError(f"'{output_yml_path}' did not parse to a mapping (got {type(output_yml).__name__}).")

    frequency_scale_factor = _cross_validate_frequency_scale_factor(
        input_value=literals['frequencyScaleFactor'], output_yml=output_yml, output_yml_path=output_yml_path)
    atom_energies = _resolve_atom_energies(
        use_atom_corrections=literals['useAtomCorrections'],
        input_atom_energies=literals['atomEnergies'],
        output_yml=output_yml,
        output_yml_path=output_yml_path)
    bond_additivity_corrections = _resolve_bond_additivity_corrections(
        use_bond_corrections=use_bond_corrections,
        bond_correction_type=bond_correction_type,
        output_yml=output_yml,
        output_yml_path=output_yml_path)

    return {
        'model_chemistry': model_chemistry,
        'frequency_scale_factor': frequency_scale_factor,
        'use_hindered_rotors': literals['useHinderedRotors'],
        'use_atom_corrections': literals['useAtomCorrections'],
        'atom_energies': atom_energies,
        'use_bond_corrections': use_bond_corrections,
        'bond_correction_type': bond_correction_type,
        'bond_additivity_corrections': bond_additivity_corrections,
        'source_paths': {
            'input_py': input_py_path,
            'output_yml': output_yml_path,
        },
    }


# The exact type every key of a FROZEN energy-settings block must have, mapped to whether None is
# also allowed for it. This is the single authoritative definition of the frozen block's shape:
# read_arc_energy_settings produces it, the capture manifest stores it, and verify_capture checks a
# manifest read back off disk against it -- a manifest is a file, so it can be hand-edited into any
# shape at all between those two points, which is exactly what verify_capture exists to catch.
_FROZEN_FIELD_TYPES = {
    'model_chemistry': ((str,), False),
    'frequency_scale_factor': ((int, float), True),
    'use_hindered_rotors': ((bool,), False),
    'use_atom_corrections': ((bool,), False),
    'atom_energies': ((dict,), True),
    'use_bond_corrections': ((bool,), False),
    'bond_correction_type': ((str,), True),
    'bond_additivity_corrections': ((dict,), True),
}


def validate_frozen_energy_settings(energy_settings: dict, context: str) -> None:
    """
    Validate a frozen energy-settings block's keys and value types.

    Args:
        energy_settings (dict): The frozen block to validate.
        context (str): What is being validated, for error messages (e.g. "the capture manifest at
                       '<path>'").

    Raises:
        ValueError: If ``energy_settings`` is not a mapping, is missing a required key, or holds a
                   value whose exact type is not the one that key requires.
    """
    if not isinstance(energy_settings, dict):
        raise ValueError(f"The frozen 'energy_settings' block in {context} is not a mapping "
                         f"(got {type(energy_settings).__name__}).")
    for field_name, (expected_types, none_allowed) in _FROZEN_FIELD_TYPES.items():
        if field_name not in energy_settings:
            raise ValueError(f"The frozen 'energy_settings' block in {context} is missing the required key "
                             f"'{field_name}'. Every block is written whole by read_arc_energy_settings, so a "
                             f"missing key means it was truncated or hand-edited.")
        value = energy_settings[field_name]
        expected_names = ' or '.join(t.__name__ for t in expected_types)
        if value is None:
            if none_allowed:
                continue
            raise ValueError(f"The frozen 'energy_settings' block in {context} has '{field_name}' set to None; "
                             f"it must be a {expected_names}.")
        if type(value) not in expected_types:
            raise ValueError(f"The frozen 'energy_settings' block in {context} has '{field_name}' set to a "
                             f"{type(value).__name__} ({value!r}); expected a {expected_names}. A wrongly typed "
                             f"frozen setting is refused rather than coerced: a non-empty string such as 'False' "
                             f"is truthy, so it would slip past every truthiness guard downstream while rendering "
                             f"into the generated Arkane input as a bare, falsy token.")
    if not energy_settings['model_chemistry'].strip():
        raise ValueError(f"The frozen 'energy_settings' block in {context} has a blank 'model_chemistry'. Without "
                         f"it, Arkane cannot apply atom energy corrections, so a QM'd transition state's E0 would "
                         f"not be on the same energy reference scale as the RMG wells around it.")


def _validate_literal_type(field_name: str, value, input_py_path: str) -> None:
    """
    Enforce that a literal directive read from ``input.py`` has exactly the type
    ``_LITERAL_FIELD_TYPES`` requires.

    Type is checked by exact ``type(value)`` membership rather than ``isinstance``, deliberately.
    ``isinstance(True, int)`` is ``True`` in Python's numeric tower, so an ``isinstance``-based
    check would let ``frequencyScaleFactor = True`` through and silently scale every frequency by
    1.0; conversely a bool subclass is not a thing ARC ever writes. See ``_LITERAL_FIELD_TYPES``
    for why a wrongly typed value is more dangerous than a missing one.

    Args:
        field_name (str): The directive name, as written in ``input.py``.
        value: The value ``ast.literal_eval`` produced for it.
        input_py_path (str): Used in the error message only.

    Raises:
        ValueError: If ``value``'s exact type is not one of the types required for ``field_name``.
    """
    expected_types = _LITERAL_FIELD_TYPES[field_name]
    if type(value) not in expected_types:
        expected_names = ' or '.join(t.__name__ for t in expected_types)
        raise ValueError(f"'{input_py_path}' assigns '{field_name}' a {type(value).__name__} ({value!r}); ARC "
                         f"always writes a {expected_names} here. Refusing to freeze a wrongly typed directive: a "
                         f"non-empty string such as 'False' is truthy, so every downstream guard would read it as "
                         f"the opposite of what this directive says, while the generated Arkane input would render "
                         f"it back as a bare (and correctly falsy) token.")


def _collect_module_level_assigns(tree: ast.Module) -> dict:
    """
    Collect every module-level ``name = <value>`` assignment in ``tree``, keyed by name.

    Only simple ``ast.Assign`` nodes with a single ``ast.Name`` target are considered;
    ``ast.AnnAssign``/``ast.AugAssign``/tuple-unpacking targets are ignored, since none of the
    fields this module reads are ever written that way by ARC.

    Args:
        tree (ast.Module): The parsed ``input.py`` module.

    Returns:
        dict: Field name -> list of ``ast.expr`` value nodes assigned to that name (more than one
             entry means the field was assigned more than once, which callers treat as ambiguous).
    """
    assigns_by_name = dict()
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name):
                assigns_by_name.setdefault(target.id, list()).append(node.value)
    return assigns_by_name


def _cross_validate_frequency_scale_factor(input_value: float | None, output_yml: dict, output_yml_path: str
                                            ) -> float | None:
    """
    Resolve the frozen ``frequency_scale_factor``, cross-validating ``input.py``'s value (if any)
    against ``output.yml``'s ``freq_scale_factor`` (if any).

    This is one of the few fields safe to cross-validate directly: it is a property of the
    frequency level of theory, not a per-subdir on/off flag, so unlike the atom/bond correction
    dicts a real disagreement here cannot be explained away by a sibling statmech subdir.

    Args:
        input_value (float | None): ``input.py``'s ``frequencyScaleFactor`` value, or ``None`` if
                                    unset.
        output_yml (dict): The parsed ``output.yml`` mapping.
        output_yml_path (str): Used in error messages only.

    Raises:
        ValueError: If ``output.yml``'s ``freq_scale_factor`` is present but not numeric, or if both
                   values are present and numerically disagree.

    Returns:
        float | None: ``input_value`` if set; else ``output.yml``'s value if set; else ``None``.
    """
    output_value = output_yml.get('freq_scale_factor')
    if output_value is not None and not isinstance(output_value, (int, float)):
        raise ValueError(f"'{output_yml_path}' has a non-numeric 'freq_scale_factor': {output_value!r}.")
    if input_value is not None and output_value is not None and not math.isclose(input_value, output_value):
        raise ValueError(f"Frequency scale factor mismatch: input.py sets frequencyScaleFactor = {input_value!r} "
                         f"but '{output_yml_path}' records freq_scale_factor = {output_value!r}.")
    return input_value if input_value is not None else (float(output_value) if output_value is not None else None)


def _resolve_atom_energies(use_atom_corrections: bool, input_atom_energies: dict | None, output_yml: dict,
                            output_yml_path: str) -> dict | None:
    """
    Resolve the frozen ``atom_energies`` mapping.

    "Null means unknown, not off": when ``use_atom_corrections`` is ``True`` but neither
    ``input.py`` nor ``output.yml`` records the actual numeric corrections, this returns ``None``
    without raising -- that is a gap in the recorded data, not evidence corrections were disabled.
    When ``use_atom_corrections`` is ``False``, the corrections are irrelevant to this pass and
    ``None`` is returned regardless of what ``output.yml`` happens to contain (see the module
    docstring's "output.yml is project-global" design note for why a non-null dict here is NOT
    treated as a contradiction).

    Args:
        use_atom_corrections (bool): This subdir's ``useAtomCorrections`` flag.
        input_atom_energies (dict | None): ``input.py``'s own ``atomEnergies`` literal, if any.
        output_yml (dict): The parsed ``output.yml`` mapping.
        output_yml_path (str): Used in error messages only.

    Raises:
        ValueError: If ``output.yml``'s ``atom_energy_corrections`` is present but not a mapping.

    Returns:
        dict | None: ``input_atom_energies`` if set; else ``output.yml``'s
             ``atom_energy_corrections`` when ``use_atom_corrections`` is ``True``; else ``None``.
    """
    if input_atom_energies is not None:
        return input_atom_energies
    if not use_atom_corrections:
        return None
    output_value = output_yml.get('atom_energy_corrections')
    if output_value is not None and not isinstance(output_value, dict):
        raise ValueError(f"'{output_yml_path}' has a non-mapping 'atom_energy_corrections': "
                         f"{type(output_value).__name__}.")
    return output_value


def _resolve_bond_additivity_corrections(use_bond_corrections: bool, bond_correction_type: str | None,
                                          output_yml: dict, output_yml_path: str) -> dict | None:
    """
    Resolve the frozen ``bond_additivity_corrections`` mapping, cross-validating
    ``bond_correction_type`` against ``output.yml``'s ``bac_type`` when this subdir explicitly
    enables bond corrections with a specific type.

    Args:
        use_bond_corrections (bool): This subdir's ``useBondCorrections`` flag.
        bond_correction_type (str | None): This subdir's ``bondCorrectionType``, if any.
        output_yml (dict): The parsed ``output.yml`` mapping.
        output_yml_path (str): Used in error messages only.

    Raises:
        ValueError: If ``output.yml``'s ``bond_additivity_corrections`` is present but not a
                   mapping; or if ``bond_correction_type`` is set and ``output.yml``'s ``bac_type``
                   is also set but they disagree.

    Returns:
        dict | None: ``output.yml``'s ``bond_additivity_corrections`` when ``use_bond_corrections``
             is ``True``; else ``None`` (see the module docstring's "output.yml is project-global"
             design note for why a non-null dict here, when ``use_bond_corrections`` is ``False``,
             is NOT treated as a contradiction).
    """
    output_bac_type = output_yml.get('bac_type')
    if bond_correction_type is not None and output_bac_type is not None and bond_correction_type != output_bac_type:
        raise ValueError(f"Bond correction type mismatch: input.py sets bondCorrectionType = "
                         f"{bond_correction_type!r} but '{output_yml_path}' records bac_type = "
                         f"{output_bac_type!r}.")
    if not use_bond_corrections:
        return None
    output_value = output_yml.get('bond_additivity_corrections')
    if output_value is not None and not isinstance(output_value, dict):
        raise ValueError(f"'{output_yml_path}' has a non-mapping 'bond_additivity_corrections': "
                         f"{type(output_value).__name__}.")
    return output_value
