"""
t3 pdep hybrid module

Writer for a *hybrid* Arkane P-dep network input file: QM statmech (RRKM) for a caller-selected
subset of transition states (TSs), RMG/ILT for all the rest. This is the counterpart of
``t3.utils.writer.write_arkane_network_input_file`` (which produces an all-ILT network) for the
case where some path reactions' transition states have already been refined with QM.

Two invariants this module exists to protect (see the module-level docstrings of
``t3/pdep/selector.py`` and ``t3/pdep/mesolver/arkane.py`` for the surrounding design):

1. **Energy reference.** A raw QM energy is not, by itself, on the same enthalpy-of-formation
   scale as an RMG well's ``E0``. Arkane's atom energy corrections (configured via
   ``modelChemistry``/``atomEnergies``) are what make a QM energy consistent with that scale
   (``arkane/encorr/corr.py``). An RMG-generated network file carries none of these directives,
   so this module injects them; it never emits ``useAtomCorrections = False``, since that
   directive silently disables the correction rather than loudly failing without it.
2. **Fail closed on a bad TS.** If a QM'd path reaction kept its ``kinetics = Arrhenius(...)``
   entry, Arkane would silently fall back to ILT for it if the statmech solve failed, and a
   caller would believe it got QM when it did not. So every QM'd path reaction's ``kinetics =
   ...`` entry is dropped (replaced with ``tunneling = '...'`` when requested).

This module never ``exec``s or ``import``s the network/artifact files it manipulates: it reads
them as text and inspects them with ``ast.parse`` only, mirroring the house pattern in
``t3/pdep/parser.py``.
"""

import ast
import hashlib
import logging
import os
import re
import shutil
import tempfile
import warnings
from dataclasses import dataclass, field
from dataclasses import fields as dataclass_fields

from t3.pdep.parser import (NetworkTextUnparseable, TGridClampRecord, format_skipped_species,
                            network_thermo_t_max, parse_pdep_network_file)
from t3.utils.writer import METHOD_LINE_CANDIDATE_RE, format_clamped_t_max, rewrite_arkane_method_line

logger = logging.getLogger(__name__)

# Arkane's transitionState() DSL call: transitionState(label, path) is the by-reference form;
# transitionState(label, E0=..., ...) is the inline form. Mixing the two on one call raises
# InputError ("can only link a quantum job or directly input information, not both"), which is
# exactly why this module never emits inline kwargs alongside a path.
_TOP_LEVEL_CALL_NAMES = ('species', 'transitionState', 'reaction')

# Characters that would let a crafted QMEnergySettings string value break out of (or corrupt) the
# quoted literal it is interpolated into in the generated, executable Arkane input; see
# ``_validate_no_injection_chars``.
_INJECTION_CHARS = ('\n', '\r', "'", '"', '\\')


def _validate_no_injection_chars(field_name: str, value: str) -> None:
    """
    Reject a string destined for raw interpolation into the generated Arkane input if it contains
    a newline, quote character or backslash, any of which could inject or corrupt directives in
    that (executable, Arkane-parsed) file.

    Args:
        field_name (str): The name of the field ``value`` came from (used in the error message).
        value (str): The string to validate.

    Raises:
        ValueError: If ``value`` contains a newline, quote character or backslash.
    """
    if any(char in value for char in _INJECTION_CHARS):
        raise ValueError(f"QMEnergySettings.{field_name} must not contain a newline, quote character or "
                         f"backslash, got {value!r}.")


# The two Arkane DSL call forms a modelChemistry directive may legitimately use in place of a
# plain string label; see arkane/encorr/models.py. Only these two bare (unquoted) object
# expressions are ever allowed through _validate_model_chemistry_expression -- anything else that
# happens to parse as an ast.Call is rejected just like any other non-object-expression string.
_MODEL_CHEMISTRY_CALL_NAMES = ('LevelOfTheory', 'CompositeLevelOfTheory')

# The real Arkane LevelOfTheory constructor schema (arkane/modelchem.py, class LevelOfTheory,
# lines 98-121 as of this writing: ``method: str``, ``basis: str = None``,
# ``auxiliary_basis: str = None``, ``cabs: str = None``, ``software: str = None``,
# ``software_version: Union[int, float, str] = None``, ``solvent: str = None``,
# ``solvation_method: str = None``, ``args: Union[str, Iterable[str]] = None``). This maps each
# real field name to the exact Python type(s) its literal value is allowed to have: every field is
# str-only EXCEPT ``software_version`` (int/float/str). ``args`` genuinely also accepts an iterable
# of str (e.g. a tuple), but that collection form is deliberately NOT accepted here (kept str-only):
# T3/ARC never emit anything but a scalar for it, and there is no concrete real-world case to widen
# the allowlist for -- narrowing the accepted surface below the real schema is always safe, only
# widening it is not. 'method' has no default in the dataclass, i.e. it is required.
_LEVEL_OF_THEORY_REQUIRED_FIELDS = ('method',)
_LEVEL_OF_THEORY_FIELD_TYPES = {
    'method': (str,),
    'basis': (str,),
    'auxiliary_basis': (str,),
    'cabs': (str,),
    'software': (str,),
    'software_version': (int, float, str),
    'solvent': (str,),
    'solvation_method': (str,),
    'args': (str,),
}

# The real Arkane CompositeLevelOfTheory constructor schema (arkane/modelchem.py, class
# CompositeLevelOfTheory, lines 174-191 as of this writing): exactly two fields, ``freq`` and
# ``energy``, both of type ``LevelOfTheory`` and both required (no defaults) -- a
# CompositeLevelOfTheory with either one missing is vacuous/malformed.
_COMPOSITE_LEVEL_OF_THEORY_REQUIRED_FIELDS = ('freq', 'energy')


def _model_chemistry_ast_call(value: str) -> ast.Call | None:
    """
    Classify a ``model_chemistry`` value as either a plain string label or a
    ``LevelOfTheory(...)``/``CompositeLevelOfTheory(...)`` bare object-expression.

    ``ast.parse(value, mode='eval')`` is used purely as a classifier here: a plain label like
    ``'wb97xd/def2tzvp'`` also happens to parse (as an ``ast.BinOp``, dividing two ``ast.Name``
    nodes), so the discriminator is specifically "is the parsed expression an ``ast.Call`` to one
    of the two known constructor names", not "did it parse at all". A value that fails to parse at
    all (e.g. it contains a newline) is treated as a plain string label, deferring to
    ``_validate_no_injection_chars`` for that case.

    Args:
        value (str): The ``model_chemistry`` value to classify.

    Returns:
        ast.Call | None: The parsed call node if ``value`` is a ``LevelOfTheory(...)`` or
                         ``CompositeLevelOfTheory(...)`` object-expression, else ``None``.
    """
    try:
        node = ast.parse(value, mode='eval').body
    except SyntaxError:
        return None
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in _MODEL_CHEMISTRY_CALL_NAMES:
        return node
    return None


def _validate_level_of_theory_call(call: ast.Call, value: str) -> None:
    """
    Structurally validate a ``LevelOfTheory(...)`` call node against the real Arkane
    ``LevelOfTheory`` constructor schema (see ``_LEVEL_OF_THEORY_FIELD_TYPES``): only keyword
    arguments (no positional args, no ``**kwargs``, no duplicates), each keyword restricted to a
    real ``LevelOfTheory`` field name, each value a literal ``ast.Constant`` of the exact Python
    type that field accepts, and the required ``method`` keyword present. An AST-structural check
    alone (any keyword name, any ``ast.Constant`` type) is not arbitrary-code-safe as a validator:
    it would let a keyword like ``__class__`` or a value like ``method=True``/``method=...``
    through, since those only fail once Arkane itself executes and constructs the object.

    Args:
        call (ast.Call): The ``LevelOfTheory(...)`` call node to validate.
        value (str): The original ``model_chemistry`` string (used in error messages).

    Raises:
        ValueError: If ``call`` has positional args, ``**kwargs``, a duplicate keyword, a keyword
                   that is not a real ``LevelOfTheory`` field, a keyword value that is not a
                   literal constant of the type that field accepts, or is missing the required
                   ``method`` keyword.
    """
    if call.args:
        raise ValueError(f"QMEnergySettings.model_chemistry's LevelOfTheory(...) call must not have positional "
                         f"arguments, got {value!r}.")
    seen_keywords = set()
    for kw in call.keywords:
        if kw.arg is None:
            raise ValueError(f"QMEnergySettings.model_chemistry's LevelOfTheory(...) call must not use "
                             f"**kwargs, got {value!r}.")
        if kw.arg in seen_keywords:
            raise ValueError(f"QMEnergySettings.model_chemistry's LevelOfTheory(...) call must not repeat "
                             f"keyword {kw.arg!r}, got {value!r}.")
        seen_keywords.add(kw.arg)
        allowed_types = _LEVEL_OF_THEORY_FIELD_TYPES.get(kw.arg)
        if allowed_types is None:
            raise ValueError(f"QMEnergySettings.model_chemistry's LevelOfTheory(...) call has unknown keyword "
                             f"{kw.arg!r}, got {value!r}.")
        # Compare exact type(), not isinstance(): isinstance(True, int) is True in Python's numeric
        # tower, so an isinstance check would silently accept e.g. method=True as if it were a str.
        if not isinstance(kw.value, ast.Constant) or type(kw.value.value) not in allowed_types:
            raise ValueError(f"QMEnergySettings.model_chemistry's LevelOfTheory(...) keyword {kw.arg!r} must be a "
                             f"literal of type {allowed_types!r}, got {value!r}.")
    missing = [field for field in _LEVEL_OF_THEORY_REQUIRED_FIELDS if field not in seen_keywords]
    if missing:
        raise ValueError(f"QMEnergySettings.model_chemistry's LevelOfTheory(...) call is missing required "
                         f"keyword(s) {missing!r}, got {value!r}.")


def _validate_composite_level_of_theory_call(call: ast.Call, value: str) -> None:
    """
    Structurally validate a ``CompositeLevelOfTheory(...)`` call node against the real Arkane
    ``CompositeLevelOfTheory`` constructor schema (see
    ``_COMPOSITE_LEVEL_OF_THEORY_REQUIRED_FIELDS``): only ``freq=``/``energy=`` keywords (no
    positional args, no ``**kwargs``, no duplicates), each a nested, itself-valid
    ``LevelOfTheory(...)`` call, and BOTH required keywords present -- ``freq`` and ``energy`` have
    no defaults in Arkane's dataclass, so e.g. a bare ``CompositeLevelOfTheory()`` with neither is
    vacuous and must be rejected, not silently accepted as structurally valid.

    Args:
        call (ast.Call): The ``CompositeLevelOfTheory(...)`` call node to validate.
        value (str): The original ``model_chemistry`` string (used in error messages).

    Raises:
        ValueError: If ``call`` has positional args, ``**kwargs``, a duplicate keyword, a keyword
                   other than ``freq``/``energy``, a keyword value that is not a valid
                   ``LevelOfTheory(...)`` call, or is missing ``freq`` and/or ``energy``.
    """
    if call.args:
        raise ValueError(f"QMEnergySettings.model_chemistry's CompositeLevelOfTheory(...) call must not have "
                         f"positional arguments, got {value!r}.")
    seen_keywords = set()
    for kw in call.keywords:
        if kw.arg not in _COMPOSITE_LEVEL_OF_THEORY_REQUIRED_FIELDS:
            raise ValueError(f"QMEnergySettings.model_chemistry's CompositeLevelOfTheory(...) call only accepts "
                             f"'freq'/'energy' keywords, got {kw.arg!r} in {value!r}.")
        if kw.arg in seen_keywords:
            raise ValueError(f"QMEnergySettings.model_chemistry's CompositeLevelOfTheory(...) call must not "
                             f"repeat keyword {kw.arg!r}, got {value!r}.")
        seen_keywords.add(kw.arg)
        if not (isinstance(kw.value, ast.Call) and isinstance(kw.value.func, ast.Name)
                and kw.value.func.id == 'LevelOfTheory'):
            raise ValueError(f"QMEnergySettings.model_chemistry's CompositeLevelOfTheory(...) keyword {kw.arg!r} "
                             f"must be a LevelOfTheory(...) call, got {value!r}.")
        _validate_level_of_theory_call(kw.value, value)
    missing = [field for field in _COMPOSITE_LEVEL_OF_THEORY_REQUIRED_FIELDS if field not in seen_keywords]
    if missing:
        raise ValueError(f"QMEnergySettings.model_chemistry's CompositeLevelOfTheory(...) call is missing "
                         f"required keyword(s) {missing!r}, got {value!r}.")


def _validate_model_chemistry_expression(field_name: str, value: str) -> None:
    """
    Validate a ``model_chemistry`` value destined for the generated (executable) Arkane input,
    accepting either a plain string label (subject to ``_validate_no_injection_chars``) or a
    structurally-validated ``LevelOfTheory(...)``/``CompositeLevelOfTheory(...)`` bare object
    expression, which ARC's real ``modelChemistry`` directives commonly use and which
    ``_validate_no_injection_chars`` alone would incorrectly reject (these expressions legitimately
    contain quote characters).

    Args:
        field_name (str): The name of the field ``value`` came from (used in error messages).
        value (str): The ``model_chemistry`` value to validate.

    Raises:
        ValueError: If ``value`` is a plain string label containing a newline, quote character or
                   backslash; or if it is a ``LevelOfTheory(...)``/``CompositeLevelOfTheory(...)``
                   call that fails structural validation (positional args, ``**kwargs``,
                   unexpected keyword, or a non-literal/non-nested-call keyword value).
    """
    call = _model_chemistry_ast_call(value)
    if call is None:
        _validate_no_injection_chars(field_name, value)
        return
    if call.func.id == 'LevelOfTheory':
        _validate_level_of_theory_call(call, value)
    else:
        _validate_composite_level_of_theory_call(call, value)


# Every QMEnergySettings field carries its own expected-type metadata (see the ``field(...)``
# calls in the dataclass below) rather than being validated against a second, hand-maintained
# table here. QMEnergySettings is a plain dataclass, so its annotations enforce NOTHING at
# runtime, and the values reaching it come from a YAML-loaded capture manifest that can carry any
# type at all. A wrongly typed value here is more dangerous than a missing one: a non-empty string
# such as 'False' is TRUTHY, so `if not energy_settings.use_atom_corrections` never fires -- while
# _build_energy_header's f-string renders it straight back out as a bare `useAtomCorrections =
# False` that Arkane evaluates as the boolean. The guard is bypassed and the directive inverted in
# a single step, and the resulting hybrid network is silently on the wrong energy scale. Types are
# compared by exact type() membership, not isinstance: isinstance(True, int) is True in Python's
# numeric tower, so an isinstance check would accept frequency_scale_factor=True and silently
# scale every frequency by 1.0.
def _validate_energy_settings_types(energy_settings) -> None:
    """
    Enforce that every ``QMEnergySettings`` field has exactly the type it claims, before any other
    guard reads it.

    This runs FIRST, ahead of every truthiness-based guard in
    ``write_hybrid_network_input_file``, because those guards are precisely what a wrongly typed
    value defeats. The type table consulted here is ``QMEnergySettings.field_types()`` -- derived
    from each field's own ``metadata=...`` -- so it is the SAME table
    ``t3.pdep.energy_settings.validate_frozen_energy_settings`` consults for the frozen capture
    manifest; there is no second, hand-maintained copy to drift out of sync.

    Args:
        energy_settings (QMEnergySettings): The settings to type-check.

    Raises:
        ValueError: If any field's exact type is not one this dataclass permits for it.
    """
    for field_name, (expected_types, none_allowed) in QMEnergySettings.field_types().items():
        value = getattr(energy_settings, field_name)
        if value is None:
            if none_allowed:
                continue
            raise ValueError(f"QMEnergySettings.{field_name} must not be None; it must be a "
                             f"{' or '.join(t.__name__ for t in expected_types)}.")
        if type(value) not in expected_types:
            raise ValueError(f"QMEnergySettings.{field_name} must be a "
                             f"{' or '.join(t.__name__ for t in expected_types)}, got a "
                             f"{type(value).__name__} ({value!r}). A wrongly typed setting is refused rather than "
                             f"coerced: a non-empty string such as 'False' is truthy, so it would slip past the "
                             f"guards below while rendering into the generated Arkane input as a bare, falsy token.")


@dataclass(frozen=True)
class QMEnergySettings:
    """
    The energy-reference and tunneling settings to inject into a hybrid network input file.

    Args:
        model_chemistry (str): The Arkane ``modelChemistry`` label (e.g. ``'wb97xd/def2tzvp'``).
                               Required and must not be blank: this is what lets Arkane apply
                               atom energy corrections so a QM'd TS's E0 lands on the same
                               enthalpy-of-formation scale as the RMG wells around it.
        frequency_scale_factor (float, optional): Arkane's ``frequencyScaleFactor``. Omitted from
                                                  the header entirely when ``None``.
        use_hindered_rotors (bool, optional): Arkane's ``useHinderedRotors``. Default: ``True``.
        use_bond_corrections (bool, optional): Arkane's ``useBondCorrections``. Default: ``False``.
                                               Must not be ``True``: Arkane's input DSL has no
                                               directive to pin bond additivity correction VALUES
                                               (``useBondCorrections``/``bondCorrectionType`` only
                                               select whether/which BAC SCHEME applies; the actual
                                               numeric corrections come from Arkane's own database
                                               at solve time). ``write_hybrid_network_input_file``
                                               raises ``ValueError`` rather than honoring ``True``.
        bond_correction_type (str, optional): Arkane's ``bondCorrectionType`` ('p' for
                                              Petersson-type or 'm' for Melius-type BAC; see
                                              ``arkane/input.py``). Must be ``None`` when
                                              ``use_bond_corrections`` is ``False`` (its only
                                              allowed value, since ``use_bond_corrections=True`` is
                                              itself always refused). Default: ``None``.
        atom_energies (dict, optional): Arkane's ``atomEnergies`` mapping -- the only directive in
                                        Arkane's input DSL that pins atom energy correction VALUES
                                        (as opposed to merely turning corrections on/off). Must not
                                        be ``None`` when ``use_atom_corrections`` is ``True``:
                                        ``write_hybrid_network_input_file`` raises ``ValueError``
                                        rather than emitting ``useAtomCorrections = True`` with no
                                        ``atomEnergies = ...`` to back it.
        use_atom_corrections (bool, optional): Arkane's ``useAtomCorrections``. Default: ``True``.
                                               Must not be ``False``: that directive silently
                                               disables the atom energy corrections that put a
                                               QM'd TS's E0 on the same enthalpy-of-formation scale
                                               as the RMG wells around it (see the module
                                               docstring); ``write_hybrid_network_input_file``
                                               raises ``ValueError`` rather than honoring
                                               ``False``.
        bond_additivity_corrections (dict, optional): The BAC VALUES ARC's ``output.yml`` records
                                                      (as opposed to ``bond_correction_type``,
                                                      which merely names a SCHEME). PROVENANCE
                                                      ONLY: Arkane's input DSL has no directive
                                                      that pins bond additivity correction values
                                                      the way ``atomEnergies`` pins atom energy
                                                      corrections, so this field is never rendered
                                                      by ``_build_energy_header`` and
                                                      ``write_hybrid_network_input_file`` refuses
                                                      ``use_bond_corrections=True`` unconditionally
                                                      regardless of what this holds. Carried through
                                                      purely so a frozen capture manifest does not
                                                      silently drop a value ARC computed. Default:
                                                      ``None``.
        tunneling (str, optional): The tunneling model name (e.g. ``'Eckart'``) written onto every
                                   QM'd reaction's ``reaction(...)`` block. Set to ``None`` to
                                   omit ``tunneling = ...`` entirely. Default: ``'Eckart'``. Left
                                   at ``'Eckart'`` by default (rather than omitted) intentionally:
                                   Arkane's Eckart tunneling raises loudly on a negative barrier
                                   height (``rmgpy/kinetics/tunneling.pyx``) instead of silently
                                   proceeding without a needed correction. Set this to ``None`` as
                                   an explicit escape hatch if that loud failure is not wanted.
                                   Unlike every other field above, this one is write-time-only: it
                                   is supplied directly by the caller of
                                   ``write_hybrid_network_input_file`` and never appears in a
                                   frozen capture manifest's ``energy_settings`` block (see
                                   ``field_types(frozen=True)``).
    """
    model_chemistry: str = field(
        metadata={'types': (str,), 'none_allowed': False, 'in_frozen_manifest': True})
    frequency_scale_factor: float | None = field(
        default=None, metadata={'types': (int, float), 'none_allowed': True, 'in_frozen_manifest': True})
    use_hindered_rotors: bool = field(
        default=True, metadata={'types': (bool,), 'none_allowed': False, 'in_frozen_manifest': True})
    use_bond_corrections: bool = field(
        default=False, metadata={'types': (bool,), 'none_allowed': False, 'in_frozen_manifest': True})
    bond_correction_type: str | None = field(
        default=None, metadata={'types': (str,), 'none_allowed': True, 'in_frozen_manifest': True})
    atom_energies: dict | None = field(
        default=None, metadata={'types': (dict,), 'none_allowed': True, 'in_frozen_manifest': True})
    use_atom_corrections: bool = field(
        default=True, metadata={'types': (bool,), 'none_allowed': False, 'in_frozen_manifest': True})
    bond_additivity_corrections: dict | None = field(
        default=None, metadata={'types': (dict,), 'none_allowed': True, 'in_frozen_manifest': True})
    tunneling: str | None = field(
        default='Eckart', metadata={'types': (str,), 'none_allowed': True, 'in_frozen_manifest': False})

    @classmethod
    def field_types(cls, *, frozen: bool = False) -> dict:
        """
        The single authoritative type table for every field this dataclass declares.

        This is the ONE definition of the settings' shape. Both the write-time guard
        (``_validate_energy_settings_types``, run first inside
        ``write_hybrid_network_input_file``) and the frozen-capture-manifest guard
        (``t3.pdep.energy_settings.validate_frozen_energy_settings``, run by ``verify_capture``)
        read this same table -- derived from each field's own ``metadata=...`` above -- instead of
        each maintaining its own hand-written copy. Adding a field to this dataclass (with its
        ``metadata``) is therefore the only edit needed to validate it at both boundaries.

        Args:
            frozen (bool): If True, omit fields whose ``metadata['in_frozen_manifest']`` is False
                          -- i.e. fields that are write-time-only (supplied directly by the caller
                          of ``write_hybrid_network_input_file``, such as ``tunneling``) and never
                          appear in the frozen block a capture manifest stores.

        Returns:
            dict: Field name -> (expected types tuple, whether ``None`` is allowed).
        """
        return {f.name: (f.metadata['types'], f.metadata['none_allowed'])
                for f in dataclass_fields(cls)
                if not frozen or f.metadata.get('in_frozen_manifest', True)}

    @classmethod
    def from_frozen(cls, frozen: dict, *, tunneling: str | None = 'Eckart') -> 'QMEnergySettings':
        """
        Build a validated ``QMEnergySettings`` from a frozen energy-settings block.

        This is the ONE place a frozen block -- as returned by
        ``t3.pdep.energy_settings.read_arc_energy_settings`` and stored verbatim under a capture
        manifest's ``'energy_settings'`` key -- is turned into the object
        ``write_hybrid_network_input_file`` actually consumes. Every frozen key this dataclass
        recognizes (``field_types(frozen=True)``) is threaded through structurally, so a caller
        never has to remember to map each field by hand (and cannot forget one, the way a manual
        ``QMEnergySettings(model_chemistry=..., atom_energies=..., ...)`` call site could).

        Args:
            frozen (dict): The frozen energy-settings block. Must contain every key
                          ``field_types(frozen=True)`` names. Any other key (e.g. ``source_paths``)
                          is ignored: it is manifest provenance, not a ``QMEnergySettings`` field.
            tunneling (str, optional): The tunneling model to use. Write-time-only, never part of
                                       a frozen block -- see the class docstring. Default:
                                       ``'Eckart'``.

        Returns:
            QMEnergySettings: Validated against ``field_types()``.

        Raises:
            ValueError: If ``frozen`` is missing a required key, or any value's exact type is not
                       one ``field_types()`` permits for it.
        """
        kwargs = {}
        for name in cls.field_types(frozen=True):
            if name not in frozen:
                raise ValueError(f"The frozen energy-settings block is missing the required key '{name}'.")
            kwargs[name] = frozen[name]
        instance = cls(tunneling=tunneling, **kwargs)
        _validate_energy_settings_types(instance)
        return instance


@dataclass(frozen=True)
class HybridNetworkResult:
    """
    The outcome of ``write_hybrid_network_input_file``.

    Args:
        dest_path (str): The path the hybrid Arkane input file was written to.
        qm_ts_labels (tuple): The transition state labels switched to QM/RRKM.
        ilt_ts_labels (tuple): The transition state labels left as RMG/ILT.
        vendored_files (tuple): The paths (relative to ``os.path.dirname(dest_path)``) written
                                under the ``qm/`` subdirectory: one ``qm/<label>.py`` per QM'd TS,
                                plus one ``qm/logs/<label>/<basename>`` per ``Log(...)`` file it
                                references.
        warnings (tuple): Human-readable warnings (e.g. a QM'd path reaction that carried no
                          ``kinetics = ...`` entry to drop in the first place).
        t_grid_clamp (TGridClampRecord): Provenance for whether the T grid written to
                                        ``dest_path`` was clamped down from what the source
                                        network requested. See ``TGridClampRecord``'s docstring
                                        for the full three-state design rationale (mirrors
                                        ``t3.utils.writer.ArkaneNetworkWriteResult``'s identical
                                        field).
    """
    dest_path: str
    qm_ts_labels: tuple
    ilt_ts_labels: tuple
    vendored_files: tuple
    warnings: tuple = field(default_factory=tuple)
    t_grid_clamp: TGridClampRecord | None = None


def write_hybrid_network_input_file(source_path: str,
                                    dest_path: str,
                                    method: str,
                                    qm_transition_states: dict,
                                    energy_settings: QMEnergySettings,
                                    sensitivity: bool = False,
                                    qm_artifacts_root: str | None = None,
                                    ) -> HybridNetworkResult:
    """
    Write a hybrid Arkane P-dep network input file: QM/RRKM for the requested TSs, RMG/ILT for
    the rest.

    Args:
        source_path (str): The path to the original RMG P-dep network file to build from.
        dest_path (str): The path to write the resulting hybrid Arkane input file to. The parent
                         directory is created if it does not already exist. A ``qm/`` subdirectory
                         next to it is populated with the vendored QM artifacts (see
                         ``HybridNetworkResult.vendored_files``).
        method (str): 'CSE', 'MSC' or 'RS' (see ``t3.utils.writer.METHOD_MAP``).
        qm_transition_states (dict): Network transition state label -> path to the ARC-written
                                    statmech artifact file for that TS (Arkane's
                                    ``species_input_template`` shape; see
                                    ``arc/statmech/arkane.py``).
        energy_settings (QMEnergySettings): The energy-reference/tunneling settings to inject.
        sensitivity (bool, optional): Whether to inject a ``sensitivity_conditions`` directive
                                     spanning the network's T/P extrema, mirroring
                                     ``write_arkane_network_input_file``. Default: ``False``.
        qm_artifacts_root (str, optional): The directory every artifact's resolved ``Log(...)``
                                          path must live under -- typically the ARC project
                                          directory the artifacts came from. When ``None``, each
                                          artifact is confined to its own parent directory (the
                                          narrowest, fail-closed default). See
                                          ``_read_qm_artifact`` for why confinement is
                                          non-negotiable: the vendoring step COPIES these files.

    Raises:
        ValueError: If ``energy_settings.model_chemistry`` is missing/blank or fails validation; if
                   ``energy_settings.use_atom_corrections`` is ``False``; if
                   ``energy_settings.use_atom_corrections`` is ``True`` and ``atom_energies`` is
                   ``None``; if ``energy_settings.use_bond_corrections`` is ``True``; if
                   ``qm_transition_states`` is empty; if a key of ``qm_transition_states`` is not
                   a transition state label in the source network; if a QM artifact file (or a
                   ``Log(...)`` file it references) does not exist; if a ``Log(...)`` path
                   resolves outside the allowed root, or its argument is not a literal string;
                   or if the source network cannot be parsed.
        RuntimeError: If the generated network file text, or any vendored QM artifact's text,
                     fails its own ``ast.parse(...)`` self-check. This should never happen in
                     practice (every edit is computed structurally from the AST), but a failure
                     here is caught loudly, before anything is written to ``dest_path``, rather
                     than silently handing back a broken hybrid network.

    Returns:
        HybridNetworkResult: The outcome of the write.
    """
    _validate_energy_settings_types(energy_settings)
    if not energy_settings.model_chemistry or not energy_settings.model_chemistry.strip():
        raise ValueError("QMEnergySettings.model_chemistry is required and must not be blank: without it, Arkane "
                         "cannot apply atom energy corrections, so a QM'd transition state's E0 would not be on "
                         "the same energy reference scale as the RMG wells around it.")
    _validate_model_chemistry_expression('model_chemistry', energy_settings.model_chemistry)
    if not energy_settings.use_atom_corrections:
        raise ValueError("QMEnergySettings.use_atom_corrections is False: this directive silently disables "
                         "Arkane's atom energy corrections, so a QM'd transition state's E0 would not be on the "
                         "same energy reference scale as the RMG wells around it. Leave use_atom_corrections=True "
                         "(the default).")
    if energy_settings.atom_energies is None:
        raise ValueError("QMEnergySettings.use_atom_corrections is True but atom_energies is None: Arkane's "
                         "atomEnergies directive is the ONLY place in its input DSL where atom energy correction "
                         "VALUES are pinned (useAtomCorrections merely turns the mechanism on/off). Leaving "
                         "atom_energies unset would emit useAtomCorrections = True with no atomEnergies = ... to "
                         "back it, silently falling through to whatever corrections Arkane's own database happens "
                         "to have for this model chemistry -- exactly the drift that can put a QM'd transition "
                         "state's E0 off the same enthalpy-of-formation scale as the RMG wells around it. Provide "
                         "the ARC-computed atom_energies dict for this model chemistry.")
    if energy_settings.tunneling is not None:
        _validate_no_injection_chars('tunneling', energy_settings.tunneling)
    if energy_settings.use_bond_corrections:
        raise ValueError("QMEnergySettings.use_bond_corrections is True: Arkane's input DSL has no directive to "
                         "pin bond additivity correction VALUES -- useBondCorrections only turns the mechanism "
                         "on, and bondCorrectionType only selects WHICH scheme ('p' Petersson-type or 'm' "
                         "Melius-type) to apply; the actual numeric corrections are looked up from Arkane's own "
                         "database at solve time, with no way for this file to pin (or even record) what was "
                         "actually applied. That is exactly the drift this hybrid writer exists to prevent for a "
                         "QM'd transition state's E0. Leave use_bond_corrections=False (the default) -- the escape "
                         "hatch this guard names.")
    if not energy_settings.use_bond_corrections and energy_settings.bond_correction_type is not None:
        raise ValueError("QMEnergySettings.bond_correction_type is set but use_bond_corrections is False, so it "
                         "would silently have no effect; either set use_bond_corrections=True or leave "
                         "bond_correction_type unset.")
    if not qm_transition_states:
        raise ValueError('qm_transition_states is empty; use t3.utils.writer.write_arkane_network_input_file(...) '
                         'instead to write an all-ILT (non-hybrid) Arkane network input file.')

    network = parse_pdep_network_file(source_path)
    known_labels = set(network.transition_state_labels)
    unknown_labels = sorted(set(qm_transition_states) - known_labels)
    if unknown_labels:
        raise ValueError(f"qm_transition_states references unknown transition state label(s) "
                         f"{unknown_labels} for network '{network.network_id}'; known transition state labels are "
                         f"{sorted(known_labels)}.")
    for ts_label, artifact_path in qm_transition_states.items():
        if not os.path.isfile(artifact_path):
            raise ValueError(f"The QM statmech artifact for transition state '{ts_label}' was not found at "
                             f"'{artifact_path}'.")

    # Pre-flight: parse every QM artifact and resolve every Log(...) path it references BEFORE
    # writing anything to disk, so a missing log file can never result in a partially-written,
    # dangling hybrid network.
    artifact_infos = {
        ts_label: _read_qm_artifact(
            ts_label=ts_label,
            artifact_path=artifact_path,
            allowed_log_root=qm_artifacts_root if qm_artifacts_root is not None
            else os.path.dirname(os.path.abspath(artifact_path)),
        )
        for ts_label, artifact_path in qm_transition_states.items()}

    with open(source_path, 'r') as f:
        text = f.read()

    tree = _parse_as_ast(text=text, path=source_path)
    line_starts = _line_start_offsets(text)

    edits = list()
    warnings_list = list()
    first_species_lineno = None

    for node in tree.body:
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
            continue
        call = node.value
        call_name = _get_call_name(call)
        if call_name not in _TOP_LEVEL_CALL_NAMES:
            continue

        # NOTE on the two ``if kw.arg is not None`` comprehensions below: they discard a ``**kwargs``
        # unpacking, which by itself would be the fail-open this module refuses everywhere else --
        # ``transitionState(**payload)`` would read as label None (so the QM artifact is never vendored
        # in) and ``reaction(**payload)`` as a reaction with no transitionState (so its ILT
        # ``kinetics = ...`` entry is never dropped), in both cases writing a file T3 believes carries
        # QM results while it carries RMG estimates.
        #
        # No guard is added here, deliberately, because one would be dead code:
        # ``parse_pdep_network_file(source_path)`` above runs unconditionally before this walk and
        # refuses any ``**kwargs`` unpacking itself, and its RECOGNIZED_TOP_LEVEL_CALLS is a strict
        # SUPERSET of this module's _TOP_LEVEL_CALL_NAMES ({'species', 'transitionState', 'reaction'}),
        # so every call reached here has already been cleared. The correctness of these two
        # comprehensions therefore rests on that upstream refusal rather than on themselves -- which is
        # exactly why it is written down here. If that parse is ever moved, made conditional, or its
        # recognized set narrowed, these two lines become a live fail-open again.
        # Pinned end-to-end by TestKwargsUnpackingIsRefusedByTheRewriterToo in
        # tests/test_pdep/test_hybrid.py, which asserts the composition premise directly.
        if call_name == 'species':
            if first_species_lineno is None:
                first_species_lineno = node.lineno

        elif call_name == 'transitionState':
            kwargs = {kw.arg: kw.value for kw in call.keywords if kw.arg is not None}
            label = _literal_or_none(kwargs.get('label'))
            if label in qm_transition_states:
                start = _pos(text, line_starts, node.lineno, 0)
                end = _consume_trailing_newline(text, _pos(text, line_starts, node.end_lineno, node.end_col_offset))
                edits.append((start, end, f"transitionState('{label}', 'qm/{label}.py')\n"))

        elif call_name == 'reaction':
            keyword_nodes = {kw.arg: kw for kw in call.keywords if kw.arg is not None}
            transition_state_kw = keyword_nodes.get('transitionState')
            if transition_state_kw is None:
                ts_label = None
            else:
                ts_label = _literal_or_none(transition_state_kw.value)
                if ts_label is None and not _is_literal_none(transition_state_kw.value):
                    # Refuse rather than degrade: a transitionState=... value that IS present but
                    # is not a literal (e.g. a bare name or a nested call) cannot be checked
                    # against qm_transition_states, so treating it as absent would silently skip
                    # dropping this reaction's 'kinetics = ...' entry even when the reaction's real
                    # (unevaluable) transition state is actually one of the QM'd ones -- the same
                    # fail-open bug already fixed for 'products='/'isomers='/'reactants=' in
                    # t3/pdep/parser.py and for Log(...) above in this module.
                    unevaluable = ast.get_source_segment(text, transition_state_kw.value)
                    raise ValueError(f"The pdep network file at '{source_path}' declares a "
                                     f"'transitionState' keyword that could not be evaluated "
                                     f"literally: transitionState={unevaluable}. Refusing to treat "
                                     f"it as having no transition state in its place.")
            if ts_label in qm_transition_states:
                kinetics_kw = keyword_nodes.get('kinetics')
                if kinetics_kw is None:
                    warnings_list.append(f"Path reaction referencing QM'd transition state '{ts_label}' carries "
                                         f"no 'kinetics = ...' entry to drop.")
                else:
                    edits.append(_kinetics_removal_edit(text=text,
                                                        line_starts=line_starts,
                                                        kinetics_kw=kinetics_kw,
                                                        call=call,
                                                        tunneling=energy_settings.tunneling))

    if first_species_lineno is not None:
        start = _pos(text, line_starts, first_species_lineno, 0)
        edits.append((start, start, _build_energy_header(energy_settings)))

    # Apply edits back-to-front (by descending start offset) so earlier, not-yet-applied offsets
    # (all computed against the pristine, pre-edit text) stay valid.
    for start, end, replacement in sorted(edits, key=lambda edit: edit[0], reverse=True):
        text = text[:start] + replacement + text[end:]

    text, t_grid_clamp = _rewrite_method_and_sensitivity(text=text, method=method, sensitivity=sensitivity,
                                                         source_path=source_path)

    # Self-check: the generated network file's own text must itself be valid Python before it is
    # ever written to disk. A splice bug elsewhere in this module could otherwise hand back a
    # broken (or subtly wrong) hybrid network without any indication something went wrong.
    try:
        ast.parse(text)
    except SyntaxError as e:
        raise RuntimeError(f"Refusing to write '{dest_path}': the generated network file text failed "
                           f"its own self-check (it is not valid Python): {e}.") from e

    dest_dir = os.path.dirname(dest_path)
    if dest_dir and not os.path.isdir(dest_dir):
        os.makedirs(dest_dir)

    # Vendor the QM artifacts (and their referenced logs) BEFORE writing dest_path. Vendoring
    # touches multiple files (one qm/<label>.py plus one or more qm/logs/<label>/<basename> per
    # TS) and can fail partway through (e.g. a log file copy raising); writing dest_path first
    # would leave behind an input.py whose 'qm/<label>.py' transitionState(...) references point
    # at a missing or incomplete qm/ tree. Vendoring first means a failure here leaves no
    # dest_path at all, rather than a dangling one.
    qm_dir = os.path.join(dest_dir, 'qm')
    vendored_files = _vendor_qm_artifacts(artifact_infos=artifact_infos, qm_dir=qm_dir, dest_dir=dest_dir)

    # Self-check: every vendored TS artifact must also be valid Python, for the same reason.
    for ts_label in artifact_infos:
        vendored_ts_path = os.path.join(qm_dir, f'{ts_label}.py')
        with open(vendored_ts_path, 'r') as f:
            vendored_text = f.read()
        try:
            ast.parse(vendored_text)
        except SyntaxError as e:
            raise RuntimeError(f"Refusing to write '{dest_path}': the vendored QM artifact for transition "
                               f"state '{ts_label}' at '{vendored_ts_path}' failed its own self-check "
                               f"(it is not valid Python): {e}.") from e

    # Written atomically -- staged into a temp file in the SAME directory as dest_path, fsync'd,
    # then os.replace'd into place, mirroring t3.main.T3._mark_arc_finalization_complete (which
    # exists precisely so a crash cannot leave a marker that certifies work that did not finish).
    # Without this, a plain open(dest_path, 'w'); f.write(text) can be interrupted mid-write (disk
    # full, OOM-kill, power loss), leaving a TORN, partially-written input.py at dest_path; on the
    # very next pass, _mark_arc_finalization_complete then runs and marks the iteration fully
    # finalized, so check_arc_finalization_complete() would report success while a downstream
    # consumer reads a truncated Arkane input as if T3 stood behind it -- the marker becoming MORE
    # durable than the very output it is meant to certify. The temp file is staged inside dest_dir
    # (the network's own output directory), never directly under 'PDep hybrid', so it can never be
    # mistaken for one of the top-level per-network entries T3._prune_stale_pdep_hybrid_outputs
    # inspects (and raises on if not a plain directory) -- it lives one level deeper.
    staging_dir = dest_dir or '.'
    fd, staged_path = tempfile.mkstemp(prefix='.hybrid-input-', suffix='.py', dir=staging_dir)
    try:
        with os.fdopen(fd, 'w') as f:
            f.write(text)
            f.flush()
            os.fsync(f.fileno())
        os.replace(staged_path, dest_path)
    except Exception:
        if os.path.isfile(staged_path):
            os.remove(staged_path)
        raise
    else:
        # fsync the containing directory too, so the rename itself (not just the file's bytes) is
        # durable before _mark_arc_finalization_complete() can run and certify this write as done.
        dir_fd = os.open(staging_dir, os.O_RDONLY)
        try:
            os.fsync(dir_fd)
        finally:
            os.close(dir_fd)

    return HybridNetworkResult(
        dest_path=dest_path,
        qm_ts_labels=tuple(sorted(qm_transition_states)),
        ilt_ts_labels=tuple(sorted(known_labels - set(qm_transition_states))),
        vendored_files=vendored_files,
        warnings=tuple(warnings_list),
        t_grid_clamp=t_grid_clamp,
    )


def _read_qm_artifact(ts_label: str, artifact_path: str, allowed_log_root: str) -> dict:
    """
    Parse a QM statmech artifact file and resolve every ``Log(...)`` path it references.

    Every resolved ``Log(...)`` path is CONFINED to ``allowed_log_root``: an artifact is
    caller-supplied data, and its ``Log(...)`` arguments name files this module's vendoring step
    will later COPY into the project -- an unconfined path (absolute, or a relative one
    traversing out via ``..``/a symlink) could therefore exfiltrate any file the process can
    read (``/etc/passwd``, an SSH key) into the run directory. The check mirrors the
    ``os.path.realpath`` + ``os.path.commonpath`` confinement in
    ``t3.runners.rmg_runner.run_arkane_job`` (and ``t3.pdep.discovery._confine_to_project``): a
    path that escapes the root RAISES; it is never silently skipped.

    Args:
        ts_label (str): The transition state label the artifact belongs to (used in error text).
        artifact_path (str): The path to the ARC-written statmech artifact file.
        allowed_log_root (str): The directory every resolved ``Log(...)`` path must live under --
            the artifact's own ARC project directory. Always supplied explicitly by the caller,
            never inferred here.

    Raises:
        ValueError: If the artifact cannot be parsed as Python; references a ``Log(...)`` whose
                   argument is not a literal string (see below -- refused, never skipped);
                   references a ``Log(...)`` file that resolves outside ``allowed_log_root``; or
                   references a ``Log(...)`` file that does not exist.

    Returns:
        dict: ``{'content': str, 'sha256': str, 'logs': [(original_quoted_path, resolved_abs_path),
        ...], 'log_occurrences': [{'original_path': str, 'lineno': int, 'col_offset': int,
        'end_lineno': int, 'end_col_offset': int}, ...]}``. ``'sha256'`` is the hex digest of the
        EXACT bytes this call read from ``artifact_path`` -- computed from the same read used to
        parse ``'content'``, never from a second, separate read of the file afterward. A caller
        that instead re-opened ``artifact_path`` later to hash it (as this module's own
        ``t3.pdep.capture`` used to) would be exposed to a TOCTOU: if the pointer file changed on
        disk between the two reads, the recorded hash would describe bytes that were never the ones
        parsed or vendored, silently making the manifest lie about provenance. ``'logs'`` is
        de-duplicated by ``original_quoted_path`` (the same log file is often referenced by more
        than one keyword, e.g. ``geometry`` and ``frequencies``) and is used to decide which files
        to copy. ``'log_occurrences'`` is NOT de-duplicated -- it carries one entry (with the exact
        source position of that ``Log(...)`` call's string-literal argument) per occurrence, so
        every reference gets rewritten even when two keywords point at the same original path.
    """
    with open(artifact_path, 'rb') as f:
        raw_content = f.read()
    content_sha256 = hashlib.sha256(raw_content).hexdigest()
    content = raw_content.decode('utf-8')
    tree = _parse_as_ast(text=content, path=artifact_path)
    artifact_dir = os.path.dirname(os.path.abspath(artifact_path))
    resolved_root = os.path.realpath(allowed_log_root)

    logs = list()
    log_occurrences = list()
    seen_original_paths = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or _get_call_name(node) != 'Log':
            continue
        if not node.args or not isinstance(node.args[0], ast.Constant) or not isinstance(node.args[0].value, str):
            # Refuse rather than degrade: a Log(...) whose argument is not a literal string
            # (``Log(some_var)``, ``Log(os.path.join(...))``) cannot be checked for existence or
            # confinement, so skipping it would silently break this function's guarantee that
            # every referenced log exists -- the same fail-open bug fixed for
            # ``products=``/``isomers=``/``reactants=`` in ``t3/pdep/parser.py``.
            unevaluable = ast.get_source_segment(content, node.args[0]) if node.args else '<no argument>'
            raise ValueError(f"Transition state '{ts_label}''s QM artifact at '{artifact_path}' contains a "
                             f"Log(...) call (line {node.lineno}) whose argument is not a literal string and "
                             f"so cannot be checked for existence or confinement: Log({unevaluable}). "
                             f"Refusing to read it.")
        arg_node = node.args[0]
        original_path = arg_node.value
        log_occurrences.append({'original_path': original_path,
                                'lineno': arg_node.lineno,
                                'col_offset': arg_node.col_offset,
                                'end_lineno': arg_node.end_lineno,
                                'end_col_offset': arg_node.end_col_offset})
        if original_path in seen_original_paths:
            continue
        seen_original_paths.add(original_path)
        resolved_path = original_path if os.path.isabs(original_path) else os.path.join(artifact_dir, original_path)
        resolved_real_path = os.path.realpath(resolved_path)
        if resolved_real_path == resolved_root \
                or os.path.commonpath([resolved_root, resolved_real_path]) != resolved_root:
            raise ValueError(f"Transition state '{ts_label}''s QM artifact at '{artifact_path}' references a "
                             f"Log(...) path that resolves outside its allowed root: '{original_path}' resolves "
                             f"to '{resolved_real_path}', which is not under '{allowed_log_root}' (resolved to "
                             f"'{resolved_root}'). Refusing to read (or later copy) it.")
        if not os.path.isfile(resolved_path):
            raise ValueError(f"Transition state '{ts_label}''s QM artifact at '{artifact_path}' references a "
                             f"Log(...) file that does not exist: '{original_path}' (resolved to "
                             f"'{resolved_path}').")
        logs.append((original_path, resolved_path))

    return {'content': content, 'sha256': content_sha256, 'logs': logs, 'log_occurrences': log_occurrences}


def _vendored_log_names(resolved_paths: list) -> dict:
    """
    Determine collision-free file names under which to vendor a transition state's referenced logs.

    Flattening by basename alone is NOT safe: ARC deliberately renames every job's output file to
    ``output.out`` "for consistency between software" (``arc/job/local.py:356``,
    ``arc/job/adapter.py:870``), so a single transition state's sp, freq and rotor-scan logs are all
    named ``output.out``, differing only in their job directory. Copying them by basename would
    silently collapse them onto one file (last write wins) and rewrite every ``Log(...)`` reference
    to point at it, so the frequencies could be read from the sp job's output. That failure is
    invisible: Arkane parses the result happily and returns plausible but wrong kinetics.

    Names are therefore disambiguated by prepending the job directory, but ONLY for the basenames
    that actually collide, so the common case keeps readable file names. The decision is made over
    the whole set rather than in encounter order, so it does not depend on the order in which the
    ``Log(...)`` calls happen to appear.

    Args:
        resolved_paths (list): The absolute paths of the log files referenced by one transition
                              state's QM artifact.

    Returns:
        dict: Absolute log path -> the file name to vendor it under.
    """
    basename_counts = dict()
    for resolved_path in resolved_paths:
        basename = os.path.basename(resolved_path)
        basename_counts[basename] = basename_counts.get(basename, 0) + 1

    names, taken_names = dict(), set()
    for resolved_path in resolved_paths:
        basename = os.path.basename(resolved_path)
        if basename_counts[basename] > 1:
            job_dir = os.path.basename(os.path.dirname(resolved_path))
            name = f'{job_dir}_{basename}' if job_dir else basename
        else:
            name = basename
        if name in taken_names:
            # Pathological: two colliding basenames whose job directories also match.
            stem, extension = os.path.splitext(name)
            index = 1
            while f'{stem}_{index}{extension}' in taken_names:
                index += 1
            name = f'{stem}_{index}{extension}'
        taken_names.add(name)
        names[resolved_path] = name
    return names


def _vendor_qm_artifacts(artifact_infos: dict, qm_dir: str, dest_dir: str) -> tuple:
    """
    Copy every QM artifact (and every log file it references) into ``qm_dir``, rewriting the
    artifacts' ``Log(...)`` paths to be relative to their new location.

    ARC recreates its statmech directory with ``delete_existing_subdir=True``, so a log file
    referenced in place could be deleted out from under a later solve; vendoring makes the hybrid
    network self-contained.

    Each ``Log(...)`` path is rewritten by splicing at that specific occurrence's own string
    ``ast.Constant`` node position (``_read_qm_artifact``'s ``log_occurrences``), not by an
    exact-string text search-and-replace over the whole artifact. A text-based rewrite would (a)
    rewrite every textual occurrence of the original path, including ones outside a ``Log(...)``
    call (e.g. a mention in a comment), and (b) miss a path built via implicit string concatenation
    (``Log('a' 'b')``), since that never appears as one contiguous quoted literal in the source.
    Splicing by node position is exact regardless of either case, and (mirroring
    ``write_hybrid_network_input_file``'s own edit-application pattern) is applied back-to-front,
    by descending start offset, so earlier positions -- all computed against the pristine,
    pre-edit content -- stay valid as edits are applied.

    Re-solving into the same output directory must not leave orphaned files behind from a prior
    solve's ``qm_dir`` -- e.g. a transition state that was QM'd previously but is left as ILT (or
    renamed/removed from ``qm_transition_states``) this time around would otherwise leave its
    stale ``qm/<label>.py`` and ``qm/logs/<label>/`` lingering, undetected, alongside the fresh
    output. So if ``qm_dir`` already exists, it is replaced by the freshly staged one. This
    replacement is guarded tightly -- only when ``qm_dir``'s basename is literally ``'qm'`` and it
    is a direct child of ``dest_dir`` (i.e. exactly the path this function itself constructs and
    populates) -- so a caller-supplied ``qm_dir`` that does not match that exact shape is never
    touched (it is merged into instead; see below).

    Every artifact and log file is staged into a sibling temporary directory FIRST, and ``qm_dir``
    is only replaced once every single one of them has been copied/written successfully. Without
    this, a mid-copy failure (e.g. a disk error partway through copying the third of five log
    files) would have already destroyed the previous, good ``qm_dir`` before discovering the
    failure, permanently destroying a prior solve's vendored artifacts and leaving ``qm_dir`` in a
    half-written state -- even though the overall write is correctly aborted (no ``dest_path`` is
    ever written). Staging first means any exception raised while populating the stage leaves the
    previous ``qm_dir`` (if any) completely untouched.

    The replacement of an existing, managed ``qm_dir`` is itself done as an atomic SWAP, not a
    destroy-then-replace: the existing ``qm_dir`` is first renamed aside to a sibling temp path
    (a same-filesystem rename, not a copy), then the staged directory is ``os.replace``'d into the
    now-vacant ``qm_dir`` path (atomic, since the target is absent), and only then is the
    renamed-aside old tree removed. This closes the gap a naive ``rmtree`` then ``replace`` leaves
    open: a crash, kill, or exception between the two steps would otherwise destroy the previous
    capture and leave nothing in its place, even though ``os.replace`` itself is atomic. If
    ``os.replace`` fails, the old tree is renamed straight back so ``qm_dir`` is left exactly as it
    was; ``qm_dir`` is never observed empty or partially populated at any point in this sequence.

    Args:
        artifact_infos (dict): Transition state label -> the dict returned by ``_read_qm_artifact``.
        qm_dir (str): The ``qm/`` directory to populate. Atomically replaced if it already exists
                     (see above), only once staging succeeds in full.
        dest_dir (str): The hybrid network input file's directory, used to render
                       ``vendored_files`` as paths relative to it, and to guard the removal above.

    Returns:
        tuple: The vendored file paths, relative to ``dest_dir``.
    """
    is_managed_qm_dir = (os.path.basename(os.path.normpath(qm_dir)) == 'qm'
                         and os.path.dirname(os.path.normpath(qm_dir)) == os.path.normpath(dest_dir))

    staging_dir = tempfile.mkdtemp(prefix='.qm-staging-', dir=dest_dir)
    try:
        vendored_files = list()
        for ts_label, info in artifact_infos.items():
            content = info['content']
            log_dir = os.path.join(staging_dir, 'logs', ts_label)
            log_names = _vendored_log_names(resolved_paths=[resolved for _, resolved in info['logs']])

            new_relative_paths = dict()  # original_path -> the path to rewrite it to
            for original_path, resolved_path in info['logs']:
                if not os.path.isdir(log_dir):
                    os.makedirs(log_dir)
                basename = log_names[resolved_path]
                dest_log_path = os.path.join(log_dir, basename)
                shutil.copyfile(resolved_path, dest_log_path)
                vendored_files.append(os.path.relpath(dest_log_path, staging_dir))
                new_relative_paths[original_path] = os.path.join('logs', ts_label, basename)

            line_starts = _line_start_offsets(content)
            edits = list()
            for occurrence in info['log_occurrences']:
                start = _pos(content, line_starts, occurrence['lineno'], occurrence['col_offset'])
                end = _pos(content, line_starts, occurrence['end_lineno'], occurrence['end_col_offset'])
                prefix, quote = _string_literal_prefix_and_quote(content[start:end], ts_label,
                                                                  occurrence['original_path'])
                new_relative_path = new_relative_paths[occurrence['original_path']]
                edits.append((start, end, f'{prefix}{quote}{new_relative_path}{quote}'))
            for start, end, replacement in sorted(edits, key=lambda edit: edit[0], reverse=True):
                content = content[:start] + replacement + content[end:]

            dest_ts_path = os.path.join(staging_dir, f'{ts_label}.py')
            with open(dest_ts_path, 'w') as f:
                f.write(content)
            vendored_files.append(os.path.relpath(dest_ts_path, staging_dir))

        # Every artifact and log file for this solve has now been staged successfully, outside
        # qm_dir. Only now -- with nothing left in this vendoring pass that can fail -- replace
        # qm_dir's previous contents (if any) with what was staged, so an exception anywhere above
        # leaves a prior solve's qm_dir completely untouched rather than half-cleared.
        if os.path.isdir(qm_dir) and not is_managed_qm_dir:
            # Never remove a directory this function does not itself own (see docstring); merge
            # the staged files into it in place instead, preserving this defensive branch's
            # original never-touch-a-foreign-qm_dir behavior.
            for root, _dirs, files in os.walk(staging_dir):
                rel_root = os.path.relpath(root, staging_dir)
                dest_root = qm_dir if rel_root == '.' else os.path.join(qm_dir, rel_root)
                os.makedirs(dest_root, exist_ok=True)
                for name in files:
                    shutil.copyfile(os.path.join(root, name), os.path.join(dest_root, name))
            shutil.rmtree(staging_dir)
        else:
            # Never destroy the previous qm_dir before the new one is safely in place: a
            # destroy-then-replace sequence has a gap (crash, kill, or an exception from
            # os.replace itself) in which qm_dir is simply gone -- with nothing in it, while an
            # older manifest elsewhere may still claim the capture is valid. Instead, rename the
            # existing qm_dir aside first (a same-filesystem, near-instant rename, not a copy),
            # then os.replace the staged directory into the now-vacant qm_dir path (atomic, since
            # the target is absent), and only rmtree the renamed-aside old tree once the new one
            # is confirmed in place. If os.replace itself fails, the old tree is renamed straight
            # back so qm_dir is left exactly as it was.
            old_qm_dir = None
            if os.path.isdir(qm_dir):
                old_qm_dir = tempfile.mkdtemp(prefix='.qm-old-', dir=dest_dir)
                os.rmdir(old_qm_dir)
                os.rename(qm_dir, old_qm_dir)
            try:
                os.replace(staging_dir, qm_dir)
            except BaseException:
                if old_qm_dir is not None:
                    os.rename(old_qm_dir, qm_dir)
                raise
            if old_qm_dir is not None:
                shutil.rmtree(old_qm_dir)
    finally:
        if os.path.isdir(staging_dir):
            shutil.rmtree(staging_dir, ignore_errors=True)

    qm_relpath = os.path.relpath(qm_dir, dest_dir)
    return tuple(os.path.normpath(os.path.join(qm_relpath, vendored_file)) for vendored_file in vendored_files)


def _build_energy_header(energy_settings: QMEnergySettings) -> str:
    """
    Render the energy-reference header injected ahead of the first ``species(...)`` block.

    Never emits ``useAtomCorrections = False``: that directive silently disables Arkane's atom
    energy corrections rather than loudly failing without them, which is exactly the dangerous
    case this module exists to avoid (``write_hybrid_network_input_file`` already raises before
    this function is ever reached if ``energy_settings.use_atom_corrections`` is ``False``, so
    this function always renders ``useAtomCorrections = True``).

    ``modelChemistry`` is rendered bare (unquoted) when ``energy_settings.model_chemistry`` is a
    ``LevelOfTheory(...)``/``CompositeLevelOfTheory(...)`` object expression (it must stay
    executable Arkane syntax, not become a quoted string containing that syntax as text), and
    quoted otherwise (the plain string label case).

    Args:
        energy_settings (QMEnergySettings): The settings to render.

    Returns:
        str: The header text, terminated with a blank line.
    """
    if _model_chemistry_ast_call(energy_settings.model_chemistry) is not None:
        lines = [f'modelChemistry = {energy_settings.model_chemistry}']
    else:
        lines = [f'modelChemistry = "{energy_settings.model_chemistry}"']
    if energy_settings.frequency_scale_factor is not None:
        lines.append(f'frequencyScaleFactor = {energy_settings.frequency_scale_factor}')
    lines.append(f'useHinderedRotors = {energy_settings.use_hindered_rotors}')
    lines.append(f'useBondCorrections = {energy_settings.use_bond_corrections}')
    if energy_settings.use_bond_corrections:
        lines.append(f'bondCorrectionType = {energy_settings.bond_correction_type!r}')
    if energy_settings.atom_energies is not None:
        lines.append(f'atomEnergies = {energy_settings.atom_energies!r}')
    lines.append(f'useAtomCorrections = {energy_settings.use_atom_corrections}')
    lines.append('')
    return '\n'.join(lines) + '\n'


def _kinetics_removal_edit(text: str, line_starts: list, kinetics_kw: ast.keyword, call: ast.Call,
                            tunneling: str | None) -> tuple:
    """
    Build the edit that drops a ``reaction(...)`` block's ``kinetics = ...`` entry, optionally
    replacing it with ``tunneling = '...'``.

    The removal span's end boundary is located structurally from the AST -- the next keyword's own
    position if 'kinetics' is not the last keyword, or the enclosing call's closing ')' position if
    it is -- rather than by pattern-matching for a trailing comma in the source text. A regex-based
    trailing-comma search cannot tell 'kinetics' apart from other keywords' text, requires an
    immediate newline right after the (optional) comma to fire, and is defeated by a trailing
    inline '# comment': any of these produces a dangling comma (or a double comma, if 'kinetics' is
    not last) that is invalid Python. Using the neighboring node's own boundary instead is correct
    regardless of where 'kinetics' sits among the call's keywords, whether the call spans one line
    or many, and whether a comment trails the removed entry.

    Args:
        text (str): The full source text the offsets are computed against.
        line_starts (list): Per-line character offsets, as returned by ``_line_start_offsets``.
        kinetics_kw (ast.keyword): The ``kinetics=`` keyword node to remove.
        call (ast.Call): The enclosing ``reaction(...)`` call node, used to locate the keyword
                         immediately following ``kinetics_kw`` (or the call's own closing paren, if
                         ``kinetics_kw`` is the last keyword).
        tunneling (str, optional): The tunneling model name to write in its place, or ``None`` to
                                   drop ``kinetics = ...`` without replacing it.

    Returns:
        tuple: ``(start_offset, end_offset, replacement_text)``.
    """
    start = _pos(text, line_starts, kinetics_kw.lineno, kinetics_kw.col_offset)

    kw_index = call.keywords.index(kinetics_kw)
    if kw_index + 1 < len(call.keywords):
        boundary_node = call.keywords[kw_index + 1]
        end = _pos(text, line_starts, boundary_node.lineno, boundary_node.col_offset)
    else:
        # 'kinetics' is the last keyword: the boundary is the call's own closing ')'.
        end = _pos(text, line_starts, call.end_lineno, call.end_col_offset - 1)

    # Back the start up to the beginning of its line if only whitespace precedes 'kinetics' on
    # it, so the removed line's own indentation is not left behind as a blank/whitespace line.
    line_start = text.rfind('\n', 0, start) + 1
    indent = text[line_start:start]
    if indent.strip() == '':
        start = line_start
    else:
        indent = ''

    # Symmetrically, back the end boundary up to the start of its own line if only whitespace
    # (or nothing) precedes it there, so the shared ', ' separator glue and any trailing comment on
    # 'kinetics'' own line(s) are swept away too, without eating the next keyword's (or the closing
    # paren's) own indentation.
    end_line_start = text.rfind('\n', 0, end) + 1
    if text[end_line_start:end].strip() == '':
        end = end_line_start

    replacement = f"{indent}tunneling = {tunneling!r},\n" if tunneling else ''
    return start, end, replacement


def _rewrite_method_and_sensitivity(text: str, method: str, sensitivity: bool,
                                    source_path: str) -> tuple[str, TGridClampRecord]:
    """
    Rewrite the ``method = '...'`` line and (optionally) inject a ``sensitivity_conditions``
    directive spanning the network's T/P extrema.

    Mirrors ``t3.utils.writer.write_arkane_network_input_file``'s line-scan approach (reusing
    ``rewrite_arkane_method_line`` for the method line itself, and the same quote-/whitespace-
    tolerant ``METHOD_LINE_CANDIDATE_RE`` guard) so the two writers stay consistent.

    Args:
        text (str): The network input file text (after TS/reaction/header edits).
        method (str): 'CSE', 'MSC' or 'RS'.
        sensitivity (bool): Whether to inject a ``sensitivity_conditions`` directive.
        source_path (str): The source network path, used only in the error message below.

    Raises:
        ValueError: If ``sensitivity`` is ``True`` and the T/P extrema could not be parsed; if
                   ``text`` does not contain exactly one ``method = ...`` line (zero means no
                   method line was found to rewrite; more than one is ambiguous); or if the
                   network's pdep ``Tmax`` exceeds the narrowest species thermo (NASA) ceiling in
                   the file AND clamping to that ceiling would put ``Tmax`` at or below ``Tmin``
                   (no valid temperature point left to solve at all). Each case is refused rather
                   than silently handing back a network solved with the wrong -- or an unrunnable
                   -- grid.

    Returns:
        tuple: The rewritten text (str), and the ``TGridClampRecord`` provenance for whether the
              T grid written was clamped down from what the source network requested. See
              ``TGridClampRecord``'s docstring for the three-state ("clamped" / "not clamped" /
              caller never asked -- N/A here, this function always computes and returns it) design.
    """
    # Computed once, from the SAME text this function rewrites (species(...) blocks are untouched
    # by the TS/reaction/header edits applied before this is called): this is the ceiling every
    # species in the network has valid NASA thermo up to, mirroring
    # t3.utils.writer.write_arkane_network_input_file's clamp so both writers hand Arkane a grid
    # it can actually solve rather than one that fails with "No valid NASA polynomial at
    # temperature ... K." See network_thermo_t_max's docstring for why None (no species
    # contributed a determinable ceiling) means "clamp nothing".
    #
    # If text does not even parse as Python, that is NOT this function's problem to report: it
    # means an earlier edit in write_hybrid_network_input_file's pipeline (TS/reaction/header
    # splicing) corrupted the text, and that pipeline already has its own dedicated self-check
    # (an ast.parse(text) right after this function returns) whose entire job is to catch exactly
    # that and raise a RuntimeError naming dest_path. Letting network_thermo_t_max's
    # NetworkTextUnparseable propagate from here would pre-empt that self-check with a less
    # specific error raised for the wrong reason (this function was only trying to compute a
    # clamp ceiling, not validate the whole file), so ONLY that specific, unparseable-text case
    # is caught here and treated the same as "no species contributed a determinable ceiling":
    # clamp nothing, and let the pipeline's own self-check be the one that catches and reports
    # the corruption. Any OTHER ValueError network_thermo_t_max might raise (e.g. a malformed
    # **kwargs unpacking in a species(...) call) is a real defect in the network text and must
    # propagate rather than being silently treated the same as "text is fine, just no ceiling".
    try:
        ceiling = network_thermo_t_max(text)
    except NetworkTextUnparseable:
        ceiling = None
    thermo_t_max = ceiling.t_max if ceiling is not None else None
    if ceiling is not None and ceiling.skipped:
        logger.warning(
            f"Network '{source_path}': the NASA thermo ceiling used to clamp the hybrid network "
            f"({thermo_t_max} K) could not account for {len(ceiling.skipped)} species whose "
            f"thermo could not be read: {format_skipped_species(ceiling.skipped)}. The true "
            f"network-wide ceiling may therefore be lower than {thermo_t_max} K.")
    lines = text.splitlines(keepends=True)
    new_lines = list()
    t_min, t_max, t_count, p_min, p_max = None, None, None, None, None
    parse_tp = False
    method_rewrite_count = 0
    # T-grid clamp provenance (mirrors t3.utils.writer.write_arkane_network_input_file's
    # identical tracking; see TGridClampRecord's docstring for why "clamped" is a three-state
    # design and why this is tracked at all): populated as the clamp decisions below are
    # actually made, not inferred after the fact from the final t_max/t_count.
    requested_t_max, clamped = None, False
    tlist_dropped, tlist_original_highest = False, None
    for line in lines:
        skip_line = False
        if 'pressureDependence(' in line:
            parse_tp = True
        if parse_tp:
            if 'Tmin' in line and '(' in line:
                t_min = line.split('(')[1].split(',')[0]
            elif 'Tmax' in line and '(' in line:
                t_max = line.split('(')[1].split(',')[0]
                requested_t_max = float(t_max)
                if thermo_t_max is not None and float(t_max) > thermo_t_max:
                    if t_min is None or thermo_t_max <= float(t_min):
                        raise ValueError(
                            f"Cannot write a hybrid network from '{source_path}': the pdep "
                            f"network's requested Tmax ({t_max} K) exceeds the narrowest species "
                            f"thermo ceiling in the network ({thermo_t_max} K), and clamping to "
                            f"that ceiling would leave a Tmax ({thermo_t_max} K) at or below Tmin "
                            f"({t_min} K), i.e. no valid temperature point to solve at all.")
                    logger.warning(
                        f"Network '{source_path}': requested pdep Tmax of {t_max} K exceeds the "
                        f"narrowest species thermo (NASA) ceiling in the file, {thermo_t_max} K; "
                        f"clamping the hybrid Arkane grid's Tmax to {thermo_t_max} K so standalone "
                        f"Arkane does not refuse the solve with 'No valid NASA polynomial at "
                        f"temperature ... K.'")
                    formatted_t_max = format_clamped_t_max(thermo_t_max)
                    paren_index = line.index('(')
                    comma_index = line.index(',', paren_index)
                    line = line[:paren_index + 1] + formatted_t_max + line[comma_index:]
                    t_max = formatted_t_max
                    clamped = True
            elif 'Tcount' in line and '=' in line:
                #     Tcount = 8,
                t_count = line.split('=', 1)[1].split(',')[0].strip()
            elif 'Tlist' in line and '(' in line:
                #     Tlist = ([3200,2290.91,...,700],'K'),
                # RMG's writer always precomputes an explicit Tlist alongside Tmin/Tmax/Tcount, and
                # Arkane's own pdep.py only calls generate_T_list() (which would build the grid from
                # Tmin/Tmax/Tcount) when self.Tlist is None. If Tlist is present, it wins outright
                # and Tmin/Tmax are ignored -- so clamping the Tmax line above is a no-op unless this
                # stale, out-of-range Tlist is also dropped so Arkane regenerates the grid over the
                # clamped range instead of solving at the network's original (too-high) grid points.
                if thermo_t_max is not None:
                    rhs = line.split('=', 1)[1].strip()
                    if rhs.endswith(','):
                        rhs = rhs[:-1]
                    try:
                        parsed = ast.literal_eval(rhs)
                    except (ValueError, SyntaxError) as e:
                        raise ValueError(
                            f"Cannot write a hybrid network from '{source_path}': a clamp is in "
                            f"play (the narrowest species thermo ceiling is {thermo_t_max} K), "
                            f"so the pdep network's Tlist line must be understood to decide "
                            f"whether to drop it, but it could not be parsed as a Python "
                            f"literal: {e}")
                    if (not isinstance(parsed, tuple) or len(parsed) != 2
                            or not isinstance(parsed[0], (list, tuple))
                            or not all(isinstance(v, (int, float)) and not isinstance(v, bool)
                                       for v in parsed[0])
                            or parsed[1] != 'K'):
                        raise ValueError(
                            f"Cannot write a hybrid network from '{source_path}': a clamp is in "
                            f"play (the narrowest species thermo ceiling is {thermo_t_max} K), "
                            f"so the pdep network's Tlist line must be understood to decide "
                            f"whether to drop it, but it is not a 2-element (sequence-of-numbers, "
                            f"unit) tuple with unit 'K'; got: {parsed!r}.")
                    entries = parsed[0]
                    if any(entry > thermo_t_max for entry in entries):
                        missing = [name for name, val in (('Tmin', t_min), ('Tmax', t_max),
                                                          ('Tcount', t_count)) if val is None]
                        if missing:
                            raise ValueError(
                                f"Cannot write a hybrid network from '{source_path}': the pdep "
                                f"network's Tlist contains an entry above the narrowest species "
                                f"thermo ceiling ({thermo_t_max} K) and must be dropped so Arkane "
                                f"regenerates the T grid itself, but {', '.join(missing)} could not "
                                f"be parsed from the pressureDependence(...) block to regenerate "
                                f"from.")
                        highest = max(entries)
                        logger.warning(
                            f"Network '{source_path}': dropping the pressureDependence(...) block's "
                            f"explicit Tlist line when writing the hybrid network -- its highest "
                            f"entry ({highest} K) exceeds the narrowest species thermo (NASA) "
                            f"ceiling in the file, {thermo_t_max} K; Arkane will regenerate the T "
                            f"grid itself from the clamped Tmin ({t_min} K) / Tmax ({t_max} K) / "
                            f"Tcount ({t_count}) instead of using this network's original (too-high) "
                            f"grid.")
                        skip_line = True
                        tlist_dropped = True
                        tlist_original_highest = float(highest)
            elif 'Pmin' in line and '(' in line:
                p_min = line.split('(')[1].split(',')[0]
            elif 'Pmax' in line and '(' in line:
                p_max = line.split('(')[1].split(',')[0]
        if METHOD_LINE_CANDIDATE_RE.match(line):
            new_lines.append(rewrite_arkane_method_line(line=line, method=method))
            method_rewrite_count += 1
        elif 'rmgmode' in line:
            new_lines.append(line)
            if sensitivity:
                if any(param is None for param in [t_min, t_max, p_min, p_max]):
                    raise ValueError(f'Could not parse all T/P parameters, got:\n'
                                     f'T min = {t_min}, T max = {t_max}, P min = {p_min}, P max = {p_max}.')
                new_lines.append(f"""    sensitivity_conditions = [[({t_min}, 'K'), ({p_min}, 'bar')],
                              [({t_max}, 'K'), ({p_min}, 'bar')],
                              [({t_min}, 'K'), ({p_max}, 'bar')],
                              [({t_max}, 'K'), ({p_max}, 'bar')]],""")
        elif not skip_line:
            new_lines.append(line)

    if method_rewrite_count != 1:
        raise ValueError(f"Expected exactly one 'method = ...' line to rewrite in '{source_path}', found "
                         f"{method_rewrite_count}. Refusing to write a hybrid network with an unresolved (zero) "
                         f"or ambiguous (more than one) method line, rather than silently solving with the "
                         f"source file's original method.")

    t_grid_clamp = TGridClampRecord(
        clamped=clamped,
        requested_t_max=requested_t_max,
        thermo_ceiling=thermo_t_max,
        written_t_max=float(t_max) if t_max is not None else None,
        tlist_dropped=tlist_dropped,
        tlist_original_highest=tlist_original_highest,
        skipped_species=tuple(ceiling.skipped) if ceiling is not None else tuple(),
    )
    return ''.join(new_lines), t_grid_clamp


def _parse_as_ast(text: str, path: str) -> ast.Module:
    """
    Parse Python text into an AST, never ``exec``ing or ``import``ing it.

    Args:
        text (str): The Python source text.
        path (str): The path it was read from (used only in the error message).

    Raises:
        ValueError: If the text cannot be parsed as Python.

    Returns:
        ast.Module: The parsed tree.
    """
    try:
        # RMG/ARC-generated files sometimes contain comment text with invalid escape sequences
        # (e.g. ``\H``) that trigger a spurious SyntaxWarning under ``ast.parse``; see the
        # matching note in ``t3.pdep.parser``.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', SyntaxWarning)
            return ast.parse(text)
    except SyntaxError as e:
        raise ValueError(f"Could not parse '{path}' as Python: {e}")


def _line_start_offsets(text: str) -> list:
    """
    Compute the character offset each line starts at, for converting an AST node's
    (1-indexed lineno, col_offset) into a flat character index into ``text``.

    Args:
        text (str): The source text.

    Returns:
        list: ``line_starts[i]`` is the character offset line ``i + 1`` starts at.
    """
    line_starts = list()
    offset = 0
    for line in text.splitlines(keepends=True):
        line_starts.append(offset)
        offset += len(line)
    return line_starts


def _pos(text: str, line_starts: list, lineno: int, col_offset: int) -> int:
    """
    Convert an AST node's (1-indexed lineno, col_offset) into a flat character index.

    Per the ``ast`` module's documented behavior, ``col_offset`` is a UTF-8 BYTE offset into its
    line, not a character index -- applying it directly as a character index (i.e.
    ``line_starts[lineno - 1] + col_offset``) silently misaligns every position on a line that has
    any non-ASCII (multi-byte-in-UTF-8) character before it. To convert correctly, this re-encodes
    the target line to UTF-8, slices it to the byte offset, and decodes back to count how many
    characters that many bytes actually span.

    Args:
        text (str): The full source text ``line_starts`` was computed against.
        line_starts (list): Per-line character offsets, as returned by ``_line_start_offsets``.
        lineno (int): The 1-indexed line number.
        col_offset (int): The 0-indexed UTF-8 byte offset within that line.

    Returns:
        int: The character offset into the source text.
    """
    line_start = line_starts[lineno - 1]
    line_end = line_starts[lineno] if lineno < len(line_starts) else len(text)
    line_text = text[line_start:line_end]
    char_offset = len(line_text.encode('utf-8')[:col_offset].decode('utf-8'))
    return line_start + char_offset


# The only string-literal prefixes that can produce an ``ast.Constant`` whose ``.value`` is a
# ``str`` (as opposed to ``bytes``, or an ``ast.JoinedStr`` for f-strings, both of which
# ``_read_qm_artifact`` already refuses before a log occurrence is ever recorded): no prefix, or a
# bare ``r``/``R`` or ``u``/``U``. ``b``/``rb``/``br`` (any case) parse to a ``bytes`` value, and
# any ``f``/``rf``/``fr`` prefix parses to ``ast.JoinedStr`` rather than ``ast.Constant`` -- neither
# ever reaches this function. This was verified directly against the ``ast`` module rather than
# assumed.
_STRING_LITERAL_PREFIX_RE = re.compile(r'^[rRuU]?')


def _string_literal_prefix_and_quote(source_segment: str, ts_label: str, original_path: str) -> tuple:
    """
    Determine the prefix and quote characters that bound a string-literal source segment (the
    exact source text of an ``ast.Constant`` string node), so a ``Log(...)`` occurrence can be
    rewritten without assuming the literal's first character is always the opening quote.

    That assumption is false for a prefixed string literal (e.g. ``Log(r'foo/bar.log')``): its
    source segment starts with the prefix character(s), not a quote, and splicing at the wrong
    offset produces invalid Python text (the closing quote survives with an extra prefix character
    stuck in front of it). This computes the actual bounding prefix and quote (single or triple,
    single- or double-quote style) instead, and refuses -- rather than guessing -- if the segment's
    shape cannot be unambiguously determined, since a wrong guess would silently corrupt the
    generated file.

    Args:
        source_segment (str): The exact source text of the string-literal node (``content[start:end]``).
        ts_label (str): The transition state the occurrence belongs to (used in the error message).
        original_path (str): The occurrence's ``original_path`` (used in the error message).

    Raises:
        ValueError: If the segment's bounding prefix/quote cannot be determined.

    Returns:
        tuple: ``(prefix, quote)``, e.g. ``('', "'")`` or ``('r', '"')`` or ``('', '\"\"\"')``.
    """
    prefix_match = _STRING_LITERAL_PREFIX_RE.match(source_segment)
    prefix = prefix_match.group(0)
    rest = source_segment[len(prefix):]
    for quote in ('"""', "'''", '"', "'"):
        if (rest.startswith(quote) and source_segment.endswith(quote)
                and len(source_segment) >= len(prefix) + 2 * len(quote)):
            return prefix, quote
    raise ValueError(f"Transition state '{ts_label}''s QM artifact has a Log(...) argument "
                     f"'{original_path}' whose source text {source_segment!r} does not have a "
                     f"recognizable string-literal prefix/quote shape. Refusing to rewrite it "
                     f"rather than risk splicing at the wrong offset and producing invalid Python.")


def _consume_trailing_newline(text: str, index: int) -> int:
    """
    Advance an offset past a single trailing newline (``\\n`` or ``\\r\\n``), if present.

    Args:
        text (str): The text ``index`` indexes into.
        index (int): The offset to advance.

    Returns:
        int: ``index``, advanced past a trailing newline if one was there.
    """
    if text[index:index + 2] == '\r\n':
        return index + 2
    if text[index:index + 1] == '\n':
        return index + 1
    return index


def _get_call_name(call: ast.Call) -> str | None:
    """
    Get the callee name of an ``ast.Call`` node (e.g., ``'Log'`` for ``Log(...)``).

    Only a bare ``ast.Name`` callee yields a name, for the same reason as
    ``t3.pdep.parser._get_call_name``: an attribute call such as ``foo.reaction(...)`` used to report
    ``'reaction'``, so every caller dispatching on this name treated it as the Arkane/RMG DSL
    directive it merely resembles. Arkane's loader binds these names in a namespace with no builtins
    and no imports, so a directive can only ever BE a bare name; anything with a dot in front of it is
    not the directive it is named after.

    Args:
        call (ast.Call): The call node.

    Returns:
        str: The callee name, or ``None`` if the callee is not a bare name.
    """
    if isinstance(call.func, ast.Name):
        return call.func.id
    return None


def _literal_or_none(node):
    """
    Safely evaluate an AST node to a Python literal, if possible.

    Args:
        node: An AST node (or ``None``).

    Returns:
        The evaluated literal, or ``None`` if ``node`` is ``None`` or is not a literal.
    """
    if node is None:
        return None
    try:
        return ast.literal_eval(node)
    except (ValueError, TypeError):
        return None


def _is_literal_none(node) -> bool:
    """
    Check whether an AST node is literally the ``None`` constant (as opposed to a node that simply
    failed to evaluate as a literal, for which ``_literal_or_none`` also returns ``None``).

    Args:
        node: An AST node (or ``None``).

    Returns:
        bool: True if ``node`` is an ``ast.Constant`` whose value is ``None``.
    """
    return isinstance(node, ast.Constant) and node.value is None
