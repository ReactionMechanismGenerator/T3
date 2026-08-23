"""
t3 pdep explorer input_file module.

Writes a valid Arkane PES-explorer input file from an RMG P-dep network file (or a T3 hybrid
Arkane network input; both are plain Arkane-DSL Python source), per
``docs/t3-pdep-qm-explorer-design.md`` section 4.

An RMG-written network file is NOT a valid Arkane explorer input on its own: it carries a
``network(...)`` block Arkane's explorer input parser does not recognize, has no ``database(...)``
directive (required for the explorer to look up thermo/kinetics for newly discovered species), and
every path reaction typically carries an explicit ``kinetics = ...`` entry that would otherwise
shadow the very kinetics Arkane's explorer is meant to estimate for it. This module rewrites the
file structurally (via ``ast.parse`` splicing, mirroring ``t3.pdep.hybrid.write_hybrid_network_input_file``)
rather than by text-scanning for ``network(``, so a coincidental match inside a comment or string
can never be mistaken for the block to drop.
"""

import ast
import io
import keyword
import math
import os
import re
import tempfile
import tokenize
import warnings
from dataclasses import dataclass, field

from arc.molecule.molecule import Molecule

from t3.pdep.hashing import hash_bytes
from t3.pdep.hybrid import _MODEL_CHEMISTRY_CALL_NAMES, _validate_model_chemistry_expression
from t3.utils.writer import METHOD_LINE_CANDIDATE_RE, METHOD_MAP, rewrite_arkane_method_line

# The only top-level Arkane DSL calls this module needs to recognize while walking the source's
# AST. Anything else (e.g. a stray helper call) is left untouched, exactly as
# ``t3.pdep.hybrid._TOP_LEVEL_CALL_NAMES`` does for the calls it cares about.
_TOP_LEVEL_CALL_NAMES = ('species', 'transitionState', 'reaction', 'network', 'pressureDependence')

# A minimal, Arkane-runnable default ``database(...)`` directive, modeled directly on the real
# fixture at RMG-Py/examples/arkane/explorer/methoxy/input.py. Every key here can be overridden
# (or added to) via the ``database_kwargs`` argument.
DEFAULT_DATABASE_KWARGS = {
    'thermoLibraries': ['primaryThermoLibrary'],
    'reactionLibraries': [],
    'kineticsDepositories': ['training'],
    'kineticsFamilies': 'default',
    'kineticsEstimator': 'rate rules',
}

# The exact keyword order the real Arkane explorer fixture uses; keeping this order is cosmetic
# (Arkane's input parser reads keyword arguments, not positional order) but makes a generated file
# diff cleanly against a hand-written one.
_EXPLORER_KWARG_ORDER = ('source', 'explore_tol', 'energy_tol', 'flux_tol', 'bathGas',
                        'maximumRadicalElectrons')

# Maps this module's snake_case argument names to the Arkane explorer(...) keyword names Arkane
# itself expects (Arkane's own DSL mixes camelCase and snake_case across directives).
_EXPLORER_KWARG_NAMES = {
    'source': 'source',
    'explore_tol': 'explore_tol',
    'energy_tol': 'energy_tol',
    'flux_tol': 'flux_tol',
    'bathGas': 'bathGas',
    'maximumRadicalElectrons': 'maximumRadicalElectrons',
}

# Characters that would let a crafted reaction label break out of (or corrupt) the quoted string
# literal it is rendered into for a ``kinetics('<label>')`` job directive in the generated,
# executable Arkane input; mirrors ``t3.pdep.hybrid._INJECTION_CHARS``. Rendering the label with
# ``!r`` (below) already makes this belt-and-suspenders rather than the only defense, but a label
# is read verbatim out of the source network file, so it is untrusted input either way.
_KINETICS_LABEL_INJECTION_CHARS = ('\n', '\r', "'", '"', '\\')

# Sentinel distinguishing "the species() block carries no 'reactive' keyword at all" from every
# value the keyword could literally hold (True/False, or None for a non-literal expression, which
# ``_literal_or_none`` cannot tell apart from a literal ``None``). The distinction is load-bearing
# for the bath-gas rule: an ABSENT keyword on a declared bath-gas species is the shape every
# RMG-written network file has (``arkane/pdep.py:654`` ``save_input_file`` never emits ``reactive``
# for any species) and gets ``reactive = False`` injected into the generated file, while an
# explicit ``reactive=True`` is a contradiction that is refused (see ``_validate_bath_gas``).
_REACTIVE_ABSENT = object()

# The subset of Arkane's own input-file namespace (``local_context`` at RMG-Py
# ``arkane/input.py:607``; T3 must not import arkane, hence a transcription rather than a reference)
# that this module will splice verbatim out of an untrusted source and into a NEW file Arkane will
# ``exec``. It is deliberately NOT the whole namespace. Arkane's namespace is grouped by its own
# comments, and only one of those groups is a set of inert value constructors -- these. Everything
# Arkane labels ``# Jobs`` (``kinetics``, ``statmech``, ``thermo``, ``pressureDependence``,
# ``explorer``, ``bac``, ``ae``) APPENDS A JOB when called, and everything it labels ``# Functions``
# (``species``, ``transitionState``, ``reaction``, ``network``, ``database``) mutates the loader's
# species/reaction/database state. Splicing a call to any of those out of an untrusted source means
# letting the source add jobs to, or reconfigure, the explorer run this module is generating -- a
# top-level ``database(...)`` in the source lands AFTER the generated one and silently overrides it.
# The five ``# Functions`` names are still legitimate as whole top-level statements (that is what a
# network source IS), which is what ``_TOP_LEVEL_CALL_NAMES`` permits; they are refused everywhere
# else, i.e. nested inside another expression, where they would only ever be a smuggled side effect.
# ``range`` is excluded as well: it is not a value constructor, and an unbounded ``range(10**9)``
# nested in a spliced literal is a memory bomb Arkane would evaluate at load time.
#
# Measured against every real Arkane/RMG file in RMG-Py's examples (33 RMG-written network*.py and
# 28 hand-written arkane input.py): the nested calls those files actually contain are drawn ONLY from
# the set below -- zero job directives, zero ``# Functions`` names, zero ``range``, in either group.
# Narrowing therefore costs no legitimate network source. If a future Arkane grows a new value
# constructor, a source using it is refused with an actionable message naming this constant.
_ARKANE_VALUE_NAMES = frozenset({
    'TransportData',
    'SingleExponentialDown',
    'Arrhenius',
    'RateUncertainty',
    'IdealGasTranslation',
    'LinearRotor',
    'NonlinearRotor',
    'KRotor',
    'SphericalTopRotor',
    'HarmonicOscillator',
    'HinderedRotor',
    'FreeRotor',
    'ThermoData',
    'Wilhoit',
    'NASA',
    'NASAPolynomial',
    'SMILES',
    'adjacencyList',
    'InChI',
    'LevelOfTheory',
    'CompositeLevelOfTheory',
})

# The names Arkane's ``# Jobs`` and ``# Functions`` groups define, i.e. the callables whose whole
# purpose is a side effect on the loader's state. Kept as a named set purely so a source that calls
# one can be refused with a message that says WHY that specific name is refused, rather than the
# generic "not a recognized Arkane DSL function" -- the distinction between "Arkane has never heard
# of this name" and "Arkane defines this name and that is exactly the problem" is the difference
# between an actionable refusal and a confusing one.
_ARKANE_SIDE_EFFECT_NAMES = frozenset({
    'kinetics',
    'statmech',
    'thermo',
    'pressureDependence',
    'explorer',
    'bac',
    'ae',
    'species',
    'transitionState',
    'reaction',
    'network',
    'database',
})


@dataclass(frozen=True)
class ExplorerInputSummary:
    """
    The outcome of ``write_arkane_explorer_input_file``.

    Args:
        dest_path (str): The path the Arkane explorer input file was written to.
        seed_species (tuple): The seed (source) species labels the explorer was configured with.
        bath_gas (dict): The bath gas composition recorded in the ``explorer(...)`` block.
        kinetics_labels_emitted (tuple): The reaction labels a ``kinetics('<label>')`` job
                                        directive was emitted for (reactions that had no explicit
                                        ``kinetics=`` keyword of their own).
        reactive_false_injected (tuple): The bath-gas species labels whose ``species()`` blocks had
                                        no ``reactive`` keyword in the source and had
                                        ``reactive = False`` injected into the generated file
                                        (issue #183: an RMG-written source never carries the
                                        keyword, so this records exactly what the writer added on
                                        the source's behalf).
        warnings (tuple): Human-readable warnings (e.g. that a multi-species bath gas's fractions
                          are recorded here but will not be honored by Arkane; see P16).
        spin_multiplicity_corrected (tuple): One ``(label, old, new)`` per species whose
                          ``spinMultiplicity`` contradicted its own adjacency-list ``multiplicity``
                          header and was rewritten to match it (see the correction block below).
                          RMG's explorer emits ``spinMultiplicity = 1`` for triplet atomic O
                          (adjacency ``multiplicity 3``, u2, thermo label ``O(T)``); re-feeding that
                          file makes Arkane round-trip an inconsistent species and crash, so any
                          RMG-explored network carrying O(3P) is otherwise un-re-feedable.
    """
    dest_path: str
    seed_species: tuple
    bath_gas: dict
    kinetics_labels_emitted: tuple = field(default_factory=tuple)
    reactive_false_injected: tuple = field(default_factory=tuple)
    warnings: tuple = field(default_factory=tuple)
    spin_multiplicity_corrected: tuple = field(default_factory=tuple)


def write_arkane_explorer_input_file(source_path: str,
                                     dest_path: str,
                                     seed_species,
                                     method: str,
                                     bath_gas: dict,
                                     explore_tol: float = None,
                                     energy_tol: float = None,
                                     flux_tol: float = None,
                                     maximum_radical_electrons: int = None,
                                     database_kwargs: dict = None,
                                     expected_source_hash: str = None,
                                     ) -> ExplorerInputSummary:
    """
    Write a valid Arkane PES-explorer input file from an RMG P-dep network (or hybrid Arkane
    network) input file.

    Args:
        source_path (str): The path to the source RMG P-dep network (or hybrid Arkane network)
                           file to build from.
        dest_path (str): The path to write the resulting Arkane explorer input file to. The parent
                         directory is created if it does not already exist.
        seed_species (list | tuple): The seed (source) species labels to explore from. Must
                                     contain exactly 1 or 2 labels, per Arkane's own limitation
                                     (P2, P3): ``explorer(...)``'s ``source`` is resolved from
                                     ``species_dict`` only (never a transition state), and Arkane's
                                     explorer accepts a unimolecular or bimolecular source only.
        method (str): 'CSE', 'MSC' or 'RS' (see ``t3.utils.writer.METHOD_MAP``), used to rewrite
                     the kept ``pressureDependence(...)`` block's ``method = ...`` line.
        bath_gas (dict): The bath gas composition, mapping species labels to mole fractions. Every
                        label must be a ``species()`` block in the source. P16: Arkane/RMG identify
                        the bath gas as every unreactive core species, not by name, so the
                        generated file must carry ``reactive=False`` on each of these blocks --
                        and since no RMG-written source ever carries the keyword
                        (``arkane/pdep.py:654`` never emits it; issue #183), this writer INJECTS
                        ``reactive = False`` into a declared bath-gas species block that has no
                        ``reactive`` keyword, and refuses one carrying an explicit literal
                        ``reactive=True`` (a contradiction) or a non-literal value (unverifiable).
                        See ``_validate_bath_gas``.
        explore_tol (float, optional): The energy tolerance for exploring new isomers/reactions.
        energy_tol (float, optional): The energy tolerance for including a well/transition state in
                                      the output network.
        flux_tol (float, optional): The flux tolerance for including a well/transition state in the
                                    output network.
        maximum_radical_electrons (int, optional): The maximum number of radical electrons allowed
                                                   in an explored species.
        database_kwargs (dict, optional): Keyword arguments overriding/extending
                                          ``DEFAULT_DATABASE_KWARGS`` for the prepended
                                          ``database(...)`` directive.
        expected_source_hash (str, optional): If given, the content hash (``t3.pdep.hashing``
                                              format) that ``source_path``'s bytes must match at
                                              the moment they are read here. A caller upstream may
                                              have already checked the file's hash before deciding
                                              to explore it, but that check and this read are two
                                              separate opens of the same path; a rewrite or symlink
                                              swap in between would let the earlier check pass on
                                              the old bytes while this function builds the explorer
                                              input from different ones. Passing the hash the
                                              caller validated against closes that window, because
                                              it is checked against the exact bytes this function
                                              reads and uses -- not a second, independent read.

    Raises:
        ValueError: If the source cannot be parsed; if ``seed_species`` does not contain exactly 1
                   or 2 labels; if a seed label is not a ``species()`` block in the source (e.g. it
                   is a ``transitionState()`` label, or does not exist at all, or is defined as
                   BOTH a ``species()`` and a ``transitionState()``); if a ``bath_gas`` label is
                   not a ``species()`` block, or is a ``species()`` block carrying an explicit
                   literal ``reactive=True`` or a non-literal ``reactive`` value (an absent
                   keyword is NOT refused -- ``reactive = False`` is injected into the generated
                   file instead); if the source does not contain exactly one
                   ``pressureDependence(...)`` block; if a ``reaction()`` block has neither an
                   explicit ``kinetics=`` keyword nor a literal ``label`` to target a
                   ``kinetics('<label>')`` job directive at; or if ``expected_source_hash`` is
                   given and does not match the content hash of the bytes actually read from
                   ``source_path``.
        RuntimeError: If the generated input file text fails its own ``ast.parse(...)``
                     self-check. This should never happen in practice (every edit is computed
                     structurally from the AST), but a failure here is caught loudly, before
                     anything is written to ``dest_path``, rather than silently handing back a
                     broken explorer input.

    Returns:
        ExplorerInputSummary: The outcome of the write.
    """
    # Same trap as in ``t3.pdep.explorer.adapter.validate_explorer_seed``, and checked here too
    # because this writer is reachable without going through an adapter at all. ``tuple('OH')`` is
    # ``('O', 'H')``, which passes the 1-or-2 rule below as a bimolecular seed.
    if isinstance(seed_species, (str, bytes)):
        raise ValueError(f"'seed_species' must be a sequence of label strings, not a bare "
                         f"{type(seed_species).__name__} ({seed_species!r}); a string is itself a sequence, "
                         f"so it would be read character by character as {tuple(seed_species)!r}.")
    seed_species = tuple(seed_species)
    for label in seed_species:
        if not isinstance(label, str) or not label:
            raise ValueError(f"Every entry of 'seed_species' must be a non-empty label string, got {label!r} "
                             f"of type {type(label).__name__} in {seed_species!r}.")
    if len(seed_species) not in (1, 2):
        raise ValueError(f"seed_species must contain exactly 1 or 2 labels (Arkane's explorer 'source' accepts a "
                         f"unimolecular or bimolecular seed only, see P3), got {len(seed_species)}: {seed_species}.")

    # Read the bytes ONCE, exactly as ``t3.pdep.parser.parse_pdep_network_file`` does, rather than
    # opening in text mode here and hashing a separate read elsewhere. Two reads would mean this
    # function could build its explorer input from content that never got hash-checked at all: if
    # the file changed between an earlier check and this open, the input written here would explore
    # bytes nobody approved.
    with open(source_path, 'rb') as f:
        data = f.read()
    if expected_source_hash is not None:
        actual_hash = hash_bytes(data)
        if actual_hash != expected_source_hash:
            raise ValueError(f"Source file '{source_path}' changed after it was validated: expected content hash "
                             f"{expected_source_hash!r} but the bytes read just now hash to {actual_hash!r}. "
                             f"Writing an explorer input from this file would explore bytes nobody approved.")
    # Decoded the way Python itself decodes source, via the PEP 263 encoding cookie, rather than with
    # a hard-coded 'utf-8' or the locale encoding open(path, 'r') would have used. See
    # ``t3.pdep.parser.parse_pdep_network_file`` for the full rationale; this mirrors it exactly so
    # the same file decodes identically in both places.
    encoding, _ = tokenize.detect_encoding(io.BytesIO(data).readline)
    text = data.decode(encoding)

    tree = _parse_as_ast(text=text, path=source_path)
    _validate_source_statements(tree=tree, text=text, source_path=source_path)
    line_starts = _line_start_offsets(text)

    species_nodes = dict()
    species_reactive = dict()
    species_structure_nodes = dict()
    species_spin_mult_nodes = dict()
    ts_nodes = dict()
    network_nodes = list()
    pdep_nodes = list()
    kinetics_labels_emitted = list()
    edits = list()

    for node in tree.body:
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
            continue
        call = node.value
        call_name = _get_call_name(call)
        if call_name not in _TOP_LEVEL_CALL_NAMES:
            continue
        kwargs = {kw.arg: kw.value for kw in call.keywords if kw.arg is not None}

        if call_name == 'species':
            label = _literal_or_none(kwargs.get('label'))
            if label is not None:
                species_nodes[label] = node
                reactive_node = kwargs.get('reactive')
                species_reactive[label] = _REACTIVE_ABSENT if reactive_node is None \
                    else _literal_or_none(reactive_node)
                species_structure_nodes[label] = kwargs.get('structure')
                species_spin_mult_nodes[label] = kwargs.get('spinMultiplicity')

        elif call_name == 'transitionState':
            label = _literal_or_none(kwargs.get('label'))
            if label is not None:
                ts_nodes[label] = node

        elif call_name == 'network':
            network_nodes.append(node)

        elif call_name == 'pressureDependence':
            pdep_nodes.append(node)

        elif call_name == 'reaction':
            if 'kinetics' in kwargs:
                continue
            label = _literal_or_none(kwargs.get('label'))
            if label is None:
                unevaluable = ast.get_source_segment(text, kwargs['label']) if 'label' in kwargs else None
                raise ValueError(f"The source at '{source_path}' declares a reaction() block with no explicit "
                                 f"'kinetics=' keyword, and its 'label' keyword {'could not be evaluated '
                                 'literally: ' + repr(unevaluable) if unevaluable is not None else 'is missing'}. "
                                 f"Refusing to write an explorer input that would end up with neither an explicit "
                                 f"kinetics entry nor a 'kinetics(...)' job directive for this reaction.")
            _validate_kinetics_label(label=label, source_path=source_path)
            kinetics_labels_emitted.append(label)

    if len(pdep_nodes) != 1:
        raise ValueError(f"Expected exactly one 'pressureDependence(...)' block in '{source_path}', found "
                         f"{len(pdep_nodes)}. Refusing to write an explorer input built on an unresolved (zero) "
                         f"or ambiguous (more than one) pressureDependence block.")
    pdep_node = pdep_nodes[0]

    # Resolve the configured seed(s) to the label THIS network actually uses BEFORE validating them.
    # The seed is configured once, against the original source network, where it is a literal
    # species() label ('[O]C=O'). From round 1 on the loop re-explores a HYBRID network whose species
    # Arkane relabelled with network indices ('[O]C=O(1)'), so the configured label is no longer a
    # literal label here and _validate_seed_species (by-label, like Arkane's own source resolution)
    # would refuse it. Resolution is by structure, never by loosening that validator -- a seed that
    # cannot be resolved is passed through unchanged, so the validator still raises its clear error.
    seed_species, seed_resolution_warnings = _resolve_seed_labels(
        seed_species=seed_species, species_nodes=species_nodes,
        species_structure_nodes=species_structure_nodes, source_path=source_path)
    _validate_seed_species(seed_species=seed_species, species_nodes=species_nodes, ts_nodes=ts_nodes)
    # Rebound from the return value, not merely called for its side effect: an integer-only field
    # coerces an integral float (2.0 -> 2), and rendering the un-coerced argument would write a count
    # as a non-count into the generated file.
    field_values = validate_explorer_field_values(explore_tol=explore_tol, energy_tol=energy_tol,
                                                  flux_tol=flux_tol,
                                                  maximum_radical_electrons=maximum_radical_electrons)
    explore_tol, energy_tol = field_values['explore_tol'], field_values['energy_tol']
    flux_tol, maximum_radical_electrons = field_values['flux_tol'], field_values['maximum_radical_electrons']
    # Checked here rather than left to ``METHOD_MAP[method]`` inside
    # ``t3.utils.writer.rewrite_arkane_method_line``, where an unknown or None method raises a bare
    # KeyError with no field name, no valid set, and no indication that the caller chose it -- and it
    # raises only after the source has been read, validated and edited. Same reasoning as the numeric
    # contracts above: the field has a small closed set of legal values, so it is knowable up front.
    if method not in METHOD_MAP:
        raise ValueError(f"The master-equation 'method' must be one of {sorted(METHOD_MAP)}, got "
                         f"{method!r}. Arkane's own name for it is written into the kept "
                         f"pressureDependence(...) block, so an unrecognized method cannot be rendered.")
    warnings_list, reactive_injection_labels = _validate_bath_gas(
        bath_gas=bath_gas, species_nodes=species_nodes,
        species_reactive=species_reactive, source_path=source_path)
    # Seed resolution runs before warnings_list exists, so fold its notes in now; they belong in the
    # same ExplorerInputSummary.warnings the caller surfaces.
    warnings_list = list(seed_resolution_warnings) + warnings_list

    # Inject 'reactive = False,' into each declared bath-gas species() block that carries no
    # 'reactive' keyword in the source -- see _validate_bath_gas's docstring for why this is a
    # faithful translation of the caller's declaration rather than a rewrite of user input, and why
    # no RMG-written source ever carries the keyword itself (issue #183). The insertion lands at
    # the first keyword's own position (a keyword argument's placement inside the call is
    # semantically free), so the edit is computed structurally from the AST like every other edit
    # here, never by text-scanning for a species block.
    for label in reactive_injection_labels:
        species_call = species_nodes[label].value
        first_kw = species_call.keywords[0]
        insert_pos = _pos(text, line_starts, first_kw.lineno, first_kw.col_offset)
        edits.append((insert_pos, insert_pos, f"reactive = False,\n{' ' * first_kw.col_offset}"))

    # Correct any spinMultiplicity that contradicts its species' own adjacency-list multiplicity.
    # RMG's explorer writes triplet atomic O ('[O]', thermo label "O(T)") with spinMultiplicity = 1
    # while its adjacency-list header says 'multiplicity 3' (u2, two unpaired electrons). Splicing
    # that block verbatim into the next round's explorer input makes Arkane re-read an
    # internally-inconsistent species and crash on Hund's rule -- so any RMG-explored network
    # carrying O(3P) is un-re-feedable until this is fixed. This is a CONSISTENCY correction, not
    # invented statmech: the adjacency list is the species' authoritative electronic structure, so
    # when an explicit 'multiplicity N' header disagrees with spinMultiplicity, the header wins and
    # spinMultiplicity is rewritten to N (energies/frequencies untouched). SMILES-only structures
    # carry no explicit header and are left alone (Arkane derives their multiplicity itself). Placed
    # in T3 rather than in RMG-Py deliberately: T3 must be robust to the network files stock RMG
    # already writes, and cannot assume a patched RMG on every box -- see the ticket's placement
    # argument. The correction is recorded (and warned) so it is never silent.
    spin_multiplicity_corrected = list()
    for label, structure_node in species_structure_nodes.items():
        adj_multiplicity = _adjacency_list_multiplicity(structure_node)
        if adj_multiplicity is None:
            continue
        spin_mult_node = species_spin_mult_nodes.get(label)
        spin_mult_value = _literal_or_none(spin_mult_node) if spin_mult_node is not None else None
        if not isinstance(spin_mult_value, int) or isinstance(spin_mult_value, bool):
            continue
        if spin_mult_value == adj_multiplicity:
            continue
        start = _pos(text, line_starts, spin_mult_node.lineno, spin_mult_node.col_offset)
        end = _pos(text, line_starts, spin_mult_node.end_lineno, spin_mult_node.end_col_offset)
        edits.append((start, end, str(adj_multiplicity)))
        spin_multiplicity_corrected.append((label, spin_mult_value, adj_multiplicity))
        message = (f"Species {label!r}: spinMultiplicity {spin_mult_value} contradicts its "
                   f"adjacency-list 'multiplicity {adj_multiplicity}' header; corrected to "
                   f"{adj_multiplicity} so Arkane can re-read the species (RMG's explorer writes "
                   f"triplet atomic O as a singlet). Energies and frequencies are untouched.")
        warnings_list.append(message)

    # Remove every network(...) block: Arkane's explorer input parser does not recognize it, and
    # keeping it would (at best) be dead text and (at worst) confuse a human re-reading the file.
    for node in network_nodes:
        start = _pos(text, line_starts, node.lineno, 0)
        end = _consume_trailing_newline(text, _pos(text, line_starts, node.end_lineno, node.end_col_offset))
        edits.append((start, end, ''))

    # Insert the kinetics(...) job directives immediately before the pressureDependence(...) block,
    # in the order their reactions were encountered.
    if kinetics_labels_emitted:
        insert_pos = _pos(text, line_starts, pdep_node.lineno, 0)
        # Rendered with !r (not an f-string interpolated into a quoted literal) so the label is
        # always a well-formed Python string literal; _validate_kinetics_label above is the
        # belt-and-suspenders check that also refuses a label containing a character that could
        # break out of (or corrupt) that literal in the first place.
        kinetics_text = ''.join(f'kinetics({label!r})\n' for label in kinetics_labels_emitted)
        edits.append((insert_pos, insert_pos, kinetics_text))

    # Apply edits back-to-front (by descending start offset) so earlier, not-yet-applied offsets
    # (all computed against the pristine, pre-edit text) stay valid.
    for start, end, replacement in sorted(edits, key=lambda edit: edit[0], reverse=True):
        text = text[:start] + replacement + text[end:]

    text = _rewrite_method_line(text=text, method=method, source_path=source_path)
    text = _build_database_block(database_kwargs) + text
    text = text + _build_explorer_block(seed_species=seed_species, bath_gas=bath_gas, explore_tol=explore_tol,
                                        energy_tol=energy_tol, flux_tol=flux_tol,
                                        maximum_radical_electrons=maximum_radical_electrons)

    # Self-check, part 1: the generated input file's own text must itself be valid Python before it
    # is ever written to disk. A splice bug elsewhere in this module could otherwise hand back a
    # broken (or subtly wrong) explorer input without any indication something went wrong.
    try:
        generated_tree = ast.parse(text)
    except SyntaxError as e:
        raise RuntimeError(f"Refusing to write '{dest_path}': the generated explorer input file text failed its "
                           f"own self-check (it is not valid Python): {e}.") from e

    # Self-check, part 2: ast.parse succeeding only proves the file PARSES, not that Arkane can
    # LOAD it -- ``inf``/``nan`` and a constructor-style repr() (e.g. pathlib.Path(...)) are both
    # valid syntax but undefined names at load time (see _render_literal). This re-walks only the
    # blocks THIS MODULE generated -- the prepended database(...) and appended explorer(...)
    # directives -- and asserts every keyword-argument value is ast.literal_eval-able, i.e. an
    # actual Python literal rather than a name/call/attribute that would raise when Arkane loads
    # the file. The spliced-in source network text (species/reaction/pressureDependence/...) is
    # deliberately NOT checked this way: it legitimately contains non-literal expressions (e.g.
    # SMILES(...), Arrhenius(...)) that ast.literal_eval correctly refuses.
    #
    # Part 2a, FIRST: the generated file's top-level statement list as a whole. Checking only the two
    # generated directives' keyword values was a blind spot, and it was exploited: a crafted
    # ``database_kwargs`` KEY (interpolated outside any quoting, unlike the values) closed the
    # generated ``database(...)``, injected a top-level statement, and reopened a ``database(`` that
    # the rendered tail completed. The result parsed, and the loop below skipped the injected
    # statement because it was not a ``database``/``explorer`` call. ``_build_database_block`` now
    # refuses a non-identifier key, which closes that specific route -- this asks the more general
    # question the self-check should have been asking all along: is EVERY statement in the file about
    # to be written one this module meant to be there? A splice bug that manufactures a statement is
    # exactly what a self-check is for, and there is no reason to enumerate the ways it might.
    _validate_generated_statements(tree=generated_tree, text=text, dest_path=dest_path)
    for node in generated_tree.body:
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
            continue
        if _get_call_name(node.value) not in ('database', 'explorer'):
            continue
        for kw in node.value.keywords:
            if kw.arg is None:
                continue
            try:
                ast.literal_eval(kw.value)
            except (ValueError, TypeError) as e:
                unevaluable = ast.get_source_segment(text, kw.value)
                raise RuntimeError(f"Refusing to write '{dest_path}': the generated '{_get_call_name(node.value)}"
                                   f"(...)' directive's '{kw.arg}' keyword value {unevaluable!r} failed its own "
                                   f"self-check (it is not a Python literal Arkane can load): {e}.") from e

    dest_dir = os.path.dirname(dest_path)
    if dest_dir and not os.path.isdir(dest_dir):
        os.makedirs(dest_dir)

    # Written atomically -- staged into a temp file in the SAME directory as dest_path, fsync'd,
    # then os.replace'd into place, mirroring ``t3.pdep.hybrid.write_hybrid_network_input_file``
    # (which exists precisely so a crash cannot leave a torn, partially-written input.py at
    # dest_path).
    staging_dir = dest_dir or '.'
    fd, staged_path = tempfile.mkstemp(prefix='.explorer-input-', suffix='.py', dir=staging_dir)
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
        # durable.
        dir_fd = os.open(staging_dir, os.O_RDONLY)
        try:
            os.fsync(dir_fd)
        finally:
            os.close(dir_fd)

    return ExplorerInputSummary(
        dest_path=dest_path,
        seed_species=seed_species,
        bath_gas=dict(bath_gas),
        kinetics_labels_emitted=tuple(kinetics_labels_emitted),
        reactive_false_injected=tuple(reactive_injection_labels),
        warnings=tuple(warnings_list),
        spin_multiplicity_corrected=tuple(spin_multiplicity_corrected),
    )


def _validate_seed_species(seed_species: tuple, species_nodes: dict, ts_nodes: dict) -> None:
    """
    Validate the explorer seed against Arkane's own source-resolution limitation (P2, P3).

    Args:
        seed_species (tuple): The seed (source) species labels.
        species_nodes (dict): species() label -> AST node, as collected from the source.
        ts_nodes (dict): transitionState() label -> AST node, as collected from the source.

    Raises:
        ValueError: If a seed label is defined as both a species() and a transitionState() block
                   (ambiguous); if a seed label is a transitionState() label only (Arkane's
                   'source' is resolved from species_dict only -- a transitionState() label there
                   raises a KeyError); or if a seed label is not defined at all.
    """
    for label in seed_species:
        is_species = label in species_nodes
        is_ts = label in ts_nodes
        if is_species and is_ts:
            raise ValueError(f"Seed species label '{label}' is defined as BOTH a species() and a "
                             f"transitionState() block; this is ambiguous and Arkane's explorer 'source' can only "
                             f"resolve a label against species_dict, never a transition state.")
        if is_ts and not is_species:
            raise ValueError(f"Seed species label '{label}' is defined as a transitionState() block, not a "
                             f"species() block. Arkane's explorer 'source' is resolved from species_dict only (see "
                             f"P2); a transitionState() label there raises a KeyError in Arkane, so this is refused "
                             f"here first.")
        if not is_species:
            raise ValueError(f"Seed species label '{label}' is not defined as a species() block anywhere in the "
                             f"source. Known species labels are {sorted(species_nodes)}.")


def _molecule_identity(molecule) -> str:
    """
    The canonical structural identity of an ARC ``Molecule``: its canonical SMILES joined to its spin
    multiplicity.

    Mirrors ``t3.pdep.pes_rounds._canonical_structure`` so a seed built from a SMILES and a network
    species built from an ``adjacencyList`` that are the same physical structure (same connectivity,
    same multiplicity) produce the same string here.
    """
    return f'{molecule.to_smiles()}|m{molecule.multiplicity}'


def _seed_structure_identity(label: str) -> str | None:
    """
    The canonical structure identity of a configured seed LABEL read as a SMILES.

    T3's PES seeds are SMILES labels (e.g. ``'[O]C=O'``), so a seed that is no longer a literal
    ``species()`` label in a re-explored (hybrid) network can still be matched to the label that
    network uses, by structure. Returns ``None`` when the label is not a parseable SMILES -- an
    index-suffixed label such as ``'O=C=O(5)'`` or any non-structural name -- in which case
    ``_resolve_seed_labels`` leaves the seed untouched for ``_validate_seed_species`` to judge by
    string.

    Args:
        label (str): The configured seed (source) label.

    Returns:
        Optional[str]: The canonical structure identity, or ``None`` if the label is not a SMILES.
    """
    try:
        molecule = Molecule().from_smiles(label)
    except Exception:
        # Deliberately broad, and the breadth is the point. This function's contract is "None when
        # the label is not a parseable SMILES", and every caller depends on that failing CLOSED --
        # an unresolved seed is passed through untouched so `_validate_seed_species` can refuse it
        # by name, with a message naming the labels the network actually defines. A parser
        # exception escaping here would replace that diagnostic with a raw backend error from
        # somewhere inside openbabel. The backend raises more than ValueError/KeyError in practice
        # (a non-str label surfaces as TypeError), and it is third-party code whose exception set
        # is not ours to enumerate. `pes_rounds._canonical_structure` guards its own parse the same
        # way, for the same reason.
        return None
    return _molecule_identity(molecule)


def _species_structure_identity(structure_node) -> str | None:
    """
    The canonical structure identity of a ``species()`` block's ``structure`` keyword value.

    Mirrors ``_adjacency_list_multiplicity``'s extraction -- the same two ``structure`` call forms
    real network files use, ``adjacencyList(...)`` and ``SMILES(...)`` -- but resolves the whole
    structure to the identity seed resolution compares against, not merely its multiplicity.

    Args:
        structure_node: The AST node of the species' ``structure`` keyword value, or ``None``.

    Returns:
        Optional[str]: The canonical structure identity, or ``None`` when the structure is absent or
        is not a parseable ``adjacencyList(...)``/``SMILES(...)`` string literal.
    """
    if not isinstance(structure_node, ast.Call) or not structure_node.args:
        return None
    call_name = _get_call_name(structure_node)
    arg = structure_node.args[0]
    text = arg.value if isinstance(arg, ast.Constant) and isinstance(arg.value, str) else None
    if text is None:
        return None
    try:
        if call_name == 'adjacencyList':
            molecule = Molecule().from_adjacency_list(text)
        elif call_name == 'SMILES':
            molecule = Molecule().from_smiles(text)
        else:
            return None
    except Exception:
        # Broad for the same reason as `_seed_structure_identity`: an unparseable `structure`
        # keyword must leave this species unmatched rather than abort the whole write.
        return None
    return _molecule_identity(molecule)


def _resolve_seed_labels(seed_species: tuple, species_nodes: dict, species_structure_nodes: dict,
                         source_path: str) -> tuple:
    """
    Resolve each configured seed (source) label to the label THIS network actually uses, matching by
    molecular structure whenever the configured label is not itself a ``species()`` label here.

    T3's PES loop configures its seed once, against the original source network, where the seed is a
    literal ``species()`` label (``'[O]C=O'``). From round 1 on, the loop re-explores a HYBRID network
    whose species Arkane's explorer relabelled with network indices (``'[O]C=O(1)'``), so the
    configured label is no longer a literal label in the file and Arkane's by-label ``source``
    resolution -- which ``_validate_seed_species`` refuses first, with a clearer message -- would
    raise. The information to bridge the two labels is in the file: every species carries
    ``structure = adjacencyList(...)``, so the configured seed can be matched to the round's actual
    label by canonical structure rather than by string.

    Resolution is per label, and never loosens ``_validate_seed_species``:

    * **Literal match** -- the label is already a ``species()`` label here (round 0 against the real
      source): keep it, with no structural work at all, so round-0 behaviour is byte-for-byte
      unchanged.
    * **Structural match** -- otherwise, if the configured label parses as a species structure (a
      SMILES; the campaign's seeds are SMILES labels), match its identity against every
      ``species()`` block's structure. Exactly one match resolves to that label. **Zero** matches
      leave the configured label unchanged, so ``_validate_seed_species`` raises its own clear
      "not defined ... known labels are [...]" error rather than this function inventing one.
      **More than one** match is refused (ambiguous): a source configuration must name exactly one
      physical species, and seeding from an arbitrary one of several structurally identical species
      could explore from the wrong copy.
    * A label that neither matches literally nor parses as a structure likewise falls through to
      ``_validate_seed_species`` unchanged.

    Args:
        seed_species (tuple): The configured seed (source) labels.
        species_nodes (dict): ``species()`` label -> AST node, as collected from the source.
        species_structure_nodes (dict): ``species()`` label -> ``structure`` keyword AST node.
        source_path (str): The network path the seed is being resolved against (for messages).

    Returns:
        tuple: ``(resolved_seed_species, warnings)`` -- the resolved labels in input order, and any
        human-readable notes about labels that were remapped.

    Raises:
        ValueError: If a seed label matches more than one ``species()`` block by structure.
    """
    resolved = list()
    resolution_warnings = list()
    network_identities = None  # computed lazily, only if some label actually needs structural work
    for label in seed_species:
        if label in species_nodes:
            resolved.append(label)
            continue
        seed_identity = _seed_structure_identity(label)
        if seed_identity is None:
            resolved.append(label)
            continue
        if network_identities is None:
            network_identities = {net_label: _species_structure_identity(node)
                                  for net_label, node in species_structure_nodes.items()}
        matches = sorted(net_label for net_label, identity in network_identities.items()
                         if identity is not None and identity == seed_identity)
        if len(matches) == 1:
            resolution_warnings.append(
                f"Seed species label {label!r} is not a species() block in '{source_path}'; resolved "
                f"it by structure ({seed_identity}) to {matches[0]!r}, the label this network uses for "
                f"that species.")
            resolved.append(matches[0])
        elif len(matches) > 1:
            raise ValueError(
                f"Seed species label {label!r} matches {len(matches)} species() blocks in "
                f"'{source_path}' by structure ({seed_identity}): {matches}. A source configuration "
                f"must identify exactly one species; refusing to seed the explorer from an arbitrary "
                f"one of several structurally identical species.")
        else:
            resolved.append(label)
    return tuple(resolved), resolution_warnings


def _validate_kinetics_label(label: str, source_path: str) -> None:
    """
    Refuse a reaction label that could break out of (or corrupt) the quoted string literal it is
    rendered into for a ``kinetics('<label>')`` job directive.

    The label is read verbatim out of ``source_path`` (an untrusted RMG/hybrid network file), and
    then interpolated into a NEW, executable Arkane input file this module writes. Rendering it
    with ``!r`` (see the call site in ``write_arkane_explorer_input_file``) already guarantees a
    well-formed Python string literal for any ``str`` value, but this check is the belt-and-
    suspenders guard against a label containing a quote, newline or backslash, mirroring
    ``t3.pdep.hybrid._validate_no_injection_chars``.

    Args:
        label (str): The reaction label to validate.
        source_path (str): The source network path the label was read from (used only in the
                           error message).

    Raises:
        ValueError: If ``label`` contains a newline, quote character or backslash.
    """
    if any(char in label for char in _KINETICS_LABEL_INJECTION_CHARS):
        raise ValueError(f"Reaction label {label!r} (from '{source_path}') must not contain a newline, quote "
                         f"character or backslash: it would be rendered into a kinetics('<label>') job directive "
                         f"in the generated, executable Arkane explorer input file, and such a character there "
                         f"could inject or corrupt directives in that file.")


def _validate_number(field_name: str, value, minimum, minimum_inclusive: bool, integer_only: bool,
                     rationale: str):
    """
    Check one numeric value against ITS OWN field's contract, not merely against "is a literal".

    ``_render_literal`` asks whether a value can be written as a Python literal that will load. This
    asks the caller's other question -- whether the value is a plausible one for this particular
    field -- and the two have different answers: ``True`` and ``'2'`` are both perfectly renderable
    and both are nonsense as ``maximumRadicalElectrons``. Neither check subsumes the other, so both
    run.

    ``bool`` is refused first and explicitly. It is the trap that makes a type check look like it is
    working when it is not: ``bool`` is an ``int`` subclass, so ``isinstance(True, int)`` is True and
    ``math.isfinite(True)`` is True, and ``True`` renders as the valid literal ``True`` -- which
    Arkane then uses in arithmetic as 1. Nothing anywhere raises. The same reasoning appears in
    ``_is_numeric_expression`` for the AST side of this module; it is written out in both places
    rather than shared, because the two operate on different kinds of thing (a value vs. a node).

    A non-finite float is diagnosed as non-finite even for an integer-only field, because that is the
    more specific and more useful message: ``repr(float('inf'))`` is the bare name ``inf``, valid
    SYNTAX that this module's ``ast.parse`` self-check accepts but an undefined NAME the moment
    Arkane loads the file.

    An integer-only field accepts a float that IS a whole number and returns it coerced to ``int``.
    "Integer-only" is a constraint on the VALUE, not on how the caller spelled it, and the spelling
    is not the caller's choice to begin with: this value arrives from T3's YAML config, where ``2.0``
    is a float. Refusing it would refuse a config asking for precisely what an accepted ``2`` asks
    for. The coercion is what makes this a widening rather than a hole -- ``2.5`` is still refused,
    and so is ``True``, which is admitted by ``float(value).is_integer()`` and would otherwise render
    a silent count of 1 (hence the bool check above running first).

    Args:
        field_name (str): The field the value came from, used in the error message.
        value: The value to check. Never None (callers skip None, which omits the keyword).
        minimum: The smallest permitted value.
        minimum_inclusive (bool): Whether ``minimum`` itself is permitted.
        integer_only (bool): Whether the field is a count, and so admits no fractional value.
        rationale (str): Why this field's bound is what it is, quoted into the error message so the
                        refusal explains itself with evidence rather than an assertion.

    Returns:
        The validated value, coerced to ``int`` if the field is integer-only and the value was an
        integral float. Callers must render THIS value, not the one they passed in.

    Raises:
        ValueError: If the value is a bool, a non-finite float, not a number of the required kind, or
                   outside the field's permitted range.
    """
    if isinstance(value, bool):
        raise ValueError(f"'{field_name}' must be a number, got the bool {value!r}. A bool is an int subclass, so "
                         f"it passes every numeric type check and every finiteness check, renders as the valid "
                         f"literal {value!r}, and is then silently used by Arkane as {int(value)} -- with nothing "
                         f"raising anywhere. {rationale}")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"'{field_name}' must be a finite number, got {value!r}. Arkane's input file cannot carry "
                         f"a non-finite literal (it would be written as the undefined name 'inf'/'nan' and raise "
                         f"when Arkane loads the file). To apply no filter, omit '{field_name}' entirely (None) "
                         f"rather than passing infinity, which leaves the keyword out and lets Arkane apply its "
                         f"own default.")
    if integer_only:
        if isinstance(value, float):
            # Safe to call: a non-finite float was already refused above, and .is_integer() is the
            # exact question being asked -- 2.0 is a count spelled as a float, 2.5 is not a count.
            if not value.is_integer():
                raise ValueError(f"'{field_name}' must be a whole number, got {value!r}. {rationale}")
            value = int(value)
        elif not isinstance(value, int):
            raise ValueError(f"'{field_name}' must be an int, got {value!r} of type {type(value).__name__}. "
                             f"{rationale}")
    elif not isinstance(value, (int, float)):
        raise ValueError(f"'{field_name}' must be a number, got {value!r} of type {type(value).__name__}. A str "
                         f"renders as a perfectly valid quoted literal and reaches Arkane, which then raises a "
                         f"TypeError in arithmetic deep inside the job, after the exploration has been set up. "
                         f"{rationale}")
    if value < minimum or (value == minimum and not minimum_inclusive):
        boundary = f'>= {minimum}' if minimum_inclusive else f'> {minimum}'
        raise ValueError(f"'{field_name}' must be {boundary}, got {value!r}. {rationale}")
    return value


# Per-field contracts for the numeric values this writer renders into ``explorer(...)``. The bounds
# are read off what Arkane DOES with each value, not invented: see the rationale on each entry. The
# observed corpus ranges (examples/arkane/explorer/, arkane/data/methoxy_explore.py) all satisfy
# these, and ``TestFieldValuesAreValidatedPerFieldNotJustPerLiteralType`` pins every one of them so a
# future tightening cannot quietly reject the inputs Arkane itself ships.
_EXPLORER_NUMBER_CONTRACTS = {
    'explore_tol': {
        'minimum': 0.0, 'minimum_inclusive': True, 'integer_only': False,
        'rationale': "Arkane explores every isomer whose leak exceeds 'explore_tol * kchar' "
                     "(arkane/explorer.py:235). A negative tolerance makes that threshold negative, so "
                     "a leak coefficient of exactly zero exceeds it and a network with no leak at all "
                     "is explored -- not a looser tolerance but a different comparison. Zero is "
                     "permitted: it means 'explore on any leak at all', which is effectively unbounded "
                     "but well defined, and refusing it would guard nothing while an equally unbounded "
                     "1e-300 is accepted. Bounding that cost is a caller-side budget question, not "
                     "something a lower bound on this field can express.",
    },
    'energy_tol': {
        'minimum': 0.0, 'minimum_inclusive': True, 'integer_only': False,
        'rationale': "'energy_tol' is a dimensionless multiple of RT, not an energy: Arkane forms "
                     "'dE = tol * R * T' (rmgpy/rmg/pdep.py:331), collects every reaction with "
                     "'E0 - E0source > dE' (:353), and REMOVES that set (arkane/explorer.py:336). A "
                     "negative tolerance therefore removes channels lying BELOW the source energy as "
                     "well, which is not a tighter version of the documented filter ('greater in Free "
                     "Energy than tol*R*T + Gf_source') but the opposite one. Zero is permitted: it is "
                     "the exact cutoff, Arkane gates on 'energy_tol != np.inf' (:303,:309) and so "
                     "passes it through, and refusing it while accepting an equally destructive 1e-300 "
                     "would refuse one spelling rather than guard anything.",
    },
    'flux_tol': {
        'minimum': 0.0, 'minimum_inclusive': True, 'integer_only': False,
        'rationale': "Zero is permitted because it is Arkane's OWN default (arkane/input.py:496) and "
                     "its documented spelling of 'apply no flux filter' -- the call site guards on "
                     "'if self.flux_tol != 0.0' (arkane/explorer.py:317). A negative tolerance, by "
                     "contrast, can never exclude anything. There is deliberately no upper bound: the "
                     "primary filter compares against an absolute k*c (rmgpy/rmg/pdep.py:380), not "
                     "against a normalized fraction, so '<= 1' would be a bound on the wrong quantity.",
    },
    'maximum_radical_electrons': {
        'minimum': 0, 'minimum_inclusive': True, 'integer_only': True,
        'rationale': "This is a count of electrons, compared with '>' against an actual radical count "
                     "(rmgpy/constraints.py:144-147). A fractional maximum is not a count, and a "
                     "negative one forbids every species including closed-shell ones. Zero is allowed: "
                     "'no radicals at all' is a coherent request.",
    },
}


def validate_explorer_field_values(explore_tol, energy_tol, flux_tol, maximum_radical_electrons) -> dict:
    """
    Check every numeric ``explorer(...)`` field against its own contract in ``_EXPLORER_NUMBER_CONTRACTS``.

    Public (no leading underscore): it is the single numeric-contract rule shared by TWO entry
    points -- this module's own ``write_arkane_explorer_input_file``, and
    ``t3.pdep.explorer.config.PDepExplorerConfig``, which delegates to it in ``__post_init__`` so a
    config built directly (never passed through the writer) is refused at construction time rather
    than only when it is later rendered. Duplicating the per-field contracts in ``config.py`` instead
    would let the two entry points' rules drift apart, which is exactly the antipattern
    ``t3.pdep.explorer.adapter.validate_explorer_seed``'s docstring warns about for the seed rules.
    So this is module-level public API now, not an internal helper.

    A None value is skipped rather than checked: the writer omits the keyword entirely for None, which
    is how a caller asks for Arkane's own default, and validating an absent field would refuse the one
    correct way to spell "no filter".

    Args:
        explore_tol: The exploration tolerance, or None.
        energy_tol: The energy filter tolerance, or None.
        flux_tol: The flux filter tolerance, or None.
        maximum_radical_electrons: The maximum radical electron count, or None.

    Returns:
        dict: The same four fields keyed by name, each either None or its validated value -- which is
              not always the value passed in, since an integer-only field coerces an integral float.
              The caller must render from this mapping, or the coercion is silently discarded.

    Raises:
        ValueError: If any given value violates its field's contract.
    """
    values = {'explore_tol': explore_tol, 'energy_tol': energy_tol, 'flux_tol': flux_tol,
              'maximum_radical_electrons': maximum_radical_electrons}
    for field_name, contract in _EXPLORER_NUMBER_CONTRACTS.items():
        if values[field_name] is not None:
            values[field_name] = _validate_number(field_name=field_name, value=values[field_name], **contract)
    return values


def _validate_bath_gas(bath_gas: dict, species_nodes: dict, species_reactive: dict,
                       source_path: str) -> tuple:
    """
    Validate the bath gas composition against Arkane/RMG's bath-gas identification rule (P16), and
    decide which species blocks need ``reactive = False`` injected into the generated file.

    Arkane/RMG identify the bath gas as every UNREACTIVE core species, never by name
    (``rmgpy/rmg/pdep.py:856-863``), and the literal keyword in the ``species()`` block is the only
    channel that survives: the explorer's own ``make_new_species(spec, reactive=False)``
    (``arkane/explorer.py:141``) is overridden when the species enters the core, because
    ``rmgpy/rmg/model.py:475`` re-reads ``reactive = object.reactive`` from the loaded species
    object. Yet no RMG-written network file ever carries the keyword -- ``arkane/pdep.py:654``
    (``save_input_file``) emits species blocks with no ``reactive`` at all (issue #183). So an
    ABSENT keyword on a declared bath-gas species is not a refusal case but the normal case, and
    since this module GENERATES the destination file, it emits the fact the caller declared
    (``bath_gas`` names this species as the collider) as ``reactive = False`` in the generated
    blocks: a translation between two encodings of one fact, not a rewrite of user input. Only a
    species that EXPLICITLY contradicts the declaration (a literal ``reactive=True``), or whose
    ``reactive`` value cannot be literally evaluated at all, is refused.

    Args:
        bath_gas (dict): Bath gas species label -> mole fraction.
        species_nodes (dict): species() label -> AST node, as collected from the source.
        species_reactive (dict): species() label -> the literal value of its 'reactive' keyword,
                                 or ``_REACTIVE_ABSENT`` when the block carries no such keyword
                                 (``None`` means the keyword exists but its value is not a
                                 literal).
        source_path (str): The source network path, used in error messages.

    Raises:
        ValueError: If the bath gas is empty; if a fraction violates its contract or the
                   composition's sum/equality contracts; if a bath gas label is not a species()
                   block in the source; if a bath gas species carries an explicit literal
                   ``reactive=True`` (a contradiction this writer refuses to resolve either way);
                   or if its ``reactive`` value is not a literal ``True``/``False``.

    Returns:
        tuple: A 2-tuple of
              - list: Warnings, e.g. that a multi-species bath gas's requested fractions are
                recorded but will not be honored by Arkane/RMG (P16).
              - tuple: The bath-gas labels whose species() blocks need ``reactive = False``
                injected (those carrying no ``reactive`` keyword in the source).
    """
    # An EMPTY bath gas is refused here rather than passed through. Every check below is a loop
    # over ``bath_gas``, so an empty dict satisfies all of them vacuously -- the guard would report
    # "valid" precisely when there is nothing to validate. It is not merely unvalidated but
    # unusable: the writer removes the source's ``network(...)`` block, which is the only other
    # place Arkane looks for a bath gas, so ``bathGas={}`` leaves Arkane to fail deep inside the
    # run with ``InputError('bathGas not specified in explorer block')`` (arkane/explorer.py:76)
    # after the exploration has already been set up.
    if not bath_gas:
        raise ValueError(f"A bath gas is required to write an Arkane explorer input file from "
                         f"'{source_path}', got {bath_gas!r}. The source's network(...) block -- "
                         f"Arkane's only other source of a bath gas -- is removed by this writer, "
                         f"so an empty bath gas would fail inside the Arkane run itself.")
    for label, fraction in bath_gas.items():
        # Deliberately STRICTER than Arkane, which checks nothing here: ``arkane/input.py:508`` is a
        # pure passthrough, and no code between it and the collision-frequency calculation normalizes
        # the fractions. ``rmgpy/pdep/configuration.pyx:190-197`` then uses them raw -- as a linear
        # weight for sigma and molecular weight, and as an EXPONENT for epsilon. So a fraction set
        # that does not sum to 1 does not fail: it produces a physically wrong bath gas and reports
        # success. That is precisely the class of error a validator restricted to "what would Arkane
        # itself refuse?" cannot catch, which is why the bound is asserted here instead.
        _validate_number(field_name=f"bath_gas['{label}']", value=fraction, minimum=0.0,
                         minimum_inclusive=False, integer_only=False,
                         rationale="A bath gas mole fraction must be greater than zero and no greater "
                                   "than one; a zero-fraction species is not present at all, and the "
                                   "fractions are used unnormalized as weights and as an exponent "
                                   "(rmgpy/pdep/configuration.pyx:190-197).")
        # There is deliberately NO separate "fraction <= 1" check. It cannot fire: every fraction is
        # already required to be > 0 and the composition is required to sum to 1, so no individual
        # fraction can exceed 1 -- except inside the sum check's rounding tolerance, where refusing it
        # would be a false positive against a composition the caller wrote correctly. A guard that is
        # unreachable except where it is wrong is worse than no guard: it reads as a defense while
        # defending nothing, and it would have to be kept correct forever.
        if label not in species_nodes:
            raise ValueError(f"Bath gas label '{label}' is not defined as a species() block anywhere in the "
                             f"source. Known species labels are {sorted(species_nodes)}.")

    # Decided per label AFTER the membership loop above so a mixed error (one unknown label, one
    # conflicted one) always reports the missing block first -- the more fundamental problem.
    inject_labels = list()
    for label in bath_gas:
        reactive_value = species_reactive[label]
        if reactive_value is _REACTIVE_ABSENT:
            # The RMG-written shape (arkane/pdep.py:654 emits no 'reactive' for any species):
            # marked for injection into the generated file, see this function's docstring.
            inject_labels.append(label)
        elif reactive_value is True:
            raise ValueError(f"Bath gas label '{label}' is a species() block carrying an explicit literal "
                             f"'reactive=True', which contradicts its declaration as bath gas. Per P16 "
                             f"(rmgpy/rmg/pdep.py:856-863), Arkane/RMG identify the bath gas as every "
                             f"unreactive core species, and the species block's own keyword is the only "
                             f"channel that survives into the core (rmgpy/rmg/model.py:475 re-reads "
                             f"'reactive' from the loaded species, overriding arkane/explorer.py:141's "
                             f"make_new_species(..., reactive=False)). Honouring the bath-gas declaration "
                             f"would silently rewrite the source's own explicit claim, and honouring the "
                             f"claim would have Arkane generate statmech for the collider -- so this "
                             f"contradiction is refused rather than resolved either way. Fix the source, "
                             f"or drop '{label}' from bath_gas.")
        elif reactive_value is not False:
            raise ValueError(f"Bath gas label '{label}' is a species() block whose 'reactive' keyword is not "
                             f"a literal True or False, so it cannot be verified either way: injecting "
                             f"'reactive = False' beside it could contradict whatever the expression "
                             f"evaluates to when Arkane loads the file, and leaving it is a bath gas "
                             f"Arkane may never recognize as one (P16, rmgpy/rmg/pdep.py:856-863). Make it "
                             f"a literal, or remove it so this writer can inject the declared fact itself.")

    # Checked once over the whole composition rather than per species, because it is the only property
    # here that no individual fraction can violate on its own. ``math.isclose`` with an absolute
    # tolerance, not ``== 1``: {'He': 0.3, 'Ar': 0.3, 'N2': 0.4} does not sum to exactly 1.0 in binary
    # floating point, and refusing a composition a chemist wrote correctly would be a false positive
    # the caller cannot fix.
    fraction_sum = math.fsum(bath_gas.values())
    if not math.isclose(fraction_sum, 1.0, abs_tol=1e-6):
        raise ValueError(f"The bath gas mole fractions must sum to 1, got {fraction_sum!r} for {bath_gas!r}. "
                         f"Arkane does not check or normalize this: the fractions reach "
                         f"rmgpy/pdep/configuration.pyx:190-197 as raw weights (and, for epsilon, as an "
                         f"exponent), so a composition that does not sum to 1 yields a physically wrong bath gas "
                         f"sigma, epsilon and molecular weight -- and the run reports success.")

    # Refused, not warned about. Per P16 the requested fractions of a MULTI-species bath gas are
    # discarded: ``arkane/explorer.py:187`` assigns them onto the network and then ``:274``/``:302``/
    # ``:346`` call ``network.update(...)``, whose body at ``rmgpy/rmg/pdep.py:857-863``
    # unconditionally does ``self.bath_gas = {}`` and then ``1.0 / len(bath_gas)`` over every
    # unreactive CORE species -- with no rmgmode guard on that block. So {'He': 0.9, 'Ar': 0.1} sums
    # to 1, passes every check above, is written into the input file, is recorded in the manifest as
    # what ran, and the run uses 0.5/0.5. A warning is the wrong instrument for a request that cannot
    # be satisfied at all: the file and the provenance would both still claim 0.9/0.1.
    #
    # Equal fractions are accepted because they are what the run will actually apply. Compared with a
    # tolerance rather than ``==``: 1/3 + 1/3 + 1/3 is exactly equal-weighted and not exactly 1.0 in
    # binary floating point, and refusing a composition a chemist wrote correctly is a false positive
    # the caller cannot act on.
    if len(bath_gas) > 1:
        equal_fraction = 1.0 / len(bath_gas)
        unequal = {label: fraction for label, fraction in bath_gas.items()
                   if not math.isclose(fraction, equal_fraction, abs_tol=1e-6)}
        if unequal:
            raise ValueError(
                f"A multi-species bath gas must specify EQUAL mole fractions "
                f"({equal_fraction} each for {len(bath_gas)} species), got {bath_gas!r} -- "
                f"{sorted(unequal)} differ. Per P16 the requested fractions would not be honoured: "
                f"Arkane assigns them onto the network (arkane/explorer.py:187) and then calls "
                f"network.update(...) (arkane/explorer.py:274), which unconditionally overwrites the "
                f"whole composition with equal weights over every unreactive core species "
                f"(rmgpy/rmg/pdep.py:857-863). Writing these fractions would put a number in the input "
                f"file, and in this run's provenance, that the run does not use. Ask for equal "
                f"fractions, or use a single-species bath gas.")

    warnings_list = list()
    if len(bath_gas) > 1:
        # The fractions are no longer the risk -- unequal ones are refused above, and equal ones are
        # what the run applies. What survives is the SET: the equal weights are spread over every
        # unreactive species in the CORE, which is not necessarily the set named here. So a two-species
        # request can still be honoured as a three-species one, at 1/3 each, with nothing raising.
        warnings_list.append(f"bath_gas specifies {len(bath_gas)} species {sorted(bath_gas)} with equal "
                             f"fractions. Per P16, PDepNetwork.update (rmgpy/rmg/pdep.py:857-863) does not "
                             f"read this composition at all: it assigns equal weights over every UNREACTIVE "
                             f"CORE species. The fractions therefore match what will run, but the species SET "
                             f"may not -- if the core carries other unreactive species, they are included and "
                             f"the per-species weight changes accordingly.")
    return warnings_list, tuple(inject_labels)


def _rewrite_method_line(text: str, method: str, source_path: str) -> str:
    """
    Rewrite the kept ``pressureDependence(...)`` block's ``method = '...'`` line.

    Reuses ``t3.utils.writer.rewrite_arkane_method_line`` and its ``METHOD_LINE_CANDIDATE_RE``
    guard, mirroring ``t3.pdep.hybrid._rewrite_method_and_sensitivity``'s approach (minus the
    sensitivity-conditions injection, which has no place in an explorer input).

    The candidate lines are restricted to the ``pressureDependence(...)`` block's own line span,
    found by re-parsing the (already edited) text rather than by scanning it. Scanning every line was
    wrong in a way the "exactly one match" guard hid rather than caught: if the real block carries no
    ``method`` line while some line elsewhere in the file matches the pattern -- a module docstring
    containing one, a comment, a string value -- the count is still exactly 1, so the rewrite lands
    outside the block and Arkane runs with an unresolved method. The count guard only catches the case
    where the stray line is an EXTRA one. Anchoring to the node makes the two cases the same case.

    Args:
        text (str): The input file text (after network(...) removal / kinetics(...) insertion).
        method (str): 'CSE', 'MSC' or 'RS'.
        source_path (str): The source network path, used only in the error message below.

    Raises:
        ValueError: If the ``pressureDependence(...)`` block does not contain exactly one
                   ``method = ...`` line.
        RuntimeError: If the edited text no longer parses, or no longer contains exactly one
                     ``pressureDependence(...)`` block -- either is a splice defect in this module.

    Returns:
        str: The rewritten text.
    """
    try:
        tree = ast.parse(text)
    except SyntaxError as e:
        raise RuntimeError(f"Refusing to write an explorer input for '{source_path}': the text failed its own "
                           f"self-check before the method line could be rewritten (it is not valid Python): "
                           f"{e}. This is a splice defect in this module.") from e
    pdep_nodes = [node for node in tree.body
                  if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                  and _get_call_name(node.value) == 'pressureDependence']
    if len(pdep_nodes) != 1:
        raise RuntimeError(f"Refusing to write an explorer input for '{source_path}': the edited text contains "
                           f"{len(pdep_nodes)} pressureDependence(...) blocks, not 1. The source was already "
                           f"required to have exactly one, so this is a splice defect in this module.")
    # 1-indexed, inclusive, as ast reports them.
    first_line, last_line = pdep_nodes[0].lineno, pdep_nodes[0].end_lineno

    lines = text.splitlines(keepends=True)
    new_lines = list()
    method_rewrite_count = 0
    for line_number, line in enumerate(lines, start=1):
        if first_line <= line_number <= last_line and METHOD_LINE_CANDIDATE_RE.match(line):
            new_lines.append(rewrite_arkane_method_line(line=line, method=method))
            method_rewrite_count += 1
        else:
            new_lines.append(line)

    if method_rewrite_count != 1:
        raise ValueError(f"Expected exactly one 'method = ...' line to rewrite inside the "
                         f"pressureDependence(...) block of '{source_path}' (lines {first_line}-{last_line}), "
                         f"found {method_rewrite_count}. Refusing to write an explorer input with an unresolved "
                         f"(zero) or ambiguous (more than one) method line.")

    return ''.join(new_lines)


def _render_literal(field_name: str, value) -> str:
    """
    Render ``value`` as a Python literal that is guaranteed to LOAD, not merely to parse.

    ``repr(value)`` is not a safe way to render an arbitrary Python object into a generated Arkane
    input file: it merely has to be valid Python SYNTAX to pass this module's ``ast.parse``
    self-check, but Arkane must actually be able to *load* it. Two classes of value defeat plain
    ``!r`` this way: a non-finite float (``repr(float('inf'))`` is the bare name ``inf``, an
    undefined NAME at load time -- see ``_validate_number``), and any object whose ``repr()``
    is a constructor call (e.g. ``pathlib.Path('/x')``, a numpy scalar/array) rather than a literal
    -- valid syntax, but a call to a name (``PosixPath``, ``array``, ...) this file never imports.

    Only ``None``, ``bool``, ``int``, a finite ``float``, ``str``, and ``list``/``tuple``/``dict``
    composed (recursively) of those are accepted; everything else is refused. ``bool`` is checked
    before ``int`` because it is a subclass of ``int``.

    Args:
        field_name (str): The name of the field ``value`` came from (used in the error message).
        value: The value to render.

    Raises:
        ValueError: If ``value`` (or something nested inside it) is not one of the accepted types,
                   is a non-finite float, or (for a dict) has a non-``str`` key.

    Returns:
        str: ``repr(value)``, once ``value`` has been proven to be a safe literal.
    """
    if value is None or isinstance(value, (bool, str)):
        return repr(value)
    if isinstance(value, int):
        return repr(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"'{field_name}' must be a finite number, got {value!r}. A non-finite float renders "
                             f"as the undefined name 'inf'/'nan' and raises when Arkane loads the generated input "
                             f"file.")
        return repr(value)
    if isinstance(value, (list, tuple)):
        rendered = ', '.join(_render_literal(field_name, item) for item in value)
        return f'[{rendered}]' if isinstance(value, list) else f'({rendered},)' if len(value) == 1 \
            else f'({rendered})'
    if isinstance(value, dict):
        parts = list()
        for key, item in value.items():
            if not isinstance(key, str):
                raise ValueError(f"'{field_name}' has a dict with a non-string key {key!r}; refusing to render "
                                 f"it into the generated Arkane input file.")
            parts.append(f'{key!r}: {_render_literal(field_name, item)}')
        return '{' + ', '.join(parts) + '}'
    raise ValueError(f"'{field_name}' has value {value!r} of type {type(value).__name__}, which cannot be "
                     f"rendered as a Python literal guaranteed to load in the generated Arkane input file "
                     f"(its repr() may be a constructor call rather than a literal). Only None, bool, int, "
                     f"finite float, str, and list/tuple/dict composed of those are accepted.")


# The exact signature of Arkane's own ``database()`` (arkane/input.py:83-84), mapped to what each
# keyword accepts THERE -- not to what is renderable. A key can be a flawless Python identifier and
# still be one ``database()`` does not take, in which case the generated file is valid Python that
# raises ``TypeError: database() got an unexpected keyword argument`` the instant Arkane loads it,
# after the whole run has been set up. That is knowable here, so it is refused here.
#
# 'str_or_list_of_str' exists because ``kineticsFamilies`` and ``kineticsDepositories`` genuinely
# accept either a magic string or a list, and validate it themselves (arkane/input.py:92-109); giving
# them the plain list contract of the library keywords would refuse Arkane's own default of
# 'default'. ``frequenciesLibraries`` and ``kineticsEstimator`` are accepted by ``database()`` and
# then never read (verified: each name appears only in the signature) -- refusing them would be
# stricter than Arkane for no benefit, since a caller passing one is asking for a documented no-op,
# not a broken run.
_DATABASE_KWARG_CONTRACTS = {
    'thermoLibraries': 'list_of_str',
    'transportLibraries': 'list_of_str',
    'reactionLibraries': 'list_of_str',
    'frequenciesLibraries': 'list_of_str',
    'kineticsFamilies': 'str_or_list_of_str',
    'kineticsDepositories': 'str_or_list_of_str',
    'kineticsEstimator': 'str',
}

# The magic strings each of the two dual-typed keywords accepts, from their own validation
# (arkane/input.py:92-100 for depositories, :102-109 for families). Anything else must be a list, so a
# typo'd 'defualt' is refused here rather than reaching Arkane's InputError mid-run.
_DATABASE_MAGIC_STRINGS = {
    'kineticsFamilies': ('default', 'all', 'none'),
    'kineticsDepositories': ('default', 'all'),
}


def _validate_database_kwarg(key: str, value) -> None:
    """
    Check one ``database(...)`` keyword against Arkane's real signature and that keyword's own type.

    Args:
        key (str): The keyword name, already proven to be a plain non-dunder identifier by the caller.
        value: The value to be rendered for it.

    Raises:
        ValueError: If ``key`` is not a keyword Arkane's ``database()`` accepts, or ``value`` is not a
                   type that keyword accepts.
    """
    if key not in _DATABASE_KWARG_CONTRACTS:
        raise ValueError(f"Refusing to render a database(...) keyword named {key!r}: Arkane's database() accepts "
                         f"only {sorted(_DATABASE_KWARG_CONTRACTS)} (arkane/input.py:83-84). An unknown keyword "
                         f"is valid Python, so the generated file would parse and pass this module's self-check, "
                         f"and then raise 'TypeError: database() got an unexpected keyword argument' the moment "
                         f"Arkane loads it -- after the run has already been set up.")
    contract = _DATABASE_KWARG_CONTRACTS[key]
    if contract == 'str':
        if not isinstance(value, str):
            raise ValueError(f"The database(...) keyword {key!r} must be a str, got {value!r} of type "
                             f"{type(value).__name__}.")
        return
    if contract == 'str_or_list_of_str' and isinstance(value, str):
        if value not in _DATABASE_MAGIC_STRINGS[key]:
            raise ValueError(f"The database(...) keyword {key!r} accepts the strings "
                             f"{list(_DATABASE_MAGIC_STRINGS[key])} or a list of names, got {value!r}. Arkane "
                             f"checks this itself (arkane/input.py:92-109), so a typo here becomes an InputError "
                             f"once the run has started rather than a refusal now.")
        return
    # Both remaining contracts require a list (of str) at this point. None is permitted because it is
    # Arkane's own default for every library keyword (arkane/input.py:83).
    if value is None:
        return
    if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
        expected = 'a list of str' if contract == 'list_of_str' \
            else f"one of {list(_DATABASE_MAGIC_STRINGS[key])} or a list of str"
        raise ValueError(f"The database(...) keyword {key!r} must be {expected}, got {value!r} of type "
                         f"{type(value).__name__}. A bare str where a list is meant is the dangerous case: "
                         f"Arkane's as_list() would silently treat it as a one-element list, so a run "
                         f"misconfigured this way succeeds while loading something other than what was asked for.")


def _build_database_block(database_kwargs: dict) -> str:
    """
    Build the ``database(...)`` directive to prepend, merging caller overrides over
    ``DEFAULT_DATABASE_KWARGS``.

    Args:
        database_kwargs (dict, optional): Keyword arguments overriding/extending
                                          ``DEFAULT_DATABASE_KWARGS``.

    Every VALUE here goes through ``_render_literal``, but the KEYS were interpolated raw, and a key
    is not a value: it lands outside any quoting, in a position where a newline ends the statement.
    A key of ``"thermoLibraries=['x'])\\n__PWNED__ = <gadget>\\ndatabase(zzz"`` closed the generated
    ``database(...)`` call, inserted a top-level statement of the attacker's choosing, and reopened a
    ``database(`` whose keyword the rendered ``=[...],\\n)`` tail completed -- so the result was
    valid Python, passed the generated-file self-check (which only inspects the keyword values of
    ``database``/``explorer`` calls and therefore never looked at the injected statement), and
    executed when Arkane loaded the file. Verified before this guard existed.

    Keys are therefore required to be plain Python identifiers. That is not a heuristic about
    dangerous characters: an identifier cannot contain a newline, a quote, a parenthesis or an
    ``=``, so it cannot leave the keyword position it is rendered into, whatever it says.

    Args:
        database_kwargs (dict, optional): Keyword arguments overriding/extending
                                          ``DEFAULT_DATABASE_KWARGS``.

    Returns:
        str: The rendered ``database(...)`` directive, terminated with a blank line.

    Raises:
        ValueError: If any key is not a plain Python identifier, or is a dunder.
    """
    merged = dict(DEFAULT_DATABASE_KWARGS)
    if database_kwargs:
        merged.update(database_kwargs)
    for key in merged:
        if not isinstance(key, str) or not key.isidentifier() or keyword.iskeyword(key):
            raise ValueError(f"Refusing to render a database(...) keyword named {key!r}: a keyword name is "
                             f"interpolated into the generated Arkane input file OUTSIDE any quoting, so only "
                             f"a plain Python identifier can be guaranteed not to escape the keyword position "
                             f"it is written into and inject a statement of its own.")
        if key.startswith('__') and key.endswith('__'):
            raise ValueError(f"Refusing to render a database(...) keyword named {key!r}: a dunder is never an "
                             f"Arkane database(...) keyword, and naming one is a sign the value came from "
                             f"somewhere it should not have.")
        _validate_database_kwarg(key=key, value=merged[key])
    lines = [f'    {key}={_render_literal(key, value)},' for key, value in merged.items()]
    return 'database(\n' + '\n'.join(lines) + '\n)\n\n'


def _build_explorer_block(seed_species: tuple, bath_gas: dict, explore_tol, energy_tol, flux_tol,
                          maximum_radical_electrons) -> str:
    """
    Build the ``explorer(...)`` directive to append, matching the keyword order/shape of the real
    Arkane fixture (RMG-Py/examples/arkane/explorer/methoxy/input.py).

    Args:
        seed_species (tuple): The seed (source) species labels.
        bath_gas (dict): Bath gas species label -> mole fraction.
        explore_tol (float, optional): The energy tolerance for exploring new isomers/reactions.
        energy_tol (float, optional): The energy tolerance for including a well/transition state in
                                      the output network.
        flux_tol (float, optional): The flux tolerance for including a well/transition state in the
                                    output network.
        maximum_radical_electrons (int, optional): The maximum number of radical electrons allowed
                                                   in an explored species.

    Returns:
        str: The rendered ``explorer(...)`` directive.
    """
    values = {
        'source': list(seed_species),
        'explore_tol': explore_tol,
        'energy_tol': energy_tol,
        'flux_tol': flux_tol,
        'bathGas': dict(bath_gas),
        'maximumRadicalElectrons': maximum_radical_electrons,
    }
    parts = [f'{_EXPLORER_KWARG_NAMES[key]}={_render_literal(key, values[key])}' for key in _EXPLORER_KWARG_ORDER
            if values[key] is not None]
    return 'explorer(' + ', '.join(parts) + ')\n'


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
        # (e.g. ``\H``) that trigger a spurious SyntaxWarning under ``ast.parse``; see the matching
        # note in ``t3.pdep.parser`` and ``t3.pdep.hybrid``.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', SyntaxWarning)
            return ast.parse(text)
    except SyntaxError as e:
        raise ValueError(f"Could not parse '{path}' as Python: {e}")


def _source_snippet(node: ast.AST, text: str, limit: int = 80) -> str:
    """
    Render a short, single-line description of ``node`` for an error message.

    Args:
        node (ast.AST): The offending node.
        text (str): The full source text ``node`` came from.
        limit (int): The maximum snippet length before truncation.

    Returns:
        str: ``node``'s source text (whitespace-collapsed to one line), truncated to ``limit``
            characters, or a bracketed node-type name if the source segment cannot be recovered.
    """
    segment = ast.get_source_segment(text, node)
    if segment is None:
        segment = f'<{type(node).__name__}>'
    segment = ' '.join(segment.split())
    if len(segment) > limit:
        segment = segment[:limit] + '...'
    return segment


def _is_numeric_expression(node: ast.AST) -> bool:
    """
    Whether ``node`` is an arithmetic expression over numeric literals only.

    Called on operands of a ``+``/``-``/``*`` ``ast.BinOp`` to keep the multiplication that a real
    input file uses (``0.5 * 4.184``) while refusing the repetition that shares its syntax
    (``'x' * 10 ** 9``, ``[0] * 10 ** 8``). ``bool`` is deliberately NOT numeric here even though
    Python treats it as an ``int`` subclass: arithmetic on ``True`` is never what an input file
    means, so refusing it costs nothing and keeps the rule easy to state.

    Args:
        node (ast.AST): The operand to classify. Assumed already whitelisted by
                       ``_validate_source_expression``, so only the node kinds it permits are
                       considered here.

    Returns:
        bool: ``True`` if the operand is numeric-literal arithmetic, ``False`` otherwise.
    """
    if isinstance(node, ast.Constant):
        return isinstance(node.value, (int, float, complex)) and not isinstance(node.value, bool)
    if isinstance(node, ast.UnaryOp):
        return _is_numeric_expression(node.operand)
    if isinstance(node, ast.BinOp):
        return _is_numeric_expression(node.left) and _is_numeric_expression(node.right)
    return False


def _validate_source_expression(node: ast.AST, text: str, source_path: str) -> None:
    """
    Recursively refuse anything inside a to-be-spliced top-level statement that is not a literal,
    numeric arithmetic, or a call to an inert Arkane value constructor.

    This is the core of the fix for the remote-code-execution hole: ``write_arkane_explorer_input_file``
    splices the untrusted source's text verbatim into a NEW file that Arkane later ``exec``s, so every
    expression that survives into that file must be provably inert. Permitted node kinds: ``ast.Constant``,
    ``ast.List``/``ast.Tuple``/``ast.Set``, ``ast.Dict``, ``ast.UnaryOp`` (``+``/``-`` only), ``ast.BinOp``
    (``+``/``-``/``*`` only, and only over numeric operands -- see ``_is_numeric_expression``), and
    ``ast.Call`` whose ``func`` is a bare ``ast.Name`` in
    ``_ARKANE_VALUE_NAMES`` (recursing into its ``args``/``keywords``; ``**kwargs`` unpacking -- a
    ``keyword`` with ``arg is None`` -- is refused). Everything else (``ast.Attribute``, ``ast.Subscript``,
    a loaded ``ast.Name``, comprehensions, ``ast.Lambda``, ``ast.IfExp``, ``ast.Compare``, ``ast.BoolOp``,
    ``ast.Starred``, ``ast.JoinedStr``, ``ast.Await``, ``ast.NamedExpr``, etc.) is refused.

    Args:
        node (ast.AST): The expression node to validate.
        text (str): The full source text (used only to render error-message snippets).
        source_path (str): The source path (used only in error messages).

    Raises:
        ValueError: If ``node`` (or anything nested inside it) is not one of the permitted forms.
    """
    if isinstance(node, ast.Constant):
        return
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        for elt in node.elts:
            _validate_source_expression(elt, text, source_path)
        return
    if isinstance(node, ast.Dict):
        for key, value in zip(node.keys, node.values):
            if key is not None:
                _validate_source_expression(key, text, source_path)
            _validate_source_expression(value, text, source_path)
        return
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        _validate_source_expression(node.operand, text, source_path)
        return
    if isinstance(node, ast.BinOp):
        if not isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
            raise ValueError(
                f"Refusing to use '{source_path}' as an Arkane explorer/network source: line {node.lineno} "
                f"uses the arithmetic operator {type(node.op).__name__} ({_source_snippet(node, text)!r}). "
                f"Only +, -, * and / are permitted in a spliced expression. '**' in particular is refused "
                f"because it turns a short literal into an arbitrarily large one that Arkane evaluates at "
                f"load time -- 10 ** 9 is three characters of source and a memory/CPU bomb -- and no real "
                f"Arkane or RMG input file uses it (measured across all 100 example input files in RMG-Py: "
                f"only '-' and '*' appear at all).")
        for operand in (node.left, node.right):
            _validate_source_expression(operand, text, source_path)
            if not _is_numeric_expression(operand):
                raise ValueError(
                    f"Refusing to use '{source_path}' as an Arkane explorer/network source: line "
                    f"{node.lineno} applies arithmetic to a non-numeric operand "
                    f"({_source_snippet(node, text)!r}). Arithmetic on a string or a container is repetition, "
                    f"not arithmetic: 'x' * 1000000000 and [0] * 100000000 are short to write and are "
                    f"evaluated -- at full size -- when Arkane execs the file this text is spliced into. Only "
                    f"numeric operands are permitted; real input files only ever do arithmetic on numbers.")
        return
    if isinstance(node, ast.Call):
        func = node.func
        if isinstance(func, ast.Name) and func.id in _ARKANE_SIDE_EFFECT_NAMES:
            raise ValueError(
                f"Refusing to use '{source_path}' as an Arkane explorer/network source: line "
                f"{node.lineno} calls '{func.id}' from inside an expression "
                f"({_source_snippet(node, text)!r}). Arkane defines '{func.id}', and that is exactly the "
                f"problem: it is one of Arkane's job/state directives, whose only effect is to add a job "
                f"to -- or reconfigure -- the run. This source's text is spliced verbatim into a NEW file "
                f"Arkane will exec, so a job directive reached from inside a spliced value would run in, "
                f"and change, the explorer job this module is generating. Only the inert value "
                f"constructors in _ARKANE_VALUE_NAMES may appear inside a spliced expression.")
        if isinstance(func, ast.Name) and func.id in _ARKANE_VALUE_NAMES:
            for arg in node.args:
                _validate_source_expression(arg, text, source_path)
            for kw in node.keywords:
                if kw.arg is None:
                    raise ValueError(
                        f"Refusing to use '{source_path}' as an Arkane explorer/network source: line "
                        f"{node.lineno} uses '**...' keyword-argument unpacking in a call "
                        f"({_source_snippet(node, text)!r}). This source's text is spliced verbatim into a "
                        f"NEW file Arkane will exec, so an unpacked mapping (which could smuggle an "
                        f"arbitrary keyword past this whitelist) is refused.")
                _validate_source_expression(kw.value, text, source_path)
            return
        if isinstance(func, ast.Name):
            raise ValueError(
                f"Refusing to use '{source_path}' as an Arkane explorer/network source: line {node.lineno} "
                f"calls '{func.id}' ({_source_snippet(node, text)!r}), which is not a recognized Arkane value "
                f"constructor. This source's text is spliced verbatim into a NEW file Arkane will exec, so "
                f"only calls to Arkane's own inert value constructors are permitted. If '{func.id}' IS a "
                f"legitimate Arkane value constructor, add it to _ARKANE_VALUE_NAMES.")
        raise ValueError(
            f"Refusing to use '{source_path}' as an Arkane explorer/network source: line {node.lineno} calls "
            f"a non-plain-name callee ({_source_snippet(node, text)!r}), e.g. an attribute call such as "
            f"'foo.species(...)'. This source's text is spliced verbatim into a NEW file Arkane will exec, "
            f"so only bare calls to Arkane's own DSL functions (by name) are permitted -- an attribute "
            f"access here is never treated as a recognized DSL callee, however it is named.")
    raise ValueError(
        f"Refusing to use '{source_path}' as an Arkane explorer/network source: line {node.lineno} contains a "
        f"{type(node).__name__} expression ({_source_snippet(node, text)!r}), which is not one of the "
        f"literal, arithmetic, or Arkane-DSL-call forms permitted inside a spliced statement. This source's "
        f"text is spliced verbatim into a NEW file Arkane will exec, so anything else here is refused.")


_GENERATED_CALL_NAMES = tuple(_TOP_LEVEL_CALL_NAMES) + ('database', 'explorer', 'kinetics')


def _validate_generated_statements(tree: ast.Module, text: str, dest_path: str) -> None:
    """
    Self-check the top-level statement list of the file this module is about to WRITE.

    ``_validate_source_statements`` judges the untrusted input; this judges this module's own output,
    which is not the same question. Everything here got in either by being spliced from an
    already-validated source or by being generated below, so a statement that is neither is evidence
    of a bug in the splicing or the rendering -- and that bug is how the ``database_kwargs`` key
    injection worked: it manufactured a top-level statement that no validator was looking at, because
    the only self-check inspected the keyword values of two named calls.

    The permitted statements are the source's shapes plus the three directives this module generates
    (``database``, ``explorer``, ``kinetics``). Refusing here raises ``RuntimeError``, not
    ``ValueError``: a rejected source is a bad input, but a rejected OUTPUT is this module's own
    defect, and the two should not arrive at a caller looking alike.

    Args:
        tree (ast.Module): The parsed generated text.
        text (str): The generated text (used only to render error-message snippets).
        dest_path (str): The path about to be written (used only in error messages).

    Raises:
        RuntimeError: If the generated text contains a top-level statement this module never meant
                     to emit.
    """
    for node in tree.body:
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
            continue
        if (isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                and _get_call_name(node.value) in _GENERATED_CALL_NAMES):
            continue
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                ast.literal_eval(node.value)
            except (ValueError, TypeError, SyntaxError, MemoryError, RecursionError):
                pass
            else:
                continue
        raise RuntimeError(
            f"Refusing to write '{dest_path}': the generated explorer input file text failed its own "
            f"self-check. Line {node.lineno} is a top-level {type(node).__name__} statement "
            f"({_source_snippet(node, text)!r}) that this module never emits and no validated source "
            f"statement can become -- only calls to {sorted(_GENERATED_CALL_NAMES)}, bare literals and "
            f"literal assignments belong here. Something spliced or rendered a statement into this file, "
            f"which is a defect in this module (or an injection through one of its arguments), not a "
            f"problem with the source file.")


def _validate_source_statements(tree: ast.Module, text: str, source_path: str) -> None:
    """
    Refuse anything at module top level that is not a recognized Arkane DSL call or a single-target
    literal assignment, closing the remote-code-execution hole this module previously had: any
    top-level statement it did not otherwise recognize (an ``import``, a ``def``, an assignment
    computed from a call it never inspected, ...) was silently spliced verbatim into a NEW file that
    Arkane later ``exec``s.

    Only three top-level statement shapes are permitted:
      - ``ast.Expr`` whose ``.value`` is an ``ast.Call`` to one of ``_TOP_LEVEL_CALL_NAMES`` -- the
        five network-describing directives this module actually understands and rewrites. Its
        arguments are validated recursively by ``_validate_source_expression``, but the callee itself
        is not routed through that function, because at top level these five names are the whole
        point of a network source while inside an expression they would only ever be a smuggled side
        effect. Arkane's other directives (``kinetics``, ``thermo``, ``explorer``, ``database``, ...)
        are refused here: each one would add a job to, or reconfigure, the explorer run this module
        is generating, and a spliced ``database(...)`` lands after the generated one and overrides it.
      - ``ast.Assign`` with exactly one target, that target a bare ``ast.Name``, whose value is a true
        literal -- ``ast.literal_eval``-able (e.g. ``title = '...'``, ``useHinderedRotors = True``).
        Anything computed is refused, which is what the surrounding comments always claimed and what
        the recursive expression whitelist did not actually enforce.
      - ``ast.Expr`` whose ``.value`` is an ``ast.Constant`` -- a module docstring or other bare
        literal. Evaluating a constant and discarding it cannot do anything, so refusing it would buy
        nothing.

    Args:
        tree (ast.Module): The parsed source tree.
        text (str): The full source text (used only to render error-message snippets).
        source_path (str): The source path (used only in error messages).

    Raises:
        ValueError: If any top-level statement (or anything nested inside a permitted one) is not one
                   of the three permitted shapes above.
    """
    for node in tree.body:
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
            continue
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            if not isinstance(node.value.func, ast.Name):
                # Delegated so an attribute callee ('foo.species(...)') is diagnosed as the attribute
                # access it is, rather than as a merely unrecognized directive name. Both refuse; only
                # one of them tells the reader that no amount of renaming will help.
                _validate_source_expression(node.value, text, source_path)
            call_name = _get_call_name(node.value)
            if call_name not in _TOP_LEVEL_CALL_NAMES:
                named = (f"Arkane defines '{call_name}', and that is exactly the problem: calling it adds a "
                         f"job to -- or reconfigures -- the run, and a spliced 'database(...)' lands after "
                         f"the generated one and silently overrides it. "
                         if call_name in _ARKANE_SIDE_EFFECT_NAMES else '')
                raise ValueError(
                    f"Refusing to use '{source_path}' as an Arkane explorer/network source: line "
                    f"{node.lineno} is a top-level call that is not one of the network-describing "
                    f"directives this module understands ({sorted(_TOP_LEVEL_CALL_NAMES)}): "
                    f"{_source_snippet(node, text)!r}. {named}This source's text is spliced verbatim into "
                    f"a NEW file Arkane will exec, so a top-level call it does not understand is refused "
                    f"rather than carried along.")
            for arg in node.value.args:
                _validate_source_expression(arg, text, source_path)
            for kw in node.value.keywords:
                if kw.arg is None:
                    raise ValueError(
                        f"Refusing to use '{source_path}' as an Arkane explorer/network source: line "
                        f"{node.lineno} uses '**...' keyword-argument unpacking in a call "
                        f"({_source_snippet(node, text)!r}). This source's text is spliced verbatim into a "
                        f"NEW file Arkane will exec, so an unpacked mapping (which could smuggle an "
                        f"arbitrary keyword past this whitelist) is refused.")
                _validate_source_expression(kw.value, text, source_path)
            continue
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            target = node.targets[0].id
            # A literal value is inert; a literal assignment is not, because the NAME it binds may be
            # one Arkane put in the namespace. 'explorer = None' is spliced above the generated
            # explorer(...) directive and turns it into a TypeError at load time; 'species = None'
            # does the same to every species block after it. The value being harmless says nothing
            # about the binding.
            if target in _ARKANE_VALUE_NAMES or target in _ARKANE_SIDE_EFFECT_NAMES or target == 'range':
                raise ValueError(
                    f"Refusing to use '{source_path}' as an Arkane explorer/network source: line "
                    f"{node.lineno} assigns to '{target}' ({_source_snippet(node, text)!r}), which is a name "
                    f"Arkane itself defines in the namespace it loads an input file in. Rebinding it shadows "
                    f"the real one for every statement that follows -- including the directives this module "
                    f"generates and appends -- so a source may not assign to it, however harmless the value.")
            # ``modelChemistry`` is the one directive whose real ARC/hybrid value is a bare
            # ``LevelOfTheory(...)``/``CompositeLevelOfTheory(...)`` call (which Arkane execs at load
            # time), not an ``ast.literal_eval``-able literal. Rather than grow a second model-chemistry
            # allowlist here, route it through the SAME structural checker T3 uses when it EMITS this
            # directive (``t3.pdep.hybrid._validate_model_chemistry_expression``); the AST is what
            # decides exec semantics and the verbatim splice reproduces it, so validating the
            # ``ast.unparse``-round-tripped node -- which that string-taking checker re-parses -- is
            # equivalent to validating the spliced text. The gate is deliberately on a genuine call
            # node to one of the two known names: the checker accepts any NON-call string as a plain
            # label (injection-char check only), so a computed ``modelChemistry`` value must keep
            # falling through to the refusal below rather than being handed over as a bare label.
            if (target == 'modelChemistry' and isinstance(node.value, ast.Call)
                    and isinstance(node.value.func, ast.Name)
                    and node.value.func.id in _MODEL_CHEMISTRY_CALL_NAMES):
                try:
                    _validate_model_chemistry_expression(target, ast.unparse(node.value))
                except ValueError as e:
                    raise ValueError(
                        f"Refusing to use '{source_path}' as an Arkane explorer/network source: line "
                        f"{node.lineno} assigns a 'modelChemistry' value that fails structural validation "
                        f"({_source_snippet(node, text)!r}): {e}. This source's text is spliced verbatim into "
                        f"a NEW file Arkane will exec, so a malformed model-chemistry call is refused.") from e
                continue
            try:
                ast.literal_eval(node.value)
            except (ValueError, TypeError, SyntaxError, MemoryError, RecursionError) as e:
                raise ValueError(
                    f"Refusing to use '{source_path}' as an Arkane explorer/network source: line "
                    f"{node.lineno} assigns something that is not a literal "
                    f"({_source_snippet(node, text)!r}): {e}. This source's text is spliced verbatim into a "
                    f"NEW file Arkane will exec, so a top-level assignment must be a plain literal value; "
                    f"anything computed at load time is refused.") from e
            continue
        raise ValueError(
            f"Refusing to use '{source_path}' as an Arkane explorer/network source: line {node.lineno} is a "
            f"top-level {type(node).__name__} statement ({_source_snippet(node, text)!r}), which is none of "
            f"a call to one of {sorted(_TOP_LEVEL_CALL_NAMES)}, a single-target literal assignment (e.g. "
            f"'title = ...'), or a bare literal. This source's text is spliced verbatim into a NEW file that "
            f"Arkane will exec, so only those statement shapes are permitted here; everything else -- "
            f"including 'import', 'def', 'class', 'if', 'for' and multi-target assignments -- is refused.")


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
    line, not a character index -- applying it directly as a character index would silently
    misalign every position on a line that has any non-ASCII (multi-byte-in-UTF-8) character
    before it. This re-encodes the target line to UTF-8, slices it to the byte offset, and decodes
    back to count how many characters that many bytes actually span.

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
    Get the callee name of an ``ast.Call`` node whose ``func`` is a bare ``ast.Name`` (e.g.,
    ``'species'`` for ``species(...)``).

    Deliberately returns ``None`` (rather than the trailing attribute) for anything else, most
    importantly an ``ast.Attribute`` callee such as ``foo.species(...)``: treating that as
    ``'species'`` would let an attribute call impersonate a recognized top-level Arkane DSL call.

    Args:
        call (ast.Call): The call node.

    Returns:
        str: The callee name, or ``None`` if ``call.func`` is not a bare ``ast.Name``.
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


# The 'multiplicity N' header a species' adjacencyList string carries, e.g. the first line of
# 'multiplicity 3\n1 O u2 p2 c0\n'. Anchored to a line start so it never matches the word inside a
# comment or a species label.
_ADJ_MULTIPLICITY_RE = re.compile(r'^\s*multiplicity\s+(\d+)\s*$', re.MULTILINE)


def _adjacency_list_multiplicity(structure_node) -> int | None:
    """
    The explicit ``multiplicity`` a species' ``adjacencyList(...)`` structure declares.

    Returns ``None`` -- meaning "no correction to make here" -- when the structure is not an
    ``adjacencyList(...)`` call (e.g. ``SMILES(...)``, whose multiplicity Arkane derives itself),
    when its argument is not a plain string literal, or when the adjacency string carries no
    explicit ``multiplicity`` header (a closed-shell species omits it, defaulting to 1, which needs
    no correction).

    Args:
        structure_node: The AST node of the species' ``structure`` keyword value, or ``None``.

    Returns:
        Optional[int]: The declared multiplicity, or ``None`` if there is nothing to correct.
    """
    if not isinstance(structure_node, ast.Call) or _get_call_name(structure_node) != 'adjacencyList':
        return None
    if not structure_node.args:
        return None
    arg = structure_node.args[0]
    adj_text = arg.value if isinstance(arg, ast.Constant) and isinstance(arg.value, str) else None
    if adj_text is None:
        return None
    match = _ADJ_MULTIPLICITY_RE.search(adj_text)
    return int(match.group(1)) if match else None
