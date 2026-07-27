"""
t3 pdep parser module

An ``ast``-based parser for RMG-generated pressure-dependent network files (e.g.,
``RMG/pdep/network4_2.py``). These files are executable Arkane DSL, but they must
NEVER be executed (via ``exec``, ``eval``, or ``import``): they may originate from
untrusted or malformed RMG runs, and executing them is unnecessary since only their
declarative structure (species, transition states, and path reactions) is needed here.
This module therefore parses the file into a syntax tree with ``ast.parse`` only.
"""

import ast
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

RECOGNIZED_TOP_LEVEL_CALLS = {'species', 'transitionState', 'reaction', 'network', 'pressureDependence'}

# ``pdepreaction(...)`` is deliberately NOT part of ``RECOGNIZED_TOP_LEVEL_CALLS`` above: it is
# not part of the Arkane *input* DSL that RMG pdep network files use. It is a write-only call
# that Arkane itself emits into a pdep job's ``output.py`` to report a solved kinetics result.
ARKANE_PDEP_OUTPUT_TOP_LEVEL_CALL = 'pdepreaction'

# Keywords of a kinetics call that describe the validity bounds or metadata of a rate expression
# rather than the rate itself. These are always finite in Arkane's output (they echo the job's
# T/P ranges), so they must never be allowed to satisfy a "the kinetics parameters are finite"
# check on their own: a ``Chebyshev(coeffs=[], Tmin=..., Tmax=..., Pmin=..., Pmax=...)`` carries
# NO rate information despite having four finite numeric leaves. ``T0`` (the modified-Arrhenius
# reference temperature, a fixed 1 K in RMG/Arkane output) is a reference quantity like ``Tref``
# and ``Pref``, not a fitted rate coefficient. The keys observed in real Arkane pdep output
# fixtures (``tests/data/pdep_me/``) are Tmin/Tmax/Pmin/Pmax/kunits.
KINETICS_BOUNDS_AND_METADATA_KEYS = ('Tmin', 'Tmax', 'Pmin', 'Pmax', 'Tref', 'Pref', 'T0',
                                     'kunits', 'comment', 'efficiencies')

# Kinetics call names that may legitimately appear nested inside a ``kinetics=`` keyword of an
# Arkane pdep output entry (e.g., ``PDepArrhenius(arrhenius=[Arrhenius(...), ...])``). Payload
# extraction recurses into these; any OTHER call in payload position (``float('nan')``,
# ``array(...)``) is unverifiable and is surfaced as a ``None`` leaf so the ME-success gate
# rejects it rather than silently skipping it.
NESTED_KINETICS_CALL_NAMES = ('Arrhenius', 'MultiArrhenius', 'PDepArrhenius', 'MultiPDepArrhenius',
                              'Chebyshev', 'ThirdBody', 'Lindemann', 'Troe')


@dataclass(frozen=True)
class PDepPathReaction:
    """
    A single ``reaction(...)`` entry (a path reaction) parsed from an RMG pdep network file.

    Args:
        label (str): The reaction label, e.g. ``'reaction2'``.
        reactants (tuple): The reactant species labels.
        products (tuple): The product species labels.
        transition_state (str, optional): The associated transition state label, if any.
        kinetics_type (str, optional): The kinetics callee name, e.g. ``'Arrhenius'``, if any.
        kinetics_comment (str): The kinetics ``comment`` keyword text, or ``''`` if absent.
    """
    label: str
    reactants: tuple
    products: tuple
    transition_state: str | None
    kinetics_type: str | None
    kinetics_comment: str


@dataclass(frozen=True)
class PDepNetwork:
    """
    A parsed RMG pdep network file.

    Args:
        network_id (str): The file stem, e.g. ``'network4_2'`` (used as a join key downstream).
        path (str): The absolute path that was parsed.
        label (str, optional): The ``network(...)`` label, if any.
        species_labels (tuple): All ``species(...)`` labels, in file order.
        transition_state_labels (tuple): All ``transitionState(...)`` labels, in file order.
        path_reactions (tuple): The parsed ``PDepPathReaction`` entries, in file order.
        isomers (tuple): The ``network(...)`` isomer labels.
        reactant_channels (tuple): The ``network(...)`` reactant channels, each a tuple of str.
        product_channels (tuple): The ``network(...)`` product channels, each a tuple of str.
    """
    network_id: str
    path: str
    label: str | None
    species_labels: tuple
    transition_state_labels: tuple
    path_reactions: tuple
    isomers: tuple
    reactant_channels: tuple
    product_channels: tuple

    def path_reactions_by_ts(self) -> dict:
        """
        Map each transition state label to the path reactions that use it.

        Several path reactions can legitimately share one transition state label, so the
        value is always a tuple and callers must handle that case explicitly. Path reactions
        with ``transition_state`` set to ``None`` are excluded from this map.

        Returns:
            dict: Keys are transition state labels (str), values are tuples of ``PDepPathReaction``.
        """
        by_ts = dict()
        for path_reaction in self.path_reactions:
            if path_reaction.transition_state is None:
                continue
            by_ts.setdefault(path_reaction.transition_state, list()).append(path_reaction)
        return {ts_label: tuple(reactions) for ts_label, reactions in by_ts.items()}


@dataclass(frozen=True)
class PDepArkaneReaction:
    """
    A single ``pdepreaction(...)`` entry parsed from an Arkane pdep job's ``output.py``.

    Args:
        reactants (tuple): The reactant species labels.
        products (tuple): The product species labels.
        kinetics_type (str, optional): The kinetics callee name, e.g. ``'Chebyshev'``, if any.
        kinetics_params (dict): The kinetics call's keyword arguments, keyword name -> parsed
            literal value (e.g., ``'coeffs'`` -> a nested list of numbers/``None``, ``'Tmin'`` ->
            a ``(300, 'K')`` tuple). Keywords whose value could not be evaluated as a literal are
            omitted.
        numeric_values (tuple): Every numeric leaf found anywhere inside ``kinetics_params``
            (recursing into nested lists/tuples), flattened into one tuple, in traversal order.
            ``None`` entries are preserved (not skipped or coerced) since detecting them is the
            whole point: a rejected chemically-significant-eigenvalues solve still writes a
            syntactically valid ``Chebyshev(coeffs=[[None, None, ...], ...])`` block. NOTE: this
            deliberately mixes the rate payload with the (always finite) T/P bounds/metadata
            numbers; success gating must use ``rate_payload_numeric_values`` instead.
        rate_payload_numeric_values (tuple): The numeric (or ``None``) leaves of the RATE PAYLOAD
            only: every kinetics keyword that is not in ``KINETICS_BOUNDS_AND_METADATA_KEYS``
            (e.g. ``coeffs`` for Chebyshev, ``A``/``n``/``Ea`` for Arrhenius, the nested
            ``arrhenius=[...]`` list for PDepArrhenius/MultiArrhenius). Extraction recurses into
            nested kinetics calls; anything unverifiable in payload position (a bare name such as
            ``nan``, an unrecognized call such as ``float('nan')`` or ``array(...)``) surfaces as
            a ``None`` leaf rather than being skipped.
        missing_kinetics_keys (tuple): The kinetics keywords that were omitted from
            ``kinetics_params`` because their values could not be evaluated as literals, so an
            omitted ``coeffs`` is distinguishable from a legitimately absent one.
    """
    reactants: tuple
    products: tuple
    kinetics_type: str | None
    kinetics_params: dict
    numeric_values: tuple
    rate_payload_numeric_values: tuple = tuple()
    missing_kinetics_keys: tuple = tuple()


def parse_arkane_pdep_output_file(path: str) -> tuple:
    """
    Parse an Arkane pdep job's ``output.py`` file from disk.

    Args:
        path (str): The path to the Arkane pdep ``output.py`` file.

    Raises:
        ValueError: If the file cannot be parsed as Python.

    Returns:
        tuple: The parsed ``PDepArkaneReaction`` entries, in file order. Empty if the file
        contains no ``pdepreaction(...)`` calls (including an empty or comment-only file);
        unlike ``parse_pdep_network_file``, this deliberately does NOT raise in that case, since
        emptiness is itself an informative outcome for an ME-success gate built on top of this.
    """
    with open(path, 'r') as f:
        text = f.read()
    return parse_arkane_pdep_output_text(text=text, path=path)


def parse_arkane_pdep_output_text(text: str, path: str = '') -> tuple:
    """
    Parse Arkane pdep ``output.py`` text into a tuple of ``PDepArkaneReaction``.

    Args:
        text (str): The Arkane pdep ``output.py`` file content.
        path (str, optional): The path the text was read from, if any (used only in error
            messages).

    Raises:
        ValueError: If the text cannot be parsed as Python.

    Returns:
        tuple: The parsed ``PDepArkaneReaction`` entries, in file order. Empty if no
        ``pdepreaction(...)`` calls are present; this is not an error (see
        ``parse_arkane_pdep_output_file``).
    """
    try:
        # Arkane output files can contain the same RMG-generated comments (with invalid escape
        # sequences) that trigger a spurious SyntaxWarning under ``ast.parse``; see the matching
        # note in ``parse_pdep_network_text`` above.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', SyntaxWarning)
            tree = ast.parse(text)
    except SyntaxError as e:
        raise ValueError(f"Could not parse the Arkane pdep output file at '{path}' as Python: {e}")

    reactions = []
    for node in tree.body:
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
            continue
        call = node.value
        if _get_call_name(call) != ARKANE_PDEP_OUTPUT_TOP_LEVEL_CALL:
            continue
        reactions.append(_parse_arkane_pdep_reaction(_call_keywords(call)))
    return tuple(reactions)


def _parse_arkane_pdep_reaction(kwargs: dict) -> PDepArkaneReaction:
    """
    Build a ``PDepArkaneReaction`` from a ``pdepreaction(...)`` call's keyword arguments.

    Args:
        kwargs (dict): A mapping of keyword name to its (still-AST) value node.

    Returns:
        PDepArkaneReaction: The parsed Arkane pdep reaction.
    """
    reactants = tuple(_literal_or_none(kwargs.get('reactants')) or ())
    products = tuple(_literal_or_none(kwargs.get('products')) or ())

    kinetics_type = None
    kinetics_params = dict()
    numeric_values = []
    rate_payload_numeric_values = []
    missing_kinetics_keys = []
    kinetics_node = kwargs.get('kinetics')
    if isinstance(kinetics_node, ast.Call):
        kinetics_type = _get_call_name(kinetics_node)
        for key, value_node in _call_keywords(kinetics_node).items():
            value = _literal_or_none(value_node)
            if value is None and not _is_literal_none(value_node):
                # The value could not be evaluated as a literal at all (e.g., itself a nested
                # Call or a bare name such as ``nan``); omit it from ``kinetics_params`` rather
                # than silently recording a fabricated ``None``, but RECORD the omission so an
                # unparseable ``coeffs`` is distinguishable from a legitimately absent one.
                missing_kinetics_keys.append(key)
            else:
                kinetics_params[key] = value
                numeric_values.extend(_extract_numeric_leaves(value))
            if key not in KINETICS_BOUNDS_AND_METADATA_KEYS:
                # Rate payload (coeffs / A / n / Ea / arrhenius=[...] / ...): extract its numeric
                # leaves from the AST node itself, so nested kinetics calls contribute their
                # payload and unverifiable nodes surface as None leaves.
                rate_payload_numeric_values.extend(_extract_payload_numeric_leaves(value_node))

    return PDepArkaneReaction(
        reactants=reactants,
        products=products,
        kinetics_type=kinetics_type,
        kinetics_params=kinetics_params,
        numeric_values=tuple(numeric_values),
        rate_payload_numeric_values=tuple(rate_payload_numeric_values),
        missing_kinetics_keys=tuple(missing_kinetics_keys),
    )


def _is_literal_none(node) -> bool:
    """
    Check whether an AST node is literally the ``None`` constant (as opposed to a node that
    simply failed to evaluate as a literal, for which ``_literal_or_none`` also returns ``None``).

    Args:
        node: An AST node (or ``None``).

    Returns:
        bool: True if ``node`` is an ``ast.Constant`` whose value is ``None``.
    """
    return isinstance(node, ast.Constant) and node.value is None


def _extract_payload_numeric_leaves(node) -> list:
    """
    Recursively collect the numeric (or ``None``) leaves of a rate-payload AST node.

    Unlike ``_extract_numeric_leaves`` (which walks an already literal-evaluated value), this
    walks the AST node directly, so it can recurse into nested kinetics calls such as
    ``PDepArrhenius(arrhenius=[Arrhenius(...), ...])`` whose values ``ast.literal_eval`` cannot
    handle. Within a recognized nested kinetics call, keywords listed in
    ``KINETICS_BOUNDS_AND_METADATA_KEYS`` are skipped (they are bounds/metadata, not payload).
    Any node in payload position that cannot be verified as a finite literal -- a bare name such
    as ``nan``, an unrecognized call such as ``float('nan')`` or ``array(...)`` -- yields a
    ``None`` leaf, so the ME-success gate rejects it instead of silently skipping it. String
    constants (unit labels inside quantity tuples like ``(1e13, 's^-1')``) are skipped.

    Args:
        node: An AST node from a rate-payload keyword's value.

    Returns:
        list: The numeric/``None`` leaves found, in traversal order.
    """
    leaves = []
    if isinstance(node, ast.Constant):
        if node.value is None:
            leaves.append(None)
        elif isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
            leaves.append(node.value)
        # Strings and other constants (unit labels, comments) carry no rate information: skip.
    elif isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        for element in node.elts:
            leaves.extend(_extract_payload_numeric_leaves(element))
    elif isinstance(node, ast.Dict):
        for value_node in node.values:
            leaves.extend(_extract_payload_numeric_leaves(value_node))
    elif isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        inner = _extract_payload_numeric_leaves(node.operand)
        if isinstance(node.op, ast.USub):
            inner = [-leaf if leaf is not None else None for leaf in inner]
        leaves.extend(inner)
    elif isinstance(node, ast.Call) and _get_call_name(node) in NESTED_KINETICS_CALL_NAMES:
        for key, value_node in _call_keywords(node).items():
            if key in KINETICS_BOUNDS_AND_METADATA_KEYS:
                continue
            leaves.extend(_extract_payload_numeric_leaves(value_node))
    else:
        # Unverifiable payload content (bare Name like ``nan``, unrecognized call like
        # ``float('nan')`` or ``array(...)``, arbitrary expression): surface as None so the
        # gate treats it as non-finite rather than pretending it is absent.
        leaves.append(None)
    return leaves


def _extract_numeric_leaves(value) -> list:
    """
    Recursively collect every numeric (or ``None``) leaf out of a parsed literal value.

    Handles nested lists/tuples such as ``Chebyshev(coeffs=[[...], [...]])`` and quantity
    tuples such as ``Tmin=(300, 'K')``: numbers are collected, non-numeric leaves (e.g., unit
    strings like ``'K'`` or ``'s^-1'``) are silently skipped, and ``None`` leaves are preserved
    as ``None`` rather than skipped or coerced, since detecting them is the whole point of this
    module's caller (the ME-success gate).

    Args:
        value: A parsed literal value (number, ``None``, str, list, tuple, dict, or a nesting
            of these).

    Returns:
        list: The numeric/``None`` leaves found, in traversal order.
    """
    leaves = []
    if isinstance(value, (list, tuple)):
        for item in value:
            leaves.extend(_extract_numeric_leaves(item))
    elif isinstance(value, dict):
        for item in value.values():
            leaves.extend(_extract_numeric_leaves(item))
    elif value is None:
        leaves.append(None)
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        leaves.append(value)
    return leaves


def parse_pdep_network_file(path: str) -> PDepNetwork:
    """
    Parse an RMG pdep network file from disk.

    Args:
        path (str): The path to the RMG pdep network file (e.g., ``.../pdep/network4_2.py``).

    Raises:
        ValueError: If the file cannot be parsed as Python, or contains no path reactions.

    Returns:
        PDepNetwork: The parsed network.
    """
    with open(path, 'r') as f:
        text = f.read()
    network_id = Path(path).stem
    return parse_pdep_network_text(text=text, network_id=network_id, path=path)


def parse_pdep_network_text(text: str, network_id: str, path: str = '') -> PDepNetwork:
    """
    Parse RMG pdep network file text into a ``PDepNetwork``.

    Args:
        text (str): The RMG pdep network file content.
        network_id (str): An identifier for this network (e.g., the file stem).
        path (str, optional): The path the text was read from, if any (stored for reference).

    Raises:
        ValueError: If the text cannot be parsed as Python, or contains no path reactions.

    Returns:
        PDepNetwork: The parsed network.
    """
    try:
        # Some RMG-generated comments contain invalid escape sequences (e.g., ``\H``) that
        # trigger a SyntaxWarning under ``ast.parse`` even though the file is otherwise valid
        # Python. This narrowly suppresses only that warning around the parse call itself.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', SyntaxWarning)
            tree = ast.parse(text)
    except SyntaxError as e:
        raise ValueError(f"Could not parse the pdep network file at '{path}' as Python: {e}")

    species_labels = list()
    transition_state_labels = list()
    path_reactions = list()
    network_label = None
    isomers = tuple()
    reactant_channels = tuple()
    product_channels = tuple()

    for node in tree.body:
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
            continue
        call = node.value
        call_name = _get_call_name(call)
        if call_name not in RECOGNIZED_TOP_LEVEL_CALLS:
            continue
        kwargs = _call_keywords(call, path=path, call_name=call_name)

        if call_name == 'species':
            label = _literal_or_none(kwargs.get('label'))
            if label is not None:
                species_labels.append(label)

        elif call_name == 'transitionState':
            label = _literal_or_none(kwargs.get('label'))
            if label is not None:
                transition_state_labels.append(label)

        elif call_name == 'reaction':
            path_reactions.append(_parse_reaction(kwargs))

        elif call_name == 'network':
            network_label = _literal_or_none(kwargs.get('label'))
            isomers = tuple(_literal_or_none(kwargs.get('isomers')) or ())
            reactant_channels = tuple(tuple(channel) for channel in (_literal_or_none(kwargs.get('reactants')) or ()))
            product_channels = tuple(tuple(channel) for channel in (_literal_or_none(kwargs.get('products')) or ()))

        # ``pressureDependence(...)`` data is not part of the public API here; it is
        # recognized (so it is not silently ignored as "unrecognized") but otherwise unused.

    if not path_reactions:
        raise ValueError(f"The pdep network file at '{path}' contains no reaction(...) entries; nothing to parse.")

    return PDepNetwork(
        network_id=network_id,
        path=path,
        label=network_label,
        species_labels=tuple(species_labels),
        transition_state_labels=tuple(transition_state_labels),
        path_reactions=tuple(path_reactions),
        isomers=isomers,
        reactant_channels=reactant_channels,
        product_channels=product_channels,
    )


def _parse_reaction(kwargs: dict) -> PDepPathReaction:
    """
    Build a ``PDepPathReaction`` from a ``reaction(...)`` call's keyword arguments.

    Args:
        kwargs (dict): A mapping of keyword name to its (still-AST) value node.

    Returns:
        PDepPathReaction: The parsed path reaction.
    """
    label = _literal_or_none(kwargs.get('label'))
    reactants = tuple(_literal_or_none(kwargs.get('reactants')) or ())
    products = tuple(_literal_or_none(kwargs.get('products')) or ())
    transition_state = _literal_or_none(kwargs.get('transitionState'))

    kinetics_type = None
    kinetics_comment = ''
    kinetics_node = kwargs.get('kinetics')
    if isinstance(kinetics_node, ast.Call):
        kinetics_type = _get_call_name(kinetics_node)
        kinetics_kwargs = _call_keywords(kinetics_node)
        comment = _literal_or_none(kinetics_kwargs.get('comment'))
        if isinstance(comment, str):
            kinetics_comment = comment
        else:
            # Composite kinetics (e.g. MultiArrhenius(arrhenius=[Arrhenius(..., comment='...'), ...]))
            # carry their comments on the nested Arrhenius calls rather than on the top-level call.
            # Walk the whole call tree and join every nested `comment=` string found, so this
            # provenance is not silently dropped for multi-component kinetics.
            kinetics_comment = _nested_kinetics_comment(kinetics_node)

    return PDepPathReaction(
        label=label,
        reactants=reactants,
        products=products,
        transition_state=transition_state,
        kinetics_type=kinetics_type,
        kinetics_comment=kinetics_comment,
    )


def _nested_kinetics_comment(kinetics_node: ast.Call) -> str:
    """
    Collect every ``comment=`` string keyword found on any call nested inside a kinetics call.

    Composite kinetics such as ``MultiArrhenius(arrhenius=[Arrhenius(..., comment='...'), ...])``
    carry their per-component comments on the nested calls rather than on the top-level call, so a
    lookup that only checks the top-level ``comment`` keyword silently drops them. This walks the
    entire call subtree (via ``ast.walk``) and joins every nested ``comment=`` string constant it
    finds, in traversal order, with newlines.

    Args:
        kinetics_node (ast.Call): The top-level ``kinetics=`` call node.

    Returns:
        str: The joined nested comments, or ``''`` if none were found.
    """
    comments = []
    for node in ast.walk(kinetics_node):
        if not isinstance(node, ast.Call):
            continue
        nested_comment = _literal_or_none(_call_keywords(node).get('comment'))
        if isinstance(nested_comment, str) and nested_comment:
            comments.append(nested_comment)
    return '\n'.join(comments)


def _get_call_name(call: ast.Call) -> str | None:
    """
    Get the callee name of an ``ast.Call`` node (e.g., ``'Arrhenius'`` for ``Arrhenius(...)``).

    Args:
        call (ast.Call): The call node.

    Returns:
        str: The callee name, or ``None`` if it could not be determined.
    """
    if isinstance(call.func, ast.Name):
        return call.func.id
    if isinstance(call.func, ast.Attribute):
        return call.func.attr
    return None


def _call_keywords(call: ast.Call, path: str = '', call_name: Optional[str] = None) -> dict:
    """
    Map a call's keyword argument names to their (still-AST) value nodes.

    A ``**kwargs`` unpacking (e.g. ``reaction(**payload)``) has ``ast.keyword.arg is None``: this
    module never executes the file, so it has no ``payload`` to expand and cannot know which
    keywords (if any) such an unpacking would actually supply. Previously this function silently
    dropped those keyword entries, which let a top-level DSL call such as ``reaction(**payload)``
    or ``species(**payload)`` parse "successfully" with every field absent instead of being
    refused -- the parser then reported a network that did not describe what the file actually
    declares (fail-open). When ``call_name`` is given (the caller is walking a RECOGNIZED
    top-level DSL call -- see ``RECOGNIZED_TOP_LEVEL_CALLS``), a ``**kwargs`` unpacking now raises
    instead, matching how the rest of this module already refuses a file it cannot resolve
    without executing it (e.g. ``_literal_or_raise``). The check runs unconditionally on every
    keyword, including when other, explicit keywords are also present (e.g.
    ``reaction(reactants=['a'], **payload)``), so a naive "kwargs ended up empty" test would not
    catch it.

    The refusal is NOT scoped to the recognized top-level calls, though only those were reported.
    Every one of this function's six call sites reads keywords and then decides something from
    whether a keyword is present -- the nested ``kinetics=`` value, an Arkane ``pdepreaction(...)``,
    a ``comment=``. Answering "absent" for a keyword that an unpacking may well supply is the same
    fail-open wherever it happens, so gating the refusal on a ``call_name`` argument would leave the
    guard open by default and make its correctness depend on every future call site remembering to
    opt in. Refusing everywhere costs nothing measurable: across the 156 parseable files under
    RMG-Py's ``examples/`` and ``arkane/data/``, ZERO calls use a ``**kwargs`` unpacking at all --
    these files are machine-written with explicit keywords, and a hand-written one that uses ``**``
    is exactly the case this module cannot resolve.

    Args:
        call (ast.Call): The call node.
        path (str): The pdep network file path being parsed (used only in the error message).
        call_name (Optional[str]): The call name being parsed (e.g. ``'reaction'``), used only to make
            the error message name the offending call. Its absence does not weaken the check.

    Raises:
        ValueError: If ``call`` contains a ``**kwargs`` unpacking.

    Returns:
        dict: Keyword name -> AST value node.
    """
    keywords = dict()
    # ``call.args`` is the sibling hole to ``**kwargs``, and the unpacking fix walked past it: this
    # function reads ``call.keywords`` and never looked at ``call.args`` at all. Arkane's DSL functions
    # take positional parameters (``reaction(reactants, products, ...)`` at arkane/input.py), so
    # ``reaction(['A'], ['B'], transitionState='TS1')`` is a legal call that parsed as a reaction with
    # NO reactants and NO products. Same measured justification as the unpacking: zero of the 156
    # parseable files under RMG-Py's examples/ and arkane/data/ pass a positional argument to any of
    # these calls. Mapping positions to names would mean duplicating Arkane's signatures here and
    # keeping them in step with it forever, which is a bigger promise than refusing a spelling nothing
    # writes.
    if call.args:
        raise ValueError(f"The pdep network file at '{path}' calls "
                         f"'{call_name or _get_call_name(call) or '<unnamed>'}(...)' with "
                         f"{len(call.args)} POSITIONAL argument(s); refusing to parse it because this module "
                         f"reads keywords by name and would silently record the positional values as absent.")
    for kw in call.keywords:
        if kw.arg is None:
            raise ValueError(f"The pdep network file at '{path}' calls "
                             f"'{call_name or _get_call_name(call) or '<unnamed>'}(...)' with a "
                             f"'**kwargs' unpacking; refusing to parse it because this module never "
                             f"executes the file and so cannot resolve which keywords the unpacking "
                             f"would actually supply.")
        keywords[kw.arg] = kw.value
    return keywords


def _literal_or_none(node):
    """
    Safely evaluate an AST node to a Python literal, if possible.

    Args:
        node: An AST node (or ``None``).

    Returns:
        The evaluated literal, or ``None`` if ``node`` is ``None`` or is not a literal
        (e.g., it is itself a ``Call`` such as ``Arrhenius(...)``).
    """
    if node is None:
        return None
    try:
        return ast.literal_eval(node)
    except (ValueError, TypeError):
        return None
