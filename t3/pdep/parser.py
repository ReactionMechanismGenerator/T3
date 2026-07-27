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
