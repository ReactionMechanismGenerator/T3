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
import io
import tokenize
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from t3.pdep.hashing import hash_bytes
# network_thermo_t_max is a network-file parsing concern callers reasonably look for in this
# module, but it lives in t3.utils so that t3.utils.writer can use it without a circular import
# (t3.pdep -> t3.utils.writer is the existing dependency direction; several t3.pdep submodules
# import from t3.utils.writer, so t3.utils.writer cannot import from t3.pdep in turn).
from t3.utils.network_thermo import (NetworkTextUnparseable,  # noqa: F401
                                     TGridClampRecord,  # noqa: F401
                                     format_skipped_species,  # noqa: F401
                                     network_thermo_t_max)  # noqa: F401

# Declared because four of the names above are re-exports: used by importers of THIS module rather
# than by anything in it. The per-line suppressions state that to ruff, but only to ruff -- CodeQL's
# py/unused-import asks solely whether a name is used within its own file, follows no cross-module
# re-export, and reads no linter pragma, so it reported all of them (alert 150). ``__all__`` is the
# language's own way to say "this is the public surface", and the query honours it.
#
# The list is not hand-guessed: it is every name any other module imports from ``t3.pdep.parser``,
# taken by parsing the repo's import statements. That sweep also found the alert was not purely a
# false positive as first reported -- ``NetworkThermoCeiling`` was re-exported here and imported
# from here by nobody, so it is gone rather than declared. Importers take it from
# ``t3.utils.network_thermo``, which is where it lives.
__all__ = [
    'NetworkTextUnparseable',
    'PDepArkaneReaction',
    'PDepNetwork',
    'PDepNetworkE0',
    'PDepPathReaction',
    'TGridClampRecord',
    'canonical_channel_pair',
    'format_skipped_species',
    'network_thermo_t_max',
    'parse_arkane_pdep_output_file',
    'parse_arkane_pdep_output_text',
    'parse_pdep_network_e0_file',
    'parse_pdep_network_e0_text',
    'parse_pdep_network_file',
    'parse_pdep_network_text',
    'to_json_safe',
]

RECOGNIZED_TOP_LEVEL_CALLS = {'species', 'transitionState', 'reaction', 'network', 'pressureDependence'}

# kJ/mol per one of each energy unit RMG itself accepts for an E0 quantity: the five named
# ``rmgpy.quantity.Energy`` units plus its two extra-dimensionality conversions ('K' via R and
# 'cm^-1' via h*c*100*Na -- see rmgpy/quantity.py:747-756). Factors are computed from CODATA 2018
# exact constants rather than typed in rounded, so they match rmgpy's own conversions to full
# precision. A unit outside this table is refused, never guessed: a wrong energy scale would
# produce a diagram that looks plausible and is quantitatively wrong everywhere.
_AVOGADRO = 6.02214076e23           # mol^-1
_PLANCK = 6.62607015e-34            # J*s
_SPEED_OF_LIGHT = 299792458.0       # m/s
_ELEMENTARY_CHARGE = 1.602176634e-19  # C
_GAS_CONSTANT = 8.31446261815324    # J/(mol*K)
E0_UNIT_TO_KJ_PER_MOL = {
    'J/mol': 1e-3,
    'kJ/mol': 1.0,
    'cal/mol': 4.184e-3,
    'kcal/mol': 4.184,
    'eV/molecule': _ELEMENTARY_CHARGE * _AVOGADRO * 1e-3,
    'cm^-1': _PLANCK * _SPEED_OF_LIGHT * 100.0 * _AVOGADRO * 1e-3,
    'K': _GAS_CONSTANT * 1e-3,
}

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

# Kinetics constructors that may appear NESTED inside another kinetics call's rate payload, and
# whose own keywords are therefore rate data to be walked rather than an opaque expression. The
# common case is ``PDepArrhenius(pressures=..., arrhenius=[Arrhenius(...), ...])``, which is what
# Arkane emits for ``interpolationModel = ('pdeparrhenius',)``; ``Troe``/``Lindemann``/``ThirdBody``
# nest their high- and low-pressure limits the same way.
#
# Membership is a whitelist on purpose. A call NOT listed here is left to the fail-closed ``else``
# branch, which surfaces it as a ``None`` leaf: recursing into an arbitrary call would let
# ``array([[1.0, 2.0]])`` or ``float('nan')`` donate their arguments to the ME-success gate as if
# those were fitted rate coefficients. Being wrong in the "unrecognized" direction costs a
# needless regeneration; being wrong the other way passes a network whose rates are not real.
NESTED_KINETICS_CALL_NAMES = ('Arrhenius', 'PDepArrhenius', 'MultiArrhenius', 'MultiPDepArrhenius',
                              'Chebyshev', 'Troe', 'Lindemann', 'ThirdBody')


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


def _canonical_channel(labels) -> tuple:
    """
    Canonically order the species labels of a single network channel.

    Arkane sorts the species of every channel it reads or derives, so membership comparisons
    between channels are order-insensitive on its side. T3 holds labels rather than species
    objects, so it sorts the labels instead. Both sides of every comparison must be canonicalized
    the same way or a membership test silently misses and a channel is counted twice.

    Args:
        labels (iterable): The species labels of one channel.

    Returns:
        tuple: The labels, sorted.
    """
    return tuple(sorted(labels))


def canonical_channel_pair(reactant_labels, product_labels) -> tuple:
    """
    Canonically order the two sides of one net reaction into a direction-insensitive pair key.

    This is the identity of a net reaction as far as the network's topology is concerned: the two
    channels it connects, with no direction. Both the prediction
    (``PDepNetwork.expected_net_reaction_channel_pairs``) and the observation (a parsed
    ``pdepreaction(...)`` entry) must be keyed through THIS function -- a check that canonicalizes
    the two sides differently compares two things that were never the same shape, and passes.

    Either side may be passed already canonical: ``_canonical_channel`` sorts, so applying it to an
    already-sorted channel is a no-op. That idempotence is why the enumeration in
    ``PDepNetwork.expected_net_reaction_channel_pairs`` can feed its own canonical channels straight
    back in rather than needing a second, subtly-different helper.

    Args:
        reactant_labels (iterable): The species labels on one side of the net reaction.
        product_labels (iterable): The species labels on the other side.

    Returns:
        tuple: A 2-tuple of canonical channels, itself in a canonical order.
    """
    return tuple(sorted((_canonical_channel(reactant_labels), _canonical_channel(product_labels))))


def _derive_product_channels(path_reactions: tuple,
                             isomers: tuple,
                             reactant_channels: tuple,
                             ) -> tuple:
    """
    Derive a network's product channels from its path reactions.

    RMG-generated pdep network files never write the ``products`` keyword of the ``network(...)``
    call, so Arkane derives the product channels itself: any path reaction side that is neither a
    declared unimolecular isomer nor a declared reactant channel becomes a product channel. This
    mirrors that derivation so T3 sees the same channels Arkane will, without importing or
    executing any RMG code. Reading the keyword literally instead yields no product channels at
    all for every real file, which under-counts the net reactions Arkane writes.

    The two arity branches deliberately test different memberships, exactly as Arkane does: a
    unimolecular side is tested against the isomers, a multi-species side against the declared
    reactant channels. Reactant sides are considered before product sides, per path reaction, in
    file order, so the resulting order matches Arkane's.

    Args:
        path_reactions (tuple): The parsed ``PDepPathReaction`` entries, in file order.
        isomers (tuple): The declared isomer labels.
        reactant_channels (tuple): The declared reactant channels, canonically ordered.

    Returns:
        tuple: The derived product channels, each a canonically ordered tuple of species labels.
    """
    product_channels = list()
    for path_reaction in path_reactions:
        for side in (_canonical_channel(path_reaction.reactants), _canonical_channel(path_reaction.products)):
            if len(side) == 1:
                if side[0] not in isomers and side not in product_channels:
                    product_channels.append(side)
            elif len(side) > 1:
                if side not in reactant_channels and side not in product_channels:
                    product_channels.append(side)
    return tuple(product_channels)


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
        reactant_channels (tuple): The ``network(...)`` reactant channels, each a canonically
                                   ordered tuple of str.
        product_channels (tuple): The RESOLVED product channels, each a canonically ordered tuple
                                  of str. RMG-generated files omit the ``products`` keyword and
                                  let Arkane derive these from the path reactions, so they are
                                  derived here too whenever the keyword is absent.
        product_channels_declared (bool): Whether ``product_channels`` were read from an explicit
                                          ``products`` keyword (``True``) or derived from the path
                                          reactions (``False``).
        source_hash (str, optional): The ``'sha256:<hex>'`` content hash of the exact bytes this
                                     network was parsed from, or ``None`` when there were no file
                                     bytes to hash (i.e. ``parse_pdep_network_text``). ``network_id``
                                     is only a FILE STEM, so it identifies a name, not a content:
                                     two different revisions of ``network4_2.py`` share it. This is
                                     what lets a decision made about one revision be refused as a
                                     gate for a run against another (see
                                     ``t3.pdep.api.explore_pdep_network``).
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
    product_channels_declared: bool = False
    source_hash: str | None = None

    def expected_net_reaction_count(self) -> int:
        """
        The number of live ``pdepreaction(...)`` entries Arkane writes for this network.

        This is deliberately just the length of ``expected_net_reaction_channel_pairs()`` rather
        than a second, independent enumeration. The count and the pair set are two views of one
        prediction, and a consumer that checks both against the same output must not be able to
        get contradictory answers from them.

        Returns:
            int: The expected number of live ``pdepreaction(...)`` entries.
        """
        return len(self.expected_net_reaction_channel_pairs())

    def expected_net_reaction_channel_pairs(self) -> tuple:
        """
        The channel pairs Arkane writes a live ``pdepreaction(...)`` entry for, as a set.

        Arkane enumerates net reactions in ``PressureDependenceJob.save()`` as a double loop whose
        source index runs over the isomers and reactant channels only, and whose destination index
        runs over those plus the product channels, skipping the diagonal. Within the
        isomer/reactant block each unordered pair is therefore visited twice, and the second
        occurrence is written out commented, contributing no parsed entry; every
        isomer-or-reactant-channel to product-channel pair is visited exactly once and is always
        written live.

        No branch of that loop omits an entry -- a rate fit that fails raises, and a rejected fit
        still writes its entry with empty coefficients -- so this count is exact rather than an
        upper bound. That is what makes it usable to detect a truncated ``output.py``.

        The duplicate suppression is deliberately simulated rather than reduced to a closed form.
        Arkane compares each candidate against EVERY previously printed reaction, by label and in
        either direction, not merely against the ones inside the isomer/reactant block. A single
        configuration may legitimately appear in both the source list and the product list -- a
        declared unimolecular reactant channel is also derived as a product channel, since the
        derivation's unimolecular branch tests only against the isomers -- and the entries reaching
        it from the other sources are then commented out as duplicates of the entries leaving it.
        The closed form ``S*(S-1)/2 + S*P`` silently over-counts exactly those, which would reject
        a complete solve. Measured: a network with one such overlapping channel writes 15 live
        entries where the closed form predicts 18.

        Pairs are returned canonically unordered (``_canonical_channel_pair``) because the
        enumeration's own duplicate suppression is direction-insensitive -- Arkane compares by label
        in either direction -- so which side of a live entry is the reactant side is a detail of the
        solver's net-reaction ordering, not part of the network's identity. That canonicalization
        cannot collapse two live entries into one: the suppression already guarantees that if
        ``(a, b)`` was printed then ``(b, a)`` never is, so the returned tuple has exactly one entry
        per live ``pdepreaction(...)``, which is what lets ``expected_net_reaction_count()`` be its
        length.

        Returns:
            tuple: The expected live channel pairs, each a canonically ordered 2-tuple of channels.
        """
        sources = self.isomers_as_channels() + self.reactant_channels
        destinations = sources + self.product_channels
        printed = list()
        for destination_index, destination in enumerate(destinations):
            for source_index, source in enumerate(sources):
                if source_index == destination_index:
                    continue
                if any((source, destination) in ((other_source, other_destination),
                                                 (other_destination, other_source))
                       for other_source, other_destination in printed):
                    continue
                printed.append((source, destination))
        return tuple(canonical_channel_pair(source, destination) for source, destination in printed)

    def isomers_as_channels(self) -> tuple:
        """
        The isomer labels expressed as single-species channels.

        Isomers are stored as bare labels while reactant and product channels are tuples of
        labels. Net reaction enumeration compares the three against one another, so they must
        share a shape.

        Returns:
            tuple: One single-element tuple per isomer.
        """
        return tuple((isomer,) for isomer in self.isomers)

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


def to_json_safe(value):
    """
    Recursively convert a value into plain YAML-safe types (and JSON-safe apart from non-finite
    floats -- see below; the shorthand "JSON/YAML-safe" is used throughout this module for that,
    since the one exception is deliberate and documented rather than an oversight).

    Tuples and lists become lists; dicts are copied with their values recursively converted;
    everything else (including ``None``) is returned unchanged. ``None`` is deliberately passed
    through untouched rather than being treated as "nothing to convert", since a ``None`` leaf
    inside ``kinetics_params`` (e.g. a rejected CSE solve's all-``None`` Chebyshev ``coeffs``) is
    meaningful and must survive the conversion exactly.

    Anything that is not a container and not one of the scalar types below is REFUSED rather than
    passed through. Passing it through is what a converter like this normally does, and it is the
    fail-open shape: the object survives the "conversion" unchanged, the returned structure still
    claims to be JSON/YAML-safe, and the failure surfaces much later at ``yaml.safe_dump`` -- or,
    worse, does not surface at all, because ``yaml.dump`` will happily emit a ``!!python/object``
    tag that only explodes on the way back in. The refusal keeps the function's promise checkable
    at the point the wrong type enters, and a caller with a dataclass converts it explicitly (as
    ``PDepExplorationResult.as_dict`` does for its selection and k(T,P) entries).

    Non-finite floats (``nan``, ``inf``) are deliberately NOT refused, though they are representable
    in YAML and not in strict JSON. A ``nan`` here is real evidence rather than a mistake: a
    degenerate or rejected solve genuinely produces one, and it reaches this function only by having
    been parsed out of a solver's own output. Refusing it would discard the record of a run whose
    numbers went wrong -- precisely when the record matters most -- so the caller writing strict JSON
    is the one that must decide what to do about it.

    Args:
        value: The value to convert, of arbitrary nesting depth.

    Raises:
        TypeError: If a leaf is not a str, bool, int, float, or None, or a mapping key is not a str.

    Returns:
        The converted value, containing no tuples and no non-scalar leaves anywhere.
    """
    if isinstance(value, (tuple, list)):
        return [to_json_safe(item) for item in value]
    if isinstance(value, dict):
        # Keys are checked, not just values. A non-str key round-trips WRONG rather than loudly:
        # yaml.safe_dump happily writes an int or tuple key, and json.dump silently stringifies it,
        # so {1: 'a'} reloads as {'1': 'a'} and a later lookup by the original key misses. Refusing
        # here keeps the corruption at the point it enters instead of at some future read.
        for key in value:
            if not isinstance(key, str):
                raise TypeError(f'Cannot render a mapping keyed by {key!r} of type '
                                f'{type(key).__name__} as a JSON/YAML-safe value; keys must be '
                                f'strings, or the key changes type across a save/load round trip.')
        return {key: to_json_safe(val) for key, val in value.items()}
    # bool before int is unnecessary here (both are allowed) but the tuple order is kept explicit
    # as a reminder that bool IS an int subclass, which matters wherever these two are told apart.
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f'Cannot render {value!r} of type {type(value).__name__} as a JSON/YAML-safe '
                    f'value. Convert it explicitly at the call site; passing it through unchanged '
                    f'would produce a structure that only claims to be serializable.')


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
            omitted. NOTE (deliberate, documented lossiness): ``as_dict()`` renders every tuple
            nested inside ``kinetics_params`` -- including a ``(300, 'K')`` bound pair -- as a
            plain list for YAML-safety, so a save/load round trip through ``as_dict()`` returns a
            list where a tuple went in. This is not a bug to fix: ``t3.pdep.yaml_safe`` deliberately
            rejects Python type tags (``!!python/tuple`` and friends) on any file read through a
            public entrypoint, since a caller-supplied path fully controls what such a tag would
            construct (see that module's docstring), so tuples cannot be made to survive the round
            trip without reopening that hole. The identical lossiness is documented for
            ``PDepExplorationResult.manifest`` in ``t3.pdep.explorer.result``.
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

    def as_dict(self) -> dict:
        """
        Render this entry as plain JSON/YAML-safe types.

        ``kinetics_params`` is copied recursively rather than shallow-copied: its values may nest
        tuples/lists/dicts several levels deep (e.g. a ``PDepArrhenius`` payload's ``arrhenius=[...]``
        list of per-component dicts), and every tuple anywhere in that structure -- including a
        ``(300, 'K')`` bound pair -- must become a list so the whole result stays free of dataclass
        instances and tuples, matching the house style set by ``PDepNetworkSelection.as_dict()``.
        ``None`` leaves (e.g. a rejected CSE solve's all-``None`` ``coeffs``) are preserved exactly:
        detecting them is the whole point of carrying them (see the class docstring), so ``as_dict()``
        must not drop or coerce them the way a naive "skip falsy values" copy would.

        This is a deliberate, documented round-trip lossiness (see ``kinetics_params`` in the class
        docstring): a save/load round trip through this method returns a list where a tuple went in.

        Returns:
            dict: A plain dict, containing no dataclass instances or tuples anywhere, including
                nested inside ``kinetics_params``.
        """
        return {
            'reactants': list(self.reactants),
            'products': list(self.products),
            'kinetics_type': self.kinetics_type,
            'kinetics_params': to_json_safe(self.kinetics_params),
            'numeric_values': list(self.numeric_values),
            'rate_payload_numeric_values': list(self.rate_payload_numeric_values),
            'missing_kinetics_keys': list(self.missing_kinetics_keys),
        }


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
        PDepNetwork: The parsed network, with ``source_hash`` set to the content hash of the exact
            bytes that were parsed.
    """
    # Read the bytes ONCE and hash those bytes, rather than parsing the file and then calling
    # t3.pdep.hashing.hash_file(path) to hash it again. Two reads means the recorded hash can
    # describe content that was never parsed: if the file is replaced between the two, the returned
    # network is the old revision while its provenance claims the new one -- and this hash exists
    # precisely to catch a network file changing underneath a decision.
    with open(path, 'rb') as f:
        data = f.read()
    # Decoded the way Python itself decodes source, via the PEP 263 encoding cookie, rather than with
    # a hard-coded 'utf-8' or the locale encoding open(path, 'r') would have used. These files ARE
    # Python source: RMG writes them, a chemist may hand-edit one, and a '# -*- coding: latin-1 -*-'
    # header or a UTF-8 BOM is legal in both cases. Hard-coded utf-8 rejects the first and turns the
    # second into a leading \ufeff that ast.parse then refuses; the locale encoding makes the parse
    # result depend on the machine's LANG. detect_encoding handles both and strips the BOM.
    encoding, _ = tokenize.detect_encoding(io.BytesIO(data).readline)
    text = data.decode(encoding)
    network_id = Path(path).stem
    return parse_pdep_network_text(text=text, network_id=network_id, path=path,
                                   source_hash=hash_bytes(data))


def parse_pdep_network_text(text: str, network_id: str, path: str = '',
                            source_hash: str | None = None) -> PDepNetwork:
    """
    Parse RMG pdep network file text into a ``PDepNetwork``.

    Args:
        text (str): The RMG pdep network file content.
        network_id (str): An identifier for this network (e.g., the file stem).
        path (str, optional): The path the text was read from, if any (stored for reference).
        source_hash (str, optional): The ``'sha256:<hex>'`` hash of the FILE BYTES ``text`` was
            decoded from, when there was a file. Deliberately not computed from ``text`` here: a
            hash of the decoded string would not match ``t3.pdep.hashing.hash_file`` for any file
            whose bytes differ from its decoded form (CRLF line endings, a BOM, a non-UTF-8 encoding
            cookie), and a provenance hash that silently disagrees with the one the SA cache sidecar
            records is worse than no hash at all. Callers with no file leave it ``None``.

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
    product_channels_declared = False

    for node in tree.body:
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
            continue
        call = node.value
        call_name = _get_call_name(call)
        if call_name not in RECOGNIZED_TOP_LEVEL_CALLS:
            continue
        kwargs = _call_keywords(call, path=path, call_name=call_name)

        if call_name == 'species':
            label = _literal_or_raise(kwargs.get('label'), path=path, keyword='label',
                                      action="omit it from the species labels")
            if label is not None:
                species_labels.append(label)

        elif call_name == 'transitionState':
            label = _literal_or_raise(kwargs.get('label'), path=path, keyword='label',
                                      action="omit it from the transition state labels")
            if label is not None:
                transition_state_labels.append(label)

        elif call_name == 'reaction':
            path_reactions.append(_parse_reaction(kwargs, path=path))

        elif call_name == 'network':
            network_label = _literal_or_none(kwargs.get('label'))
            # ``isomers`` and ``reactants`` fail closed for the same reason ``products`` does below,
            # and the cost of degrading is higher: they are the SOURCE side of the net reaction loop
            # and the memberships every path reaction side is tested against when product channels
            # are derived. Silently reading an unevaluable value as "no isomers" would both
            # under-count the sources and reclassify every unimolecular side as a product channel,
            # so the expected net reaction count would describe a topology the file never declared.
            # An explicit ``None`` literal is refused here too, rather than treated as absent: the
            # keyword was written, so its intent is not "omitted" and guessing is not ours to do.
            for keyword in ('isomers', 'reactants'):
                if keyword in kwargs and _literal_or_none(kwargs[keyword]) is None:
                    raise ValueError(f"The pdep network file at '{path}' declares a '{keyword}' keyword "
                                     f"that could not be evaluated literally; refusing to read it as an "
                                     f"empty channel list in its place.")
            isomers = tuple(_literal_or_none(kwargs.get('isomers')) or ())
            reactant_channels = tuple(_canonical_channel(channel)
                                      for channel in (_literal_or_none(kwargs.get('reactants')) or ()))
            # A ``products`` keyword that is present but not literal-evaluable (a name, a call, a
            # comprehension) must not be treated as absent: deriving the channels instead would
            # silently answer a different question than the file asks, and the resulting net
            # reaction count would reject complete solves. Refuse the file rather than guess.
            declared_product_channels = _literal_or_none(kwargs.get('products'))
            if 'products' in kwargs and declared_product_channels is None:
                raise ValueError(f"The pdep network file at '{path}' declares a 'products' keyword that "
                                 f"could not be evaluated literally; refusing to derive product channels "
                                 f"in its place.")
            product_channels_declared = declared_product_channels is not None
            product_channels = tuple(_canonical_channel(channel) for channel in (declared_product_channels or ()))

        # ``pressureDependence(...)`` data is not part of the public API here; it is
        # recognized (so it is not silently ignored as "unrecognized") but otherwise unused.

    if not path_reactions:
        raise ValueError(f"The pdep network file at '{path}' contains no reaction(...) entries; nothing to parse.")

    if not product_channels_declared:
        product_channels = _derive_product_channels(path_reactions=tuple(path_reactions),
                                                    isomers=isomers,
                                                    reactant_channels=reactant_channels)
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
        product_channels_declared=product_channels_declared,
        source_hash=source_hash,
    )


@dataclass(frozen=True)
class PDepNetworkE0:
    """
    The E0 values a pdep network file declares, converted to kJ/mol.

    A mapping deliberately SEPARATE from ``PDepNetwork``: the network dataclass is a frozen,
    hashable-by-field topology record that several persistence paths serialize, while these are
    plain dicts (unhashable) that only the PES-diagram path needs. Species without an ``E0``
    keyword are ABSENT from the mapping, never recorded as 0.0 -- a consumer that needs a missing
    species' energy must fail closed itself rather than silently place the species at the origin.

    Args:
        species (dict): Species label -> E0 (kJ/mol).
        transition_states (dict): Transition state label -> E0 (kJ/mol).
    """
    species: dict
    transition_states: dict


class _NotAConstantNumber(Exception):
    """
    Internal signal that an AST node is not a constant arithmetic numeric expression.

    Carries a human-readable description of the offending construct (e.g. ``"a name ('x')"``) so
    the catching call site can build an error message that names what it refused, without every
    call site re-deriving that description itself.
    """

    def __init__(self, description: str):
        self.description = description
        super().__init__(description)


def _describe_non_numeric_node(node) -> str:
    """
    Describe an AST node that ``_fold_constant_number`` refused, for use in an error message.

    Args:
        node: The offending AST node.

    Returns:
        str: A short, human-readable description of the node (e.g. ``"a name ('x')"``).
    """
    if isinstance(node, ast.Name):
        return f"a name ({node.id!r})"
    if isinstance(node, ast.Call):
        return f"a call ({_get_call_name(node) or '<unnamed>'}(...))"
    if isinstance(node, ast.Attribute):
        return f"an attribute access ('.{node.attr}')"
    if isinstance(node, ast.Subscript):
        return "a subscript"
    if isinstance(node, ast.Compare):
        return "a comparison"
    if isinstance(node, ast.BinOp):
        return f"a {type(node.op).__name__} operator"
    return f"a {type(node).__name__} node"


def _fold_constant_number(node) -> float:
    """
    Fold an AST expression node to a float, IF AND ONLY IF it is built entirely from numeric
    literals combined with unary +/- and binary +, -, *, / (parenthesised nesting included).

    This exists because RMG writes a TS's ``E0`` as an arithmetic expression rather than a plain
    literal whenever an energy correction was applied (e.g. ``E0 = (264.497 - 411.735,
    "kJ/mol")``), and ``ast.literal_eval`` -- correctly -- refuses any ``BinOp``. The expression is
    still fully constant, so it can be evaluated safely by hand: this function walks the tree
    itself and combines only numbers with +, -, *, /, so nothing beyond that arithmetic can ever
    run. It NEVER calls ``eval``/``exec``, and a ``Name`` is refused outright rather than looked up
    in any namespace -- there is no namespace here, so there is nothing a name could safely mean.

    Args:
        node: An AST expression node.

    Raises:
        _NotAConstantNumber: If ``node`` (or any subnode) is not built solely from numeric
            literals and +, -, *, / arithmetic -- e.g. a name, a call, an attribute access, a
            subscript, a comparison, or a non-arithmetic operator such as ``**``.
        ZeroDivisionError: If the expression divides by a folded zero.

    Returns:
        float: The folded numeric value.
    """
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) \
            and not isinstance(node.value, bool):
        return float(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        operand = _fold_constant_number(node.operand)
        return -operand if isinstance(node.op, ast.USub) else operand
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
        left = _fold_constant_number(node.left)
        right = _fold_constant_number(node.right)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if right == 0:
            raise ZeroDivisionError('division by zero in a constant numeric expression')
        return left / right
    raise _NotAConstantNumber(_describe_non_numeric_node(node))


def _e0_kj_per_mol(node, *, path: str, call_name: str, label: str) -> float | None:
    """
    Evaluate a ``species(...)``/``transitionState(...)`` ``E0=`` keyword node to kJ/mol.

    The numeric half of the pair is folded via ``_fold_constant_number`` rather than
    ``ast.literal_eval``, so a constant arithmetic expression such as
    ``E0 = (264.497 - 411.735, "kJ/mol")`` -- what RMG writes whenever an energy correction was
    applied to a transition state -- is read exactly rather than refused. The unit half keeps the
    plain-literal reading it always had: it is a bare string, never an expression.

    Args:
        node: The keyword's (still-AST) value node, or ``None`` if the keyword is absent.
        path (str): The file path being parsed (for error messages only).
        call_name (str): ``'species'`` or ``'transitionState'`` (for error messages only).
        label (str): The declaring call's label (for error messages only).

    Raises:
        ValueError: If ``E0`` is present but is not a 2-element tuple, its numeric half is not a
            constant numeric literal/expression (e.g. a name, a call, an attribute access), its
            numeric half divides by zero, or its unit is not one RMG's own Energy quantity accepts
            (see ``E0_UNIT_TO_KJ_PER_MOL``).

    Returns:
        Optional[float]: The E0 in kJ/mol, or ``None`` if the keyword is genuinely absent.
    """
    if node is None:
        return None
    # The unit half must be a bare string constant; there is no arithmetic to fold there. The
    # numeric half is folded separately below so a BinOp there (unreadable by ast.literal_eval
    # alone) does not doom the whole pair.
    if not (isinstance(node, ast.Tuple) and len(node.elts) == 2
            and isinstance(node.elts[1], ast.Constant) and isinstance(node.elts[1].value, str)):
        raise ValueError(f"The pdep network file at '{path}' declares an 'E0' keyword on "
                         f"{call_name}(label={label!r}) that is not a literal (number, "
                         f"unit-string) pair; refusing to read it as absent in its place.")
    magnitude_node, unit_node = node.elts
    try:
        magnitude = _fold_constant_number(magnitude_node)
    except _NotAConstantNumber as e:
        # A written-but-unfoldable E0 (a call such as Quantity(...), a bare name) must raise
        # rather than read as absent: a consumer would then report "no E0 declared" about a
        # species whose file plainly gives one -- the same fail-open shape _literal_or_raise
        # exists to refuse.
        raise ValueError(f"The pdep network file at '{path}' declares an 'E0' keyword on "
                         f"{call_name}(label={label!r}) whose numeric half is {e.description}, "
                         f"not a constant number or arithmetic expression; refusing to read it "
                         f"as absent in its place.") from e
    except ZeroDivisionError as e:
        raise ValueError(f"The pdep network file at '{path}' declares an 'E0' keyword on "
                         f"{call_name}(label={label!r}) whose numeric half divides by zero; "
                         f"refusing to read it as absent in its place.") from e
    unit = unit_node.value
    if unit not in E0_UNIT_TO_KJ_PER_MOL:
        raise ValueError(f"The pdep network file at '{path}' declares E0=({magnitude!r}, "
                         f"{unit!r}) on "
                         f"{call_name}(label={label!r}) with the unrecognized energy unit "
                         f"{unit!r}; the recognized units are "
                         f"{sorted(E0_UNIT_TO_KJ_PER_MOL)}. Refusing to guess an energy scale.")
    return magnitude * E0_UNIT_TO_KJ_PER_MOL[unit]


def parse_pdep_network_e0_file(path: str) -> PDepNetworkE0:
    """
    Read the E0 values (kJ/mol) a pdep network file on disk declares.

    Args:
        path (str): The path to the RMG pdep network file.

    Raises:
        ValueError: As ``parse_pdep_network_e0_text``.

    Returns:
        PDepNetworkE0: The declared species/transition-state E0 values, in kJ/mol.
    """
    with open(path, 'rb') as f:
        data = f.read()
    # Decoded via the PEP 263 cookie for the same reason parse_pdep_network_file is: these files
    # are Python source, and a hard-coded 'utf-8' or the locale encoding mis-reads legal ones.
    encoding, _ = tokenize.detect_encoding(io.BytesIO(data).readline)
    return parse_pdep_network_e0_text(text=data.decode(encoding), path=path)


def parse_pdep_network_e0_text(text: str, path: str = '') -> PDepNetworkE0:
    """
    Read the E0 values (kJ/mol) pdep network file text declares.

    This is the energy companion to ``parse_pdep_network_text``: the same never-execute ast walk,
    reading only the ``E0 = (value, unit)`` keyword of every ``species(...)`` and
    ``transitionState(...)`` call. It deliberately does NOT require any ``reaction(...)`` entries:
    it answers "what energies does this file declare", not "is this a well-formed network" -- the
    topology parser owns the latter.

    Args:
        text (str): The RMG pdep network file content.
        path (str, optional): The path the text was read from, if any (error messages only).

    Raises:
        ValueError: If the text cannot be parsed as Python; if any ``E0`` is present but not a
            literal ``(number, unit_string)`` pair, or carries an unrecognized unit; or if the
            same label is declared twice with CONFLICTING E0 values (last-wins would silently pick
            one of two contradictory statements about one stationary point; an identical duplicate
            carries no ambiguity and is tolerated).

    Returns:
        PDepNetworkE0: The declared species/transition-state E0 values, in kJ/mol. Labels whose
            declaration carries no ``E0`` keyword are absent from the mappings.
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', SyntaxWarning)
            tree = ast.parse(text)
    except SyntaxError as e:
        raise ValueError(f"Could not parse the pdep network file at '{path}' as Python: {e}")

    e0_by_call_name = {'species': dict(), 'transitionState': dict()}
    for node in tree.body:
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
            continue
        call = node.value
        call_name = _get_call_name(call)
        if call_name not in e0_by_call_name:
            continue
        kwargs = _call_keywords(call, path=path, call_name=call_name)
        label = _literal_or_raise(kwargs.get('label'), path=path, keyword='label',
                                  action="omit it from the E0 mapping")
        if label is None:
            continue
        e0 = _e0_kj_per_mol(kwargs.get('E0'), path=path, call_name=call_name, label=label)
        if e0 is None:
            continue
        recorded = e0_by_call_name[call_name]
        if label in recorded and recorded[label] != e0:
            raise ValueError(f"The pdep network file at '{path}' declares "
                             f"{call_name}(label={label!r}) twice with conflicting E0 values "
                             f"({recorded[label]} vs. {e0} kJ/mol); refusing to pick one of two "
                             f"contradictory statements about the same stationary point.")
        recorded[label] = e0
    return PDepNetworkE0(species=e0_by_call_name['species'],
                         transition_states=e0_by_call_name['transitionState'])


def _parse_reaction(kwargs: dict, path: str = '') -> PDepPathReaction:
    """
    Build a ``PDepPathReaction`` from a ``reaction(...)`` call's keyword arguments.

    Args:
        kwargs (dict): A mapping of keyword name to its (still-AST) value node.
        path (str, optional): The file path being parsed, for error messages only.

    Raises:
        ValueError: If ``reactants``, ``products``, or ``transitionState`` is present but not
            literal-evaluable.

    Returns:
        PDepPathReaction: The parsed path reaction.
    """
    label = _literal_or_none(kwargs.get('label'))
    # A present-but-unevaluable ``reactants``/``products`` must not collapse to an empty tuple as
    # if the keyword were absent: ``_derive_product_channels`` has no zero-length branch, so an
    # empty side silently contributes NO product channel, under-counting
    # ``expected_net_reaction_count()``; and callers outside this module (e.g. ``t3.main``) build a
    # TS-search reaction from ``reactants + products``, so an empty tuple there queues a bogus,
    # species-less reaction instead of failing loudly.
    reactants = tuple(_literal_or_raise(kwargs.get('reactants'), path=path, keyword='reactants',
                                        action="read it as an empty reactant list") or ())
    products = tuple(_literal_or_raise(kwargs.get('products'), path=path, keyword='products',
                                       action="read it as an empty product list") or ())
    # An ABSENT ``transitionState`` is legitimate (``path_reactions_by_ts`` relies on it to exclude
    # such reactions from its map), but one that is present and unevaluable must not collapse to
    # the same ``None`` -- that would silently and wrongly exclude a path reaction the file DID
    # associate with a transition state.
    transition_state = _literal_or_raise(kwargs.get('transitionState'), path=path, keyword='transitionState',
                                         action="treat it as having no transition state")

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

    Only a bare ``ast.Name`` callee yields a name. An attribute call such as ``foo.reaction(...)``
    used to report ``'reaction'``, which made every caller that dispatches on this name treat it as
    the Arkane DSL directive it merely resembles -- so a file containing ``foo.network(...)`` or
    ``foo.pdepreaction(...)`` parsed as a network or a net reaction, and the resulting
    ``PDepNetwork``/reaction list then served as evidence in the explorer's artifact-identity checks.
    Arkane's own loader binds these names in a namespace with no builtins and no imports, so a
    directive can only ever BE a bare name; anything with a dot in front of it is, by construction,
    not the directive it is named after.

    Args:
        call (ast.Call): The call node.

    Returns:
        str: The callee name, or ``None`` if the callee is not a bare name.
    """
    if isinstance(call.func, ast.Name):
        return call.func.id
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


def _literal_or_raise(node, *, path: str, keyword: str, action: str):
    """
    Evaluate an AST keyword value node as a literal, failing closed if it is present but unusable.

    ``_literal_or_none`` alone cannot distinguish "the keyword was never written" (the legitimate,
    absent case) from "the keyword was written but its value could not be evaluated as a literal"
    (e.g. a bare name or a nested call) -- both collapse to ``None``. Treating the latter as if it
    were the former is the fail-open defect this module has already had to fix multiple times: an
    unevaluable value must RAISE rather than be silently read as absent, since a caller that reads
    it as absent answers a different question than the file asked and can under-count channels or
    silently drop entries.

    An EXPLICITLY WRITTEN ``None`` literal (e.g. ``reactants=None``) is refused for the same
    reason, consistent with the ``network(...)`` ``isomers``/``reactants``/``products`` guard in
    ``parse_pdep_network_text``: the keyword was written, so its intent is not "omitted" and
    guessing is not ours to do. A real RMG-generated ``network.py`` never emits an explicit
    ``None`` for these keywords, so refusing it costs nothing on legitimate input and closes a
    hole on hand-edited/unusual input.

    Args:
        node: The keyword's (still-AST) value node, or ``None`` if the keyword is absent.
        path (str): The file path being parsed (for the error message only).
        keyword (str): The keyword name (for the error message only).
        action (str): What would otherwise silently happen if this were treated as absent, phrased
            to complete "refusing to ... in its place." (e.g. ``"read it as an empty reactant list"``).

    Raises:
        ValueError: If ``node`` is present but not literal-evaluable, or is literally ``None``.

    Returns:
        The literal-evaluated value, or ``None`` if ``node`` is ``None``, i.e. genuinely absent.
    """
    if node is None:
        return None
    value = _literal_or_none(node)
    if value is None:
        raise ValueError(f"The pdep network file at '{path}' declares a '{keyword}' keyword that "
                         f"could not be evaluated literally; refusing to {action} in its place.")
    return value
