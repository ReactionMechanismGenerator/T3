"""
t3 pdep explorer config module.

Defines ``PDepExplorerConfig``, a frozen, validated description of one PES-exploration request.
It is a plain data holder: it does not run an exploration itself (that is
``t3.pdep.explorer.factory.explorer_factory`` and the concrete ``PESExplorerAdapter``), and it is
deliberately validated at a NARROWER scope than "everything that could ever be wrong with this
request" -- see the class docstring's "what is validated WHERE" section for the full split and why
it is split that way rather than centralized here.
"""

import copy
import os
import types
from collections.abc import Mapping
from dataclasses import dataclass, field

from t3.pdep.explorer.adapter import refuse_bare_string_seed
from t3.pdep.explorer.input_file import validate_explorer_field_values

# There is no dedicated method-validator function or valid-method constant set anywhere in
# t3/pdep/mesolver/ or t3/pdep/explorer/input_file.py to delegate to (checked by grep before writing
# this module): every existing call site either compares directly against
# ``t3.utils.writer.METHOD_MAP`` (``input_file.py:311``, ``t3/pdep/mesolver/arkane.py:104``) or
# nothing at all. ``METHOD_MAP`` IS exactly the closed CSE/MSC/RS-style set the brief for this module
# asked to delegate to if one existed, so 'method' is validated against it here rather than merely
# checked for non-emptiness -- delegating to the same source of truth ``write_arkane_explorer_input_file``
# itself uses, so the two entry points cannot drift apart on what a valid method is.
from t3.utils.writer import METHOD_MAP


def _deep_freeze(value):
    """
    Recursively convert a value into a deeply read-only equivalent.

    ``types.MappingProxyType`` alone freezes only the OUTER mapping: a nested list value fished out
    of the proxy is still the same mutable list, so ``frozen.database_kwargs['thermoLibraries']
    .append(...)`` would mutate a supposedly frozen config after validation. This helper closes
    that hole by rebuilding every container as its immutable counterpart: mappings become read-only
    views over rebuilt dicts, sequences become tuples, sets become frozensets. A non-container leaf
    is deep-copied, so an exotic mutable leaf object at least never aliases the caller's original.

    Args:
        value: The value to freeze.

    Returns:
        The deeply frozen equivalent of ``value``.
    """
    if isinstance(value, Mapping):
        return types.MappingProxyType({key: _deep_freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_deep_freeze(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_deep_freeze(item) for item in value)
    return copy.deepcopy(value)


def deep_thaw(value):
    """
    Recursively convert a ``_deep_freeze``d value back into plain mutable Python containers.

    The consumption boundary needs this: ``t3.pdep.explorer.input_file._validate_database_kwarg``
    deliberately requires a real ``list`` (not any Sequence) for the ``database(...)`` library
    keywords, so a frozen config's tuple values must be thawed back to lists before they are handed
    to the writer. Thawing produces a fresh, independent copy every time, so the frozen config a
    result was made from cannot be edited through anything this function returned.

    Args:
        value: The value to thaw (``None`` passes through as ``None``).

    Returns:
        A plain mutable equivalent of ``value``: dicts, lists, and sets, with deep-copied leaves.
    """
    if isinstance(value, Mapping):
        return {key: deep_thaw(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [deep_thaw(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return {deep_thaw(item) for item in value}
    return copy.deepcopy(value)


@dataclass(frozen=True)
class PDepExplorerConfig:
    """
    A frozen, validated description of one PES-exploration request.

    This is a data holder, not an executor: constructing one does not read or write the filesystem,
    look up the explorer registry, or apply the seed/capability rules
    (``t3.pdep.explorer.adapter.validate_explorer_seed``) that need a concrete adapter CLASS to mean
    anything. What IS validated here is only what can be known from the arguments themselves, at
    construction time, with no other information -- everything else is validated later, by the code
    that actually has the information needed to validate it correctly. Splitting it this way is
    deliberate, not an oversight:

    - Construction time (HERE, ``__post_init__``): types; that ``trusted_output_root`` and
      ``output_directory`` are absolute and that ``output_directory`` is strictly contained in
      ``trusted_output_root``; the numeric ``explorer(...)`` field contracts (via
      ``t3.pdep.explorer.input_file.validate_explorer_field_values``). None of this can go stale
      between construction and use, and none of it needs anything this object does not already carry.
    - Factory/adapter time (``t3.pdep.explorer.factory.explorer_factory`` /
      ``t3.pdep.explorer.adapter.validate_explorer_seed``): the seed rules
      (``max_source_species``, ``supports_transition_state_seeds``). These need the concrete
      explorer CLASS's capabilities, which this config does not know -- ``explorer`` here is just a
      registry-lookup key, not the class itself. Re-implementing the seed rules here would mean
      either duplicating the registry lookup (drifting from the one in ``factory.py``, exactly the
      antipattern ``validate_explorer_seed``'s own docstring warns about) or guessing at capability
      values that might be wrong for whichever class ``explorer`` actually names.
    - API-call time (the function that consumes this config to actually run an exploration, landing
      in a later piece of this feature): filesystem state -- whether ``trusted_output_root`` and
      ``output_directory`` exist, are directories, are writable, etc. Checking this HERE would be
      validating a fact that can already be stale by the time it matters (TOCTOU): nothing stops the
      filesystem changing between this object's construction and its actual use.
    - Write time (``t3.pdep.explorer.input_file.write_arkane_explorer_input_file``): bath-gas species
      existence against the network file, reactive=False marking, renderability of every value as an
      Arkane-loadable literal. None of that is knowable from a config in isolation -- it depends on
      the actual network source file being explored, which this config also does not carry.

    'explorer' is checked only for being a non-empty str, NOT against the explorer registry
    (``t3.pdep.explorer.factory._registered_explorer_adapters``): the registry is the single source
    of truth for which names are actually registered, and a duplicate "is it registered" check here
    would inevitably drift from it (e.g. an adapter registered after this config type was written, or
    only within a particular process). The factory already refuses an unregistered name with a clear
    error; this config only guarantees it received A name.

    'method' is checked against ``t3.utils.writer.METHOD_MAP`` (see the module-level comment above
    this class for why that, and not a bespoke check, is the delegate).

    'bath_gas' and 'database_kwargs', when given, are stored DEEPLY frozen (``_deep_freeze``):
    every nested mapping becomes a read-only ``types.MappingProxyType`` view over a rebuilt dict,
    every nested sequence a tuple, every set a frozenset. A top-level-only
    ``MappingProxyType(deepcopy(...))`` would decouple the stored copy from the caller's original,
    but a nested list fished out of the proxy would still be the config's own mutable list --
    ``config.database_kwargs['thermoLibraries'].append(...)`` would then mutate a supposedly frozen
    config after validation, silently changing the Arkane input eventually generated from it.
    The flip side is that the writer must not be fed these frozen values directly:
    ``t3.pdep.explorer.input_file._validate_database_kwarg`` requires a real ``list`` for the
    ``database(...)`` library keywords, so the consumption boundary
    (``t3.pdep.api.explore_pdep_network``) passes them through ``deep_thaw`` first.

    'seed_species' and 'transition_state_seeds' are only coerced to ``tuple`` here, NOT run through
    ``validate_explorer_seed``. That rule needs ``max_source_species`` and
    ``supports_transition_state_seeds``, which are attributes of the concrete explorer CLASS, not of
    this config -- the same reasoning as 'explorer' above. The factory/adapter apply that rule with
    the correct capability values at the point they actually have them.

    No filesystem path (``trusted_output_root``, ``output_directory``) is checked for existence,
    directory-ness, or writability here -- see the "API-call time" bullet above.

    Args:
        explorer (str): The explorer adapter registry name (e.g. ``'Arkane'``). Not checked against
                        the registry itself; see the class docstring.
        trusted_output_root (str): The absolute path to the root directory explorer output is
                                   confined to.
        output_directory (str): The absolute path to the directory this exploration's files will be
                                written to. Must resolve to a location strictly inside
                                ``trusted_output_root``.
        seed_species (tuple): The source (seed) species labels to explore from, coerced to a tuple.
                              Not otherwise validated here; see the class docstring.
        method (str): The master-equation method, one of ``t3.utils.writer.METHOD_MAP``
                     (e.g. 'CSE', 'MSC', 'RS').
        bath_gas (Mapping, optional): The bath gas composition, mapping species labels to mole
                                      fractions. Stored deeply frozen; see the class docstring.
        explore_tol (float, optional): The energy tolerance for exploring new isomers/reactions.
        energy_tol (float, optional): The energy tolerance for including a well/transition state in
                                      the output network.
        flux_tol (float, optional): The flux tolerance for including a well/transition state in the
                                    output network.
        maximum_radical_electrons (int, optional): The maximum number of radical electrons allowed
                                                   in an explored species. An integral float (e.g.
                                                   ``2.0``) is coerced to ``int`` (``2``).
        transition_state_seeds (tuple): Transition state label(s) offered as additional seeds,
                                        coerced to a tuple. Not otherwise validated here; see the
                                        class docstring. Defaults to an empty tuple.
        database_kwargs (Mapping, optional): Keyword arguments describing the RMG database settings
                                             to use for the exploration. Stored deeply frozen (list
                                             values become tuples); see the class docstring.

    Raises:
        ValueError: If ``explorer`` or ``method`` is not a non-empty str, or ``method`` is not one of
                   ``t3.utils.writer.METHOD_MAP``; if ``trusted_output_root`` or ``output_directory``
                   is not a non-empty, absolute str path, or ``output_directory`` does not resolve
                   strictly inside ``trusted_output_root``; if ``bath_gas`` or ``database_kwargs`` is
                   given and is not a ``Mapping``; or if any numeric field fails its contract in
                   ``t3.pdep.explorer.input_file._EXPLORER_NUMBER_CONTRACTS``.
    """
    explorer: str
    trusted_output_root: str
    output_directory: str
    seed_species: tuple
    method: str
    bath_gas: Mapping | None = None
    explore_tol: float | None = None
    energy_tol: float | None = None
    flux_tol: float | None = None
    maximum_radical_electrons: int | None = None
    transition_state_seeds: tuple = field(default_factory=tuple)
    database_kwargs: Mapping | None = None

    def __post_init__(self):
        """
        Enforce every construction-time rule documented on the class, and nothing more.

        Raises:
            ValueError: See the class docstring.
        """
        for field_name in ('explorer', 'method'):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value:
                raise ValueError(f"'{field_name}' must be a non-empty str, got {value!r} of type "
                                 f"{type(value).__name__}.")
        if self.method not in METHOD_MAP:
            raise ValueError(f"'method' must be one of {sorted(METHOD_MAP)}, got {self.method!r}. Arkane's own "
                             f"name for it is written into the pressureDependence(...) block by "
                             f"t3.pdep.explorer.input_file.write_arkane_explorer_input_file, so an unrecognized "
                             f"method cannot be rendered.")

        for field_name in ('trusted_output_root', 'output_directory'):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value:
                raise ValueError(f"'{field_name}' must be a non-empty str, got {value!r} of type "
                                 f"{type(value).__name__}.")
            if not os.path.isabs(value):
                raise ValueError(f"'{field_name}' must be an absolute path, got {value!r}.")

        # Confinement, mirroring the realpath+commonpath idiom used across t3/pdep/ (e.g.
        # t3.pdep.discovery._confine_to_project): resolved via os.path.realpath on BOTH sides (so a
        # traversal component or a symlink cannot escape it), never str.startswith, which would
        # wrongly accept a sibling directory -- output_directory='/tmp/rootX' against
        # trusted_output_root='/tmp/root' shares the '/tmp/root' PREFIX as a string but is not inside
        # it as a path, and startswith('/tmp/root') cannot tell the difference.
        resolved_root = os.path.realpath(self.trusted_output_root)
        resolved_output_directory = os.path.realpath(self.output_directory)
        if resolved_output_directory == resolved_root \
                or os.path.commonpath([resolved_root, resolved_output_directory]) != resolved_root:
            raise ValueError(f"'output_directory' ({self.output_directory!r}, resolved to "
                             f"{resolved_output_directory!r}) must resolve to a location strictly inside "
                             f"'trusted_output_root' ({self.trusted_output_root!r}, resolved to "
                             f"{resolved_root!r}), and must not be the root itself.")

        # Delegates to the SAME numeric-contract function write_arkane_explorer_input_file uses, so
        # the two entry points cannot drift apart on what a valid explore_tol/energy_tol/flux_tol/
        # maximum_radical_electrons is (see validate_explorer_field_values's own docstring for why it
        # is public now). The COERCED return values are written back via object.__setattr__ --
        # validating and discarding them would let an integral float like maximum_radical_electrons=2.0
        # pass validation but still be rendered later as a non-count float.
        field_values = validate_explorer_field_values(
            explore_tol=self.explore_tol, energy_tol=self.energy_tol, flux_tol=self.flux_tol,
            maximum_radical_electrons=self.maximum_radical_electrons)
        for field_name, coerced_value in field_values.items():
            object.__setattr__(self, field_name, coerced_value)

        # The seed RULES (how many source species, whether transition states may seed at all) are
        # deliberately NOT applied here -- see the class docstring. But normalizing to a tuple is not
        # the neutral act it looks like: tuple('OH') is ('O', 'H'), so coercing a bare string here
        # would convert a plausible single-label typo into a well-formed bimolecular source channel,
        # and validate_explorer_seed's own refusal of that mistake could never fire afterwards --
        # by then there is nothing left to distinguish it from a caller who genuinely asked for
        # O + H. The refusal is therefore shared, and runs at BOTH entry points, before either
        # tuple() call.
        for field_name in ('seed_species', 'transition_state_seeds'):
            value = getattr(self, field_name)
            refuse_bare_string_seed(field_name=field_name, value=value)
            object.__setattr__(self, field_name, tuple(value or tuple()))

        for field_name in ('bath_gas', 'database_kwargs'):
            value = getattr(self, field_name)
            if value is None:
                continue
            if not isinstance(value, Mapping):
                raise ValueError(f"'{field_name}' must be a Mapping (e.g. a dict), got {value!r} of type "
                                 f"{type(value).__name__}.")
            # Deeply frozen (see _deep_freeze): neither mutating the caller's original dict after
            # construction, nor obtaining this attribute and mutating anything reachable through it,
            # can change what this config holds. MappingProxyType(deepcopy(...)) alone would only
            # break the aliasing -- a nested list value would still be THIS config's own mutable
            # list, appendable through the read-only view after validation.
            object.__setattr__(self, field_name, _deep_freeze(value))
