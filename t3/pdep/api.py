"""
t3 pdep api module

The public, standalone entrypoint for deciding whether a pressure-dependent (PDep) reaction
network deserves QM refinement, WITHOUT running T3. This is the module an external client (e.g.
Carmel) imports directly.

It never reimplements the selection logic: every decision is produced by calling
``t3.pdep.selector.select_from_sa_dict``, the same pure core function T3's in-run path
(``T3.determine_species_from_pdep_network`` in ``t3/main.py``) calls -- so a standalone caller and
a live T3 run always agree on the same network.

Scope note (deliberate, not an oversight): a prior design note mentioned accepting "a
T3-standardized YAML and/or an RMG-style network.py". There is no T3-standardized network YAML
schema anywhere in this repository -- it was floated once as an aspiration and never built. This
module therefore only supports the RMG ``network.py`` path (parsed via
``t3.pdep.parser.parse_pdep_network_file``), plus accepting already-parsed/loaded Python objects
directly (a ``PDepNetwork`` and/or a loaded sensitivity dict) so callers who already have these in
memory can skip the file-parsing step. No YAML network schema is invented here.

``explore``/``how`` are trigger-channel arguments for a PES explorer that does not exist yet (it
lands in a future commit). The signature accepts them now, on purpose, so it stays stable once the
explorer is implemented; passing ``explore=True`` today raises ``NotImplementedError`` rather than
silently doing nothing.
"""

from pathlib import Path

from arc.common import save_yaml_file

from t3.pdep.cache import validate_sa_cache
from t3.pdep.parser import PDepNetwork, parse_pdep_network_file
from t3.pdep.selector import (CACHE_STATUS_CACHED_REJECTED,
                              CACHE_STATUS_UNVALIDATED,
                              E0_PERTURBATION_J_PER_MOL,
                              EVALUATION_STATUS_NOT_EVALUATED,
                              SELECTOR_VERSION,
                              STRUCTURES_KEY,
                              PDepNetworkSelection,
                              select_from_sa_dict,
                              validate_selection_thresholds,
                              )
from t3.pdep.yaml_safe import read_sa_yaml_file
from t3.schema import T3Sensitivity

# Derived from the schema's own field defaults (rather than copied literals) so this module's
# defaults cannot silently diverge from the defaults T3 itself uses for a live run.
DEFAULT_RELATIVE_THRESHOLD = T3Sensitivity.model_fields['pdep_SA_threshold'].default
DEFAULT_MIN_DELTA_LN_K = T3Sensitivity.model_fields['pdep_min_delta_ln_k'].default


def select_pdep_network(network: str | PDepNetwork,
                        sa_path: str | None = None,
                        sa_dict: dict | None = None,
                        network_reaction: str | None = None,
                        relative_threshold: float = DEFAULT_RELATIVE_THRESHOLD,
                        min_delta_ln_k: float = DEFAULT_MIN_DELTA_LN_K,
                        perturbation: float = E0_PERTURBATION_J_PER_MOL,
                        method: str | None = None,
                        validate_cache: bool = True,
                        explore: bool = False,
                        how: str | None = None,
                        ) -> PDepNetworkSelection:
    """
    Decide whether a PDep network deserves QM refinement, as a standalone call (no T3 run needed).

    Exactly one of ``sa_path``/``sa_dict`` must be given. When ``network_reaction`` is ``None``
    (the default), every network reaction key present in the SA data (except the non-reaction
    ``'structures'`` metadata entry) is evaluated and combined into one aggregate decision
    answering "does this network deserve QM refinement AT ALL" (see
    ``PDepNetworkSelection.combine``). When a specific ``network_reaction`` is given, only that one
    reaction is evaluated -- mirroring the narrower question T3's in-run path
    (``T3.determine_species_from_pdep_network``) asks per call.

    Cache handling (deliberate safety choice, not a bug): when reading from ``sa_path`` and
    ``validate_cache`` is ``True`` (the default), ``validate_sa_cache`` is run first. A cache that
    is rejected must NOT be silently used to produce a qualifying decision -- unlike T3's in-run
    path, a standalone caller (e.g. Carmel) has no way to regenerate the SA data itself. So a
    rejected cache short-circuits to an unqualified decision that carries the rejection warning(s)
    as evidence, rather than reading and evaluating the (untrusted) SA data anyway. This check is
    skipped when ``sa_dict`` is passed directly, since there is no cache to validate against in
    that case.

    Args:
        network (str | PDepNetwork): Either the path to an RMG ``network*.py`` file (parsed via
            ``parse_pdep_network_file``), or an already-parsed ``PDepNetwork`` (e.g. if the caller
            already has one in memory and wants to skip re-parsing). A ``str`` is always treated as
            a path.
        sa_path (str, optional): The path to an Arkane sensitivity YAML, read with
            ``t3.pdep.yaml_safe.read_sa_yaml_file`` -- a restricted loader that understands the
            ``!!python/tuple``-tagged condition keys these files use but refuses every other
            non-safe YAML tag, since this is a public entrypoint reading a caller-supplied path.
            Exactly one of ``sa_path``/``sa_dict`` must be given.
        sa_dict (dict, optional): An already-loaded Arkane sensitivity dictionary. Exactly one of
            ``sa_path``/``sa_dict`` must be given.
        network_reaction (str, optional): A specific network reaction to evaluate, as
            ``'A + B <=> C'``. ``None`` (the default) evaluates every reaction key in the SA data
            and combines the decisions.
        relative_threshold (float, optional): The relative sensitivity gate; see
            ``select_from_sa_dict``.
        min_delta_ln_k (float, optional): The absolute sensitivity gate; see ``select_from_sa_dict``.
        perturbation (float, optional): The E0 perturbation Arkane applied, in J/mol.
        method (str, optional): The master-equation method the SA was (or should have been)
            generated with, recorded on the decision and checked against the cache when
            ``validate_cache`` is ``True``.
        validate_cache (bool, optional): Whether to validate ``sa_path`` against its T3 sidecar
            before trusting it. Ignored when ``sa_dict`` is given directly. Defaults to ``True``.
        explore (bool, optional): Reserved for the upcoming PES explorer adapter (Commit 4, not yet
            implemented). Must be left ``False``.
        how (str, optional): Reserved alongside ``explore`` for the upcoming PES explorer adapter.
            Only meaningful when ``explore`` is truthy.

    Raises:
        ValueError: If neither or both of ``sa_path``/``sa_dict`` are given, or if ``how`` is given
            without ``explore`` being truthy.
        NotImplementedError: If ``explore`` is truthy -- the PES explorer adapter lands in a future
            commit (Commit 4) and does not exist yet.

    Returns:
        PDepNetworkSelection: The decision (a combined one, if ``network_reaction`` is ``None``).
    """
    if explore:
        raise NotImplementedError(
            'The PES explorer adapter (explore/how) is not implemented yet; it lands in a future commit '
            '(Commit 4, PES explorer). Call select_pdep_network() without explore=True for now.')
    if how is not None:
        raise ValueError(f"'how'={how!r} was given without explore=True; 'how' only has meaning for "
                         f"the PES explorer adapter, which is only engaged when explore=True.")
    if (sa_path is None) == (sa_dict is None):
        sa_dict_repr = '<dict>' if sa_dict is not None else None
        raise ValueError(f'Exactly one of sa_path or sa_dict must be given, '
                         f'got sa_path={sa_path!r} and sa_dict={sa_dict_repr!r}.')
    # Validate the thresholds before any cache lookup or parsing: the CACHE_STATUS_CACHED_REJECTED
    # early-return below never reaches select_from_sa_dict()/coefficient_floor(), so a bad
    # threshold would otherwise silently produce a NOT_EVALUATED placeholder instead of raising.
    validate_selection_thresholds(relative_threshold=relative_threshold,
                                  min_delta_ln_k=min_delta_ln_k, perturbation=perturbation)

    network_path = network if isinstance(network, str) else network.path
    parsed_network = parse_pdep_network_file(path=network) if isinstance(network, str) else network

    cache_status = None
    cache_warnings: list = list()
    if sa_dict is None:
        if validate_cache:
            cache_status, cache_warnings = validate_sa_cache(
                sa_path=sa_path, network_path=network_path, perturbation=perturbation,
                min_delta_ln_k=min_delta_ln_k, method=method)
            if cache_status == CACHE_STATUS_CACHED_REJECTED:
                selection = PDepNetworkSelection(
                    network_id=parsed_network.network_id,
                    network_reaction=network_reaction,
                    method=method,
                    sa_path=sa_path,
                    cache_status=cache_status,
                    thresholds={'relative_threshold': relative_threshold,
                                'min_delta_ln_k': min_delta_ln_k,
                                'perturbation': perturbation,
                                },
                )
                # This decision was never actually computed against the SA data -- the cache was
                # rejected outright -- so it must not be reported as an evaluated (un)qualified
                # decision; see FIX C / evaluation_status.
                selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
                selection.warnings.extend(cache_warnings)
                return selection
        else:
            # validate_cache=False means the caller asked us to trust sa_path without checking its
            # provenance -- that must be recorded as 'unvalidated', not reported as though a real
            # cache-validity check ('generated'/'cached_valid') had actually taken place.
            cache_status = CACHE_STATUS_UNVALIDATED
        sa_dict = read_sa_yaml_file(path=sa_path)

    if network_reaction is not None:
        return select_from_sa_dict(sa_dict=sa_dict, network=parsed_network, network_reaction=network_reaction,
                                   relative_threshold=relative_threshold, min_delta_ln_k=min_delta_ln_k,
                                   perturbation=perturbation, method=method, sa_path=sa_path,
                                   cache_status=cache_status)

    reaction_keys = [key for key in sa_dict.keys() if key != STRUCTURES_KEY and isinstance(key, str)] \
        if isinstance(sa_dict, dict) else list()
    decisions = [select_from_sa_dict(sa_dict=sa_dict, network=parsed_network, network_reaction=key,
                                     relative_threshold=relative_threshold, min_delta_ln_k=min_delta_ln_k,
                                     perturbation=perturbation, method=method, sa_path=sa_path,
                                     cache_status=cache_status)
                 for key in reaction_keys]
    if not decisions:
        selection = PDepNetworkSelection(
            network_id=parsed_network.network_id,
            method=method,
            sa_path=sa_path,
            cache_status=cache_status,
            thresholds={'relative_threshold': relative_threshold,
                        'min_delta_ln_k': min_delta_ln_k,
                        'perturbation': perturbation,
                        },
        )
        selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
        selection.warnings.append(
            f'The sensitivity data for network {parsed_network.network_id} has no reaction keys to evaluate.')
        return selection
    return PDepNetworkSelection.combine(decisions)


def rank_pdep_networks(networks,
                       relative_threshold: float = DEFAULT_RELATIVE_THRESHOLD,
                       min_delta_ln_k: float = DEFAULT_MIN_DELTA_LN_K,
                       perturbation: float = E0_PERTURBATION_J_PER_MOL,
                       validate_cache: bool = True,
                       ) -> list:
    """
    Evaluate several PDep networks and rank them by how strongly they deserve QM refinement.

    Each entry of ``networks`` is evaluated with ``select_pdep_network(network_reaction=None, ...)``
    (the "does this network deserve QM at all" question). Every decision is returned -- qualifying,
    not-evaluated, and unqualified alike -- so a caller can see and audit what happened to every
    network it asked about, rather than having failures vanish silently. The returned list is
    ordered in three tiers:

        1. Qualified decisions, ranked most-deserving-first by the largest ``delta_ln_k`` among
           that decision's ``uncertain_path_reactions`` (descending).
        2. Decisions with ``evaluation_status == 'not_evaluated'`` (the network's files were
           missing/unparseable, its SA data could not be obtained, or its cache was rejected -- see
           ``PDepNetworkSelection.evaluation_status``).
        3. Unqualified (but actually evaluated) decisions.

    Within each tier, ties are broken by ``network_id`` (ascending) for determinism.

    A network whose files are missing or unparseable does not abort the whole call: the error is
    caught, recorded as a warning on a ``not_evaluated`` placeholder decision for that network, and
    evaluation continues with the rest.

    Args:
        networks: An iterable of dicts or tuples/lists, each identifying one network to evaluate.
            Mapping form: a dict with keys ``network_path``, ``sa_path``, and optionally
            ``sa_dict``/``method``. Positional form: ``(network_path, sa_path[, method])``.
        relative_threshold (float, optional): Passed through to every ``select_pdep_network`` call.
        min_delta_ln_k (float, optional): Passed through to every ``select_pdep_network`` call.
        perturbation (float, optional): Passed through to every ``select_pdep_network`` call.
        validate_cache (bool, optional): Passed through to every ``select_pdep_network`` call.

    Raises:
        ValueError: If ``relative_threshold``, ``min_delta_ln_k``, or ``perturbation`` is
            non-finite or out of its allowed range. Raised once, upfront, rather than per-network:
            a bad threshold is a config error affecting the whole call, not a per-network failure,
            so it must not be caught by the per-network ``except Exception`` below and turned into
            N misleading ``not_evaluated`` placeholders.

    Returns:
        list: Every ``PDepNetworkSelection`` decision, ordered qualified, then not-evaluated, then
            unqualified (see above).
    """
    validate_selection_thresholds(relative_threshold=relative_threshold,
                                  min_delta_ln_k=min_delta_ln_k, perturbation=perturbation)
    all_selections = list()
    for entry in networks:
        network_path = sa_path = network_id_hint = None
        try:
            network_path, sa_path, sa_dict, method, network_id_hint = _unpack_network_entry(entry)
            selection = select_pdep_network(network=network_path, sa_path=sa_path, sa_dict=sa_dict,
                                            network_reaction=None, relative_threshold=relative_threshold,
                                            min_delta_ln_k=min_delta_ln_k, perturbation=perturbation,
                                            method=method, validate_cache=validate_cache)
        except Exception as e:
            selection = PDepNetworkSelection(network_id=network_id_hint)
            selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
            selection.warnings.append(f'Could not evaluate network entry {entry!r} '
                                      f'(network_path={network_path!r}, sa_path={sa_path!r}): {e}')
        all_selections.append(selection)

    def _rank_key(selection: PDepNetworkSelection):
        if selection.qualified:
            tier = 0
        elif selection.evaluation_status == EVALUATION_STATUS_NOT_EVALUATED:
            tier = 1
        else:
            tier = 2
        max_delta_ln_k = max((entry.delta_ln_k for entry in selection.uncertain_path_reactions), default=0.0)
        return tier, -max_delta_ln_k, selection.network_id or ''

    return sorted(all_selections, key=_rank_key)


def _unpack_network_entry(entry) -> tuple:
    """
    Normalize one ``rank_pdep_networks`` network entry into its parts.

    Args:
        entry: A mapping (dict) or a positional tuple/list identifying one network.

    Raises:
        ValueError: If ``entry`` is not a mapping or a tuple/list of at least
            ``(network_path, sa_path)``, or a mapping is missing ``network_path``.

    Returns:
        tuple: (network_path, sa_path, sa_dict, method, network_id_hint).
    """
    if isinstance(entry, dict):
        network_path = entry.get('network_path')
        if network_path is None:
            raise ValueError(f"Network entry {entry!r} is missing a 'network_path' key.")
        sa_path = entry.get('sa_path')
        sa_dict = entry.get('sa_dict')
        method = entry.get('method')
    elif isinstance(entry, (tuple, list)):
        if len(entry) < 2:
            raise ValueError(f'Network entry {entry!r} must have at least (network_path, sa_path); '
                             f'got {len(entry)} element(s).')
        network_path = entry[0]
        sa_path = entry[1]
        method = entry[2] if len(entry) > 2 else None
        sa_dict = None
    else:
        raise ValueError(f"Network entry {entry!r} is neither a mapping nor a tuple/list of "
                         f"(network_path, sa_path[, method]); cannot unpack it.")
    network_id_hint = Path(network_path).stem if network_path else None
    return network_path, sa_path, sa_dict, method, network_id_hint


def save_pdep_network_selections(path: str, selections: list) -> str:
    """
    Save a list of PDep network selection decisions to a YAML file.

    The file is a mapping with a ``selector_version`` marker (``t3.pdep.selector.SELECTOR_VERSION``)
    and a ``selections`` list, rather than a bare list, so the on-disk format can evolve (e.g. gain
    new top-level keys) without becoming ambiguous with an old-format file.

    Args:
        path (str): The path to write the YAML file to.
        selections (list): The ``PDepNetworkSelection`` decisions to save.

    Returns:
        str: ``path``, as a string, so callers can chain it.
    """
    content = {'selector_version': SELECTOR_VERSION,
              'selections': [selection.as_dict() for selection in selections],
              }
    save_yaml_file(path=path, content=content)
    return str(path)
