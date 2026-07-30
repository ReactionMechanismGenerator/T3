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

``explore``/``how`` used to be reserved trigger-channel arguments for a PES explorer that did not
exist yet. That explorer now exists as ``explore_pdep_network()`` (below), which is the real entry
point for running one; ``select_pdep_network()`` itself never explores -- passing ``explore=True``
or ``how=<anything>`` to it raises ``ValueError`` redirecting the caller to
``explore_pdep_network()``. The parameters are kept on the signature (deprecated) purely so a
caller still passing them gets that redirect rather than an opaque ``TypeError``.
"""

import copy
import os
from pathlib import Path

from arc.common import save_yaml_file

from t3.pdep.cache import validate_sa_cache
from t3.pdep.explorer.config import PDepExplorerConfig, deep_thaw
from t3.pdep.explorer.factory import explorer_factory
from t3.pdep.explorer.result import (EXPLORATION_RESULT_SCHEMA_VERSION,
                                     EXPLORATION_STATUS_FAILED,
                                     EXPLORATION_STATUS_SKIPPED,
                                     EXPLORATION_STATUS_SUCCEEDED,
                                     PDepExplorationResult,
                                     )
from t3.pdep.parser import PDepNetwork, parse_pdep_network_file
from t3.pdep.selector import (CACHE_STATUS_CACHED_REJECTED,
                              CACHE_STATUS_UNVALIDATED,
                              E0_PERTURBATION_J_PER_MOL,
                              EVALUATION_STATUS_EVALUATED,
                              EVALUATION_STATUS_NOT_EVALUATED,
                              SELECTION_SCHEMA_VERSION,
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
        explore (bool, optional): Deprecated in favour of ``explore_pdep_network()``. This function
            never explores; passing ``explore=True`` raises ``ValueError`` naming
            ``explore_pdep_network()`` as the real entry point. Kept on the signature only so a
            caller still passing it gets that redirect rather than an opaque ``TypeError``.
        how (str, optional): Deprecated in favour of ``PDepExplorerConfig.explorer``, which replaces
            it. Passing any value raises ``ValueError``. Kept on the signature only so a caller
            still passing it gets that redirect rather than an opaque ``TypeError``.

    Raises:
        ValueError: If neither or both of ``sa_path``/``sa_dict`` are given; if ``explore`` is
            truthy; or if ``how`` is given at all. The latter two redirect the caller to
            ``explore_pdep_network()`` / ``PDepExplorerConfig.explorer``, which replace them.

    Returns:
        PDepNetworkSelection: The decision (a combined one, if ``network_reaction`` is ``None``).
    """
    if explore:
        raise ValueError(
            "select_pdep_network() never explores a PES; 'explore'/'how' are deprecated no-ops here. "
            'Call explore_pdep_network() instead -- it is the real entry point for running a PES '
            "explorer, and 'how' is replaced there by PDepExplorerConfig.explorer.")
    if how is not None:
        raise ValueError(
            f"'how'={how!r} was given, but 'how' is deprecated and has no meaning on "
            f"select_pdep_network() (it never explores). Use explore_pdep_network() with a "
            f"PDepExplorerConfig(explorer=...) instead -- 'explorer' replaces 'how'.")
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
                    network_source_hash=parsed_network.source_hash,
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
            network_source_hash=parsed_network.source_hash,
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


def explore_pdep_network(network_path: str,
                         config: PDepExplorerConfig,
                         selection: PDepNetworkSelection | None = None,
                         logger=None,
                         ) -> PDepExplorationResult:
    """
    Run a PES explorer against a PDep network, as a standalone call, optionally gated by a budget.

    Asymmetric with ``select_pdep_network()`` ON PURPOSE: that function accepts ``str |
    PDepNetwork``, this one accepts ONLY a path string. The explorer input writer
    (``t3.pdep.explorer.input_file.write_arkane_explorer_input_file``, via its ``source_path``
    parameter) needs the network's SOURCE TEXT read from a file on disk -- and an already-parsed
    ``PDepNetwork`` is not that: its ``.path`` may be ``None``, may be stale, or may point at a file
    that has since changed on disk. A parsed object is refused here rather than silently re-read
    from a ``.path`` that might not mean what it used to.

    Budget gate (only when ``selection`` is given -- see the ``selection`` Args entry below):
    select-first is a BUDGET POLICY, never a correctness claim. ``select_pdep_network()`` only sees
    reactions already present in the SA dict and the network file, so a network missing the very
    channel that would dominate its sensitivity can never qualify -- and would then never be
    explored, even though exploring it might reveal exactly that missing channel. Passing
    ``selection=None`` is the deliberate "explore regardless of budget" path for a caller who wants
    that.

    Filesystem state (``config.trusted_output_root``, ``config.output_directory``) is checked HERE,
    not in ``PDepExplorerConfig.__post_init__``: the config deliberately checks none of it (see its
    class docstring's "what is validated WHERE" section) because a filesystem fact checked at
    construction time can already be stale by the time this function runs.

    Args:
        network_path (str): The path to an existing RMG P-dep network file. Unlike
            ``select_pdep_network()``'s ``network`` argument, a parsed ``PDepNetwork`` is NOT
            accepted here; see above.
        config (PDepExplorerConfig): The validated description of the exploration request.
        selection (PDepNetworkSelection, optional): The budget-gate decision for this network, e.g.
            from ``select_pdep_network()``. When given, ``selection.method`` must equal
            ``config.method``, ``selection.network_id`` must equal the parsed network's
            ``network_id``, and ``selection.network_source_hash`` must be present and equal the
            content hash of ``network_path`` (all checked here; see Raises). This check is a
            STALE-SELECTION guard, not a run-input-integrity one: it proves only that
            ``network_path`` held the approved bytes at the moment this function parsed it, one
            read among several between here and Arkane actually consuming the file. What binds the
            bytes the explorer input writer consumes is the SEPARATE ``expected_source_hash`` this
            function forwards to ``explorer_factory()`` below, checked against the writer's own
            single read of the same path -- and even that does not close every window: the input
            file Arkane later reads is the one written into ``config.output_directory``, inside
            ``config.trusted_output_root``, and a change to THAT file between the write and
            Arkane's own read is a different exposure this function does not address. When
            ``selection.qualified`` is
            ``False``, the exploration is skipped -- no adapter is ever constructed -- and a
            ``'skipped'`` result carrying ``selection.reason()`` is returned. When ``None`` (the
            default), the exploration runs unconditionally.
        logger (Logger, optional): The current T3 Logger instance, passed through to the explorer
            adapter.

    Raises:
        ValueError: If ``network_path`` is not a ``str`` (see above); if ``selection`` is given and
            its ``method`` does not equal ``config.method``, or its ``network_id`` does not equal
            the parsed network's ``network_id``; if ``selection.network_source_hash`` is ``None`` or
            differs from the content hash of ``network_path`` (a decision bound to no content, or to
            content that has since changed, is not a gate for this run); if
            ``config.trusted_output_root`` does not exist,
            is not a directory, or is a symlink; or if ``config.output_directory`` does not resolve
            strictly inside the realpath'd ``config.trusted_output_root``.
        RuntimeError: Propagated verbatim from the explorer adapter (e.g. a rule-0 run-directory
            claim collision surfaced as an exception, rather than as a recorded 'failed' result).
            This is deliberate -- see the module-level note in the function body just above the
            ``adapter.explore()`` call.

    Returns:
        PDepExplorationResult: The outcome -- 'skipped' (budget gate declined it), 'failed' (the
            explorer ran and did not succeed), or 'succeeded'.
    """
    if not isinstance(network_path, str):
        raise ValueError(
            f'explore_pdep_network() requires network_path to be a path string to an existing '
            f'network file, not a parsed object (got {type(network_path).__name__}). Unlike '
            f"select_pdep_network(), which accepts str | PDepNetwork, this function does not accept "
            f"an already-parsed PDepNetwork: the explorer input writer needs the network's SOURCE "
            f"TEXT read from a file on disk (see write_arkane_explorer_input_file's source_path "
            f"parameter), and a parsed object is not that -- its .path may be None, stale, or point "
            f"at a file that has since changed. Pass the path string instead.")

    parsed_network = parse_pdep_network_file(path=network_path)

    if selection is not None:
        if selection.method != config.method:
            raise ValueError(
                f"selection.method ({selection.method!r}) does not match config.method "
                f"({config.method!r}). Gating a decision made under one master-equation method and "
                f"then exploring under another is silent provenance corruption: the recorded "
                f"decision would not be about the run that actually happened.")
        if selection.network_id != parsed_network.network_id:
            raise ValueError(
                f"selection.network_id ({selection.network_id!r}) does not match the parsed "
                f"network's network_id ({parsed_network.network_id!r}). A decision about a "
                f"different network used as this one's budget gate fabricates confidence in a "
                f"result nothing actually justifies -- the same hazard "
                f"PDepNetworkSelection.combine() already refuses for the same reason.")
        # network_id is a FILE STEM, so the check above only proves the decision was about a file
        # with this NAME. RMG rewrites pdep/network4_2.py on every iteration that touches the
        # network, and a selection is routinely made in one process and acted on in another, so
        # "same stem" and "same network" are not the same statement. Bind to the content instead.
        if selection.network_source_hash is None:
            raise ValueError(
                f"selection.network_source_hash is None, so this decision carries no binding to the "
                f"content it was made about, and network_id ({selection.network_id!r}) is only a file "
                f"stem -- it matches every revision of that file. Re-run the selection against the "
                f"network file (select_pdep_network records the hash whenever it parses one), or pass "
                f"selection=None to explore without a budget gate.")
        if selection.network_source_hash != parsed_network.source_hash:
            raise ValueError(
                f"selection.network_source_hash ({selection.network_source_hash!r}) does not match the "
                f"content hash of {network_path!r} ({parsed_network.source_hash!r}): the network file "
                f"has changed since the decision was made. The decision describes a network that is no "
                f"longer the one about to be explored -- its sensitivities, its channels, and the "
                f"transition states it named may all have moved. Re-run the selection against the "
                f"current file.")
        # `qualified` carries NO signal unless the decision was actually evaluated. When
        # select_pdep_network() rejects a stale SA cache it returns qualified=False AND
        # evaluation_status='not_evaluated' (api.py:174, :208), and reason() nonetheless renders the
        # full "does not qualify: no transition state the network is sensitive to ..." sentence.
        # Gating on `qualified` alone therefore reads a missing evaluation as a negative verdict:
        # the exploration is silently skipped and the caller is handed a decision that was never
        # made, phrased as though it had been. Raise rather than pick a side -- exploring anyway and
        # skipping are both guesses about what the caller meant, and the caller can express either
        # deliberately (re-run the selection, or pass selection=None).
        #
        # `qualified` is the exception, and it is why this is `and not selection.qualified` rather
        # than a bare status check. PDepNetworkSelection.combine() marks an aggregate 'not_evaluated'
        # whenever ANY component was not evaluated, because that is the truth about coverage -- but a
        # partially evaluated aggregate that DID qualify is backed by whichever evaluated component
        # qualified, and combine() only counts evaluated components' votes. Refusing it would refuse
        # a positive verdict that real evidence supports: the over-refusal failure mode, not the
        # fail-open one. The asymmetry is the point -- a partial 'no' is unsupported (an unevaluated
        # component might have been the one that qualified), a partial 'yes' is not.
        if selection.evaluation_status != EVALUATION_STATUS_EVALUATED and not selection.qualified:
            raise ValueError(
                f"selection.evaluation_status is {selection.evaluation_status!r}, so its "
                f"'qualified' field ({selection.qualified!r}) carries no verdict and cannot be used "
                f"as a budget gate -- a decision that was never evaluated is not a decision to not "
                f"explore. Re-run the selection against usable SA data, or pass selection=None to "
                f"explore without gating.")
        if not selection.qualified:
            # Nothing runs: no adapter is constructed, no filesystem state below is touched.
            return PDepExplorationResult(
                network_id=parsed_network.network_id,
                status=EXPLORATION_STATUS_SKIPPED,
                reasons=(selection.reason(),),
                selection=selection,
            )

    # Filesystem state, checked HERE and not in PDepExplorerConfig -- see the function docstring.
    # T3 must never CREATE the trusted root: creating it would mean claiming a path whose ownership
    # the caller never demonstrated (it only ever demonstrated a path STRING).
    if not os.path.isdir(config.trusted_output_root) or os.path.islink(config.trusted_output_root):
        raise ValueError(
            f"config.trusted_output_root ({config.trusted_output_root!r}) must already exist as a "
            f"real (non-symlink) directory; explore_pdep_network() never creates it.")

    # Re-verified rather than trusted from PDepExplorerConfig.__post_init__: that check compared
    # path STRINGS (via realpath) at construction time, when trusted_output_root/output_directory
    # need not have existed at all. Only now, with the root confirmed to exist, can realpath resolve
    # any symlink that actually exists on disk -- so this re-check can catch something the
    # construction-time one structurally could not.
    resolved_root = os.path.realpath(config.trusted_output_root)
    resolved_output_directory = os.path.realpath(config.output_directory)
    if resolved_output_directory == resolved_root \
            or os.path.commonpath([resolved_root, resolved_output_directory]) != resolved_root:
        raise ValueError(
            f"config.output_directory ({config.output_directory!r}, resolved to "
            f"{resolved_output_directory!r}) must resolve to a location strictly inside "
            f"config.trusted_output_root ({config.trusted_output_root!r}, resolved to "
            f"{resolved_root!r}).")

    # Intermediate directories BETWEEN the root and output_directory are created here -- permitted
    # precisely because they are inside a root the caller has already vouched for (checked above).
    # config.output_directory ITSELF is deliberately NOT created: the adapter's rule-0 os.mkdir
    # (t3.pdep.explorer.arkane.ArkaneExplorerAdapter._claim_run_directory) must stay the SOLE creator
    # of that leaf directory -- that mkdir winning IS the atomic claim of ownership over the run, and
    # pre-creating the directory here would make every run look like it collided with a previous one
    # (rule 0 refuses a pre-existing directory, even an empty one; see its docstring). Do not
    # "helpfully" add an os.makedirs(config.output_directory, exist_ok=True) below -- that is exactly
    # the line that would break rule 0.
    os.makedirs(os.path.dirname(config.output_directory), exist_ok=True)

    adapter = explorer_factory(
        explorer=config.explorer,
        seed_species=config.seed_species,
        output_directory=config.output_directory,
        network_path=network_path,
        method=config.method,
        # Thawed at this boundary, not passed frozen: the config stores bath_gas/database_kwargs
        # DEEPLY frozen (nested lists as tuples, nested mappings as read-only views), and the input
        # writer's _validate_database_kwarg deliberately requires a real list for the database(...)
        # library keywords. deep_thaw hands the adapter a fresh plain-dict/list copy, so the writer's
        # contract is met and nothing downstream can mutate the frozen config through what it got.
        bath_gas=deep_thaw(config.bath_gas),
        explore_tol=config.explore_tol,
        energy_tol=config.energy_tol,
        flux_tol=config.flux_tol,
        maximum_radical_electrons=config.maximum_radical_electrons,
        logger=logger,
        transition_state_seeds=config.transition_state_seeds,
        database_kwargs=deep_thaw(config.database_kwargs),
        # The check above (selection.network_source_hash vs. parsed_network.source_hash) proves only
        # what the bytes at network_path WERE at the moment this function parsed them. It says
        # nothing about what write_arkane_explorer_input_file() will later read, because that is a
        # separate open() of the same path. Forwarding the hash here is what binds the two reads: the
        # writer refuses if its own read does not match, rather than silently exploring whatever
        # bytes happen to be there by the time it gets to them.
        expected_source_hash=parsed_network.source_hash,
    )

    # adapter.explore() calls its own set_up() internally (see ArkaneExplorerAdapter.set_up's
    # docstring on why that ordering is load-bearing); explore_pdep_network() must not call set_up()
    # itself.
    #
    # No catch-all here, deliberately. An ordinary Arkane run failure is a RECORDED result
    # (status='failed', via adapter.reasons) -- that is adapter.explore() returning False, not
    # raising. Everything else -- a bad config, an unreadable network, an unknown explorer name from
    # the factory, a RuntimeError out of the adapter -- must stay loud and propagate. Wrapping this
    # call in `except Exception` would relabel T3 bugs and invalid input as "the explorer returned a
    # negative result", which is a wrong statement of fact a caller could act on.
    succeeded = adapter.explore()

    if succeeded:
        return PDepExplorationResult(
            network_id=parsed_network.network_id,
            status=EXPLORATION_STATUS_SUCCEEDED,
            # get_networks() rather than the adapter's own final_network_paths attribute: the former
            # is the ABC's contract (and enforces the "not before a successful explore()" rule for
            # us), the latter is one concrete adapter's internal state that no other explorer is
            # obliged to have.
            network_paths=tuple(adapter.get_networks()),
            output_paths=tuple(adapter.output_paths),
            k_tp_as_written=tuple(adapter.get_k_tp()),
            # Copied, not aliased. The result advertises itself as a frozen record of what happened;
            # sharing the adapter's live dict would let anything still holding the adapter rewrite
            # the provenance of a run that has already been reported.
            manifest=copy.deepcopy(adapter.manifest),
            selection=selection,
        )
    # An adapter that fails without saying why is violating the contract documented on
    # PESExplorerAdapter.reasons. Say exactly that, rather than letting PDepExplorationResult's
    # empty-reasons guard raise a ValueError that would read as though this function had built a
    # malformed result -- the fact worth reporting is WHICH adapter broke the contract, and the
    # underlying failure must not be swallowed on the way past.
    reasons = tuple(adapter.reasons) or (
        f'The {type(adapter).__name__} explorer reported failure without recording any reason. '
        f'This is an adapter contract violation (see PESExplorerAdapter.reasons); the exploration '
        f'genuinely did fail, but no diagnosis is available from it.',)
    return PDepExplorationResult(
        network_id=parsed_network.network_id,
        status=EXPLORATION_STATUS_FAILED,
        reasons=reasons,
        # A failed run's artifacts are exactly what someone diagnosing the failure needs, and this
        # is the only record they get -- 'reasons' says what T3 concluded, while output_paths say
        # where Arkane's own logs and partial output are. Dropping them forced a human to rediscover
        # the run directory by hand. network_paths and k(T,P) are deliberately NOT carried: those
        # would assert a usable result exists, which is the claim the failure denies.
        output_paths=tuple(adapter.output_paths),
        manifest=copy.deepcopy(adapter.manifest),
        selection=selection,
    )


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

    The file is a mapping with a ``selection_schema_version`` marker
    (``t3.pdep.selector.SELECTION_SCHEMA_VERSION``) and a ``selections`` list, rather than a bare
    list, so the on-disk format can evolve (e.g. gain new top-level keys) without becoming
    ambiguous with an old-format file. This marker describes the SHAPE of the envelope and of each
    selection record, not the decision logic that produced them -- see
    ``t3.pdep.selector.SELECTION_SCHEMA_VERSION``'s own comment for why that distinction matters.

    Args:
        path (str): The path to write the YAML file to.
        selections (list): The ``PDepNetworkSelection`` decisions to save.

    Returns:
        str: ``path``, as a string, so callers can chain it.
    """
    content = {'selection_schema_version': SELECTION_SCHEMA_VERSION,
              'selections': [selection.as_dict() for selection in selections],
              }
    save_yaml_file(path=path, content=content)
    return str(path)


def save_pdep_exploration_results(path: str, results: list) -> str:
    """
    Save a list of PDep exploration results to a YAML file.

    The file is a mapping with an ``exploration_result_schema_version`` marker
    (``t3.pdep.explorer.result.EXPLORATION_RESULT_SCHEMA_VERSION``) and a ``results`` list, rather
    than a bare list, so the on-disk format can evolve (e.g. gain new top-level keys) without
    becoming ambiguous with an old-format file -- mirroring ``save_pdep_network_selections`` above.

    The envelope deliberately does NOT carry a selection schema version, even though each result
    nests a serialized selection (``PDepExplorationResult.as_dict()['selection']``). It does not
    need to: ``selection_schema_version`` and ``selection_algorithm_version`` are FIELDS on
    ``PDepNetworkSelection`` (not envelope keys), so every nested selection record already
    self-describes. A version key on this envelope would be a second, redundant source of truth
    that could disagree with the records it describes -- so it is deliberately omitted.

    Args:
        path (str): The path to write the YAML file to.
        results (list): The ``PDepExplorationResult`` records to save.

    Returns:
        str: ``path``, as a string, so callers can chain it.
    """
    content = {'exploration_result_schema_version': EXPLORATION_RESULT_SCHEMA_VERSION,
              'results': [result.as_dict() for result in results],
              }
    save_yaml_file(path=path, content=content)
    return str(path)
