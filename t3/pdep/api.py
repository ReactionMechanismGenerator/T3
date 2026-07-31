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

import yaml

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
from t3.pdep.parser import PDepArkaneReaction, PDepNetwork, parse_pdep_network_file
from t3.pdep.selector import (CACHE_STATUS_CACHED_REJECTED,
                              CACHE_STATUS_CACHED_VALID,
                              CACHE_STATUS_GENERATED,
                              CACHE_STATUS_UNVALIDATED,
                              E0_PERTURBATION_J_PER_MOL,
                              EVALUATION_STATUS_EVALUATED,
                              EVALUATION_STATUS_NOT_EVALUATED,
                              SELECTION_ALGORITHM_VERSION,
                              SELECTION_SCHEMA_VERSION,
                              STRUCTURES_KEY,
                              PDepNetworkSelection,
                              SensitiveTransitionState,
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
        # method, network_id, the hash diagnosis (no-hash vs. hash-mismatch), and the evaluation-
        # status diagnosis below are each independently sufficient reasons this selection cannot gate
        # an exploration, and each is established by different evidence. Accumulate every applicable
        # diagnosis into ONE raise rather than stopping at the first: a caller who fixes only the
        # first-reported problem, re-runs, and is then handed a second, previously unmentioned
        # problem has been actively misled about how much was wrong.
        reasons = []
        if selection.method != config.method:
            reasons.append(
                f"selection.method ({selection.method!r}) does not match config.method "
                f"({config.method!r}). Gating a decision made under one master-equation method and "
                f"then exploring under another is silent provenance corruption: the recorded "
                f"decision would not be about the run that actually happened.")
        # network_id is a FILE STEM, so a match only proves the decision was about a file with this
        # NAME. RMG rewrites pdep/network4_2.py on every iteration that touches the network, and a
        # selection is routinely made in one process and acted on in another, so "same stem" and
        # "same network" are not the same statement. Bind to the content instead -- but only when
        # network_id itself already matches: a network_id mismatch already proves this is a decision
        # about a DIFFERENT network, so the hash-mismatch wording ("the network file has changed
        # since the decision was made") would be misleading and is short-circuited (skipped) below.
        if selection.network_id != parsed_network.network_id:
            reasons.append(
                f"selection.network_id ({selection.network_id!r}) does not match the parsed "
                f"network's network_id ({parsed_network.network_id!r}). A decision about a "
                f"different network used as this one's budget gate fabricates confidence in a "
                f"result nothing actually justifies -- the same hazard "
                f"PDepNetworkSelection.combine() already refuses for the same reason.")
        else:
            if selection.network_source_hash is None:
                reasons.append(
                    f"selection.network_source_hash is None, so this decision carries no binding to "
                    f"the content it was made about, and network_id ({selection.network_id!r}) is "
                    f"only a file stem -- it matches every revision of that file. Re-run the "
                    f"selection against the network file (select_pdep_network records the hash "
                    f"whenever it parses one), or pass selection=None to explore without a budget "
                    f"gate.")
            elif selection.network_source_hash != parsed_network.source_hash:
                reasons.append(
                    f"selection.network_source_hash ({selection.network_source_hash!r}) does not "
                    f"match the content hash of {network_path!r} ({parsed_network.source_hash!r}): "
                    f"the network file has changed since the decision was made. The decision "
                    f"describes a network that is no longer the one about to be explored -- its "
                    f"sensitivities, its channels, and the transition states it named may all have "
                    f"moved. Re-run the selection against the current file.")
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
            reasons.append(
                f"selection.evaluation_status is {selection.evaluation_status!r}, so its "
                f"'qualified' field ({selection.qualified!r}) carries no verdict and cannot be used "
                f"as a budget gate -- a decision that was never evaluated is not a decision to not "
                f"explore. Re-run the selection against usable SA data, or pass selection=None to "
                f"explore without gating.")
        if reasons:
            if len(reasons) == 1:
                raise ValueError(reasons[0])
            raise ValueError(
                f"This selection cannot gate an exploration, for {len(reasons)} independent reasons: "
                + ' '.join(f"({i}) {reason}" for i, reason in enumerate(reasons, start=1)))
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


# --- Loaders -------------------------------------------------------------------------------------
#
# These are the read side of ``save_pdep_network_selections``/``save_pdep_exploration_results``
# above. Both are STRICT and carry no migration/fallback path for an older on-disk shape: as of
# this writing ``t3/pdep`` has zero files on ``official/main``, no ``t3_sa_cache.yml`` or selections
# YAML exists anywhere on disk, and neither writer had a caller until these loaders were added -- so
# no file of an older shape can exist to migrate FROM, and a migration branch here would be
# untestable dead code (there is no fixture that could ever legitimately exercise it). If a future
# schema bump needs a migration path, add it then, driven by an actual old file that needs reading.

def _require_record_field(record: dict, key: str, *, path: str, context: str):
    """
    Fetch ``record[key]``, or raise a diagnostic ``ValueError`` if it is absent.

    Every T3-written record always carries every key its ``as_dict()`` renders, so a missing key
    means the record was not written by the matching ``save_*`` function (or was hand-edited/
    truncated). Without this, a missing key would surface as a raw ``KeyError`` out of a public API
    -- naming neither the file, nor which record, nor which field -- instead of a diagnostic
    ``ValueError`` matching the style every other refusal path in these loaders already uses.

    Args:
        record (dict): The record to read from.
        key (str): The field name to fetch.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file (e.g.
            ``'selections[3]'`` or ``'results[2].selection'``), for diagnostics only.

    Raises:
        ValueError: If ``key`` is not present in ``record``.

    Returns:
        The value at ``record[key]``.
    """
    if key not in record:
        raise ValueError(f"PDep file {path!r} has {context} missing required field {key!r}.")
    return record[key]


def _require_list_field(record: dict, key: str, *, path: str, context: str) -> list:
    """
    Fetch ``record[key]`` and require it to already be list-like (a ``list`` or ``tuple``), never a
    bare ``str``/``bytes``.

    ``list(...)``/``tuple(...)`` succeed silently on any iterable, including a string -- so
    ``list(record['warnings'])`` on ``warnings: 'AB'`` would silently coerce it to ``['A', 'B']``,
    character by character, instead of refusing an obviously malformed field. Every record field
    that this module coerces via a bare ``list(...)``/``tuple(...)`` call goes through this guard
    first so that trap is refused instead of silently "succeeding".

    Args:
        record (dict): The record to read from.
        key (str): The field name to fetch.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for diagnostics
            only (see ``_require_record_field``).

    Raises:
        ValueError: If ``key`` is missing, or its value is a ``str``/``bytes`` or not a list/tuple
            at all.

    Returns:
        list: ``record[key]``, as a fresh ``list``.
    """
    value = _require_record_field(record, key, path=path, context=context)
    if isinstance(value, (str, bytes)) or not isinstance(value, (list, tuple)):
        raise ValueError(f"PDep file {path!r} has {context} field {key!r} that must be a list "
                         f"(got {type(value).__name__}: {value!r}); a string is not accepted even "
                         f"though it is iterable, since coercing it would silently shred it "
                         f"character by character.")
    return list(value)


def _require_bool_field(record: dict, key: str, *, path: str, context: str) -> bool:
    """
    Fetch ``record[key]`` and require it to already be a real ``bool``.

    A persisted record is untrusted input (see ``_read_persisted_yaml_file``'s docstring): a
    reconstructed field feeds directly into gates like ``explore_pdep_network``'s
    ``not selection.qualified`` check, so a truthy non-bool (a non-empty string such as ``'yes'``,
    or an int such as ``1``) must be refused rather than silently accepted as "true". ``bool`` is a
    subclass of ``int`` in Python, so ``isinstance(value, bool)`` is checked BEFORE any numeric
    check would otherwise wrongly accept it -- the reverse order would let plain ints through.

    Args:
        record (dict): The record to read from.
        key (str): The field name to fetch.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for diagnostics
            only (see ``_require_record_field``).

    Raises:
        ValueError: If ``key`` is missing, or its value is not a real ``bool``.

    Returns:
        bool: ``record[key]``.
    """
    value = _require_record_field(record, key, path=path, context=context)
    if not isinstance(value, bool):
        raise ValueError(f"PDep file {path!r} has {context} field {key!r} that must be a bool "
                         f"(got {type(value).__name__}: {value!r}).")
    return value


def _require_optional_bool_field(record: dict, key: str, *, path: str, context: str) -> bool | None:
    """
    Fetch ``record[key]`` and require it to be a real ``bool`` or ``None``.

    See ``_require_bool_field`` for why a truthy non-bool must be refused rather than accepted.

    Args:
        record (dict): The record to read from.
        key (str): The field name to fetch.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for diagnostics
            only (see ``_require_record_field``).

    Raises:
        ValueError: If ``key`` is missing, or its value is neither a real ``bool`` nor ``None``.

    Returns:
        bool, optional: ``record[key]``.
    """
    value = _require_record_field(record, key, path=path, context=context)
    if value is not None and not isinstance(value, bool):
        raise ValueError(f"PDep file {path!r} has {context} field {key!r} that must be a bool or "
                         f"null (got {type(value).__name__}: {value!r}).")
    return value


def _require_optional_str_field(record: dict, key: str, *, path: str, context: str) -> str | None:
    """
    Fetch ``record[key]`` and require it to be a ``str`` or ``None``.

    Args:
        record (dict): The record to read from.
        key (str): The field name to fetch.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for diagnostics
            only (see ``_require_record_field``).

    Raises:
        ValueError: If ``key`` is missing, or its value is neither a ``str`` nor ``None``.

    Returns:
        str, optional: ``record[key]``.
    """
    value = _require_record_field(record, key, path=path, context=context)
    if value is not None and not isinstance(value, str):
        raise ValueError(f"PDep file {path!r} has {context} field {key!r} that must be a string or "
                         f"null (got {type(value).__name__}: {value!r}).")
    return value


def _require_str_field(record: dict, key: str, *, path: str, context: str) -> str:
    """
    Fetch ``record[key]`` and require it to be a non-``None`` ``str``.

    Args:
        record (dict): The record to read from.
        key (str): The field name to fetch.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for diagnostics
            only (see ``_require_record_field``).

    Raises:
        ValueError: If ``key`` is missing, or its value is not a ``str``.

    Returns:
        str: ``record[key]``.
    """
    value = _require_record_field(record, key, path=path, context=context)
    if not isinstance(value, str):
        raise ValueError(f"PDep file {path!r} has {context} field {key!r} that must be a string "
                         f"(got {type(value).__name__}: {value!r}).")
    return value


def _require_int_field(record: dict, key: str, *, path: str, context: str) -> int:
    """
    Fetch ``record[key]`` and require it to be a real ``int``, not a ``bool``.

    ``bool`` is a subclass of ``int``, so ``isinstance(value, bool)`` is checked first and refused
    even though ``isinstance(value, int)`` would otherwise accept ``True``/``False`` as ``0``/``1``.

    Args:
        record (dict): The record to read from.
        key (str): The field name to fetch.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for diagnostics
            only (see ``_require_record_field``).

    Raises:
        ValueError: If ``key`` is missing, or its value is not a real ``int``.

    Returns:
        int: ``record[key]``.
    """
    value = _require_record_field(record, key, path=path, context=context)
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"PDep file {path!r} has {context} field {key!r} that must be an int "
                         f"(got {type(value).__name__}: {value!r}).")
    return value


def _require_numeric_field(value, *, path: str, context: str, key: str):
    """
    Require ``value`` to be a real number (``int`` or ``float``), not a ``bool``.

    Used for ``thresholds`` dict values, which are already fetched by the time this is called (the
    dict itself is fetched via ``_require_record_field``), so this takes the value directly rather
    than a ``(record, key)`` pair.

    Args:
        value: The value to check.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for diagnostics
            only (see ``_require_record_field``).
        key (str): The ``thresholds`` sub-key this value came from, for diagnostics only.

    Raises:
        ValueError: If ``value`` is not a real ``int``/``float``.

    Returns:
        The numeric value, unchanged.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"PDep file {path!r} has {context} field 'thresholds' entry {key!r} that "
                         f"must be numeric (got {type(value).__name__}: {value!r}).")
    return value


def _require_enum_field(record: dict, key: str, *, path: str, context: str, allowed: tuple) -> str:
    """
    Fetch ``record[key]`` and require it to be one of ``allowed``.

    Args:
        record (dict): The record to read from.
        key (str): The field name to fetch.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for diagnostics
            only (see ``_require_record_field``).
        allowed (tuple): The values ``record[key]`` is allowed to hold.

    Raises:
        ValueError: If ``key`` is missing, or its value is not one of ``allowed``.

    Returns:
        str: ``record[key]``.
    """
    value = _require_record_field(record, key, path=path, context=context)
    if value not in allowed:
        raise ValueError(f"PDep file {path!r} has {context} field {key!r} with value {value!r}, "
                         f"which is not one of the recognized values {allowed!r}.")
    return value


def _require_optional_enum_field(record: dict, key: str, *, path: str, context: str,
                                 allowed: tuple) -> str | None:
    """
    Fetch ``record[key]`` and require it to be ``None`` or one of ``allowed``.

    Args:
        record (dict): The record to read from.
        key (str): The field name to fetch.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for diagnostics
            only (see ``_require_record_field``).
        allowed (tuple): The values ``record[key]`` is allowed to hold, besides ``None``.

    Raises:
        ValueError: If ``key`` is missing, or its value is neither ``None`` nor one of ``allowed``.

    Returns:
        str, optional: ``record[key]``.
    """
    value = _require_record_field(record, key, path=path, context=context)
    if value is not None and value not in allowed:
        raise ValueError(f"PDep file {path!r} has {context} field {key!r} with value {value!r}, "
                         f"which is neither null nor one of the recognized values {allowed!r}.")
    return value


def _require_thresholds_field(record: dict, key: str, *, path: str, context: str) -> dict:
    """
    Fetch ``record[key]`` and require it to be a mapping whose values are all numeric.

    Deliberately does NOT enforce a fixed/exact key set: ``select_from_sa_dict``'s malformed-
    sa_dict branch records a ``coefficient_floor`` key that is absent from every other thresholds
    dict, so pinning an exact key set here would refuse a legitimately-produced record.

    Args:
        record (dict): The record to read from.
        key (str): The field name to fetch.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for diagnostics
            only (see ``_require_record_field``).

    Raises:
        ValueError: If ``key`` is missing, its value is not a mapping, or any of its values is not
            numeric.

    Returns:
        dict: ``record[key]``, as a fresh ``dict``.
    """
    value = _require_record_field(record, key, path=path, context=context)
    if not isinstance(value, dict):
        raise ValueError(f"PDep file {path!r} has {context} field {key!r} that must be a mapping "
                         f"(got {type(value).__name__}: {value!r}).")
    return {sub_key: _require_numeric_field(sub_value, path=path, context=context, key=sub_key)
           for sub_key, sub_value in value.items()}


def _sensitive_transition_state_from_dict(record: dict, *, path: str,
                                          context: str) -> SensitiveTransitionState:
    """
    Reconstruct one ``SensitiveTransitionState`` from its ``as_dict()`` rendering.

    ``SensitiveTransitionState.as_dict()`` renders the ``condition`` tuple as a
    ``{'T': ..., 'T_unit': ..., 'P': ..., 'P_unit': ...}`` dict when it has the expected 4-element
    shape (see that method's docstring); this restores it to the ``(T, 'K', P, 'bar')`` tuple the
    dataclass field expects. A condition that was instead rendered as a plain list (the malformed-
    condition fallback in ``as_dict()``) is restored as a tuple of its elements, matching whatever
    shape it actually has.

    Args:
        record (dict): One entry of a ``selected_ts``/``uncertain_path_reactions`` list, as rendered
            by ``SensitiveTransitionState.as_dict()``.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file (e.g.
            ``'selections[3].selected_ts[0]'``), for diagnostics only (see
            ``_require_record_field``).

    Raises:
        ValueError: If ``record`` is not a mapping, its ``condition`` field is a string (a non-dict
            ``condition`` is coerced via ``tuple(...)``, which would otherwise silently shred a
            string character by character -- see ``_require_list_field``'s docstring), or any other
            field is missing or the wrong type.

    Returns:
        SensitiveTransitionState: The reconstructed record.
    """
    if not isinstance(record, dict):
        raise ValueError(f"PDep file {path!r} has {context} that is not a mapping "
                         f"(got {type(record).__name__}: {record!r}).")
    condition = _require_record_field(record, 'condition', path=path, context=context)
    if isinstance(condition, dict):
        condition = (condition.get('T'), condition.get('T_unit'), condition.get('P'),
                    condition.get('P_unit'))
    else:
        if isinstance(condition, (str, bytes)) or not isinstance(condition, (list, tuple)):
            raise ValueError(f"PDep file {path!r} has {context} field 'condition' that must be a "
                             f"mapping or a list (got {type(condition).__name__}: {condition!r}); a "
                             f"string is not accepted even though it is iterable, since coercing it "
                             f"would silently shred it character by character.")
        condition = tuple(condition)
    return SensitiveTransitionState(
        ts_label=_require_str_field(record, 'ts_label', path=path, context=context),
        coefficient=_require_numeric_field(
           _require_record_field(record, 'coefficient', path=path, context=context),
           path=path, context=context, key='coefficient'),
        condition=condition,
        path_reaction_label=_require_optional_str_field(record, 'path_reaction_label', path=path,
                                                        context=context),
        path_reaction_str=_require_optional_str_field(record, 'path_reaction_str', path=path,
                                                      context=context),
        kinetics_comment=_require_str_field(record, 'kinetics_comment', path=path, context=context),
        uncertain=_require_optional_bool_field(record, 'uncertain', path=path, context=context),
        delta_ln_k=_require_numeric_field(
           _require_record_field(record, 'delta_ln_k', path=path, context=context),
           path=path, context=context, key='delta_ln_k'),
    )


def _selection_from_dict(record, *, path: str, context: str,
                         allow_none: bool = False) -> PDepNetworkSelection | None:
    """
    Reconstruct one ``PDepNetworkSelection`` from its ``as_dict()`` rendering.

    Shared by ``load_pdep_network_selections`` (for each entry of the ``selections`` list) and
    ``load_pdep_exploration_results`` (for the nested ``selection`` key of each result), so the
    per-record shape and version checks live in exactly one place.

    Args:
        record: One selection record, as rendered by ``PDepNetworkSelection.as_dict()``, or
            ``None``. Only meaningful when ``allow_none`` is true (the nested-selection case, where
            an explicit ``null`` means "no gating decision"); a top-level ``selections`` list entry
            has no such meaning for ``None`` and must never be passed with ``allow_none=True``.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file (e.g.
            ``'selections[3]'`` or ``'results[2].selection'``), for diagnostics only.
        allow_none (bool): Whether ``record is None`` is a legitimate value that should pass through
            as ``None`` rather than being refused. True only for the nested ``selection`` key of an
            exploration result, where ``PDepExplorationResult.as_dict()`` always writes this key
            (rendering it as ``None`` when there is no gating decision -- see that method's
            docstring) -- so an explicit ``null`` here is meaningful, but an ABSENT key is not the
            same thing and must be refused by the caller before this function is ever invoked with
            ``allow_none=True`` (see ``load_pdep_exploration_results``).

    Raises:
        ValueError: If ``record`` is ``None`` and ``allow_none`` is false; if ``record`` is neither
            ``None`` nor a mapping; if it carries a ``selection_schema_version`` this code does not
            understand; if it carries a ``selection_algorithm_version`` this code does not
            understand; or if any required field is missing or the wrong type.

    Returns:
        PDepNetworkSelection, optional: The reconstructed selection, or ``None`` if ``record`` was
            ``None`` and ``allow_none`` is true.
    """
    if record is None:
        if allow_none:
            return None
        raise ValueError(f"PDep file {path!r} has {context} that is null; a selection record "
                         f"cannot be reconstructed from nothing.")
    if not isinstance(record, dict):
        raise ValueError(f"PDep file {path!r} contains a selection record that is not a mapping "
                         f"(got {type(record).__name__}: {record!r}) at {context}.")
    record_version = record.get('selection_schema_version')
    if record_version != SELECTION_SCHEMA_VERSION:
        raise ValueError(f"PDep file {path!r} contains a selection record with "
                         f"selection_schema_version={record_version!r}, which this code does not "
                         f"understand (only {SELECTION_SCHEMA_VERSION} is supported), at {context}.")
    # selection_algorithm_version describes the decision SEMANTICS (which gates were applied and
    # how), not the on-disk shape -- unlike selection_schema_version, a newer algorithm version could
    # plausibly change what 'qualified'/'selected_ts' MEAN for the same shape, so a schema match does
    # not imply an algorithm match. t3.pdep has shipped no release and this constant has never been
    # bumped, so there is no concrete evidence that a future bump would be safe to read forward
    # (additively or otherwise) rather than a genuine semantic break; consistent with this package's
    # existing strict no-migration stance (see the module comment above
    # _sensitive_transition_state_from_dict), refuse outright rather than silently trusting a
    # decision produced by logic this code was never shown to agree with. If a future bump turns out
    # to be safely additive, that should be a deliberate, tested decision made at THAT time, not a
    # default assumed here.
    record_algorithm_version = record.get('selection_algorithm_version')
    if record_algorithm_version != SELECTION_ALGORITHM_VERSION:
        raise ValueError(f"PDep file {path!r} contains a selection record with "
                         f"selection_algorithm_version={record_algorithm_version!r}, which this "
                         f"code cannot interpret (only {SELECTION_ALGORITHM_VERSION} is supported), "
                         f"at {context}.")
    selection = PDepNetworkSelection(
        network_id=_require_str_field(record, 'network_id', path=path, context=context),
        network_source_hash=_require_optional_str_field(record, 'network_source_hash', path=path,
                                                        context=context),
        qualified=_require_bool_field(record, 'qualified', path=path, context=context),
        network_reaction=_require_optional_str_field(record, 'network_reaction', path=path,
                                                      context=context),
        direction_key=_require_optional_str_field(record, 'direction_key', path=path,
                                                  context=context),
        direction_keys=_require_list_field(record, 'direction_keys', path=path, context=context),
        direction_ambiguous=_require_bool_field(record, 'direction_ambiguous', path=path,
                                                context=context),
        method=_require_optional_str_field(record, 'method', path=path, context=context),
        sa_path=_require_optional_str_field(record, 'sa_path', path=path, context=context),
        cache_status=_require_optional_enum_field(
           record, 'cache_status', path=path, context=context,
           allowed=(CACHE_STATUS_GENERATED, CACHE_STATUS_CACHED_VALID, CACHE_STATUS_CACHED_REJECTED,
                   CACHE_STATUS_UNVALIDATED)),
        thresholds=_require_thresholds_field(record, 'thresholds', path=path, context=context),
        selected_ts=[_sensitive_transition_state_from_dict(entry, path=path,
                                                           context=f'{context}.selected_ts[{i}]')
                    for i, entry in enumerate(_require_list_field(record, 'selected_ts', path=path,
                                                                  context=context))],
        uncertain_path_reactions=[
           _sensitive_transition_state_from_dict(
              entry, path=path, context=f'{context}.uncertain_path_reactions[{i}]')
           for i, entry in enumerate(_require_list_field(record, 'uncertain_path_reactions',
                                                         path=path, context=context))],
        warnings=_require_list_field(record, 'warnings', path=path, context=context),
        network_reactions_examined=_require_int_field(record, 'network_reactions_examined',
                                                       path=path, context=context),
        evaluation_status=_require_enum_field(
           record, 'evaluation_status', path=path, context=context,
           allowed=(EVALUATION_STATUS_EVALUATED, EVALUATION_STATUS_NOT_EVALUATED)),
        selection_schema_version=record['selection_schema_version'],
        selection_algorithm_version=record['selection_algorithm_version'],
    )
    _validate_selection_cross_field_invariants(selection, path=path, context=context)
    return selection


def _validate_selection_cross_field_invariants(selection: PDepNetworkSelection, *, path: str,
                                               context: str) -> None:
    """
    Reject a reconstructed ``PDepNetworkSelection`` whose fields, though individually well-typed,
    could not have been produced by any constructor in this package.

    ``_selection_from_dict`` validates each field's TYPE independently; it does not (and, per-field,
    cannot) validate the relationships BETWEEN fields that ``select_from_sa_dict``/``combine()``
    always maintain. A hand-edited record can satisfy every per-field check while still fabricating
    a positive verdict (``qualified=True``) with no evidence behind it, bypassing
    ``explore_pdep_network``'s budget gate. This closes exactly that hole, no more.

    The reachable ``(qualified, evaluation_status, uncertain_path_reactions, selected_ts)``
    combinations, read off ``select_from_sa_dict`` (single decision) and ``combine()``
    (aggregate), are:

    1. Single decision, not evaluated (malformed sa_dict, unresolved direction, malformed SA
       entry, or none/only-unusable/below-floor TS rows): ``qualified=False``,
       ``uncertain_path_reactions=[]``, ``selected_ts=[]`` (the loop that would populate
       ``selected_ts`` either never ran or ran and found nothing usable).
    2. Single decision, evaluated: ``select_from_sa_dict`` sets
       ``uncertain_path_reactions = [e for e in selected_ts if e.uncertain]`` and then
       ``qualified = bool(uncertain_path_reactions)`` -- so for a single decision ``qualified``
       and ``bool(uncertain_path_reactions)`` are always exactly equal, and every entry of
       ``uncertain_path_reactions`` is both an element of ``selected_ts`` and has
       ``uncertain is True``.
    3. Combined (``combine()``): ``uncertain_path_reactions``/``selected_ts`` are each unioned
       (order-preserving, de-duplicated) across ALL components regardless of evaluation status,
       while ``qualified`` is unioned (``any(...)``) over EVALUATED components only. This is why
       ``qualified=False`` with a non-empty ``uncertain_path_reactions`` IS reachable (a
       not-evaluated component's evidence rides along on a negative aggregate) -- so the exact
       equality from case 2 must NOT be enforced on a loaded record. What DOES still hold after
       ``combine()``: (a) if ANY unioned component qualified, that component's own non-empty
       ``uncertain_path_reactions`` are part of the union, so ``qualified=True`` still implies
       ``uncertain_path_reactions`` is non-empty; (b) the union operation only selects EXISTING
       elements, so every unioned ``uncertain_path_reactions`` entry is still an element of the
       unioned ``selected_ts`` with ``uncertain is True``, unchanged from case 2.

    So exactly two invariants survive every reachable path and are enforced here:

    - ``qualified=True`` requires ``uncertain_path_reactions`` non-empty.
    - Every entry of ``uncertain_path_reactions`` has ``uncertain is True`` and is also present in
      ``selected_ts``.

    Deliberately NOT enforced (each would refuse a legitimate, reachable record):

    - ``qualified == bool(uncertain_path_reactions)`` exactly -- broken by case 3's
      ``qualified=False`` with non-empty evidence.
    - Any relationship between ``qualified``/``uncertain_path_reactions`` and
      ``evaluation_status`` -- ``combine()`` produces ``qualified=True`` with
      ``evaluation_status='not_evaluated'`` as a legitimate partial-yes (see
      ``PDepNetworkSelection.reason()``'s docstring and ``explore_pdep_network``'s gate, which
      accepts exactly this state).

    Args:
        selection (PDepNetworkSelection): The freshly reconstructed decision to validate.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for
            diagnostics only (see ``_require_record_field``).

    Raises:
        ValueError: If ``qualified`` is true with an empty ``uncertain_path_reactions``, or if any
            entry of ``uncertain_path_reactions`` is not marked ``uncertain is True``, or is not
            also present in ``selected_ts``.
    """
    if selection.qualified and not selection.uncertain_path_reactions:
        raise ValueError(
            f"PDep file {path!r} has {context} with qualified=True but an empty "
            f"uncertain_path_reactions for network {selection.network_id!r}: no code path in "
            f"t3.pdep.selector produces a positive verdict without at least one uncertain "
            f"transition state behind it, so this record cannot have come from a real selection "
            f"decision and cannot be trusted to gate QM exploration.")
    for index, entry in enumerate(selection.uncertain_path_reactions):
        if entry.uncertain is not True:
            raise ValueError(
                f"PDep file {path!r} has {context}.uncertain_path_reactions[{index}] "
                f"(ts_label={entry.ts_label!r}) with uncertain={entry.uncertain!r} for network "
                f"{selection.network_id!r}: every entry of uncertain_path_reactions is, by "
                f"construction, an entry of selected_ts with uncertain is True, so this record "
                f"cannot have come from a real selection decision.")
        if entry not in selection.selected_ts:
            raise ValueError(
                f"PDep file {path!r} has {context}.uncertain_path_reactions[{index}] "
                f"(ts_label={entry.ts_label!r}) that is not among {context}.selected_ts for "
                f"network {selection.network_id!r}: uncertain_path_reactions is always a subset of "
                f"selected_ts, so evidence absent from selected_ts is fabricated, not real.")


def _read_persisted_yaml_file(path: str):
    """
    Read a T3-written PDep YAML file with ``yaml.safe_load``, NOT ``arc.common.read_yaml_file``.

    ``read_yaml_file`` uses ``yaml.FullLoader``, which constructs Python objects from tags in the
    file. ``t3.pdep.yaml_safe``'s own module docstring already spells out why that is unacceptable
    here: ``t3.pdep.api`` is a PUBLIC entrypoint reading a CALLER-SUPPLIED path, so anyone able to
    influence that path fully controls what a FullLoader read constructs. The loaders below are
    exactly that shape, and they were written against the wrong primitive.

    Plain ``yaml.safe_load`` suffices -- strictly safer even than ``yaml_safe.read_sa_yaml_file``,
    which exists only because Arkane's SA files legitimately need ``!!python/tuple``. Nothing T3
    writes here does: ``as_dict()`` renders tuples as plain lists, which is what the round-trip
    reconstruction then restores. So no tag support is needed at all, and any file containing one
    is not a file T3 wrote.

    Args:
        path (str): The path to the YAML file to read.

    Returns:
        The parsed content, using only plain Python types.

    Raises:
        ValueError: If the file cannot be parsed as safe YAML (e.g. it carries a Python object tag).
    """
    with open(path, 'r') as f:
        try:
            return yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ValueError(f'Could not parse {path!r} as plain YAML: {e}. A T3-written PDep file '
                             f'contains only plain types; a Python object tag here means the file '
                             f'was not written by T3.') from e


def load_pdep_network_selections(path: str) -> list:
    """
    Load a list of PDep network selection decisions from a YAML file written by
    ``save_pdep_network_selections``.

    This is STRICT and carries no migration/fallback path for an older on-disk shape -- see the
    module comment above ``_sensitive_transition_state_from_dict`` for why that is safe. An
    unversioned file, a file whose version this code does not recognize, a malformed envelope, or a
    malformed record are all refused outright rather than guessed at.

    Args:
        path (str): The path to read the YAML file from.

    Raises:
        ValueError: If the file's top level is not a mapping; if the envelope carries no
            ``selection_schema_version`` key (an unversioned file is not "version 1", it is of
            unknown shape); if that version is not the one this code understands; if the
            ``selections`` key is absent or not a list; or if any entry (identified by its index in
            the list) is ``None``, is not a mapping, is missing a required field, has a required
            field of the wrong type (including a list-typed field given a bare string), carries a
            ``selection_schema_version`` this code does not understand, or carries a
            ``selection_algorithm_version`` this code does not understand.

    Returns:
        list: The reconstructed ``PDepNetworkSelection`` decisions, in file order.
    """
    content = _read_persisted_yaml_file(path=path)
    if not isinstance(content, dict):
        raise ValueError(f"PDep network selections file {path!r} does not contain a mapping at its "
                         f"top level (got {type(content).__name__}); cannot read it as a selections "
                         f"envelope.")
    if 'selection_schema_version' not in content:
        raise ValueError(f"PDep network selections file {path!r} has no 'selection_schema_version' "
                         f"key: an unversioned file is not version {SELECTION_SCHEMA_VERSION}, it is "
                         f"of unknown shape and cannot be trusted.")
    envelope_version = content['selection_schema_version']
    if envelope_version != SELECTION_SCHEMA_VERSION:
        raise ValueError(f"PDep network selections file {path!r} has "
                         f"selection_schema_version={envelope_version!r}, but this code only "
                         f"understands version {SELECTION_SCHEMA_VERSION}.")
    records = content.get('selections')
    if not isinstance(records, list):
        raise ValueError(f"PDep network selections file {path!r} has no 'selections' list "
                         f"(got {type(records).__name__ if 'selections' in content else 'missing'}: "
                         f"{records!r}).")
    return [_selection_from_dict(record, path=path, context=f'selections[{index}]',
                                 allow_none=False)
           for index, record in enumerate(records)]


def _arkane_reaction_from_dict(record: dict, *, path: str, context: str) -> PDepArkaneReaction:
    """
    Reconstruct one ``PDepArkaneReaction`` from its ``as_dict()`` rendering.

    Args:
        record (dict): One entry of a ``k_tp_as_written`` list, as rendered by
            ``PDepArkaneReaction.as_dict()``.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file (e.g.
            ``'results[2].k_tp_as_written[0]'``), for diagnostics only.

    Raises:
        ValueError: If ``record`` is not a mapping, is missing a required field, has a list-typed
            field given a bare string (see ``_require_list_field``), ``kinetics_type`` is neither a
            string nor null, or ``kinetics_params`` is not a mapping.

    Returns:
        PDepArkaneReaction: The reconstructed record.
    """
    if not isinstance(record, dict):
        raise ValueError(f"PDep exploration results file {path!r} contains a k_tp_as_written entry "
                         f"that is not a mapping (got {type(record).__name__}: {record!r}) "
                         f"at {context}.")
    kinetics_params = _require_record_field(record, 'kinetics_params', path=path, context=context)
    if not isinstance(kinetics_params, dict):
        raise ValueError(f"PDep file {path!r} has {context} field 'kinetics_params' that must be a "
                         f"mapping (got {type(kinetics_params).__name__}: {kinetics_params!r}).")
    return PDepArkaneReaction(
        reactants=tuple(_require_list_field(record, 'reactants', path=path, context=context)),
        products=tuple(_require_list_field(record, 'products', path=path, context=context)),
        kinetics_type=_require_optional_str_field(record, 'kinetics_type', path=path,
                                                  context=context),
        kinetics_params=kinetics_params,
        numeric_values=tuple(_require_list_field(record, 'numeric_values', path=path,
                                                 context=context)),
        rate_payload_numeric_values=tuple(_require_list_field(record, 'rate_payload_numeric_values',
                                                              path=path, context=context)),
        missing_kinetics_keys=tuple(_require_list_field(record, 'missing_kinetics_keys', path=path,
                                                        context=context)),
    )


def load_pdep_exploration_results(path: str) -> list:
    """
    Load a list of PDep exploration results from a YAML file written by
    ``save_pdep_exploration_results``.

    This is STRICT and carries no migration/fallback path for an older on-disk shape -- see the
    module comment above ``_sensitive_transition_state_from_dict`` for why that is safe. Each
    result's nested ``selection`` (if any) is reconstructed via ``_selection_from_dict``, which
    applies the same per-record checks ``load_pdep_network_selections`` does, so a result carrying a
    selection this code does not understand is refused for exactly the same reason a bare selection
    file would be.

    Args:
        path (str): The path to read the YAML file from.

    Raises:
        ValueError: If the file's top level is not a mapping; if the envelope carries no
            ``exploration_result_schema_version`` key (an unversioned file is not "version 1", it is
            of unknown shape); if that version is not the one this code understands; if the
            ``results`` key is absent or not a list; or if any entry (identified by its index in the
            list) is not a mapping, is missing a required field (including the nested ``'selection'``
            key -- a T3-written record always carries it, as ``None`` when there is no gating
            decision, so an ABSENT key is refused rather than conflated with an explicit ``null``),
            has a required field of the wrong type (including a list-typed field given a bare
            string), or whose nested selection is not a mapping, carries a
            ``selection_schema_version`` this code does not understand, or carries a
            ``selection_algorithm_version`` this code does not understand.

    Returns:
        list: The reconstructed ``PDepExplorationResult`` records, in file order.
    """
    content = _read_persisted_yaml_file(path=path)
    if not isinstance(content, dict):
        raise ValueError(f"PDep exploration results file {path!r} does not contain a mapping at its "
                         f"top level (got {type(content).__name__}); cannot read it as an exploration "
                         f"results envelope.")
    if 'exploration_result_schema_version' not in content:
        raise ValueError(f"PDep exploration results file {path!r} has no "
                         f"'exploration_result_schema_version' key: an unversioned file is not "
                         f"version {EXPLORATION_RESULT_SCHEMA_VERSION}, it is of unknown shape and "
                         f"cannot be trusted.")
    envelope_version = content['exploration_result_schema_version']
    if envelope_version != EXPLORATION_RESULT_SCHEMA_VERSION:
        raise ValueError(f"PDep exploration results file {path!r} has "
                         f"exploration_result_schema_version={envelope_version!r}, but this code "
                         f"only understands version {EXPLORATION_RESULT_SCHEMA_VERSION}.")
    records = content.get('results')
    if not isinstance(records, list):
        raise ValueError(f"PDep exploration results file {path!r} has no 'results' list "
                         f"(got {type(records).__name__ if 'results' in content else 'missing'}: "
                         f"{records!r}).")
    results = list()
    for index, record in enumerate(records):
        context = f'results[{index}]'
        if record is None or not isinstance(record, dict):
            raise ValueError(f"PDep exploration results file {path!r} contains a result record that "
                             f"is not a mapping (got {type(record).__name__}: {record!r}) "
                             f"at {context}.")
        if 'selection' not in record:
            # PDepExplorationResult.as_dict() ALWAYS writes the 'selection' key (as None when there
            # is no gating decision -- see that method's docstring), so an absent key here cannot
            # come from a file save_pdep_exploration_results wrote; conflating "key absent" with
            # "explicit null" would let a malformed/truncated record silently pass as "explored
            # without a gating decision", which is a real, distinct state elsewhere in this code
            # (see PDepExplorationResult's docstring) and must not be guessed at.
            raise ValueError(f"PDep exploration results file {path!r} has {context} missing the "
                             f"required 'selection' key (a T3-written record always carries it, as "
                             f"null when there is no gating decision).")
        results.append(PDepExplorationResult(
            network_id=_require_record_field(record, 'network_id', path=path, context=context),
            status=_require_record_field(record, 'status', path=path, context=context),
            reasons=tuple(_require_list_field(record, 'reasons', path=path, context=context)),
            network_paths=tuple(_require_list_field(record, 'network_paths', path=path,
                                                     context=context)),
            output_paths=tuple(_require_list_field(record, 'output_paths', path=path,
                                                    context=context)),
            k_tp_as_written=tuple(
                _arkane_reaction_from_dict(entry, path=path,
                                           context=f'{context}.k_tp_as_written[{entry_index}]')
                for entry_index, entry in enumerate(
                    _require_list_field(record, 'k_tp_as_written', path=path, context=context))),
            manifest=_require_record_field(record, 'manifest', path=path, context=context),
            selection=_selection_from_dict(record['selection'], path=path,
                                           context=f'{context}.selection', allow_none=True),
        ))
    return results
