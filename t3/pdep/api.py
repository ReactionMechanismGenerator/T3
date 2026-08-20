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
import shutil
import tempfile
from pathlib import Path

import yaml

from arc.common import save_yaml_file, to_yaml

from t3.pdep.budget import (BUDGET_ALGORITHM_VERSION,
                            BUDGET_RECORD_SCHEMA_VERSION,
                            PDepBudgetNetworkOutcome,
                            PDepBudgetRecord,
                            VALID_BUDGET_OUTCOMES,
                            VALID_BUDGET_SKIP_REASON_CODES,
                            )
from t3.pdep.assessment import (ASSESSMENT_ENVELOPE_SCHEMA_VERSION,
                                ASSESSMENT_RECORD_SCHEMA_VERSION,
                                PDepNetworkAssessment,
                                )
from t3.pdep.cache import read_t_grid_clamp_record, validate_sa_cache
from t3.pdep.explorer.config import PDepExplorerConfig, deep_thaw
from t3.pdep.explorer.factory import explorer_factory
from t3.pdep.explorer.result import (ADMISSION_POLICY_CALLER_ADMITTED,
                                     ADMISSION_POLICY_QUALIFIED_SELECTION,
                                     ADMISSION_POLICY_UNGATED,
                                     EXPLORATION_RESULT_SCHEMA_VERSION,
                                     EXPLORATION_STATUS_FAILED,
                                     EXPLORATION_STATUS_SKIPPED,
                                     EXPLORATION_STATUS_SUCCEEDED,
                                     PDepExplorationResult,
                                     )
from t3.pdep.parser import PDepArkaneReaction, PDepNetwork, parse_pdep_network_file
from t3.pdep.reason_codes import VALID_ASSESSMENT_REASON_CODES, VALID_ASSESSMENT_STATUSES
from t3.pdep.selector import (CACHE_STATUS_CACHED_REJECTED,
                              CACHE_STATUS_UNVALIDATED,
                              E0_PERTURBATION_J_PER_MOL,
                              EVALUATION_STATUS_EVALUATED,
                              EVALUATION_STATUS_NOT_EVALUATED,
                              SELECTION_ALGORITHM_VERSION,
                              SELECTION_SCHEMA_VERSION,
                              STRUCTURES_KEY,
                              VALID_CACHE_STATUSES,
                              PDepNetworkSelection,
                              SensitiveTransitionState,
                              select_from_sa_dict,
                              selection_rank_key,
                              validate_selection_thresholds,
                              )
from t3.pdep.yaml_safe import read_sa_yaml_file
from t3.schema import T3Sensitivity
from t3.utils.network_thermo import t_grid_clamp_shape_error

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
    t_grid_clamp = None
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
                    # Deliberately no t_grid_clamp: the field means "the T grid this decision rests
                    # on", and this decision rests on nothing -- the SA was never read. The only
                    # value available here would come from the sidecar validate_sa_cache has just
                    # refused to trust, so reporting it would launder provenance out of a source
                    # declared untrustworthy one line earlier.
                    t_grid_clamp=None,
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
        # Read only now, once the cache this provenance describes has either been vouched for or
        # been explicitly trusted by the caller via validate_cache=False. Reading it earlier meant
        # the cache-rejected placeholder above carried a value from a sidecar that had just been
        # declared untrustworthy. Still best-effort: sa_path may point at a sidecar predating this
        # field, or at SA data produced outside T3 entirely, in which case read_t_grid_clamp_record
        # returns None (unknown provenance) rather than raising -- it must never gate
        # evaluation_status.
        t_grid_clamp = read_t_grid_clamp_record(sa_path)
        sa_dict = read_sa_yaml_file(path=sa_path)

    if network_reaction is not None:
        return select_from_sa_dict(sa_dict=sa_dict, network=parsed_network, network_reaction=network_reaction,
                                   relative_threshold=relative_threshold, min_delta_ln_k=min_delta_ln_k,
                                   perturbation=perturbation, method=method, sa_path=sa_path,
                                   cache_status=cache_status, t_grid_clamp=t_grid_clamp)

    reaction_keys = [key for key in sa_dict.keys() if key != STRUCTURES_KEY and isinstance(key, str)] \
        if isinstance(sa_dict, dict) else list()
    decisions = [select_from_sa_dict(sa_dict=sa_dict, network=parsed_network, network_reaction=key,
                                     relative_threshold=relative_threshold, min_delta_ln_k=min_delta_ln_k,
                                     perturbation=perturbation, method=method, sa_path=sa_path,
                                     cache_status=cache_status, t_grid_clamp=t_grid_clamp)
                 for key in reaction_keys]
    if not decisions:
        selection = PDepNetworkSelection(
            network_id=parsed_network.network_id,
            network_source_hash=parsed_network.source_hash,
            method=method,
            sa_path=sa_path,
            cache_status=cache_status,
            t_grid_clamp=t_grid_clamp,
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


# The admission policies a CALLER may request. A strict subset of
# ``t3.pdep.explorer.result.VALID_ADMISSION_POLICIES``, which additionally contains 'ungated' --
# that one is not a choice, it is what passing no selection means, so explore_pdep_network() derives
# it and refuses it as an argument. Deriving rather than accepting it also makes the two facts
# impossible to contradict: a caller cannot claim 'ungated' while passing a selection.
REQUESTABLE_ADMISSION_POLICIES = (ADMISSION_POLICY_QUALIFIED_SELECTION, ADMISSION_POLICY_CALLER_ADMITTED)


def explore_pdep_network(network_path: str,
                         config: PDepExplorerConfig,
                         selection: PDepNetworkSelection | None = None,
                         logger=None,
                         *,
                         admission_policy: str = ADMISSION_POLICY_QUALIFIED_SELECTION,
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

    ``admission_policy`` addresses the same pressure from the other side, and the two are NOT
    interchangeable -- passing no selection and admitting a network yourself are different claims,
    which is why 'caller_admitted' without a selection is refused rather than treated as a synonym.
    It
    exists because qualification stopped being the only admission authority: a caller can rank the
    whole field of networks and decide what to spend on itself, in which case ``selection.qualified``
    is a TIER that informed the ranking, not a veto over its outcome. (T3's own in-run path ranks
    like this -- ``t3.main`` hands the field to ``t3.pdep.budget.apply_pdep_qm_budget`` -- but still
    offers only QUALIFIED selections to the budget, so it does not itself use this policy. The
    policy is here for callers driving the public API, and for T3 if that ever changes.) Under
    ``admission_policy='caller_admitted'`` the QUALIFICATION checks stand aside -- the unqualified
    skip and the not-evaluated refusal both -- and nothing else does: ``method``, ``network_id`` and
    ``network_source_hash`` are provenance, and no budget decision makes a stale selection current.
    That is the whole difference from ``selection=None``, which drops the binding along with the
    gate; a caller that has admitted a network still wants its run bound to the content the decision
    was made about. The policy is recorded on the returned result, because a ``'succeeded'`` result
    carrying an unqualified selection is otherwise indistinguishable from a bypassed gate.

    A caveat for callers assembling their own admission decision: ``apply_pdep_qm_budget()`` is not
    a drop-in oracle for this parameter. It refuses an all-``not_evaluated`` network outright, and
    an evaluated-but-unqualified selection names no uncertain transition states -- so it costs
    nothing, always "fits", and would be admitted for free. T3's own path never meets that case
    because it only ever offers QUALIFIED selections to the budget; a caller piping the full output
    of ``rank_pdep_networks()`` into it would.

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
            ``False`` (and ``admission_policy`` is the default), the exploration is skipped -- no
            adapter is ever constructed -- and a ``'skipped'`` result carrying ``selection.reason()``
            is returned. When ``None`` (the default), the exploration runs unconditionally.
        logger (Logger, optional): The current T3 Logger instance, passed through to the explorer
            adapter.
        admission_policy (str, optional): Keyword-only. What admits this exploration: one of
            ``ADMISSION_POLICY_QUALIFIED_SELECTION`` (the default -- ``selection`` is the admission
            authority, and an unqualified or unevaluated one declines the run) or
            ``ADMISSION_POLICY_CALLER_ADMITTED`` (the caller already decided to spend on this
            network; ``selection`` is kept for its binding to the network content, not for its
            verdict). See the discussion above. Recorded on the returned result. Keyword-only so
            that it cannot be passed where ``logger`` is expected: ``logger`` is the fourth
            positional parameter and existing callers may pass it positionally.

    Raises:
        ValueError: If ``network_path`` is not a ``str`` (see above); if ``admission_policy`` is not
            one of the ``ADMISSION_POLICY_*`` constants, or is ``'caller_admitted'`` with no
            ``selection`` to bind (that combination is inert, and the likeliest way to write it is
            by forgetting the selection); if ``selection`` is given and
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
        PDepExplorationResult: The outcome -- 'skipped' (the qualification gate declined it), 'failed' (the
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
    # Checked before anything is parsed or opened: a misspelled policy that fell through to the
    # default would silently REINSTATE the qualification gate for a caller who explicitly asked it to
    # stand aside, so a network they ranked and chose to spend on would be skipped and reported as
    # declined -- a wrong statement about a decision the caller actually made.
    if admission_policy not in REQUESTABLE_ADMISSION_POLICIES:
        raise ValueError(
            f'explore_pdep_network() admission_policy must be one of '
            f'{REQUESTABLE_ADMISSION_POLICIES}, got {admission_policy!r}. (The recorded field on '
            f'PDepExplorationResult has one further value, {ADMISSION_POLICY_UNGATED!r}, but that '
            f'is not a policy a caller chooses -- it is what passing no selection MEANS, and is '
            f'derived here rather than requested.)')
    if admission_policy == ADMISSION_POLICY_CALLER_ADMITTED and selection is None:
        # Not an over-refusal: this rejects an incoherent CALL, not a run whose data supports it.
        # 'caller_admitted' exists to keep the provenance binding while overriding the qualification
        # verdict, and with no selection there is no binding to keep, so the argument would do
        # nothing whatsoever. The likeliest way to write this call is by forgetting the selection --
        # which loses exactly the binding the argument was reaching for -- so it is refused rather
        # than honoured as a no-op.
        raise ValueError(
            "explore_pdep_network() was given admission_policy='caller_admitted' but no selection. "
            "That policy only means anything alongside a selection: it says the caller made the "
            "spend decision itself and is passing the selection to bind this run to the content the "
            "decision was made about. With selection=None there is nothing to bind and nothing to "
            "override, so this call would behave identically to the default. Pass the selection, or "
            "drop the admission_policy argument.")
    # What gets RECORDED is derived, not echoed. With no selection there is no decision behind this
    # run at all, so recording the argument's default would have the result assert that a qualified
    # selection admitted an exploration for which no selection ever existed -- a false provenance
    # claim, and precisely the one an auditor of an expensive QM run would lean on.
    recorded_admission_policy = admission_policy if selection is not None else ADMISSION_POLICY_UNGATED

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
                f"({config.method!r}). Binding a decision made under one master-equation method to an "
                f"exploration run under another is silent provenance corruption: the recorded "
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
                f"different network used as this one's recorded decision fabricates confidence in a "
                f"result nothing actually justifies -- the same hazard "
                f"PDepNetworkSelection.combine() already refuses for the same reason.")
        else:
            if selection.network_source_hash is None:
                reasons.append(
                    f"selection.network_source_hash is None, so this decision carries no binding to "
                    f"the content it was made about, and network_id ({selection.network_id!r}) is "
                    f"only a file stem -- it matches every revision of that file. Re-run the "
                    f"selection against the network file (select_pdep_network records the hash "
                    f"whenever it parses one), or pass selection=None to explore with no recorded "
                    f"decision at all. Note that admission_policy='caller_admitted' does NOT "
                    f"stand this check down -- it overrides the qualification verdict, not the "
                    f"binding to content.")
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
        #
        # Both this check and the `not selection.qualified` skip below are QUALIFICATION policy, not
        # integrity, so both stand aside under 'caller_admitted'. The reasoning above is entirely
        # about using `qualified` AS A GATE -- and a caller that has declared it admitted this
        # network elsewhere is not using it as one, so the raise would be answering a question nobody
        # asked. Worse, it would force that caller back to selection=None, throwing away the method /
        # network_id / hash binding just checked, which is the whole reason to pass a selection at
        # all. Under the DEFAULT policy the raise stays: with no external admission there is no
        # positive evidence and no spend decision, so a missing evaluation is still a missing
        # verdict, and reading it as "does not qualify" is still silent corruption.
        # Unconditional, and deliberately separate from the policy-gated verdict check below: an
        # evaluation_status outside the known set is MALFORMED DATA, not a verdict this function may
        # decline to consult. Under 'caller_admitted' the verdict check stands aside, so without this
        # a hand-built selection carrying a typo'd or invented status would reach the explorer and be
        # recorded, unexamined, as the provenance of an expensive run. (The loader already refuses
        # such a record on the way in -- see _require_enum_field on 'evaluation_status' -- so this
        # closes the same hole for a selection built in memory rather than read from disk.)
        if selection.evaluation_status not in (EVALUATION_STATUS_EVALUATED,
                                               EVALUATION_STATUS_NOT_EVALUATED):
            reasons.append(
                f"selection.evaluation_status ({selection.evaluation_status!r}) is not one of "
                f"{(EVALUATION_STATUS_EVALUATED, EVALUATION_STATUS_NOT_EVALUATED)}. A decision whose "
                f"coverage cannot be read is malformed, and no admission policy makes it readable.")
        if admission_policy == ADMISSION_POLICY_QUALIFIED_SELECTION \
                and selection.evaluation_status != EVALUATION_STATUS_EVALUATED \
                and not selection.qualified:
            reasons.append(
                f"selection.evaluation_status is {selection.evaluation_status!r}, so its "
                f"'qualified' field ({selection.qualified!r}) carries no verdict and cannot be used "
                f"as an admission gate -- a decision that was never evaluated is not a decision to not "
                f"explore. Re-run the selection against usable SA data; or, if you made the spend "
                f"decision yourself, pass admission_policy='caller_admitted' to keep this "
                f"selection's binding to the network content while overriding its verdict; or pass "
                f"selection=None to explore with no recorded decision at all.")
        if reasons:
            if len(reasons) == 1:
                raise ValueError(reasons[0])
            raise ValueError(
                f"This selection cannot bind this exploration, for {len(reasons)} independent reasons: "
                + ' '.join(f"({i}) {reason}" for i, reason in enumerate(reasons, start=1)))
        if admission_policy == ADMISSION_POLICY_QUALIFIED_SELECTION and not selection.qualified:
            # Nothing runs: no adapter is constructed, no filesystem state below is touched.
            return PDepExplorationResult(
                network_id=parsed_network.network_id,
                status=EXPLORATION_STATUS_SKIPPED,
                reasons=(selection.reason(),),
                selection=selection,
                admission_policy=recorded_admission_policy,
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
        timeout=config.timeout,
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
            admission_policy=recorded_admission_policy,
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
        admission_policy=recorded_admission_policy,
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

    return sorted(all_selections, key=selection_rank_key)


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


def _refuse_content_that_would_not_parse_back(content, path: str) -> None:
    """
    Refuse to write bytes that ``_read_persisted_yaml_file`` could not even parse.

    That function's docstring states the invariant: "Nothing T3 writes here [needs tag support] ...
    any file containing one is not a file T3 wrote." It is an assertion ABOUT the writers, made on
    the read side, and nothing on the write side enforced it. ``arc.common.save_yaml_file`` renders
    with ``yaml.dump`` and a full representer, so one non-plain value anywhere in a record -- a
    ``pathlib.Path`` handed to a ``str``-annotated field is the realistic way in -- is written out as
    a ``!!python/object/apply:`` tag. The write reports success, and every loader here refuses the
    file from that moment on.

    What this guarantees, and what it does not, is worth being exact about, because the name of the
    failure it prevents is TOTAL LOSS rather than a bad record. ``yaml.safe_load`` fails before any
    per-record check runs, so one bad nested field costs the whole file -- an entire iteration's
    assessments -- and it does so with a parser error that names a YAML tag rather than the field
    that caused it. It does NOT guarantee the loaders will ACCEPT the file: a ``str`` where a list
    belongs is perfectly good YAML, and is refused later, per record, with a message naming the
    field. Every record type written here now type-checks its own TYPED fields in ``__post_init__``,
    so that second class is mostly unreachable through this package's own constructors -- what is
    left, and what this check exists for, is the CONTENTS of the dict fields that carry provenance
    rather than schema (``PDepNetworkSelection.t_grid_clamp`` is the live example), plus every field
    added after this comment was written.

    The check renders through ``to_yaml`` -- the SAME function ``save_yaml_file`` writes with, custom
    string representer included -- and parses THAT back. Checking with ``yaml.safe_dump`` instead
    would be checking bytes other than the ones written, which is the shape of the very bug this is
    here to prevent.

    Args:
        content: The content about to be written.
        path (str): The destination, named in the error so the caller knows which write was refused.

    Raises:
        ValueError: If the rendered YAML could not be parsed back by ``yaml.safe_load``.
    """
    # `yaml.YAMLError` only, deliberately. A `RecursionError` from a self-referential record, or a
    # `TypeError` from a broken `__repr__`, is a defect in the code that built the record rather than
    # a fact about the data -- converting one into a refusal here would label a bug as a data problem
    # and hide it in a provenance record. The same DATA-vs-CODE line the rest of this module draws.
    try:
        yaml.safe_load(to_yaml(py_content=content))
    except yaml.YAMLError as e:
        raise ValueError(f'Refusing to write {path!r}: the rendered YAML does not parse back as plain '
                         f'YAML, so this module could never read the file it is about to write ({e}). '
                         f'A T3-written PDep file contains only plain types -- a str-annotated field '
                         f'holding something else (a pathlib.Path is the usual culprit) is rendered as '
                         f'a Python object tag, and the WHOLE file, not just that field, is refused on '
                         f'the way back in.') from e


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
    _refuse_content_that_would_not_parse_back(content=content, path=path)
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
    _refuse_content_that_would_not_parse_back(content=content, path=path)
    save_yaml_file(path=path, content=content)
    return str(path)


def save_pdep_budget_record(path: str, record: PDepBudgetRecord) -> str:
    """
    Save one PDep QM budget record to a YAML file.

    Unlike ``save_pdep_network_selections``/``save_pdep_exploration_results`` above, the content
    written here is exactly ``record.as_dict()``, with no separate envelope wrapping it. Those two
    savers wrap a LIST in an envelope carrying its own top-level version marker, because a list can
    legitimately be empty -- an empty list has no record of its own to carry a version, so the
    envelope has to. ``PDepBudgetRecord`` is different on both counts: it is saved as a single
    object, never a list, and it already carries its own ``schema_version``/``algorithm_version``
    fields (see ``PDepBudgetRecord.__post_init__``, which pins them to
    ``BUDGET_RECORD_SCHEMA_VERSION``/``BUDGET_ALGORITHM_VERSION``). A version that must survive
    nesting inside another record belongs on the dataclass; a version merely describing the outer
    file's envelope belongs on the envelope. Here there is no outer file distinct from the record
    itself, so adding a second, envelope-level version key would just be a redundant source of truth
    that could disagree with the one the record already carries -- exactly the reasoning
    ``save_pdep_exploration_results`` gives for omitting a redundant selection-version envelope key,
    taken one step further since there is no list to wrap at all.

    Args:
        path (str): The path to write the YAML file to.
        record (PDepBudgetRecord): The budget record to save.

    Returns:
        str: ``path``, as a string, so callers can chain it.
    """
    content = record.as_dict()
    # Checked BEFORE staging, so a refusal cannot leave a `.pdep-budget-*` dropping in the iteration
    # directory.
    _refuse_content_that_would_not_parse_back(content=content, path=path)
    # Write atomically: stage the YAML in the same directory, then ``os.replace`` it onto ``path`` in
    # one filesystem operation, so a crash or full disk mid-write can never leave a truncated,
    # partial-but-still-parseable record in place that a later reader would trust as authoritative.
    # ``os.replace`` also closes the symlink write-through vector: if ``path`` is a pre-existing
    # symlink, it replaces the link itself rather than following it and overwriting whatever the link
    # points to. Mirrors ``t3.pdep.capture._write_manifest``'s mkstemp+save_yaml_file+os.replace
    # idiom, the existing precedent for atomically writing a YAML sidecar in this codebase.
    directory = os.path.dirname(path) or '.'
    fd, staged_path = tempfile.mkstemp(prefix='.pdep-budget-', dir=directory)
    os.close(fd)
    try:
        save_yaml_file(path=staged_path, content=content)
        os.replace(staged_path, path)
    finally:
        if os.path.isfile(staged_path):
            os.remove(staged_path)
    return str(path)


def _fsync_path(path: str) -> None:
    """
    Flush one path to the storage device, so what was written to it survives a power loss.

    Used for both halves of a durable replace: the staged FILE, so its bytes are on disk before
    anything points at them, and then the destination DIRECTORY, so the rename that points at them
    is on disk too. A directory is opened read-only and fsynced exactly like a file on POSIX; the two
    cases differ only in what the kernel flushes.

    Args:
        path (str): The file or directory to flush.
    """
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def save_pdep_network_assessments(path: str, assessments: list, *, complete: bool) -> str:
    """
    Save the PDep network assessment records for one T3 iteration to a YAML file.

    Like ``save_pdep_network_selections``/``save_pdep_exploration_results`` and unlike
    ``save_pdep_budget_record``, the list is wrapped in an envelope. The reason is the same one those
    two give: a list can legitimately be empty, and an iteration in which T3 found no P-dep networks
    at all is an ordinary outcome rather than an error -- an empty list has no record of its own to
    carry a version, so the envelope has to.

    The envelope's version is its OWN constant under its own key, ``assessment_envelope_schema_version``,
    not a second copy of the per-record ``assessment_record_schema_version``. Both are 1 today, and
    ``save_pdep_network_selections`` does reuse one number for both roles, but that shortcut makes
    each version a hostage of the other: adding a field to a record would force the envelope to claim
    a change it never underwent, and renaming the list key would force every record ever written to
    be re-stamped. ``load_pdep_network_assessments`` checks the two separately, so a file cannot
    claim one shape at the top level and hold another underneath.

    The write is atomic AND durable, unlike ``save_pdep_network_selections``'s. That is not a
    stylistic difference: T3's funnel rewrites this file once per network as an iteration progresses,
    so a crash or a full disk part-way through a write is a real scenario rather than a theoretical
    one, and this is precisely the file whose job is to survive the failure that interrupted it. A
    truncated but still-parseable record would be believed by the next reader with the full authority
    of provenance -- and would under-report exactly the networks whose absence this whole increment
    exists to fix.

    That same rewriting is why the envelope states whether the list is FINISHED. A file holding four
    of an iteration's twelve networks is not a smaller version of the truth, it is a different claim,
    and one that reads as authoritative: "four networks were assessed" and "four networks had been
    assessed when T3 died" are the same bytes without this flag. ``complete`` is required rather than
    defaulted precisely because the safe value depends on the caller and a forgotten argument would
    quietly assert the stronger of the two.

    Args:
        path (str): The path to write the YAML file to.
        assessments (list): The ``PDepNetworkAssessment`` records to save.
        complete (bool): Whether these are ALL the assessments of the iteration. ``False`` for the
            incremental writes T3 makes as it goes, ``True`` for the final one once every network
            has had its turn.

    Raises:
        ValueError: If ``complete`` is not a ``bool``. A truthy stand-in would be written out as
            itself and refused on the way back in, which is a confusing place to find out.

    Returns:
        str: ``path``, as a string, so callers can chain it.
    """
    if not isinstance(complete, bool):
        raise ValueError(f'The PDep network assessments `complete` flag must be a bool, got '
                         f'{complete!r} ({type(complete).__name__}). This flag is the difference '
                         f'between a finished record and a crash scene; it cannot be inferred.')
    content = {'assessment_envelope_schema_version': ASSESSMENT_ENVELOPE_SCHEMA_VERSION,
              'complete': complete,
              'assessments': [assessment.as_dict() for assessment in assessments],
              }
    # Checked before staging, for the same reason as `save_pdep_budget_record`: a refusal must not
    # leave a staging directory behind. It matters more here, where this file is rewritten once per
    # network -- a per-network dropping would accumulate across the whole iteration.
    _refuse_content_that_would_not_parse_back(content=content, path=path)
    # Stage in a private directory alongside the destination, then ``os.replace`` in one filesystem
    # operation. Staging in a fresh 0700 directory rather than beside the target (as
    # ``save_pdep_budget_record`` does) closes the gap between creating the staged file and reopening
    # it by path to write: nothing else can substitute a file at a path only this process can reach.
    # ``os.replace`` onto the destination is what makes the swap atomic, and it also means a symlink
    # sitting at ``path`` is replaced rather than written through.
    directory = os.path.dirname(path) or '.'
    staging_directory = tempfile.mkdtemp(prefix='.pdep-assessments-', dir=directory)
    staged_path = os.path.join(staging_directory, 'assessments.yml')
    try:
        save_yaml_file(path=staged_path, content=content)
        # Atomic is not the same as durable: ``os.replace`` orders the rename, but says nothing about
        # whether the bytes behind it ever reached the disk. Without these, a power loss can leave
        # the rename applied and the contents lost -- the file present, and empty or torn.
        _fsync_path(staged_path)
        os.replace(staged_path, path)
        _fsync_path(directory)
    finally:
        shutil.rmtree(staging_directory, ignore_errors=True)
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


def _require_optional_int_field(record: dict, key: str, *, path: str, context: str) -> int | None:
    """
    Fetch ``record[key]`` and require it to be a real ``int`` or ``None``.

    ``bool`` is a subclass of ``int``, so ``isinstance(value, bool)`` is checked first and refused
    even though ``isinstance(value, int)`` would otherwise accept ``True``/``False`` as ``0``/``1``
    -- mirroring ``_require_int_field`` above.

    Args:
        record (dict): The record to read from.
        key (str): The field name to fetch.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for diagnostics
            only (see ``_require_record_field``).

    Raises:
        ValueError: If ``key`` is missing, or its value is neither a real ``int`` nor ``None``.

    Returns:
        int, optional: ``record[key]``.
    """
    value = _require_record_field(record, key, path=path, context=context)
    if value is not None and (isinstance(value, bool) or not isinstance(value, int)):
        raise ValueError(f"PDep file {path!r} has {context} field {key!r} that must be an int or "
                         f"null (got {type(value).__name__}: {value!r}).")
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


def _require_optional_dict_field(record: dict, key: str, *, path: str, context: str,
                                 shape_error=None) -> dict | None:
    """
    Fetch ``record.get(key)`` and require it to be a mapping or ``None`` -- and, unlike every other
    ``_require_*_field`` helper, tolerate the key being ABSENT entirely rather than refusing.

    This is deliberately looser than ``_require_thresholds_field``/``_require_record_field``: it
    backs ``PDepNetworkSelection.t_grid_clamp``, whose whole point is a three-state design (see that
    field's docstring) where "the key was never written" (an old sidecar/selection predating this
    field) must read as unknown provenance (``None``), the SAME as an explicit ``null`` -- never as a
    refusal. A record that DOES carry the key but with a non-mapping, non-null value is still a
    genuine malformation and is refused, since ``TGridClampRecord.as_dict()`` never renders anything
    else.

    That last argument reaches further than "is it a mapping", which is why ``shape_error`` exists:
    ``as_dict()`` does not render a mapping that fails to say whether a clamp happened either. The
    caller supplies the shape rather than this helper knowing it, so the helper stays a helper. The
    record's own ``validate()`` refuses the same thing a moment later when it is constructed, so
    checking here is about diagnostics, not about whether the refusal happens: the constructor knows
    only the network, while a human fixing a hand-edited file needs the file and the position in it.

    Args:
        record (dict): The record to read from.
        key (str): The field name to fetch.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file, for diagnostics
            only (see ``_require_record_field``).
        shape_error (callable, optional): A function taking the mapping and returning a reason string
            when it is not the shape this field is documented to hold, or ``None`` when it is.
            Omitted for a field whose contents genuinely are free-form.

    Raises:
        ValueError: If the key is present with a value that is neither a mapping nor ``None``, or
            with a mapping ``shape_error`` rejects.

    Returns:
        dict, optional: ``record[key]`` as a fresh ``dict``, or ``None`` if the key is absent or
            explicitly ``null``.
    """
    if key not in record:
        return None
    value = record[key]
    if value is not None and not isinstance(value, dict):
        raise ValueError(f"PDep file {path!r} has {context} field {key!r} that must be a mapping or "
                         f"null (got {type(value).__name__}: {value!r}).")
    if value is not None and shape_error is not None:
        reason = shape_error(value)
        if reason is not None:
            raise ValueError(f"PDep file {path!r} has {context} field {key!r} that is a mapping but "
                             f"not the one this field is documented to hold: {reason}.")
    return dict(value) if value is not None else None


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
    # ``isinstance(True, int)`` and ``True == 1``, so a ``selection_schema_version: true`` file would
    # otherwise sail through this equality check and be read as correctly versioned.
    if isinstance(record_version, bool) or record_version != SELECTION_SCHEMA_VERSION:
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
    if isinstance(record_algorithm_version, bool) or record_algorithm_version != SELECTION_ALGORITHM_VERSION:
        raise ValueError(f"PDep file {path!r} contains a selection record with "
                         f"selection_algorithm_version={record_algorithm_version!r}, which this "
                         f"code cannot interpret (only {SELECTION_ALGORITHM_VERSION} is supported), "
                         f"at {context}.")
    selection = PDepNetworkSelection(
        # Optional, deliberately, and closing a real defect rather than a loosening for its own
        # sake: `rank_pdep_networks` records a decision for an entry too malformed to name a network
        # (`PDepNetworkSelection(network_id=network_id_hint)` with a `None` hint), and
        # `t3.pdep.budget` counts two such records as two DISTINCT networks rather than collapsing
        # them. Requiring a string here meant that decision could be written by one public API
        # function and refused by another on the way back in -- taking the whole file with it, not
        # just the unnamed record. `PDepNetworkAssessment.network_id` stays required: an assessment
        # is a statement ABOUT a named network, and there is no unnamed case to record.
        network_id=_require_optional_str_field(record, 'network_id', path=path, context=context),
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
        # `VALID_CACHE_STATUSES`, not a tuple respelled here: the constructor checks the same set,
        # and two hand-written copies of it are exactly how a status added to one half goes missing
        # from the other.
        cache_status=_require_optional_enum_field(record, 'cache_status', path=path, context=context,
                                                  allowed=VALID_CACHE_STATUSES),
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
        t_grid_clamp=_require_optional_dict_field(record, 't_grid_clamp', path=path, context=context,
                                                  shape_error=t_grid_clamp_shape_error),
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
    ``explore_pdep_network``'s qualification gate. This closes exactly that hole, no more.

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
      ``PDepNetworkSelection.reason()``'s docstring and ``explore_pdep_network``'s qualification
      gate, which accepts exactly this state).

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
            # Absent key -> DERIVED, not defaulted and not refused. Every record written before this
            # field existed predates 'caller_admitted' entirely, so for those records the basis is
            # not a guess: a record carrying a selection was necessarily admitted by the
            # qualification gate, and one carrying none was necessarily ungated. That is the same
            # derivation explore_pdep_network() itself performs, applied to the same evidence, so it
            # reconstructs the true value rather than inventing a plausible one.
            #
            # Refusing instead would make every such file unloadable while still calling it
            # exploration_result_schema_version 1 -- two incompatible shapes under one version
            # number, which is precisely what that version exists to prevent. A blanket
            # ``.get('admission_policy', ADMISSION_POLICY_QUALIFIED_SELECTION)`` would be worse than
            # either: it would relabel selection-less records as gate-admitted, the exact false
            # provenance claim this field was added to stop.
            admission_policy=record['admission_policy'] if 'admission_policy' in record
            else (ADMISSION_POLICY_QUALIFIED_SELECTION if record['selection'] is not None
                  else ADMISSION_POLICY_UNGATED),
        ))
    return results


def _budget_network_outcome_from_dict(record: dict, *, path: str, context: str) -> PDepBudgetNetworkOutcome:
    """
    Reconstruct one ``PDepBudgetNetworkOutcome`` from its ``as_dict()`` rendering.

    There is no per-outcome schema/algorithm version to check here: unlike a top-level
    ``PDepBudgetRecord`` or a nested ``PDepNetworkSelection``, ``PDepBudgetNetworkOutcome`` carries
    no version fields of its own -- it is only ever produced and consumed alongside the
    ``PDepBudgetRecord`` that owns it, whose ``schema_version`` already governs the shape of every
    outcome nested inside it.

    Args:
        record (dict): One entry of a ``network_outcomes`` list, as rendered by
            ``PDepBudgetNetworkOutcome.as_dict()``.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file (e.g.
            ``'network_outcomes[0]'``), for diagnostics only.

    Raises:
        ValueError: If ``record`` is not a mapping, is missing a required field, has a required field
            of the wrong type, ``outcome`` is not one of ``VALID_BUDGET_OUTCOMES``, or ``reason_code``
            is neither ``None`` nor one of ``VALID_BUDGET_SKIP_REASON_CODES``.

    Returns:
        PDepBudgetNetworkOutcome: The reconstructed outcome.
    """
    if not isinstance(record, dict):
        raise ValueError(f"PDep budget record file {path!r} contains a network outcome that is not "
                         f"a mapping (got {type(record).__name__}: {record!r}) at {context}.")
    return PDepBudgetNetworkOutcome(
        network_id=_require_optional_str_field(record, 'network_id', path=path, context=context),
        outcome=_require_enum_field(record, 'outcome', path=path, context=context,
                                    allowed=VALID_BUDGET_OUTCOMES),
        cost=_require_int_field(record, 'cost', path=path, context=context),
        network_source_hash=_require_optional_str_field(record, 'network_source_hash', path=path,
                                                         context=context),
        method=_require_optional_str_field(record, 'method', path=path, context=context),
        reason_code=_require_optional_enum_field(record, 'reason_code', path=path, context=context,
                                                 allowed=VALID_BUDGET_SKIP_REASON_CODES),
        reason=_require_optional_str_field(record, 'reason', path=path, context=context),
        remaining_transition_states=_require_optional_int_field(record, 'remaining_transition_states',
                                                                path=path, context=context),
        rank=_require_int_field(record, 'rank', path=path, context=context),
        unnamed_offer_index=_require_optional_int_field(record, 'unnamed_offer_index', path=path,
                                                        context=context),
    )


def load_pdep_budget_record(path: str) -> PDepBudgetRecord:
    """
    Load one PDep QM budget record from a YAML file written by ``save_pdep_budget_record``.

    This is STRICT and carries no migration/fallback path for an older on-disk shape -- see the
    module comment above ``_sensitive_transition_state_from_dict`` for why that is safe. An
    unversioned file, a file whose version this code does not recognize, a malformed envelope, or a
    malformed record are all refused outright rather than guessed at.

    Unlike ``load_pdep_network_selections``/``load_pdep_exploration_results``, there is no separate
    envelope to check: ``save_pdep_budget_record`` writes ``record.as_dict()`` directly, so the
    file's top level IS the record, and ``budget_record_schema_version``/``budget_algorithm_version``
    are read as ordinary record fields rather than envelope keys.

    Args:
        path (str): The path to read the YAML file from.

    Raises:
        ValueError: If the file's top level is not a mapping; if it carries no
            ``budget_record_schema_version`` key (an unversioned file is not "version 1", it is of
            unknown shape); if that version is not the one this code understands; if it carries a
            ``budget_algorithm_version`` this code does not understand; if the ``network_outcomes``
            key is absent or not a list; or if any entry (identified by its index in the list) is not
            a mapping, is missing a required field, or has a required field of the wrong type.

    Returns:
        PDepBudgetRecord: The reconstructed budget record.
    """
    content = _read_persisted_yaml_file(path=path)
    if not isinstance(content, dict):
        raise ValueError(f"PDep budget record file {path!r} does not contain a mapping at its top "
                         f"level (got {type(content).__name__}); cannot read it as a budget record.")
    if 'budget_record_schema_version' not in content:
        raise ValueError(f"PDep budget record file {path!r} has no 'budget_record_schema_version' "
                         f"key: an unversioned file is not version {BUDGET_RECORD_SCHEMA_VERSION}, it "
                         f"is of unknown shape and cannot be trusted.")
    schema_version = content['budget_record_schema_version']
    if isinstance(schema_version, bool) or schema_version != BUDGET_RECORD_SCHEMA_VERSION:
        raise ValueError(f"PDep budget record file {path!r} has "
                         f"budget_record_schema_version={schema_version!r}, but this code only "
                         f"understands version {BUDGET_RECORD_SCHEMA_VERSION}.")
    algorithm_version = content.get('budget_algorithm_version')
    if isinstance(algorithm_version, bool) or algorithm_version != BUDGET_ALGORITHM_VERSION:
        raise ValueError(f"PDep budget record file {path!r} has "
                         f"budget_algorithm_version={algorithm_version!r}, which this code cannot "
                         f"interpret (only {BUDGET_ALGORITHM_VERSION} is supported).")
    outcomes = content.get('network_outcomes')
    if not isinstance(outcomes, list):
        raise ValueError(f"PDep budget record file {path!r} has no 'network_outcomes' list "
                         f"(got {type(outcomes).__name__ if 'network_outcomes' in content else 'missing'}: "
                         f"{outcomes!r}).")
    return PDepBudgetRecord(
        iteration=_require_int_field(content, 'iteration', path=path, context='<top level>'),
        max_transition_states=_require_optional_int_field(content, 'max_transition_states', path=path,
                                                          context='<top level>'),
        max_networks=_require_optional_int_field(content, 'max_networks', path=path,
                                                 context='<top level>'),
        total_cost=_require_int_field(content, 'total_cost', path=path, context='<top level>'),
        network_outcomes=tuple(
            _budget_network_outcome_from_dict(entry, path=path, context=f'network_outcomes[{index}]')
            for index, entry in enumerate(outcomes)),
        schema_version=schema_version,
        algorithm_version=algorithm_version,
    )


def _assessment_from_dict(record, *, path: str, context: str) -> PDepNetworkAssessment:
    """
    Reconstruct one ``PDepNetworkAssessment`` from its ``as_dict()`` rendering.

    Unlike ``_selection_from_dict``, this needs no companion cross-field validator: every invariant
    that makes an assessment record trustworthy (a status agreeing with its reason code, a nested
    selection being present exactly when the reason code implies one, a nested selection's own
    ``qualified`` agreeing with the outer verdict) is enforced by ``PDepNetworkAssessment``'s
    ``__post_init__``, so simply constructing one re-checks the whole record. That is the point of
    having built the invariants into the type: a hand-edited file gets refused by the same rule that
    refused the impossible combination at the site that would have written it.

    The construction is wrapped only so those refusals name the FILE and the entry. The type's own
    message names the network, which is enough at the write site but not at the read site: on a run
    with a dozen iterations on disk, "network4_2 carries status qualified with reason
    sa_all_methods_failed" is the difference between a fixable report and a hunt. The per-field
    ``_require_*`` helpers already name both, so they are deliberately left outside the wrapper
    rather than being given a second, redundant prefix.

    Args:
        record: One assessment record, as rendered by ``PDepNetworkAssessment.as_dict()``.
        path (str): The file this record was read from, for diagnostics only.
        context (str): A short description of where this record sits in the file (e.g.
            ``'assessments[3]'``), for diagnostics only.

    Raises:
        ValueError: If ``record`` is not a mapping; if it carries an
            ``assessment_record_schema_version`` this code does not understand; if the ``selection``
            key is absent; if a required field is missing or of the wrong type; if a secondary reason
            code is not recognized; or if the record violates any of the type's own invariants.

    Returns:
        PDepNetworkAssessment: The reconstructed record.
    """
    if not isinstance(record, dict):
        raise ValueError(f"PDep network assessments file {path!r} has {context} that is not a "
                         f"mapping (got {type(record).__name__}: {record!r}); an assessment record "
                         f"cannot be reconstructed from it.")
    record_version = record.get('assessment_record_schema_version')
    if isinstance(record_version, bool) or record_version != ASSESSMENT_RECORD_SCHEMA_VERSION:
        raise ValueError(f"PDep network assessments file {path!r} has {context} with "
                         f"assessment_record_schema_version={record_version!r}, which this code does "
                         f"not understand (only {ASSESSMENT_RECORD_SCHEMA_VERSION} is supported).")
    if 'selection' not in record:
        # PDepNetworkAssessment.as_dict() ALWAYS writes this key, rendering it null for the reason
        # codes that forbid a nested selection. An absent key is therefore not the same thing as an
        # explicit null: null means "there was no selection", absent means this record was not
        # written by T3 -- and guessing null for it would silently convert a corrupt record into a
        # plausible "never evaluated" one, which is the exact misreading this file exists to prevent.
        raise ValueError(f"PDep network assessments file {path!r} has {context} missing its required "
                         f"'selection' key (a T3-written record always carries it, as null when the "
                         f"reason code forbids a nested selection); an absent key is not the same as "
                         f"an explicit null and is not guessed at.")
    secondary_reason_codes = _require_list_field(record, 'secondary_reason_codes', path=path,
                                                 context=context)
    for index, code in enumerate(secondary_reason_codes):
        if code not in VALID_ASSESSMENT_REASON_CODES:
            raise ValueError(f"PDep network assessments file {path!r} has {context} field "
                             f"'secondary_reason_codes'[{index}] with value {code!r}, which is not "
                             f"one of the recognized reason codes "
                             f"{VALID_ASSESSMENT_REASON_CODES!r}.")
    fields = dict(
        network_id=_require_str_field(record, 'network_id', path=path, context=context),
        iteration=_require_int_field(record, 'iteration', path=path, context=context),
        status=_require_enum_field(record, 'status', path=path, context=context,
                                   allowed=VALID_ASSESSMENT_STATUSES),
        reason_code=_require_enum_field(record, 'reason_code', path=path, context=context,
                                        allowed=VALID_ASSESSMENT_REASON_CODES),
        secondary_reason_codes=secondary_reason_codes,
        network_path=_require_optional_str_field(record, 'network_path', path=path, context=context),
        network_source_hash=_require_optional_str_field(record, 'network_source_hash', path=path,
                                                        context=context),
        observable_label=_require_optional_str_field(record, 'observable_label', path=path,
                                                     context=context),
        sa_rank_index=_require_optional_int_field(record, 'sa_rank_index', path=path,
                                                  context=context),
        chemkin_reaction=_require_optional_str_field(record, 'chemkin_reaction', path=path,
                                                     context=context),
        network_reaction=_require_optional_str_field(record, 'network_reaction', path=path,
                                                     context=context),
        requested_me_methods=_require_list_field(record, 'requested_me_methods', path=path,
                                                 context=context),
        final_method=_require_optional_str_field(record, 'final_method', path=path, context=context),
        sa_path=_require_optional_str_field(record, 'sa_path', path=path, context=context),
        # Read as a free optional string rather than through _require_optional_enum_field, even
        # though the selection loader enum-restricts its own cache_status. PDepNetworkAssessment
        # accepts any string here, and a loader stricter than the constructor it feeds would make a
        # record T3 legitimately wrote unreadable by T3 -- the one failure a durable record must
        # never have. Strictness that is not shared with the writer is not safety, it is a bug with a
        # delay on it.
        cache_status=_require_optional_str_field(record, 'cache_status', path=path, context=context),
        warnings=_require_list_field(record, 'warnings', path=path, context=context),
        selection=_selection_from_dict(record['selection'], path=path,
                                       context=f'{context}.selection', allow_none=True),
        schema_version=record_version,
    )
    try:
        return PDepNetworkAssessment(**fields)
    except ValueError as e:
        raise ValueError(f"PDep network assessments file {path!r} has {context} that is not a valid "
                         f"assessment record: {e}") from e


def load_pdep_network_assessments(path: str, *, allow_incomplete: bool = False) -> list:
    """
    Load the PDep network assessment records for one T3 iteration from a YAML file written by
    ``save_pdep_network_assessments``.

    This is STRICT and carries no migration/fallback path for an older on-disk shape -- see the
    module comment above ``_sensitive_transition_state_from_dict`` for why that is safe. An
    unversioned file, a file whose version this code does not recognize, a malformed envelope, or a
    malformed record are all refused outright rather than guessed at.

    An UNFINISHED file is refused too, by default. T3 rewrites this file after every network, so a
    file left behind by a crash is well-formed, parseable, and short -- the one failure this loader
    cannot detect by looking at the records, and the one most likely to be believed. Refusing it
    means a caller counting "how many networks could not be evaluated" cannot silently answer with a
    number that stopped where the crash did. Pass ``allow_incomplete=True`` to read the partial
    record deliberately, which is exactly what an operator investigating that crash wants.

    Args:
        path (str): The path to read the YAML file from.
        allow_incomplete (bool, optional): Whether to accept a file whose envelope says the
            iteration never finished writing it. Must be an exact ``bool``, for the same reason
            ``complete`` must be: an escape hatch that opens on truthiness is not a decision.

    Raises:
        ValueError: If ``allow_incomplete`` is not a ``bool``; if the file's top level is not a
            mapping; if the envelope carries no
            ``assessment_record_schema_version`` key (an unversioned file is not "version 1", it is
            of unknown shape); if that version is not the one this code understands; if the
            ``complete`` flag is absent or not a bool; if it is ``False`` and ``allow_incomplete``
            was not asked for; if the ``assessments`` key is absent or not a list; or if any entry
            (identified by its index in the list) is malformed in any of the ways
            ``_assessment_from_dict`` refuses.

    Returns:
        list: The reconstructed ``PDepNetworkAssessment`` records, in file order.
    """
    content = _read_persisted_yaml_file(path=path)
    if not isinstance(content, dict):
        raise ValueError(f"PDep network assessments file {path!r} does not contain a mapping at its "
                         f"top level (got {type(content).__name__}); cannot read it as an "
                         f"assessments envelope.")
    if 'assessment_envelope_schema_version' not in content:
        raise ValueError(f"PDep network assessments file {path!r} has no "
                         f"'assessment_envelope_schema_version' key: an unversioned file is not "
                         f"version {ASSESSMENT_ENVELOPE_SCHEMA_VERSION}, it is of unknown shape and "
                         f"cannot be trusted.")
    envelope_version = content['assessment_envelope_schema_version']
    if isinstance(envelope_version, bool) or envelope_version != ASSESSMENT_ENVELOPE_SCHEMA_VERSION:
        raise ValueError(f"PDep network assessments file {path!r} has "
                         f"assessment_envelope_schema_version={envelope_version!r}, but this code "
                         f"only understands version {ASSESSMENT_ENVELOPE_SCHEMA_VERSION}.")
    if not isinstance(allow_incomplete, bool):
        raise ValueError(f'The PDep network assessments `allow_incomplete` flag must be a bool, got '
                         f'{allow_incomplete!r} ({type(allow_incomplete).__name__}). It is held to the '
                         f'same standard as the flag it overrides: an escape hatch that opens on '
                         f'truthiness is not a decision, it is an accident.')
    if 'complete' not in content:
        raise ValueError(f"PDep network assessments file {path!r} has no 'complete' key, so there is "
                         f"no way to tell a finished iteration's record from one a crash cut short. "
                         f"Both are well-formed and only one is the whole field.")
    complete = content['complete']
    if not isinstance(complete, bool):
        raise ValueError(f"PDep network assessments file {path!r} has complete={complete!r} "
                         f"({type(complete).__name__}); this must be a bool, since anything else "
                         f"would be read for its truthiness and quietly promote an unfinished record.")
    # Checked BEFORE the completeness refusal below, which reports how many records the partial file
    # got to: counting them first would mean a hand-written `assessments: 1` raised a bare TypeError
    # from inside the message of the error it was about to raise.
    records = content.get('assessments')
    if not isinstance(records, list):
        raise ValueError(f"PDep network assessments file {path!r} has no 'assessments' list "
                         f"(got {type(records).__name__ if 'assessments' in content else 'missing'}: "
                         f"{records!r}).")
    if not complete and not allow_incomplete:
        raise ValueError(f"PDep network assessments file {path!r} says complete=False: T3 was still "
                         f"writing it when the iteration ended, so the {len(records)} record(s) in it "
                         f"are where the run stopped, not every network it looked at. Pass "
                         f"allow_incomplete=True to read it as the partial record it is.")
    return [_assessment_from_dict(record, path=path, context=f'assessments[{index}]')
           for index, record in enumerate(records)]
