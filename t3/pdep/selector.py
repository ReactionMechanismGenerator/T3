"""
t3 pdep selector module

The pure selection core deciding whether a pressure-dependent (PDep) reaction network
qualifies for expensive QM refinement (PES exploration and/or master-equation refinement).

A network qualifies iff both:

(a) **Sensitivity** -- some run-defined observable is sensitive to a reaction in that network.
    This criterion is established by the *caller* (T3's top-SA machinery); this module is only
    ever asked about a network reaction that already satisfies it.
(b) **Uncertainty** -- among the transition states that the network's k(T,P) is in turn most
    sensitive to, at least one belongs to a path reaction whose kinetics is an estimate rather
    than a library value or a training reaction. Decided by the shared provenance predicate in
    :mod:`t3.utils.uncertainty`.

Everything here is pure: it takes an already-loaded Arkane sensitivity dictionary and an already
-parsed network, and returns a decision object. It performs no I/O, holds no T3 state, and knows
nothing about paths, Arkane invocation, or species queueing -- so both T3's in-run path and the
public standalone API can call it and get the *same* decision.

Two notes on the numerics, both load-bearing:

* Coefficients are dln(k)/dE0 in mol/J, produced by perturbing E0 by Arkane's default
  2 kcal/mol. The relative gate alone cannot separate a real derivative from a denormal
  ~1e-18 structural zero, so an absolute floor is applied too. That floor is expressed as a
  **dimensionless rate response** (``min_delta_ln_k``, the smallest change in ln(k) worth
  caring about) and converted to a coefficient using the perturbation size, rather than being
  hard-coded as a raw coefficient -- so it keeps its meaning if the perturbation ever changes.
* The relative gate's denominator is taken over **transition-state entries only**. Well rows
  and TS rows are different derivative coordinates (a well energy versus a barrier), so mixing
  them would compare unlike quantities and would make the selected TS set change whenever an
  unrelated thermo row changed.
"""

import copy
import math
from dataclasses import dataclass, field

from t3.pdep.parser import PDepNetwork
from t3.utils.uncertainty import is_this_kinetics_comment_uncertain

# The Arkane sensitivity YAML marks transition-state rows with this prefix (``'(TS) ' + label``),
# and carries a non-reaction ``structures`` mapping alongside the network reaction keys.
TS_ENTRY_PREFIX = '(TS) '
STRUCTURES_KEY = 'structures'

# Arkane's default PDep sensitivity perturbation is 2 kcal/mol, applied to E0.
E0_PERTURBATION_J_PER_MOL = 8368.0

# The smallest ln(k) response worth calling "sensitive", used to derive the absolute floor.
DEFAULT_MIN_DELTA_LN_K = 1e-3

# SELECTION_SCHEMA_VERSION and SELECTION_ALGORITHM_VERSION used to be a single SELECTOR_VERSION
# doing three unrelated jobs at once -- SA-cache usability, on-disk envelope schema, and selection
# provenance -- so a change meant for one job forced a bump that wrongly invalidated the other two.
# Each constant below now has exactly one job; the fourth job (SA-cache usability) lives as
# ``t3.pdep.cache.SA_CACHE_CONTRACT_VERSION``, in cache.py, since that is the only module that
# reads or writes it.

# The SHAPE of one PDepNetworkSelection.as_dict() record: its set of keys and their types. Bump
# this when a field is added, removed, or renamed, or a field's rendered type changes, so a
# consumer reading an old saved-selection YAML from disk can tell it needs to migrate. Version 1
# means "the shape as of first ship"; pre-ship development churn on this shape is deliberately NOT
# versioned, because t3.pdep has never shipped a release -- no saved record has ever left this repo
# under a version number that would need to keep meaning the same thing. Do NOT bump this for a
# change to the decision LOGIC (see SELECTION_ALGORITHM_VERSION below) or to SA-cache usability
# (see t3.pdep.cache.SA_CACHE_CONTRACT_VERSION) -- neither of those changes the shape of a record.
SELECTION_SCHEMA_VERSION = 1

# The SEMANTICS of the decision: which gates are applied and how (relative/absolute thresholds,
# TS-only denominator, direction resolution, the uncertainty provenance predicate, ...). Bump this
# when the selection LOGIC changes in a way that could flip a past decision, so a saved selection
# can be told apart from one a newer selector would make given the same inputs. Do NOT bump this
# for a change to the on-disk SHAPE of a record (see SELECTION_SCHEMA_VERSION above) or to SA-cache
# usability (see t3.pdep.cache.SA_CACHE_CONTRACT_VERSION) -- neither of those changes what the
# decision means.
SELECTION_ALGORITHM_VERSION = 1

CACHE_STATUS_GENERATED = 'generated'
CACHE_STATUS_CACHED_VALID = 'cached_valid'
CACHE_STATUS_CACHED_REJECTED = 'cached_rejected'
# A caller passed validate_cache=False: the file is merely trusted, not actually validated against
# its T3 sidecar, and is not freshly generated either -- distinct provenance from both of the above.
CACHE_STATUS_UNVALIDATED = 'unvalidated'

# Evaluation status values for PDepNetworkSelection.evaluation_status (FIX C): 'evaluated' means the
# decision was actually computed from real SA data; 'not_evaluated' means it could not be (unreadable
# / unparseable network, missing SA data, or a rejected cache), so qualified/selected_ts carry no signal.
EVALUATION_STATUS_EVALUATED = 'evaluated'
EVALUATION_STATUS_NOT_EVALUATED = 'not_evaluated'


@dataclass(frozen=True)
class SensitiveTransitionState:
    """
    Evidence for one transition state that passed the sensitivity gates at one condition.

    Note that a path reaction's ``label`` is NOT unique within an RMG network file (several
    ``reaction(...)`` blocks routinely share e.g. ``'reaction1'``), so ``path_reaction_str`` is
    carried alongside it to identify the reaction unambiguously.

    Args:
        ts_label (str): The transition state label, e.g. ``'TS3'`` (without the ``'(TS) '`` prefix).
        coefficient (float): The signed dln(k)/dE0 sensitivity coefficient, in mol/J.
        condition (tuple): The (T, 'K', P, 'bar') condition at which it passed.
        path_reaction_label (str, optional): The label of the joined path reaction, if the join succeeded.
        path_reaction_str (str, optional): The joined path reaction as ``'A + B <=> C'``, if the join succeeded.
        kinetics_comment (str): The joined path reaction's RMG kinetics comment ( ``''`` if unjoined).
        uncertain (bool, optional): The provenance verdict, or ``None`` if the join failed.
        delta_ln_k (float): The dimensionless rate response this coefficient corresponds to,
            ``abs(coefficient) * perturbation``, i.e. the actual ``ln(k)`` change the E0
            perturbation produced. This is what the absolute floor (``min_delta_ln_k``) is
            expressed in, so it is carried alongside the raw coefficient for inspection.
    """
    ts_label: str
    coefficient: float
    condition: tuple
    path_reaction_label: str | None
    path_reaction_str: str | None
    kinetics_comment: str
    uncertain: bool | None
    delta_ln_k: float

    def as_dict(self) -> dict:
        """
        Render this record as plain JSON/YAML-safe types.

        The condition tuple ``(T, 'K', P, 'bar')`` is rendered as a dict with ``T``/``T_unit``/
        ``P``/``P_unit`` keys (rather than as a list) so that downstream consumers do not need to
        remember positional ordering.

        Returns:
            dict: A plain dict, containing no dataclass instances or tuples. ``condition`` is
                rendered as ``{'T': ..., 'T_unit': ..., 'P': ..., 'P_unit': ...}`` when it has the
                expected 4-element ``(T, T_unit, P, P_unit)`` shape; otherwise (a malformed but
                selected condition) it falls back to a plain list so serialization cannot crash.
        """
        # Deep-copied, not aliased: a condition is not guaranteed to be a flat tuple of scalars.
        # `_sensitive_transition_state_from_dict` coerces whatever a persisted record carries, so a
        # nested container can reach this record -- and handing one out would let a caller rewrite
        # the evidence behind a decision that has already been reported. Same defect class as the
        # manifest/selection/thresholds leaks already closed elsewhere in this package.
        condition = copy.deepcopy(self.condition)
        if isinstance(condition, (tuple, list)) and len(condition) == 4:
            temperature, temperature_unit, pressure, pressure_unit = condition
            condition = {'T': temperature, 'T_unit': temperature_unit,
                        'P': pressure, 'P_unit': pressure_unit}
        elif isinstance(condition, tuple):
            condition = list(condition)
        return {
            'ts_label': self.ts_label,
            'coefficient': self.coefficient,
            'condition': condition,
            'path_reaction_label': self.path_reaction_label,
            'path_reaction_str': self.path_reaction_str,
            'kinetics_comment': self.kinetics_comment,
            'uncertain': self.uncertain,
            'delta_ln_k': self.delta_ln_k,
        }


@dataclass
class PDepNetworkSelection:
    """
    The decision for one PDep network, returned by both the in-run path and the public API.

    Args:
        network_id (str): The network file stem, e.g. ``'network4_2'``. This is the join key
            downstream consumers (e.g. Carmel) use to identify the network.
        network_source_hash (str, optional): The ``'sha256:<hex>'`` content hash of the network file
            this decision was computed from, as ``t3.pdep.parser.hash_bytes`` /
            ``t3.pdep.cache.hash_file`` spell it, or ``None`` when the network was not parsed from a
            file. ``network_id`` is a FILE STEM and therefore names a file, not a content: a decision
            about one revision of ``network4_2.py`` matches every later revision of it by
            ``network_id`` alone. This field is what binds the decision to the bytes it was actually
            made about, so a consumer acting on it much later (e.g.
            ``t3.pdep.api.explore_pdep_network``, possibly in another process) can refuse a network
            that has changed since.
        qualified (bool): Whether criterion (b) holds, i.e. at least one sensitive transition
            state belongs to a path reaction with uncertain kinetics. The caller is responsible
            for criterion (a).
        network_reaction (str, optional): The network reaction whose sensitivities were examined,
            as requested by the caller.
        direction_key (str, optional): The Arkane SA key actually used, which may differ from
            ``network_reaction`` by label legalization or by direction.
        direction_keys (list): The per-component ``direction_key``s, in input order, recorded by
            ``combine()``; empty on a single-reaction decision.
        direction_ambiguous (bool): Whether the requested direction could not be resolved exactly
            and a different one was used instead. The decision is still reported, but it answers a
            slightly different question than the one asked.
        method (str, optional): The master-equation method whose SA was used, e.g. ``'MSC'``.
        sa_path (str, optional): The path to the Arkane sensitivity YAML that was read.
        cache_status (str, optional): One of ``'generated'``, ``'cached_valid'``, ``'cached_rejected'``.
        thresholds (dict): The gates actually applied, for reproducibility.
        selected_ts (list): ``SensitiveTransitionState`` entries that passed the gates.
        uncertain_path_reactions (list): The subset of ``selected_ts`` whose provenance is uncertain,
            i.e. the evidence for ``qualified``.
        warnings (list): Human-readable warnings (stale cache, failed TS join, non-finite rows,
            ambiguous direction, shared TS label, structurally dead TS rows).
        network_reactions_examined (int): The number of network reaction keys examined to produce
            this decision. For a single-reaction decision (``select_from_sa_dict``), this is ``1``
            once a direction key was resolved and ``0`` if it was not (the requested reaction could
            not be located in the SA output at all); set to the count of combined decisions by
            ``combine()``.
        evaluation_status (str): One of ``'evaluated'`` (the decision was actually computed from
            real SA data) or ``'not_evaluated'`` (it could not be: unreadable/unparseable network,
            missing SA data, or a rejected cache) -- in the latter case, ``qualified`` and
            ``selected_ts`` carry no signal and must not be read as "does not qualify".
        selection_schema_version (int): The on-disk SHAPE this record was built under; see
            ``SELECTION_SCHEMA_VERSION``. This is a FIELD, not a ``thresholds`` dict key, for two
            reasons: (a) ``PDepExplorationResult.as_dict()`` nests a serialized selection
            (``self.selection.as_dict()``); a version living only in the enclosing envelope cannot
            describe a nested record, but a field survives nesting. (b) ``thresholds`` is built by
            hand at four call sites, and a dict key can silently be omitted at one of them (as
            happened here) -- a dataclass field with a default cannot be.
        selection_algorithm_version (int): The decision SEMANTICS this record was produced by; see
            ``SELECTION_ALGORITHM_VERSION``. A field for the same two reasons as
            ``selection_schema_version`` above.
        t_grid_clamp (dict, optional): ``TGridClampRecord.as_dict()`` provenance for whether the
            Arkane SA T grid this decision rests on was clamped down from the network's original
            grid (see ``t3.utils.network_thermo.TGridClampRecord``). THREE states, not two:
            ``None`` means unknown provenance -- an old sidecar written before this field existed,
            or SA data produced outside T3 entirely -- and must never be read as "not clamped";
            a dict with ``'clamped': False`` means clamping was considered and explicitly did not
            apply; a dict with ``'clamped': True`` means the SA evidence backing this decision was
            computed over a narrower T range than the network's own thermo would otherwise allow.
            Missing/unknown provenance must never cause ``evaluation_status`` to become
            ``'not_evaluated'`` or any other refusal -- this field is purely descriptive.
    """
    network_id: str
    network_source_hash: str | None = None
    qualified: bool = False
    network_reaction: str | None = None
    direction_key: str | None = None
    direction_keys: list = field(default_factory=list)
    direction_ambiguous: bool = False
    method: str | None = None
    sa_path: str | None = None
    cache_status: str | None = None
    thresholds: dict = field(default_factory=dict)
    selected_ts: list = field(default_factory=list)
    uncertain_path_reactions: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    network_reactions_examined: int = 0
    evaluation_status: str = EVALUATION_STATUS_EVALUATED
    selection_schema_version: int = SELECTION_SCHEMA_VERSION
    selection_algorithm_version: int = SELECTION_ALGORITHM_VERSION
    t_grid_clamp: dict | None = None

    def uncertain_ts_labels(self) -> list:
        """
        The distinct transition state labels providing the evidence for qualification.

        Returns:
            list: Sorted distinct TS labels of the uncertain path reactions.
        """
        return sorted({entry.ts_label for entry in self.uncertain_path_reactions})

    def reason(self) -> str:
        """
        A human-readable rationale, in the style T3 already uses to record provenance.

        This is the only place a decision is rendered for a human, so it is the place a decision that
        was never made must not read as a negative one. The space is THREE cases, not two:

        - ``qualified`` -- a yes, and a yes stands even on an incomplete evaluation: one sensitive,
          uncertain transition state qualifies the network however much other evidence was unreadable
          (the module docstring's asymmetry -- a partial yes is supported, a partial no is not).
          ``combine()`` constructs exactly this record when one component qualified and another could
          not be evaluated, so it is a reachable state, not a defensive one.
        - not ``qualified`` and ``evaluation_status == 'evaluated'`` -- a real negative, computed over
          a complete accounting of the evidence.
        - not ``qualified`` and ``evaluation_status == 'not_evaluated'`` -- NO VERDICT. Here
          ``qualified`` is only its ``False`` dataclass default and carries no signal at all (see the
          field's docstring), so rendering "does not qualify" would attribute to T3 a decision it
          never made -- and this string is what gets written into the provenance record a human later
          reads to understand why a network was skipped. The diagnosis of WHY the criterion could not
          be computed already lives in ``warnings``; this method surfaces it rather than restating it.

        Returns:
            str: The rationale for this decision.
        """
        if self.qualified:
            labels = ', '.join(self.uncertain_ts_labels())
            partial = ('' if self.evaluation_status == EVALUATION_STATUS_EVALUATED
                       else ' Part of the evidence could not be evaluated; that does not weaken this '
                            'positive verdict, which rests on the transition state(s) named above.')
            if not labels:
                # Unreachable from ``select_from_sa_dict``/``combine()`` (both derive ``qualified``
                # from a non-empty ``uncertain_path_reactions``), but reachable from a persisted
                # record loaded off disk, whose fields are validated for TYPE and not for this
                # cross-field invariant. Naming no transition state is better than rendering
                # "sensitive to , whose kinetics ...".
                return (f'PDep network {self.network_id} is recorded as qualifying for QM refinement, but the '
                        f'record names no uncertain transition state as the evidence for that.{partial}')
            return (f'PDep network {self.network_id} qualifies for QM refinement: its rate coefficient is '
                    f'sensitive to {labels}, whose kinetics are estimates rather than library or training '
                    f'values.{partial}')
        if self.evaluation_status != EVALUATION_STATUS_EVALUATED:
            verdict = (f'PDep network {self.network_id} was NOT evaluated for QM refinement: the sensitivity '
                       f'criterion could not be computed, so no verdict on qualification was reached here -- '
                       f'this network was neither accepted nor rejected.')
            if self.warnings:
                return f'{verdict} Why it could not be computed: {" ".join(self.warnings)}'
            return f'{verdict} No diagnosis was recorded.'
        return (f'PDep network {self.network_id} does not qualify for QM refinement: no transition state '
                f'the network is sensitive to belongs to a path reaction with uncertain kinetics.')

    def as_dict(self) -> dict:
        """
        Render this decision as plain JSON/YAML-safe types.

        Returns:
            dict: A plain dict, containing no dataclass instances or tuples. ``selected_ts`` and
                ``uncertain_path_reactions`` are rendered via ``SensitiveTransitionState.as_dict()``.
                Every container field is deep-copied so that a caller mutating the returned dict
                (including anything nested inside it) can never reach back into this object's own
                live containers -- ``as_dict()`` is a reported, read-only snapshot.
        """
        return {
            'network_id': self.network_id,
            'network_source_hash': self.network_source_hash,
            'qualified': self.qualified,
            'network_reaction': self.network_reaction,
            'direction_key': self.direction_key,
            'direction_keys': copy.deepcopy(self.direction_keys),
            'direction_ambiguous': self.direction_ambiguous,
            'method': self.method,
            'sa_path': self.sa_path,
            'cache_status': self.cache_status,
            'thresholds': copy.deepcopy(self.thresholds),
            'selected_ts': [entry.as_dict() for entry in self.selected_ts],
            'uncertain_path_reactions': [entry.as_dict() for entry in self.uncertain_path_reactions],
            'warnings': copy.deepcopy(self.warnings),
            'network_reactions_examined': self.network_reactions_examined,
            'evaluation_status': self.evaluation_status,
            'selection_schema_version': self.selection_schema_version,
            'selection_algorithm_version': self.selection_algorithm_version,
            't_grid_clamp': copy.deepcopy(self.t_grid_clamp),
        }

    @classmethod
    def combine(cls, decisions: list) -> 'PDepNetworkSelection':
        """
        Combine several per-reaction decisions for the same network into one aggregate decision.

        Answers "does this network deserve QM refinement AT ALL" from the per-reaction decisions
        answering it one network reaction at a time. ``selected_ts``, ``uncertain_path_reactions``,
        and ``warnings`` are unioned across all inputs, de-duplicated while preserving first-seen
        order (not sorted arbitrarily). ``qualified`` is ``True`` iff any input decision qualified.
        ``direction_ambiguous`` is ``True`` iff any input decision had it ``True``. The per-component
        ``direction_key``s are recorded, in input order, as ``direction_keys`` on the combined
        decision (the combined decision's own ``direction_key`` is left ``None``, since it no longer
        refers to a single reaction).

        ``network_id`` is taken from the first decision, and all decisions must agree on it: since
        ``network_id`` is the join key downstream consumers use to identify the network, silently
        combining decisions for different networks would fabricate confidence in a result that does
        not actually describe one network. ``method`` and ``cache_status`` are also expected to agree
        across inputs (they come from the same SA source), but a disagreement there is a weaker
        signal (e.g. a partially re-generated cache) -- so instead of raising, the combined field is
        set to ``None`` and a warning is recorded, rather than silently adopting the first value.
        ``sa_path`` and ``thresholds`` are taken from the first decision unconditionally.

        ``network_source_hash`` is treated as identity, like ``network_id``: two components carrying
        two DIFFERENT real hashes were computed from different revisions of the network, so combining
        them raises. A hash present on some components and missing on others is only "not recorded",
        not "different", so the combined field is set to ``None`` with a warning rather than adopting
        the one real hash -- which would assert a content binding part of the aggregate does not have.

        ``evaluation_status`` is ``'not_evaluated'`` if ANY component was not evaluated, and a warning
        records how many. This is a statement of coverage, not of usability: whether a partially
        evaluated aggregate may still be acted on is a policy question, answered by the consumer
        (``t3.pdep.api.explore_pdep_network`` accepts one that qualified and refuses one that did not).

        ``selection_schema_version`` and ``selection_algorithm_version`` are treated as identity,
        like ``network_id``, and for a stronger reason than ``network_source_hash``: decisions
        combined here are always produced by ONE run of ONE process, so every component necessarily
        went through the same import of ``t3.pdep.selector`` and therefore the same versions.
        Disagreement is not a weaker signal to average away with a warning -- it means something is
        deeply wrong (e.g. components were collected across a code upgrade, or fabricated by a
        caller), so it is refused outright, exactly as a ``network_id`` disagreement is.

        ``qualified`` is unioned over the EVALUATED components only. This class is a mutable dataclass
        with no invariant enforcement, so a component can carry ``qualified=True`` alongside
        ``evaluation_status='not_evaluated'``; counting its vote would let a flag its own status
        declares meaningless decide the aggregate.

        ``network_reaction`` is set to ``None`` on the result, since it now represents the whole
        network rather than one reaction, and ``network_reactions_examined`` is set to the number
        of decisions combined. A single-component call is equivalent to that component, modulo the
        fields documented above (``network_reaction`` becomes ``None`` and ``direction_keys`` is
        populated).

        Args:
            decisions (list): The per-reaction ``PDepNetworkSelection`` decisions to combine.

        Raises:
            ValueError: If ``decisions`` is empty, if the decisions disagree on ``network_id``, if
                they carry two different non-``None`` ``network_source_hash`` values, or if they
                disagree on ``selection_schema_version`` or ``selection_algorithm_version``.

        Returns:
            PDepNetworkSelection: The combined decision.
        """
        if not decisions:
            raise ValueError('Cannot combine an empty list of decisions.')
        first = decisions[0]
        network_ids = {decision.network_id for decision in decisions}
        if len(network_ids) > 1:
            raise ValueError(f'Cannot combine decisions for different networks: {sorted(network_ids)}.')

        # ``selection_schema_version``/``selection_algorithm_version`` disagreement is refused
        # outright, like ``network_id`` -- see the docstring section above. Unlike
        # ``network_source_hash``, there is no "not recorded" case to tolerate: both fields always
        # carry a concrete int (their dataclass defaults), so any disagreement is a real one.
        schema_versions = {decision.selection_schema_version for decision in decisions}
        if len(schema_versions) > 1:
            raise ValueError(f'Cannot combine decisions with different selection_schema_version values '
                             f'for network {first.network_id}: {sorted(schema_versions)}.')
        algorithm_versions = {decision.selection_algorithm_version for decision in decisions}
        if len(algorithm_versions) > 1:
            raise ValueError(f'Cannot combine decisions with different selection_algorithm_version values '
                             f'for network {first.network_id}: {sorted(algorithm_versions)}.')

        # ``network_source_hash`` is an IDENTITY field, not a provenance nicety, so two distinct real
        # hashes are refused outright for the same reason two distinct network_ids are: the
        # components describe different content and any aggregate over them describes nothing that
        # exists. A hash missing on some components is weaker evidence -- it means "not recorded",
        # not "different" -- so that case records None and warns, exactly as a method disagreement
        # does. Adopting the one real hash there would be the fail-open: it would assert that the
        # aggregate is bound to those bytes when part of it demonstrably is not.
        source_hashes = {decision.network_source_hash for decision in decisions}
        known_source_hashes = {value for value in source_hashes if value is not None}
        if len(known_source_hashes) > 1:
            raise ValueError(f'Cannot combine decisions computed from different revisions of network '
                             f'{first.network_id}: {sorted(known_source_hashes)}.')

        methods = {decision.method for decision in decisions}
        cache_statuses = {decision.cache_status for decision in decisions}
        warnings = list()
        network_source_hash = first.network_source_hash
        if len(source_hashes) > 1:
            network_source_hash = None
            warnings.append(f'Only some combined decisions recorded a network_source_hash; recording None '
                            f'rather than binding the aggregate to bytes part of it was not computed from.')
        method = first.method
        if len(methods) > 1:
            method = None
            warnings.append(f'Combined decisions disagree on method ({sorted(str(m) for m in methods)}); '
                            f'recording None rather than silently adopting one.')
        cache_status = first.cache_status
        if len(cache_statuses) > 1:
            cache_status = None
            warnings.append(f'Combined decisions disagree on cache_status '
                            f'({sorted(str(s) for s in cache_statuses)}); recording None rather than '
                            f'silently adopting one.')

        # Only EVALUATED components get a vote on qualification. ``PDepNetworkSelection`` is a mutable
        # dataclass with no invariant enforcement, so `qualified=True, evaluation_status='not_evaluated'`
        # is a constructible state; ``any(decision.qualified ...)`` over every component would let a
        # component that was never evaluated carry the aggregate to `qualified=True` on a flag that,
        # by that component's own status, means nothing.
        not_evaluated = [decision for decision in decisions
                         if decision.evaluation_status != EVALUATION_STATUS_EVALUATED]
        qualified = any(decision.qualified for decision in decisions
                        if decision.evaluation_status == EVALUATION_STATUS_EVALUATED)

        combined = cls(
            network_id=first.network_id,
            network_source_hash=network_source_hash,
            qualified=qualified,
            network_reaction=None,
            direction_ambiguous=any(decision.direction_ambiguous for decision in decisions),
            method=method,
            sa_path=first.sa_path,
            cache_status=cache_status,
            thresholds=dict(first.thresholds),
            network_reactions_examined=len(decisions),
            selection_schema_version=first.selection_schema_version,
            selection_algorithm_version=first.selection_algorithm_version,
            t_grid_clamp=copy.deepcopy(first.t_grid_clamp),
        )
        # ``evaluation_status`` is NOT allowed to fall back to the dataclass default here: a fresh
        # PDepNetworkSelection is 'evaluated', so combining components that were never evaluated used
        # to produce a confidently 'evaluated', qualified=False aggregate -- precisely the "reports a
        # decision nobody made" failure ``evaluation_status`` exists to prevent, reintroduced one
        # level up.
        #
        # The status states a FACT about coverage ("was every component evaluated"), unconditionally.
        # It deliberately does NOT also encode the policy question "is this aggregate usable as a
        # gate" -- an earlier version exempted aggregates that qualified, on the reasoning that the
        # positive claim is backed by whichever component qualified. That reasoning is sound about
        # USABILITY and false as a statement of coverage: it made a partially evaluated aggregate
        # report that it was fully evaluated. The usability judgement now lives where it belongs, in
        # ``t3.pdep.api.explore_pdep_network``, which accepts a not_evaluated aggregate that
        # qualified and refuses one that did not.
        if not_evaluated:
            warnings.append(f'{len(not_evaluated)} of {len(decisions)} combined decisions were not '
                            f'evaluated; the aggregate covers only the evaluated ones.')
            combined.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
        combined.direction_keys = [decision.direction_key for decision in decisions]
        combined.selected_ts = _union_preserving_order(decision.selected_ts for decision in decisions)
        combined.uncertain_path_reactions = _union_preserving_order(
            decision.uncertain_path_reactions for decision in decisions)
        combined.warnings = _union_preserving_order([*(decision.warnings for decision in decisions), warnings])
        return combined


def _union_preserving_order(iterables) -> list:
    """
    Union several iterables of hashable items, de-duplicated, in first-seen order.

    Args:
        iterables: An iterable of iterables (e.g. a list of lists) to union.

    Returns:
        list: The de-duplicated union, in first-seen order.
    """
    seen = set()
    result = list()
    for iterable in iterables:
        for item in iterable:
            if item not in seen:
                seen.add(item)
                result.append(item)
    return result


def coefficient_floor(min_delta_ln_k: float = DEFAULT_MIN_DELTA_LN_K,
                      perturbation: float = E0_PERTURBATION_J_PER_MOL,
                      ) -> float:
    """
    Convert a minimum meaningful ln(k) response into an absolute sensitivity coefficient floor.

    Args:
        min_delta_ln_k (float, optional): The smallest ln(k) change worth calling sensitive.
        perturbation (float, optional): The E0 perturbation Arkane applied, in J/mol.

    Raises:
        ValueError: If ``min_delta_ln_k`` or ``perturbation`` is not a finite, positive number, or
            if both are individually valid but the DERIVED floor ``min_delta_ln_k / perturbation``
            is not a finite, positive number. Validating only the two inputs is not enough: finite,
            positive inputs can still divide to an overflowed inf (the absolute gate would then
            reject nothing) or an underflowed 0.0 (the absolute gate would then reject only exact
            zeros, silently disabling the denormal-structural-zero protection it exists for) --
            both fail-opens, not deliberate choices, so the derived value is checked too rather
            than trusting that valid inputs imply a valid result.

    Returns:
        float: The corresponding floor on abs(dln(k)/dE0), in mol/J.
    """
    if not _is_finite(min_delta_ln_k) or min_delta_ln_k <= 0:
        raise ValueError(f'min_delta_ln_k must be a finite, positive number, got {min_delta_ln_k}.')
    if not _is_finite(perturbation) or perturbation <= 0:
        raise ValueError(f'The perturbation must be a finite, positive number, got {perturbation}.')
    floor = min_delta_ln_k / perturbation
    if not _is_finite(floor) or floor <= 0:
        raise ValueError(f'The derived coefficient floor (min_delta_ln_k / perturbation) must be a '
                         f'finite, positive number, got {floor} from min_delta_ln_k={min_delta_ln_k} '
                         f'and perturbation={perturbation}.')
    return floor


def _validate_relative_threshold(relative_threshold: float) -> None:
    """
    Validate that a relative-threshold gate is a usable, non-negative fraction.

    A NaN relative_threshold makes every comparison against it False, so ``cutoff = max(nan * x,
    floor)`` would end up as whatever ``floor`` is -- silently discarding the relative gate rather
    than raising, which is a fail-open, not a deliberate "no relative gate" choice.

    Args:
        relative_threshold (float): The relative gate, as a fraction of the largest abs
            coefficient at the same condition.

    Raises:
        ValueError: If ``relative_threshold`` is not finite or is negative.
    """
    if not _is_finite(relative_threshold) or relative_threshold < 0:
        raise ValueError(f'relative_threshold must be a finite, non-negative number, '
                         f'got {relative_threshold}.')


def _bounded_cutoff(relative_threshold: float, max_abs: float, floor: float, direction_key: str) -> float:
    """
    Compute ``max(relative_threshold * max_abs, floor)``, guarding the DERIVED cutoff itself.

    ``relative_threshold`` and ``max_abs`` are each finite (validated at the call site / read from
    finite sensitivity data), but their product can still overflow to infinity. An infinite cutoff
    would silently fail every row closed -- a wrong, unreported negative verdict -- rather than
    surfacing the degenerate threshold that produced it, so it is refused instead.

    Args:
        relative_threshold (float): The relative gate, already validated finite and non-negative.
        max_abs (float): The largest absolute coefficient at this condition.
        floor (float): The absolute floor, already validated finite and positive.
        direction_key (str): The direction key this cutoff is being computed for, for the message.

    Raises:
        ValueError: If the derived cutoff is not finite.

    Returns:
        float: The cutoff, guaranteed finite.
    """
    cutoff = max(relative_threshold * max_abs, floor)
    if not _is_finite(cutoff):
        raise ValueError(f'The derived cutoff (relative_threshold * max_abs) for {direction_key} is not '
                         f'finite (got {cutoff} from relative_threshold={relative_threshold} and '
                         f'max_abs={max_abs}); refusing to silently reject every row against it.')
    return cutoff


def validate_selection_thresholds(relative_threshold: float,
                                  min_delta_ln_k: float = DEFAULT_MIN_DELTA_LN_K,
                                  perturbation: float = E0_PERTURBATION_J_PER_MOL,
                                  ) -> None:
    """
    Validate all three selection thresholds together, at a single entry point.

    Every caller-facing function that turns these three numbers into a decision
    (:func:`select_from_sa_dict`, :func:`select_sensitive_wells`,
    :func:`t3.pdep.api.select_pdep_network`, :func:`t3.pdep.api.rank_pdep_networks`) calls this
    before doing any work, so a bad threshold raises immediately at the caller rather than
    surfacing later as a silently wrong (over-permissive or over-refusing) selection.

    Args:
        relative_threshold (float): The relative gate; must be finite and >= 0.
        min_delta_ln_k (float, optional): The absolute gate's numerator; must be finite and > 0.
        perturbation (float, optional): The E0 perturbation, in J/mol; must be finite and > 0.

    Raises:
        ValueError: If any of the three is non-finite or out of its allowed range.
    """
    _validate_relative_threshold(relative_threshold)
    coefficient_floor(min_delta_ln_k=min_delta_ln_k, perturbation=perturbation)


def select_from_sa_dict(sa_dict: dict,
                        network: PDepNetwork,
                        network_reaction: str,
                        relative_threshold: float,
                        min_delta_ln_k: float = DEFAULT_MIN_DELTA_LN_K,
                        perturbation: float = E0_PERTURBATION_J_PER_MOL,
                        method: str | None = None,
                        sa_path: str | None = None,
                        cache_status: str | None = None,
                        t_grid_clamp: dict | None = None,
                        ) -> PDepNetworkSelection:
    """
    Decide whether a PDep network qualifies for QM refinement (criterion (b)).

    Args:
        sa_dict (dict): A loaded Arkane PDep sensitivity dictionary. Keys are network reaction
            strings (plus a ``'structures'`` entry); values map a (T, 'K', P, 'bar') condition to
            a mapping of entry label to sensitivity coefficient.
        network (PDepNetwork): The parsed network file, supplying the TS -> path reaction join
            and each path reaction's RMG kinetics comment.
        network_reaction (str): The network reaction that criterion (a) flagged as
            observable-sensitive, as ``'A + B <=> C'`` in network species labels.
        relative_threshold (float): The relative gate, as a fraction of the largest abs TS
            coefficient at the same condition (T3's ``pdep_SA_threshold``).
        min_delta_ln_k (float, optional): The absolute gate, as a minimum ln(k) response.
        perturbation (float, optional): The E0 perturbation Arkane applied, in J/mol.
        method (str, optional): The master-equation method used, recorded on the decision.
        sa_path (str, optional): The path the SA dictionary was read from, recorded on the decision.
        cache_status (str, optional): How the SA data was obtained, recorded on the decision.
        t_grid_clamp (dict, optional): ``TGridClampRecord.as_dict()`` provenance for whether the
            SA evidence's T grid was clamped, recorded verbatim on the decision. ``None`` means
            unknown provenance, not "not clamped" -- see ``PDepNetworkSelection``'s docstring.

    Raises:
        ValueError: If ``relative_threshold``, ``min_delta_ln_k``, or ``perturbation`` is
            non-finite or out of its allowed range (see :func:`validate_selection_thresholds`).

    Returns:
        PDepNetworkSelection: The decision, including the evidence behind it.
    """
    validate_selection_thresholds(relative_threshold=relative_threshold,
                                  min_delta_ln_k=min_delta_ln_k, perturbation=perturbation)
    if not isinstance(sa_dict, dict):
        selection = PDepNetworkSelection(
            network_id=network.network_id,
            network_source_hash=network.source_hash,
            network_reaction=network_reaction,
            method=method,
            sa_path=sa_path,
            cache_status=cache_status,
            t_grid_clamp=t_grid_clamp,
            thresholds={'relative_threshold': relative_threshold,
                        'min_delta_ln_k': min_delta_ln_k,
                        'perturbation': perturbation,
                        'coefficient_floor': coefficient_floor(min_delta_ln_k=min_delta_ln_k,
                                                               perturbation=perturbation),
                        },
        )
        selection.warnings.append(
            f'The sensitivity data for network {network.network_id} is malformed (expected a dict, got '
            f'{type(sa_dict).__name__}); cannot evaluate criterion (b).')
        selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
        return selection

    floor = coefficient_floor(min_delta_ln_k=min_delta_ln_k, perturbation=perturbation)
    selection = PDepNetworkSelection(
        network_id=network.network_id,
        network_source_hash=network.source_hash,
        network_reaction=network_reaction,
        method=method,
        sa_path=sa_path,
        cache_status=cache_status,
        t_grid_clamp=t_grid_clamp,
        thresholds={'relative_threshold': relative_threshold,
                    'min_delta_ln_k': min_delta_ln_k,
                    'perturbation': perturbation,
                    'coefficient_floor': floor,
                    },
    )

    direction_key, direction_ambiguous, direction_warnings = resolve_direction_key(
        sa_dict=sa_dict, network_reaction=network_reaction)
    selection.direction_key = direction_key
    selection.direction_ambiguous = direction_ambiguous
    selection.warnings.extend(direction_warnings)
    if direction_key is None:
        selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
        return selection

    # A direction key was resolved, so exactly one network reaction key was actually examined.
    selection.network_reactions_examined = 1

    resolved_entry = sa_dict.get(direction_key)
    if not isinstance(resolved_entry, dict):
        selection.warnings.append(
            f'The sensitivity data for SA key {direction_key} of network {selection.network_id} is malformed '
            f'(expected a dict of conditions, got {type(resolved_entry).__name__}); cannot evaluate criterion (b).')
        selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
        return selection

    by_ts = network.path_reactions_by_ts()
    non_finite, malformed, max_abs_ts_overall = 0, 0, 0.0
    ts_rows_seen, ts_rows_usable = 0, 0
    well_rows_finite, well_max_abs_overall = 0, 0.0
    for condition, entries in resolved_entry.items():
        if not isinstance(entries, dict):
            malformed += 1
            continue
        ts_entries = dict()
        for entry, coefficient in entries.items():
            if not isinstance(entry, str) or not entry.startswith(TS_ENTRY_PREFIX):
                # Not a TS row. Tracked separately (not counted in `malformed`/`non_finite`, which
                # are about the TS channel only) so that case (c) below can report whether the
                # non-TS (well) rows for this direction key carried any signal.
                if isinstance(entry, str) and _is_finite(coefficient):
                    well_rows_finite += 1
                    well_max_abs_overall = max(well_max_abs_overall, abs(float(coefficient)))
                continue
            ts_rows_seen += 1
            if not _is_finite(coefficient):
                non_finite += 1
                continue
            ts_rows_usable += 1
            ts_entries[entry[len(TS_ENTRY_PREFIX):]] = float(coefficient)
        if not ts_entries:
            continue
        max_abs_ts = max(abs(coefficient) for coefficient in ts_entries.values())
        max_abs_ts_overall = max(max_abs_ts_overall, max_abs_ts)
        cutoff = _bounded_cutoff(relative_threshold=relative_threshold, max_abs=max_abs_ts, floor=floor,
                                 direction_key=direction_key)
        for ts_label, coefficient in sorted(ts_entries.items()):
            if abs(coefficient) < cutoff:
                continue
            selection.selected_ts.extend(_evidence_for_ts(ts_label=ts_label,
                                                          coefficient=coefficient,
                                                          condition=condition,
                                                          by_ts=by_ts,
                                                          selection=selection,
                                                          perturbation=perturbation))

    if malformed:
        selection.warnings.append(f'Ignored {malformed} malformed condition entr{"y" if malformed == 1 else "ies"} '
                                  f'(expected a dict of entry label to coefficient) for {direction_key}.')

    if non_finite:
        selection.warnings.append(f'Ignored {non_finite} non-finite transition state sensitivity '
                                  f'coefficient(s) for {direction_key}.')
    # The below-floor space is four mutually exclusive, exhaustive situations over
    # (ts_rows_seen, malformed, ts_rows_usable, max_abs_ts_overall), not one. `evaluation_status`
    # states whether criterion (b) was actually evaluated; a negative verdict from data that could
    # not answer the question is exactly the fail-open this branch closes (see the module
    # docstring's asymmetry: a partial yes is supported, a partial no is not). Cases (c) and (d)
    # were once separate (an exact-zero TS response was treated as a dead instrument, while a
    # below-floor-but-nonzero response was treated as a considered, fully-evaluated negative); the
    # first PDep-QM trial run showed that distinction does not survive contact with real data (see
    # the merged case below), so they are now one case.
    if ts_rows_seen == 0 and not malformed:
        # (a) No condition block was malformed and none contained a TS-prefixed row: there was
        # genuinely nothing to measure, so this is not a below-floor response -- there is no
        # response to speak of.
        selection.warnings.append(
            f'No transition state sensitivity rows were found for {direction_key}; criterion (b) has not '
            f'been evaluated (there is nothing to measure here, not a below-floor response).')
        selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
    elif ts_rows_seen == 0:
        # (a2) malformed > 0 and no TS-prefixed row was seen in any READABLE block: unlike case (a),
        # it is not known whether TS rows existed for this direction_key -- an unreadable block may
        # have hidden the only one. This is a distinct diagnosis from case (a) ("there were
        # genuinely none") and from case (b) below ("rows were seen, none usable"): here, no rows
        # were even seen, but that is because reading failed, not because there was nothing there.
        selection.warnings.append(
            f'{malformed} condition entr{"y was" if malformed == 1 else "ies were"} malformed for '
            f'{direction_key}, and no transition state row was seen in any readable condition block; '
            f'whether transition state sensitivity rows exist here cannot be determined, so criterion '
            f'(b) has not been evaluated.')
        selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
    elif ts_rows_usable == 0:
        # (b) TS-prefixed rows were seen, but every one of them was non-finite (or otherwise
        # unusable), so none could answer the question either.
        selection.warnings.append(
            f'{ts_rows_seen} transition state sensitivity row(s) were found for {direction_key}, but none '
            f'were usable (all were non-finite); criterion (b) has not been evaluated.')
        selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
    elif max_abs_ts_overall < floor:
        # (c) Usable TS rows exist and none reaches the absolute floor -- whether every one is
        # EXACTLY zero (the structural-zero signature of an Arkane that perturbs only a
        # synthesized TS E0 without reaching the ILT rate expression) or merely a real, finite,
        # measured response smaller than the smallest response worth acting on (a denormal like
        # 1e-18), both signatures mean the same thing: a live transition state sensitivity would
        # sit orders of magnitude above this floor. The first PDep-QM trial run measured this
        # directly: TS coefficients of ~1e-16-1e-14 mol/J against a floor of ~1.2e-7 mol/J (four to
        # six orders of magnitude short), while well coefficients for the very same network sat at
        # ~1e-7-1e-6 mol/J, at or above the floor, and Arkane's own log showed it had perturbed the
        # TS by the same magnitude as the wells. The TS channel was dead, not merely small.
        # This case used to be split in two: an exact zero was a dead instrument (not_evaluated),
        # while a below-floor nonzero was treated as a considered, fully-evaluated negative
        # (evaluated, qualified=False). That distinction does not survive contact with real data --
        # the trial's below-floor-but-nonzero TS coefficients would have been reported by the old
        # case (d) as a confident "does not qualify", exactly the fail-open this case exists to
        # prevent (see the module docstring's asymmetry: a partial yes is supported, a partial no
        # is not). Both signatures are now one case: a below-floor reading is a dead instrument,
        # not a small answer, so criterion (b) cannot be evaluated either way.
        if well_rows_finite:
            well_note = (f' {well_rows_finite} finite non-transition-state (well) row(s) for the same '
                        f'direction key had a largest absolute coefficient of {well_max_abs_overall:.3e} '
                        f'mol/J, so the wells did respond to the perturbation even though the transition '
                        f'state did not.')
        else:
            well_note = ' No finite non-transition-state (well) rows were found for this direction key either.'
        exact_zero_note = ' (every usable coefficient was exactly zero)' if max_abs_ts_overall == 0.0 else ''
        selection.warnings.append(
            f'The largest transition state sensitivity coefficient for {direction_key} '
            f'({max_abs_ts_overall:.3e} mol/J) did not reach the absolute floor ({floor:.3e} mol/J)'
            f'{exact_zero_note}. A live transition state sensitivity would sit orders of magnitude above '
            f'the floor, so this reading means the transition state channel carries no signal -- not a '
            f'small-but-real answer -- and criterion (b) cannot be evaluated; this usually means the '
            f'sensitivity analysis was run with an Arkane that does not propagate a transition state '
            f'perturbation into an ILT-based rate, which RMG-Py PR #2990 '
            f'(https://github.com/ReactionMechanismGenerator/RMG-Py/pull/2990) fixes -- check that the '
            f'Arkane on PYTHONPATH carries it before treating this network as uninteresting.' + well_note)
        selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED

    selection.uncertain_path_reactions = [entry for entry in selection.selected_ts if entry.uncertain]
    selection.qualified = bool(selection.uncertain_path_reactions)
    if not selection.qualified:
        # A negative verdict is only supported if it was computed over a complete accounting of the
        # evidence. If this function is still reporting 'evaluated' at this point (i.e. none of cases
        # (a)/(a2)/(b)/(c) above already downgraded it) but rows were discarded as malformed or
        # non-finite along the way, the negative rests on an incomplete picture: one of the discarded
        # rows might have been the sensitive, uncertain one that would have flipped this verdict. This
        # is the same asymmetry as cases (a)-(c) above, applied uniformly to every negative verdict
        # rather than only to the ones already caught by a specific case.
        discarded = malformed + non_finite
        if discarded and selection.evaluation_status == EVALUATION_STATUS_EVALUATED:
            selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
            selection.warnings.append(
                f'The negative verdict for {direction_key} rests on incomplete data: {malformed} malformed '
                f'condition entr{"y" if malformed == 1 else "ies"} and {non_finite} non-finite transition '
                f'state coefficient(s) were discarded before reaching it; criterion (b) has not been '
                f'evaluated (the discarded data might have qualified this network).')
        # A partial "no" is also only supported if every selected TS was actually assessed. An entry
        # with ``uncertain is None`` (its provenance could not be joined to a path reaction, see
        # ``_evidence_for_ts``) is a non-answer, not a "certain" vote; letting it stand in for a real
        # negative would be the same fail-open closed above -- a confident False built from data that
        # never answered the question.
        unassessed_labels = sorted({entry.ts_label for entry in selection.selected_ts
                                    if entry.uncertain is None})
        if unassessed_labels:
            selection.evaluation_status = EVALUATION_STATUS_NOT_EVALUATED
            selection.warnings.append(
                f'The qualification verdict for {direction_key} rests on transition state(s) whose '
                f"provenance could not be assessed ({', '.join(unassessed_labels)}); criterion (b) "
                f'has not been evaluated.')
    return selection


def resolve_direction_key(sa_dict: dict,
                          network_reaction: str,
                          ) -> tuple:
    """
    Find the Arkane SA key corresponding to a requested network reaction.

    The two directions of a reaction appear as separate keys in the SA output, so the requested
    direction is preferred. Matching falls back through label legalization (ARC rewrites
    ``'('``/``')'`` as ``'['``/``']'``) and then to the opposite direction: Arkane computes its
    sensitivity coefficients from the network's (direction-independent) rate parameters, so the
    forward- and reverse-direction coefficients agree to solver tolerance, but using the opposite
    entry is still flagged ambiguous since it was not the direction actually requested.

    Args:
        sa_dict (dict): The loaded Arkane sensitivity dictionary.
        network_reaction (str): The requested reaction, as ``'A + B <=> C'``.

    Returns:
        tuple: (key or ``None``, direction_ambiguous (bool), warnings (list)).
    """
    warnings_list = list()
    if not isinstance(sa_dict, dict):
        warnings_list.append(f'The sensitivity data is malformed (expected a dict, got '
                             f'{type(sa_dict).__name__}); cannot locate reaction {network_reaction}.')
        return None, False, warnings_list
    candidates = [key for key in sa_dict.keys() if key != STRUCTURES_KEY and isinstance(key, str)]
    if network_reaction in candidates:
        return network_reaction, False, warnings_list

    canonical = _canonical_reaction(network_reaction)
    matches = [key for key in candidates if _canonical_reaction(key) == canonical]
    if len(matches) == 1:
        warnings_list.append(f'Matched network reaction {network_reaction} to SA key {matches[0]} '
                             f'after label canonicalization.')
        return matches[0], False, warnings_list
    if len(matches) > 1:
        warnings_list.append(f'Network reaction {network_reaction} matches {len(matches)} SA keys '
                             f'({", ".join(matches)}); using {matches[0]}.')
        return matches[0], True, warnings_list

    reverse = _canonical_reaction(_reverse_reaction(network_reaction))
    reverse_matches = [key for key in candidates if _canonical_reaction(key) == reverse]
    if reverse_matches:
        warnings_list.append(f'The SA output has no entry for the requested direction {network_reaction}; '
                             f'using the opposite-direction entry ({reverse_matches[0]}) instead. Arkane\'s '
                             f'sensitivity coefficients are equal for the two directions to solver tolerance.')
        return reverse_matches[0], True, warnings_list

    warnings_list.append(f'Could not locate reaction {network_reaction} in the SA output.')
    return None, False, warnings_list


def select_sensitive_wells(entries_by_condition: dict,
                           relative_threshold: float,
                           min_delta_ln_k: float = DEFAULT_MIN_DELTA_LN_K,
                           perturbation: float = E0_PERTURBATION_J_PER_MOL,
                           ) -> dict:
    """
    Identify the wells and channels a network reaction's rate coefficient is sensitive to.

    Ranking is by absolute value: sensitivity coefficients are signed, and a strongly *negative*
    coefficient is just as much a sensitivity as a positive one. The same absolute floor used for
    transition states is applied, which keeps denormal structural zeros out of the result.

    Args:
        entries_by_condition (dict): The SA dictionary entry for one network reaction, mapping a
            condition to a mapping of entry label to coefficient.
        relative_threshold (float): The relative gate, as a fraction of the largest abs well
            coefficient at the same condition.
        min_delta_ln_k (float, optional): The absolute gate, as a minimum ln(k) response.
        perturbation (float, optional): The E0 perturbation Arkane applied, in J/mol.

    Raises:
        ValueError: If ``relative_threshold``, ``min_delta_ln_k``, or ``perturbation`` is
            non-finite or out of its allowed range (see :func:`validate_selection_thresholds`).

    Returns:
        dict: Keys are well labels, values are lists of the conditions at which they were sensitive.
    """
    validate_selection_thresholds(relative_threshold=relative_threshold,
                                  min_delta_ln_k=min_delta_ln_k, perturbation=perturbation)
    floor = coefficient_floor(min_delta_ln_k=min_delta_ln_k, perturbation=perturbation)
    sensitive_wells = dict()
    for condition, entries in entries_by_condition.items():
        wells = {entry: float(coefficient) for entry, coefficient in entries.items()
                 if isinstance(entry, str) and not entry.startswith(TS_ENTRY_PREFIX)
                 and _is_finite(coefficient)}
        if not wells:
            continue
        max_abs_well = max(abs(coefficient) for coefficient in wells.values())
        cutoff = _bounded_cutoff(relative_threshold=relative_threshold, max_abs=max_abs_well, floor=floor,
                                 direction_key=str(condition))
        for well, coefficient in wells.items():
            if abs(coefficient) >= cutoff:
                sensitive_wells.setdefault(well, list()).append(condition)
    return sensitive_wells


def _evidence_for_ts(ts_label: str,
                     coefficient: float,
                     condition: tuple,
                     by_ts: dict,
                     selection: PDepNetworkSelection,
                     perturbation: float = E0_PERTURBATION_J_PER_MOL,
                     ) -> list:
    """
    Build the evidence records for one sensitive transition state, joining it to its path reaction.

    Args:
        ts_label (str): The transition state label.
        coefficient (float): Its sensitivity coefficient.
        condition (tuple): The condition at which it passed the gates.
        by_ts (dict): The network's TS label -> path reactions map.
        selection (PDepNetworkSelection): The decision being built, appended to on join problems.
        perturbation (float, optional): The E0 perturbation Arkane applied, in J/mol, used to
            compute ``delta_ln_k`` on the resulting records.

    Returns:
        list: The ``SensitiveTransitionState`` records for this transition state.
    """
    delta_ln_k = abs(coefficient) * perturbation
    path_reactions = by_ts.get(ts_label, tuple())
    if not path_reactions:
        message = (f'Could not join sensitive transition state {ts_label} to a path reaction of network '
                   f'{selection.network_id}; its provenance could not be assessed.')
        if message not in selection.warnings:
            selection.warnings.append(message)
        return [SensitiveTransitionState(ts_label=ts_label,
                                         coefficient=coefficient,
                                         condition=condition,
                                         path_reaction_label=None,
                                         path_reaction_str=None,
                                         kinetics_comment='',
                                         uncertain=None,
                                         delta_ln_k=delta_ln_k,
                                         )]
    if len(path_reactions) > 1:
        message = (f'Transition state {ts_label} of network {selection.network_id} is shared by '
                   f'{len(path_reactions)} path reactions; all of them were assessed.')
        if message not in selection.warnings:
            selection.warnings.append(message)
    return [SensitiveTransitionState(
        ts_label=ts_label,
        coefficient=coefficient,
        condition=condition,
        path_reaction_label=path_reaction.label,
        path_reaction_str=_path_reaction_str(path_reaction),
        kinetics_comment=path_reaction.kinetics_comment,
        uncertain=is_this_kinetics_comment_uncertain(kinetics_comment=path_reaction.kinetics_comment),
        delta_ln_k=delta_ln_k,
    ) for path_reaction in path_reactions]


def _path_reaction_str(path_reaction) -> str:
    """
    Render a path reaction as ``'A + B <=> C'``.

    Args:
        path_reaction (PDepPathReaction): The path reaction.

    Returns:
        str: The rendered reaction string.
    """
    return f"{' + '.join(path_reaction.reactants)} <=> {' + '.join(path_reaction.products)}"


def _canonical_reaction(reaction_str: str) -> tuple:
    """
    Canonicalize a reaction string so equivalent labelings compare equal.

    ARC legalizes species labels by rewriting ``'('`` and ``')'`` as ``'['`` and ``']'``, and the
    order of species within a side carries no meaning, so both are normalized away. The two sides
    are NOT sorted against each other: direction is meaningful here.

    Args:
        reaction_str (str): The reaction string, as ``'A + B <=> C'``.

    Returns:
        tuple: A canonical (reactants, products) pair of sorted label tuples.
    """
    sides = reaction_str.split(' <=> ')
    canonical_sides = list()
    for side in sides:
        labels = [label.strip().replace('[', '(').replace(']', ')') for label in side.split(' + ')]
        canonical_sides.append(tuple(sorted(labels)))
    return tuple(canonical_sides)


def _reverse_reaction(reaction_str: str) -> str:
    """
    Swap the two sides of a reaction string.

    Args:
        reaction_str (str): The reaction string, as ``'A + B <=> C'``.

    Returns:
        str: The reversed reaction string. Returned unchanged if it has no ``' <=> '`` separator.
    """
    sides = reaction_str.split(' <=> ')
    if len(sides) != 2:
        return reaction_str
    return f'{sides[1]} <=> {sides[0]}'


def _is_finite(coefficient) -> bool:
    """
    Whether a sensitivity coefficient is a usable finite number.

    Args:
        coefficient: The value read from the SA YAML (Arkane writes ``.nan`` on repeated failure).

    Returns:
        bool: ``True`` if it is a finite real number.
    """
    return isinstance(coefficient, (int, float)) and not isinstance(coefficient, bool) \
        and math.isfinite(coefficient)
