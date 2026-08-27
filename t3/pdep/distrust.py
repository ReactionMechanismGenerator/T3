"""
The pre-QM "rank by distrust" selection criterion for the PES exploration loop (I-032).

Why this exists -- the criterion it replaces is structurally blind
-----------------------------------------------------------------
The loop's other selector ranks transition states by master-equation E0-sensitivity: perturb a
saddle's ``E0``, measure the ln(k) response, queue the movers. That criterion cannot see the very
channels a mechanism most needs. A transition state in a reduced network file carries no statmech
``modes`` (nobody has run the QM yet -- that IS the queue), so ``Reaction.can_tst()`` is False and
every k(E) comes from the inverse Laplace transform, where the saddle's ``E0`` enters only as a
discrete threshold gate (``rmgpy/pdep/reaction.pyx``). For a bimolecular association/dissociation
channel that gate never binds -- the rate is fixed by the bimolecular asymptote, the
detailed-balance clamp, and high-pressure-limit renormalization -- so perturbing the saddle's
``E0`` is an exact no-op and the sensitivity is a STRUCTURAL ZERO (``0.0`` or ~1e-18), not a
measurement of low leverage. The blindness is worst exactly for the low bimolecular entrance
channels that carry the reactive flux. You cannot use a quantity's sensitivity to decide whether to
compute that quantity when not having computed it is what makes the sensitivity meaningless.

What distrust ranks by instead
------------------------------
Every input here is present in the reduced network file BEFORE any QM runs, and none of it is a
sensitivity value -- so no structural zero can suppress a candidate, and the "invert the
sensitivity" trap (which inherits the same zero) is avoided by construction. A candidate is ranked
by how little we trust the barrier we currently hold:

* **Flat energy window (the gate).** A saddle far above the lowest one on the surface carries
  negligible flux at the temperatures of interest, so it is declined regardless of provenance. Only
  saddles within ``energy_window_kj`` of the lowest surface saddle stay eligible. This is what makes
  the criterion a screen rather than a firehose: it declines the high isomerization saddle on the
  r002 CHO2 surface while keeping the low ``[H] + O=C=O`` entrance channels.
* **Provenance (dominant rank term).** A family/rate-rule estimate is distrusted; a library or
  already-computed value is trusted and not re-queued. Read from the kinetics comment.
* **RMG's own variance (secondary rank term).** When RMG produced the rate by decision-tree node
  averaging it attaches ``RateUncertainty(var=...)`` -- its own voice saying "I guessed this one".
  Larger ``var`` -> more distrust. Optional: absent on a plain estimate, in which case it simply
  contributes nothing.
* **Barrier height (tiebreak).** Among equally-provenanced, equal-variance saddles, the lower
  barrier carries more flux and so is worth getting right first.

Nothing here is computed, estimated, or invented: the energies and comments are read verbatim from
the network file, and a missing datum leaves its term neutral rather than being guessed at.
"""

from dataclasses import dataclass, replace

from t3.pdep.parser import PDepNetworkE0, PDepPathReaction
from t3.pdep.pes_rounds import (CandidateSplit, QMCandidate, SkippedChannel,
                                SKIP_NO_EVIDENCE, SKIP_OUTSIDE_WINDOW, SKIP_TRUSTED_PROVENANCE,
                                SKIP_UNMEASURABLE)

# The flat energy window half-width, in kJ/mol: a saddle more than this above the lowest saddle on
# the surface is declined. ~30 kJ/mol is the figure discussed for the r002 CHO2 network -- wide
# enough to hold both low entrance channels, tight enough to decline the ~114 kJ/mol-higher
# isomerization saddle.
DEFAULT_ENERGY_WINDOW_KJ = 30.0

# Provenance classes of the current barrier, read from the kinetics comment.
PROVENANCE_LIBRARY = 'library'
PROVENANCE_ESTIMATE = 'estimate'

# Case-insensitive substrings that mark a barrier as a trusted library/QM value rather than a
# family/rate-rule estimate. RMG writes "Reaction library: '<name>'" for a value taken from a
# kinetics library (which is where a completed QM job's fitted rate is stored, e.g. 'kineticsjobs').
_LIBRARY_MARKERS = ('reaction library', 'library:')

# The size of the provenance rank gap. It must exceed any realistic sum of the variance and barrier
# terms so that an estimate always outranks a library value: RMG's ln-space ``var`` is O(10) and the
# barrier tiebreak is in [0, 1), so 1e6 is never crossed.
_PROVENANCE_GAP = 1.0e6


@dataclass(frozen=True)
class DistrustParams:
    """
    The knobs of the distrust criterion.

    Attributes:
        energy_window_kj (float): The flat energy window half-width (kJ/mol). A candidate saddle
            more than this above the lowest saddle on the surface is declined.
    """
    energy_window_kj: float = DEFAULT_ENERGY_WINDOW_KJ


@dataclass(frozen=True)
class DistrustScore:
    """
    The transparent components behind one candidate's distrust ranking.

    Every field is derived from the pre-QM network file, so the ranking is auditable term by term
    rather than collapsed into an opaque scalar.

    Attributes:
        ts_label (str): The network-local transition state label.
        barrier_kj (float | None): The saddle's height above its lower-energy adjacent configuration
            (kJ/mol), or ``None`` if a side's species energy was missing from the file.
        height_above_lowest_saddle_kj (float | None): The saddle's E0 minus the lowest saddle E0 on
            the surface (kJ/mol), the flat-window quantity, or ``None`` if this saddle's E0 was
            missing.
        in_window (bool): Whether the saddle is within ``energy_window_kj`` of the lowest surface
            saddle. A saddle whose E0 could not be read is kept (``True``) rather than silently
            dropped -- we cannot place it, so we do not rule it out.
        provenance (str): ``PROVENANCE_LIBRARY`` (trusted) or ``PROVENANCE_ESTIMATE`` (distrusted).
        kinetics_var (float | None): RMG's ``RateUncertainty(var=...)`` for this rate, or ``None``.
        score (float): The distrust score; higher means less trusted and so more worth computing.
    """
    ts_label: str
    barrier_kj: float | None
    height_above_lowest_saddle_kj: float | None
    in_window: bool
    provenance: str
    kinetics_var: float | None
    score: float


def classify_provenance(kinetics_comment: str) -> str:
    """
    Classify a barrier's provenance from its kinetics comment.

    A comment naming a kinetics library (``"Reaction library: '...'"``) -- which is where a completed
    QM job's fitted rate is stored -- is a value we already trust and would not spend QM on again. A
    comment that does not is treated as a distrusted estimate; this deliberately fails toward
    computing, since a barrier whose provenance we cannot read is not one we should trust by default.

    Args:
        kinetics_comment (str): The reaction's kinetics ``comment`` text (possibly ``''``).

    Returns:
        str: ``PROVENANCE_LIBRARY`` or ``PROVENANCE_ESTIMATE``.
    """
    lowered = (kinetics_comment or '').lower()
    if any(marker in lowered for marker in _LIBRARY_MARKERS):
        return PROVENANCE_LIBRARY
    return PROVENANCE_ESTIMATE


def _side_energy(labels: tuple, species_e0: dict) -> float | None:
    """
    Sum the E0 of one reaction side, or ``None`` if any species' energy is missing.

    Args:
        labels (tuple): The species labels on one side of the reaction.
        species_e0 (dict): Species label -> E0 (kJ/mol), from ``PDepNetworkE0.species``.

    Returns:
        float | None: The summed side energy (kJ/mol), or ``None`` if a label had no E0.
    """
    if not labels:
        return None
    total = 0.0
    for label in labels:
        if label not in species_e0:
            return None
        total += species_e0[label]
    return total


def compute_distrust(path_reaction: PDepPathReaction,
                     ts_label: str,
                     network_e0: PDepNetworkE0,
                     lowest_saddle_e0_kj: float | None,
                     params: DistrustParams) -> DistrustScore:
    """
    Score one candidate by distrust from the pre-QM network file.

    Args:
        path_reaction (PDepPathReaction): The parsed path reaction (its species and kinetics comment
            and any ``RateUncertainty(var=...)``).
        ts_label (str): The candidate's network-local transition state label.
        network_e0 (PDepNetworkE0): The declared species and transition-state E0 values (kJ/mol).
        lowest_saddle_e0_kj (float | None): The lowest transition-state E0 on the surface (kJ/mol),
            or ``None`` if the file declared no transition-state E0 at all.
        params (DistrustParams): The window and any other knobs.

    Returns:
        DistrustScore: The transparent components and the scalar score for this candidate.
    """
    ts_e0 = network_e0.transition_states.get(ts_label)
    reactant_side = _side_energy(path_reaction.reactants, network_e0.species)
    product_side = _side_energy(path_reaction.products, network_e0.species)

    # Barrier over the lower-energy adjacent configuration: direction-robust and non-negative for a
    # real saddle. Either side may be unreadable; use whichever is present, and leave the barrier
    # None only when neither is.
    lower_side = None
    for side in (reactant_side, product_side):
        if side is not None and (lower_side is None or side < lower_side):
            lower_side = side
    barrier_kj = None
    if ts_e0 is not None and lower_side is not None:
        barrier_kj = ts_e0 - lower_side

    height_kj = None
    if ts_e0 is not None and lowest_saddle_e0_kj is not None:
        height_kj = ts_e0 - lowest_saddle_e0_kj
    # A saddle whose E0 could not be read is kept rather than dropped: we cannot place it on the
    # surface, so we do not rule it out of the window.
    in_window = height_kj is None or height_kj <= params.energy_window_kj

    provenance = classify_provenance(path_reaction.kinetics_comment)
    kinetics_var = path_reaction.kinetics_uncertainty_var

    provenance_term = _PROVENANCE_GAP if provenance == PROVENANCE_ESTIMATE else 0.0
    variance_term = kinetics_var if kinetics_var is not None else 0.0
    # Barrier tiebreak in [0, 1): monotone decreasing in the barrier, so a lower barrier scores
    # higher. A barrier we could not read is treated as maximally distrusted (1.0).
    flux_term = 1.0 if barrier_kj is None else 1.0 / (1.0 + max(barrier_kj, 0.0))
    score = provenance_term + variance_term + flux_term

    return DistrustScore(ts_label=ts_label, barrier_kj=barrier_kj,
                         height_above_lowest_saddle_kj=height_kj, in_window=in_window,
                         provenance=provenance, kinetics_var=kinetics_var, score=score)


def rank_candidates_by_distrust(candidates: tuple,
                                network_e0: PDepNetworkE0,
                                params: DistrustParams) -> tuple:
    """
    Score every candidate and split them into the eligible (in-window), ranked by distrust, and the
    declined (outside the window).

    This is the pure ranking core, with no evidence-stamping or ``CandidateSplit`` bookkeeping, so
    it can be replayed and unit-tested directly on a parsed network.

    Args:
        candidates (tuple): ``QMCandidate`` objects to rank.
        network_e0 (PDepNetworkE0): The declared E0 values (kJ/mol).
        params (DistrustParams): The window and any other knobs.

    Returns:
        tuple: ``(eligible, declined)`` where ``eligible`` is a list of ``(QMCandidate,
        DistrustScore)`` in descending distrust order (ties keeping input order, which is
        network-file order), and ``declined`` is a list of ``(QMCandidate, DistrustScore)`` for the
        out-of-window candidates in input order.
    """
    lowest_saddle_e0_kj = (min(network_e0.transition_states.values())
                           if network_e0.transition_states else None)
    scored = [(candidate,
               compute_distrust(candidate.path_reaction, candidate.ts_label, network_e0,
                                lowest_saddle_e0_kj, params))
              for candidate in candidates]
    eligible = [pair for pair in scored if pair[1].in_window]
    declined = [pair for pair in scored if not pair[1].in_window]
    # sorted() is stable, so equal scores keep network-file order -- the same determinism the
    # sensitivity path promises.
    eligible.sort(key=lambda pair: -pair[1].score)
    return eligible, declined


def _outside_window_reason(candidate: QMCandidate, score: DistrustScore,
                           params: DistrustParams) -> str:
    """Prose for a candidate declined because its saddle sits above the flat energy window."""
    height = ('unknown' if score.height_above_lowest_saddle_kj is None
              else f'{score.height_above_lowest_saddle_kj:.1f}')
    return (f"'{candidate.path_reaction.label}': transition state {candidate.ts_label} sits "
            f'{height} kJ/mol above the lowest saddle on the surface, outside the '
            f'{params.energy_window_kj:.1f} kJ/mol flat energy window; it carries too little '
            f'reactive flux to justify the QM spend, whatever the (dis)trust in its barrier.')


def select_by_distrust(split: CandidateSplit,
                       evidence_by_ts_label: dict,
                       network_e0: PDepNetworkE0,
                       params: DistrustParams) -> CandidateSplit:
    """
    Select and rank a network's QM candidates by distrust, in place of the sensitivity floor.

    This is the distrust counterpart to ``t3.pdep.pes_rounds.attach_sensitivity_evidence`` and is
    called from ``t3.pdep.pes_loop`` when ``qm.scope == 'distrust'``. It does three things:

    0. **Skips trusted (library / already-computed) candidates** before anything else. A barrier
       whose kinetics comment names a reaction library is a value distrust trusts and would not spend
       QM on again (the module contract); it is appended to ``skipped`` as ``SKIP_TRUSTED_PROVENANCE``
       and never ranked, so an all-library network queues nothing rather than redundantly recomputing
       what it already holds.
    1. **Stamps the measured sensitivity coefficient** (from ``evidence_by_ts_label``) onto every
       remaining candidate that has a finite row, exactly as the sensitivity path does. The coefficient is NOT
       the selection basis here -- distrust is -- but it is still carried so ``t3.pdep.capture`` (a
       fail-closed guard predating this scope) accepts the queued artifact, and so the round record
       keeps the sensitivity value visible (I-031: stop ranking on the structural zero, do not hide
       it). A candidate with no finite row is still skipped, since capture has no coefficient to
       record and inventing one is forbidden; it is classified ``SKIP_UNMEASURABLE`` when it was the
       structurally-unmeasurable kind, else ``SKIP_NO_EVIDENCE`` -- unchanged from the sensitivity
       path.
    2. **Declines out-of-window candidates**, appending each to ``skipped`` as
       ``SKIP_OUTSIDE_WINDOW`` with its height above the lowest saddle in the reason. This is the
       negative control: the criterion is a screen, not a firehose.
    3. **Ranks the survivors by descending distrust** and returns them in that order, each carrying
       its ``distrust_score``. ``t3.pdep.pes_loop._trim_candidates`` under this scope keeps that
       order and takes the first ``max_transition_states_per_round``.

    Args:
        split (CandidateSplit): The split from ``split_qm_candidates``.
        evidence_by_ts_label (dict): Network-local TS label -> ``(coefficient, delta_ln_k)``, from
            ``t3.pdep.pes_sa.run_round_me_sensitivity``.
        network_e0 (PDepNetworkE0): The declared species/TS E0 values (kJ/mol), from the same
            network file.
        params (DistrustParams): The window and any other knobs.

    Returns:
        CandidateSplit: The in-window candidates in descending distrust order, each stamped with its
        coefficient/``delta_ln_k`` and ``distrust_score``; every other candidate appended to
        ``skipped`` with its reason and classification.
    """
    skipped = list(split.skipped)
    stamped = []
    for candidate in split.candidates:
        if classify_provenance(candidate.path_reaction.kinetics_comment) == PROVENANCE_LIBRARY:
            # A trusted library/QM value: distrust does not re-queue it. Skipped before ranking so an
            # all-library network queues nothing rather than spending QM on what it already holds.
            skipped.append(SkippedChannel(
                label=candidate.path_reaction.label,
                reason=(f"'{candidate.path_reaction.label}': transition state {candidate.ts_label} "
                        f'already carries a trusted library / already-computed (QM) rate, which '
                        f'distrust does not re-queue; computing it again would be redundant.'),
                ts_label=candidate.ts_label, classification=SKIP_TRUSTED_PROVENANCE))
            continue
        pair = evidence_by_ts_label.get(candidate.ts_label)
        if pair is None:
            # No finite sensitivity row: capture would refuse the artifact, so this cannot be
            # queued regardless of distrust. Same skip the sensitivity path takes; the classification
            # only records whether the absence was the structural kind.
            classification = (SKIP_UNMEASURABLE if not candidate.e0_sensitivity_measurable
                              else SKIP_NO_EVIDENCE)
            skipped.append(SkippedChannel(
                label=candidate.path_reaction.label,
                reason=(f"'{candidate.path_reaction.label}': transition state {candidate.ts_label} "
                        f'has no finite sensitivity row in this round, and a captured artifact must '
                        f'carry the coefficient that will be recorded for it (t3.pdep.capture); not '
                        f'queueing it rather than inventing a number.'),
                ts_label=candidate.ts_label, classification=classification))
            continue
        coefficient, delta_ln_k = pair
        stamped.append(replace(candidate, coefficient=coefficient, delta_ln_k=delta_ln_k))

    eligible, declined = rank_candidates_by_distrust(tuple(stamped), network_e0, params)
    for candidate, score in declined:
        skipped.append(SkippedChannel(
            label=candidate.path_reaction.label,
            reason=_outside_window_reason(candidate, score, params),
            ts_label=candidate.ts_label, classification=SKIP_OUTSIDE_WINDOW,
            coefficient=candidate.coefficient, delta_ln_k=candidate.delta_ln_k))
    ranked = tuple(replace(candidate, distrust_score=score.score) for candidate, score in eligible)
    return CandidateSplit(candidates=ranked, skipped=tuple(skipped))
