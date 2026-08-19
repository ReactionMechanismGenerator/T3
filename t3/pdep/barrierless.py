"""
t3 pdep barrierless module

Which path reactions must never be sent to a transition-state search.

A barrierless reaction has no saddle point on its potential energy surface, so a TS search for it
cannot succeed no matter how many adapters are tried or how long they run. The group's own
reference note is unambiguous about the canonical case:

    R_Recombination reactions are barrierless -- there is no saddle point and ARC cannot compute
    them. Exclude during generation and record the reason.
    -- Vault/Code/ARC/Canonical Levels of Theory.md

This is not a theoretical concern. On the CHO2 surface T3 admitted ``HO + CHO <=> CH2O2`` -- two
radicals recombining to a closed-shell molecule -- for TS search via
``ts_adapters_for_unknown_unimolecular``, because ``ARCReaction.is_unimolecular()`` is
direction-agnostic and a bimolecular->unimolecular recombination therefore qualifies. One adapter
dutifully returned three guesses for a saddle point that does not exist, and they were optimized
on a cluster.

An automatically exploring loop makes this failure mode structural rather than incidental:
exploration generates recombination channels faster than any other kind, so a loop without this
gate spends its entire QM budget on reactions that cannot converge.

**This module fails OPEN.** An unrecognized family is reported as *not* barrierless and gets its
TS search. Wrongly skipping a real barrier silently loses physics from the surface; wrongly
attempting a barrierless one wastes one job and says so in the log. The asymmetry is deliberate,
and it is why the barrierless set is an explicit allowlist rather than a heuristic.

The zero-activation-energy signal is deliberately NOT used as a classifier. Every barrierless
Arrhenius fit in an RMG network carries ``Ea=(0,'kJ/mol')``, but so do plenty of barriered
reactions whose fitted Ea collapsed to zero, and acting on that would fail CLOSED on real
chemistry. Family is the only signal here.
"""

import re
from dataclasses import dataclass

from t3.pdep.parser import PDepPathReaction

# RMG families whose reactions proceed without a saddle point. Additions belong here rather than
# at a call site, so that every consumer of this module agrees on the set.
BARRIERLESS_FAMILIES = frozenset({
    'R_Recombination',
    'Birad_recombination',
    'Diradical_formation',
})

# RMG writes the family into the kinetics comment in exactly two forms, both seen in the CHO2
# network files this module was written against:
#   "Estimated from node Root_N-1R->H_N-1CNOS->N in family R_Recombination."
#   "family: 1,2_Insertion_CO"
# A family name is not [A-Za-z_]+ -- '1,2_Insertion_CO' and '1,3_Insertion_CO2' both start with a
# digit and contain a comma -- so the character class has to admit those or the two insertion
# families silently come back as None.
_FAMILY_CHARS = r"[A-Za-z0-9_,\-]+"
_FAMILY_IN_NODE_RE = re.compile(rf'\bin family ({_FAMILY_CHARS})')
_FAMILY_LABELLED_RE = re.compile(rf'^family:\s*({_FAMILY_CHARS})\s*$', re.MULTILINE)


@dataclass(frozen=True)
class BarrierlessVerdict:
    """
    Why a path reaction may or may not be sent to a TS search.

    Attributes:
        is_barrierless (bool): Whether the reaction has no saddle point and must not be queued.
        family (str | None): The RMG family the verdict was made from, or ``None`` if the kinetics
                             comment named none.
        reason (str): A human-readable justification, recorded alongside the decision so a skipped
                      channel is never silently skipped.
    """
    is_barrierless: bool
    family: str | None
    reason: str


def rmg_family(path_reaction: PDepPathReaction) -> str | None:
    """
    Extract the RMG family name from a path reaction's kinetics comment.

    Args:
        path_reaction (PDepPathReaction): The parsed path reaction.

    Returns:
        str | None: The family name, or ``None`` if the comment names none. Never guessed at:
                    a comment without a family yields ``None`` rather than a best effort.
    """
    comment = path_reaction.kinetics_comment or ''
    match = _FAMILY_IN_NODE_RE.search(comment)
    if match is not None:
        return match.group(1).rstrip('.')
    match = _FAMILY_LABELLED_RE.search(comment)
    if match is not None:
        return match.group(1).rstrip('.')
    return None


def classify_barrierless(path_reaction: PDepPathReaction) -> BarrierlessVerdict:
    """
    Decide whether a path reaction is barrierless and must be kept out of QM.

    Args:
        path_reaction (PDepPathReaction): The parsed path reaction.

    Returns:
        BarrierlessVerdict: The verdict, carrying the family it was made from and a reason string
                            suitable for logging and for the round record.
    """
    family = rmg_family(path_reaction)
    if family is None:
        return BarrierlessVerdict(
            is_barrierless=False, family=None,
            reason=f"'{path_reaction.label}': the kinetics comment names no RMG family, so no "
                   f'barrierless determination can be made; queueing it for TS search.')
    if family in BARRIERLESS_FAMILIES:
        return BarrierlessVerdict(
            is_barrierless=True, family=family,
            reason=f"'{path_reaction.label}': family {family} is barrierless -- it has no saddle "
                   f'point, so a TS search cannot converge. Taking RMG/ILT kinetics instead.')
    return BarrierlessVerdict(
        is_barrierless=False, family=family,
        reason=f"'{path_reaction.label}': family {family} is not a known barrierless family; "
               f'queueing it for TS search.')
