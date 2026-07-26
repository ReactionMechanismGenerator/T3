#!/usr/bin/env python3
# encoding: utf-8

"""
t3.utils.uncertainty module

A shared predicate for determining whether a reaction's kinetics are uncertain enough to
warrant refinement (e.g., a QM calculation). This logic is used by
``T3.reaction_requires_refinement`` (which additionally applies a dedup gate for
previously-considered reactions), and by callers that only have a bare kinetics comment
string, such as the RMG PDep network file parser, with no ``T3Reaction`` object to work with.

A reaction is considered certain (not uncertain) if:

- its kinetics method is one of ``LIBRARY`` or ``TRAINING_SET`` (a method-level short
  circuit, regardless of comment); or
- its comment carries a "certain" provenance marker (an exact rate-rule match, a matched
  training reaction, a library provenance statement, or a seed mechanism statement) and is
  not itself wrapped in an "estimate" qualifier (e.g., ``Estimated using ...``,
  ``Estimated from node ...``, a BM rule fitted to training reactions, or an Ea raised to
  match the endothermicity of the reaction).

Any other kinetics method (``RATE_RULES``, ``PDEP``, ``QM``, ``USER``, ``UNKNOWN``, or
``None``) falls through to comment analysis for a forced verdict. Note that ``QM`` and
``USER`` are intentionally NOT method-level short circuits here: T3's own kinetics parser
(``t3.utils.cantera_parser``) only ever assigns ``LIBRARY``, ``RATE_RULES``, or ``PDEP``, so
``QM``/``USER`` values seen here originate from user input or restart state, not from a fresh
parse. Short-circuiting them would mean a failed, placeholder, or stale QM/user reaction
could never be re-queued for refinement.

Before scanning for markers, lines that are per-path-reaction PDep debug annotations (i.e.,
that start with ``High-P limit:`` after stripping leading whitespace and an optional leading
``'! '``) are dropped from the comment. These are emitted once per path reaction of a PDep
network when RMG runs at DEBUG level, and must not decide the net PDep reaction's verdict in
either direction: a single library path reaction must not make an otherwise-uncertain PDep
net reaction certain, nor should an estimated path reaction make it uncertain.
"""

import re
from typing import TYPE_CHECKING

from t3.chem import KineticsMethod

if TYPE_CHECKING:
    from t3.chem import T3Reaction


# Kinetics methods always considered certain, regardless of the kinetics comment.
# NOTE: intentionally excludes QM and USER -- see the module docstring for why.
CERTAIN_KINETICS_METHODS = (KineticsMethod.LIBRARY, KineticsMethod.TRAINING_SET)

# Raw string aliases (normalized to lowercase, whitespace-collapsed) mapped to a ``KineticsMethod``.
# Used to normalize bare-comment callers' (e.g., the RMG PDep network file parser) raw method
# strings, which may not exactly match a ``KineticsMethod`` enum value.
#
# WARNING: this alias map (in particular the ``'training'`` alias below) is only safe when the
# caller passes an actual kinetics-method FIELD, never text pulled from a comment. RMG comment
# text contains phrases like ``From training reaction <n> used ...``, which denotes an ESTIMATE
# (a rate rule fitted from a training reaction and then applied by analogy to this reaction) and
# must stay uncertain -- as opposed to ``Matched reaction <n> ... in <Family>/training``, which
# means this reaction IS the training reaction and is certain. Do not derive ``kinetics_method``
# from comment text.
KINETICS_METHOD_ALIASES = {
    'library': KineticsMethod.LIBRARY,
    'library reaction': KineticsMethod.LIBRARY,
    'reaction library': KineticsMethod.LIBRARY,
    'training set': KineticsMethod.TRAINING_SET,
    'training': KineticsMethod.TRAINING_SET,
    'rate rules': KineticsMethod.RATE_RULES,
    'rate rule': KineticsMethod.RATE_RULES,
    'qm': KineticsMethod.QM,
    'user': KineticsMethod.USER,
    'pdep': KineticsMethod.PDEP,
    'unknown': KineticsMethod.UNKNOWN,
}

# Comment-level markers.
ESTIMATED_USING_MARKER = 'Estimated using'
ESTIMATED_FROM_NODE_MARKER = 'Estimated from node'
BM_RULE_FITTED_MARKER = 'BM rule fitted to'
ENDOTHERMICITY_MARKER = 'to match endothermicity of reaction'
EXACT_MATCH_MARKER = 'Exact match found'
LIBRARY_REACTION_MARKER = 'Library reaction:'
REACTION_LIBRARY_MARKER = 'reaction library'
SEED_MECHANISM_MARKER = 'Seed mechanism:'

# Per-path-reaction PDep debug annotation prefix (emitted once per path reaction when RMG runs
# at DEBUG level; see rmgpy/chemkin.pyx around lines 1803-1808). Lines that start with this
# (after stripping leading whitespace and an optional leading '! ') are dropped from the
# comment before marker scanning, in both directions.
HIGH_P_LIMIT_PREFIX = 'High-P limit:'

# ``Matched reaction <n> ... /training`` bound to a short span: RMG's own parser
# (rmgpy/data/kinetics/family.py, ``training_reaction_pattern``) matches this on a single line,
# but some fixture comments split ``... in\`` and ``<Family>/training`` across a literal line
# continuation, so we tolerate a bounded amount of intervening text (non-greedy, capped) rather
# than matching across an unbounded distance.
MATCHED_TRAINING_REGEX = re.compile(r'Matched reaction\s+\d+[\s\S]{0,200}?/training')


def _strip_high_p_limit_lines(kinetics_comment: str) -> str:
    """
    Remove per-path-reaction ``High-P limit: ...`` PDep debug annotation lines from
    ``kinetics_comment``, so that a single path reaction's provenance cannot decide the net
    PDep reaction's verdict in either direction.

    Args:
        kinetics_comment (str): The raw kinetics comment.

    Returns:
        str: The comment with any ``High-P limit:`` lines removed.
    """
    kept_lines = list()
    for line in kinetics_comment.splitlines():
        stripped = line.strip()
        if stripped.startswith('! '):
            stripped = stripped[2:]
        if stripped.startswith(HIGH_P_LIMIT_PREFIX):
            continue
        kept_lines.append(line)
    return '\n'.join(kept_lines)


def _normalize_kinetics_method(kinetics_method: 'KineticsMethod | str | None') -> 'KineticsMethod | None':
    """
    Normalize a ``kinetics_method`` argument (a ``KineticsMethod`` enum instance, a raw string, or
    ``None``) into a ``KineticsMethod`` enum instance, or ``None`` if it cannot be recognized.

    String matching is case-insensitive and whitespace-tolerant, and accepts common aliases (e.g.,
    ``'library'``, ``'Library reaction'``, ``'reaction library'`` all normalize to ``LIBRARY``).
    An unrecognized string does not raise; it normalizes like ``None`` (treated as unknown).

    Args:
        kinetics_method (KineticsMethod, str, optional): The reaction's kinetics method, either as a
                                                         ``KineticsMethod`` enum instance or a raw string.

    Returns:
        Optional[KineticsMethod]: The normalized ``KineticsMethod``, or ``None`` if unrecognized.
    """
    if isinstance(kinetics_method, KineticsMethod):
        return kinetics_method
    if kinetics_method is None:
        return None
    key = ' '.join(str(kinetics_method).lower().split())
    if key in KINETICS_METHOD_ALIASES:
        return KINETICS_METHOD_ALIASES[key]
    try:
        return KineticsMethod(kinetics_method)
    except ValueError:
        return None


def is_this_kinetics_comment_uncertain(kinetics_comment: str | None,
                                       kinetics_method: 'KineticsMethod | str | None' = None,
                                       ) -> bool:
    """
    Determine whether a kinetics comment (and, optionally, a kinetics method) indicates that a
    rate coefficient is uncertain and should be refined.

    A reaction is considered certain (not uncertain) if its kinetics method is one of ``LIBRARY``
    or ``TRAINING_SET`` (a method-level short circuit that does not depend on the comment), or if
    its comment carries a certain provenance marker (exact rate-rule match, matched training
    reaction, library provenance statement, or seed mechanism statement) that is not itself
    wrapped in an estimate qualifier.

    Any other kinetics method (``RATE_RULES``, ``PDEP``, ``QM``, ``USER``, ``UNKNOWN``, or
    ``None``) falls through to comment analysis for a forced verdict. See the module docstring
    for why ``QM``/``USER`` are not method-level short circuits here.

    Before marker scanning, per-path-reaction ``High-P limit: ...`` PDep debug annotation lines
    are dropped from the comment (see the module docstring).

    Args:
        kinetics_comment (str, optional): The reaction's RMG kinetics comment.
        kinetics_method (KineticsMethod, str, optional): The reaction's kinetics method, either as a
                                                         ``KineticsMethod`` enum instance or a raw
                                                         string value (e.g., ``'Library'``). Callers
                                                         with only a bare comment string (e.g., the
                                                         RMG PDep network file parser) will not have
                                                         an enum instance to pass here; a raw string
                                                         (case-insensitively and whitespace-tolerantly
                                                         matched, including common aliases) is accepted.
                                                         An unrecognized string does not raise; it falls
                                                         through to comment analysis, the same as ``None``.

    Returns:
        bool: Whether the reaction's kinetics are uncertain and should be refined. ``True`` if so.
    """
    method = _normalize_kinetics_method(kinetics_method)
    if method in CERTAIN_KINETICS_METHODS:
        return False
    kinetics_comment = kinetics_comment or ''
    kinetics_comment = _strip_high_p_limit_lines(kinetics_comment)
    has_estimate = (
        ESTIMATED_USING_MARKER in kinetics_comment
        or ESTIMATED_FROM_NODE_MARKER in kinetics_comment
        or BM_RULE_FITTED_MARKER in kinetics_comment
        or ENDOTHERMICITY_MARKER in kinetics_comment
    )
    has_certain_marker = (
        EXACT_MATCH_MARKER in kinetics_comment
        or MATCHED_TRAINING_REGEX.search(kinetics_comment) is not None
        or LIBRARY_REACTION_MARKER in kinetics_comment
        or REACTION_LIBRARY_MARKER in kinetics_comment.lower()
        or SEED_MECHANISM_MARKER in kinetics_comment
    )
    return not (has_certain_marker and not has_estimate)


def is_this_reaction_uncertain(reaction: 'T3Reaction | None') -> bool | None:
    """
    Determine whether a ``T3Reaction``'s kinetics are uncertain and should be refined.

    A thin wrapper around ``is_this_kinetics_comment_uncertain()`` that extracts the relevant
    fields from a ``T3Reaction`` instance. Does not apply any dedup gate for
    previously-considered reactions; that remains the caller's responsibility.

    Args:
        reaction (T3Reaction, optional): The reaction to query.

    Returns:
        Optional[bool]: Whether the reaction's kinetics are uncertain and should be refined.
                         ``True`` if so, ``None`` if ``reaction`` is ``None``.
    """
    if reaction is None:
        return None
    return is_this_kinetics_comment_uncertain(kinetics_comment=reaction.kinetics_comment,
                                               kinetics_method=reaction.kinetics_method,
                                               )
