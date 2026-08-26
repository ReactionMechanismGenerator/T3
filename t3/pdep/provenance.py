"""
t3 pdep provenance module

Classify WHERE the energy of a stationary point on a P-dep surface came from, so T3's PES diagram
can colour every level by the trustworthiness of its energy rather than by its type. Four
provenance classes, ordered by decreasing confidence:

1. ``PROVENANCE_QM``      -- computed by quantum chemistry in this campaign.
2. ``PROVENANCE_LIBRARY`` -- taken from a curated thermo library, or a kinetics library via ILT.
3. ``PROVENANCE_FAMILY``  -- estimated from a reaction-family template or tree node (ILT).
4. ``PROVENANCE_GAV``     -- thermo group-additivity estimation.

plus ``PROVENANCE_UNKNOWN`` for a level whose provenance genuinely cannot be read -- rendered
distinctly and named in the legend, never silently folded into one of the four.

The classifier reads the ``comment=`` text RMG/Arkane write into a species' ``thermo=NASA(...)`` or
a reaction's ``kinetics=Arrhenius(...)`` call in the network file. The one trap those comments set
is that they are **concatenated, not replaced**: an RMG estimate is written first and this run's
quantum-chemistry result is appended AFTER it and overrides it. Concretely, a QM'd well reads

    Thermo group additivity estimation: group(...)...Thermo library: thermojobs

and is QM, not GAV -- ``thermojobs`` (thermo) and ``kineticsjobs`` (kinetics) are ARC/Arkane's own
output libraries, this run's quantum chemistry written back as a library. A first-match parse would
label it GAV, the exact inversion the diagram exists to prevent. So the classifier finds every
provenance marker in the comment and lets the **last** one (largest start index) win.
"""

import re
from collections.abc import Iterable

PROVENANCE_QM = 'qm'
PROVENANCE_LIBRARY = 'library'
PROVENANCE_FAMILY = 'family'
PROVENANCE_GAV = 'gav'
PROVENANCE_UNKNOWN = 'unknown'

# Higher means more trustworthy. Used to pick the worst-link provenance of a multi-species channel
# (a bimolecular asymptote's energy is only as trustworthy as its least-trustworthy fragment) and
# to order the legend. ``PROVENANCE_UNKNOWN`` is the floor so a channel with any unreadable
# fragment is reported unknown rather than given a fragment's class it does not deserve.
PROVENANCE_CONFIDENCE = {
    PROVENANCE_QM: 4,
    PROVENANCE_LIBRARY: 3,
    PROVENANCE_FAMILY: 2,
    PROVENANCE_GAV: 1,
    PROVENANCE_UNKNOWN: 0,
}

# Human-readable class names for the diagram legend.
PROVENANCE_LABELS = {
    PROVENANCE_QM: 'QM (this campaign)',
    PROVENANCE_LIBRARY: 'Library',
    PROVENANCE_FAMILY: 'Family estimate (ILT)',
    PROVENANCE_GAV: 'GAV estimate',
    PROVENANCE_UNKNOWN: 'Unknown provenance',
}

# The library names ARC/Arkane give this run's own quantum-chemistry output when it writes it back
# into the network file as a "library". A ``Thermo library``/``Reaction library`` marker naming one
# of these is QM; any other name is a genuinely external, curated library. These are ARC's standard
# output-library names, independent of T3's configured ``library_name``; extend the set if a
# campaign is found to use a different one.
_QM_OUTPUT_LIBRARY_NAMES = frozenset({'thermojobs', 'kineticsjobs'})

# One entry per provenance marker RMG/Arkane emits into a comment: (kind, keyword regex). The
# regex matches only the marker KEYWORD, never the trailing name -- capturing the name greedily
# breaks on concatenated comments where the next marker abuts the previous name with no separator
# (e.g. ``primaryThermoLibraryThermo library: thermojobs``). The name is read separately, only for
# the marker that wins, by ``_library_name_after``.
_MARKERS = (
    (PROVENANCE_GAV, re.compile(r'group additivity estimation')),
    ('_thermo_library', re.compile(r'Thermo library:')),
    ('_reaction_library', re.compile(r'Reaction library:')),
    (PROVENANCE_FAMILY, re.compile(r'Estimated (?:using template|from node)')),
)
_LIBRARY_KINDS = ('_thermo_library', '_reaction_library')


def _library_name_after(comment: str, marker_end: int) -> str:
    """
    Read the library name that follows a ``Thermo library:``/``Reaction library:`` marker.

    The name is the first whitespace/quote-delimited token after the marker keyword; any
    surrounding quotes (``Reaction library: 'kineticsjobs'``) are stripped. This is only ever
    called for the LAST marker in the comment, so nothing recognised follows the name.

    Args:
        comment (str): The full comment text.
        marker_end (int): The index just past the marker keyword.

    Returns:
        str: The library name, or ``''`` if none could be read.
    """
    match = re.match(r"\s*['\"]?([^\s'\",]+)", comment[marker_end:])
    return match.group(1) if match else ''


def classify_provenance(comment: str | None,
                        qm_library_names: frozenset[str] = _QM_OUTPUT_LIBRARY_NAMES) -> str:
    """
    Classify the provenance of an energy from the RMG/Arkane ``comment`` that produced it.

    Works for both a species' ``thermo`` comment and a reaction's ``kinetics`` comment: the marker
    set spans both, and no thermo marker collides with a kinetics one. Classification is
    **last-writer-wins** -- every provenance marker in the (concatenated) comment is located and the
    one starting latest decides the class, so this run's quantum chemistry (appended last) overrides
    the RMG estimate it superseded. See the module docstring for the trap this defends against.

    Args:
        comment (str | None): The ``comment`` text, or ``''``/``None`` when the level carries none.
        qm_library_names (frozenset[str]): The library names that denote this run's own QM output
            (see ``_QM_OUTPUT_LIBRARY_NAMES``).

    Returns:
        str: One of the ``PROVENANCE_*`` constants. ``PROVENANCE_UNKNOWN`` when the comment is
            empty or carries no recognised marker -- an honest gap, never defaulted into a class.
    """
    if not comment or not comment.strip():
        return PROVENANCE_UNKNOWN

    winner = None  # (start, kind, end) of the latest-starting marker
    for kind, pattern in _MARKERS:
        for match in pattern.finditer(comment):
            if winner is None or match.start() > winner[0]:
                winner = (match.start(), kind, match.end())
    if winner is None:
        return PROVENANCE_UNKNOWN

    _, kind, end = winner
    if kind in _LIBRARY_KINDS:
        name = _library_name_after(comment, end)
        return PROVENANCE_QM if name in qm_library_names else PROVENANCE_LIBRARY
    return kind  # PROVENANCE_GAV or PROVENANCE_FAMILY


def combine_channel_provenance(provenances: Iterable[str]) -> str:
    """
    The provenance of a (possibly bimolecular) channel from its fragments' provenances.

    A channel's asymptote energy is the sum over its fragments, so it is only as trustworthy as its
    least-trustworthy fragment: the worst (lowest-confidence) class wins. An empty channel -- which
    should never occur -- is ``PROVENANCE_UNKNOWN``.

    Args:
        provenances (Iterable[str]): The per-fragment ``PROVENANCE_*`` classes.

    Returns:
        str: The worst-link ``PROVENANCE_*`` class.
    """
    classes = list(provenances)
    if not classes:
        return PROVENANCE_UNKNOWN
    return min(classes, key=lambda provenance: PROVENANCE_CONFIDENCE[provenance])
