"""
t3 pdep me_success module

Decides whether a completed Arkane pressure-dependence (ME, master equation) job's
``output.py`` represents a genuine ME solve, as opposed to a hard crash or a soft failure that
Arkane itself silently accepts as "successful" (exit code 0, empty stderr) while writing out a
kinetics result that is not actually usable.

This module never executes ``output.py`` (no ``exec``/``eval``/``import``): it is parsed with
``t3.pdep.parser.parse_arkane_pdep_output_file``, which uses ``ast.parse`` only, for the same
reasons documented in ``t3/pdep/parser.py``.
"""

import math
import os
import re
from dataclasses import dataclass

from t3.pdep.parser import parse_arkane_pdep_output_file

# Stderr noise that Arkane and its dependencies routinely emit that does not indicate an ME
# solve failure. This mirrors ARC's own ``ignorable_phrases`` filter (see
# ``arc/statmech/arkane.py::run_arkane``) rather than inventing a new filtering policy here.
IGNORABLE_STDERR_PHRASES = (
    'Open Babel Warning',
    'Accepted unusual valence',
    '==============================',
    'pjrt_executable.cc',
)

# A Python ``warnings``-module emission looks like ``<file>:<lineno>: SomeWarning: <message>``,
# optionally followed by an indented echo of the offending source line. The explorer's ME solve
# routinely prints scipy ``LinAlgWarning: Ill-conditioned matrix`` and RMG ``thermoengine``
# ``LinAlgWarning`` this way on a run that otherwise succeeds and writes a valid, count-consistent
# ``output.py``. A warning is non-fatal BY CONSTRUCTION (the warnings module never raises), so it is
# not evidence a run failed -- and the real failure signals remain untouched: a non-zero exit, a
# missing required artifact, a fatal ``Error:``/``Critical:`` log marker, and the payload checks. A
# traceback line (``Traceback (most recent call last):``, ``ValueError: ...``) does NOT match this
# pattern and still fails the run, so this filters noise without deleting the guard.
_WARNING_LINE_RE = re.compile(r':\d+:\s+\w*Warning:\s')


def real_stderr_lines(stderr_text: str) -> list:
    """
    Return the stderr lines that are genuine failure evidence, dropping ignorable noise.

    Drops blank lines, every ``IGNORABLE_STDERR_PHRASES`` line, and Python ``warnings``-module
    output -- both the ``file:line: SomeWarning: message`` line and the single indented source-echo
    line the warnings module prints after it. Everything else (tracebacks, error text) is kept,
    stripped, in order.

    Args:
        stderr_text (str): The raw stderr content.

    Returns:
        list: The real (non-ignorable) stderr lines, each stripped of surrounding whitespace.
    """
    real = []
    prev_was_warning = False
    for line in stderr_text.splitlines():
        stripped = line.strip()
        if not stripped:
            prev_was_warning = False
            continue
        if any(phrase in line for phrase in IGNORABLE_STDERR_PHRASES):
            prev_was_warning = False
            continue
        if _WARNING_LINE_RE.search(line):
            prev_was_warning = True
            continue
        # The warnings module echoes the offending source line, indented, immediately after the
        # warning itself. Drop that echo, but only when it directly follows a warning line -- an
        # indented line elsewhere (e.g. inside a traceback) is left as real evidence.
        if prev_was_warning and line[:1].isspace():
            prev_was_warning = False
            continue
        prev_was_warning = False
        real.append(stripped)
    return real


@dataclass(frozen=True)
class MESuccessResult:
    """
    The outcome of checking whether an Arkane pdep ME job succeeded.

    Args:
        succeeded (bool): Whether ``output.py`` represents a genuine ME solve.
        reasons (tuple): Human-readable failure reasons, in the order they were detected.
            Empty when ``succeeded`` is True.
        reactions (tuple): The parsed ``PDepArkaneReaction`` entries (from
            ``t3.pdep.parser``), empty if none could be parsed.
    """
    succeeded: bool
    reasons: tuple
    reactions: tuple


def check_arkane_me_success(output_path: str,
                             exit_code: int | None = None,
                             stderr: str | None = None,
                             expected_reactions: int | set | None = None,
                             ) -> MESuccessResult:
    """
    Decide whether an Arkane pdep ME job genuinely succeeded, and report why not when it did not.

    Checks, in order (each stops the check and returns as soon as it fails, since later checks
    are meaningless once an earlier one has already failed):

    1. ``output.py`` exists.
    2. It is non-empty. Arkane pre-creates ``output.py`` (opened ``'w'``) before running any
       job and only fills it in afterward, so bare existence proves nothing: a hard crash
       inside ``job.execute()`` leaves a valid, 0-byte ``output.py`` and a traceback on stderr.
    3. It parses to at least one ``pdepreaction(...)`` entry.
    4. Every parsed entry names non-empty reactants AND non-empty products, and carries a
       non-empty RATE PAYLOAD (``coeffs``, ``A``/``n``/``Ea``, a nested ``arrhenius=[...]``
       list, ...) whose numeric leaves are all finite (no ``None``, no ``nan``, no ``inf``).
       The payload is checked separately from the Tmin/Tmax/Pmin/Pmax bounds and other
       metadata, which are always finite in Arkane's output and must never satisfy this check
       on their own. This is the ONLY check that catches the observed soft failure: with
       ``method='chemically-significant eigenvalues'``, Arkane can reject a negative rate
       coefficient, print "Error: Negative rate coefficient generated; rejecting result." to
       *stdout* only, exit 0 with empty stderr, and still write a syntactically perfect
       ``pdepreaction(...)`` whose ``Chebyshev(coeffs=[[None, None, ...], ...])`` is entirely
       ``None``. Positivity of the coefficients themselves is deliberately NOT required: Chebyshev
       coefficients are log-space fit coefficients and are legitimately negative; only
       finiteness is checked here, never sign.
    5. (Optional, off by default) The parsed entries cover an expected reaction count or set.

    The subprocess exit status and stderr, if the caller has them, are additional evidence that
    gets folded into ``reasons`` -- never a substitute for the payload checks above, since the
    soft-failure case above has both a zero exit code and empty stderr. They are, however, still
    sufficient on their own to FAIL the job: Arkane can write a complete result for one network
    and then die later in the run, so a well-formed payload does not vouch for the run as a whole.

    Args:
        output_path (str): Path to the Arkane pdep job's ``output.py`` file.
        exit_code (int, optional): The Arkane subprocess's exit status, if available.
        stderr (str, optional): The Arkane subprocess's stderr text, if available. Benign lines
            (see ``IGNORABLE_STDERR_PHRASES``) are filtered out before being folded into reasons.
        expected_reactions (int or set, optional): If an ``int``, the exact number of
            ``pdepreaction(...)`` entries expected. If a ``set``, a set of
            ``(reactants, products)`` tuples (each itself a tuple of species labels) that must
            all be present among the parsed entries. Defaults to ``None`` (no check).

    Returns:
        MESuccessResult: The success/failure verdict, with reasons and the parsed reactions.
    """
    reasons = []

    if exit_code is not None and exit_code != 0:
        reasons.append(f'Arkane exited with a non-zero status ({exit_code}).')
    if stderr:
        real_errors = real_stderr_lines(stderr)
        if real_errors:
            reasons.append('Arkane reported stderr output: ' + ' | '.join(real_errors))

    # Check 1: existence.
    if not os.path.isfile(output_path):
        reasons.append(f"Arkane output file '{output_path}' does not exist.")
        return MESuccessResult(succeeded=False, reasons=tuple(reasons), reactions=tuple())

    # Check 2: non-empty (catches the hard-crash case; existence alone is worthless).
    if os.path.getsize(output_path) == 0:
        reasons.append(
            f"Arkane output file '{output_path}' is empty (0 bytes); Arkane pre-creates this "
            f"file before running any job, so this means the ME job crashed before writing "
            f"any result."
        )
        return MESuccessResult(succeeded=False, reasons=tuple(reasons), reactions=tuple())

    # Check 3: parses to at least one pdepreaction(...) entry.
    try:
        reactions = parse_arkane_pdep_output_file(output_path)
    except ValueError as e:
        reasons.append(str(e))
        return MESuccessResult(succeeded=False, reasons=tuple(reasons), reactions=tuple())

    if not reactions:
        reasons.append(f"Arkane output file '{output_path}' contains no pdepreaction(...) entries.")
        return MESuccessResult(succeeded=False, reasons=tuple(reasons), reactions=reactions)

    # Check 4: every entry has a valid identity (non-empty reactants AND products) and a
    # non-empty, entirely finite RATE PAYLOAD. The payload check deliberately uses
    # ``rate_payload_numeric_values`` and never the combined ``numeric_values``: the latter also
    # carries the (always finite) Tmin/Tmax/Pmin/Pmax bounds, which used to let a payload-free
    # ``Chebyshev(coeffs=[], Tmin=..., ...)`` -- or one whose unparseable ``coeffs`` the parser
    # omitted entirely -- pass on its bounds alone, with no rate coefficients at all.
    bad_entry_found = False
    for i, reaction in enumerate(reactions):
        reaction_label = f"{' + '.join(reaction.reactants)} <=> {' + '.join(reaction.products)}"
        if not reaction.reactants or not reaction.products:
            # A pdepreaction without both sides is not a valid k(T,P) result no matter how
            # finite its numbers are.
            bad_entry_found = True
            reasons.append(
                f"Reaction #{i} ({reaction_label}) has empty "
                f"{'reactants' if not reaction.reactants else 'products'}; a pdepreaction entry "
                f"must name both non-empty reactants and non-empty products."
            )
            continue
        if not reaction.rate_payload_numeric_values:
            # Guard against the finiteness check passing vacuously: no payload leaves at all
            # (an empty ``coeffs=[]``, or a kinetics form the parser had to omit) carries no
            # rate information, so it must not be read as "all of its values are finite".
            bad_entry_found = True
            missing = (f" (kinetics keywords that could not be parsed as literals: "
                       f"{list(reaction.missing_kinetics_keys)})"
                       if reaction.missing_kinetics_keys else '')
            reasons.append(
                f"Reaction #{i} ({reaction_label}) has no numeric {reaction.kinetics_type} "
                f"rate coefficients{missing}; finite T/P bounds alone do not constitute a "
                f"usable ME solve result."
            )
            continue
        bad_values = [v for v in reaction.rate_payload_numeric_values
                      if v is None or not math.isfinite(v)]
        if bad_values:
            bad_entry_found = True
            shown = bad_values[:5]
            more = f' and {len(bad_values) - 5} more' if len(bad_values) > 5 else ''
            reasons.append(
                f"Reaction #{i} ({reaction_label}) has non-finite {reaction.kinetics_type} "
                f"kinetics parameters: {shown}{more}."
            )
    if bad_entry_found:
        return MESuccessResult(succeeded=False, reasons=tuple(reasons), reactions=reactions)

    # Check 5 (optional): expected reaction coverage.
    if expected_reactions is not None:
        coverage_ok = True
        if isinstance(expected_reactions, int):
            if len(reactions) != expected_reactions:
                coverage_ok = False
                reasons.append(
                    f'Expected {expected_reactions} pdepreaction(...) entries, found {len(reactions)}.'
                )
        else:
            found = {(r.reactants, r.products) for r in reactions}
            missing = set(expected_reactions) - found
            if missing:
                coverage_ok = False
                reasons.append(f'Expected pdepreaction(...) entries not found in output: {sorted(missing)}.')
        if not coverage_ok:
            return MESuccessResult(succeeded=False, reasons=tuple(reasons), reactions=reactions)

    # Any process-level evidence collected up front (a non-zero exit status, real stderr output)
    # still fails the job even though the payload itself looks well-formed: Arkane can write a
    # complete result for one network and then die, so a well-formed payload does not vouch for
    # the run as a whole. Returning ``succeeded=True`` here while holding such a reason would
    # discard the only record that anything went wrong.
    return MESuccessResult(succeeded=not reasons, reasons=tuple(reasons), reactions=reactions)
