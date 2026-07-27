"""
t3 utils network_thermo module

This is a leaf module (it imports nothing from ``t3.*``): ``t3.utils.writer`` -- and, most of
``t3.pdep`` -- need to compute a pdep network's species thermo ceiling without pulling in
``t3.pdep`` itself. ``t3.pdep.parser`` re-exports ``network_thermo_t_max`` since it is a
network-file parsing concern callers reasonably look for there too, but the real home of this
logic is ``t3.utils`` so that the ``t3.pdep`` <-> ``t3.utils.writer`` import direction stays
one-way.
"""

import ast
import warnings
from dataclasses import dataclass


class NetworkTextUnparseable(ValueError):
    """Raised by ``network_thermo_t_max`` only when ``text`` itself does not parse as Python.

    This is deliberately a distinguishable subclass of ``ValueError`` (not a plain
    ``ValueError``) so callers that want to degrade "the text was corrupted" into "clamp
    nothing" -- deferring to a downstream self-check that owns reporting that corruption -- can
    catch exactly that case, without also silently swallowing a ``ValueError`` raised for a
    different, non-degradable reason (e.g. a malformed ``**kwargs`` unpacking in a
    ``species(...)`` call, which is a real defect in the network text that a caller should not
    paper over as "no species contributed a thermo ceiling").
    """


def _get_call_name(call: ast.Call) -> str | None:
    """
    Get the callee name of an ``ast.Call`` node (e.g., ``'NASA'`` for ``NASA(...)``).

    Only a bare ``ast.Name`` callee yields a name, mirroring
    ``t3.pdep.parser._get_call_name``: an attribute callee (e.g. ``foo.species(...)``) would
    never report as ``'species'`` on Arkane DSL input, since Arkane's own loader never
    constructs one -- so builtins with no attribute-access construction are simply not named,
    rather than guessed at.

    Args:
        call (ast.Call): The call node.

    Returns:
        str | None: The callee name, or ``None`` if the callee is not a bare name.
    """
    if isinstance(call.func, ast.Name):
        return call.func.id
    return None


def _call_keywords(call: ast.Call) -> dict:
    """
    Map a call's keyword arguments to their (still-AST) value nodes.

    A self-contained counterpart of ``t3.pdep.parser._call_keywords``, kept independent so this
    module remains a leaf module.

    Args:
        call (ast.Call): The call node.

    Raises:
        ValueError: If ``call`` has positional arguments or ``**kwargs`` unpacking.

    Returns:
        dict: Keyword name to AST value node.
    """
    if call.args:
        raise ValueError(f"Refusing to read keywords: {len(call.args)} positional argument(s) "
                          f"present; only keyword arguments are supported, since a positional "
                          f"argument's name is absent from the AST.")
    keywords = dict()
    for kw in call.keywords:
        if kw.arg is None:
            raise ValueError("Refusing to read keywords: '**kwargs' unpacking is present; the "
                              "unpacked names cannot be determined statically, so the keywords "
                              "actually supplied cannot be known.")
        keywords[kw.arg] = kw.value
    return keywords


def _literal_or_none(node):
    """
    Evaluate an AST node as a Python literal, or ``None`` if it cannot be.

    Args:
        node: An AST node (or ``None``).

    Returns:
        The literal value, or ``None`` if ``node`` is ``None`` or is not literal-evaluable
        (e.g. it is itself a ``Call``, such as ``NASA(...)``).
    """
    if node is None:
        return None
    try:
        return ast.literal_eval(node)
    except (ValueError, TypeError):
        return None


@dataclass(frozen=True)
class NetworkThermoCeiling:
    """
    The result of ``network_thermo_t_max``: the computed ceiling, plus which species (if any)
    could not be read and so did not contribute to it.

    Attributes:
        t_max (float | None): The minimum, over every ``species(...)`` call whose outer NASA
            ``Tmax`` could be positively determined in Kelvin, of that ``Tmax``; or ``None`` if
            no species contributed one at all (including when the file has no ``species(...)``
            calls).
        skipped (tuple[str, ...]): One human-readable entry per ``species(...)`` call that did
            NOT contribute to ``t_max``, naming the species (by its ``label=`` keyword when
            literal-evaluable, else a positional description) and why it was skipped. Empty
            when every species in the file contributed. A non-empty ``skipped`` means ``t_max``
            is only a ceiling among the species that COULD be read -- the true network-wide
            ceiling may be lower still, if a skipped species' real (unreadable) thermo ceiling
            is tighter than every readable one.
    """
    t_max: float | None
    skipped: tuple[str, ...]


def format_skipped_species(skipped: tuple, limit: int = 5) -> str:
    """
    Format a ``NetworkThermoCeiling.skipped`` tuple for display in a log message, capping the
    number of entries actually listed so a network with many skipped species does not flood the
    log.

    Args:
        skipped (tuple): The skip-reason entries to format (as returned by
                         ``network_thermo_t_max``'s ``.skipped``).
        limit (int, optional): The maximum number of entries to list verbatim before summarizing
                               the remainder as "and N more". Default: ``5``.

    Returns:
        str: A comma-joined, optionally-truncated rendering of ``skipped``. Empty string if
            ``skipped`` is empty.
    """
    if not skipped:
        return ''
    shown = list(skipped[:limit])
    remainder = len(skipped) - len(shown)
    if remainder > 0:
        shown.append(f'and {remainder} more')
    return '; '.join(shown)


def network_thermo_t_max(text: str) -> NetworkThermoCeiling:
    """
    Compute the network's NASA thermo Tmax ceiling among the species whose outer NASA ``Tmax``
    could be read, and report which species (if any) that ceiling does NOT account for.

    Arkane refuses to evaluate ANY NASA polynomial outside its own fitted range: a
    ``pressureDependence(...)`` block RMG wrote legitimately (RMG solves k(T,P) up to a pdep
    ``Tmax`` that exceeds a species' own thermo ceiling -- RMG tolerates the extrapolation at
    generation time) makes standalone Arkane fail outright with "No valid NASA polynomial at
    temperature ... K." when building statmech for that species. This function computes the
    TIGHTEST such ceiling -- the MINIMUM over species -- since one low-ceiling species makes the
    whole network's Arkane solve fail no matter how high another species' thermo is valid to.

    Only the OUTER ``Tmax=`` keyword on a species' ``NASA(...)`` call is read; a per-segment
    ``NASAPolynomial(..., Tmax=...)`` value nested inside ``polynomials=[...]`` is never
    consulted, since the outer value is the overall validity ceiling Arkane actually checks a
    requested temperature against (not the individual segments).

    IMPORTANT -- what this function does NOT guarantee: it does not promise a ceiling at which
    EVERY species has valid thermo. A species is silently excluded from the computed minimum,
    rather than forcing a spurious zero or a raise, whenever its outer thermo ceiling cannot be
    positively determined: its ``species(...)`` call has no ``thermo=`` keyword at all, its
    ``thermo=`` is not a ``NASA(...)`` call (e.g. ``ThermoData(...)``, which carries no
    polynomial-fit ceiling to compare), its outer ``Tmax`` is not literal-evaluable, or that
    ``Tmax``'s unit is not ``'K'`` (comparing a differently-unitted bound as if it were Kelvin
    would silently mis-clamp). This deliberately differs from ``t3.pdep.parser``'s fail-closed
    posture (e.g. ``_literal_or_raise``): an unreadable per-species thermo entry here is not
    evidence the FILE is malformed, only that one species cannot contribute a bound, and a
    caller computing a clamp ceiling would rather degrade to "no constraint from this species"
    than refuse the whole file over one exotic thermo entry. Bath-gas species are not
    special-cased; Arkane never builds thermo-dependent statmech for them, so they never need to
    contribute a bound in the first place.

    Consequently, when ``skipped`` in the returned ``NetworkThermoCeiling`` is non-empty, the
    returned ``t_max`` is only a ceiling among the species that COULD be read: the network's true
    ceiling may be lower still if a skipped species' real (unreadable) thermo is tighter. Callers
    that use ``t_max`` to clamp a solve should surface ``skipped`` (e.g. as a warning naming the
    skipped species) so this under-protection is visible rather than silent.

    Args:
        text (str): The RMG pdep network file content.

    Raises:
        NetworkTextUnparseable: If ``text`` cannot be parsed as Python.

    Returns:
        NetworkThermoCeiling: The computed ceiling (``None`` if no species contributed one) and
            the skip diagnostics described above.

    See also:
        t3.pdep.parser.parse_pdep_network_text: RMG-generated pdep network files can trigger a
            SyntaxWarning (e.g. an invalid escape sequence in a comment) from Python's own
            ast.parse; that warning is likewise suppressed here.
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', SyntaxWarning)
            tree = ast.parse(text)
    except SyntaxError as e:
        raise NetworkTextUnparseable(
            f"Could not parse pdep network text to compute network_thermo_t_max: {e}")

    t_max_values = []
    skipped = []
    for node in tree.body:
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
            continue
        call = node.value
        if _get_call_name(call) != 'species':
            continue
        kwargs = _call_keywords(call)
        label = _literal_or_none(kwargs.get('label'))
        name = repr(label) if isinstance(label, str) else '<species with unreadable label>'

        thermo_node = kwargs.get('thermo')
        if thermo_node is None:
            skipped.append(f"{name} (no thermo= keyword)")
            continue
        if not isinstance(thermo_node, ast.Call) or _get_call_name(thermo_node) != 'NASA':
            skipped.append(f"{name} (thermo is not a NASA(...) call)")
            continue
        thermo_kwargs = _call_keywords(thermo_node)
        t_max_value = _literal_or_none(thermo_kwargs.get('Tmax'))
        if not isinstance(t_max_value, tuple) or len(t_max_value) != 2:
            skipped.append(f"{name} (NASA Tmax is missing or not a literal (magnitude, unit) tuple)")
            continue
        magnitude, unit = t_max_value
        if not isinstance(magnitude, (int, float)) or isinstance(magnitude, bool):
            skipped.append(f"{name} (NASA Tmax magnitude is not numeric)")
            continue
        if unit != 'K':
            skipped.append(f"{name} (NASA Tmax unit is {unit!r}, not 'K')")
            continue
        t_max_values.append(float(magnitude))

    return NetworkThermoCeiling(
        t_max=min(t_max_values) if t_max_values else None,
        skipped=tuple(skipped),
    )
