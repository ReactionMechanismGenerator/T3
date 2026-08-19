"""
t3 pdep pes_sa module

The standalone PES loop's own source of REAL per-transition-state sensitivity evidence: a
master-equation (ME) sensitivity analysis run with Arkane on the round's freshly explored network.

Why this exists: ``t3.pdep.capture.capture_ts_artifacts`` refuses (by deliberate, fail-closed
design) to capture any transition-state QM artifact that does not carry the ``coefficient``/
``delta_ln_k`` sensitivity evidence that justified selecting it -- see ``verify_capture``. Inside
a T3 iteration that evidence comes from a dedicated Arkane ME SA run
(``t3.utils.writer.write_pdep_network_file`` -> ``t3.runners.rmg_runner.run_arkane_job`` ->
``sensitivity/sa_coefficients.yml`` -> ``t3.pdep.selector``), frozen onto each ``TSJoinRecord`` at
queue time (``t3.main.T3.queue_pdep_transition_states``) and read back off the durable join
records at capture time (``t3.main.T3._capture_pdep_ts_artifacts``). The standalone loop
(``t3.pdep.pes_loop``) had NO such source, so every capture it attempted was -- correctly --
refused. This module runs the same Arkane ME SA, through the same production seams
(``write_arkane_network_input_file`` with ``sensitivity=True``, ``run_arkane_job``,
``read_sa_yaml_file``), against the round's own explored network.

What the evidence means here: the loop has no single caller-declared observable reaction the way a
T3 iteration does (there, criterion (a) names the network reaction an observable is sensitive to,
and the selector examines exactly that direction key). The standalone loop's k(T,P) of interest is
the network itself, so each transition state's evidence is its strongest measured leverage --
the maximum-|coefficient| dln(k)/dE0 row for that transition state across EVERY network reaction
direction and every condition of the SA output. That is a real, measured Arkane derivative, never
a default or a sentinel: a transition state for which the SA reports no finite row gets NO
evidence here, and the loop refuses to queue it rather than inventing a number
(``t3.pdep.pes_rounds.attach_sensitivity_evidence``).

This module never imports ``rmgpy`` or ``arkane``; Arkane runs behind ``run_arkane_job`` (ARC's
own entry point), exactly as T3's in-iteration ME SA does.
"""

import logging
import os

from t3.pdep.selector import E0_PERTURBATION_J_PER_MOL, STRUCTURES_KEY, TS_ENTRY_PREFIX, _is_finite
from t3.pdep.yaml_safe import read_sa_yaml_file
from t3.runners.rmg_runner import run_arkane_job
from t3.utils.writer import write_arkane_network_input_file

_logger = logging.getLogger(__name__)

# Where Arkane writes the ME SA coefficients, relative to the SA run's output directory -- the
# same artifact ``run_arkane_job`` itself requires (its ``required_artifact`` default).
SA_COEFFICIENTS_RELPATH = os.path.join('sensitivity', 'sa_coefficients.yml')


def ts_sensitivity_evidence(sa_dict: dict,
                            perturbation: float = E0_PERTURBATION_J_PER_MOL,
                            ) -> dict[str, tuple[float, float]]:
    """
    Extract per-transition-state sensitivity evidence from a loaded Arkane ME SA dictionary.

    For each transition state, the evidence is the maximum-|coefficient| row across EVERY network
    reaction direction key and every (T, P) condition -- the strongest measured leverage this
    transition state's E0 has on any k(T,P) of the network (see the module docstring for why the
    aggregation is over all directions rather than one caller-declared observable direction).

    Malformed pieces (a non-dict direction entry, a non-dict condition entry, a non-finite
    coefficient, the ``structures`` mapping) are skipped, never guessed at: a transition state
    whose every row is unusable simply gets no evidence, which downstream refuses loudly.

    Args:
        sa_dict (dict): A loaded Arkane PDep sensitivity dictionary (see
            ``t3.pdep.selector._select_from_sa_dict`` for the shape). Keys are network reaction
            strings (plus a ``'structures'`` entry); values map a (T, 'K', P, 'bar') condition to
            a mapping of entry label to sensitivity coefficient.
        perturbation (float, optional): The E0 perturbation Arkane applied, in J/mol, used to
            compute each ``delta_ln_k`` exactly as ``t3.pdep.selector._evidence_for_ts`` does.

    Returns:
        dict[str, tuple[float, float]]: Network-local transition state label ->
        ``(coefficient, delta_ln_k)``.

    Raises:
        ValueError: If ``sa_dict`` is not a dict at all -- a malformed payload is a fact about the
                   SA artifact, not "a network with no sensitive transition states", and must not
                   be silently read as the latter.
    """
    if not isinstance(sa_dict, dict):
        raise ValueError(f'The sensitivity payload is malformed: expected a dict, got '
                         f'{type(sa_dict).__name__}.')
    evidence = dict()
    for direction_key, conditions in sa_dict.items():
        if direction_key == STRUCTURES_KEY or not isinstance(conditions, dict):
            continue
        for _condition, entries in conditions.items():
            if not isinstance(entries, dict):
                continue
            for entry, coefficient in entries.items():
                if not isinstance(entry, str) or not entry.startswith(TS_ENTRY_PREFIX):
                    continue
                if not _is_finite(coefficient):
                    continue
                ts_label = entry[len(TS_ENTRY_PREFIX):]
                coefficient = float(coefficient)
                current = evidence.get(ts_label)
                if current is None or abs(coefficient) > abs(current[0]):
                    evidence[ts_label] = (coefficient, abs(coefficient) * perturbation)
    return evidence


def run_round_me_sensitivity(network_path: str,
                             sa_dir: str,
                             method: str,
                             timeout: float | None = None,
                             logger=None,
                             ) -> dict[str, tuple[float, float]]:
    """
    Run one round's Arkane ME sensitivity analysis and return the per-TS evidence it measured.

    Mirrors T3's in-iteration pipeline exactly, through the same production seams: render an
    Arkane input for ``network_path`` with a ``sensitivity_conditions`` directive spanning the
    network's own T/P extrema (``write_arkane_network_input_file(..., sensitivity=True)``), run
    Arkane on it (``run_arkane_job``, which itself requires ``sensitivity/sa_coefficients.yml``
    to have been (re)written for the job to count as succeeded), then load and reduce that output
    (``read_sa_yaml_file`` -> ``ts_sensitivity_evidence``).

    Args:
        network_path (str): The network file to analyze -- the round's freshly explored network.
        sa_dir (str): The SA run's own output directory (``t3.pdep.pes_rounds.RoundPaths.sa``).
        method (str): The master-equation method, 'CSE', 'MSC' or 'RS' (``config.pes.method``).
        timeout (float, optional): Wall-clock deadline for the Arkane process, in seconds, passed
            through to ``run_arkane_job``. ``None`` means no deadline.
        logger: A T3 ``Logger``, or ``None``. Passed through to ``run_arkane_job``.

    Returns:
        dict[str, tuple[float, float]]: Network-local transition state label ->
        ``(coefficient, delta_ln_k)``. May name transition states the caller is not interested in
        (Arkane reports every TS row); may also lack a transition state whose every row was
        non-finite -- the caller decides what a missing entry means.

    Raises:
        ValueError: If the Arkane SA job failed, timed out, or produced a malformed payload --
                   the round then has no real sensitivity evidence, and the caller must refuse to
                   queue QM rather than invent any (see ``t3.pdep.pes_loop``).
        OSError: If the SA input could not be written, or the SA output could not be read.
    """
    os.makedirs(sa_dir, exist_ok=True)
    input_path = os.path.join(sa_dir, 'input.py')
    write_arkane_network_input_file(source_path=network_path, dest_path=input_path,
                                    method=method, sensitivity=True)
    job_result = run_arkane_job(input_file=input_path, output_directory=sa_dir,
                                logger=logger, timeout=timeout)
    if not job_result:
        raise ValueError(f'The master-equation sensitivity analysis for {network_path!r} failed: '
                         f'{job_result.reason}')
    sa_path = os.path.join(sa_dir, SA_COEFFICIENTS_RELPATH)
    evidence = ts_sensitivity_evidence(read_sa_yaml_file(sa_path))
    _logger.info(f'PES loop: the master-equation sensitivity analysis for {network_path!r} '
                 f'measured evidence for {len(evidence)} transition state(s): '
                 f'{sorted(evidence)}.')
    return evidence
