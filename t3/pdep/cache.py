"""
t3 pdep cache module

Validity metadata for cached Arkane PDep sensitivity output.

Re-using a sensitivity YAML that is already on disk saves a very expensive Arkane run, but doing
so naively is unsafe. Nothing in Arkane's own output identifies which Arkane produced it or which
network it came from, so this module writes a small T3-owned sidecar next to every generated
sensitivity YAML, recording what T3 needs in order to trust it later, and refuses any cache that
lacks one.

Cache *validity* (right hashes, right method, right contract version, parseable) is a different
question from data *usefulness* (does the SA output actually carry transition-state signal), and
this module answers only the first one. A sensitivity analysis whose transition-state perturbation
never reached the rate expression writes structurally dead transition-state rows that look exactly
like "nothing here is sensitive" -- but that is a judgment about whether a network qualifies for
QM, made per-reaction-key and at the correct granularity by
``t3.pdep.selector.select_from_sa_dict``, not a reason to distrust the cache itself. This module
still records ``max_abs_ts_coefficient`` in the sidecar, but purely as provenance for a human
reading the file, not as a gate here.
"""

import math
import os
import re
import subprocess


from arc.common import save_yaml_file

# Re-exported, not defined here: ``t3.pdep.cache.hash_file`` is named as the mandated hashing
# primitive throughout join.py, capture.py, and main.py, so the import path is kept while the
# format itself lives in one place alongside its bytes-input sibling ``hash_bytes``.
from t3.pdep.hashing import hash_file  # noqa: F401
from t3.pdep.yaml_safe import read_sa_yaml_file

from t3.pdep.selector import (CACHE_STATUS_CACHED_REJECTED,
                              CACHE_STATUS_CACHED_VALID,
                              DEFAULT_MIN_DELTA_LN_K,
                              E0_PERTURBATION_J_PER_MOL,
                              STRUCTURES_KEY,
                              TS_ENTRY_PREFIX,
                              coefficient_floor,
                              )

SA_CACHE_METADATA_FILE_NAME = 't3_sa_cache.yml'

# The label Arkane prints immediately before the RMG-Py commit hash in its own log.
ARKANE_LOG_RMG_PY_COMMIT_LABEL = 'current git HEAD for RMG-Py'

# Whether a T3-generated sensitivity YAML on disk may still be trusted -- this is the ONLY job this
# constant has. It is deliberately separate from ``t3.pdep.selector.SELECTION_SCHEMA_VERSION`` (the
# on-disk shape of a saved *selection* record) and ``SELECTION_ALGORITHM_VERSION`` (the selection
# *decision logic*): a cache can go stale for reasons that have nothing to do with either of those
# (e.g. this module starts recording a new provenance field the write side must supply), and bumping
# either of those two must not force every cached SA YAML to be needlessly regenerated. Bump this
# when what this module needs from the sidecar to trust a cached file changes. There is no fallback
# to a predecessor key: no ``t3_sa_cache.yml`` sidecar has ever existed on disk anywhere under the
# old combined marker, so a fallback would be untestable dead code, not a real migration path.
SA_CACHE_CONTRACT_VERSION = 1

# A mismatch of more than this (in J/mol) between the recorded and the assumed perturbation
# invalidates the cache, since the absolute floor is derived from the perturbation size.
PERTURBATION_TOLERANCE_J_PER_MOL = 1.0


def sa_cache_metadata_path(sa_path: str) -> str:
    """
    The path of the T3 sidecar belonging to a sensitivity YAML.

    Args:
        sa_path (str): The path to the Arkane ``sa_coefficients.yml``.

    Returns:
        str: The path to the sidecar metadata file.
    """
    return os.path.join(os.path.dirname(sa_path), SA_CACHE_METADATA_FILE_NAME)




def max_abs_ts_coefficient(sa_dict: dict) -> float | None:
    """
    Scan a loaded Arkane sensitivity dictionary for the largest absolute finite TS coefficient.

    The T3 sidecar records this value as provenance only, for a human reading the sidecar; it is
    not used by ``validate_sa_cache`` as a cache-validity gate. Whether an SA output whose
    transition-state rows are all structural zeros (i.e. this scan finds nothing above the floor)
    carries usable criterion-(b) signal is decided by ``t3.pdep.selector.select_from_sa_dict``.

    Args:
        sa_dict (dict): A loaded Arkane PDep sensitivity dictionary. Keys are network reaction
            strings (plus a ``'structures'`` entry); values map a condition to a mapping of entry
            label to sensitivity coefficient.

    Returns:
        Optional[float]: The largest ``abs()`` of any finite ``'(TS) '``-prefixed coefficient found
            anywhere in ``sa_dict``, or ``None`` if ``sa_dict`` is not a dict or contains no such
            finite coefficient.
    """
    if not isinstance(sa_dict, dict):
        return None
    largest = None
    for key, entries_by_condition in sa_dict.items():
        if key == STRUCTURES_KEY or not isinstance(entries_by_condition, dict):
            continue
        for entries in entries_by_condition.values():
            if not isinstance(entries, dict):
                continue
            for entry, coefficient in entries.items():
                if not isinstance(entry, str) or not entry.startswith(TS_ENTRY_PREFIX):
                    continue
                if not isinstance(coefficient, (int, float)) or isinstance(coefficient, bool) \
                        or not math.isfinite(coefficient):
                    continue
                value = abs(float(coefficient))
                if largest is None or value > largest:
                    largest = value
    return largest


def get_git_commit(path: str | None) -> str | None:
    """
    Best-effort lookup of the git commit of a checkout, recorded for provenance only.

    Args:
        path (str, optional): A path inside the git repository, e.g. an RMG-Py checkout.

    Returns:
        Optional[str]: The commit hash, or ``None`` if it could not be determined.
    """
    if not path or not os.path.isdir(path):
        return None
    try:
        result = subprocess.run(['git', '-C', path, 'rev-parse', 'HEAD'],
                                capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout.strip() or None if result.returncode == 0 else None


def read_arkane_log_rmg_py_commit(arkane_log_path: str | None) -> str | None:
    """
    Best-effort lookup of the RMG-Py commit that produced an Arkane run, from Arkane's own log.

    T3 cannot identify that commit by introspecting its own interpreter: ARC runs Arkane in a
    subprocess, in a different conda environment (``micromamba run -n rmg_env python -m arkane``),
    so nothing of the RMG-Py that did the work is ever loaded here. Arkane, however, records its
    own provenance in the log it writes into the output directory T3 already owns, which makes the
    log the only first-hand witness available.

    The commit is not on the label line but on the line after it, and the line after THAT is a
    commit date, so only the first following non-empty line is considered and it is returned only
    if it actually looks like a hash. Recording an arbitrary log line would be worse than recording
    nothing, since it would read as real provenance to whoever audits the sidecar later.

    Args:
        arkane_log_path (str, optional): The path to an Arkane ``arkane.log``.

    Returns:
        Optional[str]: The RMG-Py commit hash, or ``None`` if it could not be determined (no log,
            unreadable log, no RMG-Py stanza, or a value that is not a hash).
    """
    if not arkane_log_path or not os.path.isfile(arkane_log_path):
        return None
    try:
        with open(arkane_log_path, 'r', errors='replace') as f:
            lines = f.readlines()
    except OSError:
        return None
    for index, line in enumerate(lines):
        if ARKANE_LOG_RMG_PY_COMMIT_LABEL not in line:
            continue
        for candidate in lines[index + 1:]:
            candidate = candidate.strip()
            if not candidate:
                continue
            return candidate if re.fullmatch(r'[0-9a-fA-F]{7,40}', candidate) else None
    return None


def write_sa_cache_metadata(sa_path: str,
                            network_path: str,
                            network_id: str,
                            method: str,
                            perturbation: float = E0_PERTURBATION_J_PER_MOL,
                            rmg_py_path: str | None = None,
                            rmg_py_commit: str | None = None,
                            t_grid_clamp: dict | None = None,
                            ) -> str:
    """
    Write the T3 sidecar next to a freshly generated sensitivity YAML.

    Args:
        sa_path (str): The path to the generated ``sa_coefficients.yml``.
        network_path (str): The path to the RMG network file the analysis was run on.
        network_id (str): The network file stem, e.g. ``'network4_2'``.
        method (str): The master-equation method used, e.g. ``'MSC'``.
        perturbation (float, optional): The E0 perturbation Arkane applied, in J/mol.
        rmg_py_path (str, optional): The RMG-Py checkout used, recorded for provenance.
        rmg_py_commit (str, optional): The RMG-Py commit that actually ran, e.g. from
            ``read_arkane_log_rmg_py_commit``. Recorded verbatim; takes precedence over deriving
            a commit from ``rmg_py_path``.
        t_grid_clamp (dict, optional): ``t3.utils.network_thermo.TGridClampRecord.as_dict()`` for
            the Arkane input file this SA was run against, i.e. whether the T grid Arkane actually
            solved over was clamped down from what was originally requested. Written into the
            sidecar ONLY when supplied -- so a sidecar written by a caller that never passed this
            (an old sidecar, or SA data produced outside this code path) omits the key entirely
            rather than recording a value that would misrepresent "unknown" as "not clamped". See
            ``TGridClampRecord``'s docstring for why those two states must never be conflated.

    Returns:
        str: The path of the sidecar that was written.
    """
    try:
        sa_dict = read_sa_yaml_file(sa_path)
    except Exception:
        sa_dict = None
    metadata = {'sa_cache_contract_version': SA_CACHE_CONTRACT_VERSION,
                'network_id': network_id,
                'method': method,
                'perturbation': perturbation,
                'perturbation_units': 'J/mol',
                'network_file_hash': hash_file(network_path),
                'network_file_name': os.path.basename(network_path),
                # Binds this sidecar to the exact SA YAML content it vouches for, so the YAML
                # cannot be replaced or hand-edited afterwards and still validate.
                'sa_file_hash': hash_file(sa_path),
                # The largest abs finite TS coefficient anywhere in the SA output. Recorded as
                # provenance only, useful for a human inspecting the sidecar (e.g. to see whether
                # an un-patched Arkane failed to propagate a TS perturbation into an ILT-based
                # rate) -- validate_sa_cache() does NOT gate cache validity on this value; whether
                # the data is USEFUL for criterion (b) is decided by
                # t3.pdep.selector.select_from_sa_dict, not here.
                'max_abs_ts_coefficient': max_abs_ts_coefficient(sa_dict),
                'rmg_py_path': rmg_py_path,
                # A commit read out of Arkane's log is first-hand evidence of what actually ran and
                # is recorded as given; deriving one from a local checkout is the fallback for a
                # caller that has the checkout rather than the log.
                'rmg_py_commit': rmg_py_commit if rmg_py_commit is not None else get_git_commit(rmg_py_path),
                }
    if t_grid_clamp is not None:
        # Not gated on by validate_sa_cache (see SA_CACHE_CONTRACT_VERSION's comment and this
        # module's docstring): this is pure provenance, added without bumping the contract
        # version, since a sidecar written before this key existed is merely missing an optional
        # field, not misread by anything that consults sa_cache_contract_version.
        metadata['t_grid_clamp'] = t_grid_clamp
    metadata_path = sa_cache_metadata_path(sa_path)
    save_yaml_file(path=metadata_path, content=metadata)
    return metadata_path


def read_t_grid_clamp_record(sa_path: str) -> dict | None:
    """
    Best-effort read of the T-grid clamp provenance recorded alongside a sensitivity YAML.

    Deliberately as lenient as ``write_sa_cache_metadata`` is strict: unlike
    ``validate_sa_cache`` (which fails closed on anything it cannot fully trust), this function
    exists only to *disclose* provenance, never to gate anything. A missing sidecar, an
    unparseable sidecar, or a sidecar that predates this key all collapse to the same ``None``
    ("unknown provenance") -- exactly as absent as if this function had never been called -- so
    that a caller cannot mistake "I don't know" for "no clamp happened" (see
    ``TGridClampRecord``'s docstring), and so that unknown provenance can never, by itself, cause
    a refusal.

    Args:
        sa_path (str): The path to the ``sa_coefficients.yml`` whose sidecar should be consulted.

    Returns:
        Optional[dict]: The recorded ``TGridClampRecord.as_dict()`` mapping, or ``None`` if no
            sidecar exists, it could not be read, or it has no ``t_grid_clamp`` key.
    """
    metadata_path = sa_cache_metadata_path(sa_path)
    if not os.path.isfile(metadata_path):
        return None
    try:
        # The restricted loader, not arc's read_yaml_file (yaml.FullLoader): this sidecar sits
        # adjacent to a caller-supplied sa_path, and FullLoader constructs Python objects from tags
        # DURING the load -- before the except below could ever collapse the read to None.
        metadata = read_sa_yaml_file(metadata_path)
    except Exception:
        return None
    if not isinstance(metadata, dict):
        return None
    t_grid_clamp = metadata.get('t_grid_clamp')
    return t_grid_clamp if isinstance(t_grid_clamp, dict) else None


def validate_sa_cache(sa_path: str,
                      network_path: str,
                      perturbation: float = E0_PERTURBATION_J_PER_MOL,
                      min_delta_ln_k: float = DEFAULT_MIN_DELTA_LN_K,
                      method: str | None = None,
                      ) -> tuple:
    """
    Decide whether a sensitivity YAML already on disk may be reused.

    A cache is only accepted when T3 itself wrote the sidecar, the selection semantics have not
    changed since, the perturbation matches the one the gates assume, the network file has not
    changed, the SA YAML itself has not changed since the sidecar was written, and the ME method
    matches (when given). Anything else is rejected so the analysis is regenerated -- silently
    reusing a stale or foreign file would corrupt the decision with no visible symptom.

    This function answers only cache VALIDITY (is this the right, unmodified, T3-written file).
    Whether the SA data it validates is actually USEFUL -- i.e. whether its transition-state rows
    carry any signal above the absolute coefficient floor -- is a different question, answered per
    reaction-key by ``t3.pdep.selector.select_from_sa_dict``, not here. ``min_delta_ln_k`` and
    ``perturbation`` are still validated eagerly below (see the comment ahead of ``floor``) purely
    so a caller bug in either surfaces immediately as a ``ValueError``, regardless of which
    validity branch a given call happens to reach first.

    Args:
        sa_path (str): The path to the candidate ``sa_coefficients.yml``.
        network_path (str): The path to the RMG network file it should correspond to.
        perturbation (float, optional): The E0 perturbation the selector will assume, in J/mol.
        min_delta_ln_k (float, optional): The minimum meaningful ``ln(k)`` response; validated
            eagerly here (via ``t3.pdep.selector.coefficient_floor``) only to surface a bad value
            as a ``ValueError`` immediately -- the resulting floor is not otherwise used in this
            function.
        method (str, optional): The master-equation method this call expects the cache to have
            been generated with, e.g. ``'MSC'``. If given and it disagrees with the recorded
            ``method``, the cache is rejected.

    Returns:
        tuple: (status (str), warnings (list)). The status is ``'cached_valid'`` or ``'cached_rejected'``.
    """
    # Validate the thresholds before any early return: a bad perturbation/min_delta_ln_k is a
    # caller bug and should surface as a ValueError immediately, not be masked by whichever
    # early-return branch (missing file, missing sidecar, ...) a particular call happens to hit
    # first -- the floor derived below is only actually used once we reach the final gate, but its
    # inputs must be valid regardless of how far this call gets.
    floor = coefficient_floor(min_delta_ln_k=min_delta_ln_k, perturbation=perturbation)

    warnings_list = list()
    if not os.path.isfile(sa_path):
        return CACHE_STATUS_CACHED_REJECTED, [f'No cached sensitivity output at {sa_path}.']
    metadata_path = sa_cache_metadata_path(sa_path)
    if not os.path.isfile(metadata_path):
        return CACHE_STATUS_CACHED_REJECTED, [
            f'The cached sensitivity output at {sa_path} has no T3 metadata sidecar, so it cannot be trusted '
            f'(its transition state rows may predate the Arkane fix that makes them meaningful). Regenerating.']

    try:
        # The restricted loader, not arc's read_yaml_file (yaml.FullLoader): the sidecar sits
        # adjacent to a caller-supplied sa_path, and FullLoader constructs Python objects from tags
        # DURING the load -- before the except below could ever turn the read into a rejection.
        metadata = read_sa_yaml_file(metadata_path) or dict()
    except Exception as e:
        return CACHE_STATUS_CACHED_REJECTED, [f'Could not read the cache metadata at {metadata_path}: {e}. '
                                              f'Regenerating.']
    if not isinstance(metadata, dict):
        return CACHE_STATUS_CACHED_REJECTED, [f'The cache metadata at {metadata_path} is malformed. Regenerating.']

    # Fail-closed on an absent key too, not just a mismatched one: a sidecar written before this
    # key existed (or hand-edited to drop it) provides no evidence at all that this module's
    # current usability contract was satisfied, and ``None != SA_CACHE_CONTRACT_VERSION`` already
    # takes this branch, so there is no separate case to add here.
    recorded_version = metadata.get('sa_cache_contract_version')
    if recorded_version != SA_CACHE_CONTRACT_VERSION:
        warnings_list.append(f'The cached sensitivity output at {sa_path} was written under SA-cache '
                             f'contract version {recorded_version}, but this is version '
                             f'{SA_CACHE_CONTRACT_VERSION}. Regenerating.')
        return CACHE_STATUS_CACHED_REJECTED, warnings_list

    recorded_perturbation = metadata.get('perturbation')
    if not isinstance(recorded_perturbation, (int, float)) \
            or abs(float(recorded_perturbation) - perturbation) > PERTURBATION_TOLERANCE_J_PER_MOL:
        warnings_list.append(f'The cached sensitivity output at {sa_path} used a perturbation of '
                             f'{recorded_perturbation} J/mol, but {perturbation} J/mol is assumed. Regenerating.')
        return CACHE_STATUS_CACHED_REJECTED, warnings_list

    if not os.path.isfile(network_path):
        warnings_list.append(f'The network file {network_path} is missing, so the cached sensitivity output at '
                             f'{sa_path} cannot be validated against it. Regenerating.')
        return CACHE_STATUS_CACHED_REJECTED, warnings_list
    if metadata.get('network_file_hash') != hash_file(network_path):
        warnings_list.append(f'The network file {network_path} changed since the sensitivity output at {sa_path} '
                             f'was generated. Regenerating.')
        return CACHE_STATUS_CACHED_REJECTED, warnings_list

    if metadata.get('sa_file_hash') != hash_file(sa_path):
        warnings_list.append(f'The sensitivity output at {sa_path} changed since its T3 sidecar was written '
                             f'(content hash mismatch). Regenerating.')
        return CACHE_STATUS_CACHED_REJECTED, warnings_list

    recorded_method = metadata.get('method')
    if method is not None and recorded_method != method:
        warnings_list.append(f'The cached sensitivity output at {sa_path} was generated with method '
                             f'{recorded_method!r}, but method {method!r} was requested. Regenerating.')
        return CACHE_STATUS_CACHED_REJECTED, warnings_list

    # No max_abs_ts_coefficient-vs-floor check here: cache validity and data usefulness are
    # different questions (see module docstring). ``floor`` above exists purely to validate
    # ``min_delta_ln_k``/``perturbation`` eagerly; it is not otherwise consulted in this function.
    return CACHE_STATUS_CACHED_VALID, warnings_list
