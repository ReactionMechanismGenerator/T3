"""
t3 pdep cache module

Validity metadata for cached Arkane PDep sensitivity output.

Re-using a sensitivity YAML that is already on disk saves a very expensive Arkane run, but doing
so naively is unsafe. Nothing in Arkane's own output identifies which Arkane produced it or which
network it came from, and a sensitivity analysis whose transition-state perturbation never reached
the rate expression writes structurally dead transition-state rows that look exactly like "nothing
here is sensitive". Silently trusting such a file makes a network fail criterion (b) with no
warning -- the network is quietly never selected for QM.

This module therefore writes a small T3-owned sidecar next to every generated sensitivity YAML,
recording what T3 needs in order to trust it later, and refuses any cache that lacks one.
"""

import math
import os
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
                              SELECTOR_VERSION,
                              STRUCTURES_KEY,
                              TS_ENTRY_PREFIX,
                              coefficient_floor,
                              )

SA_CACHE_METADATA_FILE_NAME = 't3_sa_cache.yml'

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

    This is what the T3 sidecar binds itself to: an SA output whose transition-state rows are all
    structural zeros (i.e. this scan finds nothing above the floor) carries no criterion-(b) signal
    at all, no matter what the caller's decision looked like at write time.

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


def write_sa_cache_metadata(sa_path: str,
                            network_path: str,
                            network_id: str,
                            method: str,
                            perturbation: float = E0_PERTURBATION_J_PER_MOL,
                            rmg_py_path: str | None = None,
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

    Returns:
        str: The path of the sidecar that was written.
    """
    try:
        sa_dict = read_sa_yaml_file(sa_path)
    except Exception:
        sa_dict = None
    metadata = {'selector_version': SELECTOR_VERSION,
                'network_id': network_id,
                'method': method,
                'perturbation': perturbation,
                'perturbation_units': 'J/mol',
                'network_file_hash': hash_file(network_path),
                'network_file_name': os.path.basename(network_path),
                # Binds this sidecar to the exact SA YAML content it vouches for, so the YAML
                # cannot be replaced or hand-edited afterwards and still validate.
                'sa_file_hash': hash_file(sa_path),
                # The largest abs finite TS coefficient anywhere in the SA output. A cache whose
                # TS rows are all structural zeros (e.g. from an un-patched Arkane that never
                # propagates a TS perturbation into an ILT-based rate) must never be trusted as
                # "valid" -- see validate_sa_cache().
                'max_abs_ts_coefficient': max_abs_ts_coefficient(sa_dict),
                'rmg_py_path': rmg_py_path,
                'rmg_py_commit': get_git_commit(rmg_py_path),
                }
    metadata_path = sa_cache_metadata_path(sa_path)
    save_yaml_file(path=metadata_path, content=metadata)
    return metadata_path


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
    changed, the SA YAML itself has not changed since the sidecar was written, the ME method
    matches (when given), and the SA data itself still carries transition-state signal. Anything
    else is rejected so the analysis is regenerated -- silently reusing a stale, foreign, or
    signal-dead file would suppress the network from QM selection with no visible symptom.

    Args:
        sa_path (str): The path to the candidate ``sa_coefficients.yml``.
        network_path (str): The path to the RMG network file it should correspond to.
        perturbation (float, optional): The E0 perturbation the selector will assume, in J/mol.
        min_delta_ln_k (float, optional): The minimum meaningful ``ln(k)`` response used to derive
            the absolute coefficient floor that the cached data's ``max_abs_ts_coefficient`` must
            clear; see ``t3.pdep.selector.coefficient_floor``.
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

    recorded_version = metadata.get('selector_version')
    if recorded_version != SELECTOR_VERSION:
        warnings_list.append(f'The cached sensitivity output at {sa_path} was written by selector version '
                             f'{recorded_version}, but this is version {SELECTOR_VERSION}. Regenerating.')
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

    recorded_max_abs_ts_coefficient = metadata.get('max_abs_ts_coefficient')
    if not isinstance(recorded_max_abs_ts_coefficient, (int, float)) \
            or isinstance(recorded_max_abs_ts_coefficient, bool) \
            or recorded_max_abs_ts_coefficient < floor:
        warnings_list.append(
            f'The cached sensitivity output at {sa_path} has no transition-state coefficient at or above the '
            f'absolute floor ({floor:.3e} mol/J, recorded max was {recorded_max_abs_ts_coefficient!r}). Its '
            f'transition state rows carry no signal -- this usually means the SA was produced by an Arkane that '
            f'does not propagate a TS perturbation into an ILT-based rate. Regenerating.')
        return CACHE_STATUS_CACHED_REJECTED, warnings_list

    return CACHE_STATUS_CACHED_VALID, warnings_list
