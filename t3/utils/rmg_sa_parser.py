"""
Utility for parsing RMG sensitivity analysis CSV files into T3's
standardized per-condition sa_dict format.

RMG SA CSV files have the naming convention::

    sensitivity_{reactor_index}_SPC_{species_index}.csv

Column headers follow the format::

    Time(s), dln[Obs(ID)]/dln[k#], ..., dln[Obs(ID)]/dG[Spc(ID)], ...

This module converts these CSV files into the sa_dict structure used by
all T3 simulate adapters, without requiring rmgpy.
"""

import csv
import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

from t3.common import get_chem_to_rmg_rxn_index_map


# Regex for parsing SA column headers.
# Matches: dln[ObsLabel(ID)]/dln[k123]  or  dln[ObsLabel(ID)]/dG[SpcLabel(ID)]
_HEADER_RE = re.compile(
    r'^dln\[(?P<obs>[^\]]+)\]'    # numerator: dln[Obs(ID)]
    r'/d(?P<kind>ln|G)'           # denominator type: dln or dG
    r'\[(?P<param>[^\]]+)\]$'     # parameter: k123 or SpcLabel(ID)
)


def _strip_rmg_index(label: str) -> str:
    """Strip trailing RMG index ``(ID)`` from a species label.

    Examples::

        'OH(4)'   -> 'OH'
        'H2O'     -> 'H2O'
        'O(T)(5)' -> 'O(T)'
    """
    if '(' in label and label.endswith(')'):
        return label.rsplit('(', 1)[0]
    return label


def parse_rmg_sa_csv(csv_path: str) -> Tuple[List[float], Dict[str, Dict[int, List[float]]], Dict[str, Dict[str, List[float]]]]:
    """Parse a single RMG SA CSV file.

    Args:
        csv_path: Path to a ``sensitivity_*.csv`` file.

    Returns:
        A 3-tuple ``(time, kinetics, thermo)`` where:
        - *time* is a list of floats (time points in seconds).
        - *kinetics* maps ``observable_label -> {reaction_index: [values]}``.
        - *thermo* maps ``observable_label -> {species_label: [values]}``.
    """
    time_data: List[float] = []
    kinetics: Dict[str, Dict[int, List[float]]] = {}
    thermo: Dict[str, Dict[str, List[float]]] = {}

    with open(csv_path, 'r', newline='') as fh:
        reader = csv.reader(fh)
        headers = next(reader)

        # Identify the time column and SA columns
        time_col: Optional[int] = None
        sa_columns: List[Tuple[int, str, str, str]] = []  # (col_idx, obs_label, kind, param_raw)

        for col_idx, header in enumerate(headers):
            header = header.strip()
            if 'time' in header.lower():
                time_col = col_idx
                continue
            m = _HEADER_RE.match(header)
            if m:
                obs_label = _strip_rmg_index(m.group('obs'))
                kind = 'kinetics' if m.group('kind') == 'ln' else 'thermo'
                param_raw = m.group('param')
                sa_columns.append((col_idx, obs_label, kind, param_raw))

        # Read rows
        for row in reader:
            if not row or not row[0].strip():
                continue
            if time_col is not None:
                time_data.append(float(row[time_col]))
            for col_idx, obs_label, kind, param_raw in sa_columns:
                value = float(row[col_idx])
                if kind == 'kinetics':
                    # param_raw is e.g. 'k10' -> extract integer 10
                    rxn_idx = int(param_raw.lstrip('k'))
                    kinetics.setdefault(obs_label, {}).setdefault(rxn_idx, []).append(value)
                else:
                    spc_label = _strip_rmg_index(param_raw)
                    thermo.setdefault(obs_label, {}).setdefault(spc_label, []).append(value)

    return time_data, kinetics, thermo


def parse_rmg_sa_csvs(solver_dir: str) -> dict:
    """Parse all RMG SA CSV files in a solver directory.

    Args:
        solver_dir: Path to the ``solver/`` directory containing
                    ``sensitivity_*.csv`` files.

    Returns:
        A raw dict with keys ``'time'``, ``'kinetics'``, ``'thermo'``
        (flat format, not yet per-condition lists).  Kinetics keys are
        RMG 1-based reaction indices (int).
    """
    result: Dict[str, dict] = {'time': [], 'kinetics': {}, 'thermo': {}}

    if not os.path.isdir(solver_dir):
        return result

    csv_files = sorted(f for f in os.listdir(solver_dir)
                       if 'sensitivity' in f and f.endswith('.csv'))

    time_set = False
    for csv_file in csv_files:
        csv_path = os.path.join(solver_dir, csv_file)
        time_data, kinetics, thermo = parse_rmg_sa_csv(csv_path)

        if not time_set and time_data:
            result['time'] = time_data
            time_set = True

        # Merge kinetics
        for obs, params in kinetics.items():
            if obs not in result['kinetics']:
                result['kinetics'][obs] = {}
            result['kinetics'][obs].update(params)

        # Merge thermo
        for obs, params in thermo.items():
            if obs not in result['thermo']:
                result['thermo'][obs] = {}
            result['thermo'][obs].update(params)

    return result


def rmg_sa_csvs_to_sa_dict(solver_dir: str,
                            chem_annotated_path: Optional[str] = None,
                            ) -> dict:
    """Parse RMG SA CSVs and return a standardized per-condition sa_dict.

    This is the main entry point for converting RMG CSV output into the
    format used by ``T3.sa_dict`` and understood by
    ``sa_dict_to_yaml`` / ``sa_dict_from_yaml``.

    Args:
        solver_dir: Path to the ``solver/`` directory containing
                    ``sensitivity_*.csv`` files.
        chem_annotated_path: Optional path to ``chem_annotated.inp``.
                             When provided, RMG reaction indices in the
                             CSV headers are mapped to Chemkin indices.

    Returns:
        A dict with keys ``'time'``, ``'kinetics'``, ``'thermo'``,
        each value being a list with one entry (single condition).
    """
    raw = parse_rmg_sa_csvs(solver_dir)

    # Optionally remap reaction indices from RMG numbering to Chemkin numbering
    if chem_annotated_path and os.path.isfile(chem_annotated_path):
        chem_to_rmg = get_chem_to_rmg_rxn_index_map(chem_annotated_path)
        rmg_to_chem = {v: k for k, v in chem_to_rmg.items()}
        remapped_kinetics: Dict[str, Dict[int, np.ndarray]] = {}
        for obs, params in raw['kinetics'].items():
            remapped_kinetics[obs] = {}
            for rmg_idx, values in params.items():
                chem_idx = rmg_to_chem.get(rmg_idx, rmg_idx)
                remapped_kinetics[obs][chem_idx] = np.array(values)
        kinetics = remapped_kinetics
    else:
        kinetics = {obs: {idx: np.array(v) for idx, v in params.items()}
                    for obs, params in raw['kinetics'].items()}

    thermo = {obs: {spc: np.array(v) for spc, v in params.items()}
              for obs, params in raw['thermo'].items()}

    time_data = np.array(raw['time'])

    return {'time': [time_data], 'kinetics': [kinetics], 'thermo': [thermo]}
