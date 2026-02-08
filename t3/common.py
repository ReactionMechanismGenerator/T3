"""
t3 common module
"""

import re
import datetime
import os
import string
from typing import TYPE_CHECKING, Dict, List, Tuple, Union, Optional

import numpy as np
import yaml

if TYPE_CHECKING:
    from t3.chem import T3Species


VERSION = '0.1.0'

t3_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))  # absolute path to the T3 folder
DATA_BASE_PATH = os.path.join(t3_path, 'data')
TEST_DATA_BASE_PATH = os.path.join(t3_path, 'tests', 'data')
SIMULATE_TEST_DATA_BASE_PATH = os.path.join(t3_path, 'tests', 'test_simulate', 'data')
EXAMPLES_BASE_PATH = os.path.join(t3_path, 'examples')
SCRATCH_BASE_PATH = os.path.join(t3_path, 'tests', 'scratch')
IPYTHON_SIMULATOR_EXAMPLES_PATH = os.path.join(t3_path, 'ipython', 'simulator_adapter_examples')
PROJECTS_BASE_PATH = os.path.join(t3_path, 'Projects')
VALID_CHARS = "-_=.,%s%s" % (string.ascii_letters, string.digits)


def get_species_by_label(label: str,
                         species_list: list['T3Species'],
                         ) -> Optional['T3Species']:
    """
    Get a species from a list of species by its label.

    Args:
        label (str): A species label.
        species_list (list): Entries are T3 Species objects.

    Returns:
        Optional[T3Species]: The corresponding species object instance from the species_list.
                              Returns ``None`` if no species was found.
    """
    for species in species_list:
        if species.label == label or to_chemkin_label(species) == label:
            return species
    # ARC's check_label() legalizes '(' → '[' and ')' → ']'.
    # Try matching with that normalization so RMG-format labels still resolve.
    if '(' in label:
        bracket_label = label.replace('(', '[').replace(')', ']')
        for species in species_list:
            if species.label == bracket_label:
                return species
    if '(' in label and ')' in label:
        # try by the RMG species index
        for species in species_list:
            if species.index == int(label.split('(')[-1].split(')')[0]):
                return species
    return None


def to_chemkin_label(species: 'T3Species') -> str:
    """
    Return a string identifier for the provided `species` that can be used in a
    Chemkin file. Although the Chemkin format allows up to 16 characters for a
    species identifier, this function uses a maximum of 10 to ensure that all
    reaction equations fit in the maximum limit of 52 characters, compliant with the RMG convention.

    Args:
        species (T3Species): A T3 Species object.

    Returns:
        str: A Chemkin-compliant species label.
    """
    label = species.label
    if not getattr(species, 'reactive', True) and 0 < len(label) <= 10:
        return label
    index = getattr(species, 't3_index', -1)
    if index is None:
        index = getattr(species, 'index', -1)
    if index is None:
        index = -1

    if index == -1:
        if len(label) > 0 and not re.search(r'[^A-Za-z0-9\-_,\(\)\*#.:\[\]]+', label):
            if len(label) <= 16:
                return label
            else:
                return label
        elif species.mol is not None:
            return '{0}'.format(species.mol.get_formula())
    else:
        if len(label) > 0 and index >= 0 and not re.search(r'[^A-Za-z0-9\-_,\(\)\*#.:\[\]]+', label):
            name = '{0}({1:d})'.format(label, index)
            if len(name) <= 16:
                return name

        if species.mol is not None:
            name = '{0}({1:d})'.format(species.mol.get_formula(), index)
            if len(name) <= 16:
                return name
            if index >= 0:
                if 'X' in name:
                    name = 'SX({0:d})'.format(index)
                else:
                    name = 'S({0:d})'.format(index)
                if len(name) <= 16:
                    return name
    return label


def dict_to_str(dictionary: dict,
                level: int = 0,
                ) -> str:
    """
    A helper function to log dictionaries in a pretty way.

    Args:
        dictionary (dict): A general python dictionary.
        level (int): A recursion level counter, sets the visual indentation.

    Returns:
        str: A text representation for the dictionary.
    """
    message = ''
    for key, value in dictionary.items():
        if isinstance(value, dict):
            message += ' ' * level * 2 + str(key) + ':\n' + dict_to_str(value, level + 1)
        else:
            message += ' ' * level * 2 + str(key) + ': ' + str(value) + '\n'
    return message


def delete_root_rmg_log(project_directory: str) -> None:
    """
    Delete the 'RMG.log' file left in the root output directory, it's a left-over.

    Args:
        project_directory (str): The path to the root output folder.
    """
    rmg_log_path = os.path.join(project_directory, 'RMG.log')
    if os.path.isfile(rmg_log_path):
        os.remove(rmg_log_path)


def time_lapse(t0: datetime.datetime) -> datetime.timedelta:
    """
    A helper function returning the elapsed time since t0.

    Args:
        t0 (datetime.datetime): The initial time the count starts from.

    Returns: datetime.timedelta
        The time difference between now and t0.
    """
    return datetime.datetime.now() - t0


def convert_termination_time_to_seconds(termination_time: Tuple[float, str]):
    """
    Converts the termination_time tuple from the RMG reactor to seconds.
    This is necessary for the RMS adapters since the Julia solver expects
    the integration bounds to be in units of seconds.

    Args:
        termination_time (Tuple[float, str]):  Termination time for simulating in the RMG reactor. Example: [5, 'hours']

    Returns:
        t_final (float): The termination time in seconds.
    """
    unit_conversion = {'micro-s': 1e-6,
                       'ms': 1e-3,
                       's': 1,
                       'hrs': 3600,
                       'hours': 3600,
                       'days': 3600*24,
                       }
    t_final, units = termination_time
    t_final = t_final * unit_conversion[units]
    return t_final


def get_values_within_range(value_range: Union[int, float, list, tuple],
                            num: int,
                            use_log_scale: bool = False,
                            ) -> list:
    """
    Get evenly dispersed values within a given range.
    If ``value_range`` is a number, the function returns it as a single-entry list.

    Args:
         value_range (Union[list, tuple]): A range of values, represented by two numbers (min, max).
         num (int): The number of desired repetitions within the value range.
         use_log_scale (bool, optional): Whether to disperse the points logarithmically (e.g., 0.1, 1, 10, 100).
                                         If not, points are dispersed linearly (e.g., 0, 5, 10).

    Returns:
        List[Union[int, float]]: The dispersed values.
    """
    if isinstance(value_range, (int, float)) or len(value_range) == 1:
        return [value_range] if isinstance(value_range, (int, float)) else value_range
    if num <= 0:
        raise ValueError(f'num must be a positive value, got: {num}')
    min_val = min(value_range)
    if use_log_scale:
        max_val = max(value_range)
        return [min_val * 10 ** i for i in range(num) if min_val * 10 ** i <= max_val]
    interval = get_interval(value_range=value_range, num=num)
    if num <= 2:
        return [min_val + (i + 1) * interval for i in range(num)]
    return [min_val + i * interval for i in range(num)]


def get_interval(value_range: Union[list, tuple],
                 num: int,
                 ) -> int:
    """
    Get an interval within a range of values that can be repeated ``num`` times.
    In case ``num` equals 1 or 2 we do not place points at the ends of the interval. E.g.:
    1:  [     *     ]
    2:  [   *   *   ]
    3:  [*    *    *]
    4:  [*  *  *  *]

    Args:
         value_range (Union[list, tuple]): A range of values, represented by two numbers (min, max).
         num (int): The number of desired repetitions within the value range.

    Returns:
        int: The interval.
    """
    num_of_intervals = num + 1 if num <= 2 else num - 1
    return (max(value_range) - min(value_range)) / num_of_intervals


def get_chem_to_rmg_rxn_index_map(chem_annotated_path: str) -> Dict[int, int]:
    """
    Get a dictionary that maps "Chemkin" reaction indices to "RMG" reaction indices.
    A Chemkin file counts duplicate reactions by design, while RMG treats duplicate reactions as a single reaction
    with two rate coefficients (or more) that are summed up.
    When T3 reads the reactions from an RMG run, it gets N reactions, N being the number of reactions by the RMG count.
    However, when RMG performs SA, it gives the Chemkin reaction index to each rate coefficient it perturbs
    (it actually perturbs each of the duplicate reaction rate coefficients).

    An example Chemkin file that highlights the different indexing:

        ! Reaction index: Chemkin #4; RMG #4
        ! Library reaction: primaryH2O2
        ! Flux pairs: O2(2), HO2(8); H2(5), HO2(8); O2(2), HO2(8);
        O2(2)+O2(2)+H2(5)=HO2(8)+HO2(8)                     2.000000e+17 0.000     25.830

        ! Reaction index: Chemkin #5; RMG #5
        ! Library reaction: primaryH2O2
        HO2(8)+HO2(8)=O2(2)+H2O2(9)                         1.030000e+14 0.000     11.040
        DUPLICATE
        ! Reaction index: Chemkin #6; RMG #5
        ! Library reaction: primaryH2O2
        HO2(8)+HO2(8)=O2(2)+H2O2(9)                         1.940000e+11 0.000     -1.409
        DUPLICATE

        ! Reaction index: Chemkin #7; RMG #6
        ! Library reaction: primaryH2O2
        ! Flux pairs: HO2(8), O2(2); OH(10), H2O(6);
        OH(10)+HO2(8)=O2(2)+H2O(6)                          2.140000e+06 1.650     2.180

    Returns:
        Dict[int, int]: A map between Chemkin to RMG reaction indices.
    """
    rxn_map = dict()
    if os.path.isfile(chem_annotated_path):
        with open(chem_annotated_path, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if 'Reaction index:' in line:
                splits = line.split()
                rxn_map[int(splits[4].split('#')[1].split(';')[0])] = int(splits[-1].split('#')[1])
    return rxn_map


def get_observable_label_from_header(header: str) -> str:
    """
    Get the observable label from a header in an RMG SA csv file.

    Args:
        header (str): The header from an SA csv file, e.g. ``dln[OH(4)]/dln[k1]: ...``.

    Returns:
        str: The observable label.

    Raises:
        ValueError: If the header does not contain the expected ``[label]`` pattern.
    """
    try:
        return header.split('[')[1].split(']')[0]
    except IndexError:
        raise ValueError(f"Could not parse observable label from SA header: {header!r}")


def get_parameter_from_header(header: str) -> Optional[str]:
    """
    Get the parameter label from a header in an SA csv file.
    parameter extraction examples:
    for species thermo (RMG) get 'C2H4(8)' from ``dln[ethane(1)]/dG[C2H4(8)]``
    for species enthalpy (Cantera) get 'C2H4(8)' from ``dln[ethane(1)]/dH[C2H4(8)]``
    for reaction, get k8 from ``dln[ethane(1)]/dln[k8]: H(6)+ethane(1)=H2(12)+C2H5(5)``

    Args:
        header (str): The header from an SA csv file.

    Returns:
        Optional[str]: The parameter label.
    """
    for prefix in ('/dG[', '/dH[', '/dln['):
        start_pos = header.find(prefix)
        if start_pos != -1:
            start_pos += len(prefix)
            break
    else:
        return None
    bracket_count = 1
    text = []
    for i in range(start_pos, len(header)):
        if header[i] == '[':
            bracket_count += 1
        elif header[i] == ']':
            bracket_count -= 1
            if bracket_count == 0:
                break
        text.append(header[i])
    return ''.join(text)


def sa_dict_to_yaml(sa_dict: dict,
                            metadata: Optional[dict] = None,
                            ) -> dict:
    """
    Convert an in-memory SA dictionary (with numpy arrays) into a plain-Python
    dictionary that can be written to YAML via ``save_yaml_file``.

    The in-memory format stores per-condition lists::

        sa_dict['time']      = [time_array_cond0, ...]
        sa_dict['kinetics']  = [kinetics_dict_cond0, ...]
        sa_dict['thermo']    = [thermo_dict_cond0, ...]

    Args:
        sa_dict (dict): The SA dictionary produced by a simulate adapter's
                        ``get_sa_coefficients()`` method.
        metadata (dict, optional): Extra metadata to embed (adapter, iteration, etc.).

    Returns:
        dict: A YAML-safe dictionary with a ``conditions`` list.
    """
    def _to_list(obj):
        if hasattr(obj, 'tolist'):
            return obj.tolist()
        return list(obj)

    result = dict()
    if metadata:
        result['metadata'] = metadata

    conditions = []
    time_list = sa_dict.get('time', [])
    kin_list = sa_dict.get('kinetics', [])
    thermo_list = sa_dict.get('thermo', [])

    for i in range(len(time_list)):
        cond = dict()
        cond['time'] = _to_list(time_list[i])
        for section, src in [('kinetics', kin_list), ('thermo', thermo_list)]:
            cond[section] = dict()
            if i < len(src):
                for observable_label, params in src[i].items():
                    cond[section][observable_label] = dict()
                    for parameter, values in params.items():
                        cond[section][observable_label][parameter] = _to_list(values)
        conditions.append(cond)

    result['conditions'] = conditions
    return result


def sa_dict_from_yaml(raw: dict) -> dict:
    """
    Convert a YAML-loaded SA dictionary back into the in-memory format
    expected by T3 (with numpy arrays for coefficient data).

    Supports both the per-condition format (``conditions`` key) and the
    legacy flat format (``time``, ``kinetics``, ``thermo`` at top level).

    Args:
        raw (dict): The dictionary as returned by ``read_yaml_file`` on an
                    ``sa_coefficients.yml`` file.

    Returns:
        dict: An SA dictionary with per-condition lists.
    """
    sa_dict: Dict[str, list] = {'time': [], 'kinetics': [], 'thermo': []}

    if 'conditions' in raw:
        for cond in raw['conditions']:
            sa_dict['time'].append(np.array(cond.get('time', [])))
            for section in ('kinetics', 'thermo'):
                section_dict = dict()
                for observable_label, params in cond.get(section, {}).items():
                    section_dict[observable_label] = dict()
                    for parameter, values in params.items():
                        key = int(parameter) if section == 'kinetics' else parameter
                        section_dict[observable_label][key] = np.array(values)
                sa_dict[section].append(section_dict)
    else:
        # Legacy flat format — wrap in single-condition list
        sa_dict['time'].append(np.array(raw.get('time', [])))
        for section in ('kinetics', 'thermo'):
            section_dict = dict()
            for observable_label, params in raw.get(section, {}).items():
                section_dict[observable_label] = dict()
                for parameter, values in params.items():
                    key = int(parameter) if section == 'kinetics' else parameter
                    section_dict[observable_label][key] = np.array(values)
            sa_dict[section].append(section_dict)

    return sa_dict


def save_yaml_file(path: str,
                   content: Union[list, dict],
                   top_keys: Optional[List[str]] = None,
                   ) -> None:
    """
    Save a YAML file with control over key ordering.

    Args:
        path (str): The YAML file path to save.
        content (Union[list, dict]): The content to save.
        top_keys (Optional[List[str]]): Keys to write first, in order.
                                        Remaining keys follow in their original order.
                                        Only applies when *content* is a dict.
    """
    if os.path.dirname(path) and not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
    if top_keys and isinstance(content, dict):
        with open(path, 'w') as f:
            for key in top_keys:
                if key in content:
                    yaml.dump({key: content[key]}, f, default_flow_style=False, sort_keys=False)
                    f.write('\n')
            remainder = {k: v for k, v in content.items() if k not in top_keys}
            if remainder:
                yaml.dump(remainder, f, default_flow_style=False, sort_keys=False)
    else:
        yaml_str = yaml.dump(data=content, default_flow_style=False, sort_keys=False)
        with open(path, 'w') as f:
            f.write(yaml_str)
