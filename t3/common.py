"""
t3 common module
"""

import datetime
import os
import string
from typing import Dict, Optional, Tuple, Union

from rmgpy.species import Species

from arc.species.converter import molecules_from_xyz

VERSION = '0.1.0'

t3_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))  # absolute path to the T3 folder
DATA_BASE_PATH = os.path.join(t3_path, 'tests', 'data')
SIMULATE_DATA_BASE_PATH = os.path.join(t3_path, 'tests', 'test_simulate_adapters', 'data')
EXAMPLES_BASE_PATH = os.path.join(t3_path, 'examples')
SCRATCH_BASE_PATH = os.path.join(t3_path, 'tests', 'scratch')
IPYTHON_SIMULATOR_EXAMPLES_PATH = os.path.join(t3_path, 'ipython', 'simulator_adapter_examples')
PROJECTS_BASE_PATH = os.path.join(t3_path, 'Projects')
VALID_CHARS = "-_=.,%s%s" % (string.ascii_letters, string.digits)


def get_species_by_label(label: str,
                         species_list: list,
                         ) -> Optional[Species]:
    """
    Get a species from a list of species by its label.

    Args:
        label (str): A species label.
        species_list (list): Entries are RMG Species objects.

    Returns:
        Optional[Species]: The corresponding species object instance from the species_list.
                           Returns ``None`` if no species was found.
    """
    for species in species_list:
        if species.label == label or species.to_chemkin() == label:
            return species
    if '(' in label and ')' in label:
        # try by the RMG species index
        for species in species_list:
            if species.index == int(label.split('(')[-1].split(')')[0]):
                return species
    return None


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


def get_rmg_species_from_a_species_dict(species_dict: dict,
                                        raise_error: bool = False,
                                        ) -> Optional[Species]:
    """
    Get an RMG Species instance that corresponds to a species specified under the rmg.species
    section of the T3 input file (a species dictionary).

    Args:
        species_dict (dict): The species dictionary to process.
        raise_error (bool, optional): Whether to raise an error if a Species instance cannot be generated.
                                      Default: ``False``.

    Raises:
        ValueError: If the species dictionary does not have a specified structure (if ``raise_error`` is ``True``).

    Returns:
        Species: The corresponding RMG species instance.
    """
    species = None
    errored = False
    if species_dict['adjlist'] is not None:
        species = Species(label=species_dict['label']).from_adjacency_list(species_dict['adjlist'])
    elif species_dict['smiles'] is not None:
        species = Species(label=species_dict['label'], smiles=species_dict['smiles'])
    elif species_dict['inchi'] is not None:
        species = Species(label=species_dict['label'], inchi=species_dict['inchi'])
    elif species_dict['xyz'] is not None:
        for xyz in species_dict['xyz']:
            mol_bo = molecules_from_xyz(xyz=xyz)[1]
            if mol_bo is not None:
                species = Species(label=species_dict['label']).from_adjacency_list(mol_bo.to_adjacency_list())
                break
        else:
            errored = True
    else:
        errored = True
    if errored and raise_error:
        raise ValueError(f"The species corresponding to {species_dict['label']} does not have a specified structure.")
    return species


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
