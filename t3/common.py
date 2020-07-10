"""
t3 common module
"""

import os
import string
from typing import Optional

from rmgpy.species import Species

from arc.species.converter import molecules_from_xyz

VERSION = '0.1.0'

t3_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))  # absolute path to the T3 folder
DATA_BASE_PATH = os.path.join(t3_path, 'tests', 'data')
EXAMPLES_BASE_PATH = os.path.join(t3_path, 'examples')
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
        Optional[Species]: The corresponding species from the species_list.
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
    else:
        errored = True
    if errored and raise_error:
        raise ValueError(f"The species corresponding to {species_dict['label']} does not have a specified structure.")
    return species
