"""
t3 tests common module
"""

import os
from typing import Optional

from rmgpy.molecule import Molecule

from arc.common import read_yaml_file

from t3.common import EXAMPLES_BASE_PATH, PROJECTS_BASE_PATH
from t3.main import T3


def run_minimal(project: Optional[str] = None,
                project_directory: Optional[str] = None,
                iteration: Optional[int] = None,
                set_paths: bool = False,
                ) -> T3:
    """A helper function for running the minimal example"""
    minimal_input = os.path.join(EXAMPLES_BASE_PATH, 'minimal', 'input.yml')
    input_dict = read_yaml_file(path=minimal_input)
    input_dict['verbose'] = 10
    input_dict['project_directory'] = project_directory \
                                      or os.path.join(PROJECTS_BASE_PATH, 'test_minimal_delete_after_usage')
    if project is not None:
        input_dict['project'] = project
    t3 = T3(**input_dict)
    t3.iteration = iteration or 0
    if set_paths:
        t3.set_paths()
    return t3


def isomorphic_smiles(smiles_1: str,
                      smiles_2: str,
                      ) -> bool:
    """
    Check whether two SMILES strings represent isomorphic molecules.

    Args:
        smiles_1: A SMILES string.
        smiles_2: A SMILES string.

    Returns:
        bool: Whether the two SMILES strings represent isomorphic molecules.
    """
    mol_1 = Molecule(smiles=smiles_1)
    mol_2 = Molecule(smiles=smiles_2)
    return mol_1.is_isomorphic(mol_2)


def check_expected_generated_radicals(radicals: list, expected_radicals: list):
    """
    A helper function for testing the generator.
    """
    assert len(expected_radicals) == len(set([rad[0] for rad in radicals]))
    assert len(expected_radicals) == len(set([rad[1] for rad in radicals]))
    expected_labels = [expected_rad_tuple[0] for expected_rad_tuple in expected_radicals]
    for rad_tuple in radicals:
        assert rad_tuple[0] in expected_labels
        assert any(isomorphic_smiles(rad_tuple[1], expected_rad_tuple[1]) for expected_rad_tuple in expected_radicals)


def almost_equal(a: float,
                 b: float,
                 places: int = 4,
                 ratio: Optional[float] = None,
                 ) -> bool:
    """
    A helper function for testing almost equal assertions.

    Args:
        a (float): The first number.
        b (float): The second number.
        places (int, optional): The number of decimal places to round to. Default is 4.
        ratio (float, optional): The ratio between the two numbers. Default is None.
    """
    if ratio is not None:
        return abs(a - b) / abs(a) < ratio
    return round(abs(a - b), places) == 0
