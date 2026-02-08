"""
t3 tests common module
"""

from typing import Optional
import os
import shutil

from arc.common import read_yaml_file
from arc.molecule.molecule import Molecule

from t3.common import EXAMPLES_BASE_PATH, PROJECTS_BASE_PATH
from t3.main import T3


def run_minimal(project: Optional[str] = None,
                project_directory: Optional[str] = None,
                iteration: Optional[int] = None,
                set_paths: bool = False,
                ) -> T3:
    """
    A helper function for running the minimal example.

    Args:
        project (str, optional): The project name.
        project_directory (str, optional): The project directory.
        iteration (int, optional): The iteration number.
        set_paths (bool, optional): Whether to set the paths.

    Returns:
        T3: The T3 object.
    """
    minimal_input = os.path.join(EXAMPLES_BASE_PATH, 'minimal', 'input.yml')
    input_dict = read_yaml_file(path=minimal_input)
    input_dict['verbose'] = 10
    input_dict['project_directory'] = project_directory \
                                      or os.path.join(PROJECTS_BASE_PATH, 'test_minimal_delete_after_usage')
    if project is not None:
        input_dict['project'] = project
    if 't3' in input_dict and 'sensitivity' in input_dict['t3']:
        input_dict['t3']['sensitivity']['adapter'] = 'CanteraConstantTP'
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
    
    Checks that:
    1. The number of expected radicals matches the number of unique SMILES generated
    2. All generated radicals have SMILES that match at least one expected SMILES
    """
    # Check unique SMILES count matches expected count
    unique_smiles = set([rad[1] for rad in radicals])
    assert len(expected_radicals) == len(unique_smiles), \
        f"Expected {len(expected_radicals)} unique SMILES, got {len(unique_smiles)}"
    
    # Check all generated radicals have matching SMILES in expected
    for rad_tuple in radicals:
        assert any(isomorphic_smiles(rad_tuple[1], expected_rad_tuple[1]) for expected_rad_tuple in expected_radicals), \
            f"Generated SMILES {rad_tuple[1]} not found in expected radicals"


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


def copy_model(model_path: str, name: str = '') -> str:
    """
    Copy a model to a temporary location.

    Args:
        model_path (str): The path to the model to copy.
        name (str, optional): A unique name for the copied model in addition to "temp_".

    Returns:
        str: The path to the copied model.
    """
    model_path = os.path.abspath(model_path)
    model_name = os.path.basename(model_path)
    model_dir = os.path.dirname(model_path)
    name = f'{name}_' if name else ''
    new_model_path = os.path.join(model_dir, f'temp_{name}' + model_name)
    shutil.copyfile(model_path, new_model_path)
    return new_model_path
