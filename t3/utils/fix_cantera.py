"""
t3 utils fix_cantera module
A module to automatically fix issues with RMG-generated Cantera files, mainly resolving mislabeled duplicate reactions.
"""
import os

from typing import List, Optional

import shutil
import time
import traceback

import cantera as ct

from arc.common import read_yaml_file, save_yaml_file


def get_traceback(model_path: str) -> Optional[str]:
    """
    Try loading the Cantera model and return the traceback if it fails.

    Args:
        model_path (str): The path to the cantera YAML model file.

    Returns:
        Optional[str]: The traceback if the model fails to load.
    """
    tb = None
    try:
        ct.Solution(model_path)
    except ct.CanteraError:
        tb = traceback.format_exc()
    return tb


def fix_cantera(model_path: str):
    """
    Fix a Cantera model that has incorrectly marked duplicate reactions.
    Creates a backup copy of the Cantera model and fixes the content of the original file in place.

    Args:
        model_path (str): The path to the cantera YAML model file.

    Returns:
        bool: Whether the model was fixed.
    """
    shutil.copyfile(model_path, model_path + '.bak')
    done, fixed = False, False
    counter = 0
    while not done and counter < 1000:
        counter += 1
        tb = get_traceback(model_path)
        if tb is None:
            done = True
            break
        else:
            if 'Undeclared duplicate reactions detected' in tb:
                fix_undeclared_duplicate_reactions(model_path, tb)
                fixed = True
            elif 'No duplicate found for declared duplicate reaction' in tb:
                fix_no_duplicate_found(model_path, tb)
                fixed = True
            else:
                print(f'Could not fix {model_path}:\n\n{tb}')
                break
        time.sleep(1)
    if fixed:
        print(f'Fixing Cantera model {model_path} (and creating a backup copy with a .bak extension).')
    else:
        os.remove(model_path + '.bak')
    return done


def fix_undeclared_duplicate_reactions(model_path: str, tb: str):
    """
    Fix a Cantera model that has undeclared duplicate reactions.

    Args:
        model_path (str): The path to the cantera YAML model file.
        tb (str): The traceback.
    """
    content = read_yaml_file(model_path)
    rxns = get_dup_rxn_indices(tb)
    print(f'Marking reactions {", ".join([str(r) for r in rxns])} as duplicate.')
    for i in rxns:
        content['reactions'][i - 1]['duplicate'] = True
    save_yaml_file(model_path, content)


def fix_no_duplicate_found(model_path: str, tb: str):
    """
    Fix a Cantera model that has a reaction marked as duplicate by mistake with no other duplicate reaction.

    Args:
        model_path (str): The path to the cantera YAML model file.
        tb (str): The traceback.
    """
    content = read_yaml_file(model_path)
    rxns = get_mistakenly_marked_dup_rxns(tb)
    for i in rxns:
        if 'duplicate' in content['reactions'][i].keys():
            print(f'Marking reaction {i} as non-duplicate.')
            del content['reactions'][i]['duplicate']
    save_yaml_file(model_path, content)


def get_dup_rxn_indices(tb: str) -> List[int]:
    """
    Get the duplicate reactions from the traceback.

    Args:
        tb (str): The traceback.

    Returns:
        List[int]: The reactions indices.
    """
    rxns = list()
    if tb is None:
        return rxns
    lines = tb.split('\n')
    read = False
    for line in lines:
        if 'Undeclared duplicate reactions detected:' in line:
            read = True
        if '|  Line |' in line:
            break
        if read and 'Reaction' in line:
            rxns.append(int(line.split()[1].split(':')[0]))
    return rxns


def get_mistakenly_marked_dup_rxns(tb: str) -> List[int]:
    """
    Get the duplicate reactions from the traceback.

    Args:
        tb (str): The traceback.

    Returns:
        List[int]: The reactions indices.
    """
    rxns = list()
    if tb is None:
        return rxns
    lines = tb.split('\n')
    for line in lines:
        if 'No duplicate found for declared duplicate reaction number' in line:
            rxns.append(int(line.split()[8]))
    return rxns
