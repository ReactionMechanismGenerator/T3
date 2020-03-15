"""
t3 utils module
"""

import os
from typing import Optional, Union


def dict_to_str(dictionary: dict,
                level: Optional[int] = 0,
                ) -> str:
    """
    A helper function to log dictionaries in a pretty way.

    Args:
        dictionary (dict): A general python dictionary.
        level (Optional[int]): A recursion level counter, should only be used internally.

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


def combine_dicts(dict1: Union[dict, None],
                  dict2: Union[dict, None],
                  ) -> dict:
    """
    Combine two dictionaries.

    Args:
        dict1 (Union[dict, None]): Dictionary 1.
        dict2 (Union[dict, None]): Dictionary 2.

    Raises:
        TypeError: If both dictionaries are None.

    Returns:
        dict: The combined dictionary
    """
    if dict1 is None and dict2 is None:
        raise TypeError('Both dicts cannot be None')
    if dict1 is None:
        return dict2
    if dict2 is None:
        return dict1
    combination = dict1.copy()
    for key, val in dict2.items():
        combination[key] = val
    return combination


def delete_root_rmg_log(path: str) -> None:
    """
    Delete the 'RMG.log' file left in the root output directory, it's a left-over.

    Args:
        path (str): The path to the root output folder.
    """
    rmg_log_path = os.path.join(path, 'RMG.log')
    if os.path.isfile(rmg_log_path):
        os.remove(rmg_log_path)
