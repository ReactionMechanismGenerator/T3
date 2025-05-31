"""
t3 utils rmg module
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple

import os

from rmgpy.chemkin import load_chemkin_file

if TYPE_CHECKING:
    from rmgpy.reaction import Reaction


def get_reaction_kinetics(
        reaction: 'Reaction',
        chemkin_path: str,
        species_dict_path: str,
) -> Optional[Any]:
    """
    Get the reaction kinetics from a Chemkin file.

    Args:
        reaction (Reaction): The reaction to consider.
        chemkin_path (str): Path to the Chemkin file.
        species_dict_path (str): Path to the species dictionary file.

    Returns:
        Optional[Any]: The reaction kinetics object.
    """
    if not os.path.exists(chemkin_path):
        raise FileNotFoundError(f"Chemkin file not found: {chemkin_path}")
    if not os.path.exists(species_dict_path):
        raise FileNotFoundError(f"Species dictionary file not found: {species_dict_path}")

    chemkin_species_list, chemkin_reaction_list = load_chemkin_file(chemkin_path, species_dict_path)

    for chemkin_reaction in chemkin_reaction_list:
        if reaction.is_isomorphic(other=chemkin_reaction, save_order=True):
            return chemkin_reaction.kinetics
    return None
