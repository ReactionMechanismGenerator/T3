"""
T3 Cantera YAML Parser
A lightweight parser for Cantera YAML files.
"""

import logging
import os
from typing import Dict, List, Optional, Tuple

from arc.common import read_yaml_file

from t3.chem import T3Species, T3Reaction


def load_cantera_yaml_file(path: str,
                           species_dict_path: Optional[str] = None,
                           ) -> Tuple[List[T3Species], List[T3Reaction]]:
    """
    Load a Cantera YAML file and return a list of species and reactions.

    Args:
        path (str): The path to the Cantera YAML file.
        species_dict_path (Optional[str]): The path to the RMG species dictionary file.

    Returns:
        Tuple[List[T3Species], List[T3Reaction]]: The loaded species and reactions.
    """
    if not os.path.isfile(path):
        raise IOError(f'File {path} does not exist.')

    adjlists = {}
    if species_dict_path and os.path.isfile(species_dict_path):
        adjlists = parse_species_dictionary(species_dict_path)

    data = read_yaml_file(path)
    species_list, reactions_list = [], []

    species_data = data.get('species', [])
    for spc_datum in species_data:
        label = spc_datum.get('name')
        if not label:
            continue
        note = spc_datum.get('note', '')

        thermo_method = None
        thermo_source = None

        if 'Thermo Source:' in note:
            try:
                thermo_part = note.split('Thermo Source:')[1].split('|')[0].strip()
                if 'Thermo library:' in thermo_part:
                    thermo_method = 'Library'
                    thermo_source = thermo_part.split('Thermo library:')[1].strip()
                elif 'Thermo group additivity' in thermo_part:
                    thermo_method = 'GAV'
                    thermo_source = thermo_part.replace('Thermo group additivity estimation:', '').strip()
                elif 'QM' in thermo_part:
                    thermo_method = 'QM'
                    thermo_source = thermo_part.strip()
                else:
                    thermo_source = thermo_part
            except IndexError:
                pass

        adjlist = adjlists.get(label)
        t3_spc = T3Species(label=label,
                           thermo_method=thermo_method,
                           thermo_source=thermo_source,
                           thermo_comment=note,
                           adjlist=adjlist,
                           )
        species_list.append(t3_spc)
        
    reactions_data = data.get('reactions', [])
    for i, rxn_datum in enumerate(reactions_data):
        equation = rxn_datum.get('equation')
        if not equation:
            continue
        try:
            if '<=>' in equation:
                arrow = '<=>'
            elif '=>' in equation:
                arrow = '=>'
            else:
                arrow = '<=>'
            
            reactants_str, products_str = equation.split(arrow)
            reactants_labels = [r.strip() for r in reactants_str.split('+')]
            products_labels = [p.strip() for p in products_str.split('+')]
            
            reactants = [get_species_by_label(l, species_list) for l in reactants_labels]
            products = [get_species_by_label(l, species_list) for l in products_labels]
            
            reactants = [r for r in reactants if r is not None]
            products = [p for p in products if p is not None]
            
            if reactants and products:
                # Parse 'note' for metadata
                note = rxn_datum.get('note', '')
                kinetics_method = None
                kinetics_source = None

                if 'Source:' in note:
                    try:
                        source_part = note.split('Source:')[1].split('|')[0].strip()
                        if 'Library' in source_part:
                            kinetics_method = 'Library'
                            kinetics_source = source_part.replace('Library', '').strip()
                        elif 'Template family' in source_part:
                            kinetics_method = 'Rate Rules'
                            kinetics_source = source_part.replace('Template family', '').strip()
                        elif 'PDep' in source_part or 'Network' in source_part:
                            kinetics_method = 'PDep'
                            kinetics_source = source_part
                        else:
                            kinetics_source = source_part
                    except IndexError:
                        pass

                # Create reaction with labels instead of species objects
                # ARCReaction expects labels, not ARCSpecies objects
                rxn = T3Reaction(r_species=reactants,
                                 p_species=products,
                                 kinetics_method=kinetics_method,
                                 kinetics_source=kinetics_source,
                                 kinetics_comment=note)
                rxn.label = equation
                reactions_list.append(rxn)

        except Exception as e:
            logging.warning(f"Failed to parse reaction equation: {equation}. Error: {e}")
            continue

    return species_list, reactions_list


def parse_species_dictionary(path: str) -> Dict[str, str]:
    """
    Parse an RMG species dictionary file.

    Args:
        path (str): The path to the species dictionary file.

    Returns:
        Dict[str, str]: A dictionary mapping species labels to adjacency lists.
    """
    if not os.path.isfile(path):
        return {}

    with open(path, 'r') as f:
        lines = f.readlines()

    adjlists = {}
    current_label = None
    current_adjlist = []

    for line in lines:
        if line.strip() == "":
            if current_label and current_adjlist:
                adjlists[current_label] = "".join(current_adjlist).strip()
            current_label = None
            current_adjlist = []
        elif current_label is None:
            current_label = line.strip()
        else:
            current_adjlist.append(line)

    if current_label and current_adjlist:
        adjlists[current_label] = "".join(current_adjlist).strip()

    return adjlists


def get_species_by_label(label: str, species_list: List[T3Species]) -> Optional[T3Species]:
    """
    Get a species from a list by its label.

    Args:
        label (str): The label of the species to find.
        species_list (List[T3Species]): The list of species to search.

    Returns:
        Optional[T3Species]: The species with the matching label, or None if not found.
    """
    for spc in species_list:
        if spc.label == label:
            return spc
    return None
