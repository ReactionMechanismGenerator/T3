"""
T3 Cantera YAML Parser
A lightweight parser for Cantera YAML files.
"""

import os
import re
from typing import Dict, List, Optional, Tuple

from arc.common import read_yaml_file

from t3.chem import T3Species, T3Reaction
from t3.utils.rmg_shim import PDepNetwork


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

    # Build species list and a map from the *original* YAML label to the
    # T3Species object.  ARC's check_label() may legalize the label during
    # T3Species.__init__ (e.g. 'OH(4)' → 'OH[4]'), so reaction equations
    # must be looked up via the original label, not the legalized one.
    yaml_label_to_species: Dict[str, T3Species] = {}

    species_data = data.get('species', [])
    for spc_datum in species_data:
        label = spc_datum.get('name')
        if not label:
            continue
        species_note = spc_datum.get('note', '')
        thermo_datum = spc_datum.get('thermo', {})
        thermo_note = thermo_datum.get('note', '') if isinstance(thermo_datum, dict) else ''

        # Combine notes for comment, prioritizing thermo_note
        note = thermo_note if thermo_note else species_note

        thermo_method, thermo_source = parse_thermo_comment(note)

        adjlist = adjlists.get(label)
        t3_spc = T3Species(label=label,
                           thermo_method=thermo_method,
                           thermo_source=thermo_source,
                           thermo_comment=note,
                           adjlist=adjlist,
                           )
        species_list.append(t3_spc)
        yaml_label_to_species[label] = t3_spc

    reactions_data = data.get('reactions', [])
    seen_equations = set()
    for i, rxn_datum in enumerate(reactions_data):
        equation = rxn_datum.get('equation')
        if not equation:
            continue
        # Cantera YAML files may have multiple entries with the same equation
        # (from different RMG rate-rule families or marked ``duplicate: true``).
        # Keep only the first entry per equation to match the merged-duplicate
        # indexing used elsewhere in T3.
        if equation in seen_equations:
            continue
        seen_equations.add(equation)
        if '<=>' in equation:
            arrow = '<=>'
        elif '=>' in equation:
            arrow = '=>'
        elif '=' in equation:
            arrow = '='
        else:
            continue  # no recognizable arrow, skip

        # Strip pressure-dependent collider notation before splitting by '+'.
        # Patterns: '(+M)', '(+N2)', '(+N2(32))', and bare 'M' bath gas.
        clean_eq = re.sub(r'\(\+[^)]*(?:\([^)]*\))?\)', '', equation)
        parts = clean_eq.split(arrow, maxsplit=1)
        if len(parts) != 2:
            continue
        reactants_str, products_str = parts
        reactants_labels = [r.strip() for r in reactants_str.split('+') if r.strip()]
        products_labels = [p.strip() for p in products_str.split('+') if p.strip()]
        # Also strip bare M bath gas that isn't a real species
        reactants_labels = [l for l in reactants_labels if l != 'M']
        products_labels = [l for l in products_labels if l != 'M']

        # Look up via original YAML label to handle ARC label legalization
        reactants = [yaml_label_to_species.get(l) for l in reactants_labels]
        products = [yaml_label_to_species.get(l) for l in products_labels]

        reactants = [r for r in reactants if r is not None]
        products = [p for p in products if p is not None]

        if reactants and products:
            kinetics = rxn_datum.get('rate-constant', {})
            note = rxn_datum.get('note', '')
            kinetics_method = None
            kinetics_source = None
            network = None

            if note:
                for line in note.split('\n'):
                    line = line.strip()
                    if line.startswith('Library reaction:'):
                        kinetics_method = 'Library'
                        kinetics_source = line.split(':', 1)[1].strip()
                        break
                    elif line.startswith('Template reaction:'):
                        kinetics_method = 'Rate Rules'
                        kinetics_source = line.split(':', 1)[1].strip()
                        break
                    elif line.startswith('PDep reaction:'):
                        kinetics_method = 'PDep'
                        kinetics_source = line.split(':', 1)[1].strip()
                        # Parse network index from "PDepNetwork #N"
                        network_match = re.search(r'#(\d+)', kinetics_source)
                        if network_match:
                            network = PDepNetwork(index=int(network_match.group(1)))
                        break

            # Create reaction with label strings for reactants/products (ARCReaction requirement)
            # and species objects for r_species/p_species (T3 tracking)
            is_pdep = bool(re.search(r'\(\+', equation)) or rxn_datum.get('type') in (
                'pressure-dependent-Arrhenius', 'Chebyshev', 'falloff', 'chemically-activated')
            rxn = T3Reaction(r_species=reactants,
                             p_species=products,
                             kinetics=kinetics,
                             kinetics_method=kinetics_method,
                             kinetics_source=kinetics_source,
                             kinetics_comment=note,
                             index=i,
                             is_pressure_dependent=is_pdep,
                             network=network,
                             )
            if '<=>' in equation:
                rxn.label = equation
            elif '=>' in equation:
                rxn.label = equation.replace('=>', '<=>')
            elif '=' in equation:
                rxn.label = equation.replace('=', '<=>', 1)
            else:
                rxn.label = equation  # Fallback
            reactions_list.append(rxn)

    return species_list, reactions_list


def parse_thermo_comment(note: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parse a thermo comment note to extract the method and source.
    """
    if not note:
        return None, None

    # 1. Handle "Thermo Source:" wrapper (recurse)
    if 'Thermo Source:' in note:
        try:
            inner_note = note.split('Thermo Source:')[1].split('|')[0].strip()
            return parse_thermo_comment(inner_note)
        except IndexError:
            pass

    thermo_method = None
    thermo_source = note.strip()

    # 2. Determine Method
    # If it contains radical/group modifications, it is ALWAYS GAV,
    # even if it started as a library entry.
    if '+ radical(' in note or '+ group(' in note or 'group additivity' in note.lower():
        thermo_method = 'GAV'
    elif 'Thermo library' in note:
        thermo_method = 'Library'
    elif 'QM' in note:
        thermo_method = 'QM'

    # 3. Clean Source String
    # We must strip the prefix regardless of the method determined above.
    # e.g. "Thermo library: Lib + radical" -> method=GAV, source="Lib + radical"

    prefixes_to_remove = [
        'Thermo group additivity estimation:',
        'Thermo library corrected for liquid phase:',
        'Thermo library:',
    ]

    for prefix in prefixes_to_remove:
        if prefix in note:
            # Split on the prefix and take the second part
            # using split instead of replace ensures we only remove the header
            parts = note.split(prefix, 1)
            if len(parts) > 1:
                thermo_source = parts[1].split('|')[0].strip()
            break

    # Final fallback for QM or clean cleanup
    if thermo_method == 'QM':
        thermo_source = note.strip()

    return thermo_method, thermo_source


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
