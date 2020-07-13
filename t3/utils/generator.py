"""
t3 utils generator module

Used to generate specific species and reactions.
"""

from typing import List, Type

import rmgpy.molecule.element as elements
from rmgpy.molecule.molecule import Atom, Bond
from rmgpy.species import Species


def generate_radicals(species: Type[Species],
                      types: List[str],
                      react_aromatic_rings: bool = False,
                      ):
    """
    Generate all radicals for a species by radical type.

    Args:
        species (Species): The RMG Species instance to process.
                           The ``label`` attribute of ``species`` should not be empty.
        types (List[str]): Entries are types of radicals to return.
        react_aromatic_rings (bool, optional): Whether to also consider hydrogen atoms on aromatic rings.
                                               Default: ``False``.

    Returns:
        List[Tuple[str, str]]: Entries are tuples representing the generated radical species,
                               the first entry in the tuple is a label,
                               the second entry is the respective SMILES representation.
    """
    radicals, existing_radical_indices, relevant_radical_indices, output = list(), list(), list(), list()
    if species is None or len(species.molecule[0].atoms) == 1:
        return radicals

    species = species.copy(deep=True)
    species.generate_resonance_structures(keep_isomorphic=False, filter_structures=True)

    # generate all normal "radicals", whether requested or not
    for molecule in species.molecule:
        if not molecule.reactive:
            continue
        existing_radical_indices = [molecule.atoms.index(atom) for atom in molecule.atoms
                                    if atom.radical_electrons]
        for atom_1 in molecule.atoms:
            if atom_1.is_hydrogen():
                for atom_2, bond_12 in atom_1.edges.items():
                    if bond_12.is_single():
                        # skipping hydrogen bonds
                        break
                else:
                    continue
                if not react_aromatic_rings and any(bond.is_benzene() for bond in atom_2.edges.values()):
                    continue
                mol_copy = molecule.copy(deep=True)
                # We are about to change the connectivity of the atoms in the molecule,
                # which will invalidate any existing vertex connectivity information; thus we reset it.
                mol_copy.reset_connectivity_values()

                # get the corresponding bond_12 in mol_copy
                for atom_2_copy, bond_12_copy in mol_copy.atoms[molecule.atoms.index(atom_1)].edges.items():
                    if bond_12_copy.is_single():
                        # skipping hydrogen bonds
                        break
                else:
                    continue

                mol_copy.remove_bond(bond_12_copy)
                mol_splits = mol_copy.split()
                if len(mol_splits) == 2:
                    mol_1, mol_2 = mol_splits
                else:
                    # something went wrong, don't use these molecules
                    continue

                derivative_mol = mol_1 if len(mol_2.atoms) == 1 else mol_2

                radicals_added = 0
                for atom in derivative_mol.atoms:
                    theoretical_charge = elements.PeriodicSystem.valence_electrons[atom.symbol] \
                                         - atom.get_total_bond_order() \
                                         - atom.radical_electrons - \
                                         2 * atom.lone_pairs
                    if theoretical_charge == atom.charge + 1:
                        # we're missing a radical electron on this atom
                        atom.increment_radical()
                        radicals_added += 1
                if radicals_added != 1:
                    # something went wrong, don't use these molecules
                    continue
                derivative_mol.update(raise_atomtype_exception=False)
                species_from_derivative_mol = Species(molecule=[derivative_mol])
                species_from_derivative_mol.generate_resonance_structures(keep_isomorphic=False,
                                                                          filter_structures=True)

                for existing_radical in radicals:
                    species_from_existing_radical = Species(molecule=[existing_radical])
                    species_from_existing_radical.generate_resonance_structures(keep_isomorphic=False,
                                                                                filter_structures=True)
                    if species_from_derivative_mol.is_isomorphic(species_from_existing_radical):
                        break
                else:
                    radicals.append(derivative_mol)
                    index_shift = 1 if len(mol_1.atoms) == 1 else 0
                    radical_atom_index = [derivative_mol.atoms.index(atom)
                                          for atom in derivative_mol.atoms
                                          if atom.radical_electrons
                                          and derivative_mol.atoms.index(atom) + index_shift
                                          not in existing_radical_indices][0]
                    relevant_radical_indices.append(radical_atom_index)

    for i, radical_mol in enumerate(radicals):
        if 'radical' in types:
            output.append((f'{species.label}_radical_{i}', radical_mol.copy(deep=True).to_smiles()))
        if 'alkoxyl' in types:
            alkoxyl = radical_mol.copy(deep=True)
            oxygen = Atom(element='O', radical_electrons=1, charge=0, lone_pairs=2)
            alkoxyl.add_atom(oxygen)
            alkoxyl.atoms[relevant_radical_indices[i]].decrement_radical()
            new_bond = Bond(atom1=alkoxyl.atoms[relevant_radical_indices[i]], atom2=oxygen, order=1)
            alkoxyl.add_bond(new_bond)
            output.append((f'{species.label}_alkoxyl_{i}', alkoxyl.to_smiles()))
        if 'peroxyl' in types:
            peroxyl = radical_mol.copy(deep=True)
            oxygen_1 = Atom(element='O', radical_electrons=0, charge=0, lone_pairs=2)
            oxygen_2 = Atom(element='O', radical_electrons=1, charge=0, lone_pairs=2)
            peroxyl.add_atom(oxygen_1)
            peroxyl.add_atom(oxygen_2)
            peroxyl.atoms[relevant_radical_indices[i]].decrement_radical()
            new_bond_1 = Bond(atom1=peroxyl.atoms[relevant_radical_indices[i]], atom2=oxygen_1, order=1)
            new_bond_2 = Bond(atom1=oxygen_1, atom2=oxygen_2, order=1)
            peroxyl.add_bond(new_bond_1)
            peroxyl.add_bond(new_bond_2)
            output.append((f'{species.label}_peroxyl_{i}', peroxyl.to_smiles()))

    return output
