"""
t3 utils generator module

Used to generate specific species and reactions.
"""

from typing import List

from rmgpy.molecule.molecule import Atom, Bond
from rmgpy.species import Species

from arc.common import generate_resonance_structures
from arc.species import ARCSpecies


def generate_radicals(species: Species,
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
    existing_radical_indices, output, aromatic_rings, output_species = list(), list(), list(), list()
    if species is None or len(species.molecule[0].atoms) == 1 \
            or not any(atom.is_hydrogen() for atom in species.molecule[0].atoms):
        return output
    spc = ARCSpecies(label=species.label, adjlist=species.copy(deep=True).to_adjacency_list(), keep_mol=True)
    res_structures = generate_resonance_structures(spc.mol)
    spc.mol = res_structures[0] if res_structures is not None else spc.mol
    spc.mol.atoms = [a for a in spc.mol.atoms if not a.is_hydrogen()] + [a for a in spc.mol.atoms if a.is_hydrogen()]
    spc.final_xyz = spc.get_xyz(generate=True)
    spc.bdes = list()
    i = 0
    for atom_1 in spc.mol.atoms:
        if not atom_1.is_hydrogen():
            for atom_2, bond_12 in atom_1.edges.items():
                if atom_2.is_hydrogen() and bond_12.is_single():
                    # skipping hydrogen bonds
                    break
            else:
                continue
            if spc.mol.atoms.index(atom_1) in existing_radical_indices:
                continue
            if not react_aromatic_rings and any(bond.is_benzene() for bond in atom_1.edges.values()):
                continue
            existing_radical_indices.append(spc.mol.atoms.index(atom_1))
            i += 1
            spc.bdes.append((spc.mol.atoms.index(atom_1) + 1, spc.mol.atoms.index(atom_2) + 1))

    radicals = [rad for rad in spc.scissors(sort_atom_labels=True) if rad.label != 'H']

    rad_i, alkoxyl_i, peroxyl_i = 0, 0, 0
    for i, rad in enumerate(radicals):
        if 'radical' in types:
            if not any(rad.is_isomorphic(spc) for spc in output_species):
                output_species.append(rad)
                output.append((f'{species.label}_radical_{rad_i}', rad.mol.copy(deep=True).to_smiles()))
                rad_i += 1
        if 'alkoxyl' in types:
            alkoxyl = rad.copy()
            alkoxyl.mol_list = None
            oxygen = Atom(element='O', radical_electrons=1, charge=0, lone_pairs=2)
            alkoxyl.mol.add_atom(oxygen)
            alkoxyl.mol.atoms[existing_radical_indices[i]].decrement_radical()
            new_bond = Bond(atom1=alkoxyl.mol.atoms[existing_radical_indices[i]], atom2=oxygen, order=1)
            alkoxyl.mol.add_bond(new_bond)
            if not any(alkoxyl.is_isomorphic(spc) for spc in output_species):
                output_species.append(alkoxyl)
                output.append((f'{species.label}_alkoxyl_{alkoxyl_i}', alkoxyl.mol.to_smiles()))
                alkoxyl_i += 1
        if 'peroxyl' in types:
            peroxyl = rad.copy()
            peroxyl.mol_list = None
            oxygen_1 = Atom(element='O', radical_electrons=0, charge=0, lone_pairs=2)
            oxygen_2 = Atom(element='O', radical_electrons=1, charge=0, lone_pairs=2)
            peroxyl.mol.add_atom(oxygen_1)
            peroxyl.mol.add_atom(oxygen_2)
            peroxyl.mol.atoms[existing_radical_indices[i]].decrement_radical()
            new_bond_1 = Bond(atom1=peroxyl.mol.atoms[existing_radical_indices[i]], atom2=oxygen_1, order=1)
            new_bond_2 = Bond(atom1=oxygen_1, atom2=oxygen_2, order=1)
            peroxyl.mol.add_bond(new_bond_1)
            peroxyl.mol.add_bond(new_bond_2)
            if not any(peroxyl.is_isomorphic(spc) for spc in output_species):
                output_species.append(peroxyl)
                output.append((f'{species.label}_peroxyl_{peroxyl_i}', peroxyl.mol.to_smiles()))
                peroxyl_i += 1

    return output
