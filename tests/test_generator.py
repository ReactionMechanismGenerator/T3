#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_generator module
"""

from rmgpy.species import Species

from t3.utils.generator import generate_radicals
from tests.common import isomorphic_smiles, check_expected_generated_radicals


def test_generate_radicals():
    """Test generating radicals from a given species"""
    # test H, should not give any radicals
    radicals = generate_radicals(species=Species(label='H', smiles='[H]'),
                                 types=['radical', 'alkoxyl', 'peroxyl'],
                                 )
    assert radicals == list()

    # test CO, should not give any radicals
    radicals = generate_radicals(species=Species(label='CO', smiles='[C-]#[O+]'),
                                 types=['radical', 'alkoxyl', 'peroxyl'],
                                 )
    assert radicals == list()

    # test CH4, should only give 3 radicals, one of each type
    radicals = generate_radicals(species=Species(label='CH4', smiles='C'),
                                 types=['radical', 'alkoxyl', 'peroxyl'],
                                 )
    expected_radicals = [('CH4_radical_0', '[CH3]'),
                         ('CH4_alkoxyl_0', 'C[O]'),
                         ('CH4_peroxyl_0', 'CO[O]')]
    for rad_tuple, expected_rad_tuple in zip(radicals, expected_radicals):
        assert rad_tuple[0] == expected_rad_tuple[0]
        assert isomorphic_smiles(rad_tuple[1], expected_rad_tuple[1])

    # test CH4, 'radical' type only
    radicals = generate_radicals(species=Species(label='CH4', smiles='C'),
                                 types=['radical'],
                                 )
    expected_radicals = [('CH4_radical_0', '[CH3]')]
    for rad_tuple, expected_rad_tuple in zip(radicals, expected_radicals):
        assert rad_tuple[0] == expected_rad_tuple[0]
        assert isomorphic_smiles(rad_tuple[1], expected_rad_tuple[1])

    # test CH4, 'radical' type only with react_aromatic_rings set to False
    radicals = generate_radicals(species=Species(label='CH4', smiles='C'),
                                 types=['radical'],
                                 react_aromatic_rings=False,
                                 )
    expected_radicals = [('CH4_radical_0', '[CH3]')]
    for rad_tuple, expected_rad_tuple in zip(radicals, expected_radicals):
        assert rad_tuple[0] == expected_rad_tuple[0]
        assert isomorphic_smiles(rad_tuple[1], expected_rad_tuple[1])

    # test CH4, 'alkoxyl' type only
    radicals = generate_radicals(species=Species(label='CH4', smiles='C'),
                                 types=['alkoxyl'],
                                 )
    expected_radicals = [('CH4_alkoxyl_0', 'C[O]')]
    for rad_tuple, expected_rad_tuple in zip(radicals, expected_radicals):
        assert rad_tuple[0] == expected_rad_tuple[0]
        assert isomorphic_smiles(rad_tuple[1], expected_rad_tuple[1])

    # test CH4, 'peroxyl' type only
    radicals = generate_radicals(species=Species(label='CH4', smiles='C'),
                                 types=['peroxyl'],
                                 )
    expected_radicals = [('CH4_peroxyl_0', 'CO[O]')]
    for rad_tuple, expected_rad_tuple in zip(radicals, expected_radicals):
        assert rad_tuple[0] == expected_rad_tuple[0]
        assert isomorphic_smiles(rad_tuple[1], expected_rad_tuple[1])

    # test ethylamine
    radicals = generate_radicals(species=Species(label='ethylamine', smiles='NCC'),
                                 types=['radical', 'alkoxyl', 'peroxyl'],
                                 )
    expected_radicals = [('ethylamine_radical_0', 'C[CH]N'),
                         ('ethylamine_alkoxyl_0', 'CC([O])N'),
                         ('ethylamine_peroxyl_0', 'CC(O[O])N'),
                         ('ethylamine_radical_1', '[CH2]CN'),
                         ('ethylamine_alkoxyl_1', 'NCC[O]'),
                         ('ethylamine_peroxyl_1', 'NCCO[O]'),
                         ('ethylamine_radical_2', 'CC[NH]'),
                         ('ethylamine_alkoxyl_2', 'CCN[O]'),
                         ('ethylamine_peroxyl_2', 'CCNO[O]')]
    check_expected_generated_radicals(radicals, expected_radicals)

    # test benzyl alcohol
    radicals = generate_radicals(species=Species(label='benzyl_alcohol', smiles='c1ccccc1CO'),
                                 types=['radical', 'alkoxyl', 'peroxyl'],
                                 )
    expected_radicals = [('benzyl_alcohol_radical_0', '[O]CC1=CC=CC=C1'),
                         ('benzyl_alcohol_alkoxyl_0', '[O]OCC1=CC=CC=C1'),
                         ('benzyl_alcohol_peroxyl_0', '[O]OOCC1=CC=CC=C1'),
                         ('benzyl_alcohol_radical_1', 'O[CH]C1=CC=CC=C1'),
                         ('benzyl_alcohol_alkoxyl_1', '[O]C(O)C1=CC=CC=C1'),
                         ('benzyl_alcohol_peroxyl_1', '[O]OC(O)C1=CC=CC=C1')]
    check_expected_generated_radicals(radicals, expected_radicals)

    # test benzyl alcohol, react_aromatic_rings
    radicals = generate_radicals(species=Species(label='benzyl_alcohol', smiles='c1ccccc1CO'),
                                 types=['peroxyl'],
                                 react_aromatic_rings=True,
                                 )
    expected_radicals = [('benzyl_alcohol_peroxyl_0', '[O]OOCc1ccccc1'),
                         ('benzyl_alcohol_peroxyl_1', '[O]OC(O)c1ccccc1'),
                         ('benzyl_alcohol_peroxyl_2', '[O]Oc1ccccc1CO'),
                         ('benzyl_alcohol_peroxyl_4', '[O]Oc1cccc(CO)c1'),
                         ('benzyl_alcohol_peroxyl_5', '[O]Oc1ccc(CO)cc1')]
    check_expected_generated_radicals(radicals, expected_radicals)
