#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_generator module
"""

from rmgpy.species import Species

from t3.utils.generator import generate_radicals
from tests.common import isomorphic_smiles


def test_generate_radicals():
    """Test generating radicals from a given species"""
    species_list = [Species(label='H', smiles='[H]'),
                    Species(label='CH4', smiles='C'),
                    Species(label='ethylamine', smiles='NCC'),
                    Species(label='imipramine', smiles='CN(C)CCCN1C2=CC=CC=C2CCC3=CC=CC=C31'),
                    ]

    # test H, should not give any radicals
    radicals = generate_radicals(species=Species(label='H', smiles='[H]'),
                                 types=['radical', 'alkoxyl', 'peroxyl'],
                                 )
    assert radicals == []

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
    for rad_tuple, expected_rad_tuple in zip(radicals, expected_radicals):
        assert rad_tuple[0] == expected_rad_tuple[0]
        assert isomorphic_smiles(rad_tuple[1], expected_rad_tuple[1])

    # test imipramine
    radicals = generate_radicals(species=Species(label='imipramine', smiles='CN(C)CCCN1C2=CC=CC=C2CCC3=CC=CC=C31'),
                                 types=['radical', 'alkoxyl', 'peroxyl'],
                                 )
    expected_radicals = [('imipramine_radical_0', 'CN(C[CH]CN1C2=CC=CC=C2CCC2=C1C=CC=C2)C'),
                        ('imipramine_alkoxyl_0', 'CN(CC(CN1C2C=CC=CC=2CCC2C1=CC=CC=2)[O])C'),
                        ('imipramine_peroxyl_0', '[O]OC(CN(C)C)CN1C2=CC=CC=C2CCC2=C1C=CC=C2'),
                        ('imipramine_radical_1', 'CN([CH]CCN1C2C=CC=CC=2CCC2C1=CC=CC=2)C'),
                        ('imipramine_alkoxyl_1', 'CN(C(CCN1C2=CC=CC=C2CCC2=C1C=CC=C2)[O])C'),
                        ('imipramine_peroxyl_1', '[O]OC(N(C)C)CCN1C2=CC=CC=C2CCC2=C1C=CC=C2'),
                        ('imipramine_radical_2', 'CN(CC[CH]N1C2=CC=CC=C2CCC2=C1C=CC=C2)C'),
                        ('imipramine_alkoxyl_2', 'CN(CCC(N1C2=CC=CC=C2CCC2=C1C=CC=C2)[O])C'),
                        ('imipramine_peroxyl_2', '[O]OC(N1C2=CC=CC=C2CCC2=C1C=CC=C2)CCN(C)C'),
                        ('imipramine_radical_3', 'CN(CCCN1C2C=CC=CC=2[CH]CC2C1=CC=CC=2)C'),
                        ('imipramine_alkoxyl_3', 'CN(CCCN1C2=CC=CC=C2CC(C2=C1C=CC=C2)[O])C'),
                        ('imipramine_peroxyl_3', '[O]OC1CC2=CC=CC=C2N(C2=C1C=CC=C2)CCCN(C)C'),
                        ('imipramine_radical_4', 'CN(CCCN1C2=CC=CC=C2CCC2=C1C=CC=C2)[CH2]'),
                        ('imipramine_alkoxyl_4', '[O]CN(CCCN1C2=CC=CC=C2CCC2=C1C=CC=C2)C'),
                        ('imipramine_peroxyl_4', '[O]OCN(CCCN1C2=CC=CC=C2CCC2=C1C=CC=C2)C')]
    for rad_tuple, expected_rad_tuple in zip(radicals, expected_radicals):
        assert rad_tuple[0] == expected_rad_tuple[0]
        assert isomorphic_smiles(rad_tuple[1], expected_rad_tuple[1])

    # test imipramine, react_aromatic_rings
    radicals = generate_radicals(species=Species(label='imipramine', smiles='CN(C)CCCN1C2=CC=CC=C2CCC3=CC=CC=C31'),
                                 types=['peroxyl'],
                                 react_aromatic_rings=True,
                                 )
    expected_radicals = [('imipramine_peroxyl_0', '[O]OC(CN(C)C)CN1C2=CC=CC=C2CCC2=C1C=CC=C2'),
                         ('imipramine_peroxyl_1', '[O]OC(N(C)C)CCN1C2C=CC=CC=2CCC2C1=CC=CC=2'),
                         ('imipramine_peroxyl_2', '[O]OC(N1C2C=CC=CC=2CCC2C1=CC=CC=2)CCN(C)C'),
                         ('imipramine_peroxyl_3', '[O]OC1CC2=CC=CC=C2N(C2=C1C=CC=C2)CCCN(C)C'),
                         ('imipramine_peroxyl_4', '[O]OCN(CCCN1C2=CC=CC=C2CCC2=C1C=CC=C2)C'),
                         ('imipramine_peroxyl_5', '[O]OC1=CC=CC2=C1N(CCCN(C)C)C1=CC=CC=C1CC2'),
                         ('imipramine_peroxyl_6', '[O]OC1C=CC2=C(C=1)N(CCCN(C)C)C1=CC=CC=C1CC2'),
                         ('imipramine_peroxyl_7', '[O]OC1C=CC2=C(C=1)CCC1C(N2CCCN(C)C)=CC=CC=1'),
                         ('imipramine_peroxyl_8', '[O]OC1=CC=CC2=C1CCC1C(N2CCCN(C)C)=CC=CC=1')]
    for rad_tuple, expected_rad_tuple in zip(radicals, expected_radicals):
        assert rad_tuple[0] == expected_rad_tuple[0]
        assert isomorphic_smiles(rad_tuple[1], expected_rad_tuple[1])
