#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_cantera_parser module
Tests for t3/utils/cantera_parser.py
"""

import os
import pytest

from t3.common import TEST_DATA_BASE_PATH
from t3.chem import T3Species, ThermoMethod
from t3.utils.cantera_parser import load_cantera_yaml_file, get_species_by_label, parse_species_dictionary, parse_thermo_comment


# Path to test Cantera YAML file and dictionary
TEST_CANTERA_YAML_PATH = os.path.join(TEST_DATA_BASE_PATH, 'models', 'test_cantera_parser.yaml')
TEST_SPECIES_DICT_PATH = os.path.join(TEST_DATA_BASE_PATH, 'models', 'test_species_dictionary.txt')


def test_load_cantera_yaml_file():
    """Test loading a Cantera YAML file and parsing species and reactions."""
    species_list, reactions_list = load_cantera_yaml_file(TEST_CANTERA_YAML_PATH,
                                                         species_dict_path=TEST_SPECIES_DICT_PATH)

    # Check species were loaded
    assert len(species_list) == 6

    # Check species labels
    species_labels = [spc.label for spc in species_list]
    assert 'H2' in species_labels
    assert 'O2' in species_labels
    assert 'H2O' in species_labels
    assert 'OH' in species_labels
    assert 'H' in species_labels
    assert 'O' in species_labels

    # Check adjlists were loaded
    h2_spc = next(spc for spc in species_list if spc.label == 'H2')
    assert h2_spc.mol is not None
    assert '1 H u0 p0 c0 {2,S}' in h2_spc.adjlist

    # Check multiplicity was determined from adjlist
    assert h2_spc.multiplicity == 1

    o2_spc = next(spc for spc in species_list if spc.label == 'O2')
    assert o2_spc.multiplicity == 3  # Ground state O2 is triplet

    h_spc = next(spc for spc in species_list if spc.label == 'H')
    assert h_spc.multiplicity == 2

    # Check thermo parsing for H2 (Library source)
    assert h2_spc.thermo_method == 'Library'
    assert 'primaryThermoLibrary' in h2_spc.thermo_source

    # Check thermo parsing for H2O (GAV source)
    h2o_spc = next(spc for spc in species_list if spc.label == 'H2O')
    assert h2o_spc.thermo_method == 'GAV'

    # Check thermo parsing for OH (QM source)
    oh_spc = next(spc for spc in species_list if spc.label == 'OH')
    assert oh_spc.thermo_method == 'QM'

    # Check reactions were loaded
    assert len(reactions_list) == 3

    # Check first reaction
    rxn1 = reactions_list[0]
    assert rxn1.label == 'H2 + O <=> H + OH'
    # r_species/p_species contain T3Species objects
    assert len(rxn1.r_species) == 2
    assert len(rxn1.p_species) == 2
    reactant_labels = [r.label for r in rxn1.r_species]
    product_labels = [p.label for p in rxn1.p_species]
    assert 'H2' in reactant_labels
    assert 'O' in reactant_labels
    assert 'H' in product_labels
    assert 'OH' in product_labels

    # Check kinetics parsing for first reaction (Library source)
    assert rxn1.kinetics_method == 'Library'
    assert 'BurkeH2O2' in rxn1.kinetics_source

    # Check kinetics parsing for second reaction (Rate Rules source)
    rxn2 = reactions_list[1]
    assert rxn2.kinetics_method == 'Rate Rules'

    # Check kinetics parsing for third reaction (PDep source)
    rxn3 = reactions_list[2]
    assert rxn3.kinetics_method == 'PDep'


def test_load_species_thermo_comments():
    """Test that thermo comments are correctly parsed for species."""
    species_list, _ = load_cantera_yaml_file(
        path=os.path.join(TEST_DATA_BASE_PATH, 'determine_species', 'iteration_2', 'RMG', 'cantera', 'chem_annotated.yaml'),
        species_dict_path=os.path.join(TEST_DATA_BASE_PATH, 'determine_species', 'iteration_2', 'RMG', 'chemkin', 'species_dictionary.txt'))

    for species in species_list:
        if species.label == 'H2(1)':
            assert species.thermo_method == ThermoMethod.LIBRARY
            assert species.thermo_source == 'primaryThermoLibrary'
            assert species.thermo_comment == 'Thermo library: primaryThermoLibrary'
        if species.label == 'H2O(7)':
            assert species.thermo_method == ThermoMethod.LIBRARY
            assert species.thermo_source == 'primaryThermoLibrary'
            assert species.thermo_comment == 'Thermo library: primaryThermoLibrary'
        if species.label == 'HO2(6)':
            assert species.thermo_method == ThermoMethod.GAV
            assert species.thermo_source == 'group(O2s-OsH) + group(O2s-OsH) + radical(HOOJ)'
            assert species.thermo_comment == 'Thermo group additivity estimation: group(O2s-OsH) + group(O2s-OsH) + radical(HOOJ)'
        if species.label == 'OH(4)':
            assert species.thermo_method == ThermoMethod.GAV
            assert species.thermo_source == 'primaryThermoLibrary + radical(HOJ)'
            assert species.thermo_comment == 'Thermo library: primaryThermoLibrary + radical(HOJ)'


def test_load_cantera_yaml_file_nonexistent():
    """Test that loading a non-existent file raises IOError."""
    with pytest.raises(IOError, match='does not exist'):
        load_cantera_yaml_file('/nonexistent/path/to/file.yaml')


@pytest.mark.parametrize("comment, expected_method, expected_source", [
    # spc_1: Standard GAV
    (
            'Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-CO) + '
            'missing(O2d-CO) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + '
            'group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)',
            'GAV',
            'group(O2s-(Cds-Cd)H) + missing(O2d-CO) + missing(O2d-CO) + group(Cds-Cds(Cds-O2d)O2s) + '
            'group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + '
            'group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)'
    ),

    # spc_2: Library + Radical
    (
            'Thermo library: primaryThermoLibrary + radical(HOOj)',
            'GAV',
            'primaryThermoLibrary + radical(HOOj)'
    ),

    # spc_3: Standard Library
    (
            'Thermo library: primaryThermoLibrary',
            'Library',
            'primaryThermoLibrary'
    ),

    # spc_4: Liquid Library (Solute)
    # Note: If parser is strict about "Thermo library:", this might need adjustment in the parser logic.
    (
            'Thermo library corrected for liquid phase: thermo_DFT_CCSDTF12_BAC + Solvation correction '
            'with water as solvent and solute estimated using Data from solute library',
            'Library',
            'thermo_DFT_CCSDTF12_BAC + Solvation correction with water as solvent and solute estimated using Data from solute library'
    ),

    # spc_5: Liquid Library (Abraham)
    (
            'Thermo library corrected for liquid phase: DFT_QCI_thermo + Solvation correction with '
            'water as solvent and solute estimated using abraham(Oss- noncyclic) + abraham(OssH) + '
            'nonacentered(OssH) + abraham(OssH) + nonacentered(OssH) + abraham(CssH2) + radical(ROOJ)',
            'GAV',
            'DFT_QCI_thermo + Solvation correction with water as solvent and solute estimated using '
            'abraham(Oss- noncyclic) + abraham(OssH) + nonacentered(OssH) + abraham(OssH) + '
            'nonacentered(OssH) + abraham(CssH2) + radical(ROOJ)'
    ),

    # spc_6: Liquid Library + Radical
    (
            'Thermo library corrected for liquid phase: api_soup + radical(CCsJO) + Solvation '
            'correction with water as solvent and solute estimated using abraham(Oss-noncyclic) + '
            'nonacentered(OxRing) + abraham(N3t) + abraham(Css-noH) + abraham(CssH2) + '
            'abraham(CssH3) + abraham(Ct)',
            'GAV',
            'api_soup + radical(CCsJO) + Solvation correction with water as solvent and solute estimated using '
            'abraham(Oss-noncyclic) + nonacentered(OxRing) + abraham(N3t) + abraham(Css-noH) + '
            'abraham(CssH2) + abraham(CssH3) + abraham(Ct)'
    ),

    # spc_7: Liquid GAV
    (
            'Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(N3t-CtCs) + '
            'group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + missing(Ct-CsN3t) + '
            'ring(Ethylene_oxide) + Solvation correction with water as solvent and solute estimated '
            'using abraham(Oss-noncyclic) + nonacentered(OxRing) + abraham(OssH) + nonacentered(OssH) + '
            'abraham(N3t) + abraham(Css-noH) + abraham(CssH) + abraham(CssH3) + abraham(Ct)',
            'GAV',
            'group(O2s-CsCs) + group(O2s-CsH) + group(N3t-CtCs) + group(Cs-(Cds-Cds)CsCsOs) + '
            'group(Cs-CsOsOsH) + group(Cs-CsHHH) + missing(Ct-CsN3t) + ring(Ethylene_oxide) + '
            'Solvation correction with water as solvent and solute estimated using '
            'abraham(Oss-noncyclic) + nonacentered(OxRing) + abraham(OssH) + nonacentered(OssH) + '
            'abraham(N3t) + abraham(Css-noH) + abraham(CssH) + abraham(CssH3) + abraham(Ct)'
    ),

    # spc_8: JetSurf Library
    (
            "Thermo library: JetSurF2.0",
            'Library',
            "JetSurF2.0"
    ),

    # spc_9: Simple GAV
    (
            "Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + "
            "group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + "
            "group(Cds-CdsHH) + group(Cds-CdsHH)",
            'GAV',
            "group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + "
            "group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH)"
    ),

    # spc_10: JetSurf Library + Radical
    (
            "Thermo library: JetSurF2.0 + radical(RCCJ)",
            'GAV',
            "JetSurF2.0 + radical(RCCJ)"
    ),
])
def test_parse_thermo_comment(comment, expected_method, expected_source):
    """
    Test parsing the thermo comment to extract method and source.
    """
    method, source = parse_thermo_comment(comment)

    assert method == expected_method

    if source and expected_source:
        assert source.replace('\n', ' ').strip() == expected_source.replace('\n', ' ').strip()
    else:
        assert source == expected_source


def test_get_species_by_label():
    """Test getting a species by label from a list."""
    # Create some test species with SMILES to avoid multiplicity issues
    species_list = [
        T3Species(label='H2', smiles='[H][H]'),
        T3Species(label='O2', smiles='[O][O]'),
        T3Species(label='H2O', smiles='O'),
    ]

    # Test finding existing species
    h2_spc = get_species_by_label('H2', species_list)
    assert h2_spc is not None
    assert h2_spc.label == 'H2'

    o2_spc = get_species_by_label('O2', species_list)
    assert o2_spc is not None
    assert o2_spc.label == 'O2'

    h2o_spc = get_species_by_label('H2O', species_list)
    assert h2o_spc is not None
    assert h2o_spc.label == 'H2O'


def test_get_species_by_label_not_found():
    """Test that getting a non-existent species returns None."""
    species_list = [
        T3Species(label='H2', smiles='[H][H]'),
        T3Species(label='O2', smiles='[O][O]'),
    ]

    result = get_species_by_label('NonExistent', species_list)
    assert result is None


def test_get_species_by_label_empty_list():
    """Test that getting a species from an empty list returns None."""
    result = get_species_by_label('H2', [])
    assert result is None


def test_parse_species_dictionary():
    """Test parsing an RMG species dictionary file."""
    dict_path = os.path.join(TEST_DATA_BASE_PATH, 'models', 'test_species_dictionary.txt')
    adjlists = parse_species_dictionary(dict_path)

    assert len(adjlists) == 6
    assert 'H' in adjlists
    assert 'O' in adjlists
    assert 'O2' in adjlists
    assert 'OH' in adjlists
    assert 'H2' in adjlists
    assert 'H2O' in adjlists

    # Check a specific adjlist (OH)
    assert '1 O u1 p2 c0 {2,S}' in adjlists['OH']
    assert '2 H u0 p0 c0 {1,S}' in adjlists['OH']
    assert 'multiplicity 2' in adjlists['OH']

