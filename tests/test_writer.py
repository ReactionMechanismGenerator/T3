#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_writer module
"""

import os

from arc.common import read_yaml_file

from t3.common import DATA_BASE_PATH, EXAMPLES_BASE_PATH
from t3.schema import InputBase
from t3.utils.writer import to_camel_case, write_rmg_input_file


def test_to_camel_case():
    """Test converting underscore lowercase to camel case"""
    uvs = ['seed_mechanisms',
           'conditions_per_iteration',
           'tolerance_interrupt_simulation',
           'tolerance_move_edge_reaction_to_surface_interrupt']
    ccvs = [to_camel_case(uv) for uv in uvs]
    assert ccvs == ['seedMechanisms',
                    'conditionsPerIteration',
                    'toleranceInterruptSimulation',
                    'toleranceMoveEdgeReactionToSurfaceInterrupt']


def test_write_rmg_input_file():
    """Test writing an RMG input file"""
    minimal_input = os.path.join(EXAMPLES_BASE_PATH, 'minimal', 'input.yml')
    input_dict = read_yaml_file(path=minimal_input)
    minimal_input_file_path = os.path.join(DATA_BASE_PATH, 'minimal_input.py')
    schema = InputBase(project=input_dict['project'],
                       project_directory=DATA_BASE_PATH,
                       t3=input_dict['t3'],
                       rmg=input_dict['rmg'],
                       qm=input_dict['qm'],
                       verbose=20,
                       ).dict()
    write_rmg_input_file(kwargs=schema['rmg'],
                         iteration=2,
                         path=minimal_input_file_path,
                         walltime='01:00:00:00')

    # minimal input file
    with open(minimal_input_file_path, 'r') as f:
        content = f.read()
    expected_input = """database(
    thermoLibraries=['primaryThermoLibrary'],
    reactionLibraries=[],
    transportLibraries=['OneDMinN2', 'PrimaryTransportLibrary', 'NOx2018', 'GRI-Mech'],
    seedMechanisms=[],
    kineticsDepositories='default',
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)

species(
    label='H2',
    reactive=True,
    structure=SMILES('[H][H]'),
)

species(
    label='O2',
    reactive=True,
    structure=SMILES('[O][O]'),
)

species(
    label='H',
    reactive=True,
    structure=SMILES('[H]'),
)

species(
    label='OH',
    reactive=True,
    structure=SMILES('[OH]'),
)

simpleReactor(
    temperature=(1000.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        'H2': 0.67,
        'O2': 0.33,
        'H': 0,
        'OH': 0,
    },
    terminationConversion={'H2': 0.9},
    terminationTime=(1000000.0, 's'),
    nSims=12,
)

model(
    toleranceMoveToCore=0.001,
    toleranceInterruptSimulation=0.001,
    filterReactions=False,
    filterThreshold=100000000.0,
    maxNumObjsPerIter=1,
    terminateAtMaxObjects=False,
)

simulator(atol=1e-16, rtol=1e-08)
"""
    assert content == expected_input

    # other keys not in the minimal example
    schema['rmg']['pdep'] = {
        'method': 'CSE',
        'max_grain_size': 1.5,
        'max_number_of_grains': 250,
        'T': [300, 2500, 10],
        'P': [1, 100, 10],
        'interpolation': 'Chebyshev',
        'T_basis_set': 6,
        'P_basis_set': 4,
        'max_atoms': 10,
    }
    schema['rmg']['options'] = {
        'seed_name': 'Seed',
        'save_edge': True,
        'save_html': True,
        'generate_seed_each_iteration': False,
        'save_seed_to_database': False,
        'units': 'si',
        'generate_plots': False,
        'save_simulation_profiles': False,
        'verbose_comments': False,
        'keep_irreversible': False,
        'trimolecular_product_reversible': True,
        'save_seed_modulus': -1,
    }
    schema['rmg']['species_constraints'] = {
        'allowed': ['input species', 'seed mechanisms', 'reaction libraries'],
        'max_C_atoms': 10,
        'max_O_atoms': 4,
        'max_N_atoms': 0,
        'max_Si_atoms': 0,
        'max_S_atoms': 1,
        'max_heavy_atoms': 14,
        'max_radical_electrons': 2,
        'max_singlet_carbenes': 1,
        'max_carbene_radicals': 0,
        'allow_singlet_O2': False,
    }

    write_rmg_input_file(kwargs=schema['rmg'],
                         iteration=2,
                         path=minimal_input_file_path,
                         walltime='01:00:00:00')

    # minimal input file
    with open(minimal_input_file_path, 'r') as f:
        lines = f.readlines()
    for line in [
        "pressureDependence(\n",
        "    method='chemically-significant eigenvalues',\n",
        "    maximumGrainSize=(1.5, 'kJ/mol'),\n",
        "    temperatures=(300, 2500, 'K', 10),\n",
        "    pressures=(1, 100, 'bar', 10),\n",
        "    interpolation=('Chebyshev', 6, 4),\n",
        "options(\n",
        "    generateSeedEachIteration=False,\n",
        "    generateOutputHTML=True,\n",
        "    wallTime='01:00:00:00',\n",
        "generatedSpeciesConstraints(\n",
        "    allowed=['input species', 'seed mechanisms', 'reaction libraries'],\n",
        "    maximumCarbonAtoms=10,\n",
        "    maximumOxygenAtoms=4,\n",
        "    maximumNitrogenAtoms=0,\n",
        "    maximumSulfurAtoms=1,\n",
        "    maximumHeavyAtoms=14,\n",
        "    maximumRadicalElectrons=2,\n",
        "    maximumSingletCarbenes=1,\n",
    ]:
        assert line in lines

    os.remove(minimal_input_file_path)



