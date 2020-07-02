#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_writer module
"""

import os

from arc.common import read_yaml_file

from t3.common import DATA_BASE_PATH, EXAMPLES_BASE_PATH
from t3.schema import InputBase, RMG, T3
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
    write_rmg_input_file(rmg=schema['rmg'],
                         t3=schema['t3'],
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
    terminationTime=(5.0, 's'),
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

simulator(atol=1e-16, rtol=1e-08, sens_atol=1e-06, sens_rtol=0.0001)

options(
    name='Seed',
    generateSeedEachIteration=True,
    saveSeedToDatabase=False,
    units='si',
    generateOutputHTML=False,
    generatePlots=False,
    saveSimulationProfiles=False,
    verboseComments=False,
    saveEdgeSpecies=True,
    keepIrreversible=False,
    trimolecularProductReversible=True,
    wallTime='01:00:00:00',
    saveSeedModulus=-1,
)
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

    write_rmg_input_file(rmg=schema['rmg'],
                         t3=schema['t3'],
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


def test_write_rmg_input_file_liquid():
    """Test writing an RMG input file for a liquid phase reactor"""
    rmg = {'database': {'thermo_libraries': ['BurkeH2O2', 'api_soup', 'thermo_DFT_CCSDTF12_BAC',
                                             'DFT_QCI_thermo', 'primaryThermoLibrary', 'CBS_QB3_1dHR', 'CurranPentane'],
                        'kinetics_libraries': ['BurkeH2O2inN2', 'api_soup', 'NOx2018', 'Klippenstein_Glarborg2016']},
           'species': [{'label': 'AIBN',
                        'smiles': 'CC(C)(C#N)/N=N/C(C)(C)C#N',
                        'concentration': 4.900e-6},
                       {'label': 'MeOH',
                        'smiles': 'CO',
                        'concentration': 0.0124},
                       {'label': 'water',
                        'smiles': 'O',
                        'concentration': 0.0278,
                        'solvent': True},
                       {'label': 'O2',
                        'smiles': '[O][O]',
                        'constant': True,
                        'concentration': 2.730e-7},
                       {'label': 'OHCH2OO',
                        'smiles': 'OCO[O]',
                        'SA_observable': True},
                       {'label': 'cyanoisopropylOO',
                        'smiles': 'N#CC(C)(C)O[O]',
                        'SA_observable': True},
                       {'label': 'N2',
                        'smiles': 'N#N',
                        'constant': True,
                        'concentration': 4.819e-7,
                        'reactive': False}],
           'reactors': [{'type': 'liquid batch constant T V',
                         'T': [293, 393],
                         'termination_time': [72, 'hrs']}],
           'model': {'core_tolerance': [0.001]},
           'options': {'save_edge': True, 'save_html': True},
           'species_constraints': {'max_C_atoms': 4,
                                   'max_O_atoms': 3,
                                   'max_N_atoms': 2,
                                   'max_Si_atoms': 0,
                                   'max_S_atoms': 0,
                                   'max_heavy_atoms': 10,
                                   'max_radical_electrons': 1}}

    t3 = {'sensitivity':
              {'adapter': 'RMG',
               'atol': 1e-6,
               'rtol': 1e-4,
               }
          }

    file_path = os.path.join(DATA_BASE_PATH, 'test_write_rmg_input_file_liquid.py')
    rmg_schema = RMG(**rmg).dict()  # fill in defaults
    t3_schema = T3(**t3).dict()     # fill in defaults

    write_rmg_input_file(rmg=rmg_schema,
                         t3=t3_schema,
                         iteration=1,
                         path=file_path)

    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in ["    thermoLibraries=['BurkeH2O2', 'api_soup', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo', 'primaryThermoLibrary', 'CBS_QB3_1dHR', 'CurranPentane'],\n",
                 "    kineticsDepositories='default',\n",
                 "    structure=SMILES('CC(C)(C#N)/N=N/C(C)(C)C#N'),\n",
                 "liquidReactor(\n",
                 "    temperature=[(293.0, 'K'), (393.0, 'K')],\n",
                 "    initialConcentrations={\n",
                 "        'water': (0.0278, 'mol/cm^3'),\n",
                 "        'AIBN': (4.9e-06, 'mol/cm^3'),\n",
                 "        'O2': (2.73e-07, 'mol/cm^3'),\n",
                 "        'cyanoisopropylOO': (0, 'mol/cm^3'),\n",
                 "    terminationTime=(72.0, 'hours'),\n",
                 "    nSims=12,\n",
                 "    constantSpecies=['O2', 'N2', ],\n",
                 "    toleranceMoveToCore=0.001,\n",
                 "    toleranceInterruptSimulation=0.001,\n",
                 "    filterThreshold=100000000.0,\n",
                 "simulator(atol=1e-16, rtol=1e-08, sens_atol=1e-06, sens_rtol=0.0001)\n",
                 "    generateOutputHTML=True,\n",
                 "    allowed=['input species', 'seed mechanisms', 'reaction libraries'],\n",
                 "    maximumCarbonAtoms=4,\n",
                 "solvation(solvent='water')\n",
                 ]:
        assert line in lines
    os.remove(file_path)


def test_write_rmg_input_file_seed_all_radicals():
    """Test writing an RMG input file while seeding all radicals for one species"""
    rmg = {'database': {'thermo_libraries': ['BurkeH2O2'],
                        'kinetics_libraries': ['BurkeH2O2inN2']},
           'species': [{'label': 'methylethylester',
                        'smiles': 'COCC',
                        'concentration': 1,
                        'seed_all_rads': ['radical', 'alkoxyl', 'peroxyl']},
                       {'label': 'O2',
                        'smiles': '[O][O]',
                        'concentration': 2},
                       {'label': 'N2',
                        'smiles': 'N#N',
                        'constant': True,
                        'concentration': 6,
                        'reactive': False}],
           'reactors': [{'type': 'gas batch constant T P',
                         'T': 1250,
                         'P': [1, 10],
                         'termination_time': [10, 's']}],
           'model': {'core_tolerance': [0.001]}}

    t3 = {'sensitivity':
              {'adapter': 'RMG',
               'atol': 1e-6,
               'rtol': 1e-4,
               }
          }

    file_path = os.path.join(DATA_BASE_PATH, 'test_write_rmg_input_file_seed_rads.py')
    rmg_schema = RMG(**rmg).dict()  # fill in defaults
    t3_schema = T3(**t3).dict()     # fill in defaults

    write_rmg_input_file(rmg=rmg_schema,
                         t3=t3_schema,
                         iteration=1,
                         path=file_path)

    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in ["    thermoLibraries=['BurkeH2O2'],\n",
                 "    label='methylethylester',\n",
                 "    label='methylethylester_radical_0',\n",
                 "    label='methylethylester_alkoxyl_0',\n",
                 "    label='methylethylester_peroxyl_0',\n",
                 "    label='methylethylester_radical_1',\n",
                 "    label='methylethylester_alkoxyl_1',\n",
                 "    label='methylethylester_peroxyl_1',\n",
                 "    label='methylethylester_radical_2',\n",
                 "    label='methylethylester_alkoxyl_2',\n",
                 "    label='methylethylester_peroxyl_2',\n",
                 ]:
        assert line in lines
    os.remove(file_path)
