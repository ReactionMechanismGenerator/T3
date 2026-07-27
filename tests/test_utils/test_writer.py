#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_writer module
"""

import os

import pytest

from arc.common import read_yaml_file

from t3.common import TEST_DATA_BASE_PATH, EXAMPLES_BASE_PATH
from t3.schema import InputBase, RMG, T3
from t3.utils.writer import (
    rewrite_arkane_method_line,
    to_camel_case,
    write_arkane_network_input_file,
    write_rmg_input_file,
)

NETWORK_FIXTURE = os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'iteration_1', 'RMG', 'pdep', 'network4_1.py')


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
    minimal_input_file_path = os.path.join(TEST_DATA_BASE_PATH, 'minimal_input.py')
    schema = InputBase(project=input_dict['project'],
                       project_directory=TEST_DATA_BASE_PATH,
                       t3=input_dict['t3'],
                       rmg=input_dict['rmg'],
                       qm=input_dict['qm'],
                       verbose=20,
                       ).model_dump()
    write_rmg_input_file(rmg=schema['rmg'],
                         t3=schema['t3'],
                         iteration=2,
                         path=minimal_input_file_path,
                         walltime='01:00:00:00',
                         )

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
    },
    terminationConversion={'H2': 0.9},
    terminationTime=(5.0, 's'),
    nSims=12,
)

model(
    toleranceMoveToCore=0.001,
    toleranceInterruptSimulation=0.001,
    filterReactions=True,
    filterThreshold=100000000.0,
    maxNumObjsPerIter=1,
    terminateAtMaxObjects=False,
)

simulator(atol=1e-16, rtol=1e-08)

options(
    name='Seed',
    generateSeedEachIteration=True,
    saveSeedToDatabase=False,
    units='si',
    generateOutputHTML=False,
    generatePlots=False,
    saveSimulationProfiles=False,
    verboseComments=False,
    saveEdgeSpecies=False,
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
                        'concentration': [0.0278, 0.0502],
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
                         'V': [1, 10],
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
              {'adapter': 'CanteraConstantTP',
               'atol': 1e-6,
               'rtol': 1e-4,
               }
          }

    file_path = os.path.join(TEST_DATA_BASE_PATH, 'test_write_rmg_input_file_liquid.py')
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
                 "        'water': [(0.0278, 'mol/cm^3'), (0.0502, 'mol/cm^3')],\n",
                 "        'AIBN': (4.9e-06, 'mol/cm^3'),\n",
                 "        'O2': (2.73e-07, 'mol/cm^3'),\n",
                 "    terminationTime=(72.0, 'hours'),\n",
                 "    nSims=12,\n",
                 "    constantSpecies=['O2', 'N2', ],\n",
                 "    toleranceMoveToCore=0.001,\n",
                 "    toleranceInterruptSimulation=0.001,\n",
                 "    filterThreshold=100000000.0,\n",
                 "simulator(atol=1e-16, rtol=1e-08)\n",
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
                        'concentration': [2, 2.5]},
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
              {'adapter': 'CanteraConstantTP',
               'atol': 1e-6,
               'rtol': 1e-4,
               }
          }

    file_path = os.path.join(TEST_DATA_BASE_PATH, 'test_write_rmg_input_file_seed_rads.py')
    rmg_schema = RMG(**rmg).dict()  # fill in defaults
    t3_schema = T3(**t3).dict()     # fill in defaults

    write_rmg_input_file(rmg=rmg_schema,
                         t3=t3_schema,
                         iteration=1,
                         path=file_path)

    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in ["        'O2': [2.0, 2.5],\n",
                 "    thermoLibraries=['BurkeH2O2'],\n",
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


def test_rewrite_arkane_method_line_double_quoted():
    """B4(a): a method = "..." line (double-quoted) must be rewritten, not silently returned
    unchanged. split("'") sees no single quotes at all on a double-quoted line, so the old
    len(splits) >= 3 guard never fires."""
    line = '    method = "modified strong collision",\n'
    result = rewrite_arkane_method_line(line=line, method='CSE')
    assert 'chemically-significant eigenvalues' in result
    assert 'modified strong collision' not in result


def test_rewrite_arkane_method_line_malformed_quote_raises():
    """A line that IS a 'method = ...' candidate (matches METHOD_LINE_CANDIDATE_RE, so the
    per-file rewrite count treats it as 'found') but does not have the expected quoted-assignment
    shape (e.g. an unterminated quote) must RAISE rather than silently return the line unchanged.
    Per METHOD_LINE_CANDIDATE_RE's own docstring, such a malformed candidate must fail loudly via
    this exact ValueError, not be silently skipped or passed through -- a caller that got the
    unchanged (still source-method) line back would go on to solve with the wrong method while
    downstream cache metadata records the requested one, silently mismatching them."""
    line = "    method = 'chemically-significant eigenvalues,\n"  # missing closing quote
    with pytest.raises(ValueError, match='method'):
        rewrite_arkane_method_line(line=line, method='CSE')


def test_write_arkane_network_input_file_rewrites_unspaced_method_line(tmp_path):
    """B4(b): an unspaced method='...' line (no space around '=') must still be rewritten. The
    old call-site guard `if 'method = ' in line:` never fires for this shape (it looks for the
    literal substring 'method = ' with a space), so nothing is rewritten and the ORIGINAL
    source-file method is silently re-solved under the new method's name."""
    source_text = None
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    assert "    method = 'modified strong collision',\n" in source_text
    unspaced_text = source_text.replace(
        "    method = 'modified strong collision',\n", "    method='modified strong collision',\n")

    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(unspaced_text)

    dest_path = str(tmp_path / 'input.py')
    write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')

    with open(dest_path, 'r') as f:
        dest_text = f.read()
    assert 'chemically-significant eigenvalues' in dest_text
    assert 'modified strong collision' not in dest_text


def test_write_arkane_network_input_file_raises_if_no_method_line(tmp_path):
    """B4(c): a source file with no 'method = ...' line at all must RAISE, not silently write a
    dest file whose method was never actually rewritten (and whose written-cache metadata would
    then lie about which method was solved)."""
    with open(NETWORK_FIXTURE, 'r') as f:
        source_text = f.read()
    assert "    method = 'modified strong collision',\n" in source_text
    no_method_text = source_text.replace("    method = 'modified strong collision',\n", '')

    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(no_method_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='method'):
        write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')


CLAMP_FIXTURE_TEXT = """species(
    label = 'S1',
    structure = SMILES('[CH3]'),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.0,0,0,0,0,100,5], Tmin=(100,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.0,0,0,0,0,100,5], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(100,'K'), Tmax=(3000,'K'), E0=(1,'kJ/mol'), Cp0=(20,'J/(mol*K)'), CpInf=(30,'J/(mol*K)'), label=\"\"\"S1\"\"\", comment=\"\"\"\"\"\"),
)

species(
    label = 'S2',
    structure = SMILES('[CH4]'),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.0,0,0,0,0,100,5], Tmin=(100,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.0,0,0,0,0,100,5], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(100,'K'), Tmax=(6000,'K'), E0=(1,'kJ/mol'), Cp0=(20,'J/(mol*K)'), CpInf=(30,'J/(mol*K)'), label=\"\"\"S2\"\"\", comment=\"\"\"\"\"\"),
)

network(
    label = 'PDepNetwork #1',
    isomers = [
        'S1',
    ],
    reactants = [
        ('S1', 'S2'),
    ],
    bathGas = {
        'He': 1.0,
    },
)

pressureDependence(
    label = 'PDepNetwork #1',
    Tmin = (300,'K'),
    Tmax = (3200,'K'),
    Tcount = 8,
    Pmin = (0.1,'bar'),
    Pmax = (100,'bar'),
    Pcount = 5,
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)
"""


def test_write_arkane_network_input_file_clamps_the_written_pressure_dependence_tmax_line(tmp_path):
    """The pdep block above asks for Tmax = (3200, 'K'), but S1's outer NASA thermo is only valid
    up to 3000 K (S2's is valid to 6000 K, so the network-wide ceiling -- the MIN over species --
    is 3000 K, not 6000 K). RMG tolerated this extrapolation at generation time, but standalone
    Arkane refuses with 'No valid NASA polynomial at temperature 3200 K.' The written dest file's
    own pressureDependence(...) block -- the line Arkane's job itself reads -- must be clamped to
    3000 K, not left at the network's originally requested 3200 K."""
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(CLAMP_FIXTURE_TEXT)

    dest_path = str(tmp_path / 'input.py')
    write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')

    with open(dest_path, 'r') as f:
        dest_text = f.read()
    assert "Tmax = (3200,'K')" not in dest_text
    assert "Tmax = (3000,'K')" in dest_text


def test_write_arkane_network_input_file_clamps_sensitivity_conditions_too(tmp_path):
    """The injected sensitivity_conditions directive spans the network's T/P extrema; its high-T
    entries must use the CLAMPED Tmax (3000 K), not the network's originally requested Tmax
    (3200 K) -- otherwise the sensitivity directive would itself ask Arkane to solve at a
    temperature no species' thermo supports, defeating the point of clamping the pressureDependence
    block above it."""
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(CLAMP_FIXTURE_TEXT)

    dest_path = str(tmp_path / 'input.py')
    write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE', sensitivity=True)

    with open(dest_path, 'r') as f:
        dest_text = f.read()
    assert 'sensitivity_conditions' in dest_text
    assert "(3200, 'K')" not in dest_text
    assert "(3000, 'K')" in dest_text


def test_write_arkane_network_input_file_does_not_rewrite_a_tmax_already_within_the_thermo_ceiling(tmp_path):
    """When the pdep block's requested Tmax is already at or below the network's species thermo
    ceiling, nothing needs clamping: the written file's Tmax line must be byte-for-byte the same
    as the source's, not gratuitously rewritten (which would risk introducing a formatting
    difference, e.g. int vs. float rendering, even when no clamp was actually needed)."""
    within_ceiling_text = CLAMP_FIXTURE_TEXT.replace("Tmax = (3200,'K'),", "Tmax = (2500,'K'),")
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(within_ceiling_text)

    dest_path = str(tmp_path / 'input.py')
    write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')

    with open(dest_path, 'r') as f:
        dest_text = f.read()
    assert "Tmax = (2500,'K')" in dest_text
    assert "(2500, 'K')" in dest_text  # sensitivity_conditions also uses the un-clamped value


def test_write_arkane_network_input_file_raises_when_clamping_would_leave_tmax_at_or_below_tmin(tmp_path):
    """If the species thermo ceiling (3000 K here) is at or below the pdep block's own Tmin (also
    set to 3000 K in this fixture), clamping Tmax down to that ceiling would leave no valid
    temperature range to solve at all. This must raise a ValueError naming both numbers and the
    species-thermo ceiling, rather than silently writing a degenerate (or inverted) T range."""
    degenerate_text = CLAMP_FIXTURE_TEXT.replace("Tmin = (300,'K'),", "Tmin = (3000,'K'),")
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(degenerate_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='3000'):
        write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')


def test_write_arkane_network_input_file_does_not_clamp_when_no_species_contributes_a_thermo_ceiling(tmp_path):
    """When no species(...) block in the file has a determinable outer NASA Tmax (here, both S1
    and S2 use ThermoData(...) instead of NASA(...)), network_thermo_t_max returns None, and None
    means 'clamp nothing' -- the pdep block's requested Tmax = (3200, 'K') must be written through
    unchanged, no matter how high it is, rather than being compared against a spurious ceiling
    (e.g. 0 K) that would wrongly clamp -- or crash -- every such network."""
    no_nasa_text = CLAMP_FIXTURE_TEXT.replace(
        "thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.0,0,0,0,0,100,5], Tmin=(100,'K'), "
        "Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.0,0,0,0,0,100,5], Tmin=(1000,'K'), "
        "Tmax=(3000,'K'))], Tmin=(100,'K'), Tmax=(3000,'K'), E0=(1,'kJ/mol'), "
        "Cp0=(20,'J/(mol*K)'), CpInf=(30,'J/(mol*K)'), label=\"\"\"S1\"\"\", comment=\"\"\"\"\"\"),",
        "thermo = ThermoData(Tdata=([300,400,500],'K'), Cpdata=([20.8,20.8,20.8],'J/(mol*K)'), "
        "H298=(0,'kJ/mol'), S298=(114.7,'J/(mol*K)')),",
    ).replace(
        "thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.0,0,0,0,0,100,5], Tmin=(100,'K'), "
        "Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.0,0,0,0,0,100,5], Tmin=(1000,'K'), "
        "Tmax=(6000,'K'))], Tmin=(100,'K'), Tmax=(6000,'K'), E0=(1,'kJ/mol'), "
        "Cp0=(20,'J/(mol*K)'), CpInf=(30,'J/(mol*K)'), label=\"\"\"S2\"\"\", comment=\"\"\"\"\"\"),",
        "thermo = ThermoData(Tdata=([300,400,500],'K'), Cpdata=([20.8,20.8,20.8],'J/(mol*K)'), "
        "H298=(0,'kJ/mol'), S298=(114.7,'J/(mol*K)')),",
    )
    assert 'NASA(' not in no_nasa_text
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(no_nasa_text)

    dest_path = str(tmp_path / 'input.py')
    write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE', sensitivity=True)

    with open(dest_path, 'r') as f:
        dest_text = f.read()
    assert "Tmax = (3200,'K')" in dest_text
    assert "(3200, 'K')" in dest_text


TLIST_FIXTURE_TEXT = CLAMP_FIXTURE_TEXT.replace(
    "    Tcount = 8,\n",
    "    Tcount = 8,\n"
    "    Tlist = ([3200,2290.91,1784.07,1460.87,1236.81,1072.34,946.479,847.059,766.54,700],'K'),\n",
)


def test_write_arkane_network_input_file_drops_an_out_of_range_tlist_line(tmp_path):
    """Arkane's pdep.py only calls generate_T_list() (building the grid from Tmin/Tmax/Tcount)
    when self.Tlist is None; when an explicit Tlist is present it wins outright and Tmin/Tmax are
    ignored. RMG always writes an explicit Tlist alongside Tmin/Tmax/Tcount, so clamping only the
    Tmax line (as the earlier fix did) is a no-op: Arkane still solves at the stale Tlist's 3200 K
    entry and still raises 'No valid NASA polynomial at temperature 3200 K.' The written file's
    Tlist line must be dropped entirely so Arkane regenerates the grid from the clamped
    Tmin/Tmax/Tcount, while those three lines themselves remain."""
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(TLIST_FIXTURE_TEXT)

    dest_path = str(tmp_path / 'input.py')
    write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')

    with open(dest_path, 'r') as f:
        dest_text = f.read()
    assert 'Tlist' not in dest_text
    assert "Tmin = (300,'K')" in dest_text
    assert "Tmax = (3000,'K')" in dest_text
    assert 'Tcount = 8,' in dest_text


def test_write_arkane_network_input_file_leaves_an_in_range_tlist_line_untouched(tmp_path):
    """When every Tlist entry is already within the network's species thermo ceiling, the line
    must be left byte-for-byte as written by RMG -- dropping it would be a gratuitous rewrite that
    throws away RMG's own precomputed grid for no reason, since Arkane can solve it as-is."""
    in_range_text = TLIST_FIXTURE_TEXT.replace(
        "    Tlist = ([3200,2290.91,1784.07,1460.87,1236.81,1072.34,946.479,847.059,766.54,700],'K'),\n",
        "    Tlist = ([2900,2290.91,1784.07,1460.87,1236.81,1072.34,946.479,847.059,766.54,700],'K'),\n",
    ).replace("Tmax = (3200,'K'),", "Tmax = (2900,'K'),")
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(in_range_text)

    dest_path = str(tmp_path / 'input.py')
    write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')

    with open(dest_path, 'r') as f:
        dest_text = f.read()
    assert "Tlist = ([2900,2290.91,1784.07,1460.87,1236.81,1072.34,946.479,847.059,766.54,700],'K')" in dest_text


def test_write_arkane_network_input_file_drops_tlist_even_when_tmax_already_within_ceiling(tmp_path):
    """The drop decision must be made on the Tlist entries themselves, NOT on whether the Tmax
    clamp fired: a file can carry a Tmax already at/below the ceiling while Tlist still contains a
    stale higher entry -- exactly the shape left behind by a partial or earlier clamp -- and that
    file must still be corrected, or Arkane will silently solve at the stale out-of-range Tlist
    entry despite Tmax looking fine."""
    already_clamped_text = TLIST_FIXTURE_TEXT.replace("Tmax = (3200,'K'),", "Tmax = (3000,'K'),")
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(already_clamped_text)

    dest_path = str(tmp_path / 'input.py')
    write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')

    with open(dest_path, 'r') as f:
        dest_text = f.read()
    assert 'Tlist' not in dest_text
    assert "Tmax = (3000,'K')" in dest_text


def test_write_arkane_network_input_file_raises_on_out_of_range_tlist_without_tcount(tmp_path):
    """Dropping an out-of-range Tlist is only safe if Tmin, Tmax and Tcount were all found in the
    block, since Arkane needs all three to regenerate the grid. If Tcount is missing, silently
    leaving the stale out-of-range Tlist would reproduce the very bug being fixed, so this must
    raise a ValueError naming what's missing rather than writing a file Arkane will still fail
    on."""
    no_tcount_text = TLIST_FIXTURE_TEXT.replace('    Tcount = 8,\n', '')
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(no_tcount_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='Tcount'):
        write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')


CHEBYSHEV_FIXTURE_TEXT = CLAMP_FIXTURE_TEXT.replace(
    "Tmax=(3000,'K'), E0=(1,'kJ/mol'), Cp0=(20,'J/(mol*K)'), CpInf=(30,'J/(mol*K)'), "
    "label=\"\"\"S1\"\"\"",
    "Tmax=(2000,'K'), E0=(1,'kJ/mol'), Cp0=(20,'J/(mol*K)'), CpInf=(30,'J/(mol*K)'), "
    "label=\"\"\"S1\"\"\"",
).replace(
    "Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(100,'K'), Tmax=(3000,'K'), E0=(1,'kJ/mol'), "
    "Cp0=(20,'J/(mol*K)'), CpInf=(30,'J/(mol*K)'), label=\"\"\"S1\"\"\"",
    "Tmin=(1000,'K'), Tmax=(2000,'K'))], Tmin=(100,'K'), Tmax=(2000,'K'), E0=(1,'kJ/mol'), "
    "Cp0=(20,'J/(mol*K)'), CpInf=(30,'J/(mol*K)'), label=\"\"\"S1\"\"\"",
).replace(
    "    Tcount = 8,\n",
    "    Tcount = 8,\n"
    "    Tlist = ([302.6349,324.8037,375.649,472.2362,654.3436,1016.4935,1763.5095,2928.0671],"
    "'K'),\n",
)


def test_write_arkane_network_input_file_drops_a_genuine_chebyshev_tlist(tmp_path):
    """Finding A: the be177e1 fix drops the explicit Tlist line so Arkane's own
    generate_T_list() regenerates the grid from the clamped Tmin/Tmax/Tcount. But
    generate_T_list() does NOT always build a 1/T-linspace grid: it picks between two grid
    families depending on interpolationModel --

        - Gauss-Chebyshev nodes when interpolationModel[0] == 'Chebyshev' (the branch this test
          pins), computed as
              T_i = 2 / ((1/Tmax - 1/Tmin) * (-cos((2*i+1)*pi/(2*Tcount))) + 1/Tmax + 1/Tmin)
          for i = 0 .. Tcount-1;
        - otherwise, an evenly-spaced-in-1/T grid.

    The Tlist entries below are genuine Gauss-Chebyshev nodes for Tmin=300 K, Tmax=3200 K,
    Tcount=8 (computed with the formula above), so this fixture is what RMG would actually have
    written for a Chebyshev-interpolated network -- unlike the other Tlist tests in this module,
    whose Tlist values are arbitrary placeholders that only exercise the drop/keep comparison,
    not real grid-generation equivalence. This distinction matters because 300/300 real RMG pdep
    network files sampled use interpolationModel = ('Chebyshev', ...); the 1/T-linspace branch
    (exercised elsewhere via a pdeparrhenius network) is the UNREPRESENTATIVE case, and dropping
    the Tlist is only sound if Arkane's regeneration reproduces (or improves on) the same node
    family that produced the dropped values -- which requires the interpolationModel line itself
    to survive the rewrite untouched, since that's the only thing selecting which family
    generate_T_list() picks.

    The network's thermo ceiling here is 2000 K (S1's NASA Tmax, lower than S2's 6000 K), so of
    the 8 genuine Chebyshev nodes only the top one (2928.0671 K) exceeds it -- an out-of-range
    top entry, as required -- while the rest are in range. Asserts (1) the Tlist line is dropped
    and (2) the interpolationModel line survives byte-for-byte, since that line is the only thing
    guaranteeing Arkane's regenerated grid stays in the Chebyshev family the dropped Tlist came
    from."""
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(CHEBYSHEV_FIXTURE_TEXT)

    dest_path = str(tmp_path / 'input.py')
    write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')

    with open(dest_path, 'r') as f:
        dest_text = f.read()
    assert 'Tlist' not in dest_text
    assert "interpolationModel = ('Chebyshev', 6, 4)," in dest_text
    assert "Tmin = (300,'K')" in dest_text
    assert "Tmax = (2000,'K')" in dest_text
    assert 'Tcount = 8,' in dest_text


def test_write_arkane_network_input_file_warns_when_ceiling_skips_species(tmp_path, caplog):
    """Finding B: network_thermo_t_max's ceiling silently excludes any species whose thermo it
    could not read (no thermo=, non-NASA thermo, unreadable/non-Kelvin Tmax) -- so when that
    happens, the computed ceiling may be too high and under-protect the clamp. This must not be
    silent: the writer must log a WARNING naming the skip count and the skipped species, so a
    human can tell the clamp may be under-protective, rather than trusting a ceiling that quietly
    ignored a species. Uses a fixture where S2 has ThermoData (no NASA Tmax to read) instead of
    NASA, so S2 is skipped and only S1 (2000 K) or NASA-derived ceilings drive the clamp."""
    mixed_text = CHEBYSHEV_FIXTURE_TEXT.replace(
        "thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.0,0,0,0,0,100,5], Tmin=(100,'K'), "
        "Tmax=(1000,'K')), NASAPolynomial(coeffs=[4.0,0,0,0,0,100,5], Tmin=(1000,'K'), "
        "Tmax=(6000,'K'))], Tmin=(100,'K'), Tmax=(6000,'K'), E0=(1,'kJ/mol'), "
        "Cp0=(20,'J/(mol*K)'), CpInf=(30,'J/(mol*K)'), label=\"\"\"S2\"\"\", comment=\"\"\"\"\"\"),",
        "thermo = ThermoData(Tdata=([300,400,500],'K'), Cpdata=([20.8,20.8,20.8],'J/(mol*K)'), "
        "H298=(0,'kJ/mol'), S298=(114.7,'J/(mol*K)')),",
    )
    assert 'NASA(' in mixed_text  # S1 still uses NASA
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(mixed_text)

    dest_path = str(tmp_path / 'input.py')
    import logging
    with caplog.at_level(logging.WARNING):
        write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')

    assert any("'S2'" in record.message and 'could not account for' in record.message
               for record in caplog.records)


def test_write_arkane_network_input_file_raises_on_tlist_with_non_kelvin_unit(tmp_path):
    """Finding D: the Tlist-drop logic must validate the parsed literal's shape before comparing
    entries against the thermo ceiling. A Tlist unit other than 'K' would silently mis-compare
    (e.g. Celsius or Rankine magnitudes treated as Kelvin), so rather than crash opaquely or
    compare wrongly, this must raise a ValueError naming the file and what wasn't understood,
    since a clamp IS in play here (S1's ceiling is 2000 K, below the pdep block's Tmax)."""
    bad_unit_text = CHEBYSHEV_FIXTURE_TEXT.replace(
        "'K'),\n    Pmin", "'degC'),\n    Pmin")
    assert "'degC'" in bad_unit_text
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(bad_unit_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='Tlist'):
        write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')


def test_write_arkane_network_input_file_raises_on_tlist_with_non_2_tuple_shape(tmp_path):
    """Finding D: a Tlist RHS that literal-evaluates to a 2-tuple whose first element is not a
    sequence of numbers (here, a single scalar instead of a list) must not be indexed/iterated
    blindly -- that would either crash with a TypeError deep inside the comparison or (for a
    scalar that happens to be iterable-like) silently misread the grid. A clamp is in play (S1's
    ceiling is 2000 K), so this must raise a ValueError naming the file and what wasn't
    understood, rather than crashing opaquely."""
    bad_shape_text = CHEBYSHEV_FIXTURE_TEXT.replace(
        "Tlist = ([302.6349,324.8037,375.649,472.2362,654.3436,1016.4935,1763.5095,2928.0671],"
        "'K'),\n",
        "Tlist = (300.0,'K'),\n",
    )
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(bad_shape_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='Tlist'):
        write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')


def test_write_arkane_network_input_file_raises_on_tlist_with_non_numeric_entries(tmp_path):
    """Finding D: a Tlist sequence containing a non-numeric entry (e.g. a stray string) must not
    be compared against the numeric thermo ceiling with a bare '>' -- that would raise an opaque
    TypeError deep inside the comparison instead of a clear, file-naming ValueError. A clamp is
    in play (S1's ceiling is 2000 K), so this must raise a ValueError naming the file and what
    wasn't understood."""
    bad_entry_text = CHEBYSHEV_FIXTURE_TEXT.replace(
        "302.6349,324.8037", "'not-a-number',324.8037")
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(bad_entry_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='Tlist'):
        write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')


def test_write_arkane_network_input_file_raises_on_multiline_tlist(tmp_path):
    """Finding D: the Tlist-drop logic's line-scan parsing assumes a single-line Tlist RHS. A
    multiline Tlist (0/300 real files sampled, but not impossible) must be refused clearly with a
    ValueError naming the file, rather than the line-scan silently parsing only the first line
    (producing a truncated, wrong literal) or crashing opaquely elsewhere. A clamp is in play
    (S1's ceiling is 2000 K)."""
    multiline_text = CHEBYSHEV_FIXTURE_TEXT.replace(
        "    Tlist = ([302.6349,324.8037,375.649,472.2362,654.3436,1016.4935,1763.5095,2928.0671],"
        "'K'),\n",
        "    Tlist = ([302.6349,324.8037,375.649,472.2362,654.3436,1016.4935,1763.5095,\n"
        "        2928.0671],'K'),\n",
    )
    source_path = str(tmp_path / 'source_network.py')
    with open(source_path, 'w') as f:
        f.write(multiline_text)

    dest_path = str(tmp_path / 'input.py')
    with pytest.raises(ValueError, match='Tlist'):
        write_arkane_network_input_file(source_path=source_path, dest_path=dest_path, method='CSE')


def teardown_module():
    """Safety net: remove any generated input files left under tests/data."""
    for name in ('minimal_input.py',
                 'test_write_rmg_input_file_liquid.py',
                 'test_write_rmg_input_file_seed_rads.py'):
        path = os.path.join(TEST_DATA_BASE_PATH, name)
        if os.path.isfile(path):
            os.remove(path)
