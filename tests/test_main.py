"""
t3 tests test_tandem module
"""

import os
import shutil

from rmgpy import settings
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import NASAPolynomial, NASA, ThermoData, Wilhoit
from rmgpy.thermo.model import HeatCapacityModel

import tandem.main as main


BASE_PATH = os.path.join(os.path.dirname(__file__), 'data')

spc1 = Species().from_smiles('CC')
spc1.thermo = HeatCapacityModel()
spc1.thermo.comment = 'Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH)' \
                      ' + group(Cs-CsHHH) + radical(RCCJ)'

spc2 = Species().from_smiles('CC')
spc2.thermo = HeatCapacityModel()
spc2.thermo.comment = 'Thermo library: primaryThermoLibrary + radical(CH3)'

spc3 = Species().from_smiles('CCO')
spc3.thermo = HeatCapacityModel()
spc3.thermo.comment = 'Thermo library: primaryThermoLibrary'

# arc_input_file_path = os.path.join(BASE_PATH, 'tandem_1.yml')
# rmg_input_file_path = os.path.join(BASE_PATH, 'rmg', 'input.py')
# 
# spc1 = Species().from_smiles('CC')
# spc1.thermo = HeatCapacityModel()
# spc1.thermo.comment = 'Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) ' \
#                       '+ group(Cs-CsHHH) + radical(RCCJ)'
# 
# spc2 = Species().from_smiles('CC')
# spc2.thermo = HeatCapacityModel()
# spc2.thermo.comment = 'Thermo library: primaryThermoLibrary + radical(CH3)'
# 
# spc3 = Species().from_smiles('CCO')
# spc3.thermo = HeatCapacityModel()
# spc3.thermo.comment = 'Thermo library: primaryThermoLibrary'
# 



def test_run_rmg():
    """Test the ability to run RMG from T3"""
    run_rmg_path = os.path.join(BASE_PATH, 'rmg', 'run_rmg')
    if os.path.isdir(run_rmg_path):
        shutil.rmtree(run_rmg_path)
    os.makedirs(run_rmg_path)
    rmg_input_file_path = os.path.join(run_rmg_path, 'input.py')
    shutil.copyfile(src=os.path.join(BASE_PATH, 'rmg', 'input.py'), dst=rmg_input_file_path)
    main.run_rmg(input_file=rmg_input_file_path,
                 output_directory=run_rmg_path,
                 kwargs={'restart': '',
                         'walltime': '00:00:00:00',
                         'maxproc': 1,
                         'kineticsdatastore': False},
                 arguments={'max RMG walltime': '00:00:01:00',
                            'max RMG exceptions allowed': 0},
                 tolerance=0.01,
                 thermo_library=None,
                 verbose=False)
    with open(os.path.join(run_rmg_path, 'RMG.log'), 'r') as f:
        line = f.readline()
    assert 'RMG execution initiated' in line
    shutil.rmtree(run_rmg_path)


def test_run_arc():
    """Test the ability to run ARC from T3"""
    main.run_arc(input_dict={'level_of_theory': 'b3lyp/6-31g**//b3lyp/6-31g**',
                             'calc_freq_factor': False},
                 run_directory=BASE_PATH,
                 species_to_calc=list(),
                 verbose=False)
    with open(os.path.join(BASE_PATH, 'ARC', 'arc.log'), 'r') as f:
        line = f.readline()
    assert 'ARC execution initiated' in line
    shutil.rmtree(os.path.join(BASE_PATH, 'ARC'))


def test_parse_arc_input_file():
    """Test parsing an ARC input file and extracting T3 parameters"""
    arguments, input_dict = main.parse_arc_input_file(os.path.join(BASE_PATH, 'tandem_1.yml'),
                                                      has_sa=True,
                                                      has_pdep=False,
                                                      verbose=False)
    expected_arguments = {'SA observables': [{'label': 'OH', 'smiles': '[OH]'}],
                          'SA method': 'RMG',
                          'SA threshold': 0.001,
                          'SA species': 10,
                          'SA reactions': 10,
                          'SA pdep threshold': 0.1,
                          'collision violators': True,
                          'all core species': True,
                          'RMG tolerances': [0.1, 0.01],
                          'max tandem iterations': 10,
                          'max RMG exceptions allowed': 10,
                          'max RMG walltime': '01:00:00:00'}
    expected_input_dict = {'level_of_theory': 'b3lyp/6-31g**//b3lyp/6-31g**',
                           'job_types': {'rotors': False, 'conformers': True, 'fine': False, 'freq': True,
                                         'opt': True, 'sp': True, 'onedmin': False, 'orbitals': False},
                           'allow_nonisomorphic_2d': True}
    assert arguments == expected_arguments
    assert input_dict == expected_input_dict


def test_set_legal_species_labels():
    """Test setting legal species labels"""
    # test two species with the same formula
    species_to_calc = [Species(label='i-C3H7', smiles='C[CH]C'),
                       Species(label='n-C3H7', smiles='[CH2]CC')]
    updated_species_to_calc = main.set_legal_species_labels(species_to_calc=species_to_calc, all_species=list())
    updated_labels = [spc.label for spc in updated_species_to_calc]
    assert updated_labels == ['C3H7_0', 'C3H7_1']

    # test having a species with the same formula in all_species
    species_to_calc = [Species(label='C4H9a', smiles='[CH2]CCC'),
                       Species(label='C4H9b', smiles='C[CH]CC'),
                       Species(label='NH3', smiles='N')]
    all_species = [Species(label='C4H9_0', smiles='C[C](C)(C)'),
                   Species(label='H2O', smiles='O')]
    updated_species_to_calc = main.set_legal_species_labels(species_to_calc=species_to_calc, all_species=all_species)
    updated_labels = [spc.label for spc in updated_species_to_calc]
    assert updated_labels == ['C4H9_1', 'C4H9_2', 'H3N_0']


def test_get_species_label_by_structure():
    """Test getting the species label from a list by its structure"""
    adj = """1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}"""
    spc_list = [Species(label='propane').from_smiles('CCC'), Species(label='ethane').from_smiles('CC')]
    label = main.get_species_label_by_structure(adj=adj, species_list=spc_list)
    assert label == 'ethane'


def test_get_species_by_label():
    """Test getting a species from a list by its label"""
    species_to_calc = [Species(label='C4H9a', smiles='[CH2]CCC'),
                       Species(label='C4H9b', smiles='C[CH]CC'),
                       Species(label='NH3', smiles='N')]
    spc = main.get_species_by_label(label='C4H9a', species_list=species_to_calc)
    print(type(spc))
    assert isinstance(spc, Species)
    print(spc.label)
    assert spc.label == 'C4H9a'

    spc = main.get_species_by_label(label='x', species_list=species_to_calc)
    print(spc)
    assert spc is None


def test_species_not_in_list():
    """Test determining whether a species is NOT in a list of species"""
    # calling set_legal_species_labels() sets the global species_labels_dict
    species_to_calc = [Species(label='C4H9a', smiles='[CH2]CCC'),
                       Species(label='C4H9b', smiles='C[CH]CC'),
                       Species(label='NH3', smiles='N')]
    updated_species_to_calc = main.set_legal_species_labels(species_to_calc=species_to_calc, all_species=list())
    is_species_in_list = main.species_not_in_list('C4H9_1', updated_species_to_calc)
    assert not is_species_in_list
    is_species_in_list = main.species_not_in_list('C4H9_5', updated_species_to_calc)
    assert is_species_in_list


def test_get_reaction_by_index():
    """Test getting a reaction from a list by its index"""
    reaction_list = [Reaction(index=0,
                              reactants=[Species().from_smiles('[CH2]CC')],
                              products=[Species().from_smiles('C[CH]C')])]
    rxn = main.get_reaction_by_index(index=0, reaction_list=reaction_list)
    assert isinstance(rxn, Reaction)
    assert rxn.index == 0
    assert str(rxn) == '[CH2]CC <=> C[CH]C'

    rxn = main.get_reaction_by_index(index=5, reaction_list=reaction_list)
    assert rxn is None


def test_calc_based_on_thermo_comment():
    """Test which species are selected for calculation based on their thermo comment"""
    assert main.calc_based_on_thermo_comment(spc1)
    assert main.calc_based_on_thermo_comment(spc2)
    assert not main.calc_based_on_thermo_comment(spc3)


def test_has_high_uncertainty():
    """Test determining whether a species thermo should be calculated"""
    should_species_be_calculated = main.has_high_uncertainty(
        species=spc1, unconverged_species=list(), species_to_calc=dict())
    assert should_species_be_calculated

    should_species_be_calculated = main.has_high_uncertainty(
        species=spc1, unconverged_species=[spc1.copy()], species_to_calc=dict())
    assert not should_species_be_calculated

    should_species_be_calculated = main.has_high_uncertainty(
        species=spc1, unconverged_species=list(), species_to_calc={'label': {'spc': spc2}})
    assert not should_species_be_calculated

    should_species_be_calculated = main.has_high_uncertainty(
        species=spc1, unconverged_species=list(), species_to_calc={'label': {'spc': spc3}})
    assert should_species_be_calculated


def test_load_species_and_reactions_from_chemkin_file():
    """Test loading species and reactions from a Chemkin file"""
    run_directory = os.path.join(BASE_PATH, 'iteration_0')
    rmg_species, rmg_reactions = main.load_species_and_reactions_from_chemkin_file(run_directory, verbose=False)
    assert len(rmg_species) == 27
    assert len(rmg_reactions) == 227
    assert isinstance(rmg_species[0], Species)
    assert isinstance(rmg_reactions[0], Reaction)

# def test_determine_species_to_calculate():
#     """Test that we correctly determine the species to be calculated from an RMG job"""
#     # TODO: add an actual case with SA and PDep and coll violators
#     pass


def test_determine_species_based_on_sensitivity():
    """Test determining species to calculate based on sensitivity analysis"""
    run_directory = os.path.join(BASE_PATH, 'iteration_1')
    arguments = main.parse_arc_input_file(input_file_path=os.path.join(BASE_PATH, 'tandem_1.yml'),
                                          has_sa=True,
                                          has_pdep=False,
                                          verbose=False)[0]
    rmg_species, rmg_reactions = main.load_species_and_reactions_from_chemkin_file(run_directory, verbose=False)
    unconverged_species = list()
    species_to_calc = main.determine_species_based_on_sensitivity(run_directory,
                                                                  arguments,
                                                                  rmg_species,
                                                                  rmg_reactions,
                                                                  unconverged_species,
                                                                  iteration=1,
                                                                  executed_networks=list(),
                                                                  verbose=False)[0]
    assert len(list(species_to_calc.values())) == 7
    species_to_calc_str = main.dict_to_str(species_to_calc)
    expected_species_to_calc_str = """ethane(1):
  spc: ethane(1)
  reason: observable
C2H4(17):
  spc: C=C(17)
  reason: (iteration 1) participates in the 2nd most sensitive reaction for ethane(1): C=C(17) + H(5) <=> C[CH2](4)
C2H5(4):
  spc: C[CH2](4)
  reason: (iteration 1) participates in the 2nd most sensitive reaction for ethane(1): C=C(17) + H(5) <=> C[CH2](4)
HO2(8):
  spc: [O]O(8)
  reason: (iteration 1) participates in the 5th most sensitive reaction for ethane(1): H(5) + O2(2) <=> [O]O(8)
C2H4O(41):
  spc: C=CO(41)
  reason: (iteration 1) participates in the 6th most sensitive reaction for ethane(1): C=C(17) + O(T)(14) <=> C=CO(41)
OO(34):
  spc: OO(34)
  reason: (iteration 1) participates in the 8th most sensitive reaction for ethane(1): OH(D)(33) + OH(D)(33) <=> OO(34)
CH3(3):
  spc: [CH3](3)
  reason: (iteration 1) the 5th most sensitive species thermo for ethane(1)
"""
    assert species_to_calc_str == expected_species_to_calc_str


#     def test_determine_species_from_pdep_network():
#         """"""
# use ests/data/tandem/sa_coefficients.yml with '+' in species names
#         ch2_adj = """1 C u0 p1 c0 {2,S} {3,S}
# 2 H u0 p0 c0 {1,S}
# 3 H u0 p0 c0 {1,S}"""
#         rmg_species = [Species(label='CO[O](9)').from_smiles('CO[O]'),
#                        Species(label='[CH2]OO(10)').from_smiles('[CH2]OO'),
#                        Species(label='O2(2)').from_smiles('[O][O]'),
#                        Species(label='CH3_0(5)').from_smiles('[CH3]'),
#                        Species(label='OH(D)(27)').from_smiles('[OH]'),
#                        Species(label='C=O(26)').from_smiles('C=O'),
#                        Species(label='[O]CO(28)').from_smiles('[O]CO'),
#                        Species(label='[O]O(8)').from_smiles('[O]O'),
#                        Species(label='CH2(S)(3)').from_adjacency_list(ch2_adj),
#                        Species(label='N2').from_smiles('N#N'),
#                        Species(label='Ar').from_smiles('[Ar]'),
#                        Species(label='He').from_smiles('[He]'),
#                        Species(label='Ne').from_smiles('[Ne]')]
#         pdep_rxn = PDepReaction(index=0,
#                                 reactants=[Species(label='[O]O(8)').from_smiles('[O]O'),
#                                            Species(label='CH2(S)(3)').from_adjacency_list(ch2_adj)],
#                                 products=[Species(label='CO[O](9)').from_smiles('CO[O]')],
#                                 network=PDepNetwork(index=27))
#         pdep_rxns_to_explore = [(pdep_rxn, 1, 'CH2(S)(3)')]
#         species_to_calc, executed_networks = \
#             main.determine_species_from_pdep_network(
#                 run_directory=os.path.join(BASE_PATH, 'iteration_0'),
#                 pdep_rxns_to_explore=pdep_rxns_to_explore,
#                 unconverged_species=list(),
#                 species_to_calc=dict(),
#                 iteration=1,
#                 threshold=0.1,
#                 executed_networks=list(),
#                 rmg_species=rmg_species,
#                 verbose=False)
#         # print(species_to_calc)
#         # print(executed_networks)
#         # raise


def test_modify_pdep_network_file():
    """Test modifying an Arkane P-dep network input file"""
    pdep_sa_path = os.path.join(BASE_PATH, 'iteration_0', 'pdep_sa')
    if os.path.isdir(pdep_sa_path):
        shutil.rmtree(pdep_sa_path)

    input_file_path, output_file_path, isomer_labels = \
        main.modify_pdep_network_file(run_directory=os.path.join(BASE_PATH, 'iteration_0'),
                                      network_name='network27_2',
                                      method='CSE')
    assert 'iteration_0/pdep_sa/network27_2/CSE/input.py' in input_file_path
    assert 'iteration_0/pdep_sa/network27_2/CSE/sensitivity/sa_coefficients.yml' in output_file_path
    assert isomer_labels == ('CO[O](9)', '[CH2]OO(10)')

    with open(input_file_path, 'r') as f:
        lines = f.readlines()
    sensitivity_conditions, cse = False, False
    for line in lines:
        if "    sensitivity_conditions = [[(300, 'K'), (0.01, 'bar')]," in line:
            sensitivity_conditions = True
        elif "    method = 'chemically-significant eigenvalues'," in line:
            cse = True
    assert sensitivity_conditions
    assert cse

    input_file_path, output_file_path, isomer_labels = \
        main.modify_pdep_network_file(run_directory=os.path.join(BASE_PATH, 'iteration_0'),
                                      network_name='network3_1',
                                      method='MSC')
    assert 'iteration_0/pdep_sa/network3_1/MSC/input.py' in input_file_path
    assert 'iteration_0/pdep_sa/network3_1/MSC/sensitivity/sa_coefficients.yml' in output_file_path
    assert isomer_labels == ('ethane(1)',)

    with open(input_file_path, 'r') as f:
        lines = f.readlines()
    sensitivity_conditions, msc, wrong_msc = False, False, False
    for line in lines:
        if "    sensitivity_conditions = [[(300, 'K'), (0.01, 'bar')]," in line:
            sensitivity_conditions = True
        elif "    method = 'modified strong collision'," in line:
            msc = True
        elif "modified strong collision" in line:
            # this string must not appear in any other line
            wrong_msc = True
    assert sensitivity_conditions
    assert msc
    assert not wrong_msc

    shutil.rmtree(pdep_sa_path)


def test_determine_species_based_on_collision_violators():
    """Test determining species to calculate based on collision rate violating reactions"""
    run_directory = os.path.join(BASE_PATH, 'iteration_2')
    rmg_species, rmg_reactions = main.load_species_and_reactions_from_chemkin_file(run_directory, verbose=False)
    unconverged_species = list()
    species_to_calc = main.determine_species_based_on_collision_violators(run_directory,
                                                                          rmg_species,
                                                                          unconverged_species,
                                                                          verbose=False)
    assert len(list(species_to_calc.values())) == 1
    species_to_calc_str = main.dict_to_str(species_to_calc)
    expected_species_to_calc_str = """C3H6(19):
  spc: [CH2]C[CH2](19)
  reason: species participates in a collision rate violating reaction, C3H6(19)+H(4)=C3H7_0(15)
"""
    assert species_to_calc_str == expected_species_to_calc_str


def test_add_rmg_libraries():
    """Test adding an RMG library to the RMG database repository"""
    rmg_thermo_lib_1_path = os.path.join(settings['database.directory'], 'thermo', 'libraries', 't3_thermo.py')
    rmg_thermo_lib_2_path = os.path.join(settings['database.directory'], 'thermo', 'libraries', 't3_thermo_0.py')
    libraries_path = os.path.join(BASE_PATH, 'iteration_0')
    local_context = {'ThermoData': ThermoData, 'Wilhoit': Wilhoit, 'NASAPolynomial': NASAPolynomial, 'NASA': NASA}

    if os.path.isfile(rmg_thermo_lib_1_path):
        os.remove(rmg_thermo_lib_1_path)
    if os.path.isfile(rmg_thermo_lib_2_path):
        os.remove(rmg_thermo_lib_2_path)

    # test adding a library for the fist time
    library_name = main.add_rmg_libraries(run_directory=libraries_path, library_name=None, verbose=False)
    assert library_name == 't3_thermo'
    thermo_lib = ThermoLibrary()
    thermo_lib.load(path=rmg_thermo_lib_1_path, local_context=local_context, global_context=dict())
    assert len(list(thermo_lib.entries.values())) == 10

    # test adding a library for the fist time with a different name ('t3_thermo' is occupied)
    library_name = main.add_rmg_libraries(run_directory=libraries_path, library_name=None, verbose=False)
    assert library_name == 't3_thermo_0'
    thermo_lib = ThermoLibrary()
    thermo_lib.load(path=rmg_thermo_lib_1_path, local_context=local_context, global_context=dict())
    assert len(list(thermo_lib.entries.values())) == 10

    # test appending entries to an existing library
    libraries_path = os.path.join(BASE_PATH, 'iteration_1')
    library_name = main.add_rmg_libraries(run_directory=libraries_path, library_name='t3_thermo', verbose=False)
    assert library_name == 't3_thermo'
    thermo_lib = ThermoLibrary()
    thermo_lib.load(path=rmg_thermo_lib_1_path, local_context=local_context, global_context=dict())
    assert len(list(thermo_lib.entries.values())) == 11  # extended with one additional entry

    os.remove(rmg_thermo_lib_1_path)
    os.remove(rmg_thermo_lib_2_path)


def test_get_unconverged_species():
    """Test attaining a list of unconverged species from the ARC project info file"""
    labels = ['C3H8_0', 'C3H7_1', 'C3H6_0', 'C3H6_1', 'C3H5_1',
              'C3H5_2', 'C3H6_2', 'C3H4_0', 'C3H3_0', 'C3H4_1', 'C3H4_2']
    all_species = [Species(label=label) for label in labels]
    unconverged_species = main.get_unconverged_species(run_directory=os.path.join(BASE_PATH, 'iteration_0'),
                                                       all_species=all_species,
                                                       log_species=False,
                                                       verbose=False)
    assert len(unconverged_species) == 1
    assert unconverged_species[0].label == 'C3H4_1'
