#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_main module
"""

import datetime
import os
import shutil
import re

from rmgpy.kinetics import Arrhenius
from rmgpy.reaction import Reaction
from rmgpy.rmg.pdep import PDepNetwork, PDepReaction
from rmgpy.species import Species
from rmgpy.thermo import NASA, ThermoData

from arc.common import read_yaml_file
from arc.species import ARCSpecies

from t3.common import DATA_BASE_PATH, EXAMPLES_BASE_PATH, PROJECTS_BASE_PATH
from tests.common import run_minimal
from t3.main import (T3,
                     legalize_species_label,
                     get_reaction_by_index,
                     get_species_label_by_structure)
from t3.simulate.factory import simulate_factory
from t3.utils.writer import write_rmg_input_file


test_minimal_project_directory = os.path.join(PROJECTS_BASE_PATH, 'test_minimal_delete_after_usage')

t3_minimal = {'options': {'all_core_reactions': False,
                          'all_core_species': False,
                          'collision_violators_thermo': False,
                          'collision_violators_rates': False,
                          'external_library_path': None,
                          'fit_missing_GAV': False,
                          'flux_adapter': 'RMG',
                          'library_name': 'T3lib',
                          'max_RMG_exceptions_allowed': 10,
                          'max_RMG_walltime': '00:00:05:00',
                          'max_T3_iterations': 2,
                          'max_T3_walltime': None,
                          'max_rmg_iterations': None,
                          'max_rmg_processes': None,
                          'modify_concentration_ranges_in_reverse': False,
                          'modify_concentration_ranges_together': True,
                          'num_sa_per_concentration_range': 3,
                          'num_sa_per_pressure_range': 3,
                          'num_sa_per_temperature_range': 3,
                          'num_sa_per_volume_range': 3,
                          'profiles_adapter': 'RMG',
                          'shared_library_name': None,
                          },
              'sensitivity': {'ME_methods': ['CSE', 'MSC'],
                              'SA_threshold': 0.01,
                              'adapter': 'RMGConstantTP',
                              'atol': 1e-06,
                              'global_observables': None,
                              'pdep_SA_threshold': 0.001,
                              'rtol': 0.0001,
                              'P_list': None,
                              'T_list': None,
                              'top_SA_reactions': 10,
                              'top_SA_species': 10},
              'uncertainty': None,
              }

rmg_minimal = {'memory': None,
               'cpus': None,
               'rmg_execution_type': None,
               'database': {'kinetics_depositories': 'default',
                            'kinetics_estimator': 'rate rules',
                            'kinetics_families': 'default',
                            'kinetics_libraries': [],
                            'seed_mechanisms': [],
                            'thermo_libraries': ['primaryThermoLibrary'],
                            'transport_libraries': ['OneDMinN2',
                                                    'PrimaryTransportLibrary',
                                                    'NOx2018',
                                                    'GRI-Mech']},
               'model': {'atol': 1e-16,
                         'branching_index': None,
                         'branching_ratio_max': None,
                         'core_tolerance': [0.01, 0.001],
                         'dynamics_time_scale': None,
                         'filter_reactions': True,
                         'filter_threshold': 100000000.0,
                         'ignore_overall_flux_criterion': None,
                         'max_num_objs_per_iter': 1,
                         'max_num_species': None,
                         'maximum_edge_species': None,
                         'min_core_size_for_prune': None,
                         'min_species_exist_iterations_for_prune': None,
                         'rtol': 1e-08,
                         'terminate_at_max_objects': False,
                         'tolerance_branch_reaction_to_core': None,
                         'tolerance_interrupt_simulation': [0.01, 0.001],
                         'tolerance_keep_in_edge': None,
                         'tolerance_move_edge_reaction_to_core': None,
                         'tolerance_move_edge_reaction_to_core_interrupt': None,
                         'tolerance_move_edge_reaction_to_surface': None,
                         'tolerance_move_edge_reaction_to_surface_interrupt': None,
                         'tolerance_move_surface_reaction_to_core': None,
                         'tolerance_move_surface_species_to_core': None,
                         'tolerance_thermo_keep_species_in_edge': None},
               'options': None,
               'pdep': None,
               'reactors': [{'P': 1.0,
                             'T': 1000.0,
                             'V': None,
                             'conditions_per_iteration': 12,
                             'termination_conversion': {'H2': 0.9},
                             'termination_rate_ratio': None,
                             'termination_time': (5.0, 's'),
                             'type': 'gas batch constant T P'}],
               'species': [{'SA_observable': False,
                            'UA_observable': False,
                            'adjlist': None,
                            'balance': False,
                            'concentration': 0.67,
                            'constant': False,
                            'inchi': None,
                            'label': 'H2',
                            'observable': False,
                            'reactive': True,
                            'smiles': '[H][H]',
                            'xyz': None,
                            'seed_all_rads': None,
                            'solvent': False,
                            },
                           {'SA_observable': False,
                            'UA_observable': False,
                            'adjlist': None,
                            'balance': False,
                            'concentration': 0.33,
                            'constant': False,
                            'inchi': None,
                            'label': 'O2',
                            'observable': False,
                            'reactive': True,
                            'smiles': '[O][O]',
                            'xyz': None,
                            'seed_all_rads': None,
                            'solvent': False,
                            },
                           {'SA_observable': True,
                            'UA_observable': False,
                            'adjlist': None,
                            'balance': False,
                            'concentration': 0,
                            'constant': False,
                            'inchi': None,
                            'label': 'H',
                            'observable': False,
                            'reactive': True,
                            'smiles': '[H]',
                            'xyz': None,
                            'seed_all_rads': None,
                            'solvent': False,
                            },
                           {'SA_observable': True,
                            'UA_observable': False,
                            'adjlist': None,
                            'balance': False,
                            'concentration': 0,
                            'constant': False,
                            'inchi': None,
                            'label': 'OH',
                            'observable': False,
                            'reactive': True,
                            'smiles': '[OH]',
                            'xyz': None,
                            'seed_all_rads': None,
                            'solvent': False,
                            },],
               'species_constraints': None,
               }
rmg_minimal_defaults = rmg_minimal.copy()
rmg_minimal_defaults['options'] = {'seed_name': 'Seed',
                                   'save_edge': True,
                                   'save_html': False,
                                   'generate_seed_each_iteration': True,
                                   'save_seed_to_database': False,
                                   'units': 'si',
                                   'generate_plots': False,
                                   'save_simulation_profiles': False,
                                   'verbose_comments': False,
                                   'keep_irreversible': False,
                                   'trimolecular_product_reversible': True,
                                   'save_seed_modulus': -1
                                   }
qm_minimal = {'adapter': 'ARC',
              'job_types': {'conformers': True,
                            'fine': False,
                            'freq': True,
                            'opt': True,
                            'rotors': False,
                            'sp': True},
              'level_of_theory': 'b3lyp/6-31g(d,p)',
              'reactions': [],
              'species': [],
              }

restart_base_path = os.path.join(DATA_BASE_PATH, 'restart')
dump_species_path = os.path.join(DATA_BASE_PATH, 'test_dump_species')


def setup_module():
    """
    Setup.
    Useful for rerunning these tests after a failed test during development.
    """
    if os.path.isdir(test_minimal_project_directory):
        shutil.rmtree(test_minimal_project_directory, ignore_errors=True)


def test_args_and_attributes():
    """Test passing args and assigning attributes in T3"""
    run_minimal()
    assert os.path.isfile(os.path.join(test_minimal_project_directory, 't3.log'))
    assert not os.path.isdir(os.path.join(test_minimal_project_directory, 'log_archive'))

    t3 = run_minimal()
    assert os.path.isfile(os.path.join(test_minimal_project_directory, 't3.log'))
    assert os.path.isdir(os.path.join(test_minimal_project_directory, 'log_archive'))

    assert t3.project == 'T3_minimal_example'
    assert t3.project_directory == os.path.join(PROJECTS_BASE_PATH, 'test_minimal_delete_after_usage')
    assert t3.verbose == 10

    assert t3.rmg_exceptions_counter == 0
    assert t3.iteration == 0
    assert t3.executed_networks == list()
    assert t3.t3 == t3_minimal
    assert t3.rmg == rmg_minimal_defaults
    assert t3.qm == qm_minimal
    shutil.rmtree(test_minimal_project_directory, ignore_errors=True)


def test_as_dict():
    """Test T3.as_dict()"""
    t3 = run_minimal()
    assert t3.as_dict() == {'project': 'T3_minimal_example',
                            'project_directory': test_minimal_project_directory,
                            'qm': qm_minimal,
                            'rmg': rmg_minimal_defaults,
                            't3': t3_minimal,
                            'verbose': 10}
    shutil.rmtree(test_minimal_project_directory, ignore_errors=True)


def test_write_t3_input_file():
    """Test automatically writing a T3 input file"""
    t3 = run_minimal()
    t3.write_t3_input_file()
    assert os.path.isfile(os.path.join(test_minimal_project_directory, 'T3_auto_saved_input.yml'))
    with open(os.path.join(test_minimal_project_directory, 'T3_auto_saved_input.yml'), 'r') as f:
        assert f.readline() == 'project: T3_minimal_example\n'
    shutil.rmtree(test_minimal_project_directory, ignore_errors=True)


def test_set_paths():
    """Test updating self.paths"""
    t3 = run_minimal(iteration=1, set_paths=True)
    paths = {'ARC': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC',
             'ARC info': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC/T3_minimal_example_info.yml',
             'ARC input': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC/input.yml',
             'ARC kinetics lib': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC/output/RMG '
                                 'libraries/kinetics',
             'ARC log': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC/arc.log',
             'ARC restart': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC/restart.yml',
             'ARC thermo lib': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/ARC/output/RMG '
                               'libraries/thermo/T3_minimal_example.py',
             'PDep SA': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/PDep_SA',
             'RMG': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG',
             'RMG PDep': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/pdep',
             'RMG coll vio': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/collision_rate_violators.log',
             'RMG input': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/input.py',
             'RMG log': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/RMG.log',
             'RMG job log': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/job.log',
             'RMS': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/rms',
             'SA': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA',
             'SA input': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA/input.py',
             'SA solver': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA/solver',
             'cantera annotated': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/cantera/chem_annotated.cti',
             'chem annotated': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/chemkin/chem_annotated.inp',
             'iteration': 'T3/Projects/test_minimal_delete_after_usage/iteration_1',
             'species dict': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/chemkin/'
                             'species_dictionary.txt',
             'T3 thermo lib': 'test_minimal_delete_after_usage/Libraries/T3lib.py',
             'T3 kinetics lib': 'test_minimal_delete_after_usage/Libraries/T3',
             'shared T3 thermo lib': None,
             'shared T3 kinetics lib': None,
             }
    for key, path in t3.paths.items():
        if path is None:
            assert paths[key] is None
        else:
            assert paths[key] in path


def test_restart():
    """Test that the restart() method deduces the correct status of a project"""
    # empty folders are not saved in git, add them if they don't already exist
    empty_dirs = [os.path.join(restart_base_path, 'r0'),
                  os.path.join(restart_base_path, 'r1', 'iteration_1')]
    for empty_dir in empty_dirs:
        if not os.path.isdir(empty_dir):
            os.makedirs(empty_dir)

    # empty project directory
    # results in iteration=0, run_rmg=True
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r0'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (0, True)

    # empty 'iteration_1' folder in project directory
    # results in iteration=1, run_rmg=True
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r1'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (1, True)

    # 'iteration_2' folder with an 'RMG.log' indicating a non-converged job
    # results in iteration=2, run_rmg=True
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r2'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (2, True)

    # 'iteration_3' folder with an 'RMG.log' indicating a converged job
    # results in iteration=3, run_rmg=False
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r3'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (3, False)

    # 'iteration_4' folder with an 'RMG.log' indicating a converged job and an 'arc.log' indicating a non-converged job
    # results in iteration=4, run_rmg=False
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r4'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (4, False)

    # 'iteration_5' folder with an 'RMG.log' indicating a converged job and an 'arc.log' indicating a non-converged job
    # results in iteration=5, run_rmg=False
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r5'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (5, False)

    # 'iteration_6' folder with an 'RMG.log' indicating a converged job and an 'arc.log' indicating a non-converged job
    # and an ARC 'restart.yml' file
    # results in a complete ARC run, iteration=6+1=7, run_rmg=True
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r6'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    t3.species = {0: {'RMG label': 'Imipramine_1_peroxy',
                      'Chemkin label': 'Imipramine_1_peroxy',
                      'QM label': 'Imipramine_1_peroxy_0',
                      'object': Species(smiles='C'),
                      'reasons': ['reason'],
                      'converged': None,
                      'iteration': 2,
                      }}
    t3.dump_species_and_reactions()
    assert t3.restart() == (7, True)
    t3.process_arc_run()
    assert t3.species[0]['converged'] is True
    with open(os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'arc.log'), 'r') as f:
        lines = f.readlines()
        assert 'Starting project T3_ARC_restart_test\n' in lines
        assert 'All jobs terminated. Summary for project T3_ARC_restart_test:\n' in lines

    # 'iteration_7' folder with an 'RMG.log' indicating a converged job and an 'arc.log' indicating a converged job
    # results in iteration=7+1=8, run_rmg=True
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r7'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (8, True)

    # restore r6 log file
    with open(os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'arc.log'), 'w') as f:
        f.writelines("""Dummy ARC log file\n
Starting project T3_ARC_restart_test\n
All jobs terminated. Summary for project T3_ARC_restart_test:\n
Total execution time: 00:00:00\n
ARC execution terminated on Sun Dec  4 11:50:29 2022""")


def test_check_arc_args():
    """Test the check_arc_args() method"""
    minimal_input = os.path.join(EXAMPLES_BASE_PATH, 'minimal', 'input.yml')
    input_dict = read_yaml_file(path=minimal_input)
    input_dict['verbose'] = 10
    input_dict['project_directory'] = test_minimal_project_directory
    input_dict['qm'] = {'adapter': 'ARC',
                        'unsupported_ARC_arg': 'value',
                        'bac_type': 'm',
                        }
    t3 = T3(**input_dict)
    assert t3.qm['adapter'] == 'ARC'
    assert t3.qm['bac_type'] == 'm'
    assert 'unsupported_ARC_arg' not in t3.qm


def test_run_arc():
    """Test executing ARC"""
    t3 = run_minimal(iteration=1, set_paths=True)
    t3.run_arc(arc_kwargs=t3.qm)
    with open(t3.paths['ARC log'], 'r') as f:
        lines = f.readlines()
    for line in ['Starting project T3_minimal_example\n',
                 'Geometry optimization: b3lyp/6-31g(d,p), software: gaussian (dft)\n',
                 'All jobs terminated. Summary for project T3_minimal_example:\n',
                 'Total execution time: 00:00:00\n',
                 ]:
        assert line in lines
    assert os.path.isfile(t3.paths['ARC input'])
    shutil.rmtree(test_minimal_project_directory, ignore_errors=True)


def test_process_arc_run():
    """Tests processing an ARC run and copying over a thermo library to the RMG-database repository"""
    t3 = run_minimal(project='T3',
                     project_directory=os.path.join(DATA_BASE_PATH, 'process_arc'),
                     iteration=2,
                     set_paths=True,
                     )
    t3.species = {0: {'RMG label': 'imipramine_ol_2_ket_4',
                      'Chemkin label': 'imipramine_ol_2_ket_4',
                      'QM label': 'imipramine_ol_2_ket_4',
                      'object': Species(smiles='C'),
                      'reasons': ['reason 1', 'reason 2'],
                      'converged': None,
                      'iteration': 1},
                  1: {'RMG label': 'imipramine_ol_2_ket_5',
                      'Chemkin label': 'imipramine_ol_2_ket_5',
                      'QM label': 'imipramine_ol_2_ket_5',
                      'object': Species(smiles='CC'),
                      'reasons': ['reason 3'],
                      'converged': None,
                      'iteration': 1},
                  }
    t3.process_arc_run()
    assert t3.species[0]['converged'] is True
    assert t3.species[1]['converged'] is False
    assert os.path.isfile(t3.paths['T3 thermo lib'])
    with open(t3.paths['T3 thermo lib'], 'r') as f:
        lines = f.readlines()
    for line in ['name = "T3"\n',
                 "Species imipramine_ol_2_ket_4 (run time: 1 day, 8:24:38)\n",
                 '    label = "imipramine_ol_2_ket_4",\n',
                 "        E0 = (-171.078,'kJ/mol'),\n"]:
        assert line in lines
    os.remove(t3.paths['T3 thermo lib'])


def test_get_current_rmg_tol():
    """Test getting the correct RMG tolerances"""
    t3 = run_minimal()
    t3.rmg['model']['core_tolerance'] = [0.1, 0.05, 0.01, 0.001]
    t3.iteration = 1
    assert t3.get_current_rmg_tol() == 0.1
    t3.iteration = 2
    assert t3.get_current_rmg_tol() == 0.05
    t3.iteration = 3
    assert t3.get_current_rmg_tol() == 0.01
    t3.iteration = 4
    assert t3.get_current_rmg_tol() == 0.001
    t3.iteration = 5
    assert t3.get_current_rmg_tol() == 0.001
    t3.iteration = 6
    assert t3.get_current_rmg_tol() == 0.001
    t3.iteration = 238
    assert t3.get_current_rmg_tol() == 0.001


def test_run_rmg():
    """Test the ability to run RMG from T3"""
    t3 = run_minimal(iteration=1, set_paths=True)
    t3.rmg['rmg_execution_type'] = 'incore'
    write_rmg_input_file(
        rmg=t3.rmg,
        t3=t3.t3,
        iteration=t3.iteration,
        path=t3.paths['RMG input'],
        walltime=t3.t3['options']['max_RMG_walltime'],
    )
    t3.run_rmg()
    with open(t3.paths['RMG input'], 'r') as f:
        lines = f.readlines()
    for line in ["    thermoLibraries=['primaryThermoLibrary'],\n",
                 "simulator(atol=1e-16, rtol=1e-08, sens_atol=1e-06, sens_rtol=0.0001)\n",
                 ]:
        assert line in lines
    with open(t3.paths['RMG log'], 'r') as f:
        lines = f.readlines()
    for line in ["    thermoLibraries=['primaryThermoLibrary'],\n",
                 "simulator(atol=1e-16, rtol=1e-08, sens_atol=1e-06, sens_rtol=0.0001)\n",
                 "No collision rate violators found in the model's core.\n",
                 "MODEL GENERATION COMPLETED\n",
                 "The final model core has 13 species and 20 reactions\n",
                 ]:
        assert line in lines
    assert os.path.isfile(t3.paths['chem annotated'])
    assert os.path.isfile(t3.paths['species dict'])
    shutil.rmtree(test_minimal_project_directory, ignore_errors=True)


def test_determine_species_to_calculate():
    """Test determining the species to be calculated"""

    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'determine_species'))

    # 1. no calculations required
    t3.iteration = 1
    t3.set_paths()
    t3.t3['options']['all_core_species'] = True
    additional_calcs_required = t3.determine_species_and_reactions_to_calculate()
    assert not additional_calcs_required

    # 2. All core species
    t3.iteration = 2
    t3.set_paths()
    t3.t3['options']['all_core_species'] = True
    additional_calcs_required = t3.determine_species_and_reactions_to_calculate()
    assert additional_calcs_required
    assert len(list(t3.species.keys())) == 3
    assert all([species_dict['reasons'] == ['(i 2) All core species'] for species_dict in t3.species.values()])
    assert all([species_dict['RMG label'] in ['OH', 'HO2', 'H2O2'] for species_dict in t3.species.values()])
    assert all([species_dict['QM label'] in ['s0_OH', 's1_HO2', 's2_H2O2'] for species_dict in t3.species.values()])

    # 3. collision violators
    t3.iteration = 3
    t3.set_paths()
    t3.species = dict()
    t3.t3['options']['all_core_species'] = False
    t3.t3['options']['collision_violators_thermo'] = True
    additional_calcs_required = t3.determine_species_and_reactions_to_calculate()
    assert additional_calcs_required
    assert len(list(t3.species.keys())) == 38
    assert all(['Species participates in collision rate violating reaction:' in species_dict['reasons'][0]
                or 'Participates in a reaction for which a rate coefficient is computed' in species_dict['reasons'][0]
                for species_dict in t3.species.values() if species_dict['RMG label'] not in ['H', 'OH']])

    # 4. SA observables
    assert t3.species[0]['RMG label'] == 'H'
    assert t3.species[0]['reasons'] == \
           ['(i 3) Participates in a reaction for which a rate coefficient is computed.']
    assert t3.species[3]['RMG label'] == 'CC=[C]CCCC'
    assert t3.species[3]['reasons'] == \
           ['(i 3) Species participates in collision rate violating reaction: H(3)+C7H13(920)=C7H14(323)']
    assert t3.species[10]['RMG label'] == '[CH2]CC(=C)[C]=C'
    assert t3.species[10]['reasons'] == \
           ['(i 3) Species participates in collision rate violating reaction: C6H8(2027)=C2H4(21)+C4H4(2531)']


def test_species_requires_refinement():
    """Test properly identifying the thermo comment of a species to determine whether it requires refinement"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'determine_species'))
    spc_1 = Species(smiles='C')
    spc_1.thermo = NASA()
    spc_1.thermo.comment = 'Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-CO) + ' \
                           'missing(O2d-CO) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + ' \
                           'group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)'
    spc_2 = Species(smiles='C')
    spc_2.thermo = NASA()
    spc_2.thermo.comment = 'Thermo library: primaryThermoLibrary + radical(HOOj)'
    spc_3 = Species(smiles='C')
    spc_3.thermo = NASA()
    spc_3.thermo.comment = 'Thermo library: primaryThermoLibrary'
    spc_4 = Species(smiles='CO')
    spc_4.thermo = NASA()
    spc_4.thermo.comment = 'Thermo library corrected for liquid phase: thermo_DFT_CCSDTF12_BAC + Solvation correction ' \
                           'with water as solvent and solute estimated using Data from solute library'
    spc_5 = Species(smiles='OCO[O]')
    spc_5.thermo = NASA()
    spc_5.thermo.comment = 'Thermo library corrected for liquid phase: DFT_QCI_thermo + Solvation correction with ' \
                           'water as solvent and solute estimated using abraham(Oss- noncyclic) + abraham(OssH) + ' \
                           'nonacentered(OssH) + abraham(OssH) + nonacentered(OssH) + abraham(CssH2) + radical(ROOJ)'
    spc_6 = Species(smiles='OCO[O]')
    spc_6.thermo = NASA()
    spc_6.thermo.comment = 'Thermo library corrected for liquid phase: api_soup + radical(CCsJO) + Solvation ' \
                           'correction with water as solvent and solute estimated using abraham(Oss-noncyclic) + ' \
                           'nonacentered(OxRing) + abraham(N3t) + abraham(Css-noH) + abraham(CssH2) + ' \
                           'abraham(CssH3) + abraham(Ct)'
    spc_7 = Species(smiles='OCO[O]')
    spc_7.thermo = NASA()
    spc_7.thermo.comment = 'Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(N3t-CtCs) + ' \
                           'group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + missing(Ct-CsN3t) + ' \
                           'ring(Ethylene_oxide) + Solvation correction with water as solvent and solute estimated ' \
                           'using abraham(Oss-noncyclic) + nonacentered(OxRing) + abraham(OssH) + nonacentered(OssH) + ' \
                           'abraham(N3t) + abraham(Css-noH) + abraham(CssH) + abraham(CssH3) + abraham(Ct)'
    spc_8 = Species(label='CH4', smiles='C')
    spc_8.thermo = ThermoData()
    spc_8.thermo.comment = "Thermo library: JetSurF2.0"
    spc_9 = Species(label='CH4', smiles='C')
    spc_9.thermo = ThermoData()
    spc_9.thermo.comment = "Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)(Cds-Cds)) + " \
                           "group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + " \
                           "group(Cds-CdsHH) + group(Cds-CdsHH)"
    spc_10 = Species(label='CH4', smiles='C')
    spc_10.thermo = ThermoData()
    spc_10.thermo.comment = "Thermo library: JetSurF2.0 + radical(RCCJ)"

    assert t3.species_requires_refinement(spc_1) is True
    assert t3.species_requires_refinement(spc_2) is True
    assert t3.species_requires_refinement(spc_3) is False
    assert t3.species_requires_refinement(spc_4) is False
    assert t3.species_requires_refinement(spc_5) is False
    assert t3.species_requires_refinement(spc_6) is True
    assert t3.species_requires_refinement(spc_7) is True
    assert t3.species_requires_refinement(spc_8) is False
    assert t3.species_requires_refinement(spc_9) is True
    assert t3.species_requires_refinement(spc_10) is True


def test_reaction_requires_refinement():
    """Test properly identifying the kinetic comment of a reaction to determine whether it requires refinement"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'determine_reactions'),
                     iteration=1,
                     set_paths=True,
                     )
    reactions = t3.load_species_and_reactions_from_chemkin_file()[1]
    rxn_100_kinetic_comment = """Estimated using an average for rate rule [C/H2/NonDeC;C_rad/H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: H_Abstraction"""
    assert rxn_100_kinetic_comment == reactions[100].kinetics.comment


def test_determine_species_based_on_sa():
    """Test determining species to calculate based on sensitivity analysis"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    sa_observables = ['H2', 'OH']
    simulate_adapter = simulate_factory(simulate_method=t3.t3['sensitivity']['adapter'],
                                        t3=t3.t3,
                                        rmg=t3.rmg,
                                        paths=t3.paths,
                                        logger=t3.logger,
                                        atol=t3.rmg['model']['atol'],
                                        rtol=t3.rmg['model']['rtol'],
                                        observable_list=sa_observables,
                                        sa_atol=t3.t3['sensitivity']['atol'],
                                        sa_rtol=t3.t3['sensitivity']['rtol'],
                                        )
    simulate_adapter.simulate()
    # return the dictionary containing all SA coefficients for these species
    t3.sa_dict = simulate_adapter.get_sa_coefficients()
    species_keys = t3.determine_species_based_on_sa()
    assert species_keys == [0, 1]
    # remove directories created when performing SA
    dirs = [t3.paths['SA']]
    for dir_ in dirs:
        if os.path.isdir(dir_):
            shutil.rmtree(dir_, ignore_errors=True)
    t3_log = os.path.join(DATA_BASE_PATH, 'minimal_data', 't3.log')
    if os.path.isfile(t3_log):
        os.remove(t3_log)


def test_determine_species_from_pdep_network():
    """Test determining species from pdep network"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'pdep_network'),
                     iteration=1,
                     set_paths=True,
                     )
    t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    # focus on reaction 2 in network4_2.py, whose species correspond to the indices below
    # reactants = ['H(34)', 'C4ene(26)'],
    # products = ['C4rad(5)'],
    pdep_rxn = PDepReaction(index=1,
                            reactants=[t3.rmg_species[35],
                                       t3.rmg_species[27]],
                            products=[t3.rmg_species[6]],
                            network=PDepNetwork(index=4))
    pdep_rxns_to_explore = [(pdep_rxn, 2, t3.rmg_species[6].label)]
    species_keys = t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)
    assert len(species_keys) == 1
    shutil.rmtree(t3.paths['PDep SA'], ignore_errors=True)


def test_determine_species_based_on_collision_violators():
    """Test determining species to calculate based on collision rate violating reactions"""
    t3 = run_minimal()
    t3.paths['RMG coll vio'] = os.path.join(DATA_BASE_PATH, 'collision_rate_violators', 'collision_rate_violators.log')
    t3.paths['chem annotated'] = os.path.join(DATA_BASE_PATH, 'collision_rate_violators', 'chem_annotated.inp')
    t3.paths['species dict'] = os.path.join(DATA_BASE_PATH, 'collision_rate_violators', 'species_dictionary.txt')
    t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    species_to_calc = t3.determine_species_and_reactions_based_on_collision_violators()[0]
    assert len(species_to_calc) == 18
    expected_species_to_calc = [
        'C7H13(920)',
        'C6H9(1933)',
        'C6H8(2025)',
        'C6H9(2035)',
        'C6H8(2027)',
        'S(1752)',
        'S(11767)',
        'S(11972)',
        'S(17233)',
        'S(16488)',
        'S(16530)',
        'S(25139)',
        'S(16448)',
        'S(16972)',
        'S(13229)',
        'C6H8(8657)',
        'S(26357)',
        'S(25149)'
    ]
    expected_numeric_identifiers = [int(re.findall(r'\((\d+)\)', species)[0]) for species in expected_species_to_calc]
    numeric_identifiers = [int(re.findall(r'\((\d+)\)', t3.species[index]['Chemkin label'])[0]) for index in species_to_calc]

    # Assert that the numeric identifiers are the same
    # We do this because 'S' is a generic label for a species, and the numeric identifier is what distinguishes them
    assert expected_numeric_identifiers == numeric_identifiers


def test_trsh_rmg_tol():
    """Test troubleshooting the RMG tolerance"""
    t3 = run_minimal()
    t3.t3['options']['max_T3_iterations'] = 10

    t3.rmg['model']['core_tolerance'] = [0.1, 0.1, 0.001, 0.0001]
    t3.iteration = 1
    t3.trsh_rmg_tol()
    assert t3.rmg['model']['core_tolerance'] == [0.1, 0.05, 0.001, 0.0001]

    t3.rmg['model']['core_tolerance'] = [0.1, 0.1, 0.001, 0.0001]
    t3.iteration = 2
    t3.trsh_rmg_tol()
    assert t3.rmg['model']['core_tolerance'] == [0.1, 0.1, 0.001, 0.0001]

    t3.rmg['model']['core_tolerance'] = [0.1, 0.1, 0.001, 0.0001]
    t3.iteration = 6
    t3.trsh_rmg_tol()
    assert t3.rmg['model']['core_tolerance'] == [0.1, 0.1, 0.001, 0.0001, 0.00005, 0.00005, 0.00005]

    t3.rmg['model']['core_tolerance'] = [0.1, 0.1, 0.001, 0.0001]
    t3.iteration = 12
    t3.trsh_rmg_tol()
    assert t3.rmg['model']['core_tolerance'] == [0.1, 0.1, 0.001, 0.0001]


def test_get_species_key():
    """Test checking whether a species already exists in self.species and getting its key"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'determine_species'),
                     iteration=2,
                     set_paths=True,
                     )
    t3.t3['options']['all_core_species'] = True
    t3.determine_species_and_reactions_to_calculate()

    # 1. by species
    assert t3.get_species_key(species=Species(smiles='[OH]')) == 0
    assert t3.get_species_key(species=Species(smiles='O[O]')) == 1
    assert t3.get_species_key(species=Species(smiles='OO')) == 2
    assert t3.get_species_key(species=Species(smiles='O')) is None

    # 2. by label
    t3.species = {5: {'QM label': 'O2'}}
    key = t3.get_species_key(label='O2')
    assert key == 5


def test_load_species_and_reactions_from_chemkin_file():
    """Test loading RMG species and reactions from a Chemkin file"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'determine_species'),
                     iteration=2,
                     set_paths=True,
                     )
    rmg_species, rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    assert len(rmg_species) == 12
    assert len(rmg_reactions) == 18
    assert rmg_species[0].label == 'Ar'
    assert rmg_species[10].label == 'H2O'
    assert str(rmg_reactions[0]) == 'H(3) + H(3) <=> H2(1)'
    assert str(rmg_reactions[10]) == 'OH(4) + H2O2(9) <=> HO2(6) + H2O(7)'


def test_add_species():
    """Test adding a species to self.species and to self.qm['species']"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'determine_species'),
                     iteration=2,
                     set_paths=True,
                     )
    t3.t3['options']['all_core_species'] = True
    t3.determine_species_and_reactions_to_calculate()
    spc_1 = Species(label='OH', smiles='[OH]')
    spc_2 = Species(label='hydrazine', smiles='NN')
    spc_3 = Species(label='H2', smiles='[H][H]')
    for spc in [spc_1, spc_2, spc_3]:
        spc.thermo = ThermoData()

    assert t3.get_species_key(species=spc_1) == 0
    assert t3.species[0]['RMG label'] == 'OH'
    assert t3.species[0]['reasons'] == ['(i 2) All core species']

    t3.add_species(species=spc_1, reasons='Some other reason')
    assert t3.get_species_key(species=spc_1) == 0
    assert t3.species[0]['RMG label'] == 'OH'
    assert t3.species[0]['reasons'] == ['(i 2) All core species', 'Some other reason']

    assert t3.get_species_key(species=spc_2) is None

    t3.add_species(species=spc_2, reasons=['R1', 'R2'])
    assert t3.get_species_key(species=spc_2) == 3
    assert t3.species[3]['RMG label'] == 'hydrazine'
    assert t3.species[3]['reasons'] == ['R1', 'R2']

    h2_xyz = """H  0.0000000  0.0000000  0.3736550
H  0.0000000  0.0000000 -0.3736550"""
    for i, rmg_species in enumerate(t3.rmg['species']):
        if rmg_species['label'] == 'H2':
            rmg_species['xyz'] = [h2_xyz]
    t3.add_species(species=spc_3, reasons=['R3'])
    assert t3.get_species_key(species=spc_3) == 4
    assert t3.species[4]['RMG label'] == 'H2'
    assert t3.species[4]['reasons'] == ['R3']

    found_h2 = False
    for qm_species in t3.qm['species']:
        if qm_species.label == 's4_H2':
            found_h2 = True
            assert isinstance(qm_species, ARCSpecies)
            assert qm_species.conformers == [{'symbols': ('H', 'H'),
                                              'isotopes': (1, 1),
                                              'coords': ((0.0, 0.0, 0.373655),
                                                         (0.0, 0.0, -0.373655)),
                                              }]
    assert found_h2


def test_add_reaction():
    """Test adding a reaction to self.reactions and to self.qm['reactions']"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'determine_reactions'),
                     iteration=1,
                     set_paths=True,
                     )
    rmg_species, rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    t3.add_reaction(reaction=rmg_reactions[342], reasons='reason 1')
    t3.add_reaction(reaction=rmg_reactions[100], reasons='reason 2')
    t3.add_reaction(reaction=rmg_reactions[14], reasons=['reason 3a', 'reason 3b'])

    assert t3.get_reaction_key(reaction=rmg_reactions[342]) == 0
    assert t3.reactions[0]['RMG label'] == 's0_H + s1_CC=CCCC <=> s2_S2XC6H13'
    assert 'H(2)+C6H12(1229)<=>S2XC6H13(794)' in t3.reactions[0]['Chemkin label']
    assert t3.reactions[0]['QM label'] == 's0_H + s1_CC=CCCC <=> s2_S2XC6H13'
    assert t3.reactions[0]['SMILES label'] == '[H] + CC=CCCC <=> CC[CH]CCC'
    assert isinstance(t3.reactions[0]['object'], Reaction)
    assert t3.reactions[0]['reasons'] == ['reason 1']
    assert t3.reactions[0]['converged'] is None
    assert t3.reactions[0]['iteration'] == 1

    assert t3.get_reaction_key(reaction=rmg_reactions[100]) == 1
    assert t3.reactions[1]['RMG label'] == 's3_S2XC12H25 + s4_fuel <=> s5_S3XC12H25 + s4_fuel'
    assert 'S2XC12H25(839)+fuel(1)<=>S3XC12H25(840)+fuel(1)' in t3.reactions[1]['Chemkin label']
    assert t3.reactions[1]['QM label'] == 's3_S2XC12H25 + s4_fuel <=> s5_S3XC12H25 + s4_fuel'
    assert t3.reactions[1]['SMILES label'] == 'CC[CH]CCCCCCCCC + CCCCCCCCCCCC <=> CCC[CH]CCCCCCCC + CCCCCCCCCCCC'
    assert isinstance(t3.reactions[1]['object'], Reaction)
    assert t3.reactions[1]['reasons'] == ['reason 2']
    assert t3.reactions[1]['converged'] is None
    assert t3.reactions[1]['iteration'] == 1

    assert t3.get_reaction_key(reaction=rmg_reactions[14]) == 2
    assert t3.reactions[2]['RMG label'] == 's6_PC4H9 <=> s7_C2H4 + s8_C2H5'
    assert 'PC4H9(191)<=>C2H4(22)+C2H5(52)' in t3.reactions[2]['Chemkin label']
    assert t3.reactions[2]['QM label'] == 's6_PC4H9 <=> s7_C2H4 + s8_C2H5'
    assert t3.reactions[2]['SMILES label'] == '[CH2]CCC <=> C=C + C[CH2]'
    assert isinstance(t3.reactions[2]['object'], Reaction)
    assert t3.reactions[2]['reasons'] == ['reason 3a', 'reason 3b']
    assert t3.reactions[2]['converged'] is None
    assert t3.reactions[2]['iteration'] == 1

    # check that reactant and product labels of an RMG reaction are set correctly when adding a reaction
    rmg_rxn_1 = Reaction(label='[N-]=[N+](N=O)[O] + HON <=> [O-][N+](=N)N=O + NO',
                         reactants=[Species(label='[N-]=[N+](N=O)[O]', smiles='[N-]=[N+](N=O)[O]',
                                            thermo=ThermoData(comment='comment 1')),
                                    Species(label='HON', smiles='[N-]=[OH+]',
                                            thermo=ThermoData(comment='comment 2'))],
                         products=[Species(label='[O-][N+](=N)N=O', smiles='[O-][N+](=N)N=O',
                                           thermo=ThermoData(comment='comment 3')),
                                   Species(label='NO', smiles='[N]=O',
                                           thermo=ThermoData(comment='comment 4'))],
                         kinetics=Arrhenius(A=(1, 'cm^3/(mol*s)'), n=0, Ea=(0, 'kJ/mol'), comment='kinetic comment 0'))
    t3.add_reaction(reaction=rmg_rxn_1, reasons='reason 4')
    assert t3.get_reaction_key(reaction=rmg_rxn_1) == 3
    assert t3.reactions[3]['RMG label'] == 's9_N3O2 + s10_HON <=> s11_HN3O2 + s12_NO'
    assert t3.reactions[3]['Chemkin label'] == '! kinetic comment 0\n' + \
           'N3O2+HON<=>HN3O2+NO                                 1.000000e+00 0.000     ' + \
           '0.000    \n'
    assert t3.reactions[3]['QM label'] == 's9_N3O2 + s10_HON <=> s11_HN3O2 + s12_NO'
    assert t3.reactions[3]['SMILES label'] == '[N-]=[N+](N=O)[O] + [N-]=[OH+] <=> [O-][N+](=N)N=O + [N]=O'
    assert isinstance(t3.reactions[3]['object'], Reaction)
    assert t3.reactions[3]['reasons'] == ['reason 4']
    assert t3.reactions[3]['converged'] is None
    assert t3.reactions[3]['iteration'] == 1
    assert rmg_rxn_1.reactants[0].label == 's9_N3O2'
    assert rmg_rxn_1.reactants[1].label == 's10_HON'
    assert rmg_rxn_1.products[0].label == 's11_HN3O2'
    assert rmg_rxn_1.products[1].label == 's12_NO'


def test_dump_species():
    """Test dump species for restart purposes"""
    # create an empty `iteration_5` directory
    if not os.path.isdir(os.path.join(dump_species_path, 'iteration_5')):
        os.makedirs(os.path.join(dump_species_path, 'iteration_5'))
    t3 = T3(project='test_dump_species',
            project_directory=dump_species_path,
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    t3.species = {0: {
        'RMG label': 'Imipramine_1_peroxy',
        'Chemkin label': 'Imipramine_1_peroxy',
        'QM label': 'Imipramine_1_peroxy_0',
        'object': Species(smiles='C'),
        'reasons': ['reason'],
        'converged': None,
        'iteration': 2,
    }}
    t3.dump_species_and_reactions()
    assert os.path.isfile(os.path.join(dump_species_path, 't3.log'))
    assert os.path.isfile(os.path.join(dump_species_path, 'species.yml'))
    assert t3.restart() == (5, True)


def test_load_species():
    """Test loading the dumped species dictionary from `test_dump_species()` above"""
    t3 = T3(project='test_dump_species',
            project_directory=dump_species_path,
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    t3.load_species_and_reactions()
    assert t3.species[0]['Chemkin label'] == 'Imipramine_1_peroxy'
    assert t3.species[0]['QM label'] == 'Imipramine_1_peroxy_0'


# main functions:

def test_get_reaction_by_index():
    """Test getting reaction by index"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    rmg_species, rmg_reactions = t3.load_species_and_reactions_from_chemkin_file()
    index = 5
    reaction = get_reaction_by_index(index, rmg_reactions)
    assert reaction.reactants[0].label == 'H'
    assert reaction.reactants[1].label == '[O]O'
    assert reaction.products[0].label == 'OO'


def test_legalize_species_label():
    """Test the legalize_species_label() function"""
    species = Species(smiles='C', label='CH4')
    legalize_species_label(species=species)
    assert species.label == 'CH4'

    species = Species(smiles='C#C', label='C#C')
    legalize_species_label(species=species)
    assert species.label == 'C2H2'

    species = Species(smiles='C=CC', label='S(2398)')
    legalize_species_label(species=species)
    assert species.label == 'C3H6'


def test_get_species_label_by_structure():
    """Test getting the species label from a list by its structure"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    rmg_species = t3.load_species_and_reactions_from_chemkin_file()[0]
    adj = """1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}"""
    label = get_species_label_by_structure(adj, rmg_species)
    assert label == 'H2O'

    spc_1 = Species(label='CH2')
    adj_1 = """CH2
multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}"""
    spc_1.from_adjacency_list(adj_1)
    spc_2 = Species(label='CH2(S)')
    adj_2 = """CH2(S)
multiplicity 1
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}"""
    spc_2.from_adjacency_list(adj_2)
    species_list = [spc_1, spc_2]
    label_1 = get_species_label_by_structure(adj_1, species_list)
    assert label_1 == 'CH2'
    label_2 = get_species_label_by_structure(adj_2, species_list)
    assert label_2 == 'CH2(S)'


def test_check_overtime():
    """Test checking overtime"""
    t3 = run_minimal(project_directory=os.path.join(DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    t3.t3['options']['max_T3_walltime'] = '01:00:00:00'
    t3.t0 = datetime.datetime.today()
    assert t3.check_overtime() is False
    t3.t0 = t3.t0.replace(year=t3.t0.year - 1)
    assert t3.check_overtime() is True


def teardown_module():
    """teardown any state that was previously set up."""

    # delete log files
    for i in range(10):
        directory = os.path.join(restart_base_path, f'r{i}')
        if os.path.isdir(directory):
            log_file = os.path.join(directory, 't3.log')
            if os.path.isfile(log_file):
                os.remove(log_file)
            log_archive = os.path.join(directory, 'log_archive')
            if os.path.isdir(log_archive):
                shutil.rmtree(log_archive, ignore_errors=True)

    # other files to delete
    files = [os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'T3_ARC_restart_test.info'),
             os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'input.yml'),
             os.path.join(restart_base_path, 'r6', 'species.yml'),
             os.path.join(DATA_BASE_PATH, 'process_arc', 'species.yml'),
             ]
    for file in files:
        if os.path.isfile(file):
            os.remove(file)

    # delete folders
    for directory in [test_minimal_project_directory,
                      dump_species_path,
                      os.path.join(DATA_BASE_PATH, 'minimal_data', 'log_archive'),
                      os.path.join(DATA_BASE_PATH, 'determine_species', 'log_archive'),
                      os.path.join(DATA_BASE_PATH, 'pdep_network', 'log_archive'),
                      os.path.join(DATA_BASE_PATH, 'process_arc', 'log_archive'),
                      os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'output'),
                      os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'log_and_restart_archive'),
                      ]:
        if os.path.isdir(directory):
            shutil.rmtree(directory, ignore_errors=True)
