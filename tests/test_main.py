#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_main module
"""

import datetime
import os
import shutil
import re


from arc.common import read_yaml_file

from t3.chem import T3Reaction, T3Species, T3Status
from t3.common import TEST_DATA_BASE_PATH, EXAMPLES_BASE_PATH, PROJECTS_BASE_PATH
from tests.common import run_minimal
from t3.main import (T3,
                     auto_complete_rmg_libraries,
                     legalize_species_label,
                     get_reaction_by_index,
                     get_species_label_by_structure)
from t3.simulate.factory import simulate_factory
from t3.utils.rmg_shim import Arrhenius, PDepNetwork, PDepReaction, ThermoData
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
                              'adapter': 'CanteraConstantTP',
                              'adaptive_perturbation': False,
                              'atol': 1e-06,
                              'delta_h': 0.1,
                              'delta_k': 0.05,
                              'experimental_idt_path': None,
                              'global_observables': None,
                              'idt_criterion': 'max_dOHdt',
                              'idt_sa_method': 'brute_force',
                              'max_sa_workers': 24,
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
                             'idt_mode': 'matrix',
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
                            'role': None,
                            'equivalence_ratios': None,
                            'oxidizer_fraction': None,
                            'diluent_to_oxidizer_ratio': None,
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
                            'role': None,
                            'equivalence_ratios': None,
                            'oxidizer_fraction': None,
                            'diluent_to_oxidizer_ratio': None,
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
                            'role': None,
                            'equivalence_ratios': None,
                            'oxidizer_fraction': None,
                            'diluent_to_oxidizer_ratio': None,
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
                            'role': None,
                            'equivalence_ratios': None,
                            'oxidizer_fraction': None,
                            'diluent_to_oxidizer_ratio': None,
                            }],
               'species_constraints': None,
               }
rmg_minimal_defaults = rmg_minimal.copy()
rmg_minimal_defaults['options'] = {'seed_name': 'Seed',
                                   'save_edge': False,
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
              'job_types': {'conf_opt': True,
                            'fine': False,
                            'freq': True,
                            'opt': True,
                            'rotors': False,
                            'sp': True},
              'level_of_theory': 'b3lyp/6-31g(d,p)',
              'reactions': [],
              'species': [],
              }

restart_base_path = os.path.join(TEST_DATA_BASE_PATH, 'restart')
dump_species_path = os.path.join(TEST_DATA_BASE_PATH, 'test_dump_species')


def test_thermodata_no_duplicate_fields():
    """Test that ThermoData dataclass doesn't have duplicate field definitions in source."""
    import inspect
    src = inspect.getsource(ThermoData)
    # Count field definitions (lines that look like 'fieldname: Type = ...')
    field_lines = [line.strip().split(':')[0] for line in src.splitlines()
                   if ':' in line and '=' in line and not line.strip().startswith(('#', 'def', 'class', '"'))]
    duplicates = [name for name in field_lines if field_lines.count(name) > 1]
    assert not duplicates, f"ThermoData has duplicate field definitions: {set(duplicates)}"


def setup_module():
    """
    Setup.
    Useful for rerunning these tests after a failed test during development.
    """
    if os.path.isdir(test_minimal_project_directory):
        shutil.rmtree(test_minimal_project_directory, ignore_errors=True)


def test_args_and_attributes():
    """Test passing args and assigning attributes in T3"""
    try:
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
    finally:
        shutil.rmtree(test_minimal_project_directory, ignore_errors=True)


def test_as_dict():
    """Test T3.as_dict()"""
    try:
        t3 = run_minimal()
        assert t3.as_dict() == {'project': 'T3_minimal_example',
                                'project_directory': test_minimal_project_directory,
                                'qm': qm_minimal,
                                'rmg': rmg_minimal_defaults,
                                't3': t3_minimal,
                                'verbose': 10}
    finally:
        shutil.rmtree(test_minimal_project_directory, ignore_errors=True)


def test_write_t3_input_file():
    """Test automatically writing a T3 input file"""
    try:
        t3 = run_minimal()
        t3.write_t3_input_file()
        assert os.path.isfile(os.path.join(test_minimal_project_directory, 'T3_auto_saved_input.yml'))
        with open(os.path.join(test_minimal_project_directory, 'T3_auto_saved_input.yml'), 'r') as f:
            assert f.readline() == 'project: T3_minimal_example\n'
    finally:
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
             'figs': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/Figures',
             'SA': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA',
             'SA input': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA/input.py',
             'SA coefficients': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA/sa_coefficients.yml',
             'SA dict': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA/sa.yaml',
             'SA IDT dict': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA/sa_idt.yaml',
             'SA IDT dict top X': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA/sa_idt_top_x.yaml',
             'SA solver': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/SA/solver',
             'cantera annotated': 'T3/Projects/test_minimal_delete_after_usage/iteration_1/RMG/cantera/chem_annotated.yaml',
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
    # results in iteration=0, run_rmg=True, restart_rmg=False
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r0'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (0, True, False)

    # empty 'iteration_1' folder in project directory
    # results in iteration=1, run_rmg=True, restart_rmg=False
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r1'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (1, True, False)

    # 'iteration_2' folder with an 'RMG.log' indicating a non-converged job
    # results in iteration=2, run_rmg=False, restart_rmg=True (RMG started but didn't terminate)
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r2'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (2, False, True)

    # 'iteration_3' folder with an 'RMG.log' indicating a converged job
    # results in iteration=3, run_rmg=False, restart_rmg=False
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r3'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (3, False, False)

    # 'iteration_4' folder with an 'RMG.log' indicating a converged job and an 'arc.log' indicating a non-converged job
    # results in iteration=4, run_rmg=False, restart_rmg=False
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r4'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (4, False, False)

    # 'iteration_5' folder with an 'RMG.log' indicating a converged job and an 'arc.log' indicating a non-converged job
    # results in iteration=5, run_rmg=False, restart_rmg=False
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r5'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    assert t3.restart() == (5, False, False)

    # 'iteration_6' folder with an 'RMG.log' indicating a converged job, an 'arc.log' indicating
    # a converged ARC run, and an ARC 'restart.yml' file (left over from a previous interrupted run).
    # Because arc.log already contains the termination line, the state machine sees ARC as
    # terminated and returns (6, False, False); restart_arc is NOT triggered.
    t3 = T3(project='test_restart',
            project_directory=os.path.join(restart_base_path, 'r6'),
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    t3.species = {0: T3Species(label='Imipramine_1_peroxy',
                               qm_label='Imipramine_1_peroxy_0',
                               smiles='C',
                               reasons=['reason'],
                               t3_status=T3Status.PENDING,
                               t3_index=0,
                               created_at_iteration=2)}
    t3.dump_species_and_reactions()
    assert t3.restart() == (6, False, False)
    t3.process_arc_run()
    assert t3.species[0].is_converged is True
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
    assert t3.restart() == (7, False, False)

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
    try:
        t3 = run_minimal(iteration=1, set_paths=True)
        t3.run_arc(arc_kwargs=t3.qm)
        with open(t3.paths['ARC log'], 'r') as f:
            lines = f.readlines()
        for line in ['Starting project T3_minimal_example\n',
                     'Geometry optimization: b3lyp/6-31g(d,p), software: gaussian\n',
                     'All jobs terminated. Summary for project T3_minimal_example:\n',
                     'Total execution time: 00:00:00\n',
                     ]:
            assert line in lines
        assert os.path.isfile(t3.paths['ARC input'])
    finally:
        shutil.rmtree(test_minimal_project_directory, ignore_errors=True)


def test_process_arc_run():
    """Tests processing an ARC run and copying over a thermo library to the RMG-database repository"""
    t3 = run_minimal(project='T3',
                     project_directory=os.path.join(TEST_DATA_BASE_PATH, 'process_arc'),
                     iteration=2,
                     set_paths=True,
                     )
    t3.species = {0: T3Species(label='imipramine_ol_2_ket_4',
                               qm_label='imipramine_ol_2_ket_4',
                               smiles='C',
                               reasons=['reason 1', 'reason 2'],
                               t3_status=T3Status.PENDING,
                               t3_index=0,
                               created_at_iteration=1),
                  1: T3Species(label='imipramine_ol_2_ket_5',
                               qm_label='imipramine_ol_2_ket_5',
                               smiles='CC',
                               reasons=['reason 3'],
                               t3_status=T3Status.PENDING,
                               t3_index=1,
                               created_at_iteration=1),
                  }
    try:
        t3.process_arc_run()
        assert t3.species[0].is_converged is True
        assert t3.species[1].is_converged is False
        assert os.path.isfile(t3.paths['T3 thermo lib'])
        with open(t3.paths['T3 thermo lib'], 'r') as f:
            lines = f.readlines()
        for line in ['name = "T3"\n',
                     "Species imipramine_ol_2_ket_4 (run time: 1 day, 8:24:38)\n",
                     '    label = "imipramine_ol_2_ket_4",\n',
                     "        E0 = (-171.078,'kJ/mol'),\n"]:
            assert line in lines
    finally:
        if os.path.isfile(t3.paths['T3 thermo lib']):
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
    try:
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
                     "simulator(atol=1e-16, rtol=1e-08)\n",
                     ]:
            assert line in lines
        with open(t3.paths['RMG log'], 'r') as f:
            lines = f.readlines()
        for line in ["    thermoLibraries=['primaryThermoLibrary'],\n",
                     "simulator(atol=1e-16, rtol=1e-08)\n",
                     "MODEL GENERATION COMPLETED\n",
                     ]:
            assert line in lines
        assert os.path.isfile(t3.paths['chem annotated'])
        assert os.path.isfile(t3.paths['species dict'])
    finally:
        shutil.rmtree(test_minimal_project_directory, ignore_errors=True)


def test_determine_species_to_calculate():
    """Test determining the species to be calculated"""

    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'determine_species'))

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
    assert all([species.reasons == ['(i 2) All core species'] for species in t3.species.values()])
    assert all([species.label in ['OH[4]', 'HO2[6]', 'H2O2[9]'] for species in t3.species.values()])
    print([species.qm_label for species in t3.species.values()])
    assert all([species.qm_label in ['s0_HO', 's1_HO2', 's2_H2O2'] for species in t3.species.values()])

    # 3. collision violators
    t3.iteration = 3
    t3.set_paths()
    t3.species = dict()
    t3.t3['options']['all_core_species'] = False
    t3.t3['options']['collision_violators_thermo'] = True
    additional_calcs_required = t3.determine_species_and_reactions_to_calculate()
    assert additional_calcs_required
    assert len(list(t3.species.keys())) == 38
    assert all(['Species participates in collision rate violating reaction:' in species.reasons[0]
                or 'Participates in a reaction for which a rate coefficient is computed' in species.reasons[0]
                for species in t3.species.values() if species.label not in ['H', 'OH']])

    # 4. SA observables
    # Species labels are legalized by ARC: '(' → '[', ')' → ']'
    # Reason strings retain the original RMG format from the collision violators file.
    assert t3.species[0].label == 'H[3]'
    assert t3.species[0].reasons == \
           ['(i 3) Participates in a reaction for which a rate coefficient is computed.']
    assert t3.species[3].label == 'C7H13[920]'
    assert t3.species[3].reasons == \
           ['(i 3) Species participates in collision rate violating reaction: H(3)+C7H13(920)=C7H14(323)']
    assert t3.species[10].label == 'C6H8[2027]'
    assert t3.species[10].reasons == \
           ['(i 3) Species participates in collision rate violating reaction: C6H8(2027)=C2H4(21)+C4H4(2531)']


def test_reaction_requires_refinement():
    """Test properly identifying the kinetic comment of a reaction to determine whether it requires refinement"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'determine_reactions'),
                     iteration=1,
                     set_paths=True,
                     )
    reactions = t3.load_species_and_reactions_from_yaml_file()[1]

    # Reaction list indices are 0-based after deduplicating DUPLICATE entries.
    # Chemkin DUPLICATE reactions share an equation; only the first is kept.
    assert reactions[24].kinetics_method.value == 'Library'
    assert reactions[24].kinetics_source == 'JetSurF2.0'
    assert reactions[24].kinetics_comment.strip() == """Reaction index: Chemkin #30; RMG #4278
Library reaction: JetSurF2.0
Flux pairs: PC4H9(191), C4H8(197); CH3(23), CH4(31);"""

    assert reactions[67].kinetics_method.value == 'Library'
    assert reactions[67].kinetics_source == 'JetSurF2.0'
    assert reactions[67].kinetics_comment.strip() == """Reaction index: Chemkin #76; RMG #4237
Library reaction: JetSurF2.0
Flux pairs: IC3H7(102), C3H6(104); H(2), H2(4);"""

    assert reactions[102].kinetics_method.value == 'Rate Rules'
    assert reactions[102].kinetics_source == 'Disproportionation'
    comment_102 = '\n'.join(line.rstrip() for line in reactions[102].kinetics_comment.strip().splitlines())
    assert comment_102 == """Reaction index: Chemkin #114; RMG #6134
Template reaction: Disproportionation
Flux pairs: S(842), fuel(1); C2H5(52), C2H4(22);
Estimated from node Root_N-4R->H_4CNOS-u1_N-1R!H->O_N-4CNOS->O_Ext-4CNS-R_N-Sp-5R!H#4CCCNNNSSS_N-2R!H->S_N-5R!H->O_Sp-5CS-4CCNSS_Ext-4CNS-R
Multiplied by reaction path degeneracy 3.0"""

    assert reactions[126].kinetics_method.value == 'Library'
    assert reactions[126].kinetics_source == 'JetSurF2.0'
    assert reactions[126].kinetics_comment.strip() == """Reaction index: Chemkin #141; RMG #4977
Library reaction: JetSurF2.0
Flux pairs: fuel(1), S(838); H(2), H2(4);"""

    assert reactions[146].kinetics_method.value == 'Library'
    assert reactions[146].kinetics_source == 'JetSurF2.0'
    assert reactions[146].kinetics_comment.strip() == """Reaction index: Chemkin #163; RMG #4413
Library reaction: JetSurF2.0
Flux pairs: C5H11(428), C5H10(431); H(2), H2(4);"""


def test_determine_species_based_on_sa():
    """Test determining species to calculate based on sensitivity analysis"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    try:
        t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_yaml_file()
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
    finally:
        # remove directories created when performing SA
        shutil.rmtree(t3.paths['SA'], ignore_errors=True)
        t3_log = os.path.join(TEST_DATA_BASE_PATH, 'minimal_data', 't3.log')
        if os.path.isfile(t3_log):
            os.remove(t3_log)


def test_determine_reactions_based_on_sa_cantera():
    """Test determining reactions to calculate based on SA using the CanteraConstantTP adapter"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    try:
        t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_yaml_file()
        sa_observables = ['H2', 'OH']
        simulate_adapter = simulate_factory(simulate_method='CanteraConstantTP',
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
        t3.sa_dict = simulate_adapter.get_sa_coefficients()
        reaction_keys = t3.determine_reactions_based_on_sa()
        assert isinstance(reaction_keys, list)
        # all returned keys should correspond to reactions stored in t3.reactions
        for key in reaction_keys:
            assert key in t3.reactions
            assert isinstance(t3.reactions[key], T3Reaction)
            assert len(t3.reactions[key].reasons) > 0
    finally:
        shutil.rmtree(t3.paths['SA'], ignore_errors=True)
        t3_log = os.path.join(TEST_DATA_BASE_PATH, 'minimal_data', 't3.log')
        if os.path.isfile(t3_log):
            os.remove(t3_log)


def test_determine_reactions_based_on_sa_rmg():
    """Test determining reactions to calculate based on SA using the RMGConstantTP adapter"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    try:
        t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_yaml_file()
        sa_observables = ['H2', 'OH']
        simulate_adapter = simulate_factory(simulate_method='RMGConstantTP',
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
        t3.sa_dict = simulate_adapter.get_sa_coefficients()
        reaction_keys = t3.determine_reactions_based_on_sa()
        assert isinstance(reaction_keys, list)
        for key in reaction_keys:
            assert key in t3.reactions
            assert isinstance(t3.reactions[key], T3Reaction)
            assert len(t3.reactions[key].reasons) > 0
    finally:
        shutil.rmtree(t3.paths['SA'], ignore_errors=True)
        t3_log = os.path.join(TEST_DATA_BASE_PATH, 'minimal_data', 't3.log')
        if os.path.isfile(t3_log):
            os.remove(t3_log)


def test_determine_species_from_pdep_network():
    """Test determining species from pdep network"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'pdep_network'),
                     iteration=1,
                     set_paths=True,
                     )
    try:
        t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_yaml_file()
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

        # Also test with T3Reaction (the production path)
        t3.species = dict()
        t3_rxn = T3Reaction(
            r_species=[t3.rmg_species[35], t3.rmg_species[27]],
            p_species=[t3.rmg_species[6]],
            is_pressure_dependent=True,
            network=PDepNetwork(index=4),
        )
        pdep_rxns_to_explore = [(t3_rxn, 2, t3.rmg_species[6].label)]
        species_keys = t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)
        assert len(species_keys) == 1
    finally:
        shutil.rmtree(t3.paths['PDep SA'], ignore_errors=True)


def test_determine_species_based_on_collision_violators():
    """Test determining species to calculate based on collision rate violating reactions"""
    t3 = run_minimal()
    t3.paths['RMG coll vio'] = os.path.join(TEST_DATA_BASE_PATH, 'collision_rate_violators', 'collision_rate_violators.log')
    t3.paths['cantera annotated'] = os.path.join(TEST_DATA_BASE_PATH, 'collision_rate_violators', 'cantera', 'chem_annotated.yaml')
    t3.paths['species dict'] = os.path.join(TEST_DATA_BASE_PATH, 'collision_rate_violators', 'species_dictionary.txt')
    t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_yaml_file()
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
    # ARC legalizes '(' → '[', ')' → ']' in species labels
    numeric_identifiers = [int(re.findall(r'[\[\(](\d+)[\]\)]', t3.species[index].label)[0]) for index in species_to_calc]

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
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'determine_species'),
                     iteration=2,
                     set_paths=True,
                     )
    t3.t3['options']['all_core_species'] = True
    t3.determine_species_and_reactions_to_calculate()

    # 1. by species
    assert t3.get_species_key(species=T3Species(label='OH', smiles='[OH]')) == 0
    assert t3.get_species_key(species=T3Species(label='O2', smiles='O[O]')) == 1
    assert t3.get_species_key(species=T3Species(label='O2', smiles='OO')) == 2
    assert t3.get_species_key(species=T3Species(label='O', smiles='O')) is None

    # 2. by label
    t3.species = {5: T3Species(label='O2', smiles='O[O]', reasons=['reason'], qm_label='s5_O2')}
    key = t3.get_species_key(label='s5_O2')
    assert key == 5
    t3.species = {6: T3Species(label='O2', smiles='O[O]', reasons=['reason'], qm_label='s5_O2')}
    key = t3.get_species_key(label='O2')
    assert key == 6


def test_get_species_key_rmg_label_not_shadowed():
    """Test that get_species_key with RMG label_type can find species
    whose RMG label differs from their ARC-legalized label."""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'determine_species'),
                     iteration=2,
                     set_paths=True,
                     )
    # Add a species whose ARC label differs from its RMG label
    spc = T3Species(label='OH(4)', smiles='[OH]', reasons=['test'])
    # ARC legalizes 'OH(4)' to 'OH[4]', so spc.label == 'OH[4]'
    assert spc.label == 'OH[4]'
    t3.species = {0: spc}
    # Searching for 'OH(4)' with RMG label_type should still find the species
    # since the RMG label check should compare against the original RMG label,
    # not the legalized ARC label. Currently line 1113 is unreachable because
    # line 1109 catches all label == t3_species.label matches first.
    # With 'OH(4)' != 'OH[4]' (legalized), the generic check at line 1109 won't match,
    # but RMG check at line 1113 also compares against t3_species.label (the legalized one),
    # so neither will match. This test exposes the bug.
    key = t3.get_species_key(label='OH(4)', label_type='RMG')
    assert key == 0


def test_get_reaction_key_smiles():
    """Test that get_reaction_key with SMILES label_type finds the right reaction."""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'determine_species'),
                     iteration=2,
                     set_paths=True,
                     )
    t3.rmg_species, t3.rmg_reactions = t3.load_species_and_reactions_from_yaml_file()
    # Add a reaction to t3.reactions
    rxn = t3.rmg_reactions[0]
    t3.reactions = {0: rxn}
    smiles_label = rxn.get_reaction_smiles_label()
    # get_reaction_key with SMILES label_type should find the reaction and return its key
    key = t3.get_reaction_key(label=smiles_label, label_type='SMILES')
    assert key == 0


def test_paths_shared_lib_with_none_external_path():
    """Test that set_paths doesn't crash when shared_library_name is set
    but external_library_path is None."""
    t3 = run_minimal(iteration=1)
    t3.t3['options']['shared_library_name'] = 'my_shared_lib'
    t3.t3['options']['external_library_path'] = None
    # This should not raise a TypeError from os.path.join(None, ...)
    try:
        t3.set_paths()
    except TypeError as e:
        if 'expected str' in str(e) or 'NoneType' in str(e):
            raise AssertionError(
                f"set_paths crashed with TypeError when external_library_path is None: {e}"
            )
        raise


def test_load_species_and_reactions_from_yaml_file():
    """Test loading RMG species and reactions from a Chemkin file"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'determine_species'),
                     iteration=2,
                     set_paths=True,
                     )
    rmg_species, rmg_reactions = t3.load_species_and_reactions_from_yaml_file()
    assert len(rmg_species) == 12
    assert len(rmg_reactions) == 17
    assert rmg_species[0].label == 'Ar'
    assert rmg_species[10].label == 'H2O[7]'
    assert 'H(3) + H(3) <=> H2(1)' in str(rmg_reactions[0])
    assert 'OH(4) + H2O2(9) <=> HO2(6) + H2O(7)' in str(rmg_reactions[10])


def test_add_species():
    """Test adding a species to self.species and to self.qm['species']"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'determine_species'),
                     iteration=2,
                     set_paths=True,
                     )
    t3.t3['options']['all_core_species'] = True
    t3.determine_species_and_reactions_to_calculate()
    spc_1 = T3Species(label='OH', smiles='[OH]')
    spc_2 = T3Species(label='hydrazine', smiles='NN')
    spc_3 = T3Species(label='H2', smiles='[H][H]')
    for spc in [spc_1, spc_2, spc_3]:
        spc.thermo = ThermoData()

    assert t3.get_species_key(species=spc_1) == 0
    assert t3.species[0].label == 'OH[4]'
    assert t3.species[0].reasons == ['(i 2) All core species']

    t3.add_species(species=spc_1, reasons='Some other reason')
    assert t3.get_species_key(species=spc_1) == 0
    assert t3.species[0].label == 'OH[4]'
    assert t3.species[0].reasons == ['(i 2) All core species', 'Some other reason']

    assert t3.get_species_key(species=spc_2) is None

    t3.add_species(species=spc_2, reasons=['R1', 'R2'])
    assert t3.get_species_key(species=spc_2) == 3
    assert t3.species[3].label == 'hydrazine'
    assert t3.species[3].reasons == ['R1', 'R2']

    h2_xyz = """H  0.0000000  0.0000000  0.3736550
H  0.0000000  0.0000000 -0.3736550"""
    for i, rmg_species in enumerate(t3.rmg['species']):
        if rmg_species['label'] == 'H2':
            rmg_species['xyz'] = [h2_xyz]
    t3.add_species(species=spc_3, reasons=['R3'])
    assert t3.get_species_key(species=spc_3) == 4
    assert t3.species[4].label == 'H2'
    assert t3.species[4].reasons == ['R3']

    found_h2 = False
    for qm_species in t3.qm['species']:
        if qm_species.label == 's4_H2':
            found_h2 = True
            assert isinstance(qm_species, T3Species)
            assert qm_species.conformers == [{'symbols': ('H', 'H'),
                                              'isotopes': (1, 1),
                                              'coords': ((0.0, 0.0, 0.373655),
                                                         (0.0, 0.0, -0.373655)),
                                              }]
    assert found_h2


def test_add_reaction():
    """Test adding a reaction to self.reactions and to self.qm['reactions']"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'determine_reactions'),
                     iteration=1,
                     set_paths=True,
                     )
    rmg_species, rmg_reactions = t3.load_species_and_reactions_from_yaml_file()

    # Filter for valid/balanced reactions to use in test
    valid_reactions = []
    for r in rmg_reactions:
        # Explicitly skip known invalid reaction from test data
        if "fuel(1)" in r.label:
            continue
        try:
            if r.check_atom_balance(raise_error=False):
                valid_reactions.append(r)
        except Exception:
            pass

    assert len(valid_reactions) > 0, "No valid balanced reactions found in the test dataset!"

    # Use valid reactions for adding
    t3.add_reaction(reaction=valid_reactions[0], reasons='reason 1')
    t3.add_reaction(reaction=valid_reactions[min(100, len(valid_reactions)-1)], reasons='reason 2')
    t3.add_reaction(reaction=valid_reactions[min(14, len(valid_reactions)-1)], reasons=['reason 3a', 'reason 3b'])

    assert t3.get_reaction_key(reaction=valid_reactions[0]) == 0
    # assert t3.reactions[0].rmg_label == 's0_H + s1_CC=CCCC <=> s2_S2XC6H13'
    # The first valid reaction in the filtered set is different
    assert t3.reactions[0].label == valid_reactions[0].label
    assert t3.reactions[0].qm_label == valid_reactions[0].label

    assert isinstance(t3.reactions[0], T3Reaction)
    assert t3.reactions[0].reasons == ['reason 1']
    assert t3.reactions[0].is_converged is False
    assert t3.reactions[0].created_at_iteration == 1

    assert t3.get_reaction_key(reaction=valid_reactions[min(100, len(valid_reactions)-1)]) == 1
    assert isinstance(t3.reactions[1], T3Reaction)
    assert t3.reactions[1].reasons == ['reason 2']
    assert t3.reactions[1].is_converged is False
    assert t3.reactions[1].created_at_iteration == 1

    assert t3.get_reaction_key(reaction=valid_reactions[min(14, len(valid_reactions)-1)]) == 2
    assert isinstance(t3.reactions[2], T3Reaction)
    assert t3.reactions[2].reasons == ['reason 3a', 'reason 3b']
    assert t3.reactions[2].is_converged is False
    assert t3.reactions[2].created_at_iteration == 1

    # check that reactant and product labels of an RMG reaction are set correctly when adding a reaction
    h_species = T3Species(label='H', smiles='[H]', thermo=ThermoData(comment='comment 1'))
    ch4_species = T3Species(label='CH4', smiles='C', thermo=ThermoData(comment='comment 2'))
    ch3_species = T3Species(label='CH3', smiles='[CH3]', thermo=ThermoData(comment='comment 3'))
    h2_species = T3Species(label='H2', smiles='[H][H]', thermo=ThermoData(comment='comment 4'))

    rmg_rxn_1 = T3Reaction(label='H + CH4 <=> CH3 + H2',
                           r_species=[h_species, ch4_species],
                           p_species=[ch3_species, h2_species],
                           kinetics=Arrhenius(A=(1, 'cm^3/(mol*s)'), n=0, Ea=(0, 'kJ/mol'), comment='kinetic comment 0'))
    t3.add_reaction(reaction=rmg_rxn_1, reasons='reason 4')
    assert t3.get_reaction_key(reaction=rmg_rxn_1) == 3
    assert t3.reactions[3].rmg_label == 'H + CH4 <=> CH3 + H2'
    # qm_label is built from species labels (these species have simple labels, no RMG index)
    assert t3.reactions[3].qm_label == 'H + CH4 <=> CH3 + H2'
    assert isinstance(t3.reactions[3], T3Reaction)
    assert t3.reactions[3].reasons == ['reason 4']
    assert t3.reactions[3].is_converged is False
    assert t3.reactions[3].created_at_iteration == 1


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
    t3.species = {0: T3Species(label='Imipramine_1_peroxy',
                               qm_label='Imipramine_1_peroxy_0',
                               smiles='C',
                               reasons=['reason'],
                               t3_status=T3Status.PENDING,
                               t3_index=0,
                               created_at_iteration=2)}
    t3.dump_species_and_reactions()
    assert os.path.isfile(os.path.join(dump_species_path, 't3.log'))
    assert os.path.isfile(os.path.join(dump_species_path, 'species.yml'))
    assert t3.restart() == (5, True, False)


def test_load_species():
    """Test loading the dumped species dictionary from `test_dump_species()` above"""
    t3 = T3(project='test_dump_species',
            project_directory=dump_species_path,
            t3=t3_minimal,
            rmg=rmg_minimal,
            qm=qm_minimal,
            )
    t3.load_species_and_reactions()
    assert t3.species[0].label == 'Imipramine_1_peroxy'
    assert t3.species[0].qm_label == 'Imipramine_1_peroxy_0'


# main functions:

def test_get_reaction_by_index():
    """Test getting reaction by index"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    rmg_species, rmg_reactions = t3.load_species_and_reactions_from_yaml_file()
    index = 5
    reaction = get_reaction_by_index(index, rmg_reactions)
    # Species labels are ARC-legalized: '(' → '[', ')' → ']'
    assert reaction.r_species[0].label == 'H[3]'
    assert reaction.r_species[1].label == 'HO2[6]'
    assert reaction.p_species[0].label == 'OO[9]'


def test_legalize_species_label():
    """Test the legalize_species_label() function"""
    species = T3Species(smiles='C', label='CH4')
    legalize_species_label(species=species)
    assert species.label == 'CH4'

    species = T3Species(smiles='C#C', label='C#C')
    legalize_species_label(species=species)
    # ARC legalizes 'C#C' → 'CtC' during ARCSpecies.__init__; all chars in CtC are valid
    assert species.label == 'CtC'

    species = T3Species(smiles='C=CC', label='S(2398)')
    legalize_species_label(species=species)
    # S(2398) matches the S(...) pattern → replaced with formula, then ARC legalizes parentheses
    assert species.label == 'C3H6'


def test_get_species_label_by_structure():
    """Test getting the species label from a list by its structure"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    rmg_species = t3.load_species_and_reactions_from_yaml_file()[0]
    adj = """1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}"""
    label = get_species_label_by_structure(adj, rmg_species)
    assert label == 'H2O[7]'

    adj_1 = """CH2
multiplicity 3
1 C u2 p0 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}"""
    spc_1 = T3Species(label='CH2', adjlist=adj_1)
    adj_2 = """CH2[S]
multiplicity 1
1 C u0 p1 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}"""
    spc_2 = T3Species(label='CH2[S]', adjlist=adj_2)
    species_list = [spc_1, spc_2]
    label_1 = get_species_label_by_structure(adj_1, species_list)
    assert label_1 == 'CH2'
    label_2 = get_species_label_by_structure(adj_2, species_list)
    assert label_2 == 'CH2[S]'


def test_check_overtime():
    """Test checking overtime"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    t3.t3['options']['max_T3_walltime'] = '01:00:00:00'
    t3.t0 = datetime.datetime.today()
    assert t3.check_overtime() is False
    t3.t0 = t3.t0.replace(year=t3.t0.year - 1)
    assert t3.check_overtime() is True


def test_auto_complete_rmg_libraries():
    """Test auto completing RMG libraries"""
    t3 = run_minimal(project_directory=os.path.join(TEST_DATA_BASE_PATH, 'minimal_data'),
                     iteration=1,
                     set_paths=True,
                     )
    assert t3.rmg['database']['thermo_libraries'] == ['primaryThermoLibrary']
    assert t3.rmg['database']['kinetics_libraries'] == []
    database_1 = t3.rmg['database'].copy()
    database_1['chemistry_sets'] = None  # auto_complete_rmg_libraries() is called upon initiation, these get deleted
    database_1['use_low_credence_libraries'] = False
    database_1 = auto_complete_rmg_libraries(database_1)
    assert database_1['thermo_libraries'] == ['primaryThermoLibrary']
    assert database_1['kinetics_libraries'] == []
    assert 'chemistry_sets' not in database_1

    database_2 = t3.rmg['database'].copy()
    database_2['chemistry_sets'] = ['primary', 'nitrogen']
    database_2['use_low_credence_libraries'] = False
    database_2 = auto_complete_rmg_libraries(database_2)
    assert database_2['thermo_libraries'] == ['primaryThermoLibrary', 'BurkeH2O2', 'Spiekermann_refining_elementary_reactions',
                                              'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo', 'CBS_QB3_1dHR', 'NH3', 'NitrogenCurran',
                                              'CHON_G4', 'CN', 'NOx2018']
    assert database_2['kinetics_libraries'] == ['primaryNitrogenLibrary', 'HydrazinePDep', 'Ethylamine']
    assert database_2['seed_mechanisms'] == ['primaryH2O2']
    assert 'chemistry_sets' not in database_2

    database_3 = t3.rmg['database'].copy()
    database_3['chemistry_sets'] = ['primary', 'nitrogen']
    database_3['use_low_credence_libraries'] = True
    database_3 = auto_complete_rmg_libraries(database_3)
    assert database_3['thermo_libraries'] == ['primaryThermoLibrary', 'BurkeH2O2', 'Spiekermann_refining_elementary_reactions',
                                              'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo', 'CBS_QB3_1dHR', 'NH3', 'NitrogenCurran',
                                              'CHON_G4', 'CN', 'NOx2018', 'primaryNS', 'CHN', 'CHON', 'BurcatNS']
    assert database_3['kinetics_libraries'] == ['primaryNitrogenLibrary', 'HydrazinePDep', 'Ethylamine']
    assert database_2['seed_mechanisms'] == ['primaryH2O2']
    assert 'chemistry_sets' not in database_3


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
             os.path.join(TEST_DATA_BASE_PATH, 'process_arc', 'species.yml'),
             ]
    for file in files:
        if os.path.isfile(file):
            os.remove(file)

    # delete folders
    for directory in [
        test_minimal_project_directory,
        dump_species_path,
        os.path.join(TEST_DATA_BASE_PATH, 'minimal_data', 'log_archive'),
        os.path.join(TEST_DATA_BASE_PATH, 'determine_species', 'log_archive'),
        os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'log_archive'),
        os.path.join(TEST_DATA_BASE_PATH, 'process_arc', 'log_archive'),
        os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'output'),
        os.path.join(restart_base_path, 'r6', 'iteration_6', 'ARC', 'log_and_restart_archive'),
    ]:
        if os.path.isdir(directory):
            shutil.rmtree(directory, ignore_errors=True)
