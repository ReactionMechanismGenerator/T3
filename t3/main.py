"""
The tandem tool (T3) for iterative kinetic model generation and refinement

Todo:
    - generate n generations of reactions (circles around the reactants) and get all thermo right
    - parse errors and warnings from the ARC output/status.yml file for failed species
    - modify tolerance interrupt simulation according to tolerance move to core
    - implement Cantera and RMS SA, including global observables
    - determine whether a species should be forbidden or not
    - scan pdep networks and the core, mark non-physical species
    - utilize the uncertainty analysis script
    - write an RMG input file for the user to later use.
"""

import inspect
import os
import pandas as pd
import re
import shutil
import time
import traceback
from typing import List, Optional, Tuple, Union

from pydas.daspk import DASPKError

from arkane import Arkane
from rmgpy import settings as rmg_settings
from rmgpy.chemkin import load_chemkin_file
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.exceptions import (ChemicallySignificantEigenvaluesError,
                              ChemkinError,
                              CollisionError,
                              CoreError,
                              ILPSolutionError,
                              InputError,
                              InvalidMicrocanonicalRateError,
                              KineticsError,
                              ModifiedStrongCollisionError,
                              NetworkError,
                              PressureDependenceError,
                              ReactionError,
                              ReservoirStateError,
                              StatmechError,
                              StatmechFitError,
                              )
from rmgpy.reaction import Reaction
from rmgpy.rmg.main import initialize_log as initialize_rmg_log
from rmgpy.rmg.main import RMG
from rmgpy.rmg.pdep import PDepReaction
from rmgpy.solver.simple import SimpleReactor
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.species import Species
from rmgpy.thermo import NASAPolynomial, NASA, ThermoData, Wilhoit
from rmgpy.tools.loader import load_rmg_py_job
from rmgpy.tools.simulate import simulate

from arc.common import get_ordinal_indicator, key_by_val, read_yaml_file, save_yaml_file, time_lapse
from arc.main import ARC
from arc.settings import valid_chars

from t3.common import PROJECTS_BASE_PATH, delete_root_rmg_log
from t3.logger import Logger
from t3.schema import InputBase
from t3.utils.writer import write_pdep_network_file, write_rmg_input_file

PDEP_SA_ME_METHODS = ['CSE', 'MSC']
RMG_THERMO_LIB_BASE_PATH = os.path.join(rmg_settings['database.directory'], 'thermo', 'libraries')


class T3(object):
    """
    The main T3 class.

    Dictionary structures::

        species = {
            <int: T3_spc_index>: {
                'RMG label': <str: RMG label>,
                'Chemkin label': <str: Chemkin label>,
                'QM label': <str: The label used for the QM calc>,
                'object': <Species: RMG Species object>,
                'reasons': <List[str]: Reasons for calculating this species>,
                'converged': <Optional[bool]: whether thermo was successfully calculated, ``None`` if pending>,
                'iteration': <int: The iteration this species was originally added on>,
            },
        }

    Args:
        project (str): The project name.
        project_directory (str, optional): The project directory. Required through the API, optional through an input
                                           file (will be set to the directory of the input file if not specified).
        verbose (int, optional): The logging level, optional. 10 - debug, 20 - info, 30 - warning, default: 20.
        clean_dir (bool, optional): Whether to delete all existing files and folders in the project directory prior to
                                    execution. If set to ``True``, the restart feature will not be triggered.
                                    Default: ``False``.
        t3 (dict, optional): T3 directives.
        rmg (dict): RMG directives.
        qm (dict, optional): QM directive.

    Attributes:
        project (str): The project name.
        project_directory (str): The project directory. Required through the API, optional through an input file.
        t3 (dict): T3 directives.
        rmg (dict): RMG directives.
        qm (dict, optional): QM directive.
        verbose (int): The logging level.
        t0 (float): Initial time when the project was spawned.
        logger (Logger): The Logger class instance.
        rmg_exceptions_counter (int): Number of times RMG crashed.
        iteration (int): The current iteration number. Iteration 0 is reserved for pre-running ARC,
                         normal iteration numbers begin at 1.
        thermo_lib_base_path (str): The path to the thermo libraries folder in RMG-database.
        kinetics_lib_base_path (str): The path to the kinetics libraries folder in RMG-database.
        species (Dict[int, dict]: The T3 species dictionary. Keys are T3 species indices, values are dictionaries.
        paths (dict): Various directory and file paths.
        executed_networks (list): PDep networks for which SA was already executed. Entries are tuples of isomer labels.
        rmg_species (List[Species]): Entries are RMG species objects in the model core for a certain T3 iteration.
        rmg_reactions (List[Reaction]): Entries are RMG reaction objects in the model core for a certain T3 iteration.
    """

    def __init__(self,
                 project: str,
                 rmg: dict,
                 t3: Optional[dict] = None,
                 qm: Optional[dict] = None,
                 project_directory: Optional[str] = None,
                 verbose: int = 20,
                 clean_dir: bool = False,
                 ):

        self.t0 = time.time()  # initialize the timer

        project_directory = project_directory or os.path.join(PROJECTS_BASE_PATH, project)

        self.schema = InputBase(project=project,
                                project_directory=project_directory,
                                t3=t3,
                                rmg=rmg,
                                qm=qm,
                                verbose=verbose,
                                ).dict()

        self.schema_exclude_unset = InputBase(project=project,
                                              project_directory=project_directory,
                                              t3=t3,
                                              rmg=rmg,
                                              qm=qm,
                                              verbose=verbose,
                                              ).dict(exclude_unset=True)

        self.project = self.schema['project']
        self.project_directory = self.schema['project_directory']
        self.t3 = self.schema['t3']
        self.rmg = self.schema['rmg']
        self.qm = self.schema['qm']
        self.verbose = self.schema['verbose']

        if clean_dir and os.path.isdir(self.project_directory):
            shutil.rmtree(self.project_directory)
        if not os.path.isdir(self.project_directory):
            os.makedirs(self.project_directory)

        # initialize the logger
        self.logger = Logger(project=self.project,
                             project_directory=self.project_directory,
                             verbose=self.verbose,
                             t0=self.t0,
                             )
        self.logger.log_args(schema=self.schema_exclude_unset)

        self.rmg_exceptions_counter = 0
        self.iteration = 0
        self.thermo_lib_base_path = os.path.join(rmg_settings['database.directory'], 'thermo', 'libraries')
        self.kinetics_lib_base_path = os.path.join(rmg_settings['database.directory'], 'kinetics', 'libraries')
        self.species, self.reactions, self.paths = dict(), dict(), dict()
        self.rmg_species, self.rmg_reactions, self.executed_networks = list(), list(), list()

        if self.qm:
            # check args
            self.check_arc_args()

        if any(self.rmg['model']['core_tolerance'][i+1] >= self.rmg['model']['core_tolerance'][i]
               for i in range(len(self.rmg['model']['core_tolerance']) - 1)):
            self.logger.warning(f'The RMG tolerances are not in descending order.')
            self.logger.info(f'Got: {self.rmg["model"]["core_tolerance"]}')

    def as_dict(self) -> dict:
        """
        Generate a dictionary representation of the object's arguments.

        Returns:
            dict: The dictionary representation.
        """
        result = dict()
        result['project'] = self.project
        result['project_directory'] = self.project_directory
        result['verbose'] = self.verbose
        result['t3'] = self.t3
        result['rmg'] = self.rmg
        result['qm'] = self.qm
        return result

    def write_t3_input_file(self,
                            path: Optional[str] = None,
                            all_args: bool = False,
                            ) -> None:
        """
        Save the current **arguments** (not all attributes) as a T3 input file.
        Useful for creating an input file using the API.

        Args:
             path (str, optional): The full path for the generated input file,
                                   or to the folder where this file will be saved under a default name.
                                   If ``None``, the input file will be saved to the project directory.
             all_args (bool, optional): Whether to save all arguments in the generated input file
                                        including all default values). Default: ``False``.
        """
        if path is None:
            path = os.path.join(self.project_directory, 'T3_auto_saved_input.yml')
        if os.path.isdir(path):
            path += '/' if path[-1] != '/' else ''
            path += 'T3_auto_saved_input.yml'
        base_path = os.path.dirname(path)
        if not os.path.isdir(base_path):
            os.makedirs(base_path)
        self.logger.info(f'\n\nWriting input file to {path}')
        save_yaml_file(path=path, content=self.schema if all_args else self.schema_exclude_unset)

    def execute(self):
        """
        Execute T3.
        """
        # check whether T3 should be restarted in the project directory
        iteration_start, run_rmg_at_start = self.restart()

        if iteration_start == 0 \
                and self.qm \
                and self.qm['adapter'] == 'ARC' \
                and (len(self.qm['species']) or len(self.qm['reactions'])):
            self.set_paths(iteration=iteration_start)
            self.run_arc(arc_kwargs=self.qm)
            self.process_arc_run()
            # don't request these species and reactions again
            iteration_start += 1
        # ARC species and reactions will be loaded again if restarting and they were already sent to ARC, set to list()
        self.qm['species'], self.qm['reactions'] = list(), list()

        additional_calcs_required = False
        iteration_start = iteration_start or 1

        # main T3 loop
        max_t3_iterations = self.t3['options']['max_T3_iterations']
        for self.iteration in range(iteration_start, max_t3_iterations):

            self.logger.info(f'\n\n\nT3 iteration {self.iteration}:\n'
                             f'---------------\n')
            self.set_paths()

            # RMG
            if self.iteration > iteration_start or self.iteration == iteration_start and run_rmg_at_start:
                write_rmg_input_file(
                    kwargs=self.rmg,
                    iteration=self.iteration,
                    path=self.paths['RMG input'],
                    walltime=self.t3['options']['max_RMG_walltime'],
                )
                self.run_rmg()

            # SA
            if self.t3['sensitivity'] is not None:
                sa_success = self.run_sa()
                if not sa_success:
                    self.logger.error(f"Could not complete the sensitivity analysis using "
                                      f"{self.t3['sensitivity']['adapter']}.")
                    self.trsh_rmg_tol()

            # determine what needs to be calculated
            additional_calcs_required = self.determine_species_to_calculate()

            # ARC
            if additional_calcs_required:
                if self.qm is None:
                    self.logger.error('Could not refine the model without any QM arguments.')
                    additional_calcs_required = None
                else:
                    self.run_arc(arc_kwargs=self.qm)
                    self.process_arc_run()
            if not additional_calcs_required and self.iteration >= len(self.rmg['model']['core_tolerance']):
                # T3 iterated through all of the user requested tolerances, and there are no more calculations required
                break

        if additional_calcs_required:
            # The main T3 loop terminated, but the latest calculations were not included in the model.
            # Run RMG for the last time.
            self.iteration += 1
            self.logger.info(f'\n\n\nT3 iteration {self.iteration} (just generating a model using RMG):\n'
                             f'------------------------------------------------------\n')
            self.set_paths()
            write_rmg_input_file(
                kwargs=self.rmg,
                iteration=self.iteration,
                path=self.paths['RMG input'],
                walltime=self.t3['options']['max_RMG_walltime'],
            )
            self.run_rmg()

        self.logger.log_species_summary(species_dict=self.species)
        self.logger.log_footer()
        delete_root_rmg_log(project_directory=self.project_directory)

    def set_paths(self,
                  iteration: Optional[int] = None,
                  project_directory: Optional[str] = None,
                  ):
        """
        Set various file and folder paths (but don't create the folders).

        Args:
            iteration (Optional[int]): The iteration number. If ``None``, self.iteration will be used instead.
            project_directory (Optional[str]): The base folder to use, will use ``self.project_directory``
                                               if this argument is left empty.
        """
        iteration = iteration or self.iteration
        iteration_path = os.path.join(project_directory or self.project_directory, f'iteration_{iteration}')
        self.paths = {
            'iteration': iteration_path,
            'RMG': os.path.join(iteration_path, 'RMG'),
            'RMG PDep': os.path.join(iteration_path, 'RMG', 'pdep'),
            'RMG input': os.path.join(iteration_path, 'RMG', 'input.py'),
            'RMG log': os.path.join(iteration_path, 'RMG', 'RMG.log'),
            'RMG coll vio': os.path.join(iteration_path, 'RMG', 'collision_rate_violators.log'),
            'chem annotated': os.path.join(iteration_path, 'RMG', 'chemkin', 'chem_annotated.inp'),
            'species dict': os.path.join(iteration_path, 'RMG', 'chemkin', 'species_dictionary.txt'),
            'SA': os.path.join(iteration_path, 'SA'),
            'SA solver': os.path.join(iteration_path, 'SA', 'solver'),
            'SA input': os.path.join(iteration_path, 'SA', 'input.py'),
            'PDep SA': os.path.join(iteration_path, 'PDep_SA'),
            'ARC': os.path.join(iteration_path, 'ARC'),
            'ARC input': os.path.join(iteration_path, 'ARC', 'input.yml'),
            'ARC restart': os.path.join(iteration_path, 'ARC', 'restart.yml'),
            'ARC log': os.path.join(iteration_path, 'ARC', 'arc.log'),
            'ARC info': os.path.join(iteration_path, 'ARC',
                                     f"{self.qm['project'] if 'project' in self.qm else 'T3'}.info"),
            'ARC thermo lib': os.path.join(iteration_path, 'ARC', 'output', 'RMG libraries', 'thermo',
                                           f"{self.qm['project'] if 'project' in self.qm else 'T3'}.py"),
            'ARC kinetics lib': os.path.join(iteration_path, 'ARC', 'output', 'RMG libraries', 'kinetics'),
        }

    def restart(self) -> Tuple[int, bool]:
        """
        Restart T3 by looking for existing iteration folders.
        Restarts ARC if it ran and did not terminate.

        Returns:
            Tuple[int, bool]:
                - The current iteration number.
                - Whether to run RMG for this iteration.
        """
        # set default values
        i_max = 0
        run_rmg_i, restart_arc_i = True, False

        folders = tuple(os.walk(self.project_directory))[0][1]  # returns a 3-tuple: (dirpath, dirnames, filenames)
        iteration_folders = [folder for folder in folders if 'iteration_' in folder]

        if len(iteration_folders):
            self.load_species()
            i_max = max([int(folder.split('_')[1]) for folder in iteration_folders])  # get the latest iteration number
            self.set_paths(iteration=i_max)
            if i_max != 0 and os.path.isfile(self.paths['RMG log']):
                # iteration 0 is reserved for ARC only if needed
                with open(self.paths['RMG log'], 'r') as f:
                    lines = f.readlines()
                    for line in lines[::-1]:
                        if 'MODEL GENERATION COMPLETED' in line:
                            # RMG terminated, no need to regenerate the model
                            run_rmg_i = False
                            break
            if os.path.isfile(self.paths['ARC log']) and (not run_rmg_i or i_max == 0):
                # The ARC log file exists, and no need to run RMG (converged) or this is iteration 0
                with open(self.paths['ARC log'], 'r') as f:
                    lines = f.readlines()
                    for line in lines[::-1]:
                        if 'ARC execution terminated on' in line:
                            # ARC terminated as well, continue to the next iteration
                            i_max += 1
                            run_rmg_i = True
                            break
                    else:
                        # ARC did not terminate, see if the restart file was generated
                        if os.path.isfile(self.paths['ARC restart']):
                            restart_arc_i = True
            if i_max or not run_rmg_i or restart_arc_i:
                rmg_text = ', using the completed RMG run from this iteration' if not run_rmg_i \
                    else ', re-running RMG for this iteration'
                arc_text = ', restarting the previous ARC run in this iteration' if restart_arc_i else ''
                self.logger.log(f'\nRestarting T3 from iteration {i_max}{rmg_text}{arc_text}.\n')
            if restart_arc_i:
                self.run_arc(input_file_path=self.paths['ARC restart'])
                self.process_arc_run()
                i_max += 1
                run_rmg_i = True

        return i_max, run_rmg_i

    def check_arc_args(self):
        """
        Check that all arguments in the ARC input dictionary are legal.

        Returns:
            dict: Updated ARC arguments.
        """
        if self.qm and self.qm['adapter'] == 'ARC':
            allowed_non_arc_args = ['adapter']
            for key in list(self.qm.keys()):
                if key not in inspect.getfullargspec(ARC.__init__).args and key not in allowed_non_arc_args:
                    # This argument was not extracted above, and it's not an ARC argument, remove so ARC doesn't crush
                    self.logger.error(f'Argument "{key}" passed to ARC is not allowed. NOT using it.\n'
                                      f'(if this was meant to be a T3 argument, it was not recognized).')
                    del self.qm[key]

    def run_arc(self,
                input_file_path: Optional[str] = None,
                arc_kwargs: Optional[dict] = None,
                ):
        """
        Run ARC.
        Either ``input_file_path`` or ``arc_kwargs`` must be specified.

        Args:
            input_file_path (Optional[str]): A path to an ARC input file.
            arc_kwargs (Optional[dict]): ARC's arguments as a keyword argument dictionary.
                                         No need to remove the ``adapter`` key if present.

        Raises:
            ValueError: If neither or both ``input_file`` and ``kwargs`` were specified,
                        or an ARC input file was specified but could not be found.
            Various ARC exceptions: If ARC crashes.
        """
        if input_file_path is None and arc_kwargs is None:
            raise ValueError('Either input_file or kwargs must be given to run ARC, got neither.')
        if input_file_path is not None and arc_kwargs is not None:
            raise ValueError('Either input_file or kwargs must be given to run ARC, not both.')
        if input_file_path is not None and not os.path.isfile(input_file_path):
            raise ValueError(f'The ARC input file {input_file_path} could not be found.')

        self.logger.info('\nRunning ARC...')

        if input_file_path is not None:
            arc_kwargs = read_yaml_file(input_file_path)

        arc_kwargs = arc_kwargs.copy()
        if 'adapter' in arc_kwargs:
            del arc_kwargs['adapter']

        arc_kwargs['project_directory'] = self.paths['ARC']
        if not os.path.isdir(arc_kwargs['project_directory']):
            os.makedirs(arc_kwargs['project_directory'])

        if 'project' not in arc_kwargs:
            arc_kwargs['project'] = 'T3'
        tic = time.time()
        arc = ARC(**arc_kwargs)
        try:
            arc.execute()
        except Exception as e:
            self.logger.error(f'ARC crashed with {e.__class__}. Got the following error message:\n{e}')
            raise
        elapsed_time = time_lapse(tic)
        self.logger.info(f'ARC terminated, execution time: {elapsed_time}')

    def process_arc_run(self):
        """
        Process an ARC run.
        Sets the self.species[<key>]['converged'] parameter.

        Todo:
            - Check for non-physical species in unconverged species.
        """
        unconverged_keys, converged_keys = list(), list()
        if os.path.isfile(self.paths['ARC info']):
            with open(self.paths['ARC info'], 'r') as f:
                read = False
                for line in f:
                    if read:
                        if 'Species' in line:
                            # e.g.:
                            # "Species Imipramine_1_peroxy (run time: 1 day, 17:28:32)"
                            # "Species Imipramine_1_peroxy (Failed!) (run time: 1 day, 17:28:32)"
                            key = self.get_species_key(label=line.split()[1])
                            if key is not None:
                                if '(Failed!)' in line:
                                    unconverged_keys.append(key)
                                    self.species[key]['converged'] = False
                                else:
                                    converged_keys.append(key)
                                    self.species[key]['converged'] = True
                    if 'Considered the following species' in line:
                        read = True
                    if 'Overall time since project initiation' in line:
                        read = False
        else:
            raise ValueError(f'ARC did not save a project.info file, something must be wrong.')
        self.logger.log_unconverged_species(
            species_keys=unconverged_keys,
            species_dict=self.species,
        )
        if len(converged_keys):
            # we calculated something, add to thermo library
            self.add_to_rmg_library()
        # clear the calculated objects from self.qm:
        self.qm['species'], self.qm['reactions'] = list(), list()

    def get_current_rmg_tol(self) -> float:
        """
        Get the current RMG tolerance.

        Returns:
            float: The current RMG move to core tolerance.
        """
        # self.iteration is 1-indexed
        return self.rmg['model']['core_tolerance'][self.iteration - 1] \
            if len(self.rmg['model']['core_tolerance']) >= self.iteration \
            else self.rmg['model']['core_tolerance'][-1]

    def run_rmg(self):
        """
        Run RMG.

        Raises:
            Various RMG Exceptions: if RMG crushed too many times.
        """
        if not os.path.isfile(self.paths['RMG input']):
            raise ValueError(f"The RMG input file {self.paths['RMG input']} could not be found.")

        self.logger.info(f'Running RMG (tolerance = {self.get_current_rmg_tol()})...')
        tic = time.time()

        # setup RMG
        initialize_rmg_log(
            verbose=self.verbose,
            log_file_name=self.paths['RMG log'],
        )
        rmg = RMG(input_file=self.paths['RMG input'], output_directory=self.paths['RMG'])
        rmg_kwargs = dict()
        if self.t3['options']['max_rmg_processes'] is not None:
            rmg_kwargs['maxproc'] = self.t3['options']['max_rmg_processes']
        if self.t3['options']['max_rmg_iterations'] is not None:
            rmg_kwargs['max_iterations'] = self.t3['options']['max_rmg_iterations']

        max_rmg_exceptions_allowed = self.t3['options']['max_RMG_exceptions_allowed']
        try:
            rmg.execute(initialize=True, **rmg_kwargs)

        except (ChemicallySignificantEigenvaluesError,
                ChemkinError,
                CollisionError,
                CoreError,
                DASPKError,
                ILPSolutionError,
                InvalidMicrocanonicalRateError,
                KineticsError,
                ModifiedStrongCollisionError,
                NetworkError,
                PressureDependenceError,
                ReactionError,
                ReservoirStateError,
                StatmechError,
                StatmechFitError) as e:
            self.logger.error(f'RMG Errored with {e.__class__}. Got the following trace:')
            self.logger.info(traceback.format_exc())

            if self.rmg_exceptions_counter > max_rmg_exceptions_allowed:
                self.logger.error(f'This is the {self.rmg_exceptions_counter} exception raised by RMG.\n'
                                  f'Not allowing additional exceptions, terminating.')
                raise
            else:
                self.logger.warning(f'This is the {self.rmg_exceptions_counter} exception raised by RMG.\n'
                                    f'The maximum number of exceptions allowed is {max_rmg_exceptions_allowed}.')
            self.rmg_exceptions_counter += 1

        except InputError:
            # this exception should not be raised if the schema is functioning properly
            self.logger.error('Something seems to be wrong with the RMG input file, please check your input.')
            raise

        elapsed_time = time_lapse(tic)
        self.logger.info(f'RMG terminated, execution time: {elapsed_time}')

    def run_sa(self) -> bool:
        """
        Run a sensitivity analysis.

        Returns:
            bool: Whether the SA ran successfully. ``True`` if it did.

        Todo:
            - Run SA per SA condition for batch instead of averaging Trange/Prange (extend the schema as needed)
        """
        sa_observables = list()
        for species in self.rmg['species']:
            if species['observable'] or species['SA_observable']:
                sa_observables.append(species['label'])

        method = self.t3['sensitivity']['adapter']
        self.logger.info(f'Running SA using {method}...')
        if method == 'RMG':
            if not os.path.isdir(self.paths['SA']):
                os.mkdir(self.paths['SA'])
            if os.path.isfile(self.paths['SA input']):
                os.remove(self.paths['SA input'])
            shutil.copyfile(src=self.paths['RMG input'], dst=self.paths['SA input'])

            rmg = load_rmg_py_job(
                input_file=self.paths['SA input'],
                chemkin_file=self.paths['chem annotated'],
                species_dict=self.paths['species dict'],
                generate_images=True,
                use_chemkin_names=False,
                check_duplicates=False,
            )

            rmg_observable_species = [species for species in rmg.reaction_model.core.species
                                      if species.label in sa_observables]

            for reaction_system in rmg.reaction_systems:
                reaction_system.sensitive_species = rmg_observable_species
                reaction_system.sensitivity_threshold = self.t3['sensitivity']['SA_threshold']
                if hasattr(reaction_system, 'Trange') and reaction_system.Trange is not None:
                    temperature = sum([t.value_si for t in reaction_system.Trange]) / len(reaction_system.Trange)
                else:
                    temperature = reaction_system.T.value_si
                reaction_system.sens_conditions['T'] = temperature
                if isinstance(reaction_system, SimpleReactor):
                    if hasattr(reaction_system, 'Prange') and reaction_system.Prange is not None:
                        pressure = sum([p.value_si for p in reaction_system.Prange]) / len(reaction_system.Prange)
                    else:
                        pressure = reaction_system.P.value_si
                    reaction_system.sens_conditions['P'] = pressure
                elif isinstance(reaction_system, LiquidReactor):
                    if hasattr(reaction_system, 'Vrange') and reaction_system.Vrange is not None:
                        volume = sum([v for v in reaction_system.Vrange]) / len(reaction_system.Vrange)
                    else:
                        volume = reaction_system.V
                    reaction_system.sens_conditions['V'] = volume
                else:
                    raise NotImplementedError(f'RMG SA not implemented for Reactor type {type(reaction_system)}.')
            try:
                simulate(rmg)
            except FileNotFoundError:
                return False
        else:
            raise NotImplementedError(f'Currently only RMG is implemented as a SA method.\nGot: {method}')
        return True

    def determine_species_to_calculate(self) -> bool:
        """
        Determine which species in the executed RMG job should be calculated.
        Species which were previously attempted to be calculated but did not converge
        will not be reconsidered.

        Updates:
            self.rmg_species
            self.rmg_reactions
            self.species - via self.add_species()
            self.qm['species'] - via self.add_species()
            self.executed_networks - via self.determine_species_from_pdep_network()

        Returns:
            bool: Whether additional calculations are required.
        """
        species_keys = list()

        self.rmg_species, self.rmg_reactions = self.load_species_and_reactions_from_chemkin_file()
        self.logger.info(f'This RMG model has {len(self.rmg_species)} species '
                         f'and {len(self.rmg_reactions)} reactions in its core.')

        if self.t3['options']['all_core_species']:
            for species in self.rmg_species:
                if self.species_requires_refinement(species=species):
                    species_keys.append(self.add_species(species=species, reasons=['All core species']))
        else:
            # 1. SA observables
            sa_obserables_exist = False
            for input_species in self.rmg['species']:
                if input_species['observable'] or input_species['SA_observable']:
                    sa_obserables_exist = True
                    species_keys.append(self.add_species(
                        species=get_species_by_label(input_species['label'], self.rmg_species),
                        reasons=['SA observable'],
                    ))
            # 2. SA
            if sa_obserables_exist:
                species_keys.extend(self.determine_species_based_on_sa())
            # 3. collision violators
            if self.t3['options']['collision_violators_thermo']:
                species_keys.extend(self.determine_species_based_on_collision_violators())

        species_keys = list(set([key for key in species_keys if key is not None]))

        additional_calcs_required = bool(len(species_keys))
        self.logger.info(f'Additional calculations required: {additional_calcs_required}\n')
        if additional_calcs_required:
            self.logger.log_species_to_calculate(species_keys, self.species)
        return additional_calcs_required

    def determine_species_based_on_sa(self) -> List[int]:
        """
        Determine species to calculate based on sensitivity analysis.

        Returns:
            List[int]: Entries are T3 species indices of species determined to be calculated based on SA.

        Adapters notes:
            - This method should remain and still return a list of indices.
            - Now we don't need to check whether the observables themselves should be refined,
              it's done in determine_species_to_calculate()
            - Perhaps the get_species_by_label() and get_reaction_by_index() functions
              should be relocated to the RMG adapter.
            - Consider passing self.logger to the adapters instead if instantiating a new Logger
              (there should only be a single instance of it).
            - Need to pass self.t3['sensitivity'] to the adapters instead of ``arguments``.
        """
        species_keys, pdep_rxns_to_explore = list(), list()
        if not os.path.isdir(self.paths['SA solver']):
            self.logger.error("Could not find the path to RMG's SA solver output folder.\n"
                              "Not performing refinement based on sensitivity analysis!")
            return species_keys

        sa_files = list()
        for file_ in os.listdir(self.paths['SA solver']):
            if 'sensitivity' in file_ and file_.endswith(".csv"):
                sa_files.append(file_)

        for sa_file in sa_files:
            # iterate through all SA .csv files in the solver folder
            df = pd.read_csv(os.path.join(self.paths['SA solver'], sa_file))
            sa_dict = {'rxn': dict(), 'spc': dict()}
            for header in df.columns:
                # iterate through all headers in the SA .csv file, but skip the `Time (s)` column
                sa_type = None
                if 'dln[k' in header and self.t3['sensitivity']['top_SA_reactions']:
                    sa_type = 'rxn'
                elif 'dG' in header and self.t3['sensitivity']['top_SA_species']:
                    sa_type = 'spc'
                if sa_type is not None:
                    # proceed only if we care about this column
                    entry = dict()
                    # check whether the observable requires calculations:
                    observable_label = header.split('[')[1].split(']')[0]
                    observable = get_species_by_label(observable_label, self.rmg_species)
                    if observable is None:
                        self.logger.error(f'Could not identify observable species {observable_label}!')
                    # continue with the parameter this column represents
                    observable_label = observable.to_chemkin()
                    if observable_label not in sa_dict[sa_type]:
                        sa_dict[sa_type][observable_label] = list()
                    # parameter extraction examples:
                    # for species get 'C2H4(8)' from `dln[ethane(1)]/dG[C2H4(8)]`
                    # for reaction, get 8 from `dln[ethane(1)]/dln[k8]: H(6)+ethane(1)=H2(12)+C2H5(5)`
                    parameter = header.split('[')[2].split(']')[0]
                    if sa_type == 'rxn':
                        parameter = int(parameter[1:])
                    entry['parameter'] = parameter  # rxn number or spc label
                    entry['max_sa'] = max(df[header].max(), abs(df[header].min()))  # the coefficient could be negative
                    sa_dict[sa_type][observable_label].append(entry)

            # get the top X entries from the SA
            for observable_label, sa_list in sa_dict['rxn'].items():
                sa_list_sorted = sorted(sa_list, key=lambda item: item['max_sa'], reverse=True)
                for i in range(min(self.t3['sensitivity']['top_SA_reactions'], len(sa_list_sorted))):
                    reaction = get_reaction_by_index(sa_list_sorted[i]['parameter'] - 1, self.rmg_reactions)
                    for species in reaction.reactants + reaction.products:
                        if self.species_requires_refinement(species=species):
                            num = f'{i+1}{get_ordinal_indicator(i+1)} ' if i else ''
                            reason = f'(i {self.iteration}) participates in the {num}most sensitive reaction ' \
                                     f'for {observable_label}: {reaction}'
                            species_keys.append(self.add_species(species=species, reasons=reason))
                    if reaction.kinetics.is_pressure_dependent() \
                            and reaction not in [rxn_tup[0] for rxn_tup in pdep_rxns_to_explore] \
                            and self.t3['sensitivity']['pdep_SA_threshold'] is not None:
                        pdep_rxns_to_explore.append((reaction, i, observable_label))
            for observable_label, sa_list in sa_dict['spc'].items():
                sa_list_sorted = sorted(sa_list, key=lambda item: item['max_sa'], reverse=True)
                for i in range(min(self.t3['sensitivity']['top_SA_species'], len(sa_list_sorted))):
                    species = get_species_by_label(sa_list_sorted[i]['parameter'], self.rmg_species)
                    if species is None:
                        self.logger.error(f"Could not identify species {sa_list_sorted[i]['parameter']}!")
                    if self.species_requires_refinement(species=species):
                        num = f'{i+1}{get_ordinal_indicator(i+1)} ' if i else ''
                        reason = f'(i {self.iteration}) the {num}most sensitive species thermo for {observable_label}'
                        species_keys.append(self.add_species(species=species, reasons=reason))

        species_keys.extend(self.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore))

        return species_keys

    def determine_species_from_pdep_network(self,
                                            pdep_rxns_to_explore: List[Tuple[Union[Reaction, PDepReaction], int, str]],
                                            ) -> List[int]:
        """
        Determine species to calculate based on a pressure dependent network
        by spawning network sensitivity analyses.
        First tries the CSE method, if unsuccessful tries the MSC method.

        Args:
            pdep_rxns_to_explore (List[Tuple[Reaction, int, str]]):
                Entries are tuples of (Reaction, SA rank index, observable_label).

        Returns:
            List[int]: Entries are T3 species indices of species determined to be calculated based on SA.
        """
        species_keys = list()
        if self.t3['sensitivity']['pdep_SA_threshold'] is None:
            return species_keys
        if not os.path.isdir(self.paths['PDep SA']):
            os.mkdir(self.paths['PDep SA'])

        for reaction_tuple in pdep_rxns_to_explore:
            reaction = reaction_tuple[0]

            if not isinstance(reaction, PDepReaction):
                continue

            # identify the network name and file name
            network_file_names = list()
            for (_, _, files) in os.walk(self.paths['RMG PDep']):
                network_file_names.extend(files)
                break  # don't continue to explore subdirectories
            network_file_names = [network_file_name for network_file_name in network_file_names
                                  if f'network{reaction.network.index}_' in network_file_name]
            if not network_file_names:
                # this PDepReaction did not stem from a network file, it is probably a library reaction
                if hasattr(reaction, 'library'):
                    self.logger.info(f'Not exploring library reaction {reaction} with PES SA.')
                else:
                    self.logger.warning(f'Not exploring reaction {reaction} with PES SA '
                                        f'since it does not have a `network` attribute.')
                continue
            network_version = max([int(network_file_name.split('.')[0].split('_')[1])
                                   for network_file_name in network_file_names])
            network_name = f'network{reaction.network.index}_{network_version}'  # w/o the '.py' extension

            # try running this network using either the CSE or the MSC method
            sa_coefficients_path, arkane = None, None
            errors = list()
            for method in PDEP_SA_ME_METHODS:
                isomer_labels = write_pdep_network_file(
                    network_name=network_name,
                    method=method,
                    pdep_sa_path=self.paths['PDep SA'],
                    rmg_pdep_path=self.paths['RMG PDep'],
                )
                arkane = Arkane(input_file=os.path.join(self.paths['PDep SA'], network_name, method, 'input.py'),
                                output_directory=os.path.join(self.paths['PDep SA'], network_name, method))
                arkane.plot = True
                self.logger.info(f'\nRuning Pdep SA for network {network_name} using the {method} method...')
                try:
                    arkane.execute()
                except (ChemicallySignificantEigenvaluesError,
                        ModifiedStrongCollisionError,
                        AttributeError,
                        ValueError,
                        TypeError) as e:
                    errors.append(e)
                else:
                    # network execution was successful, mark network as executed and don't run the next method
                    self.logger.info(f'Successfully executed a Pdep SA for network {network_name} '
                                     f'using the {method} method.\n')
                    self.executed_networks.append(isomer_labels)
                    sa_coefficients_path = os.path.join(self.paths['PDep SA'], network_name, method,
                                                        'sensitivity', 'sa_coefficients.yml')
                    break
            else:
                self.logger.error(f'Could not execute a Pdep SA for network {network_name} using '
                                  f'{PDEP_SA_ME_METHODS}.\nGot the following errors:')
                for method, e in zip(PDEP_SA_ME_METHODS, errors):
                    self.logger.info(f'{e.__class__} for method {method}:\n{e}\n')

            if sa_coefficients_path is not None:
                sa_dict = read_yaml_file(sa_coefficients_path)
                reactants_label = ' + '.join([reactant.to_chemkin() for reactant in reaction.reactants])
                products_label = ' + '.join([product.to_chemkin() for product in reaction.products])
                chemkin_reaction_str = f'{reactants_label} <=> {products_label}'
                labels_map = dict()  # keys are network species labels, values are Chemkin labels of the RMG species
                for network_label, adj in sa_dict['structures'].items():
                    labels_map[network_label] = get_species_label_by_structure(adj=adj, species_list=self.rmg_species)

                reactants_label = ' + '.join([key_by_val(labels_map, reactant.label) for reactant in reaction.reactants])
                products_label = ' + '.join([key_by_val(labels_map, product.label) for product in reaction.products])
                network_reaction_str = f'{reactants_label} <=> {products_label}'
                if network_reaction_str not in sa_dict:
                    self.logger.error(f'Could not locate reaction {network_reaction_str} '
                                      f'in SA output for network {network_name}.')
                else:
                    # identify wells in this network this reaction is sensitive to
                    sensitive_wells_dict = dict()  # keys are wells labels, values are lists of sensitive conditions
                    for condition, sa_data in sa_dict[network_reaction_str].items():
                        max_sa_coeff = max([sa_coeff for sa_coeff in sa_data.values()])
                        for entry, sa_coeff in sa_data.items():
                            if '(TS)' not in entry \
                                    and sa_coeff > max_sa_coeff * self.t3['sensitivity']['pdep_SA_threshold']:
                                if entry not in sensitive_wells_dict:
                                    sensitive_wells_dict[entry] = [condition]
                                else:
                                    sensitive_wells_dict[entry].append(condition)
                    if sensitive_wells_dict:
                        # extract species from wells and add to species_to_calc if thermo is uncertain
                        for well, conditions in sensitive_wells_dict.items():
                            species_list = list()
                            for label in well.split(' + '):
                                spc_label = labels_map[label]
                                species = None
                                if spc_label is not None:
                                    species = get_species_by_label(label=labels_map[label],
                                                                   species_list=self.rmg_species)
                                elif arkane is not None:
                                    # this is an Edge species which is missing from the Core rmg_species list
                                    species = get_species_by_label(label=label,
                                                                   species_list=arkane.species_dict.values())
                                if species is not None:
                                    species_list.append(species)
                            for species in species_list:
                                if self.species_requires_refinement(species=species):
                                    num = f'{reaction_tuple[1] + 1}{get_ordinal_indicator(reaction_tuple[1] + 1)} ' \
                                        if reaction_tuple[1] else ''
                                    reason = f'(i {self.iteration}) a sensitive well in PDep ' \
                                        f'network {network_name} from which {chemkin_reaction_str} was ' \
                                        f'derived, which is the {num}most sensitive reaction for observable ' \
                                        f'{reaction_tuple[2]}, at the {conditions}.'
                                    species_keys.append(self.add_species(species=species, reasons=reason))

        return species_keys

    def determine_species_based_on_collision_violators(self) -> List[int]:
        """
        Determine species to calculate based on collision rate violating reactions.

        Returns:
            List[int]: Entries are T3 species indices of species determined to be calculated based on SA.
        """
        species_keys = list()
        if not os.path.isfile(self.paths['RMG coll vio']):
            self.logger.info('No collision rate violating reactions identified in this model.')
            return species_keys

        with open(self.paths['RMG coll vio'], 'r') as f:
            lines = f.readlines()

        for line in lines:
            if line.count('=') == 1 and ('e+' in line or 'e-' in line) and '!' not in line:
                # `line` might look like one of these:
                # C2H3O(66)+O(T)(14)=C2H2O(60)+OH(D)(33)              1.500000e+09  1.500  -0.890
                # C2H2O(60)+CH2(T)(9)(+M)=C3H4O(383)(+M)              1.000e+00     0.000   0.000
                # C2H2O(60)+CH2(T)(9)(+N2)=C3H4O(383)(+N2)            1.000e+00     0.000   0.000
                # C2H2O(60)+CH2(T)(9)(+N2(32))=C3H4O(383)(+N2(32))    1.000e+00     0.000   0.000
                rxn_to_log = line.split()[0]
                collider = re.search(r'\(\+[^)]+\)', line)
                modified_line = line
                if collider is not None:
                    collider = collider.group(0)
                    if collider.count('(') == 2 and collider.count(')') == 1:
                        collider += ')'
                    modified_line = line.replace(collider, '')
                modified_rxn_string = modified_line.split()[0].replace('+M', '')
                labels = modified_rxn_string.split('=')
                reactants = labels[0].split('+') if '+' in labels[0] else [labels[0]]
                products = labels[1].split('+') if '+' in labels[1] else [labels[1]]
                labels = reactants + products
                for label in labels:
                    species = get_species_by_label(label, self.rmg_species)
                    if species is None:
                        self.logger.error(f'Could not identify species {label}!')
                    if self.species_requires_refinement(species=species):
                        reason = f'Species participates in collision rate violating reaction: {rxn_to_log}'
                        species_keys.append(self.add_species(species=species, reasons=reason))
        return species_keys

    def trsh_rmg_tol(self, factor: float = 0.5):
        """
        Troubleshoot by lowering the RMG tolerance used by the next iteration.

        Args:
            factor (float, optional): A factor by which to multiply the current RMG tolerance move to core.
        """
        if 0 < self.iteration < self.t3['options']['max_T3_iterations']:
            # there's at least one more iteration allowed
            # try troubleshooting the SA by reducing tolerance_move_to_core
            core_tolerance = self.rmg['model']['core_tolerance']
            report_trsh = False
            if len(core_tolerance) <= self.iteration:
                core_tolerance.extend([factor * core_tolerance[-1]] * (self.iteration + 1 - len(core_tolerance)))
                report_trsh = True
            elif core_tolerance[self.iteration] >= core_tolerance[self.iteration - 1]:
                core_tolerance[self.iteration] = factor * core_tolerance[self.iteration - 1]
                report_trsh = True
            if report_trsh:
                self.logger.info(f'Regenerating the RMG model with a tolerance move to core '
                                 f'of {factor * core_tolerance[self.iteration]}.')

    def species_requires_refinement(self, species: Species) -> bool:
        """
        Determine whether a species thermochemical properties
        should be calculated based on their uncertainty.
        First check that this species was not previously considered.

        Args:
            species (Species): The species for which the query is performed.

        Returns:
            bool: Whether the species thermochemical properties should be calculated. ``True`` if they should be.
        """
        if self.get_species_key(species=species) is None \
                and ('group additivity' in species.thermo.comment or '+ radical(' in species.thermo.comment):
            return True
        return False

    def get_species_key(self,
                        species: Optional[Species] = None,
                        label: Optional[str] = None,
                        label_type: str = 'QM',
                        ) -> Optional[int]:
        """
        Get a species key (the T3 species index) if the species exists in self.species.
        Either ``species`` or ``label`` must be given.

        Args:
            species (Species, optional): The species for which the query is performed.
            label (str, optional): The species label.
            label_type (str, optional): The label type, either 'RMG', 'Chemkin', or 'QM'.

        Returns:
            Optional[int]: The species T3 index if it exists, ``None`` if it does not.
        """
        if species is None and label is None:
            raise ValueError('Either species or label must be specified, got neither.')
        if label_type not in ['RMG', 'Chemkin', 'QM']:
            raise ValueError(f"label type must be either 'RMG', 'Chemkin' or 'QM', got: '{label_type}'.")
        for key, species_dict in self.species.items():
            if species is not None and species.is_isomorphic(species_dict['object']):
                return key
            if label is not None and label == species_dict[f'{label_type} label']:
                return key
        return None

    def load_species_and_reactions_from_chemkin_file(self) -> Tuple[List[Species], List[Reaction]]:
        """
        Load RMG Species and Reaction objects from the annotated Chemkin file.

        Raises:
            ChemkinError: If the Chemkin file could not be read.

        Returns:
            Tuple[List[Species], List[Reaction]]: The loaded RMG Species and Reaction objects.
        """
        try:
            rmg_species, rmg_reactions_ = load_chemkin_file(
                self.paths['chem annotated'],
                self.paths['species dict'],
                check_duplicates=True,
            )
        except ChemkinError:
            self.logger.error(f"Could not read the Chemkin file {self.paths['chem annotated']}!\n"
                              f"Trying to read it without checking for duplicate reactions...")
            try:
                rmg_species, rmg_reactions_ = load_chemkin_file(
                    self.paths['chem annotated'],
                    self.paths['species dict'],
                    check_duplicates=False,
                )
            except ChemkinError:
                self.logger.error(f"Could still not read the Chemkin file {self.paths['chem annotated']}!")
                raise
            else:
                self.logger.warning(f"Read the Chemkin file\n{self.paths['chem annotated']}\n"
                                    f"without checking for duplicate reactions.\nSA results might be inaccurate.")
        rmg_reactions = list()
        for i, reaction in enumerate(rmg_reactions_):
            # renumber, since duplicate reactions are removed and leave gaps in the index
            # also, now the index should match the RMG index - 1 rather than the Chemkin index - 1 (it is 0-indexed)
            reaction.index = i
            rmg_reactions.append(reaction)
        return rmg_species, rmg_reactions

    def add_species(self,
                    species: Species,
                    reasons: Union[List[str], str],
                    ) -> Optional[int]:
        """
        Add a species to self.species and to self.qm['species'].
        If the species already exists in self.species, only the reasons
        will be extended, and the species will not be considered
        in self.qm['species'].

        Args:
            species (Species): The species to consider.
            reasons (Union[List[str], str]): Reasons for calculating this species.

        Returns:
            Optional[int]: The T3 species index (the respective self.species key) if the species was just added,
                           ``None`` if the species already exists.
        """
        if isinstance(reasons, str):
            reasons = [reasons]
        key = self.get_species_key(species=species)
        if key is None:
            key = len(list(self.species.keys()))
            qm_species = species.copy(deep=False)
            legalize_species_label(species=qm_species)
            qm_species.label += f'_{key}'
            self.species[key] = {'RMG label': species.label,
                                 'Chemkin label': species.to_chemkin(),
                                 'QM label': qm_species.label,
                                 'object': species,
                                 'reasons': reasons,
                                 'converged': None,
                                 'iteration': self.iteration,
                                 }
            self.qm['species'].append(qm_species)
            return key

        # species already exists, extend reasons
        for reason in reasons:
            if reason not in self.species[key]['reasons']:
                self.species[key]['reasons'].append(reason)
        return None

    def add_to_rmg_library(self):
        """
        Creates RMG libraries in the RMG database repository
        if they don't already exist, and appends with the
        respective entries from the libraries generated by ARC.
        """
        arc_thermo_lib_path = self.paths['ARC thermo lib']
        rmg_thermo_lib_path = os.path.join(RMG_THERMO_LIB_BASE_PATH, f"{self.t3['options']['library_name']}.py")
        local_context = {
            'ThermoData': ThermoData,
            'Wilhoit': Wilhoit,
            'NASAPolynomial': NASAPolynomial,
            'NASA': NASA,
        }
        if os.path.isfile(arc_thermo_lib_path):
            if os.path.isfile(rmg_thermo_lib_path):
                # this thermo library already exists in the RMG database: Load it, append new entries, and save.
                rmg_thermo_lib, arc_thermo_lib = ThermoLibrary(), ThermoLibrary()
                rmg_thermo_lib.load(path=rmg_thermo_lib_path, local_context=local_context, global_context=dict())
                arc_thermo_lib.load(path=arc_thermo_lib_path, local_context=local_context, global_context=dict())
                arc_description = arc_thermo_lib.long_desc
                description_to_append = '\n'
                append = False
                for line in arc_description.splitlines():
                    if 'Overall time since project initiation' in line:
                        append = False
                    if append:
                        description_to_append += line + '\n'
                    if 'Considered the following' in line:
                        append = True
                rmg_thermo_lib.long_desc += description_to_append
                for entry in arc_thermo_lib.entries.values():
                    rmg_thermo_lib.entries[entry.label] = entry
                rmg_thermo_lib.save(path=rmg_thermo_lib_path)
            else:
                # this thermo library doesn't exist in the RMG database: Just copy the library generated by ARC.
                shutil.copy(arc_thermo_lib_path, rmg_thermo_lib_path)
            # Add the new thermo library to the RMG database to be used in the next iteration
            if self.t3['options']['library_name'] not in self.rmg['database']['thermo_libraries']:
                self.rmg['database']['thermo_libraries'].append(self.t3['options']['library_name'])

    def dump_species(self):
        """
        Dump self.species in case T3 needs to be restarted.
        """
        species = dict()
        for key, spc_dict in self.species.items():
            mod_spc_dict = {k: v for k, v in spc_dict.items() if k != 'object'}
            mod_spc_dict['adjlist'] = spc_dict['object'].molecule[0].to_adjacency_list()
            species[key] = mod_spc_dict
        save_yaml_file(path=os.path.join(self.project_directory, 'species.yml'), content=species)

    def load_species(self):
        """
        Load the dumped species dictionary into self.species.
        """
        if os.path.isfile(os.path.join(self.project_directory, 'species.yml')):
            species = read_yaml_file(path=os.path.join(self.project_directory, 'species.yml'))
            for key, spc_dict in species.items():
                mod_spc_dict = {k: v for k, v in spc_dict.items() if k != 'adjlist'}
                mod_spc_dict['object'] = Species().from_adjacency_list(spc_dict['adjlist'])
                self.species[key] = mod_spc_dict


def get_species_by_label(label: str,
                         species_list: list,
                         ) -> Optional[Species]:
    """
    Get a species from a list of species by its label.

    Args:
        label (str): A species label.
        species_list (list): Entries are RMG Species objects.

    Returns:
        Optional[Species]: The corresponding species from the species_list.
                           Returns ``None`` if no species was found.
    """
    for species in species_list:
        if species.label == label or species.to_chemkin() == label:
            return species
    if '(' in label and ')' in label:
        # try by the RMG species index
        for species in species_list:
            if species.index == int(label.split('(')[-1].split(')')[0]):
                return species
    return None


def get_reaction_by_index(index: int,
                          reactions: list,
                          ) -> Optional[Reaction]:
    """
    Get a reaction from a list of reactions by its index.

    Args:
        index (int): The reaction index attribute, 0-indexed.
        reactions (list): Entries are RMG Reaction objects.

    Returns:
        Optional[Reaction]: The corresponding reaction from the reaction list.
                            Returns ``None`` if no reaction was found.
    """
    for reaction in reactions:
        if reaction.index == index:
            return reaction
    return None


def legalize_species_label(species: Species):
    """
    ARC uses the species label as the folder name on the server and the local machine.
    Make sure a label is legal, correct it if it's not.

    Args:
        species (Species): A species object.
    """
    for char in species.label:
        if char not in valid_chars:
            species.label = species.molecule[0].get_formula()
            break
    else:
        if species.label[:2] == 'S(' and species.label[-1] == ')' \
                and all([char.isdigit() for char in species.label[2:-1]]):
            species.label = species.molecule[0].get_formula()


def get_species_label_by_structure(adj: str,
                                   species_list: list,
                                   ) -> Union[str, None]:
    """
    Get a species from a list of species by its structure (adjacency list).

    Args:
        adj (str): The species adjacency list.
        species_list (list): Entries are RMG Species objects.

    Returns:
        Union[str, None]: The corresponding species label attribute from the species_list.
                          Returns ``None`` if no species was found.
    """
    new_spc = Species().from_adjacency_list(adj)
    for spc in species_list:
        if spc.is_isomorphic(new_spc):
            return spc.label
    return None
