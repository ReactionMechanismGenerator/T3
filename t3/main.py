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
    - SA dict: dump and load for restart, inc. pdep SA
"""

import datetime
import inspect
import os
import re
import shutil
from typing import List, Optional, Tuple, Union

from arkane import Arkane
from rmgpy import settings as rmg_settings
from rmgpy.chemkin import load_chemkin_file
from rmgpy.data.kinetics import KineticsLibrary
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.exceptions import (ChemicallySignificantEigenvaluesError,
                              ChemkinError,
                              InputError,
                              ModifiedStrongCollisionError,
                              NetworkError,
                              )
from rmgpy.kinetics import Arrhenius, KineticsData
from rmgpy.reaction import Reaction
from rmgpy.rmg.pdep import PDepReaction
from rmgpy.species import Species
from rmgpy.thermo import NASAPolynomial, NASA, ThermoData, Wilhoit

from arc.common import get_ordinal_indicator, key_by_val, read_yaml_file, save_yaml_file
from arc.exceptions import ConverterError
from arc.main import ARC
from arc.species.species import ARCSpecies, check_label
from arc.species.converter import check_xyz_dict

from t3.common import PROJECTS_BASE_PATH, VALID_CHARS, delete_root_rmg_log, get_species_by_label, time_lapse
from t3.logger import Logger
from t3.runners.rmg_runner import rmg_runner
from t3.schema import InputBase
from t3.simulate.factory import simulate_factory
from t3.utils.writer import write_pdep_network_file, write_rmg_input_file

RMG_THERMO_LIB_BASE_PATH = os.path.join(rmg_settings['database.directory'], 'thermo', 'libraries')
RMG_KINETICS_LIB_BASE_PATH = os.path.join(rmg_settings['database.directory'], 'kinetics', 'libraries')


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
                'reasons': <List[str]: Reasons for calculating thermodynamic data for this species>,
                'converged': <Optional[bool]: whether thermo was successfully calculated, ``None`` if pending>,
                'iteration': <int: The iteration number in which this species was originally considered>,
            },
        }

        reactions = {
            <int: T3_reaction_index>: {
                'RMG label': <str: RMG label>,
                'Chemkin label': <str: Chemkin label>,
                'QM label': <str: The label used for the QM calc>,
                'SMILES label': <str: A reaction label that consists of the reactants/products SMILES>,
                'object': <Reaction: RMG Reaction object>,
                'reactant_keys': <List[int]: Keys of species that participate in this reaction as reactants>,
                'product_keys': <List[int]: Keys of species that participate in this reaction as products>,
                'reasons': <List[str]: Reasons for calculating the rate coefficient for this reaction>,
                'converged': <Optional[bool]: whether the rate coeff was successfully calculated, ``None`` if pending>,
                'iteration': <int: The iteration number in which this reaction was originally considered>,
            },
        }

        sa_dict = {
            'thermo': {
                <str: SA observable>: {
                    <str: RMG species label>: <array: 1D array with one entry per time point. Each entry is
                                     dLn(observable_1) / dG_species_1 in mol / kcal at the respective time>,
                }
            }
            'kinetcs': {
                <str: SA observable>: {
                    <int: reaction number>: <array: 1D array with one entry per time point. Each entry is
                                             dLn(observable_1) / dLn(k_1) at the respective time>,
                }
            }
            'time': <array: 1D array of time points in seconds>
        }

    Args:
        project (str): The project name.
        project_directory (str, optional): The project directory. Required through the API, optional through an input
                                           file (will be set to the directory of the input file if not specified).
        verbose (int, optional): The logging level, optional. 10 - debug, 20 - info, 30 - warning, default: 20.
        clean_dir (bool, optional): Whether to delete all existing files (other than input and submit files)
                                    and folders in the project directory prior to execution.
                                    If set to ``True``, the restart feature will not be triggered. Default: ``False``.
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
        t0 (datetime.datetime): Initial time when the project was spawned, stored as a datetime object.
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
        sa_observables (list): Entries are RMG species labels for the SA observables.
        sa_dict (dict): Dictionary with keys of `kinetics`, `thermo`, and `time`.
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

        self.sa_dict = None
        self.sa_observables = list()
        self.t0 = datetime.datetime.now()  # initialize the timer as datetime object

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
            self.cleanup()
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

        if any(self.rmg['model']['core_tolerance'][i + 1] > self.rmg['model']['core_tolerance'][i]
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
                                        including all default values. Default: ``False``.
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
        for self.iteration in range(iteration_start, max_t3_iterations + 1):

            self.logger.info(f'\n\n\nT3 iteration {self.iteration}:\n'
                             f'---------------\n')
            self.set_paths()

            # RMG
            if self.iteration > iteration_start or self.iteration == iteration_start and run_rmg_at_start:
                self.run_rmg()

            # SA
            if self.t3['sensitivity'] is not None:
                # determine species to run SA for
                if not self.sa_observables:
                    for species in self.rmg['species']:
                        if species['observable'] or species['SA_observable']:
                            self.sa_observables.append(species['label'])

                # use simulate_factory to create the corresponding adapter for SA
                simulate_adapter = simulate_factory(simulate_method=self.t3['sensitivity']['adapter'],
                                                    t3=self.t3,
                                                    rmg=self.rmg,
                                                    paths=self.paths,
                                                    logger=self.logger,
                                                    atol=self.rmg['model']['atol'],
                                                    rtol=self.rmg['model']['rtol'],
                                                    observable_list=self.sa_observables,
                                                    sa_atol=self.t3['sensitivity']['atol'],
                                                    sa_rtol=self.t3['sensitivity']['rtol'],
                                                    global_observables=None,
                                                    )
                simulate_adapter.simulate()
                # obtain the dictionary containing all SA coefficients for these species
                self.sa_dict = simulate_adapter.get_sa_coefficients()

            # determine what needs to be calculated
            additional_calcs_required = self.determine_species_and_reactions_to_calculate()

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

            if self.check_overtime():
                self.logger.log_max_time_reached(max_time=self.t3['options']['max_T3_walltime'])
                break

        if additional_calcs_required:
            # The main T3 loop terminated, but the latest calculations were not included in the model.
            # Run RMG for the last time.
            self.iteration += 1
            self.logger.info(f'\n\n\nT3 iteration {self.iteration} (just generating a model using RMG):\n'
                             f'------------------------------------------------------\n')
            self.set_paths()
            self.run_rmg()

        self.logger.log_species_summary(species_dict=self.species)
        self.logger.log_reactions_summary(reactions_dict=self.reactions)
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
            'RMS': os.path.join(iteration_path, 'RMG', 'rms'),
            'cantera annotated': os.path.join(iteration_path, 'RMG', 'cantera', 'chem_annotated.cti'),
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
                                     f"{self.qm['project'] if 'project' in self.qm else 'T3'}_info.yml"),
            'ARC thermo lib': os.path.join(iteration_path, 'ARC', 'output', 'RMG libraries', 'thermo',
                                           f"{self.qm['project'] if 'project' in self.qm else 'T3'}.py"),
            'ARC kinetics lib': os.path.join(iteration_path, 'ARC', 'output', 'RMG libraries', 'kinetics'),
            'RMG T3 thermo lib': os.path.join(RMG_THERMO_LIB_BASE_PATH, f"{self.t3['options']['library_name']}.py"),
            'RMG T3 kinetics lib': os.path.join(RMG_KINETICS_LIB_BASE_PATH, f"{self.t3['options']['library_name']}"),
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
            self.load_species_and_reactions()
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

        self.dump_species_and_reactions()

        arc_kwargs = arc_kwargs.copy()
        if 'adapter' in arc_kwargs:
            del arc_kwargs['adapter']

        arc_kwargs['project_directory'] = self.paths['ARC']
        if not os.path.isdir(arc_kwargs['project_directory']):
            os.makedirs(arc_kwargs['project_directory'])

        if 'project' not in arc_kwargs:
            arc_kwargs['project'] = 'T3'
        if 'species' in arc_kwargs.keys() and arc_kwargs['species']:
            species = list()
            for spc_ in arc_kwargs['species']:
                spc = ARCSpecies(species_dict=spc_) if isinstance(spc_, dict) else spc_
                key = self.get_species_key(species=spc)
                if key is not None \
                        and all('no need to compute thermo' in reason for reason in self.species[key]['reasons']):
                    species.append(ARCSpecies(rmg_species=spc, compute_thermo=False) if isinstance(spc, Species) else spc)
                else:
                    species.append(spc)
            arc_kwargs['species'] = species

        tic = datetime.datetime.now()
        arc = ARC(**arc_kwargs)
        if not os.path.isfile(self.paths['ARC input']):
            save_yaml_file(path=self.paths['ARC input'], content=arc.as_dict())

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
        Sets the self.species[<key>]['converged'] and the self.rxns[<key>]['converged'] parameters.

        Todo:
            - Check for non-physical species in unconverged species.
        """
        unconverged_spc_keys, converged_spc_keys = list(), list()
        unconverged_rxn_keys, converged_rxn_keys = list(), list()
        if os.path.isfile(self.paths['ARC info']):
            content = read_yaml_file(path=self.paths['ARC info'])
            for species in content['species']:
                key = self.get_species_key(label=species['label'])
                if key is not None:
                    if species['success']:
                        converged_spc_keys.append(key)
                        self.species[key]['converged'] = True
                    else:
                        unconverged_spc_keys.append(key)
                        self.species[key]['converged'] = False
            for reaction in content['reactions']:
                key = self.get_reaction_key(label=reaction['label'])
                if key is not None:
                    if reaction['success']:
                        converged_rxn_keys.append(key)
                        self.reactions[key]['converged'] = True
                    else:
                        unconverged_rxn_keys.append(key)
                        self.reactions[key]['converged'] = False
        else:
            raise ValueError(f'ARC did not save a project_info.yml file (expected to find it in {self.paths["ARC info"]}, '
                             f'something must be wrong.')
        self.logger.log_unconverged_species_and_reactions(
            species_keys=unconverged_spc_keys,
            species_dict=self.species,
            reaction_keys=unconverged_rxn_keys,
            reaction_dict=self.reactions,
        )
        if len(converged_spc_keys) or len(converged_rxn_keys):
            # we calculated something, add to thermo/kinetic library
            self.add_to_rmg_libraries()
        # clear the calculated objects from self.qm:
        self.qm['species'], self.qm['reactions'] = list(), list()
        self.dump_species_and_reactions()

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
        self.logger.info(f'Running RMG (tolerance = {self.get_current_rmg_tol()})...')

        # Use the RMG T3 libraries only if they exist and not already in use.
        # 1. thermo
        if self.t3['options']['library_name'] not in self.rmg['database']['thermo_libraries'] \
                and os.path.isfile(self.paths['RMG T3 thermo lib']):
            self.rmg['database']['thermo_libraries'] = [self.t3['options']['library_name']] + \
                                                       self.rmg['database']['thermo_libraries']
        elif self.t3['options']['library_name'] in self.rmg['database']['thermo_libraries'] \
                and not os.path.isfile(self.paths['RMG T3 thermo lib']):
            self.rmg['database']['thermo_libraries'].pop(self.rmg['database']['thermo_libraries'].index(
                self.t3['options']['library_name']))
        # 2. kinetics
        if self.t3['options']['library_name'] not in self.rmg['database']['kinetics_libraries'] \
                and os.path.isdir(self.paths['RMG T3 kinetics lib']):
            self.rmg['database']['kinetics_libraries'] = self.t3['options']['library_name'] + \
                                                         self.rmg['database']['kinetics_libraries']
        elif self.t3['options']['library_name'] in self.rmg['database']['kinetics_libraries'] \
                and not os.path.isdir(self.paths['RMG T3 kinetics lib']):
            self.rmg['database']['kinetics_libraries'].pop(self.rmg['database']['kinetics_libraries'].index(
                self.t3['options']['library_name']))

        write_rmg_input_file(
            rmg=self.rmg,
            t3=self.t3,
            iteration=self.iteration,
            path=self.paths['RMG input'],
            walltime=self.t3['options']['max_RMG_walltime'],
        )
        if not os.path.isfile(self.paths['RMG input']):
            raise ValueError(f"The RMG input file {self.paths['RMG input']} could not be written.")
        tic = datetime.datetime.now()

        max_rmg_exceptions_allowed = self.t3['options']['max_RMG_exceptions_allowed']
        rmg_exception_encountered = rmg_runner(rmg_input_file_path=self.paths['RMG input'],
                                               logger=self.logger,
                                               max_iterations=self.t3['options']['max_rmg_iterations'],
                                               verbose=self.verbose,
                                               )
        if rmg_exception_encountered:
            self.rmg_exceptions_counter += 1
            if self.rmg_exceptions_counter > max_rmg_exceptions_allowed:
                self.logger.error(f'This is the {self.rmg_exceptions_counter} exception raised by RMG.\n'
                                  f'Cannot allow more than {max_rmg_exceptions_allowed} RMG exceptions during a T3 run.\n'
                                  f'Not allowing additional exceptions, terminating.')
                raise ValueError('Terminating due to RMG exceptions.')
            else:
                self.logger.warning(f'RMG did not converge. '
                                    f'This is the {self.rmg_exceptions_counter} exception raised by RMG.\n'
                                    f'The maximum number of exceptions allowed is {max_rmg_exceptions_allowed}.')

        elapsed_time = time_lapse(tic)
        self.logger.info(f'RMG terminated, execution time: {elapsed_time}')

    def determine_species_and_reactions_to_calculate(self) -> bool:
        """
        Determine which species and reactions in the executed RMG job should be calculated.
        Species/reactions which were previously attempted to be calculated but did not converge
        will not be reconsidered.

        Updates:
            self.rmg_species
            self.rmg_reactions
            self.species - via self.add_species()
            self.qm['species'] - via self.add_species()
            self.reactions - via self.add_reaction()
            self.qm['reactions'] - via self.add_reaction()
            self.executed_networks - via self.determine_species_from_pdep_network()

        Returns:
            bool: Whether additional calculations are required.
        """
        species_keys, reaction_keys, coll_vio_spc_keys, coll_vio_rxn_keys = list(), list(), list(), list()

        self.rmg_species, self.rmg_reactions = self.load_species_and_reactions_from_chemkin_file()
        self.logger.info(f'This RMG model has {len(self.rmg_species)} species '
                         f'and {len(self.rmg_reactions)} reactions in its core.')

        sa_observables_exist = False
        for input_species in self.rmg['species']:
            if input_species['observable'] or input_species['SA_observable']:
                sa_observables_exist = True
                break

        if self.t3['options']['collision_violators_thermo'] or self.t3['options']['collision_violators_rates']:
            coll_vio_spc_keys, coll_vio_rxn_keys = self.determine_species_and_reactions_based_on_collision_violators()

        if self.t3['options']['all_core_species']:
            for species in self.rmg_species:
                if self.species_requires_refinement(species=species):
                    key = self.add_species(species=species, reasons=[f'(i {self.iteration}) All core species'])
                    if key is not None:
                        species_keys.append(key)

        # 1. Species
        else:
            # 1.1. SA observables
            for input_species in self.rmg['species']:
                if input_species['observable'] or input_species['SA_observable']:
                    if self.species_requires_refinement(species=get_species_by_label(input_species['label'],
                                                                                     self.rmg_species)):
                        key = self.add_species(species=get_species_by_label(input_species['label'], self.rmg_species),
                                               reasons=['SA observable'])
                        if key is not None:
                            species_keys.append(key)
            # 1.2. SA
            if sa_observables_exist:
                species_keys.extend(self.determine_species_based_on_sa())
            # 1.3. collision violators
            if self.t3['options']['collision_violators_thermo']:
                species_keys.extend(coll_vio_spc_keys)

        # 2. Reactions
        # 2.1. SA
        if sa_observables_exist:
            reaction_keys.extend(self.determine_reactions_based_on_sa())
        # 2.2. collision violators
        if self.t3['options']['collision_violators_rates']:
            reaction_keys.extend(coll_vio_rxn_keys)

        species_keys = list(set([key for key in species_keys if key is not None]))
        reaction_keys = list(set([key for key in reaction_keys if key is not None]))

        additional_calcs_required = bool(len(species_keys)) or bool(len(reaction_keys)) \
            or any(spc['converged'] is None for spc in self.species.values())

        self.logger.info(f'Additional calculations required: {additional_calcs_required}\n')
        if len(species_keys):
            self.logger.log_species_to_calculate(species_keys, self.species)
        if len(reaction_keys):
            self.logger.log_reactions_to_calculate(reaction_keys, self.reactions)
        return additional_calcs_required

    def determine_species_based_on_sa(self) -> List[int]:
        """
        Determine species to calculate based on sensitivity analysis.

        Returns:
            List[int]: Entries are T3 species indices of species determined to be calculated based on SA.
        """
        species_keys, pdep_rxns_to_explore = list(), list()
        if self.sa_dict is None:
            self.logger.error(f"T3's sa_dict was None. Please check that the input file contains a proper 'sensitivity' "
                              f"block and/or that SA was run successfully.\n"
                              f"Not performing refinement based on sensitivity analysis!")
            return species_keys

        # Create a new dictionary that stores the absolute maximum value from each SA.
        sa_dict_max = {'kinetics': dict(), 'thermo': dict()}
        for key in ['kinetics', 'thermo']:
            for observable_label in self.sa_dict[key].keys():
                if observable_label not in sa_dict_max[key]:
                    sa_dict_max[key][observable_label] = list()
                for parameter in self.sa_dict[key][observable_label].keys():
                    entry = dict()
                    entry['parameter'] = parameter  # rxn number as int or spc label as str
                    entry['max_sa'] = max(self.sa_dict[key][observable_label][parameter].max(), abs(
                        self.sa_dict[key][observable_label][parameter].min()))  # the coefficient could be negative
                    sa_dict_max[key][observable_label].append(entry)

            # Get the top X entries from the SA.
            for observable_label, sa_list in sa_dict_max['kinetics'].items():
                sa_list_sorted = sorted(sa_list, key=lambda item: item['max_sa'], reverse=True)
                for i in range(min(self.t3['sensitivity']['top_SA_reactions'], len(sa_list_sorted))):
                    reaction = get_reaction_by_index(sa_list_sorted[i]['parameter'] - 1, self.rmg_reactions)
                    for species in reaction.reactants + reaction.products:
                        if self.species_requires_refinement(species=species):
                            num = f'{i+1}{get_ordinal_indicator(i+1)} ' if i else ''
                            reason = f'(i {self.iteration}) participates in the {num}most sensitive reaction ' \
                                     f'for {observable_label}: {reaction}'
                            key = self.add_species(species=species, reasons=reason)
                            if key is not None:
                                species_keys.append(key)
                    if reaction.kinetics.is_pressure_dependent() \
                            and reaction not in [rxn_tup[0] for rxn_tup in pdep_rxns_to_explore] \
                            and self.t3['sensitivity']['pdep_SA_threshold'] is not None:
                        pdep_rxns_to_explore.append((reaction, i, observable_label))
            for observable_label, sa_list in sa_dict_max['thermo'].items():
                sa_list_sorted = sorted(sa_list, key=lambda item: item['max_sa'], reverse=True)
                for i in range(min(self.t3['sensitivity']['top_SA_species'], len(sa_list_sorted))):
                    species = get_species_by_label(sa_list_sorted[i]['parameter'], self.rmg_species)
                    if species is None:
                        self.logger.error(f"Could not identify species {sa_list_sorted[i]['parameter']}!")
                    if self.species_requires_refinement(species=species):
                        num = f'{i+1}{get_ordinal_indicator(i+1)} ' if i else ''
                        reason = f'(i {self.iteration}) the {num}most sensitive species thermo for {observable_label}'
                        key = self.add_species(species=species, reasons=reason)
                        if key is not None:
                            species_keys.append(key)

        species_keys.extend(self.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore))

        return species_keys

    def determine_reactions_based_on_sa(self) -> List[int]:
        """
        Determine reaction rate coefficients to calculate based on sensitivity analysis.

        Returns:
            List[int]: Entries are T3 reaction indices of reactions determined to be calculated based on SA.
        """
        reaction_keys, pdep_rxns_to_explore = list(), list()
        if not os.path.isdir(self.paths['SA solver']):
            self.logger.error("Could not find the path to RMG's SA solver output folder.\n"  # todo: it doesn't have to be RMG
                              "Not performing refinement based on sensitivity analysis!")
            return reaction_keys

        # Create a new dictionary that stores the absolute maximum value from each SA.
        sa_dict_max = dict()
        for observable_label in self.sa_dict['kinetics'].keys():
            if observable_label not in sa_dict_max:
                sa_dict_max[observable_label] = list()
            for parameter in self.sa_dict['kinetics'][observable_label].keys():
                entry = dict()
                entry['parameter'] = parameter  # rxn number as int or spc label as str
                entry['max_sa'] = max(self.sa_dict['kinetics'][observable_label][parameter].max(),
                                      abs(self.sa_dict['kinetics'][observable_label][parameter].min()))
                sa_dict_max[observable_label].append(entry)

        # Get the top X entries from the SA.
        for observable_label, sa_list in sa_dict_max.items():
            sa_list_sorted = sorted(sa_list, key=lambda item: item['max_sa'], reverse=True)
            for i in range(min(self.t3['sensitivity']['top_SA_reactions'], len(sa_list_sorted))):
                reaction = get_reaction_by_index(sa_list_sorted[i]['parameter'] - 1, self.rmg_reactions)
                if self.reaction_requires_refinement(reaction=reaction):
                    num = f'{i+1}{get_ordinal_indicator(i+1)} ' if i else ''
                    reason = f'(i {self.iteration}) the {num}most sensitive reaction for {observable_label}'
                    key = self.add_reaction(reaction=reaction, reasons=reason)
                    if key is not None:
                        reaction_keys.append(key)
                # if reaction.kinetics.is_pressure_dependent() \
                #         and reaction not in [rxn_tup[0] for rxn_tup in pdep_rxns_to_explore] \
                #         and self.t3['sensitivity']['pdep_SA_threshold'] is not None:
                #     pdep_rxns_to_explore.append((reaction, i, observable_label))

        # reaction_keys.extend(self.determine_reactions_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore))  # Todo: make sure rate coefficients are perturbed in Arkane PDep SA, consider these rxns here

        return reaction_keys

    def determine_species_from_pdep_network(self,
                                            pdep_rxns_to_explore: List[Tuple[Union[Reaction, PDepReaction], int, str]],
                                            ) -> List[int]:
        """
        Determine species to calculate based on a pressure dependent network
        by spawning network sensitivity analyses.

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

            # Try running this network using user-specified methods by order.
            sa_coefficients_path, arkane = None, None
            errors = list()
            for method in self.t3['sensitivity']['ME_methods']:
                isomer_labels = write_pdep_network_file(
                    network_name=network_name,
                    method=method,
                    pdep_sa_path=self.paths['PDep SA'],
                    rmg_pdep_path=self.paths['RMG PDep'],
                )
                arkane = Arkane(input_file=os.path.join(self.paths['PDep SA'], network_name, method, 'input.py'),
                                output_directory=os.path.join(self.paths['PDep SA'], network_name, method))
                arkane.plot = True
                self.logger.info(f'\nRunning PDep SA for network {network_name} using the {method} method...')
                try:
                    arkane.execute()
                except (AttributeError,
                        ChemicallySignificantEigenvaluesError,
                        ModifiedStrongCollisionError,
                        NetworkError,
                        TypeError,
                        UnboundLocalError,
                        ValueError,
                        ) as e:
                    errors.append(e)
                else:
                    # Network execution was successful, mark network as executed and don't run the next method.
                    self.logger.info(f'Successfully executed a PDep SA for network {network_name} '
                                     f'using the {method} method.\n')
                    self.executed_networks.append(isomer_labels)
                    sa_coefficients_path = os.path.join(self.paths['PDep SA'], network_name, method,
                                                        'sensitivity', 'sa_coefficients.yml')
                    break
            else:
                self.logger.error(f"Could not execute a PDep SA for network {network_name} using "
                                  f"{self.t3['sensitivity']['ME_methods']}.\nGot the following errors:")
                for method, e in zip(self.t3['sensitivity']['ME_methods'], errors):
                    self.logger.info(f'{e.__class__} for method {method}:\n{e}\n')

            if sa_coefficients_path is not None:
                sa_dict = read_yaml_file(sa_coefficients_path)
                reactants_label = ' + '.join([reactant.to_chemkin() for reactant in reaction.reactants])
                products_label = ' + '.join([product.to_chemkin() for product in reaction.products])
                chemkin_reaction_str = f'{reactants_label} <=> {products_label}'
                labels_map = dict()  # Keys are network species labels, values are Chemkin labels of the RMG species.
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
                                        f'{reaction_tuple[2]}, at {conditions}.'
                                    key = self.add_species(species=species, reasons=reason)
                                    if key is not None:
                                        species_keys.append(key)

        return species_keys

    def determine_species_and_reactions_based_on_collision_violators(self) -> Tuple[List[int], List[int]]:
        """
        Determine species to calculate based on collision rate violating reactions.

        Returns:
            Tuple[List[int], List[int]]:
                - Entries are T3 species indices of species determined to be calculated based on SA.
                - Entries are T3 reaction indices of reactions determined to be calculated based on SA.
        """
        species_keys, reaction_keys = list(), list()
        if not os.path.isfile(self.paths['RMG coll vio']):
            self.logger.info('No collision rate violating reactions identified in this model.')
            return species_keys, reaction_keys

        with open(self.paths['RMG coll vio'], 'r') as f:
            lines = f.readlines()

        for line in lines:
            if line.count('=') == 1 and ('e+' in line or 'e-' in line) and '!' not in line:
                # ``line`` might look like one of these:
                # C2H3O(66)+O(T)(14)=C2H2O(60)+OH(D)(33)              1.500000e+09  1.500  -0.890
                # C2H2O(60)+CH2(T)(9)+M=C3H4O(383)+M                  1.000e+10     0.000   5.000
                # C2H2O(60)+CH2(T)(9)(+M)=C3H4O(383)(+M)              1.000e+00     0.000   0.000
                # C2H2O(60)+CH2(T)(9)(+N2)=C3H4O(383)(+N2)            1.000e+00     0.000   0.000
                # C2H2O(60)+CH2(T)(9)(+N2(32))=C3H4O(383)(+N2(32))    1.000e+00     0.000   0.000
                rxn_to_log = line.split()[0]
                collider = re.search(r'\(\+[^)]+\)', line)
                modified_line = line
                if collider is not None:
                    # Todo: Not considering colliders here. Change this when automating termolecular reactions.
                    collider = collider.group(0)
                    if collider.count('(') == 2 and collider.count(')') == 1:
                        collider += ')'
                    modified_line = line.replace(collider, '')
                modified_rxn_string = modified_line.split()[0].replace('+M', '')
                labels = modified_rxn_string.split('=')
                reactant_labels = labels[0].split('+') if '+' in labels[0] else [labels[0]]
                product_labels = labels[1].split('+') if '+' in labels[1] else [labels[1]]

                # 1. Species
                labels = reactant_labels + product_labels
                for label in labels:
                    species = get_species_by_label(label, self.rmg_species)
                    if species is None:
                        self.logger.error(f'Could not identify species {label}!')
                    if self.species_requires_refinement(species=species):
                        reason = f'(i {self.iteration}) Species participates in collision rate violating ' \
                                 f'reaction: {rxn_to_log}'
                        key = self.add_species(species=species, reasons=reason)
                        if key is not None:
                            species_keys.append(key)

                # 2. Reactions (not considering colliders for now)
                reactants = [get_species_by_label(label, self.rmg_species) for label in reactant_labels]
                products = [get_species_by_label(label, self.rmg_species) for label in product_labels]
                if not len(reactants) or not len(products):
                    self.logger.error(f'Could not identify reaction {rxn_to_log}!')
                reaction = Reaction(reactants=reactants, products=products)
                if self.reaction_requires_refinement(reaction) \
                        and not any(self.species_requires_refinement(species=spc) for spc in reactants + products):
                    # only consider a rate violating reaction if all the thermo was first fixed
                    reason = f'(i {self.iteration}) Reaction rate coefficient violates the collision rate.'
                    key = self.add_reaction(reaction=reaction, reasons=reason)
                    if key is not None:
                        reaction_keys.append(key)

        return species_keys, reaction_keys

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

    def species_requires_refinement(self, species: Optional[Species]) -> bool:
        """
        Determine whether a species thermochemical properties
        should be calculated based on the data uncertainty.
        First check that this species was not previously considered.

        Args:
            species (Species): The species for which the query is performed.

        Returns:
            bool: Whether the species thermochemical properties should be calculated. ``True`` if they should be.
        """
        if species is None:
            return False
        thermo_comment = species.thermo.comment.split('Solvation')[0]
        if (self.get_species_key(species=species) is None
            or self.species[self.get_species_key(species=species)]['converged'] is None) \
                and ('group additivity' in thermo_comment or '+ radical(' in thermo_comment):
            return True
        return False

    def reaction_requires_refinement(self, reaction: Reaction) -> bool:
        """
        Determine whether a reaction rate coefficient
        should be calculated based on the data uncertainty.
        First check that this reaction was not previously considered.
        Args:
            reaction (Reaction): The reaction for which the query is performed.
        Returns:
            bool: Whether the reaction rate coefficient should be calculated. ``True`` if it should be.
        Todo:
            Consider cases such as:
            #
            ! BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
            !     Total Standard Deviation in ln(k): 2.14551182899
            ! Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
            ! Euclidian distance = 0
            ! family: R_Recombination
            H(5)+S(752)=C12H26(1)                               1.766370e+13 0.153     0.000
            #
            Also, sometimes a rate rule is used exactly, but the tree might be too generic
            and the reaction should in fact be computed.
            #
            PDep reactions: need to check whether the reaction is elementary, then a high pressure limit of it could be
            computed by ARC. But if it's a network reaction, then don't send it directly into ARC, make sure to check the wells and path reactions
            via Arkane SA.
        """
        if isinstance(reaction, LibraryReaction):
            return False
        kinetics_comment = reaction.kinetics.comment if reaction.kinetics is not None else ''
        if self.get_reaction_key(reaction=reaction) is None \
                and 'Exact match found for rate rule' not in kinetics_comment \
                and ('Estimated using an average for rate rule' in kinetics_comment
                     or ('Estimated using template' in kinetics_comment and 'for rate rule' in kinetics_comment)
                     or ('Estimated using average of templates' in kinetics_comment
                         and 'for rate rule' in kinetics_comment)
                     or '' in kinetics_comment):
            return True
        return False

    def get_species_key(self,
                        species: Optional[Union[Species, ARCSpecies]] = None,
                        label: Optional[str] = None,
                        label_type: str = 'QM',
                        ) -> Optional[int]:
        """
        Get a species key (the T3 species index) if the species exists in self.species.
        Either ``species`` or ``label`` must be given.

        Args:
            species (Union[Species, ARCSpecies], optional): The species for which the query is performed.
            label (str, optional): The species label.
            label_type (str, optional): The label type, either 'RMG', 'Chemkin', or 'QM' (default: 'QM').

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

    def get_reaction_key(self,
                         reaction: Optional[Reaction] = None,
                         label: Optional[str] = None,
                         label_type: str = 'QM',
                         ) -> Optional[int]:
        """
        Get a reaction key (the T3 reaction index) if the reaction exists in self.reactions.
        Either ``reaction`` or ``label`` must be given.

        Args:
            reaction (Reaction, optional): The reaction for which the query is performed.
            label (str, optional): The reaction label.
            label_type (str, optional): The label type, either 'RMG', 'Chemkin', 'QM', or 'SMILES' (default: 'QM').

        Returns:
            Optional[int]: The reaction T3 index if it exists, ``None`` if it does not.
        """
        if reaction is None and label is None:
            raise ValueError('Either reaction or label must be specified, got neither.')
        if label_type not in ['RMG', 'Chemkin', 'QM', 'SMILES']:
            raise ValueError(f"label type must be either 'RMG', 'Chemkin' or 'QM', got: '{label_type}'.")
        for key, reaction_dict in self.reactions.items():
            if reaction is not None and reaction.is_isomorphic(reaction_dict['object']):
                return key
            if label is not None and label == reaction_dict[f'{label_type} label']:
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
                self.logger.error(f"Still could not read the Chemkin file {self.paths['chem annotated']}!")
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
        will be updated (extended), and the species will not be considered
        in self.qm['species'].

        Args:
            species (Species): The species to consider.
            reasons (Union[List[str], str]): Reasons for calculating this species.

        Returns:
            Optional[int]: The T3 species index (the respective self.species key) if the species was just added,
                           ``None`` if the species already exists.
        """
        reasons = [reasons] if isinstance(reasons, str) else reasons
        key = self.get_species_key(species=species)
        if key is None:
            key = len(list(self.species.keys()))
            qm_species = get_species_with_qm_label(species=species, key=key)
            self.species[key] = {'RMG label': species.label,
                                 'Chemkin label': species.to_chemkin(),
                                 'QM label': qm_species.label,
                                 'object': species,
                                 'reasons': reasons,
                                 'converged': None,
                                 'iteration': self.iteration,
                                 }

            # Check whether T3 has xyz information for this species, if so process it.
            for rmg_species in self.rmg['species']:
                if rmg_species['label'] == species.label and rmg_species['xyz'] is not None:
                    xyzs = list()
                    if not isinstance(rmg_species['xyz'], list):
                        rmg_species['xyz'] = [rmg_species['xyz']]
                    for xyz in rmg_species['xyz']:
                        # Only pass valid xyzs to ARC.
                        try:
                            xyz_dict = check_xyz_dict(xyz)
                        except ConverterError:
                            pass
                        else:
                            xyzs.append(xyz_dict)
                    if len(xyzs):
                        if self.qm['adapter'] == 'ARC':
                            # Make qm_species and ARCSpecies instance to consider the xyz information
                            qm_species = ARCSpecies(label=qm_species.label,
                                                    rmg_species=qm_species,
                                                    xyz=xyzs,
                                                    )
                        else:
                            raise NotImplementedError(f"Passing XYZ information to {self.qm['adapter']} "
                                                      f"is not yet implemented.")
            self.qm['species'].append(qm_species)
            return key

        # The species already exists, extend reasons.
        for reason in reasons:
            if reason not in self.species[key]['reasons']:
                self.species[key]['reasons'].append(reason)
        return None

    def add_reaction(self,
                     reaction: Reaction,
                     reasons: Union[List[str], str],
                     ) -> Optional[int]:
        """
        Add a species to self.species and to self.qm['species'].
        If the species already exists in self.species, only the reasons
        will be updated (extended), and the species will not be considered
        in self.qm['species'].
        Args:
            reaction (Reaction): The reaction to consider.
            reasons (Union[List[str], str]): Reasons for calculating this reaction.
        Returns:
            Optional[int]: The T3 reaction index (the respective self.reactions key) if the reaction was just added,
                           ``None`` if the reaction already exists.
        Todo:
            Add tests.
        """
        reasons = [reasons] if isinstance(reasons, str) else reasons
        rxn_key = self.get_reaction_key(reaction=reaction)
        if rxn_key is None:
            chemkin_label = ''
            try:
                chemkin_label = reaction.to_chemkin()
            except (AttributeError, ChemkinError, TypeError) as e:
                self.logger.debug(f'Could not generate a Chemkin label for reaction {reaction}. Got:\n{e}')

            rxn_key = len(list(self.reactions.keys()))
            for spc in reaction.reactants + reaction.products:
                if self.get_species_key(species=spc) is None:
                    self.add_species(species=spc, reasons=f'(i {self.iteration}) Participates in a reaction for which '
                                                          f'a rate coefficient is computed.')

            reaction.reactants = [get_species_with_qm_label(species=spc, key=self.get_species_key(species=spc))
                                  for spc in reaction.reactants]
            reaction.products = [get_species_with_qm_label(species=spc, key=self.get_species_key(species=spc))
                                 for spc in reaction.products]

            qm_label = ' <=> '.join([' + '.join([spc.label for spc in species_list])
                                     for species_list in [reaction.reactants, reaction.products]])
            smiles_label = ' <=> '.join([' + '.join([spc.molecule[0].to_smiles() for spc in species_list])
                                         for species_list in [reaction.reactants, reaction.products]])
            reaction.label = qm_label

            self.reactions[rxn_key] = {'RMG label': reaction.label or str(reaction),
                                       'Chemkin label': chemkin_label,
                                       'QM label': qm_label,
                                       'SMILES label': smiles_label,
                                       'object': reaction,
                                       'reasons': reasons,
                                       'converged': None,
                                       'iteration': self.iteration,
                                       'reactant_keys': [self.get_species_key(species=spc) for spc in reaction.reactants],
                                       'product_keys': [self.get_species_key(species=spc) for spc in reaction.products],
                                       }
            qm_reaction = reaction.copy()
            qm_reaction.label = qm_label
            self.qm['reactions'].append(qm_reaction)
            return rxn_key

        # The reaction already exists, extend reasons.
        for reason in reasons:
            if reason not in self.reactions[rxn_key]['reasons']:
                self.reactions[rxn_key]['reasons'].append(reason)
        return None

    def add_to_rmg_libraries(self):
        """
        Creates RMG libraries in the RMG database repository if they don't already exist,
        and appends with the respective entries from the libraries generated by ARC.

        Todo:
            Tests kinetics libraries.
        """
        # 1. Thermo:
        arc_thermo_lib_path = self.paths['ARC thermo lib']
        rmg_t3_thermo_lib_path = self.paths['RMG T3 thermo lib']
        local_context = {
            'ThermoData': ThermoData,
            'Wilhoit': Wilhoit,
            'NASAPolynomial': NASAPolynomial,
            'NASA': NASA,
        }
        if os.path.isfile(arc_thermo_lib_path) and os.path.isfile(rmg_t3_thermo_lib_path):
            # This thermo library already exists in the RMG database: Load it, append new entries, and save.
            rmg_thermo_lib, arc_thermo_lib = ThermoLibrary(), ThermoLibrary()
            rmg_thermo_lib.load(path=rmg_t3_thermo_lib_path, local_context=local_context, global_context=dict())
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
                entry_species = Species(molecule=[entry.item])
                entry_species.generate_resonance_structures(keep_isomorphic=False, filter_structures=True)
                for existing_entry in rmg_thermo_lib.entries.values():
                    if entry_species.is_isomorphic(existing_entry.item):
                        self.logger.warning(f"Not adding species {entry.label} to the "
                                            f"{self.t3['options']['library_name']} thermo library, "
                                            f"the species seems to already exist under the label "
                                            f"{existing_entry.label}.")
                        break
                rmg_thermo_lib.entries[entry.label] = entry
            rmg_thermo_lib.save(path=rmg_t3_thermo_lib_path)
        else:
            # This thermo library doesn't exist in the RMG database: Just copy the library generated by ARC.
            if os.path.isfile(arc_thermo_lib_path):
                shutil.copy(arc_thermo_lib_path, rmg_t3_thermo_lib_path)

        # 2. Kinetics:
        arc_kinetics_lib_path = self.paths['ARC kinetics lib']
        rmg_t3_kinetics_lib_path = self.paths['RMG T3 kinetics lib']
        local_context = {
            'KineticsData': KineticsData,
            'Arrhenius': Arrhenius,
        }
        if os.path.isfile(arc_kinetics_lib_path) and os.path.isfile(rmg_t3_kinetics_lib_path):
            # This kinetics library already exists in the RMG database: Load it, append new entries, and save.
            rmg_kinetics_lib, arc_kinetics_lib = KineticsLibrary(), KineticsLibrary()
            rmg_kinetics_lib.load(path=rmg_t3_kinetics_lib_path, local_context=local_context, global_context=dict())
            arc_kinetics_lib.load(path=arc_kinetics_lib_path, local_context=local_context, global_context=dict())
            arc_description = arc_kinetics_lib.long_desc
            description_to_append = '\n'
            append = False
            for line in arc_description.splitlines():
                if 'Overall time since project initiation' in line:
                    append = False
                if append:
                    description_to_append += line + '\n'
                if 'Considered the following' in line:
                    append = True
            rmg_kinetics_lib.long_desc += description_to_append
            for entry in arc_kinetics_lib.entries.values():
                entry_reaction = Reaction(reactants=entry.item.reactants[:],
                                          products=entry.item.products[:],
                                          specific_collider=entry.item.specific_collider,
                                          kinetics=entry.data,
                                          duplicate=entry.item.duplicate,
                                          reversible=entry.item.reversible,
                                          allow_pdep_route=entry.item.allow_pdep_route,
                                          elementary_high_p=entry.item.elementary_high_p,
                                          )
                for existing_entry in rmg_kinetics_lib.entries.values():
                    if entry_reaction.is_isomorphic(existing_entry.item):
                        self.logger.warning(f"Not adding reaction {entry.label} to the "
                                            f"{self.t3['options']['library_name']} kinetics library, "
                                            f"the reaction seems to already exist under the label "
                                            f"{existing_entry.label}.")
                        break
                rmg_kinetics_lib.entries[entry.label] = entry
            rmg_kinetics_lib.save(path=rmg_t3_kinetics_lib_path)
        else:
            # This kinetics library doesn't exist in the RMG database: Just copy the library generated by ARC.
            if os.path.isfile(arc_kinetics_lib_path):
                shutil.copy(arc_kinetics_lib_path, rmg_t3_kinetics_lib_path)

    def dump_species_and_reactions(self):
        """
        Dump self.species and self.reactions in case T3 needs to be restarted.
        """
        species, reactions = dict(), dict()
        for key, spc_dict in self.species.items():
            mod_spc_dict = {k: v for k, v in spc_dict.items() if k != 'object'}
            mod_spc_dict['adjlist'] = spc_dict['object'].molecule[0].to_adjacency_list()
            species[key] = mod_spc_dict
        save_yaml_file(path=os.path.join(self.project_directory, 'species.yml'), content=species)
        for key, rxn_dict in self.reactions.items():
            mod_rxn_dict = {k: v for k, v in rxn_dict.items() if k != 'object'}
            reactions[key] = mod_rxn_dict
        save_yaml_file(path=os.path.join(self.project_directory, 'reactions.yml'), content=reactions)

    def load_species_and_reactions(self):
        """
        Load the dumped species and reactions dictionaries into self.species and self.reactions, respectively.
        Resurrect the Species and Reaction objects.
        """
        if os.path.isfile(os.path.join(self.project_directory, 'species.yml')):
            species = read_yaml_file(path=os.path.join(self.project_directory, 'species.yml'))
            for key, spc_dict in species.items():
                mod_spc_dict = {k: v for k, v in spc_dict.items() if k != 'adjlist'}
                mod_spc_dict['object'] = Species().from_adjacency_list(spc_dict['adjlist'])
                self.species[key] = mod_spc_dict
        if os.path.isfile(os.path.join(self.project_directory, 'reactions.yml')):
            reactions = read_yaml_file(path=os.path.join(self.project_directory, 'reactions.yml'))
            for key, rxn_dict in reactions.items():
                mod_rxn_dict = rxn_dict.copy()
                mod_rxn_dict['object'] = Reaction(reactants=[spc_dict['object']
                                                             for key, spc_dict in self.species.items()
                                                             if key in mod_rxn_dict['reactant_keys']],
                                                  products=[spc_dict['object']
                                                            for key, spc_dict in self.species.items()
                                                            if key in mod_rxn_dict['product_keys']],
                                                  )
                self.reactions[key] = mod_rxn_dict

    def check_overtime(self) -> bool:
        """
        Check that the timer hasn't run out.

        Returns:
            bool: Whether T3 is running over time.
        """
        if self.t3['options']['max_T3_walltime'] is not None:
            delta = time_lapse(self.t0)
            splits = self.t3['options']['max_T3_walltime'].split(':')  # 01:00:00:00
            max_delta = datetime.timedelta(days=int(splits[0]),
                                           hours=int(splits[1]),
                                           minutes=int(splits[2]),
                                           seconds=int(splits[3]),
                                           )
            if delta > max_delta:
                return True
        return False

    def cleanup(self):
        """
        Clean the working directory other than the input and submit files and the log archive folder.
        """
        for root, dirs, files in os.walk(self.project_directory, topdown=True):
            for file_ in files:
                if 'input' not in file_ and 'submit' not in file_:
                    os.remove(os.path.join(root, file_))
            for folder in dirs:
                if folder != 'log_archive':
                    shutil.rmtree(os.path.join(root, folder), ignore_errors=True)


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


def legalize_species_label(species: Species,
                           return_label: bool = False,
                           check_arc_label: bool = True,
                           ) -> Optional[str]:
    """
    ARC uses the species label as the folder name on the server and the local machine.
    Make sure a label is legal, correct it if it's not.

    Args:
        species (Species): A species object.
        return_label (bool, optional): Whether to return the new label.
        check_arc_label (bool, optional): Whether to also check the label with ARC.

    Returns:
        str: The legalized species label.
    """
    for char in species.label:
        if char not in VALID_CHARS:
            species.label = species.molecule[0].get_formula()
            break
    else:
        if species.label[:2] == 'S(' and species.label[-1] == ')' \
                and all([char.isdigit() for char in species.label[2:-1]]):
            species.label = species.molecule[0].get_formula()
    if check_arc_label:
        species.label = check_label(species.label, verbose=False)[0]
    if return_label:
        return species.label
    return None


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


def get_species_with_qm_label(species: Species,
                              key: int,
                              ) -> Species:
    """
    Get a copy of the species with an updated QM label.
    This also adds the species to self.species if it's not already there.

    Args:
         species (Species): The species to consider.
         key (int): The respective species key, if exists.

    Returns:
        Species: A copy of the original species with a formatted QM species label.

    Todo:
        Add tests.
    """
    qm_species = species.copy(deep=False)
    legalize_species_label(species=qm_species)
    qm_species.label = f's{key}_{qm_species.label}'
    return qm_species
