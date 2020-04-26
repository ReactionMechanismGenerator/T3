"""
The RMG-ARC tandem module for iterative kinetic model generation and refinement

RMG is an automatic chemical mechanism generator. It is awesomely awesome.
ARC is an automatic electronic structure computation scheduler and processor. It is magnificently magnificent.
Here we iteratively execute RMG and ARC to generate and refine a chemical kinetic model.
Run RMG and ARC in tandem by passing an additional input file to RMG with the `-a` attribute, e.g.::

    python $rmgpy/rmg.py input.py -t3 arc.yml

where ``input.py`` is the legacy RMG input file, and ``arc.yml`` is an ARC input file augmented with T3 directives.

Do not specify species nor reactions in ARC's input file, they are automatically selected from RMG by this tool
according to the T3 directives in ARC's input file.

Directives by which species are selected for thermodynamic data calculations are all optional,
and may include the following (values in parentheses are defaults)::

    - SA observables: A list of observables for the SA. Entries are dictionaries, keys are 'label' and either 'smiles'
                      or 'adj' bearing the structure.
    - SA method (optional): The software to use for running SA: Either 'RMG' (default), 'RMS', or 'Cantera'.
    - SA threshold (optional): The SA threshold to use, 0.01 by default.
    - SA species (10): The top X species each observable is sensitive to for which thermo will be calculated.
    - SA reactions (10): The top X reactions each observable is sensitive to for which participating species thermo
                         will be calculated.
    - SA pdep threshold (0.001): The ratio of the sensitivity coefficient of a well relative to the maximum sensitivity
                                coefficient of all wells/TSs at each condition in a pressure-dependent network of a
                                reaction identified by ``SA reactions``, for which thermodynamic data will be calculated.
    - collision violators (``True``): Whether to calculate thermodynamic data for all species participating in
                                      collision rate violating reactions.
    - all core species (``False``): Whether to calculate thermodynamic data for all core species.

The above criteria are always matched with uncertainty analysis, currently in the form of a simple check of whether the
thermodynamic data was derived from group additivity or not.

Additional allowed directives are::

    - RMG tolerances ([0.05, 0.01, 0.001]): A list of RMG toleranceMoveToCore to use.
    - max tandem iterations (10): Maximum iteration T3 will make.
    - max RMG exceptions allowed (10): Maximum number of times RMG is allowed to crush.
    - max RMG walltime (01:00:00:00): The maximal wall time allocated per RMG execution.

T3 has a restart feature that automatically kicks in if it finds previous iteration directories in the output path.
To run T3 **afresh** in the same path, maker sure to first delete all iteration directories in that path.

Todo:
    - generate n generations of reactions (circles around the reactants) and get all thermo right
    - parse errors and warnings from the ARC output/status.yml file for failed species
    - modify tolerance interrupt simulation according to tolerance move to core
    - implement Cantera and RMS SA, including global observables
    - determine whether a species should be forbidden or not
    - scan pdep networks and the core, mark non-physical species
    - utilize the uncertainty analysis script
"""

import inspect
import os
import pandas as pd
import re
import shutil
import time
import traceback
from typing import TYPE_CHECKING, Optional, Tuple, Union

from pydas.daspk import DASPKError

from arkane import Arkane
from rmgpy import settings
from rmgpy.chemkin import load_chemkin_file
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.exceptions import \
    (ChemicallySignificantEigenvaluesError,
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
from rmgpy.rmg.main import initialize_log as initialize_rmg_log
from rmgpy.rmg.main import RMG
from rmgpy.rmg.pdep import PDepReaction
from rmgpy.solver.simple import SimpleReactor
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.species import Species
from rmgpy.thermo import NASAPolynomial, NASA, ThermoData, Wilhoit
from rmgpy.tools.loader import load_rmg_py_job
from rmgpy.tools.simulate import simulate

from arc.common import get_ordinal_indicator, key_by_val, read_yaml_file, time_lapse
from arc.main import ARC

from .config import species_labels_dict
from .logger import \
    (initialize_tandem_log,
     log,
     log_footer,
     log_species_summary,
     log_species_to_calculate,
     log_unconverged_species,
     )
from .utils import combine_dicts, delete_root_rmg_log, dict_to_str

if TYPE_CHECKING:
    from rmgpy.reaction import Reaction


ME_METHODS = ['CSE', 'MSC']  # master equation methods, could be any combination of 'CSE', 'RS', 'MSC'
MAX_ITERATIONS = 10  # this default value will be overridden by the input file's "max tandem iterations" argument
MAX_EXCEPTIONS = 10  # this default value will be overridden by the input file's "max RMG exceptions allowed" argument

t0 = None
rmg_exceptions_counter = 0
rmg_thermo_lib_base_path = os.path.join(settings['database.directory'], 'thermo', 'libraries')


def execute(args,
            kwargs,
            verbose: Optional[bool] = True,
            ) -> None:
    """
    Execute the tandem RMG and ARC model generation and refinement module.
    Note that ``args`` and ``kwargs`` should not be defined with the common `*` and `**` here.

    Args:
        args (parser.parse_args): The command line parsed arguments.
        kwargs (dict): Keyword arguments to be passed when initializing the RMG object.
        verbose (Optional[bool]): Whether or not to log to file.

    Raises:
        Input error if the number of RMG tolerances is greater than the total number of T3 iterations.
    """
    cwd = os.getcwd()
    # If the output directory is not yet set, then set it to the same directory as the input files by default
    if not args.output_directory:
        args.output_directory = cwd

    thermo_library = initialize_tandem_log(output_directory=args.output_directory, verbose=verbose)

    unconverged_species = list()    # entries are Species objects of species for which thermo calculations failed once.
    all_species = list()            # stores species_to_calc_in_iteration from all iterations, includes unconverged species
    species_dict = dict()           # keys are labels, values are dicts of {'spc': <RMG Species objects>, 'reason': ``str``}
    executed_networks = list()      # P-dep networks for which SA was already executed. Entries are tuples of isomer labels.

    rmg_input_file = os.path.join(cwd, args.rmg[0])
    with open(rmg_input_file, 'r') as f:
        content = f.read()
        has_sa = 'sensitivity=[' in content
        has_pdep = 'pressureDependence(' in content
    arc_input_file = os.path.join(cwd, args.arc[0])
    arguments, arc_input_dict = parse_arc_input_file(input_file_path=arc_input_file, has_sa=has_sa, has_pdep=has_pdep)

    log(f'\nUsing the following arguments:\n'
        f'{dict_to_str(arguments, level=1)}',
        verbose=verbose)

    max_iterations = arguments['max tandem iterations']  # default is 10, set in parse_arc_input_file()
    tolerances = arguments['RMG tolerances']
    if any(tolerances[i+1] >= tolerances[i] for i in range(len(tolerances) - 1)):
        log(f'The RMG tolerances are not in descending order, got: {tolerances}', 'warning', verbose=verbose)
    if len(tolerances) > max_iterations:
        raise InputError(f'The number of RMG tolerances ({len(tolerances)}) specified ({tolerances}) '
                         f'cannot be greater than the total number of T3 iterations requested ({max_iterations})')

    i_0, run_rmg_i_0, thermo_library = restart_t3(path=args.output_directory, thermo_library=thermo_library)
    i = -1
    additional_calcs_required = False

    # main RMG-ARC Tandem Tool loop
    for i in range(i_0, max_iterations):
        tolerance = tolerances[i] if len(tolerances) > i else tolerances[-1]
        log(f'\n\n\nT3 iteration {i}:\n'
            f'---------------\n',
            verbose=verbose)
        run_directory = os.path.join(args.output_directory, f'iteration_{i}')
        if not os.path.exists(run_directory):
            os.mkdir(run_directory)

        if i == i_0 and run_rmg_i_0 or i > i_0:
            run_rmg(input_file=rmg_input_file, output_directory=run_directory, kwargs=kwargs, arguments=arguments,
                    tolerance=tolerance, thermo_library=thermo_library)

        sa_success = run_sa(method=arguments['SA method'],
                            observable_list=arguments['SA observables'],
                            input_file=rmg_input_file,
                            run_directory=run_directory,
                            threshold=arguments['SA threshold'],
                            )
        if not sa_success:
            log(f'Could not complete the sensitivity analysis using {arguments["SA method"]}', 'error', verbose=verbose)
            if len(tolerances) > i + 1:
                # the next iteration will utilize a tighter tolerance
                continue
            elif i + 1 < max_iterations:
                log(f'Decreasing lowest tolerance to {0.5 * tolerances[-1]} and trying again', verbose=verbose)
                tolerances = tolerances + [0.5 * tolerances[-1]] * (max_iterations - len(tolerances))

        species_dict_in_iteration, additional_calcs_required, executed_networks = determine_species_to_calculate(
            run_directory=run_directory,
            arguments=arguments,
            unconverged_species=unconverged_species,
            all_species=all_species,
            iteration=i,
            executed_networks=executed_networks,
        )

        species_to_calc_in_iteration = [spc['spc'] for spc in species_dict_in_iteration.values()]
        all_species.extend(species_to_calc_in_iteration)
        species_dict.update(species_dict_in_iteration)

        if additional_calcs_required:
            run_arc(arc_input_dict, run_directory, species_to_calc_in_iteration)
            # add the calculated RMG libraries to the database and input file
            unconverged_species_in_iteration = get_unconverged_species(run_directory, all_species)
            unconverged_species.extend(unconverged_species_in_iteration)
            if len(species_to_calc_in_iteration) > len(unconverged_species_in_iteration):
                # we calculated something, add to library
                thermo_library = add_rmg_libraries(run_directory, thermo_library)
            else:
                additional_calcs_required = False
        if not additional_calcs_required and i >= len(tolerances):
            # T3 iterated through all of the user requested tolerances, and there are no more calculations required
            break

    if additional_calcs_required:
        # run RMG for the last time
        i += 1
        log(f'\n\n\nT3 iteration {i} (just generating a model using RMG):\n'
            f'---------------\n',
            verbose=verbose)
        run_directory = os.path.join(args.output_directory, f'iteration_{i}')
        if not os.path.exists(run_directory):
            os.mkdir(run_directory)
        run_rmg(input_file=rmg_input_file, output_directory=run_directory, kwargs=kwargs, arguments=arguments,
                tolerance=tolerances[i] if len(tolerances) > i else tolerances[-1], thermo_library=thermo_library)

    log_species_summary(species_dict, unconverged_species, verbose=verbose)
    log_footer(verbose=verbose)
    delete_root_rmg_log(path=args.output_directory)


def run_rmg(input_file: str,
            output_directory: str,
            kwargs: dict,
            arguments: dict,
            tolerance: float,
            thermo_library: Optional[str] = None,
            verbose: Optional[bool] = True,
            ) -> None:
    """
    Run RMG.

    Args:
        input_file (str): The path to the legacy RMG input file.
        output_directory (str): The path to the output directory for the current iteration.
        kwargs (dict): Keyword arguments to be passed when initializing the RMG object.
        arguments (dict): Various arguments for The RMG-ARC Tandem Tool.
        tolerance (float): The toleranceMoveToCore value to use.
        thermo_library (Optional[str]): The thermo library name created by the tandem tool.
                                        Consistent between iterations
                                        (it is None only in iteration 0, where it is being set).
        verbose (Optional[bool]): Whether or not to log to file.

    Raises:
        InputError: If something is wrong with the RMG input file.
        Various RMG Exceptions: if RMG crushed too many times.
    """
    global rmg_exceptions_counter
    global rmg_thermo_lib_base_path
    log(f'Running RMG (tolerance = {tolerance})...', verbose=verbose)
    # Initialize the logging system (resets the RMG.log file)
    initialize_rmg_log(verbose=20, log_file_name=os.path.join(output_directory, 'RMG.log'))  # 20 is `logging.INFO`
    max_num_exceptions_allowed = arguments['max RMG exceptions allowed']
    tic = time.time()
    rmg = RMG(input_file=input_file, output_directory=output_directory)
    if thermo_library is not None:
        rmg.thermo_libraries = [thermo_library]
    rmg.initialize(**kwargs)
    rmg.wallTime = arguments['max RMG walltime']
    rmg.model_settings_list[0].tol_move_to_core = tolerance
    if not rmg.model_settings_list[0].tol_keep_in_edge:
        # if not pruning (toleranceKeepInEdge != 0), have toleranceMoveToCore equal to toleranceInterruptSimulation
        rmg.model_settings_list[0].tol_interrupt_simulation = tolerance
        log(f'Running RMG (toleranceMoveToCore = {tolerance}, toleranceInterruptSimulation = {tolerance})...',
            verbose=verbose)
    else:
        # if pruning, then only change toleranceMoveToCore
        if rmg.model_settings_list[0].tol_interrupt_simulation < 1e8:
            log(f'Expected toleranceInterruptSimulation to be set higher than 1e8 when prunning in RMG, '
                f'got {rmg.model_settings_list[0].tol_interrupt_simulation}', 'warning', verbose=verbose)
        log(f'Running RMG (toleranceMoveToCore = {tolerance})...', verbose=verbose)

    try:
        rmg.execute(initialize=False)

    except (ChemicallySignificantEigenvaluesError, ChemkinError, CollisionError, CoreError, DASPKError, ILPSolutionError,
            InvalidMicrocanonicalRateError, KineticsError, ModifiedStrongCollisionError, NetworkError,
            PressureDependenceError, ReactionError, ReservoirStateError, StatmechError, StatmechFitError) as e:
        log(f'RMG Errored with {e.__class__}. Got the following trace:', level='error', verbose=verbose)
        log(traceback.format_exc(), verbose=verbose)

        if rmg_exceptions_counter > max_num_exceptions_allowed:
            log(f'This is the {rmg_exceptions_counter} exception raised by RMG. Not allowing additional exceptions. '
                f'Terminating the process', level='error', verbose=verbose)
            raise
        else:
            log(f'This is the {rmg_exceptions_counter} exception raised by RMG '
                f'(maximum number of allowed exceptions is {max_num_exceptions_allowed})',
                level='warning', verbose=verbose)
        rmg_exceptions_counter += 1

    except InputError:
        log('Something seems to be wrong with the RMG input file:', level='error', verbose=verbose)
        log('\n\n', verbose=verbose)
        log(input_file, verbose=verbose)
        raise

    elapsed_time = time_lapse(tic)
    if verbose:
        log(f'RMG terminated. Overall execution time: {elapsed_time}', verbose=verbose)


def run_arc(input_dict: Union[str, dict],
            run_directory: str,
            species_to_calc: Optional[list] = None,
            verbose: Optional[bool] = True,
            ) -> None:
    """
    Run ARC.

    Args:
        input_dict (str, dict): A dictionary containing directives to execute ARC, or the path to the yml file.
        run_directory (str): A path to the RMG-ARC iteration directory.
        species_to_calc (Optional[list]): Entries are RMG Species for which ARC will calculate thermo properties.
        verbose (Optional[bool]): Whether or not to log to file.

    Raises:
        TypeError: If ``input_dict`` is of wrong type.
        Various ARC Exceptions: if ARC crashed.
    """
    if isinstance(input_dict, str):
        input_dict = read_yaml_file(input_dict)
    if not isinstance(input_dict, dict):
        raise TypeError(f'The input dictionary must be a dictionary or a path to an ARC input file.\n'
                        f'Got: {input_dict} which is a {type(input_dict)}.')
    species_to_calc = species_to_calc or list()
    log('\nRunning ARC...', verbose=verbose)
    input_dict['project_directory'] = os.path.join(run_directory, 'ARC')
    if not os.path.exists(input_dict['project_directory']):
        os.mkdir(input_dict['project_directory'])
    if len(species_to_calc):
        input_dict['arc_species_list'] = species_to_calc
    if 'project' not in input_dict:
        input_dict['project'] = 't3'
    tic = time.time()
    arc0 = ARC(**input_dict)
    try:
        arc0.execute()
    except Exception as e:
        log(f'ARC crashed with {e.__class__}. Got the following message:\n{e}', level='error', verbose=verbose)
        raise
    elapsed_time = time_lapse(tic)
    log(f'ARC terminated. Overall execution time: {elapsed_time}', verbose=verbose)


def run_sa(method: str,
           observable_list: list,
           run_directory: str,
           input_file: str,
           threshold: Optional[float] = 0.001,
           verbose: Optional[bool] = True,
           ) -> bool:
    """
    Run a sensitivity analysis.

    Args:
        method (str): The software to use, either RMG, RMS, or Cantera.
        observable_list (list): Entries are dictionaries of 'label' and structure (either 'smiles' or 'adj').
        run_directory (str): A path to the RMG-ARC iteration directory.
        input_file (str): The path to the legacy RMG input file.
        threshold (Optional[float]): The sensitivity threshold to use.
        verbose (Optional[bool]): Whether or not to log to file.

    Raises:
        InputError: If ``method`` has a wrong value.

    Returns:
        bool: Whether SA ran successfully. ``True`` if it did.
    """
    if method.lower() not in ['rmg', 'rms', 'cantera']:
        raise InputError(f'The "SA method" argument must equal to either "RMG", "RMS", or "Cantera". Got: {method}')
    if method.lower() == 'rmg':
        log('Running SA using RMG...', verbose=verbose)
        model = os.path.join(run_directory, 'chemkin', 'chem_annotated.inp')
        species_dict = os.path.join(run_directory, 'chemkin', 'species_dictionary.txt')
        sa_path = os.path.join(run_directory, 'sa')
        if not os.path.isdir(sa_path):
            os.mkdir(sa_path)
        rmg_sa_input_file = os.path.join(sa_path, 'input.py')
        if os.path.isfile(rmg_sa_input_file):
            os.remove(rmg_sa_input_file)
        shutil.copyfile(src=input_file, dst=rmg_sa_input_file)

        rmg = load_rmg_py_job(input_file=rmg_sa_input_file,
                              chemkin_file=model,
                              species_dict=species_dict,
                              generate_images=True,
                              use_chemkin_names=False,
                              check_duplicates=False)

        rmg_species = rmg.reaction_model.core.species
        rmg_observable_species = list()
        for observable in observable_list:
            for rmg_spc in rmg_species:
                if 'adj' in observable:
                    observable_spc = Species(label=observable['label']).from_adjacency_list(observable['adj'])
                elif 'smiles' in observable:
                    observable_spc = Species(label=observable['label']).from_smiles(observable['smiles'])
                else:
                    raise InputError(f'All SA observables must have structure (smiles or adj), got: {observable}')
                if observable_spc.label == rmg_spc.label or observable_spc.is_isomorphic(rmg_spc):
                    rmg_observable_species.append(rmg_spc)
                    break
            else:
                raise InputError(f'Could not find the observable species {observable["label"]} '
                                 f'in the RMG species list.')

        for reaction_system in rmg.reaction_systems:
            if isinstance(reaction_system, SimpleReactor):
                reaction_system.sensitive_species = rmg_observable_species
                reaction_system.sensitivity_threshold = threshold
                if hasattr(reaction_system, 'Trange') and reaction_system.Trange is not None:
                    temperature = sum([t.value_si for t in reaction_system.Trange]) / len(reaction_system.Trange)
                else:
                    temperature = reaction_system.T.value_si
                reaction_system.sens_conditions['T'] = temperature
                if hasattr(reaction_system, 'Prange') and reaction_system.Prange is not None:
                    pressure = sum([p.value_si for p in reaction_system.Prange]) / len(reaction_system.Prange)
                else:
                    pressure = reaction_system.P.value_si
                reaction_system.sens_conditions['P'] = pressure
            elif isinstance(reaction_system, LiquidReactor):
                reaction_system.sensitive_species = rmg_observable_species
                reaction_system.sensitivity_threshold = threshold
                if hasattr(reaction_system, 'Trange') and reaction_system.Trange is not None:
                    temperature = sum([t.value_si for t in reaction_system.Trange]) / len(reaction_system.Trange)
                else:
                    temperature = reaction_system.T.value_si
                reaction_system.sens_conditions['T'] = temperature
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
        raise NotImplementedError('Currently only RMG is implemented as an SA method')  # temp
    return True


def restart_t3(path: str,
               thermo_library: Optional[str] = None,
               verbose: Optional[bool] = True,
               ) -> Tuple[int, bool, str]:
    """
    Restart T3 by looking for existing iteration folders.

    Args:
        path (str): The folder path of the current T3 run, where the 'iteration_i' folders are located.
        thermo_library (Optional[str]): The name of the thermo library parsed from the previous T3 log file.
        verbose (Optional[bool]): Whether or not to log to file.

    Returns:
        Tuple[int, bool, str]:
            - The current iteration number.
            - Whether to run RMG for this iteration.
            - The thermo library name to use.
    """
    # default values if not restarting:
    i_max = 0
    run_rmg_i, restart_arc_i = True, False
    folders = tuple(os.walk(path))[0][1]  # returns a 3-tuple: (dirpath, dirnames, filenames)
    folders = [folder for folder in folders if 'iteration_' in folder]
    if folders:
        i_max = max(int(folder.split('_')[1]) for folder in folders)  # get the recent iteration number
        run_directory = f'iteration_{i_max}'
        rmg_log_path = os.path.join(run_directory, 'RMG.log')
        arc_log_path = os.path.join(run_directory, 'ARC', 'arc.log')
        arc_restart_path = os.path.join(run_directory, 'ARC', 'restart.yml')
        if os.path.isfile(rmg_log_path):
            with open(rmg_log_path, 'r') as f:
                lines = f.readlines()
                for line in lines[::-1]:
                    if 'MODEL GENERATION COMPLETED' in line:
                        # RMG terminated, no need to regenerate the model
                        run_rmg_i = False
                        break
        if os.path.isfile(arc_log_path):
            with open(arc_log_path, 'r') as f:
                lines = f.readlines()
                for line in lines[::-1]:
                    if 'All jobs terminated. Summary for project' in line:
                        # ARC terminated as well, continue to the next iteration
                        i_max += 1
                        break
                else:
                    # ARC did not terminate, see if the restart file was generated
                    if os.path.isfile(arc_restart_path):
                        restart_arc_i = True
        elif i_max == 0:
            # This is the first iteration, and ARC did not converge, make sure thermo_library is None
            thermo_library = None
        if i_max or not run_rmg_i or restart_arc_i:
            rmg_text = ', using the completed RMG run from this iteration' if not run_rmg_i \
                else ', re-running RMG for this iteration'
            arc_text = ', restarting the previous ARC run in this iteration' if restart_arc_i else ''
            log(f'\n\n\nRestarting T3 from iteration {i_max}{rmg_text}{arc_text}.\n', verbose=verbose)
        if restart_arc_i:
            run_arc(input_dict=arc_restart_path,
                    run_directory=run_directory,
                    species_to_calc=list())

    else:
        # no iteration folders exist, make sure thermo_library is None
        thermo_library = None
    return i_max, run_rmg_i, thermo_library


def parse_arc_input_file(input_file_path: Union[str, list],
                         has_sa: bool,
                         has_pdep: bool,
                         verbose: Optional[bool] = True,
                         ) -> Tuple[dict, dict]:
    """
    Split the ARC input file into relevant arguments for The RMG-ARC Tandem Tool, and the legacy ARC input.
    This function also checks the ARC input file's validity.

    All Tandem arguments are optional, and may include the following (values in parentheses are defaults):
    'SA observables' (``True``), 'SA species' (10), 'SA reactions' (10), 'SA pdep threshold' (10%),
    'collision violators' (``True``), 'all core species' (``False``), 'RMG tolerances' ([0.1, 0.01, 0.001]),
    'max tandem iterations' (10), 'max RMG exceptions allowed' (10), 'max RMG walltime' (01:00:00:00).

    Args:
        input_file_path (str, list): A path to ARC's input file.
        has_sa (bool): Whether the RMG input file contains sensitivity analysis directives, ``True`` if it does.
        has_pdep (bool): Whether the RMG input file contains a P-dep block, ``True`` if it does.
        verbose (Optional[bool]): Whether or not to log to file.

    Raises:
        InputError: If the ARC input file contains either species or reactions.
        ValueError: If arguments['SA pdep threshold'] has a value outside of the range [0, 1].
        TypeError: If arguments have wrong type.

    Returns:
        Tuple[dict, dict]:
            - The T3 arguments.
            - ARC's input.
    """
    input_file_path = input_file_path[0] if isinstance(input_file_path, list) else input_file_path
    input_dict = read_yaml_file(input_file_path)
    arguments = dict()

    arguments['SA observables'] = input_dict['SA observables'] if 'SA observables' in input_dict else list()
    arguments['SA method'] = input_dict['SA method'] if 'SA method' in input_dict else 'RMG'
    arguments['SA threshold'] = input_dict['SA threshold'] if 'SA threshold' in input_dict else 0.001
    arguments['SA species'] = input_dict['SA species'] if 'SA species' in input_dict else 10 and has_sa
    arguments['SA reactions'] = input_dict['SA reactions'] if 'SA reactions' in input_dict else 10 and has_sa
    arguments['SA pdep threshold'] = input_dict['SA pdep threshold'] \
        if 'SA pdep threshold' in input_dict else 0.01 if has_pdep else 1.0
    arguments['collision violators'] = input_dict['collision violators'] \
        if 'collision violators' in input_dict else True
    arguments['all core species'] = input_dict['all core species'] if 'all core species' in input_dict else False
    arguments['RMG tolerances'] = input_dict['RMG tolerances'] if 'RMG tolerances' in input_dict else [0.05, 0.01, 0.001]
    arguments['max tandem iterations'] = input_dict['max tandem iterations'] \
        if 'max tandem iterations' in input_dict else MAX_ITERATIONS
    arguments['max RMG exceptions allowed'] = input_dict['max RMG exceptions allowed'] \
        if 'max RMG exceptions allowed' in input_dict else MAX_EXCEPTIONS
    arguments['max RMG walltime'] = input_dict['max RMG walltime'] \
        if 'max RMG walltime' in input_dict else '01:00:00:00'

    # check argument types:
    if not isinstance(arguments['SA observables'], list):
        raise TypeError(f'The SA observables argument must be a list, got {arguments["SA observables"]} '
                        f'which is a {type(arguments["SA observables"])}')
    if not isinstance(arguments['SA method'], str):
        raise TypeError(f'The SA method argument must be a string, got {arguments["SA method"]} '
                        f'which is a {type(arguments["SA method"])}')
    if not isinstance(arguments['SA threshold'], float):
        raise TypeError(f'The SA threshold argument must be a float, got {arguments["SA threshold"]} '
                        f'which is a {type(arguments["SA threshold"])}')
    if not isinstance(arguments['SA species'], int):
        raise TypeError(f'The SA species argument must be an integer, got {arguments["SA species"]} '
                        f'which is a {type(arguments["SA species"])}')
    if not isinstance(arguments['SA pdep threshold'], float):
        if isinstance(arguments['SA pdep threshold'], int) and 0 <= arguments['SA pdep threshold'] <= 1:
            arguments['SA pdep threshold'] = float(arguments['SA pdep threshold'])
        else:
            raise TypeError(f'The SA pdep threshold argument must be a float, got {arguments["SA pdep threshold"]} '
                            f'which is a {type(arguments["SA pdep threshold"])}')
    if arguments['SA pdep threshold'] > 1 or arguments['SA pdep threshold'] < 0:
        raise ValueError(f'The SA pdep threshold argument cannot be negative or greater than one. '
                         f'Got {arguments["SA pdep threshold"]}')
    if not isinstance(arguments['collision violators'], bool):
        raise TypeError(f'The collision violators argument must be a boolean, got {arguments["collision violators"]} '
                        f'which is a {type(arguments["collision violators"])}')
    if not isinstance(arguments['all core species'], bool):
        raise TypeError(f'The all core species argument must be a boolean, got {arguments["all core species"]} '
                        f'which is a {type(arguments["all core species"])}')
    if not isinstance(arguments['RMG tolerances'], (list, tuple)):
        raise TypeError(f'The RMG tolerances argument must be a list, got {arguments["RMG tolerances"]} '
                        f'which is a {type(arguments["RMG tolerances"])}')
    if not isinstance(arguments['max tandem iterations'], int):
        raise TypeError(f'The max tandem iterations argument must be an integer, '
                        f'got {arguments["max tandem iterations"]} '
                        f'which is a {type(arguments["max tandem iterations"])}')
    if not isinstance(arguments['max RMG exceptions allowed'], int):
        raise TypeError(f'The max RMG exceptions allowed argument must be an integer, '
                        f'got {arguments["max RMG exceptions allowed"]} '
                        f'which is a {type(arguments["max RMG exceptions allowed"])}')
    if not isinstance(arguments['SA method'], str):
        raise TypeError(f'The max RMG walltime argument must be a string (e.g., "00:05:00:00" for 5 hrs), '
                        f'got {arguments["max RMG walltime"]} '
                        f'which is a {type(arguments["max RMG walltime"])}')

    if 'RMG tolerances' not in input_dict:
        log(f'Using the following default RMG toleranceMoveToCore tolerances: {arguments["RMG tolerances"]}. '
            f'Specify the "RMG tolerances" list in the ARC input file to override this.',
            level='warning', verbose=verbose)

    # remove all T3 arguments from ARC's input
    for argument in arguments.keys():
        if argument in input_dict:
            del input_dict[argument]

    # check that the resulting ARC input dictionary is valid
    if 'species' in input_dict:
        raise InputError("The 'species' dictionary cannot be present in the ARC input file when running in tandem "
                         "with RMG. Correct ARC's input file and run again.")
    if 'reactions' in input_dict:
        raise InputError("The 'reactions' dictionary cannot be present in the ARC input file when running in tandem "
                         "with RMG. Correct ARC's input file and run again.")
    for key in list(input_dict.keys()):
        if key not in inspect.getfullargspec(ARC.__init__).args:
            # This argument was not extracted above, and it's not an ARC argument, remove so ARC doesn't crush
            log(f'Argument "{key}" passed to ARC is not allowed. Not using it.\n'
                f'(if this was meant to be a T3 argument, it was not recognized)',
                level='error', verbose=verbose)
            del input_dict[key]
    return arguments, input_dict


def set_legal_species_labels(species_to_calc: list,
                             all_species: list,
                             ) -> list:
    """
    ARC uses the species label as folder names on the servers and the local machine.
    Make sure each species has a legal and unique label (which consists of the molecular formula, underscore, and index)
    Store the new labels (keys) and original Chemkin labels (values) in the ``species_labels_dict`` dictionary.

    Args:
        species_to_calc (list): Entries are RMG Species objects for which thermo should be calculated in this iteration.
        all_species (list): Entries are RMG Species object of previously considered species,
                            used here to determine unique labels.

    Returns:
        list: Species to be calculated with modified labels.
    """
    global species_labels_dict
    for i in range(len(species_to_calc)):
        # we're changing arguments (labels) within the `species_to_calc` list, so iterate by index
        formula = species_to_calc[i].molecule[0].get_formula()
        existing_indices = list()
        for spc in all_species + species_to_calc:
            # all_species does not include the entries in species_to_calc at this point
            if '_' in spc.label:
                splits = spc.label.split('_')
                if len(splits) == 2 and splits[0] == formula and all([char.isdigit() for char in splits[1]]):
                    existing_indices.append(int(splits[1]))
        index = max(existing_indices) + 1 if len(existing_indices) else 0
        new_label = f'{formula}_{index}'
        species_labels_dict[new_label] = species_to_calc[i].to_chemkin()  # store the original Chemkin label
        species_to_calc[i].label = new_label  # reset the label
    return species_to_calc


def species_not_in_list(species: Union[Species, str],
                        species_list: list,
                        ) -> bool:
    """
    Check whether a species is NOT in species_list. Use the unique Chemkin species labels.
    Note that the labels in species_list were changed, therefore the global ``species_labels_dict`` is used.

    Args:
        species (Species, str): Either an RMG Species object or a Chemkin label thereof.
        species_list (list): Entries are RMG Species objects.

    Raises:
        TypeError: If ``species`` is of wrong type.

    Returns:
        bool: Whether the species is in the list. True if it is.
    """
    global species_labels_dict
    if isinstance(species, str):
        label = species
    elif isinstance(species, Species):
        label = species.to_chemkin()
    else:
        raise TypeError(f'species can be either an RMG Species object or a Chemkin label thereof.\n'
                        f'Got {species} which is a {type(species)}')
    for spc in species_list:
        if label == spc.label:
            return False
        elif spc.label in list(species_labels_dict.keys()) and label == species_labels_dict[spc.label]:
            return False
    return True


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


def get_species_by_label(label: str,
                         species_list: list,
                         ) -> Union[Species, None]:
    """
    Get a species from a list of species by its label.

    Args:
        label (str): A species' label (species.label).
        species_list (list): Entries are RMG Species objects.

    Raises:
        InputError: If ``label`` is None.

    Returns:
        Union[Species, None]: The corresponding species from the species_list.
                              Returns ``None`` if no species was found.
    """
    if label is None:
        raise InputError('Got None as label input')
    for spc in species_list:
        if spc.label == label or spc.to_chemkin() == label:
            return spc
    if '(' in label and ')' in label:
        # try by the RMG species index
        for spc in species_list:
            if spc.index == int(label.split('(')[-1].split(')')[0]):
                return spc
    return None


def get_reaction_by_index(index: int,
                          reaction_list: list,
                          ) -> Union['Reaction', None]:
    """
    Get a reaction from a list of reactions by its index.

    Args:
        index (int): The reaction index attribute, 0-indexed.
        reaction_list (list): Entries are RMG Reaction objects.

    Returns:
        Union['Reaction', None]: The corresponding reaction from the reaction_list.
                               Returns ``None`` if no reaction was found.
    """
    for rxn in reaction_list:
        if rxn.index == index:
            return rxn
    return None


def calc_based_on_thermo_comment(species: Union[Species, None]) -> bool:
    """
    A helper function for reading the species thermo comment and determining whether to calculate thermodynamic
    properties for it using ARC.

    Args:
        species (Union[Species, None]): The considered species.

    Raises:
        InputError: If ``species`` is None.

    Returns:
        bool: Whether thermodynamic properties should be calculated. ``True`` if they should.
    """
    if species is None:
        raise InputError('Got None as species input')
    if 'group additivity' in species.thermo.comment or '+ radical(' in species.thermo.comment:
        return True
    return False


def has_high_uncertainty(species: Union[Species, None],
                         unconverged_species: list,
                         species_to_calc: dict,
                         ) -> bool:
    """
    Determine whether a species' thermochemical properties should be calculated based on the uncertainty in its
    thermodynamic properties (currently only considering GAV).
    Also check that it is not already selected or previously attempted to be calculated but did not converge.

    Todo:
        Currently only checks whether the species was calculated using GAV, should be elaborated to use uncertainty.

    Args:
        species (Union[Species, None]): The species for which the query is performed.
        unconverged_species (list): Entries are RMG Species objects for which thermo calculations did not converge
                                    in ARC throughout the iterations.
        species_to_calc (dict): RMG Species to calculate thermo for in the current iteration.
                                Keys are labels, values are dicts of {'spc': RMG Species objects, 'reason': ``str``}.

    Raises:
        InputError: If ``species`` is None.

    Returns:
        bool: Whether thermochemical properties for the species should be calculated. ``True`` if they should be.
    """
    if species is None:
        raise InputError('Got None as species input')
    if calc_based_on_thermo_comment(species) \
            and not any(species.is_isomorphic(other) for other in unconverged_species +
                        [spc['spc'] for spc in species_to_calc.values()]):
        return True
    return False


def load_species_and_reactions_from_chemkin_file(run_directory: str,
                                                 verbose: Optional[bool] = True,
                                                 ) -> Tuple[list, list]:
    """
    Load RMG Species and Reaction objects from the annotated Chemkin file.

    Args:
        run_directory (str): The current RMG-ARC Tandem Tool iteration directory.
        verbose (Optional[bool]): Whether or not to log to file.

    Raises:
        ChemkinError: If the Chemkin file could not be read.

    Returns:
        Tuple[list, list]:
            - The Species objects.
            - The Reactions objects.
    """
    chemkin_path = os.path.join(run_directory, 'chemkin', 'chem_annotated.inp')
    dictionary_path = os.path.join(run_directory, 'chemkin', 'species_dictionary.txt')
    try:
        rmg_species, rmg_reactions_ = load_chemkin_file(chemkin_path, dictionary_path, check_duplicates=True)
    except ChemkinError:
        log(f'Could not read the Chemkin file {chemkin_path}!\n'
            f'Trying to read it without checking for duplicate reactions...',
            level='error', verbose=verbose)
        try:
            rmg_species, rmg_reactions_ = load_chemkin_file(chemkin_path, dictionary_path, check_duplicates=False)
        except ChemkinError:
            log(f'Could still not read the Chemkin file {chemkin_path}!',
                level='error', verbose=verbose)
            raise
        else:
            log(f'Read the Chemkin file\n{chemkin_path}\nwithout checking for duplicate reactions, '
                f'SA results might be inaccurate.', verbose=verbose)
    rmg_reactions = list()
    for i, reaction in enumerate(rmg_reactions_):
        # renumber, since duplicate reactions are removed and leave gaps in the index
        # also, now the index should match the RMG index - 1 rather than the Chemkin index - 1 (it is 0-indexed)
        reaction.index = i
        rmg_reactions.append(reaction)
    return rmg_species, rmg_reactions


def determine_species_to_calculate(run_directory: str,
                                   arguments: dict,
                                   unconverged_species: list,
                                   all_species: list,
                                   iteration: int,
                                   executed_networks: list,
                                   verbose: Optional[bool] = True,
                                   ) -> Tuple[dict, bool, list]:
    """
    Determine which species in the executed RMG job should be calculated by ARC.
    Species which were previously attempted to be calculated but did not converge will not be selected again.

    Args:
        run_directory (str): The current RMG-ARC Tandem Tool iteration directory.
        arguments (dict): User directives for determining which species thermodynamic properties to calculate.
        unconverged_species (list): Entries are RMG Species for which thermo calculations did not converge
                                    in ARC throughout the iterations.
        all_species (list): Entries are RMG Species identified so far to be calculated. Used for setting legal names.
        iteration (int): The current T3 iteration number.
        executed_networks (list): Entries are tuples of isomers labels from networks which were already executed.
        verbose (Optional[bool]): Whether or not to log to file.

    Returns:
        Tuple[dict, bool, list]:
            - Keys are species Chemkin labels, values are dictionaries.
              In the value dictionaries, keys are 'spc' and 'reason', values are RMG Species for which thermodynamic
              properties should be calculated, and the reason they were selected for the calculation, respectively.

            - Whether additional calculations are required. ``True`` if they are.

            - The updated executed_networks list.
    """
    species_to_calc = dict()  # keys are labels, values are dicts of {'spc': RMG Species objects, 'reason': ``str``}
    species_to_calc_sa, species_to_calc_coll = None, None

    rmg_species, rmg_reactions = load_species_and_reactions_from_chemkin_file(run_directory)
    log(f'This RMG model has {len(rmg_species)} species and {len(rmg_reactions)} reactions in its core.',
        verbose=verbose)

    if arguments['all core species']:
        # calculate all core species if needed based on their thermo comment (don't bother checking other directives)
        for species in rmg_species:
            if has_high_uncertainty(species, unconverged_species, species_to_calc):
                species_to_calc[species.to_chemkin()] = {'spc': species, 'reason': 'All core species'}
    else:

        # SA
        if arguments['SA observables'] or arguments['SA species'] or arguments['SA reactions'] \
                or arguments['SA pdep threshold'] < 1:
            species_to_calc_sa, executed_networks = determine_species_based_on_sensitivity(
                run_directory, arguments, rmg_species, rmg_reactions, unconverged_species, iteration, executed_networks)

        # collision violators
        if arguments['collision violators']:
            species_to_calc_coll = determine_species_based_on_collision_violators(
                run_directory, rmg_species, unconverged_species)

        species_to_calc = combine_dicts(species_to_calc_sa, species_to_calc_coll)

    additional_calcs_required = bool(len(species_to_calc.values()))
    log(f'Additional calculations required: {additional_calcs_required}\n',
        verbose=verbose)
    if additional_calcs_required:
        log_species_to_calculate(species_to_calc, verbose=verbose)
        set_legal_species_labels([spc['spc'] for spc in species_to_calc.values()], all_species)
    return species_to_calc, additional_calcs_required, executed_networks


def determine_species_based_on_sensitivity(run_directory: str,
                                           arguments: dict,
                                           rmg_species: list,
                                           rmg_reactions: list,
                                           unconverged_species: list,
                                           iteration: int,
                                           executed_networks: list,
                                           verbose: Optional[bool] = True,
                                           ) -> Tuple[dict, list]:
    """
    Determine species to calculate based on sensitivity analysis.

    Args:
        run_directory (str): The current RMG-ARC Tandem Tool iteration directory.
        arguments (dict): User directives for determining which species thermodynamic properties to calculate.
        rmg_species (list): All RMG Species parsed from the Chemkin file.
        rmg_reactions (list): All RMG Reactions parsed from the Chemkin file.
        unconverged_species (list): Entries are RMG Species for which thermo calculations did not converge
                                    in ARC throughout the iterations.
        iteration (int): The current T3 iteration number.
        executed_networks (list): Entries are tuples of isomers labels from networks which were already executed.
        verbose (Optional[bool]): Whether or not to log to file.

    Returns:
        Tuple[dict, list]:
            - Keys are Chemkin labels, values are dicts of {'spc': RMG Species objects, 'reason': ``str``}

            - The updated executed_networks list.
    """
    species_to_calc = dict()  # keys are labels, values are dicts of {'spc': RMG Species objects, 'reason': ``str``}
    pdep_rxns_to_explore = list()  # Entries are (pressure dependent reaction, i, observable_label) tuples
    sa_path = os.path.join(run_directory, 'sa', 'solver')
    if not os.path.exists(sa_path):
        log("Could not find the path to RMG's solver output folder. Not executing "
            "calculations in ARC based on sensitivity analysis!",
            level='error', verbose=verbose)
        return dict(), list()

    sa_files = list()
    for file_ in os.listdir(sa_path):
        if 'sensitivity' in file_ and file_.endswith(".csv"):
            sa_files.append(file_)

    for sa_file in sa_files:
        # iterate through all SA .csv files in the solver folder
        df = pd.read_csv(os.path.join(sa_path, sa_file))
        sa_dict = {'rxn': dict(), 'spc': dict()}
        for header in df.columns:
            # iterate through all headers in the SA .csv file, but skip the `Time (s)` column
            sa_type = None
            if 'dln[k' in header and arguments['SA reactions']:
                sa_type = 'rxn'
            elif 'dG' in header and arguments['SA species']:
                sa_type = 'spc'
            if sa_type is not None:
                # proceed only if we care about this column
                entry = dict()
                # check whether the observable requires calculations:
                observable_label = header.split('[')[1].split(']')[0]
                observable = get_species_by_label(observable_label, rmg_species)
                if observable is None:
                    log(f'Could not identify observable species {observable_label}!',
                        level='error', verbose=verbose)
                if arguments['SA observables'] \
                        and has_high_uncertainty(observable, unconverged_species, species_to_calc):
                    species_to_calc[observable.to_chemkin()] = {'spc': observable, 'reason': 'observable'}
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
            for i in range(min(arguments['SA reactions'], len(sa_list_sorted))):
                reaction = get_reaction_by_index(sa_list_sorted[i]['parameter'] - 1, rmg_reactions)
                for species in reaction.reactants + reaction.products:
                    if has_high_uncertainty(species, unconverged_species, species_to_calc):
                        num = f'{i+1}{get_ordinal_indicator(i+1)} ' if i else ''
                        reason = f'(iteration {iteration}) participates in the {num}most sensitive reaction ' \
                                 f'for {observable_label}: {reaction}'
                        species_to_calc[species.to_chemkin()] = {'spc': species, 'reason': reason}
                if reaction.kinetics.is_pressure_dependent() \
                        and reaction not in [rxn_tup[0] for rxn_tup in pdep_rxns_to_explore] \
                        and arguments['SA pdep threshold'] < 1:
                    pdep_rxns_to_explore.append((reaction, i, observable_label))
        for observable_label, sa_list in sa_dict['spc'].items():
            sa_list_sorted = sorted(sa_list, key=lambda item: item['max_sa'], reverse=True)
            for i in range(min(arguments['SA species'], len(sa_list_sorted))):
                species = get_species_by_label(sa_list_sorted[i]['parameter'], rmg_species)
                if species is None:
                    log(f'Could not identify species {sa_list_sorted[i]["parameter"]}!',
                        level='error', verbose=verbose)
                if has_high_uncertainty(species, unconverged_species, species_to_calc):
                    num = f'{i+1}{get_ordinal_indicator(i+1)} ' if i else ''
                    reason = f'(iteration {iteration}) the {num}most sensitive species thermo for {observable_label}'
                    species_to_calc[species.to_chemkin()] = {'spc': species, 'reason': reason}

    species_to_calc, executed_networks = determine_species_from_pdep_network(run_directory=run_directory,
                                                                             pdep_rxns_to_explore=pdep_rxns_to_explore,
                                                                             unconverged_species=unconverged_species,
                                                                             species_to_calc=species_to_calc,
                                                                             iteration=iteration,
                                                                             threshold=arguments['SA pdep threshold'],
                                                                             executed_networks=executed_networks,
                                                                             rmg_species=rmg_species)

    return species_to_calc, executed_networks


def determine_species_from_pdep_network(run_directory: str,
                                        pdep_rxns_to_explore: list,
                                        unconverged_species: list,
                                        species_to_calc: dict,
                                        iteration: int,
                                        threshold: float,
                                        executed_networks: list,
                                        rmg_species: list,
                                        verbose: Optional[bool] = True,
                                        ) -> Tuple[dict, list]:
    """
    Determine species to calculate based on a pressure dependent network by spawning a network sensitivity analysis.
    First tries to use the CSE method, if unsuccessful tries the MSC method.

    Args:
        run_directory (str): The current RMG-ARC Tandem Tool iteration directory.
        pdep_rxns_to_explore (list): Entries are tuples of (Reaction, SA rank index, observable_label).
        unconverged_species (list): Entries are RMG Species for which thermo calculations did not converge
                                    in ARC throughout the iterations.
        species_to_calc (dict): Keys are Chemkin labels, values are dictionaries of {'spc': species, 'reason': reason}.
        iteration (int): The T3 iteration number.
        threshold (float): The relative SA pdep threshold above which sensitivity coefficients are considered high.
        executed_networks (list): Entries are tuples of isomers labels from networks which were already executed.
        rmg_species (list): All RMG Species parsed from the Chemkin file.
        verbose (Optional[bool]): Whether or not to log to file.

    Returns:
        Tuple[dict, list]:
            - The updated species_to_calc dictionary.

            - The updated executed_networks list.
    """
    if not os.path.isdir(os.path.join(run_directory, 'pdep_sa')):
        os.mkdir(os.path.join(run_directory, 'pdep_sa'))
    for reaction_tuple in pdep_rxns_to_explore:
        reaction = reaction_tuple[0]
        networks_path = os.path.join(run_directory, 'pdep')
        if isinstance(reaction, PDepReaction) and 0 < threshold < 1:
            # consider wells on the PES to which this rate coefficient is sensitive

            # identify the network name and file name
            network_file_names = list()
            for (_, _, files) in os.walk(networks_path):
                network_file_names.extend(files)
                break  # don't continue to explore subdirectories
            network_file_names = [network_file_name for network_file_name in network_file_names
                                  if f'network{reaction.network.index}_' in network_file_name]
            if not network_file_names:
                # this PDepReaction did not stem from a network file, it is probably a library reaction
                continue
            network_version = max([int(network_file_name.split('.')[0].split('_')[1])
                                   for network_file_name in network_file_names])
            network_name = f'network{reaction.network.index}_{network_version}'  # w/o the '.py' extension

            # try running this network using either the CSE or the MSC method
            sa_coefficients_path = None
            errors = list()
            arkane = None
            for method in ME_METHODS:
                input_file_path, output_file_path, isomer_labels = modify_pdep_network_file(
                    run_directory, network_name, method=method)
                arkane = Arkane(input_file=input_file_path, output_directory=os.path.dirname(input_file_path))
                arkane.plot = True
                log(f'\nRuning Pdep SA for network {network_name} using the {method} method...', verbose=verbose)
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
                    log(f'Successfully executed a Pdep SA for network {network_name} using the {method} method.\n',
                        verbose=verbose)
                    executed_networks.append(isomer_labels)
                    sa_coefficients_path = output_file_path
                    break
            else:
                log(f'Could not execute a Pdep SA for network {network_name} using either of the {ME_METHODS} '
                    f'methods. Got the following errors:',
                    level='error', verbose=verbose)
                for method, e in zip(ME_METHODS, errors):
                    log(f'{e.__class__} for method {method}: {e}\n', verbose=verbose)

            if sa_coefficients_path is not None:
                sa_dict = read_yaml_file(sa_coefficients_path)
                # get_species_label_by_structure(adj, species_list=all_species)
                reactants_label = ' + '.join([reactant.to_chemkin() for reactant in reaction.reactants])
                products_label = ' + '.join([product.to_chemkin() for product in reaction.products])
                chemkin_reaction_str = f'{reactants_label} <=> {products_label}'
                labels_map = dict()  # keys are network species labels, values are chemkin labels of the RMG species
                for network_label, adj in sa_dict['structures'].items():
                    labels_map[network_label] = get_species_label_by_structure(adj=adj, species_list=rmg_species)

                reactants_label = ' + '.join([key_by_val(labels_map, reactant.label) for reactant in reaction.reactants])
                products_label = ' + '.join([key_by_val(labels_map, product.label) for product in reaction.products])
                network_reaction_str = f'{reactants_label} <=> {products_label}'
                if network_reaction_str not in sa_dict:
                    log(f'Could not locate reaction {network_reaction_str} in SA output for network {network_name}.',
                        level='error', verbose=verbose)
                else:
                    # identify wells in this network this reaction is sensitive to
                    sensitive_wells_dict = dict()  # keys are wells labels, values are lists of sensitive conditions
                    for condition, sa_data in sa_dict[network_reaction_str].items():
                        max_sa_coeff = max([sa_coeff for sa_coeff in sa_data.values()])
                        for entry, sa_coeff in sa_data.items():
                            if '(TS)' not in entry and sa_coeff > max_sa_coeff * threshold:
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
                                spc = None
                                if spc_label is not None:
                                    spc = get_species_by_label(label=labels_map[label], species_list=rmg_species)
                                elif arkane is not None:
                                    # this is an Edge species which is missing from the Core rmg_species list
                                    spc = get_species_by_label(label=label, species_list=arkane.species_dict.values())
                                if spc is not None:
                                    species_list.append(spc)
                            for species in species_list:
                                if has_high_uncertainty(species, unconverged_species, species_to_calc):
                                    num = f'{reaction_tuple[1] + 1}{get_ordinal_indicator(reaction_tuple[1] + 1)} ' \
                                        if reaction_tuple[1] else ''
                                    reason = f'(iteration {iteration}) species participates in pressure dependent ' \
                                        f'network {network_name} from which reaction {chemkin_reaction_str} was ' \
                                        f'derived, which is the {num}most sensitive reaction for observable ' \
                                        f'{reaction_tuple[2]}, and the species is part of a sensitive well in this ' \
                                        f'network at the conditions {conditions}'
                                    species_to_calc[species.to_chemkin()] = {'spc': species, 'reason': reason}

    return species_to_calc, executed_networks


def modify_pdep_network_file(run_directory: str,
                             network_name: str,
                             method: str,
                             ) -> Tuple[str, str, tuple]:
    """
    A helper function for adding a P-dep SA directive to an Arkane network input file.

    Args:
        run_directory (str): The current RMG-ARC Tandem Tool iteration directory.
        network_name (str): The name of the original network file, e.g. 'network32_1'.
        method (str): 'CSE', 'MSC' or 'RS'.

    Raises:
        InputError: If ``method`` is not identified, or T/P ranges could not be parsed from the file.

    Returns:
        Tuple[str, str, tuple]:
            - The path to the new Arkane network input file.
            - The path to the Arkane output file.
            - Isomer labels of the current network.
    """
    if method not in ['CSE', 'MSC', 'RS']:
        raise InputError(f"ME method must be either 'CSE', 'RS', or 'MSC', got {method}")

    method_dict = {'CSE': 'chemically-significant eigenvalues',
                   'RS': 'reservoir state',
                   'MSC': 'modified strong collision',
                   }

    # copy the network file into a new folder (and rename it to input.py)
    sa_pdep_path = os.path.join(run_directory, 'pdep_sa', network_name, method)
    if not os.path.isdir(sa_pdep_path):
        os.makedirs(sa_pdep_path)
    input_file_path = os.path.join(sa_pdep_path, 'input.py')
    output_file_path = os.path.join(sa_pdep_path, 'sensitivity', 'sa_coefficients.yml')
    shutil.copyfile(src=os.path.join(run_directory, 'pdep', network_name + '.py'),
                    dst=input_file_path)

    with open(input_file_path, 'r') as f:
        lines = f.readlines()
        new_lines, isomer_labels = list(), list()
        t_min, t_max, p_min, p_max = None, None, None, None
        parse_tp, parse_isomers = False, (False, False)
        for line in lines:
            if 'pressureDependence(' in line:
                parse_tp = True
            if 'network(' in line:
                parse_isomers = (True, False)
            if parse_isomers[0] and 'isomers =' in line:
                parse_isomers = (True, True)
            if 'reactants =' in line:
                parse_isomers = (False, False)
            if parse_tp:
                if 'Tmin' in line:
                    #     Tmin = (300,'K'),
                    t_min = line.split('(')[1].split(',')[0]
                elif 'Tmax' in line:
                    #     Tmax = (2200,'K'),
                    t_max = line.split('(')[1].split(',')[0]
                elif 'Pmin' in line:
                    #     Pmin = (0.01,'bar'),
                    p_min = line.split('(')[1].split(',')[0]
                elif 'Pmax' in line:
                    #     Pmax = (100,'bar'),
                    p_max = line.split('(')[1].split(',')[0]
            if all(parse_isomers) and "'," in line:
                #         'C=O(26)',
                isomer_labels.append(line.split("'")[1])
            if 'method = ' in line:
                #     method = 'chemically-significant eigenvalues',
                splits = line.split("'")
                new_lines.append(f"{splits[0]}'{method_dict[method]}'{splits[2]}")
            elif 'rmgmode' in line:
                new_lines.append(line)
                if any(param is None for param in [t_min, t_max, p_min, p_max]):
                    raise InputError(f'Could not parse all T/P parameters, got:\n'
                                     f'T min = {t_min}, T max = {t_max}, P min = {p_min}, P max = {p_max}.')
                sa_conditions = f"""    sensitivity_conditions = [[({t_min}, 'K'), ({p_min}, 'bar')],
                              [({t_max}, 'K'), ({p_min}, 'bar')],
                              [({t_min}, 'K'), ({p_max}, 'bar')],
                              [({t_max}, 'K'), ({p_max}, 'bar')]],"""
                new_lines.append(sa_conditions)
            else:
                new_lines.append(line)

    with open(input_file_path, 'w') as f:
        f.writelines(new_lines)

    return input_file_path, output_file_path, tuple(isomer_labels)


def determine_species_based_on_collision_violators(run_directory: str,
                                                   rmg_species: list,
                                                   unconverged_species: list,
                                                   verbose: Optional[bool] = True,
                                                   ) -> dict:
    """
    Determine species to calculate based on collision rate violating reactions.

    Args:
        run_directory (str): The current RMG-ARC Tandem Tool iteration directory.
        rmg_species (list): All RMG Species parsed from the Chemkin file.
        unconverged_species (list): Entries are RMG Species for which thermo calculations did not converge
                                    in ARC throughout the iterations.
        verbose (Optional[bool]): Whether or not to log to file.

    Returns:
        dict: Keys are Chemkin labels, values are dicts of {'spc': RMG Species objects, 'reason': ``str``}
    """
    species_to_calc = dict()  # keys are labels, values are dicts of {'spc': RMG Species objects, 'reason': ``str``}
    coll_violators_path = os.path.join(run_directory, 'collision_rate_violators.log')
    if not os.path.isfile(coll_violators_path):
        log('No collision rate violating reactions identified in this model.', verbose=verbose)
        return dict()

    with open(coll_violators_path, 'r') as f:
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
                species = get_species_by_label(label, rmg_species)
                if species is None:
                    log(f'Could not identify species {label}!',
                        level='error', verbose=verbose)
                if has_high_uncertainty(species, unconverged_species, species_to_calc):
                    reason = f'species participates in a collision rate violating reaction, {rxn_to_log}'
                    species_to_calc[species.to_chemkin()] = {'spc': species, 'reason': reason}
    return species_to_calc


def add_rmg_libraries(run_directory: str,
                      library_name: Optional[str] = None,
                      verbose: Optional[bool] = True,
                      ) -> str:
    """
    Creates a "t3_thermo" library in the correct place in the RMG database repository if it doesn't already exist,
    and appends the respective entries from the library generated by ARC.
    The libraries generated by ARC are located in ``run_directory``/ARC/output/RMG libraries/

    Args:
        run_directory (str): The directory of the recently terminated RMG-ARC iteration.
        library_name (Optional[str]): The library name. None in the first iteration (where it is created),
                                      other iterations should pass this argument for consistency.
        verbose (Optional[bool]): Whether or not to log to file.

    Returns:
         str: The (created) library_name.
    """
    global rmg_thermo_lib_base_path
    arc_thermo_lib_path = os.path.join(run_directory, 'ARC', 'output', 'RMG libraries', 'thermo', 't3.py')
    rmg_thermo_lib_path = ''
    if library_name is None:
        # This is the first iteration, come up with a unique library_name
        unique_library_name = False
        j = 0
        library_name = 't3_thermo'
        while not unique_library_name:
            # make sure the new library name is unique and used throughout the project
            rmg_thermo_lib_path = os.path.join(rmg_thermo_lib_base_path, f'{library_name}.py')
            if not os.path.isfile(rmg_thermo_lib_path):
                unique_library_name = True
            else:
                library_name = 't3_thermo' + '_' + str(j)
                j += 1
        log(f'Created the RMG Thermo library {library_name}', verbose=verbose)
    else:
        # This is not the first iteration, use the provided library_name
        rmg_thermo_lib_path = os.path.join(settings['database.directory'], 'thermo', 'libraries', f'{library_name}.py')
    local_context = {'ThermoData': ThermoData, 'Wilhoit': Wilhoit, 'NASAPolynomial': NASAPolynomial, 'NASA': NASA}
    if os.path.isfile(rmg_thermo_lib_path) and os.path.isfile(arc_thermo_lib_path):
        # the t3 thermo library already exists in the RMG database: load it, append new entries, and save.
        rmg_thermo_lib, arc_thermo_lib = ThermoLibrary(), ThermoLibrary()
        rmg_thermo_lib.load(path=rmg_thermo_lib_path, local_context=local_context, global_context=dict())
        arc_thermo_lib.load(path=arc_thermo_lib_path, local_context=local_context, global_context=dict())
        description = arc_thermo_lib.long_desc
        description_to_append = '\n'
        append = False
        for line in description.splitlines():
            if 'Considered the following' in line:
                append = True
            if append:
                description_to_append += line + '\n'
        rmg_thermo_lib.long_desc += description_to_append
        for entry in arc_thermo_lib.entries.values():
            unique_species_name, j = False, 0
            label = entry.label
            while not unique_species_name:
                # make sure each entry has a unique label in the unified library
                for existing_entry in rmg_thermo_lib.entries.values():
                    if label == existing_entry.label:
                        label = entry.label + '_' + str(j)
                        j += 1
                        break
                else:
                    unique_species_name = True
            entry.label = label
            rmg_thermo_lib.entries[entry.label] = entry
        rmg_thermo_lib.save(path=rmg_thermo_lib_path)
    elif not os.path.isfile(rmg_thermo_lib_path) and os.path.isfile(arc_thermo_lib_path):
        # the t3 thermo library doesn't exist in the RMG database: just copy the library generated by ARC.
        shutil.copy(arc_thermo_lib_path, rmg_thermo_lib_path)
    return library_name


def get_unconverged_species(run_directory: str,
                            all_species: list,
                            log_species: Optional[bool] = True,
                            verbose: Optional[bool] = True,
                            ) -> list:
    """
    Get the labels of unconverged species from the present ARC iteration.

    Args:
        run_directory (str): The directory name of the present iteration.
        all_species (list): Entries are RMG Species objects sent to calculation so far.
        log_species (Optional[bool]): Whether to log unconverged species, ``True`` to log, used for testing.
        verbose (Optional[bool]): Whether or not to log to file.

    Returns:
        list: Entries are RMG Species objects which could not be calculated in the present iteration.
    """
    unconverged_species_labels, unconverged_species = list(), list()
    info_path = os.path.join(run_directory, 'ARC', 't3.info')
    if os.path.isfile(info_path):
        with open(info_path, 'r') as f:
            read = False
            for line in f:
                if read:
                    if line == '\n':
                        break
                    if 'Species' in line and '(Failed!)' in line:
                        unconverged_species_labels.append(line.split()[1])
                if 'Considered the following species' in line:
                    read = True
    for label in unconverged_species_labels:
        for spc in all_species:
            if label == spc.label:
                unconverged_species.append(spc)
    if log_species:
        log_unconverged_species(unconverged_species, verbose=verbose)
    return unconverged_species
