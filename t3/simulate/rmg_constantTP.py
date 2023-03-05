"""
RMG Simulator Adapter module
Used to run mechanism analysis with RMG
"""

import datetime
import itertools
import os
import pandas as pd
import shutil
from typing import List, Optional

from rmgpy.kinetics.diffusionLimited import diffusion_limiter
from rmgpy.rmg.listener import SimulationProfilePlotter, SimulationProfileWriter
from rmgpy.rmg.settings import ModelSettings
from rmgpy.solver.liquid import LiquidReactor
from rmgpy.solver.simple import SimpleReactor
from rmgpy.tools.loader import load_rmg_py_job
from rmgpy.tools.plot import plot_sensitivity

from t3.common import get_chem_to_rmg_rxn_index_map, get_species_by_label, get_values_within_range, time_lapse
from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter
from t3.utils.writer import write_rmg_input_file


class RMGConstantTP(SimulateAdapter):
    """
    RMGConstantTP is an adapter for the abstract class SimulateAdapter that simulates mechanisms in a homogenous,
    isothermal, isobaric batch reactor under ideal gas or ideal liquid conditions.

    Args:
        t3 (dict): The T3.t3 attribute, which is a dictionary containing the t3 block from the input yaml or API.
        rmg (dict): The T3.rmg attribute, which is a dictionary containing the rmg block from the input yaml or API.
        paths (dict): The T3.paths attribute, which is a dictionary containing relevant paths.
        logger (Logger): Instance of T3's Logger class.
        atol (float): The absolute tolerance used when integrating.
        rtol (float): The relative tolerance used when integrating.
        observable_list (Optional[list]): Species used for SA. Entries are species labels as strings. Example: ['OH']
        sa_atol (float, optional): The absolute tolerance used when performing sensitivity analysis.
        sa_atol (float, optional): The relative tolerance used when performing sensitivity analysis.
        global_observables (Optional[List[str]]): List of global observables ['IgD', 'ESR', 'SL'] used by Cantera adapters.

    Attributes:
        atol (float): The absolute tolerance used when integrating during an RMG iteration.
        global_observables (List[str]): List of global observables ['IgD', 'ESR', 'SL'] used by Cantera adapters.
        logger (Logger): Instance of T3's Logger class.
        observable_list (list): Species used for SA. Entries are species labels as strings. Example: ['OH']
        observable_species (list): Species object representations of the species used for SA.
        paths (dict): The T3.paths attribute, which is a dictionary containing relevant paths.
        rmg (dict): The T3.rmg attribute, which is a dictionary containing the rmg block from the input yaml or API.
        rmg_model (RMG class): Representation of an RMG job.
        rtol (float): The relative tolerance used when integrating during an RMG iteration.
        sa_atol (float): The absolute tolerance used when performing sensitivity analysis.
        sa_rtol (float): The relative tolerance used when performing sensitivity analysis.
        t3 (dict): The T3.t3 attribute, which is a dictionary containing the t3 block from the input yaml or API.
    """

    def __init__(self,
                 t3: dict,
                 rmg: dict,
                 paths: dict,
                 logger: Logger,
                 atol: float = 1e-16,
                 rtol: float = 1e-8,
                 observable_list: Optional[list] = None,
                 sa_atol: float = 1e-6,
                 sa_rtol: float = 1e-4,
                 global_observables: Optional[List[str]] = None,
                 ):

        self.t3 = t3
        self.rmg = rmg
        self.paths = paths
        self.logger = logger
        self.observable_list = observable_list or list()

        # for this adapter, tolerances are read in from the rmg input file created by writer.py
        self.atol = atol
        self.rtol = rtol
        self.sa_atol = sa_atol
        self.sa_rtol = sa_rtol

        # cantera argument not used by this adapter
        self.global_observables = global_observables

        # initialize other attributes
        self.observable_species = list()
        self.rmg_model = None
        self.rmg_input_file = None

        self.set_up()

    def set_up(self):
        """
        Read in the chemkin file, species dictionary, and RMG input file.
        If the user requested SA, simulate the job with SA and generate the output csv files.
        If SA is not requested, only simulate the job to obtain species profiles.

        Raises:
            FileNotFoundError: If the RMG adapter does not properly read the rmg input file.
            ValueError: If RMG SA is not implemented for the given reactor type.
        """
        if len(self.observable_list):
            self.rmg_input_file = self.paths['SA input']
            if not os.path.isdir(self.paths['SA']):
                os.mkdir(self.paths['SA'])
            if not os.path.isdir(self.paths['SA solver']):
                os.mkdir(self.paths['SA solver'])
            if os.path.isfile(self.paths['SA input']):
                os.remove(self.paths['SA input'])
        else:
            self.rmg_input_file = self.paths['RMG input']

        write_rmg_input_file(
            rmg=self.generate_rmg_reactors_for_simulation(),
            t3=self.t3,
            iteration=1,  # Does not matter for simulating or computing SA.
            path=self.rmg_input_file,
            walltime=self.t3['options']['max_RMG_walltime'],
        )

        with open(self.rmg_input_file, 'r') as f:
            lines = f.readlines()
        restart_string = "restartFromSeed(path='seed')"
        new_lines = [line for line in lines if restart_string not in line]
        with open(self.rmg_input_file, 'w') as f:
            f.writelines(new_lines)

        self.rmg_model = load_rmg_py_job(input_file=self.rmg_input_file,
                                         chemkin_file=self.paths['chem annotated'],
                                         species_dict=self.paths['species dict'],
                                         generate_images=True,
                                         use_chemkin_names=False,
                                         check_duplicates=False,
                                         )
        if self.rmg_model is None:
            raise FileNotFoundError('The RMG adapter did not properly read the rmg input file.')

        if len(self.observable_list):
            self.logger.info(f'Running a simulation with SA using RMGConstantTP for '
                             f'{len(self.rmg_model.reaction_systems)} conditions...')
            self.observable_species = [species for species in self.rmg_model.reaction_model.core.species
                                       if species.label in self.observable_list]
            for i, reaction_system in enumerate(self.rmg_model.reaction_systems):
                reaction_system.sensitive_species = self.observable_species
                reaction_system.sensitivity_threshold = self.t3['sensitivity']['SA_threshold']
                reaction_system.sens_conditions['T'] = reaction_system.T.value_si
                if isinstance(reaction_system, SimpleReactor):
                    reaction_system.sens_conditions['P'] = reaction_system.P.value_si
                elif isinstance(reaction_system, LiquidReactor):
                    reaction_system.sens_conditions['V'] = reaction_system.V
                else:
                    raise NotImplementedError(f'RMG SA is not implemented for Reactor type {type(reaction_system)}.')
        else:
            self.logger.info('Running a simulation using RMGConstantTP...')

    def simulate(self):
        """
        Simulates the model using the RMG constant T and P simulator adapter.
        """
        solver_path = os.path.join(self.rmg_model.output_directory, 'solver')
        if os.path.exists(solver_path):
            shutil.rmtree(solver_path, ignore_errors=True)
        if not os.path.exists(solver_path):
            os.mkdir(solver_path)

        tic = datetime.datetime.now()
        for index, reaction_system in enumerate(self.rmg_model.reaction_systems):
            if reaction_system.sensitive_species and reaction_system.sensitive_species == ['all']:
                reaction_system.sensitive_species = self.rmg_model.reaction_model.core.species
            reaction_system.attach(SimulationProfileWriter(output_directory=self.rmg_model.output_directory,
                                                           reaction_sys_index=index,
                                                           core_species=self.rmg_model.reaction_model.core.species))
            reaction_system.attach(SimulationProfilePlotter(output_directory=self.rmg_model.output_directory,
                                                            reaction_sys_index=index,
                                                            core_species=self.rmg_model.reaction_model.core.species))

            sens_worksheet, pdep_networks = list(), list()
            for spc in reaction_system.sensitive_species:
                csv_path = os.path.join(self.rmg_model.output_directory, 'solver', f'sensitivity_{index + 1}_SPC_{spc.index}.csv')
                sens_worksheet.append(csv_path)
            for source, networks in self.rmg_model.reaction_model.network_dict.items():
                pdep_networks.extend(networks)

            model_settings = ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1, tol_interrupt_simulation=1)
            simulator_settings = self.rmg_model.simulator_settings_list[-1]

            if isinstance(reaction_system, LiquidReactor):  # Relocate to a liquid reactor RMG adapter once created.
                self.rmg_model.load_database()
                solvent_data = self.rmg_model.database.solvation.get_solvent_data(self.rmg_model.solvent)
                diffusion_limiter.enable(solvent_data, self.rmg_model.database.solvation)
            elif self.rmg_model.uncertainty is not None:
                self.rmg_model.verbose_comments = True
                self.rmg_model.load_database()

            if reaction_system.const_spc_names is not None:
                reaction_system.get_const_spc_indices(self.rmg_model.reaction_model.core.species)

            try:
                reaction_system.simulate(
                    core_species=self.rmg_model.reaction_model.core.species,
                    core_reactions=self.rmg_model.reaction_model.core.reactions,
                    edge_species=self.rmg_model.reaction_model.edge.species,
                    edge_reactions=self.rmg_model.reaction_model.edge.reactions,
                    surface_species=[],
                    surface_reactions=[],
                    pdep_networks=pdep_networks,
                    sensitivity=True if reaction_system.sensitive_species else False,
                    sens_worksheet=sens_worksheet,
                    model_settings=model_settings,
                    simulator_settings=simulator_settings,
                    conditions={'T': reaction_system.sens_conditions['T'], 'P': reaction_system.sens_conditions['P']},
                    prune=False,
                )
            except ZeroDivisionError as e:
                self.logger.warning(f'Cannot simulate reaction system, got:\n{e}')

            if reaction_system.sensitive_species:
                try:
                    plot_sensitivity(self.rmg_model.output_directory, index, reaction_system.sensitive_species)
                except FileNotFoundError as e:
                    self.logger.warning(f'Cannot plot sensitivity, got:\n{e}')

            if self.rmg_model.uncertainty is not None:
                self.rmg_model.run_uncertainty_analysis()

        self.logger.info(f'Simulation via RMG completed, execution time: {time_lapse(tic)}')

    def get_sa_coefficients(self) -> Optional[dict]:
        """
        Obtain the SA coefficients.

        Returns:
             sa_dict (Optional[dict]): An SA dictionary, structure is given in the docstring for T3/t3/main.py
        """
        chem_to_rmg_rxn_index_map = get_chem_to_rmg_rxn_index_map(chem_annotated_path=self.paths['chem annotated'])
        solver_path = os.path.join(self.paths['SA'], 'solver')
        if not os.path.exists(solver_path):
            self.logger.error("Could not find the path to RMG's SA solver output folder.")
            return None
        sa_files = list()
        for file_ in os.listdir(solver_path):
            if 'sensitivity' in file_ and file_.endswith(".csv"):
                sa_files.append(file_)
        sa_dict = {'kinetics': dict(), 'thermo': dict(), 'time': list()}
        for sa_file in sa_files:
            df = pd.read_csv(os.path.join(solver_path, sa_file))
            for header in df.columns:
                sa_type = None
                if 'Time' in header:
                    sa_dict['time'] = df[header].values
                elif '/dln[k' in header:
                    sa_type = 'kinetics'
                elif '/dG[' in header:
                    sa_type = 'thermo'
                if sa_type is not None:
                    observable_label = header.split('[')[1].split(']')[0]
                    observable = get_species_by_label(observable_label, self.rmg_model.reaction_model.core.species)
                    if observable is None:
                        self.logger.error(f'Could not identify observable species for label: {observable_label}')
                    observable_label = observable.to_chemkin()
                    if observable_label not in sa_dict[sa_type].keys():
                        sa_dict[sa_type][observable_label] = dict()
                    # parameter extraction examples:
                    # for species get 'C2H4(8)' from `dln[ethane(1)]/dG[C2H4(8)]`
                    # for reaction, get 8 from `dln[ethane(1)]/dln[k8]: H(6)+ethane(1)=H2(12)+C2H5(5)`
                    parameter = header.split('[')[2].split(']')[0]
                    if sa_type == 'kinetics':
                        parameter = parameter[1:]
                        parameter = chem_to_rmg_rxn_index_map[int(parameter)] \
                            if all(c.isdigit() for c in parameter) else parameter
                    sa_dict[sa_type][observable_label][parameter] = df[header].values
        return sa_dict

    def get_idt_by_T(self) -> dict:
        """
        Finds the ignition point by approximating dT/dt as a first order forward difference and then finds
        the point of maximum slope. However, the RMG reactors only simulate at constant T, so this implementation
        only returns a dictionary whose values are empty lists.

        Returns:
            idt_dict (dict): Dictionary whose keys include 'idt' and 'idt_index' and whose values are lists of
                             the ignition delay time in seconds and index at which idt occurs respectively.
        """
        idt_dict = {'idt': list(),
                    'idt_index': list(),
                    }
        return idt_dict

    def generate_rmg_reactors_for_simulation(self) -> dict:
        """
        Turn all RMG ranged reactors into individual reactors with specific species concentrations,
        temperature, pressure (for gas phase simulations), and volume (for liquid phase simulations).

        Returns:
            dict: A modified dictionary representing the self.rmg dictionary with individual reactor representations.
        """
        ranged_species = any(isinstance(spc['concentration'], list) for spc in self.rmg['species'])
        ranged_t = any(isinstance(reactor['T'], list) for reactor in self.rmg['reactors'])
        ranged_p = any('P' in reactor.keys() and isinstance(reactor['P'], list) for reactor in self.rmg['reactors'])
        ranged_v = any('V' in reactor.keys() and isinstance(reactor['V'], list) for reactor in self.rmg['reactors'])
        if not any([ranged_species, ranged_t, ranged_p, ranged_v]):
            return self.rmg

        mod_rmg = self.rmg.copy()
        mod_rmg['reactors'] = list()

        species_lists = self.get_species_concentration_lists_from_ranged_params()
        for reactor in self.rmg['reactors']:
            t_vals = get_values_within_range(value_range=reactor['T'],
                                             num=self.t3['options']['num_sa_per_temperature_range'])
            p_vals = v_vals = None
            if 'P' in reactor.keys():
                p_vals = get_values_within_range(value_range=reactor['P'],
                                                 num=self.t3['options']['num_sa_per_pressure_range'],
                                                 use_log_scale=True)
            elif 'V' in reactor.keys():
                v_vals = get_values_within_range(value_range=reactor['V'],
                                                 num=self.t3['options']['num_sa_per_volume_range'])
            for t_val in t_vals:
                for param in p_vals if ranged_p else v_vals:
                    for species_list in species_lists:
                        new_reactor = {k: v for k, v in reactor.items() if k not in ['T', 'P', 'V']}
                        new_reactor['T'] = t_val
                        new_reactor['P' if ranged_p else 'V'] = param
                        new_reactor['species_list'] = species_list
                        mod_rmg['reactors'].append(new_reactor)

        return mod_rmg

    def get_species_concentration_lists_from_ranged_params(self) -> List[List[dict]]:
        """
        Get a list of species concentrations with specific values from a range of concentrations.

        Returns:
            List[List[dict]]: Entries are lists of dictionaries describing an RMG species concentration by label.
        """
        species_lists = list()
        spc_indices_w_ranges = [i for i, spc in enumerate(self.rmg['species'])
                                if isinstance(spc['concentration'], (list, tuple))]
        species_list = [{'label': spc['label'], 'concentration': spc['concentration']} for spc in self.rmg['species']
                        if (isinstance(spc['concentration'], (float, int)) and spc['concentration'] > 0)
                        or spc['balance'] or not spc['reactive']]
        species_vals = [get_values_within_range(value_range=self.rmg['species'][spc_indices_w_ranges[species_index]]['concentration'],
                                                num=self.t3['options']['num_sa_per_concentration_range'])
                        for species_index in spc_indices_w_ranges]

        # 1. No ranged concentrations
        if len(spc_indices_w_ranges) == 0:
            return [[{'label': spc['label'], 'concentration': spc['concentration']} for spc in self.rmg['species']
                     if spc['concentration'] > 0 or spc['balance'] or not spc['reactive']]]

        # 2. Only two ranged concentrations and modify_concentration_ranges_in_reverse is True
        if len(spc_indices_w_ranges) == 2 and self.t3['options']['modify_concentration_ranges_in_reverse']:
            spc_0_vals = get_values_within_range(value_range=self.rmg['species'][spc_indices_w_ranges[0]]['concentration'],
                                                 num=self.t3['options']['num_sa_per_concentration_range'])
            spc_1_vals = get_values_within_range(value_range=self.rmg['species'][spc_indices_w_ranges[1]]['concentration'],
                                                 num=self.t3['options']['num_sa_per_concentration_range'])
            for val_0, val_1 in zip(spc_0_vals, spc_1_vals[::-1]):
                new_species_list = species_list
                new_species_list.append({'label': self.rmg['species'][spc_indices_w_ranges[0]]['label'],
                                         'concentration': val_0})
                new_species_list.append({'label': self.rmg['species'][spc_indices_w_ranges[1]]['label'],
                                         'concentration': val_1})
                species_lists.append(new_species_list)

        # 3. No combinations, modify_concentration_ranges_together is True
        elif self.t3['options']['modify_concentration_ranges_together']:
            for point_number in range(self.t3['options']['num_sa_per_concentration_range']):
                new_species_list = species_list
                for i, spc_index in enumerate(spc_indices_w_ranges):
                    new_species_list.append({'label': self.rmg['species'][spc_index]['label'],
                                             'concentration': species_vals[i][point_number]})
                species_lists.append(new_species_list)

        # 4. Combinations (products)
        else:
            for vals in itertools.product(*species_vals):
                new_species_list = species_list
                for i, val in enumerate(vals):
                    new_species_list.append({'label': self.rmg['species'][spc_indices_w_ranges[i]]['label'],
                                             'concentration': val})
                species_lists.append(new_species_list)

        for species_list in species_lists:
            species_list.sort(key=lambda spc: spc['concentration'][0] if isinstance(spc['concentration'], (tuple, list))
                              else spc['concentration'], reverse=True)
        return species_lists


register_simulate_adapter("RMGConstantTP", RMGConstantTP)
