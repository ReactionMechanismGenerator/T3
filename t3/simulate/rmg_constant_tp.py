"""
T3's RMG Simulator Adapter module
Used to run mechanism analysis with RMG via subprocess execution.
"""

import datetime
import itertools
import logging
import os
import numpy as np
from typing import List, Optional, TYPE_CHECKING

from arc.common import read_yaml_file

from t3.common import get_chem_to_rmg_rxn_index_map, get_values_within_range, time_lapse
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter
from t3.runners.rmg_runner import run_rmg_sa_incore
from t3.utils.writer import write_rmg_input_file

if TYPE_CHECKING:
    from t3.logger import Logger


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
        sa_rtol (float, optional): The relative tolerance used when performing sensitivity analysis.
        global_observables (Optional[List[str]]): List of global observables ['IDT', 'ESR', 'SL'] used by Cantera adapters.
    """

    def __init__(self,
                 t3: dict,
                 rmg: dict,
                 paths: dict,
                 logger: 'Logger',
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
        self.atol = atol
        self.rtol = rtol
        self.sa_atol = sa_atol
        self.sa_rtol = sa_rtol
        self.global_observables = global_observables

        self.rmg_input_file = None
        self.sa_output_file = None

        self.set_up()

    def set_up(self):
        """
        Prepare directories and input files.
        """
        if len(self.observable_list):
            self.rmg_input_file = self.paths['SA input']
            self.sa_output_file = os.path.join(self.paths['SA'], 'sa_output.yml')
            if not os.path.isdir(self.paths['SA']):
                os.mkdir(self.paths['SA'])
            # Clean up old files
            if os.path.isfile(self.paths['SA input']):
                os.remove(self.paths['SA input'])
            if os.path.isfile(self.sa_output_file):
                os.remove(self.sa_output_file)
        else:
            self.rmg_input_file = self.paths['RMG input']
            self.sa_output_file = None

        # 1. Generate Reactors
        rmg_config = self.generate_rmg_reactors_for_simulation()

        num_conditions = len(rmg_config['reactors'])
        if self.observable_list:
            self.logger.info(f'Running a simulation with SA using RMGConstantTP for {num_conditions} conditions...')
        else:
            self.logger.info('Running a simulation using RMGConstantTP...')

        # 2. Inject Sensitivity Settings into Reactors
        # This ensures the writer generates: simpleReactor(..., sensitivity=['A'], sensitivityThreshold=0.01)
        if self.observable_list:
            for reactor in rmg_config['reactors']:
                reactor['sensitivity'] = self.observable_list
                reactor['sensitivityThreshold'] = self.t3['sensitivity']['SA_threshold']

        # 3. Write the RMG input file
        write_rmg_input_file(
            rmg=rmg_config,
            t3=self.t3,
            iteration=1,
            path=self.rmg_input_file,
            walltime=self.t3['options']['max_RMG_walltime'],
        )

    def simulate(self):
        """
        Simulates the model using the RMG constant T and P simulator adapter via a subprocess.
        """
        tic = datetime.datetime.now()

        if self.observable_list:
            self.logger.info(f'Running RMG SA subprocess for observables: {self.observable_list}...')

            success, error_msg = run_rmg_sa_incore(
                rmg_input_file_path=self.rmg_input_file,
                chemkin_file_path=self.paths['chem annotated'],
                species_dict_path=self.paths['species dict'],
                output_path=self.sa_output_file,
                observables=self.observable_list,
                threshold=self.t3['sensitivity']['SA_threshold']
            )

            if not success:
                self.logger.error(f"RMG SA subprocess failed.\nDetails: {error_msg}")
            else:
                self.logger.info("RMG SA subprocess completed successfully.")
        else:
            self.logger.info('No observables listed, skipping SA subprocess.')

        self.logger.info(f'Simulation via RMG completed, execution time: {time_lapse(tic)}')

    def get_sa_coefficients(self) -> Optional[dict]:
        """
        Obtain the SA coefficients by parsing the YAML generated by the script.
        """
        if not self.sa_output_file or not os.path.isfile(self.sa_output_file):
            self.logger.error("Could not find the SA output YAML file.")
            return None

        raw_sa_dict = read_yaml_file(self.sa_output_file)

        if not raw_sa_dict:
            return None

        chem_to_rmg_map = get_chem_to_rmg_rxn_index_map(chem_annotated_path=self.paths['chem annotated'])
        rmg_to_chem_map = {v: k for k, v in chem_to_rmg_map.items()}

        kinetics = dict()
        thermo = dict()
        time_data = np.array(raw_sa_dict.get('time', []))

        if 'kinetics' in raw_sa_dict:
            for obs_label, param_dict in raw_sa_dict['kinetics'].items():
                if obs_label not in kinetics:
                    kinetics[obs_label] = dict()

                for rmg_index, values in param_dict.items():
                    arr_values = np.array(values)
                    try:
                        rmg_idx_int = int(rmg_index)
                        chem_index = rmg_to_chem_map.get(rmg_idx_int)
                        if chem_index is not None:
                            kinetics[obs_label][chem_index] = arr_values
                        else:
                            kinetics[obs_label][rmg_index] = arr_values
                    except ValueError:
                        kinetics[obs_label][rmg_index] = arr_values

        if 'thermo' in raw_sa_dict:
            for obs_label, param_dict in raw_sa_dict['thermo'].items():
                if obs_label not in thermo:
                    thermo[obs_label] = dict()
                for spc_label, values in param_dict.items():
                    thermo[obs_label][spc_label] = np.array(values)

        return {'time': [time_data], 'kinetics': [kinetics], 'thermo': [thermo]}

    def get_idt_by_T(self) -> dict:
        return {'idt': list(), 'idt_index': list()}

    def generate_rmg_reactors_for_simulation(self) -> dict:
        """
        Turn all RMG ranged reactors into individual reactors.
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
        """
        species_lists = list()
        spc_indices_w_ranges = [i for i, spc in enumerate(self.rmg['species'])
                                if isinstance(spc['concentration'], (list, tuple))]
        species_list = [{'label': spc['label'], 'concentration': spc['concentration']} for spc in self.rmg['species']
                        if (isinstance(spc['concentration'], (float, int)) and spc['concentration'] > 0)
                        or spc['balance'] or not spc['reactive']]
        species_vals = [get_values_within_range(
            value_range=self.rmg['species'][spc_indices_w_ranges.index(species_index)]['concentration'],
            num=self.t3['options']['num_sa_per_concentration_range'])
                        for species_index in spc_indices_w_ranges]

        # 1. No ranged concentrations
        if len(spc_indices_w_ranges) == 0:
            return [[{'label': spc['label'], 'concentration': spc['concentration']} for spc in self.rmg['species']
                     if spc['concentration'] > 0 or spc['balance'] or not spc['reactive']]]

        # 2. Only two ranged concentrations and modify_concentration_ranges_in_reverse is True
        if len(spc_indices_w_ranges) == 2 and self.t3['options']['modify_concentration_ranges_in_reverse']:
            spc_0_vals = get_values_within_range(
                value_range=self.rmg['species'][spc_indices_w_ranges[0]]['concentration'],
                num=self.t3['options']['num_sa_per_concentration_range'])
            spc_1_vals = get_values_within_range(
                value_range=self.rmg['species'][spc_indices_w_ranges[1]]['concentration'],
                num=self.t3['options']['num_sa_per_concentration_range'])
            for val_0, val_1 in zip(spc_0_vals, spc_1_vals[::-1]):
                new_species_list = list(species_list)
                new_species_list.append({'label': self.rmg['species'][spc_indices_w_ranges[0]]['label'],
                                         'concentration': val_0})
                new_species_list.append({'label': self.rmg['species'][spc_indices_w_ranges[1]]['label'],
                                         'concentration': val_1})
                species_lists.append(new_species_list)

        # 3. No combinations, modify_concentration_ranges_together is True
        elif self.t3['options']['modify_concentration_ranges_together']:
            for point_number in range(self.t3['options']['num_sa_per_concentration_range']):
                new_species_list = list(species_list)
                for i, spc_index in enumerate(spc_indices_w_ranges):
                    new_species_list.append({'label': self.rmg['species'][spc_index]['label'],
                                             'concentration': species_vals[i][point_number]})
                species_lists.append(new_species_list)

        # 4. Combinations (products)
        else:
            for vals in itertools.product(*species_vals):
                new_species_list = list(species_list)
                for i, val in enumerate(vals):
                    new_species_list.append({'label': self.rmg['species'][spc_indices_w_ranges[i]]['label'],
                                             'concentration': val})
                species_lists.append(new_species_list)

        for species_list in species_lists:
            total = sum(s['concentration'] for s in species_list
                        if isinstance(s['concentration'], (int, float)))
            if total > 1.0 or total <= 0:
                logging.getLogger(__name__).warning(
                    f'Species concentrations sum to {total:.4f}, expected (0, 1]. '
                    f'Cantera will normalize mole fractions.')
            species_list.sort(key=lambda spc: spc['concentration'][0] if isinstance(spc['concentration'], (tuple, list))
            else spc['concentration'], reverse=True)
        return species_lists


register_simulate_adapter("RMGConstantTP", RMGConstantTP)
