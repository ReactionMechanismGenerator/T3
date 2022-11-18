"""
RMG Simulator Adapter module
Used to run mechanism analysis with RMG
"""

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

from t3.common import get_species_by_label
from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter


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
                 global_observables: Optional[List[str]] = None
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

        if not os.path.isfile(self.rmg_input_file):
            shutil.copyfile(src=self.paths['RMG input'], dst=self.rmg_input_file)

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
            self.logger.info('Running a simulation with SA using RMGConstantTP...')
            self.observable_species = [species for species in self.rmg_model.reaction_model.core.species
                                       if species.label in self.observable_list]
            for i, reaction_system in enumerate(self.rmg_model.reaction_systems):
                reaction_system.sensitive_species = self.observable_species
                reaction_system.sensitivity_threshold = self.t3['sensitivity']['SA_threshold']
                if hasattr(reaction_system, 'Trange') and reaction_system.Trange is not None:
                    temperature = sum([t.value_si for t in reaction_system.Trange]) / len(reaction_system.Trange)
                    self.logger.info(f'An average temperature of {temperature:.1f} K is taken for the SA'
                                     f' for RMG reactor {i}')
                else:
                    temperature = reaction_system.T.value_si
                reaction_system.sens_conditions['T'] = temperature
                if isinstance(reaction_system, SimpleReactor):
                    if hasattr(reaction_system, 'Prange') and reaction_system.Prange is not None:
                        pressure = sum([p.value_si for p in reaction_system.Prange]) / len(reaction_system.Prange)
                        self.logger.info(f'An average pressure of {pressure * 1e-5:.2f} bar is taken for the SA'
                                         f' for RMG reactor {i}')
                    else:
                        pressure = reaction_system.P.value_si
                    reaction_system.sens_conditions['P'] = pressure
                elif isinstance(reaction_system, LiquidReactor):
                    if hasattr(reaction_system, 'Vrange') and reaction_system.Vrange is not None:
                        volume = sum([v for v in reaction_system.Vrange]) / len(reaction_system.Vrange)
                        self.logger.info(f'An average volume of {volume:.2f} m^3 is taken for the SA'
                                         f' for RMG reactor {i}')
                    else:
                        volume = reaction_system.V
                    reaction_system.sens_conditions['V'] = volume
                else:
                    raise NotImplementedError(f'RMG SA not implemented for Reactor type {type(reaction_system)}.')
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

            if reaction_system.sensitive_species:
                plot_sensitivity(self.rmg_model.output_directory, index, reaction_system.sensitive_species)

            if self.rmg_model.uncertainty is not None:
                self.rmg_model.run_uncertainty_analysis()

    def get_sa_coefficients(self) -> Optional[dict]:
        """
        Obtain the SA coefficients.

        Returns:
             sa_dict (Optional[dict]): An SA dictionary, structure is given in the docstring for T3/t3/main.py
        """
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
                        parameter = int(parameter) if all(c.isdigit() for c in parameter) else parameter
                    sa_dict[sa_type][observable_label][parameter] = df[header].values
        return sa_dict

    def get_idt_by_T(self) -> dict:
        """
        Finds the ignition point by approximating dT/dt as a first order forward difference and then finds
        the point of maximum slope. However, the RMG reactors only simulate at constant T, so this method
        returns a dictionary whose values are empty lists

        Returns:
            idt_dict (dict): Dictionary whose keys include 'idt' and 'idt_index' and whose values are lists of
                             the ignition delay time in seconds and index at which idt occurs respectively.
        """
        idt_dict = {'idt': list(),
                    'idt_index': list(),
                    }
        return idt_dict


register_simulate_adapter("RMGConstantTP", RMGConstantTP)
