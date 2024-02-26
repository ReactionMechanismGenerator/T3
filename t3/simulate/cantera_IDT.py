"""
Cantera Simulator Adapter module for ignition delay time (IDT) calculations.
Used to run mechanism analysis with Cantera for IDT.
"""

from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import os

import cantera as ct
import math
import numpy as np

from arc.common import save_yaml_file

from rmgpy.chemkin import load_chemkin_file

from t3.common import determine_concentrations_by_equivalence_ratios
from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter


class CanteraIDT(SimulateAdapter):
    """
    CanteraIDT is an adapter for the abstract class SimulateAdapter that simulates ignition delay time.

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


    Attributes:
        atol (float): The absolute tolerance used when integrating.
        cantera_simulation (ct.ReactorNet): Cantera reactor net object.
        inert_list (list): List of possible inert species in the model
        inert_index_list (list): List of indices corresponding to the inert species present in the model.
        logger (Logger): Instance of T3's Logger class.
        model (ct.Solution): Cantera solution object for the mechanism.
        num_ct_reactions (int): Number of reactions in the model.
        num_ct_species (int): Number of species in the model.
        paths (dict): The T3.paths attribute, which is a dictionary containing relevant paths.
        rmg (dict): The T3.rmg attribute, which is a dictionary containing the rmg block from the input yaml or API.
        rtol (float): The relative tolerance used when integrating.
        rxn_identifier_lookup (dict): Keys are reactions (str). Values are index in the model.
        sa_atol (float): The absolute tolerance used when performing sensitivity analysis.
        sa_rtol (float): The relative tolerance used when performing sensitivity analysis.
        radical_label (str): The label of the prominent ignition radical: OH if present, else H (e.g., in hydrazine).
        spc_identifier_lookup (dict): Keys are species (str). Values are index in the model.
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
                 ):

        self.t3 = t3
        self.rmg = rmg
        self.paths = paths
        self.logger = logger
        self.atol = atol
        self.rtol = rtol
        self.sa_atol = sa_atol
        self.sa_rtol = sa_rtol

        self.model = None
        self.cantera_simulation = None
        self.inert_list = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'N2']
        self.inert_index_list = list()
        self.spc_identifier_lookup, self.rxn_identifier_lookup = dict(), dict()
        self.num_ct_reactions = None
        self.num_ct_species = None
        self.T_list, self.P_list, self.reaction_time_list = list(), list(), list()
        self.idt_dict = dict()

        self.set_up()
        self.radical_label = self.determine_radical_label()
        print(f'radical label: {self.radical_label}')

    def set_up(self):
        """
        Read in the Cantera input file and set up attributes.
        """
        self.model = ct.Solution(infile=self.paths['cantera annotated'])
        self.num_ct_reactions = len(self.model.reactions())
        self.num_ct_species = len(self.model.species())

        for i, species in enumerate(self.model.species()):
            if species.name in self.inert_list:
                self.inert_index_list.append(i)
        for i, spc in enumerate(self.model.species()):
            self.spc_identifier_lookup[spc.name] = i
        for i, rxn in enumerate(self.model.reactions()):
            self.rxn_identifier_lookup[rxn.equation] = i
        self.species_names_without_indices = [self.model.species()[i].name.split('(')[0] for i in range(self.num_ct_species)]  # Todo: what about HNO(T)(21)? make more robust

        self.T_list = ([self.rmg['reactors'][0]['T']], 'K')
        self.P_list = ([self.rmg['reactors'][0]['P']], 'bar')

    def simulate(self):
        """
        Simulate the mechanism and compute ignition delay times.
        """
        self.logger.info('Running a simulation using CanteraIDT...')

        equivalence_ratios, concentration_combinations = self.get_concentration_combinations()
        reactor_idt_dict = dict()
        for r, reactor in enumerate(self.rmg['reactors']):
            T_list, P_list = get_t_and_p_lists(reactor)
            if equivalence_ratios is not None and concentration_combinations is not None:
                for i, X in enumerate(concentration_combinations):
                    for P in P_list:
                        for T in T_list:
                            self.model.TPX = T, P * 1e5, X
                            self.idt_dict[(equivalence_ratios[i], P, T)] = \
                                self.simulate_idt(fig_name=f'R{r}_{equivalence_ratios[i]}_{round(P,2)}_bar_{round(T,2)}_K.png')
            else:
                X = {spc['label']: spc['concentration'] for spc in self.rmg['species'] if spc['concentration']}
                for P in P_list:
                    for T in T_list:
                        self.model.TPX = T, P * 1e5, X
                        self.idt_dict[(0, P, T)] = \
                            self.simulate_idt(fig_name=f'R{r}_{round(P,2)}_bar_{round(T,2)}_K.png')
            if len(T_list) >= 3:
                plot_idt_vs_temperature(self.idt_dict, figs_path=self.paths['figs'], reactor_index=r)
            reactor_idt_dict[r] = self.idt_dict
        save_yaml_file(os.path.join(self.paths['figs'], 'idt_dict.yaml'), reactor_idt_dict)

    def simulate_idt(self, energy: str = 'on',
                     fig_name: str = 'idt.png',
                     ) -> Optional[float]:
        """
        Simulate an IdealGasReactor to find the IDT value.

        Returns:
            Optional[float]: The IDT in seconds.
        """
        reactor = ct.IdealGasReactor(contents=self.model, energy=energy)
        net = ct.ReactorNet([reactor])

        # for i in range(self.model.n_reactions):
        #     reactor.add_sensitivity_reaction(i)

        net.atol = self.atol
        net.rtol = self.rtol
        net.atol_sensitivity = self.sa_atol
        net.rtol_sensitivity = self.sa_rtol

        time_history = ct.SolutionArray(self.model, extra='t')
        t, counter, max_idt = 0, 0, 1000

        while t < max_idt:
            t = net.step()
            time_history.append(reactor.thermo.state, t=t)
            if counter % 100 == 0:
                concentrations = np.asarray([x[0] for x in time_history(self.radical_label).X], dtype=np.float32)
                max_c_idx = np.argmax(concentrations)
                if concentrations[max_c_idx] > concentrations[-1] * 1.2 \
                        and len(concentrations) > max_c_idx * 1.1 \
                        and time_history.t[-1] > time_history.t[max_c_idx] * 1.2:
                    break
            counter += 1

        idt = compute_idt(time_history, self.radical_label, figs_path=self.paths['figs'], fig_name=fig_name)
        return idt

    def determine_radical_label(self) -> str:
        """
        Determine the label of the prominent ignition radical: OH if present, else H (e.g., in hydrazine).

        Returns:
            str: The label of the prominent ignition radical.
        """
        h, oh = None, None
        for i, species in enumerate(self.species_names_without_indices):
            if species.lower() == 'oh':
                oh = self.model.species()[i].name
            if species.lower() == 'h':
                h = self.model.species()[i].name
            if oh is not None:
                break
        return oh or h

    def get_cantera_species_label(self, rmg_label: str) -> Optional[str]:
        """
        Determine the corresponding label in the Cantera model for an RMG input species.

        Returns:
            Optional[str]: The label of the corresponding Cantera species.
        """
        for i, label in enumerate(self.species_names_without_indices):
            if label == rmg_label:
                return self.model.species()[i].name
        return None

    def get_concentration_combinations(self) -> Tuple[Optional[List[float]], Optional[List[dict]]]:
        """
        Get concentration combinations according to the equivalence ratios of the fuel.

        Returns: Tuple[Optional[List[float]], Optional[List[dict]]]
            - List[float]: List of equivalence ratios.
            - List[dict]: List of dictionaries, each is a combination of concentrations for an IDT simulation.
        """
        objects = determine_concentrations_by_equivalence_ratios(species=self.rmg['species'])
        fuel_idx = None
        for i, spc in enumerate(self.rmg['species']):
            if 'role' in spc.keys() and spc['role'] == 'fuel':
                fuel_idx = i
                break
        if fuel_idx is None:
            return None, None
        equivalence_ratios = self.rmg['species'][fuel_idx]['equivalence_ratios']
        concentration_combinations = list()
        for i, _ in enumerate(equivalence_ratios):
            concentration_dict = dict()
            for spc in self.rmg['species']:
                if spc['role'] is None:
                    cantera_label = self.get_cantera_species_label(spc['label'])
                    if cantera_label is not None:
                        concentration_dict[cantera_label] = spc['concentration']
                if spc['role'] == 'fuel':
                    cantera_label = self.get_cantera_species_label(spc['label'])
                    if cantera_label is not None:
                        concentration_dict[cantera_label] = objects['fuel']['concentration']
                elif spc['role'] in ['oxygen', 'nitrogen'] and len(objects[spc['role']]['concentration']) == len(equivalence_ratios):
                    cantera_label = self.get_cantera_species_label(spc['label'])
                    if cantera_label is not None:
                        concentration_dict[cantera_label] = objects[spc['role']]['concentration'][i]
            concentration_combinations.append(concentration_dict)
        return equivalence_ratios, concentration_combinations

    def get_sa_coefficients(self):
        """
        Obtain the SA coefficients.

        Returns:
             sa_dict (dict): a SA dictionary, whose structure is given in the docstring for T3/t3/main.py
        """
        pass
        # sa_dict = {'kinetics': dict(), 'thermo': dict(), 'time': list()}
        #
        # for condition_data in self.all_data:
        #     time, data_list, reaction_sensitivity_data, thermodynamic_sensitivity_data = condition_data
        #     sa_dict['time'] = time.data
        #
        #     # extract kinetic SA
        #     for rxn in reaction_sensitivity_data:
        #         # for kinetics, get `ethane(1)` from `dln[ethane(1)]/dln[k8]: H(6)+ethane(1)=H2(12)+C2H5(5)`
        #         observable_label = rxn.label.split('[')[1].split(']')[0]
        #         if observable_label not in sa_dict['kinetics']:
        #             sa_dict['kinetics'][observable_label] = dict()
        #         # for kinetics, get k8 from `dln[ethane(1)]/dln[k8]: H(6)+ethane(1)=H2(12)+C2H5(5)` then only extract 8
        #         parameter = rxn.label.split('[')[2].split(']')[0]
        #         parameter = int(parameter[1:])
        #         sa_dict['kinetics'][observable_label][parameter] = rxn.data
        #
        #     # extract thermo SA
        #     for spc in thermodynamic_sensitivity_data:
        #         # for thermo get 'C2H4(8)' from `dln[ethane(1)]/dH[C2H4(8)]`
        #         observable_label = spc.label.split('[')[1].split(']')[0]
        #         if observable_label not in sa_dict['thermo']:
        #             sa_dict['thermo'][observable_label] = dict()
        #         # for thermo get 'C2H4(8)' from `dln[ethane(1)]/dH[C2H4(8)]`
        #         parameter = spc.label.split('[')[2].split(']')[0]
        #         sa_dict['thermo'][observable_label][parameter] = spc.data
        #
        # return sa_dict


def compute_idt(time_history: ct.SolutionArray,
                radical_label: str,
                figs_path: str,
                fig_name: str,
                ) -> Optional[float]:
    """
    Finds the ignition point by approximating dT/dt as a first order forward difference
    and then finds the point of maximum slope.

    Args:
        time_history (ct.SolutionArray): Cantera solution array.
        radical_label (str): The label of the prominent ignition radical.
        figs_path (str): The path to the figures' directory.
        fig_name (str): The name of the figure to save.

    Returns:
        Optional[float]: The IDT in seconds.

    Todo:
        - solve possible noice in dc/dt. ** add IDT(1000/T) figure.
    """
    figs_path = os.path.join(figs_path, 'IDTs')
    if not os.path.isdir(figs_path):
        os.makedirs(figs_path)
    times = time_history.t
    concentration = np.asarray([x[0] for x in time_history(radical_label).X], dtype=np.float32)
    dc_dt = np.diff(concentration) / np.diff(times)
    idt_index_dc_dt = np.argmax(dc_dt)
    idt_index_c = np.argmax(concentration)
    idt = times[idt_index_dc_dt]
    if idt_index_dc_dt > len(times) - 10 or idt < 1e-12 or max(concentration) < concentration[0] * 100:
        return None
    try:
        plt.plot(times, concentration)
        plt.plot(times[idt_index_dc_dt], concentration[idt_index_dc_dt], 'o')
        plt.xlabel('Time (s)')
        plt.ylabel(f'[{radical_label}]')
        plt.title(f'IDT = {idt:.2e} s')
        if times[idt_index_c] > times[idt_index_dc_dt] * 2.5:
            x_min = min(times[idt_index_dc_dt] * 0.8, times[np.argmax(concentration)] * 0.1)
            x_max = max(times[idt_index_dc_dt] * 1.2, times[np.argmax(concentration)] * 1.5)
            if times[idt_index_dc_dt] < x_max * 0.05:
                x_min = 0
        else:
            x_min, x_max = times[idt_index_dc_dt] * 0.8, times[idt_index_dc_dt] * 1.2
        plt.xlim(x_min, x_max)
        plt.savefig(os.path.join(figs_path, fig_name))
        plt.close()
    except (AttributeError, ValueError):
        pass
    return idt


def get_t_and_p_lists(reactor: dict) -> Tuple[List[float], List[float]]:
    """
    Get temperature and pressure lists for a single RMG reactor.

    Args:
        reactor (dict): RMG reactor dictionary.

    Returns: Tuple[List[float], List[float]]
        - List[float]: List of temperatures.
        - List[float]: List of pressures.
    """
    if isinstance(reactor['T'], (int, float)):
        T_list = [reactor['T']]
    else:
        inverse_ts = np.linspace(1 / reactor['T'][1],
                                 1 / reactor['T'][0], num=100)  # 15 inverse T points
        T_list = [1 / inverse_t for inverse_t in inverse_ts[::-1]]
    if isinstance(reactor['P'], (int, float)):
        P_list = [reactor['P']]
    else:
        base = 10
        log_p = np.linspace(math.log(reactor['P'][0], base),
                            math.log(reactor['P'][1], base), num=3) # 3 pressure in log10 space
        P_list = [base ** p for p in log_p]
    return T_list, P_list


def plot_idt_vs_temperature(idt_dict: dict,
                            figs_path: str,
                            reactor_index: int = 0,
                            exp_data: tuple = None,
                            ) -> None:
    """
    Plot IDT vs. 1000/T per phi and P condition combination.

    idt_dict[(equivalence_ratios[i], P, T)]

    Args:
        idt_dict (dict): A dictionary containing IDT values.
        figs_path (str): The path to the figures' directory.
        reactor_index (int, optional): The reactor index.
        exp_data (tuple): Experimental data in lists for IDT [sec] and 1000/T [K] 
    """
    figs_path = os.path.join(figs_path, 'IDT_vs_T')
    if not os.path.isdir(figs_path):
        os.makedirs(figs_path)
    data = get_idt_per_phi_p_condition(idt_dict)
    
    for phi, phi_data in data.items():
        for p, phi_p_data in phi_data.items():
            fig_name = f'R{reactor_index}_{phi}_{round(p,2)}_bar.png'
            try:
                fig, ax = plt.subplots()
                ax.set_xlabel('1000/T (1/K)')
                ax.set_ylabel('IDT (s)')
                ax.set_title(f'IDT vs. 1000/T, phi = {phi}, P = {p} bar')
                ax.scatter(phi_p_data.keys(), phi_p_data.values(),label='simulation', color='blue',marker="o")
                ax.set_yscale('log')
                #Add experimental data specified by the user
                if exp_data is not None:
                    ax.scatter(exp_data[0], exp_data[1], label='experiment', color='orange',marker="D")
                    ax.set_yscale('log')
                ax.legend(loc='lower right')
                fig.savefig(os.path.join(figs_path, fig_name))
            except (AttributeError, ValueError):
                pass


def get_idt_per_phi_p_condition(idt_dict: dict) -> dict:
    """
    Get IDT values per phi, P, and 1000/T condition combination.

    Args:
        idt_dict (dict): A dictionary containing IDT values.

    Returns:
        dict: A dictionary containing IDT values per phi and P condition combination corresponding to the 1000/T list.
    """
    data = dict()
    for triple_key in idt_dict.keys():
        phi, p, _ = triple_key
        if phi not in data.keys():
            data[phi] = dict()
        if p not in data[phi].keys():
            data[phi][p] = dict()
    for triple_key, value in idt_dict.items():
        data[triple_key[0]][triple_key[1]][1000 / triple_key[2]] = value
    return data


register_simulate_adapter("CanteraIDT", CanteraIDT)
