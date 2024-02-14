"""
Cantera Simulator Adapter module for ignition delay time (IDT) calculations.
Used to run mechanism analysis with Cantera for IDT.
"""

from typing import List, Optional, Tuple
import concurrent.futures as cf
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import os
import traceback

import cantera as ct
import math
import numpy as np

from arc.common import read_yaml_file, save_yaml_file

from t3.common import determine_concentrations_by_equivalence_ratios, remove_numeric_parentheses
from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter


DELTA_H = 0.1  # a +|- 0.1 kJ/mol perturbation
DELTA_K = 0.05  # a *|/ (1 + 5%) perturbation
R = 8.31446261815324  # J/(mol*K)
EA_UNIT_CONVERSION = {'J/mol': 1, 'kJ/mol': 1e+3, 'cal/mol': 4.184, 'kcal/mol': 4.184e+3}
P_UNIT_CONVERSION = {'bar': 1, 'atm': 1.01325, 'Pa': 1e-5}


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
        self.reactor_idt_dict = None
        self.inert_list = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'N2']
        self.inert_index_list = list()
        self.spc_identifier_lookup, self.rxn_identifier_lookup = dict(), dict()
        self.num_ct_reactions = None
        self.num_ct_species = None
        self.T_list, self.P_list, self.reaction_time_list, self.species_names_without_indices = list(), list(), list(), list()
        self.idt_sa_dict = dict()

        self.set_up()
        self.radical_label = self.determine_radical_label()

    def set_up(self):
        """
        Read in the Cantera input file and set up attributes.
        """
        self.model = ct.Solution(infile=self.paths['cantera annotated'])
        self.num_ct_reactions = len(self.model.reactions())
        self.num_ct_species = len(self.model.species())
        self.inert_index_list = [i for i, species in enumerate(self.model.species()) if species.name in self.inert_list]
        self.spc_identifier_lookup = {spc.name: i for i, spc in enumerate(self.model.species())}
        self.rxn_identifier_lookup = {rxn.equation: i for i, rxn in enumerate(self.model.reactions())}
        self.species_names_without_indices = [remove_numeric_parentheses(self.model.species()[i].name)
                                              for i in range(self.num_ct_species)]

    def simulate(self):
        """
        Simulate the mechanism and compute ignition delay times.
        """
        if self.logger is not None:
            self.logger.info('Running a simulation using CanteraIDT...')
        self.reactor_idt_dict = self.simulate_idt_for_all_reactors(save_yaml=True)

    def simulate_idt_for_all_reactors(self,
                                      save_yaml: bool = True,
                                      save_fig: bool = True,
                                      energy: str = 'on',
                                      max_idt: float = 1.0,
                                      ) -> dict:
        """
        Simulate an IdealGasReactor to find the IDT value for all working points with parallelization.

        Args:
            save_yaml (bool, optional): Save the IDT dictionary to a YAML file.
            save_fig (bool, optional): Whether to save the figures.
            energy (str, optional): The energy mode to use. Options are 'off', 'on', 'off after ignition'.
            max_idt (int, optional): Maximum IDT in seconds.

        Returns:
            dict: The IDT dictionary.
        """
        idt_dict = dict()
        infile = self.paths['cantera annotated']
        equivalence_ratios, concentration_combinations = self.get_concentration_combinations()
        reactor_idt_dict = dict()
        for r, reactor in enumerate(self.rmg['reactors']):
            T_list, P_list = get_t_and_p_lists(reactor)
            if equivalence_ratios is not None and concentration_combinations is not None:
                combinations = [
                    (r, T, P, X, equivalence_ratios[idx], infile, save_fig, energy, max_idt)
                    for T in T_list
                    for P in P_list
                    for idx, X in enumerate(concentration_combinations)
                ]
            else:
                X = {spc['label']: spc['concentration'] for spc in self.rmg['species'] if spc['concentration']}
                phi = None if equivalence_ratios is None else equivalence_ratios[0]
                combinations = [
                    (r, T, P, X, phi, infile, save_fig, energy, max_idt)
                    for T in T_list for P in P_list
                ]
            results = list()
            for args in combinations:
                results.append(self.simulate_idt_for_a_point(*args))
            for i, (_, T, P, X, _, _, _, _, _) in enumerate(combinations):
                if results[i] is not None:
                    if equivalence_ratios is not None:
                        eq_ratio = equivalence_ratios[i % len(equivalence_ratios)]
                        idt_dict.setdefault(eq_ratio, {}).setdefault(P, {})[T] = results[i]
                    else:
                        idt_dict.setdefault(0, {}).setdefault(P, {})[T] = results[i]

            if len(T_list) >= 3 and save_fig:
                plot_idt_vs_temperature(idt_dict, figs_path=self.paths['figs'], reactor_index=r)
            reactor_idt_dict[r] = idt_dict
        if save_yaml:
            save_yaml_file(os.path.join(self.paths['figs'], 'idt_dict.yaml'), reactor_idt_dict)
        return reactor_idt_dict

    def simulate_idt_for_a_point(self,
                                 r: int,
                                 t: float,
                                 p: float,
                                 x: dict,
                                 phi: Optional[float],
                                 infile: str,
                                 save_fig: bool = True,
                                 energy: str = 'on',
                                 max_idt: int = 1,
                                 ) -> Optional[float]:
        """
        Simulate an IdealGasReactor to find the IDT value for a specific working point.

        Args:
            r (int): The reactor index.
            t (float): Temperature in K.
            p (float): Pressure in bar.
            x (dict): Concentrations dictionary.
            phi (Optional[float]): Equivalence ratio.
            infile (str): The path to the Cantera input file.
            save_fig (bool, optional): Save the figure.
            energy (str, optional): The energy mode to use. Options are 'off', 'on', 'off after ignition'.
            max_idt (int, optional): Maximum IDT in seconds.

        Returns:
            Optional[float]: The IDT in seconds.
        """
        fig_name = f'R{r}_{phi}_{p:.2f}_bar_{t:.2f}_K.png' if phi is not None else f'R{r}_{p:.2f}_bar_{t:.2f}_K.png'
        model = ct.Solution(infile=infile)
        model.TPX = t, p * 1e5, x
        reactor = ct.IdealGasReactor(contents=model, energy=energy)
        net = ct.ReactorNet([reactor])
        net.atol, net.rtol, net.atol_sensitivity, net.rtol_sensitivity = self.atol, self.rtol, self.sa_atol, self.sa_rtol
        time_history = ct.SolutionArray(model, extra='t')
        t_, counter = 0, 0
        while t_ < max_idt:
            t_ = net.step()
            time_history.append(reactor.thermo.state, t=t_)
            if counter % 100 == 0:
                concentrations = np.asarray([x[0] for x in time_history(self.radical_label).X], dtype=np.float32)
                max_c_idx = np.argmax(concentrations)
                if concentrations[max_c_idx] > concentrations[-1] * 1.2 \
                        and len(concentrations) > max_c_idx * 1.1 \
                        and time_history.t[-1] > time_history.t[max_c_idx] * 1.2:
                    break
            counter += 1
        idt = compute_idt(time_history=time_history,
                          radical_label=self.radical_label,
                          figs_path=self.paths['figs'],
                          fig_name=fig_name if save_fig else None,
                          )
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
            if label == remove_numeric_parentheses(rmg_label):
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
                    if cantera_label is not None and spc['concentration'] != 0:
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

    def get_sa_coefficients(self,
                            top_SA_species: int = 10,
                            top_SA_reactions: int = 10,
                            max_workers: int = 24,
                            save_yaml: bool = True,
                            ) -> Optional[dict]:
        """
        Obtain the SA coefficients for IDT using brute force.

        Args:
            top_SA_species (int, optional): The number of top sensitive species to return.
            top_SA_reactions (int, optional): The number of top sensitive reactions to return.
            max_workers (int, optional): The maximal number of workers to use for parallel processing.
            save_yaml (bool, optional): Save the SA dictionary to a YAML file.

        Returns:
             sa_dict (Optional[dict]): a SA dictionary, whose structure is given in the docstring for T3/t3/main.py
        """
        if not self.reactor_idt_dict:
            self.simulate()
        sa_dict = {'thermo': {'IDT': dict()}, 'kinetics': {'IDT': dict()}}
        tasks = [('thermo', i) for i in range(self.num_ct_species)] + \
                [('kinetics', i) for i in range(self.num_ct_reactions)]
        with cf.ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {executor.submit(worker,
                                              task,
                                              self.paths['cantera annotated'],
                                              self.paths['SA'],
                                              self.t3,
                                              self.paths,
                                              self.rmg,
                                              self.logger,
                                              max_idt=1.0,
                                              ): task for task in tasks}
            for future in cf.as_completed(future_to_task):
                task = future_to_task[future]
                try:
                    idt_dict = future.result()
                    sa_dict[task[0]]['IDT'][task[1]] = idt_dict
                except Exception as e:
                    self.logger.error(f"Task {task} generated an exception:\n{e}\n{traceback.format_exc()}")

        idt_sa_dict = compute_idt_sa(reactor_idt_dict=self.reactor_idt_dict,
                                     perturbed_idt_dict=sa_dict,
                                     )
        self.idt_sa_dict = get_top_sa_coefficients(idt_sa_dict=idt_sa_dict,
                                                   top_species=top_SA_species,
                                                   top_reactions=top_SA_reactions)
        if save_yaml:
            save_yaml_file(path=self.paths['SA IDT dict'], content=idt_sa_dict)
            save_yaml_file(path=self.paths['SA IDT dict top X'], content=self.idt_sa_dict)
        return self.idt_sa_dict


def worker(task: tuple,
           model_path: str,
           work_dir: str,
           t3: dict,
           paths: dict,
           rmg: dict,
           logger: Logger,
           max_idt: float = 1.0,
           ) -> dict:
    """
    Worker function to perturb a model parameter and run IDT calculation for SA computation.

    Args:
        task (tuple): Tuple containing ('thermo', index) or ('kinetics', index).
        model_path (str): Path to the Cantera input file.
        work_dir (str): The path to the working directory.
        t3 (dict): The T3.t3 attribute, which is a dictionary containing the t3 directives.
        paths (dict): The T3.paths attribute, which is a dictionary containing relevant paths.
        rmg (dict): The T3.rmg attribute, which is a dictionary containing the rmg block from the input yaml or API.
        logger (Logger): Instance of T3's Logger class.
        max_idt (float, optional): Maximum IDT in seconds.

    Returns:
        dict: The IDT dictionary.
    """
    if not isinstance(task, tuple) or len(task) != 2:
        raise ValueError(f'Task must be a tuple of length 2, got: {task}')
    perturbed_model_path = os.path.join(work_dir, f'{task[0]}_{task[1]}.yaml')
    kind, index = task
    try:
        if kind == 'thermo':
            success = perturb_enthalpy(original_path=model_path,
                                       perturbed_path=perturbed_model_path,
                                       species_index=index,
                                       logger=logger,
                                       )
        elif kind == 'kinetics':
            success = perturb_reaction_rate_coefficient(original_path=model_path,
                                                        perturbed_path=perturbed_model_path,
                                                        reaction_index=index,
                                                        logger=logger,
                                                        )
        else:
            raise ValueError(f'Unknown task type {kind}')
        if success:
            new_paths = {**paths, 'cantera annotated': perturbed_model_path}
            ct_idt_adapter = CanteraIDT(t3=t3, paths=new_paths, rmg=rmg, logger=logger)
            return ct_idt_adapter.simulate_idt_for_all_reactors(save_yaml=False, save_fig=False, energy='on', max_idt=max_idt)
        else:
            return dict()
    finally:
        os.remove(perturbed_model_path) if os.path.exists(perturbed_model_path) else None


def compute_idt(time_history: ct.SolutionArray,
                radical_label: str,
                figs_path: Optional[str] = None,
                fig_name: Optional[str] = None,
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
        - solve possible noise in dc/dt. ** add IDT(1000/T) figure.
    """
    figs_path = os.path.join(figs_path, 'IDTs')
    if not os.path.isdir(figs_path):
        os.makedirs(figs_path, exist_ok=True)
    times = time_history.t
    concentration = np.asarray([x[0] for x in time_history(radical_label).X], dtype=np.float32)
    if all(c == 0 for c in concentration):
        return None
    dc_dt = np.diff(concentration) / np.diff(times)
    idt_index_dc_dt = np.argmax(dc_dt)
    idt_index_c = np.argmax(concentration)
    idt = float(times[idt_index_dc_dt])
    if idt_index_dc_dt > len(times) - 10 or idt < 1e-12 or max(concentration) < concentration[0] * 100:
        return None
    if figs_path is not None and fig_name is not None:
        try:
            plt.plot(times, concentration)
            plt.plot(times[idt_index_dc_dt], concentration[idt_index_dc_dt], 'o')
            plt.xlabel('Time (s)')
            plt.ylabel(f'[{radical_label}]')
            plt.title(f'IDT = {idt:.2e} s')
            ax = plt.gca()
            ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.xaxis.get_major_formatter().set_scientific(True)
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


def get_t_and_p_lists(reactor: dict,
                      num_t_points: int = 25,
                      ) -> Tuple[List[float], List[float]]:
    """
    Get temperature and pressure lists for a single RMG reactor.

    Args:
        reactor (dict): RMG reactor dictionary.
        num_t_points (int, optional): Number of temperature points to consider.

    Returns: Tuple[List[float], List[float]]
        - List[float]: List of temperatures.
        - List[float]: List of pressures.
    """
    if isinstance(reactor['T'], (int, float)):
        T_list = [reactor['T']]
    elif len(reactor['T']) == 2:
        inverse_ts = np.linspace(1 / reactor['T'][1],
                                 1 / reactor['T'][0], num=min(int(abs(reactor['T'][1] - reactor['T'][0]) / 10), num_t_points))
        T_list = [1 / inverse_t for inverse_t in inverse_ts[::-1]]
    else:
        T_list = [float(t) for t in reactor['T']]
    if isinstance(reactor['P'], (int, float)):
        P_list = [reactor['P']]
    elif len(reactor['P']) == 2:
        base = 10
        log_p = np.linspace(math.log(reactor['P'][0], base), math.log(reactor['P'][1], base), num=3) # 3 pressures in log10 space
        P_list = [base ** p for p in log_p]
    else:
        P_list = [float(p) for p in reactor['P']]
    T_list = [float(t) for t in T_list]
    P_list = [float(p) for p in P_list]
    return T_list, P_list


def plot_idt_vs_temperature(idt_dict: dict,
                            figs_path: str,
                            reactor_index: int = 0,
                            exp_data: dict = None,
                            ) -> None:
    """
    Plot IDT vs. 1000/T per phi and P condition combination.
    If exp_data is provided, plot the experimental data as well and only consider phi and P values that appear in it.

    idt_dict[equivalence_ratios[i]][P][T]

    Args:
        idt_dict (dict): A dictionary containing IDT values.
        figs_path (str): The path to the figures' directory.
        reactor_index (int, optional): The reactor index.
        exp_data (dict): Experimental data in the same format as idt_dict, IDT units are in s.
    """
    figs_path = os.path.join(figs_path, 'IDT_vs_T')
    if not os.path.isdir(figs_path):
        os.makedirs(figs_path)
    for phi, phi_data in idt_dict.items():
        if exp_data is not None and phi not in exp_data:
            continue
        for p, phi_p_data in phi_data.items():
            if exp_data is not None and p not in exp_data[phi]:
                continue
            fig_name = f'R{reactor_index}_{phi}_{p:.2f}_bar.png'
            try:
                fig, ax = plt.subplots()
                ax.set_xlabel('1000/T (1/K)')
                ax.set_ylabel('IDT (s)')
                ax.set_title(f'IDT vs. 1000/T, phi = {phi}, P = {p:.2f} bar')
                ax.scatter([1000 / t for t in phi_p_data.keys()], phi_p_data.values(), label='simulation', color='blue',
                           marker='o', linestyle='-')
                ax.set_yscale('log')
                if exp_data is not None:
                    ax.scatter([1000 / t for t in exp_data[phi][p].keys()], [e * 1e-6 for e in exp_data[phi][p].values()],
                               label='experiment', color='orange', marker="D")
                    ax.set_yscale('log')
                ax.legend(loc='lower right')
                fig.savefig(os.path.join(figs_path, fig_name))
            except (AttributeError, ValueError):
                pass


def perturb_enthalpy(original_path: str,
                     perturbed_path: str,
                     species_index: int,
                     logger: Optional[Logger] = None,
                     ) -> bool:
    """
    Perturb the enthalpy of a single species in the model.

    Args:
        original_path (str): Path to the original Cantera input file.
        perturbed_path (str): Path to the perturbed Cantera input file.
        species_index (int): The index of the species to perturb.
        logger (Logger, optional): Instance of T3's Logger class.

    Returns:
        bool: ``True`` if the perturbation was successful, ``False`` otherwise.
    """
    content = read_yaml_file(original_path)
    for i in range(len(content['species'])):
        if i == species_index:
            if content['species'][i]['thermo']['model'] != 'NASA7':
                if logger is not None:
                    logger.warning(f"Species '{content['species'][i]['name']}' does not use the expected NASA "
                                   f"polynomials format for thermo, not perturbing it.")
                return False
            content['species'][i]['thermo']['data'][0][5] += DELTA_H * 1e3 / R
            if len(content['species'][i]['thermo']['data']) == 2:
                content['species'][i]['thermo']['data'][1][5] += DELTA_H * 1e3 / R
            break
    save_yaml_file(perturbed_path, content)
    return True


def perturb_reaction_rate_coefficient(original_path: str,
                                      perturbed_path: str,
                                      reaction_index: int,
                                      logger: Optional[Logger] = None,
                                      ) -> bool:
    """
    Perturb the rate coefficient of a single reaction in the model by modifying the pre-exponential factor.

    Args:
        original_path (str): Path to the original Cantera input file.
        perturbed_path (str): Path to the perturbed Cantera input file.
        reaction_index (int): The index of the reaction to perturb.
        logger (Logger, optional): Instance of T3's Logger class.

    Returns:
        bool: ``True`` if the perturbation was successful, ``False`` otherwise.
    """
    content = read_yaml_file(original_path)
    for i in range(len(content['reactions'])):
        if i == reaction_index:
            if 'rate-constant' in content['reactions'][i] and 'A' in content['reactions'][i]['rate-constant']:
                # relevant for Arrhenius and three-body rate expressions
                content['reactions'][i]['rate-constant']['A'] *= (1 + DELTA_K)
            elif 'type' in content['reactions'][i] and content['reactions'][i]['type'] == 'falloff':
                content['reactions'][i]['low-P-rate-constant']['A'] *= (1 + DELTA_K)
                content['reactions'][i]['high-P-rate-constant']['A'] *= (1 + DELTA_K)
            elif 'type' in content['reactions'][i] and content['reactions'][i]['type'] == 'pressure-dependent-Arrhenius':
                for j in range(len(content['reactions'][i]['rate-constants'])):
                    content['reactions'][i]['rate-constants'][j]['A'] *= (1 + DELTA_K)
            elif 'type' in content['reactions'][i] and content['reactions'][i]['type'] == 'Chebyshev':
                content['reactions'][i]['data'][0][0] += math.log10(1 + DELTA_K)
            else:
                if logger is not None:
                    logger.warning(f"Reaction '{content['reactions'][i]['equation']}' does not use the expected "
                                   f"Arrhenius rate expression, not perturbing it.")
                return False
            break
    save_yaml_file(perturbed_path, content)
    return True


def compute_idt_sa(reactor_idt_dict: dict,
                   perturbed_idt_dict: dict,
                   ) -> dict:
    """
    Compute the sensitivity analysis coefficients for IDT.
    for kinetics use: dln(IDT)/dln(k_i) = [k_i * d(IDT)] / [IDT * d(k_i)], dimensionless
    for thermo use: dln(IDT)/d(H298) = d(IDT) / [IDT * d(H298)], units: mol/kJ

    Args:
        reactor_idt_dict (dict): The IDT dictionary.
        perturbed_idt_dict (dict): The perturbed IDT dictionary.

    Returns:
        dict: The IDT SA dictionary structured as idt_sa_dict[kind]['IDT'][r][phi][p][t].
    """
    idt_sa_dict = dict()
    for token in ['thermo', 'kinetics']:
        idt_sa_dict[token] = {'IDT': dict()}
        for r, reactor_idt_data in reactor_idt_dict.items():
            idt_sa_dict[token]['IDT'][r] = dict()
            for phi, phi_data in reactor_idt_data.items():
                idt_sa_dict[token]['IDT'][r][phi] = dict()
                for p, p_data in phi_data.items():
                    idt_sa_dict[token]['IDT'][r][phi][p] = dict()
                    for t, idt in p_data.items():
                        idt_sa_dict[token]['IDT'][r][phi][p][t] = dict()
                        if idt is None:
                            continue
                        for index, perturbed_idt_data in perturbed_idt_dict[token]['IDT'].items():
                            perturbed_idt_value = perturbed_idt_data[r][phi][p][t]
                            if perturbed_idt_value is None:
                                continue
                            delta_idt = perturbed_idt_value - idt
                            if token == 'kinetics':
                                sa_coeff = delta_idt / (idt * DELTA_K)  # [k_i * d(IDT)] / [IDT * d(k_i)]; k_i / d(k_i) = 1 / DELTA_K
                            elif token == 'thermo':
                                sa_coeff = delta_idt / (idt * DELTA_H)  # d(IDT) / [IDT * d(H298)]; in mol/kJ
                            else:
                                raise ValueError(f'Unknown token: {token}')
                            idt_sa_dict[token]['IDT'][r][phi][p][t][index] = sa_coeff
    return idt_sa_dict


def get_top_sa_coefficients(idt_sa_dict: dict,
                            top_species: int = 10,
                            top_reactions: int = 10,
                            ) -> dict:
    """
    Get the top sensitivity coefficients for IDT based on absolute values.

    Args:
        idt_sa_dict (dict): The IDT SA dictionary.
        top_species (int, optional): The number of top sensitive species to return.
        top_reactions (int, optional): The number of top sensitive reactions to return.

    Returns:
        dict: The top sensitivity coefficients dictionary.
    """
    top_sa_dict = dict()
    for token in ['thermo', 'kinetics']:
        n = top_species if token == 'thermo' else top_reactions
        top_sa_dict[token] = {'IDT': dict()}
        for r, reactor_idt_data in idt_sa_dict[token]['IDT'].items():
            top_sa_dict[token]['IDT'][r] = dict()
            for phi, phi_data in reactor_idt_data.items():
                top_sa_dict[token]['IDT'][r][phi] = dict()
                for p, p_data in phi_data.items():
                    top_sa_dict[token]['IDT'][r][phi][p] = dict()
                    for t, idt_data in p_data.items():
                        top_sa_dict[token]['IDT'][r][phi][p][t] = dict()
                        top_species_indices = sorted(idt_data, key=lambda x: abs(idt_data[x]), reverse=True)[:n]
                        for index in top_species_indices:
                            top_sa_dict[token]['IDT'][r][phi][p][t][index] = idt_data[index]
    return top_sa_dict







def get_h298(model: ct.Solution, species_index: int) -> float:
    """
    Get the enthalpy of formation at 298 K for a species in kJ/mol.

    Args:
        model (ct.Solution): Cantera solution object for the mechanism.
        species_index (int): The index of the species.

    Returns:
        float: The enthalpy of formation at 298 K in kJ/mol.
    """
    model.TP = 298, 1e5
    return model.standard_enthalpies_RT[species_index] * ct.gas_constant * 298 / 1e6  # J/kmol to kJ/mol


def calculate_arrhenius_rate_coefficient(A: float, n: float, Ea: float, T: float, Ea_units: str) -> float:
    """
    Calculate the Arrhenius rate coefficient.

    Args:
        A (float): Pre-exponential factor.
        n (float): Temperature exponent.
        Ea (float): Activation energy in J/mol.
        T (float): Temperature in Kelvin.
        Ea_units (str): Units of the rate coefficient.

    Returns:
        float: The rate coefficient at the specified temperature.
    """
    if Ea_units not in EA_UNIT_CONVERSION:
        raise ValueError(f"Unsupported Ea units: {Ea_units}")
    return A * (T ** n) * math.exp(-1 * (Ea * EA_UNIT_CONVERSION[Ea_units]) / (R * T))


def calculate_troe_rate_coefficient(reaction_data: dict, T: float, P: float, Ea_units: str) -> float:
    """
    Calculate the Troe (or Lindemann) rate coefficient.

    Args:
        reaction_data (dict): Dictionary containing reaction parameters (the Cantera YAML dict of a single reaction).
        T (float): Temperature in Kelvin.
        P (float): Pressure in bar.
        Ea_units (str): Units of the activation energy.

    Returns:
        float: The rate coefficient at the specified temperature and pressure.
    """
    high_params = reaction_data['high-P-rate-constant']
    low_params = reaction_data['low-P-rate-constant']
    k0 = calculate_arrhenius_rate_coefficient(A=low_params['A'], n=low_params['b'], Ea=low_params['Ea'], T=T, Ea_units=Ea_units)
    kinf = calculate_arrhenius_rate_coefficient(A=high_params['A'], n=high_params['b'], Ea=high_params['Ea'], T=T, Ea_units=Ea_units)
    C = P / (R * T * 10)  # bath gas concentration in mol/cm^3
    Pr = k0 * C / kinf
    F = 1
    if 'Troe' in reaction_data:
        troe_params = reaction_data['Troe']
        alpha = troe_params.get('A', 0.0)
        T1, T2, T3 = troe_params.get('T1', 1e+30), troe_params.get('T2', 1e+30), troe_params.get('T3', 1e-30)
        Fcent = (1 - alpha) * math.exp(-T / T3) + alpha * math.exp(-T / T1)
        Fcent += math.exp(-T2 / T) if T2 != 0.0 else 0.0
        c = -0.4 - 0.67 * math.log10(Fcent)
        n = 0.75 - 1.27 * math.log10(Fcent)
        d = 0.14
        F = 10.0 ** (math.log10(Fcent) / (1 + ((math.log10(Pr) + c) / (n - d * (math.log10(Pr) + c))) ** 2))
    return kinf * (Pr / (1 + Pr)) * F


def calculate_plog_rate_coefficient(reaction_data: dict, T: float, P: float, Ea_units: str) -> float:
    """
    Calculate the PLOG rate coefficient.

    Args:
        reaction_data (dict): Dictionary containing reaction parameters (the Cantera YAML dict of a single reaction).
        T (float): Temperature in Kelvin.
        P (float): Pressure in bar.
        Ea_units (str): Units of the activation energy.

    Returns:
        float: The rate coefficient at the specified temperature and pressure.
    """
    i, p_low, p_high = 0, 0.0, 0.0
    for i in range(len(reaction_data['rate-constants']) - 1):
        p_low = get_pressure_from_cantera(reaction_data['rate-constants'][i]['P'])
        p_high = get_pressure_from_cantera(reaction_data['rate-constants'][i + 1]['P'])
        if p_low <= P <= p_high:
            break
    k_low = calculate_arrhenius_rate_coefficient(A=reaction_data['rate-constants'][i]['A'],
                                                    n=reaction_data['rate-constants'][i]['b'],
                                                    Ea=reaction_data['rate-constants'][i]['Ea'],
                                                    T=T, Ea_units=Ea_units)
    k_high = calculate_arrhenius_rate_coefficient(A=reaction_data['rate-constants'][i + 1]['A'],
                                                    n=reaction_data['rate-constants'][i + 1]['b'],
                                                    Ea=reaction_data['rate-constants'][i + 1]['Ea'],
                                                    T=T, Ea_units=Ea_units)
    if P == p_low:
        return k_low
    if P == p_high:
        return k_high
    return k_low * 10 ** (math.log10(P / p_low) / math.log10(p_high / p_low) * math.log10(k_high / k_low))


def calculate_chebyshev_rate_coefficient(reaction_data: dict, T: float, P: float) -> float:
    """
    Calculate the Chebyshev rate coefficient.

    Args:
        reaction_data (dict): Dictionary containing reaction parameters (the Cantera YAML dict of a single reaction).
        T (float): Temperature in Kelvin.
        P (float): Pressure in bar.

    Returns:
        float: The rate coefficient at the specified temperature and pressure.
    """
    T_min, T_max = reaction_data['temperature-range'][0], reaction_data['temperature-range'][1]
    P_min = get_pressure_from_cantera(reaction_data['pressure-range'][0])
    P_max = get_pressure_from_cantera(reaction_data['pressure-range'][1])
    Tr = (2.0 / T - 1.0 / T_min - 1.0 / T_max) / (1.0 / T_max - 1.0 / T_min)
    Pr = (2.0 * math.log10(P) - math.log10(P_min) - math.log10(P_max)) / (math.log10(P_max) - math.log10(P_min))
    n_T = len(reaction_data['data'])
    n_P = len(reaction_data['data'][0])
    k = 0
    for t in range(n_T):
        for p in range(n_P):
            k += reaction_data['data'][t][p] * chebyshev(t, Tr) * chebyshev(p, Pr)
    return 10 ** k


def chebyshev(n: int, x: float) -> float:
    """
    Calculate the Chebyshev polynomial of the first kind.

    Args:
        n (int): Polynomial degree.
        x (float): Value at which to evaluate the polynomial.

    Returns:
        float: The interpolated Chebyshev polynomial value.
    """
    if n == 0:
        return 1
    if n == 1:
        return x
    t, t0, t1 = 0.0, 1, x
    for i in range(1, n):
        t = 2 * x * t1 - t0
        t0 = t1
        t1 = t
    return t


def get_pressure_from_cantera(p: str) -> float:
    """
    Get the pressure in bar from a Cantera pressure string.

    Args:
        p (str): Cantera pressure string.

    Returns:
        float: Pressure in bar.
    """
    p_units = p.split()[1]
    p_factor = P_UNIT_CONVERSION[p_units]
    return float(p.split()[0]) * p_factor


def get_rate_coefficient(reaction_data, T: float, P: Optional[float] = 1, Ea_units: str = 'J/mol') -> float:
    """
    Compute the rate coefficient based on the reaction type.
    Supports Arrhenius, ThreeBody, and Falloff types.

    Args:
        reaction_data (dict): Dictionary containing reaction parameters (the Cantera YAML dict of a single reaction).
        T (float): Temperature in Kelvin.
        P (float, optional): Pressure in bar (needed for some reaction types).
        Ea_units (str): Units of the activation energy.

    Returns:
        float: The rate coefficient at the specified temperature (and pressure if needed).
    """
    rate_info = reaction_data.get('rate-constant', {})
    if ('A' in rate_info and 'b' in rate_info and 'Ea' in rate_info and len(rate_info) == 3) \
            or ('type' in reaction_data and reaction_data['type'] == 'three-body'):
        return calculate_arrhenius_rate_coefficient(A=rate_info['A'], n=rate_info['b'], Ea=rate_info['Ea'], T=T, Ea_units=Ea_units)
    elif "type" in reaction_data and reaction_data['type'] == 'falloff':
        return calculate_troe_rate_coefficient(reaction_data, T, P, Ea_units)
    elif 'type' in reaction_data and reaction_data['type'] == 'pressure-dependent-Arrhenius':
        return calculate_plog_rate_coefficient(reaction_data, T, P, Ea_units)
    elif 'type' in reaction_data and reaction_data['type'] == 'Chebyshev':
        return calculate_chebyshev_rate_coefficient(reaction_data, T, P)
    else:
        raise ValueError("Unsupported reaction type or missing parameters.")


def get_Ea_units(path: str) -> str:
    """
    Get the units of the activation energy from a Cantera YAML file.

    Args:
        path (str): Path to the Cantera YAML file.

    Returns:
        str: Units of the activation energy.
    """
    content = read_yaml_file(path)
    units = content.get('units', {})
    Ea_units = units.get('activation-energy', 'kcal/mol')
    if Ea_units not in ['J/mol', 'kJ/mol', 'cal/mol', 'kcal/mol']:
        raise ValueError(f"Unsupported Ea units: {Ea_units} read from file {path}")
    return Ea_units


register_simulate_adapter("CanteraIDT", CanteraIDT)
