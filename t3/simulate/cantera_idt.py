"""
Cantera Simulator Adapter module for ignition delay time (IDT) calculations.

Models the post-reflected-shock state of a shock tube as a closed, adiabatic,
constant-volume batch reactor (``IdealGasReactor``). Also serves as the base
class for :class:`CanteraRCM`, which overrides the reactor type to
constant-pressure for rapid compression machine experiments,
preserving other IDT capabilities developed here.
"""

import concurrent.futures as cf
import logging
import math
import os
import traceback
from typing import Dict, List, Optional, Tuple

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

from arc.common import read_yaml_file, save_yaml_file
from arc.constants import R

from t3.common import determine_concentrations_by_equivalence_ratios, remove_numeric_parentheses
from t3.logger import Logger
from t3.simulate.adapter import SimulateAdapter
from t3.simulate.factory import register_simulate_adapter

py_logger = logging.getLogger(__name__)

DELTA_H = 0.1  # +/- 0.1 kJ/mol enthalpy perturbation for thermo brute-force SA (default, overridable via schema)
DELTA_K = 0.05  # *(1 + 5%) rate-coefficient perturbation for kinetics brute-force SA (default, overridable via schema)
EA_UNIT_CONVERSION = {'J/mol': 1, 'kJ/mol': 1e+3, 'cal/mol': 4.184, 'kcal/mol': 4.184e+3}
P_UNIT_CONVERSION = {'bar': 1, 'atm': 1.01325, 'Pa': 1e-5}


class CanteraIDT(SimulateAdapter):
    """
    CanteraIDT is a SimulateAdapter that runs ignition-delay-time (IDT) simulations
    mimicking a shock tube experiment. It uses Cantera's ``IdealGasReactor``
    (closed, adiabatic, constant-volume) to model the post-reflected-shock state.
    Optionally followed by a brute-force or adjoint sensitivity analysis
    that perturbs each species' enthalpy and each reaction's rate coefficient.

    Args:
        t3 (dict): The ``T3.t3`` block (input yaml or API).
        rmg (dict): The ``T3.rmg`` block (input yaml or API).
        paths (dict): The ``T3.paths`` dictionary.
        logger (Logger): A T3 ``Logger`` instance.
        atol (float): Absolute integration tolerance.
        rtol (float): Relative integration tolerance.
        observable_list (Optional[list]): Currently unused; the IDT itself is the observable.
        sa_atol (float): Absolute SA tolerance.
        sa_rtol (float): Relative SA tolerance.
        global_observables (Optional[List[str]]): Forwarded global observables (e.g. ``['IDT']``).
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
        self.atol = atol
        self.rtol = rtol
        self.sa_atol = sa_atol
        self.sa_rtol = sa_rtol
        self.observable_list = observable_list or list()
        self.global_observables = global_observables

        self.model: Optional[ct.Solution] = None
        self.cantera_simulation = None
        self.reactor_idt_dict: Optional[dict] = None
        self.inert_list = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'N2']
        self.inert_index_list: List[int] = list()
        self.spc_identifier_lookup: Dict[str, int] = dict()
        self.rxn_identifier_lookup: Dict[str, int] = dict()
        self.num_ct_reactions: Optional[int] = None
        self.num_ct_species: Optional[int] = None
        self.species_names_without_indices: List[str] = list()
        self.idt_sa_dict: dict = dict()

        sensitivity = self.t3.get('sensitivity', {})
        self.idt_criterion = sensitivity.get('idt_criterion', 'max_dOHdt')
        self.idt_sa_method = sensitivity.get('idt_sa_method', 'brute_force')
        self.delta_h = sensitivity.get('delta_h', DELTA_H)
        self.delta_k = sensitivity.get('delta_k', DELTA_K)
        self.adaptive_perturbation = sensitivity.get('adaptive_perturbation', False)
        self.experimental_idt_path = sensitivity.get('experimental_idt_path', None)

        self.set_up()
        self.radical_label = self.determine_radical_label()

    def set_up(self):
        """
        Read the Cantera input file and populate the model-derived attributes.
        """
        self.model = ct.Solution(infile=self.paths['cantera annotated'])
        self.num_ct_reactions = len(self.model.reactions())
        self.num_ct_species = len(self.model.species())
        self.inert_index_list = [i for i, s in enumerate(self.model.species()) if s.name in self.inert_list]
        self.spc_identifier_lookup = {spc.name: i for i, spc in enumerate(self.model.species())}
        self.rxn_identifier_lookup = {rxn.equation: i for i, rxn in enumerate(self.model.reactions())}
        self.species_names_without_indices = [remove_numeric_parentheses(self.model.species()[i].name)
                                              for i in range(self.num_ct_species)]

    def simulate(self):
        """
        Simulate the mechanism and compute ignition delay times for every working point.
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
        Simulate every (reactor, T, P, φ) working point and return a nested IDT dictionary.

        The expansion of T / P / φ depends on the reactor's ``idt_mode``:

        - ``'matrix'`` (default): full Cartesian product T_list × P_list × φ_list.
        - ``'row'``: index-aligned tuples (T[i], P[i], φ[i]). All three lists must
          share the same length.

        Args:
            save_yaml (bool): Whether to dump the resulting dict to ``figs/idt_dict.yaml``.
            save_fig (bool): Whether to draw IDT(1000/T) plots.
            energy (str): Cantera reactor energy mode (``'on' / 'off' / 'off after ignition'``).
            max_idt (float): Hard upper bound on simulation time, seconds.

        Returns:
            dict: ``reactor_idt_dict[r][φ][P][T] = idt_seconds`` (or ``None`` on failure).
        """
        infile = self.paths['cantera annotated']
        equivalence_ratios, concentration_combinations = self.get_concentration_combinations()
        reactor_idt_dict: dict = dict()
        for r, reactor in enumerate(self.rmg['reactors']):
            T_list, P_list = get_t_and_p_lists(reactor)
            mode = reactor.get('idt_mode', 'matrix')
            combinations = self._build_combinations(
                r=r,
                T_list=T_list,
                P_list=P_list,
                equivalence_ratios=equivalence_ratios,
                concentration_combinations=concentration_combinations,
                infile=infile,
                save_fig=save_fig,
                energy=energy,
                max_idt=max_idt,
                mode=mode,
            )
            results = [self.simulate_idt_for_a_point(*args) for args in combinations]
            idt_dict: dict = dict()
            for (_, T, P, _, phi, *_), idt in zip(combinations, results):
                key_phi = phi if phi is not None else 0
                idt_dict.setdefault(key_phi, {}).setdefault(P, {})[T] = idt
            if len(T_list) >= 3 and save_fig:
                plot_idt_vs_temperature(idt_dict, figs_path=self.paths['figs'], reactor_index=r)
            reactor_idt_dict[r] = idt_dict
        if save_yaml:
            os.makedirs(self.paths['figs'], exist_ok=True)
            save_yaml_file(os.path.join(self.paths['figs'], 'idt_dict.yaml'), reactor_idt_dict)
        return reactor_idt_dict

    def _build_combinations(self,
                            r: int,
                            T_list: List[float],
                            P_list: List[float],
                            equivalence_ratios: Optional[List[float]],
                            concentration_combinations: Optional[List[dict]],
                            infile: str,
                            save_fig: bool,
                            energy: str,
                            max_idt: float,
                            mode: str,
                            ) -> List[tuple]:
        """
        Expand T / P / φ into the list of working points to simulate.
        """
        if equivalence_ratios is None or concentration_combinations is None:
            x = {spc['label']: spc['concentration'] for spc in self.rmg['species'] if spc['concentration']}
            return [(r, T, P, x, None, infile, save_fig, energy, max_idt)
                    for T in T_list for P in P_list]

        if mode == 'row':
            n = len(equivalence_ratios)
            if not (len(T_list) == len(P_list) == n):
                raise ValueError(
                    f"Reactor {r} uses idt_mode='row' but T_list (len {len(T_list)}), "
                    f"P_list (len {len(P_list)}) and equivalence_ratios (len {n}) "
                    f"do not all share the same length.")
            return [(r, T_list[i], P_list[i], concentration_combinations[i],
                     equivalence_ratios[i], infile, save_fig, energy, max_idt)
                    for i in range(n)]

        # default: full cartesian
        return [(r, T, P, X, equivalence_ratios[idx], infile, save_fig, energy, max_idt)
                for T in T_list
                for P in P_list
                for idx, X in enumerate(concentration_combinations)]

    def simulate_idt_for_a_point(self,
                                 r: int,
                                 t: float,
                                 p: float,
                                 x: dict,
                                 phi: Optional[float],
                                 infile: str,
                                 save_fig: bool = True,
                                 energy: str = 'on',
                                 max_idt: float = 1.0,
                                 ) -> Optional[float]:
        """
        Simulate one IdealGasReactor at fixed (T, P, X) and return the IDT in seconds.
        """
        fig_name = (f'R{r}_{phi}_{p:.2f}_bar_{t:.2f}_K.png'
                    if phi is not None else f'R{r}_{p:.2f}_bar_{t:.2f}_K.png')
        model = ct.Solution(infile=infile)
        model.TPX = t, p * 1e5, x
        reactor = self._create_reactor(model, energy)
        net = ct.ReactorNet([reactor])
        net.atol, net.rtol = self.atol, self.rtol
        net.atol_sensitivity, net.rtol_sensitivity = self.sa_atol, self.sa_rtol
        time_history = ct.SolutionArray(model, extra='t')
        t_, counter = 0.0, 0
        use_radical_early_exit = self.idt_criterion != 'max_dTdt' and self.radical_label is not None
        while t_ < max_idt:
            t_ = net.step()
            time_history.append(reactor.thermo.state, t=t_)
            if use_radical_early_exit and counter % 100 == 0:
                concentrations = np.asarray([row[0] for row in time_history(self.radical_label).X], dtype=np.float32)
                max_c_idx = int(np.argmax(concentrations))
                if (concentrations[max_c_idx] > concentrations[-1] * 1.2
                        and len(concentrations) > max_c_idx * 1.1
                        and time_history.t[-1] > time_history.t[max_c_idx] * 1.2):
                    break
            counter += 1
        radical = self.radical_label
        if self.idt_criterion == 'max_radical_dt':
            radical = self._find_best_radical(time_history) or radical
        return compute_idt(time_history=time_history,
                           radical_label=radical,
                           criterion=self.idt_criterion,
                           figs_path=self.paths['figs'],
                           fig_name=fig_name if save_fig else None,
                           )

    def _create_reactor(self, model: ct.Solution, energy: str = 'on'):
        """Create the Cantera reactor for IDT simulation. Subclasses override this."""
        return ct.IdealGasReactor(model, energy=energy)

    def determine_radical_label(self) -> Optional[str]:
        """
        Pick a radical species label used for IDT detection at init time.

        - ``max_dOHdt`` (default): prefer OH, fall back to H if OH is
          absent from the mechanism.
        - ``max_dTdt``: not radical-based, returns None.
        - ``max_radical_dt``: returns the same OH/H fallback as a safety
          net; the actual per-simulation selection of the
          highest-peak-concentration radical (among non-inert species
          with ≤3 heavy atoms) happens later in :meth:`_find_best_radical`.
        """
        if self.idt_criterion == 'max_dTdt':
            return None
        h, oh = None, None
        for i, species in enumerate(self.species_names_without_indices):
            if species.lower() == 'oh':
                oh = self.model.species()[i].name
            if species.lower() == 'h':
                h = self.model.species()[i].name
            if oh is not None:
                break
        if self.idt_criterion == 'max_dOHdt':
            return oh or h
        # max_radical_dt: defer to simulation-time selection, but set a default fallback
        return oh or h

    def _find_best_radical(self, time_history: ct.SolutionArray) -> Optional[str]:
        """
        Scan all species in the simulation history and return the label of the
        radical with the highest peak concentration. Candidates are non-inert
        species with ≤3 heavy atoms (C, N, O, S).
        """
        best_label, best_peak = None, 0.0
        for i, spc in enumerate(self.model.species()):
            if i in self.inert_index_list:
                continue
            heavy = sum(v for el, v in spc.composition.items() if el not in ('H', 'He', 'Ar', 'Ne', 'Kr', 'Xe'))
            if heavy > 3:
                continue
            conc = np.asarray([row[0] for row in time_history(spc.name).X], dtype=np.float64)
            peak = float(np.max(conc))
            if peak > best_peak:
                best_peak = peak
                best_label = spc.name
        if best_label is not None:
            py_logger.info(f'Auto-selected IDT radical: {best_label} (peak conc: {best_peak:.2e})')
        return best_label

    def _prescan_radical(self,
                         infile: str,
                         T: float,
                         P: float,
                         X: dict,
                         energy: str = 'on',
                         max_idt: float = 1.0,
                         ) -> Optional[str]:
        """
        Run a quick simulation (no SA) at a single condition to identify the
        best radical for ``max_radical_dt``. Returns the radical label, or
        ``None`` if no suitable radical is found.
        """
        model = ct.Solution(infile=infile)
        model.TPX = T, P * 1e5, X
        reactor = self._create_reactor(model, energy)
        net = ct.ReactorNet([reactor])
        net.atol, net.rtol = self.atol, self.rtol
        time_history = ct.SolutionArray(model, extra='t')
        t_ = 0.0
        while t_ < max_idt:
            t_ = net.step()
            time_history.append(reactor.thermo.state, t=t_)
        if len(time_history.t) < 2:
            return self.radical_label
        radical = self._find_best_radical(time_history)
        return radical or self.radical_label

    def get_cantera_species_label(self, rmg_label: str) -> Optional[str]:
        """
        Look up the Cantera species name (which usually carries an ``(N)`` suffix) for an RMG label.
        """
        for i, label in enumerate(self.species_names_without_indices):
            if label == rmg_label:
                return self.model.species()[i].name
            if label == remove_numeric_parentheses(rmg_label):
                return self.model.species()[i].name
        return None

    def get_concentration_combinations(self) -> Tuple[Optional[List[float]], Optional[List[dict]]]:
        """
        Build the per-φ concentration dictionaries that drive each IDT simulation.

        Uses the new ``determine_concentrations_by_equivalence_ratios`` API which returns::

            {'equivalence_ratios': [φ1, φ2, ...],
             'concentrations': {<species_label>: [c_at_φ1, c_at_φ2, ...], ...}}

        Returns:
            (equivalence_ratios, concentration_combinations) — both ``None`` if no fuel
            with ``equivalence_ratios`` is declared. Each entry of
            ``concentration_combinations`` is a dict keyed by **Cantera** species name
            (with the ``(N)`` suffix) suitable for assigning to ``model.TPX``.
        """
        objects = determine_concentrations_by_equivalence_ratios(species=self.rmg['species'])
        if objects is None:
            return None, None
        equivalence_ratios = objects['equivalence_ratios']
        per_species_columns = objects['concentrations']
        concentration_combinations: List[dict] = list()
        for i in range(len(equivalence_ratios)):
            concentration_dict: Dict[str, float] = dict()
            for label, column in per_species_columns.items():
                cantera_label = self.get_cantera_species_label(label)
                if cantera_label is not None and column[i]:
                    concentration_dict[cantera_label] = column[i]
            # Inerts that are NOT in the role-driven mixture (e.g. trace bath gas) still
            # need to be included verbatim.
            role_labels = set(per_species_columns.keys())
            for spc in self.rmg['species']:
                if spc['label'] in role_labels:
                    continue
                if spc.get('role') is not None:
                    continue
                if spc.get('concentration'):
                    cantera_label = self.get_cantera_species_label(spc['label'])
                    if cantera_label is not None:
                        concentration_dict[cantera_label] = spc['concentration']
            concentration_combinations.append(concentration_dict)
        return equivalence_ratios, concentration_combinations

    def get_sa_coefficients(self) -> Optional[dict]:
        """
        Run the IDT sensitivity analysis. Dispatches to adjoint or brute-force based on
        ``self.idt_sa_method``. Tuning knobs (``top_SA_species``, ``top_SA_reactions``,
        ``max_sa_workers``, and whether to dump YAML) are read from
        ``self.t3['sensitivity']`` — no kwargs.

        The full SA result is dumped to ``paths['SA IDT dict']`` and the pruned top-X result
        to ``paths['SA IDT dict top X']`` whenever ``save_yaml`` (from the sensitivity block,
        defaults to True) is truthy.
        """
        if self.idt_sa_method == 'adjoint':
            return self._get_sa_coefficients_adjoint()
        return self._get_sa_coefficients_brute_force()

    def _sa_opts(self) -> Tuple[int, int, int, bool]:
        """Read the SA tuning knobs from ``self.t3['sensitivity']``."""
        sens = self.t3.get('sensitivity', {}) or {}
        return (int(sens.get('top_SA_species', 10)),
                int(sens.get('top_SA_reactions', 10)),
                int(sens.get('max_sa_workers', 24)),
                bool(sens.get('save_sa_yaml', True)))

    def _get_sa_coefficients_brute_force(self) -> Optional[dict]:
        """
        Brute-force IDT sensitivity analysis: for every species enthalpy and every
        reaction rate coefficient, perturb a copy of the mechanism and re-run all working
        points in a worker process. Then collapse to ``dlnIDT / dX`` coefficients and keep
        only the top ``top_SA_*`` of each.
        """
        top_SA_species, top_SA_reactions, max_workers, save_yaml = self._sa_opts()
        if self.logger is not None:
            self.logger.info(f'Running brute-force IDT SA using {max_workers} workers...')
        if not self.reactor_idt_dict:
            self.simulate()

        # Adaptive perturbation: scale delta_h relative to |H298|
        delta_h = self.delta_h
        delta_k = self.delta_k
        if self.adaptive_perturbation:
            py_logger.info('Adaptive perturbation sizing enabled')

        sa_dict: Dict[str, Dict[str, dict]] = {'thermo': {'IDT': dict()}, 'kinetics': {'IDT': dict()}}
        tasks = ([('thermo', i) for i in range(self.num_ct_species)]
                 + [('kinetics', i) for i in range(self.num_ct_reactions)])
        os.makedirs(self.paths['SA'], exist_ok=True)

        total_tasks = len(tasks)
        succeeded, failed, retried = 0, 0, 0

        with cf.ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Compute per-species adaptive delta_h if enabled
            adaptive_delta_h_map: Dict[int, float] = dict()
            if self.adaptive_perturbation:
                for i in range(self.num_ct_species):
                    h298 = get_h298(self.model, i)
                    adaptive_delta_h_map[i] = max(delta_h, 0.001 * abs(h298))

            def _submit(task_, dh_, dk_):
                return executor.submit(worker,
                                       task_,
                                       self.paths['cantera annotated'],
                                       self.paths['SA'],
                                       self.t3,
                                       self.paths,
                                       self.rmg,
                                       self.logger,
                                       1.0,
                                       dh_,
                                       dk_,
                                       )

            future_to_task = dict()
            for task in tasks:
                kind, index = task
                task_delta_h = adaptive_delta_h_map.get(index, delta_h) if kind == 'thermo' else delta_h
                task_delta_k = delta_k
                future_to_task[_submit(task, task_delta_h, task_delta_k)] = (task, task_delta_h, task_delta_k)

            # First pass: collect results, defer retries so the pool stays saturated.
            to_retry: List[tuple] = list()
            for future in cf.as_completed(future_to_task):
                task, task_delta_h, task_delta_k = future_to_task[future]
                try:
                    sa_dict[task[0]]['IDT'][task[1]] = future.result()
                    succeeded += 1
                except Exception as e:
                    py_logger.warning(f"Task {task} failed: {e}. Scheduling retry with halved perturbation.")
                    to_retry.append((task, task_delta_h, task_delta_k))

            # Second pass: re-submit all failures at once, then harvest as they complete.
            retry_future_to_task: Dict = dict()
            for task, task_delta_h, task_delta_k in to_retry:
                kind, _ = task
                retry_delta_h = task_delta_h / 2 if kind == 'thermo' else task_delta_h
                retry_delta_k = task_delta_k / 2 if kind == 'kinetics' else task_delta_k
                retry_future_to_task[_submit(task, retry_delta_h, retry_delta_k)] = task
                retried += 1

            for future in cf.as_completed(retry_future_to_task):
                task = retry_future_to_task[future]
                try:
                    sa_dict[task[0]]['IDT'][task[1]] = future.result()
                    succeeded += 1
                except Exception as e2:
                    failed += 1
                    if self.logger is not None:
                        self.logger.error(f"Task {task} permanently failed after retry:\n{e2}\n"
                                          f"{traceback.format_exc()}")

        py_logger.info(f"SA completed: {succeeded}/{total_tasks} succeeded "
                       f"({failed} failed, {retried} retried)")
        if total_tasks > 0 and failed / total_tasks > 0.1:
            py_logger.warning(f"SA results may be unreliable (>10% failures: "
                              f"{failed}/{total_tasks} = {100 * failed / total_tasks:.1f}%)")

        idt_sa_dict_all = compute_idt_sa(reactor_idt_dict=self.reactor_idt_dict,
                                         perturbed_idt_dict=sa_dict,
                                         delta_h=delta_h,
                                         delta_k=delta_k,
                                         )
        self.idt_sa_dict = get_top_sa_coefficients(idt_sa_dict=idt_sa_dict_all,
                                                   top_species=top_SA_species,
                                                   top_reactions=top_SA_reactions)
        if save_yaml:
            save_yaml_file(path=self.paths['SA IDT dict'], content=idt_sa_dict_all)
            save_yaml_file(path=self.paths['SA IDT dict top X'], content=self.idt_sa_dict)
        return self.idt_sa_dict

    def _simulate_idt_adjoint_sa(self,
                                 r: int,
                                 t: float,
                                 p: float,
                                 x: dict,
                                 phi: Optional[float],
                                 infile: str,
                                 energy: str = 'on',
                                 max_idt: float = 1.0,
                                 sa_observable: Optional[str] = None,
                                 ) -> Tuple[Optional[float], dict, dict]:
        """
        Run a single IDT simulation with Cantera's built-in adjoint sensitivity analysis.

        Args:
            sa_observable: The species or ``'temperature'`` whose sensitivity is
                extracted at each time step. If ``None``, determined automatically
                from ``self.idt_criterion`` and ``self.radical_label``.

        Returns:
            (idt, kinetics_sa, thermo_sa) where kinetics_sa and thermo_sa are dicts
            mapping parameter index to SA coefficient at t=IDT.
        """
        model = ct.Solution(infile=infile)
        model.TPX = t, p * 1e5, x
        reactor = self._create_reactor(model, energy)
        net = ct.ReactorNet([reactor])
        net.atol, net.rtol = self.atol, self.rtol
        net.atol_sensitivity, net.rtol_sensitivity = self.sa_atol, self.sa_rtol

        n_reactions = len(model.reactions())
        n_species = len(model.species())

        # Register sensitivity parameters: all reactions and all species enthalpies
        for i in range(n_reactions):
            reactor.add_sensitivity_reaction(i)
        for i in range(n_species):
            reactor.add_sensitivity_species_enthalpy(i)

        time_history = ct.SolutionArray(model, extra='t')
        # Store sensitivity matrix at each time step
        sa_history: List[np.ndarray] = list()

        # Determine the observable for adjoint SA extraction.
        if sa_observable is not None:
            obs_name = sa_observable
        elif self.idt_criterion == 'max_dTdt':
            obs_name = 'temperature'
        else:
            obs_name = self.radical_label or 'temperature'

        time_ = 0.0
        while time_ < max_idt:
            time_ = net.step()
            time_history.append(reactor.thermo.state, t=time_)
            # Extract SA for the observable at this time step
            n_params = n_reactions + n_species
            sa_row = np.zeros(n_params)
            for j in range(n_params):
                try:
                    sa_row[j] = net.sensitivity(obs_name, j)
                except (ct.CanteraError, IndexError):
                    sa_row[j] = 0.0
            sa_history.append(sa_row)

        if len(time_history.t) < 2:
            return None, dict(), dict()

        # Determine IDT index
        times = time_history.t
        if self.idt_criterion == 'max_dTdt':
            T_data = time_history.T
            dt = np.diff(times)
            dt = np.where(dt == 0, np.finfo(float).tiny, dt)
            dTdt = np.diff(T_data) / dt
            idt_idx = int(np.argmax(dTdt))
            if idt_idx > len(times) - 10 or times[idt_idx] < 1e-12:
                return None, dict(), dict()
            if T_data[idt_idx] - T_data[0] < 50:
                return None, dict(), dict()
        else:
            # For max_radical_dt, find best radical first
            radical = self.radical_label
            if self.idt_criterion == 'max_radical_dt':
                radical = self._find_best_radical(time_history) or radical
            if radical is None:
                return None, dict(), dict()
            conc = np.asarray([row[0] for row in time_history(radical).X], dtype=np.float64)
            if all(c == 0 for c in conc):
                return None, dict(), dict()
            dc_dt = np.diff(conc) / np.diff(times)
            idt_idx = int(np.argmax(dc_dt))
            if idt_idx > len(times) - 10 or times[idt_idx] < 1e-12:
                return None, dict(), dict()

        idt = float(times[idt_idx])
        sa_at_idt = sa_history[idt_idx]

        kinetics_sa: Dict[int, float] = dict()
        thermo_sa: Dict[int, float] = dict()
        for j in range(n_reactions):
            if abs(sa_at_idt[j]) > 1e-20:
                kinetics_sa[j] = float(sa_at_idt[j])
        for j in range(n_species):
            if abs(sa_at_idt[n_reactions + j]) > 1e-20:
                thermo_sa[j] = float(sa_at_idt[n_reactions + j])

        return idt, kinetics_sa, thermo_sa

    def _get_sa_coefficients_adjoint(self) -> Optional[dict]:
        """
        Run adjoint (Cantera built-in) IDT sensitivity analysis. A single simulation per
        condition extracts dln[observable]/dln[k_i] at t=IDT.
        """
        top_SA_species, top_SA_reactions, _, save_yaml = self._sa_opts()
        if self.logger is not None:
            self.logger.info('Running adjoint IDT SA...')
        if not self.reactor_idt_dict:
            self.simulate()

        infile = self.paths['cantera annotated']
        equivalence_ratios, concentration_combinations = self.get_concentration_combinations()

        idt_sa_dict: dict = {'thermo': {'IDT': dict()}, 'kinetics': {'IDT': dict()}}

        for r, reactor in enumerate(self.rmg['reactors']):
            T_list, P_list = get_t_and_p_lists(reactor)
            mode = reactor.get('idt_mode', 'matrix')

            idt_sa_dict['thermo']['IDT'][r] = dict()
            idt_sa_dict['kinetics']['IDT'][r] = dict()

            # Build working-point list (same as simulation)
            if equivalence_ratios is None or concentration_combinations is None:
                x = {spc['label']: spc['concentration'] for spc in self.rmg['species'] if spc['concentration']}
                points = [(T, P, x, None) for T in T_list for P in P_list]
            elif mode == 'row':
                points = [(T_list[i], P_list[i], concentration_combinations[i],
                           equivalence_ratios[i]) for i in range(len(equivalence_ratios))]
            else:
                points = [(T, P, X, equivalence_ratios[idx])
                          for T in T_list for P in P_list
                          for idx, X in enumerate(concentration_combinations)]

            # For max_radical_dt, run a quick prescan simulation at the median
            # condition to identify the best radical before the SA loop.
            sa_observable = None
            if self.idt_criterion == 'max_radical_dt' and points:
                mid = len(points) // 2
                T_mid, P_mid, X_mid, _ = points[mid]
                sa_observable = self._prescan_radical(infile=infile, T=T_mid, P=P_mid, X=X_mid)
                py_logger.info(f'Adjoint SA: prescan selected radical "{sa_observable}" '
                               f'at T={T_mid} K, P={P_mid} bar')

            for T, P, X, phi in points:
                key_phi = phi if phi is not None else 0
                idt, kin_sa, thermo_sa = self._simulate_idt_adjoint_sa(
                    r=r, t=T, p=P, x=X, phi=phi, infile=infile,
                    sa_observable=sa_observable,
                )
                if idt is None:
                    continue

                for token, sa_data in [('kinetics', kin_sa), ('thermo', thermo_sa)]:
                    sa_node = idt_sa_dict[token]['IDT'][r]
                    sa_node.setdefault(key_phi, {}).setdefault(P, {}).setdefault(T, {})
                    sa_node[key_phi][P][T].update(sa_data)

        self.idt_sa_dict = get_top_sa_coefficients(idt_sa_dict=idt_sa_dict,
                                                   top_species=top_SA_species,
                                                   top_reactions=top_SA_reactions)
        if save_yaml:
            os.makedirs(self.paths['SA'], exist_ok=True)
            save_yaml_file(path=self.paths['SA IDT dict'], content=idt_sa_dict)
            save_yaml_file(path=self.paths['SA IDT dict top X'], content=self.idt_sa_dict)
        return self.idt_sa_dict

    def compare_with_experiment(self, exp_data_path: str) -> dict:
        """
        Load experimental IDT data from a YAML file and compare with simulated values.

        The experimental YAML format::

            citation: "Author et al., Journal (Year)"
            data:
              - {T: 1000, P: 10, phi: 1.0, idt: 2.3e-3}
              - ...

        T in K, P in bar, idt in seconds.

        Returns:
            dict with keys ``'citation'``, ``'points'`` (list of per-point comparisons),
            ``'rmse_log'`` (RMSE of log10 errors), ``'n_points'``, ``'n_matched'``.
        """
        if not self.reactor_idt_dict:
            self.simulate()
        exp = read_yaml_file(exp_data_path)
        citation = exp.get('citation', 'unknown')
        points = exp.get('data', [])

        comparisons: list = list()
        log_errors: list = list()

        for pt in points:
            exp_T, exp_P = float(pt['T']), float(pt['P'])
            exp_phi = float(pt.get('phi', 0))
            exp_idt = float(pt['idt'])

            # Find nearest simulated point
            best_sim_idt = None
            best_dist = float('inf')
            for reactor_data in self.reactor_idt_dict.values():
                for phi, phi_data in reactor_data.items():
                    for P, p_data in phi_data.items():
                        for T, idt in p_data.items():
                            if idt is None:
                                continue
                            dist = ((T - exp_T) / exp_T) ** 2 + ((P - exp_P) / exp_P) ** 2
                            if exp_phi > 0:
                                dist += ((phi - exp_phi) / exp_phi) ** 2
                            if dist < best_dist:
                                best_dist = dist
                                best_sim_idt = idt

            if best_sim_idt is not None and exp_idt > 0 and best_sim_idt > 0:
                log_err = math.log10(best_sim_idt / exp_idt)
                log_errors.append(log_err)
                comparisons.append({
                    'T': exp_T, 'P': exp_P, 'phi': exp_phi,
                    'idt_exp': exp_idt, 'idt_sim': best_sim_idt,
                    'log10_error': log_err,
                })
            else:
                comparisons.append({
                    'T': exp_T, 'P': exp_P, 'phi': exp_phi,
                    'idt_exp': exp_idt, 'idt_sim': None,
                    'log10_error': None,
                })

        rmse_log = float(np.sqrt(np.mean(np.array(log_errors) ** 2))) if log_errors else None
        n_matched = len(log_errors)

        result = {
            'citation': citation,
            'n_points': len(points),
            'n_matched': n_matched,
            'rmse_log': rmse_log,
            'points': comparisons,
        }

        if rmse_log is not None:
            py_logger.info(f"Experimental comparison ({citation}): "
                           f"{n_matched}/{len(points)} points matched, "
                           f"RMSE(log10) = {rmse_log:.3f}")
        else:
            py_logger.info(f"Experimental comparison ({citation}): no matched points")

        return result

    def get_idt_by_T(self) -> dict:
        """
        Return the IDT(T) curve already computed by ``simulate()``.
        Format: ``{'idt': [...], 'idt_index': [...]}`` with one entry per condition.
        """
        if not self.reactor_idt_dict:
            return {'idt': list(), 'idt_index': list()}
        idts: List[float] = list()
        for reactor_data in self.reactor_idt_dict.values():
            for phi_data in reactor_data.values():
                for p_data in phi_data.values():
                    for _, idt in p_data.items():
                        if idt is not None:
                            idts.append(idt)
        return {'idt': idts, 'idt_index': list(range(len(idts)))}


def worker(task: tuple,
           model_path: str,
           work_dir: str,
           t3: dict,
           paths: dict,
           rmg: dict,
           logger: Logger,
           max_idt: float = 1.0,
           delta_h: float = DELTA_H,
           delta_k: float = DELTA_K,
           ) -> dict:
    """
    Worker function (must be picklable for ProcessPoolExecutor) that perturbs a single
    parameter, builds a temporary CanteraIDT adapter against the perturbed model, and
    returns the IDT dictionary for that perturbation.
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
                                       delta_h=delta_h,
                                       )
        elif kind == 'kinetics':
            success = perturb_reaction_rate_coefficient(original_path=model_path,
                                                        perturbed_path=perturbed_model_path,
                                                        reaction_index=index,
                                                        logger=logger,
                                                        delta_k=delta_k,
                                                        )
        else:
            raise ValueError(f'Unknown task type {kind}')
        if success:
            new_paths = {**paths, 'cantera annotated': perturbed_model_path}
            ct_idt_adapter = CanteraIDT(t3=t3, paths=new_paths, rmg=rmg, logger=logger)
            return ct_idt_adapter.simulate_idt_for_all_reactors(save_yaml=False,
                                                                save_fig=False,
                                                                energy='on',
                                                                max_idt=max_idt)
        return dict()
    finally:
        if os.path.exists(perturbed_model_path):
            os.remove(perturbed_model_path)


def compute_idt(time_history: ct.SolutionArray,
                radical_label: Optional[str] = None,
                criterion: str = 'max_dOHdt',
                figs_path: Optional[str] = None,
                fig_name: Optional[str] = None,
                ) -> Optional[float]:
    """
    Find the IDT using the specified criterion. Returns ``None`` if the trajectory
    looks degenerate (no real ignition).

    Criteria:
        - ``max_dOHdt``: time of maximum d[OH]/dt; falls back to d[H]/dt
          if OH is absent from the mechanism.
        - ``max_radical_dt``: auto-detect the highest-concentration small
          radical, use max d[radical]/dt.
        - ``max_dTdt``: time of maximum dT/dt.
    """
    if figs_path is not None:
        figs_path = os.path.join(figs_path, 'IDTs')
        os.makedirs(figs_path, exist_ok=True)
    times = time_history.t

    if criterion == 'max_dTdt':
        T_data = time_history.T
        if len(T_data) < 2:
            return None
        dt = np.diff(times)
        dt = np.where(dt == 0, np.finfo(float).tiny, dt)
        dTdt = np.diff(T_data) / dt
        idt_index = int(np.argmax(dTdt))
        idt = float(times[idt_index])
        if idt_index > len(times) - 10 or idt < 1e-12:
            return None
        T_rise = T_data[idt_index] - T_data[0]
        if T_rise < 50:  # less than 50K rise → no real ignition
            return None
        if figs_path is not None and fig_name is not None:
            fig, ax = plt.subplots()
            try:
                ax.plot(times, T_data)
                ax.plot(times[idt_index], T_data[idt_index], 'o')
                ax.set_xlabel('Time (s)')
                ax.set_ylabel('Temperature (K)')
                ax.set_title(f'IDT = {idt:.2e} s (max dT/dt)')
                ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                ax.xaxis.get_major_formatter().set_scientific(True)
                fig.savefig(os.path.join(figs_path, fig_name))
            except (AttributeError, ValueError) as e:
                py_logger.debug(f'IDT plot failed: {e}')
            finally:
                plt.close(fig)
        return idt

    # max_dOHdt or max_radical_dt: use radical concentration
    if radical_label is None:
        return None
    concentration = np.asarray([row[0] for row in time_history(radical_label).X], dtype=np.float32)
    if all(c == 0 for c in concentration):
        return None
    dt = np.diff(times)
    dt = np.where(dt == 0, np.finfo(float).tiny, dt)
    dc_dt = np.diff(concentration) / dt
    idt_index_dc_dt = int(np.argmax(dc_dt))
    idt_index_c = int(np.argmax(concentration))
    idt = float(times[idt_index_dc_dt])
    if (idt_index_dc_dt > len(times) - 10
            or idt < 1e-12
            or max(concentration) < concentration[0] * 100):
        return None
    if figs_path is not None and fig_name is not None:
        fig, ax = plt.subplots()
        try:
            ax.plot(times, concentration)
            ax.plot(times[idt_index_dc_dt], concentration[idt_index_dc_dt], 'o')
            ax.set_xlabel('Time (s)')
            ax.set_ylabel(f'[{radical_label}]')
            ax.set_title(f'IDT = {idt:.2e} s')
            ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.xaxis.get_major_formatter().set_scientific(True)
            if times[idt_index_c] > times[idt_index_dc_dt] * 2.5:
                x_min = min(times[idt_index_dc_dt] * 0.8, times[idt_index_c] * 0.1)
                x_max = max(times[idt_index_dc_dt] * 1.2, times[idt_index_c] * 1.5)
                if times[idt_index_dc_dt] < x_max * 0.05:
                    x_min = 0
            else:
                x_min, x_max = times[idt_index_dc_dt] * 0.8, times[idt_index_dc_dt] * 1.2
            ax.set_xlim(x_min, x_max)
            fig.savefig(os.path.join(figs_path, fig_name))
        except (AttributeError, ValueError) as e:
            py_logger.debug(f'IDT plot failed: {e}')
        finally:
            plt.close(fig)
    return idt


def get_t_and_p_lists(reactor: dict,
                      num_t_points: int = 25,
                      ) -> Tuple[List[float], List[float]]:
    """
    Expand the reactor's ``T`` and ``P`` fields into explicit numeric lists.

    A scalar stays a 1-element list. A 2-element list is treated as a min/max range and
    expanded (T linearly in 1/T, P logarithmically). A longer list is taken verbatim,
    which is what ``idt_mode='row'`` relies on.
    """
    if isinstance(reactor['T'], (int, float)):
        T_list = [float(reactor['T'])]
    elif len(reactor['T']) == 2:
        n = min(int(abs(reactor['T'][1] - reactor['T'][0]) / 10), num_t_points)
        n = max(n, 2)
        inverse_ts = np.linspace(1 / reactor['T'][1], 1 / reactor['T'][0], num=n)
        T_list = [1 / inverse_t for inverse_t in inverse_ts[::-1]]
    else:
        T_list = [float(t) for t in reactor['T']]
    if isinstance(reactor['P'], (int, float)):
        P_list = [float(reactor['P'])]
    elif len(reactor['P']) == 2:
        log_p = np.linspace(math.log10(reactor['P'][0]), math.log10(reactor['P'][1]), num=3)
        P_list = [10 ** p for p in log_p]
    else:
        P_list = [float(p) for p in reactor['P']]
    return [float(t) for t in T_list], [float(p) for p in P_list]


def plot_idt_vs_temperature(idt_dict: dict,
                            figs_path: str,
                            reactor_index: int = 0,
                            exp_data: Optional[dict] = None,
                            ) -> None:
    """
    Plot ``IDT(1000/T)`` per φ and pressure condition. If ``exp_data`` is provided,
    overlay it and only plot the φ/P combinations that appear in it.
    ``idt_dict`` is the nested dict ``{φ: {P: {T: idt}}}`` produced by ``simulate_idt_for_all_reactors``.
    """
    figs_path = os.path.join(figs_path, 'IDT_vs_T')
    os.makedirs(figs_path, exist_ok=True)
    for phi, phi_data in idt_dict.items():
        if exp_data is not None and phi not in exp_data:
            continue
        for p, phi_p_data in phi_data.items():
            if exp_data is not None and p not in exp_data[phi]:
                continue
            fig_name = f'R{reactor_index}_{phi}_{p:.2f}_bar.png'
            fig, ax = plt.subplots()
            try:
                sim_data = {t: v for t, v in phi_p_data.items() if v is not None}
                ax.set_xlabel('1000/T (1/K)')
                ax.set_ylabel('IDT (s)')
                ax.set_title(f'IDT vs. 1000/T, phi = {phi}, P = {p:.2f} bar')
                ax.scatter([1000 / t for t in sim_data.keys()], list(sim_data.values()),
                           label='simulation', color='blue', marker='o', linestyle='-')
                ax.set_yscale('log')
                if exp_data is not None:
                    ax.scatter([1000 / t for t in exp_data[phi][p].keys()],
                               list(exp_data[phi][p].values()),
                               label='experiment', color='orange', marker='D')
                ax.legend(loc='lower right')
                fig.savefig(os.path.join(figs_path, fig_name))
            except (AttributeError, ValueError) as e:
                py_logger.debug(f'IDT vs T plot failed: {e}')
            finally:
                plt.close(fig)


def perturb_enthalpy(original_path: str,
                     perturbed_path: str,
                     species_index: int,
                     logger: Optional[Logger] = None,
                     delta_h: float = DELTA_H,
                     ) -> bool:
    """
    Add ``delta_h`` (kJ/mol) to the constant term of one species' NASA7 polynomial(s).
    Returns ``False`` if the species' thermo is not in NASA7 form.

    NASA7 polynomials store enthalpy in dimensionless form: ``H(T)/(R*T)``.
    The 6th coefficient (``a6``) is the integration constant with units of K
    (effectively ``H_ref / R``). To shift enthalpy by ``delta_h`` kJ/mol we
    convert to J/mol (``* 1e3``) then divide by ``R`` to match the
    dimensionless convention: ``a6 += delta_h * 1e3 / R``.
    """
    content = read_yaml_file(original_path)
    for i in range(len(content['species'])):
        if i != species_index:
            continue
        if content['species'][i]['thermo']['model'] != 'NASA7':
            if logger is not None:
                logger.warning(f"Species '{content['species'][i]['name']}' does not use the expected NASA "
                               f"polynomials format for thermo, not perturbing it.")
            return False
        content['species'][i]['thermo']['data'][0][5] += delta_h * 1e3 / R
        if len(content['species'][i]['thermo']['data']) == 2:
            content['species'][i]['thermo']['data'][1][5] += delta_h * 1e3 / R
        break
    save_yaml_file(perturbed_path, content)
    return True


def perturb_reaction_rate_coefficient(original_path: str,
                                      perturbed_path: str,
                                      reaction_index: int,
                                      logger: Optional[Logger] = None,
                                      delta_k: float = DELTA_K,
                                      ) -> bool:
    """
    Multiply the pre-exponential factor of one reaction's rate expression by ``(1 + delta_k)``.
    Handles Arrhenius / three-body, falloff, PLOG, and Chebyshev forms.
    """
    content = read_yaml_file(original_path)
    for i in range(len(content['reactions'])):
        if i != reaction_index:
            continue
        rxn = content['reactions'][i]
        if 'rate-constant' in rxn and 'A' in rxn['rate-constant']:
            rxn['rate-constant']['A'] *= (1 + delta_k)
        elif rxn.get('type') == 'falloff':
            rxn['low-P-rate-constant']['A'] *= (1 + delta_k)
            rxn['high-P-rate-constant']['A'] *= (1 + delta_k)
        elif rxn.get('type') == 'pressure-dependent-Arrhenius':
            for j in range(len(rxn['rate-constants'])):
                rxn['rate-constants'][j]['A'] *= (1 + delta_k)
        elif rxn.get('type') == 'Chebyshev':
            rxn['data'][0][0] += math.log10(1 + delta_k)
        else:
            if logger is not None:
                logger.warning(f"Reaction '{rxn['equation']}' does not use a recognized "
                               f"rate expression, not perturbing it.")
            return False
        break
    save_yaml_file(perturbed_path, content)
    return True


def compute_idt_sa(reactor_idt_dict: dict,
                   perturbed_idt_dict: dict,
                   delta_h: float = DELTA_H,
                   delta_k: float = DELTA_K,
                   ) -> dict:
    """
    Convert raw perturbed-IDT samples into normalized SA coefficients:

    - kinetics: ``dln(IDT)/dln(k_i) = (IDT_pert - IDT) / (IDT * delta_k)`` (dimensionless)
    - thermo:   ``dln(IDT)/dH = (IDT_pert - IDT) / (IDT * delta_h)`` (mol/kJ)

    Returns a dict shaped ``{kind: {'IDT': {r: {phi: {P: {T: {param_index: coeff}}}}}}}``.
    """
    idt_sa_dict: dict = dict()
    for token in ('thermo', 'kinetics'):
        idt_sa_dict[token] = {'IDT': dict()}
        for r, reactor_idt_data in reactor_idt_dict.items():
            idt_sa_dict[token]['IDT'][r] = dict()
            for phi, phi_data in reactor_idt_data.items():
                idt_sa_dict[token]['IDT'][r][phi] = dict()
                for p, p_data in phi_data.items():
                    idt_sa_dict[token]['IDT'][r][phi][p] = dict()
                    for t, idt in p_data.items():
                        idt_sa_dict[token]['IDT'][r][phi][p][t] = dict()
                        if idt is None or idt == 0:
                            continue
                        for index, perturbed_idt_data in perturbed_idt_dict[token]['IDT'].items():
                            try:
                                perturbed_idt_value = perturbed_idt_data[r][phi][p][t]
                            except (KeyError, TypeError):
                                continue
                            if perturbed_idt_value is None:
                                continue
                            delta_idt = perturbed_idt_value - idt
                            if token == 'kinetics':
                                sa_coeff = delta_idt / (idt * delta_k)
                            else:
                                sa_coeff = delta_idt / (idt * delta_h)
                            idt_sa_dict[token]['IDT'][r][phi][p][t][index] = sa_coeff
    return idt_sa_dict


def get_top_sa_coefficients(idt_sa_dict: dict,
                            top_species: int = 10,
                            top_reactions: int = 10,
                            ) -> dict:
    """
    Keep only the top-N (by absolute value) sensitivity coefficients per
    (kind, reactor, phi, P, T) tuple.
    """
    top_sa_dict: dict = dict()
    for token in ('thermo', 'kinetics'):
        n = top_species if token == 'thermo' else top_reactions
        top_sa_dict[token] = {'IDT': dict()}
        for r, reactor_idt_data in idt_sa_dict[token]['IDT'].items():
            top_sa_dict[token]['IDT'][r] = dict()
            for phi, phi_data in reactor_idt_data.items():
                top_sa_dict[token]['IDT'][r][phi] = dict()
                for p, p_data in phi_data.items():
                    top_sa_dict[token]['IDT'][r][phi][p] = dict()
                    for t, idt_data in p_data.items():
                        top_indices = sorted(idt_data, key=lambda x: abs(idt_data[x]), reverse=True)[:n]
                        top_sa_dict[token]['IDT'][r][phi][p][t] = {idx: idt_data[idx] for idx in top_indices}
    return top_sa_dict


def get_h298(model: ct.Solution, species_index: int) -> float:
    """
    Standard enthalpy of formation at 298 K for one species, in kJ/mol.
    """
    model.TP = 298, 1e5
    return model.standard_enthalpies_RT[species_index] * ct.gas_constant * 298 / 1e6


def calculate_arrhenius_rate_coefficient(A: float, n: float, Ea: float, T: float, Ea_units: str) -> float:
    """
    Standard ``k(T) = A * T^n * exp(-Ea / (R T))``. ``Ea`` is given in ``Ea_units``.
    """
    if Ea_units not in EA_UNIT_CONVERSION:
        raise ValueError(f"Unsupported Ea units: {Ea_units}")
    return A * (T ** n) * math.exp(-1 * (Ea * EA_UNIT_CONVERSION[Ea_units]) / (R * T))


def calculate_troe_rate_coefficient(reaction_data: dict, T: float, P: float, Ea_units: str) -> float:
    """
    Troe (or Lindemann if no Troe params) falloff rate coefficient.
    """
    high_params = reaction_data['high-P-rate-constant']
    low_params = reaction_data['low-P-rate-constant']
    k0 = calculate_arrhenius_rate_coefficient(A=low_params['A'], n=low_params['b'],
                                              Ea=low_params['Ea'], T=T, Ea_units=Ea_units)
    kinf = calculate_arrhenius_rate_coefficient(A=high_params['A'], n=high_params['b'],
                                                Ea=high_params['Ea'], T=T, Ea_units=Ea_units)
    C = P / (R * T * 10)  # bath gas concentration in mol/cm^3
    Pr = k0 * C / kinf
    F = 1.0
    if 'Troe' in reaction_data:
        troe_params = reaction_data['Troe']
        alpha = troe_params.get('A', 0.0)
        T1 = troe_params.get('T1', 1e+30)
        T2 = troe_params.get('T2', 1e+30)
        T3 = troe_params.get('T3', 1e-30)
        Fcent = (1 - alpha) * math.exp(-T / T3) + alpha * math.exp(-T / T1)
        Fcent += math.exp(-T2 / T) if T2 != 0.0 else 0.0
        c = -0.4 - 0.67 * math.log10(Fcent)
        n_param = 0.75 - 1.27 * math.log10(Fcent)
        d = 0.14
        F = 10.0 ** (math.log10(Fcent) /
                     (1 + ((math.log10(Pr) + c) / (n_param - d * (math.log10(Pr) + c))) ** 2))
    return kinf * (Pr / (1 + Pr)) * F


def calculate_plog_rate_coefficient(reaction_data: dict, T: float, P: float, Ea_units: str) -> float:
    """
    PLOG (pressure-dependent piecewise Arrhenius) rate coefficient at ``(T, P)``.
    Clamps to the boundary rate if ``P`` falls outside the tabulated range.
    """
    rate_constants = reaction_data['rate-constants']
    n_entries = len(rate_constants)
    # Boundary clamping: return lowest/highest tabulated rate if P is out of range
    p_first = get_pressure_from_cantera(rate_constants[0]['P'])
    p_last = get_pressure_from_cantera(rate_constants[-1]['P'])
    if P <= p_first:
        return calculate_arrhenius_rate_coefficient(A=rate_constants[0]['A'], n=rate_constants[0]['b'],
                                                    Ea=rate_constants[0]['Ea'], T=T, Ea_units=Ea_units)
    if P >= p_last:
        return calculate_arrhenius_rate_coefficient(A=rate_constants[-1]['A'], n=rate_constants[-1]['b'],
                                                    Ea=rate_constants[-1]['Ea'], T=T, Ea_units=Ea_units)
    # Find the bracketing pressure interval
    i = 0
    for i in range(n_entries - 1):
        p_low = get_pressure_from_cantera(rate_constants[i]['P'])
        p_high = get_pressure_from_cantera(rate_constants[i + 1]['P'])
        if p_low <= P <= p_high:
            break
    k_low = calculate_arrhenius_rate_coefficient(A=rate_constants[i]['A'],
                                                 n=rate_constants[i]['b'],
                                                 Ea=rate_constants[i]['Ea'],
                                                 T=T, Ea_units=Ea_units)
    k_high = calculate_arrhenius_rate_coefficient(A=rate_constants[i + 1]['A'],
                                                  n=rate_constants[i + 1]['b'],
                                                  Ea=rate_constants[i + 1]['Ea'],
                                                  T=T, Ea_units=Ea_units)
    if P == p_low:
        return k_low
    if P == p_high:
        return k_high
    return k_low * 10 ** (math.log10(P / p_low) / math.log10(p_high / p_low) * math.log10(k_high / k_low))


def calculate_chebyshev_rate_coefficient(reaction_data: dict, T: float, P: float) -> float:
    """
    Chebyshev rate coefficient at ``(T, P)``.
    """
    T_min, T_max = reaction_data['temperature-range'][0], reaction_data['temperature-range'][1]
    P_min = get_pressure_from_cantera(reaction_data['pressure-range'][0])
    P_max = get_pressure_from_cantera(reaction_data['pressure-range'][1])
    Tr = (2.0 / T - 1.0 / T_min - 1.0 / T_max) / (1.0 / T_max - 1.0 / T_min)
    Pr = (2.0 * math.log10(P) - math.log10(P_min) - math.log10(P_max)) / (math.log10(P_max) - math.log10(P_min))
    n_T = len(reaction_data['data'])
    n_P = len(reaction_data['data'][0])
    k = 0.0
    for t in range(n_T):
        for p in range(n_P):
            k += reaction_data['data'][t][p] * chebyshev(t, Tr) * chebyshev(p, Pr)
    return 10 ** k


def chebyshev(n: int, x: float) -> float:
    """
    Chebyshev polynomial of the first kind, ``T_n(x)``.
    """
    if n == 0:
        return 1
    if n == 1:
        return x
    t, t0, t1 = 0.0, 1, x
    for _ in range(1, n):
        t = 2 * x * t1 - t0
        t0 = t1
        t1 = t
    return t


def get_pressure_from_cantera(p: str) -> float:
    """
    Parse a Cantera pressure string (e.g. ``"1 atm"``) and return the value in bar.
    """
    p_units = p.split()[1]
    p_factor = P_UNIT_CONVERSION[p_units]
    return float(p.split()[0]) * p_factor


def get_rate_coefficient(reaction_data: dict,
                         T: float,
                         P: Optional[float] = 1,
                         Ea_units: str = 'J/mol',
                         ) -> float:
    """
    Compute the rate coefficient at ``(T, P)`` for any supported reaction expression
    (Arrhenius / three-body, falloff, PLOG, Chebyshev).
    """
    rate_info = reaction_data.get('rate-constant', {})
    if (('A' in rate_info and 'b' in rate_info and 'Ea' in rate_info and len(rate_info) == 3)
            or reaction_data.get('type') == 'three-body'):
        return calculate_arrhenius_rate_coefficient(A=rate_info['A'], n=rate_info['b'],
                                                    Ea=rate_info['Ea'], T=T, Ea_units=Ea_units)
    if reaction_data.get('type') == 'falloff':
        return calculate_troe_rate_coefficient(reaction_data, T, P, Ea_units)
    if reaction_data.get('type') == 'pressure-dependent-Arrhenius':
        return calculate_plog_rate_coefficient(reaction_data, T, P, Ea_units)
    if reaction_data.get('type') == 'Chebyshev':
        return calculate_chebyshev_rate_coefficient(reaction_data, T, P)
    raise ValueError("Unsupported reaction type or missing parameters.")


def get_Ea_units(path: str) -> str:
    """
    Read the activation-energy units declared at the top of a Cantera YAML file
    (defaults to ``'kcal/mol'`` if the file omits the ``units`` block).
    """
    content = read_yaml_file(path)
    units = content.get('units', {})
    Ea_units = units.get('activation-energy', 'kcal/mol')
    if Ea_units not in ('J/mol', 'kJ/mol', 'cal/mol', 'kcal/mol'):
        raise ValueError(f"Unsupported Ea units: {Ea_units} read from file {path}")
    return Ea_units


register_simulate_adapter("CanteraIDT", CanteraIDT)
