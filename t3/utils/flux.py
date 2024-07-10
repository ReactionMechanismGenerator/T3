"""
t3 utils flux module
"""

from typing import Dict, List, Optional, Set, Tuple

import os

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pydot

from t3.utils.fix_cantera import fix_cantera


def generate_flux(model_path: str,
                  folder_path: str,
                  observables: List[str],
                  times: List[float],
                  composition: Dict[str, float],
                  T: float,
                  P: float,
                  V: Optional[float] = None,
                  reactor_type: str = 'JSR',
                  energy: bool = False,
                  a_tol: float = 1e-16,
                  r_tol: float = 1e-10,
                  explore_tol: float = 0.95,
                  dead_end_tol: float = 0.10,
                  generate_separate_diagrams_per_observable: bool = False,
                  report_flux_ratio: bool = True,
                  report_actual_flux: bool = False,
                  display_concentrations: bool = True,
                  display_r_n_p: bool = True,
                  scaling: Optional[float] = None,
                  fix_cantera_model: bool = True,
                  allowed_nodes: Optional[List[str]] = None,
                  max_chemical_generations: Optional[int] = None,
                  ):
    """
    Generate a flux diagram for a given model and composition.
    The flux diagram is generated for a specific time using the given reactor type.

    Args:
        model_path (str): The path to the cantera YAML model file.
        folder_path (str): The path to the folder in which to save the flux diagrams and accompanied data.
        observables (List[str]): The species to generate ROP diagrams for.
        times (List[float]): The times to generate flux diagrams for, in seconds.
                             For a JSR, these times represent residence times.
        composition (Dict[str, float]): The composition of the mixture. Keys are species labels as appears in the model,
                                        values are mole fractions.
        T (float): The temperature of the mixture, in Kelvin.
        P (float): The pressure of the mixture, in bar.
        V (Optional[float], optional): The reactor volume in cm^3, if relevant.
        reactor_type (str, optional): The reactor type. Supported reactor types are:
                                      'JSR': Jet stirred reactor, which is a CSTR with constant T/P/V
                                      'BatchP': An ideal gas constant pressure and constant volume batch reactor
        energy (bool, optional): Whether to turn energy equations on (False means adiabatic conditions).
        a_tol (float, optional): The absolute tolerance for the simulation.
        r_tol (float, optional): The relative tolerance for the simulation.
        explore_tol (float, optional): The minimal flux to capture of each species consumption pathway.
        dead_end_tol (float, optional): A flux exploration termination criterion.
                                        Don't explore further consumption is lower than this tolerance
                                        times the net rate of production.
        generate_separate_diagrams_per_observable (bool, optional): Whether to generate a separate flux diagram for each observable.
        report_flux_ratio (bool, optional): Whether to display the flux ratio.
        report_actual_flux (bool, optional): Whether to report the actual flux values rather than the relative flux.
        display_concentrations (bool, optional): Whether to display the concentrations.
        display_r_n_p (bool, optional): Whether to display the other reactants and products on each arrow.
        scaling (Optional[float], optional): The scaling of the final image, 100 means no scaling.
        fix_cantera_model (bool, optional): Whether to fix the Cantera model before running the simulation.
        allowed_nodes (Optional[List[str]], optional): A list of nodes to consider.
                                                       any node outside this list will not appear in the flux diagram.
        max_chemical_generations (Optional[int], optional): The maximal number of chemical generations to consider.

    Structures:
        profiles: {<time in s>: {'P': <pressure in bar>,
                                 'T': <temperature in K>,
                                 'X': {<species label>: <mole fraction>},
                                 'ROPs': {<species label>: {<reaction string>: <ROP>}}},
                   <time in s>: {...},
                  }
    """
    if not os.path.isdir(folder_path):
        os.makedirs(folder_path)
    if fix_cantera_model:
        fix_cantera(model_path=model_path)
    profiles = get_profiles_from_simulation(model_path=model_path,
                                            reactor_type=reactor_type,
                                            times=times,
                                            composition=composition,
                                            T=T,
                                            P=P,
                                            V=V,
                                            a_tol=a_tol,
                                            r_tol=r_tol,
                                            energy=energy,
                                            )

    generate_top_rop_bar_figs(profiles=profiles, observables=observables, folder_path=folder_path)

    folder_path = os.path.join(folder_path, 'flux_diagrams')
    if generate_separate_diagrams_per_observable:
        for observable in observables:
            folder_path_observable = os.path.join(folder_path, observable)
            if not os.path.isdir(folder_path):
                os.makedirs(folder_path_observable)
            generate_flux_diagrams(profiles=profiles,
                                   observables=[observable],
                                   folder_path=folder_path_observable,
                                   explore_tol=explore_tol,
                                   dead_end_tol=dead_end_tol,
                                   display_concentrations=display_concentrations,
                                   report_flux_ratio=report_flux_ratio,
                                   report_actual_flux=report_actual_flux,
                                   display_r_n_p=display_r_n_p,
                                   scaling=scaling,
                                   allowed_nodes=allowed_nodes,
                                   max_chemical_generations=max_chemical_generations,
                                   )
    else:
        generate_flux_diagrams(profiles=profiles,
                               observables=observables,
                               folder_path=folder_path,
                               explore_tol=explore_tol,
                               dead_end_tol=dead_end_tol,
                               display_concentrations=display_concentrations,
                               report_flux_ratio=report_flux_ratio,
                               report_actual_flux=report_actual_flux,
                               display_r_n_p=display_r_n_p,
                               scaling=scaling,
                               allowed_nodes=allowed_nodes,
                               max_chemical_generations=max_chemical_generations,
                               )


def get_profiles_from_simulation(model_path: str,
                                 reactor_type: str,
                                 times: List[float],
                                 composition: Dict[str, float],
                                 T: float,
                                 P: float,
                                 V: Optional[float] = None,
                                 a_tol: float = 1e-16,
                                 r_tol: float = 1e-10,
                                 energy: bool = False,
                                 ) -> dict:
    """
    Get the profiles of the simulation.

    Args:
        model_path (str): The path to the cantera YAML model file.
        times (List[float]): The times to generate flux diagrams for, in seconds.
                             For a JSR, these times represent residence times.
        composition (Dict[str, float]): The composition of the mixture. Keys are species labels as appears in the model,
                                        values are mole fractions.
        T (float): The temperature of the mixture, in Kelvin.
        P (float): The pressure of the mixture, in bar.
        V (Optional[float], optional): The reactor volume in cm^3, if relevant.
        reactor_type (str, optional): The reactor type. Supported reactor types are:
                                      'JSR': Jet stirred reactor, which is a CSTR with constant T/P/V
                                      'BatchP': An ideal gas constant pressure and constant volume batch reactor
        a_tol (float, optional): The absolute tolerance for the simulation.
        r_tol (float, optional): The relative tolerance for the simulation.
        energy (bool, optional): Whether to turn energy equations on (False means adiabatic conditions).

    Returns:
        dict: The profiles of the simulation.
    """
    gas = ct.Solution(model_path)
    mark_rxn_duplicity_correctly(gas)
    if reactor_type == 'BatchP':
        profiles = run_batch_p(gas=gas,
                               times=times,
                               composition=composition,
                               T=T,
                               P=P,
                               energy=energy,
                               max_steps=1e5,
                               a_tol=a_tol,
                               r_tol=r_tol,
                               )
    elif reactor_type == 'JSR':
        profiles = run_jsr(gas=gas,
                           times=times,
                           composition=composition,
                           T=T,
                           P=P,
                           V=V,
                           a_tol=a_tol,
                           r_tol=r_tol,
                           )
    else:
        raise NotImplementedError(f'Reactor type {reactor_type} is not yet supported')
    return profiles


def get_rxn_stoichiometry(gas: ct.Solution) -> Dict[str, List[float]]:
    """
    Get the stoichiometry of each species in all reactions in the model.

    Args:
        gas (ct.Solution): The cantera Solution object.

    Returns:
        Dict[str, List[float]]: Keys are species labels, values are stoichiometry coefficients per reaction in the model.
    """
    stoichiometry = {spc.name: list() for spc in gas.species()}
    for spc in gas.species():
        for i, rxn in enumerate(gas.reactions()):
            nu = rxn.products.get(spc.name, 0) - rxn.reactants.get(spc.name, 0)
            stoichiometry[spc.name].append(nu)
    return stoichiometry


def mark_rxn_duplicity_correctly(gas: ct.Solution):
    """
    Mark duplicate reactions correctly in the cantera Solution object.

    Args:
        gas (ct.Solution): The cantera Solution object.
    """
    rxn_equations = [r.equation for r in gas.reactions()]
    for rxn in gas.reactions():
        if rxn_equations.count(rxn.equation) > 1:
            rxn.duplicate = True


def set_jsr(gas: ct.Solution,
            time: float,
            composition: Dict[str, float],
            T: float,
            P: float,
            V: float,
            a_tol: float = 1e-16,
            r_tol: float = 1e-10,
            ) -> Tuple[ct.ReactorNet, ct.IdealGasReactor]:
    """
    Define a jet stirred reactor.

    Args:
        gas (ct.Solution): The cantera Solution object.
        time (float): The residence time in s.
        composition (Dict[str, float]): The composition of the mixture.
        T (float): The temperature in K.
        P (float): The pressure in Pa.
        V (float): The reactor volume in cm^3.
        a_tol (float, optional): The absolute tolerance for the simulation.
        r_tol (float, optional): The relative tolerance for the simulation.

    Returns:
        Tuple[ct.ReactorNet, ct.IdealGasReactor]:
            - The reactor network (consisting of the JSR)
            - The JSR
    """
    gas.TPX = T, P * 1e5, composition
    inlet = ct.Reservoir(gas)
    exhaust = ct.Reservoir(gas)
    stirred_reactor = ct.IdealGasReactor(gas, energy="off", volume=V * 1e-6)
    mfc = ct.MassFlowController(
        upstream=inlet,
        downstream=stirred_reactor,
        mdot=stirred_reactor.mass / time,
    )
    mfc.mass_flow_coeff = 1.0
    pr = ct.PressureController(
        upstream=stirred_reactor,
        downstream=exhaust,
        master=mfc,
    )
    pr.pressure_coeff = 0.01
    network = ct.ReactorNet([stirred_reactor])
    network.atol = a_tol
    network.rtol = r_tol
    return network, stirred_reactor


def set_batch_p(gas: ct.Solution,
                composition: Dict[str, float],
                T: float,
                P: float,
                energy: bool,
                a_tol: float = 1e-16,
                r_tol: float = 1e-10,
                ) -> Tuple[ct.ReactorNet, ct.IdealGasConstPressureReactor]:
    """
    Define a batch reactor with constant pressure and constant volume.

    Args:
        gas (ct.Solution): The cantera Solution object.
        composition (Dict[str, float]): The composition of the mixture.
        T (float): The temperature in K.
        P (float): The pressure in bar.
        energy (bool): Whether to turn energy equations on (False means adiabatic conditions).
        a_tol (float, optional): The absolute tolerance for the simulation.
        r_tol (float, optional): The relative tolerance for the simulation.

    Returns:
        Tuple[ct.ReactorNet, ct.IdealGasConstPressureReactor]:
            - The reactor network (consisting of the batch reactor)
            - The batch reactor
    """
    gas.TPX = T, P * 1e5, composition
    reactor = ct.IdealGasConstPressureReactor(gas, energy="on" if energy else "off")
    network = ct.ReactorNet([reactor])
    network.atol = a_tol
    network.rtol = r_tol
    return network, reactor


def run_jsr(gas: ct.Solution,
            times: List[float],
            composition: Dict[str, float],
            T: float,
            P: float,
            V: Optional[float] = None,
            a_tol: float = 1e-16,
            r_tol: float = 1e-10,
            ) -> Dict[float, dict]:
    """
    Run a JSR reactor with constant pressure, constant temperature, and constant volume.

    Args:
        gas (ct.Solution): The cantera Solution object.
        times (List[float]): The residence times to generate flux diagrams for, in seconds.
        composition (Dict[str, float]): The composition of the mixture.
        T (float): The temperature in K.
        P (float): The pressure in bar.
        V (float, optional): The reactor volume in cm^3.
        a_tol (float, optional): The absolute tolerance for the simulation.
        r_tol (float, optional): The relative tolerance for the simulation.

    Returns:
        dict: The T, P, X, and ROP profiles (values) at a specific time.
    """
    V = V or 100
    profiles = dict()
    stoichiometry = get_rxn_stoichiometry(gas)
    for tau in times:
        network, reactor = set_jsr(gas=gas, time=tau, composition=composition, T=T, P=P, V=V, a_tol=a_tol, r_tol=r_tol)
        t = 0
        while t < tau:
            rops = {spc.name: dict() for spc in gas.species()}
            t = network.step()
            if t > tau:
                network.advance(tau)
                t = tau
            cantera_reaction_rops = gas.net_rates_of_progress
            for spc in gas.species():
                for i, rxn in enumerate(gas.reactions()):
                    if stoichiometry[spc.name][i]:
                        if rxn.equation not in rops[spc.name].keys():
                            rops[spc.name][rxn.equation] = 0
                        rops[spc.name][rxn.equation] += cantera_reaction_rops[i] * stoichiometry[spc.name][i]
            profile = {'P': gas.P, 'T': gas.T, 'X': {s.name: x for s, x in zip(gas.species(), gas.X)}, 'ROPs': rops}
            if t == tau:
                profiles[t] = profile
    return profiles


def run_batch_p(gas: ct.Solution,
                times: List[float],
                composition: Dict[str, float],
                T: float,
                P: float,
                energy: bool,
                max_steps: float = 1e5,
                a_tol: float = 1e-16,
                r_tol: float = 1e-10,
                ) -> Dict[float, dict]:
    """
    Run a batch reactor with constant pressure and constant volume.

    Args:
        gas (ct.Solution): The cantera Solution object.
        T (float): The temperature in K.
        P (float): The pressure in bar.
        energy (bool): Whether to turn energy equations on (False means adiabatic conditions).
        composition (Dict[str, float]): The composition of the mixture. Keys are species labels as appears in the model.
        a_tol (float, optional): The absolute tolerance for the simulation.
        r_tol (float, optional): The relative tolerance for the simulation.
        times (List[float]): The times to generate flux diagrams for, in seconds.
        max_steps (float, optional): The maximal number of steps to take.

    Returns:
        Dict[float, dict]: The T, P, X, and ROP profiles (values) at a specific time.
    """
    network, reactor = set_batch_p(gas=gas, composition=composition, T=T, P=P, energy=energy, a_tol=a_tol, r_tol=r_tol)
    rops = {spc.name: dict() for spc in gas.species()}
    stoichiometry = get_rxn_stoichiometry(gas)
    profiles = dict()
    tau_i = 0
    step_counter = 0
    while tau_i < len(times):
        t = network.step()
        if t > times[tau_i]:
            network.advance(times[tau_i])
            t = times[tau_i]
            tau_i += 1
        step_counter += 1
        if step_counter > max_steps:
            break
        cantera_reaction_rops = gas.net_rates_of_progress
        for spc in gas.species():
            dups = list()
            for i, rxn in enumerate(gas.reactions()):
                if stoichiometry[spc.name][i]:
                    if rxn.equation not in rops[spc.name].keys():
                        rops[spc.name][rxn.equation] = list()
                    if rxn.duplicate and rxn.equation in dups:
                        rops[spc.name][rxn.equation][-1] += cantera_reaction_rops[i] * stoichiometry[spc.name][i]
                    else:
                        rops[spc.name][rxn.equation].append(cantera_reaction_rops[i] * stoichiometry[spc.name][i])
                    if rxn.duplicate and rxn.equation not in dups:
                        dups.append(rxn.equation)
        profile = {'P': gas.P, 'T': gas.T, 'X': {s.name: x for s, x in zip(gas.species(), gas.X)}, 'ROPs': rops}
        profiles[t] = profile
    return profiles


def get_top_rops(profiles: dict,
                 observables: List[str],
                 top_rops_to_plot: int = 10,
                 orders_of_magnitude_to_consider: int = 3,
                 ) -> dict:
    """
    Get the top ROPs for the observables.

    Args:
        profiles (dict): The T, P, X, and ROP profiles (values) at a specific time.
        observables (List[str]): The species to generate ROP diagrams for.
        top_rops_to_plot (int, optional): The number of top reaction ROPs to plot per figure.
        orders_of_magnitude_to_consider (int, optional): The number of orders of magnitude to consider for the ROPs.

    Returns:
        dict: The top ROPs for the observables.
    """
    top_rops = dict()
    for time, profile in profiles.items():
        top_rops[time] = {observable: dict() for observable in observables}
        for spc in observables:
            max_rop_value_per_rxn = {rxn: max(rops) if isinstance(rops, list) else rops for rxn, rops in profile['ROPs'][spc].items()}
            max_rop_rxns_sorted = [k for k, v in sorted(max_rop_value_per_rxn.items(), key=lambda item: abs(item[1]),
                                                        reverse=True)]
            max_rop_val = max([abs(v) for v in max_rop_value_per_rxn.values()])
            for i in range(top_rops_to_plot):
                if len(max_rop_rxns_sorted) <= i:
                    break
                if abs(profile['ROPs'][spc][max_rop_rxns_sorted[i]]) < max_rop_val / 10 ** orders_of_magnitude_to_consider:
                    continue
                top_rops[time][spc][max_rop_rxns_sorted[i]] = profile['ROPs'][spc][max_rop_rxns_sorted[i]]
    return top_rops


def generate_top_rop_bar_figs(profiles: dict,
                              observables: List[str],
                              folder_path: str,
                              ) -> None:
    """
    Generate top ROP figures for the observables.

    Args:
        profiles (dict): The T, P, X, and ROP profiles (values) at a specific time.
        observables (List[str]): The species to generate ROP diagrams for.
        folder_path (str): The path to the folder in which to save the flux diagrams and accompanied data.
    """
    base_path = os.path.join(folder_path, 'bar_ROPs')
    if not os.path.isdir(base_path):
        os.makedirs(base_path)
    top_rops = get_top_rops(profiles=profiles, observables=observables)
    for time, rop_dict in top_rops.items():
        for observable, rxn_dict in rop_dict.items():
            fig, ax = plt.subplots(figsize=(12, 8))
            ax.set_title(f'{observable} at {time} s')
            ax.set_xlabel('Reactions')
            ax.set_ylabel('ROP, mol/cm^3/s')
            ax.bar(rxn_dict.keys(), rxn_dict.values())
            ax.bar('Total ROP', sum(rxn_dict.values()))
            ax.barlogy = True
            plt.xticks(rotation=60)
            plt.tight_layout()
            fig_path = os.path.join(base_path, f'{observable}_{time}_s.png')
            plt.savefig(fig_path, dpi=300, facecolor='w', edgecolor='w', orientation='portrait', format='png',
                        transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)
            plt.close()


def generate_flux_diagrams(profiles: dict,
                           observables: List[str],
                           folder_path: str,
                           explore_tol: float = 0.95,
                           dead_end_tol: float = 0.10,
                           display_concentrations: bool = True,
                           report_flux_ratio: bool = True,
                           report_actual_flux: bool = False,
                           display_r_n_p: bool = True,
                           scaling: Optional[float] = None,
                           allowed_nodes: Optional[List[str]] = None,
                           max_chemical_generations: Optional[int] = None,
                           ):
    """
    Generate flux diagrams.

    Args:
        profiles (dict): The T, P, X, and ROP profiles (values) at a specific time.
        observables (List[str]): The species to start the flux diagram with.
        folder_path (str): The path to the folder in which to save the flux diagrams and accompanied data.
        explore_tol (float, optional): The minimal flux to capture of each species consumption pathway.
        dead_end_tol (float, optional): A flux exploration termination criterion.
                                        Don't explore further consumption is lower than this tolerance
                                        times the net rate of production.
        display_concentrations (bool, optional): Whether to display the concentrations.
        report_flux_ratio (bool, optional): Whether to display the flux ratio.
        report_actual_flux (bool, optional): Whether to report the actual flux values rather than the relative flux.
        display_r_n_p (bool, optional): Whether to display the other reactants and products on each arrow.
        scaling (Optional[float], optional): The scaling of the final image.
        allowed_nodes (Optional[List[str]], optional): A list of nodes to consider.
                                                       any node outside this list will not appear in the flux diagram.
        max_chemical_generations (Optional[int], optional): The maximal number of chemical generations to consider.

    Structures:
        graph: {<species1>: {'rxn1': [[<the species formed>], <rop_value>],
                            'rxn2': [[<the species formed>], <rop_value>],},
                <species2>: {...},
               }
    """
    if not os.path.isdir(folder_path):
        os.makedirs(folder_path)
    for time, profile in profiles.items():
        flux_graph, nodes_to_explore, min_rop, max_rop = get_flux_graph(profile=profile,
                                                                        observables=observables,
                                                                        explore_tol=explore_tol,
                                                                        dead_end_tol=dead_end_tol,
                                                                        max_chemical_generations=max_chemical_generations,
                                                                        )
        create_digraph(flux_graph=flux_graph,
                       profile=profile,
                       observables=observables,
                       nodes_to_explore=nodes_to_explore,
                       time=time,
                       min_rop=min_rop,
                       max_rop=max_rop,
                       folder_path=folder_path,
                       display_concentrations=display_concentrations,
                       report_flux_ratio=report_flux_ratio,
                       report_actual_flux=report_actual_flux,
                       display_r_n_p=display_r_n_p,
                       scaling=scaling,
                       allowed_nodes=allowed_nodes,
                       )


def create_digraph(flux_graph: dict,
                   profile: dict,
                   observables: List[str],
                   nodes_to_explore: Set[str],
                   time: float,
                   min_rop: float,
                   max_rop: float,
                   folder_path: str,
                   display_concentrations: bool = True,
                   report_flux_ratio: bool = True,
                   report_actual_flux: bool = False,
                   display_r_n_p: bool = True,
                   scaling: Optional[float] = None,
                   allowed_nodes: Optional[List[str]] = None,
                   ) -> None:
    """
    Create a directed graph from the flux graph and save it as a .dot file.

    Args:
        flux_graph (dict): The flux graph.
        profile (dict): The T, P, X, and ROP profiles (values) at a specific time.
        observables (List[str]): The species to start the flux diagram with.
        nodes_to_explore (Set[str]): Nodes that have additional downstream nodes (that are not dead end).
        time (float): The time in seconds.
        min_rop (float): The absolute minimal ROP value.
        max_rop (float): The absolute maximal ROP value.
        folder_path (str): The path to the folder in which to save the flux diagrams and accompanied data.
        display_concentrations (bool, optional): Whether to display the concentrations.
        report_flux_ratio (bool, optional): Whether to display the flux ratio.
        report_actual_flux (bool, optional): Whether to report the actual flux values rather than the relative flux.
        display_r_n_p (bool, optional): Whether to display the other reactants and products on each arrow.
        scaling (Optional[float], optional): The scaling of the final image.
        allowed_nodes (Optional[List[str]], optional): A list of nodes to consider.
                                                       any node outside this list will not appear in the flux diagram.
    """
    if not os.path.isdir(folder_path):
        os.makedirs(folder_path)

    graph = pydot.Dot(graph_type='digraph')
    stack, visited = [o for o in observables], set()
    species_to_consider = set(list(flux_graph.keys()))
    for rxn_dict in flux_graph.values():
        for rop_list in rxn_dict.values():
            species_to_consider.update(rop_list[0])
    xs = [v for k, v in profile['X'].items() if k in species_to_consider]
    if not len(xs):
        print(f'Could not create a flux diagram for observables {observables} at {time} s. '
              f'Could not simulate the system.')
        return
    x_max, x_min = max(xs), min(xs)
    abs_rops = [abs(values[1]) for inner_dict in flux_graph.values() for values in inner_dict.values()]
    rop_min, rop_max = min(abs_rops), max(abs_rops)
    nodes = dict()
    while len(stack):
        origin_label = stack.pop(-1)
        origin_node = get_node(graph=graph,
                               label=origin_label,
                               nodes=nodes,
                               observables=observables,
                               width=get_width(x=profile['X'][origin_label], x_min=x_min, x_max=x_max),
                               concentration=profile['X'][origin_label],
                               display_concentrations=display_concentrations,
                               )
        rxns_rop = flux_graph[origin_label] if origin_label in flux_graph.keys() else dict()
        for rxn, rop_list in rxns_rop.items():
            rxn_string = get_rxn_in_relevant_direction(rxn=rxn, spc=origin_label)
            downstream_node_labels, rop = rop_list[0], rop_list[1]
            downstream_node_labels, multipliers = unpack_stoichiometry(downstream_node_labels)
            for downstream_node_label in downstream_node_labels:
                if downstream_node_label in nodes_to_explore and downstream_node_label not in visited:
                    visited.add(downstream_node_label)
                    if allowed_nodes is None or downstream_node_label in allowed_nodes:
                        stack.append(downstream_node_label)
            downstream_nodes = [get_node(graph=graph,
                                         label=downstream_node_label,
                                         nodes=nodes,
                                         observables=observables,
                                         width=get_width(x=profile['X'][downstream_node_label], x_min=x_min, x_max=x_max),
                                         concentration=profile['X'][downstream_node_label],
                                         display_concentrations=display_concentrations,
                                         )
                                for downstream_node_label in downstream_node_labels]

            add_edges(graph=graph,
                      origin_node=origin_node,
                      origin_label=origin_label,
                      downstream_nodes=downstream_nodes,
                      downstream_node_labels=downstream_node_labels,
                      width=get_width(x=rop, x_min=rop_min, x_max=rop_max),
                      rop=rop,
                      max_rop=max_rop,
                      rxn=rxn_string,
                      multipliers=multipliers,
                      report_flux_ratio=report_flux_ratio,
                      report_actual_flux=report_actual_flux,
                      display_r_n_p=display_r_n_p,
                      allowed_nodes=allowed_nodes,
                      )
    graph.set(name='label', value=f'Flux diagram at {time} s, ROP range: [{min_rop:.2e}, {max_rop:.2e}] ' +
                                  'mol/cm\N{SUPERSCRIPT THREE}/s)')
    if scaling is not None:
        graph.set('size', f'{scaling},{scaling}')
    if allowed_nodes is not None:
        for node in graph.get_nodes():
            if node.get_name() not in allowed_nodes:
                graph.del_node(node)
    graph_dot_path = os.path.join(folder_path, f'flux_diagram_{time}_s.dot')
    graph_png_path = os.path.join(folder_path, f'flux_diagram_{time}_s.png')
    graph.write(graph_dot_path)
    try:
        graph.write_png(graph_png_path)
    except AssertionError:
        print(f'Could not create a flux diagram for observables {observables} at {time} s.')


def add_edges(graph: pydot.Dot,
              origin_node: pydot.Node,
              origin_label: str,
              downstream_nodes: List[pydot.Node],
              downstream_node_labels: List[str],
              width: float,
              rop: float,
              max_rop: float,
              rxn: Optional[str] = None,
              multipliers: Optional[List[float]] = None,
              report_flux_ratio: bool = True,
              report_actual_flux: bool = False,
              display_r_n_p: bool = True,
              allowed_nodes: Optional[List[str]] = None,
              ):
    """
    Add edges to the graph.

    Args:
        graph (pydot.Dot): The graph.
        origin_node (pydot.Node): The origin node.
        origin_label (str): The origin node label.
        downstream_nodes (List[pydot.Node]): The downstream nodes.
        downstream_node_labels (List[str]): The downstream node labels.
        rxn (str): The reaction string.
        width (float): The edge width.
        rop (float): The normalized ROP value.
        max_rop (float): The maximal ROP value.
        multipliers (List[float]): The stoichiometric multipliers.
        report_flux_ratio (bool, optional): Whether to display the flux ratio.
        report_actual_flux (bool, optional): Whether to report the actual flux values rather than the relative flux.
        display_r_n_p (bool, optional): Whether to display the other reactants and products on each arrow.
        allowed_nodes (Optional[List[str]], optional): A list of nodes to consider.
                                                       any node outside this list will not appear in the flux diagram.
    """
    for multiplier, node, node_label in zip(multipliers, downstream_nodes, downstream_node_labels):
        if allowed_nodes is None or (origin_label in allowed_nodes and node_label in allowed_nodes):
            rs, ps = get_other_reactants_and_products(rxn=rxn, spcs=[origin_label, node_label])
            edge = pydot.Edge(origin_node, node, penwidth=width + np.log10(multiplier), fontsize=8)
            label = ''
            rop = abs(rop)
            if report_flux_ratio:
                label = f'{rop:.1f}' if rop > 0.1 else f'{rop:.1e}'
            elif report_actual_flux:
                actual_rop = rop * max_rop
                label = f'{actual_rop:.1f}' if actual_rop > 0.1 else f'{actual_rop:.1e}'
            if display_r_n_p and rs:
                label += f'\n{rs}'
            if display_r_n_p and ps:
                label += f'\n{ps}'
            if label != '':
                edge.set('label', label)
            edge.set('arrowhead', 'vee')
            graph.add_edge(edge)


def get_width(x: float,
              x_min: float,
              x_max: float,
              log_scale: bool = True,
              ) -> float:
    """
    Get the width of a node or an edge.

    Args:
        x (float): The concentration or ROP value.
        x_min (float): The minimal value in the range.
        x_max (float): The maximal value in the range.
        log_scale (bool, optional): Whether to use a log scale when determining the width.

    Returns:
        float: The resulting width.
    """
    max_width, min_width = 4, 0.2
    x, x_min, x_max = abs(x), abs(x_min), abs(x_max)
    if not log_scale:
        if x == x_min == x_max:
            return 1
        return min_width + (x - x_min) * (max_width - min_width) / (x_max - x_min)
    return get_width(x=-np.log10(x_min) - np.log10(x_max) + np.log10(x),
                     x_min=-np.log10(x_max),
                     x_max=-np.log10(x_min),
                     log_scale=False,
                     )


def get_rxn_in_relevant_direction(rxn: str,
                                  spc: str) -> str:
    """
    Get the reaction string in the relevant direction with the given species as one of the reactants.

    Args:
        rxn (str): The reaction string.
        spc (str): The species label.

    Returns:
        str: The reaction string in the relevant direction.
    """
    arrow = ' <=> ' if ' <=> ' in rxn else ' => '
    wells = rxn.split(arrow)
    counts = wells[0].count(spc), wells[1].count(spc)
    i = int(counts[1] > counts[0])
    return arrow.join([wells[i], wells[not i]])


def get_node(graph: pydot.Dot,
             label: str,
             nodes: dict,
             observables: Optional[List[str]] = None,
             width: Optional[float] = None,
             concentration: Optional[float] = None,
             display_concentrations: bool = True,
             ) -> pydot.Node:
    """
    Get an existing node from the graph or create a new one.

    Args:
        graph (pydot.Dot): The graph.
        label (str): The node label.
        nodes (dict): Existing nodes in the graph.
        observables (Optional[List[str]], optional): Species that should have a blue fill color.
        width (Optional[float], optional): The node width.
        concentration (Optional[float], optional): The node concentration.
        display_concentrations (bool, optional): Whether to display the species concentrations next to its circle.

    Returns:
        pydot.Node: The node.
    """
    colors = {'blue': '#DCE5F4'}
    fontsize = 8
    if label not in nodes.keys():
        if observables is not None and label in observables:
            node = pydot.Node(label,
                              style='filled',
                              fillcolor=colors['blue'],
                              fontsize=fontsize,
                              )
        else:
            node = pydot.Node(label,
                              fontsize=fontsize,
                              )
        if display_concentrations:
            node.set('xlabel', f'{concentration:.2e}' if concentration is not None else '')
        graph.add_node(node)
        nodes[label] = node
    else:
        node = nodes[label]
    if width is not None:
        node.set_penwidth(width)
        if np.isnan(width):
            node.set_penwidth(1)
    return node


def get_flux_graph(profile: dict,
                   observables: List[str],
                   explore_tol: float = 0.95,
                   dead_end_tol: float = 0.10,
                   max_chemical_generations: Optional[int] = None,
                   ) -> Tuple[dict, Set[str], float, float]:
    """
    Explore the ROP profiles and get the flux graph.
    Also get the list of nodes to continue exploring when drawing the graph and the min and max rop.

    Args:
        profile (dict): The T, P, X, and ROP profiles (values) at a specific time.
        observables (List[str]): The species to start the flux diagram with.
        explore_tol (float, optional): The minimal flux to capture of each species consumption pathway.
        dead_end_tol (float, optional): A flux exploration termination criterion.
                                        Don't explore further consumption is lower than this tolerance
                                        times the net rate of production.
        max_chemical_generations (Optional[int], optional): The maximal number of chemical generations to consider.

    Returns:
        Tuple[dict, Set[str], float, float]: The flux graph and the maximal flux.
    """
    normalized_fluxes, max_rop = get_normalized_flux_profile(profile=profile)
    min_rop = None
    graph, nodes_to_explore = dict(), set()
    stack, visited = [obs for obs in observables], [obs for obs in observables]
    node_generation_dict = {obs: 0 for obs in observables}
    while len(stack):
        node = stack.pop(-1)
        node = node.split()[-1]
        rxns_rop = normalized_fluxes[node]
        rxns_rop_sorted = [(k, v) for k, v in sorted(normalized_fluxes[node].items(), key=lambda item: item[1], reverse=False)]
        consumption = 0
        total_consumption = sum([v for k, v in rxns_rop_sorted if v < 0])
        reduced_rxns_rop_sorted = list()
        for k, v in rxns_rop_sorted:
            if v < 0:
                reduced_rxns_rop_sorted.append((k, v))
                consumption += v
                if consumption < total_consumption * explore_tol:
                    break
        for rxn_rop in reduced_rxns_rop_sorted:
            rxn, rop = rxn_rop[0], rxn_rop[1]
            if min_rop is None or abs(rop) < min_rop:
                min_rop = abs(rop)
            opposite_rxn_species = get_opposite_rxn_species(rxn=rxn, spc=node)
            for spc in opposite_rxn_species:
                if spc not in visited \
                        and (max_chemical_generations is None
                             or node_generation_dict[node] < max_chemical_generations - 1):
                    if continue_exploring(rops=rxns_rop, dead_end_tol=dead_end_tol):
                        stack.append(spc)
                        nodes_to_explore.add(spc)
                        if max_chemical_generations is not None:
                            node_generation_dict[spc] = node_generation_dict[node] + 1
                    visited.append(spc)
            if node not in graph.keys():
                graph[node] = dict()
            if rxn in graph[node].keys():
                graph[node][rxn][1] += rop
            else:
                graph[node][rxn] = [opposite_rxn_species, rop]
    min_rop = min_rop * max_rop if min_rop is not None else 0
    return graph, nodes_to_explore, min_rop, max_rop


def get_normalized_flux_profile(profile: dict) -> Tuple[dict, float]:
    """
    Get the normalized flux profile.

    Args:
        profile (dict): The T, P, X, and ROP profiles (values) at a specific time.

    Returns:
        Tuple[dict, float]: The normalized flux profile and the maximal flux.
    """
    spc_max_flux = list()
    for spc, flux_dict in profile['ROPs'].items():
        if flux_dict:
            spc_max_flux.append(max([abs(v) for v in flux_dict.values()]))
    max_flux = max(spc_max_flux)
    normalized_fluxes = {spc: {rxn: flux / max_flux for rxn, flux in flux_dict.items()}
                         for spc, flux_dict in profile['ROPs'].items()}
    return normalized_fluxes, max_flux


def get_other_reactants_and_products(rxn: str,
                                     spcs: List[str],
                                     ) -> Tuple[str, str]:
    """
    Get the reactants and products in a reaction except for the given species,
    the first is considered a reactant, and the second a product.

    Args:
        rxn (str): The reaction string.
        spcs (List[str]): A length 2 list with the labels of a reactant and a product.

    Returns:
        Tuple[str, str]: The reactants and products.
    """
    arrow = ' <=> ' if ' <=> ' in rxn else ' => '
    wells = rxn.split(arrow)
    i = wells[1].count(spcs[0]) > wells[0].count(spcs[0])
    rs, ps = list(), list()
    for well, j in zip([rs, ps], [i, not i]):
        for w in wells[j].split(' + '):
            for token in [' + M', ' (+M)', 'M', ' + ']:
                w = w.replace(token, '')
            if spcs[0] not in w and spcs[1] not in w and w != '':
                well.append(w)
    rs = '+ ' + ' + '.join(rs) if len(rs) else ''
    ps = '- ' + ' - '.join(ps) if len(ps) else ''
    return rs, ps


def get_opposite_rxn_species(rxn: str, spc: str) -> List[str]:
    """
    Get the species in a reaction opposite to the given species
    (if the given species is one of the reactants, get the products, and vice versa).

    Args:
        rxn (str): The reaction string.
        spc (str): The species label.

    Returns:
        List[str]: The species in a reaction opposite to the given species.
    """
    arrow = ' <=> ' if ' <=> ' in rxn else ' => '
    wells = rxn.split(arrow)
    counts = count_species_in_well(well=wells[0], spc=spc), count_species_in_well(well=wells[1], spc=spc)
    i = int(counts[0] > counts[1])
    species = wells[i].split(' + ')
    for token in [' + M', ' (+M)', 'M', ' + ']:
        species = [s.replace(token, '') for s in species]
    return [s for s in species if s != '']


def count_species_in_well(well: str,
                          spc: str,
                          ) -> int:
    """
    Count the number of times a species appears in a well.

    Args:
        well (str): The well string.
        spc (str): The species label.

    Returns:
        int: The number of times a species appears in the well.
    """
    count = 0
    for token in [' + M', ' (+M)', 'M']:
        well = well.replace(token, '')
    splits = well.split(' + ')
    for s in splits:
        if s == spc:
            count += 1
    return count


def unpack_stoichiometry(labels: List[str]) -> Tuple[List[str], List[int]]:
    """
    Unpack stoichiometry.

    Args:
        labels (List[str]): The labels of the species in the reaction.

    Returns:
        Tuple[List[str], List[int]]: The unpacked labels and multipliers.
    """
    new_labels, multipliers = list(), list()
    for label in labels:
        splits = label.split()
        if len(splits) > 1:
            new_labels.append(splits[1])
            multipliers.append(int(splits[0]))
        else:
            new_labels.append(label)
            multipliers.append(1)
    return new_labels, multipliers


def continue_exploring(rops: Dict[str, float],
                       dead_end_tol: float = 0.10
                       ) -> bool:
    """
    Determine whether to continue exploring a species, or whether it is a dead end.

    Args:
        rops (Dict[str, float]): The net rate of production of the analyzed species. Keys are reaction strings,
                                 values are the net rate of production at the analysis time.
        dead_end_tol (float, optional): A flux exploration termination criterion.
    """
    production, consumption = 0, 0
    for rop in rops.values():
        production += rop if rop > 0 else 0
        consumption += rop if rop < 0 else 0
    if abs(consumption) < production * dead_end_tol:
        return False
    return True
