#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_flux module
"""

import os
import shutil

import cantera as ct
import numpy as np
import pydot

from t3.common import DATA_BASE_PATH, SCRATCH_BASE_PATH
import t3.utils.flux as flux
from tests.common import almost_equal


def test_generate_flux():
    """Test generating flux diagrams."""
    model_path = os.path.join(DATA_BASE_PATH, 'models', 'HOCHO.yaml')
    folder_path = os.path.join(SCRATCH_BASE_PATH, 'test_generate_flux_HOCHO')
    observables = ['H(3)', 'HOCO(10)']
    flux.generate_flux(model_path=model_path,
                       folder_path=folder_path,
                       observables=observables,
                       times=[0.001],
                       composition={'HOCHO(1)': 0.2, 'N2': 0.98},
                       T=1200,
                       P=1,
                       reactor_type='JSR',
                       explore_tol=0.95,
                       dead_end_tol=0.10,
                       generate_separate_diagrams_per_observable=False,
                       display_flux_ratio=True,
                       display_concentrations=True,
                       display_r_n_p=True,
                       fix_cantera_model=False,
                       )
    for fig in ['H(3)_0.001_s.png', 'HOCO(10)_0.001_s.png']:
        assert os.path.isfile(os.path.join(folder_path, 'bar_ROPs', fig))
    for ex in ['dot', 'png']:
        assert os.path.isfile(os.path.join(folder_path, 'flux_diagrams', f'flux_diagram_0.001_s.{ex}'))

    flux.generate_flux(model_path=model_path,
                       folder_path=folder_path,
                       observables=observables,
                       times=[0.001],
                       composition={'HOCHO(1)': 0.2, 'N2': 0.98},
                       T=1200,
                       P=1,
                       reactor_type='JSR',
                       explore_tol=0.95,
                       dead_end_tol=0.10,
                       generate_separate_diagrams_per_observable=True,
                       display_flux_ratio=True,
                       display_concentrations=True,
                       display_r_n_p=True,
                       fix_cantera_model=False,
                       )
    for observable in observables:
        for ex in ['dot', 'png']:
            assert os.path.isfile(os.path.join(folder_path, 'flux_diagrams', observable, f'flux_diagram_0.001_s.{ex}'))


def test_get_profiles_from_simulation():
    """Test getting profiles from a simulation."""
    profiles = flux.get_profiles_from_simulation(model_path=os.path.join(DATA_BASE_PATH, 'models', 'HOCHO.yaml'),
                                                 reactor_type='JSR',
                                                 times=[0.5, 5],
                                                 composition={'HOCHO(1)': 1.0},
                                                 T=1000,
                                                 P=1,
                                                 )
    assert len(profiles.keys()) == 2


def test_get_rxn_stoichiometry():
    """Test getting the stoichiometry of all species in all reactions"""
    model_path = os.path.join(DATA_BASE_PATH, 'models', 'HOCHO.yaml')
    gas = ct.Solution(model_path)
    stoichiometry = flux.get_rxn_stoichiometry(gas)
    assert all(s == 0 for s in stoichiometry['N2'])
    assert len(stoichiometry['H2(4)']) == 3933
    assert stoichiometry['H2(4)'][0] == 1.0
    assert stoichiometry['H2(4)'].count(1.0) == 134
    assert stoichiometry['H2(4)'].count(0) == 3783
    assert stoichiometry['H2(4)'].count(-1.0) == 16
    assert len(stoichiometry['CO(8)']) == 3933
    assert stoichiometry['CO(8)'][0] == 0
    assert stoichiometry['CO(8)'].count(1.0) == 272
    assert stoichiometry['CO(8)'].count(0) == 3484
    assert stoichiometry['CO(8)'].count(-1.0) == 169


def test_run_jsr():
    """Test getting ROPs from a JSR reactor"""
    model_path = os.path.join(DATA_BASE_PATH, 'models', 'HOCHO.yaml')
    gas = ct.Solution(model_path)
    profiles = flux.run_jsr(gas=gas, times=[0.5, 5], composition={'HOCHO(1)': 1.0}, T=1000, P=1)
    keys = list(profiles.keys())
    assert len(keys) == 2
    assert 0.01 < keys[0] < 2.0
    assert 2.0 < keys[1] < 20
    assert almost_equal(profiles[keys[0]]['P'], 1.0e5, places=2)
    assert almost_equal(profiles[keys[1]]['P'], 1.0e5, places=2)
    assert profiles[keys[0]]['T'] == 1000.0
    assert isinstance(profiles[keys[0]]['X'], dict)
    assert almost_equal(profiles[keys[0]]['X']['CO(8)'], 7.53994e-6)
    assert almost_equal(profiles[keys[0]]['X']['HOCHO(1)'], 0.99998358)
    assert len(profiles[keys[0]]['X']) == 152
    assert len(profiles[keys[0]]['ROPs']) == 152
    assert len(profiles[keys[0]]['ROPs']['H(3)']) == 842


def test_get_top_rops():
    """Test getting the top ROPs for the observables"""
    obervables = ['spc1', 'spc3']
    profiles = {0.5: {'ROPs': {'spc1': {'R1': 0.1, 'R2': 0.0002, 'R3': 0.3, 'R4': 0.00006, 'R5': -1e-16,
                                        'R6': -9, 'R7': 1e-16, 'R8': 0.6, 'R9': -0.9, 'R10': 2,
                                        'R11': 1e-22, 'R12': 0.5, 'R13': 0.4, 'R14': 0.33,
                                        'R15': -1e-10, 'R16': -0.008, 'R17': 1e-5, },
                               'spc2': {'R1': 0.1, 'R2': 0.004, 'R3': 0.3, 'R4': 6e-5, 'R5': -1e-11},
                               'spc3': {'R1': 0.1, 'R2': 0.0002, 'R3': 0.3, 'R4': 6e-5, 'R5': -1e-9}}},
                5: {'ROPs': {'spc1': {'R1': 0.5, 'R2': 0.008, 'R3': 0.3, 'R4': 6e-5, 'R5': -1e-16,
                                      'R6': -6e-5, 'R7': -4, 'R8': 0.6, 'R9': -0.9, 'R10': 2,
                                      'R11': 1e-9, 'R12': 0.5, 'R13': 0.4, 'R14': 0.33,
                                      'R15': -1e-8, 'R16': -0.008, 'R17': 0.00001, },
                             'spc2': {'R1': 0.1, 'R2': 4e-3, 'R3': 0.3, 'R4': 7e-5, 'R5': -1e-7},
                             'spc3': {'R1': 0.1, 'R2': 6e-4, 'R3': 0.3, 'R4': 6e-5, 'R5': -1e-7}}}}
    top_rops = flux.get_top_rops(profiles=profiles, observables=obervables)
    assert top_rops == {0.5: {'spc1': {'R6': -9, 'R10': 2, 'R9': -0.9, 'R8': 0.6, 'R12': 0.5,
                                       'R13': 0.4, 'R14': 0.33, 'R3': 0.3, 'R1': 0.1},
                              'spc3': {'R3': 0.3, 'R1': 0.1}},
                        5: {'spc1': {'R7': -4, 'R10': 2, 'R9': -0.9, 'R8': 0.6, 'R1': 0.5,
                                     'R12': 0.5, 'R13': 0.4, 'R14': 0.33, 'R3': 0.3, 'R2': 0.008},
                            'spc3': {'R3': 0.3, 'R1': 0.1, 'R2': 0.0006}}}


def test_generate_top_rop_bar_figs_1():
    """Test generating ROP bar figures for a dummy profile"""
    profiles = {0.5: {'ROPs': {'spc1': {'R1': 0.1, 'R2': 0.0002, 'R3': 0.3, 'R4': 0.00006, 'R5': -1e-16,
                                        'R6': -9, 'R7': 1e-16, 'R8': 0.6, 'R9': -0.9, 'R10': 2,
                                        'R11': 1e-22, 'R12': 0.5, 'R13': 0.4, 'R14': 0.33,
                                        'R15': -1e-10, 'R16': -0.008, 'R17': 1e-5, },
                               'spc2': {'R1': 0.1, 'R2': 0.004, 'R3': 0.3, 'R4': 6e-5, 'R5': -1e-11},
                               'spc3': {'HO2(5) + CO(7) <=> OH(4) + CO2(8)': 0.1, 'OH(4) + CO(7) <=> HOCO(9)': 0.0002,
                                        'OH(4) + CO(7) <=> OCHO(10)': 0.3, 'CH2O(11) (+M) <=> H2(3) + CO(7) (+M)': 6e-5,
                                        'H(2) + CH2O(11) <=> H(2) + H2(3) + CO(7)': -1e-9}}},
                5: {'ROPs': {'spc1': {'R1': 0.5, 'R2': 0.008, 'R3': 0.3, 'R4': 6e-5, 'R5': -1e-16,
                                      'R6': -6e-5, 'R7': -4, 'R8': 0.6, 'R9': -0.9, 'R10': 2,
                                      'R11': 1e-9, 'R12': 0.5, 'R13': 0.4, 'R14': 0.33,
                                      'R15': -1e-8, 'R16': -0.008, 'R17': 0.00001, },
                             'spc2': {'R1': 0.1, 'R2': 4e-3, 'R3': 0.3, 'R4': 7e-5, 'R5': -1e-7},
                             'spc3': {'R1': 0.1, 'R2': 6e-4, 'R3': 0.3, 'R4': 6e-5, 'R5': -1e-7}}}}
    observables = ['spc1', 'spc3']
    folder_path = os.path.join(SCRATCH_BASE_PATH, 'test_generate_top_rop_bar_figs_dummy')
    flux.generate_top_rop_bar_figs(profiles=profiles, observables=observables, folder_path=folder_path)
    for fig in ['spc1_0.5_s.png', 'spc1_5_s.png', 'spc3_0.5_s.png', 'spc3_5_s.png']:
        assert os.path.isfile(os.path.join(folder_path, 'bar_ROPs', fig))


def test_generate_top_rop_bar_figs_2():
    """Test generating ROP bar figures for formic acid pyrolysis"""
    observables = ['HOCHO(1)', 'CO2(9)']
    times = [0.001]
    profiles = flux.get_profiles_from_simulation(model_path=os.path.join(DATA_BASE_PATH, 'models', 'HOCHO.yaml'),
                                                 reactor_type='JSR',
                                                 times=times,
                                                 composition={'HOCHO(1)': 1.0},
                                                 T=1000,
                                                 P=1,
                                                 )
    folder_path = os.path.join(SCRATCH_BASE_PATH, 'test_generate_top_rop_bar_figs_HOCHO_pyrolysis')
    flux.generate_top_rop_bar_figs(profiles=profiles, observables=observables, folder_path=folder_path)
    for fig in ['HOCHO(1)_0.001_s.png', 'CO2(9)_0.001_s.png']:
        assert os.path.isfile(os.path.join(folder_path, 'bar_ROPs', fig))


def test_create_digraph_HOCHO_pyrolysis():
    """Test creating a flux diagram of formic acid pyrolysis"""
    observables = ['HOCHO(1)', 'CO2(9)']
    times = [0.001, 0.5]
    profiles = flux.get_profiles_from_simulation(model_path=os.path.join(DATA_BASE_PATH, 'models', 'HOCHO.yaml'),
                                                 reactor_type='JSR',
                                                 times=times,
                                                 composition={'HOCHO(1)': 1.0},
                                                 T=1000,
                                                 P=1,
                                                 )
    for i, (time, profile) in enumerate(profiles.items()):
        flux_graph, nodes_to_explore, min_rop, max_rop = flux.get_flux_graph(profile=profile, observables=observables)
        folder_path = os.path.join(SCRATCH_BASE_PATH, 'test_create_digraph_HOCHO_pyrolysis')
        flux.create_digraph(flux_graph=flux_graph,
                            profile=profile,
                            observables=observables,
                            nodes_to_explore=nodes_to_explore,
                            time=time,
                            min_rop=min_rop,
                            max_rop=max_rop,
                            folder_path=folder_path,
                            )
        assert os.path.isfile(os.path.join(folder_path, f'flux_diagram_{times[i]}_s.dot'))
        assert os.path.isfile(os.path.join(folder_path, f'flux_diagram_{times[i]}_s.png'))


def test_create_digraph_NH3():
    """Test creating a flux diagram of ammonia oxidation"""
    observables = ['NH3(1)']
    profiles = flux.get_profiles_from_simulation(model_path=os.path.join(DATA_BASE_PATH, 'models', 'NH3.yaml'),
                                                 reactor_type='JSR',
                                                 times=[1e-2],
                                                 composition={'NH3(1)': 1.0, 'O2(3)': 0.75, 'N2(2)': 5.65},
                                                 T=1500,
                                                 P=10,
                                                 )
    for time, profile in profiles.items():
        flux_graph, nodes_to_explore, min_rop, max_rop = flux.get_flux_graph(profile=profile, observables=observables)
    folder_path = os.path.join(SCRATCH_BASE_PATH, 'test_create_digraph_NH3')
    flux.create_digraph(flux_graph=flux_graph,
                        profile=profile,
                        observables=observables,
                        nodes_to_explore=nodes_to_explore,
                        time=1e-2,
                        min_rop=min_rop,
                        max_rop=max_rop,
                        folder_path=folder_path,
                        display_concentrations=False,
                        display_flux_ratio=False,
                        display_r_n_p=True,
                        )
    assert os.path.isfile(os.path.join(folder_path, 'flux_diagram_0.01_s.png'))


def test_create_digraph_N2H4():
    """Test creating a flux diagram of ammonia oxidation"""
    observables = ['H4N2(1)']
    profiles = flux.get_profiles_from_simulation(model_path=os.path.join(DATA_BASE_PATH, 'models', 'N2H4.yaml'),
                                                 reactor_type='JSR',
                                                 times=[0.05],
                                                 composition={'H4N2(1)': 1},
                                                 T=700,
                                                 P=0.001,
                                                 )
    for time, profile in profiles.items():
        flux_graph, nodes_to_explore, min_rop, max_rop = flux.get_flux_graph(profile=profile, observables=observables)
    folder_path = os.path.join(SCRATCH_BASE_PATH, 'test_create_digraph_N2H4')
    flux.create_digraph(flux_graph=flux_graph,
                        profile=profile,
                        observables=observables,
                        nodes_to_explore=nodes_to_explore,
                        time=0.005,
                        min_rop=min_rop,
                        max_rop=max_rop,
                        folder_path=folder_path,
                        display_concentrations=True,
                        display_flux_ratio=True,
                        display_r_n_p=True,
                        )
    assert os.path.isfile(os.path.join(folder_path, 'flux_diagram_0.005_s.png'))


def test_get_width():
    """Test getting the node or edge width for a given concentration or ROP value."""
    for x in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
        assert almost_equal(flux.get_width(x=x, x_min=0, x_max=10, log_scale=False), 0.38 * x + 0.2)
    assert almost_equal(flux.get_width(x=2, x_min=2, x_max=10, log_scale=False), 0.2)
    assert almost_equal(flux.get_width(x=6, x_min=2, x_max=10, log_scale=False), 2.1)
    assert almost_equal(flux.get_width(x=10, x_min=2, x_max=10, log_scale=False), 4)

    for x in [1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2]:
        assert almost_equal(flux.get_width(x=x, x_min=1e-12, x_max=1e-2, log_scale=True), 0.165 * np.log(x) + 4.76, places=2)
    assert almost_equal(flux.get_width(x=1e-8, x_min=1e-8, x_max=1e-2, log_scale=True), 0.2)
    assert almost_equal(flux.get_width(x=1e-5, x_min=1e-8, x_max=1e-2, log_scale=True), 2.1, places=1)
    assert almost_equal(flux.get_width(x=1e-2, x_min=1e-8, x_max=1e-2, log_scale=True), 4)

    assert almost_equal(flux.get_width(x=-1, x_min=2e-18, x_max=1, log_scale=True), 4)
    assert almost_equal(flux.get_width(x=-0.089, x_min=2e-18, x_max=1, log_scale=True), 3.77, places=2)


def test_get_rxn_in_relevant_direction():
    """Test getting a reaction string in the direction where the species in one of the reactants."""
    rxn = 'HCO(12) (+M) <=> H(2) + CO(7) (+M)'
    spc = 'HCO(12)'
    new_rxn = flux.get_rxn_in_relevant_direction(rxn=rxn, spc=spc)
    assert new_rxn == 'HCO(12) (+M) <=> H(2) + CO(7) (+M)'
    spc = 'H(2)'
    new_rxn = flux.get_rxn_in_relevant_direction(rxn=rxn, spc=spc)
    assert new_rxn == 'H(2) + CO(7) (+M) <=> HCO(12) (+M)'

    rxn = 'HO2(5) + CO2(8) <=> O(18) + CHO3(765)'
    spc = 'CHO3(765)'
    new_rxn = flux.get_rxn_in_relevant_direction(rxn=rxn, spc=spc)
    assert new_rxn == 'O(18) + CHO3(765) <=> HO2(5) + CO2(8)'
    spc = 'CO2(8)'
    new_rxn = flux.get_rxn_in_relevant_direction(rxn=rxn, spc=spc)
    assert new_rxn == 'HO2(5) + CO2(8) <=> O(18) + CHO3(765)'


def test_get_node():
    """Test getting a node from a pydot graph or creating one."""
    graph = pydot.Dot(graph_type='digraph')
    node_a = pydot.Node(name='Node A', style='filled', fillcolor='#DCE5F4')
    node_b = pydot.Node(name='Node B')
    node_c = pydot.Node(name='Node C')
    node_d = pydot.Node(name='Node D')
    graph.add_node(node_a)
    graph.add_node(node_b)
    graph.add_node(node_c)
    graph.add_node(node_d)
    graph.add_edge(pydot.Edge(node_a, node_b))
    graph.add_edge(pydot.Edge(node_b, node_c))
    graph.add_edge(pydot.Edge(node_c, node_d))
    nodes = {'Node A': node_a,
             'Node B': node_b,
             'Node C': node_c,
             'Node D': node_d}
    for name, node in nodes.items():
        assert flux.get_node(graph=graph, label=name, nodes=nodes) == node
    assert isinstance(flux.get_node(graph=graph, label='Node E', nodes=nodes), pydot.Node)
    node_f = flux.get_node(graph=graph, label='Node F', nodes=nodes, observables=['Node F'])
    assert isinstance(node_f, pydot.Node)
    assert node_f.to_string() == '"Node F" [fillcolor="#DCE5F4", fontsize=8, style=filled, xlabel=""];'


def test_get_flux_graph():
    """Test getting a normalized flux profile and generating a flux graph."""
    profiles = flux.get_profiles_from_simulation(model_path=os.path.join(DATA_BASE_PATH, 'models', 'HOCHO.yaml'),
                                                 reactor_type='JSR',
                                                 times=[0.001, 5],
                                                 composition={'HOCHO(1)': 1.0},
                                                 T=1000,
                                                 P=1,
                                                 )
    for i, (time, profile) in enumerate(profiles.items()):
        flux_graph, nodes_to_explore, min_rop, max_rop = flux.get_flux_graph(profile=profile, observables=['HOCHO(1)', 'CO2(9)'])
        if i == 0:
            assert nodes_to_explore == {'CO(8)', 'H2O(17)', 'H2(4)'}
            assert almost_equal(min_rop, 3.27e-21, ratio=100)
            assert almost_equal(max_rop, 0.001326, places=5)
            assert list(flux_graph.keys()) == ['CO2(9)', 'HOCHO(1)', 'H2(4)', 'H2O(17)', 'CO(8)']
            assert flux_graph['CO2(9)']['2 HOCO(10) <=> CO2(9) + HOCHO(1)'][0] == ['2 HOCO(10)']
            assert almost_equal(flux_graph['CO2(9)']['2 HOCO(10) <=> CO2(9) + HOCHO(1)'][1], -1.73922e-15, ratio=100)
        all_rop_values = list()
        for rxn_dict in flux_graph.values():
            for rxn, value in rxn_dict.items():
                all_rop_values.append(value[1])
        assert any(abs(v) == 1 for v in all_rop_values)

    profiles = flux.get_profiles_from_simulation(model_path=os.path.join(DATA_BASE_PATH, 'models', 'N2H4.yaml'),
                                                 reactor_type='JSR',
                                                 times=[0.05],
                                                 composition={'H4N2(1)': 1.0},
                                                 T=1200,
                                                 P=1,
                                                 )
    for i, (time, profile) in enumerate(profiles.items()):
        assert almost_equal(profile['P'], 1e5, places=2)  # Pa
        assert almost_equal(profile['T'], 1200)
        assert len(profile['ROPs']) == 36
        flux_graph, nodes_to_explore, min_rop, max_rop = flux.get_flux_graph(profile=profile, observables=['H4N2(1)'])
        if i == 0:
            assert nodes_to_explore == {'H2(4)', 'H3N2(6)', 'ammonia(9)', 'H(3)', 'H2N2(7)', '2 NH2(5)', 'HN2(10)', 'N2(2)', 'NH2(5)'}
            assert almost_equal(min_rop, 3.27e-21, ratio=100)
            assert almost_equal(max_rop, 20.3659, places=3)
            assert list(flux_graph.keys()) == ['H4N2(1)', 'NH2(5)', 'HN2(10)', 'H(3)', 'H2N2(7)', 'H2(4)', 'ammonia(9)', 'H3N2(6)']
            assert flux_graph['H4N2(1)']['H4N2(1) + NH2(5) <=> H3N2(6) + ammonia(9)'][0] == ['H3N2(6)', 'ammonia(9)']
            assert almost_equal(flux_graph['H4N2(1)']['H4N2(1) + NH2(5) <=> H3N2(6) + ammonia(9)'][1], -1.0)
        all_rop_values = list()
        for rxn_dict in flux_graph.values():
            for rxn, value in rxn_dict.items():
                all_rop_values.append(value[1])
        assert any(abs(v) == 1 for v in all_rop_values)


def test_get_opposite_rxn_species():
    """Test getting the species on the opposite side of a reaction."""
    rxn = 'HCO(12) (+M) <=> H(2) + CO(7) (+M)'
    assert flux.get_opposite_rxn_species(rxn=rxn, spc='HCO(12)') == ['H(2)', 'CO(7)']
    assert flux.get_opposite_rxn_species(rxn=rxn, spc='H(2)') == ['HCO(12)']
    assert flux.get_opposite_rxn_species(rxn=rxn, spc='CO(7)') == ['HCO(12)']

    rxn = 'HO2(5) + CO2(8) <=> O(18) + CHO3(765)'
    assert flux.get_opposite_rxn_species(rxn=rxn, spc='HO2(5)') == ['O(18)', 'CHO3(765)']
    assert flux.get_opposite_rxn_species(rxn=rxn, spc='O(18)') == ['HO2(5)', 'CO2(8)']
    assert flux.get_opposite_rxn_species(rxn=rxn, spc='CHO3(765)') == ['HO2(5)', 'CO2(8)']

    rxn = 'H(3) + HO2(6) + M <=> H2O2(18) + M'
    assert flux.get_opposite_rxn_species(rxn=rxn, spc='HO2(6)') == ['H2O2(18)']
    assert flux.get_opposite_rxn_species(rxn=rxn, spc='H(3)') == ['H2O2(18)']
    assert flux.get_opposite_rxn_species(rxn=rxn, spc='H2O2(18)') == ['H(3)', 'HO2(6)']


def test_get_other_reactants_and_products():
    """Test getting the reactants and products other than a given species."""
    rxn = 'HCO(12) (+M) <=> H(2) + CO(7) (+M)'
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['HCO(12)', 'H(2)'])
    assert rs == ''
    assert ps == '- CO(7)'
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['CO(7)', 'HCO(12)'])
    assert rs == '+ H(2)'
    assert ps == ''
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['H(2)', 'HCO(12)'])
    assert rs == '+ CO(7)'
    assert ps == ''

    rxn = 'OH(4) + HOCO(9) <=> H2O(13) + CO2(8)'
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['OH(4)', 'H2O(13)'])
    assert rs == '+ HOCO(9)'
    assert ps == '- CO2(8)'
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['CO2(8)', 'HOCO(9)'])
    assert rs == '+ H2O(13)'
    assert ps == '- OH(4)'

    rxn = 'H(2) + H(2) + M <=> H2(3) + M'
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['H(2)', 'H2(3)'])
    assert rs == ''
    assert ps == ''
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['H2(3)', 'H(2)'])
    assert rs == ''
    assert ps == ''

    rxn = 'H(2) + H(2) + H2(3) <=> H2(3) + H2(3)'
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['H(2)', 'H2(3)'])
    assert rs == ''
    assert ps == ''
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['H2(3)', 'H(2)'])
    assert rs == ''
    assert ps == ''

    rxn = 'O(18) + HOCHO(1) => OH(4) + OH(4) + CO(7)'
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['O(18)', 'OH(4)'])
    assert rs == '+ HOCHO(1)'
    assert ps == '- CO(7)'
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['HOCHO(1)', 'CO(7)'])
    assert rs == '+ O(18)'
    assert ps == '- OH(4) - OH(4)'
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['OH(4)', 'HOCHO(1)'])
    assert rs == '+ CO(7)'
    assert ps == '- O(18)'
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['CO(7)', 'O(18)'])
    assert rs == '+ OH(4) + OH(4)'
    assert ps == '- HOCHO(1)'

    rxn = 'CO2(8) + HOCHO(1) <=> 2 HOCO(9)'
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['CO2(8)', 'HOCO(9)'])
    assert rs == '+ HOCHO(1)'
    assert ps == ''
    rs, ps = flux.get_other_reactants_and_products(rxn=rxn, spcs=['HOCO(9)', 'HOCHO(1)'])
    assert rs == ''
    assert ps == '- CO2(8)'


def test_unpack_stoichiometry():
    """Test unpacking stoichiometric coefficients from labels."""
    labels = ['2 A', 'B', '3 C']
    labels, multipliers = flux.unpack_stoichiometry(labels=labels)
    assert labels == ['A', 'B', 'C']
    assert multipliers == [2, 1, 3]

    labels = ['2 HOCO(9)']
    labels, multipliers = flux.unpack_stoichiometry(labels=labels)
    assert labels == ['HOCO(9)']
    assert multipliers == [2]


def test_continue_exploring():
    """Test determining whether to continue exploring a species based on its ROP"""
    rops = {'R1': 1.0,
            'R2': 1e-6,
            'R3': 2.5,
            'R4': 1e-18,
            'R5': -1e-24}
    assert flux.continue_exploring(rops) is False
    rops = {'R1': 1.0,
            'R2': 1e-6,
            'R3': -2.5,
            'R4': 1e-18,
            'R5': -1e-24}
    assert flux.continue_exploring(rops) is True


def test_almost_equal():
    """Test the almost_equal help function"""
    assert almost_equal(1, 1)
    assert not almost_equal(1, 2)
    assert almost_equal(0.00010, 0.00013)
    assert almost_equal(0.0010, 0.0013, places=3)
    assert not almost_equal(0.0010, 0.0013)
    assert almost_equal(1e-15, 1e-22)
    assert not almost_equal(1e-15, 1e-22, places=16)
    assert almost_equal(1, 2, ratio=10)
    assert almost_equal(1e-5, 1e-6, ratio=100)
    assert almost_equal(1e-6, 1e-5, ratio=100)


def teardown_module():
    """teardown any state that was previously setup with a setup_module method."""
    shutil.rmtree(SCRATCH_BASE_PATH, ignore_errors=True)
