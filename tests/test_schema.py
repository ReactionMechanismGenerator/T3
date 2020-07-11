#!/usr/bin/env python3
# encoding: utf-8

"""
schema test module
"""

import pytest
from pydantic import ValidationError

from t3.schema import (T3Options,
                       T3Sensitivity,
                       T3Uncertainty,
                       RMGDatabase,
                       RMGSpecies,
                       RMGReactor,
                       RMGModel,
                       RMGOptions,
                       RMGPDep,
                       RMGSpeciesConstraints,
                       QM)

# define a long quote of length 444 characters to test constraints on string length
quote = 'Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et ' \
        'dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ' \
        'ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu ' \
        'fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt ' \
        'mollit anim id est laborum'


def test_t3_options_schema():
    """Test creating an instance of T3Options"""
    t3_options = T3Options(flux_adapter='RMG',
                           profiles_adapter='RMG',
                           collision_violators_thermo=False,
                           collision_violators_rates=False,
                           all_core_species=False,
                           all_core_reactions=False,
                           fit_missing_GAV=False,
                           max_T3_iterations=5,
                           max_RMG_exceptions_allowed=5,
                           max_RMG_walltime='00:02:00:00',
                           max_T3_walltime='01:00:00:00',
                           library_name='T3_library'
                           )
    assert t3_options.flux_adapter == 'RMG'
    assert t3_options.profiles_adapter == 'RMG'
    assert t3_options.collision_violators_thermo is False
    assert t3_options.collision_violators_rates is False
    assert t3_options.all_core_species is False
    assert t3_options.all_core_reactions is False
    assert t3_options.fit_missing_GAV is False
    assert t3_options.max_T3_iterations == 5
    assert t3_options.max_RMG_exceptions_allowed == 5
    assert t3_options.max_RMG_walltime == '00:02:00:00'
    assert t3_options.max_T3_walltime == '01:00:00:00'
    assert t3_options.library_name == 'T3_library'

    with pytest.raises(ValidationError):
        # check that flux_adapter is constrained to at most 255 characters
        T3Options(flux_adapter=quote)

    with pytest.raises(ValidationError):
        # check that profiles_adapter is constrained to at most 255 characters
        T3Options(profiles_adapter=quote)

    with pytest.raises(ValidationError):
        # check that max_T3_iterations is > 0
        T3Options(max_T3_iterations=0)

    with pytest.raises(ValidationError):
        # check that max_RMG_exceptions_allowed is >= 0
        T3Options(max_RMG_exceptions_allowed=-1)

    with pytest.raises(ValidationError):
        # check that max_RMG_walltime is in the expected format
        T3Options(max_RMG_walltime='05-02:00:00')

    with pytest.raises(ValidationError):
        # check that max_T3_walltime is in the expected format
        T3Options(max_T3_walltime='05-02:00:00')

    with pytest.raises(ValidationError):
        # check that library_name is constrained to at most 255 characters
        T3Options(library_name=quote)


def test_t3_sensitivity_schema():
    """Test creating an instance of T3Sensitivity"""
    t3_sensitivity = T3Sensitivity(adapter=None,
                                   atol=1e-6,
                                   rtol=1e-4,
                                   global_observables=None,
                                   SA_threshold=0.01,
                                   pdep_SA_threshold=0.001,
                                   ME_methods=['CSE', 'MSC'],
                                   top_SA_species=10,
                                   top_SA_reactions=10
                                   )
    assert t3_sensitivity.adapter is None
    assert t3_sensitivity.atol == 1e-6
    assert t3_sensitivity.rtol == 1e-4
    assert t3_sensitivity.global_observables is None
    assert t3_sensitivity.SA_threshold == 0.01
    assert t3_sensitivity.pdep_SA_threshold == 0.001
    assert t3_sensitivity.ME_methods == ['CSE', 'MSC']
    assert t3_sensitivity.top_SA_species == 10
    assert t3_sensitivity.top_SA_reactions == 10

    with pytest.raises(ValidationError):
        # check that adapter is constrained to at most 255 characters
        T3Sensitivity(adapter=quote)

    with pytest.raises(ValidationError):
        # check that atol is constrained to > 0
        T3Sensitivity(atol=0)

    with pytest.raises(ValidationError):
        # check that atol is constrained to < 1e-1
        T3Sensitivity(atol=1)

    with pytest.raises(ValidationError):
        # check that rtol is constrained to > 0
        T3Sensitivity(rtol=0)

    with pytest.raises(ValidationError):
        # check that rtol is constrained to < 1e-1
        T3Sensitivity(rtol=1)

    with pytest.raises(ValidationError):
        # check that entries in global_observables have at least 2 characters
        T3Sensitivity(global_observables='a')

    with pytest.raises(ValidationError):
        # check that entries in global_observables have at most 3 characters
        T3Sensitivity(global_observables='abcd')

    with pytest.raises(ValidationError):
        # check that SA_threshold is constrained to > 0
        T3Sensitivity(SA_threshold=0)

    with pytest.raises(ValidationError):
        # check that SA_threshold is constrained to < 0.5
        T3Sensitivity(SA_threshold=1)

    with pytest.raises(ValidationError):
        # check that pdep_SA_threshold is constrained to > 0
        T3Sensitivity(pdep_SA_threshold=0)

    with pytest.raises(ValidationError):
        # check that pdep_SA_threshold is constrained to < 0.5
        T3Sensitivity(pdep_SA_threshold=1)

    with pytest.raises(ValidationError):
        # check that entries in ME_methods have at least 2 characters
        T3Sensitivity(ME_methods='a')

    with pytest.raises(ValidationError):
        # check that entries in ME_methods have at most 3 characters
        T3Sensitivity(ME_methods='abcd')

    with pytest.raises(ValidationError):
        # check that top_SA_species is constrained to >= 0
        T3Sensitivity(top_SA_species=-1)

    with pytest.raises(ValidationError):
        # check that top_SA_reactions is constrained to >= 0
        T3Sensitivity(top_SA_reactions=-1)


def test_t3_uncertainty_schema():
    """Test creating an instance of T3Uncertainty"""
    t3_uncertainty = T3Uncertainty(adapter=None,
                                   local_analysis=False,
                                   global_analysis=False,
                                   correlated=True,
                                   local_number=10,
                                   global_number=5,
                                   termination_time=None,
                                   PCE_run_time=1800,
                                   PCE_error_tolerance=None,
                                   PCE_max_evals=None,
                                   logx=False
                                   )
    assert t3_uncertainty.adapter is None
    assert t3_uncertainty.local_analysis is False
    assert t3_uncertainty.global_analysis is False
    assert t3_uncertainty.correlated is True
    assert t3_uncertainty.local_number == 10
    assert t3_uncertainty.global_number == 5
    assert t3_uncertainty.termination_time is None
    assert t3_uncertainty.PCE_run_time == 1800
    assert t3_uncertainty.PCE_error_tolerance is None
    assert t3_uncertainty.PCE_max_evals is None
    assert t3_uncertainty.logx is False

    with pytest.raises(ValidationError):
        # check that adapter is constrained to at most 255 characters
        T3Uncertainty(adapter=quote)

    with pytest.raises(ValidationError):
        # check that local_number is constrained to > 0
        T3Uncertainty(local_number=0)

    with pytest.raises(ValidationError):
        # check that global_number is constrained to > 0
        T3Uncertainty(global_number=0)

    with pytest.raises(ValidationError):
        # check that termination_time is in the expected format
        T3Uncertainty(termination_time='05-02:00:00')

    with pytest.raises(ValidationError):
        # check that PCE_run_time is constrained to > 0
        T3Uncertainty(PCE_run_time=0)

    with pytest.raises(ValidationError):
        # check that PCE_max_evals is constrained to > 0
        T3Uncertainty(PCE_max_evals=0)


def test_rmg_database_schema():
    """Test creating an instance of RMGDatabase"""
    rmg_database = RMGDatabase(thermo_libraries=['BurkeH2O2', 'DFT_QCI_thermo', 'primaryThermoLibrary', 'CBS_QB3_1dHR'],
                               kinetics_libraries=['BurkeH2O2inN2', 'NOx2018', 'Klippenstein_Glarborg2016'],
                               transport_libraries=['PrimaryTransportLibrary', 'OneDMinN2', 'NOx2018', 'GRI-Mech'],
                               seed_mechanisms=list(),
                               kinetics_depositories='default',
                               kinetics_families='default',
                               kinetics_estimator='rate rules'
                               )
    assert rmg_database.thermo_libraries == ['BurkeH2O2', 'DFT_QCI_thermo', 'primaryThermoLibrary', 'CBS_QB3_1dHR']
    assert rmg_database.kinetics_libraries == ['BurkeH2O2inN2', 'NOx2018', 'Klippenstein_Glarborg2016']
    assert rmg_database.transport_libraries == ['PrimaryTransportLibrary', 'OneDMinN2', 'NOx2018', 'GRI-Mech']
    assert rmg_database.seed_mechanisms == list()
    assert rmg_database.kinetics_depositories == 'default'
    assert rmg_database.kinetics_families == 'default'
    assert rmg_database.kinetics_estimator == 'rate rules'


def test_rmg_species_schema():
    """Test creating an instance of RMGSpecies"""
    rmg_species = RMGSpecies(label='N2',
                             concentration=1,
                             smiles='N#N',
                             inchi='1S/N2/c1-2',
                             adjlist="""
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
""",
                             reactive=False,
                             observable=False,
                             SA_observable=False,
                             UA_observable=False,
                             constant=True,
                             balance=True,
                             solvent=False,
                             seed_all_rads=['radical', 'alkoxyl', 'peroxyl'],
                             )
    assert rmg_species.label == 'N2'
    assert rmg_species.concentration == 1
    assert rmg_species.smiles == 'N#N'
    assert rmg_species.inchi == '1S/N2/c1-2'
    # adjlist must be in the same column as adjlist above
    assert rmg_species.adjlist == """
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
"""
    assert rmg_species.reactive is False
    assert rmg_species.observable is False
    assert rmg_species.SA_observable is False
    assert rmg_species.UA_observable is False
    assert rmg_species.constant is True
    assert rmg_species.balance is True
    assert rmg_species.solvent is False
    assert rmg_species.seed_all_rads == ['radical', 'alkoxyl', 'peroxyl']

    with pytest.raises(ValidationError):
        # check that concentration is constrained to >= 0
        RMGSpecies(concentration=-1)

    with pytest.raises(ValidationError):
        # test RadicalTypeEnum
        RMGSpecies(label='N2',
                   concentration=1,
                   smiles='N#N',
                   inchi='1S/N2/c1-2',
                   adjlist="""
        1 N u0 p1 c0 {2,T}
        2 N u0 p1 c0 {1,T}
        """,
                   reactive=False,
                   observable=False,
                   SA_observable=False,
                   UA_observable=False,
                   constant=True,
                   balance=True,
                   solvent=False,
                   seed_all_rads=['radical', 'non-supported'],
                   )


def test_rmg_reactors_schema():
    """Test creating an instance of RMGReactor"""
    rmg_reactor = RMGReactor(type='gas batch constant T P',
                             T=[800, 1750],
                             P=[1e0, 1e1],
                             termination_conversion={'ethane': 0.2},
                             termination_time=[5, 's'],
                             termination_rate_ratio=0.01,
                             conditions_per_iteration=12
                             )
    assert rmg_reactor.type == 'gas batch constant T P'
    assert rmg_reactor.T == [800, 1750]
    assert rmg_reactor.P == [1e0, 1e1]
    assert rmg_reactor.termination_conversion == {'ethane': 0.2}
    assert rmg_reactor.termination_time == (5, 's')
    assert rmg_reactor.termination_rate_ratio == 0.01
    assert rmg_reactor.conditions_per_iteration == 12

    # test 'micro-s' units in termination_time
    rmg_reactor = RMGReactor(type='gas batch constant T P',
                             T=[800, 1750],
                             P=[1e0, 1e1],
                             termination_conversion={'ethane': 0.2},
                             termination_time=[5, 'micro-s'],
                             termination_rate_ratio=0.01,
                             conditions_per_iteration=12
                             )
    assert rmg_reactor.termination_time == (5000, 'ms')

    # test 'hrs' units in termination_time
    rmg_reactor = RMGReactor(type='gas batch constant T P',
                             T=[800, 1750],
                             P=[1e0, 1e1],
                             termination_conversion={'ethane': 0.2},
                             termination_time=[5, 'hrs'],
                             termination_rate_ratio=0.01,
                             conditions_per_iteration=12
                             )
    assert rmg_reactor.termination_time == (5, 'hours')

    with pytest.raises(ValidationError):
        # check that scalar T is constrained to > 0
        RMGSpecies(T=0)

    with pytest.raises(ValidationError):
        # check that elements in list T are constrained to > 0
        RMGSpecies(T=[0, 800])

    with pytest.raises(ValidationError):
        # check that scalar P is constrained to > 0
        RMGSpecies(P=0)

    with pytest.raises(ValidationError):
        # check that elements in list P is constrained to > 0
        RMGSpecies(P=[0, 800])

    with pytest.raises(ValidationError):
        # check that values for termination_conversion are constrained to > 0
        RMGSpecies(termination_conversion={'ethane': 0})

    with pytest.raises(ValidationError):
        # check that values for termination_conversion are constrained to < 1
        RMGSpecies(termination_conversion={'ethane': 1})

    with pytest.raises(ValidationError):
        # check wrong units for termination_time
        RMGSpecies(termination_time=[5, 'wrong'])

    with pytest.raises(ValidationError):
        # check that termination_time is constrained to > 0
        RMGSpecies(termination_time=[0, 's'])

    with pytest.raises(ValidationError):
        # check that termination_rate_ratio is constrained to > 0
        RMGSpecies(termination_rate_ratio=0)

    with pytest.raises(ValidationError):
        # check that termination_rate_ratio is constrained to < 1
        RMGSpecies(termination_rate_ratio=1)

    with pytest.raises(ValidationError):
        # check that conditions_per_iteration is constrained to > 0
        RMGSpecies(conditions_per_iteration=0)


def test_rmg_model_schema():
    """Test creating an instance of RMGModel"""
    rmg_model = RMGModel(core_tolerance=[0.05, 0.01],
                         atol=1e-16,
                         rtol=1e-8,
                         filter_reactions=False,
                         filter_threshold=1e8,
                         tolerance_interrupt_simulation=[0.05, 0.01],
                         min_core_size_for_prune=50,
                         min_species_exist_iterations_for_prune=2,
                         tolerance_keep_in_edge=0.02,
                         maximum_edge_species=1e6,
                         tolerance_thermo_keep_species_in_edge=100,
                         max_num_species=None,
                         tolerance_move_edge_reaction_to_core=None,
                         tolerance_move_edge_reaction_to_core_interrupt=None,
                         dynamics_time_scale=(0.0, 'sec'),
                         max_num_objs_per_iter=1,
                         terminate_at_max_objects=False,
                         ignore_overall_flux_criterion=None,
                         tolerance_branch_reaction_to_core=0.001,
                         branching_index=0.5,
                         branching_ratio_max=1
                         )
    assert rmg_model.core_tolerance == [0.05, 0.01]
    assert rmg_model.atol == 1e-16
    assert rmg_model.rtol == 1e-8
    assert rmg_model.filter_reactions is False
    assert rmg_model.filter_threshold == 1e8
    assert rmg_model.tolerance_interrupt_simulation == [0.05, 0.01]
    assert rmg_model.min_core_size_for_prune == 50
    assert rmg_model.min_species_exist_iterations_for_prune == 2
    assert rmg_model.tolerance_keep_in_edge == 0.02
    assert rmg_model.maximum_edge_species == 1e6
    assert rmg_model.tolerance_thermo_keep_species_in_edge == 100
    assert rmg_model.max_num_species is None
    assert rmg_model.tolerance_move_edge_reaction_to_core is None
    assert rmg_model.tolerance_move_edge_reaction_to_core_interrupt is None
    assert rmg_model.dynamics_time_scale == (0.0, 'sec')
    assert rmg_model.max_num_objs_per_iter == 1
    assert rmg_model.terminate_at_max_objects is False
    assert rmg_model.ignore_overall_flux_criterion is None
    assert rmg_model.tolerance_branch_reaction_to_core == 0.001
    assert rmg_model.branching_index == 0.5
    assert rmg_model.branching_ratio_max == 1

    with pytest.raises(ValidationError):
        # check that scalar core_tolerance is constrained to > 0
        RMGSpecies(core_tolerance=0)

    with pytest.raises(ValidationError):
        # check that scalar core_tolerance is constrained to < 1
        RMGSpecies(core_tolerance=1)

    with pytest.raises(ValidationError):
        # check that elements in list core_tolerance are constrained to > 0
        RMGSpecies(core_tolerance=[0, 0.01])

    with pytest.raises(ValidationError):
        # check that elements in list core_tolerance are constrained to < 1
        RMGSpecies(core_tolerance=[0.05, 1])

    with pytest.raises(ValidationError):
        # check that atol is constrained to > 0
        RMGSpecies(atol=0)

    with pytest.raises(ValidationError):
        # check that atol is constrained to < 1e-1
        RMGSpecies(atol=1)

    with pytest.raises(ValidationError):
        # check that rtol is constrained to > 0
        RMGSpecies(rtol=0)

    with pytest.raises(ValidationError):
        # check that rtol is constrained to < 1e-1
        RMGSpecies(rtol=1)

    with pytest.raises(ValidationError):
        # check that filter_threshold is constrained to > 0
        RMGSpecies(filter_threshold=0)

    with pytest.raises(ValidationError):
        # check that scalar tolerance_interrupt_simulation is constrained to > 0
        RMGSpecies(tolerance_interrupt_simulation=0)

    with pytest.raises(ValidationError):
        # check that elements in list tolerance_interrupt_simulation are constrained to > 0
        RMGSpecies(tolerance_interrupt_simulation=[0, 0.01])

    with pytest.raises(ValidationError):
        # check that min_core_size_for_prune is constrained to > 0
        RMGSpecies(min_core_size_for_prune=0)

    with pytest.raises(ValidationError):
        # check that min_species_exist_iterations_for_prune is constrained to > 0
        RMGSpecies(min_species_exist_iterations_for_prune=0)

    with pytest.raises(ValidationError):
        # check that tolerance_keep_in_edge is constrained to > 0
        RMGSpecies(tolerance_keep_in_edge=0)

    with pytest.raises(ValidationError):
        # check that maximum_edge_species is constrained to > 0
        RMGSpecies(maximum_edge_species=0)

    with pytest.raises(ValidationError):
        # check that tolerance_thermo_keep_species_in_edge is constrained to > 0
        RMGSpecies(tolerance_thermo_keep_species_in_edge=0)

    with pytest.raises(ValidationError):
        # check that max_num_species is constrained to > 0
        RMGSpecies(max_num_species=0)

    with pytest.raises(ValidationError):
        # check that tolerance_move_edge_reaction_to_core is constrained to > 0
        RMGSpecies(tolerance_move_edge_reaction_to_core=0)

    with pytest.raises(ValidationError):
        # check that tolerance_move_edge_reaction_to_core_interrupt is constrained to > 0
        RMGSpecies(tolerance_move_edge_reaction_to_core_interrupt=0)

    with pytest.raises(ValidationError):
        # check that max_num_objs_per_iter is constrained to > 0
        RMGSpecies(max_num_objs_per_iter=0)

    with pytest.raises(ValidationError):
        # check that tolerance_branch_reaction_to_core is constrained to > 0
        RMGSpecies(tolerance_branch_reaction_to_core=0)

    with pytest.raises(ValidationError):
        # check that branching_index is constrained to > 0
        RMGSpecies(branching_index=0)

    with pytest.raises(ValidationError):
        # check that branching_ratio_max is constrained to > 0
        RMGSpecies(branching_ratio_max=0)


def test_rmg_options_schema():
    """Test creating an instance of RMGOptions"""
    rmg_options = RMGOptions(seed_name='Seed',
                             save_edge=True,
                             save_html=False,
                             generate_seed_each_iteration=True,
                             save_seed_to_database=False,
                             units='si',
                             generate_plots=False,
                             save_simulation_profiles=False,
                             verbose_comments=False,
                             keep_irreversible=False,
                             trimolecular_product_reversible=True,
                             save_seed_modulus=-1
                             )
    assert rmg_options.seed_name == 'Seed'
    assert rmg_options.save_edge is True
    assert rmg_options.save_html is False
    assert rmg_options.generate_seed_each_iteration is True
    assert rmg_options.save_seed_to_database is False
    assert rmg_options.units == 'si'
    assert rmg_options.generate_plots is False
    assert rmg_options.save_simulation_profiles is False
    assert rmg_options.verbose_comments is False
    assert rmg_options.keep_irreversible is False
    assert rmg_options.trimolecular_product_reversible is True
    assert rmg_options.save_seed_modulus == -1

    with pytest.raises(ValidationError):
        # check that save_seed_modulus is constrained to > -1
        RMGOptions(save_seed_modulus=-2)


def test_rmg_pdep_schema():
    """Test creating an instance of RMGPDep"""
    rmg_pdep = RMGPDep(method='MSC',
                       max_grain_size=2,
                       max_number_of_grains=250,
                       T=[300, 2500, 10],
                       P=[0.01, 100, 10],
                       interpolation='Chebyshev',
                       T_basis_set=6,
                       P_basis_set=4,
                       max_atoms=16
                       )
    assert rmg_pdep.method == 'MSC'
    assert rmg_pdep.max_grain_size == 2
    assert rmg_pdep.max_number_of_grains == 250
    assert rmg_pdep.T == [300, 2500, 10]
    assert rmg_pdep.P == [0.01, 100, 10]
    assert rmg_pdep.interpolation == 'Chebyshev'
    assert rmg_pdep.T_basis_set == 6
    assert rmg_pdep.P_basis_set == 4
    assert rmg_pdep.max_atoms == 16

    with pytest.raises(ValidationError):
        # check that method has at least 2 characters
        RMGPDep(method='a')

    with pytest.raises(ValidationError):
        # check that method has at most 3 characters
        RMGPDep(method='abcd')

    with pytest.raises(ValidationError):
        # check that max_grain_size is constrained to > 0
        RMGPDep(max_grain_size=0)

    with pytest.raises(ValidationError):
        # check that max_number_of_grains is constrained to > 0
        RMGPDep(max_number_of_grains=0)

    with pytest.raises(ValidationError):
        # check that elements in list T are constrained to > 0
        RMGSpecies(T=[0, 0.01])

    with pytest.raises(ValidationError):
        # check that elements in list P are constrained to > 0
        RMGSpecies(P=[0, 0.01])

    with pytest.raises(ValidationError):
        # check that T_basis_set is constrained to > 0
        RMGPDep(T_basis_set=0)

    with pytest.raises(ValidationError):
        # check that P_basis_set is constrained to > 0
        RMGPDep(P_basis_set=0)

    with pytest.raises(ValidationError):
        # check that max_atoms is constrained to > 0
        RMGPDep(max_atoms=0)


def test_rmg_species_constraints_schema():
    """Test creating an instance of RMGSpeciesConstraints"""
    rmg_pdep = RMGSpeciesConstraints(allowed=['input species', 'seed mechanisms', 'reaction libraries'],
                                     max_C_atoms=30,
                                     max_O_atoms=10,
                                     max_N_atoms=10,
                                     max_Si_atoms=10,
                                     max_S_atoms=10,
                                     max_heavy_atoms=10,
                                     max_radical_electrons=2,
                                     max_singlet_carbenes=1,
                                     max_carbene_radicals=0,
                                     allow_singlet_O2=True
                                     )
    assert rmg_pdep.allowed == ['input species', 'seed mechanisms', 'reaction libraries']
    assert rmg_pdep.max_C_atoms == 30
    assert rmg_pdep.max_O_atoms == 10
    assert rmg_pdep.max_N_atoms == 10
    assert rmg_pdep.max_Si_atoms == 10
    assert rmg_pdep.max_S_atoms == 10
    assert rmg_pdep.max_heavy_atoms == 10
    assert rmg_pdep.max_radical_electrons == 2
    assert rmg_pdep.max_singlet_carbenes == 1
    assert rmg_pdep.max_carbene_radicals == 0
    assert rmg_pdep.allow_singlet_O2 is True

    with pytest.raises(ValidationError):
        # check that max_C_atoms is constrained to >= 0
        RMGPDep(max_C_atoms=-1)

    with pytest.raises(ValidationError):
        # check that max_O_atoms is constrained to >= 0
        RMGPDep(max_O_atoms=-1)

    with pytest.raises(ValidationError):
        # check that max_N_atoms is constrained to >= 0
        RMGPDep(max_N_atoms=-1)

    with pytest.raises(ValidationError):
        # check that max_Si_atoms is constrained to >= 0
        RMGPDep(max_Si_atoms=-1)

    with pytest.raises(ValidationError):
        # check that max_S_atoms is constrained to >= 0
        RMGPDep(max_S_atoms=-1)

    with pytest.raises(ValidationError):
        # check that max_heavy_atoms is constrained to >= 0
        RMGPDep(max_heavy_atoms=-1)

    with pytest.raises(ValidationError):
        # check that max_radical_electrons is constrained to >= 0
        RMGPDep(max_radical_electrons=-1)

    with pytest.raises(ValidationError):
        # check that max_singlet_carbenes is constrained to >= 0
        RMGPDep(max_singlet_carbenes=-1)

    with pytest.raises(ValidationError):
        # check that max_carbene_radicals is constrained to >= 0
        RMGPDep(max_carbene_radicals=-1)


def test_qm_schema():
    """Test The QM schema"""
    qm = {'adapter': 'ARC'}
    qm = QM(**qm)
    assert qm.adapter == 'ARC'
    assert qm.species == qm.reactions == list()
