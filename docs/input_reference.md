# Input file reference

This page documents **every** available T3 input argument.

This is not a runnable input file. It combines gas-phase and liquid-phase options
together for reference purposes.

```yaml
# A commented version of T3 input file
project: project_name  # required
project_directory: project_name  # optional

verbose: 20  # The logging level, optional. 10 - debug, 20 - info, 30 - warning, default: 20

# arguments related to T3 (optional)
t3:

  # options (optional block, all arguments have defaults)
  options:
    flux_adapter: RMG  # optional, can use any implemented simulation adapter, default: 'RMG'
    profiles_adapter: RMG  # optional, can use any implemented simulation adapter, default: 'RMG'
    collision_violators_thermo: false  # optional, whether to calculate thermo of species participating
                                      # in collision violating reactions, default: False
    collision_violators_rates: false  # optional, whether to calculate rates of core collision violating
                                     # reactions, default: False. If True, will force
                                     # collision_violators_thermo to be True if it's not
    all_core_species: false  # optional, whether to calculate thermo for all core species, default: False
    all_core_reactions: false  # optional, whether to calculate rates for all core species, default: False
    fit_missing_GAV: false  # optional, whether to capture wrong thermo groups of species estimated by
                            # RMG and attempt to calculate them, default: False
    max_T3_iterations: 10  # optional, maximum T3 iterations, default: 10
    max_RMG_exceptions_allowed: 10  # optional, maximum number of times RMG is allowed to crash,
                                    # default: 10
    max_RMG_walltime: '00:02:00:00'  # optional, format is 'DD:HH:MM:SS', default: '00:00:00:00'
    max_T3_walltime: '01:00:00:00'  # optional, format is 'DD:HH:MM:SS', default: None
    max_rmg_processes: null  # optional, maximum number of processes for RMG multiprocessing,
                             # default: None
    max_rmg_iterations: null  # optional, maximum number of internal RMG iterations, default: None
    library_name: T3lib  # optional, name of the RMG libraries T3 creates and saves as output,
                         # default: 'T3lib'
    shared_library_name: T3_shared  # optional, name of RMG libraries (kinetics and thermo) created
                                    # inside the respective RMG database paths that several T3
                                    # concurrent projects may share, default: None
    external_library_path: null  # optional, path to an external RMG library to use for saving shared
                                 # libraries, default: None (i.e., use the RMG database path)
    num_sa_per_temperature_range: 3  # optional, if a range of temperatures is given this argument
                                     # specifies how many equally distributed points to generate for
                                     # local SA runs, default: 3
    num_sa_per_pressure_range: 3  # optional, same as above for pressures, default: 3
    num_sa_per_volume_range: 3  # optional, same as above for volumes, default: 3
    num_sa_per_concentration_range: 3  # optional, same as above for species concentrations, default: 3
    modify_concentration_ranges_together: true  # optional, if species concentration is given as a
                                                # range, whether to vary them together (True) or create
                                                # combinations (False). Warning: False may result in
                                                # many slow SA runs. default: True
    modify_concentration_ranges_in_reverse: false  # optional, if concentration is given as a range
                                                   # for two species, whether to vary them inversely
                                                   # (e.g., fuel and oxidizer). default: False

  # sensitivity analysis (optional block, T3 can run w/o SA)
  sensitivity:
    adapter: CanteraConstantTP  # *required* (this is how SA is requested)
                                # Available adapters: 'CanteraConstantTP', 'CanteraConstantHP',
                                # 'CanteraConstantUV', 'CanteraPFR', 'CanteraPFRTProfile',
                                # 'CanteraJSR', 'RMGConstantTP'
    atol: 1e-6  # optional, default: 1e-6
    rtol: 1e-4  # optional, default: 1e-4
    global_observables: ['IDT', 'ESR', 'SL']  # optional, only implemented in Cantera adapters,
                                               # default: None
                                               # IDT = ignition delay time
                                               # ESR = extinction strain rate
                                               # SL  = laminar flame speed
    SA_threshold: 0.01  # optional, default: 0.01
    pdep_SA_threshold: 0.001  # optional, used to determine wells and reactions to calculate thermo
                              # and rates for from a PES of a sensitive reaction, default: 0.001
                              # Pass None to skip PES SA.
    ME_methods: ['CSE', 'MSC']  # master equation methods for PES SA,
                                # any combination of 'CSE', 'RS', 'MSC', default: ['CSE', 'MSC']
    top_SA_species: 10  # optional, used per observable to determine thermo to calculate, default: 10
    top_SA_reactions: 10  # optional, used per observable to determine rates to calculate as well as
                          # thermo of species participating in these reactions, default: 10
    T_list: [800, 1000, 1200]  # optional, explicit temperature list for SA runs (overrides
                               # range-based generation), Units: K, default: None
    P_list: [1, 10]  # optional, explicit pressure list for SA runs (overrides range-based
                     # generation), Units: bar, default: None

  # uncertainty analysis (optional block, T3 can run w/o UA)
  # either local or global UA type must be specified to execute an UA
  uncertainty:
    adapter: RMG  # *required* (this is how UA is requested), currently only implemented in the
                  # RMG adapter
    local_analysis: false  # optional, local UA using first-order sensitivity coefficients,
                           # default: False
    global_analysis: false  # optional, global UA varying largest sensitive input parameters,
                            # default: False
    correlated: true  # optional, whether to treat input parameters as correlated, default: True
    local_number: 10  # optional, number of reported parameters in local UA, default: 10
    global_number: 5  # optional, number of sensitive input parameters to vary in global UA per
                      # observable, applies independently to species and reactions, default: 5
    termination_time: null  # only necessary for global UA if a termination time wasn't specified in
                            # the RMG reactor
    PCE_run_time: 1800  # optional, time limit for adapting PCE to output, default: 1800, Units: s
    PCE_error_tolerance: null  # optional, target L2 error for PCE, default: None
    PCE_max_evals: null  # optional, max model evaluations for PCE adaptation, default: None
    logx: false  # optional, toggles output between mole fractions and log mole fractions,
                 # default: False

# arguments related to RMG, required
rmg:

  # general
  rmg_execution_type: incore  # optional, 'incore' or 'local', default: None (set from settings.py)
  memory: 15  # optional, RMG memory in GB (for server execution), default: None
  cpus: 16  # optional, number of processes for RMG, default: None

  # database (a required block)
  database:
    thermo_libraries: ['BurkeH2O2', 'DFT_QCI_thermo', 'primaryThermoLibrary', 'CBS_QB3_1dHR']
        # Can be None for auto-completion via chemistry_sets
    kinetics_libraries: ['BurkeH2O2inN2', 'NOx2018', 'Klippenstein_Glarborg2016']
        # Can be None for auto-completion via chemistry_sets
    chemistry_sets: ['primary', 'nitrogen', 'combustion']
        # Chemistry systems for which libraries will be auto-loaded. Can be None.
    use_low_credence_libraries: false  # Whether to use low credence libraries during
                                       # auto-completion, default: False
    transport_libraries: ['OneDMinN2', 'PrimaryTransportLibrary', 'NOx2018', 'GRI-Mech']
        # optional, default: ['OneDMinN2', 'PrimaryTransportLibrary', 'NOx2018', 'GRI-Mech']
    seed_mechanisms: []  # optional, default: []
    kinetics_depositories: default  # optional, default: 'default'
    kinetics_families: default  # optional, default: 'default'
    kinetics_estimator: rate rules  # optional, default: 'rate rules'

  # species (initial mixture) (a required block)
  # concentration units are mole fraction for gas phase and mol/cm3 for liquid phase
  # must specify either 'smiles', 'inchi', or 'adjlist'
  # not specifying 'concentration' is allowed and will result in a 0 initial concentration
  species:
    - label: ethane
      smiles: CC
      concentration: [1, 1.75]  # a concentration range can be defined (a length-2 list)
      reactive: true  # optional, default: True
      xyz: [ethane.log]  # each entry could be a string/dict XYZ format or a file path
      seed_all_rads: ['radical', 'alkoxyl', 'peroxyl']  # radical derivatives for the RMG input
    - label: OH
      smiles: '[OH]'
      observable: true  # optional, will be used as SA observable, default: False
      SA_observable: true  # optional, will be used as an SA observable, default: False
    - label: O2
      smiles: '[O][O]'
      concentration: 3.5
    - label: N2
      adjlist: |
        1 N u0 p1 c0 {2,T}
        2 N u0 p1 c0 {1,T}
      concentration: 13.16  # mole fraction (e.g., 3.5 * 3.76 for air)
      balance: true  # optional, only for gas phase simulations, default: False
      # note: 'constant' is only allowed in liquid phase simulations
      # note: 'solvent' is only allowed in liquid phase simulations

  # reactors (List[dict]) (a required block)
  # reactor type: 'gas batch constant T P' or 'liquid batch constant T V'
  # at least one termination criterion must be given per reactor
  # having a termination_time is recommended (also used for SA simulations)
  # users may specify multiple reactors, but all must be the same phase
  reactors:
    - type: gas batch constant T P
      T: [800, 1750]  # float (single T) or list (T range), Units: K
      P: 1e0  # float (single P) or list (P range), Units: bar
      termination_conversion:
        ethane: 0.2
      termination_time: [5, 's']  # allowed units: 'micro-s', 'ms', 's', 'hours', 'days'
      termination_rate_ratio: 0.01
      conditions_per_iteration: 12  # optional, default: 12 (nSims in RMG)

  # model (a required block)
  model:
    # primary tolerances:
    core_tolerance: [0.05, 0.01]  # *required*, float or list (toleranceMoveToCore)
    atol: 1e-16  # optional, default: 1e-16
    rtol: 1e-8  # optional, default: 1e-8
    # filtering:
    filter_reactions: true  # optional, default: True
    filter_threshold: 1e+8  # optional, filtering reactions threshold
    # pruning:
    tolerance_interrupt_simulation: [0.05, 0.01]  # optional, set equal to core_tolerance if omitted
    min_core_size_for_prune: 50  # optional, pruning
    min_species_exist_iterations_for_prune: 2  # optional, pruning
    tolerance_keep_in_edge: 0.02  # optional, pruning
    maximum_edge_species: 1000000  # optional, pruning
    tolerance_thermo_keep_species_in_edge:  # optional, thermo pruning
    # staging:
    max_num_species: null  # optional, staging
    # dynamics:
    tolerance_move_edge_reaction_to_core:  # optional, dynamics criterion
    tolerance_move_edge_reaction_to_core_interrupt: 5.0  # optional, dynamics criterion
    dynamics_time_scale: [0.0, 'sec']  # optional, dynamics criterion
    # multiple objects:
    max_num_objs_per_iter: 1  # optional, multiple objects
    terminate_at_max_objects: false  # optional, multiple objects
    # misc:
    ignore_overall_flux_criterion: false  # optional
    tolerance_branch_reaction_to_core: 0.001  # optional
    branching_index: 0.5  # optional
    branching_ratio_max: 1.0  # optional
    # surface algorithm:
    tolerance_move_edge_reaction_to_surface: null
    tolerance_move_surface_species_to_core: null
    tolerance_move_surface_reaction_to_core: null
    tolerance_move_edge_reaction_to_surface_interrupt: null

  # pressure dependence (optional block, gas phase only)
  pdep:
    method: MSC  # *required*, options: 'CSE', 'RS', 'MSC'
    max_grain_size: 2  # optional, kJ/mol, default: 2
    max_number_of_grains: 250  # optional, default: 250
    T: [300, 2500, 10]  # optional, [T_min, T_max, n_points], K, default: [300, 2500, 10]
    P: [0.01, 100, 10]  # optional, [P_min, P_max, n_points], bar, default: [0.01, 100, 10]
    interpolation: Chebyshev  # optional, 'PDepArrhenius' or 'Chebyshev', default: 'Chebyshev'
    T_basis_set: 6  # optional, Chebyshev only, default: 6
    P_basis_set: 4  # optional, Chebyshev only, default: 4
    max_atoms: 16  # optional, default: 16

  # options (optional block)
  options:
    seed_name: Seed  # optional, default: 'Seed'
    save_edge: false  # optional, default: False (saveEdgeSpecies)
    save_html: false  # optional, default: False
    generate_seed_each_iteration: true  # optional, default: True
    save_seed_to_database: false  # optional, default: False
    units: si  # optional, currently only 'si' is supported
    generate_plots: false  # optional, RMG statistics plots, default: False
    save_simulation_profiles: false  # optional, save RMG .csv profiles, default: False
    verbose_comments: false  # optional, verbose chemkin comments, default: False
    keep_irreversible: false  # optional, default: False
    trimolecular_product_reversible: true  # optional, default: True
    save_seed_modulus: -1  # optional, save seed every n iterations (-1 = last only), default: -1

  # species constraints (optional block)
  species_constraints:
    allowed: ['input species', 'seed mechanisms', 'reaction libraries']  # optional
    max_C_atoms: 10  # required
    max_O_atoms: 10  # required
    max_N_atoms: 10  # required
    max_Si_atoms: 10  # required
    max_S_atoms: 10  # required
    max_heavy_atoms: 10  # required
    max_radical_electrons: 2  # required
    max_singlet_carbenes: 1  # optional, default: 1
    max_carbene_radicals: 0  # optional, default: 0
    allow_singlet_O2: true  # optional, default: True in T3 (False in RMG)

# arguments related to QM calcs, required to run QM calcs, otherwise T3 will only spawn RMG
qm:
  # currently only ARC is supported
  # All legal ARC arguments are allowed here
  # If 'species' or 'reactions' are specified, ARC will be spawned prior to RMG (iteration 0)
  adapter: ARC
  # any legal ARC argument can come here, for example:
  adaptive_levels:
    (1, 6):
      opt_level: wb97xd/wb97xd/def2tzvp
      sp: ccsd(t)-f12/aug-cc-pvtz-f12
    (7, 30):
      conformer_level:
        method: wb97xd
        basis: def2svp
      opt_level:
        method: wb97xd
        basis: def2tzvp
      sp_level:
        method: dlpno-ccsd(T)
        basis: def2-tzvp
        auxiliary_basis: def2-tzvp/c
    (31, 'inf'):
      opt_level: wb97xd/wb97xd/def2tzvp
  species:
    - label: vinoxy
      smiles: C=C[O]
```
