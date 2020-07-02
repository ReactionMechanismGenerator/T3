"""
t3 schema module
used for input validation

Todo: add "live" validators for actually implemented adapters, need to access the respective factory register dicts
"""

from enum import Enum
from typing import Dict, List, Optional, Union

from pydantic import BaseModel, conint, confloat, constr, root_validator, validator


class TerminationTimeEnum(str, Enum):
    """
    The supported termination type units in an RMG reactor.
    """
    micro_s = 'micro-s'
    ms = 'ms'
    s = 's'
    hrs = 'hrs'
    hours = 'hours'
    days = 'days'


class T3Options(BaseModel):
    """
    A class for validating input.T3.options arguments
    """
    flux_adapter: constr(max_length=255) = 'RMG'
    profiles_adapter: constr(max_length=255) = 'RMG'
    collision_violators_thermo: bool = False
    collision_violators_rates: bool = False
    all_core_species: bool = False
    all_core_reactions: bool = False
    fit_missing_GAV: bool = False
    max_T3_iterations: conint(gt=0) = 10
    max_RMG_exceptions_allowed: Optional[conint(ge=0)] = 10
    max_RMG_walltime: constr(regex=r'\d+:\d\d:\d\d:\d\d') = '00:00:00:00'
    max_T3_walltime: Optional[constr(regex=r'\d+:\d\d:\d\d:\d\d')] = None
    library_name: constr(max_length=255) = 'T3'
    max_rmg_processes: Optional[conint(ge=1)] = None
    max_rmg_iterations: Optional[conint(ge=1)] = None

    class Config:
        extra = "forbid"

    @validator('collision_violators_rates')
    def check_collision_violators_rates(cls, value, values):
        """T3Options.collision_violators_rates validator"""
        if 'collision_violators_thermo' in values and value:
            values['collision_violators_thermo'] = True
        return value


class T3Sensitivity(BaseModel):
    """
    A class for validating input.T3.sensitivity arguments
    """
    adapter: Optional[constr(max_length=255)] = None
    atol: confloat(gt=0, lt=1e-1) = 1e-6
    rtol: confloat(gt=0, lt=1e-1) = 1e-4
    global_observables: Optional[List[constr(min_length=2, max_length=3)]] = None
    SA_threshold: confloat(gt=0, lt=0.5) = 0.01
    pdep_SA_threshold: Optional[confloat(gt=0, lt=0.5)] = 0.001
    ME_methods: List[constr(min_length=2, max_length=3)] = ['CSE', 'MSC']
    top_SA_species: conint(ge=0) = 10
    top_SA_reactions: conint(ge=0) = 10
    T_list: Optional[List[confloat(gt=0)]] = None
    P_list: Optional[List[confloat(gt=0)]] = None

    class Config:
        extra = "forbid"

    @validator('global_observables')
    def check_global_observables(cls, value):
        """T3Sensitivity.global_observables validator"""
        if value is not None:
            for i, entry in enumerate(value):
                if entry.lower() not in ['igd', 'esr', 'sl']:
                    raise ValueError(f'The global observables list must contain a combination of "IgD", "ESR", and "SL", '
                                     f'Got {entry} in {value}')
                if entry.lower() in [value[j].lower() for j in range(i)]:
                    raise ValueError(f'The global observables list must not contain repetitions, got {value}')
        return value

    @validator('ME_methods', always=True)
    def check_me_methods(cls, value):
        """T3Sensitivity.ME_methods validator"""
        if value is None or not value:
            raise ValueError(f'The ME methods argument cannot be None or empty.')
        for i, entry in enumerate(value):
            if entry.lower() not in ['cse', 'rs', 'msc']:
                raise ValueError(f'The ME methods list must contain a combination of "CSE", "RS", and "MSC", '
                                 f'Got {entry} in {value}')
            if entry.lower() in [value[j].lower() for j in range(i)]:
                raise ValueError(f'The ME methods list must not contain repetitions, got {value}')
        return value


class T3Uncertainty(BaseModel):
    """
    A class for validating input.T3.uncertainty arguments
    """
    adapter: Optional[constr(max_length=255)] = None
    local_analysis: bool = False
    global_analysis: bool = False
    correlated: bool = True
    local_number: conint(gt=0) = 10
    global_number: conint(gt=0) = 5
    termination_time: Optional[constr(regex=r'\d+:\d\d:\d\d:\d\d')] = None
    PCE_run_time: conint(gt=0) = 1800
    PCE_error_tolerance: Optional[confloat(gt=0)] = None
    PCE_max_evals: Optional[conint(gt=0)] = None
    logx: bool = False

    class Config:
        extra = "forbid"


class RMGDatabase(BaseModel):
    """
    A class for validating input.RMG.database arguments
    """
    thermo_libraries: List[str]
    kinetics_libraries: List[str]
    transport_libraries: List[str] = ['OneDMinN2', 'PrimaryTransportLibrary', 'NOx2018', 'GRI-Mech']
    seed_mechanisms: List[str] = list()
    kinetics_depositories: Union[List[str], str] = 'default'
    kinetics_families: Union[str, List[str]] = 'default'
    kinetics_estimator: str = 'rate rules'

    class Config:
        extra = "forbid"


class RadicalTypeEnum(str, Enum):
    """
    The supported radical ``types`` entries for ``generate_radicals()``.
    """
    radical = 'radical'
    alkoxyl = 'alkoxyl'
    peroxyl = 'peroxyl'


class RMGSpecies(BaseModel):
    """
    A class for validating input.RMG.species arguments
    """
    label: str
    concentration: confloat(ge=0) = 0
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    adjlist: Optional[str] = None
    reactive: bool = True
    observable: bool = False
    SA_observable: bool = False
    UA_observable: bool = False
    constant: bool = False
    balance: bool = False
    solvent: bool = False
    xyz: Optional[List[Union[dict, str]]] = None
    seed_all_rads: Optional[List[RadicalTypeEnum]] = None

    class Config:
        extra = "forbid"


class RMGReactor(BaseModel):
    """
    A class for validating input.RMG.reactors arguments
    """
    type: str
    T: Union[confloat(gt=0), List[confloat(gt=0)]]
    P: Optional[Union[confloat(gt=0), List[confloat(gt=0)]]] = None
    termination_conversion: Optional[Dict[str, confloat(gt=0, lt=1)]] = None
    termination_time: Optional[List[Union[confloat(gt=0), TerminationTimeEnum]]] = None
    termination_rate_ratio: Optional[confloat(gt=0, lt=1)] = None
    conditions_per_iteration: conint(gt=0) = 12

    class Config:
        extra = "forbid"

    @validator('type')
    def check_reactor_type(cls, value):
        """RMGReactor.type validator"""
        supported_reactors = ['gas batch constant T P', 'liquid batch constant T V']
        # all supporter reactors must contain a 'gas' or 'liquid' keyword, other schema validations depend on it
        if value not in supported_reactors:
            raise ValueError(f'Supported RMG reactors are\n{supported_reactors}\nGot: "{value}"')
        return value

    @validator('T')
    def check_t(cls, value):
        """RMGReactor.T validator"""
        if isinstance(value, list) and len(value) != 2:
            raise ValueError(f'When specifying the temperature as a list, only two values are allowed (T min, T max),\n'
                             f'got {len(value)} values: {value}.')
        return value

    @validator('P', always=True)
    def check_p(cls, value, values):
        """RMGReactor.P validator"""
        if isinstance(value, list) and len(value) != 2:
            raise ValueError(f'When specifying the pressure as a list, only two values are allowed (P min, P max),\n'
                             f'got {len(value)} values: {value}.')
        if 'type' in values and 'gas' in values['type'] and value is None:
            raise ValueError('The reactor pressure must be specified for a gas phase reactor.')
        if 'type' in values and 'liquid' in values['type'] and value is not None:
            raise ValueError('A reactor pressure cannot be specified for a liquid phase reactor.')
        return value

    @validator('termination_time')
    def check_termination_time(cls, value):
        """RMGReactor.termination_time validator"""
        if len(value) != 2 or not isinstance(value[0], float) or not isinstance(value[1], str):
            raise ValueError(f'The specified termination time must be a list of 2 entries: '
                             f'the value (a float) and the units (a string). Got: {value}')
        if value[1] == TerminationTimeEnum.micro_s:
            value[0] *= 1000
            value[1] = TerminationTimeEnum.ms
        elif value[1] == TerminationTimeEnum.hrs:
            value[1] = TerminationTimeEnum.hours
        value[1] = value[1].value  # convert the Enum class into a string
        return tuple(value)


class RMGModel(BaseModel):
    """
    A class for validating input.RMG.model arguments
    """
    # primary_tolerances:
    core_tolerance: Union[confloat(gt=0, lt=1), List[confloat(gt=0, lt=1)]]
    atol: confloat(gt=0, lt=1e-1) = 1e-16
    rtol: confloat(gt=0, lt=1e-1) = 1e-8
    # filtering:
    filter_reactions: bool = False
    filter_threshold: conint(gt=0) = 1e8
    # pruning:
    tolerance_interrupt_simulation: Optional[Union[confloat(gt=0), List[confloat(gt=0)]]] = None
    min_core_size_for_prune: Optional[conint(gt=0)] = None
    min_species_exist_iterations_for_prune: Optional[conint(gt=0)] = None
    tolerance_keep_in_edge: Optional[confloat(gt=0)] = None
    maximum_edge_species: Optional[conint(gt=0)] = None
    tolerance_thermo_keep_species_in_edge: Optional[confloat(gt=0)] = None
    # staging:
    max_num_species: Optional[conint(gt=0)] = None
    # dynamics:
    tolerance_move_edge_reaction_to_core: Optional[confloat(gt=0)] = None
    tolerance_move_edge_reaction_to_core_interrupt: Optional[confloat(gt=0)] = None
    dynamics_time_scale: Optional[tuple] = None
    # multiple_objects:
    max_num_objs_per_iter: conint(gt=0) = 1
    terminate_at_max_objects: bool = False
    # misc:
    ignore_overall_flux_criterion: Optional[bool] = None
    tolerance_branch_reaction_to_core: Optional[confloat(gt=0)] = None
    branching_index: Optional[confloat(gt=0)] = None
    branching_ratio_max: Optional[confloat(gt=0)] = None
    # surface algorithm
    tolerance_move_edge_reaction_to_surface: Optional[confloat(gt=0)] = None
    tolerance_move_surface_species_to_core: Optional[confloat(gt=0)] = None
    tolerance_move_surface_reaction_to_core: Optional[confloat(gt=0)] = None
    tolerance_move_edge_reaction_to_surface_interrupt: Optional[confloat(gt=0)] = None

    class Config:
        extra = "forbid"

    @validator('core_tolerance')
    def check_core_tolerance(cls, value):
        """
        RMGModel.core_tolerance validator
        set core_tolerance to always be a list
        """
        return [value] if isinstance(value, float) else value

    @validator('tolerance_interrupt_simulation', always=True)
    def check_tolerance_interrupt_simulation(cls, value, values):
        """
        RMGModel.tolerance_interrupt_simulation validator
        set tolerance_interrupt_simulation to match core_tolerance
        """
        if values['core_tolerance'] is not None:
            if value is None:
                return values['core_tolerance']
            else:
                if isinstance(value, float) and isinstance(values['core_tolerance'], list):
                    value = [value] * len(values['core_tolerance'])
                elif isinstance(value, list) and isinstance(values['core_tolerance'], list):
                    if len(value) < len(values['core_tolerance']):
                        value = value + values[-1] * (len(values['core_tolerance']) - len(value))
                    elif len(value) > len(values['core_tolerance']):
                        raise ValueError(f'The length of tolerance_interrupt_simulation ({len(value)}) '
                                         f'cannot be greater than the length of core_tolerance '
                                         f'({len(values["core_tolerance"])}).')
        return value


class RMGOptions(BaseModel):
    """
    A class for validating input.RMG.options arguments
    """
    seed_name: str = 'Seed'
    save_edge: bool = True
    save_html: bool = False
    generate_seed_each_iteration: bool = True
    save_seed_to_database: bool = False
    units: str = 'si'
    generate_plots: bool = False
    save_simulation_profiles: bool = False
    verbose_comments: bool = False
    keep_irreversible: bool = False
    trimolecular_product_reversible: bool = True
    save_seed_modulus: conint(ge=-1) = -1

    class Config:
        extra = "forbid"

    @validator('units', always=True)
    def check_units(cls, value):
        """RMGOptions.units validator"""
        if value.lower() is not None and value != 'si':
            raise ValueError(f'Currently RMG only supports SI units, got "{value}"')
        return value.lower()


class RMGPDep(BaseModel):
    """
    A class for validating input.RMG.pdep arguments
    """
    method: constr(min_length=2, max_length=3)
    max_grain_size: confloat(gt=0) = 2
    max_number_of_grains: conint(gt=0) = 250
    T: List[Union[conint(gt=0), confloat(gt=0)]] = [300, 2500, 10]
    P: List[Union[conint(gt=0), confloat(gt=0)]] = [0.01, 100, 10]
    interpolation: str = 'Chebyshev'
    T_basis_set: conint(gt=0) = 6
    P_basis_set: conint(gt=0) = 4
    max_atoms: conint(gt=0) = 16

    class Config:
        extra = "forbid"

    @validator('method')
    def check_method(cls, value):
        """RMGPDep.method validator"""
        if value not in ['CSE', 'RS', 'MSC']:
            raise ValueError(f"The PDep method must be either 'CSE', 'RS', or 'MSC'.\nGot: {value}")
        return value

    @validator('T')
    def check_t(cls, value):
        """RMGPDep.T validator"""
        if len(value) != 3:
            raise ValueError(f'The temperature range must be a length three list (T min, T max, T count),\n'
                             f'got a length {len(value)} list: {value}.')
        if value[1] <= value[0]:
            raise ValueError(f'The following T range (T min, T max, T count) does not make sense:\n{value}')
        if not isinstance(value[2], int):
            raise ValueError(f'T count {value[2]} must be an integer, got a {type(value[2])}')
        return value

    @validator('P')
    def check_p(cls, value):
        """RMGPDep.P validator"""
        if len(value) != 3:
            raise ValueError(f'The pressure range must be a length three list (P min, P max, P count),\n'
                             f'got a length {len(value)} list: {value}.')
        if value[1] <= value[0]:
            raise ValueError(f'The following P range (P min, P max, P count) does not make sense:\n{value}')
        if not isinstance(value[2], int):
            raise ValueError(f'P count {value[2]} must be an integer, got a {type(value[2])}')
        return value

    @validator('interpolation', always=True)
    def check_interpolation(cls, value):
        """RMGPDep.interpolation validator"""
        if value not in ['PDepArrhenius', 'Chebyshev']:
            raise ValueError(f'The RMG PDep interpolation method must be either "PDepArrhenius" or "Chebyshev" '
                             f'(recommended).\nGot {value}')
        return value

    @validator('T_basis_set', always=True)
    def check_t_basis_set(cls, value, values):
        """RMGPDep.T_basis_set validator"""
        if values['T'] is not None and value >= values['T'][2] \
                and values['interpolation'] is not None and values['interpolation'] == 'Chebyshev':
            raise ValueError(f'The T_basis_set must be lower than the number of T points.\n'
                             f'Got {value} and {values["T"]}')
        return value

    @validator('P_basis_set', always=True)
    def check_p_basis_set(cls, value, values):
        """RMGPDep.P_basis_set validator"""
        if values['P'] is not None and value >= values['P'][2] \
                and values['interpolation'] is not None and values['interpolation'] == 'Chebyshev':
            raise ValueError(f'The P_basis_set must be lower than the number of P points.\n'
                             f'Got {value} and {values["P"]}')
        return value


class RMGSpeciesConstraints(BaseModel):
    """
    A class for validating input.RMG.species_constraints arguments
    """
    allowed: List[str] = ['input species', 'seed mechanisms', 'reaction libraries']
    max_C_atoms: conint(ge=0)
    max_O_atoms: conint(ge=0)
    max_N_atoms: conint(ge=0)
    max_Si_atoms: conint(ge=0)
    max_S_atoms: conint(ge=0)
    max_heavy_atoms: conint(ge=0)
    max_radical_electrons: conint(ge=0)
    max_singlet_carbenes: conint(ge=0) = 1
    max_carbene_radicals: conint(ge=0) = 0
    allow_singlet_O2: bool = True

    class Config:
        extra = "forbid"

    @validator('allowed')
    def check_allowed(cls, value):
        """RMGSpeciesConstraints.allowed validator"""
        for val in value:
            if val not in ['input species', 'seed mechanisms', 'reaction libraries']:
                raise ValueError(f"The allowed species in the RMG species constraints list must be in\n"
                                 f"['input species', 'seed mechanisms', 'reaction libraries'].\n"
                                 f"Got: {val} in {value}")
        return value


class T3(BaseModel):
    """
    A class for validating input.T3 arguments
    """
    options: Optional[T3Options] = None
    sensitivity: Optional[T3Sensitivity] = None
    uncertainty: Optional[T3Uncertainty] = None

    class Config:
        extra = "forbid"

    @validator('options', always=True)
    def check_options(cls, value):
        """T3.options validator"""
        return value or T3Options()


class RMG(BaseModel):
    """
    A class for validating input.RMG arguments
    """
    database: RMGDatabase
    reactors: List[RMGReactor]
    species: List[RMGSpecies]
    model: RMGModel
    pdep: Optional[RMGPDep] = None
    options: Optional[RMGOptions] = None
    species_constraints: Optional[RMGSpeciesConstraints] = None

    class Config:
        extra = "forbid"

    @validator('database', always=True)
    def check_database(cls, value):
        """RMG.database validator"""
        if value is None or not value:
            raise ValueError('RMG database must be specified')
        return value

    @validator('reactors', always=True)
    def check_reactors(cls, value):
        """RMG.reactors validator"""
        if value is None or not value:
            raise ValueError('RMG reactors must be specified')
        return value

    @validator('species', always=True)
    def check_species(cls, value):
        """RMG.species validator"""
        if value is None or not value:
            raise ValueError('RMG species must be specified')
        return value

    @validator('model', always=True)
    def check_model(cls, value):
        """RMG.model validator"""
        if value is None or not value:
            raise ValueError('RMG model must be specified')
        return value

    @validator('pdep')
    def check_pdep_only_if_gas_phase(cls, value, values):
        """RMG.pdep validator"""
        if value is not None and values['reactors'] is not None:
            reactor_types = set([reactor.type for reactor in values['reactors']])
            if value is not None and not any(['gas' in reactor for reactor in reactor_types]):
                raise ValueError(f'A pdep section can only be specified for gas phase reactors, got: {reactor_types}')
        return value

    @root_validator(pre=True)
    def check_species_and_reactors(cls, values):
        if 'reactors' in values and values['reactors'] is not None \
                and 'species' in values and values['species'] is not None:
            # check species attributes (balance, solvation) make sense
            reactor_types = set([reactor['type'] for reactor in values['reactors']])
            balance_species, solvent_species = list(), list()
            for species in values['species']:
                if not balance_species and 'balance' in species and species['balance']:
                    balance_species.append(species['label'])
                if not solvent_species and 'solvent' in species and species['solvent']:
                    solvent_species.append(species['label'])
            gas_reactors, liquid_reactors = False, False
            for reactor_type in reactor_types:
                if 'gas' in reactor_type:
                    gas_reactors = True
                if 'liquid' in reactor_type:
                    liquid_reactors = True
            if liquid_reactors and len(balance_species):
                raise ValueError(f'A species ({balance_species}) cannot be set as a balance species '
                                 f'if liquid phase reactors are defined.')
            if gas_reactors and len(solvent_species):
                raise ValueError(f'A species ({solvent_species}) cannot be set as a solvent '
                                 f'if gas phase reactors are defined.')
            if len(balance_species) > 1:
                raise ValueError(f'Only a single species may be defined as balance,\ngot: {balance_species}')
            if len(solvent_species) > 1:
                raise ValueError(f'Only a single species may be defined as a solvent,\ngot: {solvent_species}')
            # check reactor termination_conversion has a corresponding species
            for reactor in values['reactors']:
                if 'termination_conversion' in reactor and reactor['termination_conversion'] is not None:
                    species_labels = [species['label'] for species in values['species']]
                    for termination_conversion_label in reactor['termination_conversion'].keys():
                        if termination_conversion_label not in species_labels:
                            raise ValueError(f'No species with label "{termination_conversion_label}" was defined. '
                                             f'Make sure all species labels defined under termination_conversion '
                                             f'have actual corresponding species.')
        return values

    @validator('options', always=True)
    def check_options(cls, value):
        """RMG.options validator"""
        return value or RMGOptions()


class QM(BaseModel):
    """
    A class for validating input.QM arguments
    """
    adapter: str = 'ARC'
    species: Optional[list] = None
    reactions: Optional[list] = None

    class Config:
        extra = "allow"

    @validator('adapter')
    def check_adapter(cls, value):
        """QM.adapter validator"""
        supported_qm_adapters = ['ARC']
        if value not in supported_qm_adapters:
            raise ValueError(f'Supported QM adapters are:\n{supported_qm_adapters}\nGot:{value}')
        return value

    @validator('species', 'reactions', always=True)
    def check_species(cls, value):
        """QM.species and QM.reactions validator"""
        value = value or list()
        return value


class InputBase(BaseModel):
    """
    An InputBase class for validating input arguments
    """
    project: constr(max_length=255)
    project_directory: Optional[constr(max_length=255)] = None
    verbose: conint(ge=10, le=30, multiple_of=10) = 20
    t3: Optional[T3] = None
    rmg: RMG
    qm: Optional[QM] = None

    class Config:
        extra = "forbid"

    @validator('t3', always=True)
    def check_t3(cls, value):
        """InputBase.t3 validator"""
        return value or T3()

    @validator('qm', always=True)
    def check_qm(cls, value):
        """InputBase.qm validator"""
        return value or dict()

    @root_validator(pre=True)
    def validate_rmg_t3(cls, values):
        """InputBase.validate_rmg_t3"""
        if 'rmg' in values and 't3' in values and values['t3']:
            # check termination time for global UA
            if 'uncertainty' in values['t3'] and values['t3']['uncertainty'] is not None:
                ua_termination_time = values['t3']['uncertainty']['termination_time']
                rmg_reactor_termination_times = [reactor['termination_time'] for reactor in values['rmg']['reactors']]
                if all([termination_time is None for termination_time in rmg_reactor_termination_times]) \
                        and ua_termination_time is None and values['t3']['uncertainty']['global_analysis']:
                    raise ValueError('If a global uncertainty analysis is requested, a termination time must be '
                                     'specified either under "t3.uncertainty.termination_time" '
                                     'or in at least one RMG reactor.')
            # check solvent for liquid phase
            reactor_types = set([reactor['type'] for reactor in values['rmg']['reactors']])
            solvents = list()
            for species in values['rmg']['species']:
                if 'solvent' in species and species['solvent']:
                    solvents.append(species['label'])
            if any(['liquid' in reactor_type for reactor_type in reactor_types]):
                if not len(solvents):
                    raise ValueError('One species must be defined as the solvent when using liquid phase reactors.')
                if len(solvents) > 1:
                    raise ValueError(f'Only one solvent can be specified, got: {solvents}')
            else:
                if len(solvents):
                    raise ValueError(f'No solvent species are allowed for gas phase reactors, got: {solvents}')
            # check core_tolerance and max_T3_iterations
            if 'model' in values['rmg'] and 'core_tolerance' in values['rmg']['model'] \
                    and 'options' in values['t3'] and 'max_T3_iterations' in values['t3']['options'] \
                    and not isinstance(values['rmg']['model']['core_tolerance'], float) \
                    and len(values['rmg']['model']['core_tolerance']) > values['t3']['options']['max_T3_iterations']:
                raise ValueError(f'The number of RMG core tolerances ({len(values["rmg"]["model"]["core_tolerance"])}) '
                                 f'cannot be greater than the max number of T3 iterations '
                                 f'({values["t3"]["options"]["max_T3_iterations"]}).')
        return values
