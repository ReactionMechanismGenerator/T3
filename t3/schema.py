"""
t3 schema module
used for input validation
"""

import os
from enum import Enum
from typing import Annotated, Dict, List, Optional, Tuple, Union

from pydantic import BaseModel, Field, ValidationInfo, field_validator, model_validator

from arc.common import read_yaml_file

from t3.common import DATA_BASE_PATH, VALID_CHARS
from t3.simulate.factory import _registered_simulate_adapters


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
    flux_adapter: Annotated[str, Field(max_length=255)] = 'RMG'
    profiles_adapter: Annotated[str, Field(max_length=255)] = 'RMG'
    collision_violators_thermo: bool = False
    collision_violators_rates: bool = False
    all_core_species: bool = False
    all_core_reactions: bool = False
    fit_missing_GAV: bool = False
    max_T3_iterations: Annotated[int, Field(gt=0)] = 10
    max_RMG_exceptions_allowed: Optional[Annotated[int, Field(ge=0)]] = 10
    max_RMG_walltime: Annotated[str, Field(pattern=r'\d+:\d\d:\d\d:\d\d')] = '00:00:00:00'
    max_T3_walltime: Optional[Annotated[str, Field(pattern=r'\d+:\d\d:\d\d:\d\d')]] = None
    max_rmg_processes: Optional[Annotated[int, Field(ge=1)]] = None
    max_rmg_iterations: Optional[Annotated[int, Field(ge=1)]] = None
    library_name: Annotated[str, Field(max_length=255)] = 'T3lib'
    shared_library_name: Optional[Annotated[str, Field(max_length=255)]] = None
    external_library_path: Optional[Annotated[str, Field(max_length=255)]] = None
    num_sa_per_temperature_range: Annotated[int, Field(ge=1)] = 3
    num_sa_per_pressure_range: Annotated[int, Field(ge=1)] = 3
    num_sa_per_volume_range: Annotated[int, Field(ge=1)] = 3
    num_sa_per_concentration_range: Annotated[int, Field(ge=1)] = 3
    modify_concentration_ranges_together: bool = True
    modify_concentration_ranges_in_reverse: bool = False

    class Config:
        extra = "forbid"

    @model_validator(mode='after')
    def enforce_collision_thermo(self) -> 'T3Options':
        """
        If collision_violators_rates is True, ensure collision_violators_thermo is also True.
        """
        if self.collision_violators_rates:
            self.collision_violators_thermo = True
        return self

    @field_validator('library_name')
    @classmethod
    def check_library_name(cls, value):
        """T3Options.library_name validator"""
        for char in value:
            if char not in VALID_CHARS:
                raise ValueError(f'The library name "{value}" contains an invalid character: {char}.\n'
                                 f'Only the following characters are allowed:\n{VALID_CHARS}')
        return value

    @field_validator('shared_library_name')
    @classmethod
    def check_shared_library_name(cls, value):
        """T3Options.shared_library_name validator"""
        if value is not None:
            for char in value:
                if char not in VALID_CHARS + '/':
                    raise ValueError(f'The shared library name "{value}" contains an invalid character: {char}.\n'
                                     f'Only the following characters are allowed:\n{VALID_CHARS}')
        return value

    @field_validator('external_library_path')
    @classmethod
    def check_external_library_path(cls, value):
        """T3Options.external_library_path validator"""
        if value is not None:
            for char in value:
                if char not in VALID_CHARS + '/':
                    raise ValueError(f'The external library path "{value}" contains an invalid character: {char}.\n'
                                     f'Only the following characters are allowed:\n{VALID_CHARS}')
        return value


class T3Sensitivity(BaseModel):
    """
    A class for validating input.T3.sensitivity arguments
    """
    adapter: Annotated[str, Field(max_length=255)] = 'RMGConstantTP'
    atol: Annotated[float, Field(gt=0, lt=1e-1)] = 1e-6
    rtol: Annotated[float, Field(gt=0, lt=1e-1)] = 1e-4
    global_observables: Optional[List[Annotated[str, Field(min_length=2, max_length=3)]]] = None
    SA_threshold: Annotated[float, Field(gt=0, lt=0.5)] = 0.01
    pdep_SA_threshold: Optional[Annotated[float, Field(gt=0, lt=0.5)]] = 0.001
    ME_methods: List[Annotated[str, Field(min_length=2, max_length=3)]] = ['CSE', 'MSC']
    top_SA_species: Annotated[int, Field(ge=0)] = 10
    top_SA_reactions: Annotated[int, Field(ge=0)] = 10
    T_list: Optional[List[Annotated[float, Field(gt=0)]]] = None
    P_list: Optional[List[Annotated[float, Field(gt=0)]]] = None

    class Config:
        extra = "forbid"

    @field_validator('adapter')
    @classmethod
    def check_adapter(cls, value):
        """T3Sensitivity.adapter validator"""
        if value not in _registered_simulate_adapters.keys():
            raise ValueError(
                f'The "T3 sensitivity adapter" argument of {value} was not present in the keys for the '
                f'_registered_simulate_adapters dictionary: {list(_registered_simulate_adapters.keys())}'
                f'\nPlease check that the simulate adapter was registered properly.')
        return value

    @field_validator('global_observables')
    @classmethod
    def check_global_observables(cls, value):
        """T3Sensitivity.global_observables validator"""
        if value is not None:
            for i, entry in enumerate(value):
                if entry.lower() not in ['idt', 'esr', 'sl']:
                    raise ValueError(f'The global observables list must contain a combination of "IDT", "ESR", and "SL", '
                                     f'Got {entry} in {value}')
                if entry.lower() in [value[j].lower() for j in range(i)]:
                    raise ValueError(f'The global observables list must not contain repetitions, got {value}')
        return value

    @field_validator('ME_methods')
    @classmethod
    def check_me_methods(cls, value):
        """T3Sensitivity.ME_methods validator"""
        if value is None or not value:
            raise ValueError('The ME methods argument cannot be None or empty.')
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
    adapter: Optional[Annotated[str, Field(max_length=255)]] = None
    local_analysis: bool = False
    global_analysis: bool = False
    correlated: bool = True
    local_number: Annotated[int, Field(gt=0)] = 10
    global_number: Annotated[int, Field(gt=0)] = 5
    termination_time: Optional[Annotated[str, Field(pattern=r'\d+:\d\d:\d\d:\d\d')]] = None
    PCE_run_time: Annotated[int, Field(gt=0)] = 1800
    PCE_error_tolerance: Optional[Annotated[float, Field(gt=0)]] = None
    PCE_max_evals: Optional[Annotated[int, Field(gt=0)]] = None
    logx: bool = False

    class Config:
        extra = "forbid"


class RMGDatabase(BaseModel):
    """
    A class for validating input.RMG.database arguments
    """
    thermo_libraries: Optional[List[str]] = None
    kinetics_libraries: Optional[List[str]] = None
    chemistry_sets: Optional[List[str]] = None
    use_low_credence_libraries: bool = False
    transport_libraries: List[str] = ['OneDMinN2', 'PrimaryTransportLibrary', 'NOx2018', 'GRI-Mech']
    seed_mechanisms: List[str] = list()
    kinetics_depositories: Union[List[str], str] = 'default'
    kinetics_families: Union[str, List[str]] = 'default'
    kinetics_estimator: str = 'rate rules'

    class Config:
        extra = "forbid"

    @field_validator('chemistry_sets')
    @classmethod
    def check_chemistry_sets(cls, value, info: ValidationInfo):
        """RMGDatabase.chemistry_sets validator"""
        libraries_dict = read_yaml_file(path=os.path.join(DATA_BASE_PATH, 'libraries.yml'))
        allowed_values = libraries_dict.keys()
        if value and any(v not in allowed_values for v in value):
            raise ValueError(f'The chemistry sets must be within of the following:\n{allowed_values}\nGot: {value}')
        if value is None and (info.data.get('thermo_libraries') is None or info.data.get('kinetics_libraries') is None):
            raise ValueError('The chemistry set must be specified if thermo or kinetics libraries are not specified.')
        return value


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
    concentration: Union[Annotated[float, Field(ge=0)], Tuple[Annotated[float, Field(ge=0)], Annotated[float, Field(ge=0)]]] = 0
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
    xyz: Optional[Union[List[Union[dict, str]], dict, str]] = None
    seed_all_rads: Optional[List[RadicalTypeEnum]] = None

    class Config:
        extra = "forbid"

    @field_validator('constant')
    @classmethod
    def check_ranged_concentration_not_constant(cls, value, info: ValidationInfo):
        """RMGSpecies.constant validator"""
        label = ' for ' + info.data.get('label', '') if 'label' in info.data else ''
        if value and isinstance(info.data.get('concentration'), tuple):
            raise ValueError(f"A constant species cannot have a concentration range.\n"
                             f"Got{label}: {info.data.get('concentration')}.")
        return value

    @field_validator('concentration')
    @classmethod
    def check_concentration_range_order(cls, value, info: ValidationInfo):
        """Make sure the concentration range is ordered from the smallest to the largest"""
        label = ' for ' + info.data.get('label', '') if 'label' in info.data else ''
        if isinstance(value, tuple):
            if value[0] == value[1]:
                raise ValueError(f"A concentration range cannot contain to identical concentrations.\n"
                                 f"Got{label}: {value}.")
            if value[0] > value[1]:
                value = (value[1], value[0])
        return value

    @field_validator('balance')
    @classmethod
    def check_concentration_of_balance_species(cls, value, info: ValidationInfo):
        """Make sure the concentration of the balance species is defined as a scalar, not a range"""
        if value and 'concentration' in info.data:
            if not isinstance(info.data.get('concentration'), (int, float)):
                raise ValueError(f"The balance species concentration cannot be defined as a range, "
                                 f"got: {info.data.get('concentration')}.")
        return value

    @model_validator(mode='after')
    def set_balance_concentration_default(self):
        """Set concentration=1 for balance species if concentration is 0 (default)"""
        if self.balance and self.concentration == 0:
            self.concentration = 1
        return self


class RMGReactor(BaseModel):
    """
    A class for validating input.RMG.reactors arguments
    """
    type: str
    T: Union[Annotated[float, Field(gt=0)], List[Annotated[float, Field(gt=0)]]]
    P: Optional[Union[Annotated[float, Field(gt=0)], List[Annotated[float, Field(gt=0)]]]] = None
    V: Optional[Union[Annotated[float, Field(gt=0)], List[Annotated[float, Field(gt=0)]]]] = None
    termination_conversion: Optional[Dict[str, Annotated[float, Field(gt=0, lt=1)]]] = None
    termination_time: Optional[Tuple[Annotated[float, Field(gt=0)], TerminationTimeEnum]] = None
    termination_rate_ratio: Optional[Annotated[float, Field(gt=0, lt=1)]] = None
    conditions_per_iteration: Annotated[int, Field(gt=0)] = 12

    class Config:
        extra = "forbid"

    @field_validator('type')
    @classmethod
    def check_reactor_type(cls, value):
        """RMGReactor.type validator"""
        supported_reactors = ['gas batch constant T P', 'liquid batch constant T V']
        # all supporter reactors must contain a 'gas' or 'liquid' keyword, other schema validations depend on it
        if value not in supported_reactors:
            raise ValueError(f'Supported RMG reactors are\n{supported_reactors}\nGot: "{value}"')
        return value

    @field_validator('T')
    @classmethod
    def check_t(cls, value):
        """RMGReactor.T validator"""
        if isinstance(value, list) and len(value) != 2:
            raise ValueError(f'When specifying the temperature as a list, only two values are allowed (T min, T max),\n'
                             f'got {len(value)} values: {value}.')
        return value

    @field_validator('P')
    @classmethod
    def check_p(cls, value, info: ValidationInfo):
        """RMGReactor.P validator"""
        if isinstance(value, list) and len(value) != 2:
            raise ValueError(f'When specifying the pressure as a list, only two values are allowed (P min, P max),\n'
                             f'got {len(value)} values: {value}.')
        reactor_type = info.data.get('type')
        if reactor_type and 'gas' in reactor_type and value is None:
            raise ValueError('The reactor pressure must be specified for a gas-phase reactor.')
        if reactor_type and 'liquid' in reactor_type and value is not None:
            raise ValueError('A reactor pressure cannot be specified for a liquid-phase reactor.')
        return value

    @field_validator('V')
    @classmethod
    def check_v(cls, value, info: ValidationInfo):
        """RMGReactor.V validator"""
        if isinstance(value, list) and len(value) != 2:
            raise ValueError(f'When specifying the volume as a list, only two values are allowed (V min, V max),\n'
                             f'got {len(value)} values: {value}.')
        reactor_type = info.data.get('type')
        if reactor_type and 'liquid' in reactor_type and value is None:
            raise ValueError('The reactor volume must be specified for a liquid-phase reactor.')
        if reactor_type and 'gas' in reactor_type and value is not None:
            raise ValueError('A reactor volume cannot be specified for a gas-phase reactor.')
        return value

    @field_validator('termination_time')
    @classmethod
    def check_termination_time(cls, value):
        """RMGReactor.termination_time validator"""
        if len(value) != 2 or not isinstance(value[0], float) or not isinstance(value[1], str):
            raise ValueError(f'The specified termination time must be a list of 2 entries: '
                             f'the value (a float) and the units (a string). Got: {value}')
        if value[1] == TerminationTimeEnum.micro_s:
            value= (value[0] * 1000, TerminationTimeEnum.ms)
        elif value[1] == TerminationTimeEnum.hrs:
            value = (value[0] ,TerminationTimeEnum.hours)
        value = (value[0], value[1].value)  # convert the Enum class into a string
        return value


class RMGModel(BaseModel):
    """
    A class for validating input.RMG.model arguments
    """
    # primary_tolerances:
    core_tolerance: Union[Annotated[float, Field(gt=0, lt=1)], List[Annotated[float, Field(gt=0, lt=1)]]]
    atol: Annotated[float, Field(gt=0, lt=1e-1)] = 1e-16
    rtol: Annotated[float, Field(gt=0, lt=1e-1)] = 1e-8
    # filtering:
    filter_reactions: bool = True
    filter_threshold: Union[Annotated[float, Field(gt=0)], Annotated[int, Field(gt=0)]] = 1e+8
    # pruning:
    tolerance_interrupt_simulation: Optional[Union[Annotated[float, Field(gt=0)], List[Annotated[float, Field(gt=0)]]]] = None
    min_core_size_for_prune: Optional[Annotated[int, Field(gt=0)]] = None
    min_species_exist_iterations_for_prune: Optional[Annotated[int, Field(gt=0)]] = None
    tolerance_keep_in_edge: Optional[Annotated[float, Field(gt=0)]] = None
    maximum_edge_species: Optional[Annotated[int, Field(gt=0)]] = None
    tolerance_thermo_keep_species_in_edge: Optional[Annotated[float, Field(gt=0)]] = None
    # staging:
    max_num_species: Optional[Annotated[int, Field(gt=0)]] = None
    # dynamics:
    tolerance_move_edge_reaction_to_core: Optional[Annotated[float, Field(gt=0)]] = None
    tolerance_move_edge_reaction_to_core_interrupt: Optional[Annotated[float, Field(gt=0)]] = None
    dynamics_time_scale: Optional[tuple] = None
    # multiple_objects:
    max_num_objs_per_iter: Annotated[int, Field(gt=0)] = 1
    terminate_at_max_objects: bool = False
    # misc:
    ignore_overall_flux_criterion: Optional[bool] = None
    tolerance_branch_reaction_to_core: Optional[Annotated[float, Field(gt=0)]] = None
    branching_index: Optional[Annotated[float, Field(gt=0)]] = None
    branching_ratio_max: Optional[Annotated[float, Field(gt=0)]] = None
    # surface algorithm
    tolerance_move_edge_reaction_to_surface: Optional[Annotated[float, Field(gt=0)]] = None
    tolerance_move_surface_species_to_core: Optional[Annotated[float, Field(gt=0)]] = None
    tolerance_move_surface_reaction_to_core: Optional[Annotated[float, Field(gt=0)]] = None
    tolerance_move_edge_reaction_to_surface_interrupt: Optional[Annotated[float, Field(gt=0)]] = None

    class Config:
        extra = "forbid"

    @field_validator('core_tolerance')
    @classmethod
    def check_core_tolerance(cls, value):
        """
        RMGModel.core_tolerance validator
        set core_tolerance to always be a list
        """
        return [value] if isinstance(value, float) else value

    @field_validator('filter_threshold')
    @classmethod
    def check_filter_threshold(cls, value):
        """
        RMGModel.filter_threshold validator
        set filter_threshold to always be an integer
        Usually it is given in scientific writing, e.g., 1e+8, which cannot be automatically parsed as an int.
        """
        return int(value) if isinstance(value, float) else value

    @model_validator(mode='after')
    def check_tolerance_interrupt_simulation(self) -> 'RMGModel':
        """
        RMGModel.tolerance_interrupt_simulation validator
        Sets tolerance_interrupt_simulation to match core_tolerance if not provided,
        and ensures length consistency.
        """
        # Access the already-validated values from self
        core_tol = self.core_tolerance
        # We need to access the raw attribute, which might be None if optional
        tol_interrupt = self.tolerance_interrupt_simulation

        if core_tol is not None:
            # 1. Default to core_tolerance if not set
            if tol_interrupt is None:
                self.tolerance_interrupt_simulation = core_tol
                return self

            # 2. Logic for broadcasting float to list
            if isinstance(tol_interrupt, float) and isinstance(core_tol, list):
                self.tolerance_interrupt_simulation = [tol_interrupt] * len(core_tol)

            # 3. Logic for list length validation/extension
            elif isinstance(tol_interrupt, list) and isinstance(core_tol, list):
                if len(tol_interrupt) < len(core_tol):
                    # Extend with the last value
                    self.tolerance_interrupt_simulation = tol_interrupt + [tol_interrupt[-1]] * (
                                len(core_tol) - len(tol_interrupt))
                elif len(tol_interrupt) > len(core_tol):
                    raise ValueError(f'The length of tolerance_interrupt_simulation ({len(tol_interrupt)}) '
                                     f'cannot be greater than the length of core_tolerance '
                                     f'({len(core_tol)}).')

        return self


class RMGOptions(BaseModel):
    """
    A class for validating input.RMG.options arguments
    """
    seed_name: str = 'Seed'
    save_edge: bool = False
    save_html: bool = False
    generate_seed_each_iteration: bool = True
    save_seed_to_database: bool = False
    units: str = 'si'
    generate_plots: bool = False
    save_simulation_profiles: bool = False
    verbose_comments: bool = False
    keep_irreversible: bool = False
    trimolecular_product_reversible: bool = True
    save_seed_modulus: Annotated[int, Field(ge=-1)] = -1

    class Config:
        extra = "forbid"

    @field_validator('units')
    @classmethod
    def check_units(cls, value):
        """RMGOptions.units validator"""
        if value.lower() != 'si':
            raise ValueError(f'Currently RMG only supports SI units, got "{value}"')
        return value.lower()


class RMGPDep(BaseModel):
    """
    A class for validating input.RMG.pdep arguments
    """
    method: Annotated[str, Field(min_length=2, max_length=3)]
    max_grain_size: Annotated[float, Field(gt=0)] = 2
    max_number_of_grains: Annotated[int, Field(gt=0)] = 250
    T: List[Union[Annotated[int, Field(gt=0)], Annotated[float, Field(gt=0)]]] = [300, 2500, 10]
    P: List[Union[Annotated[int, Field(gt=0)], Annotated[float, Field(gt=0)]]] = [0.01, 100, 10]
    interpolation: str = 'Chebyshev'
    T_basis_set: Annotated[int, Field(gt=0)] = 6
    P_basis_set: Annotated[int, Field(gt=0)] = 4
    max_atoms: Annotated[int, Field(gt=0)] = 16

    class Config:
        extra = "forbid"

    @field_validator('method')
    @classmethod
    def check_method(cls, value):
        """RMGPDep.method validator"""
        if value not in ['CSE', 'RS', 'MSC']:
            raise ValueError(f"The PDep method must be either 'CSE', 'RS', or 'MSC'.\nGot: {value}")
        return value

    @field_validator('T')
    @classmethod
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

    @field_validator('P')
    @classmethod
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

    @field_validator('interpolation')
    @classmethod
    def check_interpolation(cls, value):
        """RMGPDep.interpolation validator"""
        if value not in ['PDepArrhenius', 'Chebyshev']:
            raise ValueError(f'The RMG PDep interpolation method must be either "PDepArrhenius" or "Chebyshev" '
                             f'(recommended).\nGot {value}')
        return value

    @field_validator('T_basis_set')
    @classmethod
    def check_t_basis_set(cls, value, info: ValidationInfo):
        """RMGPDep.T_basis_set validator"""
        if info.data.get('T') is not None and value >= info.data['T'][2] \
                and info.data.get('interpolation') is not None and info.data['interpolation'] == 'Chebyshev':
            raise ValueError(f'The T_basis_set must be lower than the number of T points.\n'
                             f'Got {value} and {info.data.get("T")}')
        return value

    @field_validator('P_basis_set')
    @classmethod
    def check_p_basis_set(cls, value, info: ValidationInfo):
        """RMGPDep.P_basis_set validator"""
        if info.data.get('P') is not None and value >= info.data['P'][2] \
                and info.data.get('interpolation') is not None and info.data['interpolation'] == 'Chebyshev':
            raise ValueError(f'The P_basis_set must be lower than the number of P points.\n'
                             f'Got {value} and {info.data.get("P")}')
        return value


class RMGSpeciesConstraints(BaseModel):
    """
    A class for validating input.RMG.species_constraints arguments
    """
    allowed: List[str] = ['input species', 'seed mechanisms', 'reaction libraries']
    max_C_atoms: Annotated[int, Field(ge=0)]
    max_O_atoms: Annotated[int, Field(ge=0)]
    max_N_atoms: Annotated[int, Field(ge=0)]
    max_Si_atoms: Annotated[int, Field(ge=0)]
    max_S_atoms: Annotated[int, Field(ge=0)]
    max_heavy_atoms: Annotated[int, Field(ge=0)]
    max_radical_electrons: Annotated[int, Field(ge=0)]
    max_singlet_carbenes: Annotated[int, Field(ge=0)] = 1
    max_carbene_radicals: Annotated[int, Field(ge=0)] = 0
    allow_singlet_O2: bool = True

    class Config:
        extra = "forbid"

    @field_validator('allowed')
    @classmethod
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
    options: Optional[T3Options] = Field(default_factory=T3Options)
    sensitivity: Optional[T3Sensitivity] = None
    uncertainty: Optional[T3Uncertainty] = None

    class Config:
        extra = "forbid"


class RMG(BaseModel):
    """
    A class for validating input.RMG arguments
    """
    rmg_execution_type: Optional[str] = None
    memory: Optional[Annotated[int, Field(ge=0)]] = None
    cpus: Optional[Annotated[int, Field(gt=0)]] = None
    database: RMGDatabase
    reactors: List[RMGReactor]
    species: List[RMGSpecies]
    model: RMGModel
    pdep: Optional[RMGPDep] = None
    options: Optional[RMGOptions] = Field(default_factory=RMGOptions)
    species_constraints: Optional[RMGSpeciesConstraints] = None

    class Config:
        extra = "forbid"

    @field_validator('database')
    @classmethod
    def check_database(cls, value):
        """RMG.database validator"""
        if value is None or not value:
            raise ValueError('RMG database must be specified')
        return value

    @field_validator('reactors')
    @classmethod
    def check_reactors(cls, value):
        """RMG.reactors validator"""
        if value is None or not value:
            raise ValueError('RMG reactors must be specified')
        return value

    @field_validator('species')
    @classmethod
    def check_species(cls, value):
        """RMG.species validator"""
        if value is None or not value:
            raise ValueError('RMG species must be specified')
        return value

    @field_validator('model')
    @classmethod
    def check_model(cls, value):
        """RMG.model validator"""
        if value is None or not value:
            raise ValueError('RMG model must be specified')
        return value

    @field_validator('pdep')
    @classmethod
    def check_pdep_only_if_gas_phase(cls, value, info: ValidationInfo):
        """RMG.pdep validator"""
        if value is not None and 'reactors' in info.data and info.data['reactors'] is not None:
            reactor_types = set([reactor.type for reactor in info.data['reactors']])
            if value is not None and not any(['gas' in reactor for reactor in reactor_types]):
                raise ValueError(f'A pdep section can only be specified for gas phase reactors, got: {reactor_types}')
        return value

    @model_validator(mode='after')
    def check_species_and_reactors(self) -> 'RMG':
        if self.reactors and self.species:
            reactor_types = {reactor.type for reactor in self.reactors}
            balance_species = [s.label for s in self.species if s.balance]
            solvent_species = [s.label for s in self.species if s.solvent]
            gas_reactors = any('gas' in rt for rt in reactor_types)
            liquid_reactors = any('liquid' in rt for rt in reactor_types)
            if liquid_reactors and balance_species:
                raise ValueError(f'A species ({balance_species}) cannot be set as a balance species '
                                 f'if liquid phase reactors are defined.')
            if gas_reactors and solvent_species:
                raise ValueError(f'A species ({solvent_species}) cannot be set as a solvent '
                                 f'if gas phase reactors are defined.')
            if len(balance_species) > 1:
                raise ValueError(f'Only a single species may be defined as balance,\ngot: {balance_species}')
            if len(solvent_species) > 1:
                raise ValueError(f'Only a single species may be defined as a solvent,\ngot: {solvent_species}')
            species_labels = {s.label for s in self.species}
            for reactor in self.reactors:
                if reactor.termination_conversion:
                    for term_label in reactor.termination_conversion.keys():
                        if term_label not in species_labels:
                            raise ValueError(f'No species with label "{term_label}" was defined.')
        return self


class QM(BaseModel):
    """
    A class for validating input.QM arguments
    """
    adapter: str = 'ARC'
    species: list = Field(default_factory=list)
    reactions: list = Field(default_factory=list)

    class Config:
        extra = "allow"

    @field_validator('adapter')
    @classmethod
    def check_adapter(cls, value):
        """QM.adapter validator"""
        supported_qm_adapters = ['ARC']
        if value not in supported_qm_adapters:
            raise ValueError(f'Supported QM adapters are:\n{supported_qm_adapters}\nGot:{value}')
        return value


class InputBase(BaseModel):
    """
    An InputBase class for validating input arguments
    """
    project: Annotated[str, Field(max_length=255)]
    project_directory: Optional[Annotated[str, Field(max_length=255)]] = None
    verbose: Annotated[int, Field(ge=10, le=30, multiple_of=10)] = 20
    t3: Optional[T3] = Field(default_factory=T3)
    rmg: RMG
    qm: QM = Field(default_factory=QM)

    class Config:
        extra = "forbid"

    @model_validator(mode='after')
    def validate_rmg_t3(self) -> 'InputBase':
        """
        InputBase.validate_rmg_t3
        Validates cross-dependencies between RMG and T3 configurations.
        """
        if self.rmg and self.t3:
            if self.t3.uncertainty:
                ua_term_time = self.t3.uncertainty.termination_time
                rmg_reactor_term_times = [r.termination_time for r in self.rmg.reactors]
                if all(t is None for t in rmg_reactor_term_times) \
                        and ua_term_time is None \
                        and self.t3.uncertainty.global_analysis:
                    raise ValueError('If a global uncertainty analysis is requested, a termination time must be '
                                     'specified either under "t3.uncertainty.termination_time" '
                                     'or in at least one RMG reactor.')
            reactor_types = {r.type for r in self.rmg.reactors}
            solvents = [s.label for s in self.rmg.species if s.solvent]
            is_liquid = any('liquid' in rt for rt in reactor_types)
            if is_liquid:
                if not solvents:
                    raise ValueError('One species must be defined as the solvent when using liquid phase reactors.')
                if len(solvents) > 1:
                    raise ValueError(f'Only one solvent can be specified, got: {solvents}')
            else:
                if solvents:
                    raise ValueError(f'No solvent species are allowed for gas phase reactors, got: {solvents}')
            if self.rmg.model and self.rmg.model.core_tolerance \
                    and self.t3.options and self.t3.options.max_T3_iterations:
                core_tol = self.rmg.model.core_tolerance
                max_iter = self.t3.options.max_T3_iterations
                if isinstance(core_tol, list):
                    if len(core_tol) > max_iter:
                        raise ValueError(f'The number of RMG core tolerances ({len(core_tol)}) '
                                         f'cannot be greater than the max number of T3 iterations '
                                         f'({max_iter}).')
        return self
