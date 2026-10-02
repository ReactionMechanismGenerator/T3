"""
t3 schema module
used for input validation
"""

import math
import os
import re
from enum import Enum
from typing import Annotated, Literal

from pydantic import BaseModel, Field, ValidationInfo, field_serializer, field_validator, model_validator

from arc.common import read_yaml_file

from t3.common import (DATA_BASE_PATH, METHOD_MAP, VALID_CHARS,
                       convert_pressure_to_bar, convert_temperature_to_kelvin, convert_time_to_seconds,
                       convert_volume_to_cubic_meters)
from t3.simulate.factory import _registered_simulate_adapters


MAX_RCM_HISTORY_DURATION = 10.0


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
    generate_flux_diagrams: bool = True
    flux_diagrams_with_images: bool = True
    flux_diagram_reactors: int | list[int] | Literal['all'] | None = None
    profiles_adapter: Annotated[str, Field(max_length=255)] = 'RMG'
    collision_violators_thermo: bool = False
    collision_violators_rates: bool = False
    all_core_species: bool = False
    all_core_reactions: bool = False
    fit_missing_GAV: bool = False
    max_T3_iterations: Annotated[int, Field(gt=0)] = 10
    max_RMG_exceptions_allowed: Annotated[int, Field(ge=0)] | None = 10
    max_RMG_walltime: Annotated[str, Field(pattern=r'\d+:\d\d:\d\d:\d\d')] = '00:00:00:00'
    max_T3_walltime: Annotated[str, Field(pattern=r'\d+:\d\d:\d\d:\d\d')] | None = None
    max_rmg_processes: Annotated[int, Field(ge=1)] | None = None
    max_rmg_iterations: Annotated[int, Field(ge=1)] | None = None
    library_name: Annotated[str, Field(max_length=255)] = 'T3lib'
    shared_library_name: Annotated[str, Field(max_length=255)] | None = None
    external_library_path: Annotated[str, Field(max_length=255)] | None = None
    num_sa_per_temperature_range: Annotated[int, Field(ge=1)] = 3
    num_sa_per_pressure_range: Annotated[int, Field(ge=1)] = 3
    num_sa_per_volume_range: Annotated[int, Field(ge=1)] = 3
    num_sa_per_concentration_range: Annotated[int, Field(ge=1)] = 3
    modify_concentration_ranges_together: bool = True
    modify_concentration_ranges_in_reverse: bool = False

    class Config:
        extra = "forbid"

    @model_validator(mode='after')
    def enforce_collision_thermo(self) -> T3Options:
        """
        If collision_violators_rates is True, ensure collision_violators_thermo is also True.
        """
        if self.collision_violators_rates:
            self.collision_violators_thermo = True
        return self

    @field_validator('flux_diagram_reactors', mode='before')
    @classmethod
    def check_flux_diagram_reactors(cls, value):
        """flux_diagram_reactors: 1-based reactor number(s) or 'all'."""
        if value is None or value == 'all':
            return value
        nums = value if isinstance(value, list) else [value]
        for n in nums:
            if isinstance(n, bool) or not isinstance(n, int) or n < 1:
                raise ValueError(f'flux_diagram_reactors must be positive 1-based int(s) or '
                                 f'"all"; got {value!r}')
        return value

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


class IDTCriterionEnum(str, Enum):
    """IDT criterion for determining ignition delay time."""
    max_dOHdt = 'max_dOHdt'
    max_dTdt = 'max_dTdt'
    max_radical_dt = 'max_radical_dt'


class IDTSAMethodEnum(str, Enum):
    """SA method for IDT sensitivity analysis."""
    brute_force = 'brute_force'
    adjoint = 'adjoint'


class TemperatureUnitEnum(str, Enum):
    """Temperature units accepted by versioned experimental IDT files."""
    K = 'K'
    degC = 'degC'


class PressureUnitEnum(str, Enum):
    """Pressure units accepted by versioned experimental IDT files."""
    Pa = 'Pa'
    kPa = 'kPa'
    MPa = 'MPa'
    bar = 'bar'
    atm = 'atm'


class TimeUnitEnum(str, Enum):
    """Time units accepted by versioned experimental IDT files."""
    s = 's'
    ms = 'ms'
    us = 'us'
    micro_s = 'micro-s'


class ExperimentalApparatusEnum(str, Enum):
    """Apparatus models supported for versioned experimental IDT points."""
    shock_tube = 'shock tube'
    rapid_compression_machine = 'rapid compression machine'


class IgnitionTargetEnum(str, Enum):
    """Signals supported as source-defined experimental ignition targets."""
    pressure = 'pressure'
    temperature = 'temperature'
    OH = 'OH'
    OH_star = 'OH*'
    CH = 'CH'
    CH_star = 'CH*'


class IgnitionTypeEnum(str, Enum):
    """Source-defined methods supported for locating ignition on a target trace."""
    derivative_max = 'd/dt max'
    maximum = 'max'
    half_max = '1/2 max'
    derivative_max_extrapolated = 'd/dt max extrapolated'


class ExperimentalIDTRefusalReason(str, Enum):
    """Typed reasons why an otherwise valid experimental point was not scored."""
    unmappable_species = 'unmappable species'
    target_species_absent = 'target species absent'
    ignition_not_resolved = 'ignition not resolved'
    simulation_failed = 'simulation failed'


class ExperimentalTemperature(BaseModel):
    """A temperature with explicit units whose converted value is above 0 K."""
    value: Annotated[float, Field(allow_inf_nan=False)]
    units: TemperatureUnitEnum

    class Config:
        extra = 'forbid'

    @model_validator(mode='after')
    def validate_kelvin(self):
        """Require the temperature converted to Kelvin to be positive."""
        if convert_temperature_to_kelvin(self.value, self.units) <= 0:
            raise ValueError('temperature must be greater than 0 K (Kelvin)')
        return self


class ExperimentalPressure(BaseModel):
    """A pressure with explicit units."""
    value: Annotated[float, Field(gt=0, allow_inf_nan=False)]
    units: PressureUnitEnum

    class Config:
        extra = 'forbid'


class ExperimentalTime(BaseModel):
    """A positive time with explicit units."""
    value: Annotated[float, Field(gt=0, allow_inf_nan=False)]
    units: TimeUnitEnum

    class Config:
        extra = 'forbid'


class ExperimentalUncertainty(BaseModel):
    """A non-negative IDT uncertainty with explicit time units."""
    value: Annotated[float, Field(ge=0, allow_inf_nan=False)]
    units: TimeUnitEnum

    class Config:
        extra = 'forbid'


class ExperimentalCompositionEntry(BaseModel):
    """One SMILES-identified component of an experimental mole-fraction mixture."""
    smiles: Annotated[str, Field(min_length=1)]
    mole_fraction: Annotated[float, Field(gt=0, le=1, allow_inf_nan=False)]

    class Config:
        extra = 'forbid'


class ExperimentalIgnitionDefinition(BaseModel):
    """The source's target signal and rule for defining ignition."""
    target: IgnitionTargetEnum
    type: IgnitionTypeEnum

    class Config:
        extra = 'forbid'


class ExperimentalSourceReference(BaseModel):
    """A DOI and free-text locator for the source record."""
    doi: Annotated[str, Field(min_length=1)]
    record: Annotated[str, Field(min_length=1)]

    class Config:
        extra = 'forbid'


class ExperimentalHistoryTimes(BaseModel):
    """Explicit time samples for an RCM volume history."""
    values: list[float]
    units: Literal['s', 'ms', 'us']

    class Config:
        extra = 'forbid'


class ExperimentalHistoryVolumes(BaseModel):
    """Explicit volume samples for an RCM volume history."""
    values: list[float]
    units: Literal['m3', 'cm3', 'L']

    class Config:
        extra = 'forbid'


class ExperimentalCompressionTime(BaseModel):
    """End-of-compression time on the volume history's time axis."""
    value: float
    units: Literal['s', 'ms', 'us']

    class Config:
        extra = 'forbid'


class ExperimentalVolumeHistory(BaseModel):
    """Piecewise-linear RCM volume history spanning at most 10 s after conversion to SI."""
    time: ExperimentalHistoryTimes
    volume: ExperimentalHistoryVolumes
    compression_time: ExperimentalCompressionTime | None = None

    class Config:
        extra = 'forbid'

    def to_si(self) -> tuple[list[float], list[float]]:
        """Return time and volume samples in seconds and cubic meters."""
        return ([convert_time_to_seconds(value, self.time.units) for value in self.time.values],
                [convert_volume_to_cubic_meters(value, self.volume.units) for value in self.volume.values])

    @property
    def end_of_compression(self) -> float:
        """Use the explicit compression time, otherwise the first minimum-volume time."""
        if self.compression_time is not None:
            return convert_time_to_seconds(self.compression_time.value, self.compression_time.units)
        times, volumes = self.to_si()
        return times[volumes.index(min(volumes))]

    @model_validator(mode='after')
    def validate_history(self):
        """Reject degenerate, nonphysical or numerically unrepresentable histories."""
        times, volumes = self.to_si()
        if len(times) != len(volumes):
            raise ValueError('volume_history time and volume lists must have equal length')
        if len(times) < 2:
            raise ValueError('volume_history must contain at least two points')
        if not all(math.isfinite(value) for value in times + volumes):
            raise ValueError('volume_history values converted to SI must all be finite')
        if any(value <= 0 for value in volumes):
            raise ValueError('volume_history volumes converted to m3 must be greater than zero')
        intervals = [later - earlier for earlier, later in zip(times, times[1:])]
        if any(interval <= 0 for interval in intervals):
            raise ValueError('volume_history times converted to seconds must be strictly increasing')
        duration = times[-1] - times[0]
        if not math.isfinite(duration):
            raise ValueError('volume_history integration horizon in seconds must be finite')
        if duration > MAX_RCM_HISTORY_DURATION:
            raise ValueError('volume_history duration converted to seconds must be no greater than 10 s')
        if (not all(math.isfinite(interval) for interval in intervals)
                or any(not math.isfinite((later - earlier) / interval)
                       for earlier, later, interval in zip(volumes, volumes[1:], intervals))):
            raise ValueError('volume_history intervals and volume slopes in SI must be finite')
        compression_time = self.end_of_compression
        if not math.isfinite(compression_time):
            raise ValueError('compression_time converted to seconds must be finite')
        if not times[0] <= compression_time <= times[-1]:
            raise ValueError('compression_time must be within the history time range')
        return self


class ExperimentalIDTPoint(BaseModel):
    """One version-1 point with an IDT up to 10 s and optional RCM history spanning up to 10 s.

    History-driven integration ends no later than 10 s after end of compression,
    so its total elapsed duration is bounded by 20 s. History-free points are unchanged.
    """
    temperature: ExperimentalTemperature
    pressure: ExperimentalPressure
    composition: list[ExperimentalCompositionEntry]
    apparatus: ExperimentalApparatusEnum
    ignition_definition: ExperimentalIgnitionDefinition
    idt: ExperimentalTime
    uncertainty: ExperimentalUncertainty | None = None
    source: ExperimentalSourceReference
    volume_history: ExperimentalVolumeHistory | None = None
    initial_temperature: ExperimentalTemperature | None = None
    initial_pressure: ExperimentalPressure | None = None

    class Config:
        extra = 'forbid'

    @property
    def volume_history_horizon(self) -> float:
        """Cover the history and ten measured delays, capped at 10 s after compression."""
        times, _ = self.volume_history.to_si()
        return max(times[-1], self.volume_history.end_of_compression
                   + min(MAX_RCM_HISTORY_DURATION,
                         10.0 * convert_time_to_seconds(self.idt.value, self.idt.units)))

    @model_validator(mode='after')
    def validate_volume_history_state(self):
        """Keep existing post-compression state fields unchanged for all RCM points."""
        if self.volume_history is None:
            if self.initial_temperature is not None or self.initial_pressure is not None:
                raise ValueError('initial_temperature and initial_pressure require a volume_history')
            return self
        if self.apparatus != ExperimentalApparatusEnum.rapid_compression_machine:
            raise ValueError('volume_history is only allowed for a rapid compression machine')
        if self.initial_temperature is None or self.initial_pressure is None:
            raise ValueError('volume_history requires initial_temperature and initial_pressure')
        initial_pressure = convert_pressure_to_bar(self.initial_pressure.value, self.initial_pressure.units) * 1e5
        if not math.isfinite(initial_pressure) or initial_pressure <= 0:
            raise ValueError('initial_pressure converted to Pa must be finite and greater than zero')
        times, _ = self.volume_history.to_si()
        compression_time = self.volume_history.end_of_compression
        horizon = self.volume_history_horizon
        if not math.isfinite(horizon) or not math.isfinite(horizon - times[0]):
            raise ValueError('volume_history integration horizon in seconds must be finite')
        idt_seconds = convert_time_to_seconds(self.idt.value, self.idt.units)
        window = min(MAX_RCM_HISTORY_DURATION, 10.0 * idt_seconds)
        if (horizon - compression_time < 0.999 * window
                and times[-1] - compression_time < window):
            raise ValueError('volume_history horizon does not preserve the post-compression window')
        return self

    @model_validator(mode='after')
    def validate_idt_horizon(self):
        """Keep the per-point integration horizon finite, positive and bounded at 10 seconds.

        The lower bound is checked here, on the *converted* value, rather than on
        ``ExperimentalTime.value``: a positive subnormal survives ``Field(gt=0)`` and
        then underflows to exactly ``0.0`` once the unit factor is applied, which would
        reach the ``simulated_idt / experimental_idt`` comparison as a zero denominator.
        """
        idt_seconds = convert_time_to_seconds(self.idt.value, self.idt.units)
        if not math.isfinite(idt_seconds) or idt_seconds <= 0 or idt_seconds > 10.0:
            raise ValueError('IDT converted to seconds must be finite, greater than 0 s, '
                             'and no greater than 10 s')
        return self

    @field_validator('composition')
    @classmethod
    def validate_composition(cls, value):
        """Require a non-empty, unique, normalized mole-fraction composition."""
        if not value:
            raise ValueError('composition must contain at least one species')
        smiles = [entry.smiles for entry in value]
        if len(set(smiles)) != len(smiles):
            raise ValueError('composition SMILES entries must be unique')
        total = sum(entry.mole_fraction for entry in value)
        if abs(total - 1.0) > 1e-6:
            raise ValueError(f'composition mole fractions must sum to 1.0, got {total}')
        return value


class ExperimentalIDTFile(BaseModel):
    """Version-1 per-point experimental ignition-delay input file."""
    version: Literal[1]
    points: list[ExperimentalIDTPoint]

    class Config:
        extra = 'forbid'


class T3Sensitivity(BaseModel):
    """
    A class for validating input.T3.sensitivity arguments
    """
    adapter: Annotated[str, Field(max_length=255)] | None = 'CanteraConstantTP'
    atol: Annotated[float, Field(gt=0, lt=1e-1)] = 1e-6
    rtol: Annotated[float, Field(gt=0, lt=1e-1)] = 1e-4
    global_observables: list[Annotated[str, Field(min_length=2, max_length=3)]] | None = None
    SA_threshold: Annotated[float, Field(gt=0, lt=0.5)] = 0.01
    max_sa_workers: Annotated[int, Field(ge=1)] = 24
    pdep_SA_threshold: Annotated[float, Field(gt=0, lt=0.5)] | None = 0.001
    pdep_min_delta_ln_k: Annotated[float, Field(gt=0, lt=1)] = 1e-3
    ME_methods: list[Annotated[str, Field(min_length=2, max_length=3)]] = ['CSE', 'MSC']
    # `strict=True` on these two is deliberate, and is about `bool` rather than about types in
    # general: `bool` is a subclass of `int`, so without it `pdep_QM_max_networks: true` in a YAML
    # input validates happily as 1 and silently caps the run at a single network per iteration.
    pdep_QM_max_transition_states: Annotated[int, Field(gt=0, strict=True)] | None = None
    pdep_QM_max_networks: Annotated[int, Field(gt=0, strict=True)] | None = None
    top_SA_species: Annotated[int, Field(ge=0)] = 10
    top_SA_reactions: Annotated[int, Field(ge=0)] = 10
    T_list: list[Annotated[float, Field(gt=0)]] | None = None
    P_list: list[Annotated[float, Field(gt=0)]] | None = None
    idt_criterion: IDTCriterionEnum = IDTCriterionEnum.max_dOHdt
    idt_sa_method: IDTSAMethodEnum = IDTSAMethodEnum.brute_force
    delta_h: Annotated[float, Field(gt=0)] = 0.1
    delta_k: Annotated[float, Field(gt=0, lt=1)] = 0.05
    adaptive_perturbation: bool = False
    save_sa_yaml: bool = True
    experimental_idt_path: str | None = None

    class Config:
        extra = "forbid"

    @field_serializer('idt_criterion', 'idt_sa_method')
    @classmethod
    def serialize_enums(cls, v):
        """Serialize enum fields to plain strings for YAML and logging compatibility."""
        return v.value if isinstance(v, Enum) else v

    @field_validator('adapter')
    @classmethod
    def check_adapter(cls, value):
        """T3Sensitivity.adapter validator"""
        if value is not None and value not in _registered_simulate_adapters.keys():
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
        """T3Sensitivity.ME_methods validator.

        Accepts any casing and returns the canonical ``t3.common.METHOD_MAP`` key. The
        case-insensitive acceptance is deliberate and matches ``global_observables`` above, whose
        consumers do read case-insensitively; ``ME_methods``' consumers do not. The string is
        looked up in the uppercase-keyed ``METHOD_MAP`` when the Arkane input's ``method = ...``
        line is rewritten, and ``t3/main.py`` also uses it verbatim as a directory name, as the
        ``method`` recorded in the SA cache sidecar, and as ``requested_me_methods`` provenance.
        Accepting a spelling here without canonicalizing it therefore did not configure anything
        -- it raised a bare ``KeyError: 'cse'`` from inside the writer, which the
        ``(OSError, ValueError)`` handler around that call does not catch.
        """
        if value is None or not value:
            raise ValueError('The ME methods argument cannot be None or empty.')
        canonical_methods = {method.lower(): method for method in METHOD_MAP}
        normalized = list()
        for entry in value:
            method = canonical_methods.get(entry.lower())
            if method is None:
                raise ValueError(f'The ME methods list must contain a combination of '
                                 f'{sorted(METHOD_MAP)}, got {entry} in {value}')
            # Compared after canonicalization: ['CSE', 'cse'] is one method spelled two ways.
            if method in normalized:
                raise ValueError(f'The ME methods list must not contain repetitions, got {value}')
            normalized.append(method)
        return normalized


class T3Uncertainty(BaseModel):
    """
    A class for validating input.T3.uncertainty arguments
    """
    adapter: Annotated[str, Field(max_length=255)] | None = None
    local_analysis: bool = False
    global_analysis: bool = False
    correlated: bool = True
    local_number: Annotated[int, Field(gt=0)] = 10
    global_number: Annotated[int, Field(gt=0)] = 5
    termination_time: Annotated[str, Field(pattern=r'\d+:\d\d:\d\d:\d\d')] | None = None
    PCE_run_time: Annotated[int, Field(gt=0)] = 1800
    PCE_error_tolerance: Annotated[float, Field(gt=0)] | None = None
    PCE_max_evals: Annotated[int, Field(gt=0)] | None = None
    logx: bool = False

    class Config:
        extra = "forbid"


class RMGDatabase(BaseModel):
    """
    A class for validating input.RMG.database arguments
    """
    thermo_libraries: list[str] | None = None
    kinetics_libraries: list[str] | None = None
    chemistry_sets: list[str] | None = None
    use_low_credence_libraries: bool = False
    transport_libraries: list[str] = ['OneDMinN2', 'PrimaryTransportLibrary', 'NOx2018', 'GRI-Mech']
    seed_mechanisms: list[str] = list()
    kinetics_depositories: list[str] | str = 'default'
    kinetics_families: str | list[str] = 'default'
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


class SpeciesRoleEnum(str, Enum):
    """Allowed RMGSpecies role values for fuel/oxidizer/diluent mixtures."""
    fuel = 'fuel'
    oxidizer = 'oxidizer'
    diluent = 'diluent'


class RMGSpecies(BaseModel):
    """
    A class for validating input.RMG.species arguments
    """
    label: str
    concentration: Annotated[float, Field(ge=0)] | tuple[Annotated[float, Field(ge=0)], Annotated[float, Field(ge=0)]] = 0
    role: SpeciesRoleEnum | None = None
    equivalence_ratios: list[Annotated[float, Field(gt=0)]] | None = None
    oxidizer_fraction: Annotated[float, Field(gt=0, le=1)] | None = None
    diluent_to_oxidizer_ratio: Annotated[float, Field(gt=0)] | None = None
    smiles: str | None = None
    inchi: str | None = None
    adjlist: str | None = None
    reactive: bool = True
    observable: bool = False
    SA_observable: bool = False
    UA_observable: bool = False
    constant: bool = False
    balance: bool = False
    solvent: bool = False
    xyz: list[dict | str] | dict | str | None = None
    seed_all_rads: list[RadicalTypeEnum] | None = None

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

    @model_validator(mode='after')
    def check_role_consistency(self) -> RMGSpecies:
        """
        Cross-field role validation:
        - equivalence_ratios may only be set on a fuel species.
        - oxidizer_fraction may only be set on an oxidizer species.
        - diluent_to_oxidizer_ratio may only be set on a diluent species.
        - A fuel species must declare equivalence_ratios (the φ-driven sweep).
        """
        if self.equivalence_ratios is not None and self.role != SpeciesRoleEnum.fuel:
            raise ValueError(f"equivalence_ratios may only be set on a species with role='fuel'. "
                             f"Got role={self.role!r} for {self.label!r}.")
        if self.oxidizer_fraction is not None and self.role != SpeciesRoleEnum.oxidizer:
            raise ValueError(f"oxidizer_fraction may only be set on a species with role='oxidizer'. "
                             f"Got role={self.role!r} for {self.label!r}.")
        if self.diluent_to_oxidizer_ratio is not None and self.role != SpeciesRoleEnum.diluent:
            raise ValueError(f"diluent_to_oxidizer_ratio may only be set on a species with role='diluent'. "
                             f"Got role={self.role!r} for {self.label!r}.")
        if self.role == SpeciesRoleEnum.fuel and not self.equivalence_ratios:
            raise ValueError(f"A fuel species must declare a non-empty equivalence_ratios list. "
                             f"Got equivalence_ratios={self.equivalence_ratios!r} for {self.label!r}.")
        return self

    @field_serializer('role')
    @classmethod
    def serialize_role(cls, v):
        """Serialize the role enum to a plain string for YAML and logging compatibility."""
        return v.value if isinstance(v, Enum) else v

    @field_serializer('seed_all_rads')
    @classmethod
    def serialize_seed_all_rads(cls, v):
        """Serialize the RadicalTypeEnum list to plain strings so the schema dump stays
        yaml.safe_dump-able for write_t3_input_file."""
        return [x.value if isinstance(x, Enum) else x for x in v] if v else v

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


class IDTModeEnum(str, Enum):
    """
    Sweep mode used by the IDT adapter when several T / P / φ values are given.
    """
    matrix = 'matrix'  # full T_list × P_list × φ_list cartesian product
    row = 'row'        # zip(T_list, P_list, φ_list) — all must have the same length


class RMGReactor(BaseModel):
    """
    A class for validating input.RMG.reactors arguments
    """
    type: str
    T: Annotated[float, Field(gt=0)] | list[Annotated[float, Field(gt=0)]]
    P: Annotated[float, Field(gt=0)] | list[Annotated[float, Field(gt=0)]] | None = None
    V: Annotated[float, Field(gt=0)] | list[Annotated[float, Field(gt=0)]] | None = None
    termination_conversion: dict[str, Annotated[float, Field(gt=0, lt=1)]] | None = None
    termination_time: tuple[Annotated[float, Field(gt=0)], TerminationTimeEnum] | None = None
    termination_rate_ratio: Annotated[float, Field(gt=0, lt=1)] | None = None
    conditions_per_iteration: Annotated[int, Field(gt=0)] = 12
    idt_mode: IDTModeEnum = IDTModeEnum.matrix

    class Config:
        extra = "forbid"

    @field_serializer('idt_mode')
    @classmethod
    def serialize_idt_mode(cls, v):
        """Serialize the idt_mode enum to a plain string for YAML and logging compatibility."""
        return v.value if isinstance(v, Enum) else v

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
        """RMGReactor.T validator. Lists must have at least 2 entries (min/max range, or explicit row points)."""
        if isinstance(value, list) and len(value) < 2:
            raise ValueError(f'When specifying the temperature as a list, at least two values are required,\n'
                             f'got {len(value)} values: {value}.')
        return value

    @field_validator('P')
    @classmethod
    def check_p(cls, value, info: ValidationInfo):
        """RMGReactor.P validator"""
        if isinstance(value, list) and len(value) < 2:
            raise ValueError(f'When specifying the pressure as a list, at least two values are required,\n'
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
        if isinstance(value, list) and len(value) < 2:
            raise ValueError(f'When specifying the volume as a list, at least two values are required,\n'
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
    core_tolerance: Annotated[float, Field(gt=0, lt=1)] | list[Annotated[float, Field(gt=0, lt=1)]]
    atol: Annotated[float, Field(gt=0, lt=1e-1)] = 1e-16
    rtol: Annotated[float, Field(gt=0, lt=1e-1)] = 1e-8
    # filtering:
    filter_reactions: bool = True
    filter_threshold: Annotated[float, Field(gt=0)] | Annotated[int, Field(gt=0)] = 1e+8
    # pruning:
    tolerance_interrupt_simulation: Annotated[float, Field(gt=0)] | list[Annotated[float, Field(gt=0)]] | None = None
    min_core_size_for_prune: Annotated[int, Field(gt=0)] | None = None
    min_species_exist_iterations_for_prune: Annotated[int, Field(gt=0)] | None = None
    tolerance_keep_in_edge: Annotated[float, Field(gt=0)] | None = None
    maximum_edge_species: Annotated[int, Field(gt=0)] | None = None
    tolerance_thermo_keep_species_in_edge: Annotated[float, Field(gt=0)] | None = None
    # staging:
    max_num_species: Annotated[int, Field(gt=0)] | None = None
    # dynamics:
    tolerance_move_edge_reaction_to_core: Annotated[float, Field(gt=0)] | None = None
    tolerance_move_edge_reaction_to_core_interrupt: Annotated[float, Field(gt=0)] | None = None
    dynamics_time_scale: tuple | None = None
    # multiple_objects:
    max_num_objs_per_iter: Annotated[int, Field(gt=0)] = 1
    terminate_at_max_objects: bool = False
    # misc:
    ignore_overall_flux_criterion: bool | None = None
    tolerance_branch_reaction_to_core: Annotated[float, Field(gt=0)] | None = None
    branching_index: Annotated[float, Field(gt=0)] | None = None
    branching_ratio_max: Annotated[float, Field(gt=0)] | None = None
    # surface algorithm
    tolerance_move_edge_reaction_to_surface: Annotated[float, Field(gt=0)] | None = None
    tolerance_move_surface_species_to_core: Annotated[float, Field(gt=0)] | None = None
    tolerance_move_surface_reaction_to_core: Annotated[float, Field(gt=0)] | None = None
    tolerance_move_edge_reaction_to_surface_interrupt: Annotated[float, Field(gt=0)] | None = None

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
    def check_tolerance_interrupt_simulation(self) -> RMGModel:
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
    T: list[Annotated[int, Field(gt=0)] | Annotated[float, Field(gt=0)]] = [300, 2500, 10]
    P: list[Annotated[int, Field(gt=0)] | Annotated[float, Field(gt=0)]] = [0.01, 100, 10]
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
    allowed: list[str] = ['input species', 'seed mechanisms', 'reaction libraries']
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
    options: T3Options | None = Field(default_factory=T3Options)
    sensitivity: T3Sensitivity | None = None
    uncertainty: T3Uncertainty | None = None

    class Config:
        extra = "forbid"


class RMG(BaseModel):
    """
    A class for validating input.RMG arguments
    """
    rmg_execution_type: str | None = None
    memory: Annotated[int, Field(ge=0)] | None = None
    cpus: Annotated[int, Field(gt=0)] | None = None
    database: RMGDatabase
    reactors: list[RMGReactor]
    species: list[RMGSpecies]
    model: RMGModel
    pdep: RMGPDep | None = None
    options: RMGOptions | None = Field(default_factory=RMGOptions)
    species_constraints: RMGSpeciesConstraints | None = None

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
    def check_pdep_has_an_unreactive_species(self) -> RMG:
        """
        A pdep run requires at least one unreactive (inert) species to serve as the bath gas.

        RMG identifies the bath gas as every unreactive core species -- ``rmgpy/rmg/pdep.py:856-858``
        builds ``[spec for spec in core.species if not spec.reactive]`` and then asserts
        ``len(bath_gas) > 0`` ('No unreactive species to identify as bath gas') inside every
        pressure-dependent network update. An input with a pdep block and no ``reactive: false``
        species therefore cannot fail at input parsing: it dies deep inside network generation,
        potentially hours into the run, with an AssertionError that names nothing the user typed.
        Refusing it here converts that into an immediate, actionable input error.
        (Originally proposed in PR #60.)
        """
        if self.pdep is not None and self.species and not any(not spec.reactive for spec in self.species):
            raise ValueError(
                "A pdep section requires at least one unreactive species to serve as the bath gas: "
                "RMG identifies the bath gas as every core species declared with 'reactive: false' "
                "(rmgpy/rmg/pdep.py:856-858, which asserts that at least one exists inside every "
                "network update -- deep into the run, long after the input was accepted). Mark an "
                "inert species (e.g. N2, Ar, He) with 'reactive: false', or remove the pdep section.")
        return self

    @model_validator(mode='after')
    def check_species_and_reactors(self) -> RMG:
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


# Only def2 basis sets and the wb97xd functional must be written hyphen-free: ARC's cached
# frequency scale factors (data/freq_scale_factors.yml) are keyed without hyphens, so a dashed
# form silently falls back to a Truhlar fit, logs "Could not determine software for job type",
# and makes Gaussian reject the route line. Genuine Dunning names (cc-pvtz-f12) keep their dashes,
# so this checks only the two families that are actually affected.
# Source: Vault/Code/ARC/Canonical Levels of Theory.md
_DASHED_LEVEL_PATTERNS = (re.compile(r'\bwb97x-d\b', re.IGNORECASE),
                          re.compile(r'\bdef2-', re.IGNORECASE))


def _refuse_dashed_level(value: str, field_name: str) -> str:
    """
    Refuse a level of theory written with hyphens where ARC requires none.

    Args:
        value (str): The level of theory.
        field_name (str): The field being validated, named in the error.

    Returns:
        str: The unchanged value.

    Raises:
        ValueError: If the value contains a dashed def2 basis or a dashed wb97xd functional.
    """
    for pattern in _DASHED_LEVEL_PATTERNS:
        if pattern.search(value):
            raise ValueError(
                f"The '{field_name}' level of theory must be written undashed for def2 basis sets "
                f"and the wb97xd functional: got {value!r}. Write e.g. 'wb97xd/def2tzvp'. A dashed "
                f"form makes ARC miss the cached frequency scale factor and makes Gaussian reject "
                f"the route line. Genuine Dunning names such as 'cc-pvtz-f12' keep their dashes.")
    return value


class LevelOfTheoryMixin:
    """Shared level-of-theory validation for T3's QM sections."""

    # ``sp_level`` is intentionally exempt from the dashed guard: ARC's software-specific
    # frequency-scale keys include ``wb97xd/def2tzvp, software: gaussian``,
    # ``wb97xd3/def2-tzvp, software: qchem``, and
    # ``wb97xd3/def2-tzvpd, software: orca``. Single-point levels are where QChem/Orca
    # methods such as DLPNO are used, and no frequency is scaled for a single point.
    _LEVEL_PAIRING_FIELDS = ('opt_level', 'freq_level', 'irc_level', 'scan_level')
    _DASHED_LEVEL_FIELDS = _LEVEL_PAIRING_FIELDS + ('level_of_theory',)

    @model_validator(mode='after')
    def validate_levels(self):
        """Validate declared fields and only user-supplied extras on permissive QM models."""
        extras = self.model_extra or {}
        declared_fields = type(self).model_fields
        levels = {field_name: getattr(self, field_name)
                  for field_name in self._LEVEL_PAIRING_FIELDS if field_name in declared_fields}
        levels.update({field_name: extras[field_name]
                       for field_name in self._DASHED_LEVEL_FIELDS if field_name in extras})

        for field_name, value in levels.items():
            if isinstance(value, str):
                _refuse_dashed_level(value, field_name)

        if ('freq_level' in levels and 'opt_level' in levels
                and levels['freq_level'] != levels['opt_level']):
            raise ValueError(f"'freq_level' must equal 'opt_level' so frequencies are evaluated at "
                             f"a real minimum of the same surface. Got freq_level="
                             f"{levels['freq_level']!r}, opt_level={levels['opt_level']!r}.")
        if ('scan_level' in levels and 'freq_level' in levels
                and levels['scan_level'] != levels['freq_level']):
            raise ValueError(f"'scan_level' must equal 'freq_level' so rotors project out "
                             f"correctly. Got scan_level={levels['scan_level']!r}, freq_level="
                             f"{levels['freq_level']!r}.")
        if ('irc_level' in levels and 'opt_level' in levels
                and levels['irc_level'] != levels['opt_level']):
            raise ValueError(f"'irc_level' must equal 'opt_level'. Got irc_level="
                             f"{levels['irc_level']!r}, opt_level={levels['opt_level']!r}.")
        return self


class QM(LevelOfTheoryMixin, BaseModel):
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


class PESStrictSection(BaseModel):
    """
    Base for every section of the standalone PES exploration loop's input file.

    ``extra = "forbid"`` on ``PESLoopConfig`` alone only rejects an unknown TOP-LEVEL key: pydantic
    does not propagate a model's config into the nested models its fields declare. Without this
    base, a typo inside a section -- ``qm.max_transition_state_per_round`` for
    ``max_transition_states_per_round`` -- is silently discarded, the run proceeds on the schema
    default, and the user is never told their setting did nothing. Every section forbids extras
    here rather than repeating the class-level config four times, so a section added later cannot
    forget to.
    """

    class Config:
        extra = "forbid"


class PESSection(PESStrictSection):
    """
    A class for validating input.pes arguments of the standalone PES exploration loop.
    """
    network: Annotated[str, Field(min_length=1)]
    source: list[str]
    method: str = 'MSC'
    # Required, and required to be non-empty. There is no bath gas that is right by default for an
    # arbitrary network, and PDepExplorerConfig refuses a config without one -- but it refuses it
    # from inside run_pes_loop, by which time the CLI has already created the project directory,
    # the log file and round_0/. Refusing it here turns the likeliest first-run mistake on this
    # input file into an immediate, actionable validation error.
    bath_gas: dict
    explore_tol: Annotated[float, Field(gt=0)] | None = None
    energy_tol: Annotated[float, Field(gt=0)] | None = None
    flux_tol: Annotated[float, Field(gt=0)] | None = None
    maximum_radical_electrons: Annotated[int, Field(gt=0, strict=True)] | None = None
    # Explorer runtime is unbounded, so this is load-bearing rather than decorative: without it a
    # single pathological network parks the loop forever.
    #
    # strict=True for the same reason max_rounds uses it: bool is a subclass of int, so
    # `timeout: true` in a YAML input otherwise validates as 1.0 -- a one-second explorer budget
    # that fails every network, arrived at by a typo. allow_inf_nan=False because +inf passes
    # `gt=0` and reinstates exactly the unbounded runtime this field exists to bound. Both are
    # refused by PDepExplorerConfig too, but only from INSIDE run_pes_loop, after the CLI has
    # created the project directory, the log file and round_0/.
    timeout: Annotated[float, Field(gt=0, strict=True, allow_inf_nan=False)] = 7200.0

    @field_validator('source')
    @classmethod
    def check_source(cls, value):
        """PESSection.source validator.

        Arkane's explorer resolves ``source`` from ``species_dict`` only and accepts a
        unimolecular or bimolecular entry channel -- never three or more, and never a transition
        state. Refusing here beats a failure deep inside an Arkane run.
        """
        if not 1 <= len(value) <= 2:
            raise ValueError(
                f"The PES 'source' must name 1 or 2 species -- 1 for a unimolecular well, 2 for a "
                f"bimolecular entry channel (A + B). Arkane's explorer accepts nothing else. "
                f"Got {len(value)}: {value}.")
        return value

    @field_validator('bath_gas')
    @classmethod
    def check_bath_gas(cls, value):
        """PESSection.bath_gas validator.

        An empty mapping passes the ``dict`` type check but is exactly as useless as a missing one:
        ``t3.pdep.explorer.config.PDepExplorerConfig`` raises on both, deep inside the loop.
        """
        if not value:
            raise ValueError(
                "The PES 'bath_gas' must be a non-empty mapping of species labels to mole "
                "fractions (e.g. {'N2': 1.0}). An Arkane explorer input file cannot be written "
                "without one, and there is no default that is right for an arbitrary network.")
        return value

    @field_validator('method')
    @classmethod
    def check_method(cls, value):
        """PESSection.method validator."""
        if value not in ('CSE', 'RS', 'MSC'):
            raise ValueError(f"The PES method must be either 'CSE', 'RS', or 'MSC'.\nGot: {value}")
        return value


class PESQMSection(LevelOfTheoryMixin, PESStrictSection):
    """
    A class for validating input.qm arguments of the standalone PES exploration loop.

    Level defaults are the proven ones from Vault/Code/ARC/Canonical Levels of Theory.md, not
    ARC's repo defaults.
    """
    opt_level: str = 'wb97xd/def2tzvp'
    freq_level: str = 'wb97xd/def2tzvp'
    sp_level: str = 'wb97xd/def2tzvp'
    irc_level: str = 'wb97xd/def2tzvp'
    scan_level: str = 'wb97xd/def2tzvp'
    ts_adapters: list[str] = Field(default_factory=lambda: ['heuristics', 'linear', 'kinbot',
                                                            'goflow', 'rits'])
    rotors: bool = False
    irc: bool = False
    # 'sensitive' ranks the QM candidates by master-equation E0-sensitivity, 'all' keeps file order,
    # and 'distrust' (I-032) ranks by how little the current barrier is trusted from data present
    # BEFORE any QM runs -- a flat energy window over the E0 surface, the barrier's provenance, and
    # RMG's own RateUncertainty variance -- so it reaches the low bimolecular entrance channels whose
    # sensitivity is a structural zero the 'sensitive' screen cannot measure (see t3.pdep.distrust).
    scope: Literal['all', 'sensitive', 'distrust'] = 'sensitive'
    max_transition_states_per_round: Annotated[int, Field(gt=0, strict=True)] = 10
    # The smallest measured ln(k) response that justifies spending QM on a transition state,
    # mirroring T3's in-run ``t3.sensitivity.pdep_min_delta_ln_k`` (same default, same bounds).
    # Applies under the 'sensitive' and 'all' scopes: 'sensitive' ranks and 'all' does not, but
    # neither may queue a transition state whose measured leverage is below this floor -- its capture
    # manifest would then record a coefficient that never justified anything. The 'distrust' scope
    # does not gate on it (its whole point is that a structural-zero sensitivity is meaningless).
    min_delta_ln_k: Annotated[float, Field(gt=0, lt=1)] = 1e-3
    # The flat energy window half-width (kJ/mol) for the 'distrust' scope: a saddle more than this
    # above the lowest saddle on the surface carries negligible flux and is declined. Unused under
    # the other scopes. Default mirrors t3.pdep.distrust.DEFAULT_ENERGY_WINDOW_KJ.
    energy_window_kj: Annotated[float, Field(gt=0)] = 30.0


class PESTerminationSection(PESStrictSection):
    """
    A class for validating input.termination arguments of the standalone PES exploration loop.
    """
    # strict=True because bool is a subclass of int: without it, `max_rounds: true` in a YAML
    # input validates happily as 1 and silently caps the loop at a single round.
    max_rounds: Annotated[int, Field(gt=0, strict=True)] = 5
    stop_when_no_new_ts: bool = True


class PESReuseSection(PESStrictSection):
    """
    A class for validating input.reuse arguments of the standalone PES exploration loop.
    """
    from_t3_projects: list[str] = Field(default_factory=list)


class PESLoopConfig(PESStrictSection):
    """
    A class for validating the standalone PES exploration loop's input file.
    """
    pes: PESSection
    qm: PESQMSection = Field(default_factory=PESQMSection)
    termination: PESTerminationSection = Field(default_factory=PESTerminationSection)
    reuse: PESReuseSection = Field(default_factory=PESReuseSection)


class InputBase(BaseModel):
    """
    An InputBase class for validating input arguments
    """
    project: Annotated[str, Field(max_length=255)]
    project_directory: Annotated[str, Field(max_length=255)] | None = None
    verbose: Annotated[int, Field(ge=10, le=30, multiple_of=10)] = 20
    t3: T3 | None = Field(default_factory=T3)
    rmg: RMG
    qm: QM = Field(default_factory=QM)

    class Config:
        extra = "forbid"

    @model_validator(mode='after')
    def validate_rmg_t3(self) -> InputBase:
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
