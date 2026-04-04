# How-to guides


## Species properties

The `rmg.species` attribute is a list of dictionaries, each defines
a chemical species. The following are possible keys and corresponding
values for each species dictionary:

- `label` (str, required): The species label.
- `concentration` (Union[float, Tuple[float, float]]): Concentration units are mole fraction for gas phase
  and mol/cm3 for liquid phase. Defaults to `0`.
  A concentration range can also be specified (a length-2 tuple).
- `smiles` (str): The
  <a href="https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system" target="_blank">
  SMILES</a> representation.
- `inchi` (str): The
  <a href="https://iupac.org/who-we-are/divisions/division-details/inchi/" target="_blank">
  InChI</a> representation.
- `adjlist` (str): The
  <a href="http://reactionmechanismgenerator.github.io/RMG-Py/reference/molecule/adjlist.html" target="_blank">
  RMG adjacency list</a> representation.
- `reactive` (bool): Whether the species is treated by RMG as reactive. Default: `True`.
- `observable` (bool): Whether the species should be used as an observable for SA. Default: `False.`
- `SA_observable` (bool): Whether the species should be used as an observable for SA. Default: `False.`
- `constant` (bool): Whether the species concentration should remain constant throughout the simulation
  and model generation. Default: `False.`
- `balance` (bool): Whether this is a balance species. Default: `False.`
- `solvent` (bool): Whether this species should be used as the solvent in RMG. Can only be set to `True`
  for liquid phase simulations. Default: `False.`
- `xyz` (list): Optional 3D coordinates for a species. Entries could be either string representation,
  ARC dictionary representation, or a file from which the coordinates could be parsed
  (either an
  <a href="https://en.wikipedia.org/wiki/XYZ_file_format" target="_blank">XYZ file format</a>
  or a supported ESS input/log file).
- `seed_all_rads` (List[str]): The types of radical derivatives to add to the RMG input file
  for the species. Helpful for solving orphan radical issues early on in the model generation.
  Recommended for the main molecule undergoing oxidation/pyrolysis.
  Optional types are: 'radical' for simple R., 'alkoxyl' for RO., 'peroxyl' for ROO.

!!! Note
    Either `smiles`, `inchi`, or `adjlist` must be specified for each species.


## Reactors

T3 includes the common RMG reactors. Note that the reactor names have been changed
to explicitly represent their primary properties. The supported RMG reactors are:

- `gas batch constant T P` which corresponds to the RMG `simpleReactor`
- `liquid batch constant T V` which corresponds to the RMG `liquidReactor`


!!! Note
    The RMG reactors are only used for model generation.
    The simulation is done using the reactor specified in the corresponding
    [simulation adapter](how_to.md#simulation-adapters) which may
    be different if the simulation adapter is not `RMG`.


## Simulation adapters

T3 uses simulation adapters to simulate the generated mechanism and perform
sensitivity analysis. The following adapters are available:

| Adapter Name | Reactor Type | Description |
|---|---|---|
| `CanteraConstantTP` | Constant T, P | Isothermal, isobaric batch reactor. Energy equation disabled. |
| `CanteraConstantHP` | Constant H, P | Adiabatic, constant-pressure batch reactor. |
| `CanteraConstantUV` | Constant U, V | Adiabatic, constant-volume batch reactor. |
| `CanteraPFR` | Plug flow reactor | Isothermal, isobaric PFR (Lagrangian or chain-of-reactors method). |
| `CanteraPFRTProfile` | Plug flow reactor with prescribed `T(z)` | Non-isothermal, isobaric PFR.  A Lagrangian particle is advected through `N_CELLS` axial segments and the gas temperature at each segment is taken from a hard-coded `T(z)` profile (energy equation off). Useful for flow-tube experiments where the wall temperature has been measured. SA is supported and indexed by axial position `z`. |
| `CanteraJSR` | Jet stirred reactor | Isothermal, isobaric, steady-state JSR / PSR / CSTR. Residence time is set from `termination_time`. |
| `RMGConstantTP` | Constant T, P | Uses RMG's built-in solver (runs as a subprocess in the `rmg_env`). |

Specify the adapter in the `t3.sensitivity.adapter` field:

```yaml
t3:
  sensitivity:
    adapter: CanteraConstantTP
```

All Cantera adapters support **species-level sensitivity analysis** and
**global observables** (ignition delay time `IDT`, extinction strain rate `ESR`,
laminar flame speed `SL`).

See [Tutorial 2](tutorials/2_adding_simulate_adapter.md) for how to create
a new Cantera adapter.

### Mapping experimental setups to adapters

If you have an experimental target in mind, use this table to choose the
adapter that reproduces its idealized chemistry:

| Experimental setup | Adapter | Notes |
|---|---|---|
| Shock tube (reflected shock) | `CanteraConstantUV` | Standard treatment: closed, adiabatic, constant-volume batch reactor. Use `IDT` as a global observable for ignition delay measurements. |
| Rapid compression machine (RCM) | `CanteraConstantUV` | Same idealization as a reflected shock once the piston is at top-dead-center. Heat loss is not modeled. |
| Constant-volume bomb | `CanteraConstantUV` | Adiabatic version. |
| Flow reactor (isothermal) | `CanteraPFR` | Lagrangian PFR: a parcel of gas advected at constant T, P. |
| Flow reactor with measured wall temperature | `CanteraPFRTProfile` | Lagrangian PFR with a hard-coded axial temperature profile `T(z)` (energy equation off).  The geometry, the breakpoints, and the profile shape are set by module-level constants in `t3/simulate/cantera_pfr_t_profile.py` (defaults to an 80 cm reactor: a 20 cm raised-cosine ramp 300 → 900 K, a 40 cm 900 K plateau, and a 20 cm raised-cosine ramp 900 → 500 K).  Pressure stays at the inlet value; density evolves with `T`.  Supports SA -- the SA coefficients are returned as a function of axial position `z` (the reactor is steady-state in the lab frame), via an extra `'distance'` key in the `sa_dict`. |
| Flow reactor (adiabatic) | not yet supported | Would require an adiabatic PFR variant of `CanteraPFR`. |
| Closed isothermal/isobaric vessel | `CanteraConstantTP` | Energy equation disabled. |
| Closed adiabatic vessel at constant P | `CanteraConstantHP` | Energy equation enabled. |
| Jet-stirred reactor (JSR / PSR / CSTR) | `CanteraJSR` | Isothermal, isobaric, well-mixed open reactor. The reactor's residence time is taken from the `termination_time` of the corresponding `gas batch constant T P` reactor entry, and the mass flow rate is set so that `mdot = m_reactor / tau`. Reactor volume is controlled by the module-level `VOLUME` constant in `t3/simulate/cantera_jsr.py` (defaults to 100 cm³). |
| Laminar premixed flame | not yet supported as an adapter | Cantera's `FreeFlame` is used internally for the `SL` global observable, but there is no standalone flame adapter. |

If your target is not covered by an existing adapter, see
[Writing simulation adapters](#writing-simulation-adapters) and
[Tutorial 2](tutorials/2_adding_simulate_adapter.md).


## Sensitivity analysis

Sensitivity analysis (SA) identifies which reactions and species have the
largest influence on the observables. T3 uses SA results to decide which
thermodynamic properties and rate coefficients should be refined with QM.

### Configuring SA

SA is requested by including the `t3.sensitivity` block:

```yaml
t3:
  sensitivity:
    adapter: CanteraConstantTP
    SA_threshold: 0.01       # minimum SA coefficient to flag a parameter
    top_SA_species: 10       # max species per observable to refine
    top_SA_reactions: 10     # max reactions per observable to refine
```

### SA observables

Mark species as SA observables in the species list:

```yaml
rmg:
  species:
    - label: OH
      smiles: '[OH]'
      SA_observable: true    # SA will track sensitivity of OH
```

### Global observables

The Cantera adapters support global observables in addition to species concentrations:

```yaml
t3:
  sensitivity:
    adapter: CanteraConstantTP
    global_observables: ['IDT']  # ignition delay time
```

Available global observables:

- `IDT` - Ignition delay time
- `ESR` - Extinction strain rate
- `SL` - Laminar flame speed

### Multi-condition SA

When temperature, pressure, or concentration are specified as ranges, T3 generates
a grid of conditions and runs SA at each point:

```yaml
t3:
  options:
    num_sa_per_temperature_range: 3   # number of T points in the range
    num_sa_per_pressure_range: 3      # number of P points in the range
rmg:
  reactors:
    - type: gas batch constant T P
      T: [800, 1500]    # T range in K
      P: [1, 10]        # P range in bar
```

You can also specify explicit T and P lists for SA:

```yaml
t3:
  sensitivity:
    T_list: [800, 1000, 1200]  # explicit temperatures for SA
    P_list: [1, 10]            # explicit pressures for SA
```

See the [ranged conditions example](examples.md#ranged-conditions) for a full working example.

### The sa_dict format

SA results are stored in a dictionary with this structure:

```python
sa_dict = {
    'time': [np.array(...)],        # list of time arrays, one per condition
    'kinetics': [{                   # list of dicts, one per condition
        'OH': {                      # observable label
            1: np.array(...),        # reaction index -> SA coefficients over time
            5: np.array(...),
        },
    }],
    'thermo': [{                     # list of dicts, one per condition
        'OH': {                      # observable label
            'H2O': np.array(...),    # species label -> SA coefficients over time
        },
    }],
}
```

SA results are saved to `sa_coefficients.yml` in YAML format for persistence
and can be loaded back with `sa_dict_from_yaml()`.

!!! Note "PFR with prescribed `T(z)`: SA is indexed by distance"
    The `CanteraPFRTProfile` adapter is the one exception to the time-axis
    convention.  It is steady-state in the lab frame, so the natural
    independent variable for the SA coefficients is the **axial position
    `z`**, not the residence time of the Lagrangian particle.  Its
    `get_sa_coefficients()` therefore returns the standard `sa_dict` plus
    an extra `'distance'` key — a list of one numpy array per condition
    holding the cell-center positions [m].  `sa_dict['time']` is still
    populated (with the cumulative residence time at each cell) so all of
    T3's existing SA-consuming code keeps working unchanged; code that
    wants to plot SA against `z` should use `sa_dict['distance']`.


## Save an input file from the API

Saving an input file using the Python API may come handy
in many cases, since often it is easier to define parameters using
Python and auto-complete tools rather than hand-typing a YAML format.

To do this, define a `T3` object like you would as if running using the API.
However, instead of executing it, call `write_t3_input_file`.
Here's an example:


```Python
from t3 import T3

rmg_args = {'database': {'thermo_libraries': ['primaryThermoLibrary',
                                              'BurkeH2O2'],
                         'kinetics_libraries': ['BurkeH2O2inN2']},
            'species': [{'label': 'H2',
                         'smiles': '[H][H]',
                         'concentration': 0.67},
                        {'label': 'O2',
                         'smiles': '[O][O]',
                         'concentration': 0.33}],
            'reactors': [{'type': 'gas batch constant T P',
                          'T': 1000,
                          'P': 1,
                          'termination_conversion': {'H2': 0.9},
                          'termination_time': [5, 's']}],
            'model': {'core_tolerance': [0.01, 0.001]}}

t3_object = T3(project='T3_tutorial_1',
               rmg=rmg_args)

t3_object.write_t3_input_file()
```


The corresponding auto-generated YAML input file for the above example
would be:


```YAML
project: T3_tutorial_1
rmg:
  database:
    kinetics_libraries:
    - BurkeH2O2inN2
    thermo_libraries:
    - primaryThermoLibrary
    - BurkeH2O2
  model:
    core_tolerance:
    - 0.01
    - 0.001
  reactors:
  - P: 1.0
    T: 1000.0
    termination_conversion:
      H2: 0.9
    termination_time:
    - 5
    - 's'
    type: gas batch constant T P
  species:
  - concentration: 0.67
    label: H2
    smiles: '[H][H]'
  - concentration: 0.33
    label: O2
    smiles: '[O][O]'
```

The `write_t3_input_file` method accepts two arguments:

path (str, optional):

    The full path for the generated input file, or to the folder
    where this file will be saved under a default name.
    If ``None``, the input file will be saved to the project directory.

all_args (bool, optional):

    Whether to save all arguments in the generated input file
    including all default values). Default: ``False``.


## Writing simulation adapters

T3's Cantera-based adapters share a common base class (`CanteraBase`) that handles
all the heavy lifting: mechanism loading, condition generation, time integration,
SA extraction, and data storage. Creating a new Cantera reactor adapter requires
only defining the Cantera reactor type and how to instantiate it.

See [Tutorial 2](tutorials/2_adding_simulate_adapter.md) for a complete
worked example of adding a constant-volume isothermal reactor.

For non-Cantera adapters, the new class must inherit from the abstract
`SimulateAdapter` class in `T3/t3/simulate/adapter.py` and implement
`set_up()`, `simulate()`, `get_sa_coefficients()`, and `get_idt_by_T()`.


## Pre-QM, or: T3's iteration 0

Sometimes it is desired to conduct thermodynamic and/or rate coefficient calculations in advance,
prior to the main T3 iterations. For example, one might want to compute thermodynamic properties
for all important radicals of a fuel molecule in advance (before RMG generates a model).

This pre-QM computation is called the "iteration 0" of T3 and is achieved by specifying `species`
and/or `reactions` in the QM section of the input. The corresponding computations will be executed,
and the computed thermo-kinetic parameters will be used by RMG in all consecutive T3 iterations.

See the [pre-QM example](examples.md#pre-qm-iteration-0) for a working input file.


## The RMG core tolerances

RMG uses the `toleranceMoveToCore` argument to control the size of the generated model
(the "core"). See detailed explanation
<a href="http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/input.html#model-tolerances" target="_blank">
here</a>. The corresponding T3 argument is called `core_tolerance`, and is set under `rmg['model']`.
This `core_tolerance` argument is a list of floats, each will be used in a respective T3
iteration. For example, setting `core_tolerance = [0.02, 0.01, 0.005, 0.001]` will make T3 use a
`toleranceMoveToCore` of `0.02` for running RMG in iteration 1,
and a `toleranceMoveToCore` of `0.001` for running RMG in iteration 4.
All further iterations use the last entry in `core_tolerance`. In the above example,
iterations 5, 6, 7... will use a `toleranceMoveToCore` of `0.001` as well.

This feature is useful since at the first iteration RMG is normally executed without
system-specific knowledge which is provided as thermodynamic properties and rate coefficients
in later T3 iterations. If RMG is given a too low tolerance in the early iterations
it will likely explore unimportant chemistry. By gradually increasing the tolerance
we allow RMG to hone in on the model.
