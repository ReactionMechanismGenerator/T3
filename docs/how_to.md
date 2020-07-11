# How to guides


## Species properties

The `rmg.species` attribute is a list of dictionaries, each defines
a chemical species. The following are possible keys and corresponding
values for each species dictionary:

- `label` (str, required): The species label.
- `concentration` (float): concentration units are mole fraction for gas phase
  and mol/cm3 for liquid phase. Defaults to `0`.
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
- `observable` (bool): Whether the species should be used as an observable for both SA and UA. Default: `False.`
- `SA_observable` (bool): Whether the species should be used as an observable for SA. Default: `False.`
- `UA_observable` (bool): Whether the species should be used as an observable for UA. Default: `False.`
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
    THe simulation is done using the reactor specified in the corresponding
    [simulation adapter](how_to.md#writing-simulation-adapters) which may
    be different if the simulation adapter is not `RMG`.



## Save an input file from the API

Saving an input file using the Python documentation may come handy
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

Coming soon!


## Pre-QM, or: T3's iteration 0

Sometimes it is desired to conduct thermodynamic and/or rate coefficient calculations in advance,
prior to the main T3 iterations. For example, one might want to compute thermodynamic properties
for all important radicals of a fuel molecule in advance (before RMG generates a model).

This pre-QM computations is called the "iteration 0" of T3 and is achieved by specifying `species`
and/or `reactions` in the QM section of the input. The corresponding computations will be executed,
and the computed thermo-kinetic parameters will be used by RMG in all consecutive T3 iterations.


## The RMG core tolerances

RMG uses the `toleranceMoveToCore` argument to control the size of the generated model
(the "core"). See detailed explanation
<a href="http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/input.html#model-tolerances" target="_blank">
here</a>. The corresponding T3 argument is called `core_tolerance`, and is set under `rmg['model']`.
This `core_tolerance` argument is a list of floats, each will be used in a respective T3
iteration. For example, setting `core_tolerance = [0.02, 0.01, 0.005, 0.001]` will make T3 use a
`toleranceMoveToCore` of `0.02` for running RMG in iteration 1,
and a `toleranceMoveToCore` of `0.001` for running RMG in iteration 4.
All iterations further iterations use the last entry in `core_tolerance`. In the above example,
iterations 5, 6, 7... will use a `toleranceMoveToCore` of `0.001` as well.

This feature is useful since at the first iteration RMG is normally executed without
system-specific knowledge which is provided as thermodynamic properties and rate coefficients
in later T3 iterations. If a RMG is given a too low tolerance in the early iterations
it will likely explore unimportant chemistry. By gradually increasing the tolerance
we allow RMG to hone in on the model.
