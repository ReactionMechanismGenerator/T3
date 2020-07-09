# 1. Simple system, no QM

This tutorial will guide you through setting up T3 to generate
a detailed kinetic model for a relatively simple system,
a H<sub>2</sub>/O<sub>2</sub> flame.
The model will consist of species and reactions for which
thermodynamic properties and rate coefficients, respectively,
are already known from the literature.
It is an exercise in setting up T3 for a system for which
**no** quantum mechanical (QM) calculations are executed.

When using the API, the first directive is to import the main T3 object
(this is not needed when using an input file):


```Python hl_lines="1"
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

t3_object.execute()
```
*(This script is complete, it should run "as is")*


Next we'll define arguments related to RMG. This argument is a Python dictionary which
will be passed to the ``rmg`` argument of the ``T3`` object imported above.
You can name this argument with any legal Python name, here we use ``rmg_args``.
This is the only required primary argument -- this simple example does not
make use of the additional primary arguments, ``t3`` and ``qm``.


```Python hl_lines="3 4 5 6 7 8 9 10 11 12 13 14 15 16 17"
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

t3_object.execute()
```


Here, the ``rmg_args`` argument defines the database libraries to be used by RMG,
the chemical species with respective structure and concentration, the reactor
type, conditions and termination criteria, and the RMG core tolerances to be used
throughout the different T3 iterations.


!!! Note
    The RMG arguments in T3 *generally* correspond to arguments used by the legacy RMG input file.
    Users should refer to the
    <a href="http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/input.html" target="_blank">
    RMG documentation</a> for detailed explanations of all possible arguments.
    This tutorial, though, manages to capture all major syntax differences users should be aware of
    when defining the RMG arguments in T3:
    
    1. RMG arguments are currently in camelCase style (e.g., ``thermoLibraries``), whereas the respective
       arguments in T3 are in lower_case_underscore style (e.g., ``thermo_libraries``).
    2. Species definitions in T3 are more flexible, see the
       [how-to guides](how_to.md) for more details.
    3. Reactor types were re-named for clarity, see the
       [how-to guides](how_to.md) for more details.
    4. ``core_tolerance`` is a list of the RMG ``toleranceMoveToCore`` argument to use per T3 iteration.
       If the number of iterations is larger than the length of ``core_tolerance``, the last entry will be used.
       The `core_tolerance`` argument can also be given as flot, in which case it will be treated as a onc-entry list.
    5. Minor: RMG uses the notation ``reactionLibraries`` while T3 uses ``kinetics_libraries``.

Next, we create the ``T3`` object with all desired arguments.
Here, only the two required arguments, ``project`` and ``rmg`` are given.


```Python hl_lines="19 20"
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

t3_object.execute()
```

Finally, we call the ``execute()`` method of ``T3``.


```Python hl_lines="22"
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

t3_object.execute()
```


Note that the corresponding YAML input file for this tutorial would be:


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















