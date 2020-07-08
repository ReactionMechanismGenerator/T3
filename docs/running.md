# Running T3

## General

Executing T3 can be done via two ways:

- using an input file, or
- using the Python API

Other than subtle differences described here,
both approaches are equivalent as they can define
the same parameters which will be processed via the same routines.

!!! Note
    Most examples in these documentation pages are in a
    **Python API** format rather than a
    **<a href="https://yaml.org/" target="_blank">YAML</a>
    input file** format.

    T3 also has a function for writing a corresponding
    YAML input file after defining parameters via the API,
    see the [How-to guides](how_to.md#save-an-input-file-from-the-api)
    for more details.



## Activate the environment

To use T3, first activate the T3 environment. Type either:

    conda activate t3_env

or, if you have set up the recommended aliases, simply type:

    t3e

## Arguments

T3 has three minor arguments (``project``, ``project_directory``, and ``verbose``)
and three primary arguments (``t3``, ``rmg``, and ``qm``).

The ``project`` argument is required. It is a string representing the T3 project name.

The ``project_directory`` argument is optional. It is a string representing the
path to the local project directory where all the project files are stored.
If not specified, it will be set to the folder in which the input file is located
if T3 is being executed using an input file, or to a respective subfolder with the
project's name under the ``Projects`` folder in the T3 repository.

The ``verbose`` argument is optional. It is an integer representing the logging
level used by T3. Allowed values are: ``10``: debug level (very verbose),
``20``: info level (default), ``30``: warnings and errors only, ``40``: errors
only. Pass ``None`` to this argument to avoid saving a log file. 

The primary arguments specify various options for the different respective packages
(T3, RMG, and QM which currently only supports
<a href="https://reactionmechanismgenerator.github.io/ARC/index.html">ARC</a>).
Of these three, only the ``rmg`` argument is required. The ``qm`` argument must
be specified if QM-based model refinement is desired
(in most cases it is!). The ``t3`` argument contains optional T3-related
directives and should commonly be specified.

The RMG arguments in T3 are written in an underscore_lower_case (snake_case) syntax,
while many are in a camelCase syntax in RMG.
A few RMG arguments have different names altogether in T3. These arguments are:

- `kinetics_libraries`: In the RMG database block, the `kinetics_libraries` argument
  replaces the legacy RMG `reactionLibraries` argument.
- `core_tolerance`: In the RMG model block, the `core_tolerance` argument
  replaces the legacy RMG `toleranceMoveToCore` argument. See the
  [How-to guides](how_to.md#the-rmg-core-tolerances) for more details.
- `conditions_per_iteration`: In the RMG reactors block, the `conditions_per_iteration` argument
  replaces the legacy RMG `nSims` argument.
- Species definitions are different than in RMG, see the
  [How-to guides](how_to.md#species-properties) for more details.
- Reactors definitions are different than in RMG, see the
  [How-to guides](how_to.md#reactors) for more details.


!!! Note
    Some of the RMG default values have been changed in T3, see the
    <a href="https://github.com/ReactionMechanismGenerator/T3/blob/master/t3/schema.py" target="_blank">
    schema for details.


Use the below **reference guide** to learn more about these arguments.


## Reference guide

T3 has several types of reference guides:

- The [tutorials](tutorials/1_no_qm.md) are a great place to start with,
  and provide an excellent basic reference guide. In T3's tutorials you can find
  complete, functioning, and worked-out examples with explanations.
- The <a href="https://github.com/ReactionMechanismGenerator/T3/blob/master/examples/commented/input.yml">
  commented input file</a> in T3's examples shows all available input arguments along with a brief explanation.
- A <a href="https://github.com/ReactionMechanismGenerator/T3/blob/master/t3/schema.py" target="_blank">
  pydentic schema</a> is used to validate the input file,
  and could also be used as a reference for the various allowed arguments.


## Where next?

New users should start learning how to use T3 by reading and executing the
[tutorials](tutorials/1_no_qm.md).

For advanced features and specific examples for solving complex problems, see the
 [how-to guides](how_to.md). 
