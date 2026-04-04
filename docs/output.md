# Output

After running a project, the local project folder will contain the following directory tree
(**bold** represents folders, *italics* represents files):

- *t3.log*: Details of all project execution procedures.
- *t3_restart.yml*: T3 restart file. T3 can resume from this if interrupted.
- **iteration_x**:
    - **RMG**:
        - *input.py*: The RMG input file for iteration x (automatically written by T3).
        - *RMG.log*: The RMG log file for iteration x.
        - *chem_annotated.inp*: The Chemkin mechanism file with annotated reaction sources.
        - *species_dictionary.txt*: Species adjacency list dictionary.
        - **chemkin**: Additional Chemkin files (transport, species).
        - **solver**: RMG simulation profiles (if `save_simulation_profiles` is enabled).
    - **SA** (if sensitivity analysis is enabled):
        - *sa_coefficients.yml*: Sensitivity analysis results in YAML format.
        - **solver**: Raw SA output files (CSV format for RMG adapter).
    - **ARC** (if QM calculations are performed):
        - *restart.yml*: The ARC restart file for iteration x.
        - *arc.log*: The ARC log file for iteration x.
        - **calcs**: All spawned jobs sorted by species name (including transition states).
        - **output**: Additional ARC output files and folders.


## Key output files

### The generated mechanism

The final Chemkin mechanism is in the last iteration's RMG folder:

```
iteration_N/RMG/chem_annotated.inp
iteration_N/RMG/species_dictionary.txt
```

These files can be used directly in Cantera, CHEMKIN, or other tools.

### Sensitivity analysis results

SA results are saved in YAML format:

```
iteration_N/SA/sa_coefficients.yml
```

This file contains the SA coefficients for each observable species, organized
per condition. It can be loaded in Python:

```python
from t3.common import sa_dict_from_yaml
from arc.common import read_yaml_file

raw = read_yaml_file('iteration_0/SA/sa_coefficients.yml')
sa_dict = sa_dict_from_yaml(raw)

# Access kinetics SA for the first condition
kinetics = sa_dict['kinetics'][0]
for observable, reactions in kinetics.items():
    print(f'{observable}: {len(reactions)} sensitive reactions')
```

### T3 libraries

T3 creates RMG-compatible thermo and kinetics libraries from the QM calculations.
These are stored in the RMG database path under the name specified by `library_name`
(default: `T3lib`).


## Restarting a project

T3 automatically detects existing iteration folders and resumes from where it left off.
Simply re-run T3 in the same project directory:

```bash
python ~/Code/T3/T3.py input.yml
```

T3 will scan for `iteration_*` folders and continue from the latest completed iteration.
