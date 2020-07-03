# Output

After running a Project, the local Project folder will contain the following directory tree
(**bold** represents folders, *italics* represents files):

- *t3.log*: Details of all project execution procedures.
- **iteration_x**:
    - **RMG**:
        - *input.py*: The RMG input file for iteration x (automatically written by T3)
        - *RMG.log*: The RMG log file for iteration x.
        - Additional RMG files and folders
    - **ARC**:
        - *restart.yml*: The ARC restart file for iteration x.
        - *arc.log*: The ARC log file for iteration x.
        - **calcs**: All spawned jobs sorted by species name (including transition states).
        - **output**: Additional ARC output files and folders.
