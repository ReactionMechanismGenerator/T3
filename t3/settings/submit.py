"""
Submit scripts
"""

# Submission scripts stored as a dictionary with software as the primary key.
submit_scripts = {
    'rmg': """#!/bin/bash -l
#SBATCH -J {name}
#SBATCH -t 05-00:00:00
#SBATCH -o out.txt
#SBATCH -e err.txt
#SBATCH --ntasks={cpus}
#SBATCH --mem-per-cpu=9500


export PYTHONPATH=$PYTHONPATH:~/Code/RMG-Py/

conda activate rmg_env

touch initial_time

python-jl ~/Code/RMG-Py/rmg.py -n {cpus} input.py

touch final_time

""",
}
