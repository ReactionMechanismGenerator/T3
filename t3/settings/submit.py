"""
Submit scripts
"""

# Submission scripts stored as a dictionary with software as the primary key.
submit_scripts = {
#     'rmg': """#!/bin/bash -l
# #SBATCH -J {name}
# #SBATCH -t 05-00:00:00
# #SBATCH -o out.txt
# #SBATCH -e err.txt
# #SBATCH --ntasks={cpus}
# #SBATCH --mem-per-cpu=9500
#
#
# export PYTHONPATH=$PYTHONPATH:~/Code/RMG-Py/
#
# conda activate rmg_env
#
# touch initial_time
#
# python-jl ~/Code/RMG-Py/rmg.py -n {cpus} input.py
#
# touch final_time
#
# """,
    'rmg': """Universe      = vanilla

+JobName      = "{name}"

log           = job.log
output        = out.txt
error         = err.txt

getenv        = True

should_transfer_files = no

executable = job.sh

request_cpus  = {cpus}
request_memory = {memory}MB

queue

""",
    'rmg_job': """#!/bin/bash -l

touch initial_time

source /srv01/technion/alongd/.bashrc

conda activate rmg_env

python-jl /Local/ce_dana/Code/RMG-Py/rmg.py -n 10 input.py -v 20

touch final_time

""",
}
