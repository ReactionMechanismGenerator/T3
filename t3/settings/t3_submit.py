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
# #SBATCH --mem-per-cpu={memory / cpus}
#
#
# export PYTHONPATH=$PYTHONPATH:~/Code/RMG-Py/
#
# conda activate rmg_env
#
# touch initial_time
#
# python-jl ~/Code/RMG-Py/rmg.py -n {cpus} input.py{max_iterations}
#
# touch final_time
#
# """,
#     'rmg': """Universe      = vanilla
#
# +JobName      = "{name}"
#
# log           = job.log
# output        = out.txt
# error         = err.txt
#
# getenv        = True
#
# should_transfer_files = no
#
# executable = job.sh
#
# request_cpus  = {cpus}
# request_memory = {memory}MB
#
# queue
#
# """,
#     'rmg_job': """#!/bin/bash -l
#
# touch initial_time
#
# source /srv01/technion/$USER/.bashrc
#
# conda activate rmg_env
#
# python-jl /Local/ce_dana/Code/RMG-Py/rmg.py -n {cpus} input.py{max_iterations}
#
# touch final_time
# 
# """,
    'rmg': """#!/bin/bash -l

#PBS -N {name}
#PBS -q zeus_long_q
#PBS -l walltime=168:00:00
#PBS -l select=1:ncpus={cpus}
#PBS -o out.txt
#PBS -e err.txt

PBS_O_WORKDIR={workdir}
cd $PBS_O_WORKDIR

conda activate rmg_env

touch initial_time

python-jl $rmgpy_path/rmg.py -n {cpus} input.py{max_iterations}

touch final_time

""",
}
