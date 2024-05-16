"""
Submit scripts
"""

# Submission scripts stored as a dictionary with software as the primary key.
#submit_scripts = {
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

submit_scripts = {
    'rmg': {
        'Slurm': """#!/bin/bash -l
#SBATCH -p hpc
#SBATCH -J {name}
#SBATCH -N 1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={memory}
#SBATCH -o out.txt
#SBATCH -e err.txt

# Source 
source /home/azureuser/.bashrc
source /etc/profile.d/00-aliases.sh

# Export
export PATH=$PATH:/home/azureuser/Code/RMG-Py
export PYTHONPATH=$PYTHONPATH:/home/azureuser/Code/RMG-Py

# It appears this is needed to keep the julia files writable...
sudo chmod 0777 -R /opt/mambaforge

conda activate rmg_env

touch initial_time

python-jl /home/azureuser/Code/RMG-Py/rmg.py {max_iterations} -n {cpus} input.py 

touch final_time

""",
        'PBS': """#!/bin/bash -l
#PBS -N {name}
#PBS -q zeus_combined_q
#PBS -l select=1:ncpus={cpus}
#PBS -o out.txt
#PBS -e err.txt

PBS_O_WORKDIR={workdir}
cd $PBS_O_WORKDIR

conda activate rmg_env

touch initial_time

python-jl ~/Code/RMG-Py/rmg.py -n {cpus} input.py {max_iterations}

touch final_time

""",
    }
}