#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_rmg_runner module
"""

import os
from t3.common import TEST_DATA_BASE_PATH, EXAMPLES_BASE_PATH
from t3.runners.rmg_runner import rmg_job_converged, write_submit_script


class TestWriteSubmitScript(object):

    def test_minimal_write_submit_script(self):
        """
        Test the write_submit_script() function with minimal input.
        write_submit_script params are set as their default values
        This test will create a job.sh file in the project directory path and the assertion will check if the file exists
        and matches the expected file
        """
        project_directory_path = os.path.join(EXAMPLES_BASE_PATH, "minimal")

        write_submit_script(project_directory_path,
                            cpus=None,
                            memory=None,
                            verbose=None,
                            max_iterations=None,
                            t3_project_name=None)

        expected = f"""#!/bin/bash -l

#PBS -N None_RMG
#PBS -q zeus_long_q
#PBS -l walltime=168:00:00
#PBS -l select=1:ncpus=16
#PBS -o out.txt
#PBS -e err.txt

PBS_O_WORKDIR={project_directory_path}
cd $PBS_O_WORKDIR

conda activate rmg_env

touch initial_time

python-jl $rmgpy_path/rmg.py -n 16 input.pyNone

touch final_time

"""

        assert os.path.isfile(os.path.join(project_directory_path, "submit.sh"))
        with open(os.path.join(project_directory_path, "submit.sh"), "r") as f:
            content = f.read()
        assert content == expected

    def test_minimal_project_name_included(self):
        """Test that thew minimal project name is included in the PBS submit script."""
        project_directory_path = os.path.join(EXAMPLES_BASE_PATH, "minimal")
        t3_proj_name = "T3_test_name"
        write_submit_script(project_directory_path,
                            cpus=None,
                            memory=None,
                            verbose=None,
                            max_iterations=None,
                            t3_project_name=t3_proj_name)
        expected_submit = f"""#!/bin/bash -l

#PBS -N {t3_proj_name}_RMG
#PBS -q zeus_long_q
#PBS -l walltime=168:00:00
#PBS -l select=1:ncpus=16
#PBS -o out.txt
#PBS -e err.txt

PBS_O_WORKDIR={project_directory_path}
cd $PBS_O_WORKDIR

conda activate rmg_env

touch initial_time

python-jl $rmgpy_path/rmg.py -n 16 input.pyNone

touch final_time

"""

        assert os.path.isfile(os.path.join(project_directory_path, "submit.sh"))
        with open(os.path.join(project_directory_path, "submit.sh"), "r") as f:
            content_submit = f.read()
        assert content_submit == expected_submit

    def test_minimal_parameters_set(self):
        """Test creating a submit script for the minimal example."""
        project_directory_path = os.path.join(EXAMPLES_BASE_PATH, "minimal")
        
        # To be edited by user if required:
        t3_proj_name = "T3_test_name"
        cpus = 8
        max_iter = "-m 100"
        mem = 16000  # in MB
        #########################################
        write_submit_script(project_directory_path,
                            cpus=cpus,
                            memory=mem,  # in MB
                            verbose="-v 20",
                            max_iterations=max_iter,
                            t3_project_name=t3_proj_name)
        
        expected_submit = f"""#!/bin/bash -l

#PBS -N T3_test_name_RMG
#PBS -q zeus_long_q
#PBS -l walltime=168:00:00
#PBS -l select=1:ncpus=8
#PBS -o out.txt
#PBS -e err.txt

PBS_O_WORKDIR={project_directory_path}
cd $PBS_O_WORKDIR

conda activate rmg_env

touch initial_time

python-jl $rmgpy_path/rmg.py -n 8 input.py-m 100

touch final_time

"""
        assert os.path.isfile(os.path.join(project_directory_path, "submit.sh"))

        with open(os.path.join(project_directory_path, "submit.sh"), "r") as submit_file:
            content_submit = submit_file.read()
        assert content_submit == expected_submit

    def test_rmg_job_converged(self):
        """Test correctly identifying whether an RMG job converged ot not, and if not which error was received."""
        rmg_folder_1 = os.path.join(TEST_DATA_BASE_PATH, 'rmg_convergence', '1_frag_error')
        converged, error = rmg_job_converged(project_directory=rmg_folder_1)
        assert not converged
        assert error == "AttributeError: 'Fragment' object has no attribute 'count_internal_rotors'"

        rmg_folder_2 = os.path.join(TEST_DATA_BASE_PATH, 'rmg_convergence', '2_converged')
        converged, error = rmg_job_converged(project_directory=rmg_folder_2)
        assert converged
        assert error is None


def teardown_module():
    """teardown any state that was previously setup with a setup_module method."""
    file_paths = [os.path.join(EXAMPLES_BASE_PATH, 'minimal', 'submit.sh')]
    for file_path in file_paths:
        if os.path.isfile(file_path):
            os.remove(file_path)
