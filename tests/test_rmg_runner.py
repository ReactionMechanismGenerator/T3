#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_rmg_runner module
"""

import os
from t3.common import DATA_BASE_PATH, EXAMPLES_BASE_PATH
from t3.runners.rmg_runner import write_submit_script
import pytest
from tests.common import run_minimal
import shutil

class TestWriteSubmitScript(object):



    def test_minimial_write_submit_script(self):
        """Test the write_submit_script() function with minimal input.
            write_submit_script params are set as their default values
            This test will create a job.sh file in the project directory path and the assertion will check if the file exists
            and matches the expected file"""
        project_directory_path = os.path.join(EXAMPLES_BASE_PATH, "minimal")

        actual = write_submit_script(project_directory_path,
                                    cpus=None,
                                    memory=None,
                                    verbose=None,
                                    max_iterations=None,
                                    t3_project_name=None)
        

        expected = """#!/bin/bash -l

touch initial_time

source /srv01/technion/$USER/.bashrc

conda activate rmg_env

python-jl /Local/ce_dana/Code/RMG-Py/rmg.py -n 16 input.py

touch final_time

"""

        assert os.path.isfile(os.path.join(project_directory_path,"job.sh")) == True
        with open(os.path.join(project_directory_path,"job.sh"),"r") as bash_file:
            content = bash_file.read()
        assert content == expected

        os.remove(os.path.join(project_directory_path,"job.sh"))

    def test_minimal_project_name_included(self):
        """
        This test has been constructed to ensure that when the user provides a t3 project name, it will be displayed in the submit name whilst also
        not affecting the job file.
        """
        
        project_directory_path = os.path.join(EXAMPLES_BASE_PATH, "minimal")
        t3_proj_name = "T3_test_name"
        actual = write_submit_script(project_directory_path,
                                    cpus=None,
                                    memory=None,
                                    verbose=None,
                                    max_iterations=None,
                                    t3_project_name= t3_proj_name)
        expected_bash = """#!/bin/bash -l

touch initial_time

source /srv01/technion/$USER/.bashrc

conda activate rmg_env

python-jl /Local/ce_dana/Code/RMG-Py/rmg.py -n 16 input.py

touch final_time

"""
        expected_submit = """Universe      = vanilla

+JobName      = "{t3_project_name}"

log           = job.log
output        = out.txt
error         = err.txt

getenv        = True

should_transfer_files = no

executable = job.sh

request_cpus  = 16
request_memory = 25000MB

queue

""".format(t3_project_name = t3_proj_name + "_RMG")

        assert os.path.isfile(os.path.join(project_directory_path,"job.sh")) == True
        assert os.path.isfile(os.path.join(project_directory_path, "submit.sub")) == True


        with open(os.path.join(project_directory_path,"job.sh"),"r") as bash_file:
            content_bash = bash_file.read()
        assert content_bash == expected_bash
        with open(os.path.join(project_directory_path,"submit.sub"),"r") as submit_file:
            content_submit = submit_file.read()
        assert content_submit == expected_submit

        os.remove(os.path.join(project_directory_path,"job.sh"))
        os.remove(os.path.join(project_directory_path,"submit.sub"))


    def test_minimal_parameters_set(self):
        """
        This test has been built to ensure that chosen parameters by the user are being correctly written into the job and submit file
        """
        project_directory_path = os.path.join(EXAMPLES_BASE_PATH, "minimal")
        
        ####To be edited by user if required####
        t3_proj_name = "T3_test_name"
        cpus = 8
        max_iter = "-m 100"
        mem = 16000 #in MB
        #########################################
        actual = write_submit_script(project_directory_path,
                                    cpus=cpus,
                                    memory=mem, #in MB
                                    verbose="-v 20",
                                    max_iterations=max_iter,
                                    t3_project_name= t3_proj_name)

        expected_bash = """#!/bin/bash -l

touch initial_time

source /srv01/technion/$USER/.bashrc

conda activate rmg_env

python-jl /Local/ce_dana/Code/RMG-Py/rmg.py -n {cpus} input.py{max_iter}

touch final_time

""".format( cpus = cpus, max_iter = max_iter )
        
        expected_submit = """Universe      = vanilla

+JobName      = "{t3_project_name}"

log           = job.log
output        = out.txt
error         = err.txt

getenv        = True

should_transfer_files = no

executable = job.sh

request_cpus  = {cpus}
request_memory = {mem}MB

queue

""".format(t3_project_name = t3_proj_name+"_RMG", cpus = cpus, mem = mem )
        assert os.path.isfile(os.path.join(project_directory_path,"job.sh")) == True
        assert os.path.isfile(os.path.join(project_directory_path, "submit.sub")) == True
        with open(os.path.join(project_directory_path,"job.sh"),"r") as bash_file:
            content_bash = bash_file.read()
        assert content_bash == expected_bash

        with open(os.path.join(project_directory_path,"submit.sub"),"r") as submit_file:
            content_submit = submit_file.read()
        assert content_submit == expected_submit

        os.remove(os.path.join(project_directory_path,"job.sh"))
        os.remove(os.path.join(project_directory_path,"submit.sub"))
        
    def test_minimal_incorrect_mem(self):
        """
        This test has been constructed to raise an assertion error due to the user specifiying the incorrent memory.
        Memory selection for assertion error can either be less than 8000 or higher than 32000 (Subject to change) OR/AND
        the selection of Memory may not be divisable by 1000 (Subject to change)
        """
        project_directory_path = os.path.join(EXAMPLES_BASE_PATH, "minimal")
        
        ####To be edited by user if required####
        t3_proj_name = "T3_test_name"
        cpus = 8
        max_iter = "-m 100"
        mem = 850 #in MB - This is to be incorrect
        #########################################

        with pytest.raises(ValueError) as value_error:
                    write_submit_script(project_directory_path,
                                    cpus=cpus,
                                    memory=mem, #in MB
                                    verbose="-v 20",
                                    max_iterations=max_iter,
                                    t3_project_name= t3_proj_name)
        
        assert "Memory is in MB and it is strongly recommended to set the memory above 1000" in str(value_error.value)
        

        
    def test_run_rmg_into_write_submit_script(self):
        
        
        
        t3 = run_minimal(iteration=1, set_paths=True)
        t3.rmg['rmg_execution_type'] = 'local'
        
        rmg_base_path = os.path.join(t3.project_directory, 'iteration_1','RMG')
        
        t3.schema['t3']['options']['max_RMG_exceptions_allowed']= 0
        os.makedirs(rmg_base_path)
        with open(os.path.join(rmg_base_path,'RMG.log'), 'w') as rmg_log:
            rmg_log.write("""MODEL GENERATION COMPLETED
                          """)
        rmg_log.close()
        
        
        t3.run_rmg()

        
        expected_bash = """#!/bin/bash -l

touch initial_time

source /srv01/technion/$USER/.bashrc

conda activate rmg_env

python-jl /Local/ce_dana/Code/RMG-Py/rmg.py -n 10 input.py

touch final_time

"""
        expected_submit = """Universe      = vanilla

+JobName      = "T3_minimal_example_RMG"

log           = job.log
output        = out.txt
error         = err.txt

getenv        = True

should_transfer_files = no

executable = job.sh

request_cpus  = 10
request_memory = 25000MB

queue

"""

        assert os.path.isfile(os.path.join(rmg_base_path,"job.sh")) == True
        assert os.path.isfile(os.path.join(rmg_base_path, "submit.sub")) == True
        with open(os.path.join(rmg_base_path,"job.sh"),"r") as bash_file:
            content_bash = bash_file.read()
        assert content_bash == expected_bash

        with open(os.path.join(rmg_base_path,"submit.sub"),"r") as submit_file:
            content_submit = submit_file.read()
        assert content_submit == expected_submit
        
        shutil.rmtree(os.path.join(t3.project_directory))