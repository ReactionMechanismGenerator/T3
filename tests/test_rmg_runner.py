#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_rmg_runner module
"""

import os
from t3.common import DATA_BASE_PATH, EXAMPLES_BASE_PATH
from t3.runners.rmg_runner import write_submit_script



class TestWriteSubmitScript(object):

    #Need to think of multiple cases...

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