#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_rmg_runner module
"""

import os
import time

import pytest

from t3.common import TEST_DATA_BASE_PATH, EXAMPLES_BASE_PATH
from t3.runners.rmg_runner import rmg_job_converged, run_arkane_job, write_submit_script


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

python $rmgpy_path/rmg.py -n 16 input.pyNone

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

python $rmgpy_path/rmg.py -n 16 input.pyNone

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

python $rmgpy_path/rmg.py -n 8 input.py-m 100

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


class TestRunArkaneJob(object):
    """
    Test the run_arkane_job() success gate: it must require that
    ``sensitivity/sa_coefficients.yml`` both exists and was (re)written by the current run,
    not merely that some ``.yml`` file is sitting in ``sensitivity/`` (which could be stale
    output left over from a previous run).
    """

    @staticmethod
    def _write_input_file(tmp_path):
        input_file = tmp_path / 'input.py'
        input_file.write_text("# minimal Arkane input\n")
        return str(input_file)

    def test_fresh_sa_coefficients_file_returns_true(self, tmp_path, monkeypatch):
        """A sensitivity/sa_coefficients.yml written during this run -> True."""
        input_file = self._write_input_file(tmp_path)
        output_directory = str(tmp_path / 'output')

        def fake_run_arkane(statmech_dir):
            sensitivity_dir = os.path.join(statmech_dir, 'sensitivity')
            os.makedirs(sensitivity_dir, exist_ok=True)
            with open(os.path.join(sensitivity_dir, 'sa_coefficients.yml'), 'w') as f:
                f.write('some: data\n')

        monkeypatch.setattr('arc.statmech.arkane.run_arkane', fake_run_arkane)
        assert run_arkane_job(input_file=input_file, output_directory=output_directory) is True

    def test_no_sa_coefficients_file_returns_false(self, tmp_path, monkeypatch):
        """Arkane runs without error but never writes sa_coefficients.yml -> False."""
        input_file = self._write_input_file(tmp_path)
        output_directory = str(tmp_path / 'output')

        def fake_run_arkane(statmech_dir):
            pass  # writes nothing

        monkeypatch.setattr('arc.statmech.arkane.run_arkane', fake_run_arkane)
        assert run_arkane_job(input_file=input_file, output_directory=output_directory) is False

    def test_stale_pre_existing_sa_coefficients_file_returns_false(self, tmp_path, monkeypatch):
        """A pre-existing sa_coefficients.yml from a previous run, not rewritten -> deleted, False.

        run_arkane_job() deletes any pre-existing sa_coefficients.yml before invoking Arkane (rather
        than comparing mtimes, which is racy), so a stale file must be gone afterwards, and the job
        must report failure since Arkane never recreated it.
        """
        input_file = self._write_input_file(tmp_path)
        output_directory = str(tmp_path / 'output')
        sensitivity_dir = os.path.join(output_directory, 'sensitivity')
        os.makedirs(sensitivity_dir, exist_ok=True)
        stale_path = os.path.join(sensitivity_dir, 'sa_coefficients.yml')
        with open(stale_path, 'w') as f:
            f.write('stale: data\n')
        stale_time = time.time() - 3600
        os.utime(stale_path, (stale_time, stale_time))

        def fake_run_arkane(statmech_dir):
            pass  # does not rewrite sa_coefficients.yml

        monkeypatch.setattr('arc.statmech.arkane.run_arkane', fake_run_arkane)
        assert run_arkane_job(input_file=input_file, output_directory=output_directory) is False
        assert not os.path.isfile(stale_path), \
            'The stale sa_coefficients.yml should have been deleted before Arkane was invoked.'

    def test_traversal_required_artifact_raises_value_error(self, tmp_path, monkeypatch):
        """A required_artifact that escapes output_directory via '..' raises ValueError.

        Regression this guards: run_arkane_job() deletes the joined artifact path before running
        Arkane, so a traversal value like '../important.yml' used to DELETE a file outside the
        job directory. The path must be confined to output_directory, and the file outside must
        survive untouched.
        """
        input_file = self._write_input_file(tmp_path)
        output_directory = str(tmp_path / 'output')
        outside_file = tmp_path / 'important.yml'
        outside_file.write_text('precious: data\n')

        def fake_run_arkane(statmech_dir):
            pass  # must never be reached

        monkeypatch.setattr('arc.statmech.arkane.run_arkane', fake_run_arkane)
        traversal_artifact = os.path.join('..', 'important.yml')
        with pytest.raises(ValueError) as excinfo:
            run_arkane_job(input_file=input_file,
                           output_directory=output_directory,
                           required_artifact=traversal_artifact)
        assert 'important.yml' in str(excinfo.value)
        assert outside_file.is_file(), \
            'The file outside output_directory must not be deleted by a traversal required_artifact.'

    def test_absolute_required_artifact_raises_value_error(self, tmp_path, monkeypatch):
        """An absolute required_artifact is rejected outright with a ValueError naming it.

        Regression this guards: os.path.join(output_directory, <absolute path>) silently discards
        output_directory entirely, so an absolute value pointed the delete-then-require gate at an
        arbitrary filesystem location.
        """
        input_file = self._write_input_file(tmp_path)
        output_directory = str(tmp_path / 'output')
        absolute_target = tmp_path / 'elsewhere.yml'
        absolute_target.write_text('precious: data\n')

        def fake_run_arkane(statmech_dir):
            pass  # must never be reached

        monkeypatch.setattr('arc.statmech.arkane.run_arkane', fake_run_arkane)
        with pytest.raises(ValueError) as excinfo:
            run_arkane_job(input_file=input_file,
                           output_directory=output_directory,
                           required_artifact=str(absolute_target))
        assert 'elsewhere.yml' in str(excinfo.value)
        assert absolute_target.is_file(), \
            'The absolute-path target must not be deleted by run_arkane_job().'

    def test_relative_required_artifact_delete_then_require_preserved(self, tmp_path, monkeypatch):
        """A normal relative required_artifact keeps the delete-then-require semantics.

        A stale copy of the artifact is deleted before the run; when Arkane rewrites it the job
        succeeds, and when Arkane does not, the job fails and the stale copy is gone.
        """
        input_file = self._write_input_file(tmp_path)
        output_directory = str(tmp_path / 'output')
        os.makedirs(output_directory, exist_ok=True)
        artifact_path = os.path.join(output_directory, 'output.py')
        with open(artifact_path, 'w') as f:
            f.write('# stale artifact from a previous run\n')

        def fake_run_arkane_rewriting(statmech_dir):
            with open(os.path.join(statmech_dir, 'output.py'), 'w') as f:
                f.write('# fresh artifact\n')

        monkeypatch.setattr('arc.statmech.arkane.run_arkane', fake_run_arkane_rewriting)
        assert run_arkane_job(input_file=input_file,
                              output_directory=output_directory,
                              required_artifact='output.py') is True

        with open(artifact_path, 'w') as f:
            f.write('# stale artifact from a previous run\n')

        def fake_run_arkane_writing_nothing(statmech_dir):
            pass

        monkeypatch.setattr('arc.statmech.arkane.run_arkane', fake_run_arkane_writing_nothing)
        assert run_arkane_job(input_file=input_file,
                              output_directory=output_directory,
                              required_artifact='output.py') is False
        assert not os.path.isfile(artifact_path), \
            'The stale artifact should have been deleted before Arkane was invoked.'

    def test_differently_named_yml_present_returns_false(self, tmp_path, monkeypatch):
        """Some other .yml file in sensitivity/, but no sa_coefficients.yml -> False."""
        input_file = self._write_input_file(tmp_path)
        output_directory = str(tmp_path / 'output')

        def fake_run_arkane(statmech_dir):
            sensitivity_dir = os.path.join(statmech_dir, 'sensitivity')
            os.makedirs(sensitivity_dir, exist_ok=True)
            with open(os.path.join(sensitivity_dir, 'output.yml'), 'w') as f:
                f.write('some: other data\n')

        monkeypatch.setattr('arc.statmech.arkane.run_arkane', fake_run_arkane)
        assert run_arkane_job(input_file=input_file, output_directory=output_directory) is False


def teardown_module():
    """teardown any state that was previously setup with a setup_module method."""
    file_paths = [os.path.join(EXAMPLES_BASE_PATH, 'minimal', 'submit.sh')]
    for file_path in file_paths:
        if os.path.isfile(file_path):
            os.remove(file_path)
