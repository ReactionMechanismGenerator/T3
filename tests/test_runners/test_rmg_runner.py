#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_rmg_runner module
"""

import os
import shutil
import signal
import subprocess
import time

import pytest

from t3.common import TEST_DATA_BASE_PATH, EXAMPLES_BASE_PATH
from t3.runners import rmg_runner
from t3.runners.rmg_runner import (_kill_process_group,
                                   rmg_job_converged,
                                   run_arkane_job,
                                   run_rmg_incore,
                                   write_submit_script)


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
        result = run_arkane_job(input_file=input_file, output_directory=output_directory)
        assert bool(result) is True
        assert result.timed_out is False

    def test_no_sa_coefficients_file_returns_false(self, tmp_path, monkeypatch):
        """Arkane runs without error but never writes sa_coefficients.yml -> False."""
        input_file = self._write_input_file(tmp_path)
        output_directory = str(tmp_path / 'output')

        def fake_run_arkane(statmech_dir):
            pass  # writes nothing

        monkeypatch.setattr('arc.statmech.arkane.run_arkane', fake_run_arkane)
        result = run_arkane_job(input_file=input_file, output_directory=output_directory)
        assert bool(result) is False
        assert result.timed_out is False
        assert 'did not produce' in result.reason

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
        assert bool(run_arkane_job(input_file=input_file, output_directory=output_directory)) is False
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
        assert bool(run_arkane_job(input_file=input_file,
                                   output_directory=output_directory,
                                   required_artifact='output.py')) is True

        with open(artifact_path, 'w') as f:
            f.write('# stale artifact from a previous run\n')

        def fake_run_arkane_writing_nothing(statmech_dir):
            pass

        monkeypatch.setattr('arc.statmech.arkane.run_arkane', fake_run_arkane_writing_nothing)
        assert bool(run_arkane_job(input_file=input_file,
                                   output_directory=output_directory,
                                   required_artifact='output.py')) is False
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
        assert bool(run_arkane_job(input_file=input_file, output_directory=output_directory)) is False


class TestRunArkaneJobTimeout(object):
    """
    Test run_arkane_job()'s ``timeout`` path: with a deadline, Arkane runs in a separate process
    session that can actually be KILLED (in-process ARC calls cannot be; the subprocess ARC spawns
    exposes no timeout and no PID -- see the design note in run_arkane_job's docstring), and a
    deadline overrun surfaces as a falsy, ``timed_out`` result with a reason -- never as an
    exception, which would destroy the run record the caller is building.
    """

    @staticmethod
    def _write_input_file(tmp_path):
        input_file = tmp_path / 'input.py'
        input_file.write_text("# minimal Arkane input\n")
        return str(input_file)

    def test_overrunning_the_deadline_reports_timed_out_and_kills_the_process(self, tmp_path, monkeypatch):
        """A run overrunning ``timeout`` is killed (whole process group), reaped, and reported as a
        falsy result with ``timed_out`` set and the deadline named in the reason."""
        input_file = self._write_input_file(tmp_path)
        output_directory = str(tmp_path / 'output')
        spawned = {}

        def fake_spawn(output_directory):
            proc = subprocess.Popen(['bash', '-c', 'sleep 60'], start_new_session=True)
            spawned['proc'] = proc
            return proc

        monkeypatch.setattr('t3.runners.rmg_runner._spawn_arkane_subprocess', fake_spawn)
        result = run_arkane_job(input_file=input_file, output_directory=output_directory, timeout=0.2)

        assert bool(result) is False
        assert result.timed_out is True
        assert '0.2' in result.reason and 'timed out' in result.reason
        assert spawned['proc'].poll() is not None, \
            'The overrunning process must be dead and reaped, not abandoned to keep writing into ' \
            'a run directory already declared failed.'

    def test_finishing_within_the_deadline_succeeds(self, tmp_path, monkeypatch):
        """A run that writes the required artifact and exits before ``timeout`` reports success --
        the timeout guards the run, it must not fail a run that met it."""
        input_file = self._write_input_file(tmp_path)
        output_directory = str(tmp_path / 'output')

        def fake_spawn(output_directory):
            artifact_dir = os.path.join(output_directory, 'sensitivity')
            return subprocess.Popen(
                ['bash', '-c', f'mkdir -p "{artifact_dir}" && echo "some: data" > '
                               f'"{os.path.join(artifact_dir, "sa_coefficients.yml")}"'],
                start_new_session=True)

        monkeypatch.setattr('t3.runners.rmg_runner._spawn_arkane_subprocess', fake_spawn)
        result = run_arkane_job(input_file=input_file, output_directory=output_directory, timeout=30)

        assert bool(result) is True
        assert result.timed_out is False

    def test_nonzero_exit_within_the_deadline_reports_failure(self, tmp_path, monkeypatch):
        """A subprocess that dies (nonzero exit) before the deadline is a plain failure, mirroring
        the in-process path's exception handling -- and is NOT labelled a timeout."""
        input_file = self._write_input_file(tmp_path)
        output_directory = str(tmp_path / 'output')

        def fake_spawn(output_directory):
            return subprocess.Popen(['bash', '-c', 'exit 3'], start_new_session=True)

        monkeypatch.setattr('t3.runners.rmg_runner._spawn_arkane_subprocess', fake_spawn)
        result = run_arkane_job(input_file=input_file, output_directory=output_directory, timeout=30)

        assert bool(result) is False
        assert result.timed_out is False
        assert '3' in result.reason

    @pytest.mark.parametrize('bad_timeout', [0, -5, float('inf'), float('nan'), True, '60'])
    def test_invalid_timeout_is_refused(self, tmp_path, bad_timeout):
        """A timeout that is not a positive finite number is an argument error, refused before any
        side effect (no artifact deletion, no spawn)."""
        input_file = self._write_input_file(tmp_path)
        output_directory = str(tmp_path / 'output')
        with pytest.raises(ValueError, match='timeout'):
            run_arkane_job(input_file=input_file, output_directory=output_directory, timeout=bad_timeout)


def teardown_module():
    """teardown any state that was previously setup with a setup_module method."""
    file_paths = [os.path.join(EXAMPLES_BASE_PATH, 'minimal', 'submit.sh')]
    for file_path in file_paths:
        if os.path.isfile(file_path):
            os.remove(file_path)


def _pid_is_alive(pid: int) -> bool:
    """Whether a process with the given PID currently exists.

    Args:
        pid (int): The PID to check.

    Returns:
        bool: Whether the process exists.
    """
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        # The process exists but is not signalable by this user. That is alive, not dead --
        # reporting it dead here would let _wait_for_death() succeed against a running process,
        # which is the one way these tests could pass while the bug they guard is present.
        pass
    return True


def _wait_for_death(pid: int, timeout: float = 20.0) -> bool:
    """Wait for a PID to disappear.

    Args:
        pid (int): The PID to wait for.
        timeout (float, optional): Max seconds to wait.

    Returns:
        bool: Whether the process is gone.
    """
    deadline = time.time() + timeout
    while time.time() < deadline:
        if not _pid_is_alive(pid):
            return True
        time.sleep(0.1)
    return not _pid_is_alive(pid)


def _wait_for_pid_file(pid_file: str, timeout: float = 20.0) -> int:
    """Wait for a spawned grandchild to write its PID to a file, and return that PID.

    Args:
        pid_file (str): The path of the file the grandchild writes its PID into.
        timeout (float, optional): Max seconds to wait.

    Returns:
        int: The grandchild PID.
    """
    deadline = time.time() + timeout
    while time.time() < deadline:
        if os.path.isfile(pid_file):
            with open(pid_file) as f:
                content = f.read().strip()
            if content.isdigit():
                return int(content)
        time.sleep(0.05)
    raise AssertionError(f'The grandchild process never wrote its PID to {pid_file}')


# A shell one-liner that backgrounds a long-lived grandchild, records its PID, and then blocks.
# This mimics the RMG incore command shape (bash -> micromamba run -> python -> tee), where the
# process that must be killed is not the direct child of the runner.
_GRANDCHILD_SCRIPT = 'sleep 600 & echo $! > "{pid_file}"; wait'


@pytest.mark.skipif(os.name != 'posix' or shutil.which('bash') is None,
                    reason='needs POSIX process groups, signals and bash')
class TestKillProcessGroup(object):
    """Test that a timed-out RMG run is killed rather than orphaned."""

    def test_kill_process_group_kills_the_grandchild(self, tmp_path):
        """The whole process group dies, not just the direct child.

        Killing only the direct child (the old ``subprocess.run`` behavior) leaves the real RMG
        process running, which is what showed up in CI logs as
        ``Terminate orphan process: pid (20994) (python3.14)``.
        """
        pid_file = str(tmp_path / 'grandchild.pid')
        process = subprocess.Popen(_GRANDCHILD_SCRIPT.format(pid_file=pid_file),
                                   shell=True, executable='/bin/bash', start_new_session=True)
        try:
            grandchild_pid = _wait_for_pid_file(pid_file)
            assert _pid_is_alive(grandchild_pid)
            _kill_process_group(process)
            assert _wait_for_death(process.pid), 'the direct child survived'
            assert _wait_for_death(grandchild_pid), 'the grandchild was orphaned instead of killed'
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()

    def test_process_group_is_alive_tracks_the_group_not_the_child(self):
        """The SIGKILL guard reports a group alive while any member remains, and dead after.

        ``_kill_process_group()`` escalates to SIGKILL only when this returns True, so that it
        never signals a PGID that has already been recycled. It deliberately probes the whole
        group rather than the direct child: the child is a shell wrapper and can exit while the
        process that matters is still running, which is the orphan this module exists to prevent.
        """
        process = subprocess.Popen('sleep 600', shell=True, executable='/bin/bash',
                                   start_new_session=True)
        pgid = os.getpgid(process.pid)
        try:
            assert rmg_runner._process_group_is_alive(pgid)
        finally:
            os.killpg(pgid, signal.SIGKILL)
            process.wait()
        assert not rmg_runner._process_group_is_alive(pgid), 'an empty group reported alive'

    def test_run_rmg_incore_timeout_kills_the_grandchild(self, tmp_path, monkeypatch):
        """``run_rmg_incore()``'s timeout path kills the whole process tree it spawned.

        The command is replaced with a stand-in that backgrounds a grandchild, but every
        ``Popen`` keyword argument the production code passes is preserved, so this exercises the
        real timeout handling rather than a re-implementation of it.
        """
        pid_file = str(tmp_path / 'grandchild.pid')
        recorded_kwargs, spawned = {}, []
        real_popen = subprocess.Popen

        def fake_popen(cmd, **kwargs):
            recorded_kwargs.update(kwargs)
            process = real_popen(_GRANDCHILD_SCRIPT.format(pid_file=pid_file), **kwargs)
            spawned.append(process)
            return process

        monkeypatch.setattr(subprocess, 'Popen', fake_popen)
        # '00:00:00:02' is a 2 second walltime, so the timeout fires while the grandchild lives.
        exception_raised = run_rmg_incore(rmg_input_file_path=str(tmp_path / 'input.py'),
                                          walltime='00:00:00:02')
        monkeypatch.undo()

        assert exception_raised is True
        assert recorded_kwargs.get('start_new_session') is True, \
            'the child must lead its own process group for os.killpg() to reach the whole tree'
        grandchild_pid = _wait_for_pid_file(pid_file)
        assert _wait_for_death(spawned[0].pid), 'the direct child survived the timeout'
        assert _wait_for_death(grandchild_pid), 'the grandchild was orphaned instead of killed'
