"""
A "keep alive" runner tool for RMG on a server.
Should be executed locally on the head node using the t3 environment.
"""

import contextlib
import datetime
import logging
import math
import os
import shlex
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from typing import TYPE_CHECKING

from arc.job.local import (_determine_job_id,
                           change_mode,
                           execute_command,
                           parse_running_jobs_ids)

from t3.imports import local_t3_path, settings, submit_scripts
from t3.utils.fix_cantera import fix_cantera

if TYPE_CHECKING:
    from t3.logger import Logger


MEM = settings['rmg_initial_memory'] * 1000  # MB
SLEEP_TIME = 6  # hours
MAX_RMG_RUNS_PER_ITERATION = 5

RMG_EXECUTION_TYPE = settings['execution_type']['rmg']

# These describe the *local cluster software* (PBS/Slurm/OGE/HTCondor) and are therefore keyed on
# ``cluster_soft`` -- NOT on the RMG execution type. Gating them on ``RMG_EXECUTION_TYPE == 'local'``
# blanked all four to '' on any machine that runs RMG ``incore`` while still having a queue
# configured; the empty ``SUBMIT_FILENAME`` then made ``write_submit_script``/``submit_job`` join an
# empty component onto the project directory and open the directory itself (IsADirectoryError).
#
# These globals describe HOW to talk to a queue, so the only thing that can define them is which
# cluster software is configured -- not which execution type RMG happens to be set to. Deriving
# them from `cluster_soft` keeps them well-formed whenever a queue exists; whether a submit script
# is actually USED is a separate question, answered at each use site (`write_submit_script` /
# `submit_job`), which refuse when no cluster software is configured.
LOCAL_CLUSTER_SOFTWARE = settings['servers']['local'].get('cluster_soft') or ''
if LOCAL_CLUSTER_SOFTWARE:
    SUBMIT_COMMAND = settings['submit_command'][LOCAL_CLUSTER_SOFTWARE]
    CHECK_STATUS_COMMAND = settings['check_status_command'][LOCAL_CLUSTER_SOFTWARE]
    SUBMIT_FILENAME = settings['submit_filenames'][LOCAL_CLUSTER_SOFTWARE]
else:
    SUBMIT_COMMAND = CHECK_STATUS_COMMAND = SUBMIT_FILENAME = ''


def write_submit_script(project_directory: str,
                        cpus: int | None = None,
                        memory: int | None = None,
                        verbose: str | None = None,
                        max_iterations: str | None = None,
                        t3_project_name: str | None = None,
                        ) -> None:
    """
    Write an RMG submit script.

    Args:
        project_directory (str): The full path to the project directory.
        cpus (int, optional): The number of CPUs for an RMG parallelization.
        memory (int, optional): The memory in MB for an RMG run, defaults to ``MEM``.
        verbose (str, optional): Level of verbosity, e.g., ``-v 10``.
        max_iterations (str, optional): Max RMG iterations, e.g., ``-m 100``.
        t3_project_name (str, optional): The T3 project name, used for setting a job name on the server for the RMG run.

    Raises:
        ValueError: If no local cluster software is configured (``SUBMIT_FILENAME`` is empty). A
                    queue submit script is only meaningful for queue-based ('local') RMG execution;
                    refusing here names the real condition instead of joining an empty filename onto
                    ``project_directory`` and opening the directory itself.
    """
    global MEM
    if not SUBMIT_FILENAME:
        raise ValueError(
            "Cannot write an RMG submit script: no local cluster software is configured "
            "(settings['servers']['local']['cluster_soft'] is missing or empty). A submit script is only "
            "meaningful for queue-based ('local') RMG execution; on a queue-less machine set the "
            "RMG execution type to 'incore' instead."
        )
    cpus = cpus or 16
    submit_scripts_content = submit_scripts['rmg'].format(name=f'{t3_project_name}_RMG' or 'T3_RMG',
                                                          cpus=cpus,
                                                          memory=memory or MEM,
                                                          workdir=project_directory,
                                                          max_iterations=max_iterations,
                                                          )
    with open(os.path.join(project_directory, SUBMIT_FILENAME), 'w') as f:
        f.write(submit_scripts_content)
    if 'rmg_job' in submit_scripts.keys():
        # Write an aux submit script, e.g., as required for HTCondor.
        max_iterations = max_iterations or ''
        aux_submit_scripts_content = submit_scripts['rmg_job'].format(cpus=cpus,
                                                                      max_iterations=max_iterations,
                                                                      )
        with open(os.path.join(project_directory, 'job.sh'), 'w') as f:
            f.write(aux_submit_scripts_content)
            change_mode(mode='+x', file_name='job.sh', path=project_directory)


def submit_job(project_directory: str,
               logger: Logger,
               cluster_soft: str,
               memory: int | None = None,
               ) -> tuple[str | None, str | None]:
    """
    Submit an RMG job.

    Args:
        project_directory (str): The job (folder) name.
        logger (Logger): The T3 Logger object instance.
        cluster_soft (str): The server's cluster software.
        memory (int, optional): The memory in MB for an RMG run. Only used for reporting.

    Returns:
        Tuple[Optional[str], Optional[str]]: job_status, job_id
    """
    global MEM
    if not SUBMIT_COMMAND or not SUBMIT_FILENAME:
        raise ValueError(
            "Cannot submit an RMG job: no local cluster software is configured "
            "(settings['servers']['local']['cluster_soft'] is missing or empty). Without it the submit "
            "command would be 'cd dir;  ; cd ..'. Set the RMG execution type to 'incore' on a "
            "queue-less machine."
        )
    job_status = ''
    job_id = 0
    cmd = f"cd {project_directory}; {SUBMIT_COMMAND} {SUBMIT_FILENAME}; cd .."
    stdout, stderr = execute_command(cmd)
    if not len(stdout):
        time.sleep(10)
        stdout, stderr = execute_command(cmd)
    if not len(stdout):
        return None, None
    if len(stderr) > 0 or len(stdout) == 0:
        logger.info(f'Got the following error when trying to submit job {project_directory}:\n{stderr}.')
        job_status = 'errored'
    else:
        job_id = _determine_job_id(stdout=stdout, cluster_soft=cluster_soft.lower())
    m = memory or MEM or 0
    logger.info(f'\nSuccessfully submitted job {project_directory},\n'
                f'job ID = {job_id}, requested memory: {m / 1000:.2f} GB.')
    return job_status, job_id


def check_running_jobs_ids(cluster_soft: str) -> list[str]:
    """
    Check which jobs are still running on the server for this user.

    Args:
        cluster_soft (str): The server's cluster software.

    Returns:
        List(str): List of job IDs.
    """
    stdout = execute_command(CHECK_STATUS_COMMAND)[0]
    running_job_ids = parse_running_jobs_ids(stdout, cluster_soft=cluster_soft.lower())
    return running_job_ids


def rmg_job_converged(project_directory: str) -> tuple[bool, str | None]:
    """
    Determine whether an RMG job has converged.

    Args:
        project_directory (str): The job (folder) name.

    Returns:
        Tuple[bool, Optional[str]]:
            - bool: Whether this RMG run has converged.
            - Optional[str]: The error due to which this RMG run crashed.
    """
    rmg_converged, error = False, None
    rmg_log_path = os.path.join(project_directory, 'RMG.log')
    rmg_err_path = os.path.join(project_directory, 'err.txt')
    if os.path.isfile(rmg_log_path):
        with open(rmg_log_path) as f:
            lines = f.readlines()
            len_lines = len(lines)
            for i in range(10):
                if 'MODEL GENERATION COMPLETED' in lines[len_lines - 1 - i]:
                    rmg_converged = True
                    break
    if not rmg_converged and os.path.isfile(rmg_err_path):
        with open(rmg_err_path) as f:
            lines = f.readlines()
        for line in lines[::-1]:
            if 'Error' in line:
                error = line.strip()
                break
    return rmg_converged, error


_DEFAULT_RMG_TIMEOUT_S = 6 * 3600  # 6 hours
_KILL_GRACE_PERIOD_S = 10  # seconds to wait after SIGTERM before sending SIGKILL

logger = logging.getLogger(__name__)


def _kill_process_group(process: 'subprocess.Popen') -> None:
    """Kill a process and every descendant it spawned.

    The RMG incore command is a shell pipeline (``bash`` -> ``micromamba run`` -> ``python`` ->
    ``tee``), so killing only the direct child leaves the actual RMG process orphaned and still
    running. The child is started with ``start_new_session=True``, which makes it the leader of
    its own process group, so signalling that group reaches the whole tree.

    Args:
        process (subprocess.Popen): The process to kill, started with ``start_new_session=True``.
    """
    try:
        pgid = os.getpgid(process.pid)
    except ProcessLookupError:
        # The child is already gone. Reap it anyway: without a wait() it stays a zombie until
        # the Popen object is collected, and there is nothing left to signal.
        with contextlib.suppress(subprocess.TimeoutExpired):
            process.wait(timeout=_KILL_GRACE_PERIOD_S)
        return
    except PermissionError:
        return
    with contextlib.suppress(ProcessLookupError, PermissionError):
        os.killpg(pgid, signal.SIGTERM)
    with contextlib.suppress(subprocess.TimeoutExpired):
        # The grace period expiring is the normal path to SIGKILL below, not an error.
        process.wait(timeout=_KILL_GRACE_PERIOD_S)
    if _process_group_is_alive(pgid):
        with contextlib.suppress(ProcessLookupError, PermissionError):
            os.killpg(pgid, signal.SIGKILL)


def _process_group_is_alive(pgid: int) -> bool:
    """Whether any process remains in the given process group.

    Signal 0 performs the permission and existence checks without delivering anything. This is
    checked against the whole *group* rather than against the direct child, because the child is
    a shell wrapper: it can exit while the RMG process it spawned is still running, which is the
    orphan this module exists to prevent. Escalating to ``SIGKILL`` unconditionally would instead
    risk signalling an unrelated group, should the PGID have been recycled in the meantime.

    Args:
        pgid (int): The process group id to probe.

    Returns:
        bool: Whether the group still has at least one member.
    """
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        # The group exists but is not ours to signal -- alive, and SIGKILL will fail harmlessly.
        return True
    return True


def _parse_walltime_to_seconds(walltime: str) -> int:
    """Parse a 'DD:HH:MM:SS' walltime string to total seconds. Returns 0 for '00:00:00:00'."""
    parts = walltime.split(':')
    if len(parts) != 4:
        return 0
    try:
        days, hours, minutes, seconds = (int(p) for p in parts)
    except ValueError:
        return 0
    return days * 86400 + hours * 3600 + minutes * 60 + seconds


def run_rmg_incore(rmg_input_file_path: str,
                   verbose: int | None = None,
                   max_iterations: int | None = None,
                   walltime: str | None = None,
                   ) -> bool:
    """
    Run RMG incore under the rmg_env.

    Args:
        rmg_input_file_path (str): The path to the RMG input file.
        max_iterations (int, optional): Max RMG iterations.
        verbose (int, optional): Level of verbosity.
        walltime (str, optional): Max walltime in 'DD:HH:MM:SS' format. Defaults to 6 hours.

    Returns:
        bool: Whether an exception was raised.
    """
    timeout_s = _parse_walltime_to_seconds(walltime) if walltime else 0
    if timeout_s <= 0:
        timeout_s = _DEFAULT_RMG_TIMEOUT_S
    project_directory = os.path.abspath(os.path.dirname(rmg_input_file_path))
    verbose = f' -v {verbose}' if verbose is not None else ''
    max_iterations = f' -m {max_iterations}' if max_iterations is not None else ''
    script_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'rmg_incore_script.py')
    inner_cmd = (f'python {script_path} {rmg_input_file_path}{verbose}{max_iterations} '
                 f'> >(tee -a out.txt) 2> >(tee -a err.txt >&2)')
    shell_script = rf'''bash -lc 'set -uo pipefail
cd "{project_directory}"
if command -v micromamba >/dev/null 2>&1; then
    micromamba run -n rmg_env bash -c "{inner_cmd}"
elif command -v conda >/dev/null 2>&1 || command -v mamba >/dev/null 2>&1; then
    conda run -n rmg_env bash -c "{inner_cmd}"
else
    echo "Micromamba/Mamba/Conda required" >&2
    exit 1
fi' '''
    process = subprocess.Popen(shell_script, shell=True, executable='/bin/bash',
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                               start_new_session=True)
    try:
        _, stderr_text = process.communicate(timeout=timeout_s)
        stderr_text = stderr_text or ''
    except subprocess.TimeoutExpired:
        _kill_process_group(process)
        with contextlib.suppress(Exception):
            process.communicate(timeout=_KILL_GRACE_PERIOD_S)
        logger.error(f'RMG incore timed out after {timeout_s}s')
        return True
    if process.returncode != 0:
        logger.error(f'RMG incore exited with code {process.returncode}')
        return True
    if 'RMG threw an exception and did not converge.' in stderr_text:
        return True
    return False


def run_rmg_in_local_queue(project_directory: str,
                           logger: Logger,
                           memory: int | None = None,
                           cpus: int | None = None,
                           max_iterations: int | None = None,
                           restart_rmg: bool = False,
                           verbose: int | None = None,
                           t3_project_name: str | None = None,
                           ):
    """
    Run RMG on the queue of the local server (under the rmg_env).

    Args:
        project_directory (str): The path to the RMG folder.
        logger (Logger): The T3 Logger object instance.
        memory (int, optional): The submit script memory in MB.
        cpus (int, optional): The number of CPUs for an RMG parallelization.
        max_iterations (int, optional): Max RMG iterations.
        restart_rmg (bool, optional): Whether this RMG run should trigger a seed restart.
        verbose (int, optional): Level of verbosity.
        t3_project_name (str, optional): The T3 project name, used for setting a job name on the server for the RMG run.

    Returns:
        Optional[str]: The job ID.
    """
    verbose = f' -v {verbose}' if verbose is not None else ''
    max_iterations = f' -m {max_iterations}' if max_iterations is not None else ''
    write_submit_script(project_directory=project_directory,
                        cpus=cpus or settings['servers']['local']['cpus'],
                        memory=memory,
                        verbose=verbose,
                        max_iterations=max_iterations,
                        t3_project_name=t3_project_name,
                        )

    restart_string = "restartFromSeed(path='seed')"
    rmg_input_path = os.path.join(project_directory, 'input.py')
    with open(rmg_input_path) as f:
        content = f.read()
    seed_path = os.path.join(project_directory, 'seed')
    if restart_rmg:
        backup_rmg_files(project_directory=project_directory)
        if restart_string not in content and os.path.isdir(seed_path) and os.listdir(seed_path):
            if os.path.isfile(os.path.join(project_directory, 'restart_from_seed.py')):
                if os.path.isfile(os.path.join(project_directory, 'input.py')):
                    os.rename(src=os.path.join(project_directory, 'input.py'),
                              dst=os.path.join(project_directory, 'input.py.old'))
                os.rename(src=os.path.join(project_directory, 'restart_from_seed.py'),
                          dst=os.path.join(project_directory, 'input.py'))
            elif os.path.isfile(os.path.join(project_directory, 'input.py')):
                with open(os.path.join(project_directory, 'input.py')) as f:
                    content = f.read()
                with open(os.path.join(project_directory, 'input.py'), 'w') as f:
                    f.write("restartFromSeed(path='seed')\n\n" + content)
    job_status, job_id = submit_job(project_directory=project_directory,
                                    logger=logger,
                                    memory=memory,
                                    cluster_soft=LOCAL_CLUSTER_SOFTWARE,
                                    )
    return job_id


def rmg_runner(rmg_input_file_path: str,
               job_log_path: str,
               logger: Logger,
               memory: int | None = None,
               cpus: int | None = None,
               verbose: int | None = None,
               max_iterations: int | None = None,
               t3_project_name: str | None = None,
               rmg_execution_type: str | None = None,
               restart_rmg: bool = False,
               walltime: str | None = None,
               fix_cantera_model: bool = True,
               ) -> bool:
    """
    Run an RMG job as a subprocess under the rmg_env.

    Args:
        rmg_input_file_path (str): The path to the RMG input file.
        job_log_path (str): The path to the ``job.log`` file created on an HTCondor scheduler.
        logger (Logger): The T3 Logger object instance.
        memory (int, optional): The submit script memory in MB.
        cpus (int, optional): The number of CPUs for an RMG parallelization.
        max_iterations (int, optional): Max RMG iterations.
        verbose (int, optional): Level of verbosity.
        t3_project_name (str, optional): The T3 project name, used for setting a job name on the server for the RMG run.
        rmg_execution_type (str, optional): The RMG execution type (incore or local). Also set via settings.py.
        restart_rmg (bool, optional): Whether to restart RMG from seed.
        walltime (str, optional): Max walltime in 'DD:HH:MM:SS' format. Defaults to 6 hours.
        fix_cantera_model (bool, optional): Whether to fix the Cantera model file after the RMG run completes.

    Returns:
        bool: Whether an exception was raised.
    """
    if not os.path.isdir(local_t3_path):
        os.makedirs(local_t3_path)
    new_memory = memory

    rmg_execution_type = rmg_execution_type or RMG_EXECUTION_TYPE
    if rmg_execution_type == 'incore':
        rmg_exception_encountered = run_rmg_incore(rmg_input_file_path=rmg_input_file_path,
                                                   verbose=verbose,
                                                   max_iterations=max_iterations,
                                                   walltime=walltime,
                                                   )
        if fix_cantera_model:
            fix_cantera_model_files(rmg_path=os.path.dirname(rmg_input_file_path))
        return rmg_exception_encountered
    elif rmg_execution_type == 'local':
        runner_counter = 0
        rmg_errors = list()
        converged, run_rmg = False, True
        while run_rmg:
            runner_counter += 1
            project_directory = os.path.abspath(os.path.dirname(rmg_input_file_path))
            job_id = run_rmg_in_local_queue(project_directory=project_directory,
                                            logger=logger,
                                            memory=new_memory,
                                            cpus=cpus,
                                            verbose=verbose,
                                            max_iterations=max_iterations,
                                            restart_rmg=restart_rmg,
                                            t3_project_name=t3_project_name,
                                            )
            while job_id in check_running_jobs_ids(cluster_soft=LOCAL_CLUSTER_SOFTWARE):
                time.sleep(120)
            converged, error = rmg_job_converged(project_directory=project_directory)
            err_path = os.path.join(project_directory, 'err.txt')
            if os.path.isfile(err_path):
                os.rename(err_path, os.path.join(project_directory,
                                                 f'err_{datetime.datetime.now().strftime("%b%d_%Y_%H:%M:%S")}.txt'))
            rmg_errors.append(error)
            if not converged:
                if error is not None:
                    logger.info(f'RMG crashed with the following error:\n{error}')
                new_memory = get_new_memory_for_an_rmg_run(job_log_path,
                                                           logger=logger,
                                                           )
            run_rmg = not converged \
                      and new_memory is not None \
                      and runner_counter < MAX_RMG_RUNS_PER_ITERATION \
                      and not(len(rmg_errors) >= 2 and error is not None and error == rmg_errors[-2])
            restart_rmg = False if error is not None and 'Could not find one or more of the required files/directories ' \
                                                         'for restarting from a seed mechanism' in error else True
        if fix_cantera_model:
            fix_cantera_model_files(rmg_path=os.path.dirname(rmg_input_file_path))
        return not converged
    else:
        logger.warning(f'Expected either "incore" or "local" execution type for RMG, got {rmg_execution_type}.\n'
                       f'Not executing RMG.')
        return True


def get_new_memory_for_an_rmg_run(job_log_path: str,
                                  logger: Logger,
                                  ) -> int | None:
    """
    If an RMG job crashed due to too few or too much memory, compute a new desired memory for the run.
    Note that only on HTCondor there's a cap memory constraint rule that the job must consume at least 20%
    of the allocated memory within the first 30 min of the run, otherwise it is terminated.

    Args:
        job_log_path (str): The path to the ``job.log`` file created on an HTCondor scheduler.
        logger (Logger): The T3 Logger object instance.

    Returns:
        Optional[int]: The recommended memory value in MB.
    """
    global MEM
    new_mem = None
    if os.path.isfile(job_log_path):
        with open(job_log_path) as f:
            lines = f.readlines()
        for line in lines:
            # "Job Is Wasting Memory using less than 20 percent of requested Memory"
            if 'using less than 20 percent of requested' in line:
                mem = None
                for line_ in lines:
                    # 1852  -  MemoryUsage of job (MB)
                    if 'MemoryUsage of job (MB)' in line_:
                        mem = int(line_.split()[0])
                        new_mem = int(mem * 4.5)  # Must be less than mem * 5 to avoid the 20% rule.
                if mem is not None:
                    logger.info(f'RMG job terminated due to 20% memory rule, was {mem / 1000:.2f} GB')
                break
            # "MEMORY EXCEEDED"
            if 'memory exceeded' in line.lower():
                mem = None
                for line_ in lines:
                    # 14361  -  MemoryUsage of job (MB)
                    if 'MemoryUsage of job (MB)' in line_:
                        mem = int(line_.split()[0])
                        new_mem = int(mem * 3)
                if mem is not None:
                    logger.info(f'RMG job terminated since more memory is needed, was {mem / 1000:.2f} GB')
                break
    new_mem = min(new_mem, settings['servers']['local']['max mem'] * 1000) if new_mem is not None else MEM
    logger.info(f'Setting RMG job memory to {new_mem / 1000:.2f} GB')
    return new_mem


def backup_rmg_files(project_directory: str):
    """
    Backup the RMG files before restarting from seed.

    Args:
        project_directory (str): The path to the RMG folder.
    """
    restart_backup_dir = os.path.join(project_directory,
                                      f'restart_backup_{datetime.datetime.now().strftime("%b%d_%Y_%H-%M-%S")}')
    chemkin_folder_path = os.path.join(restart_backup_dir, 'chemkin')
    os.makedirs(chemkin_folder_path, exist_ok=True)
    files = ['RMG.log',
             os.path.join('chemkin', 'chem_annotated.inp'),
             os.path.join('chemkin', 'chem_edge_annotated.inp'),
             ]
    folders = ['pdep']
    for file in files:
        if os.path.exists(os.path.join(project_directory, file)):
            shutil.copy(src=os.path.join(project_directory, file),
                        dst=os.path.join(restart_backup_dir, file))
    for folder in folders:
        if os.path.exists(os.path.join(project_directory, folder)):
            shutil.copytree(src=os.path.join(project_directory, folder),
                            dst=os.path.join(restart_backup_dir, folder))


def fix_cantera_model_files(rmg_path: str) -> None:
    """
    Fix Cantera model files emitted by RMG (resolve mislabeled duplicates, drop
    invalid-rate-coefficient reactions, etc.).

    Args:
        rmg_path (str): The path to the RMG folder.
    """
    fix_cantera(model_path=os.path.join(rmg_path, 'cantera_from_ck', 'chem_annotated.yaml'))
    fix_cantera(model_path=os.path.join(rmg_path, 'cantera_from_ck', 'chem.yaml'))


@dataclass(frozen=True)
class ArkaneJobResult:
    """
    The outcome of ``run_arkane_job``.

    Truth-y exactly when the job succeeded (``__bool__`` below), so every caller that treated the
    old plain-``bool`` return as a truth value keeps working unchanged; what the object adds is the
    DIAGNOSIS a bool could not carry -- most importantly whether a failure was a deadline overrun
    (``timed_out``), which a caller must be able to distinguish from "Arkane ran and failed"
    without parsing a reason string.

    Args:
        succeeded (bool): Whether the job succeeded (produced its required artifact, with no
                          run/timeout failure).
        timed_out (bool): Whether the job overran its ``timeout`` deadline and was killed.
        reason (str, optional): A human-readable diagnosis when ``succeeded`` is False.
    """
    succeeded: bool
    timed_out: bool = False
    reason: str | None = None

    def __bool__(self) -> bool:
        """
        Returns:
            bool: Whether the job succeeded.
        """
        return self.succeeded


# The interpreter script the ``timeout`` path of ``run_arkane_job`` runs in its own process:
# exactly the same ARC entry point (``arc.statmech.arkane.run_arkane``) the in-process path calls,
# so the two paths cannot drift on HOW Arkane is invoked -- only on WHERE (a killable session vs.
# this process).
_ARKANE_SUBPROCESS_SCRIPT = (
    'import sys\n'
    'from arc.statmech.arkane import run_arkane\n'
    'run_arkane(statmech_dir=sys.argv[1])\n'
)

# How long a timed-out Arkane process group is given to exit on SIGTERM before it is SIGKILLed.
_ARKANE_KILL_GRACE_SECONDS = 10


def _spawn_arkane_subprocess(output_directory: str) -> subprocess.Popen:
    """
    Spawn ``arc.statmech.arkane.run_arkane`` in its own interpreter and its own process SESSION.

    ``start_new_session=True`` is the load-bearing part: ARC's ``run_arkane`` shells out through
    ``bash -lc`` and ``conda run`` (arc/statmech/arkane.py), so the actual Arkane interpreter is a
    grandchild -- killing only the direct child would orphan it, still writing into the run
    directory. A fresh session makes the child a process-group leader whose pgid covers every
    descendant, so ``_kill_process_group`` can take the whole tree down at once.

    Kept as a module-level helper (rather than inlined) so tests can substitute a controllable
    process without re-implementing the wait/kill/reap logic under test.

    Args:
        output_directory (str): The directory containing ``input.py``, passed to ``run_arkane``.

    Returns:
        subprocess.Popen: The spawned process (its pid == its pgid, per the new session).
    """
    return subprocess.Popen(
        [sys.executable, '-c', _ARKANE_SUBPROCESS_SCRIPT, output_directory],
        start_new_session=True,
    )


def _kill_process_group(process: subprocess.Popen) -> None:
    """
    Terminate a spawned process's entire process group and reap the direct child.

    SIGTERM first (Arkane/Python get a chance to flush logs), SIGKILL after
    ``_ARKANE_KILL_GRACE_SECONDS``. Every ``ProcessLookupError`` window (the group exiting between
    our decision and our signal) is tolerated: the goal state "nothing from this run is still
    running" is then already true. The final ``wait()`` reaps the zombie so the caller can observe
    ``process.poll() is not None``.

    Args:
        process (subprocess.Popen): The process to kill (spawned with ``start_new_session=True``).
    """
    try:
        pgid = os.getpgid(process.pid)
    except ProcessLookupError:
        process.wait()
        return
    try:
        os.killpg(pgid, signal.SIGTERM)
    except ProcessLookupError:
        process.wait()
        return
    try:
        process.wait(timeout=_ARKANE_KILL_GRACE_SECONDS)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            # The group exited on its own between the SIGTERM grace wait above and this SIGKILL --
            # the goal state ("nothing from this run is still running") is already true, so there is
            # nothing to signal and nothing to handle.
            pass
        process.wait()


def run_arkane_job(input_file: str,
                   output_directory: str,
                   logger: Logger | None = None,
                   required_artifact: str = os.path.join('sensitivity', 'sa_coefficients.yml'),
                   timeout: float | None = None,
                   ) -> ArkaneJobResult:
    """
    Run an Arkane job.

    With ``timeout=None`` (the default), ARC's ``run_arkane`` is called in-process, exactly as
    before. With a ``timeout``, the SAME ARC entry point runs in its own interpreter in its own
    process session instead, because an in-process call cannot be killed at all: ARC reaches
    ``subprocess.run`` (arc/job/local.py) with no ``timeout=`` and no Popen/PID exposed to any
    caller, and a Python thread cannot be killed either -- an earlier in-adapter attempt to wrap
    the call in a ThreadPoolExecutor could only relabel the outcome while Arkane kept running and
    kept writing into a run directory already declared failed. Spawning the process ourselves is
    what makes a deadline enforceable: on overrun the whole process group (the ``conda run`` bash
    wrapper AND the Arkane interpreter under it) is SIGTERMed, then SIGKILLed, then reaped.

    A deadline overrun is reported as a falsy result with ``timed_out=True`` -- never raised --
    so a caller assembling a run record gets a recordable failure, not an exception that destroys
    the record.

    Args:
        input_file (str): The path to the Arkane input file.
        output_directory (str): The path to the output directory.
        logger (Logger, optional): The logger object.
        required_artifact (str, optional): A path, relative to ``output_directory``, naming the
                                           artifact that Arkane must (re)write for the job to be
                                           considered successful. Deleted before the run (if
                                           present) and required to exist afterward. Defaults to
                                           the SA job's ``sensitivity/sa_coefficients.yml``, so
                                           this is a no-op for all existing callers.
        timeout (float, optional): The wall-clock deadline, in seconds, for the Arkane process.
                                   ``None`` (the default) preserves the historical in-process,
                                   unbounded call.

    Raises:
        ValueError: If ``required_artifact`` is an absolute path, or resolves (after ``..`` and
                    symlink resolution) to a location outside ``output_directory``. The joined
                    path is deleted before the run, so an unconfined value would delete an
                    arbitrary file elsewhere on the filesystem. Also if ``timeout`` is given and
                    is not a positive finite number of seconds -- an unenforceable deadline
                    (zero, negative, inf, nan, or a non-number) silently meaning "no deadline"
                    would be worse than refusing it.

    Returns:
        ArkaneJobResult: The outcome; truthy exactly when the job succeeded, so callers treating
                        it as the historical plain bool keep working.
    """
    from arc.statmech.arkane import run_arkane

    # Checked before ANY side effect (the artifact deletion below), like the required_artifact
    # confinement checks: an invalid deadline is an argument error about this call, not a run
    # outcome. bool is excluded explicitly because it is an int subclass, so timeout=True would
    # otherwise be accepted as a 1-second deadline nobody asked for.
    if timeout is not None:
        if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) \
                or not math.isfinite(timeout) or timeout <= 0:
            raise ValueError(f"The 'timeout' argument must be a positive finite number of seconds, or None "
                             f"for no deadline, got {timeout!r} of type {type(timeout).__name__}.")

    # Confine the artifact path to output_directory BEFORE any side effect: the joined path is
    # deleted below, so a traversal value ('../important.yml') or an absolute value (which makes
    # os.path.join discard output_directory entirely) would otherwise delete a file outside the
    # job directory.
    if os.path.isabs(required_artifact):
        raise ValueError(f"The 'required_artifact' argument must be a path relative to the output "
                         f"directory, got the absolute path '{required_artifact}'.")
    resolved_output_directory = os.path.realpath(output_directory)
    artifact_path = os.path.join(output_directory, required_artifact)
    resolved_artifact_path = os.path.realpath(artifact_path)
    if resolved_artifact_path == resolved_output_directory \
            or os.path.commonpath([resolved_output_directory, resolved_artifact_path]) != resolved_output_directory:
        raise ValueError(f"The 'required_artifact' argument must resolve to a path strictly inside "
                         f"the output directory '{output_directory}', got '{required_artifact}' "
                         f"which resolves to '{resolved_artifact_path}'.")

    # Ensure output directory exists
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # run_arkane expects the input file to be named 'input.py' inside the directory.
    target_input = os.path.join(output_directory, 'input.py')

    if os.path.abspath(input_file) != os.path.abspath(target_input):
        shutil.copyfile(input_file, target_input)

    # Check for success by looking for the required artifact, which is the actual product of a
    # successful Arkane job. Neither output.py nor a stale artifact is evidence of anything: both
    # can survive from a previous run. Rather than relying on a timestamp comparison (racy under
    # coarse filesystem mtime granularity, and meaningless if the clock or filesystem lies),
    # delete any pre-existing file before invoking Arkane and simply require that Arkane itself
    # (re)wrote it. ``artifact_path`` was computed (and confined to ``output_directory``) above.
    if os.path.isfile(artifact_path):
        os.remove(artifact_path)

    if timeout is None:
        try:
            run_arkane(statmech_dir=output_directory)
        except Exception as e:
            reason = f'Arkane run failed with error: {e}'
            if logger:
                logger.error(reason)
            return ArkaneJobResult(succeeded=False, reason=reason)
    else:
        process = _spawn_arkane_subprocess(output_directory=output_directory)
        try:
            process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            _kill_process_group(process=process)
            reason = (f'Arkane run in {output_directory} timed out after {timeout} s; its process '
                      f'group was killed (SIGTERM, then SIGKILL after {_ARKANE_KILL_GRACE_SECONDS} s) '
                      f'so nothing keeps writing into the run directory after this verdict.')
            if logger:
                logger.error(reason)
            return ArkaneJobResult(succeeded=False, timed_out=True, reason=reason)
        if process.returncode != 0:
            # Mirrors the in-process path's except-clause: a run that died is a failure regardless
            # of what artifacts it may have left behind. NOT labelled a timeout -- the deadline was
            # met, the process failed on its own.
            reason = f'Arkane run failed (exit status {process.returncode}).'
            if logger:
                logger.error(reason)
            return ArkaneJobResult(succeeded=False, reason=reason)

    if not os.path.isfile(artifact_path):
        reason = f'The Arkane job in {output_directory} did not produce {artifact_path}.'
        if logger:
            logger.error(reason)
        return ArkaneJobResult(succeeded=False, reason=reason)
    return ArkaneJobResult(succeeded=True)


def run_rmg_sa_incore(rmg_input_file_path: str,
                      chemkin_file_path: str,
                      species_dict_path: str,
                      output_path: str,
                      observables: list[str] | None = None,
                      threshold: float = 1e-3,
                      ) -> tuple[bool, str | None]:
    """
    Run RMG Sensitivity Analysis incore under the rmg_env.
    """
    project_directory = os.path.abspath(os.path.dirname(rmg_input_file_path))
    script_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'rmg_incore_sa.py')

    rmg_input_file_path = os.path.abspath(rmg_input_file_path)
    chemkin_file_path = os.path.abspath(chemkin_file_path)
    species_dict_path = os.path.abspath(species_dict_path)
    output_path = os.path.abspath(output_path)

    obs_str = ""
    if observables:
        obs_str = f"-obs {' '.join(shlex.quote(o) for o in observables)}"

    inner_cmd = (f'python {script_path} '
                 f'-i {rmg_input_file_path} '
                 f'-c {chemkin_file_path} '
                 f'-d {species_dict_path} '
                 f'-o {output_path} '
                 f'-t {threshold} '
                 f'{obs_str} '
                 f'> >(tee -a sa_out.txt) 2> >(tee -a sa_err.txt >&2)')
    shell_script = rf'''bash -lc 'set -uo pipefail
cd "{project_directory}"
if command -v micromamba >/dev/null 2>&1; then
    micromamba run -n rmg_env bash -c "{inner_cmd}"
elif command -v conda >/dev/null 2>&1 || command -v mamba >/dev/null 2>&1; then
    conda run -n rmg_env bash -c "{inner_cmd}"
else
    echo "Micromamba/Mamba/Conda required" >&2
    exit 1
fi' '''

    execute_command(shell_script, shell=True, no_fail=True, executable='/bin/bash')

    if os.path.isfile(output_path):
        return True, None

    error_msg = "Unknown error"
    err_file = os.path.join(project_directory, 'sa_err.txt')
    if os.path.isfile(err_file):
        with open(err_file) as f:
            error_msg = f.read()
    return False, error_msg
