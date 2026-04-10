"""
A "keep alive" runner tool for RMG on a server.
Should be executed locally on the head node using the t3 environment.
"""

import datetime
import os
import shlex
import shutil
import time
from typing import TYPE_CHECKING, List, Optional, Tuple

from arc.job.local import (_determine_job_id,
                           change_mode,
                           execute_command,
                           parse_running_jobs_ids)

from t3.imports import local_t3_path, settings, submit_scripts

if TYPE_CHECKING:
    from t3.logger import Logger


MEM = settings['rmg_initial_memory'] * 1000  # MB
SLEEP_TIME = 6  # hours
MAX_RMG_RUNS_PER_ITERATION = 5

RMG_EXECUTION_TYPE = settings['execution_type']['rmg']

if RMG_EXECUTION_TYPE == 'local':
    LOCAL_CLUSTER_SOFTWARE = settings['servers']['local']['cluster_soft']
    SUBMIT_COMMAND = settings['submit_command'][LOCAL_CLUSTER_SOFTWARE]
    CHECK_STATUS_COMMAND = settings['check_status_command'][LOCAL_CLUSTER_SOFTWARE]
    SUBMIT_FILENAME = settings['submit_filenames'][LOCAL_CLUSTER_SOFTWARE]
else:
    SUBMIT_COMMAND = CHECK_STATUS_COMMAND = SUBMIT_FILENAME = LOCAL_CLUSTER_SOFTWARE = ''


def write_submit_script(project_directory: str,
                        cpus: Optional[int] = None,
                        memory: Optional[int] = None,
                        verbose: Optional[str] = None,
                        max_iterations: Optional[str] = None,
                        t3_project_name: Optional[str] = None,
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
    """
    global MEM
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
               logger: 'Logger',
               cluster_soft: str,
               memory: Optional[int] = None,
               ) -> Tuple[Optional[str], Optional[str]]:
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


def check_running_jobs_ids(cluster_soft: str) -> List[str]:
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


def rmg_job_converged(project_directory: str) -> Tuple[bool, Optional[str]]:
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
        with open(rmg_log_path, 'r') as f:
            lines = f.readlines()
            len_lines = len(lines)
            for i in range(10):
                if 'MODEL GENERATION COMPLETED' in lines[len_lines - 1 - i]:
                    rmg_converged = True
                    break
    if not rmg_converged and os.path.isfile(rmg_err_path):
        with open(rmg_err_path, 'r') as f:
            lines = f.readlines()
        for line in lines[::-1]:
            if 'Error' in line:
                error = line.strip()
                break
    return rmg_converged, error


_DEFAULT_RMG_TIMEOUT_S = 6 * 3600  # 6 hours

logger = logging.getLogger(__name__)


def _parse_walltime_to_seconds(walltime: str) -> int:
    """Parse a 'DD:HH:MM:SS' walltime string to total seconds. Returns 0 for '00:00:00:00'."""
    parts = walltime.split(':')
    if len(parts) != 4:
        return 0
    days, hours, minutes, seconds = (int(p) for p in parts)
    return days * 86400 + hours * 3600 + minutes * 60 + seconds


def run_rmg_incore(rmg_input_file_path: str,
                   verbose: Optional[int] = None,
                   max_iterations: Optional[int] = None,
                   walltime: Optional[str] = None,
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
    try:
        result = subprocess.run(shell_script, shell=True, executable='/bin/bash',
                                capture_output=True, text=True, timeout=timeout_s)
        stderr_text = result.stderr or ''
    except subprocess.TimeoutExpired:
        logger.error(f'RMG incore timed out after {timeout_s}s')
        return True
    if 'RMG threw an exception and did not converge.' in stderr_text:
        return True
    return False


def run_rmg_in_local_queue(project_directory: str,
                           logger: 'Logger',
                           memory: Optional[int] = None,
                           cpus: Optional[int] = None,
                           max_iterations: Optional[int] = None,
                           restart_rmg: bool = False,
                           verbose: Optional[int] = None,
                           t3_project_name: Optional[str] = None,
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
    with open(rmg_input_path, 'r') as f:
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
                with open(os.path.join(project_directory, 'input.py'), 'r') as f:
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
               logger: 'Logger',
               memory: Optional[int] = None,
               cpus: Optional[int] = None,
               verbose: Optional[int] = None,
               max_iterations: Optional[int] = None,
               t3_project_name: Optional[str] = None,
               rmg_execution_type: Optional[str] = None,
               restart_rmg: bool = False,
               walltime: Optional[str] = None,
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
        return rmg_exception_encountered
    elif rmg_execution_type == 'local':
        runner_counter = 0
        rmg_errors = list()
        converged, restart_rmg, run_rmg = False, restart_rmg, True
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
        return not converged
    else:
        logger.warning(f'Expected wither "incore" or "local" execution type for RMG, got {rmg_execution_type}.\n'
                       f'Not executing RMG.')
        return True


def get_new_memory_for_an_rmg_run(job_log_path: str,
                                  logger: 'Logger',
                                  ) -> Optional[int]:
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
        with open(job_log_path, 'r') as f:
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
    os.mkdir(restart_backup_dir)
    os.mkdir(os.path.join(restart_backup_dir, 'chemkin'))
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


def run_arkane_job(input_file: str,
                   output_directory: str,
                   plot: bool = False,
                   logger: Optional['Logger'] = None,
                   ) -> bool:
    """
    Run an Arkane job.

    Args:
        input_file (str): The path to the Arkane input file.
        output_directory (str): The path to the output directory.
        plot (bool, optional): Whether to plot the results.
        logger (Logger, optional): The logger object.

    Returns:
        bool: Whether the job was successful.
    """
    from arc.statmech.arkane import run_arkane

    # Ensure output directory exists
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # run_arkane expects the input file to be named 'input.py' inside the directory.
    target_input = os.path.join(output_directory, 'input.py')

    if os.path.abspath(input_file) != os.path.abspath(target_input):
        shutil.copyfile(input_file, target_input)

    try:
        run_arkane(statmech_dir=output_directory)
    except Exception as e:
        if logger:
            logger.error(f'Arkane run failed with error: {e}')
        return False

    # Check for success by looking for the sensitivity output directory,
    # which is the actual product of a successful Arkane SA job.
    # Do not rely on output.py as it may pre-exist from previous runs.
    sa_dir = os.path.join(output_directory, 'sensitivity')
    if os.path.isdir(sa_dir) and any(f.endswith('.yml') or f.endswith('.yaml') for f in os.listdir(sa_dir)):
        return True
    return False


def run_rmg_sa_incore(rmg_input_file_path: str,
                      chemkin_file_path: str,
                      species_dict_path: str,
                      output_path: str,
                      observables: Optional[List[str]] = None,
                      threshold: float = 1e-3,
                      ) -> Tuple[bool, Optional[str]]:
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
        with open(err_file, 'r') as f:
            error_msg = f.read()
    return False, error_msg
