"""
A "keep alive" runner tool for RMG on a server.
Should be executed locally on the head node using the t3 environment.
"""

from typing import TYPE_CHECKING, List, Optional, Tuple

import os
import time

from arc.job.local import _determine_job_id, change_mode, execute_command, parse_running_jobs_ids, submit_job

from t3.imports import local_t3_path, settings, submit_scripts

if TYPE_CHECKING:
    from t3.logger import Logger


CPUS = 16  # A recommended value for RMG when running on a server (not incore)
MEM = 10000  # MB
SLEEP_TIME = 6  # hours

rmg_execution_type = settings['execution_type']['rmg']

if rmg_execution_type == 'local':
    LOCAL_CLUSTER_SOFTWARE = settings['servers']['local']['cluster_soft']
    SUBMIT_COMMAND = settings['submit_command'][LOCAL_CLUSTER_SOFTWARE]
    CHECK_STATUS_COMMAND = settings['check_status_command'][LOCAL_CLUSTER_SOFTWARE]
    SUBMIT_FILENAME = settings['submit_filenames'][LOCAL_CLUSTER_SOFTWARE]
else:
    SUBMIT_COMMAND = CHECK_STATUS_COMMAND = SUBMIT_FILENAME, LOCAL_CLUSTER_SOFTWARE = ''


def write_submit_script(project_directory: str,
                        cpus: Optional[int] = None,
                        memory: Optional[int] = None,
                        verbose: Optional[str] = None,
                        max_iterations: Optional[str] = None,
                        ) -> None:
    """
    Write an RMG submit script.

    Args:
        project_directory (str): The full path to the project directory.
        cpus (int, optional): The number of CPUs for an RMG parallelization, defaults to ``CPUS``.
        memory (int, optional): The memory in MB for an RMG run, defaults to ``MEM``.
        verbose (str, optional): Level of verbosity, e.g., ``-v 10``.
        max_iterations (str, optional): Max RMG iterations, e.g., ``-m 100``.
    """
    global MEM
    submit_scripts_content = submit_scripts['rmg'].format(name='T3_RMG',
                                                          cpus=cpus or CPUS,
                                                          memory=memory or MEM,
                                                          )
    with open(os.path.join(project_directory, SUBMIT_FILENAME), 'w') as f:
        f.write(submit_scripts_content)
    if 'rmg_job' in submit_scripts.keys():
        # Write an aux submit script, e.g., as required for HTCondor.
        verbose = verbose or ''
        max_iterations = max_iterations or ''
        aux_submit_scripts_content = submit_scripts['rmg_job'].format(cpus=cpus or CPUS,
                                                                      max_iterations=max_iterations,
                                                                      )
        with open(os.path.join(project_directory, 'job.sh'), 'w') as f:
            f.write(aux_submit_scripts_content)
            change_mode(mode='+x', file_name='job.sh', path=project_directory)


def submit_job(project_directory: str,
               logger: 'Logger',
               cluster_soft: str,
               ) -> Tuple[Optional[str], Optional[str]]:
    """
    Submit an RMG job.

    Args:
        project_directory (str): The job (folder) name.
        logger (Logger): The T3 Logger object instance.
        cluster_soft (str): The server's cluster software.

    Returns:
        Tuple[Optional[str], Optional[str]]: job_status, job_id
    """
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
    logger.info(f'Successfully submitted job {project_directory}, job ID = {job_id}.')
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


def rmg_job_converged(project_directory: str) -> bool:
    """
    Determine whether an RMG job has converged.

    Args:
        project_directory (str): The job (folder) name.

    Returns:
        bool: Whether this RMG run has converged.
    """
    rmg_converged = False
    rmg_log_path = os.path.join(project_directory, 'RMG.log')
    if os.path.isfile(rmg_log_path):
        with open(rmg_log_path, 'r') as f:
            lines = f.readlines()
            len_lines = len(lines)
            for i in range(10):
                if 'MODEL GENERATION COMPLETED' in lines[len_lines - 1 - i]:
                    rmg_converged = True
                    break
    return rmg_converged


def write_restart_file(name: str,  # Todo: implement this
                       logger: 'Logger',
                       ) -> None:
    """
    Convert an RMG input file into an RMG restart file.

    Args:
        name (str): The job (folder) name.
        logger (Logger): The T3 Logger object instance.
    """
    restart_string = "restartFromSeed(path='seed')"
    rmg_input_path = os.path.join(name, 'input.py')
    with open(rmg_input_path, 'r') as f:
        content = f.read()
    if restart_string not in content:
        logger.info(f'Converting the RMG input file of {name} into an RMG restart file')
        content = f'{restart_string}\n\n{content}'
        with open(rmg_input_path, 'w') as f:
            f.write(content)


# def get_names_by_sub_folders(pwd: str) -> List[str]:
#     """
#     Get the names of the runs.
#
#     Args:
#         pwd (str): The present working directory.
#
#     Returns:
#         List[str]: the names of all runs.
#     """
#     names = list()
#     for _, folders, _ in os.walk(pwd):
#         for folder in folders:
#             if folder[0] == 'x':
#                 names.append(folder)
#         # Don't continue to sub folders.
#         break
#     return sorted(names)


# def initialize_rmg_job(names: List[str],   # rewrite for a single RMG job, wait for it to finish, trsh mem if needed
#                        job_id_yml_path: str,
#                        logger: 'Logger'
#                        ) -> Tuple[Dict[str, bool], Dict[str, str], Dict[str, str]]:
#     """
#     Initialize the RMG job.
#
#     Args:
#         names (List[str]): The run names.
#         job_id_yml_path (str): The path to the job ID YAML file.
#         logger (Logger): the T3 Logger object instance.
#
#     Returns:
#         Tuple[Dict[str, bool], Dict[str, str], Dict[str, str]]: convergence, status, job_ids.
#     """
#     logger.info('\n\ninitializing RMG job...')
#     convergence, status, job_ids = dict(), dict(), dict()
#     if os.path.isfile(job_id_yml_path):
#         job_ids = read_yaml_file(job_id_yml_path)
#     server_job_ids = check_running_jobs_ids()
#     for name in names:
#         if name in job_ids.keys() and job_ids[name] in server_job_ids:
#             logger.info(f'Job {name} is already running (index {job_ids[name]}).')
#             convergence[name] = False
#             continue
#         logger.info(f'Initializing {name}')
#         if not os.path.isfile(os.path.join(name, SUBMIT_FILENAME)):
#             logger.info(f'Writing submit script for {name}')
#             write_submit_script(name)
#         if not rmg_job_converged(name):
#             job_status, job_id = submit_job(name=name, logger=logger)  # note: not writing restart file
#             convergence[name] = False
#         else:
#             logger.info(f'Job {name} already converged, not initializing it')
#             job_status, job_id = '', 0
#             convergence[name] = True
#         status[name] = job_status
#         job_ids[name] = job_id
#     save_yaml_file(job_id_yml_path, job_ids)
#     return convergence, status, job_ids


# def update_names_and_dicts(names: List[str],
#                            convergence: Dict[str, bool],
#                            status: Dict[str, str],
#                            job_ids: Dict[str, str],
#                            ) -> Tuple[Dict[str, bool], Dict[str, str], Dict[str, str]]:
#     """
#     Check whether new folders were added, and update the data dictionaries accordingly.
#
#     Args:
#         names (List[str]): Updated list of job names / folders.
#         convergence (Dict[str, bool]): convergence.
#         status (Dict[str, str]): status.
#         job_ids (Dict[str, str]): job IDs.
#
#     Returns:
#         Tuple[Dict[str, bool], Dict[str, str], Dict[str, str]]: convergence, status, job_ids.
#     """
#     for name in names:
#         if name not in convergence.keys():
#             convergence[name] = False
#         if name not in status.keys():
#             status[name] = ''
#         if name not in job_ids.keys():
#             job_ids[name] = 0
#     return convergence, status, job_ids


def run_rmg_incore(rmg_input_file_path: str,
                   verbose: Optional[int] = None,
                   max_iterations: Optional[int] = None,
                   ) -> bool:
    """
    Run RMG incore under the rmg_env.

    Args:
        rmg_input_file_path (str): The path to the RMG input file.
        max_iterations (int, optional): Max RMG iterations.
        verbose (int, optional): Level of verbosity.

    Returns:
        bool: Whether an exception was raised.
    """
    project_directory = os.path.abspath(os.path.dirname(rmg_input_file_path))
    verbose = f' -v {verbose}' if verbose is not None else ''
    max_iterations = f' -m {max_iterations}' if max_iterations is not None else ''
    script_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'rmg_incore_script.py')
    commands = ['CONDA_BASE=$(conda info --base)',
                'source $CONDA_BASE/etc/profile.d/conda.sh',
                'conda activate rmg_env',
                f'cd {project_directory}',
                f'python-jl {script_path} {rmg_input_file_path}{verbose}{max_iterations} '
                f'> >(tee -a out.txt) 2> >(tee -a err.txt >&2)',
                ]
    stdout, stderr = execute_command(commands, shell=True, no_fail=True, executable='/bin/bash')
    if 'RMG threw an exception and did not converge.\n' in stderr:
        return True
    return False


def run_rmg_in_local_queue(project_directory: str,
                           logger: 'Logger',
                           memory: Optional[int] = None,
                           verbose: Optional[int] = None,
                           max_iterations: Optional[int] = None,
                           ):
    """
    Run RMG on the queue of the local server (under the rmg_env).

    Args:
        project_directory (str): The path to the RMG folder.
        logger (Logger): The T3 Logger object instance.
        memory (int, optional): The submit script memory in MB.
        max_iterations (int, optional): Max RMG iterations.
        verbose (int, optional): Level of verbosity.

    Returns:
        Optional[str]: The job ID.
    """
    verbose = f' -v {verbose}' if verbose is not None else ''
    max_iterations = f' -m {max_iterations}' if max_iterations is not None else ''

    write_submit_script(project_directory=project_directory,
                        cpus=settings['servers']['local']['cpus'],
                        memory=memory,
                        verbose=verbose,
                        max_iterations=max_iterations,
                        )

    job_status, job_id = submit_job(project_directory=project_directory,
                                    logger=logger,
                                    cluster_soft=LOCAL_CLUSTER_SOFTWARE,
                                    )
    return job_id


def rmg_runner(rmg_input_file_path: str,
               job_log_path: str,
               logger: 'Logger',
               verbose: Optional[int] = None,
               max_iterations: Optional[int] = None,
               ) -> bool:
    """
    Run an RMG job as a subprocess under the rmg_env.

    Args:
        rmg_input_file_path (str): The path to the RMG input file.
        job_log_path (str): The path to the ``job.log`` file created on an HTCondor scheduler.
        logger (Logger): The T3 Logger object instance.
        max_iterations(int, optional): Max RMG iterations.
        verbose(int, optional): Level of verbosity.

    Returns:
        bool: Whether an exception was raised.
    """
    if not os.path.isdir(local_t3_path):
        os.makedirs(local_t3_path)
    memory_handler = False
    new_memory = None
    converged = True

    if rmg_execution_type == 'incore':
        rmg_exception_encountered = run_rmg_incore(rmg_input_file_path=rmg_input_file_path,
                                                   verbose=verbose,
                                                   max_iterations=max_iterations,
                                                   )
        return rmg_exception_encountered
    elif rmg_execution_type == 'local':
        while not memory_handler:
            project_directory = os.path.abspath(os.path.dirname(rmg_input_file_path))
            job_id = run_rmg_in_local_queue(project_directory=project_directory,
                                            logger=logger,
                                            memory=new_memory,
                                            verbose=verbose,
                                            max_iterations=max_iterations,
                                            )
            while job_id in check_running_jobs_ids(cluster_soft=LOCAL_CLUSTER_SOFTWARE):
                time.sleep(120)
            converged = rmg_job_converged(project_directory=project_directory)
            if not converged:
                new_memory = get_new_memory_for_an_rmg_run(job_log_path)
                if new_memory is None:
                    memory_handler = True

        return not converged

        # job_id_yml_path = os.path.join(local_t3_path, 'jobs.yml')
        # convergence, status, job_ids = initialize_rmg_job(names, job_id_yml_path, logger)
        #
        # while any(not conv for conv in convergence.values()):
        #     logger.info('\n\n\nlooping...')
        #     server_job_ids = check_running_jobs_ids()
        #     names = get_names_by_sub_folders(pwd)
        #     convergence, status, job_ids = update_names_and_dicts(names, convergence, status, job_ids)
        #     for name in names:
        #         job_id = job_ids[name]
        #         if job_id not in server_job_ids:
        #             logger.info(f'RMG job {name} {job_id} terminated')
        #             rmg_converged = rmg_job_converged(name)
        #             if rmg_converged:
        #                 logger.info(f'RMG job {name} {job_id} has converged!!!')
        #                 convergence[name] = True
        #                 continue
        #             if os.path.isfile(os.path.join(name, 'RMG.log')):
        #                 logger.info(f'RMG job {name} {job_id} did not converge')
        #                 write_restart_file(name)
        #                 logger.info(f'Restarting RMG job {name}')
        #             else:
        #                 logger.info(f'Running RMG job {name}')
        #             job_status, job_id = submit_job(name=name, logger=logger)
        #             status[name] = job_status
        #             job_ids[name] = job_id
        #         elif rmg_job_converged(name):
        #             convergence[name] = True
        #     currently_running_jobs = list()
        #     logger.info(f'\n\nConvergence: {convergence}\n\n')
        #     save_yaml_file(job_id_yml_path, job_ids)
        #     logger.info(f'Sleeping for {SLEEP_TIME} hours. ZZZ... ZZZ...')
        #     time.sleep(SLEEP_TIME * 60 * 60)


def get_new_memory_for_an_rmg_run(job_log_path) -> Optional[int]:
    """
    If an RMG job crashed due to too few or too much memory, compute a new desired memory for the run.
    Note that only on HTCondor there's a cap memory constraint rule that the job must consume at least 20%
    of the allocated memory within the first 30 min of the run, otherwise it is terminated.

    Args:
        job_log_path (str): The path to the ``job.log`` file created on an HTCondor scheduler.

    Returns:
        Optional[int]: The recommended memory value in MB.
    """
    new_mem = None
    if os.path.isfile(job_log_path):
        with open(job_log_path, 'r') as f:
            lines = f.readlines()
        for line in lines:
            #	Job Is Wasting Memory using less than 20 percent of requested Memory
            #	Code 26 Subcode 0
            if 'using less than 20 percent of requested' in line or 'Code 26 Subcode 0' in line:
                for line_ in lines:
                    #	1852  -  MemoryUsage of job (MB)
                    if 'MemoryUsage of job (MB)' in line_:
                        new_mem = int(int(line_.split()[0]) * 4.5)
                break
    new_mem = min(new_mem, settings['servers']['local']['max mem'] * 1000) if new_mem is not None else None
    return new_mem
