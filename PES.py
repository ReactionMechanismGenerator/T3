#!/usr/bin/env python3
# encoding: utf-8

"""
PES executable module

Drives the standalone PES exploration loop (``t3.pdep.pes_loop.run_pes_loop``): explore a
pressure-dependent network, queue its sensitive transition states to ARC for QM, re-explore with
the computed barriers, and repeat. A sibling to ``T3.py``, so the loop can be driven externally to
T3's iteration machinery.
"""

import argparse
import datetime
import os
import sys

from arc.common import read_yaml_file

from t3.logger import Logger
from t3.pdep.pes_loop import PES_LOOP_FAILED, run_pes_loop
from t3.pdep.pes_qm import arc_qm_runner
from t3.schema import PESLoopConfig
from t3.utils.dependencies import check_dependencies


def read_pes_input(path: str) -> PESLoopConfig:
    """
    Read and validate a PES loop input file.

    Args:
        path (str): The path to the YAML input file.

    Returns:
        PESLoopConfig: The validated configuration.

    Raises:
        FileNotFoundError: If ``path`` does not exist. Checked here rather than left to the YAML
                           reader so a mistyped path fails with the path in the message instead of
                           an opaque parser error.
        pydantic.ValidationError: If the input file does not validate. ``PESLoopConfig`` forbids
                                  extra top-level keys, so a typo'd section is refused loudly
                                  rather than silently ignored.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f'Could not find the PES input file {path}.')
    input_dict = read_yaml_file(path=path) or dict()
    return PESLoopConfig(**input_dict)


def parse_command_line_arguments(command_line_args=None):
    """
    Parse command-line arguments.

    Args:
        command_line_args: The command line arguments.

    Returns:
        The parsed command-line arguments by key words.
    """

    parser = argparse.ArgumentParser(description='The PES exploration loop')
    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
                        help='a file describing the PES exploration loop to execute')
    parser.add_argument('-p', '--project-directory', type=str, default=None,
                        help="the directory to run in, defaults to the input file's own directory")

    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-d', '--debug', action='store_true', help='print debug information')
    group.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')

    args = parser.parse_args(command_line_args)
    args.file = args.file[0]

    return args


def exit_code_for(status: str) -> int:
    """
    The process exit code corresponding to a ``PESLoopResult`` status.

    Only an outright failure is an error: a stall, a max_rounds stop, or a round with no candidates
    are findings to read in the log, not crashes, and a wrapper script must not treat them as such.

    Args:
        status (str): A ``t3.pdep.pes_loop`` ``PES_LOOP_*`` status constant.

    Returns:
        int: ``1`` for ``PES_LOOP_FAILED``, ``0`` for every other status.
    """
    return 1 if status == PES_LOOP_FAILED else 0


def _resolve(path: str, directory: str) -> str:
    """
    Resolve one path from the input file against the input file's own directory.

    Args:
        path (str): A path as written in the input file, absolute or relative.
        directory (str): The input file's directory.

    Returns:
        str: ``path`` if it is already absolute, otherwise its absolute form anchored at
        ``directory``.
    """
    return path if os.path.isabs(path) else os.path.abspath(os.path.join(directory, path))


def verbose_level(args) -> int:
    """
    The logging level the command-line verbosity flags ask for.

    Args:
        args: The parsed command-line arguments.

    Returns:
        int: ``10`` for ``-d/--debug``, ``30`` for ``-q/--quiet``, ``20`` otherwise -- the same
        mapping ``T3.py`` uses.
    """
    if args.debug:
        return 10
    if args.quiet:
        return 30
    return 20


def main():
    """
    The main PES executable function.

    Returns the ``PESLoopResult`` rather than calling ``sys.exit`` so it is callable (and testable)
    in-process; only the ``__main__`` guard turns the status into an exit code.

    Returns:
        PESLoopResult: The outcome of the loop.

    Raises:
        FileNotFoundError: If the input file, or the network file it names, does not exist. The
                           network is checked here rather than left to the explorer so a mistyped
                           or wrongly-anchored path fails immediately, with the resolved path in
                           the message, instead of raising from inside ``t3/pdep/parser.py`` after
                           a project directory and a log file have already been created.
    """
    args = parse_command_line_arguments()
    input_file = os.path.abspath(args.file)
    project_directory = os.path.abspath(args.project_directory) if args.project_directory \
        else os.path.abspath(os.path.dirname(input_file))

    config = read_pes_input(input_file)
    input_directory = os.path.dirname(input_file)
    # A relative network path is resolved against the input file's own directory here: left alone
    # it would only fail deep inside an Arkane run, long after the input was accepted.
    if not os.path.isabs(config.pes.network):
        config.pes.network = _resolve(config.pes.network, input_directory)
    # Every path in the input file is anchored the SAME way. A relative reuse path left unresolved
    # does not fail loudly -- it silently adopts nothing (the prior project simply is not there
    # from the loop's cwd), and the loop pays for it with a redundant round of real quantum
    # chemistry.
    config.reuse.from_t3_projects = [_resolve(path, input_directory)
                                     for path in config.reuse.from_t3_projects]
    if not os.path.isfile(config.pes.network):
        raise FileNotFoundError(
            f'Could not find the pressure-dependent network file {config.pes.network}.\n'
            f'A relative "pes.network" is resolved against the input file\'s own directory, '
            f'{input_directory}.')

    verbose = verbose_level(args)

    os.makedirs(project_directory, exist_ok=True)
    logger = Logger(project=os.path.basename(os.path.normpath(project_directory)),
                    project_directory=project_directory,
                    verbose=verbose,
                    t0=datetime.datetime.now())

    # check that ARC is available
    check_dependencies()

    result = run_pes_loop(config, project_directory, qm_runner=arc_qm_runner, logger=logger)

    logger.log(f'\n\nThe PES exploration loop terminated with status {result.status!r} '
               f'after {len(result.rounds)} round(s).', level='always')
    if result.reason:
        logger.log(f'Reason: {result.reason}', level='always')
    if result.final_network_path is not None:
        logger.log(f'Final network: {result.final_network_path}', level='always')
    if result.final_diagram_path is not None:
        logger.log(f'PES diagram: {result.final_diagram_path}', level='always')
    logger.log_footer()

    return result


if __name__ == '__main__':
    sys.exit(exit_code_for(main().status))
