#!/usr/bin/env python3
# encoding: utf-8

"""
A script for running RMG incore using the rmg_env.

Usage::

    conda activate rmg_env

    running with no optional args:
        python-jl rmg_incore_script.py full/path/to/rmg/input.py > >(tee -a out.txt) 2> >(tee -a err.txt >&2)

    or with optional args:
        python-jl rmg_incore_script.py full/path/to/rmg/input.py -v 20 -m 100 > >(tee -a out.txt) 2> >(tee -a err.txt >&2)
"""

import argparse
import os
import sys
import traceback

from pydas.daspk import DASPKError

from rmgpy.exceptions import (ChemicallySignificantEigenvaluesError,
                              ChemkinError,
                              CollisionError,
                              CoreError,
                              ILPSolutionError,
                              InputError,
                              InvalidMicrocanonicalRateError,
                              KineticsError,
                              ModifiedStrongCollisionError,
                              NetworkError,
                              PressureDependenceError,
                              ReactionError,
                              ReservoirStateError,
                              StatmechError,
                              StatmechFitError,
                              )
from rmgpy.rmg.main import RMG
from rmgpy.rmg.main import initialize_log as initialize_rmg_log


def parse_command_line_arguments(command_line_args=None):
    """
    Parse command-line arguments.

    Args:
        command_line_args: The command line arguments.

    Returns:
        The parsed command-line arguments by key words.
    """

    parser = argparse.ArgumentParser(description='The Tandem Tool (T3) RMG incore runner script')
    parser.add_argument('file', metavar='FILE', type=str, nargs=1, help='The RMG input file')
    parser.add_argument('-v', '--verbose', type=int, nargs=1, default=20, help='the logging level')
    parser.add_argument('-m', '--maxiter', type=int, nargs=1, default=0, help='max rmg iterations')
    args = parser.parse_args(command_line_args)
    args.file = args.file[0]
    return args


def main() -> None:
    """
    Run RMG incore.
    """
    args = parse_command_line_arguments()
    input_file = args.file
    project_directory = os.path.abspath(os.path.dirname(args.file))

    initialize_rmg_log(
        verbose=args.verbose[0],
        log_file_name=os.path.join(project_directory, 'RMG.log'),
    )
    rmg = RMG(input_file=input_file, output_directory=project_directory)
    rmg_kwargs = dict()
    if args.maxiter:
        rmg_kwargs['max_iterations'] = args.maxiter

    try:
        rmg.execute(initialize=True, **rmg_kwargs)
    except (ChemicallySignificantEigenvaluesError,
            ChemkinError,
            CollisionError,
            CoreError,
            DASPKError,
            ILPSolutionError,
            InvalidMicrocanonicalRateError,
            KineticsError,
            ModifiedStrongCollisionError,
            NetworkError,
            PressureDependenceError,
            ReactionError,
            ReservoirStateError,
            StatmechError,
            StatmechFitError,
            ) as e:
        sys.stderr.write(f'\n\n********\n'
                         f'RMG threw an exception and did not converge.\n'  # Keep this text unchanged.
                         f'Exception type: {e.__class__}\n'
                         f'Exception message:\n{e}\n'
                         f'********\n')
        print(f'RMG Errored with {e.__class__}. Got the following trace:')
        print(traceback.format_exc())

    except InputError:
        # This exception should not be raised.
        print('Error: Something seems to be wrong with the RMG input file, please check your input.')
        raise


if __name__ == '__main__':
    main()
