#!/usr/bin/env python3
# encoding: utf-8

"""
T3 executable module
"""

import argparse
import os

from arc.common import read_yaml_file

from t3 import T3
from t3.utils.dependencies import check_dependencies


def parse_command_line_arguments(command_line_args=None):
    """
    Parse command-line arguments.

    Args:
        command_line_args: The command line arguments.

    Returns:
        The parsed command-line arguments by key words.
    """

    parser = argparse.ArgumentParser(description='The Tandem Tool (T3)')
    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
                        help='a file describing the job to execute')

    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-d', '--debug', action='store_true', help='print debug information')
    group.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')

    args = parser.parse_args(command_line_args)
    args.file = args.file[0]

    return args


def main() -> None:
    """
    The main T3 executable function.
    """

    args = parse_command_line_arguments()
    input_file = args.file
    project_directory = os.path.abspath(os.path.dirname(args.file))
    input_dict = read_yaml_file(path=input_file)
    if 'project' not in list(input_dict.keys()):
        raise ValueError('A project name must be provided.')

    verbose = 20
    if args.debug:
        verbose = 10
    elif args.quiet:
        verbose = 30
    input_dict['verbose'] = input_dict['verbose'] if 'verbose' in input_dict else verbose
    if 'project_directory' not in input_dict or not input_dict['project_directory']:
        input_dict['project_directory'] = project_directory

    # check that RMG and ARC are available
    check_dependencies()

    t3_object = T3(**input_dict)
    t3_object.execute()


if __name__ == '__main__':
    main()
