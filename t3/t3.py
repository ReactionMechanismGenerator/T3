"""
T3 executable module
"""

import argparse

from .tandem import execute


def parse_command_line_arguments(command_line_args=None):
    """
    Parse the command-line arguments being passed to T3.

    Returns:
        The parsed command-line arguments by key words.
    """

    parser = argparse.ArgumentParser(description='The RMG-ARC Tandem Tool (T3)')

    parser.add_argument('-r', '--rmg', metavar='RMG', type=str, nargs=1,
                        help="The legacy RMG input file")

    parser.add_argument('-a', '--arc', metavar='ARC', type=str, nargs=1,
                        help="The augmented ARC input file")

    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-d', '--debug', action='store_true', help='print debug information')
    group.add_argument('-q', '--quiet', action='store_true', help='only print warnings and errors')

    args = parser.parse_args(command_line_args)

    return args


def main() -> None:
    """
    The main T3 executable function.
    """
    # Parse the command-line arguments (requires the argparse module)
    args = parse_command_line_arguments()

    kwargs = {
        'rmg': args.rmg,
        'arc': args.arc,
    }

    # Check that RMG and ARC are available
    rmg_available, arc_available = True, True
    try:
        from rmgpy.rmg.main import RMG
    except ImportError:
        rmg_available = False
    try:
        from arc import ARC
    except ImportError:
        arc_available = False

    if not rmg_available or not arc_available:
        msg = 'Error: Cannot execute T3. Missing the following critical component(s):\n'
        if not rmg_available:
            msg += '  - RMG\n'
            msg += '    (See http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/installation/index.html)\n'
        if not rmg_available:
            msg += '  - ARC\n'
            msg += '    (See https://reactionmechanismgenerator.github.io/ARC/installation.html)\n'
        msg += '\nInstall the missing packages and rerun T3.'
        return None

    execute(args, kwargs)


if __name__ == '__main__':
    main()
