"""
t3 utils dependencies module
"""


def check_dependencies():
    """
    Check that RMG and ARC are available.
    """
    arc_available = True
    try:
        from arc import ARC
    except ImportError:
        arc_available = False

    if not arc_available:
        msg = 'Error: Cannot execute T3. Missing the following critical component(s):\n'
        msg += '  - ARC\n'
        msg += '    (See https://reactionmechanismgenerator.github.io/ARC/installation.html)\n'
        msg += '\nInstall the missing packages, make sure they were added to the PYTHONPATH, and rerun T3.'
        print(msg)
        raise ValueError(f'T3 is missing core dependencies, see the above message.')
