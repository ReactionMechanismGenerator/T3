"""
t3 utils dependencies module
"""


def check_dependencies():
    """
    Check that ARC and its native deps (py_rdl) are available.
    """
    missing = []

    try:
        import py_rdl  # noqa: F401 — import IS the availability check
    except ImportError:
        missing.append(
            ('py_rdl',
             'Required by ARC. Install with `make install-pyrdl` from the T3 repo, '
             'or `pip install git+https://github.com/rareylab/RingDecomposerLib.git#subdirectory=src/Python`.')
        )

    try:
        from arc import ARC  # noqa: F401 — import IS the availability check
    except ImportError:
        missing.append(
            ('ARC',
             'See https://reactionmechanismgenerator.github.io/ARC/installation.html')
        )

    if missing:
        msg = 'Error: Cannot execute T3. Missing the following critical component(s):\n'
        for name, hint in missing:
            msg += f'  - {name}\n    ({hint})\n'
        msg += '\nInstall the missing packages, make sure they were added to the PYTHONPATH, and rerun T3.'
        print(msg)
        raise ValueError('T3 is missing core dependencies, see the above message.')
