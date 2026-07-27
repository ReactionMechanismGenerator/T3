# A hand-written pdep network fixture that DECLARES its product channels explicitly.
#
# Every RMG-generated network file omits the ``products=`` keyword and lets Arkane derive the
# product channels from the path reactions (arkane/input.py, the ``if products is None:`` branch).
# This fixture exercises the OTHER branch, which must use the declared list verbatim.
#
# It lives under ``tests/data/pdep_network_variants/`` rather than ``tests/data/pdep_network/``
# because the latter is a ``t3.paths['PDep SA']`` target that other tests ``shutil.rmtree`` during
# teardown; a fixture placed there is silently deleted by an unrelated test, making results depend
# on test-execution order.

species(
    label = 'A',
    E0 = (0.0, 'kJ/mol'),
)

species(
    label = 'B',
    E0 = (10.0, 'kJ/mol'),
)

species(
    label = 'C',
    E0 = (20.0, 'kJ/mol'),
)

species(
    label = 'D',
    E0 = (30.0, 'kJ/mol'),
)

transitionState(
    label = 'TS1',
    E0 = (100.0, 'kJ/mol'),
)

reaction(
    label = 'reaction1',
    reactants = ['A'],
    products = ['B', 'C'],
    transitionState = 'TS1',
    kinetics = Arrhenius(
        A = (1e10, 's^-1'),
        n = 0.0,
        Ea = (100.0, 'kJ/mol'),
        T0 = (1, 'K'),
        comment = 'Estimated using an average for rate rule',
    ),
)

network(
    label = 'explicit',
    isomers = ['A'],
    reactants = [],
    products = [['B', 'C'], ['D']],
    bathGas = {'A': 1.0},
)
