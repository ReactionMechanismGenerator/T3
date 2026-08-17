# A hand-written pdep network fixture that DISCRIMINATES between Arkane's two arity branches when
# deriving product channels (arkane/input.py, the ``if products is None:`` branch).
#
# Arkane's own comment there says configurations "not already defined as reactants or isomers"
# become product channels, but each branch tests only ONE of the two:
#   * a UNIMOLECULAR side is tested against the isomers only  -- so a side that is a declared
#     reactant channel still becomes a product channel too ('D' below);
#   * a MULTI-SPECIES side is tested against the reactant channels only -- so a side whose first
#     sorted label happens to be an isomer label is still a product channel (['B', 'C'] below).
# Unifying the branches, which reads like a harmless cleanup, silently drops both channels and
# under-counts the net reactions Arkane writes. This fixture exists to fail when that happens.
#
# It lives under ``tests/data/pdep_network_variants/`` rather than ``tests/data/pdep_network/``
# because the latter is a ``t3.paths['PDep SA']`` target that other tests ``shutil.rmtree``.

species(
    label = 'A',
    E0 = (0.0, 'kJ/mol'),
)

species(
    label = 'B',
    E0 = (5.0, 'kJ/mol'),
)

species(
    label = 'C',
    E0 = (10.0, 'kJ/mol'),
)

species(
    label = 'D',
    E0 = (15.0, 'kJ/mol'),
)

transitionState(
    label = 'TS1',
    E0 = (100.0, 'kJ/mol'),
)

transitionState(
    label = 'TS2',
    E0 = (120.0, 'kJ/mol'),
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

reaction(
    label = 'reaction2',
    reactants = ['D'],
    products = ['A'],
    transitionState = 'TS2',
    kinetics = Arrhenius(
        A = (1e11, 's^-1'),
        n = 0.0,
        Ea = (120.0, 'kJ/mol'),
        T0 = (1, 'K'),
        comment = 'Estimated using an average for rate rule',
    ),
)

network(
    label = 'arity',
    isomers = ['A', 'B'],
    reactants = [['D']],
    bathGas = {'A': 1.0},
)
