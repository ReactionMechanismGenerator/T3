database(
    thermoLibraries=['primaryThermoLibrary'],
    reactionLibraries=[],
    transportLibraries=['OneDMinN2', 'PrimaryTransportLibrary', 'NOx2018', 'GRI-Mech'],
    seedMechanisms=[],
    kineticsDepositories='default',
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)

species(
    label='H2',
    reactive=True,
    structure=SMILES('[H][H]'),
)

species(
    label='O2',
    reactive=True,
    structure=SMILES('[O][O]'),
)

species(
    label='H',
    reactive=True,
    structure=SMILES('[H]'),
)

species(
    label='OH',
    reactive=True,
    structure=SMILES('[OH]'),
)

simpleReactor(
    temperature=(1000.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        'H2': 0.67,
        'O2': 0.33,
        'H': 0,
        'OH': 0,
    },
    terminationConversion={'H2': 0.9},
    terminationTime=(1.0, 's'),
    nSims=12,
)

model(
    toleranceMoveToCore=0.001,
    toleranceInterruptSimulation=0.001,
    filterReactions=False,
    filterThreshold=100000000.0,
    maxNumObjsPerIter=1,
    terminateAtMaxObjects=False,
)

simulator(atol=1e-16, rtol=1e-08)
