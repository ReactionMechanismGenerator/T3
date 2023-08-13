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
    label='FA',
    reactive=True,
    structure=SMILES('OC=O'),
)

species(
    label='O2',
    reactive=True,
    structure=SMILES('[O][O]'),
)

species(
    label='N2',
    reactive=False,
    structure=SMILES('N#N'),
)

simpleReactor(
    temperature=(1000.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'FA': [1.0, 4.0],
        'O2': 1.0,
    },
    terminationRateRatio=0.1,
    nSims=12,
)

model(
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=0.01,
    filterReactions=True,
    filterThreshold=100000000.0,
    maxNumObjsPerIter=1,
    terminateAtMaxObjects=False,
)

simulator(atol=1e-16, rtol=1e-08, sens_atol=1e-06, sens_rtol=0.0001)

options(
    name='Seed',
    generateSeedEachIteration=True,
    saveSeedToDatabase=False,
    units='si',
    generateOutputHTML=False,
    generatePlots=False,
    saveSimulationProfiles=False,
    verboseComments=False,
    saveEdgeSpecies=True,
    keepIrreversible=False,
    trimolecularProductReversible=True,
    wallTime='00:00:05:00',
    saveSeedModulus=-1,
)
