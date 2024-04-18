database(
    thermoLibraries=['primaryThermoLibrary', 'BurkeH2O2', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo', 'Spiekermann_refining_elementary_reactions', 'CurranPentane', 'CBS_QB3_1dHR', 'primaryNS'],
    reactionLibraries=['primaryH2O2', 'NOx2018', '2006_Joshi_OH_CO', 'CurranPentane', 'Klippenstein_Glarborg2016', 'C2H4+O_Klipp2017', 'FFCM1(-)', 'C2H2_init', 'Narayanaswamy', 'Mebel_C6H5_C2H2', 'C10H11', 'C12H11_pdep', 'Lai_Hexylbenzene', '1989_Stewart_2CH3_to_C2H5_H', '2001_Tokmakov_H_Toluene_to_CH3_Benzene', '2003_Miller_Propargyl_Recomb_High_P', '2005_Senosiain_OH_C2H2', 'kislovB', 'c-C5H5_CH3_Sharma', 'fascella', '2006_Joshi_OH_CO', '2009_Sharma_C5H5_CH3_highP', '2015_Buras_C2H3_C4H6_highP', 'C3', 'Methylformate', 'C6H5_C4H4_Mebel', 'vinylCPD_H', 'Mebel_Naphthyl', 'Fulvene_H'],
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
    label='N2',
    reactive=False,
    structure=SMILES('N#N'),
)

species(
    label='O2',
    reactive=True,
    structure=SMILES('[O][O]'),
)

species(
    label='OH',
    reactive=True,
    structure=SMILES('[OH]'),
)

species(
    label='HO2',
    reactive=True,
    structure=SMILES('O[O]'),
)

species(
    label='H',
    reactive=True,
    structure=SMILES('[H]'),
)

species(
    label='H2',
    reactive=True,
    structure=SMILES('[H][H]'),
)

species(
    label='CH3',
    reactive=True,
    structure=SMILES('[CH3]'),
)

species(
    label='CO',
    reactive=True,
    structure=SMILES('[C-]#[O+]'),
)

species(
    label='CO2',
    reactive=True,
    structure=SMILES('O=C=O'),
)

species(
    label='HOCO',
    reactive=True,
    structure=SMILES('O[C]=O'),
)

species(
    label='OCHO',
    reactive=True,
    structure=SMILES('[O]C=O'),
)

species(
    label='CH2O',
    reactive=True,
    structure=SMILES('C=O'),
)

species(
    label='HCO',
    reactive=True,
    structure=SMILES('[CH]=O'),
)

simpleReactor(
    temperature=(500.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'O2': 1.0,
        'FA': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(500.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'FA': 2.5,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(500.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        'FA': 4.0,
        'N2': 3.76,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(500.0, 'K'),
    pressure=(10.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'O2': 1.0,
        'FA': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(500.0, 'K'),
    pressure=(10.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'FA': 2.5,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(500.0, 'K'),
    pressure=(10.0, 'bar'),
    initialMoleFractions={
        'FA': 4.0,
        'N2': 3.76,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(500.0, 'K'),
    pressure=(100.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'O2': 1.0,
        'FA': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(500.0, 'K'),
    pressure=(100.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'FA': 2.5,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(500.0, 'K'),
    pressure=(100.0, 'bar'),
    initialMoleFractions={
        'FA': 4.0,
        'N2': 3.76,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(1250.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'O2': 1.0,
        'FA': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(1250.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'FA': 2.5,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(1250.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        'FA': 4.0,
        'N2': 3.76,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(1250.0, 'K'),
    pressure=(10.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'O2': 1.0,
        'FA': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(1250.0, 'K'),
    pressure=(10.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'FA': 2.5,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(1250.0, 'K'),
    pressure=(10.0, 'bar'),
    initialMoleFractions={
        'FA': 4.0,
        'N2': 3.76,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(1250.0, 'K'),
    pressure=(100.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'O2': 1.0,
        'FA': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(1250.0, 'K'),
    pressure=(100.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'FA': 2.5,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(1250.0, 'K'),
    pressure=(100.0, 'bar'),
    initialMoleFractions={
        'FA': 4.0,
        'N2': 3.76,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(2000.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'O2': 1.0,
        'FA': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(2000.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'FA': 2.5,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(2000.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        'FA': 4.0,
        'N2': 3.76,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(2000.0, 'K'),
    pressure=(10.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'O2': 1.0,
        'FA': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(2000.0, 'K'),
    pressure=(10.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'FA': 2.5,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(2000.0, 'K'),
    pressure=(10.0, 'bar'),
    initialMoleFractions={
        'FA': 4.0,
        'N2': 3.76,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(2000.0, 'K'),
    pressure=(100.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'O2': 1.0,
        'FA': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(2000.0, 'K'),
    pressure=(100.0, 'bar'),
    initialMoleFractions={
        'N2': 3.76,
        'FA': 2.5,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)


simpleReactor(
    temperature=(2000.0, 'K'),
    pressure=(100.0, 'bar'),
    initialMoleFractions={
        'FA': 4.0,
        'N2': 3.76,
        'O2': 1.0,
    },
    terminationTime=(60.0, 's'),
    nSims=12,
)

model(
    toleranceMoveToCore=0.2,
    toleranceInterruptSimulation=0.2,
    filterReactions=True,
    filterThreshold=100000000.0,
    maxNumObjsPerIter=1,
    terminateAtMaxObjects=False,
)

simulator(atol=1e-16, rtol=1e-08, sens_atol=1e-06, sens_rtol=0.0001)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(2.0, 'kJ/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300, 2500, 'K', 10),
    pressures=(0.1, 110, 'bar', 10),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=16,
)

options(
    name='Seed',
    generateSeedEachIteration=True,
    saveSeedToDatabase=False,
    units='si',
    generateOutputHTML=True,
    generatePlots=False,
    saveSimulationProfiles=True,
    verboseComments=False,
    saveEdgeSpecies=False,
    keepIrreversible=False,
    trimolecularProductReversible=True,
    wallTime='00:00:00:00',
    saveSeedModulus=-1,
)

generatedSpeciesConstraints(
    allowed=['input species', 'seed mechanisms', 'reaction libraries'],
    maximumCarbonAtoms=2,
    maximumOxygenAtoms=4,
    maximumNitrogenAtoms=0,
    maximumSiliconAtoms=0,
    maximumSulfurAtoms=0,
    maximumHeavyAtoms=5,
    maximumRadicalElectrons=1,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2=True,
)
