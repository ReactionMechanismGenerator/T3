species(
    label = 'N2H4(357)',
    structure = adjacencyList("""1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 N u0 p1 c0 {1,S} {5,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (86.2423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180,1118.57,1118.57,1118.57,1118.57,1118.57,1118.57,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00528502,'amu*angstrom^2'), symmetry=1, barrier=(1.99137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.0456,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.70157,0.00237956,2.3518e-05,-2.95123e-08,1.08413e-11,10387.1,6.00039], Tmin=(100,'K'), Tmax=(957.493,'K')), NASAPolynomial(coeffs=[4.76806,0.00858732,-2.91176e-06,5.20713e-10,-3.70544e-14,9694.05,-1.65094], Tmin=(957.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.2423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""NH2NH2""", comment="""Thermo library: thermo_DFT_CCSDTF12_BAC"""),
)

species(
    label = 'NH3NH(403)',
    structure = adjacencyList("""1 N u0 p0 c+1 {2,S} {3,S} {4,S} {5,S}
2 N u0 p2 c-1 {1,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
    E0 = (269.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (32.0456,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40007,0.00605245,7.56422e-06,-9.63729e-09,3.18472e-12,32468.1,7.40385], Tmin=(298,'K'), Tmax=(835.206,'K')), NASAPolynomial(coeffs=[1.78544,0.0137859,-6.32582e-06,1.45068e-09,-1.34464e-13,32737.8,14.9026], Tmin=(835.206,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), E0=(269.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""NH3NH""", comment="""Thermo library: NH3"""),
)

species(
    label = 'H2(16)',
    structure = adjacencyList("""1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (-8.60349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3765.59],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (2.01594,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43536,0.000212711,-2.78627e-07,3.40269e-10,-7.76035e-14,-1031.36,-3.90842], Tmin=(100,'K'), Tmax=(1959.07,'K')), NASAPolynomial(coeffs=[2.78818,0.000587629,1.59016e-07,-5.5275e-11,4.34319e-15,-596.15,0.112677], Tmin=(1959.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.60349,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""H2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'H2NN(S)(380)',
    structure = adjacencyList("""1 N u0 p0 c+1 {2,D} {3,S} {4,S}
2 N u0 p2 c-1 {1,D}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (289.87,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([914.9,1432.15,1752.99,1753.03,2280.21,3636.69],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.0297,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3222.33,'J/mol'), sigma=(3.467,'angstroms'), dipoleMoment=(3.507,'De'), polarizability=(2.349,'angstroms^3'), rotrelaxcollnum=1.0, comment="""OneDMinN2"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.1047,-0.00370494,1.98462e-05,-1.84452e-08,5.56765e-12,34860.7,3.20243], Tmin=(100,'K'), Tmax=(1042.92,'K')), NASAPolynomial(coeffs=[2.43687,0.00745361,-3.05164e-06,5.69834e-10,-3.99393e-14,34949.6,10.0773], Tmin=(1042.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.87,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""NNH2(S)""", comment="""Thermo library: thermo_DFT_CCSDTF12_BAC"""),
)

species(
    label = 'N2H2(362)',
    structure = adjacencyList("""1 N u0 p1 c0 {2,D} {3,S}
2 N u0 p1 c0 {1,D} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (189.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1080.77,1342.51,2026.84,2026.84,2026.84,3930.43],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.0297,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2685.64,'J/mol'), sigma=(3.531,'angstroms'), dipoleMoment=(0,'De'), polarizability=(2.297,'angstroms^3'), rotrelaxcollnum=1.0, comment="""OneDMinN2"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.14011,-0.00411614,1.94636e-05,-1.728e-08,5.01298e-12,22777.4,3.13603], Tmin=(100,'K'), Tmax=(1069.33,'K')), NASAPolynomial(coeffs=[2.1084,0.00776413,-3.20566e-06,5.96337e-10,-4.15014e-14,22967.1,11.9298], Tmin=(1069.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""HNNH""", comment="""Thermo library: thermo_DFT_CCSDTF12_BAC"""),
)

species(
    label = 'Ar',
    structure = adjacencyList("""1 Ar u0 p4 c0
"""),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.8775,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (236.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (281.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (187.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['N2H4(357)'],
    products = ['NH3NH(403)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.34e+11,'s^-1'), n=0.86, Ea=(64.5,'kcal/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Reaction library: 'primaryNitrogenLibrary'"""),
)

reaction(
    label = 'reaction2',
    reactants = ['N2H4(357)'],
    products = ['H2(16)', 'H2NN(S)(380)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.38e+09,'s^-1'), n=1.255, Ea=(75.3,'kcal/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Reaction library: 'primaryNitrogenLibrary'"""),
)

reaction(
    label = 'reaction3',
    reactants = ['N2H4(357)'],
    products = ['H2(16)', 'N2H2(362)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.7e+12,'s^-1'), n=0, Ea=(52.9,'kcal/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Reaction library: 'primaryNitrogenLibrary'"""),
)

network(
    label = 'PDepNetwork #21',
    isomers = [
        'N2H4(357)',
    ],
    reactants = [
    ],
    bathGas = {
        'Ar': 1,
    },
)

pressureDependence(
    label = 'PDepNetwork #21',
    Tmin = (700,'K'),
    Tmax = (2500,'K'),
    Tcount = 10,
    Tlist = ([2500,1944.44,1590.91,1346.15,1166.67,1029.41,921.053,833.333,760.87,700],'K'),
    Pmin = (0.1,'bar'),
    Pmax = (100,'bar'),
    Pcount = 10,
    Plist = ([0.1,0.215443,0.464159,1,2.15443,4.64159,10,21.5443,46.4159,100],'bar'),
    maximumGrainSize = (2,'kJ/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('pdeparrhenius',),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

