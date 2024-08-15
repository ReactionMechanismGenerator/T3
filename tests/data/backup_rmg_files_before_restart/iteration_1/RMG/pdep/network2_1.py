species(
    label = 'HOCH2O(193)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {6,S}
2 O u1 p2 c0 {3,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {1,S}
"""),
    E0 = (-181.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,244.256],'cm^-1')),
        HinderedRotor(inertia=(0.00282506,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (47.0334,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.3898,0.00588265,2.76836e-05,-4.16177e-08,1.6744e-11,-21849.4,10.4839], Tmin=(100,'K'), Tmax=(944.176,'K')), NASAPolynomial(coeffs=[8.146,0.00487903,-1.13894e-06,2.10345e-10,-1.79653e-14,-23601,-16.7081], Tmin=(944.176,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-181.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""HOCH2O""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'OH(25)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
    E0 = (28.3945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3668.68],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (17.0074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51457,2.92734e-05,-5.32151e-07,1.01948e-09,-3.85939e-13,3414.25,2.10435], Tmin=(100,'K'), Tmax=(1145.76,'K')), NASAPolynomial(coeffs=[3.07194,0.00060402,-1.39806e-08,-2.13441e-11,2.48061e-15,3579.39,4.57801], Tmin=(1145.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.3945,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""OH(D)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'FORMALDEHYDE(16)',
    structure = adjacencyList("""1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
"""),
    E0 = (-121.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.61,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03191,-0.0019319,8.88979e-06,2.99012e-09,-7.9075e-12,-14627.3,3.4758], Tmin=(10,'K'), Tmax=(601.531,'K')), NASAPolynomial(coeffs=[1.39956,0.00949141,-4.43189e-06,9.4843e-10,-7.43389e-14,-14200.6,15.7519], Tmin=(601.531,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-121.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""P388""", comment="""Thermo library: 2FFOH_thermo"""),
)

species(
    label = 'H(40)',
    structure = adjacencyList("""multiplicity 2
1 H u1 p0 c0
"""),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-3.01681e-12,3.74582e-15,-1.50857e-18,1.86626e-22,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(4879.8,'K')), NASAPolynomial(coeffs=[4.28461,-0.00145495,4.44804e-07,-6.0436e-11,3.07922e-15,23723.1,-11.8931], Tmin=(4879.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'HOCHO(194)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {1,S}
"""),
    E0 = (-389.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.76241,-0.000943854,3.38779e-05,-4.04983e-08,1.43987e-11,-46821.2,7.63276], Tmin=(100,'K'), Tmax=(993.731,'K')), NASAPolynomial(coeffs=[5.645,0.00755182,-3.20851e-06,6.58919e-10,-5.05055e-14,-47989,-5.43104], Tmin=(993.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-389.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formic_acid""", comment="""Thermo library: thermo_DFT_CCSDTF12_BAC"""),
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
    E0 = (15.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (-28.5848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['OH(25)', 'FORMALDEHYDE(16)'],
    products = ['HOCH2O(193)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(6.3e+06,'cm^3/(mol*s)'), n=1.63, Ea=(4282,'cal/mol'), T0=(1,'K'), comment="""Reaction library: 'NOx2018'"""),
)

reaction(
    label = 'reaction2',
    reactants = ['HOCH2O(193)'],
    products = ['H(40)', 'HOCHO(194)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+14,'s^-1'), n=0, Ea=(14900,'cal/mol'), T0=(1,'K'), comment="""Reaction library: 'NOx2018'"""),
)

network(
    label = 'PDepNetwork #2',
    isomers = [
        'HOCH2O(193)',
    ],
    reactants = [
        ('OH(25)', 'FORMALDEHYDE(16)'),
    ],
    bathGas = {
        'Ar': 1,
    },
)

pressureDependence(
    label = 'PDepNetwork #2',
    Tmin = (300,'K'),
    Tmax = (2000,'K'),
    Tcount = 10,
    Tlist = ([301.578,314.572,342.653,390.652,467.665,589.953,785.229,1092.98,1528.1,1932.59],'K'),
    Pmin = (0.005,'bar'),
    Pmax = (100,'bar'),
    Pcount = 10,
    Plist = ([0.0053143,0.00857753,0.0213227,0.0746744,0.325889,1.53426,6.69574,23.4492,58.2919,94.0857],'bar'),
    maximumGrainSize = (2,'kJ/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

