species(
    label = 'C1rad(2)',
    structure = SMILES('[CH3]'),
    E0 = (136.525,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([542.511,1417.18,1419.16,2485.5,3603.52,3603.63],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96043,0.00059294,8.78575e-06,-9.88031e-09,3.63235e-12,16421.9,0.339866], Tmin=(100,'K'), Tmax=(697.711,'K')), NASAPolynomial(coeffs=[3.09511,0.0055543,-1.88159e-06,3.13335e-10,-2.05195e-14,16542.6,4.20298], Tmin=(697.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: thermo_DFT_CCSDTF12_BAC"""),
)

species(
    label = 'H(34)',
    structure = SMILES('[H]'),
    E0 = (211.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(211.8,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CH2(T)(48)',
    structure = SMILES('[CH2]'),
    E0 = (382.715,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1068.15,2790.8,3622.46],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01171,-0.00015992,3.27674e-06,-2.41782e-09,5.74073e-13,46029.4,0.537381], Tmin=(100,'K'), Tmax=(1101.44,'K')), NASAPolynomial(coeffs=[3.1529,0.00295891,-9.70593e-07,1.52925e-10,-9.4151e-15,46218.6,4.76363], Tmin=(1101.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(382.715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: thermo_DFT_CCSDTF12_BAC"""),
)

species(
    label = 'He',
    structure = SMILES('[He]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(84.8076,'J/mol'), sigma=(2.576,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""HE""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'Ne',
    structure = SMILES('[Ne]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (594.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['H(34)', 'CH2(T)(48)'],
    products = ['C1rad(2)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.13e+13,'cm^3/(mol*s)'), n=0.32, Ea=(0,'cal/mol'), T0=(1,'K'), comment="""Kinetics taken from the arrheniusHigh attribute of a Troe/Lindemann exprssion. Originally from reaction library FFCM1(-)"""),
)

network(
    label = 'PDepNetwork #1',
    isomers = [
        'C1rad(2)',
    ],
    reactants = [
    ],
    bathGas = {
        'He': 0.5,
        'Ne': 0.5,
    },
)

pressureDependence(
    label = 'PDepNetwork #1',
    Tmin = (300,'K'),
    Tmax = (2100,'K'),
    Tcount = 8,
    Tlist = ([302.491,323.355,370.585,457.988,614.983,900.017,1394.8,1985.54],'K'),
    Pmin = (0.1,'bar'),
    Pmax = (100,'bar'),
    Pcount = 5,
    Plist = ([0.118417,0.415262,3.16228,24.0812,84.4471],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

