species(
    label = 'HCOH(189)',
    structure = adjacencyList("""1 O u0 p1 c+1 {2,D} {4,S}
2 C u0 p1 c-1 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (98.8073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1005.83,1006.45,1710.52,1900.49,1901.93,1902.12],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (30.026,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09305,-0.004728,2.65298e-05,-2.62736e-08,8.32851e-12,11882.7,4.17582], Tmin=(100,'K'), Tmax=(1030.21,'K')), NASAPolynomial(coeffs=[2.9525,0.0077518,-3.36376e-06,6.57069e-10,-4.76172e-14,11690.5,7.63865], Tmin=(1030.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.8073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""HCOH(S)""", comment="""Thermo library: thermo_DFT_CCSDTF12_BAC"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.03191,-0.0019319,8.88979e-06,2.99012e-09,-7.9075e-12,-14627.3,3.4758], Tmin=(10,'K'), Tmax=(601.531,'K')), NASAPolynomial(coeffs=[1.39956,0.00949141,-4.43189e-06,9.4843e-10,-7.43389e-14,-14200.6,15.7519], Tmin=(601.531,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(-121.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""P388""", comment="""Thermo library: 2FFOH_thermo"""),
)

species(
    label = 'H2(20)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43536,0.000212707,-2.78618e-07,3.40261e-10,-7.76018e-14,-1031.36,-3.90842], Tmin=(100,'K'), Tmax=(1959.09,'K')), NASAPolynomial(coeffs=[2.78813,0.000587692,1.58987e-07,-5.52691e-11,4.34276e-15,-596.121,0.112965], Tmin=(1959.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.60349,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""H2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'CO(18)',
    structure = adjacencyList("""1 O u0 p1 c+1 {2,T}
2 C u0 p1 c-1 {1,T}
"""),
    E0 = (-118.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2193.04],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.01,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852135,2.4892e-06,-1.56333e-09,3.13602e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.62,'K')), NASAPolynomial(coeffs=[2.91303,0.00164662,-6.88638e-07,1.21042e-10,-7.84056e-15,-14180.9,6.71064], Tmin=(1571.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = 'CHO(100)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (33.4995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1110.15,1948.69,2643.94],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01437,-0.000845987,4.50981e-06,3.99869e-10,-2.79968e-12,4028.33,4.15361], Tmin=(10,'K'), Tmax=(630.566,'K')), NASAPolynomial(coeffs=[2.68857,0.00476829,-2.19444e-06,4.56083e-10,-3.40533e-14,4251.12,10.3792], Tmin=(630.566,'K'), Tmax=(3000,'K'))], Tmin=(10,'K'), Tmax=(3000,'K'), E0=(33.4989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""P387""", comment="""Thermo library: 2FFOH_thermo"""),
)

species(
    label = 'Ar',
    structure = adjacencyList("""1 Ar u0 p4 c0
"""),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.8775,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (195.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (267.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (141.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (206.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['HCOH(189)'],
    products = ['FORMALDEHYDE(16)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5.2e+13,'s^-1'), n=0, Ea=(32109,'cal/mol'), T0=(1,'K'), comment="""Reaction library: 'NOx2018'"""),
)

reaction(
    label = 'reaction2',
    reactants = ['HCOH(189)'],
    products = ['H2(20)', 'CO(18)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.6e+13,'s^-1'), n=0, Ea=(49465,'cal/mol'), T0=(1,'K'), comment="""Reaction library: 'NOx2018'"""),
)

reaction(
    label = 'reaction3',
    reactants = ['FORMALDEHYDE(16)'],
    products = ['H2(20)', 'CO(18)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.7e+13,'s^-1'), n=0, Ea=(71976,'cal/mol'), T0=(1,'K'), comment="""Kinetics taken from the arrheniusHigh attribute of a Troe/Lindemann exprssion. Originally from reaction library FFCM1(-)"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(40)', 'CHO(100)'],
    products = ['FORMALDEHYDE(16)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.913e+14,'cm^3/(mol*s)'), n=-0.033, Ea=(-142,'cal/mol'), T0=(1,'K'), comment="""Kinetics taken from the arrheniusHigh attribute of a Troe/Lindemann exprssion. Originally from reaction library FFCM1(-)"""),
)

network(
    label = 'PDepNetwork #1',
    isomers = [
        'HCOH(189)',
        'FORMALDEHYDE(16)',
    ],
    reactants = [
        ('H2(20)', 'CO(18)'),
        ('H(40)', 'CHO(100)'),
    ],
    bathGas = {
        'Ar': 1,
    },
)

pressureDependence(
    label = 'PDepNetwork #1',
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

