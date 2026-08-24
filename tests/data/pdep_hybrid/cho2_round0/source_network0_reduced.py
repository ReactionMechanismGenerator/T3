species(
    label = 'O=[C]O(8)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
4 H u0 p0 c0 {1,S}
"""),
    E0 = (-184.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,607.278,1696.14],'cm^-1')),
        HinderedRotor(inertia=(0.000457043,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.61,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.65821,0.00206077,2.03288e-05,-2.91283e-08,1.14311e-11,-22171.5,8.26011], Tmin=(100,'K'), Tmax=(964.009,'K')), NASAPolynomial(coeffs=[7.09694,0.0016118,-4.75818e-07,1.29874e-10,-1.30513e-14,-23476.6,-11.5343], Tmin=(964.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-OdOsH) + radical((O)CJOH)"""),
)

species(
    label = '[O]C=O(1)',
    structure = adjacencyList("""multiplicity 2
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (-169.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.61,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9649,-0.00518829,3.73952e-05,-4.21919e-08,1.47224e-11,-20431.8,7.33591], Tmin=(100,'K'), Tmax=(994.515,'K')), NASAPolynomial(coeffs=[5.26528,0.00541264,-2.4716e-06,5.38764e-10,-4.28526e-14,-21473.4,-2.86661], Tmin=(994.515,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-169.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-OdOsH) + radical(OJC=O)"""),
)

species(
    label = '[H](6)',
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
    label = 'O=C=O(5)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
"""),
    E0 = (-403.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([459.922,1087.69,1087.69,2296.71],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27862,0.00274142,7.16109e-06,-1.08032e-08,4.14302e-12,-48470.3,5.97934], Tmin=(100,'K'), Tmax=(988.879,'K')), NASAPolynomial(coeffs=[4.54606,0.00291919,-1.15487e-06,2.27661e-10,-1.70916e-14,-48980.3,-1.43257], Tmin=(988.879,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-403.087,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + missing(O2d-Cdd) + group(Cdd-OdOd)"""),
)

species(
    label = 'Ar',
    structure = adjacencyList("""1 Ar u0 p4 c0
"""),
    molecularWeight = (39.8775,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (21.347 - 191.282, "kJ/mol"), # removing the applied energy_correction
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (135.319 - 191.282, "kJ/mol"), # removing the applied energy_correction
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (30.706 - 191.282, "kJ/mol"), # removing the applied energy_correction
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction3',
    reactants = ['[H](6)', 'O=C=O(5)'],
    products = ['[O]C=O(1)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(7699.21,'m^3/(mol*s)'), n=1.33054, Ea=(21.3474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Od;HJ] for rate rule [CO2;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 17.5 to 21.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C=O(1)'],
    products = ['O=[C]O(8)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.3e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using template [R2H_S;O_rad_out;XH_out] for rate rule [R2H_S;O_rad_out;CO_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[H](6)', 'O=C=O(5)'],
    products = ['O=[C]O(8)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.37e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""Estimated using template [Od_R;HJ] for rate rule [Od_Cdd-O2d;HJ]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

network(
    label = 'PDepNetwork #1',
    isomers = [
        'O=[C]O(8)',
        '[O]C=O(1)',
    ],
    reactants = [
    ],
    bathGas = {
        'Ar': 1,
    },
)

pressureDependence(
    label = 'PDepNetwork #1',
    Tmin = (700,'K'),
    Tmax = (3200,'K'),
    Tcount = 10,
    Tlist = ([3200,2290.91,1784.07,1460.87,1236.81,1072.34,946.479,847.059,766.54,700],'K'),
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

