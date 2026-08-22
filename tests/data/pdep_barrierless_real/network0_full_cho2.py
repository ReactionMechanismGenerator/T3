species(
    label = 'O=[C]O(3)',
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
    label = '[O]C=O(4)',
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
    label = '[CH]1OO1(11)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u1 p0 c0 {1,S} {2,S} {4,S}
4 H u0 p0 c0 {3,S}
"""),
    E0 = (172.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,618.002,618.016,618.083,618.105],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5022,-0.00152154,5.13338e-05,-8.06025e-08,3.62137e-11,20739.7,7.02908], Tmin=(100,'K'), Tmax=(863.876,'K')), NASAPolynomial(coeffs=[13.0342,-0.011362,8.87123e-06,-1.87852e-09,1.32467e-13,17813.2,-44.9692], Tmin=(863.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-OsOsHH) + ring(dioxirane) + radical(OCJO)"""),
)

species(
    label = 'H(34)(1)',
    structure = adjacencyList("""multiplicity 2
1 H u1 p0 c0
"""),
    E0 = (211.8,'kJ/mol'),
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (1.00797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,25473.7,-0.446683], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(211.8,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: Klippenstein_Glarborg2016"""),
)

species(
    label = 'CO2(11)(2)',
    structure = adjacencyList("""1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
"""),
    E0 = (-402.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([516.743,1045.39,1045.39,2520.53],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(1622.99,'J/mol'), sigma=(3.941,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""PrimaryTransportLibrary"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28084,0.0025019,8.08178e-06,-1.20508e-08,4.66534e-12,-48400.8,6.00083], Tmin=(100,'K'), Tmax=(978.22,'K')), NASAPolynomial(coeffs=[4.67429,0.00260961,-9.85671e-07,1.95709e-10,-1.49833e-14,-48951.2,-2.11084], Tmin=(978.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-402.51,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""CO2""", comment="""Thermo library: thermo_DFT_CCSDTF12_BAC"""),
)

species(
    label = '[OH](5)',
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
    label = '[C-]#[O+](7)',
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.56838,-0.000852126,2.48918e-06,-1.56331e-09,3.13596e-13,-14284.3,3.57912], Tmin=(100,'K'), Tmax=(1571.64,'K')), NASAPolynomial(coeffs=[2.91306,0.00164658,-6.88618e-07,1.21038e-10,-7.84024e-15,-14180.9,6.71048], Tmin=(1571.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.741,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[O][C]=O(8)',
    structure = adjacencyList("""multiplicity 3
1 O u1 p2 c0 {3,S}
2 O u0 p2 c0 {3,D}
3 C u1 p0 c0 {1,S} {2,D}
"""),
    E0 = (33.3014,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.81051,-0.000257464,1.76457e-05,-2.38762e-08,9.15949e-12,4016.03,8.55809], Tmin=(100,'K'), Tmax=(975.948,'K')), NASAPolynomial(coeffs=[6.50399,-1.42575e-05,-6.91615e-08,7.04576e-11,-9.11448e-15,2952.97,-7.12367], Tmin=(975.948,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.3014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-OdOsH) + radical(OJC=O) + radical((O)CJOH)"""),
)

species(
    label = '[O](10)',
    structure = adjacencyList("""multiplicity 3
1 O u2 p2 c0
"""),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,-3.01681e-12,3.74582e-15,-1.50857e-18,1.86626e-22,29230.2,5.12616], Tmin=(100,'K'), Tmax=(4879.8,'K')), NASAPolynomial(coeffs=[4.28461,-0.00145495,4.44804e-07,-6.0436e-11,3.07922e-15,27479.1,-6.32199], Tmin=(4879.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=O(9)',
    structure = adjacencyList("""multiplicity 2
1 O u0 p2 c0 {2,D}
2 C u1 p0 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
"""),
    E0 = (33.3613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([880.033,2045.55,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.35603,-0.00347092,1.25665e-05,-9.99501e-09,2.27892e-12,3995.77,2.75111], Tmin=(100,'K'), Tmax=(1565.71,'K')), NASAPolynomial(coeffs=[4.6185,0.0050448,-4.39253e-06,9.73308e-10,-7.07456e-14,2787.59,-2.22863], Tmin=(1565.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.3613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-OdHH) + radical(HCdsJO)"""),
)

species(
    label = '[C]1OO1(13)',
    structure = adjacencyList("""multiplicity 3
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {3,S}
3 C u2 p0 c0 {1,S} {2,S}
"""),
    E0 = (444.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180,180,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0094,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.33628,-0.0119154,0.000106782,-1.61414e-07,7.09926e-11,53476.4,6.56729], Tmin=(100,'K'), Tmax=(880.815,'K')), NASAPolynomial(coeffs=[23.4559,-0.0321134,1.9978e-05,-3.98082e-09,2.7201e-13,47171.3,-103.617], Tmin=(880.815,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-OsOsHH) + ring(dioxirane) + radical(CH2_triplet)"""),
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
    E0 = (30.706 - 190.710, "kJ/mol"), # removing the applied energy_correction
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (112.916 - 190.710, "kJ/mol"), # removing the applied energy_correction
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (435.811 - 190.710, "kJ/mol"), # removing the applied energy_correction
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (134.747 - 190.710, "kJ/mol"), # removing the applied energy_correction
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (467.105 - 190.710, "kJ/mol"), # removing the applied energy_correction
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (374.912 - 190.710, "kJ/mol"), # removing the applied energy_correction
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (20.775 - 190.710, "kJ/mol"), # removing the applied energy_correction
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (435.811 - 190.710, "kJ/mol"), # removing the applied energy_correction
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (846.730 - 190.710, "kJ/mol"), # removing the applied energy_correction
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['H(34)(1)', 'CO2(11)(2)'],
    products = ['O=[C]O(3)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.37e+08,'cm^3/(mol*s)'), n=1.63, Ea=(30.7064,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1700,'K'), comment="""Estimated using template [Od_R;HJ] for rate rule [Od_Cdd-O2d;HJ]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[OH](5)', '[C-]#[O+](7)'],
    products = ['O=[C]O(3)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.41e+07,'cm^3/(mol*s)'), n=0, Ea=(12.552,'kJ/mol'), T0=(1,'K'), Tmin=(250,'K'), Tmax=(2500,'K'), comment="""Estimated using template [COm;O_rad] for rate rule [COm;O_pri_rad]
Euclidian distance = 1.0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(34)(1)', '[O][C]=O(8)'],
    products = ['O=[C]O(3)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.82057e+14,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O in family R_Recombination.
Multiplied by reaction path degeneracy 2.0"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C=O(4)'],
    products = ['O=[C]O(3)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.3e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using template [R2H_S;O_rad_out;XH_out] for rate rule [R2H_S;O_rad_out;CO_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O](10)', '[CH]=O(9)'],
    products = ['[O]C=O(4)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1667.73,'m^3/(mol*s)'), n=1.126, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [CO_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination
Ea raised from -8.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C=O(4)'],
    products = ['[CH]1OO1(11)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(6.95187e+18,'s^-1'), n=-1.62813, Ea=(354.137,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), uncertainty=RateUncertainty(mu=0.08275800860939324, var=25.73102631834909, Tref=1000.0, N=5, data_mean=0.0, correlation='Backbone0_N-2R!H-inRing_N-1R!H-inRing_Sp-2R!H-1R!H',), comment="""Estimated from node Backbone0_N-2R!H-inRing_N-1R!H-inRing_Sp-2R!H-1R!H in family Intra_R_Add_Endocyclic."""),
)

reaction(
    label = 'reaction1',
    reactants = ['H(34)(1)', 'CO2(11)(2)'],
    products = ['[O]C=O(4)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7699.21,'m^3/(mol*s)'), n=1.33054, Ea=(20.7755,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Od;HJ] for rate rule [CO2;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 17.5 to 20.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(34)(1)', '[O][C]=O(8)'],
    products = ['[O]C=O(4)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, data_mean=0.0, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""Estimated from node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O in family R_Recombination."""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(34)(1)', '[C]1OO1(13)'],
    products = ['[CH]1OO1(11)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = 'PDepNetwork #1',
    isomers = [
        'O=[C]O(3)',
        '[O]C=O(4)',
        '[CH]1OO1(11)',
    ],
    reactants = [
        ('H(34)(1)', 'CO2(11)(2)'),
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

