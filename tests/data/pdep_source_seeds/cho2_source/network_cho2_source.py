# Minimal SOURCE-ONLY CHO2 seed for the standalone PES exploration loop.
#
# This file carries ONLY the source well [O]C=O (formyloxy), the Ar bath gas, one network()
# block naming the well as the sole isomer, and one pressureDependence() block -- the minimum
# t3.pdep.explorer.input_file.write_arkane_explorer_input_file needs. There is NO reaction(),
# NO transitionState(), and NO adduct isomer: Arkane's explorer grows the CHO2 surface from this
# seed using the RMG database. Requiring reactions here would invert the data flow -- the explorer
# is what creates them (see t3/pdep/parser.py's require_reactions and t3/pdep/api.py).
#
# A UNIMOLECULAR source is deliberate: a bimolecular H + CO2 source makes Arkane's explorer
# discover TWO disjoint CHO2 networks, which ArkaneExplorerAdapter refuses (there is no single
# unambiguous result). A single well explores one connected network the loop can resolve and draw.
#
# The [O]C=O and Ar statmech blocks are lifted verbatim from a real RMG-explored CHO2 network
# (examples/pes_loop_cho2/t3_seed_run/network_cho2_clean.py in the I-006/I-007 worktree), so no
# energy is invented.

species(
    label = '[O]C=O',
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
    label = 'Ar',
    structure = adjacencyList("""1 Ar u0 p4 c0
"""),
    molecularWeight = (39.8775,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

network(
    label = 'PDepNetwork CHO2 source',
    isomers = [
        '[O]C=O',
    ],
    reactants = [
    ],
    bathGas = {
        'Ar': 1.0,
    },
)

pressureDependence(
    label = 'PDepNetwork CHO2 source',
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
