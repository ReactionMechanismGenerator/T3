#!/usr/bin/env python
# encoding: utf-8

name = "T3_test_src_thermo"
shortDesc = u"Small source thermo library for T3 merge tests"
longDesc = u"""
Subset of BurkeH2O2 thermo library, used for T3 unit tests.
Contains one structural duplicate of the destination library (H2 with a
different label) and two new species (HO2, H2O2).
Sourced from RMG-database/input/thermo/libraries/BurkeH2O2.py
"""
entry(
    index = 1,
    label = "Hydrogen",
    molecule =
"""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([6.9,6.96,7,7.02,7.07,7.21,7.73],'cal/(mol*K)'),
        H298 = (0,'kcal/mol'),
        S298 = (31.21,'cal/(mol*K)'),
    ),
    shortDesc = u"""H2 with a different label - should be detected as isomorphic duplicate""",
    longDesc = u"""""",
)

entry(
    index = 2,
    label = "HO2",
    molecule =
"""
multiplicity 2
1 O u0 p2 c0 {2,S} {3,S}
2 O u1 p2 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([8.35,8.89,9.47,10,10.77,11.38,12.48],'cal/(mol*K)'),
        H298 = (3,'kcal/mol'),
        S298 = (54.75,'cal/(mol*K)'),
    ),
    shortDesc = u"""HO2 radical""",
    longDesc = u"""""",
)

entry(
    index = 3,
    label = "H2O2",
    molecule =
"""
1 O u0 p2 c0 {2,S} {3,S}
2 O u0 p2 c0 {1,S} {4,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {2,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([10.42,11.45,12.35,13.11,14.29,15.21,16.85],'cal/(mol*K)'),
        H298 = (-32.53,'kcal/mol'),
        S298 = (55.65,'cal/(mol*K)'),
    ),
    shortDesc = u"""H2O2""",
    longDesc = u"""""",
)
