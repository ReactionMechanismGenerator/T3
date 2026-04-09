#!/usr/bin/env python
# encoding: utf-8

name = "T3_test_dest_thermo"
shortDesc = u"Small destination thermo library for T3 tests"
longDesc = u"""
Subset of BurkeH2O2 thermo library, used for T3 unit tests.
Sourced from RMG-database/input/thermo/libraries/BurkeH2O2.py
"""
entry(
    index = 0,
    label = "H",
    molecule =
"""
multiplicity 2
1 H u1 p0 c0
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([4.97,4.97,4.97,4.97,4.97,4.97,4.97],'cal/(mol*K)'),
        H298 = (52.1,'kcal/mol'),
        S298 = (27.39,'cal/(mol*K)'),
    ),
    shortDesc = u"""H atom""",
    longDesc = u"""""",
)

entry(
    index = 1,
    label = "H2",
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
    shortDesc = u"""H2""",
    longDesc = u"""""",
)

entry(
    index = 2,
    label = "OH",
    molecule =
"""
multiplicity 2
1 O u1 p2 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
    thermo = ThermoData(
        Tdata = ([300,400,500,600,800,1000,1500],'K'),
        Cpdata = ([7.16,7.08,7.05,7.06,7.15,7.34,7.87],'cal/(mol*K)'),
        H298 = (8.9,'kcal/mol'),
        S298 = (43.9,'cal/(mol*K)'),
    ),
    shortDesc = u"""OH radical""",
    longDesc = u"""""",
)
