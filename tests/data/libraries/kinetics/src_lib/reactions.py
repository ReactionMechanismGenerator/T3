#!/usr/bin/env python
# encoding: utf-8

name = "T3_test_src_kinetics"
shortDesc = u"Small source kinetics library for T3 merge tests"
longDesc = u"""
Subset of BurkeH2O2inN2 kinetics library, used for T3 unit tests.
Contains one duplicate-label reaction (should be skipped) and two new
reactions (should be merged into the destination, including new species
in the dictionary).
Sourced from RMG-database/input/kinetics/libraries/BurkeH2O2inN2/reactions.py
"""
entry(
    index = 1,
    label = "H + O2 <=> O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(9.99e+13, 'cm^3/(mol*s)'), n=0, Ea=(15000, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""Same label as the destination - should be skipped on merge.""",
    longDesc = u"""""",
)

entry(
    index = 2,
    label = "H2 + OH <=> H2O + H",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.16e+08, 'cm^3/(mol*s)'), n=1.51, Ea=(3430, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""Michael and Sutherland, J. Phys. Chem. 92:3853 (1988)""",
    longDesc = u"""""",
)

entry(
    index = 3,
    label = "H + HO2 <=> H2 + O2",
    degeneracy = 1,
    kinetics = Arrhenius(A=(2.81e+13, 'cm^3/(mol*s)'), n=0, Ea=(1068, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""Mueller et al., Int. J. Chem. Kinet., 31:113 (1999)""",
    longDesc = u"""""",
)
