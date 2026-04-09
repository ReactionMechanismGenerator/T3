#!/usr/bin/env python
# encoding: utf-8

name = "T3_test_dest_kinetics"
shortDesc = u"Small destination kinetics library for T3 tests"
longDesc = u"""
Subset of BurkeH2O2inN2 kinetics library, used for T3 unit tests.
Sourced from RMG-database/input/kinetics/libraries/BurkeH2O2inN2/reactions.py
"""
entry(
    index = 1,
    label = "H + O2 <=> O + OH",
    degeneracy = 1,
    kinetics = Arrhenius(A=(1.04e+14, 'cm^3/(mol*s)'), n=0, Ea=(15286, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""Hong et al., Proc. Comb. Inst. 33:309-316 (2011)""",
    longDesc = u"""""",
)

entry(
    index = 2,
    label = "OH + OH <=> O + H2O",
    degeneracy = 1,
    kinetics = Arrhenius(A=(33400, 'cm^3/(mol*s)'), n=2.42, Ea=(-1930, 'cal/mol'), T0=(1, 'K')),
    shortDesc = u"""Baulch et al., J. Phys. Chem. Ref. Data, 21:411 (1992)""",
    longDesc = u"""""",
)
