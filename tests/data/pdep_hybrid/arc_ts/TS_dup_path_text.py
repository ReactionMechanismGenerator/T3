#!/usr/bin/env python
# -*- coding: utf-8 -*-

# NOTE (test fixture): item 5 regression coverage. This mirrors TS1.py's shape, but adds an
# unrelated variable whose string value happens to be textually IDENTICAL to the 'energy' Log(...)
# path below. A text-based (regex) rewrite of that literal path would incorrectly rewrite this
# unrelated string too, since it cannot distinguish "this exact quoted substring, wherever it
# appears" from "the argument of this specific Log(...) call". A rewrite that splices by the
# Log(...) argument's own AST node position must leave this line untouched.

linear = False

spinMultiplicity = 2

energy = Log('logs/ts1_sp.out')

geometry = Log('logs/ts1_freq.out')

frequencies = Log('logs/ts1_freq.out')

rotors = [
    HinderedRotor(scanLog=Log('logs/ts1_scan.out'), pivots=[1, 2], top=[3, 4, 5], symmetry=1, fit='best'),
]

unrelated_note = 'logs/ts1_sp.out'
