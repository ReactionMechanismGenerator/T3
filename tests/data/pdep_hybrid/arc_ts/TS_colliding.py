#!/usr/bin/env python
# -*- coding: utf-8 -*-

# A fixture mirroring ARC's REAL on-disk layout, in which every job's output file is renamed to
# 'output.out' and only the job directory distinguishes them. The log files referenced here are
# stubs for path/vendoring assertions only; they are not Arkane-parseable.

linear = False

spinMultiplicity = 2

energy = Log('logs/sp_a1234/output.out')

geometry = Log('logs/freq_a1235/output.out')

frequencies = Log('logs/freq_a1235/output.out')

rotors = [
    HinderedRotor(
        scanLog=Log('logs/scan_a1236/output.out'),
        pivots = [1, 2],
        top = [3, 4],
        symmetry = 1,
        fit='best',
    ),
]
