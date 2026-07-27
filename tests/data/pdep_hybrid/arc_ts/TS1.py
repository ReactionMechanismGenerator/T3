#!/usr/bin/env python
# -*- coding: utf-8 -*-

# NOTE (test fixture): this mirrors the shape ARC's `arc/statmech/arkane.py::species_input_template`
# writes for a transition state (linear / spinMultiplicity / energy / geometry / frequencies /
# rotors, each pointing at a Log(...) job file), but the Log(...) targets below are placeholder
# stub text files, not real, Arkane-parseable quantum chemistry output. They exist only so this
# fixture's vendoring/rewriting logic has real files to copy and check for existence.

linear = False

spinMultiplicity = 2

energy = Log('logs/ts1_sp.out')

geometry = Log('logs/ts1_freq.out')

frequencies = Log('logs/ts1_freq.out')

rotors = [
    HinderedRotor(scanLog=Log('logs/ts1_scan.out'), pivots=[1, 2], top=[3, 4, 5], symmetry=1, fit='best'),
]
