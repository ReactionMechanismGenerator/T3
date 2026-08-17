#!/usr/bin/env python
# -*- coding: utf-8 -*-

# NOTE (test fixture): deliberately references a Log(...) file that does not exist on disk, to
# exercise write_hybrid_network_input_file's "missing referenced log" ValueError path. See the
# comment at the top of TS1.py for the general placeholder-stub caveat.

linear = False

spinMultiplicity = 2

energy = Log('logs/does_not_exist.out')

geometry = Log('logs/does_not_exist.out')

frequencies = Log('logs/does_not_exist.out')
