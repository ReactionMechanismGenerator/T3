#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_provenance module

Tests the energy-provenance classifier (``t3.pdep.provenance``): the last-writer-wins reading of
RMG/Arkane ``comment`` text that lets T3's PES diagram colour a level by where its energy came
from. The load-bearing case is the CONCATENATION TRAP -- a QM'd level's comment still carries the
superseded RMG estimate as a prefix, and a first-match parse would mislabel it (see the module
docstring). Comments here are taken verbatim from
``run/round_1/explorer/pdep/final/network0_reduced.py`` of the r002_m1-cho2-rerun run.
"""

import pytest

from t3.pdep.provenance import (
    PROVENANCE_FAMILY,
    PROVENANCE_GAV,
    PROVENANCE_LIBRARY,
    PROVENANCE_QM,
    PROVENANCE_UNKNOWN,
    classify_provenance,
    combine_channel_provenance,
)

# Real comments from the r002_m1-cho2-rerun final network file.
GAV_THEN_THERMOJOBS = ('Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + '
                       'group(Cds-OdOsH) + radical((O)CJOH)Thermo library: thermojobs')  # O=[C]O(8)
LIBRARY_THEN_THERMOJOBS = 'Thermo library: primaryThermoLibraryThermo library: thermojobs'  # [H](6)
BARE_PRIMARY_LIBRARY = 'Thermo library: primaryThermoLibrary'  # Ar
BARE_GAV = ('Thermo group additivity estimation: missing(O2d-Cdd) + missing(O2d-Cdd) + '
            'group(Cdd-OdOd)')
KINETICSJOBS = ("Fitted to 50 data points; dA = *|/ 12.9343, dn = +|- 0.336317, "
                "dEa = +|- 1.83023 kJ/molReaction library: 'kineticsjobs'")  # reaction1 / TS1
FAMILY_TEMPLATE = ('Estimated using template [Cdd_Od;HJ] for rate rule [CO2;HJ]\n'
                   'Euclidian distance = 1.0\nfamily: R_Addition_MultipleBond')  # reaction4 / TS2
FAMILY_NODE = ('Estimated from node Backbone0_1 in family Intra_R_Add_Endocyclic.')


@pytest.mark.parametrize('comment, expected', [
    # The trap: a QM'd level carries its superseded RMG estimate first, the QM library last.
    (GAV_THEN_THERMOJOBS, PROVENANCE_QM),       # GAV prefix + thermojobs -> QM (not GAV)
    (LIBRARY_THEN_THERMOJOBS, PROVENANCE_QM),   # curated library prefix + thermojobs -> QM
    (KINETICSJOBS, PROVENANCE_QM),              # kinetics fit + kineticsjobs -> QM
    # Genuinely external / estimated -- no QM library appended.
    (BARE_PRIMARY_LIBRARY, PROVENANCE_LIBRARY),
    (FAMILY_TEMPLATE, PROVENANCE_FAMILY),       # 'Estimated using template ...'
    (FAMILY_NODE, PROVENANCE_FAMILY),           # 'Estimated from node ... in family ...'
    (BARE_GAV, PROVENANCE_GAV),
    # Honest gap -- never defaulted into a class.
    ('', PROVENANCE_UNKNOWN),
    (None, PROVENANCE_UNKNOWN),
    ('   \n  ', PROVENANCE_UNKNOWN),
    ('some comment with no recognised provenance marker', PROVENANCE_UNKNOWN),
])
def test_classify_provenance(comment, expected):
    assert classify_provenance(comment) == expected


def test_reaction_library_quotes_stripped():
    """A double-quoted QM reaction library reads as QM just like the single-quoted form."""
    assert classify_provenance('...kJ/molReaction library: "kineticsjobs"') == PROVENANCE_QM


def test_curated_reaction_library_is_library_not_qm():
    """A reaction library that is NOT this run's QM output is a curated library, not QM."""
    assert classify_provenance("Matched.Reaction library: 'BurkeH2O2'") == PROVENANCE_LIBRARY


def test_last_writer_wins_over_three_concatenated_markers():
    """GAV, then a curated library, then the QM library -- the QM library (last) wins."""
    comment = ('Thermo group additivity estimation: group(x)'
               'Thermo library: primaryThermoLibrary'
               'Thermo library: thermojobs')
    assert classify_provenance(comment) == PROVENANCE_QM


def test_qm_library_names_are_configurable():
    """A caller can name a different QM output library; the default set no longer counts."""
    comment = '...kJ/molReaction library: myLocalQMlib'
    assert classify_provenance(comment) == PROVENANCE_LIBRARY
    assert classify_provenance(comment,
                               qm_library_names=frozenset({'myLocalQMlib'})) == PROVENANCE_QM


@pytest.mark.parametrize('fragments, expected', [
    ((PROVENANCE_QM, PROVENANCE_QM), PROVENANCE_QM),           # [H]+O=C=O, both QM
    ((PROVENANCE_QM, PROVENANCE_GAV), PROVENANCE_GAV),         # worst link wins
    ((PROVENANCE_QM, PROVENANCE_LIBRARY), PROVENANCE_LIBRARY),
    ((PROVENANCE_LIBRARY, PROVENANCE_FAMILY), PROVENANCE_FAMILY),
    ((PROVENANCE_QM, PROVENANCE_UNKNOWN), PROVENANCE_UNKNOWN),  # any unknown fragment -> unknown
    ((PROVENANCE_QM,), PROVENANCE_QM),                         # unimolecular
    ((), PROVENANCE_UNKNOWN),
])
def test_combine_channel_provenance(fragments, expected):
    assert combine_channel_provenance(fragments) == expected
