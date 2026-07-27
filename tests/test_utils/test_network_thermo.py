#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_network_thermo module
"""

import sys

import pytest

from t3.utils.network_thermo import (NetworkTextUnparseable, NetworkThermoCeiling,
                                     format_skipped_species, network_thermo_t_max)


# Synthetic fixtures for network_thermo_t_max: no in-repo network fixture (network1_1/network4_1/
# network4_2) has a species whose outer NASA Tmax is BELOW the pdep block's own Tmax (every real
# fixture's species thermo is valid to 5000-6000 K, well above any pressureDependence Tmax used),
# so the min-over-species / outer-vs-inner / skip behaviors this function exists for are
# necessarily exercised with small, hand-written text rather than a real network.
TWO_SPECIES_TEXT = '''
species(
    label = 'S1',
    structure = SMILES('[H]'),
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,0,0], Tmin=(100,'K'), Tmax=(2000,'K')),
            NASAPolynomial(coeffs=[2.5,0,0,0,0,0,0], Tmin=(2000,'K'), Tmax=(9000,'K')),
        ],
        Tmin = (100,'K'),
        Tmax = (3000,'K'),
    ),
)

species(
    label = 'S2',
    structure = SMILES('[H][H]'),
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,0,0], Tmin=(100,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (100,'K'),
        Tmax = (6000,'K'),
    ),
)
'''


def test_network_thermo_t_max_returns_the_minimum_over_species():
    """Test that network_thermo_t_max returns the MIN outer NASA Tmax across species (3000, from
    S1), not the max (6000, from S2) or the first species parsed (also S1's 3000, here, so this
    alone would not distinguish min from first-wins; test_..._is_not_the_first_species_parsed
    below pins that distinction separately). One species with a low ceiling makes the whole
    network's Arkane solve fail at any T above it, no matter how high another species' thermo is
    valid to, so returning anything other than the minimum would let a caller ask Arkane for a
    temperature it will refuse."""
    assert network_thermo_t_max(TWO_SPECIES_TEXT).t_max == 3000.0


def test_network_thermo_t_max_is_not_the_first_species_parsed():
    """Test that network_thermo_t_max is not simply the FIRST species' outer Tmax: with S2 (Tmax
    6000) listed before S1 (Tmax 3000), the result must still be the minimum (3000), not the first
    one encountered while walking the file (6000)."""
    reordered_text = '''
species(
    label = 'S2',
    structure = SMILES('[H][H]'),
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,0,0], Tmin=(100,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (100,'K'),
        Tmax = (6000,'K'),
    ),
)

species(
    label = 'S1',
    structure = SMILES('[H]'),
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,0,0], Tmin=(100,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (100,'K'),
        Tmax = (3000,'K'),
    ),
)
'''
    assert network_thermo_t_max(reordered_text).t_max == 3000.0


def test_network_thermo_t_max_reads_the_outer_nasa_tmax_not_the_inner_nasapolynomial_tmax():
    """Test that network_thermo_t_max reads the OUTER Tmax= on the NASA(...) call (3000 for S1),
    never a per-segment NASAPolynomial(..., Tmax=...) value nested inside polynomials=[...] (2000
    or 9000 for S1's two segments). The outer value is the overall validity ceiling Arkane
    actually checks a requested temperature against; reading an inner segment boundary instead
    would either under- or over-clamp relative to what Arkane itself enforces."""
    single_species_text = '''
species(
    label = 'S1',
    structure = SMILES('[H]'),
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,0,0], Tmin=(100,'K'), Tmax=(2000,'K')),
            NASAPolynomial(coeffs=[2.5,0,0,0,0,0,0], Tmin=(2000,'K'), Tmax=(9000,'K')),
        ],
        Tmin = (100,'K'),
        Tmax = (3000,'K'),
    ),
)
'''
    assert network_thermo_t_max(single_species_text).t_max == 3000.0


def test_network_thermo_t_max_skips_a_species_with_a_non_kelvin_tmax_unit():
    """Test that a species whose outer Tmax is not in Kelvin (e.g. degrees Rankine) contributes no
    constraint: comparing a differently-unitted bound as if it were Kelvin would silently
    mis-clamp, so it must be skipped rather than compared directly. With S1 skipped, only S2's
    6000 K remains. Also asserts (Finding B) that the skip is not silent: S1 must be named in
    ``.skipped`` so a caller can log that the returned ceiling may be higher than the network's
    true one."""
    text = TWO_SPECIES_TEXT.replace("Tmax = (3000,'K'),\n    ),\n)\n\nspecies(\n    label = 'S2'",
                                    "Tmax = (5400,'R'),\n    ),\n)\n\nspecies(\n    label = 'S2'")
    result = network_thermo_t_max(text)
    assert result.t_max == 6000.0
    assert len(result.skipped) == 1
    assert "'S1'" in result.skipped[0]


def test_network_thermo_t_max_skips_a_species_with_no_thermo_keyword():
    """Test that a species(...) call with no thermo= keyword at all contributes no constraint
    (rather than raising, or being treated as an unconstrained/zero ceiling). With S1's thermo
    keyword removed, only S2's 6000 K remains. Also asserts (Finding B) that S1 is named in
    ``.skipped`` rather than the skip being invisible."""
    text = '''
species(
    label = 'S1',
    structure = SMILES('[H]'),
)

species(
    label = 'S2',
    structure = SMILES('[H][H]'),
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,0,0], Tmin=(100,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (100,'K'),
        Tmax = (6000,'K'),
    ),
)
'''
    result = network_thermo_t_max(text)
    assert result.t_max == 6000.0
    assert len(result.skipped) == 1
    assert "'S1'" in result.skipped[0]


def test_network_thermo_t_max_skips_a_species_using_thermodata_instead_of_nasa():
    """Test that a species whose thermo= is a ThermoData(...) call (not NASA(...)) contributes no
    constraint: ThermoData carries no polynomial validity ceiling to compare, so it must be
    skipped rather than mistakenly read as if it were a NASA(...) call. With S1 using ThermoData,
    only S2's 6000 K remains. Also asserts (Finding B) that S1 is named in ``.skipped`` rather
    than the skip being invisible."""
    text = '''
species(
    label = 'S1',
    structure = SMILES('[H]'),
    thermo = ThermoData(Tdata=([300,400,500],'K'), Cpdata=([20.8,20.8,20.8],'J/(mol*K)'), H298=(0,'kJ/mol'), S298=(114.7,'J/(mol*K)')),
)

species(
    label = 'S2',
    structure = SMILES('[H][H]'),
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[2.5,0,0,0,0,0,0], Tmin=(100,'K'), Tmax=(6000,'K')),
        ],
        Tmin = (100,'K'),
        Tmax = (6000,'K'),
    ),
)
'''
    result = network_thermo_t_max(text)
    assert result.t_max == 6000.0
    assert len(result.skipped) == 1
    assert "'S1'" in result.skipped[0]


def test_network_thermo_t_max_returns_none_when_no_species_contributes_a_ceiling():
    """Test that network_thermo_t_max returns None (not zero, and not a raise) when every species
    in the file lacks a determinable outer NASA Tmax -- e.g. a file with no species(...) calls at
    all. A caller (write_arkane_network_input_file / _rewrite_method_and_sensitivity) must treat
    None as 'clamp nothing', so returning a spurious zero here would wrongly clamp every network to
    an unsolvable Tmax of 0 K. Also asserts (Finding B) that ``.skipped`` is empty here -- there
    are no species(...) calls at all, so there is nothing to name as skipped."""
    text = "reaction(label = 'reactionOnly', reactants = ['S1'], products = ['S2'])\n"
    result = network_thermo_t_max(text)
    assert result.t_max is None
    assert result.skipped == ()


def test_network_thermo_t_max_raises_network_text_unparseable_not_generic_valueerror():
    """Test that unparseable text (Finding C's groundwork) raises the distinguishable
    ``NetworkTextUnparseable`` subclass, not a bare ``ValueError``. Callers (e.g.
    ``t3.pdep.hybrid._rewrite_method_and_sensitivity``) rely on catching exactly this subclass so
    that ONLY the genuinely-unparseable-text case degrades to 'clamp nothing'; if this ever
    reverted to a plain ValueError, a narrowed `except NetworkTextUnparseable` at a call site
    would stop catching it and the pipeline's own self-check would no longer be the one reporting
    the corruption. It must also be a ValueError (so any pre-existing broad `except ValueError`
    still works)."""
    with pytest.raises(NetworkTextUnparseable):
        network_thermo_t_max('this is not valid python (')
    with pytest.raises(ValueError):
        network_thermo_t_max('this is not valid python (')


def test_format_skipped_species_caps_the_listed_entries():
    """Test that format_skipped_species (used by both writer.py and hybrid.py call sites to log a
    capped warning) lists at most `limit` entries verbatim and summarizes the rest as 'and N
    more', so a network with many skipped species does not flood the log."""
    skipped = tuple(f"'S{i}' (no thermo= keyword)" for i in range(7))
    formatted = format_skipped_species(skipped, limit=5)
    for i in range(5):
        assert f"'S{i}'" in formatted
    assert "'S5'" not in formatted
    assert "'S6'" not in formatted
    assert 'and 2 more' in formatted
    assert format_skipped_species(()) == ''


def test_importing_writer_does_not_pull_in_t3_pdep():
    """Pin the import layering that motivated moving network_thermo_t_max out of t3.pdep.parser
    and into this leaf module: t3.pdep depends on t3.utils.writer (several t3.pdep submodules do
    `from t3.utils.writer import ...`), so if t3.utils.writer ever again imports from t3.pdep
    (directly, or indirectly via t3.pdep.parser), importing t3.pdep.parser -- which eagerly runs
    t3/pdep/__init__.py, which in turn imports api/cache/hybrid/mesolver/parser/selector -- while
    t3.utils.writer is only partway through its own execution raises an ImportError ("cannot
    import name ... from partially initialized module").

    A prior version of this code "fixed" that cycle by moving `from t3.pdep.parser import
    network_thermo_t_max` to the very bottom of t3/utils/writer.py, after every other name in the
    module was already defined. That dodges the immediate ImportError, but it inverts the intended
    layering (a low-level util module depending on a higher-level package) and leaves a landmine:
    since most of t3.pdep imports FROM t3.utils.writer, any future tidy-up that moves that
    bottom-of-file import back up to the top-of-file import block (where every other import in the
    file lives, and where a reader would expect it) silently reintroduces the cycle and breaks the
    whole t3.pdep package. The real fix is this module (t3.utils.network_thermo): a leaf module
    that imports nothing from t3.* at all, so t3.utils.writer can import network_thermo_t_max from
    it, at the top of the file, with no cycle to dodge.

    This test asserts that a fresh import of t3.utils.writer never causes any t3.pdep submodule to
    be loaded, which is exactly what a regression back to importing from t3.pdep (in either the
    bottom-of-file or top-of-file style) would cause.
    """
    for module_name in list(sys.modules):
        if module_name == 't3.pdep' or module_name.startswith('t3.pdep.'):
            del sys.modules[module_name]
    if 't3.utils.writer' in sys.modules:
        del sys.modules['t3.utils.writer']

    import t3.utils.writer  # noqa: F401

    pdep_modules = [name for name in sys.modules if name == 't3.pdep' or name.startswith('t3.pdep.')]
    assert not pdep_modules, (
        f"Importing t3.utils.writer pulled in t3.pdep submodule(s) {pdep_modules}; this means "
        f"t3.utils.writer imports from t3.pdep again (directly or indirectly), reintroducing the "
        f"circular import that network_thermo_t_max's move to t3.utils.network_thermo was meant "
        f"to fix for good."
    )
