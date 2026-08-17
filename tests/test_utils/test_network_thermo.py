#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_network_thermo module
"""

import sys

import pytest

from t3.utils.network_thermo import (NetworkTextUnparseable, TGridClampRecord,
                                     format_skipped_species, network_thermo_t_max,
                                     t_grid_clamp_shape_error)


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


# --------------------------------------------------------------------------------------------
# TGridClampRecord's own field contract.
#
# This record is the SOURCE of the t_grid_clamp provenance that rides along with an SA cache
# sidecar and then with a persisted PDepNetworkSelection, and until now it type-checked nothing
# at all: every downstream contract described the DICT the record renders, while the record that
# renders it accepted anything. The same gap SensitiveTransitionState had.
# --------------------------------------------------------------------------------------------

def _valid_clamped_record() -> TGridClampRecord:
    """A record shaped exactly like the one t3.utils.writer builds when it clamps a T grid."""
    return TGridClampRecord(clamped=True, requested_t_max=2500.0, thermo_ceiling=2000.0,
                            written_t_max=2000.0, tlist_dropped=True, tlist_original_highest=2400.0,
                            skipped_species=('S3 (no thermo= keyword)',))


def _valid_unclamped_record() -> TGridClampRecord:
    """A record shaped exactly like the one t3.utils.writer builds when no clamp was needed."""
    return TGridClampRecord(clamped=False, requested_t_max=2500.0, thermo_ceiling=3000.0,
                            written_t_max=2500.0)


@pytest.mark.parametrize('record_factory', [_valid_clamped_record, _valid_unclamped_record])
def test_the_records_the_real_writers_build_still_construct_and_render(record_factory):
    """The two shapes t3.utils.writer and t3.pdep.hybrid actually produce must survive the field
    contract untouched. This is the guard against the contract being written to describe an
    idealized record rather than the real one -- if either of these ever fails, the check is wrong,
    not the writer."""
    record = record_factory()
    rendered = record.as_dict()

    assert isinstance(rendered['clamped'], bool)
    assert isinstance(rendered['skipped_species'], list)


@pytest.mark.parametrize('bad_clamped', ['yes', 1, 0, None, 'False'])
def test_a_clamp_record_refuses_a_clamped_that_is_not_an_explicit_bool(bad_clamped):
    """``clamped`` carries the entire three-state design: True/False is an EXPLICIT verdict from a
    writer that ran the clamp logic, and 'unknown' is the whole record being absent. A truthy string
    or a 1/0 int is neither -- it is a verdict nobody stated, and since nothing downstream reads this
    key back out of the dict, a wrong value would never surface as anything but a human reading a
    provenance record that quietly lies to them."""
    with pytest.raises(ValueError, match='clamped'):
        TGridClampRecord(clamped=bad_clamped)


@pytest.mark.parametrize('bad_tlist_dropped', ['yes', 1, None])
def test_a_clamp_record_refuses_a_tlist_dropped_that_is_not_an_explicit_bool(bad_tlist_dropped):
    """Same reasoning as ``clamped``: whether an explicit Tlist line was dropped from the written
    file is a statement of fact about a solve, and a non-bool is not that statement."""
    with pytest.raises(ValueError, match='tlist_dropped'):
        TGridClampRecord(clamped=True, tlist_dropped=bad_tlist_dropped)


@pytest.mark.parametrize('field_name', ['requested_t_max', 'thermo_ceiling', 'written_t_max',
                                        'tlist_original_highest'])
@pytest.mark.parametrize('bad_temperature', ['2500', True, [2500.0], {'K': 2500.0}])
def test_a_clamp_record_refuses_a_temperature_that_is_not_a_number(field_name, bad_temperature):
    """Every one of these four fields is a temperature in Kelvin. A string that looks like one is
    the dangerous case: it renders to YAML perfectly well, so no write-time guard would ever catch
    it, and it reads back as a temperature to a human while comparing against nothing. ``True`` is
    excluded separately because ``isinstance(True, int)`` is True in Python."""
    with pytest.raises(ValueError, match=field_name):
        TGridClampRecord(clamped=True, **{field_name: bad_temperature})


@pytest.mark.parametrize('field_name', ['requested_t_max', 'thermo_ceiling', 'written_t_max',
                                        'tlist_original_highest'])
def test_a_clamp_record_accepts_an_int_temperature(field_name):
    """A Tmax written as ``(3000, 'K')`` in a network file is an int, and network_thermo_t_max only
    happens to coerce it with ``float(...)``. Refusing an int here would be a contract describing
    this module's current internals rather than the quantity, and would break the first caller that
    passes a literal 3000. Over-refusal is the other way this contract can be wrong."""
    record = TGridClampRecord(clamped=True, **{field_name: 3000})
    assert record.as_dict()[field_name] == 3000


@pytest.mark.parametrize('field_name', ['requested_t_max', 'thermo_ceiling', 'written_t_max',
                                        'tlist_original_highest'])
def test_a_clamp_record_accepts_none_for_every_optional_temperature(field_name):
    """``None`` means 'not applicable / could not be read' on all four, and every one of them is
    genuinely None in some real writer path (e.g. a file with no readable Tmax line at all)."""
    assert TGridClampRecord(clamped=True, **{field_name: None}).as_dict()[field_name] is None


def test_a_clamp_record_refuses_a_bare_string_for_skipped_species():
    """``skipped_species='CH4'`` is the failure worth naming on its own: it is iterable, so
    ``as_dict()``'s ``list(...)`` silently renders it as ``['C', 'H', '4']`` -- three species that
    do not exist, in a field whose entire purpose is to warn that the computed ceiling may be looser
    than the network's true one. A corrupted warning is worse than no warning."""
    with pytest.raises(ValueError, match='skipped_species'):
        TGridClampRecord(clamped=True, skipped_species='CH4')


@pytest.mark.parametrize('bad_entry', [42, None, ['nested'], object()])
def test_a_clamp_record_refuses_a_non_string_entry_in_skipped_species(bad_entry):
    """Each entry is a human-readable sentence naming a species and why it was skipped. A non-string
    entry either renders as a Python-tagged YAML object (unreadable by every T3 loader) or reads as
    a species name that is not one."""
    with pytest.raises(ValueError, match='skipped_species'):
        TGridClampRecord(clamped=True, skipped_species=('S1 (no thermo= keyword)', bad_entry))


def test_a_clamp_record_accepts_a_list_for_skipped_species():
    """A list is accepted where the annotation says tuple, deliberately. ``as_dict()`` REBUILDS this
    container with ``list(...)`` rather than deep-copying it, so a list and a tuple render to exactly
    the same bytes -- refusing one would buy nothing but record-identity purity. This is the opposite
    call from PDepNetworkSelection's list fields, which ARE refused as tuples, and for the opposite
    reason: that record's as_dict() deep-copies, so a tuple there survives into the rendering as
    ``!!python/tuple``."""
    record = TGridClampRecord(clamped=True, skipped_species=['S1 (no thermo= keyword)'])
    assert record.as_dict()['skipped_species'] == ['S1 (no thermo= keyword)']


# --------------------------------------------------------------------------------------------
# t_grid_clamp_shape_error: the same contract, asked of a DICT rather than a record.
#
# Two callers need opposite outcomes from this one question -- t3.pdep.cache.read_t_grid_clamp_record
# must collapse a malformed sidecar dict to None (unknown provenance, never a refusal), while
# PDepNetworkSelection.validate() must refuse one -- so the question returns a reason instead of
# raising, and each caller decides what to do with it.
# --------------------------------------------------------------------------------------------

@pytest.mark.parametrize('record_factory', [_valid_clamped_record, _valid_unclamped_record])
def test_the_rendering_of_a_valid_record_always_passes_the_shape_check(record_factory):
    """The check is validated against what ``as_dict()`` ACTUALLY emits, not against a hand-written
    dict that merely looks like it. Checking a shape you wrote yourself is how a guard ends up
    describing bytes nobody produces -- the same mistake as validating with a different dumper than
    the one the write uses."""
    assert t_grid_clamp_shape_error(record_factory().as_dict()) is None


@pytest.mark.parametrize('not_a_dict', ['not-a-dict', 42, True, ['clamped'], None])
def test_the_shape_check_refuses_anything_that_is_not_a_dict(not_a_dict):
    """Including ``None``: absence is the caller's business to distinguish (it means unknown
    provenance), and this function answers only 'is this a rendering of the record'."""
    assert t_grid_clamp_shape_error(not_a_dict) is not None


def test_the_shape_check_refuses_a_dict_that_never_says_whether_a_clamp_happened():
    """A dict with no ``clamped`` key is a FOURTH state the three-state design says must not exist:
    not 'clamped', not 'not clamped', and not the absent-record 'unknown' either -- a provenance
    record that declines to answer the one question it exists to answer."""
    reason = t_grid_clamp_shape_error({'requested_t_max': 2500.0, 'written_t_max': 2000.0})
    assert reason is not None
    assert 'clamped' in reason


@pytest.mark.parametrize('key, bad_value', [('clamped', 'yes'),
                                            ('clamped', 1),
                                            ('tlist_dropped', 'no'),
                                            ('requested_t_max', '2500'),
                                            ('thermo_ceiling', True),
                                            ('written_t_max', ['2000']),
                                            ('tlist_original_highest', 'high'),
                                            ('skipped_species', 'CH4'),
                                            ('skipped_species', [42]),
                                            ('skipped_species', ('a tuple',)),
                                            ])
def test_the_shape_check_names_the_key_that_is_wrong(key, bad_value):
    """A reason that does not name the offending key sends whoever reads it back to the file to
    diff seven fields by eye. ``skipped_species`` as a tuple is refused here (unlike on the record,
    where as_dict() normalizes it) because a tuple reaching a dict has already skipped that
    normalization and would be written as ``!!python/tuple``, which no T3 loader can read back."""
    rendered = _valid_clamped_record().as_dict()
    rendered[key] = bad_value

    reason = t_grid_clamp_shape_error(rendered)
    assert reason is not None
    assert key in reason


def test_the_shape_check_tolerates_a_key_an_older_writer_never_wrote():
    """Only ``clamped`` is required. Every other key is optional, because a sidecar written by an
    older T3 -- before a field was added to the record -- is missing that key and is still perfectly
    honest provenance about the clamp it does describe."""
    assert t_grid_clamp_shape_error({'clamped': False}) is None


def test_the_shape_check_tolerates_a_key_a_newer_writer_added():
    """The mirror case, and the more dangerous one to get wrong: sidecars are files that outlive the
    version that wrote them, so a NEWER T3 adding an eighth field must not make every older T3
    silently discard the provenance in that sidecar. Unknown keys are carried, not refused --
    whether they can be persisted at all is already the write-time guard's question, not this
    one's."""
    rendered = _valid_clamped_record().as_dict()
    rendered['clamp_reason'] = 'a field a future version added'

    assert t_grid_clamp_shape_error(rendered) is None


# --- Round 61: three ways a "well-typed" record was still not one -------------------------------

@pytest.mark.parametrize('field_name', ['requested_t_max', 'thermo_ceiling', 'written_t_max',
                                        'tlist_original_highest'])
@pytest.mark.parametrize('not_finite', [float('nan'), float('inf'), float('-inf')])
def test_a_clamp_record_refuses_a_temperature_that_is_not_finite(field_name, not_finite):
    """A nan or an inf is a number and not a temperature. nan is the one that breaks this arc's
    central property outright rather than merely misinforming: it renders to YAML, reloads from YAML,
    and then compares UNEQUAL to itself -- so a record carrying one can never be shown to have
    round-tripped, and every equality assertion about it silently reports failure for a reason that
    has nothing to do with what was saved."""
    with pytest.raises(ValueError, match=field_name):
        TGridClampRecord(clamped=True, **{field_name: not_finite})


@pytest.mark.parametrize('not_finite', [float('nan'), float('inf')])
def test_the_shape_check_refuses_a_temperature_that_is_not_finite(not_finite):
    """The dict half of the same contract. A sidecar can carry `.nan` as plain YAML, so this is not
    only reachable by direct construction."""
    rendered = _valid_clamped_record().as_dict()
    rendered['written_t_max'] = not_finite

    assert t_grid_clamp_shape_error(rendered) is not None


def test_a_list_passed_for_skipped_species_is_normalized_to_a_tuple():
    """A frozen record that keeps a caller's list is not frozen in the way that matters. Without
    normalization the caller retains a live reference, appends to it after construction, and
    `as_dict()` renders entries the contract never saw -- valid at construction, invalid at the only
    moment anybody cares about. Accepting the convenient input and normalizing it keeps both."""
    entries = ['S1 (no thermo= keyword)']
    record = TGridClampRecord(clamped=True, skipped_species=entries)

    entries.append(42)

    assert isinstance(record.skipped_species, tuple)
    assert record.as_dict()['skipped_species'] == ['S1 (no thermo= keyword)']
    assert t_grid_clamp_shape_error(record.as_dict()) is None


def test_the_shape_contract_does_not_claim_to_stop_a_record_from_lying():
    """Pin the BOUNDARY that the docstring states, so nobody later reads this contract as a semantic
    guarantee. `clamped=False` alongside a written Tmax that differs from the requested one is
    incoherent provenance, and it passes: the real writers read both temperatures from the same
    parsed line and so cannot produce it, which means a cross-field rule could only ever fire on a
    record from another version -- and on the sidecar path a failed check silently DROPS provenance,
    turning a disagreement with a foreign writer into the exact loss this record exists to prevent."""
    incoherent = {'clamped': False, 'requested_t_max': 2500.0, 'written_t_max': 2000.0,
                  'tlist_dropped': True, 'tlist_original_highest': None}

    assert t_grid_clamp_shape_error(incoherent) is None
