#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_parser
"""

import ast
import os

import pytest

from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.parser import (
    _call_keywords,
    PDepArkaneReaction,
    PDepNetwork,
    PDepPathReaction,
    parse_arkane_pdep_output_file,
    parse_arkane_pdep_output_text,
    parse_pdep_network_file,
    parse_pdep_network_text,
)

PDEP_NETWORK_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'iteration_1', 'RMG', 'pdep')

# Fixtures for the Arkane pdep *output* parser live under ``tests/data/pdep_me/`` rather than
# ``tests/data/pdep_network/`` because the latter is a ``t3.paths['PDep SA']`` target that other
# tests ``shutil.rmtree`` during teardown; see the module docstring in ``test_me_success.py`` for
# the full explanation.
PDEP_ME_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_me')

# (file stem, number of species, number of transition states, number of path reactions)
NETWORK_COUNTS = [
    ('network1_1', 5, 1, 1),
    ('network4_1', 11, 5, 5),
    ('network4_2', 16, 10, 10),
]


@pytest.mark.parametrize('stem, num_species, num_ts, num_reactions', NETWORK_COUNTS)
def test_counts(stem, num_species, num_ts, num_reactions):
    """Test that species / transitionState / reaction counts match the real fixtures."""
    path = os.path.join(PDEP_NETWORK_DIR, f'{stem}.py')
    network = parse_pdep_network_file(path=path)
    assert len(network.species_labels) == num_species
    assert len(network.transition_state_labels) == num_ts
    assert len(network.path_reactions) == num_reactions


@pytest.mark.parametrize('stem, _num_species, _num_ts, _num_reactions', NETWORK_COUNTS)
def test_network_id_equals_file_stem(stem, _num_species, _num_ts, _num_reactions):
    """Test that ``network_id`` equals the file stem."""
    path = os.path.join(PDEP_NETWORK_DIR, f'{stem}.py')
    network = parse_pdep_network_file(path=path)
    assert network.network_id == stem
    assert network.path == path


def test_network4_2_reaction1():
    """Test a specific known reaction (reaction1) from network4_2.py."""
    path = os.path.join(PDEP_NETWORK_DIR, 'network4_2.py')
    network = parse_pdep_network_file(path=path)
    reaction1 = network.path_reactions[0]
    assert reaction1.label == 'reaction1'
    assert reaction1.reactants == ('CH2(S)(53)', 'C3rad(4)')
    assert reaction1.products == ('C4rad(5)',)
    assert reaction1.transition_state == 'TS1'
    assert reaction1.kinetics_type == 'Arrhenius'
    assert 'Estimated using template [carbene;C_pri]' in reaction1.kinetics_comment
    assert '\n' in reaction1.kinetics_comment


def test_network4_2_reaction2():
    """Test that reaction2 of network4_2.py has the expected kinetics comment."""
    path = os.path.join(PDEP_NETWORK_DIR, 'network4_2.py')
    network = parse_pdep_network_file(path=path)
    reaction2 = network.path_reactions[1]
    assert reaction2.label == 'reaction2'
    assert 'Exact match found for rate rule' in reaction2.kinetics_comment


def test_network4_2_path_reactions_by_ts():
    """Test ``path_reactions_by_ts`` on network4_2.py: 10 keys, each a 1-tuple."""
    path = os.path.join(PDEP_NETWORK_DIR, 'network4_2.py')
    network = parse_pdep_network_file(path=path)
    by_ts = network.path_reactions_by_ts()
    expected_keys = {f'TS{i}' for i in range(1, 11)}
    assert set(by_ts.keys()) == expected_keys
    for ts_label, reactions in by_ts.items():
        assert isinstance(reactions, tuple)
        assert len(reactions) == 1


SHARED_TS_TEXT = '''
reaction(
    label = 'reactionA',
    reactants = ['S1', 'S2'],
    products = ['S3'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""first"""),
)

reaction(
    label = 'reactionB',
    reactants = ['S3'],
    products = ['S4', 'S5'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(2.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""second"""),
)
'''


def test_shared_transition_state():
    """Test that two path reactions sharing one TS label map to a 2-tuple."""
    network = parse_pdep_network_text(text=SHARED_TS_TEXT, network_id='synthetic_shared_ts')
    by_ts = network.path_reactions_by_ts()
    assert set(by_ts.keys()) == {'TS1'}
    assert len(by_ts['TS1']) == 2
    assert {r.label for r in by_ts['TS1']} == {'reactionA', 'reactionB'}


TRICKY_COMMENT_TEXT = '''
reaction(
    label = 'reactionTricky',
    reactants = ['S1'],
    products = ['S2'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""this comment mentions reaction( inside it
and has a "quote-ish" edge case nearby"""),
)
'''


def test_comment_containing_reaction_call_substring():
    """Test that a comment containing the literal substring 'reaction(' is still parsed as one reaction."""
    network = parse_pdep_network_text(text=TRICKY_COMMENT_TEXT, network_id='synthetic_tricky')
    assert len(network.path_reactions) == 1
    reaction = network.path_reactions[0]
    assert reaction.label == 'reactionTricky'
    assert 'reaction(' in reaction.kinetics_comment
    assert 'and has a' in reaction.kinetics_comment


NO_TS_NO_COMMENT_TEXT = """
reaction(
    label = 'reactionBare',
    reactants = ['S1'],
    products = ['S2'],
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K')),
)
"""


def test_reaction_without_ts_or_comment():
    """Test a synthetic reaction with no transitionState and no comment keyword."""
    network = parse_pdep_network_text(text=NO_TS_NO_COMMENT_TEXT, network_id='synthetic_bare')
    reaction = network.path_reactions[0]
    assert reaction.transition_state is None
    assert reaction.kinetics_comment == ''
    assert reaction.kinetics_type == 'Arrhenius'


def test_parse_pdep_network_file_invalid_syntax_raises_value_error(tmp_path):
    """Test that a syntactically invalid file raises a ValueError."""
    bad_file = tmp_path / 'bad_network.py'
    bad_file.write_text('reaction(\n    label = \n')
    with pytest.raises(ValueError):
        parse_pdep_network_file(path=str(bad_file))


ANTI_EXEC_TEXT = """
raise RuntimeError('executed!')

reaction(
    label = 'reactionSafe',
    reactants = ['S1'],
    products = ['S2'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment=\"\"\"safe\"\"\"),
)
"""


def test_anti_exec_regression():
    """Test that a top-level side effect statement is never executed (only ast.parse is used)."""
    network = parse_pdep_network_text(text=ANTI_EXEC_TEXT, network_id='synthetic_anti_exec')
    assert len(network.path_reactions) == 1
    reaction = network.path_reactions[0]
    assert reaction.label == 'reactionSafe'
    assert reaction.reactants == ('S1',)
    assert reaction.products == ('S2',)
    assert reaction.kinetics_comment == 'safe'


def test_empty_network_raises_value_error():
    """Test that a file with no reaction(...) calls raises a ValueError."""
    text = "species(label = 'S1', structure = SMILES('[H]'))\n"
    with pytest.raises(ValueError):
        parse_pdep_network_text(text=text, network_id='synthetic_empty')


# FIX5: no in-repo fixture (network1_1/network4_1/network4_2) contains a MultiArrhenius
# (or other composite) kinetics call, so this case is necessarily covered with a synthetic,
# hand-written fixture rather than a real one.
MULTI_ARRHENIUS_TEXT = '''
reaction(
    label = 'reactionComposite',
    reactants = ['S1'],
    products = ['S2'],
    transitionState = 'TS1',
    kinetics = MultiArrhenius(
        arrhenius = [
            Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""low-T branch"""),
            Arrhenius(A=(2.0,'s^-1'), n=0, Ea=(10,'kJ/mol'), T0=(1,'K'), comment="""high-T branch"""),
        ],
    ),
)
'''


def test_multi_arrhenius_nested_comments_are_collected():
    """Test that comment= keywords nested inside a composite kinetics call (e.g. MultiArrhenius
    wrapping several Arrhenius components) are not silently dropped when the top-level kinetics
    call itself has no comment= keyword. This is a synthetic fixture: no in-repo network file
    (network1_1/network4_1/network4_2) contains a MultiArrhenius kinetics call.
    """
    network = parse_pdep_network_text(text=MULTI_ARRHENIUS_TEXT, network_id='synthetic_multi_arrhenius')
    reaction = network.path_reactions[0]
    assert reaction.kinetics_type == 'MultiArrhenius'
    assert 'low-T branch' in reaction.kinetics_comment
    assert 'high-T branch' in reaction.kinetics_comment
    assert '\n' in reaction.kinetics_comment


def test_parse_arkane_pdep_output_file_success_fixture():
    """Test that the real clean-ME-solve output.py parses to one PDepArkaneReaction with finite
    Chebyshev coefficients (including legitimately negative ones).
    """
    path = os.path.join(PDEP_ME_DIR, 'success', 'output.py')
    reactions = parse_arkane_pdep_output_file(path=path)
    assert isinstance(reactions, tuple)
    assert len(reactions) == 1
    reaction = reactions[0]
    assert isinstance(reaction, PDepArkaneReaction)
    assert reaction.reactants == ('C1rad(2)',)
    assert reaction.products == ('H(34)', 'CH2(T)(48)')
    assert reaction.kinetics_type == 'Chebyshev'
    assert reaction.kinetics_params['kunits'] == 's^-1'
    assert reaction.kinetics_params['Tmin'] == (300, 'K')
    coeffs = reaction.kinetics_params['coeffs']
    assert coeffs[0][0] == -31.1515
    # Legitimately negative coefficients must survive unmodified (log-space fit coefficients,
    # not k(T,P) itself).
    assert any(row[0] < 0 for row in coeffs)
    assert all(v is not None for v in reaction.numeric_values)


def test_parse_arkane_pdep_output_file_soft_failure_cse_none_survives():
    """Test that the real CSE soft-failure output.py's all-None Chebyshev coeffs survive
    parsing as literal None values, not skipped and not coerced to 0 or NaN.
    """
    path = os.path.join(PDEP_ME_DIR, 'soft_failure_cse', 'output.py')
    reactions = parse_arkane_pdep_output_file(path=path)
    assert len(reactions) == 1
    reaction = reactions[0]
    assert reaction.kinetics_type == 'Chebyshev'
    coeffs = reaction.kinetics_params['coeffs']
    assert len(coeffs) == 6
    assert all(v is None for row in coeffs for v in row)
    assert len(reaction.numeric_values) > 0
    # The coeffs entries are None (24 of them: 6x4), but numeric_values also carries the
    # perfectly finite Tmin/Tmax/Pmin/Pmax quantity numbers, so not *every* entry is None.
    assert any(v is None for v in reaction.numeric_values)
    assert sum(1 for v in reaction.numeric_values if v is None) == 24


def test_parse_arkane_pdep_output_file_hard_failure_empty_tuple():
    """Test that the real 0-byte hard-failure output.py parses to an empty tuple, not an error.

    This is the specific behavior that lets the ME-success gate distinguish "crashed before
    writing anything" from "wrote something bad": an empty file must not raise ValueError the
    way an empty RMG pdep network file does in parse_pdep_network_file.
    """
    path = os.path.join(PDEP_ME_DIR, 'hard_failure', 'output.py')
    reactions = parse_arkane_pdep_output_file(path=path)
    assert reactions == tuple()


def test_parse_arkane_pdep_output_text_empty_or_comment_only_yields_empty_tuple():
    """Test that an empty string and a comment-only string both parse to an empty tuple rather
    than raising, mirroring the hard-failure fixture case but for synthetic inputs.

    Regression this guards: a naive port of parse_pdep_network_text's "raise ValueError if no
    reactions" behavior into the Arkane-output parser would make this raise instead of returning
    an empty tuple, breaking the ME-success gate's ability to treat emptiness as informative
    rather than exceptional.
    """
    assert parse_arkane_pdep_output_text(text='') == tuple()
    assert parse_arkane_pdep_output_text(text='# just a comment\n# nothing else\n') == tuple()


ARKANE_OUTPUT_ANTI_EXEC_TEXT = """
raise RuntimeError('executed!')

pdepreaction(
    reactants = ['A'],
    products = ['B', 'C'],
    kinetics = Chebyshev(
        coeffs = [[1.0, 2.0], [None, 4.0]],
        kunits = 's^-1',
        Tmin = (300, 'K'),
        Tmax = (2100, 'K'),
    ),
)
"""


def test_parse_arkane_pdep_output_text_anti_exec_and_nested_none():
    """Test that a top-level side-effect statement is never executed (only ast.parse is used),
    and that a None nested two levels deep inside coeffs survives parsing.
    """
    reactions = parse_arkane_pdep_output_text(text=ARKANE_OUTPUT_ANTI_EXEC_TEXT)
    assert len(reactions) == 1
    reaction = reactions[0]
    assert reaction.reactants == ('A',)
    assert reaction.products == ('B', 'C')
    assert reaction.kinetics_type == 'Chebyshev'
    assert reaction.kinetics_params['coeffs'] == [[1.0, 2.0], [None, 4.0]]
    assert None in reaction.numeric_values
    assert 1.0 in reaction.numeric_values
    assert 300 in reaction.numeric_values


def test_parse_arkane_pdep_output_text_invalid_syntax_raises_value_error():
    """Test that syntactically invalid Arkane output text raises a ValueError, matching the
    RMG pdep network parser's behavior for malformed input (as opposed to merely-empty input,
    which must not raise).
    """
    with pytest.raises(ValueError):
        parse_arkane_pdep_output_text(text='pdepreaction(\n    reactants = \n')


def test_rate_payload_numeric_values_exclude_bounds_and_metadata():
    """Test that ``rate_payload_numeric_values`` carries only the rate payload's numeric leaves
    (``coeffs`` here), never the Tmin/Tmax/Pmin/Pmax bounds or metadata numbers.

    Regression this guards: the combined ``numeric_values`` field mixes the (always finite)
    T/P bounds with the actual rate coefficients, which let a payload-free kinetics call pass
    the ME-success gate on its bounds alone.
    """
    text = ("pdepreaction(\n"
            "    reactants=['A'], products=['B'],\n"
            "    kinetics=Chebyshev(coeffs=[[1.0, -2.5]], kunits='s^-1',\n"
            "                       Tmin=(300,'K'), Tmax=(2100,'K'),\n"
            "                       Pmin=(0.1,'bar'), Pmax=(100,'bar')),\n"
            ")\n")
    reactions = parse_arkane_pdep_output_text(text=text)
    assert len(reactions) == 1
    reaction = reactions[0]
    assert reaction.rate_payload_numeric_values == (1.0, -2.5)
    # The combined field keeps its original meaning (payload AND bounds).
    assert 300 in reaction.numeric_values
    assert 300 not in reaction.rate_payload_numeric_values
    assert reaction.missing_kinetics_keys == tuple()


def test_missing_kinetics_keys_records_unparseable_coeffs():
    """Test that a kinetics keyword omitted because it could not be literal-evaluated (a bare
    ``nan`` name inside ``coeffs``) is recorded in ``missing_kinetics_keys``, so an omitted
    ``coeffs`` is distinguishable from a legitimately absent one.
    """
    text = ("pdepreaction(\n"
            "    reactants=['A'], products=['B'],\n"
            "    kinetics=Chebyshev(coeffs=[[nan]], kunits='s^-1',\n"
            "                       Tmin=(300,'K'), Tmax=(2100,'K'),\n"
            "                       Pmin=(0.1,'bar'), Pmax=(100,'bar')),\n"
            ")\n")
    reactions = parse_arkane_pdep_output_text(text=text)
    assert len(reactions) == 1
    reaction = reactions[0]
    assert 'coeffs' in reaction.missing_kinetics_keys
    assert 'coeffs' not in reaction.kinetics_params
    # The unresolvable payload leaf surfaces as None (non-finite), never silently vanishes.
    assert any(v is None for v in reaction.rate_payload_numeric_values) \
        or reaction.rate_payload_numeric_values == tuple()


# Product channel derivation
#
# RMG never writes a ``products=`` keyword into a generated network file; Arkane derives the
# product channels from the path reactions instead. The expected values below are MEASURED, not
# assumed: a real Arkane MSC run on ``network4_1.py`` named these three channels in its log, and
# the ``output.py`` it produced is checked in at ``tests/data/pdep_me/success_multi/`` with
# exactly 12 live ``pdepreaction(...)`` entries.
NETWORK4_1_PRODUCT_CHANNELS = [{'CH2(S)(53)', 'C3rad(4)'}, {'CH2(T)(48)', 'C3rad(4)'}, {'butyl_2(67)'}]
NETWORK4_1_EXPECTED_NET_REACTIONS = 12

PDEP_NETWORK_VARIANTS_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_network_variants')


def test_product_channels_are_derived_when_the_keyword_is_absent():
    """Test that product channels absent from ``network(...)`` are derived from the path reactions.

    ``network4_1.py`` declares one isomer and two reactant channels and no ``products=`` keyword.
    Reading the keyword literally yields zero product channels, which would under-count the net
    reactions Arkane writes by a factor of four.
    """
    network = parse_pdep_network_file(path=os.path.join(PDEP_NETWORK_DIR, 'network4_1.py'))
    assert network.isomers == ('C4rad(5)',)
    assert len(network.reactant_channels) == 2
    assert network.product_channels_declared is False
    assert [set(channel) for channel in network.product_channels] == NETWORK4_1_PRODUCT_CHANNELS


def test_expected_net_reaction_count_matches_a_real_arkane_output():
    """Test the predicted net reaction count against a real Arkane ``output.py``.

    This is the load-bearing test of the derivation: it ties the count computed from the network
    file to the number of live ``pdepreaction(...)`` entries Arkane actually wrote for that very
    network, rather than to a number restated from the source.
    """
    network = parse_pdep_network_file(path=os.path.join(PDEP_NETWORK_DIR, 'network4_1.py'))
    assert network.expected_net_reaction_count() == NETWORK4_1_EXPECTED_NET_REACTIONS
    reactions = parse_arkane_pdep_output_file(path=os.path.join(PDEP_ME_DIR, 'success_multi', 'output.py'))
    assert len(reactions) == network.expected_net_reaction_count()


def test_expected_net_reaction_count_for_a_single_channel_network():
    """Test the degenerate one-isomer network, whose Arkane output holds a single entry."""
    network = parse_pdep_network_file(path=os.path.join(PDEP_NETWORK_DIR, 'network1_1.py'))
    assert len(network.isomers) == 1
    assert network.reactant_channels == tuple()
    assert len(network.product_channels) == 1
    assert network.expected_net_reaction_count() == 1
    reactions = parse_arkane_pdep_output_file(path=os.path.join(PDEP_ME_DIR, 'success', 'output.py'))
    assert len(reactions) == network.expected_net_reaction_count()


def test_product_channels_are_derived_for_every_real_network_fixture():
    """Test that the derivation yields at least one product channel for each real fixture.

    No real RMG-generated network file declares ``products=``, so a fixture reporting zero product
    channels means the derivation silently did nothing.
    """
    for stem, _, _, _ in NETWORK_COUNTS:
        network = parse_pdep_network_file(path=os.path.join(PDEP_NETWORK_DIR, f'{stem}.py'))
        assert network.product_channels_declared is False, stem
        assert len(network.product_channels) > 0, stem


def test_declared_product_channels_are_used_verbatim():
    """Test that an explicit ``products=`` keyword wins over the derivation."""
    network = parse_pdep_network_file(path=os.path.join(PDEP_NETWORK_VARIANTS_DIR,
                                                        'network_explicit_products.py'))
    assert network.product_channels_declared is True
    assert [set(channel) for channel in network.product_channels] == [{'B', 'C'}, {'D'}]
    # One isomer and two declared product channels, none of them overlapping.
    assert network.expected_net_reaction_count() == 2


def test_channel_labels_are_canonically_ordered():
    """Test that multi-species channels are sorted, so membership comparisons cannot silently miss."""
    network = parse_pdep_network_file(path=os.path.join(PDEP_NETWORK_DIR, 'network4_1.py'))
    for channel in network.reactant_channels + network.product_channels:
        assert list(channel) == sorted(channel)


def test_the_two_arity_branches_of_the_derivation_are_not_interchangeable():
    """Test that the derivation reproduces Arkane's asymmetric membership tests exactly.

    Arkane tests a unimolecular side against the isomers only, and a multi-species side against
    the reactant channels only -- despite its own comment claiming both are tested against both.
    Two consequences follow, and both are pinned here because unifying the branches looks like a
    harmless cleanup and silently under-counts the network:

    * ``['D']`` is a declared reactant channel, yet becomes a product channel as well, because the
      unimolecular branch never consults the reactant channels.
    * ``['B', 'C']`` survives even though ``'B'`` is an isomer, because the multi-species branch
      never consults the isomers.
    """
    network = parse_pdep_network_file(path=os.path.join(PDEP_NETWORK_VARIANTS_DIR,
                                                        'network_arity_asymmetry.py'))
    assert network.isomers == ('A', 'B')
    assert network.reactant_channels == (('D',),)
    assert network.product_channels_declared is False
    assert [set(channel) for channel in network.product_channels] == [{'B', 'C'}, {'D'}]
    # ``('D',)`` is both a source and a destination, so the three entries reaching it are
    # commented out as duplicates of the three leaving it. The closed form would say 9.
    assert network.expected_net_reaction_count() == 7


def test_a_channel_that_is_both_a_source_and_a_destination_is_not_double_counted():
    """Test the net reaction count when one configuration appears in both channel lists.

    Arkane suppresses a duplicate by comparing each candidate against EVERY previously printed
    reaction, by label and in either direction -- not only against those inside the
    isomer/reactant block. A configuration that is both a declared unimolecular reactant channel
    and a derived product channel therefore has the entries reaching it commented out as
    duplicates of the entries leaving it.

    The numbers are measured, not derived: ``tests/data/pdep_me/overlapping_channels/`` is a real
    Arkane MSC run of ``network4_1`` with ``butyl_2(67)`` additionally declared as a reactant
    channel. It holds 15 live entries and 9 commented ones. The closed form ``S*(S-1)/2 + S*P``
    predicts 18, which would reject that complete solve.
    """
    network = parse_pdep_network_file(path=os.path.join(PDEP_ME_DIR, 'overlapping_channels', 'input.py'))
    assert ('butyl_2(67)',) in network.reactant_channels
    assert ('butyl_2(67)',) in network.product_channels
    assert network.expected_net_reaction_count() == 15
    reactions = parse_arkane_pdep_output_file(path=os.path.join(PDEP_ME_DIR, 'overlapping_channels', 'output.py'))
    assert len(reactions) == network.expected_net_reaction_count()


def test_a_non_literal_products_keyword_is_refused_rather_than_derived():
    """Test that a present-but-unevaluable ``products`` keyword raises instead of being derived.

    ``_literal_or_none`` returns ``None`` both for an absent keyword and for one whose value is a
    name or a call. Conflating the two would let a file that DECLARES its product channels be
    silently treated as one that omits them, producing a net reaction count for a topology the
    file never described.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState='TS1')\n"
            "network(label='n', isomers=['A'], reactants=[], products=SOME_NAME)\n")
    with pytest.raises(ValueError, match='products'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')


def test_a_non_literal_isomers_keyword_is_refused():
    """Test that a present-but-unevaluable ``isomers`` keyword raises instead of degrading to empty.

    Degrading here is worse than for ``products``: the isomers are the source side of the net
    reaction loop AND the thing every unimolecular path reaction side is tested against when
    product channels are derived. An empty isomer list therefore both under-counts the sources
    and misclassifies every unimolecular side as a product channel, so the expected net reaction
    count describes a topology the file never declared -- and a complete solve gets rejected.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState='TS1')\n"
            "network(label='n', isomers=SOME_NAME, reactants=[], products=[])\n")
    with pytest.raises(ValueError, match='isomers'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')


def test_a_non_literal_reactants_keyword_is_refused():
    """Test that a present-but-unevaluable ``reactants`` keyword raises instead of degrading to empty.

    A comprehension is used here rather than a bare name so that the guard is shown to key off
    "``_literal_or_none`` could not evaluate it", not off any one node type.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState='TS1')\n"
            "network(label='n', isomers=['A'], reactants=[[s] for s in SOURCES], products=[])\n")
    with pytest.raises(ValueError, match='reactants'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')


def test_absent_isomers_and_reactants_keywords_still_parse():
    """Test that ABSENT ``isomers``/``reactants`` keywords remain legitimate and parse as empty.

    The refusal above must key off "present but unevaluable", not off "falsy": omitting these
    keywords is normal, so conflating the two cases would reject ordinary files.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState='TS1')\n"
            "network(label='n')\n")
    network = parse_pdep_network_text(text=text, network_id='n', path='<test>')
    assert network.isomers == tuple()
    assert network.reactant_channels == tuple()


def test_a_non_literal_reactants_keyword_in_reaction_is_refused():
    """Test that a present-but-unevaluable ``reactants`` keyword on ``reaction(...)`` raises.

    ``_parse_reaction`` used to mirror the pre-fix ``network(...)`` defect: ``tuple(_literal_or_none(
    kwargs.get('reactants')) or ())`` collapsed a name or call in ``reactants`` position to the
    empty tuple, silently contributing NO product channel for that side and letting ``t3.main``
    queue an empty ``T3Reaction`` for a TS search.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=REACTANT_LIST, products=['B'], transitionState='TS1')\n"
            "network(label='n', isomers=['A'], reactants=[], products=[])\n")
    with pytest.raises(ValueError, match='reactants'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')


def test_a_non_literal_products_keyword_in_reaction_is_refused():
    """Test that a present-but-unevaluable ``products`` keyword on ``reaction(...)`` raises.

    Same defect as the ``reactants`` case above, on the other side of the reaction.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=PRODUCT_LIST, transitionState='TS1')\n"
            "network(label='n', isomers=['A'], reactants=[], products=[])\n")
    with pytest.raises(ValueError, match='products'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')


def test_a_non_literal_transition_state_keyword_is_refused():
    """Test that a present-but-unevaluable ``transitionState`` keyword on ``reaction(...)`` raises.

    An ABSENT ``transitionState`` is legitimate (``PDepNetwork.path_reactions_by_ts`` relies on
    it to exclude such reactions from its map), but one that is present and unevaluable must not
    collapse to the same ``None`` -- that would silently and wrongly exclude a path reaction that
    the file DID associate with a transition state.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState=TS_LABELS[2])\n"
            "network(label='n', isomers=['A'], reactants=[], products=[])\n")
    with pytest.raises(ValueError, match='transitionState'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')


def test_absent_transition_state_keyword_still_parses():
    """Test that an ABSENT ``transitionState`` keyword remains legitimate and parses as ``None``.

    The refusal above must key off "present but unevaluable", not off "missing", so an ordinary
    path reaction without an associated transition state must still parse.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=['B'])\n"
            "network(label='n', isomers=['A'], reactants=[], products=[])\n")
    network = parse_pdep_network_text(text=text, network_id='n', path='<test>')
    assert network.path_reactions[0].transition_state is None


def test_a_non_literal_species_label_is_refused():
    """Test that a present-but-unevaluable ``label`` keyword on ``species(...)`` raises.

    The pre-fix handler did ``if label is not None: species_labels.append(label)``, so an
    unevaluable label silently vanished from ``species_labels`` instead of failing closed.
    """
    text = ("species(label=SOME_NAME, E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState='TS1')\n"
            "network(label='n', isomers=['A'], reactants=[], products=[])\n")
    with pytest.raises(ValueError, match='label'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')


def test_a_non_literal_transition_state_label_is_refused():
    """Test that a present-but-unevaluable ``label`` keyword on ``transitionState(...)`` raises.

    Same defect as the ``species(...)`` case above.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "transitionState(label=SOME_NAME, E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState='TS1')\n"
            "network(label='n', isomers=['A'], reactants=[], products=[])\n")
    with pytest.raises(ValueError, match='label'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')


def test_an_explicit_none_reactants_keyword_on_reaction_is_refused():
    """Test that an EXPLICITLY WRITTEN ``reactants=None`` on ``reaction(...)`` raises.

    ``network(...)``'s ``isomers``/``reactants``/``products`` guard already refuses an explicit
    ``None`` literal rather than treating it as legitimate absence, on the reasoning that the
    keyword was written so its intent is not "omitted". ``_literal_or_raise`` must be consistent
    with that: a real RMG-generated network file never emits an explicit ``None`` for
    ``reactants``, so refusing it costs nothing on legitimate input and closes a hole on
    hand-edited/unusual input.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=None, products=['B'], transitionState='TS1')\n"
            "network(label='n', isomers=['A'], reactants=[], products=[])\n")
    with pytest.raises(ValueError, match='reactants'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')


def test_an_explicit_none_products_keyword_on_reaction_is_refused():
    """Test that an EXPLICITLY WRITTEN ``products=None`` on ``reaction(...)`` raises.

    Same defect as the ``reactants`` case above, on the other side of the reaction.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=None, transitionState='TS1')\n"
            "network(label='n', isomers=['A'], reactants=[], products=[])\n")
    with pytest.raises(ValueError, match='products'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')


def test_an_explicit_none_transition_state_keyword_is_refused():
    """Test that an EXPLICITLY WRITTEN ``transitionState=None`` on ``reaction(...)`` raises.

    An ABSENT ``transitionState`` is legitimate (the keyword was never written), but one that IS
    written, even as literal ``None``, must raise rather than be read as the same absent case: the
    keyword was written, so its intent is not "omitted" and guessing is not ours to do.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState=None)\n"
            "network(label='n', isomers=['A'], reactants=[], products=[])\n")
    with pytest.raises(ValueError, match='transitionState'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')


def test_an_explicit_none_species_label_is_refused():
    """Test that an EXPLICITLY WRITTEN ``label=None`` on ``species(...)`` raises.

    Same rationale as the other three explicit-``None`` guards above.
    """
    text = ("species(label=None, E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState='TS1')\n"
            "network(label='n', isomers=['A'], reactants=[], products=[])\n")
    with pytest.raises(ValueError, match='label'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')


def test_an_explicit_none_transition_state_label_is_refused():
    """Test that an EXPLICITLY WRITTEN ``label=None`` on ``transitionState(...)`` raises.

    Same rationale as the other three explicit-``None`` guards above.
    """
    text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
            "transitionState(label=None, E0=(0.0,'kJ/mol'))\n"
            "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState='TS1')\n"
            "network(label='n', isomers=['A'], reactants=[], products=[])\n")
    with pytest.raises(ValueError, match='label'):
        parse_pdep_network_text(text=text, network_id='n', path='<test>')
class TestKwargsUnpackingOnRecognizedCallsIsRefused:
    """
    round-30 P2. ``_call_keywords`` maps an ``ast.Call``'s keyword arguments to their (still-AST)
    value nodes, but ``ast.keyword.arg`` is ``None`` for a ``**kwargs`` unpacking (e.g.
    ``reaction(**payload)``). Before this fix, ``_call_keywords`` silently dropped any keyword
    whose ``arg`` was ``None`` -- so a network file containing ``reaction(**payload)`` (or
    ``species(**payload)``, ``transitionState(**payload)``, ``network(**payload)``) "parsed"
    successfully as a recognized DSL call with every field absent, instead of being refused. This
    is a fail-open: this module never executes the file (see the module docstring), so it has no
    ``payload`` to expand and cannot know which keywords the unpacking would actually supply --
    silence here reports a network that does not describe what the file actually declares. The
    fix makes ``_call_keywords`` raise ``ValueError`` when a ``**kwargs`` unpacking appears on a
    RECOGNIZED top-level DSL call (``RECOGNIZED_TOP_LEVEL_CALLS``), fail-closed to match how the
    rest of this module already refuses unresolvable input (e.g. ``_literal_or_raise``).
    """

    def test_a_kwargs_unpacking_on_reaction_is_refused(self):
        # A bare ``**payload`` unpacking supplies zero literal keywords: a naive fix that only
        # checks "did we get zero keywords" would also catch this case, but so would the pre-fix
        # (defective) code accidentally happen to still raise elsewhere for unrelated reasons, so
        # this case alone is not sufficient evidence of a real fix -- see the mixed-call test below.
        text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
                "transitionState(label='TS1', E0=(0.0,'kJ/mol'))\n"
                "reaction(**payload)\n"
                "network(label='n', isomers=['A'], reactants=[], products=[])\n")
        with pytest.raises(ValueError, match='kwargs'):
            parse_pdep_network_text(text=text, network_id='n', path='<test>')

    def test_a_kwargs_unpacking_on_species_is_refused(self):
        text = ("species(**payload)\n"
                "transitionState(label='TS1', E0=(0.0,'kJ/mol'))\n"
                "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState='TS1')\n"
                "network(label='n', isomers=['A'], reactants=[], products=[])\n")
        with pytest.raises(ValueError, match='kwargs'):
            parse_pdep_network_text(text=text, network_id='n', path='<test>')

    def test_a_kwargs_unpacking_mixed_with_explicit_keywords_is_refused(self):
        # The case a naive fix that only checks "did ``_call_keywords`` return an empty dict"
        # would miss: real, literal-evaluable keywords ARE present alongside the unpacking, so
        # the resulting dict is non-empty even though the unpacking's own keywords are unknown.
        text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
                "transitionState(label='TS1', E0=(0.0,'kJ/mol'))\n"
                "reaction(reactants=['A'], **payload)\n"
                "network(label='n', isomers=['A'], reactants=[], products=[])\n")
        with pytest.raises(ValueError, match='kwargs'):
            parse_pdep_network_text(text=text, network_id='n', path='<test>')

    def test_an_ordinary_call_with_only_explicit_keywords_still_parses(self):
        # Over-refusal guard: a normal network with no ``**kwargs`` unpacking anywhere must still
        # parse fine -- the fix must not refuse calls that never use ``**kwargs`` at all.
        text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
                "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState='TS1')\n"
                "network(label='n', isomers=['A'], reactants=[], products=[])\n")
        network = parse_pdep_network_text(text=text, network_id='n', path='<test>')
        assert network.species_labels == ('A',)
        assert network.path_reactions[0].reactants == ('A',)

    def test_a_kwargs_unpacking_on_a_NESTED_call_is_also_refused(self):
        """
        Deliberately inverted from "is not refused". The narrow scope was the wrong contract.

        This test was first written to assert that a ``**kwargs`` unpacking on a nested call --
        ``kinetics=Arrhenius(**payload)`` -- is ACCEPTED, on the grounds that only the top-level DSL
        calls were reported. Its own rationale contained the refutation: it said the payload "is
        legitimately unresolvable without executing the file", which is the reason to refuse, not the
        reason to accept. Scoping the guard to top-level calls also required a ``call_name`` argument
        whose absence disabled it, i.e. a guard that is OPEN by default and stays correct only while
        every one of six call sites remembers to opt in.

        Refusing everywhere is free: across the 156 parseable files under RMG-Py's ``examples/`` and
        ``arkane/data/``, zero calls use a ``**kwargs`` unpacking anywhere -- these files are
        machine-written with explicit keywords. So the choice is not "strict vs. permissive on real
        inputs", it is "closed vs. open on inputs that never legitimately occur".
        """
        text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
                "reaction(label='reaction1', reactants=['A'], products=['B'], transitionState='TS1', "
                "kinetics=Arrhenius(**payload))\n"
                "network(label='n', isomers=['A'], reactants=[], products=[])\n")
        with pytest.raises(ValueError, match=r'\*\*kwargs'):
            parse_pdep_network_text(text=text, network_id='n', path='<test>')

    def test_the_refusal_does_not_depend_on_the_caller_passing_a_call_name(self):
        """
        The guard must not be openable by omitting an argument.

        ``call_name`` exists only so the message can name the offending call. If it also gated the
        check, then every future call site of ``_call_keywords`` would default to the fail-open
        behaviour, and this whole finding would silently reappear at the next one added.
        """
        tree = ast.parse('reaction(**payload)\n')
        call = tree.body[0].value
        with pytest.raises(ValueError, match=r'\*\*kwargs'):
            _call_keywords(call)


class TestPositionalArgumentsOnRecognizedCallsAreRefused:
    """
    The sibling of the ``**kwargs`` hole, in the same function, which that fix walked straight past.

    ``_call_keywords`` reads ``call.keywords`` and never looks at ``call.args``. Arkane's DSL
    functions accept positional arguments -- ``reaction(['A'], ['B'], transitionState='TS1')`` is a
    legal call to ``arkane/input.py``'s ``reaction(reactants, products, ...)`` -- so a file spelled
    that way parsed as a reaction with NO reactants and NO products. Exactly the same fail-open as an
    unpacking, one attribute over, and a demonstration that fixing a defect is not the same as fixing
    its class: the ``**kwargs`` commit widened the refusal to all six call sites and still never asked
    what else the function fails to read.

    Refused for the same measured reason: zero of the 156 parseable files under RMG-Py's ``examples/``
    and ``arkane/data/`` pass a positional argument to any of these calls -- they are machine-written
    with explicit keywords throughout.
    """

    @pytest.mark.parametrize('text,description', [
        ("species(label='A', E0=(0.0,'kJ/mol'))\n"
         "reaction(['A'], ['B'], transitionState='TS1')\n"
         "network(label='n', isomers=['A'], reactants=[], products=[])\n",
         'reactants and products passed positionally'),
        ("species('A', E0=(0.0,'kJ/mol'))\n"
         "reaction(label='r', reactants=['A'], products=['B'], transitionState='TS1')\n"
         "network(label='n', isomers=['A'], reactants=[], products=[])\n",
         'a species label passed positionally'),
        ("species(label='A', E0=(0.0,'kJ/mol'))\n"
         "reaction(label='r', reactants=['A'], products=['B'], transitionState='TS1')\n"
         "network('n', isomers=['A'], reactants=[], products=[])\n",
         'a network label passed positionally'),
    ])
    def test_a_positional_argument_on_a_dsl_call_is_refused(self, text, description):
        """Each of these silently loses the positional value today."""
        with pytest.raises(ValueError, match='positional'):
            parse_pdep_network_text(text=text, network_id='n', path='<test>')

    def test_the_refusal_does_not_depend_on_the_caller_passing_a_call_name(self):
        """Same reasoning as for the unpacking: the guard must not be openable by omission."""
        tree = ast.parse("reaction(['A'], ['B'])\n")
        with pytest.raises(ValueError, match='positional'):
            _call_keywords(tree.body[0].value)

    def test_an_all_keyword_call_is_still_accepted(self):
        """Over-refusal guard: the spelling every real RMG-written network file uses."""
        text = ("species(label='A', E0=(0.0,'kJ/mol'))\n"
                "reaction(label='r', reactants=['A'], products=['B'], transitionState='TS1')\n"
                "network(label='n', isomers=['A'], reactants=[], products=[])\n")
        network = parse_pdep_network_text(text=text, network_id='n', path='<test>')
        assert network.path_reactions[0].reactants == ('A',)
