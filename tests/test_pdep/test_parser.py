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
