#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_uncertainty module
"""

from t3.chem import KineticsMethod, T3Reaction
from t3.utils.uncertainty import is_this_kinetics_comment_uncertain, is_this_reaction_uncertain


class TestIsThisKineticsCommentUncertain:
    """Tests for the is_this_kinetics_comment_uncertain() function."""

    def test_none_comment_is_uncertain(self):
        """A ``None`` kinetics comment is uncertain (permissive default)."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment=None) is True

    def test_empty_comment_is_uncertain(self):
        """An empty kinetics comment is uncertain (permissive default)."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='') is True

    def test_library_kinetics_method_enum_is_not_uncertain(self):
        """A library reaction is never uncertain, regardless of its comment,
        when ``kinetics_method`` is given as the ``KineticsMethod`` enum."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='Library reaction: JetSurF2.0',
                                                  kinetics_method=KineticsMethod.LIBRARY) is False

    def test_library_kinetics_method_string_is_not_uncertain(self):
        """A library reaction is never uncertain when ``kinetics_method`` is given as its raw string value."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='Library reaction: JetSurF2.0',
                                                  kinetics_method='Library') is False

    def test_exact_match_comment_is_not_uncertain(self):
        """A rate-rule comment with an exact match is not uncertain."""
        comment = """Reaction index: Chemkin #10; RMG #6085
Template reaction: intra_H_migration
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0"""
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment,
                                                  kinetics_method=KineticsMethod.RATE_RULES) is False

    def test_estimated_using_template_for_rate_rule_is_uncertain(self):
        """A comment estimated using a template for a rate rule is uncertain."""
        comment = """Reaction index: Chemkin #7; RMG #6097
Template reaction: intra_H_migration
Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs)]"""
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment,
                                                  kinetics_method=KineticsMethod.RATE_RULES) is True

    def test_estimated_using_an_average_for_rate_rule_is_uncertain(self):
        """A comment estimated using an average for a rate rule is uncertain."""
        comment = """Reaction index: Chemkin #91; RMG #6629
Template reaction: H_Abstraction
Estimated using an average for rate rule [C/H2/NonDeC;C_rad/H/NonDeC]"""
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment,
                                                  kinetics_method=KineticsMethod.RATE_RULES) is True

    def test_no_kinetics_method_given_and_non_exact_match_comment_is_uncertain(self):
        """Without a kinetics method, an arbitrary comment with no certain marker is uncertain."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='some arbitrary comment') is True

    def test_kinetics_method_as_library_string_from_network_file(self):
        """A network-file caller with no T3Reaction object can still pass a raw method string."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='Library reaction: some_library',
                                                  kinetics_method='Library') is False

    def test_estimated_using_average_of_templates_is_uncertain(self):
        """A real ``Estimated using average of templates ...`` comment (from
        ``minimal_data/iteration_1/RMG/chemkin/chem_annotated.inp``) is uncertain."""
        comment = 'Estimated using average of templates [O/H/NonDeO;H_rad] + [H2O2;Y_rad] ' \
                  'for rate rule [H2O2;H_rad]'
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment) is True

    def test_library_reaction_text_in_comment_is_not_uncertain(self):
        """A comment stating ``Library reaction: ...`` is not uncertain even without a
        ``kinetics_method`` (e.g., a bare-comment caller such as an RMG pdep network file parser)."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='Library reaction: FFCM1(-)') is False

    def test_reaction_library_text_in_comment_is_not_uncertain(self):
        """A real comment from ``pdep_network/iteration_1/RMG/pdep/network1_1.py`` that only mentions
        ``Originally from reaction library ...`` (no ``kinetics_method``, no ``Library reaction:``
        marker) is not uncertain."""
        comment = 'Kinetics taken from the arrheniusHigh attribute of a Troe/Lindemann exprssion. ' \
                  'Originally from reaction library FFCM1(-)'
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment, kinetics_method=None) is False

    def test_matched_training_reaction_comment_is_not_uncertain(self):
        """A real comment from ``pdep_network/iteration_1/RMG/pdep/network4_2.py`` where the reaction
        itself IS the matched training reaction (``Matched reaction ... /training``) is not uncertain,
        even with no ``kinetics_method``."""
        comment = """Matched reaction 220 C2H4 + C2H5 <=> C4H9-2 in R_Addition_MultipleBond/training
This reaction matched rate rule [Cds-HH_Cds-HH;CsJ-CsHH]
family: R_Addition_MultipleBond"""
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment, kinetics_method=None) is False

    def test_from_training_reaction_with_exact_match_is_not_uncertain(self):
        """A real comment from ``pdep_network/iteration_1/RMG/pdep/network4_2.py`` where a rate rule
        was fit ``From training reaction ...`` is still not uncertain, because it also contains an
        ``Exact match found for rate rule`` marker -- despite mentioning "training reaction", this is
        NOT itself a matched training reaction (no ``Matched reaction`` marker)."""
        comment = """From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment, kinetics_method=None) is False

    def test_bare_training_word_without_matched_reaction_marker_is_uncertain(self):
        """A comment that merely mentions 'training' / '/training' without the ``Matched reaction``
        marker and without any other certain marker must NOT be treated as certain."""
        comment = 'Some rate rule fit from a training set entry, /training'
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment, kinetics_method=None) is True

    def test_library_method_alias_lowercase_string(self):
        """The raw string alias ``'library'`` (lowercase) is recognized as the LIBRARY method."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='some arbitrary comment',
                                                  kinetics_method='library') is False

    def test_library_method_alias_reaction_library_string(self):
        """The raw string alias ``'Reaction library'`` is recognized as the LIBRARY method."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='some arbitrary comment',
                                                  kinetics_method='Reaction library') is False

    def test_training_set_method_alias_string(self):
        """The raw string alias ``'Training Set'`` is recognized as the TRAINING_SET method."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='some arbitrary comment',
                                                  kinetics_method='Training Set') is False

    def test_training_set_method_alias_lowercase_string(self):
        """The raw string alias ``'training'`` (lowercase) is recognized as the TRAINING_SET method."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='some arbitrary comment',
                                                  kinetics_method='training') is False

    def test_qm_method_falls_through_to_comment_analysis(self):
        """A QM kinetics method is no longer a method-level short circuit (fix #6): T3's own
        kinetics parser never assigns QM, so it falls through to comment analysis like PDEP."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='some arbitrary comment',
                                                  kinetics_method=KineticsMethod.QM) is True
        assert is_this_kinetics_comment_uncertain(kinetics_comment='Exact match found for rate rule [x]',
                                                  kinetics_method=KineticsMethod.QM) is False

    def test_user_method_falls_through_to_comment_analysis(self):
        """A User kinetics method is no longer a method-level short circuit (fix #6): T3's own
        kinetics parser never assigns USER, so it falls through to comment analysis like PDEP."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='some arbitrary comment',
                                                  kinetics_method=KineticsMethod.USER) is True
        assert is_this_kinetics_comment_uncertain(kinetics_comment='Exact match found for rate rule [x]',
                                                  kinetics_method=KineticsMethod.USER) is False

    def test_pdep_method_falls_through_to_comment_analysis(self):
        """A PDep kinetics method does not short-circuit; it falls through to comment analysis."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='some arbitrary comment',
                                                  kinetics_method=KineticsMethod.PDEP) is True
        assert is_this_kinetics_comment_uncertain(kinetics_comment='Exact match found for rate rule [x]',
                                                  kinetics_method=KineticsMethod.PDEP) is False

    def test_unrecognized_method_string_falls_through_to_comment_analysis(self):
        """An unrecognized raw ``kinetics_method`` string must not raise and must not be forced to
        uncertain; it falls through to comment analysis (treated as unknown, like ``method=None``)."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='some arbitrary comment',
                                                  kinetics_method='not-a-real-method') is True
        assert is_this_kinetics_comment_uncertain(kinetics_comment='Exact match found for rate rule [x]',
                                                  kinetics_method='not-a-real-method') is False

    def test_continuation_split_exact_match_marker_is_not_uncertain(self):
        """A real comment from ``minimal_data/iteration_1/RMG/cantera_from_ck/chem_annotated.yaml``
        wraps 'Exact match found for rate rule' across a literal line continuation
        ('Exact match found for rate\\' then 'rule [...]'), so the old, longer substring
        'Exact match found for rate rule' misses it. The short marker 'Exact match found' must
        still catch it (RMG's own comment parser uses the short form)."""
        comment = ('Total Standard Deviation in ln(k): 11.5401827615\nExact match found for rate\\\n'
                  ' rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]\n'
                  'Euclidian distance = 0\nfamily: R_Recombination')
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment) is False

    def test_bm_rule_fitted_with_exact_match_is_uncertain(self):
        """A real comment from ``minimal_data/iteration_1/RMG/cantera_from_ck/chem_annotated.yaml``
        (also found in ``collision_rate_violators.log`` and ``determine_species``) has a BM rule
        fitted to training reactions at a tree node, and ALSO contains 'Exact match found'. Even
        though there's an exact match, an exact match to a fitted BM-rule node is not real data
        for this reaction, so this must be UNCERTAIN. (Every 'BM rule fitted to' occurrence in
        tests/data/ reports a 'Total Standard Deviation in ln(k)' of at least 1.089, well above 1.0,
        confirming these are never reliable rate-rule matches.)"""
        comment = ('BM rule fitted to 2 training reactions at node '
                  'Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O\n'
                  'Total Standard Deviation in ln(k): 11.5401827615\n'
                  'Exact match found for rate rule '
                  '[Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]\n'
                  'Euclidian distance = 0\nfamily: R_Recombination')
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment) is True

    def test_estimated_from_node_is_uncertain(self):
        """A real 'Estimated from node ...' comment (RMG's tree-rule estimate; e.g. from
        ``determine_species``/``minimal_data`` fixtures) is uncertain."""
        comment = 'Estimated from node Backbone0_N-2R!H-inRing_N-1R!H-inRing_Sp-2R!H-1R!H'
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment) is True

    def test_ea_raised_to_match_endothermicity_is_uncertain(self):
        """A real 'Ea raised from ... to match endothermicity of reaction.' comment is uncertain,
        even alongside a certain marker, since RMG has altered the barrier."""
        comment = ('Exact match found for rate rule [some_rule]\n'
                  'Ea raised from 0.0 to 11.2 kJ/mol to match endothermicity of reaction.')
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment) is True

    def test_ea_raised_noop_without_endothermicity_clause_is_not_uncertain(self):
        """A real no-op 'Ea raised from -0.0 to 0 kJ/mol.' comment never carries the
        endothermicity clause and must NOT be treated as an estimate; alongside a certain marker,
        this stays certain."""
        comment = 'Exact match found for rate rule [some_rule]\nEa raised from -0.0 to 0 kJ/mol.'
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment) is False

    def test_seed_mechanism_comment_is_not_uncertain(self):
        """A 'Seed mechanism: ...' comment is certain, constructed per RMG's Chemkin writer/reader
        format (rmgpy/chemkin.pyx, which treats a seed mechanism like a library reaction). There are
        zero occurrences of 'Seed mechanism:' in tests/data/, so this case is source-derived, not
        data-derived."""
        assert is_this_kinetics_comment_uncertain(kinetics_comment='Seed mechanism: BasicSeed') is False

    def test_pdep_high_p_limit_library_path_reaction_is_uncertain(self):
        """A PDep net reaction comment with a per-path-reaction 'High-P limit: ... (Library
        reaction: ...)' debug annotation (exact format from rmgpy/chemkin.pyx lines 1803-1808) must
        NOT be made certain merely because one high-P path reaction happens to be a library
        reaction; that per-path annotation must be stripped before marker scanning, leaving a bare
        'PDep reaction: ...' comment with no certain marker, which is uncertain by default."""
        comment = ('! PDep reaction: PDepNetwork #17\n'
                  '! High-P limit: C2H4 + H <=> C2H5 (Library reaction: FFCM1(-))\n'
                  '! High-P limit: C2H5 + O2 <=> C2H4 + HO2 (Template reaction: Disproportionation)')
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment) is True

    def test_pdep_high_p_limit_estimate_does_not_force_uncertain(self):
        """Symmetrically, a 'High-P limit: ...' line must not be allowed to force a PDep net
        reaction uncertain either -- it is stripped before both the certain scan AND the estimate
        scan. Here the net reaction comment itself is certain ('Exact match found'), and the only
        estimate-looking text ('Estimated using') lives inside a stripped High-P limit line, so the
        overall verdict must remain certain."""
        comment = ('! PDep reaction: PDepNetwork #17\n'
                  '! Exact match found for rate rule [some_rule]\n'
                  '! High-P limit: C2H4 + H <=> C2H5 (Template reaction: Estimated using template [x])')
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment) is False

    def test_matched_training_reaction_continuation_split_is_not_uncertain(self):
        """A real comment from ``pdep_network/iteration_1/RMG/pdep/network4_1.py`` where 'Matched
        reaction ... in <Family>/training' is split across lines (continuation) is still certain,
        thanks to the bounded, non-greedy ``MATCHED_TRAINING_REGEX``."""
        comment = ('Matched reaction 220 C2H4 + C2H5 <=> C4H9-2 in R_Addition_MultipleBond/training\n'
                  'This reaction matched rate rule [Cds-HH_Cds-HH;CsJ-CsHH]\n'
                  'family: R_Addition_MultipleBond')
        assert is_this_kinetics_comment_uncertain(kinetics_comment=comment) is False


class TestIsThisReactionUncertain:
    """Tests for the is_this_reaction_uncertain() function."""

    def test_none_reaction_returns_none(self):
        """A ``None`` reaction returns ``None``."""
        assert is_this_reaction_uncertain(reaction=None) is None

    def test_library_reaction_is_not_uncertain(self):
        """A library reaction is not uncertain."""
        reaction = T3Reaction(label='A+B <=> C+D',
                              kinetics_method=KineticsMethod.LIBRARY,
                              kinetics_comment='Library reaction: JetSurF2.0',
                              )
        assert is_this_reaction_uncertain(reaction=reaction) is False

    def test_exact_match_reaction_is_not_uncertain(self):
        """A rate-rule reaction with an exact match comment is not uncertain."""
        reaction = T3Reaction(label='A+B <=> C+D',
                              kinetics_method=KineticsMethod.RATE_RULES,
                              kinetics_comment='Exact match found for rate rule [some_rule]',
                              )
        assert is_this_reaction_uncertain(reaction=reaction) is False

    def test_rate_rule_reaction_without_exact_match_is_uncertain(self):
        """A rate-rule reaction without an exact match comment is uncertain."""
        reaction = T3Reaction(label='A+B <=> C+D',
                              kinetics_method=KineticsMethod.RATE_RULES,
                              kinetics_comment='Estimated using an average for rate rule [some_rule]',
                              )
        assert is_this_reaction_uncertain(reaction=reaction) is True

    def test_reaction_with_empty_comment_is_uncertain(self):
        """A reaction with an empty kinetics comment is uncertain (permissive default)."""
        reaction = T3Reaction(label='A+B <=> C+D',
                              kinetics_method=KineticsMethod.RATE_RULES,
                              kinetics_comment='',
                              )
        assert is_this_reaction_uncertain(reaction=reaction) is True
