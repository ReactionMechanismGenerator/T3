"""Tests for t3.pdep.barrierless."""

import pytest

from t3.pdep.barrierless import (BARRIERLESS_FAMILIES, classify_barrierless,
                                 rmg_family)
from t3.pdep.parser import PDepPathReaction


def _rxn(comment: str, transition_state: str | None = 'TS1') -> PDepPathReaction:
    """Build a path reaction carrying only the fields the classifier reads."""
    return PDepPathReaction(label='rxn', reactants=('A',), products=('B',),
                            transition_state=transition_state, kinetics_type='Arrhenius',
                            kinetics_comment=comment)


class TestRMGFamily(object):
    """Family extraction from the two forms RMG actually writes."""

    def test_node_estimated_form(self):
        """'... in family R_Recombination.' -- the node-estimated form."""
        comment = 'Estimated from node Root_N-1R->H_N-1CNOS->N in family R_Recombination.'
        assert rmg_family(_rxn(comment)) == 'R_Recombination'

    def test_rate_rule_form(self):
        """'family: 1,2_Insertion_CO' on its own line -- the rate-rule/training form."""
        comment = 'Estimated using an average for rate rule [CO;RO_H]\nEuclidian distance: 0\nfamily: 1,2_Insertion_CO'
        assert rmg_family(_rxn(comment)) == '1,2_Insertion_CO'

    def test_family_with_commas_and_digits_survives(self):
        """A family name is not [A-Za-z_]+: '1,3_Insertion_CO2' must come back whole."""
        comment = 'Matched reaction 0 H2 + CO2 <=> CH2O2 in 1,3_Insertion_CO2/training\nfamily: 1,3_Insertion_CO2'
        assert rmg_family(_rxn(comment)) == '1,3_Insertion_CO2'

    def test_family_with_a_plus_sign_survives(self):
        """A family name may contain '+': '1+2_Cycloaddition' (real family, network799_1
        reaction1) must not silently come back as None (N4)."""
        comment = ('Estimated using template [o_atom_singlet;multiplebond] for rate rule '
                  '[o_atom_singlet;mb_carbonyl]\nEuclidian distance = 1.0\nMultiplied by '
                  'reaction path degeneracy 2.0\nfamily: 1+2_Cycloaddition')
        assert rmg_family(_rxn(comment)) == '1+2_Cycloaddition'

    def test_no_family_returns_none(self):
        """A comment with no family is not guessed at."""
        assert rmg_family(_rxn('From a library.')) is None

    def test_empty_comment_returns_none(self):
        assert rmg_family(_rxn('')) is None


class TestClassifyBarrierless(object):
    """The verdict the loop gates QM on."""

    def test_r_recombination_is_barrierless(self):
        """The case that motivated this module: the live campaign spent its whole TS-search
        budget on HO + CHO <=> CH2O2, an R_Recombination, which has no saddle point at all."""
        comment = 'Estimated from node Root_1R->H in family R_Recombination.'
        verdict = classify_barrierless(_rxn(comment))
        assert verdict.is_barrierless is True
        assert verdict.family == 'R_Recombination'
        assert 'R_Recombination' in verdict.reason

    def test_insertion_family_is_not_barrierless(self):
        """1,2_Insertion_CO has a real saddle point -- this campaign converged one at -632.77 cm^-1."""
        verdict = classify_barrierless(_rxn('family: 1,2_Insertion_CO'))
        assert verdict.is_barrierless is False
        assert verdict.family == '1,2_Insertion_CO'

    def test_unknown_family_is_not_barrierless(self):
        """Fail OPEN, not closed. An unrecognized family gets its TS search: wrongly skipping a
        real barrier silently loses physics, while wrongly attempting a barrierless one merely
        wastes a job and is visible in the log."""
        verdict = classify_barrierless(_rxn('family: Some_New_Family'))
        assert verdict.is_barrierless is False
        assert 'not a known barrierless family' in verdict.reason

    def test_no_family_is_not_barrierless(self):
        verdict = classify_barrierless(_rxn('From a library.'))
        assert verdict.is_barrierless is False
        assert verdict.family is None

    def test_every_listed_family_classifies_as_barrierless(self):
        """The constant and the classifier must not drift apart."""
        for family in BARRIERLESS_FAMILIES:
            verdict = classify_barrierless(_rxn(f'family: {family}'))
            assert verdict.is_barrierless is True, f'{family} is listed but did not classify'

    def test_verdict_is_frozen(self):
        """A verdict is evidence; a caller must not edit it after the fact."""
        verdict = classify_barrierless(_rxn('family: R_Recombination'))
        with pytest.raises(Exception):
            verdict.is_barrierless = False
