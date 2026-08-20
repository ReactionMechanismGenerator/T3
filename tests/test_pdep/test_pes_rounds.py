"""Tests for t3.pdep.pes_rounds."""

import os
import pytest

from t3.pdep.parser import PDepNetwork, PDepPathReaction
from t3.pdep.pes_rounds import (CandidateSplit, PES_LOOP_DIAGRAM_FILENAME, QMCandidate,
                                SkippedChannel, attach_sensitivity_evidence, round_paths,
                                split_qm_candidates)


def _rxn(label: str, comment: str, ts: str | None = None) -> PDepPathReaction:
    return PDepPathReaction(label=label, reactants=('A',), products=('B',),
                            transition_state=ts if ts is not None else f'TS_{label}',
                            kinetics_type='Arrhenius', kinetics_comment=comment)


def _network(path_reactions) -> PDepNetwork:
    return PDepNetwork(network_id='network1_1', path='/abs/network1_1.py', label=None,
                       species_labels=(), transition_state_labels=(),
                       path_reactions=tuple(path_reactions), isomers=(),
                       reactant_channels=(), product_channels=())


class TestSplitQMCandidates(object):

    def test_barriered_reaction_without_qm_is_a_candidate(self):
        network = _network([_rxn('r1', 'family: 1,2_Insertion_CO')])
        split = split_qm_candidates(network, computed_ts_labels=frozenset())
        assert len(split.candidates) == 1
        assert split.candidates[0].ts_label == 'TS_r1'
        assert split.candidates[0].family == '1,2_Insertion_CO'
        assert split.skipped == ()

    def test_barrierless_reaction_is_skipped_with_a_reason(self):
        """The reason is recorded, never silently dropped."""
        network = _network([_rxn('r1', 'in family R_Recombination.')])
        split = split_qm_candidates(network, computed_ts_labels=frozenset())
        assert split.candidates == ()
        assert len(split.skipped) == 1
        assert 'R_Recombination' in split.skipped[0].reason
        assert 'no saddle point' in split.skipped[0].reason

    def test_already_computed_ts_is_not_a_candidate_again(self):
        """A TS with QM in hand is done; re-queueing it would spend the budget twice."""
        network = _network([_rxn('r1', 'family: 1,2_Insertion_CO')])
        split = split_qm_candidates(network, computed_ts_labels=frozenset({'TS_r1'}))
        assert split.candidates == ()
        assert len(split.skipped) == 1
        assert 'already has QM' in split.skipped[0].reason

    def test_reaction_without_a_transition_state_is_skipped(self):
        """A path reaction declaring no TS has nothing to compute."""
        network = _network([PDepPathReaction(label='r1', reactants=('A',), products=('B',),
                                             transition_state=None, kinetics_type='Arrhenius',
                                             kinetics_comment='family: 1,2_Insertion_CO')])
        split = split_qm_candidates(network, computed_ts_labels=frozenset())
        assert split.candidates == ()
        assert 'declares no transition state' in split.skipped[0].reason

    def test_mixed_network_splits_correctly(self):
        network = _network([_rxn('r1', 'family: 1,2_Insertion_CO'),
                            _rxn('r2', 'in family R_Recombination.'),
                            _rxn('r3', 'family: 1,3_Insertion_CO2'),
                            _rxn('r4', 'family: 1,2_Insertion_CO')])
        split = split_qm_candidates(network, computed_ts_labels=frozenset({'TS_r4'}))
        assert {c.ts_label for c in split.candidates} == {'TS_r1', 'TS_r3'}
        assert {s.label for s in split.skipped} == {'r2', 'r4'}

    def test_empty_network_yields_empty_split(self):
        split = split_qm_candidates(_network([]), computed_ts_labels=frozenset())
        assert split == CandidateSplit(candidates=(), skipped=())

    def test_candidate_order_is_deterministic(self):
        """Two runs over the same network must queue in the same order, or the budget admits a
        different subset each time and results stop being reproducible."""
        rxns = [_rxn(f'r{i}', 'family: 1,2_Insertion_CO') for i in range(6)]
        first = split_qm_candidates(_network(rxns), computed_ts_labels=frozenset())
        second = split_qm_candidates(_network(list(rxns)), computed_ts_labels=frozenset())
        assert [c.ts_label for c in first.candidates] == [c.ts_label for c in second.candidates]


class TestRoundLayout(object):
    """Each round is fully self-contained on disk."""

    def test_paths_are_under_the_round_root(self):
        paths = round_paths('/proj', 0)
        assert paths.root == os.path.join('/proj', 'round_0')
        for attr in ('arc_project', 'explorer_output', 'capture', 'hybrid'):
            assert getattr(paths, attr).startswith(paths.root + os.sep)

    def test_each_round_gets_its_own_arc_project(self):
        """This is what sidesteps capture's single-shot window: ARC wipes calcs/statmech/ on each
        rate pass, so a second ARC run must not share a project directory with the first."""
        assert round_paths('/proj', 0).arc_project != round_paths('/proj', 1).arc_project

    def test_capture_is_not_inside_the_arc_project(self):
        """capture_ts_artifacts refuses a capture_dir that resolves inside the ARC project."""
        paths = round_paths('/proj', 0)
        assert not paths.capture.startswith(paths.arc_project + os.sep)

    def test_diagram_sits_at_the_round_root(self):
        paths = round_paths('/proj', 2)
        assert paths.diagram == os.path.join('/proj', 'round_2', PES_LOOP_DIAGRAM_FILENAME)

    def test_negative_round_refused(self):
        with pytest.raises(ValueError, match='non-negative'):
            round_paths('/proj', -1)

    def test_relative_project_directory_refused(self):
        """An absolute root is what makes the confinement checks downstream meaningful."""
        with pytest.raises(ValueError, match='absolute'):
            round_paths('proj', 0)


class TestAttachSensitivityEvidence(object):
    """attach_sensitivity_evidence is what stands between the loop and a fabricated coefficient:
    real evidence is stamped, absent evidence removes the candidate (with its reason), and no
    default ever fills the gap."""

    def test_evidence_is_stamped_onto_candidates_in_order(self):
        network = _network([_rxn('r1', 'family: 1,2_Insertion_CO'),
                            _rxn('r2', 'family: 1,3_Insertion_CO2')])
        split = split_qm_candidates(network, computed_ts_labels=frozenset())
        stamped = attach_sensitivity_evidence(
            split, {'TS_r1': (-2.0e-5, 0.17), 'TS_r2': (1.0e-6, 0.008)})
        assert [c.ts_label for c in stamped.candidates] == ['TS_r1', 'TS_r2']
        assert stamped.candidates[0].coefficient == -2.0e-5
        assert stamped.candidates[0].delta_ln_k == 0.17
        assert stamped.candidates[1].coefficient == 1.0e-6
        assert stamped.skipped == ()

    def test_candidate_without_evidence_is_skipped_with_its_reason_never_defaulted(self):
        network = _network([_rxn('r1', 'family: 1,2_Insertion_CO'),
                            _rxn('r2', 'family: 1,3_Insertion_CO2')])
        split = split_qm_candidates(network, computed_ts_labels=frozenset())
        stamped = attach_sensitivity_evidence(split, {'TS_r2': (1.0e-6, 0.008)})
        assert [c.ts_label for c in stamped.candidates] == ['TS_r2']
        assert len(stamped.skipped) == 1
        assert stamped.skipped[0].label == 'r1'
        assert 'no finite sensitivity evidence' in stamped.skipped[0].reason
        assert 'inventing' in stamped.skipped[0].reason

    def test_prior_skips_are_preserved(self):
        network = _network([_rxn('r1', 'in family R_Recombination.'),
                            _rxn('r2', 'family: 1,2_Insertion_CO')])
        split = split_qm_candidates(network, computed_ts_labels=frozenset())
        stamped = attach_sensitivity_evidence(split, {'TS_r2': (1.0e-6, 0.008)})
        assert {s.label for s in stamped.skipped} == {'r1'}
        assert [c.ts_label for c in stamped.candidates] == ['TS_r2']

    def test_sa_dir_is_part_of_the_round_layout(self):
        paths = round_paths('/proj', 1)
        assert paths.sa == os.path.join('/proj', 'round_1', 'SA')
        assert paths.sa.startswith(paths.root + os.sep)

    def test_a_measured_response_below_the_floor_is_skipped_with_its_value(self):
        """A measured zero (or any below-floor response) is real data saying 'not worth the
        spend' -- the candidate is skipped with the measured value in its reason, never queued
        and never defaulted past the floor."""
        network = _network([_rxn('r1', 'family: 1,2_Insertion_CO'),
                            _rxn('r2', 'family: 1,3_Insertion_CO2')])
        split = split_qm_candidates(network, computed_ts_labels=frozenset())
        stamped = attach_sensitivity_evidence(
            split, {'TS_r1': (0.0, 0.0), 'TS_r2': (-2.0e-5, 0.17)}, min_delta_ln_k=1e-3)
        assert [c.ts_label for c in stamped.candidates] == ['TS_r2']
        assert len(stamped.skipped) == 1
        assert stamped.skipped[0].label == 'r1'
        assert 'below the min_delta_ln_k floor' in stamped.skipped[0].reason
        assert '0.000e+00' in stamped.skipped[0].reason

    def test_the_floor_defaults_off_and_an_at_floor_response_passes(self):
        network = _network([_rxn('r1', 'family: 1,2_Insertion_CO')])
        split = split_qm_candidates(network, computed_ts_labels=frozenset())
        assert attach_sensitivity_evidence(split, {'TS_r1': (0.0, 0.0)}).candidates
        at_floor = attach_sensitivity_evidence(split, {'TS_r1': (1.2e-7, 1e-3)},
                                               min_delta_ln_k=1e-3)
        assert [c.ts_label for c in at_floor.candidates] == ['TS_r1']
