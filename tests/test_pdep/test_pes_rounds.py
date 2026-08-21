"""Tests for t3.pdep.pes_rounds."""

import dataclasses
import os

import pytest

from arc.molecule.molecule import Molecule

from t3.pdep.parser import PDepNetwork, PDepPathReaction
from t3.pdep.pes_rounds import (CandidateSplit, PES_LOOP_DIAGRAM_FILENAME,
                                adoption_channel_keys_by_ts_label, attach_sensitivity_evidence,
                                channel_keys_by_ts_label, round_paths, split_qm_candidates,
                                structural_channel_key)


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


class TestStructuralChannelKeys(object):
    """The identity QM artifacts are carried across network-file writes on. Labels cannot be it:
    every TS label Arkane writes is purely positional (rmgpy/rmg/pdep.py:854 replaces every path
    reaction's transition state with a fresh, label-less object before each write), and reaction
    labels are positional too."""

    def test_keys_are_structural_and_direction_insensitive(self):
        adj_a = Molecule().from_smiles('C').to_adjacency_list()
        adj_b = Molecule().from_smiles('CC').to_adjacency_list()
        network = _network([_rxn('r1', 'family: 1,2_Insertion_CO')])
        network = dataclasses.replace(network, species_structures={'A_r1': adj_a, 'B_r1': adj_b})
        forward = _rxn('f', 'family: x')
        forward = dataclasses.replace(forward, reactants=('A_r1',), products=('B_r1',))
        backward = dataclasses.replace(forward, reactants=('B_r1',), products=('A_r1',))
        key_f = structural_channel_key(forward, network.species_structures)
        key_b = structural_channel_key(backward, network.species_structures)
        assert key_f == key_b == (('CC|m1',), ('C|m1',))

    def test_key_is_label_and_atom_order_independent(self):
        """The same molecule written under a different species label and a permuted atom order
        keys identically -- this is what makes cross-project adoption structural."""
        co2_a = ('1 O u0 p2 c0 {2,D}\n'
                 '2 C u0 p0 c0 {1,D} {3,D}\n'
                 '3 O u0 p2 c0 {2,D}\n')
        co2_b = ('1 C u0 p0 c0 {2,D} {3,D}\n'
                 '2 O u0 p2 c0 {1,D}\n'
                 '3 O u0 p2 c0 {1,D}\n')
        water = Molecule().from_smiles('O').to_adjacency_list()
        rxn_a = dataclasses.replace(_rxn('a', 'family: x'), reactants=('CO2(11)',),
                                    products=('W(1)',))
        rxn_b = dataclasses.replace(_rxn('b', 'family: x'), reactants=('carbon dioxide(999)',),
                                    products=('W(2)',))
        key_a = structural_channel_key(rxn_a, {'CO2(11)': co2_a, 'W(1)': water})
        key_b = structural_channel_key(rxn_b, {'carbon dioxide(999)': co2_b, 'W(2)': water})
        assert key_a == key_b

    def test_missing_or_unparseable_structure_yields_no_key(self):
        rxn = dataclasses.replace(_rxn('r', 'family: x'), reactants=('A',), products=('B',))
        assert structural_channel_key(rxn, {'A': Molecule().from_smiles('C').to_adjacency_list()}) \
            is None
        assert structural_channel_key(
            rxn, {'A': 'not an adjacency list', 'B': 'nope'}) is None

    def test_duplicate_channels_are_both_omitted_from_the_map(self):
        """Two TSs connecting structurally identical channels cannot be told apart by this key;
        carrying either would risk the wrong-saddle-point misattribution, so both are omitted."""
        adj_a = Molecule().from_smiles('C').to_adjacency_list()
        adj_b = Molecule().from_smiles('CC').to_adjacency_list()
        adj_c = Molecule().from_smiles('CCC').to_adjacency_list()
        r1 = dataclasses.replace(_rxn('r1', 'family: x', ts='TS1'), reactants=('A',),
                                 products=('B',))
        r2 = dataclasses.replace(_rxn('r2', 'family: x', ts='TS2'), reactants=('B',),
                                 products=('A',))
        r3 = dataclasses.replace(_rxn('r3', 'family: x', ts='TS3'), reactants=('B',),
                                 products=('C',))
        network = dataclasses.replace(
            _network([r1, r2, r3]),
            species_structures={'A': adj_a, 'B': adj_b, 'C': adj_c})
        keys = channel_keys_by_ts_label(network)
        assert set(keys) == {'TS3'}

    def test_a_ts_shared_by_several_path_reactions_is_omitted(self):
        adj_a = Molecule().from_smiles('C').to_adjacency_list()
        adj_b = Molecule().from_smiles('CC').to_adjacency_list()
        adj_c = Molecule().from_smiles('CCC').to_adjacency_list()
        r1 = dataclasses.replace(_rxn('r1', 'family: x', ts='TS1'), reactants=('A',),
                                 products=('B',))
        r2 = dataclasses.replace(_rxn('r2', 'family: x', ts='TS1'), reactants=('B',),
                                 products=('C',))
        network = dataclasses.replace(
            _network([r1, r2]), species_structures={'A': adj_a, 'B': adj_b, 'C': adj_c})
        assert channel_keys_by_ts_label(network) == {}

    def test_spin_states_do_not_collapse_onto_one_key(self):
        """BLOCKING regression pin: SMILES does not encode spin state -- singlet and triplet CH2
        both render '[CH2]' -- so without the explicit multiplicity suffix a prior
        CH2(S) + CO <=> CH2CO barrier would be adopted onto the CH2(T) channel: a FALSE structural
        match, i.e. exactly the wrong-saddle-point misattribution structural keying exists to
        prevent."""
        ch2_triplet = ('multiplicity 2\n'
                       '1 C u2 p0 c0 {2,S} {3,S}\n'
                       '2 H u0 p0 c0 {1,S}\n'
                       '3 H u0 p0 c0 {1,S}\n').replace('multiplicity 2', 'multiplicity 3')
        ch2_singlet = ('1 C u0 p1 c0 {2,S} {3,S}\n'
                       '2 H u0 p0 c0 {1,S}\n'
                       '3 H u0 p0 c0 {1,S}\n')
        co = Molecule().from_smiles('[C-]#[O+]').to_adjacency_list()
        ch2co = Molecule().from_smiles('C=C=O').to_adjacency_list()
        rxn_triplet = dataclasses.replace(_rxn('t', 'family: x'), reactants=('CH2(T)', 'CO'),
                                          products=('CH2CO',))
        rxn_singlet = dataclasses.replace(_rxn('s', 'family: x'), reactants=('CH2(S)', 'CO'),
                                          products=('CH2CO',))
        key_triplet = structural_channel_key(
            rxn_triplet, {'CH2(T)': ch2_triplet, 'CO': co, 'CH2CO': ch2co})
        key_singlet = structural_channel_key(
            rxn_singlet, {'CH2(S)': ch2_singlet, 'CO': co, 'CH2CO': ch2co})
        assert key_triplet is not None and key_singlet is not None
        assert key_triplet != key_singlet


class TestAdoptionChannelKeysAreNotEndpointsOnly(object):
    """The endpoints-only structural key identifies reactants and products, not a saddle point.
    channel_keys_by_ts_label's duplicate guard covers two A<=>B pathways sitting in ONE file, but
    cannot reach across two files: a prior capture holding one A<=>B pathway and this run holding
    a DIFFERENT one key identically, no duplicate is present anywhere, and the prior artifact is
    attached to a saddle point it was never computed for."""

    @staticmethod
    def _network(comments, network_id='network1_1'):
        adjlists = {'A': Molecule().from_smiles('CC').to_adjacency_list(),
                    'B': Molecule().from_smiles('C').to_adjacency_list()}
        path_reactions = tuple(
            PDepPathReaction(label=f'reaction{i}', reactants=('A',), products=('B',),
                             transition_state=f'TS{i}', kinetics_type='Arrhenius',
                             kinetics_comment=comment)
            for i, comment in enumerate(comments))
        return PDepNetwork(network_id=network_id, path=f'/abs/{network_id}.py', label=None,
                           species_labels=('A', 'B'),
                           transition_state_labels=tuple(f'TS{i}' for i in range(len(comments))),
                           path_reactions=path_reactions, isomers=(), reactant_channels=(),
                           product_channels=(), species_structures=adjlists)

    def test_two_different_pathways_between_the_same_endpoints_key_differently(self):
        """The whole defect in one assertion: endpoints-only, these two are indistinguishable."""
        prior = self._network(['family: H_Abstraction'])
        current = self._network(['family: Disproportionation'])
        assert channel_keys_by_ts_label(prior)['TS0'] == channel_keys_by_ts_label(current)['TS0'], \
            'the endpoints-only key cannot tell these apart -- that is why it is not what ' \
            'adoption matches on'
        assert adoption_channel_keys_by_ts_label(prior)['TS0'] \
            != adoption_channel_keys_by_ts_label(current)['TS0']

    def test_a_channel_with_no_discriminator_at_all_is_refused(self):
        """A refused match costs quantum chemistry that was already paid for once; a false one
        costs the correctness of the barrier. No discriminator, no adoption."""
        network = self._network([''])
        assert 'TS0' in channel_keys_by_ts_label(network)
        assert adoption_channel_keys_by_ts_label(network) == {}

    def test_a_library_sourced_channel_falls_back_to_its_kinetics_comment(self):
        """network21_1's three path reactions all come from a reaction library and name no family
        at all. Refusing them outright would make reuse dead for a whole class of real network --
        the comment names the library, which is a pathway identity the endpoints are not."""
        library = self._network(["Reaction library: 'primaryNitrogenLibrary'"])
        node = self._network(['Estimated from node Root_N-1R!H->C'])
        assert 'TS0' in adoption_channel_keys_by_ts_label(library)
        assert adoption_channel_keys_by_ts_label(library)['TS0'] \
            != adoption_channel_keys_by_ts_label(node)['TS0']

    def test_the_family_wins_over_the_comment_when_one_is_named(self):
        """The family is the coarser and more stable of the two -- it survives a degeneracy
        multiplier or rate-rule retraining that rewords the comment around it, so a family-bearing
        channel stays adoptable across RMG database versions."""
        a = self._network(['Estimated using template [x] for rate rule [y]\n'
                           'Multiplied by reaction path degeneracy 2.0\n'
                           'family: H_Abstraction'])
        b = self._network(['Estimated using template [x] for rate rule [y]\n'
                           'Multiplied by reaction path degeneracy 3.0\n'
                           'family: H_Abstraction'])
        assert adoption_channel_keys_by_ts_label(a)['TS0'] \
            == adoption_channel_keys_by_ts_label(b)['TS0']

    def test_the_same_pathway_still_matches_across_files(self):
        """The refusal must not be so broad that reuse can never fire at all."""
        prior = self._network(['family: H_Abstraction'], network_id='network9_9')
        current = self._network(['family: H_Abstraction'])
        assert adoption_channel_keys_by_ts_label(prior)['TS0'] \
            == adoption_channel_keys_by_ts_label(current)['TS0']

    def test_a_structurally_unkeyable_channel_is_still_refused(self):
        """The family is an ADDITIONAL discriminator, never a replacement: a channel with no
        parseable structure must not become adoptable just because it names a family."""
        network = self._network(['family: H_Abstraction'])
        network = dataclasses.replace(network, species_structures={'A': 'not an adjacency list'})
        assert adoption_channel_keys_by_ts_label(network) == {}
