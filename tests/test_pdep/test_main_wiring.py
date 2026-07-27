#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_main_wiring module

Acceptance test for T3.determine_species_from_pdep_network()'s in-run wiring: both cache paths
(reuse a valid cached SA, and generate + persist a new one), the cache-rejected path (an SA
output present without T3's own sidecar must not be trusted), and a cross-check that the wiring's
`qualified`/`network_id` verdict agrees with a direct, independent call into
`t3.pdep.selector.select_from_sa_dict()` for the same inputs.

Real Arkane is never invoked here: `t3.main.run_arkane_job` is monkeypatched throughout, and a
synthetic Arkane sensitivity dictionary stands in for real output. The real fixture network file
(`tests/data/pdep_network/iteration_1/RMG/pdep/network4_2.py`) is used, but the whole
`iteration_1` fixture tree is copied into `tmp_path` first, so nothing under `tests/data/` is
ever written to.
"""

import os

import pytest

import t3.main as t3_main
from t3.pdep.cache import sa_cache_metadata_path, write_sa_cache_metadata
from t3.pdep.join import (JOIN_STATUS_ALREADY_PRESENT,
                          JOIN_STATUS_NOT_QUEUED,
                          JOIN_STATUS_QUEUED,
                          TSJoinRecord,
                          arc_ts_label,
                          read_ts_join_sidecar,
                          ts_join_sidecar_path,
                          write_ts_join_sidecar,
                          )
from t3.pdep.parser import parse_pdep_network_file, parse_pdep_network_text
from t3.pdep.selector import (CACHE_STATUS_CACHED_VALID,
                              CACHE_STATUS_GENERATED,
                              PDepNetworkSelection,
                              SensitiveTransitionState,
                              select_from_sa_dict,
                              )
from arc.common import save_yaml_file

from tests.test_pdep._wiring_helpers import (CONDITION,
                                             NETWORK_NAME,
                                             NETWORK_REACTION_STR,
                                             build_pdep_rxns_to_explore as _build_pdep_rxns_to_explore,
                                             build_sa_dict as _build_sa_dict,
                                             build_t3 as _build_t3,
                                             candidate_sa_path as _candidate_sa_path,
                                             network_path as _network_path,
                                             )


# Adjacency lists for the two TS1 species the sensitivity fixture does not carry, generated with
# arc.molecule from the SMILES that network4_2.py itself declares for them
# (CH2(S)(53) -> '[CH2]', C3rad(4) -> '[CH2]CC'), so the structures are the file's own, not invented.
# With these present, TS1 (CH2(S)(53) + C3rad(4) <=> C4rad(5)) becomes queueable, which is what lets
# the in-run path be exercised end to end rather than only on its refusal branch.
TS1_STRUCTURES = {
    'CH2(S)(53)': 'multiplicity 3\n1 C u2 p0 c0 {2,S} {3,S}\n2 H u0 p0 c0 {1,S}\n3 H u0 p0 c0 {1,S}\n',
    'C3rad(4)': 'multiplicity 2\n'
                '1  C u1 p0 c0 {2,S} {4,S} {5,S}\n'
                '2  C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}\n'
                '3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}\n'
                '4  H u0 p0 c0 {1,S}\n5  H u0 p0 c0 {1,S}\n'
                '6  H u0 p0 c0 {2,S}\n7  H u0 p0 c0 {2,S}\n'
                '8  H u0 p0 c0 {3,S}\n9  H u0 p0 c0 {3,S}\n10 H u0 p0 c0 {3,S}\n',
}


def _sa_dict_with_ts1_structures(t3) -> dict:
    """The sensitivity fixture, extended so TS1's species can actually be built."""
    sa_dict = _build_sa_dict(t3)
    sa_dict['structures'] = {**sa_dict['structures'], **TS1_STRUCTURES}
    return sa_dict


def _selection(ts_labels: list) -> PDepNetworkSelection:
    """A qualified selection naming the given transition states as the uncertain ones.

    Built directly rather than via the selector so that a transition state can be put under test
    independently of whether the fixture's own sensitivity data happens to select it.
    """
    entries = [SensitiveTransitionState(ts_label=ts_label,
                                        coefficient=0.05,
                                        condition=CONDITION,
                                        path_reaction_label=None,
                                        path_reaction_str=None,
                                        kinetics_comment='',
                                        uncertain=True,
                                        delta_ln_k=0.05,
                                        ) for ts_label in ts_labels]
    return PDepNetworkSelection(network_id=NETWORK_NAME,
                                qualified=True,
                                selected_ts=list(entries),
                                uncertain_path_reactions=list(entries),
                                )


class TestDetermineSpeciesFromPdepNetworkWiring(object):

    def test_cached_valid_path_never_invokes_arkane(self, tmp_path, monkeypatch):
        """A valid, T3-vouched-for cache must be reused; Arkane must never be invoked."""
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_path = _candidate_sa_path(t3, method='CSE')
        os.makedirs(os.path.dirname(sa_path), exist_ok=True)
        save_yaml_file(sa_path, _build_sa_dict(t3))
        write_sa_cache_metadata(sa_path=sa_path,
                                network_path=_network_path(t3),
                                network_id=NETWORK_NAME,
                                method='CSE',
                                )

        def _fail_if_called(*args, **kwargs):
            pytest.fail('run_arkane_job should not be invoked when a valid cache is present.')

        monkeypatch.setattr(t3_main, 'run_arkane_job', _fail_if_called)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert len(t3.pdep_network_selections) == 1
        selection = t3.pdep_network_selections[0]
        assert selection.network_id == NETWORK_NAME
        assert selection.cache_status == CACHE_STATUS_CACHED_VALID
        assert selection.qualified is True

    def test_generate_path_invokes_arkane_and_persists_sidecar(self, tmp_path, monkeypatch):
        """No cache present: Arkane must be invoked, and a sidecar must be written on success."""
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_dict = _build_sa_dict(t3)
        calls = list()

        def _fake_run_arkane_job(input_file, output_directory, plot=True, logger=None, **kwargs):
            calls.append(output_directory)
            sensitivity_dir = os.path.join(output_directory, 'sensitivity')
            os.makedirs(sensitivity_dir, exist_ok=True)
            save_yaml_file(os.path.join(sensitivity_dir, 'sa_coefficients.yml'), sa_dict)
            return True

        monkeypatch.setattr(t3_main, 'run_arkane_job', _fake_run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert len(calls) == 1, 'run_arkane_job should have been invoked exactly once (CSE, the first ME method).'
        assert len(t3.pdep_network_selections) == 1
        selection = t3.pdep_network_selections[0]
        assert selection.network_id == NETWORK_NAME
        assert selection.cache_status == CACHE_STATUS_GENERATED
        assert selection.qualified is True

        sa_path = _candidate_sa_path(t3, method='CSE')
        assert os.path.isfile(sa_path)
        assert os.path.isfile(sa_cache_metadata_path(sa_path)), \
            'A T3 sidecar should have been written alongside the freshly generated SA output.'

    def test_sa_output_without_sidecar_is_rejected_and_arkane_is_rerun(self, tmp_path, monkeypatch):
        """An SA output present without T3's own sidecar must not be trusted; Arkane re-runs."""
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_dict = _build_sa_dict(t3)

        # Pre-existing SA output from some prior, untracked run: no t3_sa_cache.yml sidecar.
        stale_sa_path = _candidate_sa_path(t3, method='CSE')
        os.makedirs(os.path.dirname(stale_sa_path), exist_ok=True)
        save_yaml_file(stale_sa_path, {'structures': {}, 'stale reaction': {}})
        assert not os.path.isfile(sa_cache_metadata_path(stale_sa_path))

        calls = list()

        def _fake_run_arkane_job(input_file, output_directory, plot=True, logger=None, **kwargs):
            calls.append(output_directory)
            sensitivity_dir = os.path.join(output_directory, 'sensitivity')
            os.makedirs(sensitivity_dir, exist_ok=True)
            save_yaml_file(os.path.join(sensitivity_dir, 'sa_coefficients.yml'), sa_dict)
            return True

        monkeypatch.setattr(t3_main, 'run_arkane_job', _fake_run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert len(calls) == 1, 'Arkane must be re-run when an SA output lacks a T3 sidecar.'
        selection = t3.pdep_network_selections[0]
        assert selection.cache_status == CACHE_STATUS_GENERATED
        assert selection.qualified is True

    def test_qualified_and_network_id_agree_with_direct_selector_call(self, tmp_path, monkeypatch):
        """The wiring's decision must agree with select_from_sa_dict() called directly."""
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_dict = _build_sa_dict(t3)
        sa_path = _candidate_sa_path(t3, method='CSE')
        os.makedirs(os.path.dirname(sa_path), exist_ok=True)
        save_yaml_file(sa_path, sa_dict)
        write_sa_cache_metadata(sa_path=sa_path,
                                network_path=_network_path(t3),
                                network_id=NETWORK_NAME,
                                method='CSE',
                                )

        def _fail_if_called(*args, **kwargs):
            pytest.fail('run_arkane_job should not be invoked when a valid cache is present.')

        monkeypatch.setattr(t3_main, 'run_arkane_job', _fail_if_called)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)
        selection = t3.pdep_network_selections[0]

        network = parse_pdep_network_file(_network_path(t3))
        direct_selection = select_from_sa_dict(
            sa_dict=sa_dict,
            network=network,
            network_reaction=NETWORK_REACTION_STR,
            relative_threshold=t3.t3['sensitivity']['pdep_SA_threshold'],
            min_delta_ln_k=t3.t3['sensitivity']['pdep_min_delta_ln_k'],
        )

        assert selection.network_id == direct_selection.network_id == NETWORK_NAME
        assert selection.qualified == direct_selection.qualified is True


class TestQueuePdepTransitionStates(object):
    """Queueing a qualified network's uncertain transition states to ARC, and recording the join.

    ``TS2`` of ``network4_2`` is used for the queued cases because all three of its species
    (``H(34) + C4ene(26) <=> C4rad(5)``) have adjacency lists in the sensitivity fixture. Whether the
    fixture's own sensitivity data would select TS2 is a separate question, answered by the selector
    and not by this code, so the selection is built directly here.
    """

    def test_a_queued_transition_state_carries_the_chosen_arc_label_onto_the_qm_queue(self, tmp_path):
        """Test the whole point of the design: the label T3 chooses survives onto what ARC receives.

        ``add_reaction`` puts ``reaction.copy()`` on the QM queue, and that copy is a plain
        ``ARCReaction`` -- every T3-only field, ``t3_index`` included, is dropped by it. ``ts_label``
        is the one identity that survives, which is why the join is built on it. If this assertion
        ever fails, the join is broken at its root, not merely mislabeled.
        """
        t3 = _build_t3(tmp_path)
        network = parse_pdep_network_file(_network_path(t3))
        structures = _build_sa_dict(t3)['structures']

        records = t3.queue_pdep_transition_states(network=network,
                                                  selection=_selection(['TS2']),
                                                  structures=structures,
                                                  network_name=NETWORK_NAME,
                                                  )

        assert len(records) == 1
        assert records[0].status == JOIN_STATUS_QUEUED
        assert records[0].arc_ts_label == arc_ts_label(NETWORK_NAME, 'TS2') == 'T3PDep_network4_2_TS2'
        assert records[0].t3_reaction_key is not None
        assert t3.qm['reactions'], 'the reaction should have been queued for ARC'
        assert t3.qm['reactions'][-1].ts_label == 'T3PDep_network4_2_TS2'

    def test_the_queued_label_is_outside_arcs_own_index_namespace(self, tmp_path):
        """Test that ARC will not hand the same label to some other reaction it had to name itself."""
        t3 = _build_t3(tmp_path)
        network = parse_pdep_network_file(_network_path(t3))
        t3.queue_pdep_transition_states(network=network,
                                        selection=_selection(['TS2']),
                                        structures=_build_sa_dict(t3)['structures'],
                                        network_name=NETWORK_NAME,
                                        )
        label = t3.qm['reactions'][-1].ts_label
        assert not (label.startswith('TS') and label[2:].isdigit())

    def test_a_transition_state_whose_species_have_no_structures_is_recorded_not_dropped(self, tmp_path):
        """Test that an unqueueable transition state still produces a record.

        ``TS1`` is ``CH2(S)(53) + C3rad(4) <=> C4rad(5)``, and the sensitivity fixture carries no
        adjacency lists for the first two, so it cannot be sent to ARC. Silently dropping it would
        later be indistinguishable from a transition state whose quantum chemistry failed -- and
        those two call for opposite responses.
        """
        t3 = _build_t3(tmp_path)
        network = parse_pdep_network_file(_network_path(t3))

        records = t3.queue_pdep_transition_states(network=network,
                                                  selection=_selection(['TS1']),
                                                  structures=_build_sa_dict(t3)['structures'],
                                                  network_name=NETWORK_NAME,
                                                  )

        assert len(records) == 1
        assert records[0].status == JOIN_STATUS_NOT_QUEUED
        assert 'adjacency list' in records[0].reason
        assert records[0].expected_artifact_path is None
        assert not t3.qm['reactions'], 'nothing should have been queued for ARC'

    def test_an_already_known_reaction_is_recorded_as_such_with_no_expected_artifact(self, tmp_path):
        """Test the case that would otherwise fail silently.

        ``add_reaction`` returns ``None`` for a reaction T3 already knows and does NOT re-queue it,
        so no artifact will appear for it in this iteration. Claiming an expected artifact path here
        would make a never-computed transition state look merely unfinished.
        """
        t3 = _build_t3(tmp_path)
        network = parse_pdep_network_file(_network_path(t3))
        structures = _build_sa_dict(t3)['structures']
        kwargs = dict(network=network, structures=structures, network_name=NETWORK_NAME)

        first = t3.queue_pdep_transition_states(selection=_selection(['TS2']), **kwargs)
        t3.pdep_ts_join_records = list()  # a later iteration starts with an empty record set
        second = t3.queue_pdep_transition_states(selection=_selection(['TS2']), **kwargs)

        assert first[0].status == JOIN_STATUS_QUEUED
        assert second[0].status == JOIN_STATUS_ALREADY_PRESENT
        assert second[0].expected_artifact_path is None
        assert second[0].t3_reaction_key == first[0].t3_reaction_key

    def test_offering_the_same_network_twice_in_one_iteration_is_idempotent(self, tmp_path, monkeypatch):
        """Test the common path: several sensitive PDep reactions belong to the same network.

        ``determine_species_from_pdep_network`` calls this method once per sensitive reaction, and
        several such reactions can point at the same network, so the same ``ts_label`` is legitimately
        offered again within one iteration without ``self.pdep_ts_join_records`` having been reset in
        between (it is only reset at entry to ``determine_species_from_pdep_network``, not between
        these calls). Re-deciding an already-recorded transition state would re-invoke
        ``add_reaction`` for something already settled and would collide with its own record in
        ``merge_ts_join_records``, so the second call must recognize the key and skip it.
        """
        t3 = _build_t3(tmp_path)
        network = parse_pdep_network_file(_network_path(t3))
        structures = _build_sa_dict(t3)['structures']
        kwargs = dict(network=network, structures=structures, network_name=NETWORK_NAME)

        reaction_count_before_second_call = None

        first = t3.queue_pdep_transition_states(selection=_selection(['TS2']), **kwargs)
        reaction_count_before_second_call = len(t3.reactions)
        second = t3.queue_pdep_transition_states(selection=_selection(['TS2']), **kwargs)

        matching = [record for record in t3.pdep_ts_join_records if record.key == (NETWORK_NAME, 'TS2')]
        assert len(matching) == 1, 'the second call must not add a conflicting record for the same key'
        assert matching[0].status == JOIN_STATUS_QUEUED
        assert matching[0].t3_reaction_key == first[0].t3_reaction_key
        assert len(t3.reactions) == reaction_count_before_second_call, \
            'add_reaction should not be invoked again for an already-recorded transition state'
        assert second == list(), \
            'the second call should not append a duplicate record to its own return value either'

    def test_an_unsafe_network_name_is_recorded_not_silently_dropped(self, tmp_path):
        """Test that a network name `arc_ts_label` refuses still produces a record.

        Every selected transition state must produce a record, including this one: a network name
        containing a character ARC would rewrite (here, a space) makes `arc_ts_label` raise
        `ValueError`, and silently skipping it would break the documented invariant that a selected
        transition state is never dropped without a trace.
        """
        t3 = _build_t3(tmp_path)
        network = parse_pdep_network_file(_network_path(t3))
        structures = _build_sa_dict(t3)['structures']
        unsafe_network_name = 'network 4_2'

        records = t3.queue_pdep_transition_states(network=network,
                                                  selection=_selection(['TS2']),
                                                  structures=structures,
                                                  network_name=unsafe_network_name,
                                                  )

        assert len(records) == 1
        assert records[0].status == JOIN_STATUS_NOT_QUEUED
        assert records[0].arc_ts_label is None
        assert 'network_id' in records[0].reason
        assert not t3.qm['reactions'], 'nothing should have been queued for ARC'

    def test_a_species_with_a_malformed_adjlist_is_recorded_not_raised(self, tmp_path):
        """Test that a malformed adjacency list is recorded rather than crashing the whole run.

        The sensitivity output's ``structures`` mapping is generated by a separate tool (Arkane), so
        a corrupt adjacency list for one species is an external-data problem, not a T3 bug -- it must
        be handled the same way as the already-handled missing-structure case (log and drop this
        transition state's record), not propagate up through `determine_species_from_pdep_network`
        and lose every join record collected so far in this iteration.
        """
        t3 = _build_t3(tmp_path)
        network = parse_pdep_network_file(_network_path(t3))
        structures = dict(_build_sa_dict(t3)['structures'])
        structures['H(34)'] = 'this is not a valid adjacency list at all !!!'

        records = t3.queue_pdep_transition_states(network=network,
                                                  selection=_selection(['TS2']),
                                                  structures=structures,
                                                  network_name=NETWORK_NAME,
                                                  )

        assert len(records) == 1
        assert records[0].status == JOIN_STATUS_NOT_QUEUED
        assert 'H(34)' in records[0].reason
        assert not t3.qm['reactions'], 'nothing should have been queued for ARC'

    def test_a_shared_transition_state_is_refused_rather_than_guessed(self, tmp_path):
        """Test that a transition state owning several path reactions is not queued.

        There is no basis for picking one of them: computing a transition state for ``A <=> B`` and
        then using it for ``C <=> D`` is silently wrong chemistry, not an approximation.
        """
        t3 = _build_t3(tmp_path)
        structures = _build_sa_dict(t3)['structures']
        text = ("species(label='H(34)')\n"
                "species(label='C4ene(26)')\n"
                "species(label='C4rad(5)')\n"
                "transitionState(label='TS1')\n"
                "reaction(label='reaction1', reactants=['H(34)','C4ene(26)'], products=['C4rad(5)'], "
                "transitionState='TS1')\n"
                "reaction(label='reaction2', reactants=['C4rad(5)'], products=['H(34)','C4ene(26)'], "
                "transitionState='TS1')\n"
                "network(label='n', isomers=['C4rad(5)'], reactants=[['H(34)','C4ene(26)']])\n")
        network = parse_pdep_network_text(text=text, network_id=NETWORK_NAME)

        records = t3.queue_pdep_transition_states(network=network,
                                                  selection=_selection(['TS1']),
                                                  structures=structures,
                                                  network_name=NETWORK_NAME,
                                                  )

        assert len(records) == 1
        assert records[0].status == JOIN_STATUS_NOT_QUEUED
        assert 'ambiguous' in records[0].reason
        assert len(records[0].path_reaction_strs) == 2
        assert not t3.qm['reactions']

    def test_the_join_sidecar_is_written_by_the_in_run_path(self, tmp_path, monkeypatch):
        """Test that a real in-run pass writes the sidecar into the ARC project directory."""
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_path = _candidate_sa_path(t3, method='CSE')
        os.makedirs(os.path.dirname(sa_path), exist_ok=True)
        save_yaml_file(sa_path, _sa_dict_with_ts1_structures(t3))
        write_sa_cache_metadata(sa_path=sa_path,
                                network_path=_network_path(t3),
                                network_id=NETWORK_NAME,
                                method='CSE',
                                )
        monkeypatch.setattr(t3_main, 'run_arkane_job',
                            lambda *args, **kwargs: pytest.fail('Arkane should not run with a valid cache.'))

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        records = read_ts_join_sidecar(t3.paths['ARC'])
        assert records, 'the qualified network should have produced at least one join record'
        assert {record.network_ts_label for record in records} == {'TS1'}
        assert records == t3.pdep_ts_join_records
        # End to end: with TS1's structures available it is genuinely queued, and the reaction handed
        # to ARC carries the very label the sidecar promises its artifact will be filed under.
        assert records[0].status == JOIN_STATUS_QUEUED
        assert records[0].expected_artifact_path == os.path.join(
            t3.paths['ARC'], 'calcs', 'statmech', 'kinetics', 'TSs', 'T3PDep_network4_2_TS1.py')
        assert t3.qm['reactions'][-1].ts_label == records[0].arc_ts_label

    def test_a_stale_sidecar_is_removed_when_this_pass_selects_nothing(self, tmp_path, monkeypatch):
        """Test that an empty pass removes a stale sidecar left by an earlier, interrupted attempt.

        A restarted run of the same iteration can leave a `t3_pdep_ts_join.yml` sidecar on disk from a
        crashed/interrupted earlier attempt. If the new pass selects nothing, the stale sidecar must
        not survive to lie to the post-ARC discovery step about transition states this pass never
        even considered.
        """
        t3 = _build_t3(tmp_path)
        os.makedirs(t3.paths['ARC'], exist_ok=True)
        stale_records = [TSJoinRecord(network_id='network9_9',
                                      network_ts_label='TS7',
                                      arc_ts_label=arc_ts_label('network9_9', 'TS7'),
                                      status=JOIN_STATUS_QUEUED,
                                      reason='stale from a previous, interrupted attempt')]
        write_ts_join_sidecar(arc_project_directory=t3.paths['ARC'], records=stale_records)
        assert os.path.isfile(ts_join_sidecar_path(t3.paths['ARC']))

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=[])

        assert not t3.pdep_ts_join_records
        assert not os.path.isfile(ts_join_sidecar_path(t3.paths['ARC'])), \
            'a stale sidecar must not survive a pass that selects nothing'

    def test_stale_records_from_a_previous_iteration_are_discarded(self, tmp_path, monkeypatch):
        """Test that the record set is emptied at the start of each pass.

        The sidecar is written to a per-iteration ARC project directory, so a carried-over record
        would advertise a previous iteration's artifact path as if it belonged to this one. Seeding a
        stale record and checking it is gone tests the reset directly; comparing two identical passes
        would not, since the merge absorbs an identical repeat and would look correct either way.
        """
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_path = _candidate_sa_path(t3, method='CSE')
        os.makedirs(os.path.dirname(sa_path), exist_ok=True)
        save_yaml_file(sa_path, _sa_dict_with_ts1_structures(t3))
        write_sa_cache_metadata(sa_path=sa_path,
                                network_path=_network_path(t3),
                                network_id=NETWORK_NAME,
                                method='CSE',
                                )
        monkeypatch.setattr(t3_main, 'run_arkane_job', lambda *args, **kwargs: True)
        t3.pdep_ts_join_records = [TSJoinRecord(network_id='network9_9',
                                                network_ts_label='TS7',
                                                arc_ts_label=arc_ts_label('network9_9', 'TS7'),
                                                status=JOIN_STATUS_QUEUED,
                                                expected_artifact_path='/a/previous/iteration/TS7.py',
                                                )]

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert t3.pdep_ts_join_records, 'this pass should have produced its own records'
        assert 'network9_9' not in {record.network_id for record in t3.pdep_ts_join_records}
        assert all(record.network_id == NETWORK_NAME for record in t3.pdep_ts_join_records)
