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

import builtins
import os

import pytest

import t3.main as t3_main
from t3.chem import T3Species
from t3.pdep.cache import sa_cache_metadata_path, write_sa_cache_metadata
from t3.pdep.capture import CAPTURE_MANIFEST_FILE_NAME, capture_ts_artifacts, verify_capture
from t3.pdep.join import (JOIN_STATUS_ALREADY_PRESENT,
                          JOIN_STATUS_NOT_QUEUED,
                          JOIN_STATUS_QUEUED,
                          TSJoinRecord,
                          arc_ts_label,
                          expected_ts_artifact_path,
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
from arc.common import read_yaml_file, save_yaml_file

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


def _write_arc_info(t3, species=None, reactions=None) -> None:
    """Write a minimal ARC ``project_info.yml`` at ``t3.paths['ARC info']``.

    ``process_arc_run()`` raises if this file is absent, regardless of whether any PDep
    transition states were queued; empty ``species``/``reactions`` lists make its convergence
    bookkeeping loop a no-op, keeping these tests focused on the capture/marker wiring rather
    than on species/reaction convergence.
    """
    os.makedirs(os.path.dirname(t3.paths['ARC info']), exist_ok=True)
    save_yaml_file(path=t3.paths['ARC info'], content={'species': species or [], 'reactions': reactions or []})


def _mark_rmg_terminated(t3) -> None:
    """Write an ``RMG.log`` that ``check_rmg_status()`` reads as a converged RMG run."""
    os.makedirs(os.path.dirname(t3.paths['RMG log']), exist_ok=True)
    with open(t3.paths['RMG log'], 'w') as f:
        f.write('MODEL GENERATION COMPLETED\n')


def _mark_arc_terminated(t3) -> None:
    """Write an ``arc.log`` that ``check_arc_status()`` reads as a converged, terminated ARC run."""
    os.makedirs(os.path.dirname(t3.paths['ARC log']), exist_ok=True)
    with open(t3.paths['ARC log'], 'w') as f:
        f.write('ARC execution terminated on Sun Dec  4 11:50:29 2022\n')


def _write_ts_artifact(path: str) -> None:
    """Write a minimal artifact at ``path`` in the ARC statmech species-file shape, mirroring
    ``tests/test_pdep/test_capture.py::_write_artifact``."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    log_path = os.path.join(os.path.dirname(path), 'output.out')
    with open(log_path, 'w') as f:
        f.write('stub quantum chemistry log\n')
    content = """linear = False

spinMultiplicity = 2

energy = Log('output.out')

geometry = Log('output.out')

frequencies = Log('output.out')
"""
    with open(path, 'w') as f:
        f.write(content)


def _write_ts_status_yml(arc_dir: str, label: str, converged: bool) -> None:
    """Append a converged/unconverged ``status.yml`` entry for ``label``, mirroring
    ``tests/test_pdep/test_capture.py::_write_status_yml``."""
    output_dir = os.path.join(arc_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'status.yml'), 'a') as f:
        f.write(
            f"{label}:\n"
            f"  convergence: {str(converged).lower()}\n"
            "  job_types: {}\n"
            "  paths: {}\n"
            "  info: ''\n"
            "  errors: ''\n"
        )


def _queue_usable_ts(t3, network_id='network4_2', network_ts_label='TS9') -> TSJoinRecord:
    """Queue one usable PDep transition state against ``t3``'s current ARC directory: write the
    join sidecar record T3 would have written pre-ARC, plus the converged ARC artifact and status
    entry ARC would have produced, so ``process_arc_run()``'s capture step has something real to
    discover and vendor."""
    arc_dir = t3.paths['ARC']
    label = arc_ts_label(network_id, network_ts_label)
    expected_path = expected_ts_artifact_path(arc_dir, label)
    record = TSJoinRecord(network_id=network_id,
                          network_ts_label=network_ts_label,
                          status=JOIN_STATUS_QUEUED,
                          arc_ts_label=label,
                          expected_artifact_path=expected_path,
                          reason='Queued to ARC.',
                          )
    _write_ts_artifact(expected_path)
    _write_ts_status_yml(arc_dir, label, converged=True)
    write_ts_join_sidecar(arc_dir, [record])
    return record


class TestProcessArcRunFinalizationWiring(object):
    """Acceptance tests for the ARC finalization wiring added to ``process_arc_run()`` and
    ``restart()``: durable one-shot transition-state capture, and the marker-backed
    ``restart()`` branch that repairs a T3 crash occurring after ARC terminates but before
    finalization (capture + marker) completes.
    """

    def test_process_arc_run_captures_queued_transition_states(self, tmp_path):
        """The ordinary path: a queued, usable PDep transition state is captured into
        'PDep capture', and the finalization marker is written once capture succeeds."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        record = _queue_usable_ts(t3)

        t3.process_arc_run()

        capture_dir = t3.paths['PDep capture']
        assert os.path.isdir(capture_dir)
        assert os.path.isfile(t3.paths['ARC finalization marker'])
        # the capture dir must be a sibling of ARC, never nested inside it
        assert os.path.commonpath([os.path.realpath(capture_dir), os.path.realpath(t3.paths['ARC'])]) \
            != os.path.realpath(t3.paths['ARC'])
        assert record.arc_ts_label == arc_ts_label('network4_2', 'TS9')

        # Concrete structure, not merely "some file exists somewhere": the vendored pointer file
        # capture_ts_artifacts writes must be exactly where the manifest says it is, its referenced
        # log(s) must actually exist inside the capture, and both must be recorded in the manifest
        # for this transition state -- a directory walk finding "any file" would also pass on a
        # capture that vendored the wrong thing, or recorded paths that don't resolve.
        manifest = read_yaml_file(path=os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME))
        entries = [entry for entry in manifest['transition_states']
                  if entry['network_id'] == record.network_id
                  and entry['network_ts_label'] == record.network_ts_label]
        assert len(entries) == 1
        entry = entries[0]
        assert entry['captured_artifact_path'] is not None
        vendored_pointer_path = os.path.join(capture_dir, entry['captured_artifact_path'])
        assert os.path.isfile(vendored_pointer_path)
        assert os.path.basename(vendored_pointer_path) == f'{record.arc_ts_label}.py'
        assert entry['captured_log_paths'], 'a usable capture must record at least one vendored log'
        for relpath in entry['captured_log_paths'].values():
            assert os.path.isfile(os.path.join(capture_dir, relpath))

    def test_process_arc_run_with_no_join_records_is_a_noop_for_capture(self, tmp_path):
        """An iteration that queued no PDep transition states -- the common case -- must incur no
        capture side effects at all, while finalization still completes and the marker is written."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        assert not os.path.isfile(ts_join_sidecar_path(t3.paths['ARC']))

        t3.process_arc_run()

        assert not os.path.exists(t3.paths['PDep capture'])
        assert os.path.isfile(t3.paths['ARC finalization marker'])

    def test_restart_finalizes_when_arc_terminated_but_marker_absent(self, tmp_path):
        """The crash-recovery gap: ARC began and terminated for this iteration, but no PDep
        finalization marker exists yet (e.g. T3 crashed between ARC terminating and
        process_arc_run() finishing). restart() must detect this via check_arc_finalization_complete()
        and finalize (capture + mark), then advance to the next iteration requesting an RMG run --
        rather than silently falling through to 'skipping RMG' and losing the completed run's
        transition-state artifacts on ARC's next rate pass."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        _mark_rmg_terminated(t3)
        _mark_arc_terminated(t3)
        # Captured before restart() runs: handle_arc_finalization_restart() advances i_max and
        # re-calls set_paths(), so t3.paths mutates mid-call to point at the NEXT iteration.
        # Finalization itself belongs to the iteration that just terminated (iteration 1 here),
        # so the paths to assert against must be resolved now, not read back off t3.paths after
        # restart() returns.
        marker_path = t3.paths['ARC finalization marker']
        capture_dir = t3.paths['PDep capture']
        assert not os.path.isfile(marker_path)

        result = t3.restart()

        assert result == (2, True, False)
        assert os.path.isfile(marker_path)
        assert os.path.isdir(capture_dir)
        # A bare directory existing proves nothing about whether restart() actually finalized the
        # capture correctly -- verify_capture() is the real oracle: it raises on a missing/malformed/
        # torn manifest, and its counts confirm this is a genuine, non-empty, one-transition-state
        # capture rather than an empty or partially-written one that happened to leave a directory
        # behind.
        verified = verify_capture(capture_dir)
        assert verified.record_count == 1
        assert verified.captured_artifact_count == 1

    def test_restart_does_not_double_capture_when_marker_already_present(self, tmp_path, monkeypatch):
        """Once finalization has already completed for an iteration (marker present),
        restart() must not re-run it: not incrementing the iteration, and not invoking capture a
        second time (ARC deletes calcs/statmech/kinetics/ on its next rate pass, so a second
        capture attempt after that point would find nothing, or -- if it ran before ARC's next
        pass -- would be pure repeated, unnecessary work)."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _mark_rmg_terminated(t3)
        _mark_arc_terminated(t3)
        # Finalize once for real, so the marker legitimately exists going into restart().
        t3.process_arc_run()
        assert os.path.isfile(t3.paths['ARC finalization marker'])
        capture_calls = []
        monkeypatch.setattr(t3, '_capture_pdep_ts_artifacts', lambda: capture_calls.append(1))

        result = t3.restart()

        assert result == (1, False, False)
        assert capture_calls == [], 'capture must not be invoked again once the marker is present'

    def test_marker_is_not_written_if_capture_fails_partway(self, tmp_path, monkeypatch):
        """Fail-closed: if capture_ts_artifacts() raises, process_arc_run() must propagate the
        exception rather than swallow it, and -- because the marker is written only as the very
        last step -- the marker must be left absent so a subsequent restart() retries finalization
        rather than wrongly treating the failed run as done."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)

        def _boom(**kwargs):
            raise ValueError('simulated capture failure')

        monkeypatch.setattr(t3_main, 'capture_ts_artifacts', _boom)

        with pytest.raises(ValueError, match='simulated capture failure'):
            t3.process_arc_run()

        assert not os.path.isfile(t3.paths['ARC finalization marker'])
        # `_queue_usable_ts` writes the join sidecar straight to disk and never populates
        # `t3.pdep_ts_join_records` (that in-memory list is only ever filled within a single
        # `determine_species_from_pdep_network()` call), so asserting against it here would be
        # vacuously true regardless of what capture did. What must actually hold is that the durable
        # join sidecar on disk -- the only record `_capture_pdep_ts_artifacts()` reads on a retry --
        # survives a failed capture attempt untouched, still naming the transition state as queued,
        # so a subsequent retry does not silently lose it.
        join_records_on_disk = read_ts_join_sidecar(t3.paths['ARC'])
        assert len(join_records_on_disk) == 1
        assert join_records_on_disk[0].status == JOIN_STATUS_QUEUED

    def test_capture_is_skipped_when_an_existing_verified_capture_already_matches(self, tmp_path, monkeypatch):
        """Replay guard (round-20 finding 6): a verified capture already on disk, whose identity
        matches this iteration's join records, is authoritative -- re-running capture_ts_artifacts
        would be redundant at best, and dangerous if ARC has since deleted the artifacts the first
        capture read. run_arc() never touches the capture directory (only the finalization marker),
        so simulating a second finalization attempt for the same iteration means clearing only the
        marker, exactly as run_arc() does."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()
        assert os.path.isfile(t3.paths['ARC finalization marker'])

        t3._clear_arc_finalization_marker()
        assert not os.path.isfile(t3.paths['ARC finalization marker'])

        def _boom(**kwargs):
            raise AssertionError('capture_ts_artifacts must not be called again once an existing '
                                 'verified capture already matches this iteration')
        monkeypatch.setattr(t3_main, 'capture_ts_artifacts', _boom)
        t3.process_arc_run()

        assert os.path.isfile(t3.paths['ARC finalization marker'])

    def test_stale_capture_is_not_replayed_when_network_source_changes(self, tmp_path, monkeypatch):
        """Round-23 finding: the ARC project directory and the transition-state set can both still
        match while the underlying RMG network file has been regenerated (different
        ``source_sha256``) since the capture was written. That capture no longer describes the
        current network and must NOT be replayed as authoritative -- capture must be attempted
        again (and, in real use, would vendor the new file fresh)."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        record = _queue_usable_ts(t3)
        t3.process_arc_run()
        assert os.path.isfile(t3.paths['ARC finalization marker'])
        manifest_path = os.path.join(t3.paths['PDep capture'], CAPTURE_MANIFEST_FILE_NAME)
        original_sha256 = read_yaml_file(path=manifest_path)['networks']['network4_2']['source_sha256']

        t3._clear_arc_finalization_marker()
        # Regenerate the RMG network file (as a fresh RMG run would) and re-write the join sidecar
        # against it, exactly as the in-run path would for a new selection pass.
        network_path = os.path.join(t3.paths['RMG PDep'], 'network4_2.py')
        with open(network_path, 'a') as f:
            f.write("# regenerated: an extra reaction was added to this network\n")
        _write_sidecar(t3, [record])
        new_sha256 = read_ts_join_networks(t3.paths['ARC'])['network4_2']['source_sha256']
        assert new_sha256 != original_sha256, 'the fixture mutation must actually change the hash'

        def _boom(**kwargs):
            raise AssertionError('capture_ts_artifacts must be called again: the existing capture is '
                                 'stale on network identity and must not be replayed as authoritative')
        monkeypatch.setattr(t3_main, 'capture_ts_artifacts', _boom)
        with pytest.raises(AssertionError, match='must be called again'):
            t3.process_arc_run()

    def test_stale_capture_is_not_replayed_when_me_method_changes(self, tmp_path, monkeypatch):
        """Round-23 finding: the ARC project directory, the transition-state set, and the network
        source file can all still match while the master-equation method (CSE vs MSC vs RS) recorded
        for the network has changed since the capture was written. The frozen artifacts describe a
        different calculation and must NOT be replayed as authoritative."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        record = _queue_usable_ts(t3)
        t3.process_arc_run()
        assert os.path.isfile(t3.paths['ARC finalization marker'])
        manifest_path = os.path.join(t3.paths['PDep capture'], CAPTURE_MANIFEST_FILE_NAME)
        original_method = read_yaml_file(path=manifest_path)['networks']['network4_2']['method']
        assert original_method == 'MSC'

        t3._clear_arc_finalization_marker()
        # Re-write the join sidecar naming a different ME method for the same, unchanged network file
        # -- exactly what a re-run with a different writer.gen_configuration setting would produce.
        _write_sidecar(t3, [record], method='CSE')
        new_method = read_ts_join_networks(t3.paths['ARC'])['network4_2']['method']
        assert new_method != original_method, 'the fixture mutation must actually change the method'

        def _boom(**kwargs):
            raise AssertionError('capture_ts_artifacts must be called again: the existing capture is '
                                 'stale on ME method identity and must not be replayed as authoritative')
        monkeypatch.setattr(t3_main, 'capture_ts_artifacts', _boom)
        with pytest.raises(AssertionError, match='must be called again'):
            t3.process_arc_run()

    def test_replay_of_a_partial_capture_still_refuses_to_finalize(self, tmp_path):
        """The two new guards must not cancel each other out. A capture in which SOME queued
        transition states produced artifacts and others did not is refused on the first attempt by
        the missing-artifact guard -- correctly, and the marker stays absent. But that partial
        capture is left on disk, and it both verifies and matches this iteration's identity, with
        captured_artifact_count > 0. On the replay the authoritative-capture guard would therefore
        skip re-capturing and return EARLY, jumping over the missing-artifact refusal entirely, and
        finalization would complete -- writing the marker and permanently abandoning the transition
        state that never produced QM. The refusal must be evaluated on the skip path too, from the
        manifest's own frozen per-transition-state statuses, not only on the freshly-captured path."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        arc_dir = t3.paths['ARC']
        captured_label = arc_ts_label('network4_2', 'TS9')
        captured_path = expected_ts_artifact_path(arc_dir, captured_label)
        _write_ts_artifact(captured_path)
        _write_ts_status_yml(arc_dir, captured_label, converged=True)
        missing_label = arc_ts_label('network4_2', 'TS7')
        records = [TSJoinRecord(network_id='network4_2',
                                network_ts_label='TS9',
                                status=JOIN_STATUS_QUEUED,
                                arc_ts_label=captured_label,
                                expected_artifact_path=captured_path,
                                reason='Queued to ARC.',
                                ),
                   # queued, but ARC never produced this one
                   TSJoinRecord(network_id='network4_2',
                                network_ts_label='TS7',
                                status=JOIN_STATUS_QUEUED,
                                arc_ts_label=missing_label,
                                expected_artifact_path=expected_ts_artifact_path(arc_dir, missing_label),
                                reason='Queued to ARC.',
                                ),
                   ]
        write_ts_join_sidecar(arc_dir, records)

        with pytest.raises(ValueError, match='produced no artifact'):
            t3.process_arc_run()
        assert not os.path.isfile(t3.paths['ARC finalization marker'])
        # the partial capture is left behind, and it is a VALID, identity-matching, non-empty one
        verified = verify_capture(t3.paths['PDep capture'])
        assert verified.captured_artifact_count == 1

        # the replay: same iteration, marker still absent, capture still on disk
        with pytest.raises(ValueError, match='produced no artifact'):
            t3.process_arc_run()

        assert not os.path.isfile(t3.paths['ARC finalization marker']), \
            'a replay must not finalize past a transition state that never produced QM'

    def test_capture_raises_when_existing_capture_names_a_different_arc_project_directory(self, tmp_path):
        """Fail-closed identity check: a capture directory that already holds a verified capture, but
        one whose manifest names a different ARC project directory than this iteration's, must never
        be silently reused or overwritten -- doing so could attribute another ARC run's QM results to
        this iteration."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()
        capture_dir = t3.paths['PDep capture']
        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        manifest = read_yaml_file(path=manifest_path)
        manifest['arc_project_directory'] = os.path.join(str(tmp_path), 'a_different_arc_project')
        save_yaml_file(path=manifest_path, content=manifest)
        t3._clear_arc_finalization_marker()

        with pytest.raises(ValueError, match='different ARC project directory'):
            t3.process_arc_run()

    def test_capture_raises_when_existing_capture_names_a_different_transition_state_set(self, tmp_path):
        """Fail-closed identity check: a capture directory that already holds a verified capture from
        the SAME ARC project directory, but naming a different set of transition states than this
        iteration's join sidecar, is exactly the ambiguous case this module refuses to guess through
        -- reusing it could silently drop or fabricate transition state coverage."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()
        capture_dir = t3.paths['PDep capture']
        manifest_path = os.path.join(capture_dir, CAPTURE_MANIFEST_FILE_NAME)
        manifest = read_yaml_file(path=manifest_path)
        manifest['transition_states'][0]['network_ts_label'] = 'TS_SOME_OTHER_LABEL'
        save_yaml_file(path=manifest_path, content=manifest)
        t3._clear_arc_finalization_marker()

        with pytest.raises(ValueError, match='different set of transition states'):
            t3.process_arc_run()

    def test_process_arc_run_raises_when_a_queued_transition_state_has_no_artifact(self, tmp_path):
        """Defect 2: a transition state actually QUEUED to ARC, but whose artifact is missing after
        this capture pass (ARC never produced it, or it was cleaned up before capture ran), must
        refuse to finalize -- silently continuing would abandon QM work T3 believed it was waiting on,
        leaving an incomplete manifest with no trace anything is missing at all."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        arc_dir = t3.paths['ARC']
        label = arc_ts_label('network4_2', 'TS9')
        expected_path = expected_ts_artifact_path(arc_dir, label)
        record = TSJoinRecord(network_id='network4_2',
                              network_ts_label='TS9',
                              status=JOIN_STATUS_QUEUED,
                              arc_ts_label=label,
                              expected_artifact_path=expected_path,
                              reason='Queued to ARC.',
                              )
        # deliberately do NOT write the artifact or its status.yml entry: ARC either never produced
        # it, or it was removed before capture ran here.
        write_ts_join_sidecar(arc_dir, [record])

        with pytest.raises(ValueError, match='produced no artifact'):
            t3.process_arc_run()

        assert not os.path.isfile(t3.paths['ARC finalization marker'])

    def test_process_arc_run_does_not_raise_for_missing_artifacts_that_were_never_queued(self, tmp_path):
        """Defect 2's refusal must be scoped to transition states that were actually QUEUED to ARC:
        `already_present` (passthrough, no artifact ever expected from this ARC run) and `not_queued`
        (never sent to ARC at all) both legitimately produce an ARTIFACT_STATUS_MISSING discovery
        record, and neither is evidence of lost QM work, so neither may trip the refusal."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        arc_dir = t3.paths['ARC']
        already_present_record = TSJoinRecord(network_id='network4_2',
                                              network_ts_label='TS1',
                                              status=JOIN_STATUS_ALREADY_PRESENT,
                                              reason='Already present upstream; ARC never queued it.',
                                              )
        not_queued_record = TSJoinRecord(network_id='network4_2',
                                         network_ts_label='TS2',
                                         status=JOIN_STATUS_NOT_QUEUED,
                                         reason='Not selected for QM refinement.',
                                         )
        write_ts_join_sidecar(arc_dir, [already_present_record, not_queued_record])

        t3.process_arc_run()  # must not raise

        assert os.path.isfile(t3.paths['ARC finalization marker'])

    def test_check_arc_finalization_complete_treats_an_empty_marker_as_absent(self, tmp_path):
        """Guards against trusting mere existence of the marker file: an empty marker (e.g. a crash
        that created but never wrote to the staged file, or a hand-touched file) must read as
        'not finalized', or a run that never actually finalized would be skipped forever."""
        t3 = _build_t3(tmp_path)
        marker_path = t3.paths['ARC finalization marker']
        os.makedirs(os.path.dirname(marker_path), exist_ok=True)
        with open(marker_path, 'w') as f:
            f.write('')

        assert t3.check_arc_finalization_complete() is False

    def test_check_arc_finalization_complete_treats_garbage_content_as_absent(self, tmp_path):
        """Guards against a marker file that exists but does not carry the expected finalization
        text (e.g. corrupted, or written by unrelated code that happened to reuse the path): only
        content starting with ARC_FINALIZATION_MARKER_TEXT may be trusted as proof finalization
        actually ran, otherwise a run that never finalized would look done."""
        t3 = _build_t3(tmp_path)
        marker_path = t3.paths['ARC finalization marker']
        os.makedirs(os.path.dirname(marker_path), exist_ok=True)
        with open(marker_path, 'w') as f:
            f.write('this is not the marker text at all\n')

        assert t3.check_arc_finalization_complete() is False

    def test_check_arc_finalization_complete_returns_false_on_oserror(self, tmp_path, monkeypatch):
        """Guards against an unreadable marker crashing restart()/process_arc_run() instead of being
        treated as 'not finalized': a marker that exists but cannot be read (permissions, I/O error,
        a full filesystem) must be caught and read as absence rather than propagate and abort the
        whole run -- redoing finalization is recoverable, aborting the run is not.

        The marker is a real, readable-by-stat file here so that os.path.isfile() passes and the
        OSError branch is genuinely reached; failing that, this test would pass on the earlier
        isfile() guard alone and assert nothing about OSError handling at all."""
        t3 = _build_t3(tmp_path)
        marker_path = t3.paths['ARC finalization marker']
        os.makedirs(os.path.dirname(marker_path), exist_ok=True)
        with open(marker_path, 'w') as f:
            f.write(f'{t3_main.ARC_FINALIZATION_MARKER_TEXT} at some point\n')
        real_open = builtins.open

        def _open_raising_on_the_marker(file, *args, **kwargs):
            if file == marker_path:
                raise OSError('simulated unreadable marker')
            return real_open(file, *args, **kwargs)

        monkeypatch.setattr(builtins, 'open', _open_raising_on_the_marker)

        assert t3.check_arc_finalization_complete() is False

    def test_restart_finalizes_despite_a_bad_marker_when_join_sidecar_present(self, tmp_path):
        """Guards against a corrupted/garbage marker being read as 'already finalized' and getting
        restart() stuck skipping a run that was never actually captured: with a join sidecar present
        (there IS something to rescue) and ARC terminated, a marker file that exists but fails
        content validation must still be treated as absent, so restart() finalizes rather than
        silently skipping to 'skipping RMG' and losing the completed run's TS artifacts."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        _mark_rmg_terminated(t3)
        _mark_arc_terminated(t3)
        marker_path = t3.paths['ARC finalization marker']
        os.makedirs(os.path.dirname(marker_path), exist_ok=True)
        with open(marker_path, 'w') as f:
            f.write('garbage, not the real marker text')
        capture_dir = t3.paths['PDep capture']

        result = t3.restart()

        assert result == (2, True, False)
        assert os.path.isdir(capture_dir)

    def test_process_arc_run_redundant_direct_call_does_not_re_append_to_rmg_libraries(self, tmp_path, monkeypatch):
        """Guards against the idempotency short-circuit at the top of process_arc_run() being
        bypassed or removed: a second, redundant DIRECT call after the marker is already present
        (not via restart()'s own capture-monkeypatching path) must be a true no-op, in particular
        must NOT invoke append_to_rmg_libraries again -- doing so would append the same converged
        results to the RMG libraries a second time."""
        t3 = _build_t3(tmp_path)
        t3.species = {0: T3Species(label='dummy', smiles='[OH]')}
        _write_arc_info(t3, species=[{'label': 'dummy', 'success': True}])
        calls = []
        monkeypatch.setattr(t3_main, 'append_to_rmg_libraries', lambda **kwargs: calls.append(1))

        t3.process_arc_run()
        assert calls == [1], 'the first, real call should have appended to the RMG libraries once'
        assert os.path.isfile(t3.paths['ARC finalization marker'])

        t3.process_arc_run()

        assert calls == [1], \
            'a redundant direct call to process_arc_run() must not re-invoke append_to_rmg_libraries'

    def test_restart_does_not_finalize_a_legacy_project_with_no_join_sidecar(self, tmp_path, monkeypatch):
        """Guards against the legacy-project regression: a project finalized by an older T3 that
        never wrote join sidecars or markers must NOT be re-finalized just because a marker happens
        to be absent. With RMG terminated, ARC terminated, no marker, and no join sidecar, restart()
        must fall through to the pre-existing 'skipping RMG' behavior -- not call process_arc_run()
        (which would re-append convergence results to the RMG libraries and advance the iteration
        under a project that was already done)."""
        t3 = _build_t3(tmp_path)
        _mark_rmg_terminated(t3)
        _mark_arc_terminated(t3)
        assert not os.path.isfile(t3.paths['ARC finalization marker'])
        assert not os.path.isfile(ts_join_sidecar_path(t3.paths['ARC']))
        calls = []
        monkeypatch.setattr(t3, 'process_arc_run', lambda: calls.append(1))

        result = t3.restart()

        assert result == (1, False, False)
        assert calls == [], 'process_arc_run() must not be called for a legacy project with no join sidecar'
        assert not os.path.isfile(t3.paths['ARC finalization marker'])

    def test_run_arc_clears_a_stale_finalization_marker_before_executing(self, tmp_path, monkeypatch):
        """Guards against a new ARC run being treated as already finalized by a marker left over
        from an earlier run of the same iteration: run_arc() must clear the marker before invoking
        ARC, not after -- otherwise a crash during arc.execute() would leave the stale marker in
        place, and a subsequent restart would wrongly skip finalizing the new run's results."""
        t3 = _build_t3(tmp_path)
        marker_path = t3.paths['ARC finalization marker']
        os.makedirs(os.path.dirname(marker_path), exist_ok=True)
        with open(marker_path, 'w') as f:
            f.write(f'{t3_main.ARC_FINALIZATION_MARKER_TEXT} at some point\n')
        assert os.path.isfile(marker_path)
        marker_present_at_execute_time = []

        class _StubARC(object):
            def __init__(self, **kwargs):
                pass

            def execute(self):
                marker_present_at_execute_time.append(os.path.isfile(marker_path))

            def as_dict(self):
                return {}

        monkeypatch.setattr(t3_main, 'ARC', _StubARC)

        t3.run_arc(arc_kwargs={'project': 'test'})

        assert marker_present_at_execute_time == [False], \
            'the marker must already be cleared by the time arc.execute() runs'
        assert not os.path.isfile(marker_path)

    def test_capture_raises_when_orphaned_ts_artifacts_exist_without_a_join_sidecar(self, tmp_path):
        """Guards against silently abandoning completed QM work: if T3-labelled transition state
        artifacts exist under ARC's TS directory but no join sidecar names them (the sidecar that
        would have attributed them to a P-dep network was lost), process_arc_run() must raise rather
        than treat the absent sidecar as 'nothing was queued' and quietly proceed -- ARC deletes and
        recreates calcs/statmech/kinetics/ on its next rate pass, so there would be no second chance
        to capture what failed here."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        orphan_label = arc_ts_label('network4_2', 'TS9')
        orphan_path = expected_ts_artifact_path(t3.paths['ARC'], orphan_label)
        _write_ts_artifact(orphan_path)
        assert not os.path.isfile(ts_join_sidecar_path(t3.paths['ARC']))

        with pytest.raises(ValueError, match='join record was lost'):
            t3.process_arc_run()

        assert not os.path.isfile(t3.paths['ARC finalization marker'])

    def test_capture_ignores_an_unrelated_non_ts_file_without_a_join_sidecar(self, tmp_path):
        """Guards the negative case for the orphaned-artifact refusal: a .py file in ARC's TS
        directory that does NOT carry the T3PDep_ label prefix is not T3's own work and must not
        trigger the refusal -- process_arc_run() must proceed normally, proving the orphan check is
        actually scoped to T3-labelled artifacts and not to 'any .py file in that directory'."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        orphan_path = expected_ts_artifact_path(t3.paths['ARC'], 'network4_2_TS9')
        _write_ts_artifact(orphan_path)  # no T3PDep_ prefix
        assert not os.path.isfile(ts_join_sidecar_path(t3.paths['ARC']))

        t3.process_arc_run()

        assert os.path.isfile(t3.paths['ARC finalization marker'])
