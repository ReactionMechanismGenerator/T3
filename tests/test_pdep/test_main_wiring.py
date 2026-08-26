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
import hashlib
import os
import shutil

import pytest

import t3.main as t3_main
from t3.chem import T3Species
from t3.pdep.budget import (BUDGET_OUTCOME_ADMITTED,
                            BUDGET_OUTCOME_REFUSED,
                            BUDGET_SKIP_DOES_NOT_FIT_REMAINING,
                            PDepBudgetConsideration,
                            PDepBudgetDecision,
                            PDepBudgetSkip,
                            )
from t3.pdep.cache import hash_file, sa_cache_metadata_path, write_sa_cache_metadata
from t3.pdep.capture import CAPTURE_MANIFEST_FILE_NAME, verify_capture
from t3.pdep.discovery import ARTIFACT_STATUS_USABLE
from t3.pdep.join import (JOIN_STATUS_ALREADY_PRESENT,
                          JOIN_STATUS_NOT_QUEUED,
                          JOIN_STATUS_QUEUED,
                          TSJoinRecord,
                          arc_ts_label,
                          expected_ts_artifact_path,
                          read_ts_join_networks,
                          read_ts_join_sidecar,
                          ts_join_sidecar_path,
                          write_ts_join_sidecar,
                          )
from t3.pdep.api import load_pdep_budget_record
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


@pytest.fixture(autouse=True)
def _stub_extract_dof_conformers(monkeypatch):
    """Stub the one Arkane-invoking seam ``_write_pdep_hybrid_network_inputs()`` now calls
    (``extract_dof_conformers``, which shells out to ``rmg_env`` and is unavailable in tests) so the
    real hybrid writer still runs end to end against synthetic, correctly-shaped conformer data.
    ``wells`` is always ``{}`` from ``t3/main.py``, so the well half stays empty. Module-scoped
    (autouse) so every process_arc_run-driving test in this module gets it."""
    def _fake_extract(transition_states, wells, energy_settings, **kwargs):
        def _conformer(label, is_ts):
            conformer = {'label': label, 'is_ts': is_ts, 'E0_kJ_mol': -38.0,
                         'frequencies_cm_1': [500.0, 800.0, 1200.0],
                         'spin_multiplicity': 1, 'optical_isomers': 1, 'hindered_rotors': []}
            if is_ts:
                conformer['imaginary_frequency_cm_1'] = -1800.0
            return conformer
        return ({label: _conformer(label, True) for label in transition_states},
                {label: _conformer(label, False) for label in wells})
    monkeypatch.setattr(t3_main, 'extract_dof_conformers', _fake_extract)


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


def _selection_with_evidence(entries: list) -> PDepNetworkSelection:
    """A qualified selection whose ``uncertain_path_reactions`` carry the given
    ``(ts_label, coefficient, delta_ln_k)`` evidence tuples verbatim, allowing several entries to
    share one ``ts_label`` (unlike ``_selection``, which assigns a single uniform pair per label).
    """
    built = [SensitiveTransitionState(ts_label=ts_label,
                                      coefficient=coefficient,
                                      condition=CONDITION,
                                      path_reaction_label=None,
                                      path_reaction_str=None,
                                      kinetics_comment='',
                                      uncertain=True,
                                      delta_ln_k=delta_ln_k,
                                      ) for ts_label, coefficient, delta_ln_k in entries]
    return PDepNetworkSelection(network_id=NETWORK_NAME,
                                qualified=True,
                                selected_ts=list(built),
                                uncertain_path_reactions=list(built),
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

    def test_cached_valid_path_never_rewrites_the_arkane_input(self, tmp_path, monkeypatch):
        """A cache hit must leave ``input.py`` exactly as it was, because nothing consumes it here.

        No Arkane job runs on this path, so rewriting the input buys nothing and costs provenance:
        after any change to the writer, the directory would hold a freshly rendered ``input.py``
        that did NOT produce the ``sensitivity/sa_coefficients.yml`` sitting beside it, and neither
        the sidecar's ``network_file_hash`` nor its ``sa_file_hash`` covers ``input.py`` to catch
        the divergence.
        """
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
        input_path = os.path.join(t3.paths['PDep SA'], NETWORK_NAME, 'CSE', 'input.py')
        sentinel = '# the input.py that actually produced the cached SA\n'
        with open(input_path, 'w') as f:
            f.write(sentinel)

        def _fail_if_called(*args, **kwargs):
            pytest.fail('run_arkane_job should not be invoked when a valid cache is present.')

        monkeypatch.setattr(t3_main, 'run_arkane_job', _fail_if_called)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert t3.pdep_network_selections[0].cache_status == CACHE_STATUS_CACHED_VALID
        with open(input_path, 'r') as f:
            assert f.read() == sentinel, 'the cache hit rewrote an input.py that nothing was going to read'

    @pytest.mark.skipif(hasattr(os, 'geteuid') and os.geteuid() == 0,
                        reason='root bypasses directory permissions, so the read-only dir would be writable')
    def test_a_valid_cache_is_not_denied_by_an_unwritable_output_directory(self, tmp_path, monkeypatch):
        """A cache hit must not depend on being able to write into the output directory.

        The Arkane input used to be rendered before the cache was consulted, so a read-only or full
        ``<network>/<method>/`` directory raised an OSError and ended the network -- discarding a
        cached SA that was perfectly valid and that the run never needed to write anything to use.
        """
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

        method_dir = os.path.join(t3.paths['PDep SA'], NETWORK_NAME, 'CSE')
        os.chmod(method_dir, 0o555)
        try:
            t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)
        finally:
            os.chmod(method_dir, 0o755)

        assert len(t3.pdep_network_selections) == 1
        selection = t3.pdep_network_selections[0]
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

    def test_a_successful_arkane_run_records_the_rmg_py_provenance_in_the_sidecar(self, tmp_path, monkeypatch):
        """When Arkane runs successfully, the sidecar T3 writes must carry the RMG-Py commit that
        actually produced the sensitivity analysis, read out of the arkane.log Arkane leaves in the
        output directory -- a stale real run went undetected for a week because this field was
        never populated by the only production caller (t3/main.py). The log is the only witness
        available: ARC runs Arkane in a subprocess in a different conda environment, so the RMG-Py
        that did the work is never loaded into this interpreter and cannot be introspected here."""
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_dict = _build_sa_dict(t3)
        rmg_py_commit = 'e720866ae94eca51652978c15a0fb33c6827be67'

        def _fake_run_arkane_job(input_file, output_directory, plot=True, logger=None, **kwargs):
            sensitivity_dir = os.path.join(output_directory, 'sensitivity')
            os.makedirs(sensitivity_dir, exist_ok=True)
            save_yaml_file(os.path.join(sensitivity_dir, 'sa_coefficients.yml'), sa_dict)
            # The stanza Arkane really writes: the hash is on the line after the label, and the
            # line after that is a date.
            with open(os.path.join(output_directory, 'arkane.log'), 'w') as f:
                f.write(f'Arkane execution initiated at Sat Aug  1 17:28:25 2026\n\n'
                        f'The current git HEAD for RMG-Py is:\n'
                        f'\t{rmg_py_commit}\n'
                        f'\tFri Jul 31 17:25:49 2026 +0300\n')
            return True

        monkeypatch.setattr(t3_main, 'run_arkane_job', _fake_run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        selection = t3.pdep_network_selections[0]
        assert selection.cache_status == CACHE_STATUS_GENERATED
        assert selection.qualified is True

        sa_path = _candidate_sa_path(t3, method='CSE')
        assert os.path.isfile(sa_path)
        metadata = read_yaml_file(sa_cache_metadata_path(sa_path))
        assert metadata['rmg_py_commit'] == rmg_py_commit, \
            'The sidecar must record the RMG-Py commit that actually ran Arkane, not null.'

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

    def test_cached_valid_path_threads_t_grid_clamp_from_the_sidecar(self, tmp_path, monkeypatch):
        """The cached-valid wiring path reads ``t_grid_clamp`` from the sidecar via
        ``read_t_grid_clamp_record`` and threads it onto the resulting selection, so a saved
        ``PDepNetworkSelection`` remains self-describing about T-grid clamping even when Arkane
        is never re-invoked (pure cache reuse).
        """
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_path = _candidate_sa_path(t3, method='CSE')
        os.makedirs(os.path.dirname(sa_path), exist_ok=True)
        save_yaml_file(sa_path, _build_sa_dict(t3))
        clamp_record = {'clamped': True,
                         'requested_t_max': 3200.0,
                         'thermo_ceiling': 3000.0,
                         'written_t_max': 3000.0,
                         'tlist_dropped': True,
                         'tlist_original_highest': 3200.0,
                         'skipped_species': ['spcA'],
                         }
        write_sa_cache_metadata(sa_path=sa_path,
                                network_path=_network_path(t3),
                                network_id=NETWORK_NAME,
                                method='CSE',
                                t_grid_clamp=clamp_record,
                                )

        def _fail_if_called(*args, **kwargs):
            pytest.fail('run_arkane_job should not be invoked when a valid cache is present.')

        monkeypatch.setattr(t3_main, 'run_arkane_job', _fail_if_called)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert len(t3.pdep_network_selections) == 1
        selection = t3.pdep_network_selections[0]
        assert selection.cache_status == CACHE_STATUS_CACHED_VALID
        assert selection.t_grid_clamp == clamp_record

    def test_cached_valid_path_records_no_t_grid_clamp_when_the_sidecar_has_none(self, tmp_path, monkeypatch):
        """A sidecar written without a ``t_grid_clamp`` key (or an older sidecar predating this
        feature) must surface as ``None`` (unknown provenance), never as some other value, and
        must not cause the selection to be rejected or marked ``not_evaluated``.
        """
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

        selection = t3.pdep_network_selections[0]
        assert selection.cache_status == CACHE_STATUS_CACHED_VALID
        assert selection.t_grid_clamp is None
        assert selection.qualified is True

    def test_queueing_is_deferred_until_every_network_has_been_evaluated(self, tmp_path, monkeypatch):
        """No network may be queued before the last one has been evaluated.

        Whether a network is worth its QM cost is a question about the whole field of candidates, so
        the answer cannot be given while candidates are still being discovered. Two entries are
        explored here (both resolving to the same fixture network, which is all that is needed to
        order the evaluations against the queue calls), and every call into
        ``queue_pdep_transition_states`` must see BOTH selections already recorded. When the queue
        call still lived inside the loop, the first one saw only the single selection its own pass
        had just produced.
        """
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3) * 2
        sa_path = _candidate_sa_path(t3, method='CSE')
        os.makedirs(os.path.dirname(sa_path), exist_ok=True)
        save_yaml_file(sa_path, _sa_dict_with_ts1_structures(t3))
        write_sa_cache_metadata(sa_path=sa_path,
                                network_path=_network_path(t3),
                                network_id=NETWORK_NAME,
                                method='CSE',
                                )
        monkeypatch.setattr(t3_main, 'run_arkane_job',
                            lambda *args, **kwargs: pytest.fail('Arkane must not run; a valid cache is present.'))

        selections_when_queued = list()
        queue_pdep_transition_states = t3.queue_pdep_transition_states

        def _spy(**kwargs):
            selections_when_queued.append(len(t3.pdep_network_selections))
            return queue_pdep_transition_states(**kwargs)

        monkeypatch.setattr(t3, 'queue_pdep_transition_states', _spy)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert len(t3.pdep_network_selections) == 2
        assert selections_when_queued == [2, 2], \
            'Every queueing decision must be taken with the whole field of networks already evaluated.'
        # Deferring the decision must not drop it: the queueing itself still has to happen.
        assert {record.network_ts_label for record in t3.pdep_ts_join_records} == {'TS1'}
        assert NETWORK_NAME in t3.pdep_ts_join_networks

    def test_a_species_justified_by_a_sensitive_well_is_reported_even_when_a_queued_ts_also_needs_it(
            self, tmp_path, monkeypatch):
        """A species both paths reach must still be reported as determined by the well analysis.

        ``add_reaction`` cascades into ``add_species`` for a queued transition state's own species,
        and ``add_species`` returns ``None`` for a species it already knows. So whichever path
        reaches a shared species first is the one that gets a key back, and only the well path
        appends to the returned ``species_keys``. ``C4rad(5)`` is exactly such a species here: it is
        a sensitive well of this network AND the product of TS1's path reaction, which is queued.

        Deferring the queueing past the loop is what makes the well analysis reach it first, so this
        pins a deliberate consequence of that ordering rather than an incidental one: the return
        value means "species determined to be calculated based on SA", and a species justified by a
        sensitive well belongs in it however else it is also reachable.
        """
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_dict = _sa_dict_with_ts1_structures(t3)
        # A sensitive well, in the same shape the real Arkane fixture uses (a bare species label
        # alongside the '(TS) ' entries), strong enough to clear the selection thresholds.
        sa_dict[NETWORK_REACTION_STR][CONDITION]['C4rad(5)'] = 0.03
        sa_path = _candidate_sa_path(t3, method='CSE')
        os.makedirs(os.path.dirname(sa_path), exist_ok=True)
        save_yaml_file(sa_path, sa_dict)
        write_sa_cache_metadata(sa_path=sa_path,
                                network_path=_network_path(t3),
                                network_id=NETWORK_NAME,
                                method='CSE',
                                )
        monkeypatch.setattr(t3_main, 'run_arkane_job',
                            lambda *args, **kwargs: pytest.fail('Arkane must not run; a valid cache is present.'))

        species_keys = t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert species_keys, 'The sensitive well should have determined at least one species.'
        # T3 legalizes the label for ARC's benefit ('(' -> '[', ')' -> ']'), so the stored species
        # carries 'C4rad[5]' rather than the network's own 'C4rad(5)'.
        reported = {t3.species[key].label for key in species_keys}
        assert 'C4rad[5]' in reported, \
            'A species justified by a sensitive well must be reported even when a queued transition ' \
            'state also brings it onto the QM queue.'
        reasons = t3.species[[key for key in species_keys if t3.species[key].label == 'C4rad[5]'][0]].reasons
        assert any('sensitive well' in reason for reason in reasons)

    def test_the_configured_qm_budget_is_what_reaches_the_budget_decision(self, tmp_path, monkeypatch):
        """The budget knobs must be read from the run's own settings, not defaulted away.

        ``t3/pdep/budget.py`` is unit-tested on its own, so this pins the other half of the contract:
        that the numbers the user configured are the numbers the budget is actually applied with. A
        knob that is validated by the schema and then never passed on would leave every budget test
        green while the budget did nothing.
        """
        t3 = _build_t3(tmp_path)
        t3.t3['sensitivity']['pdep_QM_max_transition_states'] = 7
        t3.t3['sensitivity']['pdep_QM_max_networks'] = 3
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
                            lambda *args, **kwargs: pytest.fail('Arkane must not run; a valid cache is present.'))

        observed = list()
        apply_pdep_qm_budget = t3_main.apply_pdep_qm_budget

        def _spy(selections, **kwargs):
            observed.append(kwargs)
            return apply_pdep_qm_budget(selections, **kwargs)

        monkeypatch.setattr(t3_main, 'apply_pdep_qm_budget', _spy)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert observed == [{'max_transition_states': 7, 'max_networks': 3}]
        # The budget is generous here, so the network is still refined: the knob must not refuse by
        # merely being set.
        assert {record.network_ts_label for record in t3.pdep_ts_join_records} == {'TS1'}

    def test_a_network_the_budget_refuses_is_neither_queued_nor_silently_dropped(self, tmp_path, monkeypatch):
        """A budget refusal must stop the queueing AND be reported.

        The failure mode this guards against is a quiet one: if the refusal only shrank the work,
        the run's logs would be indistinguishable from a run that had nothing left to refine, and a
        network that qualified for QM would disappear without ever being mentioned.
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
        monkeypatch.setattr(t3_main, 'run_arkane_job',
                            lambda *args, **kwargs: pytest.fail('Arkane must not run; a valid cache is present.'))
        refusal = PDepBudgetSkip(network_id=NETWORK_NAME,
                                 cost=1,
                                 remaining_transition_states=0,
                                 reason_code=BUDGET_SKIP_DOES_NOT_FIT_REMAINING,
                                 reason='it needs 1 transition state(s) and only 0 remain(s) of the budget of 4',
                                 )
        # ``considered`` must stay consistent with ``admitted_indices``/``skipped``, because
        # ``build_pdep_budget_record`` (now wired into this code path) refuses to build a record
        # off a decision whose own bookkeeping disagrees with itself.
        consideration = PDepBudgetConsideration(identity=NETWORK_NAME,
                                                network_id=NETWORK_NAME,
                                                offer_indices=(0,),
                                                cost=refusal.cost,
                                                rank=0,
                                                remaining_before=refusal.remaining_transition_states,
                                                outcome=BUDGET_OUTCOME_REFUSED,
                                                reason_code=refusal.reason_code,
                                                reason=refusal.reason,
                                                )
        monkeypatch.setattr(t3_main, 'apply_pdep_qm_budget',
                            lambda selections, **kwargs: PDepBudgetDecision(admitted_indices=tuple(),
                                                                           skipped=(refusal,),
                                                                           total_cost=0,
                                                                           considered=(consideration,),
                                                                           ))
        warnings_logged = list()
        monkeypatch.setattr(t3.logger, 'warning', lambda message: warnings_logged.append(message))

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert not t3.pdep_ts_join_records, 'A refused network must not be queued.'
        assert not t3.pdep_ts_join_networks
        assert not os.path.isfile(ts_join_sidecar_path(t3.paths['ARC']))
        # The network still qualified, and that verdict is not rewritten by the budget.
        assert len(t3.pdep_network_selections) == 1 and t3.pdep_network_selections[0].qualified is True
        assert any(NETWORK_NAME in message and refusal.reason in message for message in warnings_logged), \
            f'The refusal must be reported; logged warnings were: {warnings_logged}'
        # The refusal must also survive as a durable record, not just a log line: a log is not
        # queryable after the run, and this is the artifact downstream tooling would read.
        budget_record_file = t3.paths['PDep QM budget']
        assert os.path.isfile(budget_record_file), \
            'The budget record must be written even though (especially because) it refused a network.'
        record = load_pdep_budget_record(budget_record_file)
        assert record.iteration == t3.iteration
        assert record.total_cost == 0
        assert len(record.network_outcomes) == 1
        outcome = record.network_outcomes[0]
        assert outcome.network_id == NETWORK_NAME
        assert outcome.outcome == BUDGET_OUTCOME_REFUSED
        assert outcome.reason_code == refusal.reason_code
        assert outcome.reason == refusal.reason

    def test_the_budget_record_is_written_even_when_nothing_is_refused(self, tmp_path, monkeypatch):
        """Absence of the budget record file must mean "the budget never ran this iteration", not
        "the budget ran and refused nothing" -- so the record has to be written unconditionally,
        including on the all-admitted path where there is nothing to complain about.
        """
        t3 = _build_t3(tmp_path)
        t3.t3['sensitivity']['pdep_QM_max_transition_states'] = 7
        t3.t3['sensitivity']['pdep_QM_max_networks'] = 3
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
                            lambda *args, **kwargs: pytest.fail('Arkane must not run; a valid cache is present.'))

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        # Sanity check that this really is the nothing-refused path (the budget is generous).
        assert {record.network_ts_label for record in t3.pdep_ts_join_records} == {'TS1'}

        budget_record_file = t3.paths['PDep QM budget']
        assert os.path.isfile(budget_record_file), \
            'The budget record must be written even when nothing was refused: its absence must ' \
            'mean the budget never ran, not that it ran and admitted everything.'
        record = load_pdep_budget_record(budget_record_file)
        assert record.iteration == t3.iteration
        assert record.max_transition_states == 7
        assert record.max_networks == 3
        assert all(outcome.outcome == BUDGET_OUTCOME_ADMITTED for outcome in record.network_outcomes)
        assert {outcome.network_id for outcome in record.network_outcomes} == {NETWORK_NAME}

    def test_a_stale_budget_record_is_cleared_when_a_rerun_of_the_same_iteration_turns_the_feature_off(
            self, tmp_path, monkeypatch):
        """Absence of the budget record file must mean "the budget never ran this iteration" -- but
        the record is only ever written once execution reaches the bottom of this method, so a
        re-run of the SAME iteration that this time exits early (``pdep_SA_threshold`` now ``None``)
        would otherwise leave the PREVIOUS run's record on disk, indistinguishable from a record
        produced by the current run. Pin that the early-return path actively clears any such
        leftover, restoring the invariant rather than abandoning it.
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
        monkeypatch.setattr(t3_main, 'run_arkane_job',
                            lambda *args, **kwargs: pytest.fail('Arkane must not run; a valid cache is present.'))

        # First run: the feature is on, so a budget record is written for this iteration.
        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)
        budget_record_file = t3.paths['PDep QM budget']
        assert os.path.isfile(budget_record_file), 'Precondition: the first run must have written a record.'

        # Second run of the SAME iteration: the feature is now off. The stale record from the first
        # run must not survive to be mistaken for a decision made by this run.
        t3.t3['sensitivity']['pdep_SA_threshold'] = None
        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert not os.path.isfile(budget_record_file), \
            'A stale budget record from a superseded run of this iteration must be cleared when a ' \
            're-run turns the feature off, so its presence cannot be mistaken for a decision made ' \
            'by the current (feature-off) run.'


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
        # The sidecar also freezes the qualified network's identity: the real network file's path,
        # its content hash at selection time (asserted against an independent hashlib computation,
        # never the primitive production used), and the ME method the selection's SA actually used.
        networks = read_ts_join_networks(t3.paths['ARC'])
        assert set(networks) == {NETWORK_NAME}
        entry = networks[NETWORK_NAME]
        assert entry['source_path'] == _network_path(t3)
        with open(_network_path(t3), 'rb') as f:
            assert entry['source_sha256'] == f'sha256:{hashlib.sha256(f.read()).hexdigest()}'
        assert entry['method'] == 'CSE'
        assert networks == t3.pdep_ts_join_networks

    def test_network_identity_freezing_is_fail_closed(self, tmp_path):
        """Test the guards of ``_record_pdep_network_identity``: a selection with no resolved ME
        method must raise (never default), a vanished network file must raise, and a conflicting
        re-record for the same network within one iteration must raise rather than last-write-win.
        The happy path is exercised end to end by the sidecar test above."""
        t3 = _build_t3(tmp_path)
        record = TSJoinRecord(network_id=NETWORK_NAME,
                              network_ts_label='TS1',
                              arc_ts_label=arc_ts_label(NETWORK_NAME, 'TS1'),
                              status=JOIN_STATUS_QUEUED,
                              reason='Queued to ARC.')
        t3.pdep_ts_join_records = [record]
        t3.pdep_ts_join_networks = dict()
        selection = _selection(ts_labels=['TS1'])
        assert selection.method is None
        with pytest.raises(ValueError, match='method'):
            t3._record_pdep_network_identity(network_name=NETWORK_NAME,
                                             network_path=_network_path(t3),
                                             selection=selection)

        selection.method = 'MSC'
        # A selection that recorded no content hash cannot freeze an identity either: the sidecar's
        # source_sha256 is supposed to be the hash of the bytes the selection EXAMINED, and there is
        # nothing honest to put there when the selection never recorded one. Checked before the
        # filesystem guard below, because it is a question about the decision, not about the path.
        assert selection.network_source_hash is None
        with pytest.raises(ValueError, match='no content hash'):
            t3._record_pdep_network_identity(network_name=NETWORK_NAME,
                                             network_path=_network_path(t3),
                                             selection=selection)

        selection.network_source_hash = hash_file(path=_network_path(t3))
        with pytest.raises(ValueError, match='gone'):
            t3._record_pdep_network_identity(network_name=NETWORK_NAME,
                                             network_path=os.path.join(t3.paths['RMG PDep'], 'nonexistent.py'),
                                             selection=selection)

        t3._record_pdep_network_identity(network_name=NETWORK_NAME,
                                         network_path=_network_path(t3),
                                         selection=selection)
        assert set(t3.pdep_ts_join_networks) == {NETWORK_NAME}
        assert t3.pdep_ts_join_networks[NETWORK_NAME]['method'] == 'MSC'
        # The recorded hash is the one the SELECTION carries, not a fresh hash of whatever is on disk
        # now. Those agree in a healthy run -- which is why the end-to-end test above cannot tell the
        # two apart -- so pin the distinction directly: a selection carrying a hash that does NOT
        # match the current file must record its own, because the whole point of the sidecar is to
        # name the bytes the decision examined.
        assert t3.pdep_ts_join_networks[NETWORK_NAME]['source_sha256'] == selection.network_source_hash
        stale = 'sha256:' + 'b' * 64
        assert stale != hash_file(path=_network_path(t3))
        selection.network_source_hash = stale
        t3.pdep_ts_join_networks = dict()
        t3._record_pdep_network_identity(network_name=NETWORK_NAME,
                                         network_path=_network_path(t3),
                                         selection=selection)
        assert t3.pdep_ts_join_networks[NETWORK_NAME]['source_sha256'] == stale

        # An identical re-record (a second reaction of the same network) is absorbed...
        t3.pdep_ts_join_networks = dict()
        selection.network_source_hash = hash_file(path=_network_path(t3))
        t3._record_pdep_network_identity(network_name=NETWORK_NAME,
                                         network_path=_network_path(t3),
                                         selection=selection)
        t3._record_pdep_network_identity(network_name=NETWORK_NAME,
                                         network_path=_network_path(t3),
                                         selection=selection)
        # ...but a conflicting one is refused. The network file changed on disk mid-iteration, and
        # the next reaction of the same network re-parsed it (main.py parses per reaction), so that
        # reaction's selection carries the NEW hash while the sidecar already holds the old one.
        # Note the conflict is now detected between two hashes each of which some selection actually
        # examined, rather than between one examined hash and one re-read at record time.
        with open(_network_path(t3), 'a') as f:
            f.write('# the network file changed mid-iteration\n')
        selection.network_source_hash = hash_file(path=_network_path(t3))
        with pytest.raises(ValueError, match='Conflicting'):
            t3._record_pdep_network_identity(network_name=NETWORK_NAME,
                                             network_path=_network_path(t3),
                                             selection=selection)

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
        write_ts_join_sidecar(arc_project_directory=t3.paths['ARC'], records=stale_records,
                              networks=_networks_for_records(t3, stale_records))
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

    def test_a_queued_record_freezes_the_sensitivity_evidence_that_justified_its_selection(self, tmp_path):
        """Test that the evidence which justified selecting a transition state is frozen onto its
        join record at queue time -- the whole point of this design: a restart cannot re-derive it
        from ``self.pdep_network_selections``, since that in-memory attribute is empty on a restart,
        so it must already live on the record itself.
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
        assert records[0].coefficient == 0.05
        assert records[0].delta_ln_k == 0.05

    def test_a_not_queued_record_still_freezes_the_sensitivity_evidence(self, tmp_path):
        """Test that a transition state recorded as NOT_QUEUED (here: missing adjacency lists) still
        carries the evidence that made it a candidate in the first place -- a partial capture must be
        able to tell, from the record alone, whether the missing transition state was in fact the
        dominant one, and that requires the evidence regardless of whether ARC ever saw the reaction.
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
        assert records[0].coefficient == 0.05
        assert records[0].delta_ln_k == 0.05

    def test_an_already_present_record_freezes_the_sensitivity_evidence(self, tmp_path):
        """Test that re-offering an already-known reaction still freezes the evidence onto its
        ALREADY_PRESENT record, not just the original QUEUED one.
        """
        t3 = _build_t3(tmp_path)
        network = parse_pdep_network_file(_network_path(t3))
        structures = _build_sa_dict(t3)['structures']
        kwargs = dict(network=network, structures=structures, network_name=NETWORK_NAME)

        t3.queue_pdep_transition_states(selection=_selection(['TS2']), **kwargs)
        t3.pdep_ts_join_records = list()  # a later iteration starts with an empty record set
        second = t3.queue_pdep_transition_states(selection=_selection(['TS2']), **kwargs)

        assert second[0].status == JOIN_STATUS_ALREADY_PRESENT
        assert second[0].coefficient == 0.05
        assert second[0].delta_ln_k == 0.05

    def test_the_dominant_entry_by_max_abs_coefficient_is_frozen_when_several_entries_share_one_label(
            self, tmp_path):
        """Test the aggregation rule: when several ``uncertain_path_reactions`` entries share one
        ``ts_label`` (several path reactions of the same network route through the same transition
        state), the single frozen ``(coefficient, delta_ln_k)`` pair must be the entry with the
        largest magnitude coefficient -- the one that actually justifies the selection, not merely
        whichever entry happened to be seen last or first.
        """
        t3 = _build_t3(tmp_path)
        network = parse_pdep_network_file(_network_path(t3))
        structures = _build_sa_dict(t3)['structures']
        selection = _selection_with_evidence([
            ('TS2', 0.02, 0.01),
            ('TS2', -0.09, 0.45),  # largest magnitude: this one must win
            ('TS2', 0.03, 0.015),
        ])

        records = t3.queue_pdep_transition_states(network=network,
                                                  selection=selection,
                                                  structures=structures,
                                                  network_name=NETWORK_NAME,
                                                  )

        assert len(records) == 1
        assert records[0].coefficient == -0.09
        assert records[0].delta_ln_k == 0.45

    def test_an_internal_invariant_raises_when_a_reported_label_has_no_matching_evidence(self, tmp_path,
                                                                                          monkeypatch):
        """Test the defensive internal-invariant guard: ``uncertain_ts_labels()`` is derived from
        ``uncertain_path_reactions`` in production ``PDepNetworkSelection``, so this mismatch cannot
        arise through ordinary use -- but if a future refactor ever broke that derivation silently, a
        transition state would be queued (or recorded) with no sensitivity evidence at all, which is
        exactly the defect this whole change exists to prevent. The guard is exercised directly here
        by monkeypatching ``uncertain_ts_labels`` to report a label absent from the evidence list.
        """
        t3 = _build_t3(tmp_path)
        network = parse_pdep_network_file(_network_path(t3))
        structures = _build_sa_dict(t3)['structures']
        selection = _selection(['TS2'])
        monkeypatch.setattr(selection, 'uncertain_ts_labels', lambda: ['TS2', 'TS_GHOST'])

        with pytest.raises(ValueError, match='TS_GHOST'):
            t3.queue_pdep_transition_states(network=network,
                                            selection=selection,
                                            structures=structures,
                                            network_name=NETWORK_NAME,
                                            )


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


def _write_energy_settings_fixture(arc_dir: str, statmech_subdir: str = 'kinetics') -> None:
    """Write a minimal, valid ARC energy-settings fixture -- ``calcs/statmech/<statmech_subdir>/
    input.py`` and ``output/output.yml`` -- under ``arc_dir``, exactly what
    ``t3.pdep.energy_settings.read_arc_energy_settings`` requires to freeze a complete
    energy-settings block. Mirrors ``tests/test_pdep/test_capture.py::_write_energy_settings_fixture``,
    including its ``useAtomCorrections = True`` + populated ``atomEnergies`` pairing: ARC computes
    ``useAtomCorrections`` as ``bool(model_chemistry or atom_energies)`` and records the correction
    values alongside it, so a fixture with the flag on but no values would be an ARC-impossible
    state that ``write_hybrid_network_input_file`` rightly refuses (its atom-energies guard) --
    leaving it out would make every capture-to-hybrid wiring test fail on that guard instead of on
    what it actually tests. ``useBondCorrections`` stays ``False`` so no cross-validation against
    ``output.yml`` is ever triggered by this fixture. Gated so it never overwrites an ``output.yml``
    a test has already written of its own (``output.yml`` is project-global, shared across statmech
    subdirs)."""
    statmech_dir = os.path.join(arc_dir, 'calcs', 'statmech', statmech_subdir)
    os.makedirs(statmech_dir, exist_ok=True)
    input_py_path = os.path.join(statmech_dir, 'input.py')
    if not os.path.isfile(input_py_path):
        with open(input_py_path, 'w') as f:
            f.write(
                "modelChemistry = 'CBS-QB3'\n\n"
                "useHinderedRotors = True\n\n"
                "useAtomCorrections = True\n\n"
                "atomEnergies = {'C': -37.844411, 'H': -0.499818, 'N': -54.581501, 'O': -75.062219}\n\n"
                "useBondCorrections = False\n"
            )
    output_dir = os.path.join(arc_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    output_yml_path = os.path.join(output_dir, 'output.yml')
    if not os.path.isfile(output_yml_path):
        with open(output_yml_path, 'w') as f:
            f.write('{}\n')


def _networks_for_records(t3, records: list, method: str = 'MSC') -> dict:
    """Build the sidecar networks block for the given records against ``t3``'s RMG pdep directory,
    writing a stub network file for any referenced network that has no real fixture file there.
    Hashes are computed with ``t3.pdep.cache.hash_file``, the same primitive production code
    records them with."""
    rmg_pdep_dir = t3.paths['RMG PDep']
    os.makedirs(rmg_pdep_dir, exist_ok=True)
    networks = dict()
    for network_id in sorted({record.network_id for record in records}):
        source_path = os.path.join(rmg_pdep_dir, f'{network_id}.py')
        if not os.path.isfile(source_path):
            with open(source_path, 'w') as f:
                f.write(f"# stub RMG network file\nnetwork(label='{network_id}')\n")
        networks[network_id] = {'source_path': source_path,
                                'source_sha256': hash_file(source_path),
                                'method': method}
    return networks


def _write_sidecar(t3, records: list, method: str = 'MSC') -> None:
    """Write the join sidecar for ``records`` with a matching networks block (see
    ``_networks_for_records``), exactly as the in-run path would have pre-ARC."""
    write_ts_join_sidecar(arc_project_directory=t3.paths['ARC'],
                          records=records,
                          networks=_networks_for_records(t3, records, method=method))


def _queue_usable_ts(t3, network_id='network4_2', network_ts_label='TS9') -> TSJoinRecord:
    """Queue one usable PDep transition state against ``t3``'s current ARC directory: write the
    join sidecar record T3 would have written pre-ARC (with its networks identity block), plus the
    converged ARC artifact and status entry ARC would have produced, so ``process_arc_run()``'s
    capture step has something real to discover and vendor. Also writes the ARC energy-settings
    fixture (``calcs/statmech/kinetics/input.py`` + ``output/output.yml``) that
    ``capture_ts_artifacts`` now requires -- without it, ``read_arc_energy_settings`` fails closed
    and every caller of this helper would raise."""
    arc_dir = t3.paths['ARC']
    label = arc_ts_label(network_id, network_ts_label)
    expected_path = expected_ts_artifact_path(arc_dir, label)
    record = TSJoinRecord(network_id=network_id,
                          network_ts_label=network_ts_label,
                          status=JOIN_STATUS_QUEUED,
                          arc_ts_label=label,
                          expected_artifact_path=expected_path,
                          reason='Queued to ARC.',
                          coefficient=0.05,
                          delta_ln_k=0.02,
                          )
    _write_ts_artifact(expected_path)
    _write_ts_status_yml(arc_dir, label, converged=True)
    _write_sidecar(t3, [record])
    _write_energy_settings_fixture(arc_dir)
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
        _write_energy_settings_fixture(arc_dir)
        missing_label = arc_ts_label('network4_2', 'TS7')
        records = [TSJoinRecord(network_id='network4_2',
                                network_ts_label='TS9',
                                status=JOIN_STATUS_QUEUED,
                                arc_ts_label=captured_label,
                                expected_artifact_path=captured_path,
                                reason='Queued to ARC.',
                                coefficient=0.05,
                                delta_ln_k=0.02,
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
        _write_sidecar(t3, records)

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
        _write_sidecar(t3, [record])

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
        _write_sidecar(t3, [already_present_record, not_queued_record])

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


def _queue_ts_set(t3, ts_specs: list) -> list:
    """Queue several PDep transition states of ``network4_2`` in ONE join sidecar (unlike
    ``_queue_usable_ts``, which writes a single-record sidecar), each with its own ARC artifact and
    a per-transition-state convergence verdict. ``ts_specs`` entries are
    ``(network_ts_label, converged)`` pairs: ``converged=False`` yields an artifact that discovery
    classifies UNUSABLE (present but explicitly unconverged), which is what the partial-hybrid
    refusal tests need to exist alongside a usable sibling."""
    arc_dir = t3.paths['ARC']
    records = list()
    for network_ts_label, converged in ts_specs:
        label = arc_ts_label('network4_2', network_ts_label)
        expected_path = expected_ts_artifact_path(arc_dir, label)
        records.append(TSJoinRecord(network_id='network4_2',
                                    network_ts_label=network_ts_label,
                                    status=JOIN_STATUS_QUEUED,
                                    arc_ts_label=label,
                                    expected_artifact_path=expected_path,
                                    reason='Queued to ARC.',
                                    coefficient=0.05,
                                    delta_ln_k=0.02,
                                    ))
        _write_ts_artifact(expected_path)
        _write_ts_status_yml(arc_dir, label, converged=converged)
    _write_sidecar(t3, records)
    _write_energy_settings_fixture(arc_dir)
    return records


class TestProcessArcRunHybridWiring(object):
    """Acceptance tests for the hybrid-network wiring added to ``process_arc_run()``: one hybrid
    Arkane P-dep input file per network, built exclusively from the durable capture (never from the
    live ARC project directory or the live RMG pdep directory), with only USABLE captured
    transition states switched to QM/RRKM and every other transition state left as RMG/ILT.
    Also covers the finalization-marker VERSIONING that keeps a marker written by an older,
    pre-hybrid finalization format from being silently treated as current.
    """

    def test_process_arc_run_writes_a_hybrid_network_input_from_the_capture(self, tmp_path):
        """The ordinary path: after capture, process_arc_run() must write one hybrid Arkane input
        for network4_2 under 'PDep hybrid', with the captured-usable TS9 switched to QM/RRKM (its
        vendored artifact referenced and present) and the other transition states left ILT."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)  # network4_2 / TS9, against the REAL network4_2.py fixture

        t3.process_arc_run()

        hybrid_network_dir = os.path.join(t3.paths['PDep hybrid'], 'network4_2')
        hybrid_input_path = os.path.join(hybrid_network_dir, 'input.py')
        assert os.path.isfile(hybrid_input_path)
        with open(hybrid_input_path, 'r') as f:
            written = f.read()
        # TS9 must be QM/RRKM: its transitionState is now spliced inline and vibration-only (its
        # E0, frequencies and imaginary mode), and the network references no external qm/ file ...
        assert "label = 'TS9'" in written
        assert "frequency = (-1800.0,'cm^-1')" in written
        assert 'qm/' not in written
        # ... while a never-selected transition state survives as ILT rather than being dropped.
        assert "'TS1'" in written
        # The frozen energy reference must have been injected from the capture's manifest.
        assert 'useAtomCorrections = True' in written
        assert os.path.isfile(t3.paths['ARC finalization marker'])

    def test_the_capture_alone_drives_the_production_hybrid_write_with_arc_absent(self, tmp_path):
        """THE acceptance property this whole capture layer exists for, on the PRODUCTION path:
        with the ARC project directory and the RMG pdep directory both deleted entirely, the
        production hybrid writer must still rebuild every hybrid network input from the capture
        alone -- vendored network source, frozen ME method, frozen energy settings, and vendored
        QM artifacts. Reaching for either deleted directory fails loudly here."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()

        # Non-vacuity: the capture must actually hold a usable TS9 artifact before anything is
        # deleted -- otherwise the assertions below would exercise an empty QM set and prove
        # nothing about the capture driving a hybrid write.
        manifest = read_yaml_file(path=os.path.join(t3.paths['PDep capture'], CAPTURE_MANIFEST_FILE_NAME))
        usable_labels = [entry['network_ts_label'] for entry in manifest['transition_states']
                         if entry['status'] == ARTIFACT_STATUS_USABLE]
        assert usable_labels == ['TS9']

        shutil.rmtree(t3.paths['ARC'])
        shutil.rmtree(t3.paths['RMG PDep'])
        shutil.rmtree(t3.paths['PDep hybrid'])  # force a full rebuild, not a stale leftover
        assert not os.path.exists(t3.paths['ARC'])
        assert not os.path.exists(t3.paths['RMG PDep'])

        results = t3._write_pdep_hybrid_network_inputs()

        assert sorted(results.keys()) == ['network4_2']
        result = results['network4_2']
        # Both sides of the hybrid must be non-empty: an empty QM set would mean the capture drove
        # nothing, and an empty ILT set would mean this is not a hybrid at all.
        assert result.qm_ts_labels == ('TS9',)
        assert len(result.ilt_ts_labels) > 0
        assert 'TS9' not in result.ilt_ts_labels
        assert os.path.isfile(result.dest_path)
        # Real containment, not a prefix check: str.startswith() is exactly what
        # T3._confine_to_pdep_hybrid_root deliberately avoids (a sibling directory like
        # '<root>-old' shares '<root>' as a string prefix but is NOT contained in it), so a test
        # asserting containment via startswith would pass on an escape the production code
        # refuses. os.path.commonpath on realpath'd paths is what the production code itself uses.
        resolved_hybrid_root = os.path.realpath(t3.paths['PDep hybrid'])
        resolved_dest_path = os.path.realpath(result.dest_path)
        assert os.path.commonpath([resolved_hybrid_root, resolved_dest_path]) == resolved_hybrid_root
        with open(result.dest_path, 'r') as f:
            written = f.read()
        assert "label = 'TS9'" in written
        assert 'qm/' not in written
        assert 'useAtomCorrections = True' in written

    def test_hybrid_write_refuses_when_no_capture_exists(self, tmp_path):
        """Fail-closed: the production hybrid writer must never treat a missing capture as 'nothing
        to do' -- called without a capture on disk it must raise (via verify_capture), because a
        silently-skipped hybrid write is indistinguishable from a successful empty one."""
        t3 = _build_t3(tmp_path)
        assert not os.path.exists(t3.paths['PDep capture'])

        with pytest.raises(ValueError, match='no capture manifest'):
            t3._write_pdep_hybrid_network_inputs()

    def test_no_hybrid_is_written_for_a_verified_zero_artifact_capture(self, tmp_path):
        """A capture that verified but holds zero captured artifacts (here: the only queued
        transition state came back explicitly unconverged, i.e. UNUSABLE and never vendored) has
        nothing to hybridize: finalization must complete (marker written) with no hybrid directory
        materialized, and the skip must be the outcome of a VERIFIED capture, not of a missing
        one."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_ts_set(t3, [('TS9', False)])

        t3.process_arc_run()

        # Non-vacuity: the capture exists, verified, and its single record is UNUSABLE.
        manifest = read_yaml_file(path=os.path.join(t3.paths['PDep capture'], CAPTURE_MANIFEST_FILE_NAME))
        assert [entry['status'] for entry in manifest['transition_states']] == ['unusable']
        assert not os.path.exists(t3.paths['PDep hybrid'])
        assert os.path.isfile(t3.paths['ARC finalization marker'])

    def test_a_partially_usable_selected_set_refuses_the_network_hybrid(self, tmp_path):
        """Strict-by-default gate (evaluate_pdep_hybrid): when only SOME of a network's selected
        transition states came back usable (TS9 converged, TS10 explicitly unconverged), the
        network must NOT get a half-QM hybrid input -- a hybrid missing QM for one of its most
        uncertain transition states is a degradation of what was asked for, accepted only by
        explicit opt-in. Finalization itself still completes: the refusal is a recorded decision
        about this network, not a crash."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_ts_set(t3, [('TS9', True), ('TS10', False)])

        t3.process_arc_run()

        # Non-vacuity: a usable TS9 artifact WAS captured, so the absent hybrid below is the
        # gate's doing, not a consequence of having nothing to write.
        manifest = read_yaml_file(path=os.path.join(t3.paths['PDep capture'], CAPTURE_MANIFEST_FILE_NAME))
        statuses = {entry['network_ts_label']: entry['status'] for entry in manifest['transition_states']}
        assert statuses == {'TS9': ARTIFACT_STATUS_USABLE, 'TS10': 'unusable'}
        assert not os.path.isfile(os.path.join(t3.paths['PDep hybrid'], 'network4_2', 'input.py'))
        assert os.path.isfile(t3.paths['ARC finalization marker'])

    def test_process_arc_run_propagates_a_hybrid_write_failure_and_leaves_marker_absent(self, tmp_path,
                                                                                        monkeypatch):
        """Fail-closed: if the hybrid write raises, process_arc_run() must propagate the exception
        rather than swallow it, and the finalization marker must be left absent so a subsequent
        restart() retries finalization instead of treating the failed run as done."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)

        def _boom(**kwargs):
            raise ValueError('simulated hybrid write failure')

        monkeypatch.setattr(t3_main, 'write_hybrid_network_input_file', _boom)

        with pytest.raises(ValueError, match='simulated hybrid write failure'):
            t3.process_arc_run()

        assert not os.path.isfile(t3.paths['ARC finalization marker'])

    def test_an_unversioned_legacy_marker_reads_as_not_finalized(self, tmp_path):
        """Marker versioning, backward direction: a marker written by the pre-versioning format
        (the bare legacy text, no version token) must read as NOT finalized -- that marker was
        written by a finalization that never produced hybrid network inputs, and trusting it would
        skip the hybrid write forever. The legacy text is pinned literally here on purpose: writing
        it via the (now versioned) constant would make this test vacuously compare the marker with
        itself."""
        t3 = _build_t3(tmp_path)
        marker_path = t3.paths['ARC finalization marker']
        os.makedirs(os.path.dirname(marker_path), exist_ok=True)
        with open(marker_path, 'w') as f:
            f.write('ARC finalization completed at 2026-01-01T00:00:00\n')

        assert t3.check_arc_finalization_complete() is False

    def test_a_future_version_marker_reads_as_not_finalized(self, tmp_path):
        """Marker versioning, forward direction: a marker carrying a NEWER version token than this
        code writes must also read as not finalized (fail closed) -- this code cannot know what a
        future finalization format guaranteed, so it must redo finalization rather than assume."""
        t3 = _build_t3(tmp_path)
        marker_path = t3.paths['ARC finalization marker']
        os.makedirs(os.path.dirname(marker_path), exist_ok=True)
        with open(marker_path, 'w') as f:
            f.write('ARC finalization completed (v99) at 2026-01-01T00:00:00\n')

        assert t3.check_arc_finalization_complete() is False

    def test_a_marker_written_by_this_code_reads_as_finalized(self, tmp_path):
        """The positive control for the two versioning tests above: a marker written by
        _mark_arc_finalization_complete() itself must read back as finalized -- without this, the
        versioning checks could 'pass' by rejecting every marker including our own."""
        t3 = _build_t3(tmp_path)

        t3._mark_arc_finalization_complete()

        assert t3.check_arc_finalization_complete() is True

    def test_process_arc_run_redoes_finalization_over_a_legacy_marker(self, tmp_path):
        """End-to-end consequence of marker versioning: with a legacy (unversioned) marker on disk
        and real queued work present, process_arc_run() must NOT take the already-finalized
        short-circuit -- it must redo finalization (capture + hybrid write) and leave a
        current-version marker behind."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        marker_path = t3.paths['ARC finalization marker']
        os.makedirs(os.path.dirname(marker_path), exist_ok=True)
        with open(marker_path, 'w') as f:
            f.write('ARC finalization completed at 2026-01-01T00:00:00\n')

        t3.process_arc_run()

        assert os.path.isdir(t3.paths['PDep capture'])
        assert os.path.isfile(os.path.join(t3.paths['PDep hybrid'], 'network4_2', 'input.py'))
        assert t3.check_arc_finalization_complete() is True


class TestWritePdepHybridNetworkInputsPruning(object):
    """Acceptance tests for the stale-output pruning and network_id path confinement added to
    ``_write_pdep_hybrid_network_inputs()``. The method's documented guarantee is that a network
    this pass REFUSES (or does not even see) stays entirely RMG/ILT -- which means it must have no
    hybrid input file at all. Without pruning, a directory this method wrote on a PRIOR call (a
    prior run of the same iteration, e.g. via the legacy-marker redo path covered above) survives a
    later call that refuses or drops that network, so a downstream consumer can read a hybrid QM
    input for a network T3 no longer stands behind.
    """

    @staticmethod
    def _manifest_path(t3):
        return os.path.join(t3.paths['PDep capture'], CAPTURE_MANIFEST_FILE_NAME)

    def test_a_refused_networks_stale_hybrid_input_is_removed_on_rerun(self, tmp_path):
        """A network accepted (and written) on one call must have its hybrid input REMOVED by a
        later call that refuses it -- here, because a second, unconverged transition state
        (TS10) was added to the manifest, tipping the strict-by-default gate to refuse the whole
        network. Before the fix, the ``continue`` on a refused network left the previous call's
        ``input.py`` on disk untouched."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()
        hybrid_input_path = os.path.join(t3.paths['PDep hybrid'], 'network4_2', 'input.py')
        assert os.path.isfile(hybrid_input_path), 'setup: the first pass must have written the hybrid input'

        manifest = read_yaml_file(path=self._manifest_path(t3))
        first_entry = manifest['transition_states'][0]
        stale_entry = dict(first_entry)
        stale_entry.update(network_ts_label='TS10', arc_ts_label='T3PDep_network4_2_TS10', status='unusable',
                            converged=False, reason='simulated unconverged transition state',
                            captured_artifact_path=None, captured_artifact_sha256=None, captured_log_paths=dict(),
                            captured_log_sha256=dict())
        manifest['transition_states'].append(stale_entry)
        save_yaml_file(path=self._manifest_path(t3), content=manifest)

        results = t3._write_pdep_hybrid_network_inputs()

        assert 'network4_2' not in results
        assert not os.path.isfile(hybrid_input_path), \
            'a network refused on THIS pass must not leave a hybrid input written by a prior pass'
        assert not os.path.isdir(os.path.dirname(hybrid_input_path)), \
            'the whole stale per-network directory must be gone, not merely its input.py'

    def test_a_stale_network_absent_from_the_current_capture_is_removed(self, tmp_path):
        """A per-network output directory left by a prior run, for a network that is not merely
        refused but entirely ABSENT from the current capture's manifest, must also be removed --
        nothing in the transition_states loop ever revisits a network it never iterates over, so
        this case needs pruning against what is actually on disk under 'PDep hybrid', not against
        the manifest."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()
        hybrid_input_path = os.path.join(t3.paths['PDep hybrid'], 'network4_2', 'input.py')
        assert os.path.isfile(hybrid_input_path), 'setup: the first pass must have written the hybrid input'

        ghost_dir = os.path.join(t3.paths['PDep hybrid'], 'ghost_network')
        os.makedirs(ghost_dir)
        with open(os.path.join(ghost_dir, 'input.py'), 'w') as f:
            f.write('# stale hybrid input from a network absent from the current capture\n')

        results = t3._write_pdep_hybrid_network_inputs()

        assert sorted(results.keys()) == ['network4_2']
        assert not os.path.exists(ghost_dir), 'a network absent from the current capture must not survive'
        # No over-deletion regression: the network THIS pass accepts must still be written.
        assert os.path.isfile(hybrid_input_path)

    def test_zero_captured_artifacts_clears_stale_hybrid_outputs(self, tmp_path):
        """The zero-captured-artifact early return is a positive verdict ('nothing to
        hybridize'), not an excuse to leave a prior pass's outputs in place: it must clear the
        'PDep hybrid' tree before returning ``dict()``, exactly like the ordinary path does."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()
        hybrid_input_path = os.path.join(t3.paths['PDep hybrid'], 'network4_2', 'input.py')
        assert os.path.isfile(hybrid_input_path), 'setup: the first pass must have written the hybrid input'

        manifest = read_yaml_file(path=self._manifest_path(t3))
        only_entry = manifest['transition_states'][0]
        only_entry.update(status='unusable', converged=False, reason='simulated unconverged transition state',
                          captured_artifact_path=None, captured_artifact_sha256=None, captured_log_paths=dict(),
                          captured_log_sha256=dict())
        save_yaml_file(path=self._manifest_path(t3), content=manifest)

        results = t3._write_pdep_hybrid_network_inputs()

        assert results == dict()
        assert not os.path.isfile(hybrid_input_path)
        assert not os.path.isdir(os.path.join(t3.paths['PDep hybrid'], 'network4_2'))

    def test_a_manifest_network_id_containing_dotdot_is_refused(self, tmp_path):
        """A corrupted or crafted manifest whose network_id contains '..' must never be allowed to
        escape 'PDep hybrid' via a raw os.path.join(): the method must refuse it with a ValueError
        rather than write outside its own output root."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()

        manifest = read_yaml_file(path=self._manifest_path(t3))
        for entry in manifest['transition_states']:
            entry['network_id'] = '../escaped'
        manifest['networks'] = {'../escaped': manifest['networks'].pop('network4_2')}
        save_yaml_file(path=self._manifest_path(t3), content=manifest)

        escaped_dir = os.path.realpath(os.path.join(t3.paths['PDep hybrid'], '..', 'escaped'))
        assert not os.path.exists(escaped_dir), 'setup: nothing must already exist at the escape target'

        # '../escaped' contains a path separator, so it is now refused by the single-safe-filename-
        # component check (added for P2 #4) before the realpath+commonpath "resolves outside" check
        # below it is ever reached -- either message proves the escape is refused.
        with pytest.raises(ValueError, match='resolves outside|filename component'):
            t3._write_pdep_hybrid_network_inputs()

        assert not os.path.exists(escaped_dir), 'the manifest must never be able to write outside PDep hybrid'

    def test_an_escaping_network_id_is_refused_before_anything_is_pruned(self, tmp_path):
        """Refusing a manifest is not enough: the refusal must come BEFORE the stale-output prune
        touches anything.

        Pruning is destructive and is driven by the manifest's own network_ids -- so running it
        against a manifest that has not been validated yet lets a corrupt or crafted manifest
        delete a legitimate prior pass's hybrid outputs on its way to being rejected. The method
        still fails closed either way, but 'validate the whole manifest, then destroy' is the
        invariant worth holding: a manifest T3 refuses to act on must not be able to act on the
        output tree at all."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()
        surviving_input_path = os.path.join(t3.paths['PDep hybrid'], 'network4_2', 'input.py')
        assert os.path.isfile(surviving_input_path), 'setup: the first pass must write a hybrid input'

        manifest = read_yaml_file(path=self._manifest_path(t3))
        for entry in manifest['transition_states']:
            entry['network_id'] = '../escaped'
        manifest['networks'] = {'../escaped': manifest['networks'].pop('network4_2')}
        save_yaml_file(path=self._manifest_path(t3), content=manifest)

        # See the sibling dotdot test above: '../escaped' is now refused by the single-safe-
        # filename-component check before the "resolves outside" check is reached.
        with pytest.raises(ValueError, match='resolves outside|filename component'):
            t3._write_pdep_hybrid_network_inputs()

        assert os.path.isfile(surviving_input_path), \
            'a manifest that is refused must not have pruned the previous pass\'s outputs first'

    def test_a_manifest_network_id_that_is_absolute_is_refused(self, tmp_path):
        """A manifest network_id that is an absolute path must also be refused: os.path.join()
        with an absolute second argument discards the first entirely, so this is a second,
        independent way the same defect is reachable."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()

        evil_abs = str(tmp_path / 'evil_abs')
        manifest = read_yaml_file(path=self._manifest_path(t3))
        for entry in manifest['transition_states']:
            entry['network_id'] = evil_abs
        manifest['networks'] = {evil_abs: manifest['networks'].pop('network4_2')}
        save_yaml_file(path=self._manifest_path(t3), content=manifest)

        assert not os.path.exists(evil_abs), 'setup: nothing must already exist at the escape target'

        # See the sibling dotdot test above: an absolute network_id is now refused by the
        # single-safe-filename-component check before the "resolves outside" check is reached.
        with pytest.raises(ValueError, match='resolves outside|filename component'):
            t3._write_pdep_hybrid_network_inputs()

        assert not os.path.exists(evil_abs), 'the manifest must never be able to write outside PDep hybrid'

    def test_a_manifest_network_id_containing_a_path_separator_is_refused(self, tmp_path):
        """A network_id like 'a/b' passes ``_confine_to_pdep_hybrid_root``'s realpath+commonpath
        check (it still resolves under the hybrid root), but ``_prune_stale_pdep_hybrid_outputs``
        compares TOP-LEVEL ``os.listdir()`` names against the raw accepted ids, so 'a/b' and a
        top-level 'a' entry would collide there. A network_id must be a single safe filename
        component -- no path separator, not '.'/'..'  -- so this collision can never arise."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()

        evil_id = 'a/b'
        manifest = read_yaml_file(path=self._manifest_path(t3))
        for entry in manifest['transition_states']:
            entry['network_id'] = evil_id
        manifest['networks'] = {evil_id: manifest['networks'].pop('network4_2')}
        save_yaml_file(path=self._manifest_path(t3), content=manifest)

        escape_target = os.path.realpath(os.path.join(t3.paths['PDep hybrid'], evil_id))
        assert not os.path.exists(escape_target), 'setup: nothing must already exist at the escape target'

        with pytest.raises(ValueError, match='filename component'):
            t3._write_pdep_hybrid_network_inputs()

        assert not os.path.exists(escape_target), \
            "a network_id containing a path separator must never be written under 'PDep hybrid'"

    def test_writer_consumes_only_the_records_verify_capture_already_verified(self, tmp_path, monkeypatch):
        """The writer used to re-read the manifest itself (via ``read_yaml_file(verified.manifest_
        path)``) AFTER ``verify_capture`` had already read and validated it, and rebuild
        TSArtifactRecords from that second, untrusted read. Proves that re-read is gone: the
        manifest is deleted entirely right after ``verify_capture`` returns (from inside a wrapper
        substituted for ``t3_main.verify_capture``), so any second read of it would raise. The
        writer must still succeed, using only the records ``verify_capture`` already handed back."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()
        hybrid_input_path = os.path.join(t3.paths['PDep hybrid'], 'network4_2', 'input.py')
        assert os.path.isfile(hybrid_input_path), 'setup: first pass wrote the hybrid input'
        os.remove(hybrid_input_path)
        assert not os.path.isfile(hybrid_input_path), 'setup: removed so a second, real write can be observed'

        real_verify_capture = t3_main.verify_capture

        def verify_then_delete_manifest(capture_dir):
            result = real_verify_capture(capture_dir)
            os.remove(result.manifest_path)
            return result

        monkeypatch.setattr(t3_main, 'verify_capture', verify_then_delete_manifest)

        # Must NOT raise: a second read of the (now-deleted) manifest would raise
        # arc.common.read_yaml_file's "could not find the YAML file" InputError.
        t3._write_pdep_hybrid_network_inputs()

        assert os.path.isfile(hybrid_input_path), \
            'the writer must rebuild the hybrid input from the records verify_capture already ' \
            'returned, never by re-reading the manifest a second time'

    def test_an_accepted_networks_destination_that_is_a_plain_file_is_refused(self, tmp_path):
        """An ACCEPTED network's destination under 'PDep hybrid' must itself be a plain directory
        (or absent) before anything is pruned. If a prior pass (or anything else) left a plain FILE
        at that path, the writer would later try to create this network's input.py underneath it
        and blow up -- but only AFTER the same call had already deleted other, unrelated stale
        directories, leaving a half-pruned tree behind. This must be caught during the preflight,
        before any deletion, exactly like an unacceptable stale entry is."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()
        hybrid_root = t3.paths['PDep hybrid']
        network_dir = os.path.join(hybrid_root, 'network4_2')
        hybrid_input_path = os.path.join(network_dir, 'input.py')
        assert os.path.isfile(hybrid_input_path), 'setup: the first pass must have written the hybrid input'

        # A genuinely stale, otherwise-legitimate directory that a same-pass validate-and-delete
        # implementation WOULD remove before ever reaching the accepted network's corrupted
        # destination below (network4_2 sorts after this name).
        stale_dir = os.path.join(hybrid_root, 'aaa_stale_network')
        os.makedirs(stale_dir)
        with open(os.path.join(stale_dir, 'input.py'), 'w') as f:
            f.write('# stale hybrid input from a prior pass, no longer among the accepted networks\n')

        # Replace the accepted network's own destination with a plain FILE -- something T3 itself
        # would never have written there.
        shutil.rmtree(network_dir)
        with open(network_dir, 'w') as f:
            f.write('not a directory T3 itself would have written\n')

        with pytest.raises(ValueError, match='not a plain directory'):
            t3._write_pdep_hybrid_network_inputs()

        assert os.path.isdir(stale_dir), (
            'a half-pruned tree must be impossible: the pre-existing stale directory must still be '
            'intact after the accepted network destination check raises'
        )

    def test_a_pre_existing_stale_directory_survives_when_a_later_entry_makes_pruning_raise(self, tmp_path,
                                                                                            monkeypatch):
        """A half-pruned tree must be impossible: ``_prune_stale_pdep_hybrid_outputs`` must judge
        every entry under the hybrid root acceptable BEFORE deleting any of them, not
        validate-and-delete in the same ``os.listdir()`` pass. Forces a deterministic iteration
        order (a stale, otherwise-legitimate directory sorts before an unacceptable top-level file)
        so that, pre-fix, the stale directory is rmtree'd before the later entry is ever reached and
        raises -- leaving a half-pruned tree behind."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()
        hybrid_root = t3.paths['PDep hybrid']

        # Sorts before 'network4_2' (the accepted network) and, more importantly, before the bad
        # entry below -- so a same-pass validate-and-delete loop reaches (and removes) it first.
        stale_dir = os.path.join(hybrid_root, 'aaa_stale_network')
        os.makedirs(stale_dir)
        with open(os.path.join(stale_dir, 'input.py'), 'w') as f:
            f.write('# stale hybrid input from a prior pass, no longer among the accepted networks\n')

        # Sorts after both of the above, so it is only reached once the stale directory has
        # already been removed by a same-pass implementation.
        bad_entry = os.path.join(hybrid_root, 'zzz_bad_entry')
        with open(bad_entry, 'w') as f:
            f.write('not a directory T3 itself would have written\n')

        real_listdir = os.listdir

        def sorted_listdir(path):
            names = real_listdir(path)
            if os.path.realpath(path) == os.path.realpath(hybrid_root):
                return sorted(names)
            return names

        monkeypatch.setattr(os, 'listdir', sorted_listdir)

        assert os.path.isdir(stale_dir), 'setup: the stale directory exists before pruning runs'

        with pytest.raises(ValueError, match='not a plain directory'):
            t3._write_pdep_hybrid_network_inputs()

        assert os.path.isdir(stale_dir), (
            'a half-pruned tree must be impossible: this pre-existing stale directory must still '
            'be intact after pruning raises on a later, unacceptable entry'
        )

    def test_a_usable_status_with_a_null_captured_artifact_path_is_refused_through_process_arc_run(self, tmp_path, monkeypatch):
        """P1 #2, driven through the real, end-to-end path (not just verify_capture in isolation):
        a manifest hand-edited (or corrupted) to claim status: usable with a null
        captured_artifact_path must be refused by process_arc_run(), and must leave the hybrid
        output tree unpruned and the finalization marker absent -- exactly as any other refused
        capture does."""
        t3 = _build_t3(tmp_path)
        _write_arc_info(t3)
        _queue_usable_ts(t3)
        t3.process_arc_run()
        hybrid_input_path = os.path.join(t3.paths['PDep hybrid'], 'network4_2', 'input.py')
        assert os.path.isfile(hybrid_input_path), 'setup: first pass wrote the hybrid input'
        assert t3.check_arc_finalization_complete(), 'setup: first pass finalized cleanly'

        t3._clear_arc_finalization_marker()
        manifest = read_yaml_file(path=self._manifest_path(t3))
        entry = manifest['transition_states'][0]
        assert entry['status'] == ARTIFACT_STATUS_USABLE
        entry['captured_artifact_path'] = None
        entry['captured_artifact_sha256'] = None
        save_yaml_file(path=self._manifest_path(t3), content=manifest)

        # _capture_pdep_ts_artifacts() would otherwise treat this corrupted manifest as "nothing
        # valid to reuse" (verify_capture raises -> _pdep_capture_is_authoritative catches that and
        # returns False), and silently HEAL the corruption by recapturing a fresh, valid manifest
        # from the still-live ARC project directory before process_arc_run ever reaches the writer
        # -- which would make this test pass for the wrong reason (recapture, not the P1 #2 guard).
        # Forcing the "an authoritative capture already exists" branch here is what actually
        # exercises _write_pdep_hybrid_network_inputs() -> verify_capture() against the hand-edited
        # manifest, which is the real target of this test.
        monkeypatch.setattr(t3, '_pdep_capture_is_authoritative', lambda capture_dir, join_records: True)

        with pytest.raises(ValueError, match='usable'):
            t3.process_arc_run()

        assert os.path.isfile(hybrid_input_path), \
            'a refused manifest must not prune the still-valid hybrid output from the prior pass'
        assert not t3.check_arc_finalization_complete(), \
            'the finalization marker must not be (re-)written when process_arc_run raises'
