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
from t3.pdep.parser import parse_pdep_network_file
from t3.pdep.selector import (CACHE_STATUS_CACHED_VALID,
                              CACHE_STATUS_GENERATED,
                              select_from_sa_dict,
                              )
from arc.common import save_yaml_file

from tests.test_pdep._wiring_helpers import (NETWORK_NAME,
                                             NETWORK_REACTION_STR,
                                             build_pdep_rxns_to_explore as _build_pdep_rxns_to_explore,
                                             build_sa_dict as _build_sa_dict,
                                             build_t3 as _build_t3,
                                             candidate_sa_path as _candidate_sa_path,
                                             network_path as _network_path,
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
