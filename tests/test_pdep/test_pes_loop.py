"""Tests for t3.pdep.pes_loop."""

import os

import pytest

from t3.pdep.explorer.result import (EXPLORATION_STATUS_FAILED, EXPLORATION_STATUS_SUCCEEDED,
                                     PDepExplorationResult)
from t3.pdep.parser import PDepNetwork, PDepPathReaction
from t3.pdep.pes_loop import (PES_LOOP_CONVERGED, PES_LOOP_MAX_ROUNDS, PES_LOOP_NO_CANDIDATES,
                              PES_LOOP_STALLED, PESLoopResult, RoundRecord, run_pes_loop)
from t3.schema import PESLoopConfig


@pytest.fixture
def config():
    """A minimal config pointing at a network fixture.

    ``bath_gas`` is required by ``PDepExplorerConfig`` (a config without one is a guaranteed
    write-time failure -- see t3/pdep/explorer/config.py, issue #183), so it must be non-empty
    here even though every test below stubs ``explore_pdep_network`` and never actually reaches
    that construction's use.
    """
    return PESLoopConfig(pes={'network': '/abs/network1_1.py', 'source': ['HOCHO'],
                              'bath_gas': {'He': 1.0}},
                         termination={'max_rounds': 3})


def _stub_explorer(monkeypatch, families, fail=False):
    """
    Patch ``explore_pdep_network``, ``parse_pdep_network_file``, and ``draw_pes_diagram`` at the
    names bound in ``t3.pdep.pes_loop``, so ``run_pes_loop`` can be tested without Arkane.

    Args:
        monkeypatch: pytest's monkeypatch fixture.
        families (list): One path reaction is fabricated per entry, carrying that family in its
            kinetics comment (recognized by ``t3.pdep.barrierless.classify_barrierless``).
        fail (bool): If True, the stubbed explorer reports a failed exploration instead of a
            succeeded one.
    """
    network_path = '/abs/network1_1.py'
    network_id = 'network1_1'

    if fail:
        result = PDepExplorationResult(network_id=network_id, status=EXPLORATION_STATUS_FAILED,
                                       reasons=('the explorer adapter crashed',))
    else:
        result = PDepExplorationResult(network_id=network_id, status=EXPLORATION_STATUS_SUCCEEDED,
                                       network_paths=(network_path,))

    path_reactions = tuple(
        PDepPathReaction(label=f'reaction{i}', reactants=('A',), products=('B',),
                         transition_state=f'TS{i}', kinetics_type='Arrhenius',
                         kinetics_comment=f'family: {family}')
        for i, family in enumerate(families))
    network = PDepNetwork(network_id=network_id, path=network_path, label=None,
                          species_labels=('A', 'B'),
                          transition_state_labels=tuple(f'TS{i}' for i in range(len(families))),
                          path_reactions=path_reactions, isomers=(), reactant_channels=(),
                          product_channels=())

    def _fake_explore(*args, **kwargs):
        return result

    def _fake_parse(*args, **kwargs):
        return network

    def _fake_draw(network_path, output_path):
        open(output_path, 'w').close()

    monkeypatch.setattr('t3.pdep.pes_loop.explore_pdep_network', _fake_explore)
    monkeypatch.setattr('t3.pdep.pes_loop.parse_pdep_network_file', _fake_parse)
    monkeypatch.setattr('t3.pdep.pes_loop.draw_pes_diagram', _fake_draw)


class TestRunPESLoop(object):

    def test_no_candidates_terminates_immediately(self, tmp_path, monkeypatch, config):
        """A network whose every channel is barrierless has nothing to compute: one explore, one
        diagram, done -- not three empty rounds."""
        _stub_explorer(monkeypatch, families=['R_Recombination', 'R_Recombination'])
        result = run_pes_loop(config, project_directory=str(tmp_path))
        assert result.status == PES_LOOP_NO_CANDIDATES
        assert len(result.rounds) == 1
        assert result.rounds[0].queued_ts_labels == ()
        assert len(result.rounds[0].skipped) == 2

    def test_max_rounds_is_honoured(self, tmp_path, monkeypatch, config):
        """A runner that never converges anything must not loop forever."""
        config = PESLoopConfig(pes={'network': '/abs/network1_1.py', 'source': ['HOCHO'],
                                    'bath_gas': {'He': 1.0}},
                               termination={'max_rounds': 3, 'stop_when_no_new_ts': False})
        _stub_explorer(monkeypatch, families=['1,2_Insertion_CO'])
        result = run_pes_loop(config, project_directory=str(tmp_path),
                              qm_runner=lambda candidates, paths, cfg, network_id: frozenset())
        assert result.status == PES_LOOP_MAX_ROUNDS
        assert len(result.rounds) == 3

    def test_a_round_that_stalls_stops_early(self, tmp_path, monkeypatch, config):
        """With the default stop_when_no_new_ts=True, a runner that converges nothing must stop
        after round 1 rather than spend the rest of the round budget -- this is what distinguishes
        'stalled' from 'max_rounds' (ruling C3)."""
        _stub_explorer(monkeypatch, families=['1,2_Insertion_CO'])
        result = run_pes_loop(config, project_directory=str(tmp_path),
                              qm_runner=lambda candidates, paths, cfg, network_id: frozenset())
        assert result.status == PES_LOOP_STALLED
        assert len(result.rounds) == 1

    def test_converges_when_every_candidate_is_computed(self, tmp_path, monkeypatch, config):
        """Once QM exists for every barriered channel, the loop stops on its own."""
        _stub_explorer(monkeypatch, families=['1,2_Insertion_CO'])
        result = run_pes_loop(config, project_directory=str(tmp_path),
                              qm_runner=lambda candidates, paths, cfg, network_id:
                                  frozenset(c.ts_label for c in candidates))
        assert result.status == PES_LOOP_CONVERGED
        assert len(result.rounds) == 2

    def test_a_diagram_is_drawn_every_round(self, tmp_path, monkeypatch, config):
        """The whole point: a picture per round, each from that round's own network."""
        _stub_explorer(monkeypatch, families=['1,2_Insertion_CO'])
        result = run_pes_loop(config, project_directory=str(tmp_path),
                              qm_runner=lambda candidates, paths, cfg, network_id:
                                  frozenset(c.ts_label for c in candidates))
        for record in result.rounds:
            assert record.diagram_path is not None
            assert os.path.basename(record.diagram_path) == 'pes_diagram.png'

    def test_barrierless_channels_are_never_queued(self, tmp_path, monkeypatch, config):
        """The F21 gate, asserted end to end rather than only at the unit level."""
        queued = []
        _stub_explorer(monkeypatch, families=['R_Recombination', '1,2_Insertion_CO'])

        def _runner(candidates, paths, cfg, network_id):
            queued.extend(c.family for c in candidates)
            return frozenset(c.ts_label for c in candidates)

        run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert 'R_Recombination' not in queued

    def test_no_qm_runner_explores_and_draws_only(self, tmp_path, monkeypatch, config):
        """The honest behaviour without a configured runner: one round, a diagram, no crash."""
        _stub_explorer(monkeypatch, families=['1,2_Insertion_CO'])
        result = run_pes_loop(config, project_directory=str(tmp_path))
        assert len(result.rounds) == 1
        assert result.rounds[0].diagram_path is not None

    def test_a_failed_explore_stops_the_loop_and_says_why(self, tmp_path, monkeypatch, config):
        """A round that cannot explore must not be papered over as convergence."""
        _stub_explorer(monkeypatch, families=['1,2_Insertion_CO'], fail=True)
        result = run_pes_loop(config, project_directory=str(tmp_path))
        assert result.status == 'failed'
        assert result.reason
