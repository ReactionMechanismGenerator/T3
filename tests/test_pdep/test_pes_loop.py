"""Tests for t3.pdep.pes_loop."""

import dataclasses
import os

import pytest

from t3.pdep.explorer.result import (EXPLORATION_STATUS_FAILED, EXPLORATION_STATUS_SUCCEEDED,
                                     PDepExplorationResult)
from t3.pdep.parser import PDepNetwork, PDepPathReaction
from t3.pdep.pes_loop import (PES_LOOP_CONVERGED, PES_LOOP_DIAGRAM_ONLY, PES_LOOP_FAILED,
                              PES_LOOP_MAX_ROUNDS, PES_LOOP_NO_CANDIDATES, PES_LOOP_STALLED,
                              PESLoopResult, RoundRecord, hybrid_network_path, run_pes_loop)
from t3.pdep.pes_rounds import round_paths
from t3.schema import PESLoopConfig


@pytest.fixture
def config():
    """A minimal config pointing at a network fixture.

    ``bath_gas`` is required by ``PDepExplorerConfig`` (a config without one is a guaranteed
    ``__post_init__`` failure -- see t3/pdep/explorer/config.py, issue #183), so it must be
    non-empty here. Every test below stubs ``explore_pdep_network`` itself, but ``run_pes_loop``
    still builds a real ``PDepExplorerConfig`` from this fixture before handing it to the stub, so
    that construction IS reached and validated -- only the actual exploration is faked out.
    """
    return PESLoopConfig(pes={'network': '/abs/network1_1.py', 'source': ['HOCHO'],
                              'bath_gas': {'He': 1.0}},
                         termination={'max_rounds': 3})


def _stub_explorer(monkeypatch, tmp_path, families, fail=False):
    """
    Patch ``explore_pdep_network``, ``parse_pdep_network_file``, and ``draw_pes_diagram`` at the
    names bound in ``t3.pdep.pes_loop``, so ``run_pes_loop`` can be tested without Arkane.

    Unlike a stub that discards its arguments, this one is faithful on two points that matter for
    catching a real "loop that never actually loops" bug:

    - The exploration stub records the ``network_path`` it was called with (in ``calls``, in call
      order), so a test can pin that round 1 explores round 0's hybrid output rather than
      re-exploring round 0's original input.
    - Each call returns a DISTINCT, real (on-disk, in a scratch temp dir) output path named
      ``explored_round{N}.py``, and the parse stub derives ``network_id`` from that path's stem
      exactly as the real ``parse_pdep_network_file`` does (``Path(path).stem`` --
      ``t3/pdep/parser.py:729``). The file must actually exist because ``run_pes_loop`` checks for
      it at the top of every round after the first (Important 4). This means a loop that ever passes
      the wrong network id to ``qm_runner`` or draws from the wrong path is directly observable,
      not masked by every round looking identical.

    Args:
        monkeypatch: pytest's monkeypatch fixture.
        tmp_path: pytest's per-test scratch directory fixture, used instead of tempfile.mkdtemp
            so this helper leaves nothing behind for pytest to clean up itself.
        families (list): One path reaction is fabricated per entry, carrying that family in its
            kinetics comment (recognized by ``t3.pdep.barrierless.classify_barrierless``).
        fail (bool): If True, every call to the stubbed explorer reports a failed exploration
            instead of a succeeded one.

    Returns:
        tuple: ``(calls, drawn)`` -- the ``network_path`` each ``explore_pdep_network`` call
        received, and the ``network_path`` each ``draw_pes_diagram`` call received, both in call
        order.
    """
    path_reactions = tuple(
        PDepPathReaction(label=f'reaction{i}', reactants=('A',), products=('B',),
                         transition_state=f'TS{i}', kinetics_type='Arrhenius',
                         kinetics_comment=f'family: {family}')
        for i, family in enumerate(families))
    base_network = PDepNetwork(network_id='network1_1', path='/abs/network1_1.py', label=None,
                               species_labels=('A', 'B'),
                               transition_state_labels=tuple(f'TS{i}' for i in
                                                              range(len(families))),
                               path_reactions=path_reactions, isomers=(), reactant_channels=(),
                               product_channels=())

    calls = []
    drawn = []
    explore_dir = os.path.join(str(tmp_path), 'explored')
    os.makedirs(explore_dir, exist_ok=True)

    def _fake_explore(*, network_path, config, logger=None):
        calls.append(network_path)
        if fail:
            return PDepExplorationResult(network_id=base_network.network_id,
                                         status=EXPLORATION_STATUS_FAILED,
                                         reasons=('the explorer adapter crashed',))
        # A real explore_pdep_network writes a real file to disk; the file-existence guard at the
        # top of every round after the first (t3.pdep.pes_loop.run_pes_loop) depends on that, so
        # this stub must too -- a bare fabricated path string is not enough once a round is
        # allowed to re-explore its own prior output (Important 4) rather than always advancing to
        # a hybrid file a separate helper (_touch_hybrid_file) creates.
        output_path = os.path.join(explore_dir, f'explored_round{len(calls) - 1}.py')
        open(output_path, 'w').close()
        return PDepExplorationResult(network_id=base_network.network_id,
                                     status=EXPLORATION_STATUS_SUCCEEDED,
                                     network_paths=(output_path,))

    def _fake_parse(path):
        network_id = os.path.splitext(os.path.basename(path))[0]
        return dataclasses.replace(base_network, network_id=network_id, path=path)

    def _fake_draw(network_path, output_path):
        drawn.append(network_path)
        open(output_path, 'w').close()

    monkeypatch.setattr('t3.pdep.pes_loop.explore_pdep_network', _fake_explore)
    monkeypatch.setattr('t3.pdep.pes_loop.parse_pdep_network_file', _fake_parse)
    monkeypatch.setattr('t3.pdep.pes_loop.draw_pes_diagram', _fake_draw)

    return calls, drawn


def _touch_hybrid_file(paths, network_id):
    """
    Create the round's hybrid network file, standing in for what a real ``qm_runner`` (Task 6)
    must write there. ``run_pes_loop`` checks for exactly this file at the top of every round after
    the first, so any test whose loop is meant to continue past round 0 must call this from its
    ``qm_runner`` stub.
    """
    path = hybrid_network_path(paths, network_id)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, 'w').close()
    return path


class TestRunPESLoop(object):

    def test_no_candidates_terminates_immediately(self, tmp_path, monkeypatch, config):
        """A network whose every channel is barrierless has nothing to compute: one explore, one
        diagram, done -- not three empty rounds."""
        _stub_explorer(monkeypatch, tmp_path, families=['R_Recombination', 'R_Recombination'])
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
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])

        def _runner(candidates, paths, cfg, network_id):
            _touch_hybrid_file(paths, network_id)
            return frozenset(), frozenset(c.ts_label for c in candidates)

        result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert result.status == PES_LOOP_MAX_ROUNDS
        assert len(result.rounds) == 3

    def test_a_round_that_stalls_stops_early(self, tmp_path, monkeypatch, config):
        """With the default stop_when_no_new_ts=True, a runner that converges nothing must stop
        after round 1 rather than spend the rest of the round budget -- this is what distinguishes
        'stalled' from 'max_rounds' (ruling C3)."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])
        result = run_pes_loop(config, project_directory=str(tmp_path),
                              qm_runner=lambda candidates, paths, cfg, network_id:
                              (frozenset(), frozenset(c.ts_label for c in candidates)))
        assert result.status == PES_LOOP_STALLED
        assert len(result.rounds) == 1

    def test_converges_when_every_candidate_is_computed(self, tmp_path, monkeypatch, config):
        """Once QM exists for every barriered channel, the loop stops on its own."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])

        def _runner(candidates, paths, cfg, network_id):
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert result.status == PES_LOOP_CONVERGED
        assert len(result.rounds) == 2

    def test_a_diagram_is_drawn_every_round(self, tmp_path, monkeypatch, config):
        """The whole point: a picture per round, each from that round's own FRESHLY EXPLORED
        network -- not the pre-explore input, and not another round's output. Drawing from the
        input instead of the explored network is exactly the bug this feature exists to remove
        (byte-identical diagrams across runs even when a saddle point converged)."""
        calls, drawn = _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])

        def _runner(candidates, paths, cfg, network_id):
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert len(result.rounds) == 2
        for record in result.rounds:
            assert record.diagram_path is not None
            assert os.path.basename(record.diagram_path) == 'pes_diagram.png'
        assert [os.path.basename(path) for path in drawn] == \
            ['explored_round0.py', 'explored_round1.py']
        assert drawn != calls

    def test_barrierless_channels_are_never_queued(self, tmp_path, monkeypatch, config):
        """The F21 gate, asserted end to end rather than only at the unit level."""
        queued = []
        _stub_explorer(monkeypatch, tmp_path, families=['R_Recombination', '1,2_Insertion_CO'])

        def _runner(candidates, paths, cfg, network_id):
            queued.extend(c.family for c in candidates)
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert 'R_Recombination' not in queued
        assert '1,2_Insertion_CO' in queued

    def test_no_qm_runner_explores_and_draws_only(self, tmp_path, monkeypatch, config):
        """The honest behaviour without a configured runner: one round, a diagram, no crash -- and
        a status distinct from an operational stall, per the ruling on concern 4."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])
        result = run_pes_loop(config, project_directory=str(tmp_path))
        assert result.status == PES_LOOP_DIAGRAM_ONLY
        assert len(result.rounds) == 1
        assert result.rounds[0].diagram_path is not None

    def test_diagram_only_branch_reports_every_offered_candidate_as_queued(self, tmp_path,
                                                                           monkeypatch, config):
        """N3 (round 3): qm_runner=None never runs a runner at all, so nothing CAN be dropped --
        the diagram-only branch's queued_ts_labels must report every candidate that was offered.
        Mutation M6 (reverting this branch to queued_ts_labels=()) must fail this test."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO', '1,2_Insertion_CO'])
        result = run_pes_loop(config, project_directory=str(tmp_path))
        assert result.rounds[0].queued_ts_labels == ('TS0', 'TS1')

    def test_queued_ts_labels_reflects_what_the_runner_actually_queued_not_what_was_offered(
            self, tmp_path, monkeypatch, config):
        """N3 (round 3): every qm_runner double elsewhere in this file returns the offered set
        verbatim as its queued half, which cannot distinguish RoundRecord.queued_ts_labels from
        the offered set the loop built before the runner ran. This runner drops one offered
        candidate from what it reports as queued (as build_arc_input legitimately can, for a
        candidate missing structure) -- proving the loop records the runner's OWN report, not what
        it offered. Mutation M5 (queued_ts_labels = tuple(offered_ts_labels) at pes_loop.py:348)
        must fail this test."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO', '1,2_Insertion_CO'])

        def _runner(candidates, paths, cfg, network_id):
            offered = tuple(c.ts_label for c in candidates)
            assert offered == ('TS0', 'TS1')
            _touch_hybrid_file(paths, network_id)
            # Converge nothing, and report only TS0 as queued -- TS1 is dropped, as build_arc_input
            # legitimately drops a candidate missing structure.
            return frozenset(), frozenset({'TS0'})

        result = run_pes_loop(config, project_directory=str(tmp_path),
                              qm_runner=_runner)
        assert result.rounds[0].queued_ts_labels == ('TS0',)
        assert 'TS1' not in result.rounds[0].queued_ts_labels

    def test_a_failed_explore_stops_the_loop_and_says_why(self, tmp_path, monkeypatch, config):
        """A round that cannot explore must not be papered over as convergence."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'], fail=True)
        result = run_pes_loop(config, project_directory=str(tmp_path))
        assert result.status == 'failed'
        assert result.reason

    def test_round_1_explores_the_hybrid_path_not_the_original_network(self, tmp_path, monkeypatch,
                                                                        config):
        """Round 1 must explore the QM-informed hybrid network qm_runner wrote, not re-explore
        round 0's original input -- otherwise this is N re-explorations of one file, never a loop
        (Critical 1)."""
        calls, _ = _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])

        def _runner(candidates, paths, cfg, network_id):
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert calls[0] == '/abs/network1_1.py'
        expected_round1_path = hybrid_network_path(round_paths(str(tmp_path), 0),
                                                    'explored_round0')
        assert calls[1] == expected_round1_path
        assert calls[1] != calls[0]

    def test_qm_runner_receives_the_explored_networks_id(self, tmp_path, monkeypatch, config):
        """qm_runner must be handed the id of the network THIS round explored, not the loop's
        static input path's id -- ruling C4 exists precisely so a wrong id cannot misnamespace ARC
        jobs (Critical 1)."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])
        received = []

        def _runner(candidates, paths, cfg, network_id):
            received.append(network_id)
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert received[0] == 'explored_round0'
        assert received[0] != 'network1_1'

    def test_adopted_ts_labels_seed_removes_that_ts_from_round_0_candidates(self, tmp_path,
                                                                             monkeypatch, config):
        """C1: TS labels adopted from an earlier T3 project must be treated as already computed
        before round 0 ever runs, not just before some later round -- there was previously zero
        coverage of this ruling at all."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO', '1,2_Insertion_CO'])
        queued = []

        def _runner(candidates, paths, cfg, network_id):
            queued.extend(c.ts_label for c in candidates)
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner,
                    adopted_ts_labels=frozenset({'TS0'}))
        assert 'TS0' not in queued
        assert 'TS1' in queued

    def test_reuse_config_calls_adopt_prior_qm_and_seeds_round_0(self, tmp_path, monkeypatch,
                                                                  config):
        """Ruling 1: ``config.reuse.from_t3_projects`` being non-empty must actually call
        ``adopt_prior_qm`` and remove what it returns from round 0's candidates -- a config knob
        that has no reachable call site is exactly the "ships dead" failure this test exists to
        catch."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO', '1,2_Insertion_CO'])
        adopt_calls = []

        def _fake_adopt(from_t3_projects, network_id, level_of_theory):
            adopt_calls.append((tuple(from_t3_projects), network_id, level_of_theory))
            return {'TS0': '/fake/prior/TS0.py'}

        monkeypatch.setattr('t3.pdep.pes_loop.adopt_prior_qm', _fake_adopt)
        reuse_config = PESLoopConfig(pes={'network': '/abs/network1_1.py', 'source': ['HOCHO'],
                                          'bath_gas': {'He': 1.0}},
                                     termination={'max_rounds': 3},
                                     reuse={'from_t3_projects': ['/prior/project']})
        queued = []

        def _runner(candidates, paths, cfg, network_id):
            queued.extend(c.ts_label for c in candidates)
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        run_pes_loop(reuse_config, project_directory=str(tmp_path), qm_runner=_runner)
        assert adopt_calls == [(('/prior/project',), 'network1_1', reuse_config.qm.sp_level)]
        assert 'TS0' not in queued
        assert 'TS1' in queued

    def test_empty_reuse_config_does_not_call_adopt_prior_qm(self, tmp_path, monkeypatch, config):
        """The default (empty) ``reuse.from_t3_projects`` must not pay for a filesystem walk that
        can never find anything to adopt."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])
        adopt_calls = []
        monkeypatch.setattr('t3.pdep.pes_loop.adopt_prior_qm',
                            lambda *a, **kw: adopt_calls.append((a, kw)) or {})

        def _runner(candidates, paths, cfg, network_id):
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert adopt_calls == []

    def test_a_missing_hybrid_file_fails_loudly_instead_of_a_bare_traceback(self, tmp_path,
                                                                             monkeypatch, config):
        """qm_runner's contract is to write the round's hybrid network file; if it doesn't, that
        must surface as a PES_LOOP_FAILED round with a stated reason, not a FileNotFoundError
        escaping from round N after round N-1 already spent real ARC time (Important 3)."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])

        def _runner(candidates, paths, cfg, network_id):
            # Deliberately does NOT write the hybrid file, breaking its contract.
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert result.status == PES_LOOP_FAILED
        assert len(result.rounds) == 2
        assert result.rounds[-1].status == PES_LOOP_FAILED
        assert 'hybrid' in result.reason.lower()

    def test_a_round_that_converges_nothing_does_not_advance_to_a_hybrid_file_never_written(
            self, tmp_path, monkeypatch, config):
        """A round whose qm_runner converges nothing correctly writes no hybrid file (its
        empty-convergence contract -- see t3.pdep.pes_qm.arc_qm_runner). The loop must re-explore
        that round's own explored network next round rather than advance current_network_path to
        a hybrid file that was never written -- otherwise the next round's file-existence guard
        blames qm_runner for breaking a contract it never had, turning an honest stall into a
        bogus PES_LOOP_FAILED (Important 4)."""
        config = PESLoopConfig(pes={'network': '/abs/network1_1.py', 'source': ['HOCHO'],
                                    'bath_gas': {'He': 1.0}},
                               termination={'max_rounds': 2, 'stop_when_no_new_ts': False})

        path_reactions = (PDepPathReaction(label='reaction0', reactants=('A',), products=('B',),
                                           transition_state='TS0', kinetics_type='Arrhenius',
                                           kinetics_comment='family: 1,2_Insertion_CO'),)
        base_network = PDepNetwork(network_id='network1_1', path='/abs/network1_1.py', label=None,
                                   species_labels=('A', 'B'), transition_state_labels=('TS0',),
                                   path_reactions=path_reactions, isomers=(), reactant_channels=(),
                                   product_channels=())

        calls = []
        explored_paths = []

        def _fake_explore(*, network_path, config, logger=None):
            calls.append(network_path)
            output_path = os.path.join(str(tmp_path), f'explored_round{len(calls) - 1}.py')
            open(output_path, 'w').close()
            explored_paths.append(output_path)
            return PDepExplorationResult(network_id=base_network.network_id,
                                         status=EXPLORATION_STATUS_SUCCEEDED,
                                         network_paths=(output_path,))

        def _fake_parse(path):
            network_id = os.path.splitext(os.path.basename(path))[0]
            return dataclasses.replace(base_network, network_id=network_id, path=path)

        def _fake_draw(network_path, output_path):
            open(output_path, 'w').close()

        monkeypatch.setattr('t3.pdep.pes_loop.explore_pdep_network', _fake_explore)
        monkeypatch.setattr('t3.pdep.pes_loop.parse_pdep_network_file', _fake_parse)
        monkeypatch.setattr('t3.pdep.pes_loop.draw_pes_diagram', _fake_draw)

        def _runner(candidates, paths, cfg, network_id):
            # Deliberately does NOT write the hybrid file -- this is the correct, contractual
            # behaviour when nothing converged, not a bug being simulated.
            return frozenset(), frozenset(c.ts_label for c in candidates)

        result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert result.status == PES_LOOP_MAX_ROUNDS
        assert len(result.rounds) == 2
        assert all(round_record.status != PES_LOOP_FAILED for round_record in result.rounds)
        # Round 1 must re-explore round 0's own explored output, not a hybrid file round 0's
        # qm_runner never wrote.
        assert calls[1] == explored_paths[0]
        never_written_hybrid = hybrid_network_path(round_paths(str(tmp_path), 0),
                                                    'explored_round0')
        assert calls[1] != never_written_hybrid
        assert not os.path.isfile(never_written_hybrid)
