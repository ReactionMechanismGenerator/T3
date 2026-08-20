"""Tests for t3.pdep.pes_loop."""

import dataclasses
import os

import pytest

from arc.molecule.molecule import Molecule

from t3.pdep.capture import captured_qm_artifact_path
from t3.pdep.explorer.result import (EXPLORATION_STATUS_FAILED, EXPLORATION_STATUS_SUCCEEDED,
                                     PDepExplorationResult)
from t3.pdep.join import arc_ts_label
from t3.pdep.parser import PDepNetwork, PDepPathReaction
from t3.pdep.pes_loop import (PES_LOOP_CONVERGED, PES_LOOP_DIAGRAM_ONLY, PES_LOOP_FAILED,
                              PES_LOOP_MAX_ROUNDS, PES_LOOP_NO_CANDIDATES, PES_LOOP_STALLED,
                              PESLoopResult, RoundRecord, hybrid_network_path, run_pes_loop)
from t3.pdep.pes_rounds import round_paths
from t3.schema import PESLoopConfig


@pytest.fixture(autouse=True)
def stub_me_sensitivity(monkeypatch):
    """Stand in for the round's Arkane ME sensitivity analysis (t3.pdep.pes_sa), which this
    module's loop-sequencing tests cannot run (Arkane stays out of scope here, exactly as the
    explorer does). Returns finite evidence for every TS label the stub networks below can
    declare, with strictly decreasing |coefficient| in file order so the default scope='sensitive'
    ranking preserves the file order these tests were written against. Tests about the SA step
    itself re-monkeypatch over this."""
    def _fake_sa(*, network_path, sa_dir, method, timeout=None, logger=None):
        return {f'TS{i}': (0.05 - 0.001 * i, (0.05 - 0.001 * i) * 8368.0) for i in range(10)}

    monkeypatch.setattr('t3.pdep.pes_loop.run_round_me_sensitivity', _fake_sa)


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


def _alkane_adjlist(n_carbons: int) -> str:
    """A real RMG adjacency list for the linear alkane with ``n_carbons`` carbons, rendered by
    ARC's own Molecule so the stub networks below carry genuinely parseable, mutually distinct
    structures."""
    return Molecule().from_smiles('C' * n_carbons).to_adjacency_list()


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
        PDepPathReaction(label=f'reaction{i}', reactants=(f'A{i}',), products=(f'B{i}',),
                         transition_state=f'TS{i}', kinetics_type='Arrhenius',
                         kinetics_comment=f'family: {family}')
        for i, family in enumerate(families))
    # Every channel gets its own pair of REAL, distinct molecules (linear alkanes of increasing
    # length): the loop carries QM across rounds by STRUCTURAL channel identity
    # (t3.pdep.pes_rounds.channel_keys_by_ts_label), so a stub whose species had no parseable
    # structures -- or whose every reaction connected the same pair -- would silently disable the
    # cross-round memory these tests exist to exercise.
    species_structures = {}
    for i in range(len(families)):
        species_structures[f'A{i}'] = _alkane_adjlist(2 * i + 1)
        species_structures[f'B{i}'] = _alkane_adjlist(2 * i + 2)
    base_network = PDepNetwork(network_id='network1_1', path='/abs/network1_1.py', label=None,
                               species_labels=tuple(species_structures),
                               transition_state_labels=tuple(f'TS{i}' for i in
                                                              range(len(families))),
                               path_reactions=path_reactions, isomers=(), reactant_channels=(),
                               product_channels=(), species_structures=species_structures)

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

        def _runner(candidates, paths, cfg, network_id, adopted=None):
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
                              qm_runner=lambda candidates, paths, cfg, network_id, adopted=None:
                              (frozenset(), frozenset(c.ts_label for c in candidates)))
        assert result.status == PES_LOOP_STALLED
        assert len(result.rounds) == 1

    def test_converges_when_every_candidate_is_computed(self, tmp_path, monkeypatch, config):
        """Once QM exists for every barriered channel, the loop stops on its own."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])

        def _runner(candidates, paths, cfg, network_id, adopted=None):
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

        def _runner(candidates, paths, cfg, network_id, adopted=None):
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

        def _runner(candidates, paths, cfg, network_id, adopted=None):
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

        def _runner(candidates, paths, cfg, network_id, adopted=None):
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

        def _runner(candidates, paths, cfg, network_id, adopted=None):
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

        def _runner(candidates, paths, cfg, network_id, adopted=None):
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

        def _runner(candidates, paths, cfg, network_id, adopted=None):
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

        def _fake_adopt(from_t3_projects, network_id, level_of_theory, path_reaction_labels_by_ts_label):
            adopt_calls.append((tuple(from_t3_projects), network_id, level_of_theory,
                                path_reaction_labels_by_ts_label))
            return {'TS0': '/fake/prior/TS0.py'}

        monkeypatch.setattr('t3.pdep.pes_loop.adopt_prior_qm', _fake_adopt)
        # Distinct sp_level from the default opt_level (mutation (e)): a fix that accidentally
        # compared against opt_level rather than sp_level must fail this test.
        reuse_config = PESLoopConfig(pes={'network': '/abs/network1_1.py', 'source': ['HOCHO'],
                                          'bath_gas': {'He': 1.0}},
                                     qm={'sp_level': 'ccsd(t)-f12/cc-pvtz-f12'},
                                     termination={'max_rounds': 3},
                                     reuse={'from_t3_projects': ['/prior/project']})
        queued = []

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            queued.extend(c.ts_label for c in candidates)
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        run_pes_loop(reuse_config, project_directory=str(tmp_path), qm_runner=_runner)
        assert len(adopt_calls) == 1
        call_projects, call_network_id, call_level, call_labels = adopt_calls[0]
        assert call_projects == ('/prior/project',)
        assert call_network_id == 'network1_1'
        assert call_level == reuse_config.qm.sp_level == 'ccsd(t)-f12/cc-pvtz-f12'
        # The 4th argument is the seed network's STRUCTURAL channel-key map (canonical species
        # structures, direction-insensitive), never a positional-label map -- reaction and TS
        # labels are both positional in Arkane-written files and cannot identify a channel.
        assert call_labels == {'TS0': (('C',), ('CC',)), 'TS1': (('CCC',), ('CCCC',))}
        assert 'TS0' not in queued
        assert 'TS1' in queued

    def test_empty_reuse_config_does_not_call_adopt_prior_qm(self, tmp_path, monkeypatch, config):
        """The default (empty) ``reuse.from_t3_projects`` must not pay for a filesystem walk that
        can never find anything to adopt."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])
        adopt_calls = []
        monkeypatch.setattr('t3.pdep.pes_loop.adopt_prior_qm',
                            lambda *a, **kw: adopt_calls.append((a, kw)) or {})

        def _runner(candidates, paths, cfg, network_id, adopted=None):
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

        def _runner(candidates, paths, cfg, network_id, adopted=None):
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

        def _runner(candidates, paths, cfg, network_id, adopted=None):
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

    def test_round_0_full_adoption_does_not_early_return_no_candidates(self, tmp_path, monkeypatch,
                                                                        config):
        """C2 (fix round 2): when EVERY round-0 candidate was satisfied by reuse,
        split.candidates is empty -- but that must not be read as "nothing to do" and short-circuit
        to PES_LOOP_NO_CANDIDATES carrying the raw, all-ILT explored network. The headline reuse
        scenario (reuse everything -> get nothing back) was silently reported as a benign status
        before this fix; qm_runner must still be called (with an empty candidates tuple) so it can
        fold the adopted artifact into a hybrid network."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])

        def _fake_adopt(from_t3_projects, network_id, level_of_theory, path_reaction_labels_by_ts_label):
            return {'TS0': '/fake/prior/TS0.py'}

        monkeypatch.setattr('t3.pdep.pes_loop.adopt_prior_qm', _fake_adopt)
        reuse_config = PESLoopConfig(pes={'network': '/abs/network1_1.py', 'source': ['HOCHO'],
                                          'bath_gas': {'He': 1.0}},
                                     termination={'max_rounds': 3},
                                     reuse={'from_t3_projects': ['/prior/project']})
        runner_calls = []

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            runner_calls.append((tuple(candidates), adopted))
            _touch_hybrid_file(paths, network_id)
            return frozenset(), frozenset()

        result = run_pes_loop(reuse_config, project_directory=str(tmp_path), qm_runner=_runner)
        assert len(runner_calls) == 1
        called_candidates, called_adopted = runner_calls[0]
        assert called_candidates == ()
        assert called_adopted == {'TS0': '/fake/prior/TS0.py'}
        assert result.status != PES_LOOP_NO_CANDIDATES

    def test_round_0_full_adoption_advances_to_the_hybrid_file_it_wrote(self, tmp_path, monkeypatch,
                                                                        config):
        """Important (fix round 2): a round-0-only-adoption round's hybrid file, once written, must
        actually be what the NEXT round explores -- qm_runner's returned newly_computed excludes
        adopted labels by contract, so advancing current_network_path on newly_computed alone would
        ignore the hybrid this round demonstrably wrote."""
        calls, _drawn = _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])

        def _fake_adopt(from_t3_projects, network_id, level_of_theory, path_reaction_labels_by_ts_label):
            return {'TS0': '/fake/prior/TS0.py'}

        monkeypatch.setattr('t3.pdep.pes_loop.adopt_prior_qm', _fake_adopt)
        reuse_config = PESLoopConfig(pes={'network': '/abs/network1_1.py', 'source': ['HOCHO'],
                                          'bath_gas': {'He': 1.0}},
                                     termination={'max_rounds': 2, 'stop_when_no_new_ts': False},
                                     reuse={'from_t3_projects': ['/prior/project']})
        written_hybrid_paths = []

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            written_hybrid_paths.append(_touch_hybrid_file(paths, network_id))
            return frozenset(), frozenset()

        result = run_pes_loop(reuse_config, project_directory=str(tmp_path), qm_runner=_runner)
        assert len(written_hybrid_paths) == 1
        # TS0 was folded into computed_ts_labels by adoption before round 0 even ran, so round 1
        # (having nothing left to queue) correctly converges rather than exhausting max_rounds --
        # that convergence, not a particular status, is the honest outcome here. What this test
        # pins is that round 1 reached that conclusion by exploring the hybrid round 0 actually
        # wrote, not the raw pre-adoption network (the fix-round-2 regression).
        assert result.status == PES_LOOP_CONVERGED
        assert len(result.rounds) == 2
        assert result.rounds[0].status != PES_LOOP_FAILED
        assert calls[1] == written_hybrid_paths[0]


class TestRoundMeSensitivityWiring(object):
    """Defect 1: the loop's own ME SA is the source of the real sensitivity evidence every queued
    transition state must carry. These tests re-monkeypatch over this module's autouse SA stub."""

    def test_candidates_reach_the_runner_stamped_with_the_sas_own_evidence(self, tmp_path,
                                                                           monkeypatch, config):
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])
        sa_calls = []

        def _fake_sa(*, network_path, sa_dir, method, timeout=None, logger=None):
            sa_calls.append({'network_path': network_path, 'sa_dir': sa_dir, 'method': method,
                             'timeout': timeout})
            return {'TS0': (-3.0e-5, 0.25)}

        monkeypatch.setattr('t3.pdep.pes_loop.run_round_me_sensitivity', _fake_sa)
        received = []

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            received.extend(candidates)
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert len(sa_calls) == 1
        # The SA runs on THIS round's freshly explored network, into this round's own SA dir,
        # with the configured ME method and the explorer's own timeout budget.
        assert os.path.basename(sa_calls[0]['network_path']) == 'explored_round0.py'
        assert sa_calls[0]['sa_dir'] == round_paths(str(tmp_path), 0).sa
        assert sa_calls[0]['method'] == config.pes.method
        assert sa_calls[0]['timeout'] == config.pes.timeout
        assert received[0].coefficient == -3.0e-5
        assert received[0].delta_ln_k == 0.25

    def test_a_failed_sa_fails_the_round_rather_than_inventing_evidence(self, tmp_path,
                                                                        monkeypatch, config):
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])

        def _fake_sa(*, network_path, sa_dir, method, timeout=None, logger=None):
            raise ValueError('the master-equation sensitivity analysis failed')

        monkeypatch.setattr('t3.pdep.pes_loop.run_round_me_sensitivity', _fake_sa)
        runner_calls = []

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            runner_calls.append(candidates)
            return frozenset(), frozenset()

        result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert result.status == PES_LOOP_FAILED
        assert 'sensitivity' in result.reason
        assert runner_calls == []

    def test_a_candidate_without_evidence_is_skipped_not_queued(self, tmp_path, monkeypatch, config):
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO', '1,2_Insertion_CO'])

        def _fake_sa(*, network_path, sa_dir, method, timeout=None, logger=None):
            return {'TS1': (-3.0e-5, 0.25)}  # nothing for TS0

        monkeypatch.setattr('t3.pdep.pes_loop.run_round_me_sensitivity', _fake_sa)
        queued = []

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            queued.extend(c.ts_label for c in candidates)
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert 'TS0' not in queued
        assert 'TS1' in queued
        skipped_labels = {s.label: s.reason for record in result.rounds for s in record.skipped}
        assert 'reaction0' in skipped_labels
        assert 'no finite sensitivity evidence' in skipped_labels['reaction0']

    def test_no_candidate_with_evidence_fails_the_round(self, tmp_path, monkeypatch, config):
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])
        monkeypatch.setattr('t3.pdep.pes_loop.run_round_me_sensitivity',
                            lambda **kwargs: dict())
        result = run_pes_loop(config, project_directory=str(tmp_path),
                              qm_runner=lambda candidates, paths, cfg, network_id, adopted=None:
                              (frozenset(), frozenset()))
        assert result.status == PES_LOOP_FAILED
        assert 'no finite sensitivity evidence' in result.reason

    def test_scope_sensitive_ranks_by_the_sas_own_evidence_not_file_order(self, tmp_path,
                                                                          monkeypatch):
        """Before this fix, scope='sensitive' degraded to file order with a warning because the
        standalone loop had no SA to rank by. It has one now: with four candidates and a budget
        of two, the runner must receive the two largest-|coefficient| ones, ordered by descending
        |coefficient| -- a full ranking, not a two-candidate argmax, and nothing like file order
        (which would have offered TS0 then TS1)."""
        config = PESLoopConfig(pes={'network': '/abs/network1_1.py', 'source': ['HOCHO'],
                                    'bath_gas': {'He': 1.0}},
                               qm={'scope': 'sensitive', 'max_transition_states_per_round': 2},
                               termination={'max_rounds': 1, 'stop_when_no_new_ts': False})
        _stub_explorer(monkeypatch, tmp_path,
                       families=['1,2_Insertion_CO'] * 4)

        def _fake_sa(*, network_path, sa_dir, method, timeout=None, logger=None):
            # Signs deliberately mixed: the ranking is by |coefficient|, not by signed value.
            return {'TS0': (1.0e-6, 0.008), 'TS1': (-9.0e-5, 0.75),
                    'TS2': (5.0e-5, 0.42), 'TS3': (-2.0e-6, 0.017)}

        monkeypatch.setattr('t3.pdep.pes_loop.run_round_me_sensitivity', _fake_sa)
        offered = []

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            offered.append(tuple(c.ts_label for c in candidates))
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert offered == [('TS1', 'TS2')]

    def test_diagram_only_runs_no_sa_at_all(self, tmp_path, monkeypatch, config):
        """qm_runner=None spends no QM and captures nothing, so there is no evidence to justify
        and no Arkane SA to pay for."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])
        sa_calls = []

        def _fake_sa(**kwargs):
            sa_calls.append(kwargs)
            return {}

        monkeypatch.setattr('t3.pdep.pes_loop.run_round_me_sensitivity', _fake_sa)
        result = run_pes_loop(config, project_directory=str(tmp_path))
        assert result.status == PES_LOOP_DIAGRAM_ONLY
        assert sa_calls == []


class TestCumulativeAdoptedPlumbing(object):
    """Defect 3's LOOP-SIDE plumbing, pinned against runner doubles: what run_pes_loop passes as
    `adopted` each round, and that a runner cannot corrupt the loop's own record. Whether the
    artifacts actually survive into a later round's hybrid is deliberately NOT claimed here --
    that end-to-end fact is pinned against the real capture and hybrid writer in
    test_pes_loop_integration.py."""

    def test_every_round_receives_the_cumulative_artifacts_as_adopted(self, tmp_path, monkeypatch):
        """Three TSs: TS0 adopted from a prior project, TS1 converges in round 0, TS2 in round 1.
        Round 0's runner must see the prior-project artifact; round 1's runner must ALSO see round
        0's converged TS1, pointing at round 0's own capture (the originating capture, whose
        verified manifest still exists), alongside the prior-project TS0."""
        config = PESLoopConfig(pes={'network': '/abs/network1_1.py', 'source': ['HOCHO'],
                                    'bath_gas': {'He': 1.0}},
                               termination={'max_rounds': 3},
                               reuse={'from_t3_projects': ['/prior/project']})
        _stub_explorer(monkeypatch, tmp_path,
                       families=['1,2_Insertion_CO', '1,2_Insertion_CO', '1,2_Insertion_CO'])
        monkeypatch.setattr('t3.pdep.pes_loop.adopt_prior_qm',
                            lambda *args: {'TS0': '/prior/capture/qm/TS0.py'})
        adopted_per_round = []

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            adopted_per_round.append(dict(adopted))
            _touch_hybrid_file(paths, network_id)
            converged = frozenset({candidates[0].ts_label}) if candidates else frozenset()
            return converged, frozenset(c.ts_label for c in candidates)

        result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert result.status == PES_LOOP_CONVERGED
        assert len(adopted_per_round) == 2
        assert adopted_per_round[0] == {'TS0': '/prior/capture/qm/TS0.py'}
        round0_ts1_path = captured_qm_artifact_path(
            round_paths(str(tmp_path), 0).capture, arc_ts_label('explored_round0', 'TS1'))
        assert adopted_per_round[1] == {'TS0': '/prior/capture/qm/TS0.py',
                                        'TS1': round0_ts1_path}

    def test_a_runner_mutating_its_adopted_dict_cannot_corrupt_the_loops_own_record(
            self, tmp_path, monkeypatch, config):
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO', '1,2_Insertion_CO'])
        adopted_per_round = []

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            adopted_per_round.append(dict(adopted))
            adopted['INJECTED'] = '/nowhere.py'  # a hostile/buggy runner
            _touch_hybrid_file(paths, network_id)
            converged = frozenset({candidates[0].ts_label}) if candidates else frozenset()
            return converged, frozenset(c.ts_label for c in candidates)

        run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert len(adopted_per_round) >= 2
        assert 'INJECTED' not in adopted_per_round[1]


class TestSensitivityFloor(object):
    """qm.min_delta_ln_k mirrors T3's in-run pdep_min_delta_ln_k: a below-floor (e.g. measured
    zero) response never buys QM, and 'every remaining channel below the floor' is a measured
    terminal outcome, not a fault."""

    def test_a_below_floor_candidate_is_never_queued(self, tmp_path, monkeypatch, config):
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO', '1,2_Insertion_CO'])

        def _fake_sa(*, network_path, sa_dir, method, timeout=None, logger=None):
            return {'TS0': (0.0, 0.0), 'TS1': (-3.0e-5, 0.25)}

        monkeypatch.setattr('t3.pdep.pes_loop.run_round_me_sensitivity', _fake_sa)
        queued = []

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            queued.extend(c.ts_label for c in candidates)
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert 'TS0' not in queued
        assert 'TS1' in queued
        reasons = {s.label: s.reason for record in result.rounds for s in record.skipped}
        assert 'below the min_delta_ln_k floor' in reasons['reaction0']

    def test_all_below_floor_at_round_0_is_no_candidates_not_failed(self, tmp_path, monkeypatch,
                                                                    config):
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO'])
        monkeypatch.setattr('t3.pdep.pes_loop.run_round_me_sensitivity',
                            lambda **kwargs: {'TS0': (0.0, 0.0)})
        runner_calls = []

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            runner_calls.append(candidates)
            return frozenset(), frozenset()

        result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert result.status == PES_LOOP_NO_CANDIDATES
        assert runner_calls == []
        assert 'below the min_delta_ln_k floor' in result.rounds[0].skipped[0].reason

    def test_all_below_floor_after_progress_is_converged(self, tmp_path, monkeypatch, config):
        """Round 0 converges the above-floor channel; round 1's only remaining channel measures
        below the floor -- a measured 'done', reported as convergence with the skip recorded."""
        _stub_explorer(monkeypatch, tmp_path, families=['1,2_Insertion_CO', '1,2_Insertion_CO'])

        def _fake_sa(*, network_path, sa_dir, method, timeout=None, logger=None):
            return {'TS0': (-3.0e-5, 0.25), 'TS1': (0.0, 0.0)}

        monkeypatch.setattr('t3.pdep.pes_loop.run_round_me_sensitivity', _fake_sa)

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            _touch_hybrid_file(paths, network_id)
            return frozenset(c.ts_label for c in candidates), frozenset(c.ts_label for c in candidates)

        result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)
        assert result.status == PES_LOOP_CONVERGED
        assert [record.status for record in result.rounds] == ['continuing', PES_LOOP_CONVERGED]
        assert sorted(result.rounds[0].queued_ts_labels) == ['TS0']


class TestRenumberedNetworksKeepQMOnTheRightChannel(object):
    """THE wrong-saddle-point pin: every TS label Arkane writes is positional
    (rmgpy/rmg/pdep.py:854), so a re-exploration can relabel channels between rounds. A
    label-keyed carry would then skip the WRONG channel as 'already has QM' and fold round 0's
    barrier onto a different reaction; the structural carry must keep the barrier on its own
    channel regardless of what the labels do."""

    def test_round_1_relabeling_does_not_misattribute_round_0s_artifact(self, tmp_path,
                                                                        monkeypatch, config):
        adj = {name: _alkane_adjlist(n)
               for name, n in (('A0', 1), ('B0', 2), ('A1', 3), ('B1', 4))}

        def _reaction(label, ts, reactants, products):
            return PDepPathReaction(label=label, reactants=reactants, products=products,
                                    transition_state=ts, kinetics_type='Arrhenius',
                                    kinetics_comment='family: 1,2_Insertion_CO')

        # Round 0's file: TS0 is the (A0, B0) channel. Round 1's re-exploration writes the SAME
        # two channels with the positions -- and therefore the labels -- swapped: 'TS0' now names
        # the (A1, B1) channel and 'TS1' names (A0, B0).
        round0_network = PDepNetwork(
            network_id='explored_round0', path='/x0.py', label=None,
            species_labels=tuple(adj), transition_state_labels=('TS0', 'TS1'),
            path_reactions=(_reaction('reaction0', 'TS0', ('A0',), ('B0',)),
                            _reaction('reaction1', 'TS1', ('A1',), ('B1',))),
            isomers=(), reactant_channels=(), product_channels=(), species_structures=adj)
        round1_network = dataclasses.replace(
            round0_network, network_id='explored_round1',
            path_reactions=(_reaction('reaction0', 'TS0', ('A1',), ('B1',)),
                            _reaction('reaction1', 'TS1', ('A0',), ('B0',))))
        networks = [round0_network, round1_network, round1_network]

        calls = []

        def _fake_explore(*, network_path, config, logger=None):
            calls.append(network_path)
            output_path = os.path.join(str(tmp_path), f'explored_round{len(calls) - 1}.py')
            open(output_path, 'w').close()
            return PDepExplorationResult(network_id='network1_1',
                                         status=EXPLORATION_STATUS_SUCCEEDED,
                                         network_paths=(output_path,))

        def _fake_parse(path):
            # The loop parses each round's freshly explored network right after exploring it, so
            # the round currently being parsed is simply the latest explore call.
            network = networks[min(len(calls) - 1, 2)] if calls else networks[0]
            return dataclasses.replace(network, path=path)

        def _fake_draw(network_path, output_path):
            open(output_path, 'w').close()

        monkeypatch.setattr('t3.pdep.pes_loop.explore_pdep_network', _fake_explore)
        monkeypatch.setattr('t3.pdep.pes_loop.parse_pdep_network_file', _fake_parse)
        monkeypatch.setattr('t3.pdep.pes_loop.draw_pes_diagram', _fake_draw)

        offered_channels_per_round = []
        adopted_per_round = []

        def _runner(candidates, paths, cfg, network_id, adopted=None):
            offered_channels_per_round.append(
                tuple((c.ts_label, c.path_reaction.reactants) for c in candidates))
            adopted_per_round.append(dict(adopted))
            _touch_hybrid_file(paths, network_id)
            # Round 0 converges only its first offered channel; round 1 converges the rest.
            converged = frozenset({candidates[0].ts_label}) if candidates else frozenset()
            return converged, frozenset(c.ts_label for c in candidates)

        result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=_runner)

        assert result.status == PES_LOOP_CONVERGED
        # Round 0: both channels offered; (A0, B0) converged under its round-0 label 'TS0'.
        assert offered_channels_per_round[0] == (('TS0', ('A0',)), ('TS1', ('A1',)))
        # Round 1 (relabeled): the computed (A0, B0) channel now carries label 'TS1' and must be
        # SKIPPED under that new label; the still-uncomputed (A1, B1) channel -- now labeled
        # 'TS0', which a label-keyed carry would have wrongly skipped as computed -- must be the
        # one offered.
        assert offered_channels_per_round[1] == (('TS0', ('A1',)),)
        # And round 0's artifact reaches round 1 as an adopted artifact under the channel's NEW
        # label 'TS1', still pointing at the ORIGINATING (round 0) capture, where it was vendored
        # under round 0's arc label for 'TS0' -- the right barrier, on the right channel.
        round0_artifact = captured_qm_artifact_path(
            round_paths(str(tmp_path), 0).capture, arc_ts_label('explored_round0', 'TS0'))
        assert adopted_per_round[1] == {'TS1': round0_artifact}
