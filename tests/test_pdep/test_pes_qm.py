"""Tests for t3.pdep.pes_qm."""

import os
import shutil

from arc.reaction.reaction import ARCReaction
from arc.species.species import ARCSpecies

import t3.pdep.pes_qm as pes_qm
from t3.pdep.capture import CaptureResult
from t3.pdep.discovery import ARTIFACT_STATUS_UNVERIFIED, ARTIFACT_STATUS_USABLE, TSArtifactRecord
from t3.pdep.hybrid import HybridNetworkResult
from t3.pdep.join import JOIN_STATUS_QUEUED, arc_ts_label
from t3.pdep.parser import PDepPathReaction
from t3.pdep.pes_loop import hybrid_network_path
from t3.pdep.pes_qm import ARC_INPUT_FILE_NAME, arc_qm_runner, build_arc_input
from t3.pdep.pes_rounds import QMCandidate, round_paths
from t3.schema import PESLoopConfig


def _config(**qm) -> PESLoopConfig:
    return PESLoopConfig(pes={'network': '/abs/n.py', 'source': ['A']}, **({'qm': qm} if qm else {}))


def _candidate(ts_label='TS1', family='1,2_Insertion_CO', reactants=('A',), products=('B',),
              label='reaction1') -> QMCandidate:
    path_reaction = PDepPathReaction(
        label=label, reactants=reactants, products=products,
        transition_state=ts_label, kinetics_type='Arrhenius', kinetics_comment='',
    )
    return QMCandidate(path_reaction=path_reaction, ts_label=ts_label, family=family)


# Adjacency-list text is opaque to build_arc_input -- any non-empty string exercises the
# species_structures plumbing without needing a chemically valid structure.
_STRUCTURES = {'A': '1 A u0 p0 c0\n', 'B': '1 B u0 p0 c0\n'}

# A real, atom-balanced RMG adjacency list (ethyl radical) used for the one test (N1) that must
# construct a genuine ARCReaction/ARCSpecies rather than merely assert on dict shape -- ARC's own
# atom-balance check (check_atom_balance) requires both wells to carry the same molecular formula.
_ETHYL_ADJLIST = ('1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}\n'
                  '2 C u1 p0 c0 {1,S} {6,S} {7,S}\n'
                  '3 H u0 p0 c0 {1,S}\n'
                  '4 H u0 p0 c0 {1,S}\n'
                  '5 H u0 p0 c0 {1,S}\n'
                  '6 H u0 p0 c0 {2,S}\n'
                  '7 H u0 p0 c0 {2,S}')
_BALANCED_STRUCTURES = {'A': _ETHYL_ADJLIST, 'B': _ETHYL_ADJLIST}


class TestBuildARCInput(object):

    def test_levels_reach_arc_undashed(self):
        arc_input = build_arc_input((), round_paths('/proj', 0), _config(), 'network1_1', {})
        assert arc_input['opt_level'] == 'wb97xd/def2tzvp'
        assert '-' not in arc_input['opt_level'].split('/')[1]

    def test_rotors_and_irc_off_by_default(self):
        arc_input = build_arc_input((), round_paths('/proj', 0), _config(), 'network1_1', {})
        assert arc_input['job_types']['rotors'] is False
        assert arc_input['job_types']['irc'] is False

    def test_rotors_and_irc_can_be_turned_on(self):
        arc_input = build_arc_input((), round_paths('/proj', 0),
                                    _config(rotors=True, irc=True), 'network1_1', {})
        assert arc_input['job_types']['rotors'] is True
        assert arc_input['job_types']['irc'] is True

    def test_ts_adapters_are_passed_through(self):
        arc_input = build_arc_input((), round_paths('/proj', 0),
                                    _config(ts_adapters=['linear', 'rits']), 'network1_1', {})
        assert arc_input['ts_adapters'] == ['linear', 'rits']

    def test_project_directory_is_the_rounds_own_arc_project(self):
        paths = round_paths('/proj', 2)
        arc_input = build_arc_input((), paths, _config(), 'network1_1', {})
        assert arc_input['project_directory'] == paths.arc_project

    def test_ts_labels_are_namespaced_not_network_local(self):
        """A network-local 'TS1' collides across networks; join.arc_ts_label namespaces it so a
        network can never be joined to another network's quantum chemistry."""
        candidate = _candidate(ts_label='TS1')
        arc_input = build_arc_input((candidate,), round_paths('/proj', 0), _config(), 'network1_1',
                                    _STRUCTURES)
        assert arc_input['reactions'][0]['ts_label'] == 'T3PDep_network1_1_TS1'

    def test_sp_scan_irc_level_and_family_reach_arc_input(self):
        candidate = _candidate(family='1,2_Insertion_CO')
        arc_input = build_arc_input(
            (candidate,), round_paths('/proj', 0),
            _config(sp_level='ccsd(t)-f12/cc-pvtz-f12', scan_level='wb97xd/def2tzvp',
                   irc_level='wb97xd/def2tzvp'),
            'network1_1', _STRUCTURES)
        assert arc_input['sp_level'] == 'ccsd(t)-f12/cc-pvtz-f12'
        assert arc_input['scan_level'] == 'wb97xd/def2tzvp'
        assert arc_input['irc_level'] == 'wb97xd/def2tzvp'
        assert arc_input['reactions'][0]['family'] == '1,2_Insertion_CO'

    def test_reaction_and_species_dicts_construct_a_real_arc_reaction(self):
        """N1: build_arc_input must not set a 'label' key on the reaction dict it emits.
        ARCReaction.from_dict only synthesizes a label from reactants/products when self.label
        starts empty (arc/reaction/reaction.py, set_label_reactants_products) -- a raw network
        label like 'reaction1' passed through as 'label' short-circuits that fallback, and
        check_atom_balance then crashes with IndexError on self.label.split('<=>')[well] because
        'reaction1' has no '<=>' in it. This constructs a REAL ARCReaction (mirroring ARC.__init__'s
        own order: ARCSpecies built first, then ARCReaction(reaction_dict=..., species_list=...))
        to prove the emitted dict actually survives ARC's construction path, not just that it has
        the right shape."""
        candidate = _candidate(ts_label='TS1', reactants=('A',), products=('B',), label='reaction1')
        arc_input = build_arc_input((candidate,), round_paths('/proj', 0), _config(), 'network1_1',
                                    _BALANCED_STRUCTURES)
        reaction = arc_input['reactions'][0]
        assert 'label' not in reaction

        species_list = [ARCSpecies(**spc) for spc in arc_input['species']]
        rxn = ARCReaction(reaction_dict=reaction, species_list=species_list)

        assert rxn.ts_label == reaction['ts_label']
        assert [spc.label for spc in rxn.r_species] == ['A']
        assert [spc.label for spc in rxn.p_species] == ['B']

    def test_compute_thermo_is_explicitly_disabled_on_every_species(self):
        """ARC defaults a bare species dict's compute_thermo to True, which would queue full
        thermo for every reactant/product reaching ARC through this path -- unbounded by
        config.qm.max_transition_states_per_round. This loop only ever needs the TS itself."""
        candidate = _candidate()
        arc_input = build_arc_input((candidate,), round_paths('/proj', 0), _config(), 'network1_1',
                                    _STRUCTURES)
        assert arc_input['species']
        assert all(spc['compute_thermo'] is False for spc in arc_input['species'])

    def test_candidate_with_no_structure_is_dropped_not_crashed(self, caplog):
        candidate = _candidate(reactants=('A',), products=('unknown_species',))
        with caplog.at_level('WARNING', logger='t3.pdep.pes_qm'):
            arc_input = build_arc_input((candidate,), round_paths('/proj', 0), _config(),
                                        'network1_1', _STRUCTURES)
        assert arc_input['reactions'] == []
        assert arc_input['species'] == []
        warnings = [r.message for r in caplog.records if r.levelname == 'WARNING']
        assert any('unknown_species' in message for message in warnings), \
            f'expected a WARNING naming the missing species, got: {warnings}'

    def test_dead_network_ts_label_key_is_not_emitted(self):
        candidate = _candidate()
        arc_input = build_arc_input((candidate,), round_paths('/proj', 0), _config(), 'network1_1',
                                    _STRUCTURES)
        assert 'network_ts_label' not in arc_input['reactions'][0]


# --- arc_qm_runner ------------------------------------------------------------------------------
#
# Real fixture network (tests/data/pdep_real_networks/network799_1/network799_1.py): 3 reactions
# (TS1/TS2/TS3), matching real adjacency-list species. arc_qm_runner is exercised against argument-
# recording doubles for every impure boundary it crosses (ARC, capture_ts_artifacts,
# write_hybrid_network_input_file, save_yaml_file) rather than canned-return stubs, so each test
# asserts on what arc_qm_runner actually HANDED those boundaries, not merely that it called them.

_NETWORK_ID = 'network799_1'
_FIXTURE_NETWORK_PATH = os.path.join(
    os.path.dirname(__file__), '..', 'data', 'pdep_real_networks', 'network799_1', 'network799_1.py')

_FROZEN_ENERGY_SETTINGS = {
    'model_chemistry': 'wb97xd/def2tzvp',
    'frequency_scale_factor': None,
    'use_hindered_rotors': True,
    'use_bond_corrections': False,
    'bond_correction_type': None,
    'atom_energies': None,
    'use_atom_corrections': True,
    'bond_additivity_corrections': None,
}


class _FakeARC(object):
    """Records every construction's kwargs and every .execute()/.as_dict() call, class-wide."""

    constructions = []
    execute_calls = 0

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        _FakeARC.constructions.append(kwargs)

    def execute(self):
        _FakeARC.execute_calls += 1

    def as_dict(self):
        # N5 (round 3): must be DERIVED from self.kwargs but not equal to it, or a test
        # asserting content == _FakeARC.constructions[0] can no longer tell 'the persisted
        # content came from arc.as_dict()' apart from 'the persisted content is just
        # arc_kwargs verbatim' -- which is exactly the I2 defect (persisting arc_kwargs
        # instead of calling as_dict()) that round 2's tightening silently stopped catching.
        return dict(self.kwargs) | {'_from_as_dict': True}


class _FakeCalls(object):
    """Argument recorder shared by the capture_ts_artifacts/write_hybrid_network_input_file/
    save_yaml_file doubles below -- each call's kwargs are appended verbatim."""

    def __init__(self):
        self.capture_ts_artifacts_calls = []
        self.write_hybrid_calls = []
        self.save_yaml_file_calls = []
        self.capture_result = None
        self.hybrid_result = None

    def fake_capture_ts_artifacts(self, **kwargs):
        self.capture_ts_artifacts_calls.append(kwargs)
        return self.capture_result

    def fake_write_hybrid_network_input_file(self, **kwargs):
        self.write_hybrid_calls.append(kwargs)
        return self.hybrid_result

    def fake_save_yaml_file(self, **kwargs):
        self.save_yaml_file_calls.append(kwargs)


def _arc_qm_runner_fixture(tmp_path, monkeypatch, capture_result):
    """Build real RoundPaths/candidates against the real network799_1 fixture, monkeypatch every
    impure boundary arc_qm_runner crosses with argument-recording doubles, and return the recorder
    plus the RoundPaths used, so a test can drive arc_qm_runner and then inspect exactly what it
    handed each boundary."""
    project_directory = str(tmp_path)
    paths = round_paths(project_directory, 0)
    final_dir = os.path.join(paths.explorer_output, 'pdep', 'final')
    os.makedirs(final_dir)
    dest_network_path = os.path.join(final_dir, f'{_NETWORK_ID}.py')
    shutil.copyfile(_FIXTURE_NETWORK_PATH, dest_network_path)

    candidate = _candidate(ts_label='TS1', family='1+2_Cycloaddition',
                           reactants=('O-2(13598)', 'CO2(11)'), products=('O=C1OO1(21648)',),
                           label='reaction1')

    _FakeARC.constructions = []
    _FakeARC.execute_calls = 0
    monkeypatch.setattr(pes_qm, 'ARC', _FakeARC)

    recorder = _FakeCalls()
    recorder.capture_result = capture_result
    monkeypatch.setattr(pes_qm, 'capture_ts_artifacts', recorder.fake_capture_ts_artifacts)
    monkeypatch.setattr(pes_qm, 'write_hybrid_network_input_file',
                        recorder.fake_write_hybrid_network_input_file)
    monkeypatch.setattr(pes_qm, 'save_yaml_file', recorder.fake_save_yaml_file)

    return paths, (candidate,), recorder, dest_network_path


def _usable_capture_result(paths, capture_dir, network_path):
    return CaptureResult(
        capture_dir=capture_dir,
        manifest_path=os.path.join(capture_dir, 'manifest.yml'),
        records=(
            TSArtifactRecord(network_id=_NETWORK_ID, network_ts_label='TS1',
                             arc_ts_label='T3PDep_network799_1_TS1', status=ARTIFACT_STATUS_USABLE,
                             artifact_path=os.path.join(capture_dir, 'qm', 'TS1.py'), converged=True),
        ),
        energy_settings=_FROZEN_ENERGY_SETTINGS,
        networks={
            _NETWORK_ID: {
                'source_path': network_path,
                'source_sha256': 'deadbeef',
                'captured_path': os.path.join('networks', f'{_NETWORK_ID}.py'),
                'method': 'MSC',
            },
        },
    )


def _empty_capture_result(capture_dir):
    return CaptureResult(
        capture_dir=capture_dir,
        manifest_path=os.path.join(capture_dir, 'manifest.yml'),
        records=(
            TSArtifactRecord(network_id=_NETWORK_ID, network_ts_label='TS1', arc_ts_label=None,
                             status=ARTIFACT_STATUS_UNVERIFIED,
                             artifact_path=os.path.join(capture_dir, 'qm', 'TS1.py')),
        ),
        energy_settings=None,
        networks=None,
    )


class TestArcQmRunner(object):

    def test_arc_is_constructed_with_the_reactions_and_species_build_arc_input_built(self, tmp_path, monkeypatch):
        """Regression for C1: reactions/species must reach ARC(**arc_kwargs), not be dropped."""
        capture_dir = os.path.join(str(tmp_path), 'capture')
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(
            tmp_path, monkeypatch, _empty_capture_result(capture_dir))
        arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert len(_FakeARC.constructions) == 1
        kwargs = _FakeARC.constructions[0]
        assert kwargs['reactions'], 'ARC received zero reactions'
        assert kwargs['species'], 'ARC received zero species'
        assert kwargs['reactions'][0]['ts_label'] == 'T3PDep_network799_1_TS1'

    def test_arc_is_constructed_with_this_rounds_own_project_and_project_directory(self, tmp_path, monkeypatch):
        """#9: project/project_directory must be asserted on the recorded ARC(**kwargs) call
        itself, not merely inferred from build_arc_input's own output."""
        capture_dir = os.path.join(str(tmp_path), 'capture')
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(
            tmp_path, monkeypatch, _empty_capture_result(capture_dir))
        arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        kwargs = _FakeARC.constructions[0]
        assert kwargs['project'] == _NETWORK_ID
        assert kwargs['project_directory'] == paths.arc_project

    def test_arc_execute_is_called_exactly_once(self, tmp_path, monkeypatch):
        capture_dir = os.path.join(str(tmp_path), 'capture')
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(
            tmp_path, monkeypatch, _empty_capture_result(capture_dir))
        arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert _FakeARC.execute_calls == 1

    def test_arc_input_yml_is_persisted_from_arc_as_dict_after_construction(self, tmp_path, monkeypatch):
        """I2: persisted only after ARC(**arc_kwargs) exists, from arc.as_dict() -- not from
        arc_kwargs itself. _FakeARC.as_dict() returns arc_kwargs PLUS a '_from_as_dict' marker
        (N5, round 3) specifically so this assertion can tell the two apart: content must carry
        the marker (proving it came from as_dict()) AND match arc_kwargs on every other key
        (proving as_dict() wasn't just handed something unrelated)."""
        capture_dir = os.path.join(str(tmp_path), 'capture')
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(
            tmp_path, monkeypatch, _empty_capture_result(capture_dir))
        arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert len(recorder.save_yaml_file_calls) == 1
        call = recorder.save_yaml_file_calls[0]
        assert call['path'] == os.path.join(paths.arc_project, ARC_INPUT_FILE_NAME)
        assert call['content'].get('_from_as_dict') is True
        assert call['content'] == dict(_FakeARC.constructions[0]) | {'_from_as_dict': True}

    def test_arc_input_yml_is_not_rewritten_when_it_already_exists(self, tmp_path, monkeypatch):
        capture_dir = os.path.join(str(tmp_path), 'capture')
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(
            tmp_path, monkeypatch, _empty_capture_result(capture_dir))
        os.makedirs(paths.arc_project)
        with open(os.path.join(paths.arc_project, ARC_INPUT_FILE_NAME), 'w') as f:
            f.write('already here')
        arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert recorder.save_yaml_file_calls == []

    def test_capture_ts_artifacts_is_called_with_this_rounds_own_project_and_capture_dirs(self, tmp_path, monkeypatch):
        capture_dir = os.path.join(str(tmp_path), 'capture')
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(
            tmp_path, monkeypatch, _empty_capture_result(capture_dir))
        arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert len(recorder.capture_ts_artifacts_calls) == 1
        call = recorder.capture_ts_artifacts_calls[0]
        assert call['arc_project_directory'] == paths.arc_project
        assert call['capture_dir'] == paths.capture
        assert call['networks'][_NETWORK_ID]['source_path'] == network_path
        assert call['networks'][_NETWORK_ID]['method'] == 'MSC'
        assert len(call['join_records']) == 1
        join_record = call['join_records'][0]
        assert join_record.network_ts_label == 'TS1'
        assert join_record.arc_ts_label == arc_ts_label(_NETWORK_ID, 'TS1')
        assert join_record.network_id == _NETWORK_ID
        assert join_record.status == JOIN_STATUS_QUEUED

    def test_no_usable_artifact_returns_empty_and_skips_hybrid_write(self, tmp_path, monkeypatch):
        capture_dir = os.path.join(str(tmp_path), 'capture')
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(
            tmp_path, monkeypatch, _empty_capture_result(capture_dir))
        converged, queued = arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert converged == frozenset()
        assert queued == frozenset({'TS1'})
        assert recorder.write_hybrid_calls == []

    def test_dropped_candidate_is_not_reported_as_queued(self, tmp_path, monkeypatch):
        """N3: build_arc_input silently drops a candidate with no adjacency list for a
        reactant/product (a warning, not a crash) -- arc_qm_runner's own queued_ts_labels must not
        claim that dropped candidate was ever sent to ARC, even though it was among the candidates
        this function was handed."""
        capture_dir = os.path.join(str(tmp_path), 'capture')
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(
            tmp_path, monkeypatch, _empty_capture_result(capture_dir))
        dropped = _candidate(ts_label='TS_missing', family='1,2_Insertion_CO',
                             reactants=('no_such_species',), products=('also_missing',),
                             label='reaction_missing')
        candidates = candidates + (dropped,)
        converged, queued = arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert queued == frozenset({'TS1'})
        assert 'TS_missing' not in queued
        # The dropped candidate must not have reached ARC either.
        assert len(_FakeARC.constructions) == 1
        ts_labels_sent_to_arc = {r['ts_label'] for r in _FakeARC.constructions[0]['reactions']}
        assert 'TS_missing' not in {label.rsplit('_', 1)[-1] for label in ts_labels_sent_to_arc}
        # And the join record capture_ts_artifacts was handed must not claim it was queued either.
        assert len(recorder.capture_ts_artifacts_calls) == 1
        join_ts_labels = {jr.network_ts_label
                          for jr in recorder.capture_ts_artifacts_calls[0]['join_records']}
        assert join_ts_labels == {'TS1'}

    def test_usable_artifact_triggers_hybrid_write_from_the_captures_own_network_copy(self, tmp_path, monkeypatch):
        """I1: source_path/method must come from CaptureResult.networks, never the live network."""
        capture_dir = os.path.join(str(tmp_path), 'capture')
        os.makedirs(os.path.join(capture_dir, 'networks'))
        shutil.copyfile(_FIXTURE_NETWORK_PATH,
                        os.path.join(capture_dir, 'networks', f'{_NETWORK_ID}.py'))
        capture_result = _usable_capture_result(None, capture_dir, None)
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(tmp_path, monkeypatch, capture_result)
        recorder.hybrid_result = HybridNetworkResult(
            dest_path=hybrid_network_path(paths, _NETWORK_ID), qm_ts_labels=('TS1',), ilt_ts_labels=(),
            vendored_files=(), warnings=())
        arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert len(recorder.write_hybrid_calls) == 1
        call = recorder.write_hybrid_calls[0]
        assert call['source_path'] == os.path.join(capture_dir, 'networks', f'{_NETWORK_ID}.py')
        assert call['method'] == 'MSC'
        assert call['qm_artifacts_root'] == capture_dir

    def test_hybrid_write_dest_path_is_the_loops_own_hybrid_network_path(self, tmp_path, monkeypatch):
        capture_dir = os.path.join(str(tmp_path), 'capture')
        os.makedirs(os.path.join(capture_dir, 'networks'))
        shutil.copyfile(_FIXTURE_NETWORK_PATH,
                        os.path.join(capture_dir, 'networks', f'{_NETWORK_ID}.py'))
        capture_result = _usable_capture_result(None, capture_dir, None)
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(tmp_path, monkeypatch, capture_result)
        recorder.hybrid_result = HybridNetworkResult(
            dest_path=hybrid_network_path(paths, _NETWORK_ID), qm_ts_labels=('TS1',), ilt_ts_labels=(),
            vendored_files=(), warnings=())
        arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert recorder.write_hybrid_calls[0]['dest_path'] == hybrid_network_path(paths, _NETWORK_ID)

    def test_missing_network_entry_in_capture_result_raises_keyerror(self, tmp_path, monkeypatch):
        """Round 2's guard: ARC converged transition states for network_id, but
        CaptureResult.networks has no entry for it (e.g. a concurrent/prior capture vendored a
        different network under this id, or capture_ts_artifacts's own bookkeeping is stale).
        Writing the hybrid network file from a nonexistent captured source would either crash
        confusingly deep inside write_hybrid_network_input_file or -- worse -- silently read the
        wrong file; the guard must fail loudly and immediately instead."""
        capture_dir = os.path.join(str(tmp_path), 'capture')
        os.makedirs(os.path.join(capture_dir, 'networks'))
        capture_result = CaptureResult(
            capture_dir=capture_dir,
            manifest_path=os.path.join(capture_dir, 'manifest.yml'),
            records=(
                TSArtifactRecord(network_id=_NETWORK_ID, network_ts_label='TS1',
                                 arc_ts_label='T3PDep_network799_1_TS1', status=ARTIFACT_STATUS_USABLE,
                                 artifact_path=os.path.join(capture_dir, 'qm', 'TS1.py'), converged=True),
            ),
            energy_settings=_FROZEN_ENERGY_SETTINGS,
            networks={},  # no entry for _NETWORK_ID even though TS1 converged for it
        )
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(tmp_path, monkeypatch, capture_result)
        try:
            arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
            assert False, 'expected a KeyError'
        except KeyError as err:
            assert _NETWORK_ID in str(err)
        assert len(recorder.write_hybrid_calls) == 0

    def test_qm_transition_states_dict_only_carries_usable_records(self, tmp_path, monkeypatch):
        capture_dir = os.path.join(str(tmp_path), 'capture')
        os.makedirs(os.path.join(capture_dir, 'networks'))
        shutil.copyfile(_FIXTURE_NETWORK_PATH,
                        os.path.join(capture_dir, 'networks', f'{_NETWORK_ID}.py'))
        capture_result = CaptureResult(
            capture_dir=capture_dir, manifest_path=os.path.join(capture_dir, 'manifest.yml'),
            records=(
                TSArtifactRecord(network_id=_NETWORK_ID, network_ts_label='TS1',
                                 arc_ts_label='T3PDep_network799_1_TS1', status=ARTIFACT_STATUS_USABLE,
                                 artifact_path=os.path.join(capture_dir, 'qm', 'TS1.py'), converged=True),
                TSArtifactRecord(network_id=_NETWORK_ID, network_ts_label='TS2', arc_ts_label=None,
                                 status=ARTIFACT_STATUS_UNVERIFIED,
                                 artifact_path=os.path.join(capture_dir, 'qm', 'TS2.py')),
            ),
            energy_settings=_FROZEN_ENERGY_SETTINGS,
            networks={_NETWORK_ID: {'source_path': '', 'source_sha256': '', 'method': 'MSC',
                                    'captured_path': os.path.join('networks', f'{_NETWORK_ID}.py')}},
        )
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(tmp_path, monkeypatch, capture_result)
        recorder.hybrid_result = HybridNetworkResult(
            dest_path=hybrid_network_path(paths, _NETWORK_ID), qm_ts_labels=('TS1',), ilt_ts_labels=(),
            vendored_files=(), warnings=())
        arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        qm_transition_states = recorder.write_hybrid_calls[0]['qm_transition_states']
        assert qm_transition_states == {'TS1': os.path.join(capture_dir, 'qm', 'TS1.py')}

    def test_return_value_is_converged_and_queued_ts_labels_as_a_tuple(self, tmp_path, monkeypatch):
        """arc_qm_runner returns (converged_ts_labels, queued_ts_labels) -- N3."""
        capture_dir = os.path.join(str(tmp_path), 'capture')
        os.makedirs(os.path.join(capture_dir, 'networks'))
        shutil.copyfile(_FIXTURE_NETWORK_PATH,
                        os.path.join(capture_dir, 'networks', f'{_NETWORK_ID}.py'))
        capture_result = _usable_capture_result(None, capture_dir, None)
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(tmp_path, monkeypatch, capture_result)
        recorder.hybrid_result = HybridNetworkResult(
            dest_path=hybrid_network_path(paths, _NETWORK_ID), qm_ts_labels=('TS1',), ilt_ts_labels=(),
            vendored_files=(), warnings=())
        converged, queued = arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert converged == frozenset({'TS1'})
        assert queued == frozenset({'TS1'})
