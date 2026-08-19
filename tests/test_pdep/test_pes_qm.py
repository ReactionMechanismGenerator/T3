"""Tests for t3.pdep.pes_qm."""

import os
import shutil

import t3.pdep.pes_qm as pes_qm
from t3.pdep.capture import CaptureResult
from t3.pdep.discovery import ARTIFACT_STATUS_UNVERIFIED, ARTIFACT_STATUS_USABLE, TSArtifactRecord
from t3.pdep.hybrid import HybridNetworkResult
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

    def test_reaction_and_species_dicts_can_construct_a_real_arc_reaction(self):
        """C2: a reaction dict alone (no top-level 'species' list) cannot construct a valid
        ARCReaction -- ARCReaction.from_dict resolves r_species/p_species by matching their
        'label' against ARC's own top-level species list. This asserts build_arc_input emits
        both halves of that contract."""
        candidate = _candidate(reactants=('A',), products=('B',))
        arc_input = build_arc_input((candidate,), round_paths('/proj', 0), _config(), 'network1_1',
                                    _STRUCTURES)
        reaction = arc_input['reactions'][0]
        assert reaction['reactants'] == ['A']
        assert reaction['products'] == ['B']
        assert reaction['r_species'] == [{'label': 'A'}]
        assert reaction['p_species'] == [{'label': 'B'}]
        species_by_label = {spc['label']: spc for spc in arc_input['species']}
        assert set(species_by_label) == {'A', 'B'}
        assert species_by_label['A']['adjlist'] == _STRUCTURES['A']
        assert species_by_label['B']['adjlist'] == _STRUCTURES['B']

    def test_compute_thermo_is_explicitly_disabled_on_every_species(self):
        """ARC defaults a bare species dict's compute_thermo to True, which would queue full
        thermo for every reactant/product reaching ARC through this path -- unbounded by
        config.qm.max_transition_states_per_round. This loop only ever needs the TS itself."""
        candidate = _candidate()
        arc_input = build_arc_input((candidate,), round_paths('/proj', 0), _config(), 'network1_1',
                                    _STRUCTURES)
        assert arc_input['species']
        assert all(spc['compute_thermo'] is False for spc in arc_input['species'])

    def test_candidate_with_no_structure_is_dropped_not_crashed(self):
        candidate = _candidate(reactants=('A',), products=('unknown_species',))
        arc_input = build_arc_input((candidate,), round_paths('/proj', 0), _config(), 'network1_1',
                                    _STRUCTURES)
        assert arc_input['reactions'] == []
        assert arc_input['species'] == []

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
        return {'fake_as_dict_of': self.kwargs.get('project')}


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

    def test_arc_execute_is_called_exactly_once(self, tmp_path, monkeypatch):
        capture_dir = os.path.join(str(tmp_path), 'capture')
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(
            tmp_path, monkeypatch, _empty_capture_result(capture_dir))
        arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert _FakeARC.execute_calls == 1

    def test_arc_input_yml_is_persisted_from_arc_as_dict_after_construction(self, tmp_path, monkeypatch):
        """I2: persisted only after ARC(**arc_kwargs) exists, from arc.as_dict()."""
        capture_dir = os.path.join(str(tmp_path), 'capture')
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(
            tmp_path, monkeypatch, _empty_capture_result(capture_dir))
        arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert len(recorder.save_yaml_file_calls) == 1
        call = recorder.save_yaml_file_calls[0]
        assert call['path'] == os.path.join(paths.arc_project, ARC_INPUT_FILE_NAME)
        assert call['content'] == {'fake_as_dict_of': _NETWORK_ID}

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
        assert call['join_records'][0].network_ts_label == 'TS1'

    def test_no_usable_artifact_returns_empty_and_skips_hybrid_write(self, tmp_path, monkeypatch):
        capture_dir = os.path.join(str(tmp_path), 'capture')
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(
            tmp_path, monkeypatch, _empty_capture_result(capture_dir))
        result = arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert result == frozenset()
        assert recorder.write_hybrid_calls == []

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

    def test_return_value_is_the_frozenset_of_usable_network_ts_labels_only(self, tmp_path, monkeypatch):
        capture_dir = os.path.join(str(tmp_path), 'capture')
        os.makedirs(os.path.join(capture_dir, 'networks'))
        shutil.copyfile(_FIXTURE_NETWORK_PATH,
                        os.path.join(capture_dir, 'networks', f'{_NETWORK_ID}.py'))
        capture_result = _usable_capture_result(None, capture_dir, None)
        paths, candidates, recorder, network_path = _arc_qm_runner_fixture(tmp_path, monkeypatch, capture_result)
        recorder.hybrid_result = HybridNetworkResult(
            dest_path=hybrid_network_path(paths, _NETWORK_ID), qm_ts_labels=('TS1',), ilt_ts_labels=(),
            vendored_files=(), warnings=())
        result = arc_qm_runner(candidates, paths, _config(), _NETWORK_ID)
        assert result == frozenset({'TS1'})
