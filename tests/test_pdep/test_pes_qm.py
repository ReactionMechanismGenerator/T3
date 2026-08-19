"""Tests for t3.pdep.pes_qm."""

from t3.pdep.pes_qm import build_arc_input
from t3.pdep.pes_rounds import QMCandidate, round_paths
from t3.schema import PESLoopConfig


def _config(**qm) -> PESLoopConfig:
    return PESLoopConfig(pes={'network': '/abs/n.py', 'source': ['A']}, **({'qm': qm} if qm else {}))


class TestBuildARCInput(object):

    def test_levels_reach_arc_undashed(self):
        arc_input = build_arc_input((), round_paths('/proj', 0), _config(), 'network1_1')
        assert arc_input['opt_level'] == 'wb97xd/def2tzvp'
        assert '-' not in arc_input['opt_level'].split('/')[1]

    def test_rotors_and_irc_off_by_default(self):
        arc_input = build_arc_input((), round_paths('/proj', 0), _config(), 'network1_1')
        assert arc_input['job_types']['rotors'] is False
        assert arc_input['job_types']['irc'] is False

    def test_rotors_and_irc_can_be_turned_on(self):
        arc_input = build_arc_input((), round_paths('/proj', 0),
                                    _config(rotors=True, irc=True), 'network1_1')
        assert arc_input['job_types']['rotors'] is True
        assert arc_input['job_types']['irc'] is True

    def test_ts_adapters_are_passed_through(self):
        arc_input = build_arc_input((), round_paths('/proj', 0),
                                    _config(ts_adapters=['linear', 'rits']), 'network1_1')
        assert arc_input['ts_adapters'] == ['linear', 'rits']

    def test_project_directory_is_the_rounds_own_arc_project(self):
        paths = round_paths('/proj', 2)
        arc_input = build_arc_input((), paths, _config(), 'network1_1')
        assert arc_input['project_directory'] == paths.arc_project

    def test_ts_labels_are_namespaced_not_network_local(self):
        """A network-local 'TS1' collides across networks; join.arc_ts_label namespaces it so a
        network can never be joined to another network's quantum chemistry."""
        candidate = QMCandidate(path_reaction=None, ts_label='TS1', family='1,2_Insertion_CO')
        arc_input = build_arc_input((candidate,), round_paths('/proj', 0), _config(), 'network1_1')
        assert arc_input['reactions'][0]['ts_label'] != 'TS1'
        assert 'TS1' in arc_input['reactions'][0]['ts_label']
