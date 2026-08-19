"""
Tests for t3.pdep.pes_sa -- the standalone PES loop's own source of real sensitivity evidence.

``test_run_round_me_sensitivity_runs_real_arkane`` runs the REAL Arkane master-equation SA (via
the real ``write_arkane_network_input_file`` and the real ``run_arkane_job``) on the real
``network799_1`` fixture -- deliberately no double at any seam, per the rule this branch keeps
re-learning: mock the thing under test and you test the mock. It costs ~8 s and is the only test
in the pdep suite that proves the loop's evidence pipeline end to end against Arkane itself.
"""

import math
import os
import shutil

import pytest

from t3.common import TEST_DATA_BASE_PATH
from t3.pdep import pes_sa
from t3.pdep.pes_sa import SA_COEFFICIENTS_RELPATH, run_round_me_sensitivity, ts_sensitivity_evidence
from t3.pdep.selector import E0_PERTURBATION_J_PER_MOL
from t3.pdep.yaml_safe import read_sa_yaml_file
from t3.runners.rmg_runner import ArkaneJobResult

_FIXTURE_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_real_networks', 'network799_1')
_FIXTURE_NETWORK_PATH = os.path.join(_FIXTURE_DIR, 'network799_1.py')
_FIXTURE_SA_PATH = os.path.join(_FIXTURE_DIR, 'sensitivity', 'sa_coefficients.yml')


class TestTsSensitivityEvidence(object):

    def test_real_recorded_sa_output_yields_max_abs_coefficient_per_ts(self):
        """Against a REAL recorded Arkane SA output. The pinned TS3 value comes from the SECOND
        direction key of the file (its first-key row is ~1.3e-08), so a mutant that only reads one
        direction key -- or the first -- fails this test."""
        evidence = ts_sensitivity_evidence(read_sa_yaml_file(_FIXTURE_SA_PATH))
        assert set(evidence) == {'TS1', 'TS2', 'TS3'}
        assert evidence['TS1'][0] == pytest.approx(-9.515582575266414e-05)
        assert evidence['TS2'][0] == pytest.approx(-8.660634384118037e-05)
        assert evidence['TS3'][0] == pytest.approx(-8.83268679433372e-05)
        for coefficient, delta_ln_k in evidence.values():
            assert delta_ln_k == pytest.approx(abs(coefficient) * E0_PERTURBATION_J_PER_MOL)

    def test_non_dict_payload_is_refused_not_read_as_no_evidence(self):
        with pytest.raises(ValueError, match='malformed'):
            ts_sensitivity_evidence(['not', 'a', 'dict'])

    def test_malformed_and_non_finite_pieces_are_skipped_never_guessed(self):
        sa_dict = {
            'structures': {'A': 'adjlist'},
            'A <=> B': {
                (300.0, 'K', 1.0, 'bar'): {'(TS) TS1': float('nan'),
                                           '(TS) TS2': 2.0e-5,
                                           'A': 1.0e-3},
                (600.0, 'K', 1.0, 'bar'): 'not a dict',
            },
            'A <=> C': 'not a dict either',
        }
        evidence = ts_sensitivity_evidence(sa_dict)
        assert evidence == {'TS2': (2.0e-5, pytest.approx(2.0e-5 * E0_PERTURBATION_J_PER_MOL))}

    def test_a_zero_coefficient_is_real_evidence_not_a_missing_row(self):
        """0.0 is a measured 'no leverage', which is finite, honest evidence -- distinct from a
        missing or NaN row, which is no evidence at all."""
        sa_dict = {'A <=> B': {(300.0, 'K', 1.0, 'bar'): {'(TS) TS1': 0.0}}}
        assert ts_sensitivity_evidence(sa_dict) == {'TS1': (0.0, 0.0)}


class TestRunRoundMeSensitivity(object):

    def test_run_round_me_sensitivity_runs_real_arkane(self, tmp_path):
        """The real pipeline, no double anywhere: real input rendering (sensitivity_conditions
        injected into a copy of the real network799_1), real Arkane ME SA run, real YAML read,
        real reduction. Every coefficient must be finite -- these are measured derivatives, and
        capture_ts_artifacts downstream refuses anything else."""
        network_path = os.path.join(str(tmp_path), 'network0_full.py')
        shutil.copyfile(_FIXTURE_NETWORK_PATH, network_path)
        sa_dir = os.path.join(str(tmp_path), 'SA')
        evidence = run_round_me_sensitivity(network_path=network_path, sa_dir=sa_dir,
                                            method='MSC', timeout=600)
        assert set(evidence) == {'TS1', 'TS2', 'TS3'}
        for coefficient, delta_ln_k in evidence.values():
            assert math.isfinite(coefficient) and math.isfinite(delta_ln_k)
            assert delta_ln_k == pytest.approx(abs(coefficient) * E0_PERTURBATION_J_PER_MOL)
        # The rendered input really carried the SA directive, and Arkane really (re)wrote the
        # coefficients artifact run_arkane_job requires.
        with open(os.path.join(sa_dir, 'input.py')) as f:
            assert 'sensitivity_conditions' in f.read()
        assert os.path.isfile(os.path.join(sa_dir, SA_COEFFICIENTS_RELPATH))

    def test_a_failed_arkane_job_raises_with_its_reason(self, tmp_path, monkeypatch):
        network_path = os.path.join(str(tmp_path), 'network0_full.py')
        shutil.copyfile(_FIXTURE_NETWORK_PATH, network_path)

        def _fake_run_arkane_job(*, input_file, output_directory, logger, timeout):
            return ArkaneJobResult(succeeded=False, reason='the solver exploded')

        monkeypatch.setattr(pes_sa, 'run_arkane_job', _fake_run_arkane_job)
        with pytest.raises(ValueError, match='the solver exploded'):
            run_round_me_sensitivity(network_path=network_path,
                                     sa_dir=os.path.join(str(tmp_path), 'SA'), method='MSC')
