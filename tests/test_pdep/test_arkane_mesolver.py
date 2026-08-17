#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_arkane_mesolver module

Tests the concrete ``ArkaneMESolverAdapter``. The Arkane subprocess itself is never invoked:
``t3.pdep.mesolver.arkane.run_arkane_job`` is monkeypatched with a fake that copies one of the
``tests/data/pdep_me/`` fixtures into place as ``output.py``, mirroring the
``monkeypatch``/``tmp_path`` idiom used in
``tests/test_runners/test_rmg_runner.py::TestRunArkaneJob``.

Fixtures live under ``tests/data/pdep_me/`` and ``tests/data/pdep_network/``; the latter is a
``t3.paths['PDep SA']``/``t3.paths['RMG PDep']`` target elsewhere in the suite and gets
``shutil.rmtree``'d in teardown by other tests (see e.g. ``tests/test_main.py``), so nothing is
ever written INTO those directories here -- only read from. All adapter output is written under
``tmp_path``.
"""

import os
import shutil

import pytest

from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.mesolver.arkane import ArkaneMESolverAdapter
from t3.pdep.mesolver.factory import _registered_mesolver_adapters, mesolver_factory
from t3.utils.writer import METHOD_MAP

PDEP_ME_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_me')
SUCCESS_OUTPUT = os.path.join(PDEP_ME_DIR, 'success', 'output.py')
SOFT_FAILURE_CSE_OUTPUT = os.path.join(PDEP_ME_DIR, 'soft_failure_cse', 'output.py')

NETWORK_FIXTURE = os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'iteration_1', 'RMG', 'pdep',
                               'network1_1.py')


def _fake_run_arkane_job_writing(fixture_output_path):
    """
    Build a fake ``run_arkane_job`` that reports success and copies ``fixture_output_path`` into
    ``output_directory`` as ``output.py``, WITHOUT writing any ``stderr.log`` -- reproducing the
    empirically-observed soft-failure case: exit 0, empty stderr, broken payload.
    """
    def _fake(input_file, output_directory, plot=False, logger=None, required_artifact='output.py'):
        os.makedirs(output_directory, exist_ok=True)
        shutil.copyfile(fixture_output_path, os.path.join(output_directory, 'output.py'))
        return True
    return _fake


class TestArkaneMESolverAdapterRegistration(object):

    def test_registered_under_arkane(self):
        """ArkaneMESolverAdapter is registered under the 'Arkane' key."""
        assert _registered_mesolver_adapters['Arkane'] is ArkaneMESolverAdapter

    def test_factory_returns_arkane_adapter_instance(self, tmp_path):
        """mesolver_factory(me_solver='Arkane', ...) returns an ArkaneMESolverAdapter instance."""
        output_directory = str(tmp_path / 'out')
        adapter = mesolver_factory(me_solver='Arkane',
                                   network_path=NETWORK_FIXTURE,
                                   output_directory=output_directory,
                                   method='CSE',
                                   )
        assert isinstance(adapter, ArkaneMESolverAdapter)

    def test_supports_ilt_complement_is_true(self):
        """Arkane does not require every network node to be QM-computed."""
        assert ArkaneMESolverAdapter.supports_ilt_complement is True

    def test_factory_allow_ilt_complement_true_succeeds_for_arkane(self, tmp_path):
        """allow_ilt_complement=True is accepted for the Arkane adapter (unlike MESS/Mesmer)."""
        output_directory = str(tmp_path / 'out')
        adapter = mesolver_factory(me_solver='Arkane',
                                   network_path=NETWORK_FIXTURE,
                                   output_directory=output_directory,
                                   method='CSE',
                                   allow_ilt_complement=True,
                                   )
        assert isinstance(adapter, ArkaneMESolverAdapter)
        assert adapter.allow_ilt_complement is True


class TestArkaneMESolverAdapterSetUp(object):

    @pytest.mark.parametrize('shorthand', ['CSE', 'RS', 'MSC'])
    def test_set_up_writes_correct_method_and_no_sensitivity_directive(self, tmp_path, shorthand):
        """set_up() writes the Arkane method string via METHOD_MAP and never injects sensitivity_conditions."""
        output_directory = str(tmp_path / 'out')
        adapter = ArkaneMESolverAdapter(network_path=NETWORK_FIXTURE,
                                        output_directory=output_directory,
                                        method=shorthand,
                                        )
        input_file_path = os.path.join(output_directory, 'input.py')
        assert os.path.isfile(input_file_path)
        with open(input_file_path, 'r') as f:
            content = f.read()
        assert f"method = '{METHOD_MAP[shorthand]}'" in content
        assert 'sensitivity_conditions' not in content
        assert adapter.isomer_labels == ('C1rad(2)',)

    def test_invalid_method_raises_value_error_naming_valid_keys(self, tmp_path):
        """An invalid method raises ValueError naming the valid METHOD_MAP keys."""
        output_directory = str(tmp_path / 'out')
        with pytest.raises(ValueError) as excinfo:
            ArkaneMESolverAdapter(network_path=NETWORK_FIXTURE,
                                 output_directory=output_directory,
                                 method='not-a-real-method',
                                 )
        message = str(excinfo.value)
        for key in ('CSE', 'RS', 'MSC'):
            assert key in message


class TestArkaneMESolverAdapterSolve(object):

    def _make_adapter(self, tmp_path):
        output_directory = str(tmp_path / 'out')
        return ArkaneMESolverAdapter(network_path=NETWORK_FIXTURE,
                                     output_directory=output_directory,
                                     method='CSE',
                                     )

    def test_solve_returns_true_for_success_fixture(self, tmp_path, monkeypatch):
        """solve() returns True when the run leaves the SUCCESS fixture's content as output.py."""
        adapter = self._make_adapter(tmp_path)
        monkeypatch.setattr('t3.pdep.mesolver.arkane.run_arkane_job',
                            _fake_run_arkane_job_writing(SUCCESS_OUTPUT))
        assert adapter.solve() is True
        assert adapter.me_success_result.succeeded is True

    def test_stale_stderr_log_does_not_poison_next_run(self, tmp_path, monkeypatch):
        """A stderr.log left behind by a previous failed run must not fail a subsequent good solve.

        Regression this guards: solve() reads ``<output_directory>/stderr.log`` after the run. A
        failed run leaves one behind; if the next run in the same directory writes a clean
        output.py and produces no new stderr, the adapter used to read the OLD stderr.log and
        fail a good solve. The stale file must be deleted before invoking the run, mirroring the
        delete-then-require handling of the artifact itself in run_arkane_job().
        """
        adapter = self._make_adapter(tmp_path)
        stale_stderr_path = os.path.join(adapter.output_directory, 'stderr.log')
        with open(stale_stderr_path, 'w') as f:
            f.write('Traceback (most recent call last):\nrmgpy.exceptions.NetworkError: boom\n')
        monkeypatch.setattr('t3.pdep.mesolver.arkane.run_arkane_job',
                            _fake_run_arkane_job_writing(SUCCESS_OUTPUT))
        assert adapter.solve() is True
        assert adapter.me_success_result.succeeded is True
        assert not os.path.isfile(stale_stderr_path), \
            'The stale stderr.log should have been deleted before the run was invoked.'

    def test_solve_returns_false_for_soft_failure_cse_fixture(self, tmp_path, monkeypatch):
        """
        The most important test in this module.

        The fake run reports success (return True) and writes NO stderr.log, exactly reproducing
        the empirically-observed Arkane soft failure: exit 0, empty stderr, a syntactically
        perfect pdepreaction(...) block whose Chebyshev coefficients are all None. If solve()
        trusted the subprocess's own report of success, this would return True and silently hand
        back broken kinetics. It must return False.
        """
        adapter = self._make_adapter(tmp_path)
        monkeypatch.setattr('t3.pdep.mesolver.arkane.run_arkane_job',
                            _fake_run_arkane_job_writing(SOFT_FAILURE_CSE_OUTPUT))
        assert adapter.solve() is False
        assert adapter.me_success_result.succeeded is False
        assert adapter.me_success_result.reasons


class TestArkaneMESolverAdapterExpectedReactions(object):
    """
    Tests for the optional ``expected_reactions`` completeness check.

    Without it, a truncated output.py containing a single finite pdepreaction is accepted as a
    "solved" network and get_k_tp() returns partial kinetics. The count is NOT derived from the
    network topology (the number of net reactions in a PDep network is not trivially the number
    of path reactions); it is plumbed through as an explicit optional parameter.
    """

    def test_solve_returns_false_when_expected_reactions_not_met(self, tmp_path, monkeypatch):
        """solve() returns False when the output has fewer pdepreactions than expected.

        The SUCCESS fixture contains exactly one pdepreaction; expecting two must fail the
        solve even though the single entry present is perfectly finite.
        """
        output_directory = str(tmp_path / 'out')
        adapter = ArkaneMESolverAdapter(network_path=NETWORK_FIXTURE,
                                        output_directory=output_directory,
                                        method='CSE',
                                        expected_reactions=2,
                                        )
        monkeypatch.setattr('t3.pdep.mesolver.arkane.run_arkane_job',
                            _fake_run_arkane_job_writing(SUCCESS_OUTPUT))
        assert adapter.solve() is False
        assert any('Expected 2' in reason for reason in adapter.me_success_result.reasons)

    def test_solve_returns_true_when_expected_reactions_met(self, tmp_path, monkeypatch):
        """solve() returns True when the output meets the expected pdepreaction count."""
        output_directory = str(tmp_path / 'out')
        adapter = ArkaneMESolverAdapter(network_path=NETWORK_FIXTURE,
                                        output_directory=output_directory,
                                        method='CSE',
                                        expected_reactions=1,
                                        )
        monkeypatch.setattr('t3.pdep.mesolver.arkane.run_arkane_job',
                            _fake_run_arkane_job_writing(SUCCESS_OUTPUT))
        assert adapter.solve() is True

    def test_factory_passes_expected_reactions_through(self, tmp_path):
        """mesolver_factory() plumbs expected_reactions through to the adapter."""
        output_directory = str(tmp_path / 'out')
        adapter = mesolver_factory(me_solver='Arkane',
                                   network_path=NETWORK_FIXTURE,
                                   output_directory=output_directory,
                                   method='CSE',
                                   expected_reactions=3,
                                   )
        assert adapter.expected_reactions == 3

    def test_default_none_skips_the_completeness_check(self, tmp_path, monkeypatch):
        """expected_reactions=None (the default) behaves exactly as before: no completeness check."""
        output_directory = str(tmp_path / 'out')
        adapter = ArkaneMESolverAdapter(network_path=NETWORK_FIXTURE,
                                        output_directory=output_directory,
                                        method='CSE',
                                        )
        assert adapter.expected_reactions is None
        monkeypatch.setattr('t3.pdep.mesolver.arkane.run_arkane_job',
                            _fake_run_arkane_job_writing(SUCCESS_OUTPUT))
        assert adapter.solve() is True


class TestArkaneMESolverAdapterGetKTP(object):

    def test_get_k_tp_raises_before_solve(self, tmp_path):
        """get_k_tp() raises RuntimeError before solve() has been called."""
        output_directory = str(tmp_path / 'out')
        adapter = ArkaneMESolverAdapter(network_path=NETWORK_FIXTURE,
                                        output_directory=output_directory,
                                        method='CSE',
                                        )
        with pytest.raises(RuntimeError):
            adapter.get_k_tp()

    def test_get_k_tp_raises_after_failed_solve(self, tmp_path, monkeypatch):
        """get_k_tp() raises RuntimeError after a failed solve(), never returning broken kinetics silently."""
        output_directory = str(tmp_path / 'out')
        adapter = ArkaneMESolverAdapter(network_path=NETWORK_FIXTURE,
                                        output_directory=output_directory,
                                        method='CSE',
                                        )
        monkeypatch.setattr('t3.pdep.mesolver.arkane.run_arkane_job',
                            _fake_run_arkane_job_writing(SOFT_FAILURE_CSE_OUTPUT))
        assert adapter.solve() is False
        with pytest.raises(RuntimeError):
            adapter.get_k_tp()

    def test_get_k_tp_returns_reactions_after_successful_solve(self, tmp_path, monkeypatch):
        """get_k_tp() returns the parsed reactions after a successful solve()."""
        output_directory = str(tmp_path / 'out')
        adapter = ArkaneMESolverAdapter(network_path=NETWORK_FIXTURE,
                                        output_directory=output_directory,
                                        method='CSE',
                                        )
        monkeypatch.setattr('t3.pdep.mesolver.arkane.run_arkane_job',
                            _fake_run_arkane_job_writing(SUCCESS_OUTPUT))
        assert adapter.solve() is True
        reactions = adapter.get_k_tp()
        assert len(reactions) > 0
