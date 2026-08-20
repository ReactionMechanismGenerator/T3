"""Tests for the PES.py CLI."""

import os
import sys
from types import SimpleNamespace

import pytest
import yaml

import PES
from PES import exit_code_for, parse_command_line_arguments, read_pes_input
from t3.pdep.pes_loop import (PES_LOOP_CONVERGED, PES_LOOP_DIAGRAM_ONLY, PES_LOOP_FAILED,
                              PES_LOOP_MAX_ROUNDS, PES_LOOP_NO_CANDIDATES, PES_LOOP_STALLED)


class TestReadPESInput(object):

    def test_reads_a_valid_input_file(self, tmp_path):
        path = tmp_path / 'input.yml'
        path.write_text(yaml.dump({'pes': {'network': '/abs/n.py', 'source': ['HOCHO']}}))
        config = read_pes_input(str(path))
        assert config.pes.source == ['HOCHO']

    def test_an_unknown_top_level_key_is_refused(self, tmp_path):
        """extra='forbid': a typo'd key must fail loudly, not be silently ignored, or a user
        thinks they configured something they did not."""
        path = tmp_path / 'input.yml'
        path.write_text(yaml.dump({'pes': {'network': '/abs/n.py', 'source': ['A']},
                                   'termnation': {'max_rounds': 2}}))
        with pytest.raises(Exception):
            read_pes_input(str(path))

    def test_a_missing_file_raises_a_clear_error(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_pes_input(str(tmp_path / 'nope.yml'))


class TestParseCommandLineArguments(object):

    def test_input_file_is_positional(self):
        args = parse_command_line_arguments(['input.yml'])
        assert args.file == 'input.yml'

    def test_project_directory_defaults_to_the_inputs_directory(self):
        # Unset at the argparse level, resolved by main() -- the design the addendum settled on.
        args = parse_command_line_arguments(['/runs/pes/input.yml'])
        assert args.project_directory is None

    def test_project_directory_can_be_given_explicitly(self):
        args = parse_command_line_arguments(['/runs/pes/input.yml', '-p', '/scratch/run'])
        assert args.project_directory == '/scratch/run'


class TestExitCodeFor(object):
    """A stall or a max_rounds stop is a finding to read in the log, not a crash: only an
    outright failure is worth a non-zero exit code."""

    @pytest.mark.parametrize('status', [PES_LOOP_CONVERGED,
                                        PES_LOOP_MAX_ROUNDS,
                                        PES_LOOP_NO_CANDIDATES,
                                        PES_LOOP_STALLED,
                                        PES_LOOP_DIAGRAM_ONLY])
    def test_every_non_failure_status_exits_zero(self, status):
        assert exit_code_for(status) == 0

    def test_a_failed_loop_exits_one(self):
        assert exit_code_for(PES_LOOP_FAILED) == 1


class _RecordingLogger(object):
    """Stands in for t3.logger.Logger at exactly the boundary main() constructs it, so the
    verbosity flags can be checked against what the Logger was actually handed."""

    instances = []

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        _RecordingLogger.instances.append(self)

    def log(self, message, level='info'):
        pass

    def log_footer(self):
        pass


class TestVerbosity(object):
    """``-d``/``-q`` are only worth anything if the level they map to reaches the Logger main()
    builds; asserting the mapping in isolation would not notice main() ignoring it."""

    @pytest.mark.parametrize('flags, expected_verbose', [([], 20), (['-d'], 10), (['-q'], 30)])
    def test_the_flags_reach_the_logger(self, tmp_path, monkeypatch, flags, expected_verbose):
        input_path = tmp_path / 'input.yml'
        input_path.write_text(yaml.dump({'pes': {'network': '/abs/n.py', 'source': ['HOCHO']}}))
        _RecordingLogger.instances = []
        monkeypatch.setattr(PES, 'Logger', _RecordingLogger)
        # The loop itself is out of scope here -- it is driven for real, end to end, by
        # tests/test_pdep/test_pes_loop_integration.py.
        monkeypatch.setattr(PES, 'run_pes_loop',
                            lambda *args, **kwargs: SimpleNamespace(
                                status=PES_LOOP_CONVERGED, reason='', rounds=(),
                                final_network_path=None, final_diagram_path=None))
        monkeypatch.setattr(sys, 'argv', ['PES.py', str(input_path)] + flags)

        PES.main()

        (logger,) = _RecordingLogger.instances
        assert logger.kwargs['verbose'] == expected_verbose
        assert logger.kwargs['project_directory'] == os.path.abspath(str(tmp_path))
