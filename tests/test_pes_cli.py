"""Tests for the PES.py CLI."""

import pytest
import yaml

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
        args = parse_command_line_arguments(['/runs/pes/input.yml'])
        assert args.project_directory in (None, '/runs/pes')

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
