#!/usr/bin/env python3

"""
t3 tests test_check_rmg_env module

Covers the pure logic of devtools/check_rmg_env.py: reading pins out of the rmg_env spec and judging a
child report against them. The child itself (a built rmg_env) is exercised in CI, by the workflow step
that runs the script; these tests need no rmg_env.
"""

import importlib.util
import os

import pytest

from t3.common import t3_path

SCRIPT_PATH = os.path.join(t3_path, 'devtools', 'check_rmg_env.py')
SPEC_PATH = os.path.join(t3_path, 'devtools', 'rmg_env.yml')

_spec = importlib.util.spec_from_file_location('check_rmg_env', SCRIPT_PATH)
check = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check)

SITE = '/envs/rmg_env/lib/python3.11/site-packages'
PINS = {'python': ('=', '3.11'), 'rmg': ('=', '4.0.0'), 'rmgdatabase': ('=', '4.0.0'),
        'numdifftools': ('>=', '0.9.41')}


def _good_report():
    """A child report that satisfies PINS in package mode."""
    return {'python': '3.11.18',
            'site_packages': [SITE],
            'installed': {'python': '3.11.18', 'rmg': '4.0.0', 'rmgdatabase': '4.0.0', 'numdifftools': '0.9.42'},
            'rmgpy_file': f'{SITE}/rmgpy/__init__.py',
            'rmgpy_version': '4.0.0',
            'arkane_file': f'{SITE}/arkane/__init__.py',
            'database': '/envs/rmg_env/share/rmgdatabase',
            'imports': {'t3/runners/rmg_incore_script.py:40': 'ok'}}


def _evaluate(report=None, decoy_report=None, rmg_path=None, rmg_db_path='/envs/rmg_env/share/rmgdatabase', rc=0):
    """Run evaluate() with the good defaults, overriding one input at a time."""
    report = report if report is not None else _good_report()
    decoy_report = decoy_report if decoy_report is not None else _good_report()
    return check.evaluate(PINS, report, decoy_report, rmg_path, rmg_db_path, rc)


def test_the_committed_spec_pins_python_and_rmg_in_a_form_the_check_can_read():
    """The committed spec parses, and every dependency is a pin the check can compare."""
    pins = check.load_spec(SPEC_PATH)
    assert pins['python'][0] == '='
    assert pins['rmg'][0] == '='
    assert all(op in ('=', '>=') for op, _ in pins.values())


def test_load_spec_reads_channel_qualified_and_floor_pins(tmp_path):
    """Channel prefixes are dropped and '>=' floors kept."""
    spec = tmp_path / 'env.yml'
    spec.write_text('name: e\nchannels: [rmg]\ndependencies:\n  - python=3.11\n  - rmg::rmg=4.0.0\n'
                    '  - conda-forge::numdifftools>=0.9.41\n')
    assert check.load_spec(str(spec)) == {'python': ('=', '3.11'), 'rmg': ('=', '4.0.0'),
                                          'numdifftools': ('>=', '0.9.41')}


@pytest.mark.parametrize('dependencies', ['  - python=3.11\n  - rmg\n',
                                          '  - python=3.11\n  - rmg>=4.0\n',
                                          '  - rmg=4.0.0\n',
                                          '  - python=3.11\n  - rmg=4.0.0\n  - numpy<2\n'])
def test_load_spec_refuses_a_spec_it_cannot_check(tmp_path, dependencies):
    """An unpinned or unparseable dependency, or a missing python/rmg pin, is refused, not skipped."""
    spec = tmp_path / 'env.yml'
    spec.write_text(f'name: e\ndependencies:\n{dependencies}')
    with pytest.raises(ValueError):
        check.load_spec(str(spec))


@pytest.mark.parametrize('op, pinned, installed, expected', [
    ('=', '3.11', '3.11.18', True),
    ('=', '3.11', '3.11', True),
    ('=', '3.1', '3.11.18', False),
    ('=', '4.0.0', '4.0.1', False),
    ('>=', '0.9.41', '0.9.42', True),
    ('>=', '0.9.41', '0.9.41', True),
    ('>=', '0.9.41', '0.9.39', False),
    ('>=', '0.9.41', '0.10.0', True),
    ('=', '4.0.0', None, False),
])
def test_pin_holds(op, pinned, installed, expected):
    """conda's '=' means V or V.*, never a longer number sharing the prefix; '>=' compares numerically."""
    assert check.pin_holds(op, pinned, installed) is expected


def test_evaluate_passes_a_report_that_matches_the_pins():
    """No failure for an environment that is the pin, in package mode."""
    assert _evaluate() == []


def test_evaluate_accepts_rmg_database_input_level():
    """RMG-database's git layout keeps the data under input/; that names the same database."""
    report = _good_report()
    report['database'] = '/envs/rmg_env/share/rmgdatabase/input'
    assert _evaluate(report=report) == []


@pytest.mark.parametrize('field, value, fragment', [
    ('python', '3.9.20', "child's Python is 3.9.20"),
    ('rmgpy_version', '3.3.0', 'rmgpy.__version__'),
    ('rmgpy_file', '/home/u/Code/RMG-Py/rmgpy/__init__.py', 'rmgpy loads from'),
    ('arkane_file', '/home/u/Code/RMG-Py/arkane/__init__.py', 'arkane loads from'),
    ('error', 'AttributeError: machar', 'importing rmgpy/arkane failed'),
    ('imports', {'t3/runners/pdep_source_seed.py:44': 'ModuleNotFoundError: x'}, 'pdep_source_seed.py:44'),
    ('imports', {}, 'no import statement'),
    ('database', '/home/u/Code/RMG-database/input', "ARC's RMG_DB_PATH"),
])
def test_evaluate_flags_each_departure_from_the_pin(field, value, fragment):
    """Each assertion fires on its own departure."""
    report = _good_report()
    report[field] = value
    failures = _evaluate(report=report)
    assert any(fragment in f for f in failures), failures


def test_evaluate_flags_an_installed_package_that_breaks_a_pin():
    """The rmg channel's numdifftools 0.9.39 is below the floor."""
    report = _good_report()
    report['installed']['numdifftools'] = '0.9.39'
    failures = _evaluate(report=report)
    assert any('numdifftools>=0.9.41' in f for f in failures), failures


def test_evaluate_flags_a_configured_rmg_path():
    """ARC with RMG_PATH set would put a source tree on the child's PYTHONPATH."""
    failures = _evaluate(rmg_path='/home/u/Code/RMG-Py')
    assert any('RMG_PATH=/home/u/Code/RMG-Py' in f for f in failures), failures


@pytest.mark.parametrize('decoy_update', [{'rmgpy_file': '/tmp/decoy_site/rmgpy/__init__.py',
                                           'rmgpy_version': 'decoy-from-parent-PYTHONPATH'},
                                          {'rmgpy_file': None}])
def test_evaluate_flags_a_parent_pythonpath_that_reaches_the_child(decoy_update):
    """The decoy run must import the same rmgpy as the plain run; a missing answer is a failure too."""
    decoy = _good_report()
    decoy.update(decoy_update)
    failures = _evaluate(decoy_report=decoy)
    assert any('PYTHONPATH set in the parent reached the child' in f for f in failures), failures


def test_evaluate_flags_a_failing_arkane_help():
    """`python -m arkane --help` must exit 0."""
    assert any('exited 1' in f for f in _evaluate(rc=1))


def test_import_snippets_collects_nested_imports_of_every_rmg_env_script():
    """Function-level imports count too (statmech_conformer_extract imports CompositeLevelOfTheory lazily)."""
    snippets = check.import_snippets()
    scripts = {label.rsplit(':', 1)[0] for label in snippets}
    assert scripts == set(check.RMG_ENV_SCRIPTS)
    assert any('CompositeLevelOfTheory' in code for code in snippets.values())
    assert any(code == 'from rmgpy.rmg.main import RMG' for code in snippets.values())
