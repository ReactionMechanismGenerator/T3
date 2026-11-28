#!/usr/bin/env python
"""
Assert that the built ``rmg_env`` is the one ``devtools/rmg_env.yml`` describes, in package mode.

Run it in ``t3_env``, with ARC importable, after ``rmg_env`` exists::

    python devtools/check_rmg_env.py [path/to/rmg_env.yml]

It reaches ``rmg_env`` the way T3 does in production -- through ARC's ``rmg_env_command`` launcher,
in a ``bash`` child -- and fails (exit 1, every failure listed) unless all of these hold:

1. Each pin in the spec matches the package installed in the child's environment: ``name=V`` means
   the installed version is ``V`` or ``V.*`` (conda's own meaning), ``name>=V`` means at least ``V``.
2. The child's Python is the spec's Python, ``rmgpy.__version__`` is the spec's ``rmg`` version,
   and ``rmgpy`` and ``arkane`` load from inside that environment's ``site-packages``.
3. ``RMG_PATH`` is unset in ARC's resolved settings, so ARC selects package mode.
4. A ``PYTHONPATH`` set in the parent does not reach the child: a decoy ``rmgpy`` placed on it is
   not the one the child imports.
5. ARC's ``RMG_DB_PATH`` and the child's ``rmgpy`` settings name the same database directory.
6. Every import statement of the scripts T3 runs inside ``rmg_env`` succeeds there.
7. ``python -m arkane --help`` exits 0.

Item 3 matters because ARC sets ``RMG_PATH`` on its own whenever it can find an RMG-Py checkout: an
``RMG-Py`` entry on ``sys.path``, ``$GITHUB_WORKSPACE/RMG-Py``, or ``~/Code/RMG-Py``. A source tree
found that way is put on the child's ``PYTHONPATH`` and shadows the package.
"""

import ast
import json
import os
import re
import subprocess
import sys
import tempfile

import yaml

from arc.imports import settings
from arc.job.env_run import rmg_env_command

DEVTOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
T3_ROOT = os.path.dirname(DEVTOOLS_DIR)
DEFAULT_SPEC_PATH = os.path.join(DEVTOOLS_DIR, 'rmg_env.yml')

# The scripts T3 executes under rmg_env's Python. Their import statements are the set rmg_env must
# satisfy; the t3_env modules that launch them import nothing from rmgpy themselves.
RMG_ENV_SCRIPTS = (
    os.path.join('t3', 'runners', 'rmg_incore_script.py'),
    os.path.join('t3', 'runners', 'pdep_source_seed.py'),
    os.path.join('t3', 'runners', 'statmech_conformer_extract.py'),
)

REPORT_MARKER = '@@T3_RMG_ENV_CHECK@@'
DECOY_VERSION = 'decoy-from-parent-PYTHONPATH'
_REQUIREMENT_RE = re.compile(r'^(?:(?P<channel>[^:\s]+)::)?(?P<name>[A-Za-z0-9_.\-]+)\s*(?P<op>>=|=)\s*(?P<version>[0-9][0-9A-Za-z.]*)$')

# Executed by the rmg_env interpreter, NOT by this one. Standard library only until the rmgpy
# imports, so a broken environment still reports what it is.
_CHILD_PROBE = r'''
import glob
import json
import os
import sys
import sysconfig

report = {
    "python": ".".join(str(part) for part in sys.version_info[:3]),
    "sys_version": sys.version,
    "prefix": sys.prefix,
    "site_packages": sorted({sysconfig.get_paths()["purelib"], sysconfig.get_paths()["platlib"]}),
    "pythonpath_set": bool(os.environ.get("PYTHONPATH")),
    "installed": {},
    "imports": {},
}
for meta in glob.glob(os.path.join(sys.prefix, "conda-meta", "*.json")):
    with open(meta) as f:
        record = json.load(f)
    report["installed"][record["name"]] = record["version"]
try:
    import rmgpy
    import arkane
    from rmgpy import settings as rmg_settings
    report["rmgpy_file"] = rmgpy.__file__
    report["rmgpy_version"] = rmgpy.__version__
    report["arkane_file"] = arkane.__file__
    report["database"] = rmg_settings["database.directory"]
except Exception as e:
    report["error"] = f"{type(e).__name__}: {e}"
for label, code in IMPORTS.items():
    try:
        exec(compile(code, label, "exec"), {})
        report["imports"][label] = "ok"
    except Exception as e:
        report["imports"][label] = f"{type(e).__name__}: {e}"
print(MARKER + json.dumps(report))
'''


def load_spec(spec_path: str) -> dict:
    """
    Read the pins out of an environment spec.

    Args:
        spec_path (str): The path to the conda environment file.

    Returns:
        dict: Keys are package names; values are ``(op, version)`` tuples, ``op`` being ``'='`` or ``'>='``.

    Raises:
        ValueError: If a dependency is not a plain ``[channel::]name=V`` or ``[channel::]name>=V`` pin,
                    or if the spec pins no ``python`` or no ``rmg``.
    """
    with open(spec_path) as f:
        spec = yaml.safe_load(f)
    pins = dict()
    for dependency in spec.get('dependencies') or list():
        match = _REQUIREMENT_RE.match(str(dependency).strip())
        if match is None:
            raise ValueError(f"{spec_path}: cannot check the dependency {dependency!r}; write it as "
                             f"'[channel::]name=V' or '[channel::]name>=V'.")
        pins[match.group('name')] = (match.group('op'), match.group('version'))
    for required in ('python', 'rmg'):
        if pins.get(required, (None,))[0] != '=':
            raise ValueError(f"{spec_path}: the spec must pin '{required}=<version>'.")
    return pins


def version_tuple(version: str) -> tuple:
    """
    Turn a dotted version into a tuple of ints for ordering, ignoring any non-numeric tail of a part.

    Args:
        version (str): A version such as ``'0.9.41'`` or ``'3.11.18'``.

    Returns:
        tuple: E.g. ``(0, 9, 41)``.
    """
    parts = list()
    for part in version.split('.'):
        digits = re.match(r'\d+', part)
        parts.append(int(digits.group()) if digits else 0)
    return tuple(parts)


def pin_holds(op: str, pinned: str, installed: str | None) -> bool:
    """
    Whether an installed version satisfies one pin, with conda's meaning of ``=``.

    Args:
        op (str): ``'='`` (``V`` or ``V.*``) or ``'>='``.
        pinned (str): The version in the spec.
        installed (str): The installed version, or ``None`` if the package is absent.

    Returns:
        bool: Whether the pin holds.
    """
    if installed is None:
        return False
    if op == '=':
        return installed == pinned or installed.startswith(pinned + '.')
    return version_tuple(installed) >= version_tuple(pinned)


def is_inside(path: str | None, directories: list) -> bool:
    """
    Whether ``path`` lies inside any of ``directories``, after resolving symlinks.

    Args:
        path (str): A file path, or ``None``.
        directories (list): Directory paths.

    Returns:
        bool: Whether it does.
    """
    if not path:
        return False
    real = os.path.realpath(path)
    return any(real.startswith(os.path.realpath(d) + os.sep) for d in directories)


def same_directory(first: str | None, second: str | None) -> bool:
    """
    Whether two directory paths name the same directory, tolerating RMG-database's ``input/`` level.

    Args:
        first (str): A path, or ``None``.
        second (str): A path, or ``None``.

    Returns:
        bool: Whether they do.
    """
    if not first or not second:
        return False
    return _database_root(first) == _database_root(second)


def _database_root(path: str) -> str:
    """The resolved database directory, one level up when ``path`` is RMG-database's ``input/``."""
    real = os.path.realpath(path)
    return os.path.dirname(real) if os.path.basename(real) == 'input' else real


def import_snippets(t3_root: str = T3_ROOT, scripts: tuple = RMG_ENV_SCRIPTS) -> dict:
    """
    Collect every import statement of the scripts that run inside rmg_env, nested ones included.

    Args:
        t3_root (str): The T3 repository root.
        scripts (tuple): Script paths relative to ``t3_root``.

    Returns:
        dict: Keys are ``'<script>:<line>'``; values are the import statement as source.
    """
    snippets = dict()
    for script in scripts:
        with open(os.path.join(t3_root, script)) as f:
            tree = ast.parse(f.read(), filename=script)
        for node in ast.walk(tree):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                snippets[f'{script}:{node.lineno}'] = ast.unparse(node)
    return snippets


def evaluate(pins: dict,
             report: dict,
             decoy_report: dict,
             arc_rmg_path: str | None,
             arc_rmg_db_path: str | None,
             arkane_help_rc: int,
             ) -> list:
    """
    Compare what the child reported against the spec and the package-mode rules.

    Args:
        pins (dict): The output of ``load_spec``.
        report (dict): The child's report, launched with the parent's environment as it is.
        decoy_report (dict): The child's report, launched with a decoy ``rmgpy`` on the parent's PYTHONPATH.
        arc_rmg_path (str): ARC's resolved ``RMG_PATH``.
        arc_rmg_db_path (str): ARC's resolved ``RMG_DB_PATH``.
        arkane_help_rc (int): The exit code of ``python -m arkane --help`` in the child.

    Returns:
        list: One message per failed assertion; empty when every assertion holds.
    """
    failures = list()
    installed = report.get('installed', dict())
    for name, (op, pinned) in sorted(pins.items()):
        if not pin_holds(op, pinned, installed.get(name)):
            failures.append(f"pin {name}{op}{pinned} does not hold: installed {name} is {installed.get(name)}")
    python_pin = pins['python'][1]
    if not pin_holds('=', python_pin, report.get('python')):
        failures.append(f"the child's Python is {report.get('python')}, the spec pins {python_pin}")
    if report.get('error'):
        failures.append(f"importing rmgpy/arkane failed in the child: {report['error']}")
    if report.get('rmgpy_version') != pins['rmg'][1]:
        failures.append(f"rmgpy.__version__ is {report.get('rmgpy_version')!r}, the spec pins rmg={pins['rmg'][1]}")
    for module in ('rmgpy', 'arkane'):
        path = report.get(f'{module}_file')
        if not is_inside(path, report.get('site_packages', list())):
            failures.append(f"{module} loads from {path}, not from {report.get('site_packages')} (package mode)")
    if arc_rmg_path:
        failures.append(f"ARC resolved RMG_PATH={arc_rmg_path}, so the child would run that source tree; package "
                        f"mode needs it unset (set RMG_PATH = None in ~/.arc/settings.py)")
    if decoy_report.get('rmgpy_file') is None or decoy_report.get('rmgpy_version') == DECOY_VERSION \
            or decoy_report.get('rmgpy_file') != report.get('rmgpy_file'):
        failures.append(f"a PYTHONPATH set in the parent reached the child: rmgpy came from "
                        f"{decoy_report.get('rmgpy_file')} (version {decoy_report.get('rmgpy_version')!r})")
    if not same_directory(arc_rmg_db_path, report.get('database')):
        failures.append(f"ARC's RMG_DB_PATH is {arc_rmg_db_path} but rmgpy in rmg_env reads {report.get('database')}")
    for label, outcome in sorted(report.get('imports', dict()).items()):
        if outcome != 'ok':
            failures.append(f"import at {label} failed in rmg_env: {outcome}")
    if not report.get('imports'):
        failures.append('no import statement of the rmg_env scripts was checked')
    if arkane_help_rc != 0:
        failures.append(f"python -m arkane --help exited {arkane_help_rc}")
    return failures


def run_in_rmg_env(py_args: list, env: dict) -> subprocess.CompletedProcess:
    """
    Run ``python <py_args>`` in rmg_env through ARC's launcher, as T3 does.

    Args:
        py_args (list): Everything after ``python``.
        env (dict): The parent environment to launch with.

    Returns:
        subprocess.CompletedProcess: The finished child.
    """
    return subprocess.run(['bash', '-c', rmg_env_command(py_args=py_args)],
                          capture_output=True, text=True, env=env, timeout=600)


def probe(probe_path: str, env: dict) -> dict:
    """
    Run the child probe and parse its report.

    Args:
        probe_path (str): The path to the written probe script.
        env (dict): The parent environment to launch with.

    Returns:
        dict: The report, or ``{'error': ...}`` if the child printed none.
    """
    proc = run_in_rmg_env([probe_path], env)
    for line in reversed(proc.stdout.splitlines()):
        if line.startswith(REPORT_MARKER):
            return json.loads(line[len(REPORT_MARKER):])
    return {'error': f'the child exited {proc.returncode} with no report; stderr tail: {proc.stderr[-2000:]}'}


def main(spec_path: str = DEFAULT_SPEC_PATH) -> int:
    """
    Run every check and print the verdict.

    Args:
        spec_path (str): The environment spec to check against.

    Returns:
        int: 0 if every assertion holds, else 1.
    """
    pins = load_spec(spec_path)
    with tempfile.TemporaryDirectory(prefix='t3_rmg_env_check_') as work_dir:
        probe_path = os.path.join(work_dir, 'probe.py')
        with open(probe_path, 'w') as f:
            f.write(f'IMPORTS = {import_snippets()!r}\nMARKER = {REPORT_MARKER!r}\n{_CHILD_PROBE}')
        decoy_dir = os.path.join(work_dir, 'decoy_site')
        os.makedirs(os.path.join(decoy_dir, 'rmgpy'))
        with open(os.path.join(decoy_dir, 'rmgpy', '__init__.py'), 'w') as f:
            f.write(f'__version__ = {DECOY_VERSION!r}\n')

        report = probe(probe_path, dict(os.environ))
        decoy_env = dict(os.environ)
        decoy_env['PYTHONPATH'] = os.pathsep.join([decoy_dir] + [p for p in [os.environ.get('PYTHONPATH')] if p])
        decoy_report = probe(probe_path, decoy_env)
        arkane_help = run_in_rmg_env(['-m', 'arkane', '--help'], dict(os.environ))

    arc_rmg_path, arc_rmg_db_path = settings.get('RMG_PATH'), settings.get('RMG_DB_PATH')
    print(f'spec:                 {spec_path}')
    print(f'pins:                 {", ".join(f"{n}{op}{v}" for n, (op, v) in sorted(pins.items()))}')
    print(f'sys.version:          {report.get("sys_version")}')
    print(f'sys.prefix:           {report.get("prefix")}')
    print(f'rmgpy.__file__:       {report.get("rmgpy_file")}')
    print(f'rmgpy.__version__:    {report.get("rmgpy_version")}')
    print(f'arkane.__file__:      {report.get("arkane_file")}')
    print(f'installed (pinned):   {", ".join(f"{n}={report.get("installed", {}).get(n)}" for n in sorted(pins))}')
    print(f'rmgpy database:       {report.get("database")}')
    print(f'ARC RMG_DB_PATH:      {arc_rmg_db_path}')
    print(f'ARC RMG_PATH:         {arc_rmg_path}')
    print(f'with decoy PYTHONPATH: rmgpy.__file__ = {decoy_report.get("rmgpy_file")}, '
          f'version = {decoy_report.get("rmgpy_version")!r}')
    print(f'imports checked:      {len(report.get("imports", {}))} '
          f'({sum(v == "ok" for v in report.get("imports", {}).values())} ok)')
    print(f'arkane --help:        exit {arkane_help.returncode}')

    failures = evaluate(pins, report, decoy_report, arc_rmg_path, arc_rmg_db_path, arkane_help.returncode)
    if failures:
        print(f'\nrmg_env does NOT match {spec_path}:', file=sys.stderr)
        for failure in failures:
            print(f'  - {failure}', file=sys.stderr)
        if arkane_help.returncode != 0:
            print(f'arkane --help stderr tail:\n{arkane_help.stderr[-2000:]}', file=sys.stderr)
        return 1
    print(f'\nrmg_env matches {spec_path} (package mode).')
    return 0


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:2]))
