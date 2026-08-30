"""
Generate a source-only PES exploration seed from a species identifier.

T3's PES loop is autonomous once seeded, but a new surface has until now required a hand-assembled
Arkane source file (one ``species()`` for the source well carrying estimated statmech/thermo/
transport, the bath gas, one ``network()`` naming the well as the sole isomer, and one
``pressureDependence()`` block). This module builds that file from a single species identifier so a
new surface costs a SMILES: it drives RMG's own estimators for the source well and bath gas -- in a
sibling ``rmg_env`` subprocess, since the RMG database imports only there, not under ``t3_env`` --
and wraps the two rendered ``species()`` blocks in the ``network()``/``pressureDependence()``
boilerplate the explorer needs. The produced file is what
``t3.pdep.explorer.input_file.write_arkane_explorer_input_file`` consumes.

**The source must be unimolecular.** A bimolecular source (e.g. ``H + CO2``, or a dotted SMILES
like ``[H].O=C=O``) makes Arkane's explorer discover two disjoint pressure-dependent networks,
which ``ArkaneExplorerAdapter`` refuses outright -- there is no single unambiguous result to fold
back, and the failure surfaces deep inside Arkane where the cause is unrecognisable. This module
refuses such an identifier HERE, at the entry point, before spawning the estimator subprocess, with
a message that names the disjoint-network reason (see :class:`BimolecularSourceError`).

This module imports no ``rmgpy``/``arkane``: the estimation runs behind ``rmg_env_command`` (the
same launcher ARC uses everywhere), exactly as ``t3.pdep.dof_conformers`` does. The only RMG-derived
import is ``arc.molecule.molecule.Molecule`` -- already a top-level import in
``t3.pdep.explorer.input_file`` -- used only to parse the identifier and count its disconnected
fragments for the unimolecular check.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

from arc.job.env_run import rmg_env_command
from arc.molecule.molecule import Molecule

# Framing markers the estimator driver writes around its JSON payload (RMG logs noise onto the same
# stream, so the payload is recovered between these rather than by parsing the whole stream).
_BEGIN_MARKER = '__SEED_JSON_BEGIN__'
_END_MARKER = '__SEED_JSON_END__'

_DRIVER_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                            'runners', 'pdep_source_seed.py')

# RMG databases the estimator loads. Kept small and matched to the worked CHO2 seed
# (``primaryThermoLibrary`` for library thermo of common species/bath gases, NOx2018 + GRI-Mech for
# transport groups). Overridable per call.
_DEFAULT_THERMO_LIBRARIES = ('primaryThermoLibrary',)
_DEFAULT_TRANSPORT_LIBRARIES = ('NOx2018', 'GRI-Mech')

# SMILES for the bath gases a caller can name by symbol; an identifier not in this map is passed to
# the estimator as a SMILES verbatim, so an exotic bath gas still works.
_BATH_GAS_SMILES = {'Ar': '[Ar]', 'He': '[He]', 'Ne': '[Ne]', 'Kr': '[Kr]', 'N2': 'N#N'}

# The loop's standard pressure-dependence configuration, identical to the worked CHO2 seed. These
# are solver settings (the T/P grids, grain size, method, interpolation and active-rotor flags),
# NOT per-species estimates, so they are the same for any system and are templated verbatim. A
# caller may override any key via ``pdep_settings=``. The comma style (no space after the comma
# inside a quantity tuple) mirrors RMG's own ``Quantity.__repr__`` so a generated file diffs cleanly
# against an RMG-written one; ``_render_pdep_block`` reproduces it.
_DEFAULT_PDEP_SETTINGS = {
    'Tmin': (700, 'K'),
    'Tmax': (3200, 'K'),
    'Tcount': 10,
    'Tlist': ([3200, 2290.91, 1784.07, 1460.87, 1236.81, 1072.34, 946.479, 847.059, 766.54, 700], 'K'),
    'Pmin': (0.1, 'bar'),
    'Pmax': (100, 'bar'),
    'Pcount': 10,
    'Plist': ([0.1, 0.215443, 0.464159, 1, 2.15443, 4.64159, 10, 21.5443, 46.4159, 100], 'bar'),
    'maximumGrainSize': (2, 'kJ/mol'),
    'minimumGrainCount': 250,
    'method': 'modified strong collision',
    'interpolationModel': ('pdeparrhenius',),
    'activeKRotor': True,
    'activeJRotor': True,
    'rmgmode': True,
}


class BimolecularSourceError(ValueError):
    """A source identifier that is not a single connected molecule.

    A bimolecular (or larger) source makes Arkane's explorer discover two disjoint pressure-
    dependent networks, which ``ArkaneExplorerAdapter`` refuses outright -- so it is refused here,
    at the entry point, rather than deep inside Arkane where the cause is unrecognisable.
    """


def generate_source_seed(identifier: str,
                         dest_path: str,
                         *,
                         source_label: str = None,
                         network_name: str = None,
                         bath_gas: str = 'Ar',
                         pdep_settings: dict = None,
                         thermo_libraries=None,
                         transport_libraries=None,
                         timeout: float = 1800.0) -> str:
    """Generate a source-only PES exploration seed from a species identifier.

    Args:
        identifier (str): The source well's SMILES. Must describe a single connected molecule; a
                          bimolecular identifier is refused (see :class:`BimolecularSourceError`).
        dest_path (str): Where to write the seed file. The parent directory is created if needed.
        source_label (str, optional): The ``species()``/``network`` isomer label. Defaults to
                                      ``identifier``.
        network_name (str, optional): The human name spliced into the network label
                                      ``'PDepNetwork <network_name> source'``. Defaults to
                                      ``source_label``.
        bath_gas (str, optional): The bath gas, by symbol (``'Ar'``, ``'He'``, ``'Ne'``, ``'Kr'``,
                                  ``'N2'``) or as a SMILES. Defaults to ``'Ar'``.
        pdep_settings (dict, optional): Overrides merged onto ``_DEFAULT_PDEP_SETTINGS`` (the loop's
                                        standard solver configuration).
        thermo_libraries (sequence, optional): RMG thermo libraries to load. Defaults to
                                               ``_DEFAULT_THERMO_LIBRARIES``.
        transport_libraries (sequence, optional): RMG transport libraries to load. Defaults to
                                                  ``_DEFAULT_TRANSPORT_LIBRARIES``.
        timeout (float): Wall-clock deadline for the estimator subprocess, in seconds.

    Returns:
        str: ``dest_path``.

    Raises:
        BimolecularSourceError: If ``identifier`` is not a single connected molecule.
        ValueError: If ``identifier`` cannot be parsed as a SMILES.
        RuntimeError: If the estimator subprocess fails, times out, emits no framed payload, emits a
                      malformed payload, or returns other than exactly two ``species()`` blocks.
    """
    _assert_unimolecular(identifier)

    source_label = source_label or identifier
    network_name = network_name or source_label
    bath_identifier = _BATH_GAS_SMILES.get(bath_gas, bath_gas)

    spec = {
        'source': {'identifier': identifier, 'label': source_label},
        'bath_gas': {'identifier': bath_identifier, 'label': bath_gas},
        'thermo_libraries': list(thermo_libraries or _DEFAULT_THERMO_LIBRARIES),
        'transport_libraries': list(transport_libraries or _DEFAULT_TRANSPORT_LIBRARIES),
    }
    source_block, bath_block = _run_driver(spec, timeout=timeout)

    text = _assemble_seed_file(source_block=source_block,
                               bath_block=bath_block,
                               source_label=source_label,
                               bath_label=bath_gas,
                               network_name=network_name,
                               pdep_settings=pdep_settings)

    parent = os.path.dirname(os.path.abspath(dest_path))
    os.makedirs(parent, exist_ok=True)
    with open(dest_path, 'w') as f:
        f.write(text)
    return dest_path


def _assert_unimolecular(identifier: str) -> Molecule:
    """Parse ``identifier`` and refuse it unless it is a single connected molecule."""
    try:
        mol = Molecule(smiles=identifier)
    except Exception as e:
        raise ValueError(f"Could not parse source identifier {identifier!r} as a SMILES: {e}") from e
    fragments = mol.split()
    if len(fragments) != 1:
        raise BimolecularSourceError(
            f"Source identifier {identifier!r} is not unimolecular: it describes {len(fragments)} "
            f"disconnected molecules. A bimolecular (or larger) source makes Arkane's explorer "
            f"discover two disjoint pressure-dependent networks, which ArkaneExplorerAdapter refuses "
            f"outright -- there is no single unambiguous result to fold back. Seed the PES loop from "
            f"one connected well instead.")
    return mol


def _run_driver(spec: dict, timeout: float) -> list:
    """Run the estimator driver under ``rmg_env`` and return its two ``species()`` blocks.

    Every failure mode is re-raised as ``RuntimeError`` carrying stdout/stderr context, honouring
    ``generate_source_seed``'s documented contract: a bare ``subprocess.TimeoutExpired``, a malformed
    payload, or a wrong-shaped result would otherwise escape as an opaque ``TimeoutExpired`` /
    ``JSONDecodeError`` / ``KeyError`` / tuple-unpack error the caller cannot diagnose."""
    fd, spec_path = tempfile.mkstemp(prefix='.seedgen-spec-', suffix='.json')
    try:
        with os.fdopen(fd, 'w') as f:
            json.dump(spec, f)
        inner = rmg_env_command(py_args=[_DRIVER_PATH, spec_path])
        try:
            proc = subprocess.run(['bash', '-c', inner], capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired as e:
            raise RuntimeError(
                f"Source-seed estimation timed out after {timeout:g} s.\n"
                f"stdout tail:\n{_as_text(e.stdout)[-2000:]}\n"
                f"stderr tail:\n{_as_text(e.stderr)[-2000:]}") from e
        if proc.returncode != 0:
            raise RuntimeError(f"Source-seed estimation failed (exit {proc.returncode}).\n"
                               f"stderr tail:\n{proc.stderr[-2000:]}")
        out = proc.stdout
        if _BEGIN_MARKER not in out or _END_MARKER not in out:
            raise RuntimeError("Source-seed estimation produced no framed JSON payload.\n"
                               f"stdout tail:\n{out[-2000:]}\nstderr tail:\n{proc.stderr[-2000:]}")
        body = out[out.index(_BEGIN_MARKER) + len(_BEGIN_MARKER):out.index(_END_MARKER)]
        try:
            payload = json.loads(body)
        except (ValueError, TypeError) as e:  # json.JSONDecodeError is a ValueError subclass
            raise RuntimeError(
                f"Source-seed estimation emitted a framed payload that is not valid JSON: {e}.\n"
                f"payload tail:\n{body[-2000:]}\nstderr tail:\n{proc.stderr[-2000:]}") from e
        if not isinstance(payload, dict) or 'species_blocks' not in payload:
            keys = list(payload) if isinstance(payload, dict) else type(payload).__name__
            raise RuntimeError(
                f"Source-seed estimation payload has no 'species_blocks' key (got: {keys}).\n"
                f"payload tail:\n{body[-2000:]}\nstderr tail:\n{proc.stderr[-2000:]}")
        blocks = payload['species_blocks']
        if not isinstance(blocks, list) or len(blocks) != 2:
            found = len(blocks) if isinstance(blocks, (list, tuple)) else f'a non-list ({type(blocks).__name__})'
            raise RuntimeError(
                f"Source-seed estimation returned {found} species block(s); exactly two are required "
                f"(the source well and the bath gas).\n"
                f"payload tail:\n{body[-2000:]}\nstderr tail:\n{proc.stderr[-2000:]}")
        return blocks
    finally:
        if os.path.isfile(spec_path):
            os.remove(spec_path)


def _as_text(stream) -> str:
    """Coerce a captured subprocess stream (``str`` / ``bytes`` / ``None``) to text for a message."""
    if stream is None:
        return ''
    if isinstance(stream, bytes):
        return stream.decode(errors='replace')
    return stream


def _assemble_seed_file(source_block: str, bath_block: str, source_label: str, bath_label: str,
                        network_name: str, pdep_settings: dict = None) -> str:
    """Assemble the full seed file: header, the two estimated ``species()`` blocks, and the
    templated ``network()`` and ``pressureDependence()`` blocks."""
    settings = dict(_DEFAULT_PDEP_SETTINGS)
    if pdep_settings:
        settings.update(pdep_settings)
    network_label = f'PDepNetwork {network_name} source'
    parts = [
        _header(source_label, bath_label),
        source_block,
        bath_block,
        _render_network_block(network_label, source_label, {bath_label: 1.0}),
        _render_pdep_block(network_label, settings),
    ]
    return '\n\n'.join(parts) + '\n'


def _header(source_label: str, bath_label: str) -> str:
    return (
        f"# Source-only PES exploration seed for {source_label!r}, generated by "
        f"t3.pdep.explorer.seed (I-038).\n"
        "#\n"
        f"# Carries ONLY the source well {source_label!r}, the {bath_label!r} bath gas, one "
        "network() block naming\n"
        "# the well as the sole isomer, and one pressureDependence() block -- the minimum\n"
        "# t3.pdep.explorer.input_file.write_arkane_explorer_input_file needs. There is NO "
        "reaction() and NO\n"
        "# transitionState(): Arkane's explorer is what creates them from the RMG database, so "
        "supplying them\n"
        "# would invert the data flow.\n"
        "#\n"
        "# Every statmech/thermo/transport field below comes from RMG's own estimators (group "
        "additivity,\n"
        "# statmech group frequencies, transport groups, library thermo) driven for one species -- "
        "no quantum\n"
        "# chemistry, no invented number.\n"
        "#\n"
        "# The source is UNIMOLECULAR by construction: a bimolecular source makes the explorer "
        "discover two\n"
        "# disjoint networks, which ArkaneExplorerAdapter refuses."
    )


def _render_network_block(network_label: str, isomer_label: str, bath_gas: dict) -> str:
    lines = ['network(', f'    label = {network_label!r},', '    isomers = [',
             f'        {isomer_label!r},', '    ],', '    reactants = [', '    ],', '    bathGas = {']
    for label, frac in bath_gas.items():
        lines.append(f'        {label!r}: {float(frac)!r},')
    lines += ['    },', ')']
    return '\n'.join(lines)


def _render_pdep_block(network_label: str, settings: dict) -> str:
    """Render the ``pressureDependence()`` block, reproducing RMG's ``Quantity.__repr__`` comma
    style (no space after the comma inside a quantity tuple) so the block diffs cleanly against an
    RMG-written one."""
    lines = ['pressureDependence(', f'    label = {network_label!r},']
    lines.append(f'    Tmin = {_q(settings["Tmin"])},')
    lines.append(f'    Tmax = {_q(settings["Tmax"])},')
    lines.append(f'    Tcount = {settings["Tcount"]:d},')
    lines.append(f'    Tlist = {_ql(settings["Tlist"])},')
    lines.append(f'    Pmin = {_q(settings["Pmin"])},')
    lines.append(f'    Pmax = {_q(settings["Pmax"])},')
    lines.append(f'    Pcount = {settings["Pcount"]:d},')
    lines.append(f'    Plist = {_ql(settings["Plist"])},')
    lines.append(f'    maximumGrainSize = {_q(settings["maximumGrainSize"])},')
    lines.append(f'    minimumGrainCount = {settings["minimumGrainCount"]:d},')
    lines.append(f'    method = {settings["method"]!r},')
    lines.append(f'    interpolationModel = {settings["interpolationModel"]!r},')
    lines.append(f'    activeKRotor = {settings["activeKRotor"]!r},')
    lines.append(f'    activeJRotor = {settings["activeJRotor"]!r},')
    if settings.get('rmgmode'):
        lines.append(f'    rmgmode = {settings["rmgmode"]!r},')
    lines.append(')')
    return '\n'.join(lines)


def _q(value_unit) -> str:
    """Render a ``(value, 'unit')`` quantity as ``(value,'unit')`` -- RMG's own comma style."""
    value, unit = value_unit
    return f'({value!r},{unit!r})'


def _ql(list_unit) -> str:
    """Render a ``([v, ...], 'unit')`` list quantity as ``([v,...],'unit')`` -- RMG's comma style."""
    values, unit = list_unit
    inner = ','.join(repr(v) for v in values)
    return f'([{inner}],{unit!r})'


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate a source-only PES exploration seed from a species identifier (SMILES).")
    parser.add_argument('identifier', help="Source well SMILES (must be a single connected molecule).")
    parser.add_argument('dest_path', help="Where to write the seed file.")
    parser.add_argument('--source-label', default=None, help="species()/network isomer label.")
    parser.add_argument('--network-name', default=None, help="Name spliced into the network label.")
    parser.add_argument('--bath-gas', default='Ar', help="Bath gas symbol or SMILES (default: Ar).")
    parser.add_argument('--timeout', type=float, default=1800.0, help="Estimator subprocess timeout (s).")
    args = parser.parse_args(argv)
    try:
        path = generate_source_seed(args.identifier, args.dest_path,
                                    source_label=args.source_label,
                                    network_name=args.network_name,
                                    bath_gas=args.bath_gas,
                                    timeout=args.timeout)
    except BimolecularSourceError as e:
        sys.stderr.write(f"REFUSED: {e}\n")
        return 2
    sys.stdout.write(f"Wrote source-only PES seed to {path}\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
