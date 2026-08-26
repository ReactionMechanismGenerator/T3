"""
Orchestrate DOF-normalized conformer extraction for the hybrid P-dep network writer.

``t3.pdep.hybrid.write_hybrid_network_input_file`` splices every QM conformer -- transition states
AND adopted wells -- inline and vibration-only (``HarmonicOscillator`` [+ ``HinderedRotor``], plus
the imaginary frequency for a TS), never as Arkane's by-reference full conformer. Turning a raw
ARC/Arkane statmech artifact (``Log(...)`` references to Gaussian ``opt``/``freq`` logs) into that
inline data means loading it through Arkane -- which imports only under the ``rmg_env`` (Python 3.9),
not T3's ``t3_env``. This module runs ``t3.runners.statmech_conformer_extract`` in that sibling env
via the same ``rmg_env_command`` launcher ARC uses everywhere else, and parses the conformer data
back for the writer.

The energy reference is the hybrid network's own header (``model_chemistry``, ``atom_energies``,
``frequency_scale_factor``): applied uniformly to every artifact, it puts each conformer's E0 on one
self-consistent scale, which is all a rate's barriers and detailed balance depend on.
"""
import json
import os
import subprocess
import tempfile

from arc.job.env_run import rmg_env_command

# Framing markers the extraction driver writes around its JSON payload (Arkane logs noise onto the
# same stream, so the payload is recovered between these rather than by parsing the whole stream).
_BEGIN_MARKER = '__CONFORMER_JSON_BEGIN__'
_END_MARKER = '__CONFORMER_JSON_END__'

_DRIVER_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            'runners', 'statmech_conformer_extract.py')


def extract_dof_conformers(transition_states: dict, wells: dict, energy_settings,
                           timeout: float = 900.0) -> tuple[dict, dict]:
    """Extract vibration-only conformer data for every QM transition state and well.

    Args:
        transition_states (dict): Network TS label -> path to its ARC statmech artifact file.
        wells (dict): Network isomer (well) label -> path to its ARC statmech artifact file.
        energy_settings (QMEnergySettings): The hybrid header; supplies the atom energies, frequency
                                            scale factor and model chemistry every artifact is loaded
                                            under, so all conformers land on one energy reference.
        timeout (float): Wall-clock deadline for the extraction subprocess, in seconds.

    Returns:
        tuple[dict, dict]: ``(ts_conformers, well_conformers)``, each label -> conformer-data dict
                           in the shape ``write_hybrid_network_input_file`` consumes.

    Raises:
        RuntimeError: If the extraction subprocess fails, times out, or emits no framed JSON payload.
    """
    # The model chemistry is passed through as ARC's own repr string and reconstructed inside the
    # driver (which has Arkane), keeping this module import-clean under t3_env.
    artifacts = ([{'label': label, 'path': path, 'is_ts': True} for label, path in transition_states.items()]
                 + [{'label': label, 'path': path, 'is_ts': False} for label, path in wells.items()])
    spec = {
        'frequency_scale_factor': energy_settings.frequency_scale_factor,
        'atom_energies': energy_settings.atom_energies,
        'model_chemistry_repr': energy_settings.model_chemistry,
        'artifacts': artifacts,
    }
    payload = _run_driver(spec, timeout=timeout)
    by_label = {c['label']: c for c in payload}
    ts_conformers = {label: by_label[label] for label in transition_states}
    well_conformers = {label: by_label[label] for label in wells}
    return ts_conformers, well_conformers


def _run_driver(spec: dict, timeout: float) -> list:
    fd, spec_path = tempfile.mkstemp(prefix='.dofconf-spec-', suffix='.json')
    try:
        with os.fdopen(fd, 'w') as f:
            json.dump(spec, f)
        inner = rmg_env_command(py_args=[_DRIVER_PATH, spec_path])
        proc = subprocess.run(['bash', '-c', inner], capture_output=True, text=True, timeout=timeout)
        if proc.returncode != 0:
            raise RuntimeError(f"DOF-conformer extraction failed (exit {proc.returncode}).\n"
                               f"stderr tail:\n{proc.stderr[-2000:]}")
        out = proc.stdout
        if _BEGIN_MARKER not in out or _END_MARKER not in out:
            raise RuntimeError("DOF-conformer extraction produced no framed JSON payload.\n"
                               f"stdout tail:\n{out[-2000:]}\nstderr tail:\n{proc.stderr[-2000:]}")
        body = out[out.index(_BEGIN_MARKER) + len(_BEGIN_MARKER):out.index(_END_MARKER)]
        return json.loads(body)
    finally:
        if os.path.isfile(spec_path):
            os.remove(spec_path)
