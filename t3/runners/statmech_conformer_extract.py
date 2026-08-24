#!/usr/bin/env python
"""
Extract DOF-normalized (vibration-only) conformer data from ARC/Arkane statmech artifacts.

Runs under the ``rmg_env`` (Python 3.9, Arkane importable), NOT under T3's ``t3_env``: Arkane is
the only thing that turns a raw ``opt``/``freq`` Gaussian log into a scaled-frequency, atom-energy-
corrected ``Conformer`` on a single self-consistent energy reference. This script loads each QM
statmech artifact through Arkane's own ``StatMechJob`` (so the E0 and frequencies are exactly what
Arkane would build for a by-reference splice), then strips the translational and external-rotational
modes and emits ONLY what a DOF-consistent RMG pressure-dependent network species carries:
``E0`` (kJ/mol), the ``HarmonicOscillator`` frequencies (cm^-1), any ``HinderedRotor`` modes, and --
for a transition state -- the single imaginary frequency (cm^-1).

Why vibration-only: an RMG p-dep network expresses every well as ``HarmonicOscillator`` (+ optional
``HinderedRotor``) with translation and overall rotation handled separately by the active-rotor /
rmgmode machinery. A QM transition state spliced in as Arkane's full conformer (translation +
rotation + vibration) sits on a network whose estimated wells are vibration-only, so
``Q_TS / Q_reactant`` never cancels the TS's translational and rotational partition functions -- a
~1e14 factor that makes ``calculate_tst_rate_coefficient`` return physically impossible rates and
the modified-strong-collision matrix go singular. Reducing every QM conformer to the same
vibration-only treatment restores the invariant that no network mixes DOF provenance.

The energy reference is whatever ``atomEnergies``/``frequencyScaleFactor``/``modelChemistry`` are
passed in (the hybrid network's own header): because a rate's forward/reverse barriers and its
detailed balance depend only on E0 *differences*, any single reference applied uniformly to every
QM conformer (wells AND the TS) is self-consistent.

Usage:
    python statmech_conformer_extract.py <spec_json_path>

where the JSON at ``spec_json_path`` is::

    {
      "frequency_scale_factor": 0.988,
      "atom_energies": {"C": ..., "H": ..., "O": ...},
      "model_chemistry": {"method": "wb97xd2023", "basis": "def2tzvp", "software": "gaussian"},
      "artifacts": [
        {"label": "TS2",        "path": "/abs/qm/TS2.py",   "is_ts": true},
        {"label": "O=[C]O(8)",  "path": "/abs/well_8.py",   "is_ts": false}
      ]
    }

The result (one entry per artifact) is written as JSON to stdout, framed between the markers
``__CONFORMER_JSON_BEGIN__`` and ``__CONFORMER_JSON_END__`` so a caller can recover it even when
Arkane logs noise onto the same stream.
"""
import json
import sys

from arkane.modelchem import LevelOfTheory
from arkane.statmech import StatMechJob
from rmgpy.species import Species, TransitionState
from rmgpy.statmech import HarmonicOscillator, HinderedRotor

BEGIN_MARKER = '__CONFORMER_JSON_BEGIN__'
END_MARKER = '__CONFORMER_JSON_END__'

# amu*angstrom^2 per kg*m^2, for RMG's SI inertia values.
_AMU_ANG2_PER_KG_M2 = 6.022140857e46


def _hindered_rotor_dict(mode) -> dict:
    """Serialize a HinderedRotor (or FreeRotor) into the kwargs an inline RMG mode needs."""
    out = {'type': type(mode).__name__}
    if getattr(mode, 'inertia', None) is not None:
        out['inertia'] = float(mode.inertia.value_si * _AMU_ANG2_PER_KG_M2)
        out['inertia_units'] = 'amu*angstrom^2'
    if getattr(mode, 'symmetry', None) is not None:
        out['symmetry'] = int(mode.symmetry)
    if isinstance(mode, HinderedRotor):
        if getattr(mode, 'fourier', None) is not None:
            out['fourier'] = [[float(x) for x in row] for row in mode.fourier.value_si.tolist()]
        elif getattr(mode, 'barrier', None) is not None:
            out['barrier'] = float(mode.barrier.value_si / 1000.0)  # J/mol -> kJ/mol
            out['barrier_units'] = 'kJ/mol'
        out['semiclassical'] = bool(getattr(mode, 'semiclassical', False))
    return out


def extract_one(label: str, path: str, is_ts: bool, scale: float, atom_energies: dict,
                lot: LevelOfTheory) -> dict:
    """Load one artifact through Arkane and return its vibration-only conformer data."""
    sp = TransitionState(label=label) if is_ts else Species(label=label)
    job = StatMechJob(species=sp, path=path)
    job.level_of_theory = lot
    job.frequencyScaleFactor = scale
    job.includeHinderedRotors = True
    job.applyBondEnergyCorrections = False
    job.applyAtomEnergyCorrections = True
    job.atomEnergies = atom_energies
    job.load(pdep=True, plot=False)
    conf = sp.conformer

    frequencies = None
    hindered_rotors = []
    for mode in conf.modes:
        if isinstance(mode, HarmonicOscillator):
            # RMG stores (and returns from value_si) harmonic frequencies already in cm^-1.
            frequencies = [float(f) for f in mode.frequencies.value_si.tolist()]
        elif isinstance(mode, HinderedRotor):
            hindered_rotors.append(_hindered_rotor_dict(mode))
        # IdealGasTranslation, LinearRotor, NonlinearRotor, KRotor, SphericalTopRotor are DROPPED --
        # that is the DOF normalization this script exists to perform.

    if frequencies is None:
        raise ValueError(f"No HarmonicOscillator mode found for '{label}' in '{path}'.")

    result = {
        'label': label,
        'is_ts': is_ts,
        'E0_kJ_mol': float(conf.E0.value_si / 1000.0),
        'frequencies_cm_1': frequencies,
        'spin_multiplicity': int(conf.spin_multiplicity),
        'optical_isomers': int(conf.optical_isomers),
        'hindered_rotors': hindered_rotors,
    }
    if is_ts:
        if sp.frequency is None:
            raise ValueError(f"Transition state '{label}' has no imaginary frequency in '{path}'.")
        result['imaginary_frequency_cm_1'] = float(sp.frequency.value_si)
    return result


def _resolve_level_of_theory(spec) -> LevelOfTheory:
    """Build the LevelOfTheory the E0 atom-energy corrections key on, from either spec form.

    ``model_chemistry`` is an explicit ``{method, basis, software}`` dict; ``model_chemistry_repr``
    is ARC's own ``LevelOfTheory(...)``/``CompositeLevelOfTheory(...)`` repr string (a composite is
    reduced to its energy component -- the level whose corrections set the reference)."""
    if 'model_chemistry' in spec:
        mc = spec['model_chemistry']
        return LevelOfTheory(method=mc['method'], basis=mc.get('basis'), software=mc.get('software'))
    from arkane.modelchem import CompositeLevelOfTheory
    obj = eval(spec['model_chemistry_repr'],  # noqa: S307 -- constrained namespace, ARC's own repr
               {'LevelOfTheory': LevelOfTheory, 'CompositeLevelOfTheory': CompositeLevelOfTheory})
    return obj.energy if isinstance(obj, CompositeLevelOfTheory) else obj


def main() -> int:
    spec = json.load(open(sys.argv[1]))
    scale_raw = spec.get('frequency_scale_factor')
    scale = 1.0 if scale_raw is None else float(scale_raw)
    atom_energies = spec['atom_energies']
    lot = _resolve_level_of_theory(spec)
    out = [extract_one(a['label'], a['path'], bool(a.get('is_ts', False)), scale, atom_energies, lot)
           for a in spec['artifacts']]
    sys.stdout.write(BEGIN_MARKER + '\n')
    sys.stdout.write(json.dumps(out, indent=2) + '\n')
    sys.stdout.write(END_MARKER + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
