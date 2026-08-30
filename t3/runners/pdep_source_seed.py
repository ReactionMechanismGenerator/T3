#!/usr/bin/env python
"""
Render the estimated ``species()`` blocks for a source-only PES exploration seed.

Runs under ``rmg_env`` (Python 3.9, the RMG database importable), NOT under T3's ``t3_env``: every
field a source seed carries -- ``E0``, the ``HarmonicOscillator`` frequencies, spin multiplicity,
optical isomers, molecular weight, ``TransportData``, ``SingleExponentialDown`` and the NASA
``thermo`` -- comes out of RMG's own databases (group-additivity thermo, statmech group
frequencies, transport groups, library thermo), driven for a single species. **None of it is
quantum chemistry** (see I-038): this driver drives the estimators the same way RMG's own model
generation does, and renders each ``species()`` block byte-for-byte the way
``arkane/pdep.py::PressureDependenceJob.save_input_file`` does, so a generated seed is
indistinguishable from a hand-lifted RMG-explored one.

The thermo goes through ``process_thermo_data`` (the exact NASA-fit step RMG's ``thermoengine``
runs), not the raw group-additivity ``ThermoData``: that is what turns the estimate into the NASA
polynomials a seed carries AND sets ``conformer.E0`` from the fitted NASA (a raw ``ThermoData``
gives an E0 ~0.1 kJ/mol off the fitted one -- see I-038's premise check). The bath gas takes
library thermo directly (already NASA) and gets NO statmech: a bath-gas block needs only collision
and thermo, and leaving ``conformer`` unset omits ``E0``/``modes``/``spinMultiplicity``/
``opticalIsomers`` exactly as an RMG-written bath-gas block does.

Usage:
    python pdep_source_seed.py <spec_json_path>

where the JSON at ``spec_json_path`` is::

    {
      "source":    {"identifier": "[O]C=O", "label": "[O]C=O"},
      "bath_gas":  {"identifier": "[Ar]",   "label": "Ar"},
      "thermo_libraries":    ["primaryThermoLibrary"],
      "transport_libraries": ["NOx2018", "GRI-Mech"]
    }

The result is written to stdout as JSON, framed between the markers ``__SEED_JSON_BEGIN__`` and
``__SEED_JSON_END__`` so a caller can recover it even when RMG logs noise onto the same stream::

    {"species_blocks": [<source species() block>, <bath species() block>]}
"""
import json
import sys

from rmgpy import settings as rmg_settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.species import Species
from rmgpy.thermo.thermoengine import process_thermo_data

BEGIN_MARKER = '__SEED_JSON_BEGIN__'
END_MARKER = '__SEED_JSON_END__'


def _load_database(thermo_libraries, transport_libraries):
    """Load only the RMG sub-databases a source seed needs: thermo (group additivity + the given
    libraries), transport, and statmech. No kinetics families are loaded -- the seed carries no
    reactions, so estimating them here would be wasted work."""
    db = RMGDatabase()
    db.load(
        path=rmg_settings['database.directory'],
        thermo_libraries=list(thermo_libraries),
        transport_libraries=list(transport_libraries),
        reaction_libraries=[],
        seed_mechanisms=[],
        kinetics_families='none',
        kinetics_depositories=[],
        statmech_libraries=None,
        depository=False,
    )
    return db


def _estimate_source(db, identifier, label):
    """Drive every estimator a source well needs, in RMG's own order."""
    spc = Species().from_smiles(identifier)
    spc.label = label
    spc.generate_resonance_structures()
    # process_thermo_data: the NASA fit RMG's thermoengine runs on the raw group-additivity data.
    spc.thermo = process_thermo_data(spc, db.thermo.get_thermo_data(spc))
    # generate_statmech reads the fitted thermo for E0 and the statmech DB for the frequencies.
    spc.generate_statmech()
    spc.molecular_weight = (spc.molecule[0].get_molecular_weight() * 1000, 'amu')
    spc.generate_transport_data()
    spc.generate_energy_transfer_model()
    return spc


def _estimate_bath_gas(db, identifier, label):
    """A bath gas needs only collision and thermo. Deliberately NO statmech: leaving ``conformer``
    unset omits E0/modes/spinMultiplicity/opticalIsomers, matching an RMG-written bath-gas block."""
    spc = Species().from_smiles(identifier)
    spc.label = label
    # Same NASA-fit step as the source (process_thermo_data): for a library-hit bath gas like Ar it
    # passes the library polynomials through unchanged but carries the library's E0 into the NASA's
    # own E0 field -- which is exactly what an RMG-written network file's bath-gas block shows. The
    # library label/comment survive it untouched.
    spc.thermo = process_thermo_data(spc, db.thermo.get_thermo_data(spc))
    spc.molecular_weight = (spc.molecule[0].get_molecular_weight() * 1000, 'amu')
    spc.generate_transport_data()
    spc.generate_energy_transfer_model()
    # process_thermo_data populates a Conformer (E0/spin/optical) as a side effect. A bath-gas block
    # carries none of that (it is a collider, not a well); dropping it omits those species() fields,
    # matching an RMG-written bath-gas block, while the thermo (now carrying its E0) is untouched.
    spc.conformer = None
    return spc


def _render_species_block(spec):
    """Render a ``species()`` block byte-for-byte as ``arkane/pdep.py::save_input_file`` does
    (RMG-Py arkane/pdep.py lines 674-697): same field order, same ``{!r}`` / ``{:d}`` formatting,
    so the RMG statmech/thermo/transport reprs (which ARE the Arkane DSL) render identically."""
    lines = ['species(']
    lines.append('    label = {0!r},'.format(str(spec)))
    if len(spec.molecule) > 0:
        lines.append('    structure = adjacencyList("""{0}"""),'.format(spec.molecule[0].to_adjacency_list()))
    if spec.conformer is not None:
        if spec.conformer.E0 is not None:
            lines.append('    E0 = {0!r},'.format(spec.conformer.E0))
        if len(spec.conformer.modes) > 0:
            lines.append('    modes = [')
            for mode in spec.conformer.modes:
                lines.append('        {0!r},'.format(mode))
            lines.append('    ],')
        lines.append('    spinMultiplicity = {0:d},'.format(spec.conformer.spin_multiplicity))
        lines.append('    opticalIsomers = {0:d},'.format(spec.conformer.optical_isomers))
    if spec.molecular_weight is not None:
        lines.append('    molecularWeight = {0!r},'.format(spec.molecular_weight))
    if spec.transport_data is not None:
        lines.append('    collisionModel = {0!r},'.format(spec.transport_data))
    if spec.energy_transfer_model is not None:
        lines.append('    energyTransferModel = {0!r},'.format(spec.energy_transfer_model))
    if spec.thermo is not None:
        lines.append('    thermo = {0!r},'.format(spec.thermo))
    lines.append(')')
    return '\n'.join(lines)


def main():
    with open(sys.argv[1]) as f:
        spec = json.load(f)

    db = _load_database(spec['thermo_libraries'], spec['transport_libraries'])

    source = _estimate_source(db, spec['source']['identifier'], spec['source']['label'])
    bath = _estimate_bath_gas(db, spec['bath_gas']['identifier'], spec['bath_gas']['label'])

    out = {'species_blocks': [_render_species_block(source), _render_species_block(bath)]}

    sys.stdout.write(BEGIN_MARKER + '\n')
    sys.stdout.write(json.dumps(out) + '\n')
    sys.stdout.write(END_MARKER + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
