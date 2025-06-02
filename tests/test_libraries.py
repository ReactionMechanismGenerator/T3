#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_libraries module
"""

import os
import shutil

from rmgpy.data.thermo import ThermoLibrary
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import ThermoData
from rmgpy.statmech import Conformer, IdealGasTranslation, NonlinearRotor, HarmonicOscillator

from t3.common import TEST_DATA_BASE_PATH
from tests.common import run_minimal
import t3.utils.libraries as libraries


def test_add_to_rmg_library():
    """Test adding thermo calculations to an existing thermo library"""
    libraries_path = os.path.join(TEST_DATA_BASE_PATH, 'libraries')
    if not os.path.isdir(libraries_path):
        os.makedirs(libraries_path)

    spc_1 = Species(
        index=1,
        label='C2H4',
        thermo=ThermoData(
            Tdata=([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], 'K'),
            Cpdata=([3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0], 'cal/(mol*K)'),
            H298=(-20.0, 'kcal/mol'),
            S298=(50.0, 'cal/(mol*K)'),
            Tmin=(300.0, 'K'),
            Tmax=(2000.0, 'K'),
        ),
        conformer=Conformer(
            E0=(0.0, 'kJ/mol'),
            modes=[
                IdealGasTranslation(mass=(28.03, 'amu')),
                NonlinearRotor(inertia=([5.6952e-47, 2.7758e-46, 3.3454e-46], 'kg*m^2'), symmetry=1),
                HarmonicOscillator(frequencies=([834.50, 973.31, 975.37, 1067.1, 1238.5, 1379.5, 1472.3, 1691.3,
                                                 3121.6, 3136.7, 3192.5, 3221.0], 'cm^-1')),
            ],
            spin_multiplicity=1,
            optical_isomers=1,
        ),
        smiles='C=C',
    )

    spc_2 = Species(
        index=2,
        label='CH4',
        thermo=ThermoData(
            Tdata=([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], 'K'),
            Cpdata=([3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0], 'cal/(mol*K)'),
            H298=(-50.0, 'kcal/mol'),
            S298=(100.0, 'cal/(mol*K)'),
            Tmin=(300.0, 'K'),
            Tmax=(2000.0, 'K'),
        ),
        conformer=Conformer(
            E0=(0.0, 'kJ/mol'),
            modes=[
                IdealGasTranslation(mass=(28.03, 'amu')),
                NonlinearRotor(inertia=([5.6952e-47, 2.7758e-46, 3.3454e-46], 'kg*m^2'), symmetry=1),
            ],
            spin_multiplicity=1,
            optical_isomers=1,
        ),
        smiles='C',
    )

    spc_3 = Species(
        index=2,
        label='C3H7',
        thermo=ThermoData(
            Tdata=([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0], 'K'),
            Cpdata=([3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0], 'cal/(mol*K)'),
            H298=(-92.0, 'kcal/mol'),  # this is different
            S298=(12.0, 'cal/(mol*K)'),  # this is different
            Tmin=(300.0, 'K'),
            Tmax=(2000.0, 'K'),
        ),
        conformer=Conformer(
            E0=(0.0, 'kJ/mol'),
            modes=[
                IdealGasTranslation(mass=(28.03, 'amu')),
                NonlinearRotor(inertia=([5.6952e-47, 2.7758e-46, 3.3454e-46], 'kg*m^2'), symmetry=1),
            ],
            spin_multiplicity=1,
            optical_isomers=1,
        ),
        smiles='[CH2]CC',
    )

    # 1. Test adding one species to an existing library.
    for lib_name, spc_list in [('RMG_library', [spc_1, spc_2]), ('ARC_library', [spc_3])]:
        thermo_library = ThermoLibrary(name=lib_name, long_desc=lib_name)
        for i, spc in enumerate(spc_list):
            thermo_library.load_entry(index=i,
                                      label=spc.label,
                                      molecule=spc.to_adjacency_list(),
                                      thermo=spc.thermo,
                                      shortDesc=spc.label,
                                      longDesc=spc.label)
        thermo_library.save(os.path.join(libraries_path, f'{lib_name}.py'))

    t3 = run_minimal()
    t3.set_paths()
    t3.paths['ARC thermo lib'] = os.path.join(libraries_path, 'ARC_library.py')
    t3.paths['T3 thermo lib'] = os.path.join(libraries_path, 'RMG_library.py')
    libraries.add_to_rmg_libraries(library_name=t3.t3['options']['library_name'],
                                   shared_library_name=t3.t3['options']['shared_library_name'],
                                   paths=t3.paths,
                                   logger=t3.logger)
    with open(t3.paths['T3 thermo lib'], 'r') as f:
        lines = f.readlines()
    for line in ["        H298 = (-92,'kcal/mol'),\n",
                 "        S298 = (12,'cal/(mol*K)'),\n",
                 ]:
        assert line in lines

    # 2. Test adding one species to an existing library when the new library has a species that also exists.
    for lib_name, spc_list in [('RMG_library', [spc_1, spc_2]), ('ARC_library', [spc_1, spc_3])]:
        thermo_library = ThermoLibrary(name=lib_name, long_desc=lib_name)
        for i, spc in enumerate(spc_list):
            thermo_library.load_entry(index=i,
                                      label=spc.label,
                                      molecule=spc.to_adjacency_list(),
                                      thermo=spc.thermo,
                                      shortDesc=spc.label,
                                      longDesc=spc.label)
        thermo_library.save(os.path.join(libraries_path, f'{lib_name}.py'))

    t3 = run_minimal()
    t3.set_paths()
    t3.paths['ARC thermo lib'] = os.path.join(libraries_path, 'ARC_library.py')
    t3.paths['T3 thermo lib'] = os.path.join(libraries_path, 'RMG_library.py')
    libraries.add_to_rmg_libraries(library_name=t3.t3['options']['library_name'],
                                   shared_library_name=t3.t3['options']['shared_library_name'],
                                   paths=t3.paths,
                                   logger=t3.logger)
    with open(t3.paths['T3 thermo lib'], 'r') as f:
        lines = f.readlines()
    count = 0
    for line in lines:
        if 'entry(' in line:
            count += 1
    assert count == 3


def test_get_rxn_composition():
    """Test getting the reaction composition from a Chemkin file."""
    rxn = Reaction(reactants=[Species(label='C2H4', smiles='C=C'), Species(label='OH', smiles='[OH]')],
                   products=[Species(label='C2H3', smiles='C=[CH]'), Species(label='H2O', smiles='O')])
    pes_composition = libraries.get_rxn_composition(reaction=rxn)
    print(pes_composition)
    assert pes_composition == {
        'reactants': {'C2H4': 1, 'OH': 1},
        'products': {'C2H3': 1, 'H2O': 1}
    }



def teardown_module():
    """teardown any state that was previously set up."""
    path = os.path.join(TEST_DATA_BASE_PATH, 'libraries')
    if os.path.isdir(path):
        shutil.rmtree(path, ignore_errors=True)
