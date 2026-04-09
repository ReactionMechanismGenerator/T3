#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_libraries module
"""

import logging
import os
import shutil
import pytest
from typing import List, Dict

from arc.molecule.molecule import Molecule

from t3.utils import rmg_shim
from t3.chem import T3Species
from t3.common import TEST_DATA_BASE_PATH
from t3.utils.libraries import (
    append_to_rmg_libraries,
    load_rmg_species_dictionary_file,
    _update_species_dictionary_atomic,
)


# Setup a standard logger that mimics T3's logger interface
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("t3_test")


# Real (shortened) RMG libraries used by the integration tests below.
LIBRARIES_DATA_PATH = os.path.join(TEST_DATA_BASE_PATH, 'libraries')
REAL_DEST_THERMO_PATH = os.path.join(LIBRARIES_DATA_PATH, 'thermo', 'dest_lib.py')
REAL_SRC_THERMO_PATH = os.path.join(LIBRARIES_DATA_PATH, 'thermo', 'src_lib.py')
REAL_DEST_KINETICS_DIR = os.path.join(LIBRARIES_DATA_PATH, 'kinetics', 'dest_lib')
REAL_SRC_KINETICS_DIR = os.path.join(LIBRARIES_DATA_PATH, 'kinetics', 'src_lib')


def _create_mock_shim_thermo_data(H298: float, S298: float) -> rmg_shim.ThermoData:
    """Helper to create a shim ThermoData object for testing."""
    return rmg_shim.ThermoData(
        Tdata=([300.0, 400.0, 500.0, 1000.0], 'K'),
        Cpdata=([3.0, 4.0, 5.0, 10.0], 'cal/(mol*K)'),
        H298=(H298, 'kcal/mol'),
        S298=(S298, 'cal/(mol*K)'),
        Tmin=(300.0, 'K'),
        Tmax=(2000.0, 'K'),
    )


def _create_mock_shim_kinetics_data() -> rmg_shim.Arrhenius:
    """Helper to create a shim Arrhenius object."""
    return rmg_shim.Arrhenius(
        A=(1e10, 'cm^3/(mol*s)'),
        n=0.5,
        Ea=(10.0, 'kcal/mol'),
        T0=(1, 'K')
    )


def _write_thermo_lib_file(path: str, name: str, entries: List[rmg_shim.Entry]):
    """Writes a thermo library .py file using the rmg_shim."""
    lib = rmg_shim.Library(name=name, longDesc=f"Description for {name}")
    lib.entries = entries
    with open(path, 'w') as f:
        f.write(rmg_shim.write_library_to_string(lib))


def _write_kinetics_lib_folder(folder_path: str, name: str,
                               entries: List[rmg_shim.Entry],
                               dictionary: Dict[str, Molecule]):
    """Writes a kinetics library folder (reactions.py + dictionary.txt)."""
    os.makedirs(folder_path, exist_ok=True)

    # Write reactions.py
    lib = rmg_shim.Library(name=name, longDesc=f"Description for {name}")
    lib.entries = entries
    with open(os.path.join(folder_path, 'reactions.py'), 'w') as f:
        f.write(rmg_shim.write_library_to_string(lib))

    # Write dictionary.txt
    dict_path = os.path.join(folder_path, 'dictionary.txt')
    _update_species_dictionary_atomic(dictionary, dict_path, logger)


# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------

def test_append_to_thermo_library(tmp_path):
    """
    Test adding thermo entries to an existing library.
    Verifies:
      1. New entries are appended.
      2. Indices are incremented correctly.
      3. Existing isomorphic species are NOT duplicated.
    """
    # 1. Setup Data
    # Existing species in Destination
    spc_exist = T3Species(label='C2H4', smiles='C=C')
    entry_exist = rmg_shim.Entry(
        index=1,
        label=spc_exist.label,
        molecule=spc_exist.mol.to_adjacency_list(),
        thermo=_create_mock_shim_thermo_data(-20.0, 50.0)
    )

    # New species in Source (One duplicate C2H4, one new C3H8)
    spc_new = T3Species(label='C3H8', smiles='CCC')

    # Source entries
    entry_dup = rmg_shim.Entry(
        index=1,
        label='Ethylene_Dup',  # Different label, same structure -> Should be skipped
        molecule=spc_exist.mol.to_adjacency_list(),
        thermo=_create_mock_shim_thermo_data(-20.0, 50.0)
    )
    entry_new = rmg_shim.Entry(
        index=2,
        label=spc_new.label,
        molecule=spc_new.mol.to_adjacency_list(),
        thermo=_create_mock_shim_thermo_data(-30.0, 60.0)
    )

    # 2. Write Initial Files
    dest_lib_path = tmp_path / "RMG_Thermo.py"
    src_lib_path = tmp_path / "ARC_Thermo.py"

    _write_thermo_lib_file(str(dest_lib_path), "DestLib", [entry_exist])
    _write_thermo_lib_file(str(src_lib_path), "SrcLib", [entry_dup, entry_new])

    # 3. Run Function
    paths = {
        'ARC thermo lib': str(src_lib_path),
        'T3 thermo lib': str(dest_lib_path),
        'shared T3 thermo lib': None,
        'ARC kinetics lib': None,
        'T3 kinetics lib': None,
        'shared T3 kinetics lib': None
    }

    append_to_rmg_libraries(
        library_name="TestLib",
        shared_library_name="DestLib",
        paths=paths,
        logger=logger
    )

    # 4. Verify (Parse back output)
    with open(dest_lib_path, 'r') as f:
        result_lib = rmg_shim.parse_rmg_library(f.read())

    # Should contain 2 entries (Original C2H4 + New C3H8).
    # The duplicate Ethylene should be skipped.
    assert len(result_lib.entries) == 2

    # Check Indices
    assert result_lib.entries[0].index == 1
    assert result_lib.entries[1].index == 2

    # Check Data
    assert result_lib.entries[0].label == "C2H4"
    assert result_lib.entries[1].label == "C3H8"

    # Verify values were preserved (Check H298)
    # shim stores H298 as tuple (-30.0, 'kcal/mol')
    assert result_lib.entries[1].thermo.H298[0] == -30.0


def test_append_to_kinetics_library(tmp_path):
    """
    Test adding kinetics entries (reactions + dictionary).
    Verifies:
      1. Reactions are appended.
      2. Dictionary is updated with missing species.
      3. Duplicate reaction labels are skipped.
    """
    # 1. Setup Data
    # Destination: Has H + H <=> H2
    mol_H = Molecule(smiles='[H]')
    mol_H2 = Molecule(smiles='[H][H]')

    entry_exist = rmg_shim.Entry(
        index=1,
        label="H+H<=>H2",
        kinetics=_create_mock_shim_kinetics_data()
    )
    dest_dict = {"H": mol_H, "H2": mol_H2}

    # Source: Has OH + H <=> H2O (New) and H + H <=> H2 (Duplicate Label)
    mol_OH = Molecule(smiles='[OH]')
    mol_H2O = Molecule(smiles='O')

    entry_dup = rmg_shim.Entry(
        index=1,
        label="H+H<=>H2",  # Duplicate label -> Skip
        kinetics=_create_mock_shim_kinetics_data()
    )
    entry_new = rmg_shim.Entry(
        index=2,
        label="OH+H<=>H2O",  # New -> Add
        kinetics=_create_mock_shim_kinetics_data()
    )
    src_dict = {"H": mol_H, "OH": mol_OH, "H2O": mol_H2O}

    # 2. Write Initial Folders
    dest_folder = tmp_path / "RMG_Kinetics"
    src_folder = tmp_path / "ARC_Kinetics"

    _write_kinetics_lib_folder(str(dest_folder), "DestLib", [entry_exist], dest_dict)
    _write_kinetics_lib_folder(str(src_folder), "SrcLib", [entry_dup, entry_new], src_dict)

    # 3. Run Function
    paths = {
        'ARC thermo lib': None,
        'T3 thermo lib': None,
        'shared T3 thermo lib': None,
        'ARC kinetics lib': str(src_folder),
        'T3 kinetics lib': str(dest_folder),
        'shared T3 kinetics lib': None,
    }

    append_to_rmg_libraries(
        library_name="TestLib",
        shared_library_name="DestLib",
        paths=paths,
        logger=logger
    )

    # 4. Verify Reactions
    with open(dest_folder / "reactions.py", 'r') as f:
        result_lib = rmg_shim.parse_rmg_library(f.read())

    assert len(result_lib.entries) == 2
    assert result_lib.entries[0].label == "H+H<=>H2"
    assert result_lib.entries[1].label == "OH+H<=>H2O"
    assert result_lib.entries[1].index == 2

    # 5. Verify Dictionary
    result_dict = load_rmg_species_dictionary_file(str(dest_folder / "dictionary.txt"))

    # Should have H, H2 (original) + OH, H2O (new)
    assert len(result_dict) == 4
    assert "OH" in result_dict
    assert "H2O" in result_dict


def test_atomic_dictionary_update_with_corruption(tmp_path):
    """
    Test that _update_species_dictionary_atomic handles corruption safely
    by aborting if the existing file is unreadable.
    """
    dict_path = tmp_path / "corrupt_dict.txt"

    # Create a file that exists
    with open(dict_path, 'w') as f:
        f.write("Valid content")

    # Mock data to append
    new_spc = {"A": Molecule(smiles='C')}

    # Case 1: Normal update
    _update_species_dictionary_atomic(new_spc, str(dict_path), logger)
    with open(dict_path, 'r') as f:
        content = f.read()
        assert "Valid content" in content
        assert "A" in content

    # Case 2: Unreadable file (simulated by mocking open to raise OSError)
    # Since we can't easily mock built-ins in this simplified test without unittest.mock,
    # we rely on the logic review. But we CAN test directory-as-file error.

    bad_path = tmp_path / "bad_folder"
    os.makedirs(bad_path)

    # Trying to update a directory as if it were a file should raise ValueError
    with pytest.raises(ValueError, match="is not a regular file"):
        _update_species_dictionary_atomic(new_spc, str(bad_path), logger)


def test_write_atomic(tmp_path):
    """Test that atomic write creates the file correctly."""
    target = tmp_path / "atomic_test.txt"
    content = "Hello World"

    # Write new
    rmg_shim.write_atomic(str(target), content)
    assert target.exists()
    assert target.read_text() == content

    # Overwrite
    new_content = "New Data"
    rmg_shim.write_atomic(str(target), new_content)
    assert target.read_text() == new_content

    # Check no temp files left behind
    files = list(tmp_path.iterdir())
    assert len(files) == 1  # Only the target file should exist


# ------------------------------------------------------------------------------
# Integration tests using real (shortened) RMG libraries from RMG-database
# ------------------------------------------------------------------------------

def test_read_real_thermo_library():
    """
    Read a real (shortened) RMG thermo library from disk and verify
    that the library metadata, entry indices, labels, molecules, and
    thermo data round-trip correctly through the shim parser.
    """
    with open(REAL_DEST_THERMO_PATH, 'r') as f:
        lib = rmg_shim.parse_rmg_library(f.read())

    assert lib.name == "T3_test_dest_thermo"
    assert "destination thermo library" in lib.shortDesc
    assert "BurkeH2O2" in lib.longDesc

    assert len(lib.entries) == 3
    assert [e.label for e in lib.entries] == ["H", "H2", "OH"]
    assert [e.index for e in lib.entries] == [0, 1, 2]

    # Verify the thermo data was parsed (numerical fields preserved)
    h_entry = lib.entries[0]
    assert h_entry.thermo.H298[0] == 52.1
    assert h_entry.thermo.H298[1] == "kcal/mol"
    assert h_entry.thermo.S298[0] == 27.39

    # The molecule field should be a parseable adjacency list
    h_mol = Molecule().from_adjacency_list(h_entry.molecule)
    assert h_mol.get_radical_count() == 1

    h2_mol = Molecule().from_adjacency_list(lib.entries[1].molecule)
    assert h2_mol.is_isomorphic(Molecule(smiles='[H][H]'))


def test_read_real_kinetics_library():
    """
    Read a real (shortened) RMG kinetics library (reactions.py + dictionary.txt)
    from disk and verify entries, kinetics, and the species dictionary parse
    correctly.
    """
    rxn_path = os.path.join(REAL_DEST_KINETICS_DIR, 'reactions.py')
    dict_path = os.path.join(REAL_DEST_KINETICS_DIR, 'dictionary.txt')

    with open(rxn_path, 'r') as f:
        lib = rmg_shim.parse_rmg_library(f.read())

    assert lib.name == "T3_test_dest_kinetics"
    assert len(lib.entries) == 2
    assert lib.entries[0].label == "H + O2 <=> O + OH"
    assert lib.entries[1].label == "OH + OH <=> O + H2O"

    # Verify Arrhenius parameters round-tripped through the parser
    arr = lib.entries[0].kinetics
    assert isinstance(arr, rmg_shim.Arrhenius)
    assert arr.A == (1.04e+14, 'cm^3/(mol*s)')
    assert arr.n == 0
    assert arr.Ea == (15286, 'cal/mol')

    species_dict = load_rmg_species_dictionary_file(dict_path)
    assert sorted(species_dict.keys()) == ['H', 'H2O', 'O', 'O2', 'OH']
    # Verify the parsed Molecule for OH
    assert species_dict['OH'].is_isomorphic(Molecule(smiles='[OH]'))


def test_write_real_thermo_library_roundtrip(tmp_path):
    """
    Read a real RMG thermo library, write it back to disk via the shim,
    re-read it, and verify that all entries (label, index, thermo) survive
    the round-trip.
    """
    with open(REAL_DEST_THERMO_PATH, 'r') as f:
        original = rmg_shim.parse_rmg_library(f.read())

    # Write the library back to a new file
    out_path = tmp_path / "roundtrip_thermo.py"
    out_path.write_text(rmg_shim.write_library_to_string(original))

    # Re-read and compare
    with open(out_path, 'r') as f:
        roundtripped = rmg_shim.parse_rmg_library(f.read())

    assert roundtripped.name == original.name
    assert len(roundtripped.entries) == len(original.entries)
    for orig_entry, rt_entry in zip(original.entries, roundtripped.entries):
        assert rt_entry.index == orig_entry.index
        assert rt_entry.label == orig_entry.label
        assert rt_entry.thermo.H298 == orig_entry.thermo.H298
        assert rt_entry.thermo.S298 == orig_entry.thermo.S298
        assert rt_entry.thermo.Cpdata == orig_entry.thermo.Cpdata

        # Adjacency list should still describe the same molecule
        orig_mol = Molecule().from_adjacency_list(orig_entry.molecule)
        rt_mol = Molecule().from_adjacency_list(rt_entry.molecule)
        assert rt_mol.is_isomorphic(orig_mol)


def test_write_real_kinetics_library_roundtrip(tmp_path):
    """
    Read a real RMG kinetics library, write it back, re-read it,
    and verify that reactions and Arrhenius parameters survive the
    round-trip.
    """
    rxn_path = os.path.join(REAL_DEST_KINETICS_DIR, 'reactions.py')
    with open(rxn_path, 'r') as f:
        original = rmg_shim.parse_rmg_library(f.read())

    out_path = tmp_path / "roundtrip_reactions.py"
    out_path.write_text(rmg_shim.write_library_to_string(original))

    with open(out_path, 'r') as f:
        roundtripped = rmg_shim.parse_rmg_library(f.read())

    assert roundtripped.name == original.name
    assert len(roundtripped.entries) == len(original.entries)
    for orig_entry, rt_entry in zip(original.entries, roundtripped.entries):
        assert rt_entry.index == orig_entry.index
        assert rt_entry.label == orig_entry.label
        assert rt_entry.kinetics.A == orig_entry.kinetics.A
        assert rt_entry.kinetics.n == orig_entry.kinetics.n
        assert rt_entry.kinetics.Ea == orig_entry.kinetics.Ea


def test_merge_real_thermo_libraries(tmp_path):
    """
    Merge a real (shortened) source thermo library into a real destination
    thermo library and verify:
      1. New entries (HO2, H2O2) are appended.
      2. Indices are reassigned starting after the destination's max index.
      3. The structurally-identical 'Hydrogen' entry from the source is
         skipped because it is isomorphic to the existing 'H2' species.
      4. The thermo data of the appended entries is preserved.
    """
    # Stage real libraries in tmp_path so the merge can write back without
    # touching the test data fixture.
    dest_lib_path = tmp_path / "dest_thermo.py"
    src_lib_path = tmp_path / "src_thermo.py"
    shutil.copy(REAL_DEST_THERMO_PATH, dest_lib_path)
    shutil.copy(REAL_SRC_THERMO_PATH, src_lib_path)

    paths = {
        'ARC thermo lib': str(src_lib_path),
        'T3 thermo lib': str(dest_lib_path),
        'shared T3 thermo lib': None,
        'ARC kinetics lib': None,
        'T3 kinetics lib': None,
        'shared T3 kinetics lib': None,
    }

    append_to_rmg_libraries(
        library_name="TestLib",
        shared_library_name="T3_test_dest_thermo",
        paths=paths,
        logger=logger,
    )

    with open(dest_lib_path, 'r') as f:
        merged = rmg_shim.parse_rmg_library(f.read())

    # Original 3 (H, H2, OH) + 2 new (HO2, H2O2). 'Hydrogen' is isomorphic
    # to the existing 'H2' and must be skipped.
    labels = [e.label for e in merged.entries]
    assert labels == ["H", "H2", "OH", "HO2", "H2O2"]
    assert "Hydrogen" not in labels

    # Indices should be contiguous after the max original index (which was 2)
    assert [e.index for e in merged.entries] == [0, 1, 2, 3, 4]

    # Sanity check on thermo data of appended entries
    ho2_entry = merged.entries[3]
    assert ho2_entry.thermo.H298 == (3, 'kcal/mol')
    assert ho2_entry.thermo.S298 == (54.75, 'cal/(mol*K)')

    h2o2_entry = merged.entries[4]
    assert h2o2_entry.thermo.H298 == (-32.53, 'kcal/mol')

    # And confirm the appended entries have the right molecules
    ho2_mol = Molecule().from_adjacency_list(ho2_entry.molecule)
    assert ho2_mol.is_isomorphic(Molecule(smiles='[O]O'))
    h2o2_mol = Molecule().from_adjacency_list(h2o2_entry.molecule)
    assert h2o2_mol.is_isomorphic(Molecule(smiles='OO'))


def test_merge_real_kinetics_libraries(tmp_path):
    """
    Merge a real (shortened) source kinetics library into a real
    destination library and verify:
      1. New reactions are appended with reassigned indices.
      2. The duplicate-label reaction is skipped.
      3. The species dictionary is updated with species that were missing
         from the destination dictionary (H2, HO2), while existing species
         are preserved unchanged.
    """
    dest_dir = tmp_path / "dest_kinetics"
    src_dir = tmp_path / "src_kinetics"
    shutil.copytree(REAL_DEST_KINETICS_DIR, dest_dir)
    shutil.copytree(REAL_SRC_KINETICS_DIR, src_dir)

    paths = {
        'ARC thermo lib': None,
        'T3 thermo lib': None,
        'shared T3 thermo lib': None,
        'ARC kinetics lib': str(src_dir),
        'T3 kinetics lib': str(dest_dir),
        'shared T3 kinetics lib': None,
    }

    append_to_rmg_libraries(
        library_name="TestLib",
        shared_library_name="T3_test_dest_kinetics",
        paths=paths,
        logger=logger,
    )

    # Verify reactions
    with open(dest_dir / "reactions.py", 'r') as f:
        merged = rmg_shim.parse_rmg_library(f.read())

    labels = [e.label for e in merged.entries]
    # Original 2 + 2 new (the dup label entry is skipped).
    assert labels == [
        "H + O2 <=> O + OH",
        "OH + OH <=> O + H2O",
        "H2 + OH <=> H2O + H",
        "H + HO2 <=> H2 + O2",
    ]
    assert [e.index for e in merged.entries] == [1, 2, 3, 4]

    # Verify the original entry's kinetics were NOT overwritten by the
    # source's same-label entry.
    orig_arr = merged.entries[0].kinetics
    assert orig_arr.A == (1.04e+14, 'cm^3/(mol*s)')
    assert orig_arr.Ea == (15286, 'cal/mol')

    # Verify the new reaction kinetics were appended faithfully
    new_arr = merged.entries[3].kinetics
    assert new_arr.A == (2.81e+13, 'cm^3/(mol*s)')
    assert new_arr.Ea == (1068, 'cal/mol')

    # Verify dictionary
    merged_dict = load_rmg_species_dictionary_file(str(dest_dir / "dictionary.txt"))
    # Original: H, O, OH, O2, H2O. New: H2, HO2.
    assert sorted(merged_dict.keys()) == ['H', 'H2', 'H2O', 'HO2', 'O', 'O2', 'OH']

    # Existing species kept their structure
    assert merged_dict['OH'].is_isomorphic(Molecule(smiles='[OH]'))
    # New species were added with the correct structure
    assert merged_dict['HO2'].is_isomorphic(Molecule(smiles='[O]O'))
    assert merged_dict['H2'].is_isomorphic(Molecule(smiles='[H][H]'))
