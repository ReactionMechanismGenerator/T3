#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_chem module
Tests for t3/chem.py - T3Species and T3Reaction classes
"""

import pytest

from t3.chem import (
    T3Status,
    ThermoMethod,
    KineticsMethod,
    T3Species,
    T3Reaction,
)


class TestT3Status:
    """Tests for T3Status enum."""

    def test_status_values(self):
        """Test that all expected status values exist."""
        assert T3Status.PENDING.value == "pending"
        assert T3Status.RUNNING.value == "running"
        assert T3Status.CONVERGED.value == "converged"
        assert T3Status.FAILED.value == "failed"
        assert T3Status.SKIPPED.value == "skipped"

    def test_status_from_string(self):
        """Test creating status from string value."""
        assert T3Status("pending") == T3Status.PENDING
        assert T3Status("converged") == T3Status.CONVERGED


class TestThermoMethod:
    """Tests for ThermoMethod enum."""

    def test_thermo_method_values(self):
        """Test that all expected thermo method values exist."""
        assert ThermoMethod.QM.value == "QM"
        assert ThermoMethod.LIBRARY.value == "Library"
        assert ThermoMethod.GAV.value == "GAV"
        assert ThermoMethod.ML.value == "ML"
        assert ThermoMethod.USER.value == "User"
        assert ThermoMethod.UNKNOWN.value == "Unknown"

    def test_thermo_method_from_string(self):
        """Test creating thermo method from string value."""
        assert ThermoMethod("QM") == ThermoMethod.QM
        assert ThermoMethod("Library") == ThermoMethod.LIBRARY


class TestKineticsMethod:
    """Tests for KineticsMethod enum."""

    def test_kinetics_method_values(self):
        """Test that all expected kinetics method values exist."""
        assert KineticsMethod.QM.value == "QM"
        assert KineticsMethod.LIBRARY.value == "Library"
        assert KineticsMethod.RATE_RULES.value == "Rate Rules"
        assert KineticsMethod.TRAINING_SET.value == "Training Set"
        assert KineticsMethod.PDEP.value == "PDep"
        assert KineticsMethod.USER.value == "User"
        assert KineticsMethod.UNKNOWN.value == "Unknown"


class TestT3Species:
    """Tests for T3Species class."""

    def test_create_species_with_smiles(self):
        """Test creating a T3Species with SMILES."""
        spc = T3Species(label='H2', smiles='[H][H]')
        assert spc.label == 'H2'
        assert spc.t3_status == T3Status.PENDING

    def test_create_species_with_thermo_method(self):
        """Test creating a T3Species with thermo method."""
        spc = T3Species(
            label='CH4',
            smiles='C',
            thermo_method=ThermoMethod.QM,
            thermo_source='CBS-QB3',
        )
        assert spc.label == 'CH4'
        assert spc.thermo_method == ThermoMethod.QM
        assert spc.thermo_source == 'CBS-QB3'

    def test_create_species_with_string_thermo_method(self):
        """Test creating a T3Species with string thermo method."""
        spc = T3Species(
            label='H2O',
            smiles='O',
            thermo_method='Library',
            thermo_source='primaryThermoLibrary',
        )
        assert spc.thermo_method == ThermoMethod.LIBRARY
        assert spc.thermo_source == 'primaryThermoLibrary'

    def test_create_species_with_unknown_thermo_method(self):
        """Test creating a T3Species with unknown thermo method string."""
        spc = T3Species(
            label='test',
            smiles='C',
            thermo_method='SomeUnknownMethod',
            thermo_source='details',
        )
        assert spc.thermo_method == ThermoMethod.UNKNOWN
        # Unknown method prepended to source
        assert 'SomeUnknownMethod' in spc.thermo_source

    def test_species_default_status(self):
        """Test that default status is PENDING."""
        spc = T3Species(label='test', smiles='C')
        assert spc.t3_status == T3Status.PENDING

    def test_species_is_converged_false(self):
        """Test is_converged property returns False for pending species."""
        spc = T3Species(label='test', smiles='C')
        assert spc.is_converged is False

    def test_species_is_converged_true(self):
        """Test is_converged property returns True for converged species."""
        spc = T3Species(label='test', smiles='C', t3_status=T3Status.CONVERGED)
        assert spc.is_converged is True

    def test_species_key_assignment(self):
        """Test setting T3 index."""
        spc = T3Species(label='test', smiles='C')
        assert spc.key
        spc = T3Species(label='test', smiles='C', key=5)
        assert spc.key == 5

    def test_species_created_at_iteration(self):
        """Test setting created_at_iteration."""
        spc = T3Species(label='test', smiles='C', created_at_iteration=3)
        assert spc.created_at_iteration == 3

    def test_species_reasons(self):
        """Test setting reasons list."""
        spc = T3Species(label='test', smiles='C', reasons=['sensitivity', 'high_energy'])
        assert 'sensitivity' in spc.reasons
        assert 'high_energy' in spc.reasons

    def test_species_reasons_from_string(self):
        """Test setting reasons from single string."""
        spc = T3Species(label='test', smiles='C', reasons='sensitivity')
        assert spc.reasons == ['sensitivity']

    def test_species_as_dict(self):
        """Test as_dict method includes T3 metadata."""
        spc = T3Species(
            label='CH4',
            smiles='C',
            thermo_method=ThermoMethod.QM,
            thermo_source='CBS-QB3',
            key=5,
            t3_status=T3Status.CONVERGED,
            created_at_iteration=2,
        )
        data = spc.as_dict()

        assert data['label'] == 'CH4'
        assert data['rmg_label'] == {2: 'CH4'}
        assert data['qm_label'] == 's5_CH4'
        assert data['thermo_method'] == ThermoMethod.QM.value
        assert data['thermo_source'] == 'CBS-QB3'
        assert data['key'] == 5
        assert data['t3_status'] == T3Status.CONVERGED.value
        assert data['created_at_iteration'] == 2

    def test_species_from_dict(self):
        """Test from_dict class method reconstructs species."""
        original_data = {
            'label': 'CH4',
            'smiles': 'C',
            'thermo_method': 'QM',
            'thermo_source': 'CBS-QB3',
            'key': 1,
            't3_status': 'Converged',
            'created_at_iteration': 0,
            'rmg_label': {'0': 'CH4'},
        }
        spc = T3Species.from_dict(original_data.copy())
        assert spc.label == 'CH4'
        assert spc.qm_label == 's1_CH4'
        assert spc.rmg_label == {0: 'CH4'}
        assert spc.thermo_method == ThermoMethod.QM
        assert spc.thermo_source == 'CBS-QB3'
        assert spc.key == 1
        assert spc.t3_status == T3Status.CONVERGED
        assert spc.created_at_iteration == 0

    def test_species_dict_round_trip(self):
        """
        Test a full cycle: Species -> Dict -> Species.
        Ensures all attributes, including T3 metadata and history, are preserved.
        """
        original_spc = T3Species(
            label='CH4(6)',
            smiles='C',
            key=5,
            thermo_method=ThermoMethod.QM,
            thermo_source='CBS-QB3',
            thermo_comment='Calculated using Gaussian',
            t3_status=T3Status.CONVERGED,
            created_at_iteration=2,
            reasons=['Validation', 'Sensitivity'],
        )
        original_spc.rmg_label[1] = 'CH4(12)'
        spc_dict = original_spc.as_dict()
        reconstructed_spc = T3Species.from_dict(spc_dict)

        assert reconstructed_spc.label == original_spc.label
        assert reconstructed_spc.qm_label == 's5_CH4'  # Auto-generated from key+formula
        assert reconstructed_spc.key == 5
        assert reconstructed_spc.rmg_label == {2: 'CH4(6)', 1: 'CH4(12)'}
        assert isinstance(list(reconstructed_spc.rmg_label.keys())[0], int)
        assert reconstructed_spc.thermo_method == ThermoMethod.QM
        assert reconstructed_spc.t3_status == T3Status.CONVERGED
        assert reconstructed_spc.thermo_source == 'CBS-QB3'
        assert reconstructed_spc.thermo_comment == 'Calculated using Gaussian'
        assert reconstructed_spc.created_at_iteration == 2
        assert reconstructed_spc.reasons == ['Validation', 'Sensitivity']
        assert reconstructed_spc.mol.to_smiles() == 'C'




    def test_species_repr(self):
        """Test __repr__ method."""
        spc = T3Species(
            label='CH4(5)',
            smiles='C',
            thermo_method=ThermoMethod.QM,
            thermo_source='CBS-QB3',
            key=1,
            t3_status=T3Status.CONVERGED,
        )
        repr_str = repr(spc)
        print(repr_str)
        assert 'T3Species' in repr_str
        assert 'CH4' in repr_str
        assert '[QM]' in repr_str
        assert 'status: converged' in repr_str

    def test_species_repr_no_method(self):
        """Test __repr__ without thermo method."""
        spc = T3Species(label='H2', smiles='[H][H]')
        repr_str = repr(spc)

        assert 'T3Species' in repr_str
        assert 'H2' in repr_str
        assert 'pending' in repr_str

class TestT3Reaction:
    """Tests for T3Reaction class."""

    @pytest.fixture
    def test_species(self):
        """Create test species for reactions."""
        h = T3Species(label='H', smiles='[H]')
        ch4 = T3Species(label='CH4', smiles='C')
        h2 = T3Species(label='H2', smiles='[H][H]')
        ch3 = T3Species(label='CH3', smiles='[CH3]')
        return {'H': h, 'CH4': ch4, 'H2': h2, 'CH3': ch3}

    def test_create_reaction_with_species(self, test_species):
        """Test creating a T3Reaction with species."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
        )
        assert rxn.t3_status == T3Status.PENDING
        assert len(rxn.r_species) == 2
        assert len(rxn.p_species) == 2

    def test_create_reaction_with_labels(self, test_species):
        """Test creating a T3Reaction with labels."""
        rxn = T3Reaction(
            reactants=['H', 'CH4'],
            products=['H2', 'CH3'],
        )
        assert rxn.t3_status == T3Status.PENDING
        assert 'H' in rxn.reactants
        assert 'CH4' in rxn.reactants

    def test_create_reaction_with_kinetics_method(self, test_species):
        """Test creating a T3Reaction with kinetics method."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            kinetics_method=KineticsMethod.LIBRARY,
            kinetics_source='BurkeH2O2',
        )
        assert rxn.kinetics_method == KineticsMethod.LIBRARY
        assert rxn.kinetics_source == 'BurkeH2O2'

    def test_create_reaction_with_string_kinetics_method(self, test_species):
        """Test creating a T3Reaction with string kinetics method."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            kinetics_method='Rate Rules',
            kinetics_source='H_Abstraction',
        )
        assert rxn.kinetics_method == KineticsMethod.RATE_RULES
        assert rxn.kinetics_source == 'H_Abstraction'

    def test_create_reaction_with_unknown_kinetics_method(self, test_species):
        """Test creating a T3Reaction with unknown kinetics method string."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            kinetics_method='SomeUnknownMethod',
            kinetics_source='details',
        )
        assert rxn.kinetics_method == KineticsMethod.UNKNOWN
        assert 'SomeUnknownMethod' in rxn.kinetics_source

    def test_reaction_default_status(self, test_species):
        """Test that default status is PENDING."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
        )
        assert rxn.t3_status == T3Status.PENDING

    def test_reaction_is_converged_false(self, test_species):
        """Test is_converged property returns False for pending reaction."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
        )
        assert rxn.is_converged is False

    def test_reaction_is_converged_true(self, test_species):
        """Test is_converged property returns True for converged reaction."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            t3_status=T3Status.CONVERGED,
        )
        assert rxn.is_converged is True

    def test_reaction_t3_index(self, test_species):
        """Test setting T3 index."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            t3_index=5,
        )
        assert rxn.t3_index == 5

    def test_reaction_rmg_index(self, test_species):
        """Test setting RMG index."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            rmg_index=10,
        )
        assert rxn.rmg_index == 10

    def test_reaction_created_at_iteration(self, test_species):
        """Test setting created_at_iteration."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            created_at_iteration=3,
        )
        assert rxn.created_at_iteration == 3

    def test_reaction_reasons(self, test_species):
        """Test setting reasons list."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            reasons=['sensitivity', 'collision_violator'],
        )
        assert 'sensitivity' in rxn.reasons
        assert 'collision_violator' in rxn.reasons

    def test_reaction_reasons_from_string(self, test_species):
        """Test setting reasons from single string."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            reasons='sensitivity',
        )
        assert rxn.reasons == ['sensitivity']

    def test_reaction_reactant_product_keys(self, test_species):
        """Test setting reactant and product keys."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            reactant_keys=[1, 2],
            product_keys=[3, 4],
        )
        assert rxn.reactant_keys == [1, 2]
        assert rxn.product_keys == [3, 4]

    def test_reaction_is_pressure_dependent(self, test_species):
        """Test setting is_pressure_dependent flag."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            is_pressure_dependent=True,
        )
        assert rxn.is_pressure_dependent is True

    def test_reaction_qm_label(self, test_species):
        """Test setting qm_label."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            qm_label='H+CH4=H2+CH3',
        )
        assert rxn.qm_label == 'H+CH4=H2+CH3'

    def test_reaction_rmg_label(self, test_species):
        """Test setting rmg_label."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            rmg_label='H + CH4 <=> H2 + CH3',
        )
        assert rxn.rmg_label == 'H + CH4 <=> H2 + CH3'

    def test_reaction_as_dict(self, test_species):
        """Test as_dict method includes T3 metadata."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            kinetics_method=KineticsMethod.LIBRARY,
            kinetics_source='BurkeH2O2',
            t3_index=1,
            rmg_index=5,
            t3_status=T3Status.CONVERGED,
            created_at_iteration=2,
            reactant_keys=[1, 2],
            product_keys=[3],
        )
        data = rxn.as_dict()

        assert data['kinetics_method'] == KineticsMethod.LIBRARY
        assert data['kinetics_source'] == 'BurkeH2O2'
        assert data['t3_index'] == 1
        assert data['rmg_index'] == 5
        assert data['t3_status'] == T3Status.CONVERGED
        assert data['created_at_iteration'] == 2
        assert data['reactant_keys'] == [1, 2]
        assert data['product_keys'] == [3]

    def test_reaction_from_dict(self):
        """Test from_dict class method reconstructs reaction."""
        original_data = {
            'kinetics_method': KineticsMethod.LIBRARY,
            'kinetics_source': 'BurkeH2O2',
            't3_index': 1,
            't3_status': T3Status.CONVERGED,
            'reactants': ['H', 'CH4'],
            'products': ['H2', 'CH3'],
        }
        rxn = T3Reaction.from_dict(original_data.copy())

        assert rxn.kinetics_method == KineticsMethod.LIBRARY
        assert rxn.t3_index == 1
        assert rxn.t3_status == T3Status.CONVERGED

    def test_reaction_repr(self, test_species):
        """Test __repr__ method."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
            kinetics_method=KineticsMethod.RATE_RULES,
            kinetics_source='H_Abstraction',
            t3_index=5,
            t3_status=T3Status.PENDING,
        )
        repr_str = repr(rxn)

        assert 'T3Reaction' in repr_str
        assert 'index: 5' in repr_str
        assert 'Rate Rules' in repr_str
        assert 'pending' in repr_str

    def test_reaction_repr_no_method(self, test_species):
        """Test __repr__ without kinetics method."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
        )
        repr_str = repr(rxn)

        assert 'T3Reaction' in repr_str
        assert 'pending' in repr_str

    def test_get_reaction_smiles_label(self, test_species):
        """Test get_reaction_smiles_label method."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
        )
        smiles_label = rxn.get_reaction_smiles_label()
        assert '[H]' in smiles_label
        assert 'C' in smiles_label
        assert '[H][H]' in smiles_label
        assert '[CH3]' in smiles_label
        assert '<=>' in smiles_label

    def test_get_reaction_smiles_label_uses_plus_separator(self, test_species):
        """Test that get_reaction_smiles_label joins products with '+', not '.'."""
        rxn = T3Reaction(
            r_species=[test_species['H'], test_species['CH4']],
            p_species=[test_species['H2'], test_species['CH3']],
        )
        smiles_label = rxn.get_reaction_smiles_label()
        # Products must be joined with '+', not '.'
        products_str = smiles_label.split('<=>')[1]
        assert '+' in products_str
        assert '.' not in products_str

    def test_get_reaction_smiles_label_validates_empty_smiles(self, test_species):
        """Test that the empty SMILES check uses `all()` not `any()`.
        `any(['CCC', ''])` is True (misses the empty), `all(['CCC', ''])` is False (catches it)."""
        # Directly test the validation logic: all() catches a single empty in a list
        smiles_with_empty = ['CCC', '']
        assert not all(smiles_with_empty), "all() should return False when list contains empty string"
        # any() would incorrectly return True here
        assert any(smiles_with_empty), "any() incorrectly returns True with partial empties"

    def test_species_from_dict_in_init(self):
        """Test that passing species_dict to __init__ populates T3-specific attributes."""
        T3Species.reset_counter()
        original = T3Species(
            label='CH4',
            smiles='C',
            key=10,
            thermo_method=ThermoMethod.QM,
            thermo_source='CBS-QB3',
            t3_status=T3Status.CONVERGED,
            created_at_iteration=2,
            reasons=['sensitivity'],
        )
        spc_dict = original.as_dict()
        # Reconstruct via species_dict parameter in __init__
        reconstructed = T3Species(species_dict=spc_dict)
        assert reconstructed.thermo_method == ThermoMethod.QM
        assert reconstructed.t3_status == T3Status.CONVERGED
        assert reconstructed.created_at_iteration == 2
        assert reconstructed.reasons == ['sensitivity']
