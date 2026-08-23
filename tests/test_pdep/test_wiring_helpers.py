#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_wiring_helpers module

Guards the memoisation in ``_wiring_helpers.build_t3()``. That helper is called by ~100 tests in
this directory, each of which used to re-parse the same 1.55 MB chem_annotated.yaml; the parse is
now done once and every caller is handed a ``deepcopy``. The entire correctness of that trade
rests on two properties, which are what this module tests:

1. Fidelity -- the cached build is indistinguishable from a real end-to-end parse.
2. Isolation -- a caller that mutates what it was handed cannot affect the next caller.

Without (2), tests in this directory would start passing or failing for reasons unrelated to what
they assert, which is far worse than a slow suite.
"""

from tests.test_pdep._wiring_helpers import (C4ENE_INDEX,
                                             C4RAD_INDEX,
                                             H_INDEX,
                                             build_t3)


def _labels(t3):
    return [species.label for species in t3.rmg_species]


def _adjacency_lists(t3):
    return [species.mol.to_adjacency_list() if species.mol is not None else None
            for species in t3.rmg_species]


def _reaction_signatures(t3):
    """A cheap structural fingerprint of every reaction.

    ``str(T3Reaction)`` costs ~0.15 s per reaction (~230 s over the 1461 reactions of this
    fixture), so comparing reactant/product labels and the pressure-dependence flag is used
    instead of stringifying the whole set.
    """
    return [(tuple(species.label for species in reaction.r_species),
             tuple(species.label for species in reaction.p_species),
             reaction.is_pressure_dependent)
            for reaction in t3.rmg_reactions]


class TestBuildT3Cache(object):
    """Test that the memoised build_t3() is faithful to, and isolated from, a real parse."""

    def test_cached_build_matches_a_real_parse(self, tmp_path):
        """The cached model is identical to one parsed end-to-end from the fixture file."""
        uncached = build_t3(tmp_path / 'uncached', use_cache=False)
        cached = build_t3(tmp_path / 'cached', use_cache=True)

        assert len(cached.rmg_species) == len(uncached.rmg_species)
        assert len(cached.rmg_reactions) == len(uncached.rmg_reactions)
        assert _labels(cached) == _labels(uncached)
        assert _adjacency_lists(cached) == _adjacency_lists(uncached)
        assert _reaction_signatures(cached) == _reaction_signatures(uncached)

    def test_each_caller_gets_distinct_objects(self, tmp_path):
        """Two builds share no species, reaction, molecule or atom objects."""
        first = build_t3(tmp_path / 'first')
        second = build_t3(tmp_path / 'second')

        assert first.rmg_species is not second.rmg_species
        assert first.rmg_reactions is not second.rmg_reactions
        for index in (H_INDEX, C4ENE_INDEX, C4RAD_INDEX):
            assert first.rmg_species[index] is not second.rmg_species[index]
            assert first.rmg_species[index].mol is not second.rmg_species[index].mol
            first_atoms = {id(atom) for atom in first.rmg_species[index].mol.atoms}
            second_atoms = {id(atom) for atom in second.rmg_species[index].mol.atoms}
            assert not first_atoms & second_atoms
        assert first.rmg_reactions[0] is not second.rmg_reactions[0]

    def test_object_sharing_within_a_build_is_preserved(self, tmp_path):
        """A species referenced by a reaction is the same object as the one in rmg_species.

        deepcopy of the (species, reactions) pair in one call preserves the internal aliasing a
        real parse produces; copying the two lists separately would not, and would silently
        change what these tests exercise.
        """
        t3 = build_t3(tmp_path)
        by_label = {species.label: species for species in t3.rmg_species}
        checked = 0
        for reaction in t3.rmg_reactions:
            for species in list(reaction.r_species) + list(reaction.p_species):
                if species.label in by_label:
                    assert species is by_label[species.label]
                    checked += 1
        assert checked, 'no reaction species were matched against the species list'

    def test_mutating_a_species_does_not_leak_to_the_next_caller(self, tmp_path):
        """Mutating labels, attributes and molecules of one build leaves the next build pristine."""
        first = build_t3(tmp_path / 'first')
        original_label = first.rmg_species[H_INDEX].label
        original_adjacency_list = first.rmg_species[H_INDEX].mol.to_adjacency_list()
        original_species_count = len(first.rmg_species)
        original_reaction_count = len(first.rmg_reactions)

        first.rmg_species[H_INDEX].label = 'MUTATED'
        first.rmg_species[H_INDEX].reasons = ['mutated reason']
        first.rmg_species[H_INDEX].thermo = 'mutated thermo'
        first.rmg_species[H_INDEX].mol.atoms[0].element = \
            first.rmg_species[C4ENE_INDEX].mol.atoms[0].element
        first.rmg_species.pop()
        first.rmg_reactions.pop()

        second = build_t3(tmp_path / 'second')
        assert second.rmg_species[H_INDEX].label == original_label
        assert second.rmg_species[H_INDEX].reasons != ['mutated reason']
        assert second.rmg_species[H_INDEX].thermo != 'mutated thermo'
        assert second.rmg_species[H_INDEX].mol.to_adjacency_list() == original_adjacency_list
        assert len(second.rmg_species) == original_species_count
        assert len(second.rmg_reactions) == original_reaction_count

    def test_mutating_a_reaction_does_not_leak_to_the_next_caller(self, tmp_path):
        """Mutating a reaction of one build leaves the next build's reactions pristine."""
        first = build_t3(tmp_path / 'first')
        original_signature = _reaction_signatures(first)[0]

        first.rmg_reactions[0].r_species = []
        first.rmg_reactions[0].is_pressure_dependent = not first.rmg_reactions[0].is_pressure_dependent

        second = build_t3(tmp_path / 'second')
        assert _reaction_signatures(second)[0] == original_signature
