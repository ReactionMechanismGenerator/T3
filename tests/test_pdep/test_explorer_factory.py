#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_explorer_factory module
"""

import pytest

from t3.pdep.explorer.adapter import PESExplorerAdapter
from t3.pdep.explorer.factory import (
    _registered_explorer_adapters,
    explorer_factory,
    register_explorer_adapter,
)


class _DummyPESExplorerAdapter(PESExplorerAdapter):
    """
    A minimal, fully-implemented PESExplorerAdapter used only for testing the factory/registry.

    Mirrors the design's ONE deliberate divergence from the mesolver precedent: the seed/
    capability rules (``max_source_species``, ``supports_transition_state_seeds``) are enforced in
    the adapter layer rather than in the factory, which makes direct construction just as safe as
    going through ``explorer_factory``.

    Note that this dummy deliberately does NOT re-implement those rules: it delegates to
    ``PESExplorerAdapter.__init__`` exactly as a real adapter does. A test double that carries its
    own copy of the logic would only ever prove that the copy works.
    """

    def __init__(self,
                 seed_species,
                 output_directory: str,
                 bath_gas: dict = None,
                 explore_tol: float = None,
                 energy_tol: float = None,
                 flux_tol: float = None,
                 maximum_radical_electrons: int = None,
                 logger=None,
                 transition_state_seeds: tuple = None,
                 database_kwargs: dict = None,
                 ):
        super().__init__(seed_species=seed_species, transition_state_seeds=transition_state_seeds)
        self.output_directory = output_directory
        self.bath_gas = bath_gas
        self.explore_tol = explore_tol
        self.energy_tol = energy_tol
        self.flux_tol = flux_tol
        self.maximum_radical_electrons = maximum_radical_electrons
        self.logger = logger
        self.transition_state_seeds = transition_state_seeds
        self.database_kwargs = database_kwargs
        self.explored = False

    def set_up(self):
        pass

    def explore(self):
        self.explored = True
        return True

    def get_networks(self):
        if not self.explored:
            raise RuntimeError('get_networks() was called before explore(); no exploration has been attempted yet.')
        return tuple()

    def get_k_tp(self):
        if not self.explored:
            raise RuntimeError('get_k_tp() was called before explore(); no exploration has been attempted yet.')
        return dict()


class _DummyTSSeedCompatibleAdapter(_DummyPESExplorerAdapter):
    """A dummy adapter that declares support for transition-state-seeded exploration."""

    supports_transition_state_seeds = True


class _NotAnAdapter(object):
    """A plain class that does not inherit from PESExplorerAdapter, used to test the TypeError guard."""

    def set_up(self):
        pass

    def explore(self):
        return True

    def get_networks(self):
        return tuple()

    def get_k_tp(self):
        return dict()


class TestPESExplorerAdapterABC(object):

    def test_cannot_instantiate_abstract_class_directly(self):
        """PESExplorerAdapter is an ABC and cannot be instantiated directly."""
        with pytest.raises(TypeError):
            PESExplorerAdapter()

    def test_max_source_species_defaults_to_two(self):
        """max_source_species defaults to 2 (Arkane's hard limit) on the base class."""
        assert PESExplorerAdapter.max_source_species == 2

    def test_supports_transition_state_seeds_defaults_to_false(self):
        """supports_transition_state_seeds defaults to False (Arkane: never) on the base class."""
        assert PESExplorerAdapter.supports_transition_state_seeds is False

    def test_dummy_subclass_can_be_instantiated(self):
        """A minimal concrete subclass implementing all abstract methods can be instantiated."""
        adapter = _DummyPESExplorerAdapter(seed_species=['OH'], output_directory='out')
        assert isinstance(adapter, PESExplorerAdapter)

    def test_an_empty_seed_is_refused(self):
        """
        An empty seed is refused by the base class, and is NOT reported as a too-large seed.

        Arkane reaches the same ``else`` branch for zero source species as for three or more, and
        reports both as "Reactant channels with greater than 2 reactants not supported"
        (``arkane/explorer.py:159-169``). That message is actively misleading for an empty seed, so
        T3 refuses it here with its own accurate one.
        """
        with pytest.raises(ValueError, match='at least one source species'):
            _DummyPESExplorerAdapter(seed_species=[], output_directory='out')

    def test_seed_rules_are_enforced_by_the_base_class_not_by_each_subclass(self):
        """
        The seed rules live in ``PESExplorerAdapter.__init__`` itself, not in each concrete adapter.

        This is the point of the test: a subclass that implements only the four abstract methods
        and re-implements NO validation of its own still gets every seed rule enforced. Without
        this, ``max_source_species``/``supports_transition_state_seeds`` would be documentation
        plus a test double -- the factory deliberately performs no capability check, so a concrete
        adapter whose author simply forgot to validate would have no enforcement anywhere.
        """
        class _BarelyImplemented(PESExplorerAdapter):
            def set_up(self):
                pass

            def explore(self):
                return True

            def get_networks(self):
                return tuple()

            def get_k_tp(self):
                return dict()

        assert _BarelyImplemented(seed_species=['OH']).seed_species == ('OH',)
        with pytest.raises(ValueError, match='at most 2 source species'):
            _BarelyImplemented(seed_species=['OH', 'H', 'O2'])
        with pytest.raises(ValueError, match='transition-state-seeded'):
            _BarelyImplemented(seed_species=['OH'], transition_state_seeds=['TS1'])


class TestRegisterExplorerAdapter(object):

    def setup_method(self):
        """Snapshot the registry so each test can restore it, keeping tests order-independent."""
        self._original_registry = dict(_registered_explorer_adapters)

    def teardown_method(self):
        """Restore the registry to its pre-test state, removing any keys added by the test."""
        _registered_explorer_adapters.clear()
        _registered_explorer_adapters.update(self._original_registry)

    def test_register_dummy_adapter(self):
        """register_explorer_adapter registers a valid PESExplorerAdapter subclass cleanly."""
        register_explorer_adapter('DummyExplorer', _DummyPESExplorerAdapter)
        assert _registered_explorer_adapters['DummyExplorer'] is _DummyPESExplorerAdapter

    def test_register_non_adapter_raises_type_error(self):
        """register_explorer_adapter raises TypeError for a class that is not a PESExplorerAdapter subclass."""
        with pytest.raises(TypeError):
            register_explorer_adapter('NotAnAdapter', _NotAnAdapter)


class TestExplorerFactory(object):

    def setup_method(self):
        """Snapshot the registry and register the dummy adapters needed by this test class."""
        self._original_registry = dict(_registered_explorer_adapters)
        register_explorer_adapter('DummyExplorer', _DummyPESExplorerAdapter)
        register_explorer_adapter('DummyTSExplorer', _DummyTSSeedCompatibleAdapter)

    def teardown_method(self):
        """Restore the registry to its pre-test state, removing any keys added by the test."""
        _registered_explorer_adapters.clear()
        _registered_explorer_adapters.update(self._original_registry)

    def test_factory_returns_instance_of_registered_class(self):
        """explorer_factory returns an instance of the registered class."""
        adapter = explorer_factory(explorer='DummyExplorer',
                                   seed_species=['OH'],
                                   output_directory='out',
                                   )
        assert isinstance(adapter, _DummyPESExplorerAdapter)

    def test_factory_lookup_is_case_sensitive(self):
        """explorer_factory lookup is case-sensitive; a wrong-case name raises ValueError."""
        with pytest.raises(ValueError):
            explorer_factory(explorer='dummyexplorer',
                             seed_species=['OH'],
                             output_directory='out',
                             )

    def test_factory_unknown_name_raises_value_error_listing_keys(self):
        """The ValueError message for an unknown explorer name mentions the registered keys."""
        with pytest.raises(ValueError) as excinfo:
            explorer_factory(explorer='NoSuchExplorer',
                             seed_species=['OH'],
                             output_directory='out',
                             )
        assert 'DummyExplorer' in str(excinfo.value)
        assert 'DummyTSExplorer' in str(excinfo.value)

    def test_factory_routed_invalid_seed_raises(self):
        """Too many source species, routed through the factory, raises ValueError (adapter-enforced)."""
        with pytest.raises(ValueError):
            explorer_factory(explorer='DummyExplorer',
                             seed_species=['OH', 'H', 'O2'],
                             output_directory='out',
                             )

    def test_direct_construction_invalid_seed_is_not_a_bypass(self):
        """
        Constructing the concrete adapter directly (bypassing explorer_factory entirely) with an
        invalid seed still raises. This proves the seed/capability rule is enforced in the
        adapter itself, not only in the factory, per the design's deliberate divergence from the
        mesolver precedent (whose allow_ilt_complement rule lives in the factory ONLY).
        """
        with pytest.raises(ValueError):
            _DummyPESExplorerAdapter(seed_species=['OH', 'H', 'O2'], output_directory='out')

    def test_direct_construction_ts_seed_unsupported_is_not_a_bypass(self):
        """Direct construction with transition_state_seeds against a non-supporting adapter still raises."""
        with pytest.raises(ValueError):
            _DummyPESExplorerAdapter(seed_species=['OH'], output_directory='out', transition_state_seeds=('TS1',))

    def test_factory_ts_seed_compatible_adapter_succeeds(self):
        """A registered adapter that supports transition-state seeds accepts them via the factory."""
        adapter = explorer_factory(explorer='DummyTSExplorer',
                                   seed_species=['OH'],
                                   output_directory='out',
                                   transition_state_seeds=('TS1',),
                                   )
        assert isinstance(adapter, _DummyTSSeedCompatibleAdapter)

    def test_get_networks_before_explore_raises_runtime_error(self):
        """get_networks() before a successful explore() raises RuntimeError, never returns empty results silently."""
        adapter = explorer_factory(explorer='DummyExplorer', seed_species=['OH'], output_directory='out')
        with pytest.raises(RuntimeError):
            adapter.get_networks()

    def test_get_k_tp_before_explore_raises_runtime_error(self):
        """get_k_tp() before a successful explore() raises RuntimeError, never returns empty results silently."""
        adapter = explorer_factory(explorer='DummyExplorer', seed_species=['OH'], output_directory='out')
        with pytest.raises(RuntimeError):
            adapter.get_k_tp()

    def test_get_networks_and_get_k_tp_after_explore_succeed(self):
        """get_networks() and get_k_tp() work after a successful explore()."""
        adapter = explorer_factory(explorer='DummyExplorer', seed_species=['OH'], output_directory='out')
        assert adapter.explore() is True
        assert adapter.get_networks() == tuple()
        assert adapter.get_k_tp() == dict()
