#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_mesolver_factory module
"""

import pytest

from t3.pdep.mesolver.adapter import MESolverAdapter
from t3.pdep.mesolver.factory import (
    _registered_mesolver_adapters,
    mesolver_factory,
    register_mesolver_adapter,
)


class _DummyMESolverAdapter(MESolverAdapter):
    """A minimal, fully-implemented MESolverAdapter used only for testing the factory/registry."""

    def __init__(self,
                 network_path: str,
                 output_directory: str,
                 method: str,
                 logger=None,
                 allow_ilt_complement: bool = False,
                 expected_reactions: int | set | None = None,
                 ):
        self.network_path = network_path
        self.output_directory = output_directory
        self.method = method
        self.logger = logger
        self.allow_ilt_complement = allow_ilt_complement
        self.expected_reactions = expected_reactions

    def set_up(self):
        pass

    def solve(self):
        return True

    def get_k_tp(self):
        return dict()


class _DummyILTCompatibleAdapter(_DummyMESolverAdapter):
    """A dummy adapter that declares support for a partially QM-computed (ILT-complement) network."""

    supports_ilt_complement = True


class _NotAnAdapter(object):
    """A plain class that does not inherit from MESolverAdapter, used to test the TypeError guard."""

    def set_up(self):
        pass

    def solve(self):
        return True

    def get_k_tp(self):
        return dict()


class TestMESolverAdapterABC(object):

    def test_cannot_instantiate_abstract_class_directly(self):
        """MESolverAdapter is an ABC and cannot be instantiated directly."""
        with pytest.raises(TypeError):
            MESolverAdapter()

    def test_supports_ilt_complement_defaults_to_false(self):
        """supports_ilt_complement defaults to False on the base class."""
        assert MESolverAdapter.supports_ilt_complement is False

    def test_dummy_subclass_can_be_instantiated(self):
        """A minimal concrete subclass implementing all abstract methods can be instantiated."""
        adapter = _DummyMESolverAdapter(network_path='network.txt', output_directory='out', method='CSE')
        assert isinstance(adapter, MESolverAdapter)


class TestRegisterMESolverAdapter(object):

    def setup_method(self):
        """Snapshot the registry so each test can restore it, keeping tests order-independent."""
        self._original_registry = dict(_registered_mesolver_adapters)

    def teardown_method(self):
        """Restore the registry to its pre-test state, removing any keys added by the test."""
        _registered_mesolver_adapters.clear()
        _registered_mesolver_adapters.update(self._original_registry)

    def test_register_dummy_adapter(self):
        """register_mesolver_adapter registers a valid MESolverAdapter subclass cleanly."""
        register_mesolver_adapter('DummyMESolver', _DummyMESolverAdapter)
        assert _registered_mesolver_adapters['DummyMESolver'] is _DummyMESolverAdapter

    def test_register_non_adapter_raises_type_error(self):
        """register_mesolver_adapter raises TypeError for a class that is not a MESolverAdapter subclass."""
        with pytest.raises(TypeError):
            register_mesolver_adapter('NotAnAdapter', _NotAnAdapter)


class TestMESolverFactory(object):

    def setup_method(self):
        """Snapshot the registry and register the dummy adapters needed by this test class."""
        self._original_registry = dict(_registered_mesolver_adapters)
        register_mesolver_adapter('DummyMESolver', _DummyMESolverAdapter)
        register_mesolver_adapter('DummyILTMESolver', _DummyILTCompatibleAdapter)

    def teardown_method(self):
        """Restore the registry to its pre-test state, removing any keys added by the test."""
        _registered_mesolver_adapters.clear()
        _registered_mesolver_adapters.update(self._original_registry)

    def test_factory_returns_instance_of_registered_class(self):
        """mesolver_factory returns an instance of the registered class."""
        adapter = mesolver_factory(me_solver='DummyMESolver',
                                   network_path='network.txt',
                                   output_directory='out',
                                   method='CSE',
                                   )
        assert isinstance(adapter, _DummyMESolverAdapter)

    def test_factory_lookup_is_case_sensitive(self):
        """mesolver_factory lookup is case-sensitive; a wrong-case name raises ValueError."""
        with pytest.raises(ValueError):
            mesolver_factory(me_solver='dummymesolver',
                             network_path='network.txt',
                             output_directory='out',
                             method='CSE',
                             )

    def test_factory_unknown_name_raises_value_error_listing_keys(self):
        """The ValueError message for an unknown me_solver name mentions the registered keys."""
        with pytest.raises(ValueError) as excinfo:
            mesolver_factory(me_solver='NoSuchSolver',
                             network_path='network.txt',
                             output_directory='out',
                             method='CSE',
                             )
        assert 'DummyMESolver' in str(excinfo.value)
        assert 'DummyILTMESolver' in str(excinfo.value)

    def test_allow_ilt_complement_true_against_unsupported_adapter_raises(self):
        """allow_ilt_complement=True against an adapter with supports_ilt_complement=False raises ValueError."""
        with pytest.raises(ValueError):
            mesolver_factory(me_solver='DummyMESolver',
                             network_path='network.txt',
                             output_directory='out',
                             method='CSE',
                             allow_ilt_complement=True,
                             )

    def test_allow_ilt_complement_true_against_supported_adapter_succeeds(self):
        """allow_ilt_complement=True against an adapter with supports_ilt_complement=True succeeds."""
        adapter = mesolver_factory(me_solver='DummyILTMESolver',
                                   network_path='network.txt',
                                   output_directory='out',
                                   method='CSE',
                                   allow_ilt_complement=True,
                                   )
        assert isinstance(adapter, _DummyILTCompatibleAdapter)
