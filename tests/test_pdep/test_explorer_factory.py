#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_explorer_factory module
"""

import pytest

from t3.pdep.explorer.adapter import PESExplorerAdapter, validate_explorer_seed
from t3.pdep.explorer.factory import (
    _registered_explorer_adapters,
    explorer_factory,
    register_explorer_adapter,
)


class _DummyPESExplorerAdapter(PESExplorerAdapter):
    """
    A minimal, fully-implemented PESExplorerAdapter used only for testing the factory/registry.

    The seed/capability rules (``max_source_species``, ``supports_transition_state_seeds``) are
    enforced in the adapter layer, via ``PESExplorerAdapter.__init__``, so that direct construction
    is as safe as going through ``explorer_factory``. That was once described as the design's ONE
    divergence from the mesolver precedent, with the rationale that adapter-layer enforcement made a
    factory-layer check unnecessary. The rationale was incomplete: enforcement in ``__init__`` is
    bypassed by any subclass that overrides ``__init__`` without calling ``super()``, which
    registration cannot detect (``issubclass`` is still True). The factory therefore re-asserts the
    same shared rule at its own entry point -- the two checks cover two different entry paths and
    neither is redundant. See ``TestRegistrationCannotBlessAnAdapterThatSkipsSeedValidation``.

    Note that this dummy deliberately does NOT re-implement those rules: it delegates to
    ``PESExplorerAdapter.__init__`` exactly as a real adapter does. A test double that carries its
    own copy of the logic would only ever prove that the copy works.
    """

    def __init__(self,
                 seed_species,
                 output_directory: str,
                 network_path: str = None,
                 method: str = None,
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
        self.network_path = network_path
        self.method = method
        self.bath_gas = bath_gas
        self.explore_tol = explore_tol
        self.energy_tol = energy_tol
        self.flux_tol = flux_tol
        self.maximum_radical_electrons = maximum_radical_electrons
        self.logger = logger
        # NOT re-assigned from the raw argument here: super().__init__ above already stored the
        # NORMALIZED tuple, and re-assigning would undo that normalization, leaving this "well-behaved"
        # dummy holding a list where the real adapter holds a tuple -- i.e. quietly not well-behaved.
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
        adapter = explorer_factory(explorer='DummyExplorer', network_path='network.py', method='CSE',
                                   seed_species=['OH'],
                                   output_directory='out',
                                   )
        assert isinstance(adapter, _DummyPESExplorerAdapter)

    def test_factory_lookup_is_case_sensitive(self):
        """explorer_factory lookup is case-sensitive; a wrong-case name raises ValueError."""
        with pytest.raises(ValueError):
            explorer_factory(explorer='dummyexplorer', network_path='network.py', method='CSE',
                             seed_species=['OH'],
                             output_directory='out',
                             )

    def test_factory_unknown_name_raises_value_error_listing_keys(self):
        """The ValueError message for an unknown explorer name mentions the registered keys."""
        with pytest.raises(ValueError) as excinfo:
            explorer_factory(explorer='NoSuchExplorer', network_path='network.py', method='CSE',
                             seed_species=['OH'],
                             output_directory='out',
                             )
        assert 'DummyExplorer' in str(excinfo.value)
        assert 'DummyTSExplorer' in str(excinfo.value)

    def test_factory_routed_invalid_seed_raises(self):
        """Too many source species, routed through the factory, raises ValueError (adapter-enforced)."""
        with pytest.raises(ValueError):
            explorer_factory(explorer='DummyExplorer', network_path='network.py', method='CSE',
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
        adapter = explorer_factory(explorer='DummyTSExplorer', network_path='network.py', method='CSE',
                                   seed_species=['OH'],
                                   output_directory='out',
                                   transition_state_seeds=('TS1',),
                                   )
        assert isinstance(adapter, _DummyTSSeedCompatibleAdapter)

    def test_get_networks_before_explore_raises_runtime_error(self):
        """get_networks() before a successful explore() raises RuntimeError, never returns empty results silently."""
        adapter = explorer_factory(explorer='DummyExplorer', network_path='network.py', method='CSE', seed_species=['OH'], output_directory='out')
        with pytest.raises(RuntimeError):
            adapter.get_networks()

    def test_get_k_tp_before_explore_raises_runtime_error(self):
        """get_k_tp() before a successful explore() raises RuntimeError, never returns empty results silently."""
        adapter = explorer_factory(explorer='DummyExplorer', network_path='network.py', method='CSE', seed_species=['OH'], output_directory='out')
        with pytest.raises(RuntimeError):
            adapter.get_k_tp()

    def test_get_networks_and_get_k_tp_after_explore_succeed(self):
        """get_networks() and get_k_tp() work after a successful explore()."""
        adapter = explorer_factory(explorer='DummyExplorer', network_path='network.py', method='CSE', seed_species=['OH'], output_directory='out')
        assert adapter.explore() is True
        assert adapter.get_networks() == tuple()
        assert adapter.get_k_tp() == dict()


class TestTheFactoryCanActuallyBuildTheRegisteredArkaneAdapter:
    """
    ``explorer_factory('Arkane', ...)`` has to be able to construct the adapter it advertises.

    ``ArkaneExplorerAdapter`` registers itself under 'Arkane' and requires ``network_path`` and
    ``method``; the factory's signature carried neither, so the one registration that exists could
    never be built through the one entry point the design says to use -- every call raised
    ``TypeError: __init__() missing 2 required positional arguments``. The registry advertised a
    capability that did not exist, which is worse than not registering it, because the failure only
    appears at the call site.

    ``t3/pdep/mesolver/factory.py::mesolver_factory`` already takes exactly these two arguments, so
    this is restoring an existing idiom rather than inventing one.
    """

    def test_the_arkane_adapter_is_constructible_through_the_factory(self, tmp_path):
        """The registered adapter must be reachable through the factory, not just by direct import."""
        from t3.pdep.explorer.arkane import ArkaneExplorerAdapter

        adapter = explorer_factory(explorer='Arkane',
                                   seed_species=['methoxy'],
                                   output_directory=str(tmp_path / 'run'),
                                   network_path=str(tmp_path / 'network.py'),
                                   method='CSE',
                                   bath_gas={'He': 1.0})
        assert isinstance(adapter, ArkaneExplorerAdapter)
        assert adapter.network_path == str(tmp_path / 'network.py')
        assert adapter.method == 'CSE'

    def test_the_factory_forwards_every_argument_the_adapter_declares(self, tmp_path):
        """
        Forwarding is asserted argument by argument, not just "it constructed".

        A factory that quietly drops an argument is the failure this class exists to catch: the
        adapter would run with Arkane's defaults instead of the caller's tolerances, produce a
        plausible network, and nothing would ever indicate that what was requested was not what ran.
        """
        adapter = explorer_factory(explorer='Arkane',
                                   seed_species=['methoxy'],
                                   output_directory=str(tmp_path / 'run'),
                                   network_path=str(tmp_path / 'network.py'),
                                   method='MSC',
                                   bath_gas={'He': 1.0},
                                   explore_tol=0.02,
                                   energy_tol=45.0,
                                   flux_tol=1e-8,
                                   maximum_radical_electrons=2,
                                   database_kwargs={'thermoLibraries': ['primaryThermoLibrary']})
        assert adapter.explore_tol == 0.02
        assert adapter.energy_tol == 45.0
        assert adapter.flux_tol == 1e-8
        assert adapter.maximum_radical_electrons == 2
        assert adapter.database_kwargs == {'thermoLibraries': ['primaryThermoLibrary']}
        assert adapter.bath_gas == {'He': 1.0}


class TestRegistrationCannotBlessAnAdapterThatSkipsSeedValidation:
    """
    ``issubclass`` is a claim about ancestry, not about whether ``__init__`` reached the base class.

    The seed and capability rules (``max_source_species``, ``supports_transition_state_seeds``) live
    in ``PESExplorerAdapter.__init__`` deliberately, so that direct construction is as safe as going
    through the factory. But a subclass that overrides ``__init__`` and forgets
    ``super().__init__(...)`` inherits those rules and enforces none of them -- and
    ``register_explorer_adapter`` blesses it anyway, because ``issubclass`` is still True. The result
    is an adapter that passes registration, passes construction, and then hands Arkane a three-species
    source channel that Arkane reports as "Reactant channels with 3 reactants are not supported".

    A malicious subclass cannot be defended against and this does not try to. The target is the
    accidental one: an author of a future adapter who overrides ``__init__`` and forgets the super
    call. The factory therefore re-asserts the invariants at the entry point, against the arguments it
    was given, using the shared rule function rather than anything the subclass could override.
    """

    @staticmethod
    def _register_forgetful(monkeypatch):
        """Register an adapter whose __init__ never reaches PESExplorerAdapter.__init__."""
        class _ForgetfulAdapter(_DummyPESExplorerAdapter):
            """A subclass that overrides __init__ and forgets super() -- the accidental case."""

            def __init__(self, seed_species, output_directory, **kwargs):
                # Deliberately no super().__init__(...): this is the defect under test.
                self.seed_species = tuple(seed_species)
                self.transition_state_seeds = kwargs.get('transition_state_seeds')
                self.output_directory = output_directory
                self.explored = False

        monkeypatch.setitem(_registered_explorer_adapters, 'ForgetfulExplorer', _ForgetfulAdapter)
        return _ForgetfulAdapter

    def test_a_seed_too_large_is_refused_even_when_the_adapter_never_validates_it(self, monkeypatch, tmp_path):
        """Three source species must be refused through the factory even with the rules bypassed."""
        forgetful = self._register_forgetful(monkeypatch)
        # The premise of this test, asserted rather than assumed: direct construction really does skip
        # the rules. If this ever stops being true the test below would pass for the wrong reason.
        bypassed = forgetful(seed_species=['A', 'B', 'C'], output_directory='out')
        assert bypassed.seed_species == ('A', 'B', 'C'), 'premise: the adapter really does skip validation'

        with pytest.raises(ValueError, match='most'):
            explorer_factory(explorer='ForgetfulExplorer',
                             seed_species=['A', 'B', 'C'],
                             output_directory=str(tmp_path / 'run'),
                             network_path=str(tmp_path / 'network.py'),
                             method='CSE')

    def test_an_empty_seed_is_refused_even_when_the_adapter_never_validates_it(self, monkeypatch, tmp_path):
        """An empty seed is the other half of the rule, and it is not implied by the too-large half."""
        self._register_forgetful(monkeypatch)
        with pytest.raises(ValueError, match='at least one'):
            explorer_factory(explorer='ForgetfulExplorer',
                             seed_species=[],
                             output_directory=str(tmp_path / 'run'),
                             network_path=str(tmp_path / 'network.py'),
                             method='CSE')

    def test_transition_state_seeds_are_refused_for_an_adapter_that_does_not_support_them(
            self, monkeypatch, tmp_path):
        """The capability rule is enforced at the entry point too, not only the seed-count rule."""
        self._register_forgetful(monkeypatch)
        with pytest.raises(ValueError, match='transition'):
            explorer_factory(explorer='ForgetfulExplorer',
                             seed_species=['OH'],
                             output_directory=str(tmp_path / 'run'),
                             network_path=str(tmp_path / 'network.py'),
                             method='CSE',
                             transition_state_seeds=('TS1',))

    def test_an_adapter_that_does_not_record_the_seed_it_was_given_is_refused(self, monkeypatch, tmp_path):
        """
        A VALID seed that the adapter fails to store is the quiet half of this defect.

        Re-running the rules only catches a forgetful adapter when the seed happens to violate one. If
        the seed is valid and the adapter still never stored it, every rule passes and the object is
        returned in a state where ``self.seed_species`` does not exist -- surfacing much later as an
        AttributeError inside explore(), far from the cause. So the factory also checks the outcome,
        not just the rules.
        """
        class _NonRecordingAdapter(_DummyPESExplorerAdapter):
            """Overrides __init__, validates nothing, records nothing."""

            def __init__(self, seed_species, output_directory, **kwargs):
                self.output_directory = output_directory
                self.explored = False

        monkeypatch.setitem(_registered_explorer_adapters, 'NonRecordingExplorer', _NonRecordingAdapter)
        with pytest.raises(TypeError, match='seed'):
            explorer_factory(explorer='NonRecordingExplorer',
                             seed_species=['OH'],
                             output_directory=str(tmp_path / 'run'),
                             network_path=str(tmp_path / 'network.py'),
                             method='CSE')

    def test_an_adapter_that_records_the_seed_but_drops_the_ts_seeds_is_refused(self, monkeypatch, tmp_path):
        """
        The seed is not the only thing a forgetful ``__init__`` can drop, and checking it alone made
        the post-condition an enumeration of one field rather than a statement about the object.

        This adapter is the awkward case the previous check waved through: it records
        ``seed_species`` faithfully, so every rule passes and the recorded-seed check passes, and it
        silently discards the transition-state seeds. The exploration then runs -- successfully, and
        reporting success -- as an ordinary unseeded exploration, which is a different and cheaper
        calculation than the one that was asked for. Nothing raises, so the result is a network the
        caller believes was TS-seeded and was not.
        """
        class _DropsTSSeedsAdapter(_DummyTSSeedCompatibleAdapter):
            """Records the seed, drops the transition-state seeds."""

            def __init__(self, seed_species, output_directory, **kwargs):
                self.seed_species = tuple(seed_species)
                self.transition_state_seeds = tuple()
                self.output_directory = output_directory
                self.explored = False

        monkeypatch.setitem(_registered_explorer_adapters, 'DropsTSSeedsExplorer', _DropsTSSeedsAdapter)
        with pytest.raises(TypeError, match='transition_state_seeds'):
            explorer_factory(explorer='DropsTSSeedsExplorer',
                             seed_species=['OH'],
                             output_directory=str(tmp_path / 'run'),
                             network_path=str(tmp_path / 'network.py'),
                             method='CSE',
                             transition_state_seeds=('TS1',))

    def test_a_list_seed_and_list_ts_seeds_satisfy_the_post_condition_and_are_stored_as_tuples(
            self, monkeypatch, tmp_path):
        """
        Over-refusal guard for the widened post-condition: passing lists must not trip it.

        This does NOT pin "compares the normalized return rather than the raw argument", which was the
        claim it was first written to make. That distinction is not observable: normalization is
        ``tuple(x or tuple())`` (adapter.py:50-51), so a version of the check that re-derived the
        expected values from the raw arguments would agree in every case, and a mutation doing exactly
        that survives this whole suite. Using the return value is a single-source-of-truth choice
        argued in the comment at the call site, not a behaviour a test can distinguish today -- and
        saying so here is cheaper than leaving a future reader to rediscover it via the same mutation.
        """
        monkeypatch.setitem(_registered_explorer_adapters, 'TSDummyExplorer', _DummyTSSeedCompatibleAdapter)
        adapter = explorer_factory(explorer='TSDummyExplorer',
                                   seed_species=['OH'],
                                   output_directory=str(tmp_path / 'run'),
                                   network_path=str(tmp_path / 'network.py'),
                                   method='CSE',
                                   transition_state_seeds=['TS1'])
        assert adapter.seed_species == ('OH',)
        assert adapter.transition_state_seeds == ('TS1',)

    def test_a_one_shot_iterable_seed_survives_being_validated_before_construction(self, monkeypatch, tmp_path):
        """
        The factory validates the seed and then constructs the adapter, so it hands the seed onward
        TWICE -- and a one-shot iterable does not survive that.

        ``validate_explorer_seed`` consumes whatever it is given into tuples (adapter.py:50-51). If
        the factory then forwards the ORIGINAL object, a generator arrives at the adapter exhausted,
        the adapter's own ``super().__init__`` re-validation sees an empty sequence, and the caller is
        told 'A PES exploration requires at least one source species' about a seed they demonstrably
        supplied. The message sends them looking for a missing seed instead of a consumed one.

        Forwarding the normalized tuples fixes it at the cause. This is not merely a nicer message:
        the two calls must agree on WHAT was validated, and they cannot if the first one destroys it.
        """
        monkeypatch.setitem(_registered_explorer_adapters, 'DummyExplorer', _DummyPESExplorerAdapter)
        adapter = explorer_factory(explorer='DummyExplorer',
                                   seed_species=(label for label in ['OH']),
                                   output_directory=str(tmp_path / 'run'),
                                   network_path=str(tmp_path / 'network.py'),
                                   method='CSE')
        assert adapter.seed_species == ('OH',)

    def test_the_adapter_receives_the_normalized_seed_not_the_raw_argument(self, monkeypatch, tmp_path):
        """
        The adapter must be constructed with the tuples validation produced, not the caller's object.

        Pinned separately from the generator case because it is the general rule and the generator is
        only its most visible symptom: the value the factory VALIDATED and the value it CONSTRUCTS
        WITH have to be the same value, or the validation vouched for something else.
        """
        seen = {}

        class _RecordsWhatItWasGiven(_DummyTSSeedCompatibleAdapter):
            def __init__(self, seed_species, output_directory, **kwargs):
                seen['seed_species'] = seed_species
                seen['transition_state_seeds'] = kwargs.get('transition_state_seeds')
                super().__init__(seed_species=seed_species, output_directory=output_directory, **kwargs)

        monkeypatch.setitem(_registered_explorer_adapters, 'RecordingExplorer', _RecordsWhatItWasGiven)
        explorer_factory(explorer='RecordingExplorer',
                         seed_species=['OH'],
                         output_directory=str(tmp_path / 'run'),
                         network_path=str(tmp_path / 'network.py'),
                         method='CSE',
                         transition_state_seeds=['TS1'])
        assert seen['seed_species'] == ('OH',), 'the adapter got the raw list, not the normalized tuple'
        assert seen['transition_state_seeds'] == ('TS1',)

    def test_a_well_behaved_adapter_is_not_refused_by_the_entry_point_check(self, monkeypatch, tmp_path):
        """Over-refusal guard: the normal adapter, which does call super(), must be unaffected."""
        monkeypatch.setitem(_registered_explorer_adapters, 'DummyExplorer', _DummyPESExplorerAdapter)
        adapter = explorer_factory(explorer='DummyExplorer',
                                   seed_species=['OH'],
                                   output_directory=str(tmp_path / 'run'),
                                   network_path=str(tmp_path / 'network.py'),
                                   method='CSE')
        assert adapter.seed_species == ('OH',)


class TestASeedGivenAsABareStringIsRefusedAtTheAdapterLayerToo:
    """
    The same ``tuple('OH') == ('O', 'H')`` trap in the shared seed rule, which the writer's copy
    cannot cover: an adapter constructed directly never reaches
    ``write_arkane_explorer_input_file``'s own check, and ``max_source_species`` is 2, so a
    two-letter string passes the count rule as a bimolecular seed.
    """

    def test_a_bare_string_seed_is_refused_by_the_shared_rule(self):
        """Refused where the rule lives, so both the adapter and the factory inherit it."""
        with pytest.raises(ValueError, match='seed_species'):
            validate_explorer_seed(seed_species='OH', transition_state_seeds=None,
                                   max_source_species=2, supports_transition_state_seeds=False)

    def test_a_bare_string_seed_is_refused_through_the_factory(self, monkeypatch, tmp_path):
        """And end-to-end, since that is the path a caller actually uses."""
        monkeypatch.setitem(_registered_explorer_adapters, 'DummyExplorer', _DummyPESExplorerAdapter)
        with pytest.raises(ValueError, match='seed_species'):
            explorer_factory(explorer='DummyExplorer', seed_species='OH',
                             output_directory=str(tmp_path / 'run'),
                             network_path='network.py', method='CSE')

    def test_a_bare_string_seed_is_refused_by_direct_construction(self):
        """The adapter layer is the other entry path, and it must not be weaker."""
        with pytest.raises(ValueError, match='seed_species'):
            _DummyPESExplorerAdapter(seed_species='OH', output_directory='out')

    @pytest.mark.parametrize('seed_species', [['OH', None], ['OH', ''], [1], [['nested']], [b'OH']])
    def test_a_non_string_seed_ELEMENT_is_refused_by_the_shared_rule(self, seed_species):
        """
        The element check needs its own tests at THIS layer, not just at the writer's.

        Deleting the per-element check here left the suite green: the writer
        (``t3.pdep.explorer.input_file``) has its own copy with its own tests, and those were doing all
        the work. But an adapter can be constructed and its seed stored without the writer ever being
        called, so the shared rule's own element check has to be pinned where it lives.
        """
        with pytest.raises(ValueError, match='seed_species'):
            validate_explorer_seed(seed_species=seed_species, transition_state_seeds=None,
                                   max_source_species=2, supports_transition_state_seeds=False)

    def test_a_non_string_transition_state_seed_element_is_refused(self):
        """The same rule applies to the TS seeds, which are a separate sequence."""
        with pytest.raises(ValueError, match='transition_state_seeds'):
            validate_explorer_seed(seed_species=['OH'], transition_state_seeds=[None],
                                   max_source_species=2, supports_transition_state_seeds=True)

    def test_a_proper_sequence_is_still_accepted_by_the_shared_rule(self):
        """Over-refusal guard, including the tuple and list spellings."""
        assert validate_explorer_seed(seed_species=['OH'], transition_state_seeds=None,
                                      max_source_species=2,
                                      supports_transition_state_seeds=False)[0] == ('OH',)
        assert validate_explorer_seed(seed_species=('OH', 'CH3'), transition_state_seeds=None,
                                      max_source_species=2,
                                      supports_transition_state_seeds=False)[0] == ('OH', 'CH3')
