#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_explorer_config module

Tests for ``t3.pdep.explorer.config.PDepExplorerConfig``: a frozen, validated description of one
PES-exploration request. Every ``__post_init__`` guard listed in the class docstring gets a
dedicated refusal test, plus tests for the coercions (``tuple``, integral-float, deep-freeze)
this config performs and for a fully valid round trip.
"""

import copy

import pytest

from t3.pdep.explorer.config import PDepExplorerConfig, deep_thaw


def _valid_kwargs(**overrides) -> dict:
    """
    A minimal set of constructor kwargs that passes every ``__post_init__`` guard.

    Args:
        **overrides: Any keyword to override in the returned dict.

    Returns:
        dict: Kwargs suitable for ``PDepExplorerConfig(**kwargs)``.
    """
    kwargs = dict(
        explorer='Arkane',
        trusted_output_root='/tmp/root',
        output_directory='/tmp/root/sub',
        seed_species=['methoxy'],
        method='CSE',
    )
    kwargs.update(overrides)
    return kwargs


def test_refuses_empty_explorer_name():
    """An empty (or non-str) 'explorer' would be an unusable registry lookup key downstream."""
    with pytest.raises(ValueError, match="'explorer' must be a non-empty str"):
        PDepExplorerConfig(**_valid_kwargs(explorer=''))


def test_refuses_non_str_explorer_name():
    """A non-str 'explorer' (e.g. None) would fail the registry lookup with a confusing error."""
    with pytest.raises(ValueError, match="'explorer' must be a non-empty str"):
        PDepExplorerConfig(**_valid_kwargs(explorer=None))


def test_refuses_empty_method():
    """An empty 'method' cannot be looked up in METHOD_MAP."""
    with pytest.raises(ValueError, match="'method' must be a non-empty str"):
        PDepExplorerConfig(**_valid_kwargs(method=''))


def test_refuses_method_not_in_method_map():
    """A method outside {'CSE', 'MSC', 'RS'} cannot be rendered by write_arkane_explorer_input_file."""
    with pytest.raises(ValueError, match="'method' must be one of"):
        PDepExplorerConfig(**_valid_kwargs(method='bogus'))


def test_refuses_relative_trusted_output_root():
    """A relative root cannot be safely compared against a resolved output_directory."""
    with pytest.raises(ValueError, match="'trusted_output_root' must be an absolute path"):
        PDepExplorerConfig(**_valid_kwargs(trusted_output_root='relative/root'))


def test_refuses_relative_output_directory():
    """A relative output_directory cannot be safely compared against a resolved root."""
    with pytest.raises(ValueError, match="'output_directory' must be an absolute path"):
        PDepExplorerConfig(**_valid_kwargs(output_directory='relative/sub'))


def test_refuses_empty_trusted_output_root():
    """An empty string is not a usable path and must be refused before the isabs check."""
    with pytest.raises(ValueError, match="'trusted_output_root' must be a non-empty str"):
        PDepExplorerConfig(**_valid_kwargs(trusted_output_root=''))


def test_refuses_empty_output_directory():
    """An empty string is not a usable path and must be refused before the isabs check."""
    with pytest.raises(ValueError, match="'output_directory' must be a non-empty str"):
        PDepExplorerConfig(**_valid_kwargs(output_directory=''))


def test_refuses_output_directory_equal_to_trusted_output_root():
    """output_directory == trusted_output_root must be refused: it is the root, not inside it."""
    with pytest.raises(ValueError, match='must not be the root itself'):
        PDepExplorerConfig(**_valid_kwargs(trusted_output_root='/tmp/root', output_directory='/tmp/root'))


def test_refuses_sibling_directory_of_trusted_output_root():
    """
    A sibling directory sharing a string prefix with the root must be refused.

    '/tmp/rootX' shares the string prefix '/tmp/root' with trusted_output_root='/tmp/root', which is
    exactly the case a naive str.startswith() containment check would WRONGLY accept. This test
    catches a regression from os.path.realpath+commonpath back to str.startswith.
    """
    with pytest.raises(ValueError, match='resolve to a location strictly inside'):
        PDepExplorerConfig(**_valid_kwargs(trusted_output_root='/tmp/root', output_directory='/tmp/rootX'))


def test_refuses_output_directory_outside_trusted_output_root():
    """An output_directory that is not even a sibling, just unrelated, must be refused."""
    with pytest.raises(ValueError, match='resolve to a location strictly inside'):
        PDepExplorerConfig(**_valid_kwargs(trusted_output_root='/tmp/root', output_directory='/var/tmp/other'))


def test_refuses_non_mapping_bath_gas():
    """A non-Mapping bath_gas (e.g. a list of pairs) would fail later with a confusing AttributeError."""
    with pytest.raises(ValueError, match="'bath_gas' must be a Mapping"):
        PDepExplorerConfig(**_valid_kwargs(bath_gas=[('He', 1.0)]))


def test_refuses_non_mapping_database_kwargs():
    """A non-Mapping database_kwargs (e.g. a list) would fail later with a confusing AttributeError."""
    with pytest.raises(ValueError, match="'database_kwargs' must be a Mapping"):
        PDepExplorerConfig(**_valid_kwargs(database_kwargs=['thermoLibraries']))


def test_refuses_non_finite_numeric_field():
    """Numeric fields are validated by delegating to validate_explorer_field_values."""
    with pytest.raises(ValueError, match='finite'):
        PDepExplorerConfig(**_valid_kwargs(energy_tol=float('inf')))


def test_refuses_negative_numeric_field():
    """Numeric fields are validated by delegating to validate_explorer_field_values."""
    with pytest.raises(ValueError, match="'flux_tol' must be"):
        PDepExplorerConfig(**_valid_kwargs(flux_tol=-1.0))


def test_integral_float_maximum_radical_electrons_is_coerced_to_int():
    """
    An integral float count must land on the attribute AS an int, not merely pass validation.

    validate_explorer_field_values returns the COERCED value; __post_init__ must write it back via
    object.__setattr__, not just call the function for its validation side effect and discard the
    return value.
    """
    config = PDepExplorerConfig(**_valid_kwargs(maximum_radical_electrons=2.0))
    assert config.maximum_radical_electrons == 2
    assert isinstance(config.maximum_radical_electrons, int)
    assert not isinstance(config.maximum_radical_electrons, bool)


def test_seed_species_is_coerced_to_tuple():
    """A list passed for seed_species must be stored as a tuple, matching the declared field type."""
    config = PDepExplorerConfig(**_valid_kwargs(seed_species=['methoxy', 'CH2O']))
    assert config.seed_species == ('methoxy', 'CH2O')
    assert isinstance(config.seed_species, tuple)


def test_transition_state_seeds_is_coerced_to_tuple():
    """A list passed for transition_state_seeds must be stored as a tuple."""
    config = PDepExplorerConfig(**_valid_kwargs(transition_state_seeds=['TS1', 'TS2']))
    assert config.transition_state_seeds == ('TS1', 'TS2')
    assert isinstance(config.transition_state_seeds, tuple)


def test_mutating_original_bath_gas_dict_after_construction_does_not_affect_config():
    """
    The config must be immune to later mutation of the dict the caller passed in.

    A shallow dict(value) copy alone would still share any nested mutable value; this bath_gas is
    flat, so the primary defect this test targets is aliasing of the top-level dict itself (removing
    the deep-copy-and-view step entirely).
    """
    original = {'He': 1.0}
    config = PDepExplorerConfig(**_valid_kwargs(bath_gas=original))
    original['He'] = 0.5
    original['Ar'] = 0.5
    assert config.bath_gas == {'He': 1.0}


def test_mutating_original_database_kwargs_dict_after_construction_does_not_affect_config():
    """The config must be immune to later mutation of the dict the caller passed in."""
    original = {'thermoLibraries': ['primaryThermoLibrary']}
    config = PDepExplorerConfig(**_valid_kwargs(database_kwargs=original))
    original['thermoLibraries'].append('injected')
    original['kineticsFamilies'] = 'all'
    assert config.database_kwargs == {'thermoLibraries': ('primaryThermoLibrary',)}


def test_stored_bath_gas_cannot_be_mutated_through_the_config_attribute():
    """The stored bath_gas is a read-only view (MappingProxyType), not a mutable dict."""
    config = PDepExplorerConfig(**_valid_kwargs(bath_gas={'He': 1.0}))
    with pytest.raises(TypeError):
        config.bath_gas['He'] = 0.5


def test_stored_database_kwargs_nested_values_cannot_be_mutated_through_the_config_attribute():
    """
    The freeze must be recursive, not just the outer mapping.

    A MappingProxyType over a deep copy already decouples the config from the caller's dict, but a
    nested LIST value fished out of the proxy would still be the config's own mutable list --
    ``config.database_kwargs['thermoLibraries'].append(...)`` would mutate a supposedly frozen
    config after validation, changing the Arkane input eventually generated from it. Deep freezing
    stores that value as a tuple, so the append is impossible.
    """
    config = PDepExplorerConfig(**_valid_kwargs(database_kwargs={'thermoLibraries': ['primaryThermoLibrary']}))
    assert isinstance(config.database_kwargs['thermoLibraries'], tuple)
    with pytest.raises(AttributeError):
        config.database_kwargs['thermoLibraries'].append('injected')


def test_deep_thaw_restores_plain_mutable_structures():
    """
    ``deep_thaw`` must hand the consumption boundary plain dicts/lists, freshly copied.

    The input writer's ``_validate_database_kwarg`` requires a real list for the ``database(...)``
    library keywords, so the frozen tuples must come back as lists -- and mutating what deep_thaw
    returned must not reach the frozen config it came from.
    """
    config = PDepExplorerConfig(**_valid_kwargs(database_kwargs={'thermoLibraries': ['primaryThermoLibrary']}))
    thawed = deep_thaw(config.database_kwargs)
    assert thawed == {'thermoLibraries': ['primaryThermoLibrary']}
    assert isinstance(thawed, dict)
    assert isinstance(thawed['thermoLibraries'], list)
    thawed['thermoLibraries'].append('injected')
    assert config.database_kwargs == {'thermoLibraries': ('primaryThermoLibrary',)}
    assert deep_thaw(None) is None


def test_valid_config_round_trips_all_constructor_values():
    """A fully valid, well-formed config must store every constructor value correctly."""
    seed_species = ['methoxy']
    transition_state_seeds = ['TS1']
    bath_gas = {'He': 1.0}
    database_kwargs = {'thermoLibraries': ['primaryThermoLibrary']}
    config = PDepExplorerConfig(
        explorer='Arkane',
        trusted_output_root='/tmp/root',
        output_directory='/tmp/root/sub/exploration',
        seed_species=seed_species,
        method='MSC',
        bath_gas=bath_gas,
        explore_tol=0.5,
        energy_tol=64.0,
        flux_tol=1e-15,
        maximum_radical_electrons=1,
        transition_state_seeds=transition_state_seeds,
        database_kwargs=database_kwargs,
    )
    assert config.explorer == 'Arkane'
    assert config.trusted_output_root == '/tmp/root'
    assert config.output_directory == '/tmp/root/sub/exploration'
    assert config.seed_species == ('methoxy',)
    assert config.method == 'MSC'
    assert config.bath_gas == bath_gas
    assert config.explore_tol == 0.5
    assert config.energy_tol == 64.0
    assert config.flux_tol == 1e-15
    assert config.maximum_radical_electrons == 1
    assert config.transition_state_seeds == ('TS1',)
    # Deeply frozen: the nested list value is stored as a tuple.
    assert config.database_kwargs == {'thermoLibraries': ('primaryThermoLibrary',)}
    # Rebuilt, so the returned collections must be equal in content, not identical, to the originals.
    assert config.bath_gas is not bath_gas
    assert config.database_kwargs is not database_kwargs


@pytest.mark.parametrize('field_name', ['seed_species', 'transition_state_seeds'])
def test_a_bare_string_seed_is_refused_here_and_not_silently_tupled(field_name):
    """
    A bare string seed must be refused BY THE CONFIG, not merely by the factory later on.

    This is not a duplicate of the factory's identical refusal, and moving it there would not do.
    The config normalizes its seed fields to tuples, and ``tuple('OH')`` is ``('O', 'H')`` -- so
    without this guard the config would hand the factory a perfectly well-formed two-species
    bimolecular source channel, and ``validate_explorer_seed``'s bare-string refusal could never
    fire, because by then nothing distinguishes the typo from a caller who genuinely meant O + H.
    The guard has to run before whichever ``tuple()`` call comes first, which is this one.
    """
    with pytest.raises(ValueError) as exc_info:
        PDepExplorerConfig(**_valid_kwargs(**{field_name: 'OH'}))
    message = str(exc_info.value)
    assert field_name in message
    assert 'bare str' in message
    # The message must show what the silent reading WOULD have been, so the reader sees the hazard
    # rather than just being told 'no'.
    assert "('O', 'H')" in message


def test_a_bytes_seed_is_refused_for_the_same_reason():
    """
    bytes is a sequence too, and iterating it yields ints -- an even more confusing failure.

    Would catch a guard narrowed to ``str`` alone.
    """
    with pytest.raises(ValueError, match='bare bytes'):
        PDepExplorerConfig(**_valid_kwargs(seed_species=b'OH'))


def test_a_genuine_two_species_seed_is_still_accepted():
    """
    The refusal must reject the SPELLING, not the bimolecular source channel itself.

    Would catch a guard that over-refuses by rejecting any multi-entry seed, which is the
    over-refusal failure mode this package has hit three times: refusing one spelling while the
    thing supposedly being prevented remains fully reachable one keystroke away.
    """
    config = PDepExplorerConfig(**_valid_kwargs(seed_species=['O', 'H']))
    assert config.seed_species == ('O', 'H')
