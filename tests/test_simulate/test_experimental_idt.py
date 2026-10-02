"""Tests for versioned, per-point experimental IDT comparisons."""

import copy
import math
import os
from pathlib import Path
from time import perf_counter

import cantera as ct
import numpy as np
import pytest
import yaml
from pydantic import ValidationError

from t3.schema import ExperimentalIDTFile
from t3.simulate.cantera_idt import (CanteraIDT, compute_source_defined_idt,
                                     source_defined_idt_is_resolved)
from t3.common import SIMULATE_TEST_DATA_BASE_PATH, convert_time_to_seconds


TEST_DIR_IDT = os.path.join(SIMULATE_TEST_DATA_BASE_PATH, 'cantera_idt_test')


def _point() -> dict:
    """Return a minimal valid version-1 experimental point."""
    return {
        'temperature': {'value': 1200.0, 'units': 'K'},
        'pressure': {'value': 10.0, 'units': 'bar'},
        'composition': [
            {'smiles': 'C', 'mole_fraction': 0.095057},
            {'smiles': '[O][O]', 'mole_fraction': 0.190114},
            {'smiles': 'N#N', 'mole_fraction': 0.714829},
        ],
        'apparatus': 'shock tube',
        'ignition_definition': {'target': 'temperature', 'type': 'd/dt max'},
        'idt': {'value': 2.5, 'units': 'ms'},
        'uncertainty': {'value': 0.2, 'units': 'ms'},
        'source': {'doi': '10.0000/example', 'record': 'Table 2, point 17'},
    }


def _adapter() -> CanteraIDT:
    """Construct the reduced methane adapter used by per-point integration tests."""
    model = os.path.join(TEST_DIR_IDT, 'iteration_4', 'RMG', 'cantera_from_ck', 'chem_annotated.yaml')
    rmg = {'species': [
        {'label': 'methane', 'smiles': 'C', 'concentration': 0, 'role': 'fuel',
         'equivalence_ratios': [1.0]},
        {'label': 'O2', 'smiles': '[O][O]', 'concentration': 0, 'role': 'oxidizer'},
        {'label': 'N2', 'smiles': 'N#N', 'concentration': 0, 'role': 'diluent'},
    ]}
    return CanteraIDT(t3={}, rmg=rmg, paths={'cantera annotated': model, 'figs': None}, logger=None)


def _rcm_history_point() -> dict:
    """Use the repository methane mechanism with a synthetic compression stroke."""
    point = _point()
    point.update({
        'apparatus': 'rapid compression machine',
        'initial_temperature': {'value': 1100.0, 'units': 'K'},
        'initial_pressure': {'value': 2.0, 'units': 'bar'},
        'temperature': {'value': 1400.0, 'units': 'K'},
        'pressure': {'value': 5.0, 'units': 'bar'},
        'idt': {'value': 10.0, 'units': 'ms'},
        'volume_history': {
            'time': {'values': [0.0, 1.0, 2.0, 20.0], 'units': 'ms'},
            'volume': {'values': [100.0, 75.0, 50.0, 50.0], 'units': 'cm3'},
        },
    })
    return point


def _parsed_rcm_point(point=None):
    """Validate a history-driven point for physics tests."""
    return ExperimentalIDTFile.model_validate({'version': 1, 'points': [point or _rcm_history_point()]}).points[0]


def test_rcm_volume_history_schema_accepts_units_and_initial_states():
    point = _rcm_history_point()
    point['volume_history']['compression_time'] = {'value': 2000.0, 'units': 'us'}
    parsed = _parsed_rcm_point(point)
    times, volumes = parsed.volume_history.to_si()
    assert times == [0.0, 0.001, 0.002, 0.02]
    assert volumes == pytest.approx([1e-4, 7.5e-5, 5e-5, 5e-5])
    assert parsed.volume_history.end_of_compression == 0.002
    assert parsed.initial_temperature.value == 1100.0
    assert parsed.temperature.value == 1400.0


@pytest.mark.parametrize('mutation, message', [
    ({'time': {'values': [0.0], 'units': 's'}, 'volume': {'values': [1.0], 'units': 'm3'}}, 'at least two'),
    ({'time': {'values': [0.0, 1.0], 'units': 's'}}, 'equal length'),
    ({'time': {'values': [0.0, 1.0, 1.0, 2.0], 'units': 's'}}, 'strictly increasing'),
    ({'time': {'values': [0.0, 2.0, 1.0, 3.0], 'units': 's'}}, 'strictly increasing'),
    ({'time': {'values': [0.0, math.inf, 2.0, 3.0], 'units': 's'}}, 'finite'),
    ({'time': {'values': [0.0, math.nan, 2.0, 3.0], 'units': 's'}}, 'finite'),
    ({'volume': {'values': [1.0, 0.0, 0.5, 0.5], 'units': 'L'}}, 'greater than zero'),
    ({'volume': {'values': [1.0, -1.0, 0.5, 0.5], 'units': 'L'}}, 'greater than zero'),
    ({'volume': {'values': [1.0, math.inf, 0.5, 0.5], 'units': 'L'}}, 'finite'),
    ({'volume': {'values': [1.0, math.nan, 0.5, 0.5], 'units': 'L'}}, 'finite'),
    ({'volume': {'values': [1.0, 5e-324, 0.5, 0.5], 'units': 'cm3'}}, 'greater than zero'),
    ({'time': {'values': [0.0, 5e-324, 1.0, 2.0], 'units': 'us'}}, 'strictly increasing'),
    ({'time': {'values': [0.0, 1.0, 2.0, 3.0], 'units': 'minutes'}}, 's'),
    ({'volume': {'values': [1.0, 0.75, 0.5, 0.5], 'units': 'mm3'}}, 'm3'),
    ({'compression_time': {'value': math.inf, 'units': 's'}}, 'finite'),
    ({'compression_time': {'value': 21.0, 'units': 'ms'}}, 'within the history'),
    ({'compression_time': {'value': -1.0, 'units': 'ms'}}, 'within the history'),
    ({'time': {'values': [-1e308, 0.0, 1.0, 1e308], 'units': 's'}}, 'horizon.*finite'),
    ({'time': {'values': [-1e308, 1e308, 1e308, 1e308], 'units': 's'}}, 'strictly increasing'),
    ({'volume': {'values': [1.0, 1e308, 0.5, 0.5], 'units': 'm3'}}, 'slopes.*finite'),
])
def test_rcm_volume_history_schema_rejects_invalid_history(mutation, message):
    point = _rcm_history_point()
    point['volume_history'].update(mutation)
    with pytest.raises(ValidationError, match=message):
        _parsed_rcm_point(point)


def test_rcm_volume_history_schema_rejects_shock_tube():
    point = _rcm_history_point()
    point['apparatus'] = 'shock tube'
    with pytest.raises(ValidationError, match='only.*rapid compression machine'):
        _parsed_rcm_point(point)


@pytest.mark.parametrize('field', ['initial_temperature', 'initial_pressure'])
def test_rcm_volume_history_schema_requires_initial_state(field):
    point = _rcm_history_point()
    del point[field]
    with pytest.raises(ValidationError, match='initial_temperature and initial_pressure'):
        _parsed_rcm_point(point)


@pytest.mark.parametrize('time_units, time_factor', [('s', 1.0), ('ms', 1000.0), ('us', 1e6)])
@pytest.mark.parametrize('volume_units, volume_factor', [('m3', 1.0), ('cm3', 1e6), ('L', 1000.0)])
def test_rcm_volume_history_schema_converts_all_units(time_units, time_factor, volume_units, volume_factor):
    point = _rcm_history_point()
    point['volume_history']['time'] = {'values': [0.0, 0.002 * time_factor], 'units': time_units}
    point['volume_history']['volume'] = {'values': [1e-4 * volume_factor, 5e-5 * volume_factor], 'units': volume_units}
    times, volumes = _parsed_rcm_point(point).volume_history.to_si()
    assert times == pytest.approx([0.0, 0.002])
    assert volumes == pytest.approx([1e-4, 5e-5])


def test_rcm_volume_history_schema_compression_time_overrides_minimum():
    point = _rcm_history_point()
    point['volume_history']['compression_time'] = {'value': 1.0, 'units': 'ms'}
    assert _parsed_rcm_point(point).volume_history.end_of_compression == 0.001
    del point['volume_history']['compression_time']
    assert _parsed_rcm_point(point).volume_history.end_of_compression == 0.002


@pytest.mark.parametrize('time_values, compression_time, accepted', [
    ([2**52, 2**52 + 2, 2**52 + 4], 2**52 + 4, False),
    ([0.0, 2.0, 4.0], 4.0, True),
])
def test_rcm_volume_history_schema_rejects_unrepresentable_horizon(
        time_values, compression_time, accepted):
    point = _rcm_history_point()
    point['idt'] = {'value': 10.0, 'units': 'ms'}
    point['volume_history']['time'] = {'values': time_values, 'units': 's'}
    point['volume_history']['volume']['values'] = [100.0, 75.0, 50.0]
    point['volume_history']['compression_time'] = {'value': compression_time, 'units': 's'}
    if accepted:
        _parsed_rcm_point(point)
    else:
        with pytest.raises(ValidationError, match='horizon.*post-compression window'):
            _parsed_rcm_point(point)


@pytest.mark.parametrize('field', ['initial_temperature', 'initial_pressure'])
def test_rcm_volume_history_schema_rejects_initial_state_without_history(field):
    point = _point()
    point[field] = _rcm_history_point()[field]
    with pytest.raises(ValidationError, match='require a volume_history'):
        ExperimentalIDTFile.model_validate({'version': 1, 'points': [point]})


@pytest.mark.parametrize('pressure', [
    {'value': 1e308, 'units': 'atm'},
    {'value': 5e-324, 'units': 'Pa'},
])
def test_rcm_volume_history_schema_checks_initial_pressure_in_si(pressure):
    point = _rcm_history_point()
    point['initial_pressure'] = pressure
    with pytest.raises(ValidationError, match='initial_pressure converted to Pa'):
        _parsed_rcm_point(point)


def test_rcm_volume_history_schema_validates_existing_versioned_files():
    root = Path(__file__).resolve().parents[2]
    validated = []
    for folder in ('examples', 'tests'):
        for path in (root / folder).rglob('*.yaml'):
            with path.open() as stream:
                data = yaml.safe_load(stream)
            if isinstance(data, dict) and data.get('version') == 1 and 'points' in data:
                ExperimentalIDTFile.model_validate(data)
                validated.append(path)
    assert root / 'examples/idt_with_experiment/experimental_idt_v1.yaml' in validated


def test_versioned_example_compositions_are_configured_in_example_input():
    """Every versioned example composition SMILES is present in its input species."""
    root = Path(__file__).resolve().parents[2]
    with (root / 'examples/idt_with_experiment/experimental_idt_v1.yaml').open() as stream:
        experimental = yaml.safe_load(stream)
    with (root / 'examples/idt_with_experiment/input.yml').open() as stream:
        configured = yaml.safe_load(stream)

    configured_smiles = {species['smiles'] for species in configured['rmg']['species']}
    example_smiles = {
        species['smiles']
        for point in experimental['points']
        for species in point['composition']
    }
    assert example_smiles <= configured_smiles


def test_rcm_volume_history_nonreactive_isentropic_compression():
    adapter = _adapter()
    point = _rcm_history_point()
    point['composition'] = [{'smiles': 'N#N', 'mole_fraction': 1.0}]
    parsed = _parsed_rcm_point(point)
    mixture, _ = adapter._map_experimental_composition(parsed.composition)
    history = adapter._simulate_rcm_volume_history(parsed, mixture, reactive=False)
    gas = ct.Solution(adapter.paths['cantera annotated'])
    gas.TPX = 1100.0, 2e5, mixture
    initial_entropy, initial_density = gas.entropy_mass, gas.density
    gas.SV = initial_entropy, 1.0 / (2.0 * initial_density)
    compression_index = int(np.flatnonzero(np.isclose(history.t, 0.002, atol=1e-14))[0])
    assert history.T[compression_index] == pytest.approx(gas.T, rel=0.005)
    assert history.P[compression_index] == pytest.approx(gas.P, rel=0.005)
    assert history.volume[compression_index] == pytest.approx(5e-5, rel=1e-6)
    assert history.T[0] == 1100.0
    assert history.P[0] == pytest.approx(2e5)
    assert history.volume == pytest.approx(np.interp(history.t, [0.0, 0.001, 0.002, 0.02],
                                                    [1e-4, 7.5e-5, 5e-5, 5e-5]), rel=1e-6)
    assert history.t[-1] == pytest.approx(0.102)
    assert history.volume[-1] == pytest.approx(5e-5, rel=1e-6)


def test_rcm_volume_history_integrates_through_long_history_end():
    adapter = _adapter()
    point = _rcm_history_point()
    point['idt'] = {'value': 1.0, 'units': 'ms'}
    point['composition'] = [{'smiles': 'N#N', 'mole_fraction': 1.0}]
    parsed = _parsed_rcm_point(point)
    mixture, _ = adapter._map_experimental_composition(parsed.composition)
    history = adapter._simulate_rcm_volume_history(parsed, mixture)
    assert history.t[-1] == 0.02


def test_rcm_volume_history_reactive_ignition_is_relative_to_compression():
    adapter = _adapter()
    parsed = _parsed_rcm_point()
    mixture, _ = adapter._map_experimental_composition(parsed.composition)
    history = adapter._simulate_rcm_volume_history(parsed, mixture)
    mask = history.t >= parsed.volume_history.end_of_compression
    expected = compute_source_defined_idt(history.t[mask] - 0.002, history.T[mask], 'd/dt max')
    result = adapter._compare_versioned_experiment({'version': 1, 'points': [_rcm_history_point()]})
    assert result['n_matched'] == 1
    assert result['points'][0]['idt_sim'] == pytest.approx(expected, rel=1e-9)
    assert expected > 0
    assert np.max(history.T) > 2000


def test_rcm_volume_history_ignition_does_not_depend_on_time_origin():
    adapter = _adapter()
    point = _rcm_history_point()
    original = adapter._compare_versioned_experiment({'version': 1, 'points': [point]})
    point['volume_history']['time']['values'] = [5.0, 6.0, 7.0, 25.0]
    shifted = adapter._compare_versioned_experiment({'version': 1, 'points': [point]})
    assert original['n_matched'] == shifted['n_matched'] == 1
    assert shifted['points'][0]['idt_sim'] == pytest.approx(original['points'][0]['idt_sim'], rel=1e-6)


def test_rcm_volume_history_long_tail_does_not_hide_fast_ignition():
    adapter = _adapter()
    point = _rcm_history_point()
    original = adapter._compare_versioned_experiment({'version': 1, 'points': [point]})
    point['volume_history']['time']['values'][-1] = 2000.0
    extended = adapter._compare_versioned_experiment({'version': 1, 'points': [point]})
    assert original['n_matched'] == extended['n_matched'] == 1
    assert extended['points'][0]['idt_sim'] == pytest.approx(original['points'][0]['idt_sim'], rel=0.005)


def _dense_rcm_ignition_delay(adapter, point, event_start, event_end):
    parsed = _parsed_rcm_point(point)
    mixture, _ = adapter._map_experimental_composition(parsed.composition)
    samples = np.unique(np.concatenate((parsed.volume_history.to_si()[0],
                                        [parsed.volume_history.end_of_compression, parsed.volume_history_horizon],
                                        np.linspace(event_start, event_end, 5001))))
    history = adapter._simulate_rcm_volume_history(parsed, mixture, sample_times=samples)
    mask = (history.t >= event_start) & (history.t <= event_end)
    return compute_source_defined_idt(history.t[mask] - parsed.volume_history.end_of_compression,
                                     history.T[mask], 'd/dt max')


def test_rcm_volume_history_resolves_ignition_between_long_delay_grid_samples():
    adapter = _adapter()
    point = _rcm_history_point()
    point['initial_temperature']['value'] = 1000.0
    point['idt'] = {'value': 10.0, 'units': 's'}
    expected = _dense_rcm_ignition_delay(adapter, point, 0.002, 0.052)
    assert 0.004 < expected < 0.008
    result = adapter._compare_versioned_experiment({'version': 1, 'points': [point]})
    assert result['n_matched'] == 1
    assert result['points'][0]['idt_sim'] == pytest.approx(expected, rel=0.05)


def test_rcm_volume_history_resolves_late_ignition_between_history_grid_samples():
    adapter = _adapter()
    point = _rcm_history_point()
    point['initial_temperature']['value'] = 800.0
    point['volume_history']['compression_time'] = {'value': 0.002, 'units': 's'}
    point['volume_history']['time'] = {'values': [0.0, 0.002, 6.0, 6.001, 9.0], 'units': 's'}
    point['volume_history']['volume']['values'] = [100.0, 99.0, 99.0, 20.0, 20.0]
    expected = _dense_rcm_ignition_delay(adapter, point, 6.001, 6.041)
    late_stroke_delay = 6.001 - 0.002
    assert 0.0002 < expected - late_stroke_delay < 0.0004
    result = adapter._compare_versioned_experiment({'version': 1, 'points': [point]})
    assert result['n_matched'] == 1
    actual = result['points'][0]['idt_sim']
    assert actual == pytest.approx(expected, rel=0.05)
    assert actual - late_stroke_delay == pytest.approx(expected - late_stroke_delay, rel=0.05)


@pytest.mark.parametrize('units, factor', [('s', 1.0), ('ms', 1000.0), ('us', 1e6)])
def test_rcm_volume_history_duration_limit_is_validated_in_si(units, factor):
    point = _rcm_history_point()
    point['volume_history']['time'] = {
        'values': [factor, 1.001 * factor, 1.002 * factor, 11.000001 * factor], 'units': units,
    }
    with pytest.raises(ValidationError, match='volume_history duration.*no greater than 10 s'):
        _parsed_rcm_point(point)
    point['volume_history']['time']['values'][-1] = 11.0 * factor
    point['idt'] = {'value': 10.0, 'units': 's'}
    parsed = _parsed_rcm_point(point)
    assert parsed.volume_history.to_si()[0][-1] == 11.0
    assert parsed.volume_history_horizon == parsed.volume_history.end_of_compression + 10.0


def test_rcm_volume_history_large_history_uses_adaptive_steps(monkeypatch):
    adapter = _adapter()
    point = _rcm_history_point()
    point['composition'] = [{'smiles': 'N#N', 'mole_fraction': 1.0}]
    point['volume_history']['time'] = {'values': np.linspace(0.0, 0.02, 10000).tolist(), 'units': 's'}
    point['volume_history']['volume'] = {'values': np.linspace(100.0, 50.0, 10000).tolist(), 'units': 'cm3'}
    parsed = _parsed_rcm_point(point)
    mixture, _ = adapter._map_experimental_composition(parsed.composition)

    def reject_fixed_grid(*args, **kwargs):
        raise AssertionError('RCM integration must sample adaptive steps, not a fixed linspace grid')

    monkeypatch.setattr('t3.simulate.cantera_idt.np.linspace', reject_fixed_grid)
    started = perf_counter()
    history = adapter._simulate_rcm_volume_history(parsed, mixture)
    reference = adapter._simulate_rcm_volume_history(parsed, mixture, reactive=False, sample_times=history.t)
    elapsed = perf_counter() - started
    assert elapsed < 30.0, f'Two 10,000-point runs exceeded the 30 s budget: {elapsed:.3f} s'
    assert np.array_equal(history.t, reference.t)
    assert np.all(np.isin(parsed.volume_history.to_si()[0], history.t))
    assert history.t[-1] == parsed.volume_history_horizon


def test_rcm_volume_history_half_max_mechanical_crossing_is_not_ignition():
    adapter = _adapter()
    point = _rcm_history_point()
    point['initial_temperature']['value'] = 900.0
    point['volume_history']['compression_time'] = {'value': 2.0, 'units': 'ms'}
    point['volume_history']['time']['values'] = [0.0, 1.0, 2.0, 2.002, 2.004, 20.0]
    point['volume_history']['volume']['values'] = [100.0, 75.0, 50.0, 30.0, 50.0, 50.0]
    point['ignition_definition'] = {'target': 'pressure', 'type': '1/2 max'}
    parsed = _parsed_rcm_point(point)
    mixture, _ = adapter._map_experimental_composition(parsed.composition)
    history = adapter._simulate_rcm_volume_history(parsed, mixture)
    assert np.max(history.T) > 2000.0
    result = adapter._compare_versioned_experiment({'version': 1, 'points': [point]})
    assert result['points'][0]['idt_sim'] is None
    assert 'no ignition within horizon' in result['points'][0]['refusal']['detail']


def test_rcm_volume_history_simulation_errors_remain_typed_failures(monkeypatch):
    adapter = _adapter()

    def fail_simulation(*args, **kwargs):
        raise ValueError('RCM integration failed')

    monkeypatch.setattr(adapter, '_simulate_rcm_volume_history', fail_simulation)
    result = adapter._compare_versioned_experiment({'version': 1, 'points': [_rcm_history_point()]})
    assert result['points'][0]['refusal'] == {'reason': 'simulation failed', 'detail': 'RCM integration failed'}
    assert result['points'][0]['idt_sim'] is None


@pytest.mark.parametrize('criterion', ['d/dt max', 'd/dt max extrapolated', 'max', '1/2 max'])
def test_rcm_volume_history_chemical_guard_uses_defining_peak(criterion):
    times = np.asarray([0.0, 0.1, 0.2, 0.3, 0.4])
    target = np.asarray([1.0, 1.0, 2.0, 10.0, 8.0])
    reference = np.full(5, 1000.0)
    heated = reference + np.asarray([0.0, 0.0, 3.0, 10.0, 8.0])
    assert CanteraIDT._rcm_has_chemical_heating(times, target, heated, reference, criterion)
    assert not CanteraIDT._rcm_has_chemical_heating(times, target, reference, reference, criterion)
    assert not CanteraIDT._rcm_has_chemical_heating(times[:2], target[:2], heated[:2], reference[:2], criterion)


def test_rcm_volume_history_mechanical_pressure_peak_is_not_ignition():
    adapter = _adapter()
    point = _rcm_history_point()
    point['composition'] = [{'smiles': 'N#N', 'mole_fraction': 1.0}]
    point['ignition_definition']['target'] = 'pressure'
    point['volume_history']['time']['values'] = [0.0, 2.0, 2.4, 2.6, 3.0, 4.0, 20.0]
    point['volume_history']['volume']['values'] = [100.0, 50.0, 70.0, 55.0, 70.0, 70.0, 70.0]
    point['volume_history']['compression_time'] = {'value': 2.0, 'units': 'ms'}
    parsed = _parsed_rcm_point(point)
    mixture, _ = adapter._map_experimental_composition(parsed.composition)
    history = adapter._simulate_rcm_volume_history(parsed, mixture)
    mask = history.t >= 0.002
    mechanical_delay = compute_source_defined_idt(history.t[mask] - 0.002, history.P[mask], 'd/dt max')
    assert 0.0004 < mechanical_delay < 0.0006
    result = adapter._compare_versioned_experiment({'version': 1, 'points': [point]})
    assert result['points'][0]['idt_sim'] is None
    assert result['points'][0]['refusal']['reason'] == 'ignition not resolved'
    assert 'no ignition within horizon' in result['points'][0]['refusal']['detail']


def test_rcm_volume_history_free_result_matches_existing_path_exactly():
    adapter = _adapter()
    point = _point()
    point['apparatus'] = 'rapid compression machine'
    parsed = ExperimentalIDTFile.model_validate({'version': 1, 'points': [point]}).points[0]
    assert parsed.volume_history is None
    mixture, _ = adapter._map_experimental_composition(parsed.composition)
    history = adapter.simulate_idt_for_a_point(
        r=0, t=1200.0, p=10.0, x=mixture, phi=None, infile=adapter.paths['cantera annotated'],
        save_fig=False, max_idt=1.0, apparatus='rapid compression machine', return_time_history=True,
    )
    expected = compute_source_defined_idt(history.t, history.T, 'd/dt max')
    result = adapter._compare_versioned_experiment({'version': 1, 'points': [point]})
    assert result['points'][0]['idt_sim'] == expected


def test_versioned_experimental_idt_schema_accepts_explicit_point():
    """The versioned schema accepts a complete, normalized mole-fraction point."""
    parsed = ExperimentalIDTFile.model_validate({'version': 1, 'points': [_point()]})

    assert parsed.version == 1
    assert parsed.points[0].idt.units == 'ms'
    assert parsed.points[0].source.record == 'Table 2, point 17'


@pytest.mark.parametrize('idt', [
    {'value': 1e308, 'units': 's'},
    {'value': 11, 'units': 's'},
])
def test_versioned_experimental_idt_schema_rejects_idt_over_ten_seconds(idt):
    """The versioned schema bounds converted experimental IDTs to ten seconds."""
    point = _point()
    point['idt'] = idt

    with pytest.raises(ValidationError, match='10 s'):
        ExperimentalIDTFile.model_validate({'version': 1, 'points': [point]})


def test_versioned_experimental_idt_schema_accepts_idt_at_ten_seconds():
    """The ten-second IDT limit is inclusive and conversion respects units."""
    point = _point()
    point['idt'] = {'value': 10000, 'units': 'ms'}

    parsed = ExperimentalIDTFile.model_validate({'version': 1, 'points': [point]})

    assert parsed.points[0].idt.value == 10000


@pytest.mark.parametrize('temperature', [
    {'value': -273.15, 'units': 'degC'},
    {'value': -300, 'units': 'degC'},
    {'value': 0, 'units': 'K'},
])
def test_versioned_experimental_idt_schema_rejects_nonpositive_kelvin(temperature):
    """Temperature validation occurs after conversion to Kelvin."""
    point = _point()
    point['temperature'] = temperature

    with pytest.raises(ValidationError, match='Kelvin'):
        ExperimentalIDTFile.model_validate({'version': 1, 'points': [point]})


def test_versioned_experimental_idt_schema_accepts_zero_celsius():
    """Zero Celsius is a valid positive Kelvin temperature."""
    point = _point()
    point['temperature'] = {'value': 0, 'units': 'degC'}

    parsed = ExperimentalIDTFile.model_validate({'version': 1, 'points': [point]})

    assert parsed.points[0].temperature.value == 0


@pytest.mark.parametrize('idt', [
    {'value': 5e-324, 'units': 'us'},
    {'value': 5e-324, 'units': 'micro-s'},
    {'value': 1e-322, 'units': 'ms'},
])
def test_versioned_experimental_idt_schema_rejects_idt_underflowing_to_zero(idt):
    """A positive subnormal IDT that converts to exactly zero seconds is rejected.

    ``Field(gt=0)`` on the raw value cannot catch this: the unit factor is applied
    afterwards, and ``5e-324 * 1e-6`` underflows to ``0.0``. Such a point used to
    validate cleanly and then divide by zero at the ``math.log10(simulated_idt /
    experimental_idt)`` comparison.
    """
    # Pin the fixture to the state the guard defends against. Without these two assertions a
    # value that stays non-zero after conversion -- 1e-320 ms is 1e-323, not 0.0 -- would make
    # this test silently exercise nothing.
    assert idt['value'] > 0, 'the fixture must be positive, or Field(gt=0) is what rejects it'
    assert convert_time_to_seconds(idt['value'], idt['units']) == 0.0, \
        'the fixture must underflow to exactly zero seconds, or it tests the wrong guard'

    point = _point()
    point['idt'] = idt

    with pytest.raises(ValidationError, match='greater than 0 s'):
        ExperimentalIDTFile.model_validate({'version': 1, 'points': [point]})


def test_versioned_experimental_idt_schema_accepts_smallest_representable_idt():
    """The lower bound is on the converted value, so a tiny-but-representable IDT still validates."""
    point = _point()
    point['idt'] = {'value': 1e-3, 'units': 'us'}

    parsed = ExperimentalIDTFile.model_validate({'version': 1, 'points': [point]})

    assert parsed.points[0].idt.value == 1e-3


def test_simulate_idt_for_a_point_rejects_nonfinite_horizon():
    """The integration horizon must be finite before reactor integration starts."""
    with pytest.raises(ValueError, match='finite integration horizon'):
        _adapter().simulate_idt_for_a_point(
            r=0,
            t=1200.0,
            p=10.0,
            x={'CH4': 0.1, 'O2': 0.2, 'N2': 0.7},
            phi=None,
            infile='not-a-real-file.yaml',
            max_idt=math.inf,
        )


@pytest.mark.parametrize(
    ('mutation', 'message'),
    [
        (lambda point: point.__setitem__('temperature', {'value': 1200.0}), 'units'),
        (lambda point: point['ignition_definition'].__setitem__('target', 'CH2O'), 'target'),
        (lambda point: point['ignition_definition'].__setitem__('type', 'onset'), 'type'),
        (lambda point: point.__setitem__('composition', []), 'at least one species'),
        (lambda point: point['composition'].append(copy.deepcopy(point['composition'][0])), 'must be unique'),
        (lambda point: point['composition'][0].__setitem__('mole_fraction', 0.5), 'must sum to 1.0'),
        (lambda point: point['temperature'].__setitem__('value', math.inf), 'finite'),
        (lambda point: point['temperature'].__setitem__('value', math.nan), 'finite'),
        (lambda point: point['pressure'].__setitem__('value', math.inf), 'finite'),
        (lambda point: point['pressure'].__setitem__('value', math.nan), 'finite'),
        (lambda point: point['idt'].__setitem__('value', math.inf), 'finite'),
        (lambda point: point['idt'].__setitem__('value', math.nan), 'finite'),
        (lambda point: point['uncertainty'].__setitem__('value', math.inf), 'finite'),
        (lambda point: point['uncertainty'].__setitem__('value', math.nan), 'finite'),
        (lambda point: point['composition'][0].__setitem__('mole_fraction', math.inf), 'finite'),
        (lambda point: point['composition'][0].__setitem__('mole_fraction', math.nan), 'finite'),
    ],
)
def test_versioned_experimental_idt_schema_rejects_invalid_fields(mutation, message):
    """Missing units and unknown source criteria fail validation, not simulation."""
    point = _point()
    mutation(point)

    with pytest.raises(ValidationError, match=message):
        ExperimentalIDTFile.model_validate({'version': 1, 'points': [point]})


@pytest.mark.parametrize(
    ('criterion_type', 'expected'),
    [
        ('d/dt max', 2.0),
        ('d/dt max extrapolated', 1.5),
    ],
)
def test_pressure_trace_derivative_criteria_have_hand_known_answers(criterion_type, expected):
    """A piecewise pressure trace has max slope at 2 s and baseline tangent at 1.5 s."""
    times = [0.0, 1.0, 2.0, 3.0, 4.0]
    pressure = [1.0, 1.0, 2.0, 5.0, 5.0]

    assert compute_source_defined_idt(times, pressure, criterion_type) == pytest.approx(expected)


@pytest.mark.parametrize(
    ('criterion_type', 'expected'),
    [
        ('max', 3.0),
        ('1/2 max', 2.0),
    ],
)
def test_temperature_trace_peak_criteria_have_hand_known_answers(criterion_type, expected):
    """The peak is 4 at 3 s and its rising-side half-value crossing is exactly 2 s."""
    times = [0.0, 1.0, 2.0, 3.0, 4.0]
    temperature = [0.0, 1.0, 2.0, 4.0, 3.0]

    assert compute_source_defined_idt(times, temperature, criterion_type) == pytest.approx(expected)


def test_source_defined_criterion_refuses_degenerate_trace_and_unknown_type():
    """Trace validation and unknown criteria fail explicitly."""
    assert compute_source_defined_idt([0.0], [1.0], 'max') is None
    with pytest.raises(ValueError, match='Unknown source-defined'):
        compute_source_defined_idt([0.0, 1.0], [0.0, 1.0], 'onset')


@pytest.mark.parametrize('criterion_type', ['d/dt max', 'd/dt max extrapolated'])
def test_derivative_criterion_requires_post_ignition_resolution(criterion_type):
    """A continuing ramp is truncated, while a trace with a decayed slope is resolved."""
    times = [0, 1, 2, 3, 4, 5, 6, 7]

    assert not source_defined_idt_is_resolved(times, [300, 301, 302, 303, 304, 305, 306, 307],
                                              criterion_type, 'temperature')
    assert source_defined_idt_is_resolved(times, [300, 300, 310, 400, 500, 500, 500, 500],
                                          criterion_type, 'temperature')


def test_adaptive_derivative_criterion_accepts_declining_pressure_after_ignition():
    """A post-ignition pressure decline does not invalidate an interior slope peak."""
    times = [0, 1, 2, 3, 4, 5, 6]
    pressure = [1.0, 1.0, 2.0, 20.0, 19.5, 19.0, 18.5]

    assert source_defined_idt_is_resolved(
        times, pressure, 'd/dt max', 'pressure', adaptive_sampling=True,
    )


@pytest.mark.parametrize('criterion_type', ['max', '1/2 max'])
def test_peak_criterion_requires_peak_or_plateau_within_trace(criterion_type):
    """A boundary peak is truncated, while an interior plateau resolves the source criterion."""
    times = [0, 1, 2, 3, 4]

    assert not source_defined_idt_is_resolved(times, [1, 2, 3, 4, 5], criterion_type, 'pressure')
    assert source_defined_idt_is_resolved(times, [1, 2, 5, 5, 5], criterion_type, 'pressure')


@pytest.mark.parametrize('criterion_type', ['max', '1/2 max'])
def test_adaptive_peak_criterion_rejects_clustered_still_rising_tail(criterion_type):
    """A clustered adaptive tail does not resolve a trace that is still rising."""
    times = [0.0, 0.2, 0.5, 0.9, 0.9998, 0.9999, 1.0]
    temperature = [600.0 + 1500.0 * time for time in times]

    assert not source_defined_idt_is_resolved(
        times, temperature, criterion_type, 'temperature', adaptive_sampling=True,
    )


def test_adaptive_peak_criterion_accepts_clustered_tail_after_plateau():
    """Clustered adaptive samples resolve when the trace plateau spans elapsed time."""
    times = [0.0, 0.2, 0.5, 0.8, 0.9, 0.9998, 0.9999, 1.0]
    temperature = [600.0, 900.0, 1300.0, 2100.0, 2100.0, 2100.0, 2100.0, 2100.0]

    assert source_defined_idt_is_resolved(
        times, temperature, 'max', 'temperature', adaptive_sampling=True,
    )


def test_fixed_grid_peak_criterion_keeps_sample_tail_behavior():
    """Fixed-grid traces retain their existing sample-count tail behavior."""
    times = [0, 1, 2, 3, 4]

    assert not source_defined_idt_is_resolved(
        times, [1, 2, 3, 4, 5], 'max', 'pressure', adaptive_sampling=False,
    )
    assert source_defined_idt_is_resolved(
        times, [1, 2, 5, 5, 5], 'max', 'pressure', adaptive_sampling=False,
    )


def test_versioned_shock_tube_comparison_scores_point_and_refuses_excited_targets(tmp_path):
    """A direct shock-tube point is finite; absent OH* and CH* are typed refusals."""
    points = [_point()]
    for target in ('pressure', 'OH'):
        matched = copy.deepcopy(points[0])
        matched['ignition_definition']['target'] = target
        points.append(matched)
    for target in ('OH*', 'CH*'):
        refused = copy.deepcopy(points[0])
        refused['ignition_definition']['target'] = target
        points.append(refused)
    path = tmp_path / 'experimental.yaml'
    path.write_text(yaml.safe_dump({'version': 1, 'points': points}))

    result = _adapter().compare_with_experiment(str(path))

    assert result['n_points'] == 5
    assert result['n_matched'] == 3
    assert result['n_refused'] == 2
    assert result['refusals_by_reason'] == {'target species absent': 2}
    assert math.isfinite(result['points'][0]['idt_sim'])
    assert result['points'][0]['criterion'] == points[0]['ignition_definition']
    assert result['points'][0]['source'] == points[0]['source']
    assert result['points'][3]['refusal']['reason'] == 'target species absent'
    assert result['points'][4]['refusal']['reason'] == 'target species absent'


def test_versioned_comparison_refuses_unmappable_species(tmp_path):
    """A structurally valid SMILES absent from configured/model species is refused per point."""
    point = _point()
    point['composition'] = [{'smiles': '[Xe]', 'mole_fraction': 1.0}]
    path = tmp_path / 'unmappable.yaml'
    path.write_text(yaml.safe_dump({'version': 1, 'points': [point, copy.deepcopy(point)]}))

    result = _adapter().compare_with_experiment(str(path))

    assert result['n_matched'] == 0
    assert result['refusals_by_reason'] == {'unmappable species': 2}
    assert result['points'][0]['refusal']['reason'] == 'unmappable species'
    assert '[Xe]' in result['points'][0]['refusal']['detail']


def test_versioned_rcm_point_uses_post_compression_constant_pressure_reactor(tmp_path):
    """An RCM point is simulated directly at its post-compression state."""
    point = _point()
    point['apparatus'] = 'rapid compression machine'
    path = tmp_path / 'rcm.yaml'
    path.write_text(yaml.safe_dump({'version': 1, 'points': [point]}))

    result = _adapter().compare_with_experiment(str(path))

    assert result['n_matched'] == 1
    assert math.isfinite(result['points'][0]['idt_sim'])
    assert result['points'][0]['apparatus'] == 'rapid compression machine'


def test_versioned_point_reports_typed_unresolved_ignition(tmp_path):
    """A constant-pressure RCM half-maximum at t=0 is not counted as matched."""
    point = _point()
    point['apparatus'] = 'rapid compression machine'
    point['ignition_definition'] = {'target': 'pressure', 'type': '1/2 max'}
    path = tmp_path / 'degenerate-rcm.yaml'
    path.write_text(yaml.safe_dump({'version': 1, 'points': [point]}))

    result = _adapter().compare_with_experiment(str(path))

    assert result['n_matched'] == 0
    assert result['refusals_by_reason'] == {'ignition not resolved': 1}
    assert result['points'][0]['refusal']['reason'] == 'ignition not resolved'


def test_unversioned_comparison_preserves_official_main_result_shape_and_values(tmp_path):
    """A golden official/main legacy result remains byte-for-byte equal as Python data."""
    path = tmp_path / 'legacy.yaml'
    path.write_text(yaml.safe_dump({
        'citation': 'Legacy et al., 2024',
        'data': [{'T': 1001, 'P': 10, 'phi': 1.0, 'idt': 0.002}],
    }))
    adapter = _adapter()
    adapter.reactor_idt_dict = {0: {1.0: {10.0: {1000.0: 0.001}}}}

    result = adapter.compare_with_experiment(str(path))

    assert result == {
        'citation': 'Legacy et al., 2024',
        'n_points': 1,
        'n_matched': 1,
        'rmse_log': 0.3010299956639812,
        'points': [{
            'T': 1001.0,
            'P': 10.0,
            'phi': 1.0,
            'idt_exp': 0.002,
            'idt_sim': 0.001,
            'log10_error': -0.3010299956639812,
        }],
    }
