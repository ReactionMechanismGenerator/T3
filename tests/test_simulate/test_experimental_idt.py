"""Tests for versioned, per-point experimental IDT comparisons."""

import copy
import math
import os

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


@pytest.mark.parametrize('criterion_type', ['max', '1/2 max'])
def test_peak_criterion_requires_peak_or_plateau_within_trace(criterion_type):
    """A boundary peak is truncated, while an interior plateau resolves the source criterion."""
    times = [0, 1, 2, 3, 4]

    assert not source_defined_idt_is_resolved(times, [1, 2, 3, 4, 5], criterion_type, 'pressure')
    assert source_defined_idt_is_resolved(times, [1, 2, 5, 5, 5], criterion_type, 'pressure')


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
