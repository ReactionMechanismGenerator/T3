#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_diagram module

Tests T3's own PES energy diagram (``t3.pdep.diagram``): the pure energy-level computation
(normalization to the lowest-lying isomer) and the matplotlib rendering. The rendering tests
assert artifact existence and format only -- layout aesthetics are eyeballed, not pinned.
"""

import os

import matplotlib
matplotlib.use('Agg')  # headless rendering; must precede any pyplot-touching import below
import matplotlib.pyplot as plt

import pytest

from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.diagram import (
    REFERENCE_LOWEST_CONFIGURATION,
    REFERENCE_LOWEST_ISOMER,
    T3_PES_DIAGRAM_FILENAME,
    _build_pes_figure,
    compute_pes_diagram_data,
    draw_pes_diagram,
    stagger_label_offsets,
)
from t3.pdep.parser import parse_pdep_network_e0_file, parse_pdep_network_e0_text, \
    parse_pdep_network_file, parse_pdep_network_text

PDEP_NETWORK_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'iteration_1', 'RMG', 'pdep')
REAL_NETWORKS_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_real_networks')

NETWORK_799_PATH = os.path.join(REAL_NETWORKS_DIR, 'network799_1', 'network799_1.py')
NETWORK_21_PATH = os.path.join(REAL_NETWORKS_DIR, 'network21_1', 'network21_1.py')
NETWORK_4_1_PATH = os.path.join(PDEP_NETWORK_DIR, 'network4_1.py')
NETWORK_4_2_PATH = os.path.join(PDEP_NETWORK_DIR, 'network4_2.py')
NETWORK_1_1_PATH = os.path.join(PDEP_NETWORK_DIR, 'network1_1.py')

ALL_FIXTURE_NETWORK_PATHS = [NETWORK_799_PATH, NETWORK_21_PATH, NETWORK_4_1_PATH,
                             NETWORK_4_2_PATH, NETWORK_1_1_PATH]


def _read(path):
    with open(path) as f:
        return f.read()


def _compute_from_text(text):
    network = parse_pdep_network_text(text=text, network_id='synthetic')
    e0 = parse_pdep_network_e0_text(text=text)
    return compute_pes_diagram_data(network=network, e0=e0)


def _config_by_labels(data, labels):
    return next(c for c in data.configurations if set(c.labels) == set(labels))


class TestNormalization:

    def test_the_lowest_isomer_sits_exactly_at_zero(self):
        network = parse_pdep_network_text(text=_read(NETWORK_799_PATH), network_id='network799_1')
        e0 = parse_pdep_network_e0_text(text=_read(NETWORK_799_PATH))
        data = compute_pes_diagram_data(network=network, e0=e0)
        assert data.reference_kind == REFERENCE_LOWEST_ISOMER
        assert data.reference_e0 == pytest.approx(-170.357)
        isomer = _config_by_labels(data, {'O=C1OO1(21648)'})
        assert isomer.kind == 'isomer'
        assert isomer.relative_e0 == 0.0

    def test_relative_spacings_are_preserved_under_the_shift(self):
        text = _read(NETWORK_799_PATH)
        data = _compute_from_text(text)
        # Raw E0 sums: isomer -170.357; O(S)+CO2 = 432.331 - 402.51 = 29.821;
        # CO3t = -60.3782; [O]O[C]=O = 106.943. All shifted by +170.357:
        assert _config_by_labels(data, {'O-2(13598)', 'CO2(11)'}).relative_e0 == pytest.approx(200.178)
        assert _config_by_labels(data, {'[O]C([O])=O(8160)'}).relative_e0 == pytest.approx(109.9788)
        assert _config_by_labels(data, {'[O]O[C]=O(8100)'}).relative_e0 == pytest.approx(277.300)
        # And the shift is uniform: pairwise differences equal the raw pairwise differences.
        for config in data.configurations:
            assert config.relative_e0 == pytest.approx(config.e0 - data.reference_e0)
        for ts in data.transition_states:
            assert ts.relative_e0 == pytest.approx(ts.e0 - data.reference_e0)

    def test_the_lowest_of_several_isomers_is_the_reference(self):
        text = _read(NETWORK_4_2_PATH)
        data = _compute_from_text(text)
        assert data.reference_kind == REFERENCE_LOWEST_ISOMER
        # Isomers: C4rad(5) at 63.0573, butyl_2(67) at 52.4618 kJ/mol -- the LOWER one is zero.
        assert data.reference_e0 == pytest.approx(52.4618)
        assert _config_by_labels(data, {'butyl_2(67)'}).relative_e0 == 0.0
        assert _config_by_labels(data, {'C4rad(5)'}).relative_e0 == pytest.approx(10.5955)

    def test_a_bimolecular_asymptote_below_the_lowest_isomer_stays_below_zero(self):
        """The spec's edge case, decided deliberately: the reference is the lowest ISOMER, never
        silently redefined to the lowest point overall -- so an exit channel lying deeper than
        every well legitimately lands at a negative relative energy."""
        text = """species(label='I', E0=(0.0, 'kJ/mol'))
species(label='A', E0=(-10.0, 'kJ/mol'))
species(label='B', E0=(-20.0, 'kJ/mol'))
transitionState(label='TS1', E0=(15.0, 'kJ/mol'))
reaction(label='r1', reactants=['I'], products=['A', 'B'], transitionState='TS1')
network(label='N', isomers=['I'], reactants=[])
"""
        data = _compute_from_text(text)
        assert data.reference_kind == REFERENCE_LOWEST_ISOMER
        assert data.reference_e0 == pytest.approx(0.0)
        assert _config_by_labels(data, {'I'}).relative_e0 == 0.0
        assert _config_by_labels(data, {'A', 'B'}).relative_e0 == pytest.approx(-30.0)

    def test_a_network_with_no_isomers_falls_back_to_the_lowest_configuration(self):
        """With no isomer the spec's reference is undefined; the deliberate fallback is the lowest
        configuration overall, and the data says so out loud via reference_kind."""
        text = """species(label='A', E0=(10.0, 'kJ/mol'))
species(label='B', E0=(20.0, 'kJ/mol'))
species(label='C', E0=(5.0, 'kJ/mol'))
transitionState(label='TS1', E0=(50.0, 'kJ/mol'))
reaction(label='r1', reactants=['A', 'B'], products=['C'], transitionState='TS1')
network(label='N', isomers=[], reactants=[('A', 'B')])
"""
        data = _compute_from_text(text)
        assert data.reference_kind == REFERENCE_LOWEST_CONFIGURATION
        assert data.reference_e0 == pytest.approx(5.0)
        assert _config_by_labels(data, {'C'}).relative_e0 == 0.0
        assert _config_by_labels(data, {'A', 'B'}).relative_e0 == pytest.approx(25.0)


class TestConfigurationsAndTransitionStates:

    def test_bath_gas_species_participate_in_no_configuration(self):
        """network799_1 declares Ar and Ne (bath gas) as species with E0 values; they are part of
        no isomer or channel and must not appear as diagram configurations."""
        data = _compute_from_text(_read(NETWORK_799_PATH))
        all_labels = {label for config in data.configurations for label in config.labels}
        assert 'Ar' not in all_labels
        assert 'Ne' not in all_labels
        assert len(data.configurations) == 4  # 1 isomer + 3 derived product channels

    def test_configuration_kinds(self):
        data = _compute_from_text(_read(NETWORK_4_2_PATH))
        kinds = {tuple(sorted(c.labels)): c.kind for c in data.configurations}
        assert kinds[('C4rad(5)',)] == 'isomer'
        assert kinds[('butyl_2(67)',)] == 'isomer'
        assert kinds[('C4ene(26)', 'H(34)')] == 'reactant_channel'

    def test_each_transition_state_connects_two_known_configurations(self):
        data = _compute_from_text(_read(NETWORK_799_PATH))
        config_keys = {tuple(sorted(c.labels)) for c in data.configurations}
        assert len(data.transition_states) == 3
        for ts in data.transition_states:
            assert ts.reactants in config_keys
            assert ts.products in config_keys
        ts1 = next(ts for ts in data.transition_states if ts.label == 'TS1')
        assert ts1.e0 == pytest.approx(31.8584)
        assert ts1.relative_e0 == pytest.approx(202.2154)

    def test_a_reaction_without_a_transition_state_is_barrierless(self):
        text = """species(label='I', E0=(0.0, 'kJ/mol'))
species(label='P', E0=(-5.0, 'kJ/mol'))
reaction(label='r1', reactants=['I'], products=['P'])
network(label='N', isomers=['I'], reactants=[])
"""
        data = _compute_from_text(text)
        (ts,) = data.transition_states
        assert ts.label is None
        assert ts.e0 is None
        assert ts.relative_e0 is None

    def test_a_labeled_transition_state_without_an_e0_keeps_its_label_and_no_energy(self):
        text = """species(label='I', E0=(0.0, 'kJ/mol'))
species(label='P', E0=(-5.0, 'kJ/mol'))
transitionState(label='TS1')
reaction(label='r1', reactants=['I'], products=['P'], transitionState='TS1')
network(label='N', isomers=['I'], reactants=[])
"""
        data = _compute_from_text(text)
        (ts,) = data.transition_states
        assert ts.label == 'TS1'
        assert ts.e0 is None
        assert ts.relative_e0 is None

    def test_a_participating_species_with_no_e0_is_refused_by_name(self):
        """A configuration whose species has no declared E0 cannot be placed on the surface;
        silently dropping it (or placing it at zero) would draw a diagram of a different network."""
        text = """species(label='I', E0=(0.0, 'kJ/mol'))
species(label='P')
reaction(label='r1', reactants=['I'], products=['P'])
network(label='N', isomers=['I'], reactants=[])
"""
        network = parse_pdep_network_text(text=text, network_id='synthetic')
        e0 = parse_pdep_network_e0_text(text=text)
        with pytest.raises(ValueError, match="'P'"):
            compute_pes_diagram_data(network=network, e0=e0)


class TestStaggerLabelOffsets:

    def test_two_near_degenerate_points_get_distinct_levels(self):
        levels = stagger_label_offsets(points=((0.0, 100.0), (0.5, 101.0)), min_dx=0.7, min_dy=5.0)
        assert levels[0] != levels[1]

    def test_well_separated_points_all_sit_at_level_zero(self):
        levels = stagger_label_offsets(points=((0.0, 0.0), (2.0, 0.0), (0.1, 300.0)),
                                       min_dx=0.7, min_dy=5.0)
        assert levels == (0, 0, 0)

    def test_a_cluster_of_three_gets_three_distinct_levels(self):
        levels = stagger_label_offsets(points=((0.0, 10.0), (0.2, 11.0), (0.4, 12.0)),
                                       min_dx=0.7, min_dy=5.0)
        assert len(set(levels)) == 3

    def test_the_result_is_deterministic(self):
        points = ((0.3, 5.0), (0.0, 5.5), (0.6, 4.5))
        assert stagger_label_offsets(points=points, min_dx=0.7, min_dy=5.0) \
            == stagger_label_offsets(points=points, min_dx=0.7, min_dy=5.0)


class TestDrawPESDiagram:

    @pytest.mark.parametrize('network_path', [NETWORK_799_PATH, NETWORK_21_PATH, NETWORK_4_2_PATH])
    def test_a_diagram_is_produced_for_a_representative_network(self, tmp_path, network_path):
        output_path = os.path.join(str(tmp_path), T3_PES_DIAGRAM_FILENAME)
        data = draw_pes_diagram(network_path=network_path, output_path=output_path)
        assert os.path.isfile(output_path)
        assert os.path.getsize(output_path) > 0
        # PNG magic bytes: the format follows the file extension.
        with open(output_path, 'rb') as f:
            assert f.read(8) == b'\x89PNG\r\n\x1a\n'
        assert data.configurations

    def test_the_format_follows_the_extension(self, tmp_path):
        output_path = os.path.join(str(tmp_path), 'pes_diagram.svg')
        draw_pes_diagram(network_path=NETWORK_799_PATH, output_path=output_path)
        with open(output_path, 'rb') as f:
            head = f.read(200)
        assert b'<svg' in head or b'<?xml' in head

    def test_drawing_a_no_isomer_network_annotates_the_fallback_reference(self, tmp_path):
        text = """species(label='A', E0=(10.0, 'kJ/mol'))
species(label='B', E0=(20.0, 'kJ/mol'))
species(label='C', E0=(5.0, 'kJ/mol'))
transitionState(label='TS1', E0=(50.0, 'kJ/mol'))
reaction(label='r1', reactants=['A', 'B'], products=['C'], transitionState='TS1')
network(label='N', isomers=[], reactants=[('A', 'B')])
"""
        network_path = os.path.join(str(tmp_path), 'network_no_isomer.py')
        with open(network_path, 'w') as f:
            f.write(text)
        output_path = os.path.join(str(tmp_path), T3_PES_DIAGRAM_FILENAME)
        data = draw_pes_diagram(network_path=network_path, output_path=output_path)
        assert data.reference_kind == REFERENCE_LOWEST_CONFIGURATION
        assert os.path.getsize(output_path) > 0


def _bbox_overlap_area(a, b):
    """Return the overlap area (px^2) of two Matplotlib Bbox objects, 0.0 if disjoint."""
    dx = min(a.x1, b.x1) - max(a.x0, b.x0)
    dy = min(a.y1, b.y1) - max(a.y0, b.y0)
    return dx * dy if dx > 0 and dy > 0 else 0.0


def _segment_crosses_box(p0, p1, box):
    """Return True if the segment p0->p1 (display coords) enters the interior of `box`.

    Liang-Barsky clipping: solve for the parametric range [t0, t1] along the segment that
    survives clipping against each of the box's four half-planes: a nonempty, non-degenerate
    range means the segment actually enters the box's interior.
    """
    x0, y0 = float(p0[0]), float(p0[1])
    dx, dy = float(p1[0]) - x0, float(p1[1]) - y0
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, x0 - box.x0), (dx, box.x1 - x0),
                 (-dy, y0 - box.y0), (dy, box.y1 - y0)):
        if p == 0.0:
            if q < 0.0:
                return False
            continue
        r = q / p
        if p < 0.0:
            t0 = max(t0, r)
        else:
            t1 = min(t1, r)
        if t0 > t1:
            return False
    return True


def _measure_rendered_geometry(data):
    """Render `data` via ``_build_pes_figure`` and return (texts, boxes, segments).

    ``texts``/``boxes`` are parallel lists of non-blank text artists and their display-space
    bounding boxes (``get_window_extent``); ``segments`` are all drawn line segments (level
    bars from ``ax.hlines`` land as ``LineCollection``s in ``ax.collections``; TS connectors
    and barrierless dashes are plain ``ax.lines``), also in display space.
    """
    fig = _build_pes_figure(data=data)
    try:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        texts = [t for ax in fig.axes for t in ax.texts if t.get_text().strip()]
        boxes = [t.get_window_extent(renderer=renderer) for t in texts]
        segments = []
        for ax in fig.axes:
            for line in ax.lines:
                pts = ax.transData.transform(line.get_xydata())
                segments.extend(zip(pts[:-1], pts[1:]))
            for coll in ax.collections:  # ax.hlines(...) -> LineCollection
                for seg in coll.get_segments():
                    pts = ax.transData.transform(seg)
                    segments.extend(zip(pts[:-1], pts[1:]))
        return texts, boxes, segments
    finally:
        plt.close(fig)


def _has_opaque_backing(text):
    """Return True if `text` carries an opaque (alpha >= 0.8) backing patch."""
    patch = text.get_bbox_patch()
    return patch is not None and (patch.get_alpha() or 1.0) >= 0.8


class TestRenderedLabelLegibility:
    """
    Verifies the ACTUAL rendered geometry, not just the data model: dense networks used to have
    text annotations struck through by the dotted TS connectors and level bars they overlap in
    z-order. Every annotation now carries an opaque backing patch (see ``_LABEL_BBOX`` in
    ``t3.pdep.diagram``), so a line crossing a label's box is only a defect when that label has
    no such patch. Separately, ``stagger_label_offsets`` is responsible for zero text-vs-text
    overlaps; both properties are checked here on every fixture network.
    """

    @pytest.mark.parametrize('network_path', ALL_FIXTURE_NETWORK_PATHS)
    def test_no_two_labels_overlap(self, network_path):
        data = compute_pes_diagram_data(network=parse_pdep_network_file(path=network_path),
                                        e0=parse_pdep_network_e0_file(path=network_path))
        texts, boxes, _ = _measure_rendered_geometry(data=data)
        overlapping = [(texts[i].get_text(), texts[j].get_text())
                       for i in range(len(boxes)) for j in range(i + 1, len(boxes))
                       if _bbox_overlap_area(boxes[i], boxes[j]) > 0.0]
        assert overlapping == []

    @pytest.mark.parametrize('network_path', ALL_FIXTURE_NETWORK_PATHS)
    def test_no_label_is_struck_through_by_a_line(self, network_path):
        data = compute_pes_diagram_data(network=parse_pdep_network_file(path=network_path),
                                        e0=parse_pdep_network_e0_file(path=network_path))
        texts, boxes, segments = _measure_rendered_geometry(data=data)
        struck = []
        for text, box in zip(texts, boxes):
            if _has_opaque_backing(text):
                continue  # an opaque backing patch makes any crossing line harmless
            if any(_segment_crosses_box(p0, p1, box) for p0, p1 in segments):
                struck.append(text.get_text())
        assert struck == []
