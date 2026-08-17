"""
t3 pdep diagram module

T3's own potential energy surface (PES) diagram for an explored P-dep network.

Arkane draws a network diagram too (``rmgpy.pdep.draw.NetworkDrawer``, run by
``arkane/explorer.py``), but T3 cannot call it -- T3 never imports ``rmgpy``/``arkane`` (see
``t3.pdep.parser``'s module docstring), and that drawer needs cairo, which the target RMG
environment does not reliably provide. This module renders the same surface from the network
FILE alone, through the safe never-execute parser, with matplotlib (a dependency T3 already
carries for ``t3.utils.flux``).

Energies are NORMALIZED: every level is shifted so the lowest-lying ISOMER sits at 0.0 kJ/mol.
The reference is deliberately the lowest isomer, never "the lowest point on the surface" -- a
bimolecular asymptote lying deeper than every well legitimately lands at a negative relative
energy, and silently re-anchoring the zero to it would hide exactly the feature (an exothermic
exit channel) worth seeing. Only a network with NO isomers falls back to the lowest
configuration overall, and the diagram data says so via ``reference_kind``.
"""

from dataclasses import dataclass

import matplotlib.pyplot as plt

from t3.pdep.parser import (PDepNetwork,
                            PDepNetworkE0,
                            parse_pdep_network_e0_file,
                            parse_pdep_network_file)

# The file name the explorer adapter writes the diagram under, inside the exploration run
# directory (next to Arkane's own input.py/output.py for that run).
T3_PES_DIAGRAM_FILENAME = 'pes_diagram.png'

REFERENCE_LOWEST_ISOMER = 'lowest_isomer'
REFERENCE_LOWEST_CONFIGURATION = 'lowest_configuration'

CONFIGURATION_KIND_ISOMER = 'isomer'
CONFIGURATION_KIND_REACTANT_CHANNEL = 'reactant_channel'
CONFIGURATION_KIND_PRODUCT_CHANNEL = 'product_channel'

# One categorical hue per configuration role plus one for transition states, assigned in fixed
# order by the job each color does (identity), all four validated together for colorblind-safe
# adjacent-pair separation on a white surface (OKLab CVD dE >= 8; the aqua's 2.7:1 surface
# contrast is relieved by every level carrying a direct text label).
_COLOR_ISOMER = '#2a78d6'             # blue
_COLOR_REACTANT_CHANNEL = '#eb6834'   # orange
_COLOR_PRODUCT_CHANNEL = '#1baf7a'    # aqua
_COLOR_TRANSITION_STATE = '#4a3aa7'   # violet
_COLOR_BARRIERLESS = '#7a7a7a'        # neutral gray, dashed: "connected, no barrier datum"
_COLOR_TEXT_PRIMARY = '#0b0b0b'
_COLOR_TEXT_SECONDARY = '#52514e'

_KIND_COLORS = {
    CONFIGURATION_KIND_ISOMER: _COLOR_ISOMER,
    CONFIGURATION_KIND_REACTANT_CHANNEL: _COLOR_REACTANT_CHANNEL,
    CONFIGURATION_KIND_PRODUCT_CHANNEL: _COLOR_PRODUCT_CHANNEL,
}

# Horizontal geometry, in column units (one column per configuration).
_CONFIG_HALF_WIDTH = 0.27
_TS_HALF_WIDTH = 0.14
_TS_SPREAD = 0.18  # extra x offset per additional TS sharing the same midpoint

# Opaque backing patch behind every text annotation: on dense networks the dotted TS connectors
# and level bars strike through labels that share their z-order space. A white halo under the
# text (zorder 4) and above the lines (zorder 2-3) guarantees legibility regardless of layout --
# cheaper and more robust than trying to prove no line ever crosses a label.
_LABEL_BBOX = dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='none', alpha=1.0)


@dataclass(frozen=True)
class PESConfiguration:
    """
    One stationary configuration of the surface: an isomer (well) or a reactant/product channel
    (an asymptote; possibly bimolecular, in which case ``e0`` is the sum over its fragments).

    Args:
        labels (tuple): The species labels, canonically (lexicographically) ordered.
        kind (str): One of the ``CONFIGURATION_KIND_*`` constants.
        e0 (float): The configuration E0 (kJ/mol) as declared in the network file.
        relative_e0 (float): The E0 after the normalization shift (kJ/mol).
    """
    labels: tuple
    kind: str
    e0: float
    relative_e0: float


@dataclass(frozen=True)
class PESTransitionState:
    """
    One path-reaction connection between two configurations.

    Args:
        label (Optional[str]): The transition state label, or ``None`` for a reaction declared
            with no ``transitionState`` keyword.
        e0 (Optional[float]): The TS E0 (kJ/mol), or ``None`` when the reaction is barrierless /
            the TS declares no E0 (drawn as a direct dashed connection, never as a fabricated
            barrier).
        relative_e0 (Optional[float]): The E0 after the normalization shift, or ``None``.
        reactants (tuple): The reactant-side configuration key (canonically ordered labels).
        products (tuple): The product-side configuration key (canonically ordered labels).
    """
    label: str | None
    e0: float | None
    relative_e0: float | None
    reactants: tuple
    products: tuple


@dataclass(frozen=True)
class PESDiagramData:
    """
    The computed, normalized energy levels of one explored PES.

    Args:
        network_id (str): The network identifier (file stem).
        reference_e0 (float): The absolute E0 (kJ/mol) subtracted from every level.
        reference_kind (str): ``REFERENCE_LOWEST_ISOMER``, or ``REFERENCE_LOWEST_CONFIGURATION``
            for the no-isomer fallback (see the module docstring).
        configurations (tuple): The ``PESConfiguration`` entries.
        transition_states (tuple): The ``PESTransitionState`` entries, one per path reaction.
    """
    network_id: str
    reference_e0: float
    reference_kind: str
    configurations: tuple
    transition_states: tuple


def compute_pes_diagram_data(network: PDepNetwork, e0: PDepNetworkE0) -> PESDiagramData:
    """
    Combine a parsed network topology with its declared E0 values into normalized diagram levels.

    Args:
        network (PDepNetwork): The parsed network topology.
        e0 (PDepNetworkE0): The declared E0 values (kJ/mol), from the same file.

    Raises:
        ValueError: If a species participating in any configuration declares no E0 (placing it at
            zero or dropping it would draw a diagram of a different network); if the network has
            no configurations at all; or if a path reaction side matches no known configuration.

    Returns:
        PESDiagramData: The normalized levels.
    """
    configurations = []
    seen_keys = dict()
    for kind, channels in (
            (CONFIGURATION_KIND_ISOMER, tuple((isomer,) for isomer in network.isomers)),
            (CONFIGURATION_KIND_REACTANT_CHANNEL, network.reactant_channels),
            (CONFIGURATION_KIND_PRODUCT_CHANNEL, network.product_channels)):
        for channel in channels:
            key = tuple(sorted(channel))
            if key in seen_keys:
                continue
            missing = sorted(label for label in key if label not in e0.species)
            if missing:
                raise ValueError(
                    f"Cannot place the {kind} configuration {key} of network "
                    f"'{network.network_id}' on the PES: species {missing} declare(s) no E0 in "
                    f"the network file. Refusing to draw the configuration at a fabricated "
                    f"energy or to silently omit it.")
            configurations.append((key, kind, sum(e0.species[label] for label in key)))
            seen_keys[key] = kind
    if not configurations:
        raise ValueError(f"Network '{network.network_id}' declares no isomers and no channels; "
                         f"there is no surface to draw.")

    isomer_energies = [energy for _, kind, energy in configurations
                       if kind == CONFIGURATION_KIND_ISOMER]
    if isomer_energies:
        reference_e0 = min(isomer_energies)
        reference_kind = REFERENCE_LOWEST_ISOMER
    else:
        reference_e0 = min(energy for _, _, energy in configurations)
        reference_kind = REFERENCE_LOWEST_CONFIGURATION

    transition_states = []
    for path_reaction in network.path_reactions:
        sides = []
        for side in (path_reaction.reactants, path_reaction.products):
            key = tuple(sorted(side))
            if key not in seen_keys:
                raise ValueError(
                    f"The path reaction '{path_reaction.label}' of network "
                    f"'{network.network_id}' has a side {key} matching no isomer or channel of "
                    f"the network; refusing to draw a connection to a configuration the network "
                    f"does not declare.")
            sides.append(key)
        ts_label = path_reaction.transition_state
        ts_e0 = e0.transition_states.get(ts_label) if ts_label is not None else None
        transition_states.append(PESTransitionState(
            label=ts_label,
            e0=ts_e0,
            relative_e0=ts_e0 - reference_e0 if ts_e0 is not None else None,
            reactants=sides[0],
            products=sides[1],
        ))

    return PESDiagramData(
        network_id=network.network_id,
        reference_e0=reference_e0,
        reference_kind=reference_kind,
        configurations=tuple(
            PESConfiguration(labels=key, kind=kind, e0=energy,
                             relative_e0=energy - reference_e0)
            for key, kind, energy in configurations),
        transition_states=tuple(transition_states),
    )


def stagger_label_offsets(points, min_dx: float, min_dy: float) -> tuple:
    """
    Resolve label collisions among near-degenerate levels, deterministically.

    Two levels closer than ``min_dx`` horizontally AND ``min_dy`` vertically would overlay their
    text labels; each such point is bumped to the smallest offset level (0, 1, 2, ...) distinct
    from every already-placed colliding neighbor, scanning in the caller's order. The caller maps
    an offset level to a vertical text displacement.

    Args:
        points: A sequence of ``(x, y)`` positions, in a deterministic order.
        min_dx (float): The horizontal distance below which two labels can collide.
        min_dy (float): The vertical distance below which two labels can collide.

    Returns:
        tuple: One integer offset level per point, in input order.
    """
    levels = []
    placed = []
    for x, y in points:
        taken = {level for (px, py, level) in placed
                 if abs(px - x) < min_dx and abs(py - y) < min_dy}
        level = 0
        while level in taken:
            level += 1
        levels.append(level)
        placed.append((x, y, level))
    return tuple(levels)


def _format_channel_label(labels: tuple) -> str:
    """
    Format a configuration's species labels for direct labeling (one species per line, joined
    with '+' -- multi-line so a bimolecular channel's label stays inside its own column).

    Args:
        labels (tuple): The configuration's species labels.

    Returns:
        str: The display label.
    """
    return ' +\n'.join(labels)


def _build_pes_figure(data: PESDiagramData):
    """
    Build (but do not save or close) the Matplotlib figure for a PES diagram.

    Split out of ``render_pes_diagram`` so tests can measure the actual rendered geometry
    (text bounding boxes, line positions) before the figure is closed.

    Layout: one column per configuration -- reactant channels, then isomers, then product
    channels, left to right -- each drawn as a horizontal level bar with a direct label
    (species names below, relative E0 above). Each path reaction with a TS energy gets a short
    TS level at the midpoint between the two configurations it connects, joined to both by thin
    dotted connectors; a reaction without a TS energy is a dashed straight line (barrierless --
    no fabricated barrier). Near-degenerate TS levels sharing a midpoint are spread apart, and
    colliding energy annotations are staggered (see ``stagger_label_offsets``).

    Args:
        data (PESDiagramData): The computed levels.

    Returns:
        matplotlib.figure.Figure: The built (unsaved, unclosed) figure.
    """
    ordered = sorted(
        data.configurations,
        key=lambda c: ({CONFIGURATION_KIND_REACTANT_CHANNEL: 0,
                        CONFIGURATION_KIND_ISOMER: 1,
                        CONFIGURATION_KIND_PRODUCT_CHANNEL: 2}[c.kind], c.labels))
    x_by_key = {config.labels: float(column) for column, config in enumerate(ordered)}

    energies = [c.relative_e0 for c in ordered] \
        + [ts.relative_e0 for ts in data.transition_states if ts.relative_e0 is not None]
    span = max(energies) - min(energies) or 1.0

    fig, ax = plt.subplots(figsize=(max(7.0, 1.9 * len(ordered) + 1.5), 6.5))

    # Configuration level bars with direct labels.
    for config in ordered:
        x = x_by_key[config.labels]
        color = _KIND_COLORS[config.kind]
        ax.hlines(y=config.relative_e0, xmin=x - _CONFIG_HALF_WIDTH, xmax=x + _CONFIG_HALF_WIDTH,
                  colors=color, linewidth=2.6, zorder=3)
        ax.annotate(_format_channel_label(config.labels),
                    xy=(x, config.relative_e0), xytext=(0, -10), textcoords='offset points',
                    ha='center', va='top', fontsize=7.5, color=_COLOR_TEXT_PRIMARY, zorder=4,
                    bbox=_LABEL_BBOX)
        ax.annotate(f'{config.relative_e0:.1f}',
                    xy=(x, config.relative_e0), xytext=(0, 4), textcoords='offset points',
                    ha='center', va='bottom', fontsize=7, color=_COLOR_TEXT_SECONDARY, zorder=4,
                    bbox=_LABEL_BBOX)

    # Transition state levels. Reactions sharing the same midpoint (e.g. two channels into the
    # same well pair) are spread horizontally, deterministically, before collision staggering.
    with_energy = [ts for ts in data.transition_states if ts.relative_e0 is not None]
    with_energy.sort(key=lambda ts: (ts.reactants, ts.products, ts.relative_e0, ts.label or ''))
    midpoint_counts = dict()
    ts_positions = []
    for ts in with_energy:
        x1, x2 = x_by_key[ts.reactants], x_by_key[ts.products]
        midpoint = (x1 + x2) / 2.0
        # A connection spanning non-adjacent columns can put its midpoint exactly on an
        # INTERMEDIATE column, where the TS level and its annotation land on top of that
        # column's own bar and labels. Nudge it half a column toward the product side: column
        # centers are integers, so a half-integer position is never on top of one.
        if midpoint not in (x1, x2) and midpoint in x_by_key.values():
            midpoint += 0.5 if x2 > x1 else -0.5
        occurrence = midpoint_counts.get(midpoint, 0)
        midpoint_counts[midpoint] = occurrence + 1
        # 0, +s, -s, +2s, ... around the shared midpoint.
        shift = ((occurrence + 1) // 2) * _TS_SPREAD * (1 if occurrence % 2 else -1)
        ts_positions.append((ts, midpoint + (shift if occurrence else 0.0), x1, x2))

    annotation_levels = stagger_label_offsets(
        points=tuple((x, ts.relative_e0) for ts, x, _, _ in ts_positions),
        min_dx=2.0 * _TS_SPREAD + 0.01, min_dy=0.05 * span)
    for (ts, x, x1, x2), annotation_level in zip(ts_positions, annotation_levels):
        ax.hlines(y=ts.relative_e0, xmin=x - _TS_HALF_WIDTH, xmax=x + _TS_HALF_WIDTH,
                  colors=_COLOR_TRANSITION_STATE, linewidth=2.2, zorder=3)
        for x_config, key in ((x1, ts.reactants), (x2, ts.products)):
            config_e0 = next(c.relative_e0 for c in ordered if c.labels == key)
            direction = 1 if x_config < x else -1
            ax.plot([x_config + direction * _CONFIG_HALF_WIDTH, x - direction * _TS_HALF_WIDTH],
                    [config_e0, ts.relative_e0],
                    color=_COLOR_TRANSITION_STATE, linewidth=0.9, linestyle=':', zorder=2)
        annotation = f'{ts.label} ({ts.relative_e0:.1f})' if ts.label else f'{ts.relative_e0:.1f}'
        ax.annotate(annotation,
                    xy=(x, ts.relative_e0), xytext=(0, 4 + 11 * annotation_level),
                    textcoords='offset points', ha='center', va='bottom', fontsize=7,
                    color=_COLOR_TEXT_SECONDARY, zorder=4, bbox=_LABEL_BBOX)

    # Barrierless connections: a direct dashed line, never a fabricated barrier.
    for ts in data.transition_states:
        if ts.relative_e0 is not None:
            continue
        x1, x2 = x_by_key[ts.reactants], x_by_key[ts.products]
        e1 = next(c.relative_e0 for c in ordered if c.labels == ts.reactants)
        e2 = next(c.relative_e0 for c in ordered if c.labels == ts.products)
        direction = 1 if x1 < x2 else -1
        ax.plot([x1 + direction * _CONFIG_HALF_WIDTH, x2 - direction * _CONFIG_HALF_WIDTH],
                [e1, e2], color=_COLOR_BARRIERLESS, linewidth=1.1, linestyle='--', zorder=2)
        if ts.label:
            ax.annotate(ts.label, xy=((x1 + x2) / 2.0, (e1 + e2) / 2.0),
                        xytext=(0, 4), textcoords='offset points', ha='center', va='bottom',
                        fontsize=7, color=_COLOR_TEXT_SECONDARY, zorder=4, bbox=_LABEL_BBOX)

    ax.axhline(y=0.0, color=_COLOR_TEXT_SECONDARY, linewidth=0.6, linestyle=(0, (6, 6)),
               alpha=0.45, zorder=1)
    reference_note = 'E0 relative to the lowest isomer' \
        if data.reference_kind == REFERENCE_LOWEST_ISOMER \
        else 'E0 relative to the lowest configuration (network declares no isomers)'
    ax.set_title(f'PES of {data.network_id} — {reference_note}',
                 fontsize=10, color=_COLOR_TEXT_PRIMARY)
    ax.set_ylabel('E0 (kJ/mol)', fontsize=9, color=_COLOR_TEXT_PRIMARY)
    ax.set_xticks([])
    ax.margins(x=0.06, y=0.18)
    for spine in ('top', 'right', 'bottom'):
        ax.spines[spine].set_visible(False)
    ax.yaxis.grid(True, linewidth=0.4, alpha=0.25)
    ax.set_axisbelow(True)
    fig.tight_layout()
    return fig


def render_pes_diagram(data: PESDiagramData, output_path: str) -> None:
    """
    Render computed PES levels to an image file (format follows the extension).

    Args:
        data (PESDiagramData): The computed levels.
        output_path (str): The image file path to write (e.g. ``.../pes_diagram.png``).
    """
    fig = _build_pes_figure(data)
    fig.savefig(output_path, dpi=300, facecolor='w', edgecolor='w')
    plt.close(fig)


def draw_pes_diagram(network_path: str, output_path: str) -> PESDiagramData:
    """
    Parse a pdep network file and render its normalized PES diagram.

    Args:
        network_path (str): The path to the network file (typically the exploration's resolved
            final network artifact, ``pdep/final/network<i>_(full|reduced).py``).
        output_path (str): The image file path to write (format follows the extension).

    Raises:
        ValueError: As ``parse_pdep_network_file``/``parse_pdep_network_e0_file``/
            ``compute_pes_diagram_data``.

    Returns:
        PESDiagramData: The computed levels the rendered diagram shows.
    """
    data = compute_pes_diagram_data(network=parse_pdep_network_file(path=network_path),
                                    e0=parse_pdep_network_e0_file(path=network_path))
    render_pes_diagram(data=data, output_path=output_path)
    return data
