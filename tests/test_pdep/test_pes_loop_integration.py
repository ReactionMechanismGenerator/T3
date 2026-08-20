"""
t3 tests test_pdep test_pes_loop_integration module

N6 (round 3 of the Task 6 review): nothing anywhere paired the REAL ``t3.pdep.pes_loop.run_pes_loop``
with the REAL ``t3.pdep.pes_qm.arc_qm_runner``. Every ``test_pes_loop.py`` test injects a
hand-written ``qm_runner`` double whose two-tuple return can never disagree with what the loop
offered (that gap is exactly what round 3's N3 fix addressed on the loop side). Every
``test_pes_qm.py`` test drives ``arc_qm_runner`` directly, with no ``run_pes_loop`` around it, so it
can never observe whether the hybrid network file ``arc_qm_runner`` writes for round N is actually
the file ``run_pes_loop`` hands its round-(N+1) explorer -- each half's tests only ever check that
convention against itself.

This module is that missing pairing. It fakes only the three boundaries ``arc_qm_runner`` itself
cannot cross in a test process -- ``ARC`` (no cluster here), ``capture_ts_artifacts`` (no real QM
job to capture from), and ``write_hybrid_network_input_file`` (no real QM artifacts to fold into a
hybrid network) -- with argument-recording doubles. ``explore_pdep_network`` is also faked, the same
way every other ``pes_loop`` test fakes it (Arkane itself is out of scope for this loop's own test
suite, and always has been); everything downstream of that -- the real ``parse_pdep_network_file``,
the real ``split_qm_candidates``, the real ``draw_pes_diagram``, the real ``build_arc_input``, and
every join/capture computation inside ``arc_qm_runner`` -- runs unfaked, against the real
``network799_1`` fixture also used by ``test_pes_qm.py``.

Wiring ``run_pes_loop``'s ``qm_runner`` argument into ``t3.main`` stays deliberately out of scope --
this loop is standalone until corroborated (see ``pes_loop.py``'s own module docstring) -- this
module only corroborates that the loop and the real runner agree with each other.
"""

import os
import re
import shutil
import sys

import yaml

import PES
from t3.pdep.capture import (CaptureResult, capture_ts_artifacts, captured_qm_artifact_path,
                             verify_capture)
from t3.pdep.discovery import ARTIFACT_STATUS_USABLE, TSArtifactRecord
from t3.pdep.explorer.result import EXPLORATION_STATUS_SUCCEEDED, PDepExplorationResult
from t3.pdep.hashing import hash_file
from t3.pdep.hybrid import HybridNetworkResult
from t3.pdep.join import (JOIN_STATUS_QUEUED, TSJoinRecord, arc_ts_label,
                          expected_ts_artifact_path)
from t3.pdep.parser import parse_pdep_network_file
import t3.pdep.pes_qm as pes_qm
from t3.pdep.pes_loop import run_pes_loop
from t3.pdep.pes_qm import _explored_network_path, arc_qm_runner
from t3.pdep.pes_rounds import hybrid_network_path, round_paths, split_qm_candidates
from t3.schema import PESLoopConfig

_FIXTURE_NETWORK_PATH = os.path.join(
    os.path.dirname(__file__), '..', 'data', 'pdep_real_networks', 'network799_1', 'network799_1.py')

# The same frozen energy-settings shape t3.pdep.pes_qm.QMEnergySettings.from_frozen expects,
# mirroring test_pes_qm.py's own _FROZEN_ENERGY_SETTINGS.
_ENERGY_SETTINGS = {
    'model_chemistry': 'wb97xd/def2tzvp',
    'frequency_scale_factor': None,
    'use_hindered_rotors': True,
    'use_bond_corrections': False,
    'bond_correction_type': None,
    'atom_energies': None,
    'use_atom_corrections': True,
    'bond_additivity_corrections': None,
}


class _FakeARC(object):
    """Stands in for arc.main.ARC: no cluster, no execute -- just records what it was built with."""

    constructions = []
    execute_calls = 0

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        _FakeARC.constructions.append(kwargs)

    def execute(self):
        _FakeARC.execute_calls += 1

    def as_dict(self):
        return dict(self.kwargs)


def _config(network_path, max_rounds=3):
    return PESLoopConfig(pes={'network': network_path, 'source': ['A'],
                              'bath_gas': {'He': 1.0}},
                         termination={'max_rounds': max_rounds, 'stop_when_no_new_ts': False})


def test_real_run_pes_loop_wires_the_real_arc_qm_runner_across_rounds(tmp_path, monkeypatch):
    """Drive the real run_pes_loop with the real arc_qm_runner injected as qm_runner, faking only
    ARC, capture_ts_artifacts, and write_hybrid_network_input_file. Assert the loop completes at
    least two rounds, and -- the actual N6 contract -- that the exact hybrid file arc_qm_runner
    writes for round 0 is the exact file run_pes_loop hands its round-1 explorer, observed end to
    end through the real loop rather than inferred from two separately-passing unit tests."""
    # The fixture file's own stem ('network799_1') follows T3's own hybrid-file convention
    # (network<digits>_<digits>), not the real Arkane explorer's
    # ('^network(\\d+)_(full|reduced)\\.py$', t3.pdep.explorer.arkane._FINAL_NETWORK_FILENAME_RE).
    # This fake explorer stands in for that real explorer, so it must hand downstream a network_id
    # shaped the way Arkane's really is -- parse from a renamed copy rather than the fixture path
    # directly, so real_network.network_id and every path derived from it looks like what Arkane
    # would actually produce.
    fixture_network_path = os.path.join(str(tmp_path), 'network0_full.py')
    shutil.copyfile(_FIXTURE_NETWORK_PATH, fixture_network_path)
    real_network = parse_pdep_network_file(path=fixture_network_path)
    real_split = split_qm_candidates(real_network, frozenset())
    assert real_split.candidates, 'fixture must offer at least one real QM candidate'
    target = real_split.candidates[0]

    project_directory = str(tmp_path)
    # What run_pes_loop actually handed this fake as `network_path` each round -- the N6 contract
    # (round 1 gets round 0's hybrid file) is about THIS list, not about wherever the fake itself
    # chooses to write.
    received_network_paths = []

    def _fake_explore(*, network_path, config, logger=None):
        received_network_paths.append(network_path)
        # Not one of the three edges N6 names -- Arkane itself stays faked here exactly as it is
        # in every other pes_loop test -- but the file this hands downstream is a real, on-disk
        # copy of network799_1, written to the SAME canonical location the real Arkane explorer
        # would use (paths.explorer_output/pdep/final/<network_id>.py, round_paths(...)'s own
        # convention), because arc_qm_runner's real, unfaked _explored_network_path(paths,
        # network_id) recomputes that exact path independently rather than trusting whatever this
        # fake returns in network_paths -- so writing anywhere else would desync the two halves.
        round_index = len(received_network_paths) - 1
        paths = round_paths(project_directory, round_index)
        dest_path = _explored_network_path(paths, real_network.network_id)
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        shutil.copyfile(fixture_network_path, dest_path)
        return PDepExplorationResult(network_id=real_network.network_id,
                                     status=EXPLORATION_STATUS_SUCCEEDED,
                                     network_paths=(dest_path,))

    monkeypatch.setattr('t3.pdep.pes_loop.explore_pdep_network', _fake_explore)

    _FakeARC.constructions = []
    _FakeARC.execute_calls = 0
    monkeypatch.setattr(pes_qm, 'ARC', _FakeARC)

    write_hybrid_calls = []

    def _fake_capture_ts_artifacts(*, join_records, arc_project_directory, capture_dir, networks,
                                   sensitivity_by_ts, sensitivity_aggregation):
        # Defect 1: the real capture_ts_artifacts refuses (via its verify_capture self-check) any
        # captured artifact whose join record carries no finite sensitivity evidence -- so even
        # this double must insist the runner actually passed it, keyed and finite, or the pairing
        # this module corroborates would still be one that cannot run against the real capture.
        assert sensitivity_by_ts, 'arc_qm_runner must pass the sensitivity evidence to capture'
        for key, (coefficient, delta_ln_k) in sensitivity_by_ts.items():
            assert coefficient is not None and delta_ln_k is not None, key
        os.makedirs(capture_dir, exist_ok=True)
        # Defect 3: the loop carries a converged TS across the round boundary as an adopted
        # artifact at the capture's own vendored path (captured_qm_artifact_path), and the REAL
        # _vendor_adopted_artifacts in round 1 fails closed on a missing file -- so this double
        # must actually write the artifact it claims was captured, not merely name it in a result
        # object (the same rule the hybrid double below already follows).
        vendored_path = captured_qm_artifact_path(
            capture_dir, arc_ts_label(real_network.network_id, target.ts_label))
        os.makedirs(os.path.dirname(vendored_path), exist_ok=True)
        with open(os.path.join(os.path.dirname(vendored_path), 'output.out'), 'w') as f:
            f.write('# stub ARC log output\n')
        with open(vendored_path, 'w') as f:
            f.write("geometry = Log('output.out')\n")
        return CaptureResult(
            capture_dir=capture_dir,
            manifest_path=os.path.join(capture_dir, 'manifest.yml'),
            records=(
                TSArtifactRecord(network_id=real_network.network_id,
                                 network_ts_label=target.ts_label,
                                 arc_ts_label=f'{real_network.network_id}_{target.ts_label}',
                                 status=ARTIFACT_STATUS_USABLE,
                                 artifact_path=os.path.join(capture_dir, f'{target.ts_label}.py'),
                                 converged=True,
                                 path_reaction_labels=(target.path_reaction.label,)),
            ),
            energy_settings=_ENERGY_SETTINGS,
            networks={
                real_network.network_id: {
                    'source_path': _FIXTURE_NETWORK_PATH,
                    'captured_path': f'{target.ts_label}.py',
                    'method': 'MSC',
                },
            },
        )

    def _fake_write_hybrid_network_input_file(*, source_path, dest_path, method,
                                              qm_transition_states, energy_settings,
                                              qm_artifacts_root):
        write_hybrid_calls.append({'source_path': source_path, 'dest_path': dest_path})
        # run_pes_loop's own round>0 guard checks for this file's existence on disk for real, so
        # the double must actually write it, not merely return a result object claiming to.
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        open(dest_path, 'w').close()
        return HybridNetworkResult(dest_path=dest_path, qm_ts_labels=tuple(qm_transition_states),
                                   ilt_ts_labels=(), vendored_files=(), warnings=())

    monkeypatch.setattr(pes_qm, 'capture_ts_artifacts', _fake_capture_ts_artifacts)
    monkeypatch.setattr(pes_qm, 'write_hybrid_network_input_file',
                        _fake_write_hybrid_network_input_file)

    config = _config(fixture_network_path, max_rounds=3)
    result = run_pes_loop(config, project_directory=str(tmp_path), qm_runner=arc_qm_runner)

    assert len(result.rounds) >= 2, (
        f'expected at least two rounds, got {len(result.rounds)}: '
        f'{[r.status for r in result.rounds]} ({result.status}: {result.reason})')
    assert len(write_hybrid_calls) >= 1
    round0_hybrid_path = write_hybrid_calls[0]['dest_path']
    assert os.path.isfile(round0_hybrid_path)
    assert len(received_network_paths) >= 2
    # The actual N6 contract: round 1's explorer must be handed EXACTLY the file arc_qm_runner
    # wrote for round 0.
    assert received_network_paths[1] == round0_hybrid_path


# ---------------------------------------------------------------------------------------------
# Real-capture pairing (this fix round). Everything below drives the REAL run_pes_loop with the
# REAL arc_qm_runner, the REAL capture_ts_artifacts (no capture double -- the gap every prior
# test on this branch shared, which is how defects 1-3 stayed green), the REAL Arkane
# master-equation SA (t3.pdep.pes_sa -- Arkane actually runs, once per round with candidates),
# the REAL adopt_prior_qm against a REAL prior capture, the REAL vendoring, and the REAL hybrid
# writer. Faked: ONLY the explorer (as everywhere in this suite) and ARC (execute() writes
# converged statmech artifacts into the round's own ARC project instead of submitting cluster
# jobs).
#
# The fake explorer RENUMBERS between rounds -- from round 1 on it swaps the fixture's
# 'TS1'<->'TS2' and 'reaction1'<->'reaction2' labels while the channels stay put -- because every
# label Arkane writes is purely positional (rmgpy/rmg/pdep.py:854) and a test whose labels are
# trivially stable is blind to the wrong-saddle-point misattribution the structural carry exists
# to prevent. Assertions are therefore by CHANNEL (reactant/product species), never by label, and
# each fake ARC artifact carries a channel marker in its (hash-verified, verbatim-vendored) log
# file so "the right barrier sits on the right channel" is checked on the bytes, not inferred.
# ---------------------------------------------------------------------------------------------

_MODEL_CHEMISTRY = "LevelOfTheory(method='wb97xd2023',basis='def2tzvp',software='gaussian')"

# The fixture's three channels, by species labels (direction-insensitively sorted, as
# _hybrid_channels below renders them).
_CH1 = ('CO2(11) + O-2(13598)', 'O=C1OO1(21648)')
_CH2 = ('O=C1OO1(21648)', '[O]C([O])=O(8160)')
_CH3 = ('O=C1OO1(21648)', '[O]O[C]=O(8100)')


def _channel_marker(reactants, products) -> str:
    """A stable, label-free identity line for one channel, embedded in fake QM log files."""
    return 'barrier-marker: ' + ' <=> '.join(
        sorted((' + '.join(sorted(reactants)), ' + '.join(sorted(products)))))


def _swap_quoted(text: str, a: str, b: str) -> str:
    """Swap every quoted occurrence of two labels in a network file's text."""
    return text.replace(f"'{a}'", "'@@SWAP@@'").replace(f"'{b}'", f"'{a}'") \
               .replace("'@@SWAP@@'", f"'{b}'")


def _renumbered_fixture_text() -> str:
    """The real fixture with 'TS1'<->'TS2' and 'reaction1'<->'reaction2' swapped: the same three
    channels, renumbered exactly the way a pruning/discovery re-exploration renumbers them."""
    with open(_FIXTURE_NETWORK_PATH) as f:
        text = f.read()
    text = _swap_quoted(text, 'TS1', 'TS2')
    return _swap_quoted(text, 'reaction1', 'reaction2')


def _renamed_prior_network_text() -> str:
    """The real fixture as a DIFFERENT project would have written it: every species label carries
    that project's own indices, and the TS labels sit at different positions ('TS1'<->'TS3').
    The adjacency lists -- the structures -- are untouched, so structural matching must still
    recognize the channels."""
    with open(_FIXTURE_NETWORK_PATH) as f:
        text = f.read()
    text = _swap_quoted(text, 'TS1', 'TS3')
    for original, renamed in (('O=C1OO1(21648)', 'O=C1OO1(9)'),
                              ('O-2(13598)', 'O-2(5)'),
                              ('CO2(11)', 'CO2(7)'),
                              ('[O]C([O])=O(8160)', '[O]C([O])=O(2)'),
                              ('[O]O[C]=O(8100)', '[O]O[C]=O(3)')):
        text = text.replace(f"'{original}'", f"'{renamed}'")
    return text


def _write_artifact(path, log_text='stub quantum chemistry log\n'):
    # One log file PER artifact, named after it: ARC writes every transition state's statmech
    # input into one shared TSs/ directory, so a shared 'output.out' would silently overwrite
    # each artifact's channel marker with the last one written.
    log_name = os.path.splitext(os.path.basename(path))[0] + '.out'
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(os.path.join(os.path.dirname(path), log_name), 'w') as f:
        f.write(log_text)
    with open(path, 'w') as f:
        f.write(f"linear = False\n\nspinMultiplicity = 2\n\nenergy = Log('{log_name}')\n\n"
                f"geometry = Log('{log_name}')\n\nfrequencies = Log('{log_name}')\n")


def _write_status(arc_dir, label, converged):
    output_dir = os.path.join(arc_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'status.yml'), 'a') as f:
        f.write(f"{label}:\n  convergence: {str(converged).lower()}\n  job_types: {{}}\n"
                f"  paths: {{}}\n  info: ''\n  errors: ''\n")


def _write_energy_settings(arc_dir):
    statmech_dir = os.path.join(arc_dir, 'calcs', 'statmech', 'kinetics')
    os.makedirs(statmech_dir, exist_ok=True)
    with open(os.path.join(statmech_dir, 'input.py'), 'w') as f:
        f.write(f"modelChemistry = {_MODEL_CHEMISTRY}\n\nuseHinderedRotors = True\n\n"
                "useAtomCorrections = True\n\n"
                "atomEnergies = {'C': -37.844411, 'H': -0.499818, 'O': -75.062219}\n\n"
                "useBondCorrections = False\n")
    output_dir = os.path.join(arc_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    output_yml = os.path.join(output_dir, 'output.yml')
    if not os.path.isfile(output_yml):
        with open(output_yml, 'w') as f:
            f.write('{}\n')


def _build_prior_capture(prior_project_dir, network_text, ts_labels):
    """A REAL prior capture (built by capture_ts_artifacts itself) whose vendored network is a
    REAL, parseable network file -- ``network_text`` -- and whose artifacts carry channel markers
    in their log files. ``ts_labels`` are labels of THAT network's own numbering."""
    arc_dir = os.path.join(prior_project_dir, 'arc_project')
    network_id = 'prior_run_network'
    networks_dir = os.path.join(prior_project_dir, 'networks')
    os.makedirs(networks_dir, exist_ok=True)
    source_path = os.path.join(networks_dir, f'{network_id}.py')
    with open(source_path, 'w') as f:
        f.write(network_text)
    prior_network = parse_pdep_network_file(path=source_path)
    reactions_by_ts = prior_network.path_reactions_by_ts()
    records = []
    for ts_label in ts_labels:
        label = arc_ts_label(network_id, ts_label)
        record = TSJoinRecord(network_id=network_id, network_ts_label=ts_label,
                              status=JOIN_STATUS_QUEUED, arc_ts_label=label,
                              expected_artifact_path=expected_ts_artifact_path(arc_dir, label),
                              reason='Queued to ARC.', coefficient=-9.5e-05, delta_ln_k=0.795)
        path_reaction = reactions_by_ts[ts_label][0]
        _write_artifact(record.expected_artifact_path,
                        log_text=_channel_marker(path_reaction.reactants,
                                                 path_reaction.products) + '\n')
        _write_status(arc_dir, label, converged=True)
        records.append(record)
    _write_energy_settings(arc_dir)
    capture_ts_artifacts(
        records, arc_dir, os.path.join(prior_project_dir, 'capture'),
        networks={network_id: {'source_path': source_path,
                               'source_sha256': hash_file(source_path), 'method': 'MSC'}},
        sensitivity_by_ts={record.key: (record.coefficient, record.delta_ln_k)
                           for record in records})


class _ArtifactWritingFakeARC(object):
    """Stands in for arc.main.ARC at exactly the cluster boundary: execute() writes the statmech
    artifacts, convergence statuses, and energy settings a real ARC run would leave in the round's
    own project directory. ``converge_plan`` names one set of CHANNELS (as _channel_marker
    strings) per execute() call -- never labels, which renumber between rounds -- and every
    converged artifact's log file records its channel marker."""

    converge_plan = ()
    constructions = []
    executions = 0

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        _ArtifactWritingFakeARC.constructions.append(kwargs)

    def as_dict(self):
        return dict(self.kwargs)

    def execute(self):
        arc_dir = self.kwargs['project_directory']
        to_converge = _ArtifactWritingFakeARC.converge_plan[_ArtifactWritingFakeARC.executions]
        _ArtifactWritingFakeARC.executions += 1
        _write_energy_settings(arc_dir)
        for reaction in self.kwargs['reactions']:
            label = reaction['ts_label']
            channel = _channel_marker(reaction['reactants'], reaction['products'])
            if channel in to_converge:
                _write_artifact(expected_ts_artifact_path(arc_dir, label),
                                log_text=channel + '\n')
                _write_status(arc_dir, label, converged=True)
            else:
                _write_status(arc_dir, label, converged=False)


def _hybrid_channels(hybrid_path):
    """Read a hybrid file's reaction blocks by CHANNEL (never by label): channel pair ->
    ``(local_ts_label, is_qm)``, where a QM'd channel's block has dropped its
    ``kinetics = Arrhenius(...)`` line. Scanned textually because the hybrid's consumer is Arkane
    (which accepts its positional ``transitionState(...)`` calls), not ``t3.pdep.parser``."""
    with open(hybrid_path) as f:
        content = f.read()
    channels = dict()
    for block in re.split(r'\breaction\(', content)[1:]:
        block = block.split('\n\n')[0]
        # Greedy-to-end-of-line bracket matching: species labels themselves contain ']' (e.g.
        # '[O]C([O])=O(8160)'), so a lazy or ]-excluding match truncates the list.
        reactants = re.search(r"reactants = \[(.*)\],", block).group(1)
        products = re.search(r"products = \[(.*)\],", block).group(1)
        ts_label = re.search(r"transitionState = '([^']*)'", block).group(1)
        sides = tuple(sorted(
            ' + '.join(sorted(part.strip().strip("'") for part in side.split(',')))
            for side in (reactants, products)))
        channels[sides] = (ts_label, 'kinetics = Arrhenius(' not in block)
    return channels


def _hybrid_marker(hybrid_path, ts_label):
    """The channel marker inside the hybrid's own vendored log for ``ts_label`` -- the bytes that
    prove WHICH barrier was folded under that label."""
    log_dir = os.path.join(os.path.dirname(hybrid_path), 'qm', 'logs', ts_label)
    (log_name,) = os.listdir(log_dir)
    with open(os.path.join(log_dir, log_name)) as f:
        return f.read().strip()


def _fixture_explorer(monkeypatch, project_directory, network_id):
    """The canonical fake explorer of this module: round 0's 'exploration' deposits the pristine
    fixture; every later round deposits the RENUMBERED variant (_renumbered_fixture_text), so any
    label-keyed carry across the round boundary is exposed rather than trivially stable."""
    explore_calls = []

    def _fake_explore(*, network_path, config, logger=None):
        explore_calls.append(network_path)
        round_index = len(explore_calls) - 1
        paths = round_paths(project_directory, round_index)
        dest_path = _explored_network_path(paths, network_id)
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        if round_index == 0:
            shutil.copyfile(_FIXTURE_NETWORK_PATH, dest_path)
        else:
            with open(dest_path, 'w') as f:
                f.write(_renumbered_fixture_text())
        return PDepExplorationResult(network_id=network_id,
                                     status=EXPLORATION_STATUS_SUCCEEDED,
                                     network_paths=(dest_path,))

    monkeypatch.setattr('t3.pdep.pes_loop.explore_pdep_network', _fake_explore)
    return explore_calls


_CH1_MARKER = _channel_marker(('O-2(13598)', 'CO2(11)'), ('O=C1OO1(21648)',))
_CH2_MARKER = _channel_marker(('O=C1OO1(21648)',), ('[O]C([O])=O(8160)',))


def test_real_loop_real_capture_keeps_qm_on_the_right_channel_across_a_renumber(tmp_path,
                                                                                monkeypatch):
    """Against the REAL capture_ts_artifacts, across a REAL renumbering: round 0 (pristine
    labels) converges channel 1 (then labeled 'TS1'); round 1 explores the RENUMBERED network, in
    which that same channel is labeled 'TS2' -- the loop must skip it under its new label, queue
    channel 2 (now labeled 'TS1', which a label-keyed carry would have wrongly skipped as
    computed), and fold round 0's barrier under the new label with the artifact bytes proving the
    channel. Channel 3 measures a ln(k) response of exactly 0.0 in the live ME SA -- below
    qm.min_delta_ln_k -- so it is never queued, and the loop converges at round 2."""
    project_directory = os.path.join(str(tmp_path), 'loop_project')
    os.makedirs(project_directory)
    seed_path = os.path.join(str(tmp_path), 'network0_full.py')
    shutil.copyfile(_FIXTURE_NETWORK_PATH, seed_path)

    _fixture_explorer(monkeypatch, project_directory, 'network0_full')
    _ArtifactWritingFakeARC.converge_plan = ({_CH1_MARKER}, {_CH2_MARKER})
    _ArtifactWritingFakeARC.constructions = []
    _ArtifactWritingFakeARC.executions = 0
    monkeypatch.setattr(pes_qm, 'ARC', _ArtifactWritingFakeARC)

    config = PESLoopConfig(
        pes={'network': seed_path, 'source': ['O-2(13598)', 'CO2(11)'],
             'bath_gas': {'Ar': 1.0}, 'method': 'MSC'},
        termination={'max_rounds': 4},
    )
    result = run_pes_loop(config, project_directory=project_directory, qm_runner=arc_qm_runner)

    assert result.status == 'converged', f'{result.status}: {result.reason}'
    assert [record.status for record in result.rounds] == ['continuing', 'continuing', 'converged']
    # Round 0 (pristine labels): channels 1 and 2 queued; channel 3 gated by its measured zero.
    assert sorted(result.rounds[0].queued_ts_labels) == ['TS1', 'TS2']
    # Round 1 (renumbered): ONLY the still-uncomputed channel 2 queued -- under its NEW label
    # 'TS1'. A label-keyed carry would have skipped 'TS1' as already computed and re-queued 'TS2'.
    assert result.rounds[1].queued_ts_labels == ('TS1',)
    for record in result.rounds:
        floor_skips = [s for s in record.skipped if 'below the min_delta_ln_k floor' in s.reason]
        assert len(floor_skips) == 1

    round0 = _hybrid_channels(hybrid_network_path(round_paths(project_directory, 0),
                                                  'network0_full'))
    round1_hybrid_path = hybrid_network_path(round_paths(project_directory, 1), 'network0_full')
    round1 = _hybrid_channels(round1_hybrid_path)
    # Round 0 hybrid (pristine labels): channel 1 is QM as 'TS1'; channels 2, 3 are ILT.
    assert round0[_CH1] == ('TS1', True)
    assert round0[_CH2] == ('TS2', False)
    assert round0[_CH3] == ('TS3', False)
    # Round 1 hybrid (renumbered labels): cumulative QM, each channel under its NEW label.
    assert round1[_CH1] == ('TS2', True)
    assert round1[_CH2] == ('TS1', True)
    assert round1[_CH3] == ('TS3', False)
    # The bytes prove the attribution: under round 1's 'TS2' sits round 0's channel-1 barrier,
    # and under 'TS1' sits this round's channel-2 barrier -- the right barrier on the right
    # channel, across the renumber.
    assert _hybrid_marker(round1_hybrid_path, 'TS2') == _CH1_MARKER
    assert _hybrid_marker(round1_hybrid_path, 'TS1') == _CH2_MARKER

    # Both rounds' captures are REAL and re-verify cleanly, sensitivity evidence included --
    # exactly what defect 1 made impossible.
    verified_artifacts = 0
    for round_index in (0, 1):
        verified = verify_capture(round_paths(project_directory, round_index).capture)
        # Loop-written manifests carry the aggregation marker, so their ungated all-directions
        # evidence can never be silently compared against T3's in-run (observable-gated) values.
        assert verified.sensitivity_aggregation == 'all_directions_max_abs'
        for record in verified.ts_records:
            if record.artifact_path is not None:
                verified_artifacts += 1
                assert record.coefficient is not None and record.delta_ln_k is not None
    assert verified_artifacts == 2


def test_real_loop_round_0_full_adoption_completes_with_real_vendoring(tmp_path, monkeypatch):
    """Defect 2 end to end, with adoption made maximally hostile to label matching: the prior
    project's vendored network carries DIFFERENT species label indices and different TS positions
    ('TS1'<->'TS3') -- only the structures match. All three channels adopt, round 0 queues
    nothing, never touches ARC or capture, and still writes a hybrid carrying all three channels
    as QM (each barrier's bytes on its own channel); round 1 then converges."""
    prior_project_dir = os.path.join(str(tmp_path), 'prior_project')
    _build_prior_capture(prior_project_dir, network_text=_renamed_prior_network_text(),
                         ts_labels=('TS1', 'TS2', 'TS3'))
    project_directory = os.path.join(str(tmp_path), 'loop_project')
    os.makedirs(project_directory)
    seed_path = os.path.join(str(tmp_path), 'network0_full.py')
    shutil.copyfile(_FIXTURE_NETWORK_PATH, seed_path)

    _fixture_explorer(monkeypatch, project_directory, 'network0_full')
    _ArtifactWritingFakeARC.converge_plan = ()
    _ArtifactWritingFakeARC.constructions = []
    _ArtifactWritingFakeARC.executions = 0
    monkeypatch.setattr(pes_qm, 'ARC', _ArtifactWritingFakeARC)

    config = PESLoopConfig(
        pes={'network': seed_path, 'source': ['O-2(13598)', 'CO2(11)'],
             'bath_gas': {'Ar': 1.0}, 'method': 'MSC'},
        termination={'max_rounds': 2},
        reuse={'from_t3_projects': [prior_project_dir]},
    )
    result = run_pes_loop(config, project_directory=project_directory, qm_runner=arc_qm_runner)

    assert result.status == 'converged', f'{result.status}: {result.reason}'
    assert [record.status for record in result.rounds] == ['continuing', 'converged']
    # Nothing was ever queued, so ARC must never even have been constructed.
    assert _ArtifactWritingFakeARC.constructions == []
    round0_hybrid_path = hybrid_network_path(round_paths(project_directory, 0), 'network0_full')
    round0 = _hybrid_channels(round0_hybrid_path)
    assert round0[_CH1] == ('TS1', True)
    assert round0[_CH2] == ('TS2', True)
    assert round0[_CH3] == ('TS3', True)
    # Each adopted barrier's bytes sit on their own channel, matched purely structurally: the
    # prior project's markers carry ITS species labels, so equality here proves the artifact came
    # from the structurally-matching prior channel, not from any label agreement.
    assert _hybrid_marker(round0_hybrid_path, 'TS1') \
        == _channel_marker(('O-2(5)', 'CO2(7)'), ('O=C1OO1(9)',))
    assert _hybrid_marker(round0_hybrid_path, 'TS2') \
        == _channel_marker(('O=C1OO1(9)',), ('[O]C([O])=O(2)',))
    assert _hybrid_marker(round0_hybrid_path, 'TS3') \
        == _channel_marker(('O=C1OO1(9)',), ('[O]O[C]=O(3)',))
    # The hybrid's energy settings came from the adopted artifacts' own prior manifest -- the
    # _adopted_energy_settings path a capture-less round is forced through.
    with open(round0_hybrid_path) as f:
        assert 'wb97xd2023' in f.read()


def test_pes_cli_main_drives_the_real_loop_end_to_end(tmp_path, monkeypatch):
    """``PES.main()`` for real, from a real input file: the only things faked are the two external
    engines (Arkane's explorer and ARC). The real ``run_pes_loop``, the real ``arc_qm_runner`` and
    the real ``capture_ts_artifacts`` all execute, and the assertions are on artifacts only
    production writes -- the hybrid network under the input file's own directory (the
    ``project_directory`` default reaching ``round_paths``), the path the explorer was handed (the
    relative ``pes.network`` resolved against the input file), and the round count
    (``termination.max_rounds`` parsed out of the YAML, 2, governing the run rather than the
    schema default of 5)."""
    project_directory = str(tmp_path)
    input_path = os.path.join(project_directory, 'input.yml')
    shutil.copyfile(_FIXTURE_NETWORK_PATH, os.path.join(project_directory, 'network0_full.py'))
    with open(input_path, 'w') as f:
        # 'network' is deliberately RELATIVE: the CLI must resolve it against the input file's
        # directory, or the explorer never sees a readable seed.
        yaml.dump({'pes': {'network': 'network0_full.py',
                           'source': ['O-2(13598)', 'CO2(11)'],
                           'bath_gas': {'Ar': 1.0},
                           'method': 'MSC'},
                   'termination': {'max_rounds': 2}}, f)

    explore_calls = _fixture_explorer(monkeypatch, project_directory, 'network0_full')
    _ArtifactWritingFakeARC.converge_plan = ({_CH1_MARKER}, {_CH2_MARKER})
    _ArtifactWritingFakeARC.constructions = []
    _ArtifactWritingFakeARC.executions = 0
    monkeypatch.setattr(pes_qm, 'ARC', _ArtifactWritingFakeARC)
    monkeypatch.setattr(sys, 'argv', ['PES.py', input_path])

    result = PES.main()

    # What the explorer was actually HANDED for round 0: the relative 'network0_full.py' of the
    # input file, resolved to an absolute path against that file's own directory. Asserted on the
    # engine's argument rather than on the run's success, because the fake explorer would have
    # succeeded on an unresolved path too -- and the real one would not.
    assert explore_calls[0] == os.path.join(project_directory, 'network0_full.py')

    # max_rounds: 2 from the input file, not the schema's default of 5.
    assert result.status == 'max_rounds', f'{result.status}: {result.reason}'
    assert len(result.rounds) == 2
    assert _ArtifactWritingFakeARC.executions == 2
    assert PES.exit_code_for(result.status) == 0

    # The fake explorer makes its own round directories under this same path, so these two say
    # only that the run got as far as two rounds -- what proves the CLI's project_directory
    # default reached round_paths is the hybrid below, which production alone writes. The absence
    # of round_2 does reflect production: it is the round budget, via the explorer call count.
    assert os.path.isdir(os.path.join(project_directory, 'round_0'))
    assert os.path.isdir(os.path.join(project_directory, 'round_1'))
    assert not os.path.isdir(os.path.join(project_directory, 'round_2'))

    # Ruling 6: the final status, its reason, and the diagram reached the log, and the run was
    # closed out -- not merely returned to a caller who might never print it.
    with open(os.path.join(project_directory, 't3.log')) as f:
        log_text = f.read()
    assert "terminated with status 'max_rounds' after 2 round(s)." in log_text
    assert f'Reason: {result.reason}' in log_text
    assert f'PES diagram: {result.final_diagram_path}' in log_text
    assert f'Final network: {result.final_network_path}' in log_text
    assert 'Total T3 execution time' in log_text

    round1_hybrid_path = hybrid_network_path(round_paths(project_directory, 1), 'network0_full')
    assert os.path.isfile(round1_hybrid_path)
    round1 = _hybrid_channels(round1_hybrid_path)
    # Round 1 explored the RENUMBERED network, so channel 1 -- QM'd in round 0 -- carries the new
    # label 'TS2' here, and this round's channel 2 carries 'TS1'.
    assert round1[_CH1] == ('TS2', True)
    assert round1[_CH2] == ('TS1', True)
    assert round1[_CH3] == ('TS3', False)
    # The bytes prove which barrier landed on which channel through the CLI's own run.
    assert _hybrid_marker(round1_hybrid_path, 'TS2') == _CH1_MARKER
    assert _hybrid_marker(round1_hybrid_path, 'TS1') == _CH2_MARKER


def test_pes_cli_project_directory_flag_moves_the_whole_run(tmp_path, monkeypatch):
    """``-p`` end to end: the run must land in the given directory and NOTHING may land beside the
    input file. Asserted on the hybrid network and the log -- artifacts only production writes --
    because the fake explorer creates round directories under whatever path it was handed and so
    cannot distinguish the flag from the default."""
    input_directory = os.path.join(str(tmp_path), 'inputs')
    run_directory = os.path.join(str(tmp_path), 'elsewhere', 'run')
    os.makedirs(input_directory)
    input_path = os.path.join(input_directory, 'input.yml')
    shutil.copyfile(_FIXTURE_NETWORK_PATH, os.path.join(input_directory, 'network0_full.py'))
    with open(input_path, 'w') as f:
        yaml.dump({'pes': {'network': 'network0_full.py',
                           'source': ['O-2(13598)', 'CO2(11)'],
                           'bath_gas': {'Ar': 1.0},
                           'method': 'MSC'},
                   'termination': {'max_rounds': 1}}, f)

    _fixture_explorer(monkeypatch, run_directory, 'network0_full')
    _ArtifactWritingFakeARC.converge_plan = ({_CH1_MARKER},)
    _ArtifactWritingFakeARC.constructions = []
    _ArtifactWritingFakeARC.executions = 0
    monkeypatch.setattr(pes_qm, 'ARC', _ArtifactWritingFakeARC)
    monkeypatch.setattr(sys, 'argv', ['PES.py', input_path, '-p', run_directory])

    result = PES.main()

    assert result.status == 'max_rounds', f'{result.status}: {result.reason}'
    # The hybrid network -- written by arc_qm_runner, at a path derived from the project_directory
    # the LOOP was given -- landed under the -p directory, which the CLI also had to create.
    round0_hybrid_path = hybrid_network_path(round_paths(run_directory, 0), 'network0_full')
    assert os.path.isfile(round0_hybrid_path)
    assert _hybrid_channels(round0_hybrid_path)[_CH1] == ('TS1', True)
    assert os.path.isfile(os.path.join(run_directory, 't3.log'))
    # Nothing whatsoever landed beside the input file.
    assert sorted(os.listdir(input_directory)) == ['input.yml', 'network0_full.py']


def test_pes_cli_resolves_a_relative_reuse_path_against_the_input_file(tmp_path, monkeypatch):
    """``reuse.from_t3_projects`` is anchored exactly like ``pes.network``: relative to the input
    file. Unresolved it fails silently rather than loudly -- the prior project is simply not found,
    nothing adopts, and the loop pays for it by queueing real quantum chemistry it already had.
    Here the prior capture holds all three channels, so a correctly resolved path means ARC is
    never even constructed."""
    prior_project_dir = os.path.join(str(tmp_path), 'prior_project')
    _build_prior_capture(prior_project_dir, network_text=_renamed_prior_network_text(),
                         ts_labels=('TS1', 'TS2', 'TS3'))
    project_directory = os.path.join(str(tmp_path), 'loop_project')
    os.makedirs(project_directory)
    input_path = os.path.join(project_directory, 'input.yml')
    shutil.copyfile(_FIXTURE_NETWORK_PATH, os.path.join(project_directory, 'network0_full.py'))
    with open(input_path, 'w') as f:
        yaml.dump({'pes': {'network': 'network0_full.py',
                           'source': ['O-2(13598)', 'CO2(11)'],
                           'bath_gas': {'Ar': 1.0},
                           'method': 'MSC'},
                   'termination': {'max_rounds': 2},
                   'reuse': {'from_t3_projects': ['../prior_project']}}, f)

    _fixture_explorer(monkeypatch, project_directory, 'network0_full')
    _ArtifactWritingFakeARC.converge_plan = ()
    _ArtifactWritingFakeARC.constructions = []
    _ArtifactWritingFakeARC.executions = 0
    monkeypatch.setattr(pes_qm, 'ARC', _ArtifactWritingFakeARC)
    monkeypatch.setattr(sys, 'argv', ['PES.py', input_path])

    result = PES.main()

    assert result.status == 'converged', f'{result.status}: {result.reason}'
    # All three channels adopted from the prior project, so nothing was ever queued: an
    # unresolved reuse path would have adopted nothing and run a real ARC round instead.
    assert _ArtifactWritingFakeARC.constructions == []
    round0_hybrid_path = hybrid_network_path(round_paths(project_directory, 0), 'network0_full')
    round0 = _hybrid_channels(round0_hybrid_path)
    assert [round0[channel][1] for channel in (_CH1, _CH2, _CH3)] == [True, True, True]
    # Ruling 3: the Logger main() built was actually handed to the loop -- this line is the
    # LOOP's own, not the CLI's, so it can only be in the file if logger= reached run_pes_loop.
    with open(os.path.join(project_directory, 't3.log')) as f:
        assert "PES loop: reusing 3 prior QM result(s)" in f.read()
