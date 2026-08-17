#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_explorer_arkane module

Tests the concrete ``ArkaneExplorerAdapter``. The Arkane subprocess itself is never invoked:
``t3.pdep.explorer.arkane.run_arkane_job`` is monkeypatched with a fake that writes whichever
``pdep/final/network<i>_(full|reduced).py`` and ``output*.py``/``stdout.log``/``arkane.log``
fixture content a given test needs, mirroring the ``monkeypatch`` idiom used in
``tests/test_pdep/test_arkane_mesolver.py``. All fixture trees are hand-built directly in this
file (mirroring ``tests/test_pdep/test_explorer_input_file.py``'s ``SOURCE_NO_KINETICS``-style
constants) rather than reused from ``tests/data/pdep_me/``, which models ME-solver output only and
has no ``pdep/final/`` layout at all.
"""

import os

import pytest

from t3.pdep.explorer.factory import _registered_explorer_adapters

# A minimal RMG P-dep network fixture, structurally identical to
# ``tests/test_pdep/test_explorer_input_file.py::SOURCE_NO_KINETICS`` (same bath gas/species/
# reaction shape), used as the ``network_path`` every adapter under test is built from.
NETWORK_SOURCE = """species(
    label='methoxy',
    structure=SMILES('[CH2]O'),
    E0=(0, 'kJ/mol'),
)
species(
    label='CH2O',
    structure=SMILES('C=O'),
    E0=(-50, 'kJ/mol'),
)
species(
    label='H',
    structure=SMILES('[H]'),
    E0=(0, 'kJ/mol'),
)
species(
    label='He',
    structure=SMILES('[He]'),
    E0=(0, 'kJ/mol'),
    reactive=False,
)
transitionState(
    label='TS1',
    E0=(50, 'kJ/mol'),
)
reaction(
    label='CH2O+H=methoxy',
    reactants=['CH2O', 'H'],
    products=['methoxy'],
    transitionState='TS1',
)
network(
    label='PDepNetwork #1',
    isomers=['methoxy'],
    reactants=[('CH2O', 'H')],
    bathGas={'He': 1.0},
)
pressureDependence(
    label='PDepNetwork #1',
    Tmin=(300, 'K'),
    Tmax=(2000, 'K'),
    Tcount=8,
    Pmin=(0.01, 'bar'),
    Pmax=(100, 'bar'),
    Pcount=5,
    maximumGrainSize=(0.5, 'kcal/mol'),
    maximumGrainCount=250,
    method = 'modified strong collision',
    activeKRotor=True,
    activeJRotor=True,
    rmgmode=False,
)
"""

# A valid, finite Chebyshev pdepreaction(...) payload -- a genuine ME success (mirrors
# ``tests/data/pdep_me/success/output.py`` in shape, hand-built here so this file has no
# dependency on that fixture's directory layout).
VALID_OUTPUT_PY = """pdepreaction(
    reactants = ['CH2O', 'H'],
    products = ['methoxy'],
    kinetics = Chebyshev(
        coeffs = [
            [-31.1515, 0.88105, -0.176824, -0.00220675],
            [33.9767, 0.317263, 0.0556114, -0.0172293],
        ],
        kunits = 's^-1',
        Tmin = (300, 'K'),
        Tmax = (2000, 'K'),
        Pmin = (0.01, 'bar'),
        Pmax = (100, 'bar'),
    ),
)
"""

# A well-formed FINAL NETWORK artifact, i.e. what Arkane actually copies into
# ``pdep/final/network<i>_(full|reduced).py``. That file is a serialized network
# (``PressureDependenceJob.save_input_file``, arkane/pdep.py:654-779), NOT an ME output file: it
# carries species(...) / transitionState(...) / reaction(...) / network(...) /
# pressureDependence(...) and never a pdepreaction(...). ``NETWORK_SOURCE`` has exactly that
# shape, so it doubles as a realistic final artifact.
VALID_NETWORK_PY = NETWORK_SOURCE

# A syntactically valid but entirely-None Chebyshev payload -- the soft-CSE-failure archetype:
# exit 0, empty stderr, yet no usable rate coefficients.
SOFT_FAILURE_OUTPUT_PY = """pdepreaction(
    reactants = ['CH2O', 'H'],
    products = ['methoxy'],
    kinetics = Chebyshev(
        coeffs = [
            [None, None],
            [None, None],
        ],
        kunits = 's^-1',
        Tmin = (300, 'K'),
        Tmax = (2000, 'K'),
        Pmin = (0.01, 'bar'),
        Pmax = (100, 'bar'),
    ),
)
"""

# A richer network fixture with two isomers ('A', 'B'), one reactant channel (('R1', 'R2')) and
# one product channel derived from the path reactions (('P1', 'P2')). Used to test rule 5b's
# net-reaction-count and channel-coverage checks, which need more than the single net reaction
# ``NETWORK_SOURCE`` produces to have anything meaningful to truncate/pad/starve.
# ``PDepNetwork.expected_net_reaction_count()`` for this network is 6 (verified with
# ``t3.pdep.parser.parse_pdep_network_text`` directly, not hand-derived): the six net reactions it
# implies are, canonically, B->A, (R1,R2)->A, (R1,R2)->B, A->(P1,P2), B->(P1,P2) and
# (R1,R2)->(P1,P2).
MULTI_NETWORK_PY = """species(
    label='A',
    structure=SMILES('[CH2]O'),
    E0=(0, 'kJ/mol'),
)
species(
    label='B',
    structure=SMILES('C=O'),
    E0=(-10, 'kJ/mol'),
)
species(
    label='R1',
    structure=SMILES('[H]'),
    E0=(0, 'kJ/mol'),
)
species(
    label='R2',
    structure=SMILES('C'),
    E0=(0, 'kJ/mol'),
)
species(
    label='P1',
    structure=SMILES('O'),
    E0=(0, 'kJ/mol'),
)
species(
    label='P2',
    structure=SMILES('N'),
    E0=(0, 'kJ/mol'),
)
species(
    label='He',
    structure=SMILES('[He]'),
    E0=(0, 'kJ/mol'),
    reactive=False,
)
transitionState(
    label='TS1',
    E0=(50, 'kJ/mol'),
)
transitionState(
    label='TS2',
    E0=(60, 'kJ/mol'),
)
transitionState(
    label='TS3',
    E0=(70, 'kJ/mol'),
)
reaction(
    label='R1+R2=A',
    reactants=['R1', 'R2'],
    products=['A'],
    transitionState='TS1',
)
reaction(
    label='A=B',
    reactants=['A'],
    products=['B'],
    transitionState='TS2',
)
reaction(
    label='B=P1+P2',
    reactants=['B'],
    products=['P1', 'P2'],
    transitionState='TS3',
)
network(
    label='PDepNetwork #2',
    isomers=['A', 'B'],
    reactants=[('R1', 'R2')],
    bathGas={'He': 1.0},
)
pressureDependence(
    label='PDepNetwork #2',
    Tmin=(300, 'K'),
    Tmax=(2000, 'K'),
    Tcount=8,
    Pmin=(0.01, 'bar'),
    Pmax=(100, 'bar'),
    Pcount=5,
    maximumGrainSize=(0.5, 'kcal/mol'),
    maximumGrainCount=250,
    method = 'modified strong collision',
    activeKRotor=True,
    activeJRotor=True,
    rmgmode=False,
)
"""


def _pdepreaction(reactants, products):
    """Build one finite-Arrhenius ``pdepreaction(...)`` text block for the given channel pair."""
    return (f"pdepreaction(\n"
            f"    reactants = {list(reactants)!r},\n"
            f"    products = {list(products)!r},\n"
            f"    kinetics = Arrhenius(\n"
            f"        A = (1.0e+13, 's^-1'),\n"
            f"        n = 0.0,\n"
            f"        Ea = (10.0, 'kJ/mol'),\n"
            f"        T0 = (1, 'K'),\n"
            f"    ),\n"
            f")\n")


# The six net reactions ``MULTI_NETWORK_PY`` implies, in the order Arkane's own enumeration would
# print them (verified against ``expected_net_reaction_count()``'s internal enumeration, not
# hand-derived).
_MULTI_NETWORK_NET_REACTIONS = (
    (('B',), ('A',)),
    (('R1', 'R2'), ('A',)),
    (('R1', 'R2'), ('B',)),
    (('A',), ('P1', 'P2')),
    (('B',), ('P1', 'P2')),
    (('R1', 'R2'), ('P1', 'P2')),
)

# Exactly the six net reactions the network implies: the over-refusal guard case that must accept.
MULTI_OUTPUT_EXACT_PY = ''.join(_pdepreaction(r, p) for r, p in _MULTI_NETWORK_NET_REACTIONS)

# Only five of the six -- a truncated output.py, still containment-clean (every species it names is
# still declared by the network) and still touching every isomer/channel, so ONLY the count check
# can catch it.
MULTI_OUTPUT_TRUNCATED_PY = ''.join(_pdepreaction(r, p) for r, p in _MULTI_NETWORK_NET_REACTIONS[:-1])

# All six plus one duplicate -- a surplus, which is just as much evidence the pairing is not from
# one run as a shortfall is.
MULTI_OUTPUT_SURPLUS_PY = MULTI_OUTPUT_EXACT_PY + _pdepreaction(*_MULTI_NETWORK_NET_REACTIONS[0])

# Six entries total (the count check alone would accept this), but none of them touch isomer 'B'
# on either side: the three entries that do are each replaced by a duplicate of one that does not,
# so the total stays 6 while 'B' goes completely untouched.
MULTI_OUTPUT_CHANNEL_GAP_PY = ''.join(
    _pdepreaction(r, p) for r, p in (_MULTI_NETWORK_NET_REACTIONS[1],
                                     _MULTI_NETWORK_NET_REACTIONS[3],
                                     _MULTI_NETWORK_NET_REACTIONS[5])
) * 2

# Five of the six net reactions, dropping index 0 (('B',), ('A',)) instead of the last entry. The
# remaining five still touch every declared channel -- A via indices 1 and 3, B via indices 2 and 4,
# (R1, R2) via indices 1, 2 and 5, (P1, P2) via indices 3, 4 and 5 -- so channel coverage alone would
# accept this; only the exact-count check (isolated from the >, mismatch-blind-to-shortfall mutant)
# can catch it.
MULTI_OUTPUT_TRUNCATED_DROP_FIRST_PY = ''.join(
    _pdepreaction(r, p) for r, p in _MULTI_NETWORK_NET_REACTIONS[1:]
)

# Six entries total (right count), each a duplicate of one of the three source-to-source net
# reactions (indices 0, 1, 2), so 'A', 'B' and (R1, R2) are all touched but (P1, P2) -- the declared
# product channel -- never is. Count check alone would accept this; only channel coverage that
# includes ``network.product_channels`` can catch it.
MULTI_OUTPUT_MISSING_PRODUCT_CHANNEL_PY = ''.join(
    _pdepreaction(r, p) for r, p in (_MULTI_NETWORK_NET_REACTIONS[0],
                                     _MULTI_NETWORK_NET_REACTIONS[1],
                                     _MULTI_NETWORK_NET_REACTIONS[2],
                                     _MULTI_NETWORK_NET_REACTIONS[0],
                                     _MULTI_NETWORK_NET_REACTIONS[1],
                                     _MULTI_NETWORK_NET_REACTIONS[2])
)

# The forgery neither the count check nor channel coverage can see: six entries (right count) in
# which the (R1, R2) <-> (P1, P2) pair is replaced by a second copy of A <-> (P1, P2). Every declared
# channel is still touched -- 'A' by indices 0, 1, 3, 'B' by 0, 2, 4, (R1, R2) by 1, 2 and (P1, P2)
# by 3, 4 -- and the total is still 6, so ONLY comparing the expected channel-PAIR set catches it.
MULTI_OUTPUT_SWAPPED_PAIR_PY = ''.join(
    _pdepreaction(r, p) for r, p in (_MULTI_NETWORK_NET_REACTIONS[0],
                                     _MULTI_NETWORK_NET_REACTIONS[1],
                                     _MULTI_NETWORK_NET_REACTIONS[2],
                                     _MULTI_NETWORK_NET_REACTIONS[3],
                                     _MULTI_NETWORK_NET_REACTIONS[4],
                                     _MULTI_NETWORK_NET_REACTIONS[3])
)

# Exactly the six expected net reactions, but one written in the opposite direction. Arkane's own
# duplicate suppression compares by label in EITHER direction, and which side the solver reports as
# the reactant side is not part of the network's identity, so this must still be accepted.
MULTI_OUTPUT_REVERSED_ENTRY_PY = ''.join(
    _pdepreaction(r, p) for r, p in (_MULTI_NETWORK_NET_REACTIONS[0],
                                     _MULTI_NETWORK_NET_REACTIONS[1],
                                     _MULTI_NETWORK_NET_REACTIONS[2],
                                     (_MULTI_NETWORK_NET_REACTIONS[3][1],
                                      _MULTI_NETWORK_NET_REACTIONS[3][0]),
                                     _MULTI_NETWORK_NET_REACTIONS[4],
                                     _MULTI_NETWORK_NET_REACTIONS[5])
)


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(text)


def _write_network_source(tmp_path):
    network_path = os.path.join(str(tmp_path), 'network_source.py')
    _write(network_path, NETWORK_SOURCE)
    return network_path


def _make_fake_run_arkane_job(*, success=True, output_files=None, final_files=None,
                              stdout_lines=None, arkane_log_lines=None, stderr_lines=None):
    """
    Build a fake ``run_arkane_job`` writing the given fixture content into ``output_directory``.

    Args:
        success (bool): The value the fake reports back (mirrors Arkane's own exit-code-derived
                        bool).
        output_files (dict, optional): Maps a file name (e.g. ``'output.py'``, ``'output0.py'``)
                                       to the text to write directly under ``output_directory``.
        final_files (dict, optional): Maps a file name (e.g. ``'network0_full.py'``) to the text
                                      to write under ``output_directory/pdep/final/``.
        stdout_lines (list, optional): Lines to write to ``output_directory/stdout.log``.
        arkane_log_lines (list, optional): Lines to write to ``output_directory/arkane.log``.
        stderr_lines (list, optional): Lines to write to ``output_directory/stderr.log``.
    """
    output_files = output_files or {}
    final_files = final_files or {}

    def _fake(input_file, output_directory, plot=False, logger=None, required_artifact='output.py',
              timeout=None):
        for name, content in output_files.items():
            _write(os.path.join(output_directory, name), content)
        for name, content in final_files.items():
            _write(os.path.join(output_directory, 'pdep', 'final', name), content)
        if stdout_lines:
            _write(os.path.join(output_directory, 'stdout.log'), '\n'.join(stdout_lines) + '\n')
        if arkane_log_lines:
            _write(os.path.join(output_directory, 'arkane.log'), '\n'.join(arkane_log_lines) + '\n')
        if stderr_lines:
            _write(os.path.join(output_directory, 'stderr.log'), '\n'.join(stderr_lines) + '\n')
        # Mirror the real run_arkane_job's own gate (t3/runners/rmg_runner.py): a caller-declared
        # required_artifact must actually exist for the run to be reported as successful, exactly
        # as the real runner requires Arkane itself to have (re)written it. Reporting ``success``
        # regardless of the artifact's presence would be a lie the real runner never tells, and
        # would make the job-failure signal that depends on it untestable.
        artifact_path = os.path.join(output_directory, required_artifact)
        return success and os.path.isfile(artifact_path)

    return _fake


def _build_adapter(tmp_path, monkeypatch, output_directory=None, **kwargs):
    from t3.pdep.explorer.arkane import ArkaneExplorerAdapter
    network_path = _write_network_source(tmp_path)
    output_directory = output_directory or os.path.join(str(tmp_path), 'explorer_run')
    defaults = dict(
        seed_species=('methoxy',),
        output_directory=output_directory,
        network_path=network_path,
        method='CSE',
        bath_gas={'He': 1.0},
    )
    defaults.update(kwargs)
    return ArkaneExplorerAdapter(**defaults)


class TestArkaneExplorerAdapterRegistration:

    def test_registered_under_arkane(self):
        from t3.pdep.explorer.arkane import ArkaneExplorerAdapter
        assert _registered_explorer_adapters.get('Arkane') is ArkaneExplorerAdapter


class TestArkaneExplorerAdapterTimeout:
    """
    The adapter's ``timeout`` is a passthrough to ``run_arkane_job``, and a deadline overrun is a
    RECORDED, distinguishable failure (issue #183: a hung Arkane used to hang the whole campaign,
    with no kill path anywhere between the adapter and ARC's un-timeboxed subprocess call).
    """

    def test_timeout_is_forwarded_to_run_arkane_job(self, tmp_path, monkeypatch):
        """The configured deadline must reach the runner that enforces it, verbatim."""
        captured = {}
        inner = _make_fake_run_arkane_job(success=True)

        def _fake(**kwargs):
            captured.update(kwargs)
            return inner(**{key: value for key, value in kwargs.items() if key != 'timeout'})

        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job', _fake)
        adapter = _build_adapter(tmp_path, monkeypatch, timeout=42.5)
        adapter.explore()
        assert captured['timeout'] == 42.5

    def test_timeout_defaults_to_none(self, tmp_path, monkeypatch):
        """With no timeout configured, the runner must receive None (the historical, unbounded
        in-process call), not some invented default deadline."""
        captured = {}
        inner = _make_fake_run_arkane_job(success=True)

        def _fake(**kwargs):
            captured.update(kwargs)
            return inner(**{key: value for key, value in kwargs.items() if key != 'timeout'})

        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job', _fake)
        adapter = _build_adapter(tmp_path, monkeypatch)
        adapter.explore()
        assert captured['timeout'] is None

    def test_timed_out_run_is_a_recorded_distinguishable_failure(self, tmp_path, monkeypatch):
        """A deadline overrun surfaces as succeeded=False with a reason naming the timeout --
        never as an exception (which would destroy the run record), and never folded into the
        generic 'Arkane reported job failure' diagnosis."""
        from t3.runners.rmg_runner import ArkaneJobResult

        timeout_reason = 'Arkane run in /x timed out after 7 s; its process group was killed.'
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            lambda **kwargs: ArkaneJobResult(succeeded=False, timed_out=True, reason=timeout_reason))
        adapter = _build_adapter(tmp_path, monkeypatch, timeout=7)

        assert adapter.explore() is False
        assert adapter.succeeded is False
        assert any('timed out' in reason for reason in adapter.reasons), adapter.reasons
        assert not any('Arkane reported job failure' in reason for reason in adapter.reasons), \
            'A timeout must be its own diagnosis, not the generic job-failure one.'


class TestArkaneExplorerAdapterExclusivity:
    """Rule 0: refuse to launch into a non-exclusive run directory. No cleanup-and-proceed."""

    def test_construction_touches_the_filesystem_not_at_all(self, tmp_path, monkeypatch):
        """
        Constructing an adapter must write nothing -- not even the run directory itself.

        Rule 0 promises an EXCLUSIVE run directory. That promise is only true if nothing is
        written before the directory has been judged, so ``set_up()`` may not run from
        ``__init__``: a constructor that writes ``input.py`` has already violated the exclusivity
        it is about to check. The two are one defect, not two: once rule 0 requires an absent or
        empty directory, a writing constructor would make every second ``explore()`` refuse the
        mess the first construction left behind.
        """
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)
        assert not os.path.exists(output_directory), \
            'Constructing an adapter must not create the run directory, let alone write into it.'

    def test_refusal_is_side_effect_free(self, tmp_path, monkeypatch):
        """A refused run directory must be left exactly as it was found."""
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        os.makedirs(output_directory)
        _write(os.path.join(output_directory, 'output0.py'), 'stale')
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)

        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job',
                            _make_fake_run_arkane_job(success=True))
        assert adapter.explore() is False
        assert sorted(os.listdir(output_directory)) == ['output0.py'], \
            'A refused run must not leave input.py (or anything else) behind in the directory.'

    def test_refuses_any_unexpected_pre_existing_entry(self, tmp_path, monkeypatch):
        """
        Rule 0 is "absent or empty", not "free of three known names".

        Enumerating conflicts (``output*.py``, ``pdep/``, ``chem.inp``) leaves ``species/``,
        ``RMG_libraries/``, ``plots/``, ``supporting_information.csv`` and a stale ``input.py``
        all silently acceptable, so "exclusive run directory" was false as written.
        """
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        os.makedirs(os.path.join(output_directory, 'RMG_libraries'))
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)

        called = {'n': 0}

        def _poison(*args, **kwargs):
            called['n'] += 1
            return True

        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job', _poison)
        assert adapter.explore() is False
        assert called['n'] == 0, 'run_arkane_job must never be invoked once rule 0 refuses the directory.'
        assert any('RMG_libraries' in reason for reason in adapter.reasons)

    def test_a_pre_existing_empty_directory_is_refused(self, tmp_path, monkeypatch):
        """
        DELIBERATE INVERSION of a previously passing test, which asserted that an existing empty
        run directory "is the normal T3 case and must be accepted".

        It was not the normal case -- nothing in T3 pre-creates this directory -- and accepting it is
        incompatible with claiming the directory atomically, which is the whole point of rule 0. An
        empty directory is indistinguishable from one a competing explorer claimed microseconds ago
        and has not written into yet, so "empty" cannot mean "free"; the only way to know a directory
        is ours is to have been the one that created it. Accepting a pre-created directory would make
        the claim a check again, restoring exactly the check-then-write window this closes.
        """
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        os.makedirs(output_directory)
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)
        called = {'n': 0}
        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job',
                            lambda *args, **kwargs: called.__setitem__('n', called['n'] + 1))
        assert adapter.explore() is False
        assert called['n'] == 0, 'Arkane must never be launched into a directory this run did not create.'
        assert any('empty' in reason and 'unowned' in reason for reason in adapter.reasons), adapter.reasons

    def test_a_failed_claim_leaves_the_claimed_flag_false(self, tmp_path, monkeypatch):
        """
        The claimed flag exists to let ``set_up()`` tell a genuine rule-0 claim apart from a
        standalone call, so it must only flip on the success path of ``_claim_run_directory()``. If
        it flipped unconditionally (e.g. at the top of the method, before the ``mkdir`` even runs),
        a refused claim would still leave ``set_up()`` willing to write into a directory this run
        never actually claimed -- silently reopening the exclusivity hole rule 0 exists to close.
        """
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        os.makedirs(output_directory)
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)
        reasons = adapter._claim_run_directory()
        assert reasons, 'A pre-existing directory must be refused.'
        assert adapter._run_directory_claimed is False, \
            'A failed claim must not leave the adapter believing it owns the directory.'

    def test_an_absent_directory_is_created_and_claimed(self, tmp_path, monkeypatch):
        """The normal case: rule 0 creates the directory it claims, and the run proceeds."""
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
            ))
        assert adapter.explore() is True, adapter.reasons
        assert os.path.isdir(output_directory)

    def test_a_missing_parent_directory_is_refused_not_raised(self, tmp_path, monkeypatch):
        """
        The claim no longer creates parent directories, only the leaf: only a single ``mkdir``
        syscall for the leaf component is atomic. If this claim also walked and created
        intermediate path components (as ``os.makedirs`` would), those intermediate
        creates/traversals would not be atomic, and an intermediate symlink could redirect where
        the "claimed" directory actually lands. So a missing parent is now the caller's problem,
        and the claim reports it as an ordinary refusal rather than raising.
        """
        output_directory = os.path.join(str(tmp_path), 'nested', 'explorer_run')
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)
        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job',
                            _make_fake_run_arkane_job(success=True))
        assert adapter.explore() is False
        assert any('could not be created' in reason for reason in adapter.reasons), adapter.reasons

    def test_set_up_called_directly_refuses_and_writes_nothing(self, tmp_path, monkeypatch):
        """
        ``set_up()`` must stay public: it is ``@abstractmethod`` on ``PESExplorerAdapter``
        (t3/pdep/explorer/adapter.py:99), so it cannot be made private to stop a caller from
        invoking it standalone. Instead it refuses via a runtime check -- only ``explore()``'s
        call to ``_claim_run_directory()`` (rule 0) can make ``set_up()`` proceed, and a direct
        call without that claim must raise before writing anything.
        """
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)
        with pytest.raises(RuntimeError, match='not been claimed'):
            adapter.set_up()
        assert not os.path.exists(output_directory), \
            'set_up() must write nothing when the run directory has not been claimed.'

    def test_a_second_adapter_cannot_claim_the_same_directory(self, tmp_path, monkeypatch):
        """
        The race the old check-then-write rule 0 could not stop, reduced to its observable core: two
        adapters pointed at one directory, and the second must refuse rather than interleave its
        artifacts with the first's. Under the old rule both saw an absent directory and both
        proceeded; now the first one's claim is what the second one trips over.
        """
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        first = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)
        second = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
            ))
        assert first.explore() is True, first.reasons
        assert second.explore() is False
        assert any('not this run' in reason for reason in second.reasons), second.reasons

    def test_an_uncreatable_run_directory_is_refused_rather_than_raising(self, tmp_path, monkeypatch):
        """
        A claim can fail for reasons other than 'already taken' -- here the parent path is a regular
        file, so ``mkdir`` raises ``NotADirectoryError``. Rule 0 must report that as a refusal like
        any other, not let the OSError escape ``explore()``: a caller that gets an exception where
        every other rule-0 failure returns False has to special-case the filesystem.
        """
        blocker = os.path.join(str(tmp_path), 'a_file')
        _write(blocker, 'not a directory')
        output_directory = os.path.join(blocker, 'explorer_run')
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)
        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job',
                            _make_fake_run_arkane_job(success=True))
        assert adapter.explore() is False
        assert any('could not be created' in reason for reason in adapter.reasons), adapter.reasons

    def test_a_symlinked_run_directory_is_refused_without_an_islink_check_deciding_it(
            self, tmp_path, monkeypatch):
        """
        A symlink at the run-directory path is refused because ``mkdir`` fails on ANY existing path,
        not because a separate ``islink`` test rejected it. That distinction is the point: an
        ``islink``-based guard is a check, with the same window as the old rule, and it has to
        enumerate what to reject. ``mkdir`` refuses a symlinked path as a side effect of being
        atomic, so writes can never follow the link somewhere unintended.
        """
        elsewhere = os.path.join(str(tmp_path), 'elsewhere')
        os.makedirs(elsewhere)
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        os.symlink(elsewhere, output_directory)
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)
        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job',
                            _make_fake_run_arkane_job(success=True))
        assert adapter.explore() is False
        assert any('symlink' in reason for reason in adapter.reasons), adapter.reasons
        assert not os.listdir(elsewhere), 'A refused claim must not write through the symlink.'

    def test_refuses_stale_output0_py(self, tmp_path, monkeypatch):
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        os.makedirs(output_directory)
        _write(os.path.join(output_directory, 'output0.py'), 'stale')
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)

        called = {'n': 0}

        def _poison(*args, **kwargs):
            called['n'] += 1
            return True

        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job', _poison)
        assert adapter.explore() is False
        assert called['n'] == 0, 'run_arkane_job must never be invoked once rule 0 refuses the directory.'
        assert any('output0.py' in reason for reason in adapter.reasons)

    def test_refuses_preexisting_pdep_dir(self, tmp_path, monkeypatch):
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        os.makedirs(os.path.join(output_directory, 'pdep'))
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)

        called = {'n': 0}

        def _poison(*args, **kwargs):
            called['n'] += 1
            return True

        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job', _poison)
        assert adapter.explore() is False
        assert called['n'] == 0, 'run_arkane_job must never be invoked once rule 0 refuses the directory.'
        assert any('pdep' in reason for reason in adapter.reasons)

    def test_refuses_preexisting_chem_inp(self, tmp_path, monkeypatch):
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        os.makedirs(output_directory)
        _write(os.path.join(output_directory, 'chem.inp'), 'stale')
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)

        called = {'n': 0}

        def _poison(*args, **kwargs):
            called['n'] += 1
            return True

        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job', _poison)
        assert adapter.explore() is False
        assert called['n'] == 0, 'run_arkane_job must never be invoked once rule 0 refuses the directory.'
        assert any('chem.inp' in reason for reason in adapter.reasons)


class TestArkaneExplorerAdapterResolution:
    """Rules 2-4: index-set resolution, exact equality, and reduced-vs-full selection."""

    def test_success_single_network_full(self, tmp_path, monkeypatch):
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
            ))
        assert adapter.explore() is True
        assert adapter.reasons == tuple()
        networks = adapter.get_networks()
        assert len(networks) == 1
        assert networks[0].endswith('network0_full.py')

    def test_missing_member_fails(self, tmp_path, monkeypatch):
        """A resolved output0.py with no matching pdep/final/network0_full.py must fail."""
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={},
            ))
        assert adapter.explore() is False
        assert any('final' in reason.lower() for reason in adapter.reasons)

    def test_stale_extra_output_fails(self, tmp_path, monkeypatch):
        """An extra, unaccounted-for output file must fail exact index-set equality."""
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY, 'output1.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
            ))
        assert adapter.explore() is False
        assert any('output1.py' in reason for reason in adapter.reasons)


class TestArkaneExplorerAdapterSetUpFailure:
    """``set_up()`` now runs inside ``explore()``, so its refusals surface there instead."""

    def test_a_refused_input_file_is_raised_but_the_state_is_recorded_first(self, tmp_path, monkeypatch):
        """
        A refused input file re-raises -- and must not leave the adapter looking like it never ran.

        Swallowing this into a plain ``return False`` would let a misconfigured caller read an
        input error as "Arkane explored and found nothing". But leaving ``succeeded`` at its None
        "never ran" sentinel is just as misleading in the other direction: explore() did run, and
        it did fail.
        """
        adapter = _build_adapter(tmp_path, monkeypatch, bath_gas={})
        called = {'n': 0}
        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job',
                            lambda *a, **k: called.__setitem__('n', called['n'] + 1))

        with pytest.raises(ValueError, match='bath gas'):
            adapter.explore()

        assert called['n'] == 0
        assert adapter.succeeded is False, \
            'succeeded must not stay at the None "explore() was never called" sentinel.'
        assert adapter.reasons != tuple()
        with pytest.raises(RuntimeError):
            adapter.get_networks()

    def test_an_expected_source_hash_that_does_not_match_network_path_fails_rather_than_writing(
            self, tmp_path, monkeypatch):
        """
        The adapter must forward ``expected_source_hash`` all the way to the writer's single read.

        This is the end-to-end chain for the TOCTOU fix: ``t3.pdep.api.explore_pdep_network``
        content-hashes the network once and hands that hash, plus only the pathname, down through
        ``explorer_factory`` and ``ArkaneExplorerAdapter.__init__`` to
        ``write_arkane_explorer_input_file``, which reopens the pathname. If a wrong hash reaches
        the adapter (standing in here for a network file that changed after the API's own check),
        ``set_up()`` (invoked via ``explore()``, per this class's convention) must raise rather than
        silently write an explorer input from the mismatched bytes, and no Arkane subprocess may
        run.
        """
        adapter = _build_adapter(tmp_path, monkeypatch, expected_source_hash='sha256:' + '0' * 64)
        called = {'n': 0}
        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job',
                            lambda *a, **k: called.__setitem__('n', called['n'] + 1))

        with pytest.raises(ValueError, match='changed after it was validated'):
            adapter.explore()

        assert called['n'] == 0
        assert adapter.succeeded is False
        assert adapter.reasons != tuple()
        assert not os.path.isfile(os.path.join(adapter.output_directory, 'input.py')), \
            'A hash mismatch must be caught before anything is written into the run directory.'

    @pytest.mark.parametrize('error', [TypeError('unsupported operand'),
                                       KeyError('method'),
                                       AttributeError('lower'),
                                       RuntimeError('unclaimed')])
    def test_any_set_up_failure_records_the_state_not_just_value_and_os_errors(
            self, tmp_path, monkeypatch, error):
        """
        The recorded state must not depend on WHICH exception type ``set_up()`` raised.

        The narrow ``except (ValueError, OSError)`` encoded an assumption about the failure
        taxonomy that the validators do not actually guarantee: a non-numeric tolerance reaches
        ``math.isfinite`` as a ``TypeError``, a bad ``method`` key a ``KeyError``, a non-string
        where ``.lower()`` is called an ``AttributeError``. Every one of those escaped with
        ``succeeded`` still at None, so an adapter whose ``explore()`` demonstrably ran and failed
        reported itself as "explore() was never called" -- a lie in exactly the direction that
        hides a real configuration bug. The re-raise is the contract; the state record must be
        unconditional.
        """
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr('t3.pdep.explorer.arkane.write_arkane_explorer_input_file',
                            lambda *a, **k: (_ for _ in ()).throw(error))
        called = {'n': 0}
        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job',
                            lambda *a, **k: called.__setitem__('n', called['n'] + 1))

        with pytest.raises(type(error)):
            adapter.explore()

        assert called['n'] == 0
        assert adapter.succeeded is False, \
            f'a {type(error).__name__} from set_up() left the None "never ran" sentinel.'
        assert adapter.reasons != tuple()
        assert any(type(error).__name__ in reason for reason in adapter.reasons), \
            f'the reason must name the failure type so the caller can tell these apart: {adapter.reasons}'
        with pytest.raises(RuntimeError):
            adapter.get_networks()

    @pytest.mark.parametrize('target', ['check_arkane_me_success', 'parse_pdep_network_file'])
    def test_a_crash_while_ANALYSING_the_run_is_also_recorded_not_just_one_while_writing(
            self, tmp_path, monkeypatch, target):
        """
        The post-run analysis is the MORE exposed half of this defect, and it was left uncovered.

        Widening only ``set_up()``'s catch fixes the half of ``explore()`` that reads T3's own
        arguments and leaves untouched the half that reads Arkane's *output files* -- a far larger
        raise surface (``open()`` on a file that vanished or is unreadable, a UnicodeDecodeError on
        binary content, a RecursionError from a pathologically nested payload, or a plain bug in a
        payload checker). Every one of those propagated out of ``explore()`` with ``succeeded``
        still None, i.e. "explore() was never called" about a run that had already spawned Arkane
        and written artifacts to disk. That is the same lie as the ``set_up()`` one, told about a
        more expensive event.
        """
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
            ))
        monkeypatch.setattr(f't3.pdep.explorer.arkane.{target}',
                            lambda *a, **k: (_ for _ in ()).throw(RecursionError('too deep')))

        with pytest.raises(RecursionError):
            adapter.explore()

        assert adapter.succeeded is False, \
            f'a RecursionError from {target}() left the None "explore() was never called" sentinel.'
        assert any('RecursionError' in reason for reason in adapter.reasons), adapter.reasons
        with pytest.raises(RuntimeError):
            adapter.get_networks()

    def test_an_analysis_crash_does_not_publish_artifacts(self, tmp_path, monkeypatch):
        """
        A run whose analysis crashed must expose no artifacts -- not even ones already resolved.

        ``succeeded = False`` is only half the record. If the manifest or the resolved paths were
        published anyway, a caller that (wrongly but plausibly) reads ``final_network_paths``
        without first consulting ``succeeded`` would build QM on a network no check ever cleared.
        """
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
            ))
        monkeypatch.setattr('t3.pdep.explorer.arkane.check_arkane_me_success',
                            lambda *a, **k: (_ for _ in ()).throw(TypeError('bad payload')))

        with pytest.raises(TypeError):
            adapter.explore()

        assert adapter.final_network_paths == tuple()
        assert adapter.output_paths == tuple()
        assert adapter.manifest == dict()

    def test_a_crash_AFTER_the_verdict_unpublishes_the_artifacts_and_the_verdict(
            self, tmp_path, monkeypatch):
        """
        The one crash site that lands with ``succeeded`` already True and the paths already published.

        Manifest recording (rule 6) is the LAST step of the analysis and runs after
        ``succeeded = not reasons`` and after the resolved paths are assigned. So a raise from
        ``_build_manifest`` is the only case where the handler's reset is doing real work rather than
        clearing already-empty fields -- every earlier crash site is upstream of publication. Without
        this test the reset lines are decorative: deleting all four leaves the suite green.
        """
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
            ))
        monkeypatch.setattr(type(adapter), '_build_manifest',
                            lambda *a, **k: (_ for _ in ()).throw(KeyError('level_of_theory')))

        with pytest.raises(KeyError):
            adapter.explore()

        assert adapter.succeeded is False, \
            'the verdict was already True when the manifest raised; it must be withdrawn.'
        assert adapter.final_network_paths == tuple(), \
            'a published network path survived a crashed analysis.'
        assert adapter.output_paths == tuple()
        assert adapter.manifest == dict()
        with pytest.raises(RuntimeError):
            adapter.get_networks()


class TestArkaneExplorerAdapterFinalNetworkPayload:
    """
    The resolved final network file must carry a real network, not merely a matching NAME.

    This is the fail-open the rest of the criterion was built to prevent and then walked straight
    past: ``_resolve_artifacts`` matched ``pdep/final/network<i>_(full|reduced).py`` by filename
    only, and the payload check (``check_arkane_me_success``) ran exclusively on ``output.py``. So
    a zero-byte ``network0_full.py`` beside a perfectly valid ``output.py`` was a SUCCESS, and
    ``get_networks()`` then handed the caller an empty file to build QM on. Design section 5 says
    completeness is a property of the artifact SET; checking one member of the set and taking the
    other on trust is not that.
    """

    @staticmethod
    def _explore_with_final(tmp_path, monkeypatch, final_content, output_content=VALID_OUTPUT_PY):
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': output_content},
                final_files={'network0_full.py': final_content},
            ))
        return adapter

    def test_a_zero_byte_final_network_is_refused(self, tmp_path, monkeypatch):
        """The canonical fail-open: a valid output.py must not vouch for an empty network file."""
        adapter = self._explore_with_final(tmp_path, monkeypatch, final_content='')
        assert adapter.explore() is False
        assert any('network0_full.py' in reason for reason in adapter.reasons), adapter.reasons

    def test_a_whitespace_only_final_network_is_refused(self, tmp_path, monkeypatch):
        """Non-zero size is not evidence of content -- a size check alone would pass this."""
        adapter = self._explore_with_final(tmp_path, monkeypatch, final_content='\n\n   \n')
        assert adapter.explore() is False

    def test_an_output_file_masquerading_as_a_final_network_is_refused(self, tmp_path, monkeypatch):
        """
        An ME ``output.py`` is not a network file, even though both are valid Arkane-DSL Python.

        The earlier tests wrote ``VALID_OUTPUT_PY`` as ``network0_full.py`` and asserted success,
        which blessed exactly this confusion and masked the fail-open above.
        """
        adapter = self._explore_with_final(tmp_path, monkeypatch, final_content=VALID_OUTPUT_PY)
        assert adapter.explore() is False
        assert any('network0_full.py' in reason for reason in adapter.reasons), adapter.reasons

    def test_a_truncated_final_network_is_refused(self, tmp_path, monkeypatch):
        """A file cut off after the species blocks parses fine but declares no reaction."""
        truncated = VALID_NETWORK_PY.split('transitionState(')[0]
        adapter = self._explore_with_final(tmp_path, monkeypatch, final_content=truncated)
        assert adapter.explore() is False

    def test_a_final_network_that_is_not_python_is_refused(self, tmp_path, monkeypatch):
        adapter = self._explore_with_final(tmp_path, monkeypatch,
                                           final_content='this is not python (((')
        assert adapter.explore() is False

    def test_a_well_formed_final_network_is_accepted(self, tmp_path, monkeypatch):
        adapter = self._explore_with_final(tmp_path, monkeypatch, final_content=VALID_NETWORK_PY)
        assert adapter.explore() is True, adapter.reasons
        assert adapter.get_networks()[0].endswith('network0_full.py')

    def test_output_and_final_network_must_describe_the_same_network(self, tmp_path, monkeypatch):
        """
        Codex's P2, which is this same gap seen from the other side.

        Nothing checked that ``output.py`` and ``network0_full.py`` describe the SAME network, so a
        valid output for network A beside a valid network file for network B was a success. Every
        species named by a net reaction in output.py must be declared by the network file; a
        reaction over species the network never heard of means the two artifacts did not come from
        one run.
        """
        foreign_output = VALID_OUTPUT_PY.replace("'CH2O'", "'not_in_this_network'")
        adapter = self._explore_with_final(tmp_path, monkeypatch,
                                           final_content=VALID_NETWORK_PY,
                                           output_content=foreign_output)
        assert adapter.explore() is False
        assert any('not_in_this_network' in reason for reason in adapter.reasons), adapter.reasons

    def test_a_truncated_output_with_fewer_net_reactions_than_the_network_implies_is_refused(
            self, tmp_path, monkeypatch):
        """
        The gap this task exists to close: containment alone is silent about count.

        ``MULTI_OUTPUT_TRUNCATED_PY`` carries only five of the six net reactions
        ``MULTI_NETWORK_PY`` implies, but every species it names is still legitimately declared by
        the network (containment is clean) and every isomer/channel is still touched by at least
        one entry -- only the exact-count check can catch this.
        """
        adapter = self._explore_with_final(tmp_path, monkeypatch,
                                           final_content=MULTI_NETWORK_PY,
                                           output_content=MULTI_OUTPUT_TRUNCATED_PY)
        assert adapter.explore() is False
        assert any('5' in reason and '6' in reason for reason in adapter.reasons), adapter.reasons

    def test_an_output_with_more_net_reactions_than_the_network_implies_is_refused(
            self, tmp_path, monkeypatch):
        """A surplus is just as much evidence of a mismatched pairing as a shortfall."""
        adapter = self._explore_with_final(tmp_path, monkeypatch,
                                           final_content=MULTI_NETWORK_PY,
                                           output_content=MULTI_OUTPUT_SURPLUS_PY)
        assert adapter.explore() is False
        assert any('7' in reason and '6' in reason for reason in adapter.reasons), adapter.reasons

    def test_an_output_with_exactly_the_expected_net_reaction_count_is_accepted(
            self, tmp_path, monkeypatch):
        """The over-refusal guard: a genuinely complete solve over a richer network must pass."""
        adapter = self._explore_with_final(tmp_path, monkeypatch,
                                           final_content=MULTI_NETWORK_PY,
                                           output_content=MULTI_OUTPUT_EXACT_PY)
        assert adapter.explore() is True, adapter.reasons

    def test_an_untouched_isomer_is_refused_even_when_the_total_count_matches(
            self, tmp_path, monkeypatch):
        """
        Count alone is foolable: ``MULTI_OUTPUT_CHANNEL_GAP_PY`` has exactly six entries, but none
        of them touch isomer 'B' on either side (three of the genuine six were each duplicated to
        pad the count back up instead of solving 'B'). Only channel coverage catches this.
        """
        adapter = self._explore_with_final(tmp_path, monkeypatch,
                                           final_content=MULTI_NETWORK_PY,
                                           output_content=MULTI_OUTPUT_CHANNEL_GAP_PY)
        assert adapter.explore() is False
        assert any("'B'" in reason or '"B"' in reason for reason in adapter.reasons), adapter.reasons

    def test_a_truncation_that_still_covers_every_channel_is_caught_by_count_alone(
            self, tmp_path, monkeypatch):
        """
        Isolate the count check from the channel-coverage check: ``MULTI_OUTPUT_TRUNCATED_DROP_FIRST_PY``
        drops one net reaction (index 0, ('B',) -> ('A',)) but every declared channel -- 'A', 'B',
        ('R1', 'R2') and ('P1', 'P2') -- is still touched by one of the remaining five entries, so
        channel coverage alone would accept it. Only the exact-count check (``!=``, not a
        surplus-only ``>``) can catch a five-of-six shortfall shaped this way.
        """
        adapter = self._explore_with_final(tmp_path, monkeypatch,
                                           final_content=MULTI_NETWORK_PY,
                                           output_content=MULTI_OUTPUT_TRUNCATED_DROP_FIRST_PY)
        assert adapter.explore() is False
        assert any('5' in reason and '6' in reason for reason in adapter.reasons), adapter.reasons

    def test_an_uncovered_product_channel_is_refused_even_when_the_total_count_matches(
            self, tmp_path, monkeypatch):
        """
        Isolate the product-channel pairs from the source-to-source ones:
        ``MULTI_OUTPUT_MISSING_PRODUCT_CHANNEL_PY`` has exactly six entries and touches 'A', 'B' and
        ('R1', 'R2'), but never touches the declared product channel ('P1', 'P2') on either side.
        This is only caught if the expected pair set is derived from ``network.product_channels`` as
        destinations too, not just from the isomers and reactant channels.
        """
        adapter = self._explore_with_final(tmp_path, monkeypatch,
                                           final_content=MULTI_NETWORK_PY,
                                           output_content=MULTI_OUTPUT_MISSING_PRODUCT_CHANNEL_PY)
        assert adapter.explore() is False
        assert any("('P1', 'P2')" in reason for reason in adapter.reasons), adapter.reasons

    def test_a_forged_output_with_the_right_count_and_full_coverage_is_still_refused(
            self, tmp_path, monkeypatch):
        """
        Codex's round-29 P1 A: count equality plus 'every channel is touched somewhere' is not
        topology equality. ``MULTI_OUTPUT_SWAPPED_PAIR_PY`` omits one expected channel pair and
        duplicates another, keeping the total at six AND leaving every declared channel touched, so
        both of the older checks accept it. Only comparing the expected channel-PAIR set refuses it.
        """
        adapter = self._explore_with_final(tmp_path, monkeypatch,
                                           final_content=MULTI_NETWORK_PY,
                                           output_content=MULTI_OUTPUT_SWAPPED_PAIR_PY)
        assert adapter.explore() is False
        reason = ' '.join(adapter.reasons)
        # The omitted pair is named as missing AND the repeated one as surplus: both halves of the
        # multiset comparison have to fire, or a forgery that only duplicates would read as clean.
        assert "never connects the channel pair(s) [((\'P1\', \'P2\'), (\'R1\', \'R2\'))]" in reason, reason
        assert "which the network file does not imply" in reason, reason
        assert reason.count("(('A',), ('P1', 'P2'))") == 1, reason

    def test_an_output_naming_a_channel_pair_the_network_never_declares_is_refused(
            self, tmp_path, monkeypatch):
        """
        The surplus half of the comparison, isolated from mere repetition: every species named is
        declared (so the containment check passes) and the count is exactly six, but one entry
        connects ('A',) to ('R1',) -- and ('R1',) alone is not a declared channel of this network,
        only ('R1', 'R2') is. A guard that asked only 'are all the expected pairs present, and are
        all species declared' would accept it.
        """
        forged = ''.join(
            _pdepreaction(r, p) for r, p in (_MULTI_NETWORK_NET_REACTIONS[0],
                                             _MULTI_NETWORK_NET_REACTIONS[1],
                                             _MULTI_NETWORK_NET_REACTIONS[2],
                                             _MULTI_NETWORK_NET_REACTIONS[3],
                                             _MULTI_NETWORK_NET_REACTIONS[4],
                                             (('A',), ('R1',))))
        adapter = self._explore_with_final(tmp_path, monkeypatch,
                                          final_content=MULTI_NETWORK_PY,
                                          output_content=forged)
        assert adapter.explore() is False
        reason = ' '.join(adapter.reasons)
        assert "(('A',), ('R1',))" in reason, reason
        assert "which the network file does not imply" in reason, reason

    def test_an_entry_written_in_the_reverse_direction_is_still_accepted(
            self, tmp_path, monkeypatch):
        """
        The pair comparison is deliberately direction-insensitive. Arkane's own duplicate
        suppression compares by label in either direction, and which side of a net reaction the
        solver reports as the reactant side is not part of the network's identity -- so flipping one
        entry must NOT be read as a forged topology. Pins that the new check did not tighten
        identity into a direction requirement the artifacts never promised.
        """
        adapter = self._explore_with_final(tmp_path, monkeypatch,
                                           final_content=MULTI_NETWORK_PY,
                                           output_content=MULTI_OUTPUT_REVERSED_ENTRY_PY)
        assert adapter.explore() is True, adapter.reasons

    def test_a_reversed_entry_reaches_get_k_tp_still_reversed(self, tmp_path, monkeypatch):
        """
        The direction the identity check tolerates is the direction a consumer receives.

        The test above pins that a reversed entry passes identity; this one pins where it then goes,
        which is the half that can produce a wrong NUMBER rather than a wrong verdict. ``get_k_tp()``
        does not canonicalize -- a rate coefficient is genuinely directional, so normalizing it would
        destroy information nothing downstream could recover -- and it does not reverse the entry back
        to whatever the caller had in mind, because it does not know what the caller had in mind.

        So a consumer that reads ``reactants`` as "the direction I asked for" gets the reverse rate,
        silently. This test exists to make that contract falsifiable at the boundary rather than only
        asserted in a docstring: the very entry the identity check waves through is observed here,
        still reversed, in the value handed to the consumer. api.py must resolve direction itself.
        """
        adapter = self._explore_with_final(tmp_path, monkeypatch,
                                           final_content=MULTI_NETWORK_PY,
                                           output_content=MULTI_OUTPUT_REVERSED_ENTRY_PY)
        assert adapter.explore() is True, adapter.reasons

        reversed_entry = _MULTI_NETWORK_NET_REACTIONS[3]
        as_written = {(entry.reactants, entry.products) for entry in adapter.get_k_tp()}
        assert (reversed_entry[1], reversed_entry[0]) in as_written, (
            'the reversed entry should reach the consumer exactly as Arkane wrote it')
        assert (reversed_entry[0], reversed_entry[1]) not in as_written, (
            'get_k_tp() must not silently canonicalize direction back to the network file order')


class TestArkaneExplorerAdapterReductionPredicate:
    """
    Which final artifact a run must produce is decided ONLY by mirroring Arkane's own gate.

    Arkane runs its reduction loop iff ``self.energy_tol != np.inf or self.flux_tol != 0.0``
    (``arkane/explorer.py:303``); ``network<i>_full.py`` is written before that loop and
    ``network<i>_reduced.py`` after it. So the gate decides which file exists on disk, and T3 must
    reproduce it verbatim -- including the cases that look like nonsense, because it is Arkane's
    behavior and not Arkane's intent that determines the artifact.

    These call ``_resolve_artifacts()`` directly rather than going through ``explore()``. The
    writer refuses a non-finite tolerance outright, so such a value cannot reach a real run through
    ``set_up()``; the point of these tests is that the resolution rule is correct ON ITS OWN, not
    merely correct because a guard in a different module happens to run first. Two guards that only
    compose correctly by accident have already produced defects on this branch.
    """

    @staticmethod
    def _stage(tmp_path, monkeypatch, kinds, **kwargs):
        """Build an adapter and stage ``pdep/final/network0_<kind>.py`` for each kind."""
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory, **kwargs)
        _write(os.path.join(output_directory, 'output.py'), VALID_OUTPUT_PY)
        for kind in kinds:
            _write(os.path.join(output_directory, 'pdep', 'final', f'network0_{kind}.py'),
                   VALID_NETWORK_PY)
        return adapter

    @pytest.mark.parametrize('kwargs, expected', [
        # Arkane's own defaults for both tolerances: the gate is False, no reduction runs.
        (dict(energy_tol=None, flux_tol=None), 'network0_full.py'),
        # inf is Arkane's OWN default for energy_tol and means "no filter".
        (dict(energy_tol=float('inf'), flux_tol=None), 'network0_full.py'),
        (dict(energy_tol=None, flux_tol=0.0), 'network0_full.py'),
        # A real filter request: the reduction loop runs and writes the reduced network.
        (dict(energy_tol=5.0, flux_tol=None), 'network0_reduced.py'),
        (dict(energy_tol=None, flux_tol=0.1), 'network0_reduced.py'),
        # The defect: Arkane gates on `flux_tol != 0.0`, NOT `flux_tol > 0`. A finite negative
        # flux_tol makes Arkane reduce, so the reduced network is the one that matches output.py.
        # Requiring the pre-reduction full network here would accept a network that does not
        # correspond to the kinetics Arkane actually wrote.
        (dict(energy_tol=None, flux_tol=-1.0), 'network0_reduced.py'),
    ])
    def test_the_required_artifact_mirrors_arkanes_gate(self, tmp_path, monkeypatch, kwargs, expected):
        adapter = self._stage(tmp_path, monkeypatch, kinds=('full', 'reduced'), **kwargs)
        reasons, final_network_paths, _ = adapter._resolve_artifacts()
        assert reasons == tuple(), reasons
        assert os.path.basename(final_network_paths[0]) == expected

    def test_a_negative_flux_tol_refuses_a_run_that_only_produced_the_full_network(self, tmp_path, monkeypatch):
        """The mirrored gate must also REFUSE, not just select -- otherwise it is decorative."""
        adapter = self._stage(tmp_path, monkeypatch, kinds=('full',), energy_tol=None, flux_tol=-1.0)
        reasons, final_network_paths, _ = adapter._resolve_artifacts()
        assert reasons != tuple()
        assert final_network_paths == tuple()
        assert any('network0_reduced.py' in reason for reason in reasons)

    def test_an_infinite_energy_tol_accepts_a_run_that_only_produced_the_full_network(self, tmp_path, monkeypatch):
        """
        The mirror image of a fail-open, and just as wrong: rejecting a perfectly good run.

        Treating a non-finite tolerance as a reduction REQUEST would demand a
        ``network0_reduced.py`` that Arkane will never produce.
        """
        adapter = self._stage(tmp_path, monkeypatch, kinds=('full',),
                              energy_tol=float('inf'), flux_tol=None)
        reasons, final_network_paths, _ = adapter._resolve_artifacts()
        assert reasons == tuple(), reasons
        assert os.path.basename(final_network_paths[0]) == 'network0_full.py'

    def test_reduction_filter_selects_reduced_not_full(self, tmp_path, monkeypatch):
        """When energy_tol is finite, completeness must be measured against network0_reduced.py."""
        adapter = _build_adapter(tmp_path, monkeypatch, energy_tol=1.0)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                # _full.py is always written first (before the reduction loop); only _reduced.py
                # should be required/selected here.
                final_files={'network0_full.py': VALID_NETWORK_PY,
                            'network0_reduced.py': VALID_NETWORK_PY},
            ))
        assert adapter.explore() is True
        networks = adapter.get_networks()
        assert networks[0].endswith('network0_reduced.py')

    def test_reduction_filter_without_reduced_file_fails(self, tmp_path, monkeypatch):
        """A reduction filter was requested but only network0_full.py exists: must fail."""
        adapter = _build_adapter(tmp_path, monkeypatch, flux_tol=0.5)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
            ))
        assert adapter.explore() is False
        assert any('reduced' in reason.lower() for reason in adapter.reasons)

    def test_multi_network_refused(self, tmp_path, monkeypatch):
        """Settled verdict: refuse a multi-network run in this increment, naming the seed."""
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY, 'output0.py': VALID_OUTPUT_PY,
                             'output1.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY,
                            'network1_full.py': VALID_NETWORK_PY},
            ))
        assert adapter.explore() is False
        assert any('methoxy' in reason for reason in adapter.reasons), \
            'The refusal reason must name the seed.'
        assert any('multi' in reason.lower() or '2' in reason for reason in adapter.reasons)


class TestArkaneExplorerAdapterFailureSignals:
    """Rule 1: four independent failure signals, any one of which fails the run."""

    def test_stdout_only_fatal_marker_fails(self, tmp_path, monkeypatch):
        """Exit 0, empty stderr, a valid payload -- but a fatal marker in stdout.log must still fail."""
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
                stdout_lines=['Error: Negative rate coefficient generated; rejecting result.'],
            ))
        assert adapter.explore() is False
        assert any('Error: ' in reason for reason in adapter.reasons)

    def test_arkane_log_fatal_marker_fails(self, tmp_path, monkeypatch):
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
                arkane_log_lines=['Critical: something fatal happened'],
            ))
        assert adapter.explore() is False
        assert any('Critical: ' in reason for reason in adapter.reasons)

    def test_nonzero_exit_fails(self, tmp_path, monkeypatch):
        """Even with an otherwise-complete, valid artifact set, a reported job failure alone
        must fail the run -- this isolates the nonzero-exit signal from the artifact-resolution
        checks, which would otherwise also fail on their own if no artifacts were written."""
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=False,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
            ))
        assert adapter.explore() is False
        assert any('job failure' in reason.lower() for reason in adapter.reasons)

    def test_fake_run_arkane_job_honours_required_artifact(self, tmp_path, monkeypatch):
        """
        The stub backing these tests must mirror the real ``run_arkane_job``
        (``t3/runners/rmg_runner.py``) exactly: reporting success without producing the artifact
        it was told is required is precisely the lie the real runner refuses to tell (it deletes
        any stale artifact up front and then requires Arkane itself to have (re)written it). A
        stub that reports ``success`` regardless of whether ``required_artifact`` exists would
        make this failure signal untestable and could paper over a genuine fail-open in the
        production adapter.

        This pins the job-failure reason distinctively: ``_resolve_artifacts`` would also refuse
        a missing ``output.py`` on its own (rule 2-3), so asserting ``explore() is False`` alone
        would not prove the stub's own gate did anything. Only the job-failure text below is
        exclusive to ``run_arkane_job``'s own report.
        """
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                final_files={'network0_full.py': VALID_NETWORK_PY},
                # Deliberately no output_files: the required artifact 'output.py' is never
                # written, exactly as if Arkane exited 0 but never produced it.
            ))
        assert adapter.explore() is False
        assert any("the required 'output.py' artifact was never created" in reason
                  for reason in adapter.reasons), adapter.reasons

    def test_real_stderr_content_fails(self, tmp_path, monkeypatch):
        """Failure signal 2 of 4: real (non-ignorable) stderr content alone must fail the run,
        even with an otherwise-complete, valid artifact set and a zero exit status."""
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
                stderr_lines=['Traceback (most recent call last):', 'ValueError: bad input'],
            ))
        assert adapter.explore() is False
        assert any('stderr' in reason.lower() for reason in adapter.reasons)

    def test_ignorable_stderr_content_does_not_fail(self, tmp_path, monkeypatch):
        """Stderr lines matching IGNORABLE_STDERR_PHRASES must not, by themselves, fail the run."""
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
                stderr_lines=['Open Babel Warning  in PerceiveBondOrders'],
            ))
        assert adapter.explore() is True

    def test_soft_failure_payload_fails(self, tmp_path, monkeypatch):
        """Rule 5: the payload check (via check_arkane_me_success) must catch an all-None result."""
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': SOFT_FAILURE_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
            ))
        assert adapter.explore() is False

    def test_a_stale_stderr_log_is_refused_rather_than_deleted(self, tmp_path, monkeypatch):
        """
        A stale stderr.log cannot poison the fatal-marker scan, because rule 0 refuses it first.

        ``ArkaneMESolverAdapter.solve()`` deletes a stale stderr.log before running, because it
        tolerates a re-used run directory. The explorer does not: rule 0 requires an absent or
        empty directory, so the deletion idiom is not merely unnecessary here, it is unreachable.
        This test pins WHICH guard does the work, so that a future weakening of rule 0 fails
        loudly instead of silently reintroducing a poisoned scan.
        """
        output_directory = os.path.join(str(tmp_path), 'explorer_run')
        os.makedirs(output_directory)
        _write(os.path.join(output_directory, 'stderr.log'), 'Some previous unrelated stale error\n')
        adapter = _build_adapter(tmp_path, monkeypatch, output_directory=output_directory)
        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job',
                            _make_fake_run_arkane_job(success=True))
        assert adapter.explore() is False
        assert any('stderr.log' in reason for reason in adapter.reasons)
        assert os.path.isfile(os.path.join(output_directory, 'stderr.log')), \
            'Rule 0 refuses; it must not delete the operator\'s file on the way out.'


class TestArkaneExplorerAdapterGetKTP:

    def test_get_networks_before_explore_raises(self, tmp_path, monkeypatch):
        adapter = _build_adapter(tmp_path, monkeypatch)
        with pytest.raises(RuntimeError):
            adapter.get_networks()

    def test_get_k_tp_before_explore_raises(self, tmp_path, monkeypatch):
        adapter = _build_adapter(tmp_path, monkeypatch)
        with pytest.raises(RuntimeError):
            adapter.get_k_tp()

    def test_get_k_tp_after_failed_explore_raises(self, tmp_path, monkeypatch):
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(success=False))
        assert adapter.explore() is False
        with pytest.raises(RuntimeError):
            adapter.get_k_tp()
        with pytest.raises(RuntimeError):
            adapter.get_networks()

    def test_get_k_tp_after_success_returns_reactions(self, tmp_path, monkeypatch):
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
            ))
        assert adapter.explore() is True
        reactions = adapter.get_k_tp()
        assert len(reactions) == 1
        assert reactions[0].reactants == ('CH2O', 'H')
        assert reactions[0].products == ('methoxy',)


class TestArkaneExplorerAdapterManifest:

    def test_manifest_recorded_on_success(self, tmp_path, monkeypatch):
        adapter = _build_adapter(tmp_path, monkeypatch)
        monkeypatch.setattr(
            't3.pdep.explorer.arkane.run_arkane_job',
            _make_fake_run_arkane_job(
                success=True,
                output_files={'output.py': VALID_OUTPUT_PY},
                final_files={'network0_full.py': VALID_NETWORK_PY},
            ))
        assert adapter.explore() is True
        assert adapter.manifest
        assert 'artifacts' in adapter.manifest
        paths_recorded = {a['path'] for a in adapter.manifest['artifacts']}
        assert any(p.endswith('output.py') for p in paths_recorded)
        assert any(p.endswith('network0_full.py') for p in paths_recorded)
        for artifact in adapter.manifest['artifacts']:
            assert artifact['sha256'].startswith('sha256:')
            assert artifact['size'] > 0


class TestTheClaimStageCannotEscapeTheStateContract:
    """
    The claim runs BEFORE both guards, so a crash there was still the "never ran" lie.

    ``explore()`` now records its verdict before re-raising from ``set_up()`` and from the post-run
    analysis, but ``_claim_run_directory()`` is called at ``explore()``'s top level, outside both. Two
    ways it raises:

    * The ``except FileExistsError`` handler inspects the offending path to say WHAT was in the way --
      ``os.path.islink``, ``os.path.isdir``, ``os.listdir``. On a directory that exists but cannot be
      read (mode 000, or one whose permissions change under us), ``os.listdir`` raises PermissionError
      from inside the handler. Its own comment says the diagnosis "decides nothing", so a failure to
      diagnose must degrade the MESSAGE, never the verdict -- the claim has already failed.
    * ``os.mkdir`` itself raises TypeError, not OSError, for a non-path ``output_directory`` (None is
      the realistic one, since api.py may pass it through unset).
    """

    def test_an_unreadable_existing_directory_is_refused_not_raised(self, tmp_path, monkeypatch):
        """A PermissionError while DIAGNOSING must not escape, and must not change the verdict."""
        adapter = _build_adapter(tmp_path, monkeypatch)
        os.makedirs(adapter.output_directory)
        real_listdir = os.listdir

        def _refusing_listdir(path, *args, **kwargs):
            if os.path.realpath(path) == os.path.realpath(adapter.output_directory):
                raise PermissionError(13, 'Permission denied')
            return real_listdir(path, *args, **kwargs)

        monkeypatch.setattr(os, 'listdir', _refusing_listdir)
        called = {'n': 0}
        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job',
                            lambda *a, **k: called.__setitem__('n', called['n'] + 1))

        assert adapter.explore() is False
        assert called['n'] == 0
        assert adapter.succeeded is False
        assert any(adapter.output_directory in reason for reason in adapter.reasons), adapter.reasons
        assert adapter._run_directory_claimed is False

    def test_an_unclaimable_directory_of_the_wrong_type_records_the_state(self, tmp_path, monkeypatch):
        """
        ``output_directory=None`` makes os.mkdir raise TypeError, which is not an OSError.

        Recorded and re-raised, matching how set_up() and the analysis handle the same class of
        failure: the caller gets the exception, and the adapter never claims it was not run.
        """
        adapter = _build_adapter(tmp_path, monkeypatch)
        adapter.output_directory = None
        called = {'n': 0}
        monkeypatch.setattr('t3.pdep.explorer.arkane.run_arkane_job',
                            lambda *a, **k: called.__setitem__('n', called['n'] + 1))

        with pytest.raises(TypeError):
            adapter.explore()

        assert called['n'] == 0
        assert adapter.succeeded is False, 'a TypeError from the claim left the None "never ran" sentinel.'
        assert any('TypeError' in reason for reason in adapter.reasons), adapter.reasons
        assert adapter._run_directory_claimed is False
        with pytest.raises(RuntimeError):
            adapter.get_networks()
