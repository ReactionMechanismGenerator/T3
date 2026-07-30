"""
t3 pdep explorer arkane module.

Concrete ``PESExplorerAdapter`` implementation that drives Arkane's PES explorer (the
``explorer(...)`` DSL block) from a seed of one or two source species and decides, per
``docs/t3-pdep-qm-explorer-design.md`` section 5, whether the exploration genuinely succeeded.

Arkane's explorer job is much easier to fool than a plain ME solve: it pre-creates ``output.py``
before any job runs, has a job-per-network naming convention that is only unambiguous once the
final network count is known, and (for a multi-network source) writes a REDUCED final network file
on top of the FULL one whenever a reduction filter was requested, with no marker distinguishing
"nothing was removed" from "the filter never ran". This module resolves each of those ambiguities
explicitly rather than accepting any single success-looking signal at face value.
"""

import importlib.metadata
import math
import os
import re
from collections import Counter
from typing import TYPE_CHECKING

from t3.pdep.cache import hash_file
from t3.pdep.explorer.adapter import PESExplorerAdapter
from t3.pdep.explorer.factory import register_explorer_adapter
from t3.pdep.explorer.input_file import write_arkane_explorer_input_file
from t3.pdep.me_success import IGNORABLE_STDERR_PHRASES, check_arkane_me_success
from t3.pdep.parser import canonical_channel_pair, parse_pdep_network_file
from t3.runners.rmg_runner import run_arkane_job

if TYPE_CHECKING:
    from t3.logger import Logger

# Arkane's own log formatter (see ``arkane/main.py::initialize_log``) maps CRITICAL/ERROR/WARNING
# to the literal line prefixes below (INFO/DEBUG carry no prefix at all), and attaches that same
# formatter to BOTH the stdout stream handler and the ``arkane.log`` file handler. A line starting
# with one of these prefixes, in either stream, is therefore a low-false-positive fatal marker.
_FATAL_LOG_MARKERS = ('Critical: ', 'Error: ')

# Any file directly under a candidate run directory matching this pattern is evidence of a prior
# (possibly stale) Arkane explorer run. Rule 0 refuses to launch into a directory containing any
# such file, or a 'pdep/' directory, or a 'chem.inp' -- with no cleanup-and-proceed, ever.
_OUTPUT_FILENAME_RE = re.compile(r'^output(\d*)\.py$')

# Matches a final network artifact under 'pdep/final/', capturing its network index and whether it
# is the FULL (always written) or REDUCED (only written when a reduction filter ran) variant.
_FINAL_NETWORK_FILENAME_RE = re.compile(r'^network(\d+)_(full|reduced)\.py$')


def _get_rmgpy_revision() -> str | None:
    """
    Best-effort lookup of the installed RMG-Py package version, never by importing it.

    Explorer output is a function of RMG database contents (thermo/kinetics libraries, rate
    rules); a run recorded with no RMG-Py provenance at all is not replayable later.

    Returns:
        Optional[str]: The installed 'rmgpy' distribution version, or None if it could not be
                      determined.
    """
    try:
        return importlib.metadata.version('rmgpy')
    except importlib.metadata.PackageNotFoundError:
        return None


class ArkaneExplorerAdapter(PESExplorerAdapter):
    """
    A PESExplorerAdapter that drives Arkane's PES explorer from a 1- or 2-species seed.

    Note:
        ``t3.pdep.explorer.factory.explorer_factory`` does not currently pass a ``network_path``
        or ``method`` through to the adapter it constructs, even though both are required to build
        an Arkane explorer input file (``write_arkane_explorer_input_file`` requires a source RMG
        P-dep network path and a method). This is a real, deliberate gap: wiring
        ``explorer_factory``/``t3/pdep/api.py`` to supply these is out of scope for this adapter
        and is left for a follow-up commit (the ``explore=True`` wiring). Until then, this adapter
        must be constructed directly with ``network_path``/``method`` rather than through the
        factory.

        This increment refuses (rather than resolves) a multi-network exploration: Arkane's
        explorer can, from a single seed, discover MORE than one distinct reaction network, in
        which case there is no single unambiguous "the" result to hand back from ``get_networks()``
        / ``get_k_tp()``. Detecting this case still goes through the same index-resolution logic
        used for the (currently supported) single-network case -- refusal rides atop resolution,
        it does not replace it.
    """

    def __init__(self,
                 seed_species,
                 output_directory: str,
                 network_path: str,
                 method: str,
                 bath_gas: dict = None,
                 explore_tol: float = None,
                 energy_tol: float = None,
                 flux_tol: float = None,
                 maximum_radical_electrons: int = None,
                 logger: 'Logger | None' = None,
                 transition_state_seeds: tuple = None,
                 database_kwargs: dict = None,
                 expected_source_hash: str = None,
                 ):
        """
        Args:
            seed_species (list | tuple): The source (seed) species labels to explore from (1 or 2
                                        labels; see ``PESExplorerAdapter.max_source_species``).
            output_directory (str): The path to the directory to run Arkane in. Refused if it
                                    already contains any run artifact of a prior exploration (rule
                                    0); created if it does not yet exist.
            network_path (str): The path to the source RMG P-dep network (or hybrid Arkane
                                network) file to build the explorer input from.
            method (str): 'CSE', 'MSC' or 'RS' (see ``t3.utils.writer.METHOD_MAP``).
            bath_gas (dict, optional): The bath gas composition, mapping species labels to mole
                                      fractions.
            explore_tol (float, optional): The energy tolerance for exploring new isomers/reactions.
            energy_tol (float, optional): The energy tolerance for including a well/transition
                                         state in the output network. A non-None value means a
                                         reduction filter was requested (rule 4).
            flux_tol (float, optional): The flux tolerance for including a well/transition state in
                                       the output network. A value greater than 0 means a reduction
                                       filter was requested (rule 4).
            maximum_radical_electrons (int, optional): The maximum number of radical electrons
                                                       allowed in an explored species.
            logger (Logger, optional): The current T3 Logger instance.
            transition_state_seeds (tuple, optional): Unused: Arkane's explorer never supports
                                                      transition-state seeds (see
                                                      ``PESExplorerAdapter.supports_transition_state_seeds``).
            database_kwargs (dict, optional): Keyword arguments describing the RMG database
                                             settings to use for the exploration.
            expected_source_hash (str, optional): The content hash (``t3.pdep.hashing`` format)
                                                  ``network_path``'s bytes must match at the moment
                                                  ``set_up()`` reads them, forwarded verbatim to
                                                  ``write_arkane_explorer_input_file``. See that
                                                  function's docstring for why this must be the
                                                  hash checked against the same read that consumes
                                                  the bytes, not an earlier, separate one.
        """
        super().__init__(seed_species=seed_species, transition_state_seeds=transition_state_seeds)
        self.output_directory = output_directory
        self.network_path = network_path
        self.method = method
        self.bath_gas = bath_gas or dict()
        self.explore_tol = explore_tol
        self.energy_tol = energy_tol
        self.flux_tol = flux_tol
        self.maximum_radical_electrons = maximum_radical_electrons
        self.logger = logger
        self.database_kwargs = database_kwargs
        self.expected_source_hash = expected_source_hash

        # Set True only by ``_claim_run_directory()`` on its success path (rule 0). ``set_up()``
        # refuses to run while this is False, so it cannot be used to bypass the claim.
        self._run_directory_claimed = False

        # ``None`` means "explore() has not been called yet"; this is distinct from a False
        # ``succeeded``, mirroring ``ArkaneMESolverAdapter``'s ``me_success_result`` sentinel so
        # that get_networks()/get_k_tp() can tell "never ran" apart from "ran and failed".
        self.succeeded = None
        self.reasons = tuple()
        self.final_network_paths = tuple()
        self.output_paths = tuple()
        self.manifest = dict()
        self._me_success_results = tuple()

    def set_up(self):
        """
        Write the Arkane explorer input file into the (already-claimed) run directory.

        Deliberately NOT called from ``__init__``. ``explore()`` calls it only after rule 0 has
        CLAIMED the run directory, because a constructor that writes ``input.py`` has already broken
        the exclusivity that rule 0 exists to guarantee -- and, since rule 0 claims the directory by
        creating it, a writing constructor would make the adapter refuse the very directory its own
        construction created. ``ArkaneMESolverAdapter`` documents the ``set_up()``-from-``__init__``
        idiom as a known sharp edge; the edge is load-bearing here, so it is not mirrored.

        The run directory is created by the claim (``_claim_run_directory()``), not by this method:
        this method refuses -- raises, before writing anything -- if that claim has not happened.
        ``set_up()`` must stay public (it is ``@abstractmethod`` on ``PESExplorerAdapter``), so a
        caller can still invoke it directly without going through ``explore()``; this runtime check
        is what stops that direct call from writing into a directory nobody claimed.

        Raises:
            RuntimeError: If the run directory has not been claimed yet.
        """
        if not self._run_directory_claimed:
            raise RuntimeError(
                f'set_up() cannot run: the run directory {self.output_directory!r} has not been '
                f'claimed. Only explore() claims the run directory (rule 0); set_up() cannot be '
                f'used standalone.')
        input_file_path = os.path.join(self.output_directory, 'input.py')
        write_arkane_explorer_input_file(source_path=self.network_path,
                                         dest_path=input_file_path,
                                         seed_species=self.seed_species,
                                         method=self.method,
                                         bath_gas=self.bath_gas,
                                         explore_tol=self.explore_tol,
                                         energy_tol=self.energy_tol,
                                         flux_tol=self.flux_tol,
                                         maximum_radical_electrons=self.maximum_radical_electrons,
                                         database_kwargs=self.database_kwargs,
                                         expected_source_hash=self.expected_source_hash,
                                         )

    def _claim_run_directory(self) -> tuple:
        """
        Rule 0: the single most important rule. CLAIM the run directory, atomically, or refuse.

        Every completeness and index-resolution check downstream reads this directory by pattern, so
        anything already in it is a potential false witness. The rule used to be "the directory must
        be absent or empty", which was the right intent expressed as a check -- and a check is the
        wrong shape for it. ``exists``/``listdir`` answered a question about the past; between that
        answer and the ``input.py`` write there was a window in which a second explorer could ask the
        same question and get the same answer, so two runs could both conclude they owned the
        directory and then interleave their artifacts. "Nothing was here before us" is not something
        a caller can observe; it is something it has to establish.

        So the directory is now CREATED here, and creating it IS the claim: ``os.mkdir`` is one
        ``mkdir`` syscall for exactly the leaf component, which either creates the directory or
        fails, atomically, with no window in between. It does NOT create parent directories -- the
        parent must already exist, and creating it is the caller's responsibility. Only the leaf
        mkdir is atomic; if this method walked and created intermediate path components too (as
        ``os.makedirs`` does), those intermediate creates/traversals would not be atomic, and an
        intermediate symlink could redirect where the "claimed" directory actually lands.

        Two consequences, both deliberate:

        - A pre-existing EMPTY directory is now REFUSED, where it used to be accepted. There is no
          way to keep both: accepting a directory somebody else created means the claim is a check
          again, and an empty directory is indistinguishable from one a competing explorer claimed
          microseconds ago and has not written into yet. Nothing in T3 pre-creates this directory
          (``t3/pdep/api.py`` does not), so accepting it bought a caller convenience nobody uses at
          the cost of the only property rule 0 exists to provide.
        - A symlink at the run-directory path is refused for free, without an ``islink`` check:
          ``mkdir`` fails with ``FileExistsError`` whether the path is a directory, a file, or a
          symlink to either, so a symlinked output directory can no longer receive writes. Path
          CONTAINMENT (is this directory inside the tree the caller intended?) is a different
          question and stays at the public API boundary, where the intended root is known.

        There is deliberately no stale-claim recovery and no liveness probe, unlike
        ``t3.pdep.capture._acquire_capture_lock``: that lock guards a shared resource whose owner
        may crash, and wedging it forever would be worse than a rare wrong steal. Here the claim IS
        the run directory of one specific run, an abandoned one still holds the artifacts that a
        forensic read may want, and stealing it means writing this run's output on top of another
        run's evidence. Fail closed; a human names a fresh directory.

        Returns:
            tuple: Human-readable reasons the directory could not be claimed; empty on success, in
                  which case the directory now exists and is owned by this adapter.
        """
        try:
            os.mkdir(self.output_directory)
        except FileExistsError:
            # Purely to say WHAT was in the way. The claim already failed; nothing below decides it.
            # Which is exactly why the diagnosis is wrapped: os.listdir/os.path.isdir raise
            # PermissionError on a directory that exists but cannot be read (mode 000, or permissions
            # changed under us), and letting that escape would turn a DECIDED refusal into an
            # exception out of explore() with ``succeeded`` still at its None "never ran" sentinel.
            # A guard whose own comment says it decides nothing must not be able to change the
            # outcome -- so a failure to diagnose degrades the MESSAGE and nothing else.
            try:
                if os.path.islink(self.output_directory) or not os.path.isdir(self.output_directory):
                    found = ("something already exists at that path and it is not a plain directory (a "
                             "file, or a symlink, which could point anywhere)")
                else:
                    existing = sorted(os.listdir(self.output_directory))
                    found = (f"it already exists and contains {existing}" if existing else
                             "it already exists and is empty -- which is not a free directory but an "
                             "unowned one: it may be a leftover, or another explorer may have claimed it "
                             "moments ago and not written into it yet")
            except OSError as e:
                found = (f"something already exists at that path and it could not even be inspected "
                         f"to say what ({type(e).__name__}: {e})")
            return (f"Refusing to launch an Arkane exploration into '{self.output_directory}': "
                    f"{found}. This adapter claims its run directory by creating it, so an existing "
                    f"path is by definition not this run's. There is no cleanup-and-proceed and no "
                    f"stale-claim recovery here -- anything already there could poison the "
                    f"index-resolution and completeness checks with evidence from another run, and "
                    f"overwriting it would destroy that run's evidence. Use a fresh "
                    f"output_directory.",)
        except OSError as e:
            return (f"Refusing to launch an Arkane exploration: the run directory "
                    f"'{self.output_directory}' could not be created: {e}.",)
        self._run_directory_claimed = True
        return tuple()

    def explore(self) -> bool:
        """
        Execute the Arkane exploration and decide whether it genuinely succeeded.

        Returns:
            bool: Whether the exploration succeeded.
        """
        self.succeeded = None
        self.reasons = tuple()
        self.final_network_paths = tuple()
        self.output_paths = tuple()
        self.manifest = dict()
        self._me_success_results = tuple()

        # Rule 0 runs BEFORE ``set_up()``, which is why ``set_up()`` is not called from
        # ``__init__``: the exclusivity guarantee is only real if nothing has been written yet. It
        # does not merely inspect the directory, it creates it -- see ``_claim_run_directory``.
        # Wrapped for the same reason set_up() and _run_and_evaluate() are. ``_claim_run_directory``
        # returns refusal reasons for the failures it anticipates, but ``os.mkdir`` raises TypeError --
        # not OSError -- for a non-path ``output_directory`` (None being the realistic one, since a
        # caller may pass it through unset), and that would leave ``succeeded`` at None about a run
        # whose explore() was called. Every one of explore()'s three stages now records before it
        # re-raises; a contract that holds in two of three places is not a contract.
        try:
            reasons = list(self._claim_run_directory())
        except Exception as e:
            self.reasons = (f'The Arkane exploration run directory could not be claimed '
                            f'({type(e).__name__}): {e}',)
            self.succeeded = False
            raise
        if reasons:
            self.reasons = tuple(reasons)
            self.succeeded = False
            return False
        try:
            self.set_up()
        except Exception as e:
            # Re-raised, not swallowed: a refused seed / bath gas / tolerance is an input or
            # configuration error, and turning it into a quiet ``return False`` would let a
            # misconfigured caller read it as "Arkane explored and found nothing". But the state is
            # recorded FIRST, because leaving ``succeeded`` at None would make get_networks() and
            # get_k_tp() report "explore() was never called" about an adapter whose explore() did
            # run and did fail.
            #
            # Deliberately ``Exception``, not ``(ValueError, OSError)``. Naming the types encoded an
            # assumption about the failure taxonomy that the validators do not guarantee: a
            # non-numeric tolerance reaches ``math.isfinite`` as a TypeError, an unexpected
            # ``database_kwargs`` key can surface as a KeyError, a non-string where ``.lower()`` is
            # called as an AttributeError. Each of those escaped with ``succeeded`` still None, so an
            # adapter that demonstrably ran and failed reported "explore() was never called" -- the
            # lie in the direction that hides a configuration bug. The type is folded into the reason
            # so widening the catch does not blur the diagnosis. BaseException is deliberately NOT
            # caught: a KeyboardInterrupt/SystemExit is not a determination about this network, and
            # the None sentinel is then the honest state.
            self.reasons = (f'The Arkane explorer input file could not be written '
                            f'({type(e).__name__}): {e}',)
            self.succeeded = False
            raise

        # Wrapped, for the same reason ``set_up()`` above is: an exception escaping the analysis
        # would leave ``succeeded`` at its None "explore() was never called" sentinel about a run
        # that had already spawned Arkane and written artifacts. This half of explore() has the
        # LARGER raise surface of the two, because it reads Arkane's output rather than T3's own
        # arguments -- an ``open()`` on a file that vanished mid-run, a UnicodeDecodeError on binary
        # content, a RecursionError from a pathologically nested payload, or a plain bug in a payload
        # checker. Re-raised rather than folded into ``return False``: an analysis that could not
        # complete has not proven anything, and a crash in our own checker must stay loud.
        try:
            self._run_and_evaluate()
        except Exception as e:
            self.reasons = (f'The Arkane exploration could not be evaluated '
                            f'({type(e).__name__}): {e}',)
            self.succeeded = False
            # Nothing resolved mid-analysis is published: a caller reading final_network_paths
            # without consulting succeeded must not receive a network no check ever cleared.
            self.final_network_paths = tuple()
            self.output_paths = tuple()
            self.manifest = dict()
            self._me_success_results = tuple()
            raise
        return self.succeeded

    def _run_and_evaluate(self) -> None:
        """
        Run Arkane and evaluate every success signal, recording the verdict on ``self``.

        Split out of ``explore()`` only so the whole post-run analysis sits under one guard; the
        caller owns the "record the state, then re-raise" contract. ``reasons`` starts empty here
        because ``explore()`` returns early when the rule-0 claim produced any.
        """
        reasons = list()
        input_file_path = os.path.join(self.output_directory, 'input.py')
        stdout_path = os.path.join(self.output_directory, 'stdout.log')
        stderr_path = os.path.join(self.output_directory, 'stderr.log')
        arkane_log_path = os.path.join(self.output_directory, 'arkane.log')
        # No stale-log deletion here, deliberately. ``ArkaneMESolverAdapter.solve()`` deletes a
        # stale stderr.log before running because it tolerates a re-used directory; rule 0 above
        # does not, so by this point the directory was absent or empty and a stale
        # stdout.log/stderr.log/arkane.log cannot exist. A deletion here would be dead code that
        # reads as a defense, and would quietly mask any future weakening of rule 0.

        # Called synchronously, exactly as ``ArkaneMESolverAdapter.solve()`` calls it. There is
        # deliberately NO timeout here, because T3 cannot honour one and a timeout that cannot be
        # honoured is worse than none. The earlier attempt wrapped this call in a
        # ThreadPoolExecutor and gave up on the future, which failed twice over:
        #   1. It did not even return early -- `finally: executor.shutdown(wait=True)` ran before
        #      the `return False` propagated and blocked until Arkane finished anyway, so the
        #      "timeout" bought nothing but a relabelled outcome.
        #   2. It could not have worked at any wait= setting, because nothing cancellable escapes:
        #      run_arkane_job (t3/runners/rmg_runner.py:543) calls ARC's run_arkane in-process,
        #      which reaches subprocess.run at ARC/arc/job/local.py:59-61 with no timeout= and no
        #      Popen/PID exposed to any caller. Python threads cannot be killed either.
        # Abandoning the future without killing the process is not a lesser version of a timeout,
        # it is a corruption hazard: Arkane would keep writing into a run directory T3 had already
        # declared failed, and rule 0 would then refuse that directory forever after. A real
        # timeout has to be built where the process is spawned -- either a `timeout=` on ARC's
        # execute_command, or a runner that spawns Arkane in its own process group and kills the
        # group -- and belongs in the runner layer, not smuggled in here.
        job_ran = run_arkane_job(input_file=input_file_path,
                                 output_directory=self.output_directory,
                                 logger=self.logger,
                                 required_artifact='output.py')

        # Failure signal 1 of 4: nonzero exit (surfaced here as run_arkane_job's own bool).
        if not job_ran:
            reasons.append('Arkane reported job failure (a non-zero exit status, or the required '
                           "'output.py' artifact was never created).")

        # Failure signal 2 of 4: real stderr content. Stderr alone is not sufficient to prove
        # success (a well-formed payload can still coexist with a real stderr problem elsewhere in
        # the run), but real (non-ignorable) stderr content alone IS sufficient to fail the run.
        stderr_text = None
        if os.path.isfile(stderr_path):
            with open(stderr_path, 'r') as f:
                stderr_text = f.read()
        if stderr_text:
            real_stderr_lines = [
                line.strip() for line in stderr_text.splitlines()
                if line.strip() and not any(phrase in line for phrase in IGNORABLE_STDERR_PHRASES)
            ]
            if real_stderr_lines:
                reasons.append('Arkane reported stderr output: ' + ' | '.join(real_stderr_lines))

        # Failure signal 3 of 4: a fatal marker in stdout OR Arkane's own log. This is the ONLY
        # signal that catches Arkane's soft-CSE-failure archetype at the process level: exit 0,
        # empty stderr, yet "Error: Negative rate coefficient generated; rejecting result." printed
        # to stdout only.
        fatal_lines = []
        for log_path in (stdout_path, arkane_log_path):
            if os.path.isfile(log_path):
                with open(log_path, 'r') as f:
                    for line in f:
                        if line.startswith(_FATAL_LOG_MARKERS):
                            fatal_lines.append(f'{os.path.basename(log_path)}: {line.strip()}')
        if fatal_lines:
            reasons.append('Arkane reported fatal log marker(s): ' + ' | '.join(fatal_lines))

        # Failure signal 4 of 4: a missing/invalid artifact, decided by index resolution (rules
        # 2-4) below.
        resolution_reasons, final_network_paths, output_paths = self._resolve_artifacts()
        reasons.extend(resolution_reasons)

        # Rule 5: per-file payload check, reusing check_arkane_me_success -- only meaningful once
        # the output set itself has been resolved unambiguously.
        me_success_results = []
        if output_paths and not resolution_reasons:
            for output_path in output_paths:
                result = check_arkane_me_success(output_path=output_path)
                me_success_results.append(result)
                if not result.succeeded:
                    reasons.extend(result.reasons)

        # Rule 5b: the final network artifact carries a real network, and the SAME network the
        # output describes. Rules 2-4 match 'network<i>_(full|reduced).py' by NAME only, and rule 5
        # checks payloads of output files alone -- so without this step a zero-byte
        # 'network0_full.py' beside a valid 'output.py' was a success, and get_networks() handed
        # the caller an empty file. Design section 5 makes completeness a property of the artifact
        # SET; checking one member and taking the other on trust is not that.
        if final_network_paths and not resolution_reasons:
            reasons.extend(self._check_final_network_payload(
                final_network_path=final_network_paths[0],
                me_success_results=tuple(me_success_results)))

        self.reasons = tuple(reasons)
        self.succeeded = not reasons
        if self.succeeded:
            self.final_network_paths = final_network_paths
            self.output_paths = output_paths
            self._me_success_results = tuple(me_success_results)
            # Rule 6: record a run manifest. Explorer output is a function of RMG database
            # contents, so a run recorded with no provenance at all is not replayable.
            self.manifest = self._build_manifest(input_file_path=input_file_path,
                                                 arkane_log_path=arkane_log_path)


    def _check_final_network_payload(self, final_network_path: str, me_success_results: tuple) -> tuple:
        """
        Rule 5b: prove the resolved final network file is a real network, and the right one.

        Four distinct questions, all previously unasked:

        1. Is it a network at all? ``t3.pdep.parser.parse_pdep_network_file`` is the existing
           parser for exactly this file shape and already fails closed on unparseable text and on
           a file declaring no ``reaction(...)`` at all, which covers the zero-byte, whitespace-
           only and truncated cases. On top of that this requires the ``network(...)`` block that
           ``PressureDependenceJob.save_input_file`` ALWAYS emits (arkane/pdep.py:741-755), so an
           ME ``output.py`` -- valid Arkane-DSL Python, but full of ``pdepreaction(...)`` and
           carrying no network block -- cannot pass as a network file.
        2. Is it the SAME network ``output.py`` describes? Nothing tied the two artifacts
           together, so a valid output for network A beside a valid network file for network B was
           a success. Every species named by a net reaction in the output must be declared by the
           network file: Arkane writes both from one in-memory network, so a net reaction over a
           species the network never declared means these two files did not come from one run.
           Containment is the sound test here, not equality -- a network legitimately declares
           species (bath gas, unreacted channels) that no surviving net reaction mentions.
        3. Does the output carry EVERY net reaction the network implies, and no more? Containment
           (question 2) is one-directional and count-blind: an ``output.py`` truncated after the
           first few ``pdepreaction(...)`` entries names only species the network legitimately
           declares, so it passes containment while silently handing a partial solve downstream.
           ``PDepNetwork.expected_net_reaction_count()`` is documented (and tested against a real
           Arkane output, ``tests/test_pdep/test_parser.py``) as EXACT -- it simulates Arkane's own
           duplicate-suppressing enumeration rather than a closed form -- so ANY mismatch, fewer or
           more, means this pairing is not a complete, single-run solve.
        4. Does the output touch every channel the network declares -- isomers, reactant channels
           AND product channels? Question 3 alone would still accept an output whose entries all
           land on the same subset of channels while an unrelated duplicate entry pads the count
           back up to the expected total. Every isomer and reactant channel is a source in Arkane's
           own net-reaction enumeration (``expected_net_reaction_count``'s ``sources``), so one that
           no parsed net reaction touches, on either side, was never solved even though the
           aggregate count matches.

           Product channels are required too, which is stronger than "every source" but still
           cannot over-refuse a complete solve. They appear only as destinations, yet the
           enumeration visits EVERY (source, product channel) pair, and the sole reason such a pair
           is not written live is that this exact pair was already written -- which happens only
           when the configuration is also a source, and leaves the channel present in the output
           either way. So a declared channel absent from every parsed net reaction means the solve
           is incomplete, never merely that Arkane deduplicated it away.

        Args:
            final_network_path (str): Path to the resolved 'network<i>_(full|reduced).py'.
            me_success_results (tuple): The per-output ``check_arkane_me_success`` results, whose
                                        parsed reactions supply the species labels to cross-check.

        Returns:
            tuple: Human-readable reasons the final network artifact is not trustworthy; empty if
                  it is.
        """
        file_name = os.path.basename(final_network_path)
        try:
            network = parse_pdep_network_file(path=final_network_path)
        except ValueError as e:
            return (f"The resolved final network file '{file_name}' is not a readable P-dep "
                    f"network: {e}",)
        if not network.isomers and not network.reactant_channels:
            return (f"The resolved final network file '{file_name}' declares no network(...) "
                    f"block with isomers or reactant channels, so it is not a complete final "
                    f"network artifact. Arkane always writes that block "
                    f"(arkane/pdep.py:741-755); its absence means the file is truncated, or is "
                    f"not a network file at all.",)

        declared = set(network.species_labels)
        undeclared = set()
        for result in me_success_results:
            for reaction in result.reactions:
                undeclared.update(set(reaction.reactants) - declared)
                undeclared.update(set(reaction.products) - declared)
        if undeclared:
            return (f"The output file(s) and the final network file '{file_name}' do not describe "
                    f"the same network: the output names species {sorted(undeclared)} which the "
                    f"network file never declares (it declares {sorted(declared)}). Arkane writes "
                    f"both artifacts from one in-memory network, so this pairing did not come "
                    f"from a single run.",)

        # Exactly one output.py reaches this function: _resolve_artifacts enforces
        # expected_outputs = {'output.py'} (see below) before either me_success_results or this
        # method are ever built, and explore() passes final_network_paths[0] as the sole final
        # network. A different count here is a programming error in the caller, not a runtime
        # condition this method should paper over with speculative multi-output aggregation.
        assert len(me_success_results) == 1, (
            f"_check_final_network_payload expects exactly one output.py result (enforced by "
            f"_resolve_artifacts' single-output requirement), got {len(me_success_results)}.")
        parsed_reactions = me_success_results[0].reactions

        expected_count = network.expected_net_reaction_count()
        actual_count = len(parsed_reactions)
        if actual_count != expected_count:
            direction = 'fewer' if actual_count < expected_count else 'more'
            return (f"The output file names {actual_count} pdepreaction(...) entries, {direction} "
                    f"than the {expected_count} the final network file '{file_name}' implies "
                    f"(PDepNetwork.expected_net_reaction_count() simulates Arkane's own "
                    f"duplicate-suppressing enumeration, so this count is exact, not an upper "
                    f"bound). Fewer means the output was truncated partway through the solve; "
                    f"more means this pairing is not from a single run either.",)

        # The count above is necessary but not sufficient, and so is asking merely whether every
        # declared channel is touched SOMEWHERE: an output can omit one expected channel pair and
        # duplicate another, keeping the total exact AND still touching every channel, while
        # describing a different topology than the network file declares (Codex, round-29 P1 A).
        # Identity therefore means set equality against the pairs Arkane's own enumeration would
        # write. That subsumes per-channel coverage -- a declared channel occurs in at least one
        # expected pair, so its absence from the output shows up as a missing pair -- which is why
        # a separate coverage pass is deliberately NOT kept here: two overlapping guards derived
        # from the same prediction can only disagree by being wrong.
        # Compared as MULTISETS, so that a repeated pair is reported as a surplus rather than
        # silently collapsing into the expected set: Arkane's duplicate suppression never writes the
        # same pair live twice, so a repeat is evidence in its own right.
        expected_pairs = Counter(network.expected_net_reaction_channel_pairs())
        actual_pairs = Counter(canonical_channel_pair(reaction.reactants, reaction.products)
                               for reaction in parsed_reactions)
        missing_pairs = sorted((expected_pairs - actual_pairs).elements())
        surplus_pairs = sorted((actual_pairs - expected_pairs).elements())
        if missing_pairs or surplus_pairs:
            faults = list()
            if missing_pairs:
                faults.append(f"it never connects the channel pair(s) {missing_pairs}")
            if surplus_pairs:
                faults.append(f"it connects the channel pair(s) {surplus_pairs}, which the network "
                              f"file does not imply (a pair listed here that the network DOES imply "
                              f"is one the output names more than once)")
            return (f"The output file and the final network file '{file_name}' describe different "
                    f"network topologies: {'; and '.join(faults)}. "
                    f"PDepNetwork.expected_net_reaction_channel_pairs() replays Arkane's own "
                    f"duplicate-suppressing enumeration, so the set of connected channel pairs is "
                    f"exact -- a total count that matches while the pairs do not means this output "
                    f"was not written from this network in a single complete run.",)
        return tuple()

    def _resolve_artifacts(self) -> tuple:
        """
        Rules 2-4: resolve the output/final-network artifact set by index, and enforce exact
        index-set equality between them.

        Returns:
            tuple: (reasons (tuple), final_network_paths (tuple), output_paths (tuple)). On any
                  failure, ``final_network_paths``/``output_paths`` are empty even if some
                  artifacts were found -- a partially-resolved set must never be handed back as
                  though it were complete.
        """
        reasons = []
        final_dir = os.path.join(self.output_directory, 'pdep', 'final')
        final_files_by_index = dict()
        if os.path.isdir(final_dir):
            for name in sorted(os.listdir(final_dir)):
                match = _FINAL_NETWORK_FILENAME_RE.match(name)
                if match:
                    index = int(match.group(1))
                    kind = match.group(2)
                    final_files_by_index.setdefault(index, dict())[kind] = os.path.join(final_dir, name)

        network_indices = sorted(final_files_by_index.keys())
        if not network_indices:
            reasons.append(f"No 'network<i>_full.py' final network artifact was found under "
                           f"'{final_dir}'; the exploration produced no final network to resolve.")
            return tuple(reasons), tuple(), tuple()

        if len(network_indices) > 1:
            # Settled verdict (design section 8.2): refuse a multi-network run in this increment.
            # Detection (network_indices, above) happens regardless -- refusal rides atop
            # resolution, it does not replace it.
            reasons.append(f"Arkane explored {len(network_indices)} distinct networks "
                           f"(indices {network_indices}) from seed {self.seed_species}; "
                           f"multi-network exploration is refused in this increment, since there is "
                           f"no single unambiguous result to resolve.")
            return tuple(reasons), tuple(), tuple()

        index = network_indices[0]
        # Rule 4: a reduction filter was requested whenever energy_tol was given (finite, as
        # opposed to this API's "unset" sentinel of None) or flux_tol was given and is > 0. In
        # that case 'network<i>_full.py' is copied BEFORE the reduction loop and is therefore not
        # evidence of completeness on its own -- 'network<i>_reduced.py' is written AFTER it,
        # unconditionally, whenever the reduction loop actually ran.
        # This test mirrors Arkane's own gate EXACTLY, verbatim:
        #     if self.energy_tol != np.inf or self.flux_tol != 0.0:   # arkane/explorer.py:303
        # The earlier version of this line used `flux_tol > 0` while its own comment claimed to
        # quote Arkane -- and the quote itself was wrong. A finite NEGATIVE flux_tol makes Arkane
        # run the reduction (it is != 0.0) while T3 would have gone on requiring the PRE-reduction
        # 'network<i>_full.py', accepting a network that does not correspond to the kinetics in
        # output.py. Mirroring means mirroring, including the cases that look like nonsense: it is
        # Arkane's behavior, not Arkane's intent, that decides which artifact exists on disk.
        # `None` is this API's "kwarg omitted" sentinel, which makes Arkane apply its own defaults
        # (energy_tol=np.inf, flux_tol=0.0, arkane/input.py:496) -- i.e. exactly the no-filter
        # values -- so None must be substituted before the comparison, not special-cased around it.
        # This is deliberately NOT left to rely on the writer's refusal of non-finite tolerances:
        # an infinite energy_tol means "no filter" to Arkane, so treating it as a request would
        # demand a 'network<i>_reduced.py' that Arkane will never write and REJECT a good run. Two
        # guards in two modules that only compose correctly by accident is how this branch has
        # produced defects before, so this one is correct on its own terms.
        energy_tol = math.inf if self.energy_tol is None else self.energy_tol
        flux_tol = 0.0 if self.flux_tol is None else self.flux_tol
        reduction_requested = energy_tol != math.inf or flux_tol != 0.0
        required_kind = 'reduced' if reduction_requested else 'full'
        final_network_path = final_files_by_index[index].get(required_kind)
        if final_network_path is None:
            reasons.append(
                f"A reduction filter was requested (energy_tol={self.energy_tol}, "
                f"flux_tol={self.flux_tol}), so completeness requires 'network{index}_reduced.py', "
                f"but only found {sorted(final_files_by_index[index].keys())} under '{final_dir}'."
                if reduction_requested else
                f"Expected 'network{index}_full.py' under '{final_dir}', but only found "
                f"{sorted(final_files_by_index[index].keys())}.")
            return tuple(reasons), tuple(), tuple()

        # Rule 2-3: resolve the expected output file set by index (single network -> plain
        # 'output.py' only; a non-empty 'output.py' is never itself evidence of a single-network
        # run, since kinetics(...) jobs append to it regardless of network count), and require
        # EXACT equality against what is actually on disk -- never "at least one".
        expected_outputs = {'output.py'}
        actual_outputs = {name for name in os.listdir(self.output_directory)
                          if _OUTPUT_FILENAME_RE.match(name)}
        if actual_outputs != expected_outputs:
            missing = expected_outputs - actual_outputs
            extra = actual_outputs - expected_outputs
            detail = []
            if missing:
                detail.append(f'missing {sorted(missing)}')
            if extra:
                detail.append(f'unexpected extra {sorted(extra)}')
            reasons.append(f"The resolved output file set {sorted(actual_outputs)} does not exactly "
                           f"match the expected single-network set {sorted(expected_outputs)}: "
                           f"{'; '.join(detail)}.")
            return tuple(reasons), tuple(), tuple()

        output_paths = tuple(sorted(os.path.join(self.output_directory, name)
                                    for name in actual_outputs))
        return tuple(), (final_network_path,), output_paths

    def _build_manifest(self, input_file_path: str, arkane_log_path: str) -> dict:
        """
        Rule 6: record path/size/sha256 for every artifact this run produced, plus the RMG-Py
        revision if obtainable. Explorer output is a function of RMG database contents, so a run
        with no provenance record is not replayable.

        Args:
            input_file_path (str): The path to the Arkane explorer input file this run used.
            arkane_log_path (str): The path to Arkane's own log for this run, if it exists.

        Returns:
            dict: The run manifest.
        """
        artifact_paths = [input_file_path, *self.output_paths, *self.final_network_paths]
        if os.path.isfile(arkane_log_path):
            artifact_paths.append(arkane_log_path)
        artifacts = [
            {'path': path, 'size': os.path.getsize(path), 'sha256': hash_file(path)}
            for path in artifact_paths
        ]
        return {'rmgpy_revision': _get_rmgpy_revision(), 'artifacts': artifacts}

    def get_networks(self) -> tuple:
        """
        Obtain the explored network(s) -- the resolved 'pdep/final/network<i>_(full|reduced).py'
        path(s).

        Raises:
            RuntimeError: If ``explore()`` has not been called yet, or did not succeed.

        Returns:
            tuple: The resolved final network file path(s).
        """
        if self.succeeded is None:
            raise RuntimeError('get_networks() was called before explore(); no exploration has '
                              'been attempted yet.')
        if not self.succeeded:
            raise RuntimeError(f'get_networks() was called after a failed explore(); the '
                              f'exploration did not succeed: {"; ".join(self.reasons)}')
        return self.final_network_paths

    def get_k_tp(self):
        """
        Obtain the parsed k(T,P) results for the explored network(s).

        DIRECTION IS ARKANE'S, NOT THE CALLER'S. Each entry carries the ``reactants``/``products`` as
        Arkane wrote them in ``output.py``, which is a detail of the solver's net-reaction ordering
        and is NOT guaranteed to match the direction of anything the caller asked about. A consumer
        that needs a rate in a particular direction must resolve the direction itself and reverse the
        entry when required; treating ``reactants`` as "the direction I requested" silently yields the
        reverse rate, which is a wrong number rather than an error.

        This is deliberate on both sides and neither side should be "fixed" into the other. The
        network-identity check is direction-INSENSITIVE by design (``canonical_channel_pair``,
        t3/pdep/parser.py:106, reasoned at :223-231), because which side of a net reaction is the
        reactant side is not part of the network's identity -- so a reversed entry legitimately passes
        identity (pinned by ``test_an_entry_written_in_the_reverse_direction_is_still_accepted``).
        Conversely the entries here are NOT canonicalized, because a rate coefficient genuinely is
        directional and discarding that would destroy information no later step can recover. What was
        missing was neither behaviour but this statement of the contract between them.

        Raises:
            RuntimeError: If ``explore()`` has not been called yet, or did not succeed.

        Returns:
            tuple: The parsed ``PDepArkaneReaction`` entries (from ``t3.pdep.parser``), each directed
                   as Arkane wrote it. See the direction note above before using ``reactants``.
        """
        if self.succeeded is None:
            raise RuntimeError('get_k_tp() was called before explore(); no exploration has been '
                              'attempted yet.')
        if not self.succeeded:
            raise RuntimeError(f'get_k_tp() was called after a failed explore(); the exploration '
                              f'did not succeed: {"; ".join(self.reasons)}')
        reactions = tuple()
        for result in self._me_success_results:
            reactions += result.reactions
        return reactions


register_explorer_adapter("Arkane", ArkaneExplorerAdapter)
