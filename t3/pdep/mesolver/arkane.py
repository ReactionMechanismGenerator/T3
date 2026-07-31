"""
t3 pdep mesolver arkane module.

Concrete ``MESolverAdapter`` implementation that drives Arkane to perform a plain (non-SA)
master-equation solve of a single P-dep network and parse its resulting k(T,P) kinetics.
"""

import os

from t3.logger import Logger
from t3.pdep.me_success import check_arkane_me_success
from t3.pdep.mesolver.adapter import MESolverAdapter
from t3.pdep.mesolver.factory import register_mesolver_adapter
from t3.runners.rmg_runner import run_arkane_job
from t3.utils.writer import METHOD_MAP, write_arkane_network_input_file


class ArkaneMESolverAdapter(MESolverAdapter):
    """
    A MESolverAdapter that solves a P-dep network's master equation using Arkane.

    Unlike MESS or Mesmer, Arkane does not require every node in the network to be
    QM-computed, so it can consume a network in which some path reactions are still RMG/ILT
    estimates (``allow_ilt_complement=True``).

    Note:
        ``set_up()`` runs from ``__init__`` (mirroring ``t3/simulate/cantera_base.py``'s
        pattern), so merely CONSTRUCTING an adapter writes ``input.py`` into
        ``output_directory``. Two adapters sharing one ``output_directory`` will clobber each
        other's input/output files; ``output_directory`` must therefore be unique per adapter
        instance. This is a known, deliberately-unfixed sharp edge of the mirrored pattern.

    Args:
        network_path (str): The path to the RMG P-dep network file to solve.
        output_directory (str): The path to the directory in which to write the Arkane
                                input/output files.
        method (str): The approximation method to use for solving the master equation.
                     One of ``'CSE'``, ``'RS'`` or ``'MSC'`` (see ``t3.utils.writer.METHOD_MAP``).
        logger (Logger, optional): The current T3 Logger instance.
        allow_ilt_complement (bool, optional): Whether the network being solved is only
                                              partially QM-computed, with the remaining path
                                              reactions estimated via RMG/ILT.
        expected_reactions (int | set, optional): The expected coverage of ``pdepreaction(...)``
                                                 entries in Arkane's output: either an exact
                                                 count (int) or a set of
                                                 ``(reactants, products)`` label-tuple pairs
                                                 that must all be present. When ``None`` (the
                                                 default), completeness of the output is NOT
                                                 verified: a truncated ``output.py`` whose
                                                 entries are all finite is still accepted.
                                                 Deriving this automatically from the network
                                                 topology is deferred future work (the number
                                                 of net reactions in a P-dep network is not
                                                 trivially the number of path reactions).

    Raises:
        ValueError: If ``method`` is not a key of ``t3.utils.writer.METHOD_MAP``.

    Attributes:
        isomer_labels (tuple): The isomer labels of the network, as scraped from the RMG
                               network file while writing the Arkane input file.
        me_success_result (MESuccessResult): The outcome of the most recent ``solve()`` call.
                                             Not set until ``solve()`` has been called at least
                                             once.
    """

    supports_ilt_complement = True

    def __init__(self,
                 network_path: str,
                 output_directory: str,
                 method: str,
                 logger: Logger | None = None,
                 allow_ilt_complement: bool = False,
                 expected_reactions: int | set | None = None,
                 ):
        self.network_path = network_path
        self.output_directory = output_directory
        self.method = method
        self.logger = logger
        self.allow_ilt_complement = allow_ilt_complement
        self.expected_reactions = expected_reactions

        self.isomer_labels = tuple()
        self.me_success_result = None

        self.set_up()

    def set_up(self):
        """
        Create the output directory (if needed) and write the Arkane ME input file.

        The input file is written WITHOUT a ``sensitivity_conditions`` directive: this is a
        plain ME solve, not a sensitivity-analysis run.

        This is called from ``__init__`` (the same pattern as
        ``t3/simulate/cantera_base.py``), so constructing the adapter already writes
        ``input.py``; see the class docstring for the resulting requirement that
        ``output_directory`` be unique per adapter instance.

        Raises:
            ValueError: If ``self.method`` is not a key of ``t3.utils.writer.METHOD_MAP``.
        """
        if self.method not in METHOD_MAP:
            raise ValueError(f"The 'method' argument must be one of {list(METHOD_MAP.keys())}, "
                             f"got '{self.method}'.")
        if not os.path.isdir(self.output_directory):
            os.makedirs(self.output_directory)
        input_file_path = os.path.join(self.output_directory, 'input.py')
        write_result = write_arkane_network_input_file(source_path=self.network_path,
                                                        dest_path=input_file_path,
                                                        method=self.method,
                                                        sensitivity=False,
                                                        )
        self.isomer_labels = write_result.isomer_labels

    def solve(self):
        """
        Run Arkane and determine whether the ME solve genuinely succeeded.

        Arkane pre-creates ``output.py`` before running any job and can exit 0 with empty
        stderr while writing a syntactically valid but numerically unusable result (e.g., a
        rejected chemically-significant-eigenvalues solve). Bare completion of the subprocess is
        therefore not evidence of success: the actual criterion is
        ``t3.pdep.me_success.check_arkane_me_success``, which parses ``output.py`` and requires
        every kinetics parameter to be finite.

        Returns:
            bool: Whether the ME solve genuinely succeeded.
        """
        input_file_path = os.path.join(self.output_directory, 'input.py')
        output_path = os.path.join(self.output_directory, 'output.py')
        stderr_path = os.path.join(self.output_directory, 'stderr.log')

        # Delete any stderr.log left behind by a previous (failed) run in this directory before
        # invoking Arkane, mirroring run_arkane_job's delete-then-require handling of the
        # artifact itself: a stale stderr.log read after a clean re-run would otherwise fail a
        # perfectly good solve.
        if os.path.isfile(stderr_path):
            os.remove(stderr_path)

        job_ran = run_arkane_job(input_file=input_file_path,
                                 output_directory=self.output_directory,
                                 logger=self.logger,
                                 required_artifact='output.py')

        stderr_text = None
        if os.path.isfile(stderr_path):
            with open(stderr_path, 'r') as f:
                stderr_text = f.read()

        # ``run_arkane_job`` reports a bool rather than a raw exit status, so it is folded in as a
        # 0/1 exit code. It is weaker evidence than the payload check below (Arkane exits 0 on the
        # silent-CSE failure), but it is not redundant either: it catches an ARC-level failure that
        # still happens to leave a plausible-looking output.py behind.
        self.me_success_result = check_arkane_me_success(output_path=output_path,
                                                         exit_code=0 if job_ran else 1,
                                                         stderr=stderr_text,
                                                         expected_reactions=self.expected_reactions,
                                                         )

        if not self.me_success_result.succeeded and self.logger is not None:
            for reason in self.me_success_result.reasons:
                self.logger.error(reason)

        return self.me_success_result.succeeded

    def get_k_tp(self):
        """
        Obtain the solved k(T,P) results.

        DIRECTION IS ARKANE'S, NOT THE CALLER'S. Each entry carries the ``reactants``/``products`` as
        Arkane wrote them in ``output.py``; that ordering is not guaranteed to match the direction of
        anything the caller asked about, so a consumer needing a specific direction must resolve it
        and reverse the entry when required. The same contract holds for the explorer adapter's
        ``get_k_tp`` (t3/pdep/explorer/arkane.py), where the reasoning is written out in full.

        Raises:
            RuntimeError: If ``solve()`` has not been called yet, or was called and did not
                         succeed. A caller must never silently receive ``None``/NaN kinetics as
                         if they were genuine results.

        Returns:
            tuple: The parsed ``PDepArkaneReaction`` entries (from ``t3.pdep.parser``).
        """
        if self.me_success_result is None:
            raise RuntimeError('get_k_tp() was called before solve(); no ME solve has been attempted yet.')
        if not self.me_success_result.succeeded:
            raise RuntimeError(f'get_k_tp() was called after a failed solve(); the ME solve did not succeed: '
                               f'{"; ".join(self.me_success_result.reasons)}')
        return self.me_success_result.reactions


register_mesolver_adapter("Arkane", ArkaneMESolverAdapter)
