# Vendored QM logs for the P-dep hybrid tests

## `TS_nC3H7-iC3H7_freq.out`

A real, completed **Gaussian** frequency-job output for a genuine transition state: the
`nC3H7 -> iC3H7` intramolecular H-migration TS (species `C3H7`, doublet, one imaginary
frequency, normal termination).

Copied verbatim from ARC's own test data, `arc/testing/freq/TS_nC3H7-iC3H7.out`. The `.out`
extension is kept (ARC's own, and matching the other tracked QM logs under
`tests/data/pdep_hybrid/arc_ts/logs/`), so the repo's `*.log` gitignore rule does not swallow it.

Why it lives here: the fabricated ARC project in
`tests/test_pdep/test_pes_qm.py::TestArcQmRunnerRealCaptureHybrid` drives the DOF-conformer
extraction (`t3.pdep.dof_conformers`), which shells into `rmg_env` and loads this artifact's
`Log(...)` reference through Arkane. Arkane's `ess_factory` must be able to identify and parse the
referenced log, so a placeholder text file no longer suffices — the fabricated project needs a log
that is faithful to what a real ARC run leaves on disk. It contains only C and H atoms, so the
CHONS atom-energy corrections the tests supply cover it.
