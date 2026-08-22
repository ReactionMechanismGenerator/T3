# pdep_barrierless_real

A frozen copy of one real, explored potential-energy-surface network, kept so the barrierless gate
can be tested against data T3 actually produced rather than a hand-written double that merely
resembles it. This project has a standing rule that at least one test per external seam drives real
data, because fully-mocked seams have hidden every serious defect found here so far.

## `network0_full_cho2.py`

- **What it is.** The final explored network for the CHO2 (H + CO₂) surface, produced by a
  standalone Arkane explorer run seeded from the real RMG database (ticket I-007, report R-003).
- **Copied from.**
  `examples/pes_loop_cho2/explorer_run/pdep/final/network0_full.py`, in a worktree checked out on
  the `i006-cho2-scope` branch, on 2026-08-22. That source is an *untracked generated artifact* in a
  worktree rather than a committed file, so it can
  be regenerated out from under this test at any time — which is exactly why it is vendored frozen
  here. sha256 of the vendored copy:
  `66460d1fb2e3b4e9e2fc7736d9c05d36c109bedd906ae707537fa7b2f8abd4bb`.
- **Why it is the right fixture.** It carries nine TS-bearing path reactions spanning six RMG
  families, two of which are barrierless-by-name (`R_Recombination`, `Birad_R_Recombination`) and
  one of which (`R_Addition_COm`) is a near-miss the gate must NOT skip. Driven through
  `split_qm_candidates`, it exercises the `Birad_R_Recombination` gap this fixture was added for:
  before the family was added to `BARRIERLESS_FAMILIES` the network yields **7** QM candidates;
  after, **5**, with the two `Birad_R_Recombination` channels moved to `skipped`.

Note: an earlier generation of this network (measured in R-003) had ten path reactions and three
`Birad_R_Recombination` channels; the regenerated file vendored here has nine and two. The gap and
the fix are unchanged — only the exact leak count differs.
