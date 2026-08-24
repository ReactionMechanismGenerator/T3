# `positional_ts/network0_reduced.py`

A byte-for-byte copy (sha256
`46c152d162e789fcb003d4a79af57fa96eec6a07ca5e0df69757855ea4f43d96`) of the hybrid network file the
PES exploration loop wrote in round 0 of the `r001_m1-cho2-pilot` run and then crashed reading in
round 1. T3's own hybrid writer (`t3/pdep/hybrid.py`) produced it; the copy is vendored here because
the run tree under `/home/alon/runs/` is scratch and will be pruned.

Its load-bearing line is:

```python
transitionState('TS2', 'qm/TS2.py')
```

TS2 was adopted from quantum chemistry, so the writer rewrites it to the two-positional
`transitionState(label, path)` spelling — the ONLY form Arkane (`arkane/input.py:241`,
`if len(args) == 1 and len(kwargs) == 0`) accepts for loading a stat-mech file from a path; there is
no keyword spelling for it. The pdep parser's `_call_keywords` refused every positional argument and
crashed on this file. `tests/test_pdep/test_parser.py::TestPositionalTransitionState` drives this
file through `parse_pdep_network_file` and asserts the parser recovers both TS2's label and its path.
Do not regenerate or "clean up" this file: its value is that it is exactly what a real run emitted.
