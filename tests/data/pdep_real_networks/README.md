# Real PDep network fixtures

Two complete, unmodified chains from a real T3 run, vendored so the PDep qualification gate can be
tested against data it did not generate. Every file here is a byte-for-byte copy of what that run
produced; none of it was hand-written, trimmed, or regenerated to make a test pass. Three of the
four files per network are RMG/Arkane output -- the fourth, `t3_sa_cache.yml`, is T3's own sidecar,
written by T3 at the end of that run.

Each network directory holds the three files `select_pdep_network()` reads, plus `arkane.log`,
which it does not read and which is here only as provenance evidence for a human or for
`test_the_vendored_arkane_log_names_the_rmg_py_that_carries_the_ilt_fixes`:

```
<network_id>/
    <network_id>.py                  the RMG PDep network file
    arkane.log                       Arkane's own log for this run (provenance, see below)
    sensitivity/
        sa_coefficients.yml          Arkane's sensitivity output
        t3_sa_cache.yml              T3's sidecar for that output
```

The sidecar must sit beside `sa_coefficients.yml` because `sa_cache_metadata_path()` derives it
from that file's directory. `validate_sa_cache()` gates on the SA-cache contract version, the
perturbation, the ME method, and **content hashes** of the network and SA files -- it compares no
absolute paths, which is why these files still validate (`cache_status == 'cached_valid'`) after
being copied out of the run directory that produced them.

## The two networks are a deliberate pair

| | `network21_1` | `network799_1` |
|---|---|---|
| ME method | CSE | MSC |
| selected transition states | 29 | 20 |
| uncertain transition states | **0** | **20 (all)** |
| strongest response, relative to the coefficient floor | **4190x** | 796x |
| `qualified` | **False** | **True** |

They are kept together because they disagree in the one way that matters: `network21_1` is the
*more* sensitive network of the two by a factor of five, and it still does not qualify, because
none of its path reactions carry uncertain kinetics. Qualification is
sensitivity **times** uncertainty, and this pair is the evidence -- a fixture built by hand would
almost certainly make sensitivity and uncertainty agree, and would therefore pass just as happily
against a gate that ignored uncertainty entirely.

`network21_1` also exercises the coefficient floor from both sides in a single record: its
strongest selected transition state is 4190x the floor while its weakest is 2.2x, so the selection
is neither trivially everything nor trivially nothing.

## Provenance, stated exactly

Both `arkane.log`s record the RMG-Py commit that produced them:

```
e720866ae94eca51652978c15a0fb33c6827be67
```

That commit contains both ILT sensitivity fixes (`ac4fbf1d5` and `e720866ae`). This matters: an
earlier T3 experiment produced transition-state coefficients nine orders of magnitude below the
floor and looked merely *negative* for a week, when in fact it had run against an Arkane predating
those fixes. The log is what distinguishes the two situations.

**The `t3_sa_cache.yml` sidecars here record `rmg_py_commit: null`.** They are not defective and
they have not been edited -- they were simply written before T3 learned to read the commit out of
Arkane's log. The `arkane.log` in each directory is the first-hand witness, and
`tests/test_pdep/test_real_networks.py` asserts the commit against it rather than against the
sidecar. That null is itself asserted by a test, so do not "fix" the sidecars by filling the field
in: their hashes are what make the cache validate, and any commit written there by hand would be
invented provenance rather than recorded provenance.

**What the log does not prove.** Nothing binds an `arkane.log` to the `sa_coefficients.yml` beside
it -- no hash, no shared run id -- and `validate_sa_cache()` never reads a log. The log is evidence
about where these files came from, and it is only as good as the fact that it was copied together
with them. It is not proof that this log produced this SA file.

## Source

`/home/alon/runs/T3-pdep-qm-trial-004/networks/`, a read-only prior experiment:

- `network21_1` <- `none_18r17_network21_1/` (CSE)
- `network799_1` <- `all_10r9_network799_1/` (MSC)
