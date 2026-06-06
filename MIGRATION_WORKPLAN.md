# T3 Migration Workplan: Removing RMG-Py API Dependencies

## 1. Current State (updated 2026-03-28)

### 1.1 Test Suite Status

**147 passed, 4 failed** (excl. flux tests; 19 deselected).

All 4 remaining failures are blocked on a single RMG-Py issue:
`rmgpy.data.libraries` is in the git stash, not in the working tree.

The `py_rdl` / import-time blocker is **resolved** — PyRDL is now installed
into `t3_env` via `devtools/install_pyrdl.sh`.

### 1.2 Remaining Failures

| Test | Root Cause | Blocked On |
|------|-----------|------------|
| `test_main::test_determine_species_from_pdep_network` | Arkane PDep SA produces negative rates for test network | RMG-Py stash (`rmgpy.data.libraries`) |
| `test_functional::test_computing_thermo` | RMG subprocess crashes at import | RMG-Py stash |
| `test_simulate::test_set_up_no_sa` | RMG SA subprocess crashes at import | RMG-Py stash |
| `test_simulate::test_get_sa_coefficients` | RMG SA subprocess crashes at import | RMG-Py stash |

**Single action to unblock all 4:**
```bash
cd /home/alon/Code/RMG-Py && git stash pop
```
This restores `rmgpy/data/libraries.py` which is imported at
`rmgpy/rmg/main.py:566` (`from rmgpy.data.libraries import auto_select_libraries`).

---

## 2. What Was Done

### Phase 0: Unblock Tests ✅

- Installed PyRDL into `t3_env` — all 170 tests now collect and run.
- `devtools/install_pyrdl.sh` created as thin wrapper → ARC's script with `t3_env`.
- `devtools/install_all.sh` cleaned: removed all Julia/RMS/pyrms code,
  calls `install_pyrdl.sh` at the end.
- Deleted stale `devtools/install_arc.sh` and `devtools/install_pyrms.sh`.
- Updated `Makefile`: removed `install-arc`/`install-pyrms`, added `install-pyrdl`.

### Phase 1: Fix T3Species ↔ ARCSpecies interface ✅

- `T3Species.__init__` strips T3-only kwargs via `_T3_ONLY_KWARGS` before `super().__init__()`
- Explicit `qm_label`/`rmg_label` passed as kwargs override computed values
- `label=None` handled: derives label from adjlist header or SMILES
- Test labels updated for ARC legalization (`OH(4)` → `OH[4]`, `C#C` → `CtC`)
- `test_common.py`: `t3_index` → `key` (correct T3Species field name)

### Phase 2: Fix integration logic ✅

- **Cantera parser overhaul** (`t3/utils/cantera_parser.py`):
  - `yaml_label_to_species` dict maps original YAML labels → T3Species (handles ARC label legalization)
  - `re.sub` strips `(+M)`/`(+N2)` collider notation before splitting reaction equations
  - Bare `M` bath gas stripped from reactant/product labels
  - `is_pressure_dependent` flag detected from `(+` pattern and YAML reaction type
  - Duplicate reaction deduplication (Chemkin DUPLICATE entries share equations)
  - `index=i` passed to T3Reaction constructor
- **Label normalization** (`t3/common.py`): `get_species_by_label` tries bracket-normalized labels
- **Pressure-dependent check** (`t3/main.py:715`): `getattr(reaction, 'is_pressure_dependent', False)`
- **`from_dict` round-trip** (`t3/chem.py`): `qm_label` now preserved through dump/reload cycle

### Phase 3: Fix test assertions & fixtures ✅

- Test fixture `test_cantera_parser.yaml` updated to use actual RMG note format
- Reaction indices updated for dedup (29→24, 75→67, 113→102, 140→126, 162→146)
- Species label assertions updated for ARC legalization (`H2O` → `H2O[7]`)
- Reaction string assertions use `in` check for ARCReaction repr format
- `qm_label` assertion corrected for plain-label species

### Phase 4: Fix RMG subprocess compat ✅

- **`rmg_incore_script.py`**: Removed `ILPSolutionError` from except block —
  no longer exists in current RMG-Py exceptions module.
- **`tests/data/functional_2_thermo/input.yml`**: Fixed project name
  (`functional_1_thermo_1_rate` → `functional_2_thermo`).
- **`test_computing_thermo`**: Reset `T3Species._index_counter` at test start;
  species convergence check uses regex for robustness against counter state.

### Installation scripts (ARC side)

- `ARC/devtools/install_all.sh`: fixed pyrdl call to use `run_devtool`.
- `ARC/devtools/install_pyrdl.sh`: accepts optional env name `$1` (default: `arc_env`).

---

## 3. Recommended Next Steps

### Step 1: Restore RMG-Py stash (5 min) — UNBLOCKS ALL 4 TESTS

```bash
cd /home/alon/Code/RMG-Py && git stash pop
```

Then verify:
```bash
conda run -n rmg_env python -c "from rmgpy.data.libraries import auto_select_libraries; print('OK')"
```

Re-run the 4 failing tests:
```bash
pytest tests/test_main.py::test_determine_species_from_pdep_network \
       tests/test_functional.py::test_computing_thermo \
       tests/test_simulate/test_rmg_constant_tp.py::test_set_up_no_sa \
       tests/test_simulate/test_rmg_constant_tp.py::test_get_sa_coefficients \
       --tb=short -q
```

### Step 2: Run flux tests (~20 min)

The flux tests were excluded during development because they're slow.
Run them to confirm no regressions:
```bash
pytest tests/test_main.py -k "flux" --tb=short -q
```

### Step 3: Clean up and commit

Once all tests pass:
1. Run `pytest tests/ --tb=short -q` for the full green suite
2. Commit the T3 migration changes on the `migration` branch
3. Open a PR against `main`

### Step 4 (optional): PDep network test data

If `test_determine_species_from_pdep_network` still fails after restoring the
stash (Arkane numerical issue with CSE/MSC for network4_2), the test data may
need refreshing — the PDep network input files in
`tests/data/pdep_network/iteration_1/RMG/pdep/` may be stale relative to
the current Arkane version.

---

## 4. Validation Criteria

Full green:
1. `python -c "import t3"` succeeds
2. `pytest tests/ --co` collects all 170 tests
3. `pytest tests/ -q` — **all 151 pass** (incl. flux)
4. No rmgpy imports outside `runners/rmg_incore_*.py` and `slim_rmg.py`

---

## 5. Architecture Notes

- **`t3/utils/slim_rmg.py`** — shim classes replacing rmgpy data structures
  (Reaction, PDepReaction, PDepNetwork, etc.) for simulation adapters and
  library I/O.
- **`t3/runners/rmg_incore_script.py`** and **`rmg_incore_sa.py`** — run in
  `rmg_env` as subprocesses. Their rmgpy imports are correct and intentional.
- **RDKit** (`rdkit >= 2025.03`) — declared dependency for molecular operations.
- **Julia/RMS/pyrms** — fully removed from T3's installation.
