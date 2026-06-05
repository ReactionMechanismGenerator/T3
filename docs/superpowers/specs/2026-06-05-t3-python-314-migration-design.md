# T3 → Python 3.14 Migration — Design / Spec

**Date:** 2026-06-05
**Branch:** `official-migration` (tracks `official/official-migration`)
**Goal:** Finalize T3's migration to Python 3.14, mirroring the already-completed
ARC migration (ARC is an in-process dependency of T3).

---

## 1. Context

This branch already did the heavy lifting of **removing RMG-Py API dependencies**
(cantera-based simulate adapters, `cantera_parser.py`, the `slim_rmg.py` shim, and
isolating remaining `rmgpy` imports into `runners/rmg_incore_*.py` subprocess
scripts that run in a separate `rmg_env`). See `MIGRATION_WORKPLAN.md`.

What remains is the **Python 3.14 runtime bump itself**, which has *not* happened
yet: `environment.yml` still pins `python >=3.12.*`, `pyproject.toml` still says
`requires-python = ">=3.8"`, and the working `t3_env` is Python 3.12.12.

## 2. Hard constraint that drives sequencing

T3 imports **ARC in-process** (same env) across ~10 modules — `ARCSpecies`,
`arc.molecule`, `arc.common`, `arc.species.perceive`, etc. ARC `main` is now
**Python 3.14-only**. Therefore `t3_env` itself must be 3.14, and there is **no
meaningful "green on 3.12" checkpoint** — the current ARC will not import there.

Consequence: the two migrations (RMG-decoupling finish + 3.14 bump) collapse into
one. We cannot cleanly isolate "is it decoupling or is it 3.14?"; we recover some
of that separation by triaging failures **by category** (3.14 syntax/stdlib
breakage reads differently from ARC-interface/decoupling gaps).

## 3. Key decisions (resolved)

| Topic | Decision |
|---|---|
| Python | Pin `=3.14` (env), `>=3.14` (pyproject/requirements), mirroring ARC |
| ARC source | `main` (3.14), installed editable in `t3_env` |
| RMG-Py source | `main` — already contains `rmgpy/data/auto_database.py` (`auto_select_libraries`); the old `rmgpy.data.libraries` stash hack is obsolete |
| RMG-Py runtime | Stays in a **separate Python-3.9 `rmg_env`**, invoked as a subprocess (unchanged) |
| Cantera | **Core T3 dep**, kept loose `>=3.2.0`. conda-forge has a `py314` build (`cantera 3.2.0 py314he8fa917_0`), so 3.14-compatible. (Unlike ARC, which dropped cantera entirely.) |
| cython / rdkit | Bump to `>=3.1` / `>=2026.03` to match ARC's 3.14 pins |
| Typing | **Full sweep** like ARC (`List→list`, `Optional[X]→X | None`, `Union→|`, drop dead `from __future__`) across ~21 files |
| Commits | Clean, logical commits with good messages; squash the 21 `TMP`/fixup commits before the PR to `official/main` |

## 4. Phased plan

Each phase ends at a stated gate.

### Phase 0 — Stand up the 3.14 environment (riskiest unknown)
- `environment.yml`: `python =3.14`, `cython >=3.1`, `rdkit >=2026.03`,
  keep `cantera >=3.2.0`, `numpy >=2.1.*`.
- `pyproject.toml`: `requires-python = ">=3.14"` (from stale `>=3.8`).
- `requirements.txt`: matching bumps.
- Create a fresh 3.14 `t3_env`; install **ARC `main`** editable + cantera;
  keep `rmg_env` (Python 3.9, RMG-Py `main`) unchanged.
- **Gate:** `import t3`, `import arc`, `import cantera` succeed;
  `pytest --co` collects all 170 tests.

### Phase 1 — Green the suite on 3.14
- Run the full suite; triage failures by category.
- The 4 previously stash-blocked tests should now pass via `auto_database`.
- Run the excluded flux tests.
- Fix genuine 3.14 runtime breakage.
- Update `MIGRATION_WORKPLAN.md`.
- **Gate:** full suite green on 3.14.

### Phase 2 — Typing modernization (full sweep, ~21 files)
- `List/Dict/Tuple/Optional/Union` → builtins / `X | None` / `X | Y`.
- Drop dead `from __future__ import`.
- Optionally add an ARC-style `check_python()` ≥3.14 gate.
- Re-run suite.

### Phase 3 — CI, packaging & commit hygiene
- `.github/workflows/ci.yml` → 3.14 env, ARC `main` + RMG-Py `main`;
  verify the `rmg_env` subprocess path works in CI.
- Squash the 21 `TMP`/fixup commits into clean logical commits before the
  PR to `official/main`.
- **Gate:** green CI.

## 5. Validation criteria
1. `python -c "import t3"` succeeds on Python 3.14.
2. `pytest tests/ --co` collects all 170 tests.
3. `pytest tests/ -q` — full suite green (incl. flux).
4. No `rmgpy` imports outside `runners/rmg_incore_*.py` and `slim_rmg.py`.
5. CI green on Python 3.14 against ARC `main` + RMG-Py `main`.
