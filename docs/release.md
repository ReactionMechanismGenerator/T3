# Release notes


## T3 0.2.0

Major changes since 0.1.0:

- Removed all runtime dependencies on the RMG-Py Python API. T3 now talks to
  RMG only via subprocesses and parses RMG outputs directly, making the
  `t3_env` independent of the RMG environment.
- Refactored the simulate adapter layer around a shared `CanteraBase` class.
  Adding a new Cantera reactor adapter now requires only a few lines (see
  [Tutorial 2](tutorials/2_adding_simulate_adapter.md)).
- New simulate adapters:
    - `CanteraConstantHP` — adiabatic, constant-pressure batch reactor.
    - `CanteraConstantUV` — adiabatic, constant-volume batch reactor (shock
      tubes, RCMs, constant-volume bombs).
    - `CanteraPFR` — Lagrangian isothermal/isobaric plug flow reactor.
    - `CanteraPFRTProfile` — non-isothermal PFR with a hard-coded `T(z)`
      profile, for flow-tube experiments with a measured wall temperature.
      Sensitivity analysis is supported and indexed by axial position `z`.
    - `CanteraJSR` — jet-stirred / PSR / CSTR adapter.
- New runnable examples for each new adapter under `examples/`.
- Documentation overhaul: rewritten how-to guide with an experimental-setup
  → adapter mapping table, expanded examples page, and a consolidated
  [input reference](input_reference.md) page (the old `examples/commented/`
  duplicate has been removed).


## T3 0.1.0

This is the first version of T3.

## Version style

T3 uses [Semantic Versioning](https://semver.org/)
