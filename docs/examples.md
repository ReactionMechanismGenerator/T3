# Examples

T3 ships with several example input files in the
[`examples/`](https://github.com/ReactionMechanismGenerator/T3/tree/main/examples) directory.
Each demonstrates a different feature or use case.

To run any example:

```bash
conda activate t3_env
cd ~/Code/T3/examples/<example_name>
python ~/Code/T3/T3.py input.yml
```


## Minimal

**Directory:** `examples/minimal/`

The simplest complete T3 input: H<sub>2</sub>/O<sub>2</sub> combustion with Cantera-based
sensitivity analysis and ARC for QM refinement. A good starting point for new users.

**Features demonstrated:**

- Basic species definition with SMILES and concentrations
- `CanteraConstantTP` adapter for sensitivity analysis
- `SA_observable` on H and OH
- ARC QM with a simple level of theory (`b3lyp/6-31g(d,p)`)
- Progressive `core_tolerance: [0.01, 0.001]`


## Cantera SA

**Directory:** `examples/cantera_sa/`

Demonstrates using a Cantera simulation adapter for sensitivity analysis
with global observables.

**Features demonstrated:**

- `CanteraConstantTP` adapter
- `global_observables: ['IDT']` for ignition delay time sensitivity
- Three-stage `core_tolerance: [0.01, 0.005, 0.001]`
- Species-level SA observables (OH, H)


## Shock tube

**Directory:** `examples/shock_tube/`

Idealized reflected-shock experiment for H<sub>2</sub>/O<sub>2</sub>/Ar.
A shock tube post-shock state is treated as a closed, adiabatic,
constant-volume batch reactor.

**Features demonstrated:**

- `CanteraConstantUV` adapter — closed adiabatic constant-V batch (same `IdealGasReactor` as a shock tube post-reflected-shock state). For IDT-specific sensitivity analysis use `CanteraIDT` instead.
- Dilute mixture in argon (`balance: true` on Ar, `reactive: false`)
- `global_observables: ['IDT']` — IDT is the primary measurement of a shock tube
- Short `termination_time` in milliseconds, matching typical IDT magnitudes
- Post-shock T = 1400 K, P = 2 bar


## Plug flow reactor (PFR)

**Directory:** `examples/pfr/`

An isothermal, isobaric flow reactor (a Lagrangian parcel of gas advected
at constant T and P) using the `CanteraPFR` adapter.

**Features demonstrated:**

- `CanteraPFR` adapter — Lagrangian PFR
- Residence time set via the reactor `termination_time`
- Reactor geometry (LENGTH, AREA, N_CELLS) controlled by module-level
  constants in `t3/simulate/cantera_pfr.py`
- Species-level SA observables (OH, H)


## PFR with prescribed T(z) profile

**Directory:** `examples/pfr_t_profile/`

A non-isothermal, isobaric flow reactor in which the gas temperature follows
a hard-coded axial profile `T(z)`.  Useful for flow-tube experiments where
the wall temperature is known a priori (e.g., a measured profile).
A single Lagrangian particle is advected through `N_CELLS` segments and the
gas temperature at each segment is overridden from the profile (energy
equation off); the velocity (and therefore the time step per segment) is
recomputed from the local density at each `T`.

**Features demonstrated:**

- `CanteraPFRTProfile` adapter — Lagrangian non-isothermal PFR
- Default profile: 80 cm reactor, raised-cosine ramp 300 → 900 K over 20 cm,
  isothermal plateau at 900 K for 40 cm, raised-cosine ramp 900 → 500 K over
  20 cm — `C¹` continuous everywhere
- Reactor geometry, breakpoints, and the `temperature_profile(z)` function
  are all controlled by module-level constants in
  `t3/simulate/cantera_pfr_t_profile.py` — advanced users either tune the
  constants in place or copy the file as the basis for a brand-new adapter
- Mass flow is constant (set from the inlet density and the
  `termination_time`); density and velocity evolve with `T`
- Sensitivity analysis is supported and indexed by axial position `z`
  (the reactor is steady-state in the lab frame).  `get_sa_coefficients()`
  adds a `'distance'` key to the standard `sa_dict` carrying the
  cell-center positions [m]


## Jet stirred reactor (JSR)

**Directory:** `examples/jsr/`

An isothermal, isobaric, well-mixed steady-state reactor with continuous
feed and outflow, using the `CanteraJSR` adapter.

**Features demonstrated:**

- `CanteraJSR` adapter — JSR / PSR / CSTR
- Inlet mass flow rate set so `mdot = m_reactor / tau`, with `tau` taken
  from the reactor `termination_time`
- Reactor volume controlled by the module-level `VOLUME` constant in
  `t3/simulate/cantera_jsr.py` (defaults to 100 cm³)
- Dilute fuel in nitrogen (`balance: true` on N<sub>2</sub>)
- Species-level SA observables (OH, H)


## Full QM workflow

**Directory:** `examples/with_qm/`

Shows the core T3 iterative workflow: RMG model generation followed by
QM refinement with ARC, repeated over multiple iterations.

**Features demonstrated:**

- Complete RMG + ARC iterative cycle
- `adaptive_levels` for selecting level of theory based on species size:
    - Small species (1-6 heavy atoms): high-level CCSD(T)-F12
    - Medium species (7-14): DLPNO-CCSD(T)
    - Large species (15+): DFT only
- Four-stage `core_tolerance` for progressive model refinement
- Species constraints to limit model size
- Balance species (N<sub>2</sub>)


## Ranged conditions

**Directory:** `examples/ranged_conditions/`

Demonstrates running sensitivity analysis across a grid of temperatures,
pressures, and concentrations.

**Features demonstrated:**

- Temperature range: `T: [800, 1500]` K
- Pressure range: `P: [1, 10]` bar
- Concentration ranges on fuel and oxidizer: `concentration: [0.10, 0.67]`
- Grid density control: `num_sa_per_temperature_range`, `num_sa_per_pressure_range`
- `modify_concentration_ranges_together: true` to vary species together
- `conditions_per_iteration: 12` for RMG


## Pressure dependence

**Directory:** `examples/pressure_dependence/`

Shows how to enable pressure-dependent network generation in RMG.

**Features demonstrated:**

- `pdep` block with `method: MSC` (modified strong collision)
- Chebyshev interpolation over a T-P grid
- `pdep_SA_threshold` for PES-level sensitivity analysis
- `ME_methods: ['CSE', 'MSC']` for master equation analysis
- `max_atoms: 16` controlling which species are included in PDep networks
- Species constraints to keep the model tractable


## Pre-QM (iteration 0)

**Directory:** `examples/pre_qm_iteration_0/`

Demonstrates computing thermodynamic properties for key species
*before* the first RMG iteration.

**Features demonstrated:**

- QM `species` block with radicals to compute upfront (n-propyl, i-propyl, propylperoxyl)
- `seed_all_rads: ['radical', 'alkoxyl', 'peroxyl']` on the fuel to auto-generate radical derivatives
- Pre-computed thermo available to RMG from iteration 1 onward
- CBS-QB3 composite method for iteration 0 species


## Liquid phase

**Directory:** `examples/liquid_phase/`

Demonstrates a liquid-phase simulation with a solvent species.

**Features demonstrated:**

- `liquid batch constant T V` reactor type
- Temperature range for the reactor
- `solvent: true` on the water species
- `reactive: false` on N<sub>2</sub> (non-reactive dissolved gas)
- `all_core_species: true` to calculate thermo for every core species
- Liquid-phase specific species constraints


## IDT with experimental comparison

**Directory:** `examples/idt_with_experiment/`

Demonstrates running IDT-focused sensitivity analysis with experimental data
comparison. T3 computes predicted IDTs and compares them against experimental
measurements, reporting per-point log-errors and RMSE.

**Features demonstrated:**

- `global_observables: ['IDT']` for IDT-based SA
- `idt_criterion: max_dOHdt` (configurable: `max_dTdt`, `max_radical_dt`)
- `idt_sa_method: adjoint` (`brute_force` or `adjoint` for Cantera built-in SA)
- `experimental_idt_path` pointing to a YAML file with experimental data
- Experimental YAML format with citation, T/P/phi/idt fields


## Input reference

For a comprehensive reference showing **every** available input argument with
inline documentation, see the [input reference](input_reference.md) page.
