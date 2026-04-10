
[![CI](https://github.com/ReactionMechanismGenerator/T3/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/ReactionMechanismGenerator/T3/actions/workflows/ci.yml)
[![Docs](https://github.com/ReactionMechanismGenerator/T3/actions/workflows/gh-pages.yml/badge.svg)](https://reactionmechanismgenerator.github.io/T3/)
[![codecov](https://codecov.io/gh/ReactionMechanismGenerator/T3/branch/main/graph/badge.svg)](https://codecov.io/gh/ReactionMechanismGenerator/T3)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
![python](https://img.shields.io/badge/Python-3.12+-blue.svg)

# The Tandem Tool (T3) for automated chemical kinetic model development

**T3** automates the development of detailed chemical kinetic models.
Given a set of initial species and conditions, it produces a validated model with
high-fidelity thermochemistry and rate parameters by iteratively combining
[RMG](https://rmg.mit.edu/) for model generation,
[Cantera](https://cantera.org/) for sensitivity analysis,
and [ARC](https://reactionmechanismgenerator.github.io/ARC/) for quantum mechanical calculations.

![T3 scheme][cycle]

[cycle]: /grf/T3-circle.gif "T3 scheme"

---

## Features

- **Iterative refinement:** runs RMG, identifies sensitive species and reactions via SA, computes
  accurate thermo/kinetics with ARC, and feeds the results back into the next iteration
- **Multiple SA adapters:** constant TP/HP/UV, PFR, JSR, and brute-force ignition delay time (IDT)
- **Equivalence-ratio sweeps:** fuel/oxidizer/diluent role taxonomy with automatic concentration
  computation across a grid of temperatures, pressures, and equivalence ratios
- **Uncertainty analysis:** local and global UA via polynomial chaos expansion
- **Restart support:** T3 can resume from any interrupted iteration

---

## Installation

### Prerequisites

- Python 3.12+
- [Conda](https://docs.conda.io/en/latest/) (we recommend [Miniforge](https://github.com/conda-forge/miniforge))

### Quick start

```bash
git clone https://github.com/ReactionMechanismGenerator/T3.git
cd T3
make install
conda activate t3_env
```

This clones and installs ARC, RMG-Py, RMG-database, and all external dependencies.
See the [full installation guide](https://reactionmechanismgenerator.github.io/T3/installation/)
for manual installation and troubleshooting.

---

## Usage

Create a YAML input file and run T3:

```bash
conda activate t3_env
python T3.py input.yml
```

A minimal input file looks like:

```yaml
project: my_project

t3:
  options:
    max_T3_iterations: 10
  sensitivity:
    adapter: CanteraConstantTP
    top_SA_species: 10
    top_SA_reactions: 10

rmg:
  database:
    thermo_libraries: ['primaryThermoLibrary', 'BurkeH2O2']
    kinetics_libraries: ['primaryH2O2']
  species:
    - label: H2
      smiles: '[H][H]'
      concentration: 0.67
    - label: O2
      smiles: '[O][O]'
      concentration: 0.33
    - label: OH
      smiles: '[OH]'
      SA_observable: true
  reactors:
    - type: gas batch constant T P
      T: 1000
      P: 1
      termination_conversion:
        H2: 0.9
      termination_time: [5, 's']
  model:
    core_tolerance: [0.05, 0.01]

qm:
  adapter: ARC
  level_of_theory: CBS-QB3
```

See the [input reference](https://reactionmechanismgenerator.github.io/T3/input_reference/)
for all available options.

---

## Documentation

See the [full documentation](https://reactionmechanismgenerator.github.io/T3/) for:

- [Installation](https://reactionmechanismgenerator.github.io/T3/installation/)
- [Running T3](https://reactionmechanismgenerator.github.io/T3/running/)
- [Tutorials](https://reactionmechanismgenerator.github.io/T3/tutorials/1_no_qm/)
- [Input reference](https://reactionmechanismgenerator.github.io/T3/input_reference/)
- [Examples](https://reactionmechanismgenerator.github.io/T3/examples/)

---

## Development

```bash
make test            # Run tests with coverage
make test-main       # Run main unit tests only
make test-functional # Run functional (integration) tests
```

---

## Contributing

We welcome contributions!

- Found a bug or have a feature request? [Open an issue](https://github.com/ReactionMechanismGenerator/T3/issues)
- See the [contributing guide](https://reactionmechanismgenerator.github.io/T3/contribute/) to get started

---

## License

T3 is released under the [MIT License](https://github.com/ReactionMechanismGenerator/T3/blob/main/LICENSE).
