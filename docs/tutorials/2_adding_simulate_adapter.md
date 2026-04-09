# 2. Adding a Cantera simulation adapter

This tutorial walks through creating a new Cantera-based simulation adapter for T3.
We will implement a **constant-volume, isothermal** batch reactor (`CanteraConstantTV`),
which is not currently part of T3 but is straightforward to add using the existing
`CanteraBase` class.

All Cantera adapters in T3 share a common base class that handles mechanism loading,
condition generation, time integration, sensitivity analysis, and data storage.
A new adapter only needs to:

1. Set the `cantera_reactor_type` class attribute.
2. Implement `create_reactor()` to return the appropriate Cantera reactor.
3. Optionally override `get_idt_by_T()` if the default behavior is not appropriate.


## Step 1: Create the adapter file

Create a new file at `T3/t3/simulate/cantera_constant_tv.py`:

```python
"""
Cantera Simulator Adapter module
Used to run mechanism analysis with Cantera as an ideal gas in a batch reactor at constant T-V
"""

import cantera as ct

from t3.simulate.cantera_base import CanteraBase
from t3.simulate.factory import register_simulate_adapter
```

These imports are the minimum required for any Cantera adapter:

- `cantera` provides the reactor objects.
- `CanteraBase` provides all shared simulation logic.
- `register_simulate_adapter` registers the adapter so T3 can find it by name.


## Step 2: Define the adapter class

```python
class CanteraConstantTV(CanteraBase):
    """
    Simulates ideal gases in a batch reactor at constant temperature and volume.
    Uses ``ct.IdealGasReactor`` with ``energy='off'`` to maintain constant T,
    and the default constant-volume behavior of ``IdealGasReactor``.
    """

    cantera_reactor_type = 'IdealGasReactor'

    def create_reactor(self):
        """Create a constant-volume reactor with the energy equation disabled (isothermal)."""
        return ct.IdealGasReactor(self.model, energy='off')

    def get_idt_by_T(self):
        """
        Since this adapter simulates at constant T, there is no ignition.

        Returns:
            idt_dict (dict): Dictionary whose values are empty lists.
        """
        return {'idt': list(), 'idt_index': list()}
```

Key points:

- **`cantera_reactor_type`**: A string identifying the Cantera reactor class. Used
  internally for logging and identification.
- **`create_reactor()`**: The only required method. It returns a Cantera reactor instance.
  The `self.model` attribute (a `ct.Solution` object) is set up by `CanteraBase` during
  initialization with the mechanism and initial conditions.
- **`get_idt_by_T()`**: The default implementation in `CanteraBase` computes ignition
  delay from the maximum `dT/dt`. Since a constant-T reactor has no temperature rise,
  we override it to return empty lists.


## Step 3: Register the adapter

At the bottom of the file, register the adapter with the factory:

```python
register_simulate_adapter("CanteraConstantTV", CanteraConstantTV)
```

The first argument is the string name that users will use in their input files.


## Complete file

Here is the complete `cantera_constant_tv.py`:

```python
"""
Cantera Simulator Adapter module
Used to run mechanism analysis with Cantera as an ideal gas in a batch reactor at constant T-V
"""

import cantera as ct

from t3.simulate.cantera_base import CanteraBase
from t3.simulate.factory import register_simulate_adapter


class CanteraConstantTV(CanteraBase):
    """
    Simulates ideal gases in a batch reactor at constant temperature and volume.
    Uses ``ct.IdealGasReactor`` with ``energy='off'`` to maintain constant T,
    and the default constant-volume behavior of ``IdealGasReactor``.
    """

    cantera_reactor_type = 'IdealGasReactor'

    def create_reactor(self):
        """Create a constant-volume reactor with the energy equation disabled (isothermal)."""
        return ct.IdealGasReactor(self.model, energy='off')

    def get_idt_by_T(self):
        """
        Since this adapter simulates at constant T, there is no ignition.

        Returns:
            idt_dict (dict): Dictionary whose values are empty lists.
        """
        return {'idt': list(), 'idt_index': list()}


register_simulate_adapter("CanteraConstantTV", CanteraConstantTV)
```


## Step 4: Register the import

Add the import to `T3/t3/simulate/__init__.py`:

```python
from t3.simulate.cantera_constant_tv import CanteraConstantTV
```

This ensures the adapter is registered when T3 starts.


## Step 5: Use the adapter

The new adapter can now be used in T3 input files:

```yaml
t3:
  sensitivity:
    adapter: CanteraConstantTV
    SA_threshold: 0.01
```


## Comparison of existing Cantera adapters

For reference, here is how the existing adapters differ:

| Adapter | Reactor Class | `energy` | Behavior |
|---|---|---|---|
| `CanteraConstantTP` | `IdealGasConstPressureReactor` | `'off'` | Constant T, constant P |
| `CanteraConstantHP` | `IdealGasConstPressureReactor` | (default: on) | Adiabatic, constant P |
| `CanteraConstantUV` | `IdealGasReactor` | (default: on) | Adiabatic, constant V |
| `CanteraPFR` | `IdealGasConstPressureReactor` | `'off'` | Isothermal PFR (Lagrangian) |
| `CanteraPFRTProfile` | `IdealGasConstPressureReactor` | `'off'` | PFR with prescribed `T(z)` (Lagrangian, no SA) |
| `CanteraJSR` | `IdealGasReactor` (in flow network) | `'off'` | Isothermal, isobaric jet stirred reactor (steady-state, open) |
| **CanteraConstantTV** | `IdealGasReactor` | `'off'` | **Constant T, constant V** |

The pattern is consistent: choose the Cantera reactor class that matches your
volume/pressure constraint, then set `energy='off'` for isothermal operation
or leave it on for adiabatic.


## Adding tests

It is recommended to add tests for the new adapter under `T3/tests/test_simulate/`.
Follow the pattern in `test_cantera_constant_tp.py`:

- Test that the adapter is registered in the factory.
- Test `set_up()` with a sample mechanism.
- Test `simulate()` produces reasonable species profiles.
- Test `get_sa_coefficients()` returns the expected SA dictionary format.
