"""
Shim module that re-implements just enough of rmgpy's API for T3 to read,
write, and merge RMG library files (thermo + kinetics) and to construct the
data objects that the simulation adapters pass to Cantera, without taking an
in-process dependency on rmgpy itself.
"""

import logging
import os
import tempfile
from dataclasses import dataclass, field
from typing import List, Tuple, Union, Optional, Any, Dict


class Quantity:
    """
    Minimal shim for RMG Quantity.
    Assumes the value passed is already in SI or the code will handle it.
    Use this if the code expects an object with .value_si
    """
    def __init__(self, value, units=None):
        self.value_si = value
        self.units = units

    def __repr__(self):
        return f"Quantity({self.value_si}, '{self.units}')"


class GenericData:
    """
    Shim for rmgpy.tools.data.GenericData
    """
    def __init__(self,
                 label: str,
                 data: Any,
                 units: Optional[str] = None,
                 species: Optional[Any] = None,
                 reaction: Optional[Any] = None,
                 index: Optional[int] = None):
        self.label = label
        self.data = data
        self.units = units
        self.species = species
        self.reaction = reaction
        self.index = index


class Reaction:
    """
    Shim for rmgpy.reaction.Reaction
    """
    def __init__(self,
                 label: str = '',
                 reactants: Optional[List[Any]] = None,
                 products: Optional[List[Any]] = None,
                 kinetics: Optional[Any] = None,
                 comment: str = '',
                 index: Optional[int] = None):
        self.label = label
        self.reactants = reactants or []
        self.products = products or []
        self.kinetics = kinetics
        self.comment = comment
        self.index = index
        self.is_pressure_dependent = False

    def __repr__(self):
        return f"<Reaction '{self.label}'>"

    def is_isomorphic(self, other_reaction) -> bool:
        """
        Check if this reaction is isomorphic to another reaction.
        """
        # Get species object lists
        other_reactants = getattr(other_reaction, 'r_species', getattr(other_reaction, 'reactants', []))
        other_products = getattr(other_reaction, 'p_species', getattr(other_reaction, 'products', []))

        # Check forward
        if self._check_species_lists(self.reactants, other_reactants) and \
           self._check_species_lists(self.products, other_products):
            return True

        # Check reverse (assuming all reversible for simplicity in shim, or check attribute?)
        # T3Reaction has arrow/reversible. Reaction shim does not strictly have arrow.
        # But let's be generous for checking duplicates.
        if self._check_species_lists(self.reactants, other_products) and \
           self._check_species_lists(self.products, other_reactants):
            return True

        return False

    def _check_species_lists(self, list1, list2) -> bool:
        if len(list1) != len(list2):
            return False
        l2 = list(list2)
        for s1 in list1:
            match = None
            for s2 in l2:
                # Use is_isomorphic if available, else equality
                if hasattr(s1, 'is_isomorphic'):
                    if s1.is_isomorphic(s2):
                        match = s2
                        break
                elif s1 == s2:
                    match = s2
                    break
            if match:
                l2.remove(match)
            else:
                return False
        return True


class PDepNetwork:
    """
    Shim for rmgpy.rmg.pdep.PDepNetwork
    """
    def __init__(self, index: Optional[int] = None):
        self.index = index


class PDepReaction(Reaction):
    """
    Shim for rmgpy.rmg.pdep.PDepReaction
    """
    def __init__(self,
                 index: Optional[int] = None,
                 reactants: Optional[List[Any]] = None,
                 products: Optional[List[Any]] = None,
                 network: Optional[PDepNetwork] = None,
                 comment: str = '',
                 **kwargs):
        super().__init__(reactants=reactants, products=products, index=index, comment=comment, **kwargs)
        self.network = network
        self.is_pressure_dependent = True


class CanteraCondition:
    """
    Shim for the condition object returned by generate_cantera_conditions.
    """
    def __init__(self, reactor_type: str, reaction_time: Quantity, mol_frac: dict,
                 T0: Optional[Quantity] = None, P0: Optional[Quantity] = None, V0: Optional[Quantity] = None):
        self.reactor_type = reactor_type
        self.reaction_time = reaction_time
        self.mol_frac = mol_frac
        self.T0 = T0
        self.P0 = P0
        self.V0 = V0


def generate_cantera_conditions(reactor_type_list: List[str],
                                reaction_time_list: tuple,
                                mol_frac_list: List[dict],
                                T0_list: Optional[tuple] = None,
                                P0_list: Optional[tuple] = None,
                                V0_list: Optional[tuple] = None,
                                ) -> List[CanteraCondition]:
    """
    Generates a list of CanteraCondition objects (shim).
    Handles the combinatorics of input lists.
    Args inputs are typically tuples of ([list_of_values], 'units').
    """
    conditions = []

    # helper to extract list and units
    def extract_values(arg_tuple):
        if arg_tuple is None:
            return [None], None
        values, units = arg_tuple
        if not isinstance(values, list):
            values = [values]
        return values, units

    r_times, r_time_units = extract_values(reaction_time_list)
    T0s, T0_units = extract_values(T0_list)
    P0s, P0_units = extract_values(P0_list)
    V0s, V0_units = extract_values(V0_list)

    # We assume reactor_type_list and mol_frac_list are lists (or single items to be broadcasted,
    # but the original code loops over them or assumes structure.
    # Looking at usage in cantera_constantTP.py:
    # reactor_type_list = [self.cantera_reactor_type] (length 1)
    # mol_frac_list = [self.initialconds] (length 1)
    # Tlist = ([...], 'K')
    # So we mainly iterate over T, P, V.

    # Iterate over cartesian product
    for reactor_type in reactor_type_list:
        for mol_frac in mol_frac_list:
            for reaction_time in r_times:
                for T0 in T0s:
                    for P0 in P0s:
                        for V0 in V0s:
                            # Create Quantity objects
                            qt_time = Quantity(reaction_time, r_time_units)
                            qt_T0 = Quantity(T0, T0_units) if T0 is not None else None
                            qt_P0 = Quantity(P0, P0_units) if P0 is not None else None
                            qt_V0 = Quantity(V0, V0_units) if V0 is not None else None

                            cond = CanteraCondition(
                                reactor_type=reactor_type,
                                reaction_time=qt_time,
                                mol_frac=mol_frac,
                                T0=qt_T0,
                                P0=qt_P0,
                                V0=qt_V0
                            )
                            conditions.append(cond)

    return conditions


# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

def fmt_val(val: Any) -> str:
    """
    Formats a value for output in an RMG library file.

    We generally use repr() to ensure valid Python syntax and preserve values
    faithfully (including strings with escapes, lists, tuples, etc.).
    """
    return repr(val)


def _quote_multiline_u(text: str) -> str:
    """
    Return a valid Python unicode string literal for potentially-multiline text.

    Goals:
      - Preserve content exactly (no added leading newline, no stripping).
      - Produce readable triple-quoted blocks when possible.
      - Avoid syntax errors if content includes triple-quote sequences.
      - Avoid syntax errors when content ends with the delimiter quote character.

    Strategy:
      1) Prefer triple double quotes: u\"\"\"...\"\"\" if:
           - content does not contain \"\"\"
           - content does not end with a double quote (")
      2) Else prefer triple single quotes: u'''...''' if:
           - content does not contain '''
           - content does not end with a single quote (')
      3) Else fall back to repr-based literal: u'...\\n...' (with escapes), preserving content
    """
    if text is None:
        text = ""

    has_triple_double = '"""' in text
    has_triple_single = "'''" in text

    ends_with_double = text.endswith('"')
    ends_with_single = text.endswith("'")

    if not has_triple_double and not ends_with_double:
        return 'u"""' + text + '"""'

    if not has_triple_single and not ends_with_single:
        return "u'''" + text + "'''"

    return "u" + repr(text)


# ------------------------------------------------------------------------------
# Thermo Classes
# ------------------------------------------------------------------------------

@dataclass
class NASAPolynomial:
    coeffs: List[float] = field(default_factory=list)
    Tmin: Optional[Tuple[float, str]] = None
    Tmax: Optional[Tuple[float, str]] = None

    def __repr__(self):
        return (
            f"NASAPolynomial(coeffs={fmt_val(self.coeffs)}, "
            f"Tmin={fmt_val(self.Tmin)}, Tmax={fmt_val(self.Tmax)})"
        )


@dataclass
class NASA:
    polynomials: List[NASAPolynomial] = field(default_factory=list)
    Tmin: Optional[Tuple[float, str]] = None
    Tmax: Optional[Tuple[float, str]] = None
    E0: Optional[Tuple[float, str]] = None
    Cp0: Optional[Tuple[float, str]] = None
    CpInf: Optional[Tuple[float, str]] = None
    comment: str = ""

    def __repr__(self):
        args = [
            f"polynomials={fmt_val(self.polynomials)}",
            f"Tmin={fmt_val(self.Tmin)}",
            f"Tmax={fmt_val(self.Tmax)}",
        ]
        if self.E0 is not None:
            args.append(f"E0={fmt_val(self.E0)}")
        if self.Cp0 is not None:
            args.append(f"Cp0={fmt_val(self.Cp0)}")
        if self.CpInf is not None:
            args.append(f"CpInf={fmt_val(self.CpInf)}")
        if self.comment:
            args.append(f"comment={self.comment!r}")
        return f"NASA({', '.join(args)})"


@dataclass
class ThermoData:
    Tdata: Optional[Tuple[List[float], str]] = None
    Cpdata: Optional[Tuple[List[float], str]] = None
    H298: Optional[Tuple[float, str]] = None
    S298: Optional[Tuple[float, str]] = None
    Tmin: Optional[Tuple[float, str]] = None
    Tmax: Optional[Tuple[float, str]] = None
    Cp0: Optional[Tuple[float, str]] = None
    CpInf: Optional[Tuple[float, str]] = None
    comment: str = ""

    def __repr__(self):
        args = [
            f"Tdata={fmt_val(self.Tdata)}",
            f"Cpdata={fmt_val(self.Cpdata)}",
            f"H298={fmt_val(self.H298)}",
            f"S298={fmt_val(self.S298)}",
        ]
        if self.Tmin is not None:
            args.append(f"Tmin={fmt_val(self.Tmin)}")
        if self.Tmax is not None:
            args.append(f"Tmax={fmt_val(self.Tmax)}")
        if self.Cp0 is not None:
            args.append(f"Cp0={fmt_val(self.Cp0)}")
        if self.comment:
            args.append(f"comment={self.comment!r}")
        if self.CpInf is not None:
            args.append(f"CpInf={fmt_val(self.CpInf)}")

        joined_args = ",\n        ".join(args)
        return f"ThermoData(\n        {joined_args},\n    )"


@dataclass
class Wilhoit:
    coeffs: List[float]
    Cp0: Tuple[float, str]
    CpInf: Tuple[float, str]
    H0: Tuple[float, str]
    S0: Tuple[float, str]
    B: Tuple[float, str]
    Tmin: Tuple[float, str]
    Tmax: Tuple[float, str]

    def __repr__(self):
        return (
            f"Wilhoit(coeffs={fmt_val(self.coeffs)}, Cp0={fmt_val(self.Cp0)}, "
            f"CpInf={fmt_val(self.CpInf)}, H0={fmt_val(self.H0)}, S0={fmt_val(self.S0)}, "
            f"B={fmt_val(self.B)}, Tmin={fmt_val(self.Tmin)}, Tmax={fmt_val(self.Tmax)})"
        )


# ------------------------------------------------------------------------------
# Kinetics Classes
# ------------------------------------------------------------------------------

@dataclass
class Arrhenius:
    A: Optional[Union[float, Tuple]] = None
    n: float = 0.0
    Ea: Optional[Union[float, Tuple]] = None
    T0: Optional[Union[float, Tuple]] = None
    Tmin: Optional[Tuple[float, str]] = None
    Tmax: Optional[Tuple[float, str]] = None
    Pmin: Optional[Tuple[float, str]] = None
    Pmax: Optional[Tuple[float, str]] = None
    comment: str = ""

    def __repr__(self):
        args = [
            f"A={fmt_val(self.A)}",
            f"n={self.n}",
            f"Ea={fmt_val(self.Ea)}",
            f"T0={fmt_val(self.T0)}",
        ]
        if self.Tmin is not None:
            args.append(f"Tmin={fmt_val(self.Tmin)}")
        if self.Tmax is not None:
            args.append(f"Tmax={fmt_val(self.Tmax)}")
        if self.Pmin is not None:
            args.append(f"Pmin={fmt_val(self.Pmin)}")
        if self.Pmax is not None:
            args.append(f"Pmax={fmt_val(self.Pmax)}")
        return f"Arrhenius({', '.join(args)})"


@dataclass
class MultiArrhenius:
    arrhenius: List[Arrhenius]

    def __repr__(self):
        arr_str = ",\n            ".join([repr(x) for x in self.arrhenius])
        return f"MultiArrhenius(\n        arrhenius=[\n            {arr_str},\n        ],\n    )"


@dataclass
class PDepArrhenius:
    pressures: Tuple[List[float], str]
    arrhenius: List[Arrhenius]

    def __repr__(self):
        arr_str = ",\n            ".join([repr(x) for x in self.arrhenius])
        return (
            f"PDepArrhenius(\n        pressures={fmt_val(self.pressures)},\n"
            f"        arrhenius=[\n            {arr_str},\n        ],\n    )"
        )


@dataclass
class Chebyshev:
    coeffs: List[List[float]]
    kunits: str
    Tmin: Tuple[float, str]
    Tmax: Tuple[float, str]
    Pmin: Tuple[float, str]
    Pmax: Tuple[float, str]
    comment: str = ""

    def __repr__(self):
        coeffs_str = "[" + ",\n            ".join([repr(row) for row in self.coeffs]) + "]"
        return (
            f"Chebyshev(\n        coeffs={coeffs_str},\n"
            f"        kunits={repr(self.kunits)}, "
            f"Tmin={fmt_val(self.Tmin)}, Tmax={fmt_val(self.Tmax)}, "
            f"Pmin={fmt_val(self.Pmin)}, Pmax={fmt_val(self.Pmax)}\n    )"
        )


@dataclass
class ThirdBody:
    arrheniusLow: Arrhenius
    efficiencies: Optional[Dict[str, float]] = None
    comment: str = ""

    def __repr__(self):
        args = [f"arrheniusLow={self.arrheniusLow}"]
        if self.efficiencies is not None:
            eff_str = "{" + ", ".join([f"{repr(k)}: {v}" for k, v in sorted(self.efficiencies.items())]) + "}"
            args.append(f"efficiencies={eff_str}")
        return f"ThirdBody({', '.join(args)})"


@dataclass
class Lindemann:
    arrheniusHigh: Arrhenius
    arrheniusLow: Arrhenius
    efficiencies: Optional[Dict[str, float]] = None
    comment: str = ""

    def __repr__(self):
        args = [
            f"arrheniusHigh={self.arrheniusHigh}",
            f"arrheniusLow={self.arrheniusLow}",
        ]
        if self.efficiencies is not None:
            eff_str = "{" + ", ".join([f"{repr(k)}: {v}" for k, v in sorted(self.efficiencies.items())]) + "}"
            args.append(f"efficiencies={eff_str}")
        return f"Lindemann(\n        {',\n        '.join(args)},\n    )"


@dataclass
class Troe:
    arrheniusHigh: Arrhenius
    arrheniusLow: Arrhenius
    alpha: float
    T3: Optional[Tuple[float, str]] = None
    T1: Optional[Tuple[float, str]] = None
    T2: Optional[Tuple[float, str]] = None
    efficiencies: Optional[Dict[str, float]] = None
    comment: str = ""

    def __repr__(self):
        args = [
            f"arrheniusHigh={self.arrheniusHigh}",
            f"arrheniusLow={self.arrheniusLow}",
            f"alpha={self.alpha}",
        ]
        if self.T3 is not None:
            args.append(f"T3={fmt_val(self.T3)}")
        if self.T1 is not None:
            args.append(f"T1={fmt_val(self.T1)}")
        if self.T2 is not None:
            args.append(f"T2={fmt_val(self.T2)}")
        if self.efficiencies is not None:
            eff_str = "{" + ", ".join([f"{repr(k)}: {v}" for k, v in sorted(self.efficiencies.items())]) + "}"
            args.append(f"efficiencies={eff_str}")
        return f"Troe(\n        {',\n        '.join(args)},\n    )"


@dataclass
class KineticsData:
    Tdata: Tuple[List[float], str]
    kdata: Tuple[List[float], str]
    Tmin: Optional[Tuple[float, str]] = None
    Tmax: Optional[Tuple[float, str]] = None

    def __repr__(self):
        return (
            f"KineticsData(Tdata={fmt_val(self.Tdata)}, kdata={fmt_val(self.kdata)}, "
            f"Tmin={fmt_val(self.Tmin)}, Tmax={fmt_val(self.Tmax)})"
        )


# ==============================================================================
# LIBRARY STRUCTURE
# ==============================================================================

@dataclass
class Entry:
    index: int
    label: str
    molecule: Optional[Union[str, List[str]]] = None
    thermo: Any = None
    kinetics: Any = None
    shortDesc: str = ""
    longDesc: str = ""
    degeneracy: float = 1.0
    duplicate: bool = False
    reversible: bool = True

    def __repr__(self):
        s = "entry(\n"
        s += f"    index = {self.index},\n"
        s += f"    label = {self.label!r},\n"

        if self.molecule is not None:
            if isinstance(self.molecule, list):
                s += f"    molecule = {self.molecule!r},\n"
            else:
                mol_str = self.molecule.rstrip()
                s += f"    molecule = {_quote_multiline_u(mol_str)},\n"

        if self.thermo is not None:
            t_str = repr(self.thermo).replace("\n", "\n    ")
            s += f"    thermo = {t_str},\n"

        if self.kinetics is not None:
            k_str = repr(self.kinetics).replace("\n", "\n    ")
            s += f"    kinetics = {k_str},\n"

        if self.degeneracy != 1.0:
            s += f"    degeneracy = {self.degeneracy},\n"

        if self.duplicate:
            s += f"    duplicate = {self.duplicate},\n"

        if not self.reversible:
            s += f"    reversible = {self.reversible},\n"

        s += f"    shortDesc = u{self.shortDesc!r},\n"

        long_d = self.longDesc.rstrip()
        if not long_d:
            s += '    longDesc = u"""\n""",\n'
        else:
            s += f"    longDesc = {_quote_multiline_u(long_d)},\n"

        s += ")"
        return s


@dataclass
class Library:
    name: str = ""
    shortDesc: str = ""
    longDesc: str = ""
    entries: List[Entry] = field(default_factory=list)


# ==============================================================================
# PARSER LOGIC
# ==============================================================================

def parse_rmg_library(content: str) -> Library:
    """
    Parses a Python file content (thermo or kinetics) into a shim Library object.
    Uses exec() with a restricted context mapped to shim classes.
    """
    lib = Library()

    def entry_collector(**kwargs):
        valid_keys = Entry.__annotations__.keys()
        filtered = {k: v for k, v in kwargs.items() if k in valid_keys}

        if "index" not in filtered or "label" not in filtered:
            raise ValueError(
                f"Entry missing required fields 'index' or 'label'. Found: {list(filtered.keys())}"
            )

        ent = Entry(**filtered)
        lib.entries.append(ent)

    safe_builtins = {
        "True": True,
        "False": False,
        "None": None,
        "list": list,
        "tuple": tuple,
        "dict": dict,
        "str": str,
        "int": int,
        "float": float,
        # Common low-risk helpers some library files might use:
        "set": set,
        "len": len,
        "range": range,
        "min": min,
        "max": max,
    }

    local_context = {
        "entry": entry_collector,

        # Thermo
        "NASA": NASA,
        "NASAPolynomial": NASAPolynomial,
        "ThermoData": ThermoData,
        "Wilhoit": Wilhoit,

        # Kinetics
        "Arrhenius": Arrhenius,
        "PDepArrhenius": PDepArrhenius,
        "MultiArrhenius": MultiArrhenius,
        "Chebyshev": Chebyshev,
        "ThirdBody": ThirdBody,
        "Lindemann": Lindemann,
        "Troe": Troe,
        "KineticsData": KineticsData,
    }

    try:
        exec(content, {"__builtins__": safe_builtins}, local_context)
    except Exception as e:
        raise ValueError(f"Error parsing RMG library content: {e}") from e

    lib.name = local_context.get("name", "")
    lib.shortDesc = local_context.get("shortDesc", "")
    lib.longDesc = local_context.get("longDesc", "")

    return lib


def write_library_to_string(lib: Library) -> str:
    """Writes the shim Library object back to a valid RMG-style python string."""
    output: List[str] = []
    output.append("#!/usr/bin/env python")
    output.append("# encoding: utf-8")
    output.append("")
    output.append(f"name = {lib.name!r}")
    output.append(f"shortDesc = u{lib.shortDesc!r}")

    if lib.longDesc:
        safe_desc = _quote_multiline_u(lib.longDesc.rstrip())
        output.append(f"longDesc = {safe_desc}")
    else:
        output.append('longDesc = u"""\n"""')

    output.append("")

    for ent in lib.entries:
        output.append(repr(ent))
        output.append("")

    return "\n".join(output)


def write_atomic(path: str, content: str) -> None:
    """
    Atomically write content to a file.
    
    Writes to a temporary file first, then renames to the target path.
    This ensures the target file is never in a partially-written state.
    
    Args:
        path: Target file path.
        content: Content to write.
    """
    dir_name = os.path.dirname(path) or '.'
    fd, tmp_path = tempfile.mkstemp(dir=dir_name, prefix='.tmp_atomic_')
    try:
        with os.fdopen(fd, 'w') as f:
            f.write(content)
        os.replace(tmp_path, path)
    except Exception:
        try:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
        except OSError:
            logging.getLogger(__name__).warning(f'Failed to clean up temp file: {tmp_path}')
        raise
