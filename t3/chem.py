"""
t3.chem for representing chemical species and reactions in T3.

Defines T3-specific data structures that extend ARC/RMG objects.
This isolates T3's workflow state (indices, sources, status) from the
underlying chemical definitions found in ARC.
"""

from __future__ import annotations
from enum import Enum
from typing import Optional, Dict, Any, Union, List

from arc.species.species import ARCSpecies
from arc.reaction.reaction import ARCReaction


class T3Status(str, Enum):
    """Standardized status flags for T3 workflow steps."""
    PENDING = "pending"
    RUNNING = "running"
    CONVERGED = "converged"
    FAILED = "failed"
    SKIPPED = "skipped"


class ThermoMethod(str, Enum):
    """Standardized methods for obtaining thermodynamic properties."""
    QM = "QM"
    LIBRARY = "Library"
    GAV = "GAV"
    ML = "ML"
    USER = "User"
    UNKNOWN = "Unknown"


class KineticsMethod(str, Enum):
    """Standardized methods for obtaining reaction kinetics."""
    QM = "QM"
    LIBRARY = "Library"
    RATE_RULES = "Rate Rules"
    TRAINING_SET = "Training Set"
    PDEP = "PDep"
    USER = "User"
    UNKNOWN = "Unknown"


class T3Species(ARCSpecies):
    """
    A wrapper for ARCSpecies that tracks T3 workflow metadata.

    Attributes:
        thermo_method (ThermoMethod): The high-level category of the thermo source.
        thermo_source (str): The specific details of the source.
                             If a non-standard method string is provided at init,
                             it is prepended here and method becomes UNKNOWN.
        thermo_comment (str): A verbose, formatted string for the final RMG library output.
        t3_status (T3Status): Current workflow status.
        t3_index (int): Permanent unique ID within this T3 run (0, 1, 2...).
        rmg_index (int): Transient ID from the most recent RMG model generation.
        created_at_iteration (int): The T3 iteration number where this species first appeared.
        reasons (List[str]): Reasons for calculating this species.
        qm_label (str): The species label used for QM calculations.
    """

    def __init__(self,
                 label: Optional[str] = None,
                 qm_label: Optional[str] = None,
                 # Source Tracking
                 thermo_method: Optional[Union[ThermoMethod, str]] = None,
                 thermo_source: Optional[str] = None,
                 thermo_comment: str = "",
                 # T3 State
                 t3_status: Union[T3Status, str] = T3Status.PENDING,
                 t3_index: Optional[int] = None,
                 rmg_index: Optional[int] = None,
                 created_at_iteration: int = 0,
                 reasons: Optional[Union[List[str], str]] = None,
                 # ARCSpecies arguments
                 *args,
                 **kwargs):

        super().__init__(label=label, *args, **kwargs)

        self.label = label
        self.qm_label = qm_label

        # Thermo Method Normalization
        if thermo_method is None:
            self.thermo_method: Optional[ThermoMethod] = None
            self.thermo_source = thermo_source
        elif isinstance(thermo_method, ThermoMethod):
            self.thermo_method = thermo_method
            self.thermo_source = thermo_source
        else:
            # Try to cast string to Enum
            try:
                self.thermo_method = ThermoMethod(thermo_method)
                self.thermo_source = thermo_source
            except ValueError:
                # Fallback: unknown string -> UNKNOWN enum, preserve raw string in source
                self.thermo_method = ThermoMethod.UNKNOWN
                self.thermo_source = f"[{thermo_method}] {thermo_source}" if thermo_source else thermo_method

        try:
            self.t3_status = T3Status(t3_status)
        except ValueError:
            self.t3_status = T3Status.PENDING

        self.thermo_comment = thermo_comment
        self.t3_index = t3_index
        self.rmg_index = rmg_index
        self.created_at_iteration = created_at_iteration
        self.reasons = [reasons] if isinstance(reasons, str) else reasons or []

    @property
    def is_converged(self) -> bool:
        """Helper to check if the species is effectively 'done'."""
        return self.t3_status == T3Status.CONVERGED

    def as_dict(self, reset_atom_ids: bool = False) -> dict:
        """
        Extended dictionary representation including T3 metadata.
        Crucial for saving T3 state to YAML/JSON between restarts.
        """
        data = super().as_dict(reset_atom_ids=reset_atom_ids)
        data.update({
            "thermo_method": self.thermo_method,
            "thermo_source": self.thermo_source,
            "thermo_comment": self.thermo_comment,
            "t3_status": self.t3_status,
            "t3_index": self.t3_index,
            "rmg_index": self.rmg_index,
            "created_at_iteration": self.created_at_iteration,
        })
        return data

    @classmethod
    def from_dict(cls, species_dict: Dict[str, Any]) -> T3Species:
        """
        Reconstruct a T3Species from a dictionary.
        Handles the separation of T3-specific args vs ARC args.
        """
        # Extract T3 specific keys
        t3_keys = {
            "thermo_method", "thermo_source", "thermo_comment",
            "t3_status", "t3_index", "rmg_index", "created_at_iteration", "reasons",
            "qm_label",
        }
        t3_kwargs = {k: species_dict.pop(k) for k in t3_keys if k in species_dict}
        
        # Everything else goes to ARCSpecies via kwargs
        return cls(**t3_kwargs, **species_dict)

    def __repr__(self) -> str:
        """
        Readable representation for debugging logs.
        Example: <T3Species 'CH4' (index: 1) [QM/CBS-QB3] status: converged>
        """
        if self.thermo_method:
            method_str = f" [{self.thermo_method.value}"
            if self.thermo_source:
                method_str += f"/{self.thermo_source}"
            method_str += "]"
        elif self.thermo_source:
            method_str = f" [Source: {self.thermo_source}]"
        else:
            method_str = ""

        index_str = f" (index: {self.t3_index})" if self.t3_index is not None else ""

        status_val = self.t3_status.value if isinstance(self.t3_status, T3Status) else self.t3_status
        
        return f"<T3Species '{self.label}'{index_str}{method_str} status: {status_val}>"


class T3Reaction(ARCReaction):
    """
    A wrapper for ARCReaction that tracks T3 workflow metadata.

    Attributes:
        kinetics_method (KineticsMethod): The high-level category of the kinetics source.
        kinetics_source (str): The specific details of the source.
                               If a non-standard method string is provided, it is prepended here.
        kinetics_comment (str): A verbose, formatted string for the final RMG library output.
        t3_status (T3Status): Current workflow status.
        t3_index (int): Permanent unique ID within this T3 run (if tracked persistently).
        rmg_index (int): Transient ID from the most recent RMG model generation (Reaction #).
        created_at_iteration (int): The T3 iteration number where this reaction first appeared.
        reasons (List[str]): Reasons for calculating this reaction.
        qm_label (str): The reaction label used for QM calculations.
        reactant_keys (List[int]): T3 species indices of the reactants.
        product_keys (List[int]): T3 species indices of the products.
        is_pressure_dependent (bool): Whether the reaction is pressure-dependent.
    """
    
    def __init__(self,
                 # Source Tracking
                 kinetics_method: Optional[Union[KineticsMethod, str]] = None,
                 kinetics_source: Optional[str] = None,
                 kinetics_comment: str = "",
                 # T3 State
                 t3_status: Union[T3Status, str] = T3Status.PENDING,
                 t3_index: Optional[int] = None,
                 rmg_index: Optional[int] = None,
                 created_at_iteration: int = 0,
                 reasons: Optional[Union[List[str], str]] = None,
                 qm_label: Optional[str] = None,
                 rmg_label: Optional[str] = None,
                 reactant_keys: Optional[List[int]] = None,
                 product_keys: Optional[List[int]] = None,
                 is_pressure_dependent: Optional[bool] = None,
                 # ARCReaction arguments
                 *args, 
                 **kwargs):

        super().__init__(*args, **kwargs)

        self.qm_label = qm_label
        self.rmg_label = rmg_label
        self.reactant_keys = reactant_keys
        self.product_keys = product_keys
        self.is_pressure_dependent = is_pressure_dependent

        # Kinetics Method Normalization
        if kinetics_method is None:
            self.kinetics_method: Optional[KineticsMethod] = None
            self.kinetics_source = kinetics_source
        elif isinstance(kinetics_method, KineticsMethod):
            self.kinetics_method = kinetics_method
            self.kinetics_source = kinetics_source
        else:
            try:
                self.kinetics_method = KineticsMethod(kinetics_method)
                self.kinetics_source = kinetics_source
            except ValueError:
                self.kinetics_method = KineticsMethod.UNKNOWN
                prefix = f"[{kinetics_method}] "
                self.kinetics_source = f"{prefix}{kinetics_source}" if kinetics_source else kinetics_method

        try:
            self.t3_status = T3Status(t3_status)
        except ValueError:
            self.t3_status = T3Status.PENDING

        self.kinetics_comment = kinetics_comment
        self.t3_index = t3_index
        self.rmg_index = rmg_index
        self.created_at_iteration = created_at_iteration
        self.reasons = [reasons] if isinstance(reasons, str) else reasons or []

    @property
    def is_converged(self) -> bool:
        """Helper to check if the reaction is effectively 'done'."""
        return self.t3_status == T3Status.CONVERGED

    def as_dict(self,
                reset_atom_ids: bool = False,
                report_family: bool = True,
                ) -> dict:
        """
        Extended dictionary representation including T3 metadata.
        """
        data = super().as_dict(reset_atom_ids=reset_atom_ids, report_family=report_family)
        data.update({
            "kinetics_method": self.kinetics_method,
            "kinetics_source": self.kinetics_source,
            "kinetics_comment": self.kinetics_comment,
            "t3_status": self.t3_status,
            "t3_index": self.t3_index,
            "rmg_index": self.rmg_index,
            "created_at_iteration": self.created_at_iteration,
            "reasons": self.reasons,
            "qm_label": self.qm_label,
            "rmg_label": self.rmg_label,
            "reactant_keys": self.reactant_keys,
            "product_keys": self.product_keys,
        })
        return data

    @classmethod
    def from_dict(cls,
                  reaction_dict: dict,
                  species_list: Optional[list] = None) -> 'T3Reaction':
        """
        Reconstruct a T3Reaction from a dictionary.

        This separates T3-specific metadata (passed to __init__) from
        standard ARC attributes (loaded via ARCReaction.from_dict).
        """
        # 1. Define T3-specific keys to extract for __init__
        t3_keys = {
            "kinetics_method", "kinetics_source", "kinetics_comment",
            "t3_status", "t3_index", "rmg_index", "created_at_iteration", "reasons",
            "qm_label", "rmg_label", "reactant_keys", "product_keys",
            "rmg_reaction_family",
        }

        # 2. Extract T3 kwargs (removing them from dict so ARC doesn't see them)
        t3_kwargs = {}
        for k in t3_keys:
            if k in reaction_dict:
                t3_kwargs[k] = reaction_dict.pop(k)

        # 3. Instantiate the T3Reaction
        # We pass T3 args here. ARC args are left for the loader below.
        obj = cls(**t3_kwargs, **reaction_dict)

        # 4. Populate standard ARC attributes
        # We explicitly call the parent class's instance method on the new object.
        # This preserves the logic for linking 'species_list' to reactants/products.
        ARCReaction.from_dict(obj, reaction_dict, species_list=species_list)

        return obj

    def __repr__(self) -> str:
        """
        Readable representation for debugging logs.
        Example: <T3Reaction 'H + CH4 <=> H2 + CH3' (index: 5) [Rate Rules/H_Abstraction] status: pending>
        """
        if self.kinetics_method:
            method_str = f" [{self.kinetics_method.value}"
            if self.kinetics_source:
                method_str += f"/{self.kinetics_source}"
            method_str += "]"
        elif self.kinetics_source:
            method_str = f" [Source: {self.kinetics_source}]"
        else:
            method_str = ""
            
        index_str = f" (index: {self.t3_index})" if self.t3_index is not None else ""
        status_val = self.t3_status.value if isinstance(self.t3_status, T3Status) else self.t3_status

        return f"<T3Reaction '{self.label}'{index_str}{method_str} status: {status_val}>"

    def get_reaction_smiles_label(self) -> str:
        """
        Get the reaction SMILES label.
        """
        reactants, products = self.get_reactants_and_products(return_copies=True)
        smiles_r = [reactant.mol.copy(deep=True).to_smiles() for reactant in reactants]
        smiles_p = [product.mol.copy(deep=True).to_smiles() for product in products]
        if not any(smiles_r) or not any(smiles_p):
            raise ValueError(f"""Could not find smiles for one or more species
                                 got: reactants: {smiles_r}
                                      products: {smiles_p}""")
        return "+".join(smiles_r)+"<=>"+".".join(smiles_p)
