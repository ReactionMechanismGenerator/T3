"""
t3.chem for representing chemical species and reactions in T3.

Defines T3-specific data structures that extend ARC/RMG objects.
This isolates T3's workflow state (indices, sources, status) from the
underlying chemical definitions found in ARC.
"""

from __future__ import annotations
from enum import Enum
from typing import Optional, Dict, Any, Union, List

from arc.species.species import ARCSpecies, rmg_mol_from_dict_repr
from arc.reaction.reaction import ARCReaction

from t3.common import to_chemkin_label


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

    Args:
        label (str): The species label as parsed from RMG.
        key (int, optional): A unique integer identifier for the species. If None, it will be auto-assigned.
        species_dict (dict, optional): A dictionary representation of the species, used for reconstruction.
        thermo_method (ThermoMethod or str, optional): The high-level category of the thermo source. Can be a
                                                       ThermoMethod enum or a string that will be cast to the enum.
        thermo_source (str, optional): The specific details of the source.
        thermo_comment (str, optional): The raw RMG thermo comment.
        t3_status (T3Status or str, optional): Current workflow status. Can be a T3Status enum or a string that will be
                                               cast to the enum if possible. Defaults to PENDING.
        created_at_iteration (int, optional): The T3 iteration number where this species first appeared. Default: 0.
        reasons (List[str] or str, optional): Reasons for calculating this species. Can  a list of strings or a single
                                              string that will be converted into a list.
        thermo (Any, optional): The thermodynamic data object.

    Attributes:
        key (int): The unique integer identifier for the species within the T3 workflow.
        qm_label (str): The standardized label used for QM calculations, formatted as 's{key}_{formula}'.
        rmg_label (Dict[int, str]): A history dictionary mapping iteration numbers to the species label
                                    assigned by RMG (e.g., {0: "CH4(10)"}).
        formula (str): The chemical formula of the species.
        thermo_method (ThermoMethod): The high-level category of the thermodynamic data source.
        thermo_source (str): The specific details of the thermodynamic data source.
                             If a non-standard method string was provided, it is prepended here.
        thermo_comment (str): The raw RMG thermo comment or metadata.
        t3_status (T3Status): The current status of the species in the T3 workflow (e.g., CONVERGED, FAILED).
        created_at_iteration (int): The T3 iteration number where this species was first identified.
        reasons (List[str]): A list of strings explaining why this species was selected for calculation.
        thermo (Any): The thermodynamic data object (e.g., an RMG ThermoData object).
    """

    _index_counter = 0

    # Fields that belong to T3Species but not to ARCSpecies.
    # They are silently consumed here so they never leak into super().__init__().
    _T3_ONLY_KWARGS = frozenset({
        'qm_label', 't3_index', 'rmg_label', 'formula',
    })

    def __init__(self,
                 label: Optional[str] = None,
                 key: Optional[int] = None,
                 species_dict: Optional[dict] = None,
                 thermo_method: Optional[Union[ThermoMethod, str]] = None,
                 thermo_source: Optional[str] = None,
                 thermo_comment: str = "",
                 t3_status: Union[T3Status, str] = T3Status.PENDING,
                 created_at_iteration: int = 0,
                 reasons: Optional[Union[List[str], str]] = None,
                 thermo: Optional[Any] = None,
                 *args,
                 **kwargs):

        # If species_dict is given, preprocess it the same way from_dict() does,
        # then merge its ARC-relevant keys into kwargs for super().__init__().
        _t3_dict = None
        if species_dict is not None:
            species_dict = species_dict.copy()
            _t3_keys = {'key', 'rmg_label', 'thermo_method', 'thermo_source',
                         'thermo_comment', 't3_status', 'created_at_iteration', 'reasons', 'qm_label'}
            _t3_dict = {k: species_dict.pop(k) for k in list(species_dict) if k in _t3_keys}
            is_ts = species_dict.get('is_ts', False)
            if 'mol' in species_dict and isinstance(species_dict['mol'], dict):
                species_dict['mol'] = rmg_mol_from_dict_repr(species_dict['mol'], is_ts=is_ts)
            species_dict['xyz'] = species_dict.get('final_xyz', None) or species_dict.get('initial_xyz', None)
            species_dict = remove_bad_arc_keys(species_dict)
            if label is None:
                label = species_dict.pop('label', None)
            else:
                species_dict.pop('label', None)
            if key is None:
                key = _t3_dict.pop('key', None)
            else:
                _t3_dict.pop('key', None)
            kwargs.update(species_dict)

        # Strip T3-only fields so ARCSpecies.__init__ doesn't choke on them.
        t3_extras = {k: kwargs.pop(k) for k in list(kwargs) if k in self._T3_ONLY_KWARGS}

        # ARCSpecies requires a non-None string label.  When the caller only
        # supplies an adjlist or SMILES, derive a temporary label so ARC's
        # check_label() doesn't raise.
        if label is None:
            if 'adjlist' in kwargs:
                first_line = kwargs['adjlist'].strip().splitlines()[0]
                if not first_line[0].isdigit():
                    label = first_line.split()[0]
            if label is None:
                label = kwargs.get('smiles') or 'unknown'

        super().__init__(label=label, *args, **kwargs)

        if key is None:
            self.key = T3Species._index_counter
            T3Species._index_counter += 1
        else:
            self.key = key
            if key >= T3Species._index_counter:
                T3Species._index_counter = key + 1

        self.created_at_iteration = created_at_iteration
        self.rmg_label : Dict[int, str] = t3_extras.get('rmg_label') or {self.created_at_iteration: label}
        self.formula = self.mol.get_formula()
        _dict_qm_label = _t3_dict.pop('qm_label', None) if _t3_dict else None
        self.qm_label = _dict_qm_label or t3_extras.get('qm_label') or f"s{self.key}_{self.formula}"
        self.thermo = thermo

        if thermo_method is None:
            self.thermo_method: Optional[ThermoMethod] = None
            self.thermo_source = thermo_source
        elif isinstance(thermo_method, ThermoMethod):
            self.thermo_method = thermo_method
            self.thermo_source = thermo_source
        else:
            try:
                self.thermo_method = ThermoMethod(thermo_method)
                self.thermo_source = thermo_source
            except ValueError:
                self.thermo_method = ThermoMethod.UNKNOWN
                self.thermo_source = f"[{thermo_method}] {thermo_source}" if thermo_source else thermo_method

        try:
            self.t3_status = T3Status(t3_status)
        except ValueError:
            try:
                self.t3_status = T3Status(t3_status.lower())
            except (ValueError, AttributeError):
                self.t3_status = T3Status.PENDING

        self.thermo_comment = thermo_comment
        self.reasons = [reasons] if isinstance(reasons, str) else reasons or []

        if _t3_dict is not None:
            # Apply T3-specific attributes from the dictionary, overriding defaults.
            if 'rmg_label' in _t3_dict and _t3_dict['rmg_label']:
                self.rmg_label = {int(k): v for k, v in _t3_dict['rmg_label'].items()}
            for attr in ('thermo_source', 'thermo_comment', 'created_at_iteration', 'reasons'):
                if attr in _t3_dict:
                    setattr(self, attr, _t3_dict[attr])
            if 'thermo_method' in _t3_dict and _t3_dict['thermo_method'] is not None:
                val = _t3_dict['thermo_method']
                self.thermo_method = val if isinstance(val, ThermoMethod) else ThermoMethod(val)
            if 't3_status' in _t3_dict and _t3_dict['t3_status'] is not None:
                val = _t3_dict['t3_status']
                self.t3_status = val if isinstance(val, T3Status) else T3Status(val)

    def __repr__(self) -> str:
        """
        Return a readable representation of the T3Species object for debugging.
        Example: <T3Species 's5_CH4' [QM] status: Converged>
        """
        label = getattr(self, 'qm_label', self.label)
        method_str = ""
        if self.thermo_method:
            val = self.thermo_method.value if hasattr(self.thermo_method, 'value') else str(self.thermo_method)
            method_str = f" [{val}]"
        status_val = "Unknown"
        if self.t3_status:
            status_val = self.t3_status.value if hasattr(self.t3_status, 'value') else str(self.t3_status)
        return f"<T3Species '{label}'{method_str} status: {status_val}>"

    @classmethod
    def reset_counter(cls):
        """Resets the global index counter."""
        cls._index_counter = 0

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
        thermo_method_val = self.thermo_method.value if self.thermo_method else None
        t3_status_val = self.t3_status.value if isinstance(self.t3_status, T3Status) else self.t3_status
        data.update({
            "key": self.key,
            "qm_label": self.qm_label,
            "rmg_label": self.rmg_label,
            "thermo_method": thermo_method_val,
            "thermo_source": self.thermo_source,
            "thermo_comment": self.thermo_comment,
            "t3_status": t3_status_val,
            "created_at_iteration": self.created_at_iteration,
            "reasons": self.reasons,
        })
        return data

    @classmethod
    def from_dict(cls, species_dict: Dict[str, Any]) -> 'T3Species':
        """
        Reconstruct a T3Species from a dictionary.
        """
        t3_keys = {
            "key",
            "rmg_label",
            "thermo_method",
            "thermo_source",
            "thermo_comment",
            "t3_status",
            "created_at_iteration",
            "reasons",
            "qm_label",
        }
        t3_kwargs = {k: species_dict.pop(k) for k in t3_keys if k in species_dict}
        rmg_label_history = t3_kwargs.pop("rmg_label", None)
        qm_label = t3_kwargs.pop("qm_label", None)
        is_ts = species_dict.get('is_ts', False)
        if 'mol' in species_dict and isinstance(species_dict['mol'], dict):
            species_dict['mol'] = rmg_mol_from_dict_repr(species_dict['mol'], is_ts=is_ts)
        species_dict['xyz'] = species_dict.get('final_xyz', None) or species_dict.get('initial_xyz', None)
        species_dict = remove_bad_arc_keys(species_dict)
        instance = cls(qm_label=qm_label, **t3_kwargs, **species_dict)
        if rmg_label_history:
            instance.rmg_label = {int(k): v for k, v in rmg_label_history.items()}
        return instance

    def copy(self):
        """
        Get a copy of this object instance.

        Returns:
            ARCSpecies: A copy of this object instance.
        """
        species_dict = self.as_dict(reset_atom_ids=True)
        return self.__class__.from_dict(species_dict)

    def to_chemkin(self) -> str:
        """
        Return a Chemkin-compliant label for the species.
        """
        return to_chemkin_label(self)

    def to_species_dictionary_entry(self) -> str:
        """
        Formats the species for the dictionary.txt file.
        Output format:
        s14_C2H4 (S(842))
        1 C u0 p0 c0 {2,D} {3,S} {4,S} ...
        """
        history_str = ", ".join([f"Iter {k}: {v}" for k, v in self.rmg_label.items()])
        header = f"{self.qm_label} ({history_str})"
        return f"{header}\n{self.adjlist}"


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

    # Fields that belong to T3Reaction but not to ARCReaction.
    _T3_ONLY_KWARGS = frozenset({'index'})

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
                 network: Optional[Any] = None,
                 # ARCReaction arguments
                 *args,
                 **kwargs):

        # Strip T3-only fields so ARCReaction.__init__ doesn't choke on them.
        t3_extras = {k: kwargs.pop(k) for k in list(kwargs) if k in self._T3_ONLY_KWARGS}

        super().__init__(*args, **kwargs)

        # Restore T3-only fields as attributes.
        self.index = t3_extras.get('index')

        self.qm_label = qm_label
        self.rmg_label = rmg_label
        self.reactant_keys = reactant_keys
        self.product_keys = product_keys
        self.is_pressure_dependent = is_pressure_dependent
        self.network = network

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

    def to_chemkin(self) -> str:
        """
        Return a Chemkin-compliant label for the reaction.
        """
        arrow = ' <=> ' if self.is_pressure_dependent is not False else ' = '
        if self.r_species and self.p_species:
            reactants_label = ' + '.join([r.to_chemkin() if hasattr(r, 'to_chemkin') else str(r) for r in self.r_species])
            products_label = ' + '.join([p.to_chemkin() if hasattr(p, 'to_chemkin') else str(p) for p in self.p_species])
            return f'{reactants_label}{arrow}{products_label}'
        return self.label

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
        Serializes Enums to string values to ensure clean YAML output.
        """
        data = super().as_dict(reset_atom_ids=reset_atom_ids, report_family=report_family)
        kinetics_method_val = self.kinetics_method.value if self.kinetics_method else None
        t3_status_val = self.t3_status.value if isinstance(self.t3_status, T3Status) else self.t3_status
        data.update({
            "kinetics_method": kinetics_method_val,
            "kinetics_source": self.kinetics_source,
            "kinetics_comment": self.kinetics_comment,
            "t3_status": t3_status_val,
            "t3_index": self.t3_index,
            "rmg_index": self.rmg_index,
            "created_at_iteration": self.created_at_iteration,
            "reasons": self.reasons,
            "qm_label": self.qm_label,
            "rmg_label": self.rmg_label,
            "reactant_keys": self.reactant_keys,
            "product_keys": self.product_keys,
            "is_pressure_dependent": self.is_pressure_dependent,
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
        t3_keys = {
            "kinetics_method", "kinetics_source", "kinetics_comment",
            "t3_status", "t3_index", "rmg_index", "created_at_iteration", "reasons",
            "qm_label", "rmg_label", "reactant_keys", "product_keys",
            "rmg_reaction_family",
        }
        t3_kwargs = {}
        for k in t3_keys:
            if k in reaction_dict:
                t3_kwargs[k] = reaction_dict.pop(k)
        obj = cls(**t3_kwargs, **reaction_dict)
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
        if not smiles_r or not smiles_p or not all(smiles_r) or not all(smiles_p):
            raise ValueError(f"""Could not find smiles for one or more species
                                 got: reactants: {smiles_r}
                                      products: {smiles_p}""")
        return "+".join(smiles_r)+"<=>"+"+".join(smiles_p)

    def is_isomorphic(self, other: 'T3Reaction') -> bool:
        """
        Check if this reaction is isomorphic to another reaction.
        Two reactions are isomorphic if they have the same reactants and products
        (considering species isomorphism), regardless of order.
        
        Args:
            other (T3Reaction): The other reaction to compare with.
            
        Returns:
            bool: True if reactions are isomorphic, False otherwise.
        """
        if not isinstance(other, T3Reaction):
            return False

        if len(self.r_species) != len(other.r_species) or len(self.p_species) != len(other.p_species):
            return False

        for r_spc in self.r_species:
            found_match = False
            for other_r_spc in other.r_species:
                if hasattr(r_spc, 'is_isomorphic') and r_spc.is_isomorphic(other_r_spc):
                    found_match = True
                    break
                elif hasattr(r_spc, 'label') and hasattr(other_r_spc, 'label') and r_spc.label == other_r_spc.label:
                    found_match = True
                    break
            if not found_match:
                return False

        for p_spc in self.p_species:
            found_match = False
            for other_p_spc in other.p_species:
                if hasattr(p_spc, 'is_isomorphic') and p_spc.is_isomorphic(other_p_spc):
                    found_match = True
                    break
                elif hasattr(p_spc, 'label') and hasattr(other_p_spc, 'label') and p_spc.label == other_p_spc.label:
                    found_match = True
                    break
            if not found_match:
                return False

        return True


def remove_bad_arc_keys(species_dict: dict) -> dict:
    """
    Remove keys from the species dict that ARC doesn't expect in __init__ and would error on.
    Includes calculated results, internal state flags, and output geometries.
    """
    bad_keys = {
        # --- Attributes that are NOT Init Arguments ---

        # Identity / Metadata
        'original_label',
        'index',
        'symmetry_number',

        # Calculated Energies & Properties
        'e_elect',
        'e0',
        't1',
        'zmat',
        'bond_corrections',  # Warning: This IS an arg, but often calculated internally. Keep if you want to force it.
        # If T3 re-calculates it, remove it. Usually safe to keep if valid dict.
        # I'll leave it out of bad_keys for now, but watch out for it.

        # Geometry & Conformer Results
        'initial_xyz',  # Mapped to 'xyz' in from_dict, so remove original key
        'final_xyz',  # Mapped to 'xyz' in from_dict, so remove original key
        'conf_is_isomorphic',
        'conformers',  # List of results
        'conformer_energies',
        'conformers_before_opt',
        'cheap_conformer',
        'most_stable_conformer',
        'recent_md_conformer',
        'rotors_dict',  # Complex internal dictionary
        'number_of_rotors',
        '_radius',
        '_is_linear',
        '_number_of_atoms',
        'mol_list',  # Derived from mol

        # Thermo / Transport Outputs
        'thermo',  # Result object
        'rmg_thermo',
        'long_thermo_description',
        'transport_data',

        # Settings / Levels (Saved state, not always args)
        'opt_level',
        'freq_level',
        'composite_level',
        'scan_res',

        # TS Specific Internal State
        'ts_guesses',
        'ts_report',
        'ts_checks',
        'ts_conf_spawned',
        'tsg_spawned',
        'ts_guesses_exhausted',
        'successful_methods',
        'unsuccessful_methods',
        'chosen_ts',
        'chosen_ts_list',
        'chosen_ts_method',
        'rxn_zone_atom_indices',

        # Files / Paths
        'arkane_file',
        'checkfile',
        'keep_mol',
        'neg_freqs_trshed'
    }

    return {k: v for k, v in species_dict.items() if k not in bad_keys}
