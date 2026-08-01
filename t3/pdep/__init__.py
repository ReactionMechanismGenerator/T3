"""
t3 pdep module

Deciding which pressure-dependent reaction networks deserve expensive QM refinement, and
parsing the RMG network files that decision is based on.
"""

from t3.pdep.api import (
    load_pdep_budget_record,
    load_pdep_exploration_results,
    load_pdep_network_selections,
    rank_pdep_networks,
    save_pdep_budget_record,
    save_pdep_exploration_results,
    save_pdep_network_selections,
    select_pdep_network,
)
from t3.pdep.budget import (
    BUDGET_ALGORITHM_VERSION,
    BUDGET_OUTCOME_ADMITTED,
    BUDGET_OUTCOME_REFUSED,
    BUDGET_RECORD_FILE_NAME,
    BUDGET_RECORD_SCHEMA_VERSION,
    BUDGET_SKIP_DOES_NOT_FIT_REMAINING,
    BUDGET_SKIP_EXCEEDS_BUDGET,
    BUDGET_SKIP_MAX_NETWORKS_REACHED,
    BUDGET_SKIP_NOT_EVALUATED,
    PDepBudgetNetworkOutcome,
    PDepBudgetRecord,
    VALID_BUDGET_OUTCOMES,
    VALID_BUDGET_SKIP_REASON_CODES,
    apply_pdep_qm_budget,
    budget_record_path,
    build_pdep_budget_record,
)
from t3.pdep.cache import (
    SA_CACHE_CONTRACT_VERSION,
    sa_cache_metadata_path,
    validate_sa_cache,
    write_sa_cache_metadata,
)
from t3.pdep.hybrid import (
    HybridNetworkResult,
    QMEnergySettings,
    write_hybrid_network_input_file,
)
from t3.pdep.mesolver import (
    ArkaneMESolverAdapter,
    MESolverAdapter,
    mesolver_factory,
    register_mesolver_adapter,
)
from t3.pdep.parser import (
    PDepNetwork,
    PDepPathReaction,
    parse_pdep_network_file,
    parse_pdep_network_text,
)
from t3.pdep.selector import (
    CACHE_STATUS_CACHED_REJECTED,
    CACHE_STATUS_CACHED_VALID,
    CACHE_STATUS_GENERATED,
    CACHE_STATUS_UNVALIDATED,
    DEFAULT_MIN_DELTA_LN_K,
    E0_PERTURBATION_J_PER_MOL,
    EVALUATION_STATUS_EVALUATED,
    EVALUATION_STATUS_NOT_EVALUATED,
    PDepNetworkSelection,
    SELECTION_ALGORITHM_VERSION,
    SELECTION_SCHEMA_VERSION,
    SensitiveTransitionState,
    coefficient_floor,
    resolve_direction_key,
    select_from_sa_dict,
    select_sensitive_wells,
)
