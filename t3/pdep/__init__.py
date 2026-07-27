"""
t3 pdep module

Deciding which pressure-dependent reaction networks deserve expensive QM refinement, and
parsing the RMG network files that decision is based on.
"""

from t3.pdep.cache import (
    sa_cache_metadata_path,
    validate_sa_cache,
    write_sa_cache_metadata,
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
    DEFAULT_MIN_DELTA_LN_K,
    E0_PERTURBATION_J_PER_MOL,
    PDepNetworkSelection,
    SELECTOR_VERSION,
    SensitiveTransitionState,
    coefficient_floor,
    resolve_direction_key,
    select_from_sa_dict,
    select_sensitive_wells,
)
