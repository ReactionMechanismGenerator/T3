#!/usr/bin/env python3
# encoding: utf-8

"""
t3.pdep.reason_codes module

The closed vocabulary for WHY T3 reached the conclusion it did about one P-dep network, and the
status each reason implies.

This module holds nothing but constants and has no ``t3`` imports, deliberately. The vocabulary is
shared by a producer (:mod:`t3.pdep.selector`, which reports which of its refusal sites fired) and a
consumer (:mod:`t3.pdep.assessment`, which persists it). Putting the constants in either one would
force a dependency between them in one direction or the other: with them in the selector the durable
record type would depend on the decision engine, and with them in the assessment the decision engine
would depend on persistence -- and that second direction is what stopped ``assessment`` from being
able to import ``PDepNetworkSelection`` to type-check its own nested payload, which let any object
at all be stored there. A neutral third module lets the dependency run
``selector -> reason_codes <- assessment``, with ``assessment -> selector`` for the type only.

**Three categories, and the category decides what a record must carry.** A code the selector
produced arrives WITH the selection that produced it, and that selection is the evidence for the
code, so it must be nested. A code from before the selector ran has no selection to nest, and
attaching one would misattribute evidence from another network or another pass. An internal error
may or may not have got far enough to have one.
"""


__all__ = ['INTERNAL_ERROR_REASON_CODES',
           'PRE_SELECTOR_REASON_CODES',
           'REASON_CODE_STATUS',
           'REASON_EVALUATED_NO_UNCERTAIN_TS',
           'REASON_INTERNAL_ERROR',
           'REASON_NETWORK_DISCOVERY_FAILED',
           'REASON_NETWORK_INPUT_WRITE_FAILED',
           'REASON_NETWORK_PARSE_FAILED',
           'REASON_QUALIFIED_UNCERTAIN_TS',
           'REASON_SA_ALL_METHODS_FAILED',
           'REASON_SA_OUTPUT_MALFORMED',
           'REASON_SA_OUTPUT_MISSING',
           'REASON_SA_OUTPUT_UNREADABLE',
           'REASON_SA_STRUCTURES_MISSING',
           'REASON_SELECTOR_DIRECTION_ENTRY_MALFORMED',
           'REASON_SELECTOR_DIRECTION_UNRESOLVED',
           'REASON_SELECTOR_MALFORMED_CONDITIONS_NO_TS_ROWS',
           'REASON_SELECTOR_NEGATIVE_INCOMPLETE_DATA',
           'REASON_SELECTOR_NO_TS_ROWS',
           'REASON_SELECTOR_NO_USABLE_TS_ROWS',
           'REASON_SELECTOR_SA_PAYLOAD_MALFORMED',
           'REASON_SELECTOR_TS_PROVENANCE_UNASSESSED',
           'REASON_SELECTOR_TS_RESPONSE_BELOW_FLOOR',
           'REASON_SPECIES_LABEL_MAPPING_FAILED',
           'SELECTION_BEARING_REASON_CODES',
           'STATUS_EVALUATED_NEGATIVE',
           'STATUS_INTERNAL_ERROR',
           'STATUS_NOT_EVALUATED',
           'STATUS_QUALIFIED',
           'VALID_ASSESSMENT_REASON_CODES',
           'VALID_ASSESSMENT_STATUSES',
           ]


# The four things a network's fate can be. Three of them are outcomes: the criterion was computed
# and came out yes, or computed and came out no, or could not be computed. The fourth is not an
# outcome at all but a BUG, and it is kept separate precisely so it is never counted alongside the
# legitimately unevaluable ones -- "7 networks could not be evaluated" means something an operator
# can act on, and quietly including a T3 crash in that number would make a defect look like the
# ordinary cost of doing business.
STATUS_QUALIFIED = 'qualified'
STATUS_EVALUATED_NEGATIVE = 'evaluated_negative'
STATUS_NOT_EVALUATED = 'not_evaluated'
STATUS_INTERNAL_ERROR = 'internal_error'
VALID_ASSESSMENT_STATUSES = (STATUS_QUALIFIED, STATUS_EVALUATED_NEGATIVE, STATUS_NOT_EVALUATED,
                             STATUS_INTERNAL_ERROR)

# --- The two verdicts ----------------------------------------------------------------------------
REASON_QUALIFIED_UNCERTAIN_TS = 'qualified_uncertain_ts'
REASON_EVALUATED_NO_UNCERTAIN_TS = 'evaluated_no_uncertain_ts'

# --- Failures BEFORE the selector could run ------------------------------------------------------
# Each names a point at which T3 gave up on a network before any sensitivity verdict was possible.
# Most of these were uncaught exceptions that killed the whole T3 run until increment 38; they are
# recorded outcomes now because one bad network must not end a multi-day RMG+ARC campaign.
REASON_NETWORK_DISCOVERY_FAILED = 'network_discovery_failed'
REASON_NETWORK_INPUT_WRITE_FAILED = 'network_input_write_failed'
REASON_SA_ALL_METHODS_FAILED = 'sa_all_methods_failed'
REASON_SA_OUTPUT_MISSING = 'sa_output_missing'
REASON_SA_OUTPUT_UNREADABLE = 'sa_output_unreadable'
# The SA payload was read but is not a mapping at all. This is the FUNNEL's code, not the selector's
# -- T3 inspects ``structures`` before it can build a network reaction string, so a non-mapping
# payload is caught there and never reaches the selector. See REASON_SELECTOR_SA_PAYLOAD_MALFORMED
# for the selector's own, separately reachable, branch.
REASON_SA_OUTPUT_MALFORMED = 'sa_output_malformed'
REASON_SA_STRUCTURES_MISSING = 'sa_structures_missing'
REASON_SPECIES_LABEL_MAPPING_FAILED = 'species_label_mapping_failed'
REASON_NETWORK_PARSE_FAILED = 'network_parse_failed'

# --- Refusals the selector itself reached --------------------------------------------------------
# These correspond one-for-one to the nine sites in t3.pdep.selector that set
# ``evaluation_status = EVALUATION_STATUS_NOT_EVALUATED``. The selector distinguishes them only in
# prose warnings, so select_from_sa_dict_with_diagnostics reports which one fired; without that,
# counting "how often did the coefficient floor bite?" would mean grepping English.
#
# The first is reachable only through the selector's PUBLIC entry points, which any caller may hand
# a non-mapping payload -- unlike T3's own funnel, which cannot get here.
REASON_SELECTOR_SA_PAYLOAD_MALFORMED = 'selector_sa_payload_malformed'
REASON_SELECTOR_DIRECTION_UNRESOLVED = 'selector_direction_unresolved'
REASON_SELECTOR_DIRECTION_ENTRY_MALFORMED = 'selector_direction_entry_malformed'
REASON_SELECTOR_NO_TS_ROWS = 'selector_no_ts_rows'
REASON_SELECTOR_MALFORMED_CONDITIONS_NO_TS_ROWS = 'selector_malformed_conditions_no_ts_rows'
REASON_SELECTOR_NO_USABLE_TS_ROWS = 'selector_no_usable_ts_rows'
REASON_SELECTOR_TS_RESPONSE_BELOW_FLOOR = 'selector_ts_response_below_floor'
REASON_SELECTOR_NEGATIVE_INCOMPLETE_DATA = 'selector_negative_incomplete_data'
REASON_SELECTOR_TS_PROVENANCE_UNASSESSED = 'selector_ts_provenance_unassessed'

# --- Not an outcome ------------------------------------------------------------------------------
REASON_INTERNAL_ERROR = 'internal_error'


SELECTION_BEARING_REASON_CODES = (REASON_QUALIFIED_UNCERTAIN_TS,
                                  REASON_EVALUATED_NO_UNCERTAIN_TS,
                                  REASON_SELECTOR_SA_PAYLOAD_MALFORMED,
                                  REASON_SELECTOR_DIRECTION_UNRESOLVED,
                                  REASON_SELECTOR_DIRECTION_ENTRY_MALFORMED,
                                  REASON_SELECTOR_NO_TS_ROWS,
                                  REASON_SELECTOR_MALFORMED_CONDITIONS_NO_TS_ROWS,
                                  REASON_SELECTOR_NO_USABLE_TS_ROWS,
                                  REASON_SELECTOR_TS_RESPONSE_BELOW_FLOOR,
                                  REASON_SELECTOR_NEGATIVE_INCOMPLETE_DATA,
                                  REASON_SELECTOR_TS_PROVENANCE_UNASSESSED,
                                  )
PRE_SELECTOR_REASON_CODES = (REASON_NETWORK_DISCOVERY_FAILED,
                             REASON_NETWORK_INPUT_WRITE_FAILED,
                             REASON_SA_ALL_METHODS_FAILED,
                             REASON_SA_OUTPUT_MISSING,
                             REASON_SA_OUTPUT_UNREADABLE,
                             REASON_SA_OUTPUT_MALFORMED,
                             REASON_SA_STRUCTURES_MISSING,
                             REASON_SPECIES_LABEL_MAPPING_FAILED,
                             REASON_NETWORK_PARSE_FAILED,
                             )
# Its own category because a crash can happen at ANY point -- before the selector ran, when there is
# no selection to keep, or after it returned, when discarding the evidence already obtained would
# throw away the most useful part of the breadcrumb. So this is the one category for which a nested
# selection is optional rather than required or forbidden.
INTERNAL_ERROR_REASON_CODES = (REASON_INTERNAL_ERROR,)

VALID_ASSESSMENT_REASON_CODES = (SELECTION_BEARING_REASON_CODES
                                 + PRE_SELECTOR_REASON_CODES
                                 + INTERNAL_ERROR_REASON_CODES)

# Every reason code implies exactly one status. This mapping IS the cross-field invariant: written
# once here rather than re-derived at each construction site, so a new code cannot be added without
# deciding what it means.
REASON_CODE_STATUS = {code: STATUS_NOT_EVALUATED for code in VALID_ASSESSMENT_REASON_CODES}
REASON_CODE_STATUS[REASON_QUALIFIED_UNCERTAIN_TS] = STATUS_QUALIFIED
REASON_CODE_STATUS[REASON_EVALUATED_NO_UNCERTAIN_TS] = STATUS_EVALUATED_NEGATIVE
REASON_CODE_STATUS[REASON_INTERNAL_ERROR] = STATUS_INTERNAL_ERROR
