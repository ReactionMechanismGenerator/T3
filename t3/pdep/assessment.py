#!/usr/bin/env python3
# encoding: utf-8

"""
t3.pdep.assessment module

The durable record of what T3 decided about every P-dep network it considered in an iteration --
**including the networks it could not decide about at all.**

Why this exists as its own type. Until this module, the only in-run record of a network's fate was
appended at exactly one site deep inside the success path. Every fail-closed exit before that site
produced no record whatsoever: the network was logged and then forgotten. On a real 12-network trial
that meant 7 networks -- the majority -- left no durable trace that they had been considered and
skipped, while the log's "3 qualified" read like an unqualified success. Fail-closed was and remains
the correct BEHAVIOUR; the defect was that refusing to answer looked exactly like never having been
asked.

Why it is not a :class:`t3.pdep.selector.PDepNetworkSelection`. A selection is RANKABLE, and
``selection_rank_key`` sorts a ``not_evaluated`` one into tier 1, ahead of evaluated negatives --
because an absent measurement genuinely might have qualified. That is right for a selection and
wrong for a placeholder: recording unevaluable networks as selections would make every future
consumer that ranks the full list treat "we never found out" as a stronger candidate than "we
checked and it does not qualify". So a selection is kept as an OPTIONAL NESTED payload here, and the
outer record carries the outcome.

**The invariant this type exists to enforce is that the record cannot disagree with itself.** A
durable file is believed; one that says ``status: qualified`` beside ``reason_code:
sa_all_methods_failed``, or claims a selector verdict for a network whose SA never parsed, would
mislead with the full authority of provenance. Every reason code therefore implies exactly one
status and one rule about the nested selection (see :mod:`t3.pdep.reason_codes`), and construction
refuses any combination that violates either -- at the site that built it, rather than at read time
months later.

**The record is frozen, and the nested selection is snapshotted rather than aliased.** Freezing the
outer record alone would not be enough: a ``PDepNetworkSelection`` is mutable, so a caller holding a
reference to the one it passed in could flip ``qualified`` afterwards and leave a validated record
quietly contradicting itself. The selection is deep-copied at construction so this record describes
the decision as it stood when the decision was made.
"""

import copy
import os
from dataclasses import dataclass

from t3.pdep.reason_codes import (INTERNAL_ERROR_REASON_CODES,
                                  PRE_SELECTOR_REASON_CODES,
                                  REASON_CODE_STATUS,
                                  REASON_EVALUATED_NO_UNCERTAIN_TS,
                                  REASON_INTERNAL_ERROR,
                                  REASON_NETWORK_DISCOVERY_FAILED,
                                  REASON_NETWORK_INPUT_WRITE_FAILED,
                                  REASON_NETWORK_PARSE_FAILED,
                                  REASON_QUALIFIED_UNCERTAIN_TS,
                                  REASON_SA_ALL_METHODS_FAILED,
                                  REASON_SA_OUTPUT_MALFORMED,
                                  REASON_SA_OUTPUT_MISSING,
                                  REASON_SA_OUTPUT_UNREADABLE,
                                  REASON_SA_STRUCTURES_MISSING,
                                  REASON_SELECTOR_DIRECTION_ENTRY_MALFORMED,
                                  REASON_SELECTOR_DIRECTION_UNRESOLVED,
                                  REASON_SELECTOR_MALFORMED_CONDITIONS_NO_TS_ROWS,
                                  REASON_SELECTOR_NEGATIVE_INCOMPLETE_DATA,
                                  REASON_SELECTOR_NO_TS_ROWS,
                                  REASON_SELECTOR_NO_USABLE_TS_ROWS,
                                  REASON_SELECTOR_SA_PAYLOAD_MALFORMED,
                                  REASON_SELECTOR_TS_PROVENANCE_UNASSESSED,
                                  REASON_SELECTOR_TS_RESPONSE_BELOW_FLOOR,
                                  REASON_SPECIES_LABEL_MAPPING_FAILED,
                                  SELECTION_BEARING_REASON_CODES,
                                  STATUS_EVALUATED_NEGATIVE,
                                  STATUS_INTERNAL_ERROR,
                                  STATUS_NOT_EVALUATED,
                                  STATUS_QUALIFIED,
                                  VALID_ASSESSMENT_REASON_CODES,
                                  VALID_ASSESSMENT_STATUSES,
                                  )
from t3.pdep.selector import PDepNetworkSelection


__all__ = ['ASSESSMENT_ENVELOPE_SCHEMA_VERSION',
           'ASSESSMENT_RECORD_FILE_NAME',
           'ASSESSMENT_RECORD_SCHEMA_VERSION',
           'PDepNetworkAssessment',
           'assessments_record_path',
           ]


# The assessment record is written under the ITERATION directory, beside the QM budget record it
# complements: the budget describes what happened to the networks that QUALIFIED, and this describes
# what happened to every network that was looked at, so the two belong at the same level.
ASSESSMENT_RECORD_FILE_NAME = 't3_pdep_network_assessments.yml'

# The on-disk SHAPE of one assessment record -- not what it means, only what it looks like. Bump on
# a field being added, removed, renamed, or changing type. The nested selection carries its own
# ``selection_schema_version``/``selection_algorithm_version`` FIELDS and is versioned independently;
# this number deliberately says nothing about it. Mirrors BUDGET_RECORD_SCHEMA_VERSION in
# t3.pdep.budget.
ASSESSMENT_RECORD_SCHEMA_VERSION = 1

# The on-disk shape of the FILE that holds a list of those records: which keys the envelope has and
# what the list is called. Deliberately a separate constant with a separate on-disk key, even though
# both are 1 today and were introduced together. Sharing one number would make each a hostage of the
# other: adding a field to a record would force the envelope to claim a change it did not undergo,
# and renaming the list key would force every record ever written to be re-stamped. The two shapes
# can and should move independently, and a version number that has to lie is worse than none.
ASSESSMENT_ENVELOPE_SCHEMA_VERSION = 1


def assessments_record_path(iteration_directory: str) -> str:
    """
    The path of the PDep network assessment record for one T3 iteration.

    Args:
        iteration_directory (str): The T3 iteration directory (``self.paths['iteration']``).

    Returns:
        str: ``iteration_directory``, joined with ``ASSESSMENT_RECORD_FILE_NAME``.
    """
    return os.path.join(iteration_directory, ASSESSMENT_RECORD_FILE_NAME)


def _validate_string_sequence(value, field_name: str, network_id: str) -> tuple:
    """
    Coerce a sequence of strings to a tuple, refusing a bare string.

    A bare ``str`` is iterable, so ``tuple('MSC')`` is ``('M', 'S', 'C')`` -- a silent corruption
    that would render three nonsense ME methods into a durable record and read back as though T3
    had genuinely tried them.

    Args:
        value: The value to coerce.
        field_name (str): The field being validated, for the error message.
        network_id (str): The network being recorded, for the error message.

    Raises:
        ValueError: If ``value`` is a bare string or is not a sequence of strings.

    Returns:
        tuple: ``value`` as a tuple of strings.
    """
    if isinstance(value, (str, bytes)):
        raise ValueError(f'PDep assessment for network {network_id} got a bare string for '
                         f'{field_name} ({value!r}); pass a sequence of strings, since iterating a '
                         f'string would silently record it one character at a time.')
    try:
        entries = tuple(value)
    except TypeError:
        raise ValueError(f'PDep assessment for network {network_id} got a non-iterable '
                         f'{field_name} ({value!r}).')
    for entry in entries:
        if not isinstance(entry, str):
            raise ValueError(f'PDep assessment for network {network_id} got a non-string entry in '
                             f'{field_name}: {entry!r} ({type(entry).__name__}).')
    return entries


def _validate_optional_string(value, field_name: str, network_id: str):
    """
    Check that a value is either ``None`` or a string.

    Args:
        value: The value to check.
        field_name (str): The field being validated, for the error message.
        network_id (str): The network being recorded, for the error message.

    Raises:
        ValueError: If ``value`` is neither ``None`` nor a string.
    """
    if value is not None and not isinstance(value, str):
        raise ValueError(f'PDep assessment for network {network_id} got a non-string {field_name}: '
                         f'{value!r} ({type(value).__name__}).')


@dataclass(frozen=True)
class PDepNetworkAssessment:
    """
    What T3 decided about one P-dep network offer, and why.

    One record per network OFFER, not per network: the same network can be reached from more than
    one sensitive observable reaction in a single iteration, and each offer asked its own question
    ("is THIS reaction's rate sensitive to an uncertain transition state?"). De-duplicating by
    network here would silently drop the answer to all but one of them.

    Args:
        network_id (str): The network file stem, e.g. ``'network4_2'``; the join key downstream
            consumers use.
        iteration (int): The T3 iteration this assessment was made in.
        status (str): One of ``VALID_ASSESSMENT_STATUSES``; determined by ``reason_code``.
        reason_code (str): One of ``VALID_ASSESSMENT_REASON_CODES``.
        secondary_reason_codes (tuple): Further codes that ALSO applied. The selector can reach more
            than one refusal in a single pass -- e.g. rows discarded as non-finite AND a selected
            transition state whose provenance never got assessed. Those are independent defects in
            the evidence, not a cause and its symptom, so reporting only the first would undercount
            the second wherever these are tallied.
        network_path (str, optional): The network file this offer was about.
        network_source_hash (str, optional): The ``'sha256:<hex>'`` content hash of that file, which
            is what binds this record to the bytes it was made about -- ``network_id`` names a file,
            not a content, and matches every later revision of it.
        observable_label (str, optional): The observable whose sensitivity produced this offer.
        sa_rank_index (int, optional): The offer's rank in that observable's sensitivity list.
        chemkin_reaction (str, optional): The offering reaction in Chemkin labels.
        network_reaction (str, optional): The offering reaction in network labels.
        requested_me_methods (tuple): The ME methods configured to be tried, in order.
        final_method (str, optional): The ME method that actually produced the SA, if any.
        sa_path (str, optional): The Arkane sensitivity YAML that was read, if any.
        cache_status (str, optional): One of ``'generated'``, ``'cached_valid'``, ``'cached_rejected'``.
        warnings (tuple): Human-readable diagnosis. The reason code says WHICH failure; these say
            what it looked like.
        selection (PDepNetworkSelection, optional): The selector's own record, deep-copied AND
            rendered to an immutable snapshot at construction, so neither the caller's alias nor
            this record's own copy can change what is later persisted. Required for a
            selector-produced code, forbidden for a code from before the selector ran, optional for
            an internal error -- see :mod:`t3.pdep.reason_codes`. Must describe the same network as
            this record, and must not contradict the provenance fields they both carry.
        schema_version (int): Pinned to ``ASSESSMENT_RECORD_SCHEMA_VERSION``.

    Raises:
        ValueError: If any field is of the wrong type, or the status and reason code disagree, or
            the presence of ``selection`` does not match what the reason code implies.
    """
    network_id: str
    iteration: int = 0
    status: str = STATUS_NOT_EVALUATED
    reason_code: str | None = None
    secondary_reason_codes: tuple = ()
    network_path: str | None = None
    network_source_hash: str | None = None
    observable_label: str | None = None
    sa_rank_index: int | None = None
    chemkin_reaction: str | None = None
    network_reaction: str | None = None
    requested_me_methods: tuple = ()
    final_method: str | None = None
    sa_path: str | None = None
    cache_status: str | None = None
    warnings: tuple = ()
    selection: PDepNetworkSelection | None = None
    schema_version: int = ASSESSMENT_RECORD_SCHEMA_VERSION

    def __post_init__(self):
        """Refuse any record that contradicts itself, at the site that built it."""
        if not isinstance(self.network_id, str) or not self.network_id:
            raise ValueError(f'A PDep assessment needs a non-empty string network_id, got '
                             f'{self.network_id!r} ({type(self.network_id).__name__}).')
        if self.reason_code not in VALID_ASSESSMENT_REASON_CODES:
            raise ValueError(f'Unrecognized PDep assessment reason_code {self.reason_code!r} for '
                             f'network {self.network_id}; expected one of '
                             f'{VALID_ASSESSMENT_REASON_CODES}.')
        if self.status not in VALID_ASSESSMENT_STATUSES:
            raise ValueError(f'Unrecognized PDep assessment status {self.status!r} for network '
                             f'{self.network_id}; expected one of {VALID_ASSESSMENT_STATUSES}.')
        implied_status = REASON_CODE_STATUS[self.reason_code]
        if self.status != implied_status:
            raise ValueError(f'PDep assessment for network {self.network_id} carries status '
                             f'{self.status!r} with reason_code {self.reason_code!r}, which implies '
                             f'status {implied_status!r}. A record whose status and reason disagree '
                             f'would misreport this network wherever it is later read.')
        # ``bool`` is a subclass of ``int``, so ``isinstance(True, int)`` is True and ``True == 1``
        # -- both checked explicitly, or ``iteration=True``/``schema_version=True`` would be stored
        # and rendered as ``true`` in YAML.
        if isinstance(self.iteration, bool) or not isinstance(self.iteration, int) or self.iteration < 0:
            raise ValueError(f'The PDep assessment iteration must be a non-negative integer, '
                             f'got {self.iteration!r} ({type(self.iteration).__name__}) for network '
                             f'{self.network_id}.')
        if isinstance(self.schema_version, bool) or self.schema_version != ASSESSMENT_RECORD_SCHEMA_VERSION:
            raise ValueError(f'PDep assessment for network {self.network_id} carries '
                             f'schema_version {self.schema_version!r}, but this code writes '
                             f'{ASSESSMENT_RECORD_SCHEMA_VERSION}. Nothing here migrates an old '
                             f'record forward; re-derive it rather than guessing.')
        if self.sa_rank_index is not None and (isinstance(self.sa_rank_index, bool)
                                               or not isinstance(self.sa_rank_index, int)
                                               or self.sa_rank_index < 0):
            raise ValueError(f'The PDep assessment sa_rank_index must be a non-negative integer or '
                             f'None, got {self.sa_rank_index!r} for network {self.network_id}.')
        for field_name in ('network_path', 'network_source_hash', 'observable_label',
                           'chemkin_reaction', 'network_reaction', 'final_method', 'sa_path',
                           'cache_status'):
            _validate_optional_string(getattr(self, field_name), field_name, self.network_id)
        self._validate_secondary_reason_codes()
        self._validate_selection()
        # ``frozen=True`` blocks ordinary assignment, so normalization goes through
        # ``object.__setattr__`` -- the same idiom PDepBudgetRecord uses.
        object.__setattr__(self, 'requested_me_methods',
                           _validate_string_sequence(self.requested_me_methods,
                                                     'requested_me_methods', self.network_id))
        object.__setattr__(self, 'warnings',
                           _validate_string_sequence(self.warnings, 'warnings', self.network_id))
        # Snapshot the selection, TWICE over, because ``frozen=True`` freezes the reference and not
        # the mutable object behind it -- so every invariant checked above can be undone the instant
        # this returns, and it is the undone state that would otherwise be serialized.
        #
        # The deep copy defends against the CALLER's alias: whoever passed the selection in still
        # holds it and could flip ``qualified`` afterwards. The rendered snapshot defends against
        # this record's OWN copy, which is handed out freely through ``record.selection`` -- that
        # second hole is the one the deep copy alone was mistakenly believed to close.
        # ``as_dict()`` serializes the rendering captured here, so no post-construction mutation can
        # reach the file, and it refuses outright if the live selection has drifted from it rather
        # than quietly persisting one story while the in-memory record tells another.
        if self.selection is not None:
            object.__setattr__(self, 'selection', copy.deepcopy(self.selection))
        object.__setattr__(self, '_selection_snapshot',
                           self.selection.as_dict() if self.selection is not None else None)

    def _validate_secondary_reason_codes(self):
        """Check that every secondary code is real, and that none of them repeats the primary."""
        codes = _validate_string_sequence(self.secondary_reason_codes, 'secondary_reason_codes',
                                          self.network_id)
        for code in codes:
            if code not in VALID_ASSESSMENT_REASON_CODES:
                raise ValueError(f'Unrecognized PDep assessment secondary reason_code {code!r} for '
                                 f'network {self.network_id}.')
            if code == self.reason_code:
                raise ValueError(f'PDep assessment for network {self.network_id} repeats its '
                                 f'primary reason_code {code!r} among its secondary codes, which '
                                 f'would double-count it.')
        if len(set(codes)) != len(codes):
            raise ValueError(f'PDep assessment for network {self.network_id} repeats a secondary '
                             f'reason_code: {codes!r}.')
        object.__setattr__(self, 'secondary_reason_codes', codes)

    def _validate_selection(self):
        """Check the nested selection against what the reason code implies about it."""
        if self.selection is not None and not isinstance(self.selection, PDepNetworkSelection):
            raise ValueError(f'PDep assessment for network {self.network_id} got a '
                             f'{type(self.selection).__name__} as its selection; only a '
                             f'PDepNetworkSelection can be nested here. Anything else would pass '
                             f'construction and fail later, when the record is serialized.')
        if self.reason_code in INTERNAL_ERROR_REASON_CODES:
            # A crash may happen before or after the selector returned, so either presence or absence
            # of a selection is coherent, and the field-by-field agreement checks below are skipped:
            # a crash part-way through populating a record can legitimately leave its two halves
            # disagreeing, and refusing it would discard the breadcrumb over the very inconsistency
            # the breadcrumb exists to report.
            #
            # The identity claim is NOT skipped. "These two halves are mid-update" is a coherent
            # thing for a crash record to say; "this is evidence about a different network" is not,
            # at any status.
            if self.selection is not None and self.selection.network_id != self.network_id:
                raise ValueError(f'PDep assessment for network {self.network_id} nests a selection '
                                 f'that is about network {self.selection.network_id!r}. A crash does '
                                 f'not make another network\'s evidence relevant to this one.')
            return
        if self.reason_code in PRE_SELECTOR_REASON_CODES:
            if self.selection is not None:
                raise ValueError(f'PDep assessment for network {self.network_id} carries a nested '
                                 f'selection with reason_code {self.reason_code!r}, which names a '
                                 f'failure that occurred before the selector ran. That selection '
                                 f'cannot be about this offer.')
            return
        if self.selection is None:
            raise ValueError(f'PDep assessment for network {self.network_id} has reason_code '
                             f'{self.reason_code!r}, which the selector produced, but carries no '
                             f'nested selection. The selection is the evidence for the code.')
        if self.status == STATUS_QUALIFIED and not self.selection.qualified:
            raise ValueError(f'PDep assessment for network {self.network_id} reports status '
                             f'{STATUS_QUALIFIED!r} but nests a selection that did not qualify.')
        if self.status == STATUS_EVALUATED_NEGATIVE and self.selection.qualified:
            raise ValueError(f'PDep assessment for network {self.network_id} reports status '
                             f'{STATUS_EVALUATED_NEGATIVE!r} but nests a selection that qualified.')
        self._validate_shared_provenance()

    def _validate_shared_provenance(self):
        """
        Check that the record and its nested selection describe the SAME network and the same work.

        Six fields are carried by both. They are duplicated for a reason -- the pre-selector codes
        forbid a nested selection precisely because there is none, and those records still need
        ``network_path``/``sa_path``/``network_id`` to be worth anything -- but every duplicated pair
        is somewhere the two halves of one record can contradict each other. A record naming one SA
        run while the evidence inside it names another is worse than a record naming neither: it
        reads as fully-provenanced and is wrong.

        ``network_id`` must match exactly, since that is the identity claim itself, and attaching
        another network's evidence is exactly the misattribution :mod:`t3.pdep.reason_codes` says the
        category rules exist to prevent. The rest are only required not to CONFLICT: absence is not
        disagreement, and the funnel legitimately knows things the selector never recorded (and the
        other way round), so demanding that both halves be populated would refuse ordinary records.
        """
        if self.selection.network_id != self.network_id:
            raise ValueError(f'PDep assessment for network {self.network_id} nests a selection that '
                             f'is about network {self.selection.network_id!r}. Evidence from another '
                             f'network cannot support a conclusion about this one.')
        for field_name, selection_field_name in (('network_source_hash', 'network_source_hash'),
                                                 ('network_reaction', 'network_reaction'),
                                                 ('sa_path', 'sa_path'),
                                                 ('cache_status', 'cache_status'),
                                                 ('final_method', 'method')):
            own = getattr(self, field_name)
            nested = getattr(self.selection, selection_field_name, None)
            if own is not None and nested is not None and own != nested:
                raise ValueError(f'PDep assessment for network {self.network_id} and its nested '
                                 f'selection disagree about {field_name}: the record says {own!r}, '
                                 f'the selection says {nested!r} (as {selection_field_name!r}). One '
                                 f'of the two is describing different work.')

    def as_dict(self) -> dict:
        """
        Render this assessment as plain JSON/YAML-safe types.

        The nested selection is rendered by its OWN ``as_dict()`` rather than flattened into this
        mapping: it carries ``selection_schema_version``/``selection_algorithm_version`` as fields
        precisely so it stays self-describing when nested, and flattening would strand those
        versions in an envelope that does not describe them.

        Returns:
            dict: The YAML-safe rendering of this record.
        """
        return {'assessment_record_schema_version': self.schema_version,
                'network_id': self.network_id,
                'iteration': self.iteration,
                'status': self.status,
                'reason_code': self.reason_code,
                'secondary_reason_codes': list(self.secondary_reason_codes),
                'network_path': self.network_path,
                'network_source_hash': self.network_source_hash,
                'observable_label': self.observable_label,
                'sa_rank_index': self.sa_rank_index,
                'chemkin_reaction': self.chemkin_reaction,
                'network_reaction': self.network_reaction,
                'requested_me_methods': list(self.requested_me_methods),
                'final_method': self.final_method,
                'sa_path': self.sa_path,
                'cache_status': self.cache_status,
                'warnings': list(self.warnings),
                # Always present, never omitted: an absent key and a null value read alike to a
                # careless consumer, so "no selection" is stated rather than implied.
                'selection': self._persisted_selection(),
                }

    def _persisted_selection(self):
        """
        Return the nested selection exactly as it was when this record was validated.

        Serializing the rendering captured at construction, rather than re-rendering the live
        object, is what makes the snapshot claim true: a caller holding ``record.selection`` cannot
        change what reaches the file. The drift check turns the remaining silent failure -- an
        in-memory record telling one story while the file it wrote tells another -- into a loud one.
        It cannot fire on the ordinary path, where nothing mutates a selection after handing it over.

        Raises:
            ValueError: If the live nested selection no longer matches the snapshot.

        Returns:
            dict, optional: The captured rendering, deep-copied so a caller editing the returned
                mapping cannot reach back into this record, or ``None`` if there is no selection.
        """
        if self._selection_snapshot is None:
            return None
        if self.selection.as_dict() != self._selection_snapshot:
            raise ValueError(f'The nested selection of the PDep assessment for network '
                             f'{self.network_id} no longer matches the one this record was validated '
                             f'against; it was mutated after construction. Persisting it would '
                             f'record evidence that no invariant here ever checked.')
        return copy.deepcopy(self._selection_snapshot)
