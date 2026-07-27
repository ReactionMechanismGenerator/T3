"""
t3 pdep join module

The network-transition-state <-> ARC-transition-state join.

A PDep network file names its transition states in RMG's namespace (``TS1``, ``TS3``, ...), local
to that one file. ARC names them in its own (``TS0``, ``TS1``, ...), assigned from the position of
the reaction in ARC's reaction list. The two namespaces are unrelated and they collide: ``TS1`` in
``network4_1.py`` is almost never ARC's ``TS1``. A later step has to find the Arkane statmech
``.py`` that ARC wrote for a selected transition state, so the two names must be tied together
explicitly and durably -- guessing from list order or from label text has no basis.

The tie is made by *choosing* ARC's name rather than by predicting it. ``ARCReaction.ts_label`` is
honored by ARC when it is already set -- it only falls back to ``f'TS{rxn.index}'`` when the label
is ``None`` (``arc/scheduler.py:379``) -- and, unlike T3's own ``t3_index``, it survives
``ARCReaction.copy()``, which is what ``T3.add_reaction`` puts on the QM queue. So T3 assigns a
deterministic, namespaced label up front and the join is known at queue time, before ARC has run
and without depending on ARC's indexing at all.

Labels are refused rather than sanitized. Rewriting an unsafe character would let two distinct
networks (or two distinct transition states) collapse onto one ARC label, which would silently
join a network to another network's quantum-chemistry result -- a far worse outcome than declining
to queue an oddly-named network.
"""

import os

from arc.common import read_yaml_file, save_yaml_file

# The sidecar is written under the ARC project directory of the iteration that queued the work,
# next to the artifacts it describes, so a project directory is self-describing after the fact.
TS_JOIN_SIDECAR_FILE_NAME = 't3_pdep_ts_join.yml'

# Prefix for every T3-assigned ARC transition state label. It exists to keep T3's labels out of
# ARC's own ``TS<index>`` namespace: ARC assigns those to any reaction whose ``ts_label`` it had to
# fill in itself, so an unprefixed label could be handed to two different reactions.
ARC_TS_LABEL_PREFIX = 'T3PDep'

# The characters ARC accepts in a species label without rewriting it
# (``arc/settings/settings.py:278`` is broader; this is the conservative subset that is also safe
# as a path component, since the label becomes a file name under ``calcs/statmech/kinetics/TSs/``).
LABEL_SAFE_CHARS = frozenset('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_')

# A transition state was queued to ARC in this iteration; an artifact is expected to appear.
JOIN_STATUS_QUEUED = 'queued'
# The reaction was already known to T3, so ``add_reaction`` did not queue it again. Whether an
# artifact exists is a question for the post-ARC discovery step, not something to assume here.
JOIN_STATUS_ALREADY_PRESENT = 'already_present'
# The transition state was selected but could not be queued (e.g. its species could not be built).
JOIN_STATUS_NOT_QUEUED = 'not_queued'

TS_JOIN_STATUSES = (JOIN_STATUS_QUEUED, JOIN_STATUS_ALREADY_PRESENT, JOIN_STATUS_NOT_QUEUED)


def arc_ts_label(network_id: str, network_ts_label: str) -> str:
    """
    Build the deterministic ARC transition state label for one network transition state.

    The label is a pure function of the network and the transition state, so the same selection
    always resolves to the same ARC label and the same artifact path, across iterations and across
    restarts.

    Args:
        network_id (str): The network file stem, e.g. ``'network4_1'``.
        network_ts_label (str): The transition state label in the network file's namespace,
            e.g. ``'TS3'``.

    Raises:
        ValueError: If either component is empty or contains a character that ARC would rewrite.
            Refusing is deliberate: sanitizing could map two distinct transition states onto one
            label and silently join a network to the wrong quantum-chemistry result.

    Returns:
        str: The ARC transition state label, e.g. ``'T3PDep_network4_1_TS3'``.
    """
    for name, value in (('network_id', network_id), ('network_ts_label', network_ts_label)):
        if not value:
            raise ValueError(f'Cannot build an ARC transition state label from an empty {name}.')
        illegal = sorted({char for char in value if char not in LABEL_SAFE_CHARS})
        if illegal:
            raise ValueError(f"Cannot build an ARC transition state label from {name} '{value}': it contains "
                             f"the character(s) {illegal}, which ARC would rewrite. Refusing rather than "
                             f"sanitizing, since two different transition states could otherwise collapse "
                             f"onto one label and be joined to each other's results.")
    return f'{ARC_TS_LABEL_PREFIX}_{network_id}_{network_ts_label}'


def expected_ts_artifact_path(arc_project_directory: str, label: str) -> str:
    """
    The path at which ARC is expected to write the Arkane statmech input for a transition state.

    ARC writes one file per transition state, named after the transition state's own label
    (``arc/statmech/arkane.py::generate_ts_files`` -> ``generate_species_file``), so choosing the
    label fixes the path. Note that ARC deletes and re-creates this directory at the start of every
    rate pass, so a restarted ARC run destroys the previous one -- the artifact must be vendored
    when it is found, not merely pointed at.

    Args:
        arc_project_directory (str): The ARC project directory of the iteration that queued the work.
        label (str): The ARC transition state label.

    Returns:
        str: The expected path of the per-transition-state Arkane statmech ``.py`` input file.
    """
    return os.path.join(arc_project_directory, 'calcs', 'statmech', 'kinetics', 'TSs', f'{label}.py')


class TSJoinRecord:
    """
    One network transition state, and the ARC transition state it was tied to.

    Args:
        network_id (str): The network file stem, e.g. ``'network4_1'``.
        network_ts_label (str): The transition state label in the network file's namespace.
        status (str): One of ``TS_JOIN_STATUSES``.
        arc_ts_label (str, optional): The ARC transition state label assigned by ``arc_ts_label()``,
            or ``None`` when no such label could be built (e.g. ``arc_ts_label()`` refused an unsafe
            component). A transition state still gets a record in that case -- silently dropping it
            would be indistinguishable from one whose quantum chemistry simply failed.
        path_reaction_labels (tuple): The labels of the path reactions this transition state owns.
            A transition state may own several, so this is a tuple rather than a single label, and
            a path reaction label is not unique within a network file.
        path_reaction_strs (tuple): Those path reactions as ``'A + B <=> C'`` strings, which -- unlike
            the labels -- identify them unambiguously.
        t3_reaction_key (int, optional): The T3 reaction key, when one is known.
        expected_artifact_path (str, optional): Where the statmech input is expected, if resolved.
        reason (str): Why this record has the status it has, in particular when it was not queued.
    """

    def __init__(self,
                 network_id: str,
                 network_ts_label: str,
                 status: str,
                 arc_ts_label: str | None = None,
                 path_reaction_labels: tuple = tuple(),
                 path_reaction_strs: tuple = tuple(),
                 t3_reaction_key: int | None = None,
                 expected_artifact_path: str | None = None,
                 reason: str = '',
                 ):
        if status not in TS_JOIN_STATUSES:
            raise ValueError(f"Unrecognized transition state join status '{status}'; "
                             f"expected one of {list(TS_JOIN_STATUSES)}.")
        self.network_id = network_id
        self.network_ts_label = network_ts_label
        self.arc_ts_label = arc_ts_label
        self.status = status
        self.path_reaction_labels = tuple(path_reaction_labels)
        self.path_reaction_strs = tuple(path_reaction_strs)
        self.t3_reaction_key = t3_reaction_key
        self.expected_artifact_path = expected_artifact_path
        self.reason = reason

    @property
    def key(self) -> tuple:
        """
        The identity of this record: the network transition state it describes.

        Returns:
            tuple: ``(network_id, network_ts_label)``.
        """
        return self.network_id, self.network_ts_label

    def as_dict(self) -> dict:
        """
        Render this record as plain YAML-safe types.

        Returns:
            dict: A plain dict containing no tuples.
        """
        return {
            'network_id': self.network_id,
            'network_ts_label': self.network_ts_label,
            'arc_ts_label': self.arc_ts_label,
            'status': self.status,
            'path_reaction_labels': list(self.path_reaction_labels),
            'path_reaction_strs': list(self.path_reaction_strs),
            't3_reaction_key': self.t3_reaction_key,
            'expected_artifact_path': self.expected_artifact_path,
            'reason': self.reason,
        }

    @classmethod
    def from_dict(cls, record_dict: dict) -> 'TSJoinRecord':
        """
        Reconstruct a record from its ``as_dict()`` rendering.

        Args:
            record_dict (dict): The rendered record.

        Raises:
            ValueError: If a required field is missing, or the status is unrecognized.

        Returns:
            TSJoinRecord: The reconstructed record.
        """
        for required in ('network_id', 'network_ts_label', 'status'):
            if not record_dict.get(required):
                raise ValueError(f"A transition state join record is missing the required field "
                                 f"'{required}': {record_dict}")
        return cls(
            network_id=record_dict['network_id'],
            network_ts_label=record_dict['network_ts_label'],
            arc_ts_label=record_dict.get('arc_ts_label'),
            status=record_dict['status'],
            path_reaction_labels=tuple(record_dict.get('path_reaction_labels') or ()),
            path_reaction_strs=tuple(record_dict.get('path_reaction_strs') or ()),
            t3_reaction_key=record_dict.get('t3_reaction_key'),
            expected_artifact_path=record_dict.get('expected_artifact_path'),
            reason=record_dict.get('reason') or '',
        )

    def __repr__(self) -> str:
        return (f'TSJoinRecord({self.network_id}/{self.network_ts_label} -> '
                f'{self.arc_ts_label}, {self.status})')

    def __eq__(self, other) -> bool:
        if not isinstance(other, TSJoinRecord):
            return NotImplemented
        return self.as_dict() == other.as_dict()


def validate_ts_join_records(records: list) -> None:
    """
    Refuse a set of join records whose mapping is not one-to-one.

    Both directions matter. Two records for the same network transition state mean T3 does not know
    which ARC result to read back, and two records sharing one ARC label mean two network
    transition states would read back the *same* result -- silently attributing one network's
    quantum chemistry to another. Neither can be resolved after the fact, so both are refused here,
    before anything is written or queued.

    Args:
        records (list): The ``TSJoinRecord`` entries to validate.

    Raises:
        ValueError: If a network transition state or an ARC label appears more than once.
    """
    seen_keys, seen_arc_labels = dict(), dict()
    for record in records:
        if record.key in seen_keys:
            raise ValueError(f'Ambiguous transition state join: {record.network_id}/'
                             f'{record.network_ts_label} is mapped to both '
                             f'{seen_keys[record.key]} and {record.arc_ts_label}.')
        seen_keys[record.key] = record.arc_ts_label
        if record.arc_ts_label:
            if record.arc_ts_label in seen_arc_labels:
                previous = seen_arc_labels[record.arc_ts_label]
                raise ValueError(f"Ambiguous transition state join: the ARC label '{record.arc_ts_label}' is "
                                 f"claimed by both {previous[0]}/{previous[1]} and "
                                 f"{record.network_id}/{record.network_ts_label}.")
            seen_arc_labels[record.arc_ts_label] = record.key


def merge_ts_join_records(existing: list, new: list) -> list:
    """
    Merge newly built join records into those already collected in this iteration.

    The in-run path assesses one network reaction at a time and several reactions may point at the
    same network, so the same transition state can legitimately be offered more than once. An
    identical repeat is absorbed; a repeat that disagrees is an ambiguity and is refused by
    ``validate_ts_join_records``. First-seen order is preserved, so the sidecar reads in the order
    the work was decided rather than in an arbitrary one.

    Args:
        existing (list): The ``TSJoinRecord`` entries collected so far.
        new (list): The ``TSJoinRecord`` entries to add.

    Raises:
        ValueError: If the merged set is not a one-to-one mapping.

    Returns:
        list: The merged records.
    """
    merged, by_key = list(existing), {record.key: record for record in existing}
    for record in new:
        previous = by_key.get(record.key)
        if previous is not None:
            if previous == record:
                continue
            raise ValueError(f'Conflicting transition state join records for {record.network_id}/'
                             f'{record.network_ts_label}: {previous.as_dict()} vs {record.as_dict()}.')
        by_key[record.key] = record
        merged.append(record)
    validate_ts_join_records(merged)
    return merged


def ts_join_sidecar_path(arc_project_directory: str) -> str:
    """
    The path of the join sidecar for an ARC project directory.

    Args:
        arc_project_directory (str): The ARC project directory.

    Returns:
        str: The sidecar path.
    """
    return os.path.join(arc_project_directory, TS_JOIN_SIDECAR_FILE_NAME)


def write_ts_join_sidecar(arc_project_directory: str, records: list) -> str:
    """
    Write the join sidecar for an ARC project directory.

    The records are validated before anything is written, so a refused join leaves no partial
    sidecar behind for a later step to misread as authoritative.

    Args:
        arc_project_directory (str): The ARC project directory.
        records (list): The ``TSJoinRecord`` entries to write.

    Raises:
        ValueError: If the records are not a one-to-one mapping.

    Returns:
        str: The path written.
    """
    validate_ts_join_records(records)
    path = ts_join_sidecar_path(arc_project_directory)
    if not os.path.isdir(arc_project_directory):
        os.makedirs(arc_project_directory)
    save_yaml_file(path=path, content={'transition_states': [record.as_dict() for record in records]})
    return path


def read_ts_join_sidecar(arc_project_directory: str) -> list:
    """
    Read the join sidecar of an ARC project directory.

    Args:
        arc_project_directory (str): The ARC project directory.

    Raises:
        ValueError: If the sidecar is malformed, or its mapping is not one-to-one.

    Returns:
        list: The ``TSJoinRecord`` entries, or an empty list if no sidecar exists (which means no
            PDep transition states were queued in that iteration, not that the join was lost).
    """
    path = ts_join_sidecar_path(arc_project_directory)
    if not os.path.isfile(path):
        return list()
    content = read_yaml_file(path) or dict()
    if not isinstance(content, dict) or 'transition_states' not in content:
        raise ValueError(f"The transition state join sidecar at '{path}' is malformed: expected a mapping "
                         f"with a 'transition_states' key, got {type(content)}.")
    records = [TSJoinRecord.from_dict(entry) for entry in (content['transition_states'] or list())]
    validate_ts_join_records(records)
    return records
