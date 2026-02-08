"""
t3 utils libraries module
for working with RMG thermo and kinetics libraries
"""

from dataclasses import replace
from typing import TYPE_CHECKING, Dict, List, Optional

import datetime
import os
import re
import shutil
import tempfile
import time

from arc.molecule.molecule import Molecule

from t3.utils import rmg_shim

if TYPE_CHECKING:
    from t3.logger import Logger


_LOCK_TIMEOUT_S = 3600  # 1 hour wait for lock acquisition
_LOCK_STALE_S = 3600    # 1 hour before breaking a stale lock
_LOCK_POLL_S = 10.0     # Poll every 10 seconds


def write_atomic(path: str, content: str) -> None:
    """
    Write content to a file atomically by writing to a temp file and os.replace()'ing.
    This ensures the file is never left in a partial/corrupted state.
    """
    dir_name = os.path.dirname(path) or "."
    os.makedirs(dir_name, exist_ok=True)

    fd, tmp_path = tempfile.mkstemp(dir=dir_name, text=True)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(content)
        os.replace(tmp_path, path)
    except Exception:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise


def _write_lock_file_atomic(race_path: str) -> bool:
    """
    Try to create the lock file atomically. Returns True if created, False otherwise.
    Raises OSError if permission/IO errors occur (distinct from FileExists).
    """
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC

    fd = os.open(race_path, flags, 0o644)
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        f.write(datetime.datetime.now(datetime.timezone.utc).isoformat())
    return True


def _read_lock_timestamp(race_path: str) -> Optional[datetime.datetime]:
    """
    Best-effort read of a lock timestamp from `race_path`.

    Returns:
        A timezone-aware datetime normalized to UTC, or None if the file can't be read
        or the timestamp can't be parsed.
    """
    try:
        with open(race_path, "r", encoding="utf-8") as f:
            content = f.read().strip()
    except (OSError, UnicodeDecodeError):
        return None

    if not content:
        return None

    # 1) Preferred: ISO-8601
    try:
        dt = datetime.datetime.fromisoformat(content)
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=datetime.timezone.utc)
        return dt.astimezone(datetime.timezone.utc)
    except ValueError:
        pass

    # 2) Backwards compatibility
    parts = content.split()
    if not parts:
        return None
    token = parts[-1]

    try:
        dt = datetime.datetime.strptime(token, "%H%M%S_%b%d_%Y")
        return dt.replace(tzinfo=datetime.timezone.utc)
    except ValueError:
        return None


def check_race_condition(race_path: str, logger: Optional["Logger"] = None) -> bool:
    """
    Acquire a filesystem lock using atomic file creation.
    Returns True if lock acquired, False if timed out or IO error.
    """
    start = time.monotonic()

    while True:
        try:
            try:
                if _write_lock_file_atomic(race_path):
                    return True
            except FileExistsError:
                pass
        except OSError as e:
            if logger:
                logger.error(f"IO Error trying to acquire lock at {race_path}: {e}")
            return False

        created = _read_lock_timestamp(race_path)

        # Fallback to mtime if file content is unreadable/corrupt
        if created is None:
            try:
                mtime = os.path.getmtime(race_path)
                created = datetime.datetime.fromtimestamp(mtime, tz=datetime.timezone.utc)
            except OSError:
                created = None

        if created is not None:
            age = (datetime.datetime.now(datetime.timezone.utc) - created).total_seconds()
            if age > _LOCK_STALE_S:
                if logger:
                    logger.warning(f"Breaking stale lock at {race_path} (age: {age:.0f}s)")
                try:
                    lift_race_condition(race_path)
                except OSError as e:
                    if logger:
                        logger.warning(f"Failed to break stale lock at {race_path}: {e}")
                    time.sleep(_LOCK_POLL_S)
                continue

        if time.monotonic() - start > _LOCK_TIMEOUT_S:
            return False

        time.sleep(_LOCK_POLL_S)


def lift_race_condition(race_path: str) -> None:
    """Remove the lock file (best-effort)."""
    try:
        if os.path.isfile(race_path):
            os.remove(race_path)
    except OSError:
        pass


def _extract_longdesc_append_block(long_desc: str) -> str:
    """
    Extracts the relevant block from an ARC long description to append.
    """
    lines_to_append: List[str] = []
    append_line = False
    for line in long_desc.splitlines():
        if "Overall time since project initiation" in line:
            append_line = False
        if "Considered the following" in line:
            append_line = True
        if append_line:
            lines_to_append.append(line)
    return "\n".join(lines_to_append).strip()


def _first_adjlist(molecule_field: object) -> Optional[str]:
    """
    Normalize shim Entry.molecule field (Optional[str | List[str]]) into a single adjacency list string.
    For List[str], returns the first element. Returns None if unavailable/invalid.
    """
    if molecule_field is None:
        return None
    if isinstance(molecule_field, str):
        return molecule_field
    if isinstance(molecule_field, list):
        if not molecule_field:
            return None
        first = molecule_field[0]
        return first if isinstance(first, str) else None
    return None


def load_rmg_species_dictionary_file(
    dict_path: str,
    logger: Optional["Logger"] = None,
    *,
    strict: bool = False,
) -> Dict[str, Molecule]:
    """
    Parses an RMG-style species dictionary file into a {label: Molecule} dictionary.

    Args:
        dict_path: Path to dictionary.txt
        strict: If True and the file exists but cannot be read, raise instead of returning {}.
    """
    if not os.path.exists(dict_path):
        return {}

    species_dict: Dict[str, Molecule] = {}
    try:
        with open(dict_path, "r", encoding="utf-8") as f:
            content = f.read()
    except IOError as e:
        if logger:
            logger.warning(f"Failed to read dictionary file {dict_path}: {e}")
        if strict:
            raise
        return {}

    blocks = re.split(r"\n\s*\n", content.strip())

    for block in blocks:
        lines = block.strip().splitlines()
        if len(lines) < 2:
            continue
        label = lines[0].strip()
        adj_list = "\n".join(lines[1:])
        try:
            mol = Molecule().from_adjacency_list(adj_list)
            species_dict[label] = mol
        except ValueError:
            if logger:
                logger.debug(f"Skipping unparseable dictionary block for label '{label}' in {dict_path}")

    return species_dict


def _update_species_dictionary_atomic(new_species: Dict[str, Molecule], path: str, logger: "Logger") -> None:
    """
    Updates the species dictionary atomically.
    Aborts immediately if the existing file cannot be read to prevent data loss.
    """
    if not new_species:
        return

    current_content = ""
    if os.path.exists(path):
        if not os.path.isfile(path):
            error_msg = f"Path exists but is not a regular file: {path}"
            logger.error(error_msg)
            raise ValueError(error_msg)

        try:
            with open(path, "r", encoding="utf-8") as f:
                current_content = f.read()
        except (OSError, UnicodeDecodeError):
            logger.error(f"Failed to read existing species dictionary {path}. Aborting update to prevent data loss.")
            raise

    clean_content = current_content.rstrip()
    parts: List[str] = []

    if clean_content:
        parts.append(clean_content + "\n\n")

    for label in sorted(new_species.keys()):
        mol = new_species[label]
        parts.append(f"{label}\n{mol.to_adjacency_list()}\n\n")

    final_content = "".join(parts)

    try:
        write_atomic(path, final_content)
    except OSError:
        logger.error(f"Failed to atomically write updated species dictionary {path}")
        raise


def _stable_lock_path(dest_path: str, lib_name: str, token: str) -> str:
    """
    Compute a stable lock path for a given library target.
    Includes 'token' to ensure thermo and kinetics locks don't collide if in same dir.
    """
    lock_dir = os.path.dirname(dest_path) or "."
    os.makedirs(lock_dir, exist_ok=True)
    return os.path.join(lock_dir, f"{lib_name}_{token}.race")


def append_to_rmg_libraries(
    library_name: str,
    shared_library_name: str,
    paths: Dict[str, str],
    logger: "Logger",
) -> None:
    """
    Creates RMG libraries in the RMG database repository if they don't already exist,
    and appends with the respective entries from the libraries generated by ARC.

    Uses top-level locking to protect against creation-time race conditions.
    """
    for token in ["thermo", "kinetics"]:
        arc_lib_path = paths.get(f"ARC {token} lib")
        t3_lib_path = paths.get(f"T3 {token} lib")
        shared_lib_path = paths.get(f"shared T3 {token} lib")

        # (Destination Path, Library Name, Should Lock?)
        targets = [
            (shared_lib_path, library_name, True),
            (t3_lib_path, shared_library_name, False),
        ]

        for dest_path, lib_name, use_lock in targets:
            if dest_path is None or arc_lib_path is None:
                continue
            if not os.path.exists(arc_lib_path):
                continue

            lock_acquired = False
            race_path: Optional[str] = None

            if use_lock:
                race_path = _stable_lock_path(dest_path=dest_path, lib_name=lib_name, token=token)
                if not check_race_condition(race_path, logger):
                    logger.error(f"Race condition failure or timeout at {race_path}. Skipping.")
                    continue
                lock_acquired = True

            try:
                if token == "thermo":
                    exists = os.path.isfile(dest_path)
                else:  # kinetics
                    exists = os.path.isdir(dest_path) and os.path.isfile(os.path.join(dest_path, "reactions.py"))

                if exists:
                    append_to_rmg_library(
                        library_name=lib_name,
                        source_lib_path=arc_lib_path,
                        destination_lib_path=dest_path,
                        lib_type=token,
                        logger=logger,
                    )
                    continue

                # Create new
                if token == "thermo":
                    if os.path.exists(dest_path) and os.path.isdir(dest_path):
                        logger.error(f"Cannot create thermo lib: Destination {dest_path} is a directory.")
                        continue
                    try:
                        dest_dir = os.path.dirname(dest_path) or "."
                        os.makedirs(dest_dir, exist_ok=True)
                        shutil.copy(arc_lib_path, dest_path)
                    except IOError as e:
                        logger.error(f"Failed to create new thermo library at {dest_path}: {e}")

                else:  # kinetics
                    try:
                        if not os.path.exists(dest_path):
                            dest_parent = os.path.dirname(dest_path) or "."
                            os.makedirs(dest_parent, exist_ok=True)
                            shutil.copytree(arc_lib_path, dest_path)
                        else:
                            # Edge case: path exists but isn't a directory
                            if not os.path.isdir(dest_path):
                                logger.error(
                                    f"Cannot create kinetics lib: Destination {dest_path} exists but is not a directory."
                                )
                                continue

                            # Folder exists but reactions.py missing -> copy missing files (best-effort)
                            src_rxn = os.path.join(arc_lib_path, "reactions.py")
                            src_dict = os.path.join(arc_lib_path, "dictionary.txt")

                            if os.path.isfile(src_rxn):
                                shutil.copy(src_rxn, os.path.join(dest_path, "reactions.py"))

                            # Do NOT overwrite existing dictionary
                            dest_dict = os.path.join(dest_path, "dictionary.txt")
                            if os.path.isfile(src_dict):
                                if not os.path.exists(dest_dict):
                                    shutil.copy(src_dict, dest_dict)
                                else:
                                    logger.warning(
                                        f"Existing dictionary found at {dest_dict}. "
                                        "Not overwriting during library creation. "
                                        "Library may be inconsistent if new reactions reference missing species."
                                    )

                    except IOError as e:
                        logger.error(f"Failed to create new kinetics library at {dest_path}: {e}")

            finally:
                if lock_acquired and race_path:
                    lift_race_condition(race_path)


def append_to_rmg_library(
    library_name: str,
    source_lib_path: str,
    destination_lib_path: str,
    lib_type: str,
    logger: "Logger",
) -> None:
    """
    Append entries from an ARC-generated library to an RMG library.
    Assumes external locking handled by caller.
    """
    if lib_type == "thermo":
        to_py_path = destination_lib_path
        from_py_path = source_lib_path
        to_dict_path = None
        from_dict_path = None
    elif lib_type == "kinetics":
        dest_is_dir = os.path.isdir(destination_lib_path)
        src_is_dir = os.path.isdir(source_lib_path)
        to_py_path = os.path.join(destination_lib_path, "reactions.py") if dest_is_dir else destination_lib_path
        from_py_path = os.path.join(source_lib_path, "reactions.py") if src_is_dir else source_lib_path

        to_dict_path = os.path.join(os.path.dirname(to_py_path), "dictionary.txt")
        from_dict_path = os.path.join(os.path.dirname(from_py_path), "dictionary.txt")
    else:
        raise ValueError(f"Unsupported lib_type: {lib_type}")

    # 1. Parse Libraries
    try:
        with open(from_py_path, "r", encoding="utf-8") as f:
            source_lib = rmg_shim.parse_rmg_library(f.read())
        with open(to_py_path, "r", encoding="utf-8") as f:
            destination_lib = rmg_shim.parse_rmg_library(f.read())
    except (IOError, ValueError) as e:
        logger.error(f"Failed to read/parse libraries for {library_name}: {e}")
        return

    # 2. Merge Descriptions
    description_to_append = _extract_longdesc_append_block(destination_lib.longDesc)
    if description_to_append and description_to_append.strip() not in destination_lib.longDesc:
        if destination_lib.longDesc.strip():
            destination_lib.longDesc += "\n"
        destination_lib.longDesc += description_to_append

    # 3. Calculate New Entries
    current_indices = [e.index for e in destination_lib.entries]
    max_index = max(current_indices) if current_indices else 0
    entries_to_add: List[rmg_shim.Entry] = []

    if lib_type == "thermo":
        existing_labels = {e.label for e in destination_lib.entries}

        # Only pre-existing species for clearer logging
        dest_species_by_label: Dict[str, Molecule] = {}
        for entry in destination_lib.entries:
            adj_list = _first_adjlist(entry.molecule)
            if not adj_list:
                continue
            try:
                dest_species_by_label[entry.label] = Molecule().from_adjacency_list(adj_list)
            except ValueError:
                pass

        batch_mols: List[Molecule] = []
        for entry in source_lib.entries:
            if entry.label in existing_labels:
                logger.warning(f"Thermo entry {entry.label} skipped (duplicate label).")
                continue

            adj_list = _first_adjlist(entry.molecule)
            if not adj_list:
                logger.warning(f"Skipping thermo entry {entry.label}: missing/invalid molecule.")
                continue

            try:
                entry_mol = Molecule().from_adjacency_list(adj_list)
            except ValueError:
                logger.warning(f"Skipping thermo entry {entry.label}: invalid molecule definition.")
                continue

            # Check isomorphic duplicates against pre-existing destination species (log with dest label)
            is_dup = False
            for dest_label, dest_mol in dest_species_by_label.items():
                if entry_mol.is_isomorphic(dest_mol):
                    logger.info(f"Thermo species {entry.label} is isomorphic to existing {dest_label}. Skipping.")
                    is_dup = True
                    break
            if is_dup:
                continue

            # Check within-batch isomorphic duplicates (avoid confusing “existing label” logging)
            for m in batch_mols:
                if entry_mol.is_isomorphic(m):
                    is_dup = True
                    break
            if is_dup:
                continue

            max_index += 1
            new_entry = replace(entry, index=max_index)
            entries_to_add.append(new_entry)
            batch_mols.append(entry_mol)
            existing_labels.add(entry.label)

        if entries_to_add:
            destination_lib.entries.extend(entries_to_add)
            content = rmg_shim.write_library_to_string(destination_lib)
            write_atomic(to_py_path, content)

    elif lib_type == "kinetics":
        existing_rxn_labels = {e.label for e in destination_lib.entries}

        for entry in source_lib.entries:
            if entry.label in existing_rxn_labels:
                logger.warning(f"Reaction {entry.label} skipped (duplicate label).")
                continue
            max_index += 1
            new_entry = replace(entry, index=max_index)
            entries_to_add.append(new_entry)

        if not entries_to_add:
            return

        if to_dict_path is None or from_dict_path is None:
            raise ValueError("Kinetics dictionary paths are not set.")

        # 1. Update Dictionary First (strict reads to avoid “append storm” on transient read failure)
        dest_species = load_rmg_species_dictionary_file(to_dict_path, logger=logger, strict=True)
        src_species = load_rmg_species_dictionary_file(from_dict_path, logger=logger, strict=True)

        species_to_append: Dict[str, Molecule] = {}
        for label, mol in src_species.items():
            if label not in dest_species:
                species_to_append[label] = mol
            else:
                existing_mol = dest_species[label]
                if not mol.is_isomorphic(existing_mol):
                    logger.warning(
                        f"Kinetics species {label} exists but structure differs. Using existing definition."
                    )

        dict_existed_before = os.path.exists(to_dict_path)
        dict_backup_content: Optional[str] = None

        if species_to_append and dict_existed_before:
            try:
                with open(to_dict_path, "r", encoding="utf-8") as f:
                    dict_backup_content = f.read()
            except IOError:
                dict_backup_content = None

            # If we cannot back up an existing dictionary, do not proceed (rollback would be impossible)
            if dict_backup_content is None:
                raise RuntimeError(
                    f"Cannot back up existing kinetics dictionary for rollback: {to_dict_path}. Aborting append."
                )

        dict_updated = False
        if species_to_append:
            _update_species_dictionary_atomic(species_to_append, to_dict_path, logger)
            dict_updated = True

        # 2. Update Reactions Second
        try:
            destination_lib.entries.extend(entries_to_add)
            content = rmg_shim.write_library_to_string(destination_lib)
            write_atomic(to_py_path, content)
        except Exception as e:
            logger.error(f"Failed to write reactions.py: {e}. Attempting dictionary rollback.")
            if dict_updated:
                if dict_existed_before:
                    if dict_backup_content is not None:
                        try:
                            write_atomic(to_dict_path, dict_backup_content)
                            logger.info("Dictionary rollback successful.")
                        except Exception as rollback_e:
                            logger.error(
                                f"Dictionary rollback FAILED: {rollback_e}. Library may be in inconsistent state."
                            )
                    else:
                        logger.error(
                            "Dictionary rollback unavailable (no backup content). Library may be inconsistent."
                        )
                else:
                    # Dictionary did not exist before; best-effort delete it.
                    try:
                        if os.path.isfile(to_dict_path):
                            os.remove(to_dict_path)
                            logger.info("Dictionary rollback successful (deleted newly created dictionary).")
                    except Exception as rollback_e:
                        logger.error(
                            f"Dictionary rollback FAILED (delete): {rollback_e}. Library may be in inconsistent state."
                        )
            raise

    else:
        raise ValueError(f"Unsupported lib_type: {lib_type}")
