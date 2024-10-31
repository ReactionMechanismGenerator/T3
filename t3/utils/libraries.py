"""
t3 utils libraries module
for working with RMG thermo and kinetics libraries
"""

import datetime
import os
import shutil
import time

from typing import Dict, TYPE_CHECKING

from rmgpy.data.kinetics import KineticsLibrary
from rmgpy.data.thermo import ThermoLibrary
from rmgpy.kinetics import Arrhenius, KineticsData
from rmgpy.reaction import Reaction
from rmgpy.thermo import NASAPolynomial, NASA, ThermoData, Wilhoit
from rmgpy.species import Species

if TYPE_CHECKING:
    from t3.logger import Logger


def add_to_rmg_libraries(library_name: str,
                         shared_library_name: str,
                         paths: Dict[str, str],
                         logger: 'Logger',
                         ):
    """
    Creates RMG libraries in the RMG database repository if they don't already exist,
    and appends with the respective entries from the libraries generated by ARC.

    Args:
        library_name (str): The name of the RMG library.
        shared_library_name (str): The name of an RMG database library shared between T3 projects.
        paths (Dict[str, str]): T3's dictionary of paths.
        logger (Logger): Instance of T3's Logger class.
    """
    for token in ['thermo', 'kinetics']:
        arc_lib_path, t3_lib_path, shared_lib_path = \
            paths[f'ARC {token} lib'], paths[f'T3 {token} lib'], paths[f'shared T3 {token} lib']
        if token == 'thermo':
            local_context = {
                'ThermoData': ThermoData,
                'Wilhoit': Wilhoit,
                'NASAPolynomial': NASAPolynomial,
                'NASA': NASA,
            }
        else:
            local_context = {
                'KineticsData': KineticsData,
                'Arrhenius': Arrhenius,
            }
        for to_lib_path, lib_name, race in zip([shared_lib_path, t3_lib_path],
                                               [library_name, shared_library_name],
                                               [True, False]):
            if os.path.isfile(arc_lib_path) and to_lib_path is not None and os.path.isfile(to_lib_path):
                append_to_rmg_library(library_name=lib_name,
                                      from_lib_path=arc_lib_path,
                                      to_lib_path=to_lib_path,
                                      local_context=local_context,
                                      lib_type=token,
                                      logger=logger,
                                      race=race,
                                      )
            else:
                # The destination library (T3's or the shared) doesn't exist. Just copy the library generated by ARC.
                if os.path.isfile(arc_lib_path) and to_lib_path is not None:
                    if not os.path.isdir(os.path.dirname(to_lib_path)):
                        os.makedirs(os.path.dirname(to_lib_path))
                    shutil.copy(arc_lib_path, to_lib_path)


def append_to_rmg_library(library_name: str,
                          from_lib_path: str,
                          to_lib_path: str,
                          local_context: dict,
                          lib_type: str,
                          logger: 'Logger',
                          race: bool = False,
                          ):
    """
    Append the entries from the ARC-generated library to an RMG library.

    Args:
        library_name (str): The name of the RMG library.
        from_lib_path (str): The path to the ARC-generated library.
        to_lib_path (str): The path to the RMG library to append to.
        local_context (dict): The local context to use when loading the libraries.
        lib_type (str): The type of the library (either 'thermo' or 'kinetics').
        logger (Logger): Instance of T3's Logger class.
        race (bool, optional): Whether to take measures to avoid a race condition when appending to the library.
    """
    race_path = os.path.join(os.path.dirname(to_lib_path), f'{library_name}.race')
    if race:
        race_free = check_race_condition(race_path)
        if not race_free:
            logger.error(f'Could not write to library {to_lib_path} due to a race condition.\n'
                         f'Check whether it is safe to delete the {race_path} file to continue.')
            return
    from_lib, to_lib = (ThermoLibrary(), ThermoLibrary()) if lib_type == 'thermo' \
        else (KineticsLibrary(), KineticsLibrary())
    to_lib.load(path=to_lib_path, local_context=local_context, global_context=dict())
    from_lib.load(path=from_lib_path, local_context=local_context, global_context=dict())
    from_description = from_lib.long_desc
    description_to_append = '\n'
    append = False
    for line in from_description.splitlines():
        if 'Overall time since project initiation' in line:
            append = False
        if append:
            description_to_append += line + '\n'
        if 'Considered the following' in line:
            append = True
    to_lib.long_desc += description_to_append
    for entry in from_lib.entries.values():
        skip_entry = False
        if lib_type == 'thermo':
            entry_species = Species(molecule=[entry.item])
            entry_species.generate_resonance_structures(keep_isomorphic=False, filter_structures=True)
            for existing_entry in to_lib.entries.values():
                if entry_species.is_isomorphic(existing_entry.item):
                    if entry.label != existing_entry.label:
                        logger.warning(f"Not adding species {entry.label} to the {library_name} thermo library, "
                                       f"the species seems to already exist under the label {existing_entry.label}.")
                    skip_entry = True
                    break
        elif lib_type == 'kinetics':
            entry_reaction = Reaction(reactants=entry.item.reactants[:],
                                      products=entry.item.products[:],
                                      specific_collider=entry.item.specific_collider,
                                      kinetics=entry.data,
                                      duplicate=entry.item.duplicate,
                                      reversible=entry.item.reversible,
                                      allow_pdep_route=entry.item.allow_pdep_route,
                                      elementary_high_p=entry.item.elementary_high_p,
                                      )
            for existing_entry in to_lib.entries.values():
                if entry_reaction.is_isomorphic(existing_entry.item):
                    logger.warning(f"Not adding reaction {entry.label} to the {library_name} kinetics library, "
                                   f"the reaction seems to already exist under the label {existing_entry.label}.")
                    skip_entry = True
                    break
        if not skip_entry:
            to_lib.entries[entry.label] = entry
    to_lib.save(path=to_lib_path)
    lift_race_condition(race_path)


def check_race_condition(race_path: str,
                         ) -> bool:
    """
    Check for a race condition and avoid one by creating a race holder file.

    Args:
        race_path (str): The path to the race file to check.

    Returns:
        bool: Whether there is no race condition and T3 may continue (True) or an unavoidable race exists (False).
    """
    counter = 0
    while os.path.isfile(race_path):
        with open(race_path, 'r') as f:
            content = f.read()
            if content:
                creation_date = content.split(' ')[-1]
                creation_datetime = datetime.datetime.strptime(creation_date, "%H%M%S_%b%d_%Y")
                time_delta = datetime.datetime.now() - creation_datetime
                if time_delta.total_seconds() > 1000:
                    lift_race_condition(race_path)
                    return True
        counter += 1
        time.sleep(10)
        if counter > 1000:
            return False
    with open(race_path, 'w') as f:
        f.write(f'Race created at {datetime.datetime.now().strftime("%H%M%S_%b%d_%Y")}')
    return True


def lift_race_condition(race_path: str) -> None:
    """
    Lift the race condition by deleting the race holder file.

    Args:
        race_path (str): The path to the race file to check.
    """
    if os.path.isfile(race_path):
        os.remove(race_path)