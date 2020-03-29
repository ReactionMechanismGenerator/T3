"""
t3 logger module

Using a custom logger to avoid intervening with RMG's / ARC's loggers.
"""

import datetime
import os
import shutil
import time
from typing import Optional

from arc.common import time_lapse

from .config import species_labels_dict
from .utils import dict_to_str


log_file = None


def log(message: str,
        level: Optional[str] = 'info',
        verbose: Optional[bool] = True,
        ) -> None:
    """
    The tandem tool logging function.
    RMG and ARC have loggers that will override a conventional logger used here imported from logging
    Hence we define our own simple logging tool here.

    Args:
        message (str): The message to be logged.
        level (str, optional): The message level. Controls the prefix and suffix to be added to the message.
                               Allowed values are: 'info' (default), 'warning', and 'error'.
        verbose (Optional[bool]): Whether or not to log to file.
    """
    global log_file
    if level not in ['info', 'warning', 'error']:
        log(f'Got an illegal level argument "{level}"', level='error')
        level = 'info'
    prefix = {'info': '', 'warning': '\nWARNING: ', 'error': '\n\n\nERROR: '}
    suffix = {'info': '', 'warning': '\n', 'error': '\n\n'}
    if isinstance(message, dict):
        message = dict_to_str(message)
    elif not isinstance(message, str):
        message = str(message)
    message = prefix[level] + message + suffix[level]
    # also print to stdout
    print(message)
    # log to file
    message += '\n'
    if verbose:
        with open(log_file, 'a') as f:
            f.write(message)


def initialize_tandem_log(output_directory: str,
                          verbose: Optional[bool] = True,
                          ) -> str:
    """
    Set up the logger.

    Args:
        output_directory (str): The name of the output directory where the log file will be saved.
        verbose (Optional[bool]): Whether or not to log to file.

    Returns:
        str: The thermo library name from the previous T3 run, used for restarting.
    """
    global log_file
    thermo_library = None  # `None` upon fist call to add_rmg_libraries()
    log_file = os.path.join(output_directory, 't3.log')
    if os.path.isfile(log_file):
        with open(log_file, 'r') as f:
            for line in f:
                if 'Created the RMG Thermo library' in line:
                    thermo_library = line.split()[-1]
                    break
        if not os.path.isdir(os.path.join(os.path.dirname(log_file), 'log_archive')):
            os.mkdir(os.path.join(os.path.dirname(log_file), 'log_archive'))
        local_time = datetime.datetime.now().strftime("%H%M%S_%b%d_%Y")
        log_backup_name = 't3.' + local_time + '.log'
        shutil.copy(log_file, os.path.join(os.path.dirname(log_file), 'log_archive', log_backup_name))
        os.remove(log_file)
    log_header(verbose=verbose)
    return thermo_library


def log_header(verbose: Optional[bool] = True) -> None:
    """
    Output a header containing identifying information about the RMG-ARC feature to the log.

    Args:
        verbose (Optional[bool]): Whether or not to log to file.
    """
    global t0
    t0 = time.time()
    log(f'********    The RMG-ARC Tandem Tool (T3)   ********\n'
        f'** for automated model generation and refinement **\n\n'
        f'T3 execution initiated on {time.asctime()}',
        verbose=verbose)


def log_footer(verbose: Optional[bool] = True) -> None:
    """
    Output a footer to the log.

    Args:
        verbose (Optional[bool]): Whether or not to log to file.
    """
    global t0
    execution_time = time_lapse(t0)
    log(f'\n\nTotal T3 execution time: {execution_time}', verbose=verbose)
    log(f'T3 execution terminated on {time.asctime()}\n', verbose=verbose)


def log_species_to_calculate(species_dict: dict,
                             verbose: Optional[bool] = True,
                             ) -> None:
    """
    Report the species to be calculated in the next iteration.

    Args:
        species_dict (dict): A dictionary of RMG Species, containing the reason for calculating them.
        verbose (Optional[bool]): Whether or not to log to file.
    """
    log('Species to calculate thermodynamic data for:', verbose=verbose)
    max_label_length = max([len(label) for label in species_dict.keys()])
    max_smiles_length = max([len(spc_dict['spc'].molecule[0].to_smiles()) for spc_dict in species_dict.values()])
    space1 = ' ' * (max_label_length - len('label') + 1)
    space2 = ' ' * (max_smiles_length - len('SMILES') + 1)
    log(f'Label{space1} SMILES{space2} Reason for calculating thermo for this species', verbose=verbose)
    log(f'-----{space1} ------{space2} ----------------------------------------------', verbose=verbose)
    for label, spc_dict in species_dict.items():
        smiles = spc_dict['spc'].molecule[0].to_smiles()
        space1 = ' ' * (max_label_length - len(label) + 1)
        space2 = ' ' * (max_smiles_length - len(smiles) + 1)
        log(f'{label}{space1} {smiles}{space2} {spc_dict["reason"]}', verbose=verbose)


def log_species_summary(species_dict: dict,
                        unconverged_species: list,
                        verbose: Optional[bool] = True,
                        ) -> None:
    """
    Report species summary.
    The report will be saved as `RMG_ARC_species.log` under the run_directory the RMG run folder.

    Args:
        species_dict (dict): Keys are species labels, values are RMG Species containing the reason for calculating them.
        unconverged_species (list): Entries are RMG Species objects for which thermo calculations did not converge
                                    in ARC throughout the iterations.
        verbose (Optional[bool]): Whether or not to log to file.
    """
    global species_labels_dict
    log('\n\n\nSPECIES SUMMARY', verbose=verbose)
    log('\nSpecies for which thermodynamic data was calculate by ARC:\n', verbose=verbose)
    max_label_length = max([len(label) for label in species_dict.keys()] + [6])
    max_smiles_length = max([len(spc_dict['spc'].molecule[0].to_smiles()) for spc_dict in species_dict.values()]
                            + [6])
    space1 = ' ' * (max_label_length - len('label') + 1)
    space2 = ' ' * (max_smiles_length - len('SMILES') + 1)
    log(f'Label{space1} SMILES{space2} Reason for calculating thermo for this species', verbose=verbose)
    log(f'-----{space1} ------{space2} ----------------------------------------------', verbose=verbose)
    for label, spc_dict in species_dict.items():
        smiles = spc_dict['spc'].molecule[0].to_smiles()
        space1 = ' ' * (max_label_length - len(label) + 1)
        space2 = ' ' * (max_smiles_length - len(smiles) + 1)
        if all([label != species_labels_dict[spc.label] for spc in unconverged_species]):
            log(f'{label}{space1} {smiles}{space2} {spc_dict["reason"]}', verbose=verbose)
    if unconverged_species:
        log('\nSpecies for which thermodynamic data did not converge:', verbose=verbose)
        space1 = ' ' * (max_label_length - len('label') + 1)
        space2 = ' ' * (max_smiles_length - len('SMILES') + 1)
        log(f'         Label{space1} SMILES{space2} Reason for attempting to calculate thermo for this species',
            verbose=verbose)
        log(f'         -----{space1} ------{space2} ----------------------------------------------------------',
            verbose=verbose)
        for uc_spc in unconverged_species:
            for label, spc_dict in species_dict.items():
                if label == species_labels_dict[uc_spc.label]:
                    smiles = spc_dict['spc'].molecule[0].to_smiles()
                    space1 = ' ' * (max_label_length - len(label) + 1)
                    space2 = ' ' * (max_smiles_length - len(smiles) + 1)
                    log(f'(FAILED) {label}{space1} {smiles}{space2} {spc_dict["reason"]}', verbose=verbose)
    else:
        log('\nAll species calculated in ARC successfully converged', verbose=verbose)


def log_unconverged_species(unconverged_species: list,
                            verbose: Optional[bool] = True,
                            ) -> None:
    """
    Report unconverged species.

    Args:
        unconverged_species (list): Entries are RMG Species objects for which thermo calculations did not converge
                                    in ARC throughout the iterations.
        verbose (Optional[bool]): Whether or not to log to file.
    """
    global species_labels_dict
    if unconverged_species:
        log('\nThermodynamic calculations for the following species did NOT converge:', verbose=verbose)
        max_label_length = max([len(species_labels_dict[spc.label]) for spc in unconverged_species] + [6])
        space1 = ' ' * (max_label_length - len('label') + 1)
        log(f'Label{space1} SMILES', verbose=verbose)
        log(f'-----{space1} ------', verbose=verbose)
        for spc in unconverged_species:
            label = species_labels_dict[spc.label]
            space1 = ' ' * (max_label_length - len(label))
            log(f'{label}{space1} {spc.molecule[0].to_smiles()}', verbose=verbose)
        log('\n', verbose=verbose)
    else:
        log('\nAll species thermodynamic calculations in this iteration successfully converged.', verbose=verbose)
