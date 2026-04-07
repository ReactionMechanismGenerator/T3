"""
t3 logger module

Using a custom logger to avoid interference with RMG's / ARC's loggers.
"""

import datetime
import os
import shutil
import time
from typing import Dict, List, Optional, Tuple

from arc.common import get_git_branch, get_git_commit

from t3.common import dict_to_str, t3_path, time_lapse
from t3.version import __version__ as VERSION
from t3.chem import T3Species, T3Reaction, T3Status


class Logger(object):
    """
    The T3 Logger class.

    Args:
        project (str): The project name.
        project_directory (str): The project directory path.
        verbose (Optional[int]): The logging level, optional. 10 - debug, 20 - info, 30 - warning.
                                 ``None`` to avoid logging to file.
        t0 (datetime.datetime): Initial time when the project was spawned, stored as a datetime object.

    Attributes:
        project (str): The project name.
        project_directory (str): The project directory path.
        verbose (Optional[int]): The logging level, optional. 10 - debug, 20 - info, 30 - warning.
                                 ``None`` to avoid logging to file.
        t0 (datetime.datetime): Initial time when the project was spawned, stored as a datetime object.
        log_file (str): The path to the log file.
    """

    def __init__(self,
                 project: str,
                 project_directory: str,
                 verbose: Optional[int],
                 t0: datetime.datetime,
                 ):

        self.project = project
        self.project_directory = project_directory
        self.verbose = verbose
        self.t0 = t0
        self.log_file = os.path.join(self.project_directory, 't3.log')

        if os.path.isfile(self.log_file):
            if not os.path.isdir(os.path.join(os.path.dirname(self.log_file), 'log_archive')):
                os.mkdir(os.path.join(os.path.dirname(self.log_file), 'log_archive'))
            local_time = datetime.datetime.now().strftime("%b%d_%Y_%H:%M:%S")
            log_backup_name = f't3.{local_time}.log'
            shutil.copy(self.log_file, os.path.join(os.path.dirname(self.log_file), 'log_archive', log_backup_name))
            os.remove(self.log_file)

        self.log_header()

    def log(self, message: str, level: str = 'info'):
        """
        Log a message.
    
        Args:
            message (str): The message to be logged.
            level (str, optional): The message level. Controls the prefix and suffix to be added to the message.
                                   Allowed values are: 'info' (default), 'warning', and 'error'.
        """
        if level not in ['info', 'debug', 'warning', 'error', 'always'] and level is not None:
            self.log(f'Got an illegal level argument "{level}"', level='error')
            level = 'info'
        prefix = {'debug': '', 'info': '', 'warning': '\nWARNING: ', 'error': '\n\n\nERROR: ', 'always': ''}
        suffix = {'debug': '', 'info': '', 'warning': '\n', 'error': '\n\n', 'always': ''}
        if isinstance(message, dict):
            message = dict_to_str(message)
        elif not isinstance(message, str):
            message = str(message)
        message = prefix[level] + message + suffix[level] if level is not None else message
        if self.verbose is not None:
            if level == 'debug' and self.verbose <= 10 \
                    or level == 'info' and self.verbose <= 20 \
                    or level == 'warning' and self.verbose <= 30 \
                    or level in ['error', 'always'] \
                    or level is None:
                # print to stdout
                print(message)

                if level is not None:
                    # log to file
                    with open(self.log_file, 'a') as f:
                        f.write(message + '\n')

    def debug(self, message: str):
        """
        Log a debug level message.
        """
        self.log(message=message, level='debug')

    def info(self, message: str):
        """
        Log an info level message.
        """
        self.log(message=message, level='info')

    def warning(self, message: str):
        """
        Log a warning level message.
        """
        self.log(message=message, level='warning')

    def error(self, message: str):
        """
        Log an error level message.
        """
        self.log(message=message, level='error')

    def log_header(self):
        """
        Output a header to the log.
        """
        self.log(f'T3 execution initiated on {time.asctime()}\n\n'
                 f'################################################################\n'
                 f'#                                                              #\n'
                 f'#                  The   Tandem   Tool   (T3)                  #\n'
                 f'#       for automated chemical kinetic model development       #\n'
                 f'#                                                              #\n'
                 f'#                        Version: {VERSION}{" " * (10 - len(VERSION))}                   #\n'
                 f'#                                                              #\n'
                 f'################################################################\n\n',
                 level='always')

        # Extract HEAD git commit from T3
        head, date = get_git_commit(path=t3_path)
        branch_name = get_git_branch(path=t3_path)
        if head != '' and date != '':
            self.log(f'The current git HEAD for T3 is:\n'
                     f'    {head}\n    {date}',
                     level='always')
        if branch_name and branch_name != 'main':
            self.log(f'    (running on the {branch_name} branch)\n',
                     level='always')
        else:
            self.log('\n', level='always')
        self.log(f'Starting project {self.project}', level='always')

    def log_max_time_reached(self, max_time: str):
        """
        Log that the maximum run time was reached.

        Args:
            max_time (str): The maximum T3 walltime.
        """
        execution_time = time_lapse(self.t0)
        self.log(f'Terminating T3 due to time limit.\n'
                 f'Max time set: {max_time}\n'
                 f'Current run time: {execution_time}\n', level='always')

    def log_footer(self):
        """
        Output a footer to the log.
        """
        execution_time = time_lapse(self.t0)
        self.log(f'\n\n\nTotal T3 execution time: {execution_time}', level='always')
        self.log(f'T3 execution terminated on {time.asctime()}\n', level='always')

    def log_species_to_calculate(self,
                                 species_keys: List[int],
                                 species_dict: Dict[int, T3Species],
                                 ):
        """
        Log the species to be calculated.

        Args:
            species_keys (List[int]): The T3 species indices to calculate.
            species_dict (Dict[int, T3Species]): The T3 species dictionary.
        """
        self.log(f'\nSpecies to be calculated ({len(species_keys)}):')
        for key in species_keys:
            self.log(f'{key}: {species_dict[key].qm_label} '
                     f'(reasons: {species_dict[key].reasons})')

    def log_reactions_to_calculate(self,
                                   reaction_keys: List[int],
                                   reaction_dict: Dict[int, T3Reaction],
                                   ):
        """
        Log the reactions to be calculated.

        Args:
            reaction_keys (List[int]): The T3 reaction indices to calculate.
            reaction_dict (Dict[int, T3Reaction]): The T3 reactions dictionary.
        """
        self.log(f'\nReactions to be calculated ({len(reaction_keys)}):')
        for key in reaction_keys:
            self.log(f'{key}: {reaction_dict[key].qm_label} '
                     f'(reasons: {reaction_dict[key].reasons})')

    def log_species_summary(self,
                            species_dict: Dict[int, T3Species],
                            ):
        """
        Log a summary of the species.

        Args:
            species_dict (Dict[int, T3Species]): The T3 species dictionary.
        """
        self.log('\n\n\nSpecies Summary:\n'
                 '----------------')
        for key, species in species_dict.items():
            smiles = f' "{species.mol.to_smiles()}"' if species.mol else ''
            self.log(f'{key}: {species.qm_label}{smiles} '
                     f'(status: {clean_t3_status(species)})')

    def log_reactions_summary(self,
                              reactions_dict: Dict[int, T3Reaction],
                              ):
        """
        Log a summary of the reactions.

        Args:
            reactions_dict (Dict[int, T3Reaction]): The T3 reactions dictionary.
        """
        self.log('\n\nReactions Summary:\n'
                 '------------------')
        for key, reaction in reactions_dict.items():
            self.log(f'{key}: {reaction.qm_label} '
                     f'(status: {reaction.t3_status})')

    def log_unconverged_species_and_reactions(self,
                                              species_keys: List[int],
                                              species_dict: Dict[int, T3Species],
                                              reaction_keys: List[int],
                                              reaction_dict: Dict[int, T3Reaction],
                                              ):
        """
        Log the unconverged species and reactions.

        Args:
            species_keys (List[int]): The T3 species indices that did not converge.
            species_dict (Dict[int, T3Species]): The T3 species dictionary.
            reaction_keys (List[int]): The T3 reaction indices that did not converge.
            reaction_dict (Dict[int, T3Reaction]): The T3 reactions dictionary.
        """
        if species_keys:
            self.log('\nThe following species did not converge:')
            for key in species_keys:
                smiles = f' "{species_dict[key].mol.to_smiles()}"' if species_dict[key].mol else ''
                self.log(f'{key}: {species_dict[key].qm_label}{smiles} '
                         f'(reasons: {species_dict[key].reasons})')
        if reaction_keys:
            self.log('\nThe following reactions did not converge:')
            for key in reaction_keys:
                self.log(f'{key}: {reaction_dict[key].qm_label} '
                         f'(reasons: {reaction_dict[key].reasons})')

    def log_args(self, schema: dict):
        """
        Log the arguments used in T3.

        Args:
            schema (dict): All non-default arguments.
        """
        verbose_map = {10: 'debug', 20: 'info', 30: 'warning', None: None}
        schema['verbose'] = verbose_map[schema.get('verbose', 20)]
        self.info(f'\n\nUsing the following arguments:\n\n'
                  f'{dict_to_str(schema)}')


    def _get_converged_and_unconverged_keys(self,
                                            mapping: dict,
                                            ) -> Tuple[List[int], List[int]]:
        """
        Get the converged and unconverged keys from a dictionary of species or reactions.

        Args:
            mapping (dict): The dictionary of species or reactions.

        Returns:
            Tuple[List[int], List[int]]: The converged and unconverged keys, respectively.
        """
        converged, unconverged = list(), list()
        for key, obj in mapping.items():
            if obj.t3_status == T3Status.CONVERGED:
                converged.append(key)
            else:
                unconverged.append(key)
        return converged, unconverged


def clean_t3_status(t3_object) -> str:
    """
    Returns a formatted status string for a T3 Species or Reaction.
    Prioritizes the 'SA observable' flag, otherwise cleans the T3Status enum to Title Case.

    Args:
        t3_object (Union[T3Species, T3Reaction]): The object to check.

    Returns:
        str: The cleaned status string (e.g., "SA observable", "Converged").
    """
    status = str(t3_object.t3_status)
    if "T3Status." in status:
        status = status.replace("T3Status.", "")
    return status.title()
