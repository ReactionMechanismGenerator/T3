"""
t3 logger module

Using a custom logger to avoid interference with RMG's / ARC's loggers.
"""

import datetime
import os
import shutil
import time
from typing import Dict, List, Optional

from arc.common import get_git_branch, get_git_commit

from t3.common import VERSION, dict_to_str, t3_path, time_lapse


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
            local_time = datetime.datetime.now().strftime("%H%M%S_%b%d_%Y")
            log_backup_name = f't3.{local_time}].log'
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
                 f'#      Automated kinetic model generation and refinement       #\n'
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
        if branch_name and branch_name != 'master':
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
                                 species_dict: Dict[int, dict]):
        """
        Report the species to be calculated in the next iteration.
        The 'QM label' is used for reporting, it is principally the
        RMG species label as entered by the user, or the label determined
        by RMG if it does not contain forbidden characters and is
        descriptive (i.e., NOT "S(1056)").
    
        Args:
            species_keys (List[int]): Entries are T3 species indices.
            species_dict (dict): The T3 species dictionary.

        Todo:
            Log the reasons one by one with line breaks and enumerate
        """
        if len(species_keys):
            self.info('Species to calculate thermodynamic data for:')
            max_label_length = max([len(spc_dict['QM label'])
                                    for key, spc_dict in species_dict.items() if key in species_keys])
            max_smiles_length = max([len(spc_dict['object'].molecule[0].to_smiles())
                                    for key, spc_dict in species_dict.items() if key in species_keys])
            space1 = ' ' * (max_label_length - len('label') + 1)
            space2 = ' ' * (max_smiles_length - len('SMILES') + 1)
            self.info(f'Label{space1} SMILES{space2} Reason for calculating thermo')
            self.info(f'-----{space1} ------{space2} -----------------------------')
            for key in species_keys:
                spc_dict = species_dict[key]
                smiles = spc_dict['object'].molecule[0].to_smiles()
                space1 = ' ' * (max_label_length - len(spc_dict['QM label']) + 1)
                space2 = ' ' * (max_smiles_length - len(smiles) + 1)
                self.info(f"{spc_dict['QM label']}{space1} {smiles}{space2} {spc_dict['reasons']}")

    def log_species_summary(self, species_dict: Dict[int, dict]):
        """
        Report species summary.
        The report will be saved as `RMG_ARC_species.log` under the run_directory the RMG run folder.
    
        Args:
            species_dict (dict): The T3 species dictionary.
        """
        if species_dict:
            self.info('\n\n\nSPECIES SUMMARY')
            converged_keys, unconverged_keys = list(), list()
            for key, spc_dict in species_dict.items():
                if spc_dict['converged']:
                    converged_keys.append(key)
                else:
                    unconverged_keys.append(key)

            max_label_length = max([len(spc_dict['QM label'])
                                    for spc_dict in species_dict.values()] + [6])
            max_smiles_length = max([len(spc_dict['object'].molecule[0].to_smiles())
                                     for spc_dict in species_dict.values()] + [6])
            if len(converged_keys):
                self.info('\nSpecies for which thermodynamic data was calculate:\n')
                space1 = ' ' * (max_label_length - len('label') + 1)
                space2 = ' ' * (max_smiles_length - len('SMILES') + 1)
                self.info(f'Label{space1} SMILES{space2} Reason for calculating thermo for this species')
                self.info(f'-----{space1} ------{space2} ----------------------------------------------')
                for key in converged_keys:
                    smiles = species_dict[key]['object'].molecule[0].to_smiles()
                    space1 = ' ' * (max_label_length - len(species_dict[key]['QM label']) + 1)
                    space2 = ' ' * (max_smiles_length - len(smiles) + 1)
                    self.info(f"{species_dict[key]['QM label']}{space1} {smiles}{space2} {species_dict[key]['reasons']}")
            else:
                self.info('\nNo species thermodynamic calculation converged!')

            if len(unconverged_keys):
                self.info('\nSpecies for which thermodynamic data did not converge:')
                space1 = ' ' * (max_label_length - len('label') + 1)
                space2 = ' ' * (max_smiles_length - len('SMILES') + 1)
                self.info(f'         Label{space1} SMILES{space2} Reason for calculating thermo for this species')
                self.info(f'         -----{space1} ------{space2} ----------------------------------------------')
                for key in unconverged_keys:
                    smiles = species_dict[key]['object'].molecule[0].to_smiles()
                    space1 = ' ' * (max_label_length - len(species_dict[key]['QM label']) + 1)
                    space2 = ' ' * (max_smiles_length - len(smiles) + 1)
                    self.info(f"(FAILED) {species_dict[key]['QM label']}{space1} "
                              f"{smiles}{space2} {species_dict[key]['reasons']}")
            else:
                self.info('\nAll species calculated in ARC successfully converged')

    def log_unconverged_species(self,
                                species_keys: List[int],
                                species_dict: Dict[int, dict]):
        """
        Report unconverged species.
    
        Args:
            species_keys (List[int]): Entries are T3 species indices.
            species_dict (dict): The T3 species dictionary.
        """
        if len(species_keys):
            self.info('\nThermodynamic calculations for the following species did NOT converge:')
            max_label_length = max([len(spc_dict['QM label'])
                                    for key, spc_dict in species_dict.items() if key in species_keys] + [6])
            space1 = ' ' * (max_label_length - len('label') + 1)
            self.info(f'Label{space1} SMILES')
            self.info(f'-----{space1} ------')
            for key in species_keys:
                label = species_dict[key]['QM label']
                space1 = ' ' * (max_label_length - len(label))
                self.info(f"{label}{space1} {species_dict[key]['object'].molecule[0].to_smiles()}")
            self.info('\n')
        else:
            self.info('\nAll species thermodynamic calculations in this iteration successfully converged.')

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
