#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_logger module
"""

import datetime
import os
import shutil
from typing import Optional

from t3.common import DATA_BASE_PATH
import t3.logger as logger


log_project_directory = os.path.join(DATA_BASE_PATH, 'log_file_testing_dir')


def init_logger(project: str = 'project_name',
                project_directory: str = log_project_directory,
                verbose: Optional[int] = 10,
                t0: Optional[datetime.datetime] = None,
                ) -> logger.Logger:
    """Initialize the logger"""
    if t0 is None:
        t0 = datetime.datetime.now()
    if not os.path.isdir(log_project_directory):
        os.mkdir(log_project_directory)
    return logger.Logger(project=project,
                         project_directory=project_directory,
                         verbose=verbose,
                         t0=t0,
                         )


def test_initialize_logger():
    """Test initializing the logger"""
    logger_object = init_logger()
    assert isinstance(logger_object, logger.Logger)
    assert logger_object.project == 'project_name'
    assert logger_object.project_directory == log_project_directory
    assert logger_object.verbose == 10
    assert isinstance(logger_object.t0, datetime.datetime)
    shutil.rmtree(log_project_directory, ignore_errors=True)


def log_messages(logger_object):
    """A helper function for logging dummy messages"""
    logger_object.log('message A')
    logger_object.log('message B', level='debug')
    logger_object.log('message C', level='info')
    logger_object.log('message D', level='warning')
    logger_object.log('message E', level='error')
    logger_object.log('message F', level='always')
    logger_object.log('not logging this line to file', level=None)


def test_log():
    """Test logging messages"""
    logger_object = init_logger(verbose=10)  # debug
    log_messages(logger_object)
    with open(os.path.join(log_project_directory, 't3.log')) as f:
        lines = f.readlines()
    for line in ['message A\n',
                 'message B\n',
                 'message C\n',
                 'WARNING: message D\n',
                 'ERROR: message E\n',
                 'message F\n',
                 ]:
        assert line in lines
    for line in ['not logging this line to file\n',
                 ]:
        assert line not in lines
    shutil.rmtree(log_project_directory, ignore_errors=True)

    logger_object = init_logger(verbose=20)  # info
    log_messages(logger_object)
    with open(os.path.join(log_project_directory, 't3.log')) as f:
        lines = f.readlines()
    for line in ['message A\n',
                 'message C\n',
                 'WARNING: message D\n',
                 'ERROR: message E\n',
                 'message F\n',
                 ]:
        print('*****', line)
        assert line in lines
    for line in ['message B\n',
                 'not logging this line to file\n',
                 ]:
        assert line not in lines
    shutil.rmtree(log_project_directory, ignore_errors=True)

    logger_object = init_logger(verbose=30)  # warning
    log_messages(logger_object)
    with open(os.path.join(log_project_directory, 't3.log')) as f:
        lines = f.readlines()
    for line in [
                 'WARNING: message D\n',
                 'ERROR: message E\n',
                 'message F\n',
                 ]:
        assert line in lines
    for line in ['message A\n',
                 'message B\n',
                 'message C\n',
                 'not logging this line to file\n',
                 ]:
        assert line not in lines
    shutil.rmtree(log_project_directory, ignore_errors=True)


def test_log_max_time_reached():
    """Test logging reaching the maximum walltime"""
    logger_object = init_logger(verbose=30)  # warning
    logger_object.log_max_time_reached(max_time='01:00:00:00')
    with open(os.path.join(log_project_directory, 't3.log')) as f:
        lines = f.readlines()
    assert 'Terminating T3 due to time limit.\n' in lines
    shutil.rmtree(log_project_directory, ignore_errors=True)
