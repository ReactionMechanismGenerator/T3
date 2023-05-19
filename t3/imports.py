"""
This module contains functionality to import user settings and fill in default values from T3's settings.
"""

import os
import sys

import t3.settings.t3_settings as t3_settings
from t3.settings.t3_submit import submit_scripts


# Common imports where the user can optionally put a modified copy of settings.py or t3_submit.py file under ~/.t3
home = os.getenv("HOME") or os.path.expanduser("~")
local_t3_path = os.path.join(home, '.t3')

local_t3_settings_path = os.path.join(local_t3_path, 't3_settings.py')
settings = {key: val for key, val in vars(t3_settings).items() if '__' not in key}
if os.path.isfile(local_t3_settings_path):
    local_settings = dict()
    if local_t3_path not in sys.path:
        sys.path.insert(0, local_t3_path)
    try:
        import t3_settings as local_settings
    except ImportError:
        print('No local .t3/t3_settings.py file found')
    if local_settings:
        local_settings_dict = {key: val for key, val in vars(local_settings).items() if '__' not in key}
        settings.update(local_settings_dict)

local_t3_submit_path = os.path.join(local_t3_path, 't3_submit.py')
if os.path.isfile(local_t3_submit_path):
    local_submit_scripts = dict()
    if local_t3_path != sys.path[0]:
        sys.path.insert(0, local_t3_path)
    try:
        from t3_submit import submit_scripts as local_submit_scripts
    except ImportError:
        print('No local .t3/t3_submit.py file found')
    if local_submit_scripts:
        submit_scripts.update(local_submit_scripts)
