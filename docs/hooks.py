"""
mkdocs build hooks for T3.

Substitutes ``{{T3_VERSION}}`` in any page's markdown with the version
defined in ``t3/version.py``, so the docs stay in sync with the package
without manually bumping the number in multiple places.
"""

import importlib.util
from pathlib import Path

_VERSION_FILE = Path(__file__).resolve().parent.parent / 't3' / 'version.py'
_spec = importlib.util.spec_from_file_location('t3_version', _VERSION_FILE)
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
T3_VERSION = _module.__version__


def on_page_markdown(markdown, **kwargs):
    return markdown.replace('{{T3_VERSION}}', T3_VERSION)
