#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_yaml_safe module

FIX D: t3.pdep.api reads caller-supplied SA/cache-metadata paths through a restricted YAML
loader (t3.pdep.yaml_safe.read_sa_yaml_file), rather than arc.common.read_yaml_file's
yaml.FullLoader, precisely because it is a public entrypoint over a path a caller controls. These
tests confirm the restricted loader (a) parses the real Arkane SA fixture identically to the
unsafe loader, including its ``!!python/tuple``-tagged condition keys, and (b) refuses any other
non-safe tag rather than silently constructing it.
"""

import os

import pytest
import yaml
from arc.common import read_yaml_file

from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.yaml_safe import read_sa_yaml_file

SA_PATH = os.path.join(TEST_DATA_BASE_PATH, 'pdep_sa', 'network4_2_MSC', 'sa_coefficients.yml')


def test_read_sa_yaml_file_matches_unsafe_loader_on_the_real_fixture():
    """Test that the restricted loader parses the real SA fixture identically to read_yaml_file,
    including reconstructing !!python/tuple-tagged condition keys as real tuples."""
    safe = read_sa_yaml_file(path=SA_PATH)
    unsafe = read_yaml_file(path=SA_PATH)
    assert safe == unsafe
    # Sanity check that the fixture actually exercises the !!python/tuple tag (a plain
    # yaml.SafeLoader would otherwise refuse this file and the parity assertion above would be
    # vacuous rather than meaningful).
    condition_keys = [key for value in safe.values() if isinstance(value, dict)
                      for key in value if isinstance(key, tuple)]
    assert condition_keys, 'Fixture is expected to contain at least one !!python/tuple condition key.'


def test_read_sa_yaml_file_refuses_other_unsafe_tags(tmp_path):
    """Test that a tag other than !!python/tuple (e.g. !!python/object) is refused, not
    constructed -- the whole point of restricting the loader to a single whitelisted tag."""
    path = str(tmp_path / 'malicious.yml')
    with open(path, 'w') as f:
        f.write("payload: !!python/object:builtins.str ['not actually safe']\n")

    with pytest.raises(yaml.YAMLError):
        read_sa_yaml_file(path=path)
