"""
t3 pdep yaml_safe module

A restricted YAML loader for reading Arkane sensitivity-analysis (SA) YAML files.

Arkane's PDep sensitivity output tags its ``(T, 'K', P, 'bar')`` condition keys with
``!!python/tuple``, so the file cannot be parsed by a plain ``yaml.SafeLoader``. T3's in-run path
reads it with ``arc.common.read_yaml_file`` (``yaml.FullLoader``), which -- unlike ``SafeLoader``
-- can construct arbitrary Python objects from other tags it happens to encounter. That is
acceptable for T3's own in-run path, where the file was just written by T3's own Arkane
invocation, but ``t3.pdep.api`` is a *public* entrypoint that reads a caller-supplied path: a
caller (or anyone who can influence the file at that path) fully controls what gets constructed
by a ``FullLoader``/``Loader`` read. This module registers a constructor for ``!!python/tuple``
(and only that tag) onto ``yaml.SafeLoader``, so condition keys parse correctly while every other
non-safe tag is refused rather than constructed.
"""

import yaml


class SaYamlSafeLoader(yaml.SafeLoader):
    """A ``yaml.SafeLoader`` that additionally understands the ``!!python/tuple`` tag, and
    nothing else beyond what ``yaml.SafeLoader`` already supports."""


def _construct_python_tuple(loader: yaml.SafeLoader, node: yaml.Node) -> tuple:
    """
    Construct a real Python tuple from a ``!!python/tuple``-tagged sequence node.

    Args:
        loader (yaml.SafeLoader): The active loader.
        node (yaml.Node): The YAML sequence node tagged ``!!python/tuple``.

    Returns:
        tuple: The constructed tuple.
    """
    return tuple(loader.construct_sequence(node))


SaYamlSafeLoader.add_constructor('tag:yaml.org,2002:python/tuple', _construct_python_tuple)


def read_sa_yaml_file(path: str) -> dict | list:
    """
    Safely read an Arkane sensitivity-analysis YAML file from a (possibly caller-supplied) path.

    Args:
        path (str): The path to the SA YAML file (e.g. ``sa_coefficients.yml``).

    Raises:
        yaml.YAMLError: If the file contains a tag other than the safe set plus
            ``!!python/tuple`` (e.g. ``!!python/object...``), rather than constructing it.

    Returns:
        dict | list: The parsed content, with ``!!python/tuple``-tagged condition keys
            constructed as real tuples and every other value restricted to YAML's safe types.
    """
    with open(path, 'r') as f:
        return yaml.load(stream=f, Loader=SaYamlSafeLoader)
