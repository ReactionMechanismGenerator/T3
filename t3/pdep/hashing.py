"""
t3 pdep hashing module

The one home for T3's content-hash format for PDep artifacts.

Everything T3 records as provenance -- the SA cache sidecar's ``network_file_hash`` and
``sa_file_hash``, the join sidecar's ``source_sha256``, a captured network's vendored digest, and a
``PDepNetworkSelection``'s ``network_source_hash`` -- has to be comparable to everything else, or
the comparisons that decide whether a cache is stale or a network changed under a decision quietly
stop meaning anything. Comparability is a property of the *format*, so the format lives in exactly
one place: this module. Two functions, differing only in whether the input is a path or bytes
already in hand.

This module deliberately imports nothing from ``t3``. It sits below ``t3.pdep.parser`` (which needs
``hash_bytes`` at parse time) and below ``t3.pdep.cache`` (which needs ``hash_file``), and
``cache`` -> ``selector`` -> ``parser`` means neither of those could import from the other.

``t3.pdep.cache.hash_file`` remains a valid import path and is named as the mandated primitive
throughout ``join.py``, ``capture.py``, and ``main.py``; it is now this module's ``hash_file``,
re-exported there rather than defined there.
"""

import hashlib

# The chunk size ``hash_file`` streams with. Network and SA files are small, but captured artifacts
# are not necessarily, and a digest primitive that loads an arbitrary file into memory is a
# footgun waiting for the first large one.
_CHUNK_SIZE = 1 << 16


def hash_bytes(data: bytes) -> str:
    """
    Content-hash bytes already in hand.

    Use this instead of ``hash_file`` whenever the bytes have ALREADY been read and the digest must
    describe exactly those bytes: re-opening the file to hash it would record the digest of content
    that was never used, if the file changed in between. That window is precisely what a provenance
    hash exists to catch, so opening it while recording one is self-defeating.

    Args:
        data (bytes): The bytes to hash.

    Returns:
        str: The prefixed ``'sha256:<hexdigest>'`` string.
    """
    return f'sha256:{hashlib.sha256(data).hexdigest()}'


def hash_file(path: str) -> str:
    """
    Content-hash a file, streaming it rather than loading it.

    Args:
        path (str): The path to the file to hash.

    Returns:
        str: The prefixed ``'sha256:<hexdigest>'`` string, identical to what ``hash_bytes`` returns
            for the same bytes.
    """
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(_CHUNK_SIZE), b''):
            digest.update(chunk)
    return f'sha256:{digest.hexdigest()}'
