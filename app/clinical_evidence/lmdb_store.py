"""LMDB-backed read-only store for clinical evidence edges.

The clinical evidence data is a static, bulk-loaded lookup table keyed by
``"{subject}_{object}"`` mapping to a JSON-encoded list of clinical KP edges.
It is served out of a memory-mapped LMDB file instead of a separate Redis
process: reads are in-process (no network round-trip) and the OS page cache is
shared across worker processes.
"""

import lmdb


def open_env(path: str) -> lmdb.Environment:
    """Open the clinical evidence LMDB environment read-only.

    ``lock=False`` is safe because the file is never written to while being
    served (it is rebuilt offline). ``subdir=False`` treats ``path`` as a single
    file rather than a directory.
    """
    return lmdb.open(
        path,
        readonly=True,
        lock=False,
        readahead=False,
        subdir=False,
        max_readers=512,
    )


class LMDBReader:
    """Adapter exposing a Redis-like ``.get()`` over an LMDB read transaction.

    Keeps :func:`compute_clinical_evidence` agnostic to the backing store: it
    accepts a ``str`` key and returns the raw ``bytes`` value (or ``None`` when
    the key is absent), matching the previous ``redis.Redis`` interface. A
    ``None`` transaction is tolerated so callers can degrade gracefully when the
    store is unavailable.
    """

    def __init__(self, txn):
        self._txn = txn

    def get(self, key):
        if self._txn is None:
            return None
        if isinstance(key, str):
            key = key.encode()
        return self._txn.get(key)
