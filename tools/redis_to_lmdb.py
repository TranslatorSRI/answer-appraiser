#!/usr/bin/env python3
"""Copy every key/value pair from a Redis database into an LMDB store.

One-off transition helper for moving the answer appraiser's static lookup
tables off Redis and onto the memory-mapped LMDB files read at request time
(see ``app/lmdb_store.py``). Both source databases hold plain string values, so
keys and values are copied verbatim as bytes; the resulting file is a drop-in
for the ``LMDBReader`` used by the scorers.

Keys are streamed with ``SCAN`` (never ``KEYS``) and fetched/written in
batches, so this is safe to run against a large production database.

Requires: ``pip install redis lmdb``

Examples
--------
Clinical evidence edges (Redis db 0) -> ./data/clinical_evidence.mdb::

    python tools/redis_to_lmdb.py --redis-db 0 \\
        --redis-password supersecretpassword \\
        --output ./data/clinical_evidence.mdb

Publication years (Redis db 1) -> ./data/publications.mdb::

    python tools/redis_to_lmdb.py --redis-db 1 \\
        --redis-password supersecretpassword \\
        --output ./data/publications.mdb
"""

import argparse

import lmdb
import redis


def copy_redis_to_lmdb(client, env, batch_size=10000, progress=None):
    """Copy all string key/value pairs from ``client`` into LMDB ``env``.

    Iterates the Redis keyspace with ``SCAN``, fetches each batch with a single
    ``MGET``, and writes it in one LMDB transaction. Keys whose value is missing
    or is not a plain string (``MGET`` returns nil) are skipped. Returns the
    number of pairs written. ``progress`` is an optional callable receiving the
    running copied-count after each batch.
    """
    copied = 0
    batch = []

    def flush():
        nonlocal copied
        if not batch:
            return
        values = client.mget(batch)
        with env.begin(write=True) as txn:
            for key, value in zip(batch, values):
                if value is None:
                    continue
                txn.put(key, value)
                copied += 1
        batch.clear()
        if progress is not None:
            progress(copied)

    for key in client.scan_iter(count=batch_size):
        batch.append(key)
        if len(batch) >= batch_size:
            flush()
    flush()
    return copied


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--redis-host", default="localhost")
    parser.add_argument("--redis-port", type=int, default=6380)
    parser.add_argument("--redis-password", default=None)
    parser.add_argument("--redis-db", type=int, default=0)
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the LMDB file to create (single-file, subdir=False)",
    )
    parser.add_argument(
        "--map-size",
        type=int,
        default=64 * 1024**3,
        help="LMDB map size in bytes; sparse, so it may exceed the dataset "
        "(default 64 GiB)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10000,
        help="Keys per SCAN/MGET/write batch (default 10000)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    client = redis.Redis(
        host=args.redis_host,
        port=args.redis_port,
        password=args.redis_password,
        db=args.redis_db,
    )
    client.ping()

    total = client.dbsize()
    print(f"Source redis db {args.redis_db} reports {total} keys")

    env = lmdb.open(args.output, subdir=False, map_size=args.map_size)
    try:
        def progress(copied):
            print(f"  copied {copied}/{total}", end="\r", flush=True)

        copied = copy_redis_to_lmdb(
            client, env, batch_size=args.batch_size, progress=progress
        )
        env.sync()
    finally:
        env.close()

    print(f"\nDone. Copied {copied} key/value pairs into {args.output}")


if __name__ == "__main__":
    main()
