"""Mock clinical evidence store."""

import json
import os
import tempfile

import lmdb

from app.clinical_evidence.lmdb_store import LMDBReader


def dbMock():
    """Build a temporary LMDB store and return a reader over it."""
    tmpdir = tempfile.mkdtemp()
    path = os.path.join(tmpdir, "clinical_evidence.mdb")
    env = lmdb.open(path, subdir=False, map_size=10 * 1024 * 1024)
    with env.begin(write=True) as txn:
        txn.put(
            b"UMLS:C0021641_MONDO:0005015",
            json.dumps(
                [
                    {
                        "log_odds_ratio": 1.5,
                        "total_sample_size": 100,
                    },
                    {
                        "log_odds_ratio": 0.2,
                        "total_sample_size": 10000,
                    },
                ]
            ).encode(),
        )
    reader = LMDBReader(env.begin())
    # Keep the environment referenced so it isn't garbage-collected while the
    # reader's transaction is still in use.
    reader._env = env
    return reader


response = {
    "query_graph": {
        "nodes": {
            "n0": {"categories": ["biolink:Drug"]},
            "n1": {"ids": ["MONDO:0005015"]},
        },
        "edges": {
            "n0n1": {
                "subject": "n0",
                "object": "n1",
                "predicates": ["biolink:treats"],
            }
        },
    },
    "knowledge_graph": {
        "nodes": {
            "MONDO:0005015": {
                "categories": ["biolink:Disease"],
                "name": "Diabetes",
            },
            "UMLS:C0021641": {
                "categories": [
                    "biolink:Drub",
                ],
                "name": "Insulin",
            },
        },
        "edges": {
            "n0n1": {
                "subject": "UMLS:C0021641",
                "object": "MONDO:0005015",
                "predicate": "biolink:treats",
                "sources": [
                    {
                        "resource_id": "infores:kp0",
                        "resource_role": "primary_knowledge_source",
                    }
                ],
            },
        },
    },
    "results": [
        {
            "node_bindings": {
                "n0": {"ids": ["UMLS:C0021641"]},
                "n1": {"ids": ["MONDO:0005015"]},
            },
            "analyses": [
                {
                    "resource_id": "kp0",
                    "edge_bindings": {
                        "n0n1": {"ids": ["n0n1"]},
                    },
                }
            ],
        },
    ],
}
