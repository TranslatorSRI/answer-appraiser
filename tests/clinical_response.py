"""Mock Redis."""

import fakeredis
import json


def redisMock():
    redis = fakeredis.FakeRedis()
    redis.set(
        "UMLS:C0021641_MONDO:0005015",
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
        ),
    )
    # set up mock function
    return redis


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
                "n0": [
                    {
                        "id": "UMLS:C0021641",
                    },
                ],
                "n1": [
                    {
                        "id": "MONDO:0005015",
                    },
                ],
            },
            "analyses": [
                {
                    "resource_id": "kp0",
                    "edge_bindings": {
                        "n0n1": [
                            {
                                "id": "n0n1",
                            },
                        ],
                    },
                }
            ],
        },
    ],
}
