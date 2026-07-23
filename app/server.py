import json
import logging
import os
import traceback
import warnings
import zstandard

from fastapi import HTTPException, status, Request
from fastapi.responses import JSONResponse, Response
from starlette.middleware.cors import CORSMiddleware
from uuid import uuid4

from .config import settings
from .lmdb_store import open_env
from .logger import setup_logger, get_logger
from .trapi import TRAPI
from .ordering_components import get_ordering_components

setup_logger()
LOGGER = logging.getLogger(__name__)

openapi_args = dict(
    title="SRI Answer Appraiser",
    version="0.8.2",
    terms_of_service="",
    description="SRI service that provides metrics for scoring and ordering of results",
    trapi="1.6.0",
    biolink_version="4.2.0",
    contact={
        "name": "Max Wang",
        "email": "max@covar.com",
        "x-id": "https://github.com/maximusunc",
        "x-role": "responsible developer",
    },
)

if settings.openapi_server_url:
    openapi_args["servers"] = [
        {
            "url": settings.openapi_server_url,
            "x-maturity": settings.openapi_server_maturity,
            "x-location": settings.openapi_server_location,
        },
    ]

APP = TRAPI(**openapi_args)

APP.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

EXAMPLE = {
    "message": {
        "query_graph": {
            "nodes": {
                "n0": {"ids": ["MESH:D008687"]},
                "n1": {"categories": ["biolink:Disease"]},
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
                "MESH:D008687": {
                    "categories": ["biolink:SmallMolecule"],
                    "name": "Metformin",
                    "attributes": [],
                },
                "MONDO:0005148": {
                    "categories": ["biolink:Disease"],
                    "name": "type 2 diabetes mellitus",
                    "attributes": [],
                },
            },
            "edges": {
                "n0n1": {
                    "subject": "MESH:D008687",
                    "object": "MONDO:0005148",
                    "predicate": "biolink:treats",
                    "sources": [
                        {
                            "resource_id": "infores:kp0",
                            "resource_role": "primary_knowledge_source",
                        }
                    ],
                    "attributes": [],
                    "knowledge_level": "knowledge_assertion",
                    "agent_type": "manual_agent",
                }
            },
        },
        "results": [
            {
                "node_bindings": {
                    "n0": {"ids": ["MESH:D008687"]},
                    "n1": {"ids": ["MONDO:0005148"]},
                },
                "analyses": [
                    {
                        "resource_id": "kp0",
                        "edge_bindings": {"n0n1": {"ids": ["n0n1"]}},
                    }
                ],
            }
        ],
    }
}


@APP.post("/get_appraisal")
async def sync_get_appraisal(request: Request):
    qid = str(uuid4())[:8]
    logger = get_logger(qid, "INFO")
    logger.info("Starting sync appraisal")
    compressed = False
    if request.headers.get("content-encoding") == "zstd":
        try:
            raw_body = await request.body()
            query = json.loads(zstandard.decompress(raw_body))
            compressed = True
        except Exception:
            return Response("Invalid request. Failed to decompress and ingest.", 400)
    else:
        try:
            query = await request.json()
        except json.JSONDecodeError:
            return Response("Invalid request. Failed to parse JSON.", 400)
    log_level = query.get("log_level") or "INFO"
    logger.setLevel(log_level)
    message = query["message"]
    if not message.get("results"):
        return JSONResponse(
            content={"status": "Rejected", "description": "No Results.", "job_id": qid},
            status_code=400,
        )
    try:
        await get_ordering_components(message, logger)
    except Exception:
        logger.error(f"Something went wrong while appraising: {traceback.format_exc()}")
    if compressed:
        query = zstandard.compress(json.dumps(query).encode())
    else:
        query = json.dumps(query)
    logger.info("Done appraising")
    return Response(query)


def _check_lmdb_ready(path):
    """Raise 503 unless the LMDB store at ``path`` is present and readable."""
    try:
        env = open_env(path)
        with env.begin() as txn:
            txn.stat()
    except Exception:
        raise HTTPException(status_code=status.HTTP_503_SERVICE_UNAVAILABLE)


@APP.get("/clinical_evidence_ready")
def check_clinical_evidence_readiness():
    """Check if the clinical evidence LMDB store is present and readable."""
    _check_lmdb_ready(settings.lmdb_path)


@APP.get("/publications_ready")
def check_publications_readiness():
    """Check if the publications LMDB store is present and readable."""
    _check_lmdb_ready(settings.publications_lmdb_path)
