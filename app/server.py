import logging
import redis
import traceback

from fastapi import Body, BackgroundTasks, HTTPException, status
from fastapi.responses import JSONResponse
import httpx
from starlette.middleware.cors import CORSMiddleware
from uuid import uuid4

from reasoner_pydantic import AsyncQuery, AsyncQueryResponse, Response, Query

from .config import settings
from .logger import setup_logger, get_logger
from .trapi import TRAPI
from .ordering_components import get_ordering_components


setup_logger()
LOGGER = logging.getLogger(__name__)

openapi_args = dict(
    title="SRI Answer Appraiser",
    version="0.3.6",
    terms_of_service="",
    description="SRI service that provides metrics for scoring and ordering of results",
    trapi="1.4.0",
    biolink_version="3.4.2",
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
                },
                "MONDO:0005148": {
                    "categories": ["biolink:Disease"],
                    "name": "type 2 diabetes mellitus",
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
                }
            },
        },
        "results": [
            {
                "node_bindings": {
                    "n0": [{"id": "MESH:D008687"}],
                    "n1": [{"id": "MONDO:0005148"}],
                },
                "analyses": [
                    {"resource_id": "kp0", "edge_bindings": {"n0n1": [{"id": "n0n1"}]}}
                ],
            }
        ],
    }
}

ASYNC_EXAMPLE = {
    "callback": "http://test",
    **EXAMPLE,
}


async def async_appraise(message, callback, logger: logging.Logger):
    try:
        await get_ordering_components(message, logger)
    except Exception:
        logger.error(f"Something went wrong while appraising: {traceback.format_exc()}")
    logger.info("Done appraising")
    try:
        logger.info(f"Posting to callback {callback}")
        async with httpx.AsyncClient(timeout=httpx.Timeout(timeout=600.0)) as client:
            res = await client.post(callback, json=message)
            logger.info(f"Posted to {callback} with code {res.status_code}")
    except Exception as e:
        logger.error(f"Unable to post to callback {callback}.")


@APP.post("/async_get_appraisal", response_model=AsyncQueryResponse)
async def get_appraisal(
    background_tasks: BackgroundTasks,
    query: AsyncQuery = Body(..., example=ASYNC_EXAMPLE),
):
    """Appraise Answers"""
    qid = str(uuid4())[:8]
    query_dict = query.dict()
    log_level = query_dict.get("log_level") or "INFO"
    logger = get_logger(qid, log_level)
    logger.info("Starting async appraisal")
    message = query_dict["message"]
    if not message.get("results"):
        logger.warning("No results given.")
        return JSONResponse(
            content={"status": "Rejected", "description": "No Results.", "job_id": qid},
            status_code=400,
        )
    callback = query_dict["callback"]
    background_tasks.add_task(async_appraise, message, callback, logger)
    return JSONResponse(
        content={
            "status": "Accepted",
            "description": f"Appraising answers. Will send response to {callback}",
            "job_id": qid,
        },
        status_code=200,
    )


@APP.post("/get_appraisal", response_model=Response)
async def sync_get_appraisal(query: Query = Body(..., example=EXAMPLE)):
    qid = str(uuid4())[:8]
    query_dict = query.dict()
    log_level = query_dict.get("log_level") or "INFO"
    logger = get_logger(qid, log_level)
    logger.info("Starting sync appraisal")
    message = query_dict["message"]
    if not message.get("results"):
        return JSONResponse(
            content={"status": "Rejected", "description": "No Results.", "job_id": qid},
            status_code=400,
        )
    try:
        await get_ordering_components(message, logger)
    except Exception:
        logger.error(f"Something went wrong while appraising: {traceback.format_exc()}")
    logger.info("Done appraising")
    return Response(message=message)


@APP.get("/redis_ready")
def check_redis_readiness():
    """Check if redis is started and ready to accept connections"""
    try:
        r = redis.Redis(
            host=settings.redis_host,
            port=settings.redis_port,
            password=settings.redis_password,
        )
        r.ping()
    except Exception:
        raise HTTPException(status_code=status.HTTP_503_SERVICE_UNAVAILABLE)
