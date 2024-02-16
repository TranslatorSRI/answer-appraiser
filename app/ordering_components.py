"""Compute scores for each result in the given message."""
import redis
from tqdm import tqdm
import traceback

from .config import settings
from .clinical_evidence.compute_clinical_evidence import compute_clinical_evidence
from .novelty.compute_novelty import compute_novelty


redis_pool = redis.ConnectionPool(
    host=settings.redis_host,
    port=settings.redis_port,
    db=0,
    password=settings.redis_password,
)


def get_confidence(result, message, logger):
    """
    This function iterates through the answers from multiple ARAs,
    It multiplies values of (1- score(ara[i])) for each ara
    Finally this product value is subtracted from 1
    """
    score_product = 1
    for analysis in result.get("analyses") or []:
        if analysis.get("score") is not None:
            score_product = score_product * (1 - analysis["score"])
    confidence_score = 1 - score_product
    return confidence_score


def get_clinical_evidence(result, message, logger, db_conn):
    return compute_clinical_evidence(result, message, logger, db_conn)


async def get_novelty(message, logger):
    novelty_df = await compute_novelty(message, logger)
    novelty_dict = novelty_df.to_dict(orient="index")
    novelty_scores = {
        node["drug"]: node["novelty_score"] for node in novelty_dict.values()
    }
    return novelty_scores


async def get_ordering_components(message, logger):
    logger.debug(f"Computing scores for {len(message['results'])} results")
    db_conn = redis.Redis(connection_pool=redis_pool)
    novelty_scores = {}
    try:
        novelty_scores = await get_novelty(message, logger)
    except Exception:
        logger.error(f"Novelty score failed: {traceback.format_exc()}")
    for result in tqdm(message.get("results") or []):
        confidence = 0.0
        try:
            confidence = get_confidence(result, message, logger)
        except Exception:
            logger.error(f"Confidence score failed: {traceback.format_exc()}")
        clinical_evidence_score = 0.0
        try:
            clinical_evidence_score = get_clinical_evidence(
                result,
                message,
                logger,
                db_conn,
            )
        except Exception:
            logger.error(f"Clinical evidence score failed: {traceback.format_exc()}")
        result["ordering_components"] = {
            "confidence": confidence,
            "clinical_evidence": clinical_evidence_score,
            "novelty": 0.0,
        }
        for node_bindings in result.get("node_bindings", {}).values():
            for node_binding in node_bindings:
                result["ordering_components"]["novelty"] = novelty_scores.get(
                    node_binding["id"], 0.0
                )
