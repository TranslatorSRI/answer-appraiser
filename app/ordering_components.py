"""Compute scores for each result in the given message."""
import os
import redis
from tqdm import tqdm

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
    This function iterates through the results from multiple ARAs,
    If only a single score is non-zero the result is thresholded to be in [0,1-eps]
    If a result has non-zero scores from multiple ARAs,
    then all the scores are added together and thresholded to be in [0,1]

    eps is set to 0.001
    """
    score_sum = 0
    non_zero_count = 0
    eps = 0.001
    for analysis in result.get("analyses") or []:
        if analysis.get("score") is not None:
            score_sum += analysis["score"]
            if analysis["score"] > 0:
                non_zero_count += 1
    if non_zero_count == 1 and score_sum > 1 - eps:
        score_sum = 1 - eps
    elif non_zero_count > 1 and score_sum > 1:
        score_sum = 1
    return score_sum


def get_clinical_evidence(result, message, logger, db_conn):
    return compute_clinical_evidence(result, message, logger, db_conn)


def get_novelty(message, logger):
    return compute_novelty(message, logger)


def get_ordering_components(message, logger):
    logger.debug(f"Computing scores for {len(message['results'])} results")
    db_conn = redis.Redis(connection_pool=redis_pool)
    novelty_scores_dict = get_novelty(message, logger).to_dict(orient="index")
    novelty_scores = {node["drug"]: node["novelty_score"] for node in novelty_scores_dict.values()}
    for result in tqdm(message.get("results") or []):
        clinical_evidence_score = get_clinical_evidence(
            result,
            message,
            logger,
            db_conn,
        )
        result["ordering_components"] = {
            "confidence": get_confidence(result, message, logger),
            "clinical_evidence": clinical_evidence_score,
            "novelty": 0,
        }
        if clinical_evidence_score == 0:
            # Only compute novelty if there is no clinical evidence
            for node_bindings in result.get("node_bindings", {}).values():
                for node_binding in node_bindings:
                    if node_binding["id"] in novelty_scores:
                        result["ordering_components"]["novelty"] = novelty_scores[node_binding["id"]]
