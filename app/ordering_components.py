"""Compute scores for each result in the given message."""

import contextlib
from tqdm import tqdm
import traceback

from .config import settings
from .clinical_evidence.compute_clinical_evidence import compute_clinical_evidence
from .clinical_evidence.lmdb_store import open_env, LMDBReader
from .novelty.compute_novelty import compute_novelty

# Lazily-opened, shared read-only LMDB environment for clinical evidence lookups.
_db_env = None


def _get_db_env():
    """Open (once) and return the clinical evidence LMDB environment."""
    global _db_env
    if _db_env is None:
        _db_env = open_env(settings.lmdb_path)
    return _db_env


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
        node["Result ID"]: node["novelty_score"] for node in novelty_dict.values()
    }
    return novelty_scores


async def get_ordering_components(message, logger):
    logger.debug(f"Computing scores for {len(message['results'])} results")
    novelty_scores = {}
    try:
        novelty_scores = await get_novelty(message, logger)
    except Exception:
        logger.error(f"Novelty score failed: {traceback.format_exc()}")

    # Open a single read transaction for the whole message so every clinical
    # evidence lookup is served from the same consistent snapshot. If the store
    # can't be opened, degrade gracefully (clinical evidence scores stay 0).
    try:
        txn_cm = _get_db_env().begin(buffers=False)
    except Exception:
        logger.error(f"Clinical evidence store unavailable: {traceback.format_exc()}")
        txn_cm = contextlib.nullcontext(None)

    with txn_cm as txn:
        db_conn = LMDBReader(txn)
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
                logger.error(
                    f"Clinical evidence score failed: {traceback.format_exc()}"
                )
            result["ordering_components"] = {
                "confidence": confidence,
                "clinical_evidence": clinical_evidence_score,
                "novelty": 0.0,
            }
            for binding in result.get("node_bindings", {}).values():
                for kg_id in binding["ids"]:
                    if kg_id in novelty_scores:
                        result["ordering_components"]["novelty"] = novelty_scores[kg_id]
