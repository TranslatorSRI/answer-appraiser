"""Clinical Evidence Scoring."""
import json
import logging
import os
import redis

REDIS_PSWD = os.getenv("REDIS_PSWD", "supersecretpassword")
r = redis.Redis(
    host="0.0.0.0",
    port=6379,
    password=REDIS_PSWD,
)

def compute_clinical_evidence(
    result: dict, message, logger: logging.Logger,
):
    """Given a result, compute the clinical evidence score,

    which is the weighted average of all the clinical kp edges that contain
    the result's subject and object.
    """
    clinical_evidence_score = 0
    found_edges = []
    # loop over all analyses in the given result and append any clinical kp edges to found_edges
    for analysis in result.get("analyses") or []:
        for edge_bindings in analysis.get("edge_bindings", {}).values():
            for edge_binding in edge_bindings:
                try:
                    kg_edge = message["knowledge_graph"]["edges"][edge_binding["id"]]
                except KeyError:
                    # this is malformed TRAPI
                    logger.error("malformed TRAPI")
                    continue
                clinical_edge_id = f"{kg_edge['subject']}_{kg_edge['object']}"
                kg_edge = r.get(clinical_edge_id)
                if kg_edge is not None:
                    found_edges.extend(json.loads(kg_edge))

    # Compute the clinical evidence score given all clinical kp edges
    # Score is computed by:
    # log_odds_ratio = OR
    # total_sample_size = N
    # weight = W = N / (N + N_2)
    # score = (W * OR + W_2 * OR_2) / (W + W_2)
    complete_sample_size = 0
    for kp_edge in found_edges:
        complete_sample_size += kp_edge["total_sample_size"]
    if complete_sample_size > 0:
        total_weights = 0
        for kp_edge in found_edges:
            total_weights += kp_edge["total_sample_size"] / complete_sample_size
            clinical_evidence_score += kp_edge["log_odds_ratio"] * (
                kp_edge["total_sample_size"] / complete_sample_size
            )
        if total_weights > 0:
            clinical_evidence_score /= total_weights
    return clinical_evidence_score
