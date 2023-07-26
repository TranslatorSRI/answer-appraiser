"""Test Clinical Evidence function."""
import logging

from app.clinical_evidence.compute_clinical_evidence import compute_clinical_evidence
from tests.clinical_response import response, redisMock

logger = logging.getLogger(__name__)


def test_clinical_evidence():
    """Test that the clinical evidence function works."""
    score = compute_clinical_evidence(
        response["results"][0],
        response,
        logger,
        redisMock(),
    )
    assert score == 0.10603553615150196
