import json
import zstandard

from fastapi.testclient import TestClient

from app.server import APP

testclient = TestClient(APP)


def test_sync_get_appraisal_400():
    """Test calling /query endpoint."""
    response = testclient.post(
        "/get_appraisal",
        json={"message": {}},
    )
    assert response.status_code == 400


def test_zstd_compression():
    """Test that server decompressed zstandard payloads."""
    message = {
        "message": {
            "query_graph": {},
            "knowledge_graph": {},
            "results": [],
        }
    }
    payload = zstandard.compress(json.dumps(message).encode())
    response = testclient.post(
        "/get_appraisal",
        headers={
            "Content-Encoding": "zstd",
            "Content-Type": "application/json",
        },
        data=payload,
    )

    assert response.status_code == 400
    err = response.json()
    # we were able to parse the zipped payload and found no results
    assert err["description"] == "No Results."
