import time

from fastapi.testclient import TestClient

from dna_app import config
from dna_app.server import app

from .conftest import simulate_genotypes, write_23andme


def _wait(client, job_id):
    for _ in range(200):
        job = client.get(f"/api/jobs/{job_id}").json()
        if job["state"] in ("done", "error"):
            return job
        time.sleep(0.05)
    raise AssertionError("job did not finish")


def test_upload_analyse_fetch_delete(configured, kg_sites):
    path = configured / "genome.txt"
    write_23andme(path, kg_sites, simulate_genotypes(kg_sites, {"EUR": 1.0}))
    client = TestClient(app)

    status = client.get("/api/status").json()
    assert status["max_upload_bytes"] == config.MAX_UPLOAD_BYTES
    assert status["databases"]["clinvar"]["records"] == "14"

    res = client.post("/api/upload", params={"filename": "../../genome.txt"}, content=path.read_bytes())
    assert res.status_code == 202
    job_id = res.json()["job_id"]
    assert _wait(client, job_id)["state"] == "done"
    assert not list(config.UPLOAD_DIR.iterdir())  # raw upload removed

    result = client.get(f"/api/jobs/{job_id}/result").json()
    assert result["file"]["name"] == "genome.txt"
    assert result["ancestry"]["populations"][0]["code"] == "EUR"

    assert client.delete(f"/api/jobs/{job_id}").json() == {"deleted": True}
    assert client.get(f"/api/jobs/{job_id}/result").status_code == 404


def test_rejects_oversized_upload(configured, monkeypatch):
    monkeypatch.setattr(config, "MAX_UPLOAD_BYTES", 1000)
    client = TestClient(app)
    res = client.post("/api/upload", content=b"x" * 1001)
    assert res.status_code == 413
    # Streamed body without Content-Length is also capped.
    res = client.post("/api/upload", content=iter([b"x" * 600, b"x" * 600]))
    assert res.status_code == 413
    assert not list(config.UPLOAD_DIR.iterdir())


def test_bad_file_reports_error(configured):
    client = TestClient(app)
    job_id = client.post("/api/upload", content=b"not a genome\n").json()["job_id"]
    job = _wait(client, job_id)
    assert job["state"] == "error" and "Could not recognise" in job["error"]


def test_invalid_job_id(configured):
    client = TestClient(app)
    assert client.get("/api/jobs/..%2F..%2Fetc").status_code == 404
    assert client.get("/").status_code == 200
