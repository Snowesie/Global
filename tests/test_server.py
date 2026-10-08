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


def test_password_gate(configured, monkeypatch):
    monkeypatch.setattr(config, "APP_PASSWORD", "s3cret")
    client = TestClient(app)
    assert client.get("/api/health").status_code == 200  # Render's health check stays open
    res = client.get("/")
    assert res.status_code == 401 and res.headers["www-authenticate"].startswith("Basic")
    assert client.get("/api/status", auth=("anyone", "wrong")).status_code == 401
    assert client.get("/api/status", auth=("anyone", "s3cret")).status_code == 200


def test_chunked_upload(configured, kg_sites, monkeypatch):
    from dna_app import server
    monkeypatch.setattr(server, "CHUNK_BYTES", 10_000)
    path = configured / "genome.txt"
    write_23andme(path, kg_sites, simulate_genotypes(kg_sites, {"EAS": 1.0}))
    data = path.read_bytes()
    client = TestClient(app)

    start = client.post("/api/uploads", params={"size": len(data), "filename": "genome.txt"}).json()
    uid = start["upload_id"]
    assert start["chunk_bytes"] == 10_000
    assert client.post(f"/api/uploads/{uid}/complete").status_code == 409  # nothing sent yet

    offset = 0
    while offset < len(data):
        res = client.put(f"/api/uploads/{uid}", params={"offset": offset}, content=data[offset:offset + 10_000])
        assert res.status_code == 200
        # A resent chunk (lost response) is refused, and the client resumes from the server's count.
        assert client.put(f"/api/uploads/{uid}", params={"offset": offset}, content=b"x").status_code == 409
        offset = client.get(f"/api/uploads/{uid}").json()["received"]
    assert client.put(f"/api/uploads/{uid}", params={"offset": 0}, content=b"x" * 20_000).status_code == 409

    res = client.post(f"/api/uploads/{uid}/complete")
    assert res.status_code == 202
    assert _wait(client, res.json()["job_id"])["state"] == "done"
    result = client.get(f"/api/jobs/{uid}/result").json()
    assert result["file"]["size"] == len(data)
    assert result["ancestry"]["populations"][0]["code"] == "EAS"


def test_chunk_larger_than_allowed_is_rejected(configured, monkeypatch):
    from dna_app import server
    monkeypatch.setattr(server, "CHUNK_BYTES", 100)
    client = TestClient(app)
    assert client.post("/api/uploads", params={"size": config.MAX_UPLOAD_BYTES + 1}).status_code == 413
    uid = client.post("/api/uploads", params={"size": 150}).json()["upload_id"]
    assert client.put(f"/api/uploads/{uid}", params={"offset": 0}, content=b"x" * 101).status_code == 413
    assert client.get(f"/api/uploads/{uid}").json()["received"] == 0
    assert (config.UPLOAD_DIR / f"{uid}.upload").stat().st_size == 0  # partial chunk rolled back
    assert client.put(f"/api/uploads/{uid}", params={"offset": 0}, content=b"x" * 100).json()["received"] == 100


def test_uploads_wait_for_reference_build(configured, monkeypatch):
    monkeypatch.setattr(config, "DATA_DIR", configured)
    (configured / "build_status.txt").write_text("Downloading ClinVar…")
    client = TestClient(app)
    assert client.get("/api/status").json()["building"] == "Downloading ClinVar…"
    assert client.post("/api/uploads", params={"size": 10}).status_code == 503
