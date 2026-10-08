"""HTTP API and static UI.

Run with:  uvicorn dna_app.server:app --host 0.0.0.0 --port 8000
Set DNA_APP_PASSWORD to require a password.
"""

from __future__ import annotations

import base64
import hmac
import json
import os
import re
import secrets
import shutil
import threading
import time
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from . import config
from .analysis import analyse, database_status
from .parsers import UnsupportedFormat

STATIC_DIR = Path(__file__).parent / "static"
JOB_ID = re.compile(r"^[A-Za-z0-9_-]{16,64}$")

app = FastAPI(title="Genome Reader", docs_url="/api/docs", openapi_url="/api/openapi.json")


class PasswordGate:
    """HTTP Basic auth on every route except the health check, when a password is configured.

    The browser shows its own login prompt; any username works.
    """

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        password = config.APP_PASSWORD
        if scope["type"] != "http" or not password or scope["path"] == "/api/health":
            return await self.app(scope, receive, send)
        auth = dict(scope["headers"]).get(b"authorization", b"").decode("latin-1")
        if auth.lower().startswith("basic "):
            try:
                _, _, given = base64.b64decode(auth[6:]).decode("utf-8").partition(":")
            except ValueError:
                given = ""
            if hmac.compare_digest(given.encode(), password.encode()):
                return await self.app(scope, receive, send)
        await send({"type": "http.response.start", "status": 401, "headers": [
            (b"www-authenticate", b'Basic realm="Genome Reader", charset="UTF-8"'),
            (b"content-type", b"text/plain; charset=utf-8"),
        ]})
        await send({"type": "http.response.body", "body": b"Password required"})


app.add_middleware(PasswordGate)
_executor = ThreadPoolExecutor(max_workers=2)
_jobs: dict[str, dict] = {}
_lock = threading.Lock()


def _job_dir(job_id: str) -> Path:
    if not JOB_ID.match(job_id):
        raise HTTPException(404, "Unknown job")
    return config.RESULTS_DIR / job_id


def _update(job_id: str, **fields) -> None:
    with _lock:
        _jobs[job_id].update(fields, updated=time.time())


def _run(job_id: str, upload: Path, filename: str, size: int, assembly: str | None) -> None:
    def progress(message: str, fraction: float) -> None:
        _update(job_id, message=message, progress=round(fraction, 3))

    _update(job_id, state="running")
    try:
        result = analyse(str(upload), _job_dir(job_id), progress, assembly)
        result["file"].update(name=filename, size=size)
        result["job_id"] = job_id
        result["created"] = _jobs[job_id]["created"]
        (_job_dir(job_id) / "result.json").write_text(json.dumps(result))
        _update(job_id, state="done", progress=1.0, message="Done")
    except (UnsupportedFormat, ValueError) as e:
        _update(job_id, state="error", error=str(e))
    except Exception as e:  # noqa: BLE001
        traceback.print_exc()
        _update(job_id, state="error", error=f"Analysis failed: {e}")
    finally:
        if config.DELETE_UPLOADS:
            upload.unlink(missing_ok=True)


@app.get("/api/health")
def health() -> dict:
    return {"ok": True}


@app.get("/api/status")
def status() -> dict:
    build_status = config.DATA_DIR / "build_status.txt"
    return {
        "max_upload_bytes": config.MAX_UPLOAD_BYTES,
        "chunk_bytes": CHUNK_BYTES,
        "databases": database_status(),
        "building": build_status.read_text() if build_status.exists() else None,
    }


def _check_ready() -> None:
    path = config.DATA_DIR / "build_status.txt"
    if path.exists():
        raise HTTPException(503, path.read_text())


def _check_assembly(assembly: str) -> str | None:
    if assembly not in ("auto", "GRCh37", "GRCh38"):
        raise HTTPException(400, "assembly must be auto, GRCh37 or GRCh38")
    return None if assembly == "auto" else assembly


def _too_large() -> HTTPException:
    return HTTPException(413, f"File exceeds the {config.MAX_UPLOAD_BYTES // 1024 ** 2} MB limit")


async def _append_body(request: Request, dest: Path, limit: int) -> int:
    """Append the streamed request body to ``dest``; returns bytes written, 413 past ``limit``."""
    written = 0
    with open(dest, "ab") as out:
        async for chunk in request.stream():
            written += len(chunk)
            if written > limit:
                raise _too_large()
            out.write(chunk)
    return written


def _queue(job_id: str, path: Path, filename: str, size: int, assembly: str | None) -> JSONResponse:
    now = time.time()
    with _lock:
        _jobs[job_id] = {"state": "queued", "progress": 0.0, "message": "Queued",
                         "created": now, "updated": now, "filename": filename, "size": size}
    _executor.submit(_run, job_id, path, filename, size, assembly)
    return JSONResponse({"job_id": job_id, "size": size}, status_code=202)


def _safe_name(filename: str) -> str:
    return Path(filename).name[:200] or "upload"


@app.post("/api/upload")
async def upload(request: Request, filename: str = "upload", assembly: str = "auto") -> JSONResponse:
    """Single-request upload: stream the raw body to disk and queue analysis."""
    _check_ready()
    asm = _check_assembly(assembly)
    declared = request.headers.get("content-length")
    if declared and int(declared) > config.MAX_UPLOAD_BYTES:
        raise _too_large()
    config.UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    job_id = secrets.token_urlsafe(18)
    dest = config.UPLOAD_DIR / f"{job_id}.upload"
    try:
        size = await _append_body(request, dest, config.MAX_UPLOAD_BYTES)
    except BaseException:
        dest.unlink(missing_ok=True)
        raise
    if size == 0:
        dest.unlink(missing_ok=True)
        raise HTTPException(400, "Empty upload")
    return _queue(job_id, dest, _safe_name(filename), size, asm)


# Chunked uploads keep each request small, so a hosting proxy's body-size
# limit never applies, and a dropped chunk can be retried on its own.
CHUNK_BYTES = int(os.environ.get("DNA_CHUNK_BYTES", 32 * 1024 ** 2))
STALE_UPLOAD_SECONDS = 6 * 3600
_uploads: dict[str, dict] = {}


def _drop_stale_uploads() -> None:
    cutoff = time.time() - STALE_UPLOAD_SECONDS
    with _lock:
        stale = [k for k, u in _uploads.items() if u["updated"] < cutoff]
        for k in stale:
            _uploads.pop(k)["path"].unlink(missing_ok=True)


def _get_upload(upload_id: str) -> dict:
    with _lock:
        u = _uploads.get(upload_id)
    if u is None:
        raise HTTPException(404, "Unknown or expired upload")
    return u


@app.post("/api/uploads")
def start_upload(size: int, filename: str = "upload", assembly: str = "auto") -> dict:
    _check_ready()
    asm = _check_assembly(assembly)
    if size <= 0:
        raise HTTPException(400, "Empty upload")
    if size > config.MAX_UPLOAD_BYTES:
        raise _too_large()
    _drop_stale_uploads()
    config.UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    upload_id = secrets.token_urlsafe(18)
    path = config.UPLOAD_DIR / f"{upload_id}.upload"
    path.touch()
    with _lock:
        _uploads[upload_id] = {"path": path, "size": size, "received": 0, "filename": _safe_name(filename),
                               "assembly": asm, "updated": time.time(), "busy": False}
    return {"upload_id": upload_id, "chunk_bytes": CHUNK_BYTES}


@app.get("/api/uploads/{upload_id}")
def upload_progress(upload_id: str) -> dict:
    u = _get_upload(upload_id)
    return {"received": u["received"], "size": u["size"]}


@app.put("/api/uploads/{upload_id}")
async def upload_chunk(upload_id: str, offset: int, request: Request) -> dict:
    u = _get_upload(upload_id)
    with _lock:
        if u["busy"]:
            raise HTTPException(409, "Another chunk is being written")
        if offset != u["received"]:
            raise HTTPException(409, f"Expected offset {u['received']}")
        u["busy"] = True
    try:
        try:
            written = await _append_body(request, u["path"], min(CHUNK_BYTES, u["size"] - offset))
        except BaseException:
            with open(u["path"], "ab") as fh:  # drop the partial chunk so it can be resent
                fh.truncate(offset)
            raise
        with _lock:
            u["received"] += written
            u["updated"] = time.time()
        return {"received": u["received"]}
    finally:
        u["busy"] = False


@app.post("/api/uploads/{upload_id}/complete")
def complete_upload(upload_id: str) -> JSONResponse:
    u = _get_upload(upload_id)
    if u["busy"] or u["received"] != u["size"]:
        raise HTTPException(409, f"Received {u['received']} of {u['size']} bytes")
    with _lock:
        _uploads.pop(upload_id, None)
    return _queue(upload_id, u["path"], u["filename"], u["size"], u["assembly"])


@app.get("/api/jobs/{job_id}")
def job_status(job_id: str) -> dict:
    with _lock:
        job = _jobs.get(job_id)
    if job is None:
        if (_job_dir(job_id) / "result.json").exists():
            return {"state": "done", "progress": 1.0, "message": "Done"}
        raise HTTPException(404, "Unknown job")
    return job


@app.get("/api/jobs/{job_id}/result")
def job_result(job_id: str) -> FileResponse:
    path = _job_dir(job_id) / "result.json"
    if not path.exists():
        raise HTTPException(404, "Result not ready")
    return FileResponse(path, media_type="application/json")


@app.delete("/api/jobs/{job_id}")
def delete_job(job_id: str) -> dict:
    d = _job_dir(job_id)
    with _lock:
        job = _jobs.get(job_id)
        if job and job["state"] in ("queued", "running"):
            raise HTTPException(409, "Job is still running")
        _jobs.pop(job_id, None)
    shutil.rmtree(d, ignore_errors=True)
    return {"deleted": True}


@app.get("/")
def index() -> FileResponse:
    return FileResponse(STATIC_DIR / "index.html")


app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
