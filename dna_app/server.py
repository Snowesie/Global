"""HTTP API and static UI.

Run with:  uvicorn dna_app.server:app --host 0.0.0.0 --port 8000
"""

from __future__ import annotations

import json
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


@app.get("/api/status")
def status() -> dict:
    return {"max_upload_bytes": config.MAX_UPLOAD_BYTES, "databases": database_status()}


@app.post("/api/upload")
async def upload(request: Request, filename: str = "upload", assembly: str = "auto") -> JSONResponse:
    """Stream the raw request body to disk (no multipart buffering) and queue analysis."""
    if assembly not in ("auto", "GRCh37", "GRCh38"):
        raise HTTPException(400, "assembly must be auto, GRCh37 or GRCh38")
    declared = request.headers.get("content-length")
    if declared and int(declared) > config.MAX_UPLOAD_BYTES:
        raise HTTPException(413, f"File exceeds the {config.MAX_UPLOAD_BYTES // 1024 ** 2} MB limit")

    config.UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    job_id = secrets.token_urlsafe(18)
    dest = config.UPLOAD_DIR / f"{job_id}.upload"
    size = 0
    try:
        with open(dest, "wb") as out:
            async for chunk in request.stream():
                size += len(chunk)
                if size > config.MAX_UPLOAD_BYTES:
                    raise HTTPException(
                        413, f"File exceeds the {config.MAX_UPLOAD_BYTES // 1024 ** 2} MB limit"
                    )
                out.write(chunk)
    except BaseException:
        dest.unlink(missing_ok=True)
        raise
    if size == 0:
        dest.unlink(missing_ok=True)
        raise HTTPException(400, "Empty upload")

    safe_name = Path(filename).name[:200] or "upload"
    now = time.time()
    with _lock:
        _jobs[job_id] = {"state": "queued", "progress": 0.0, "message": "Queued",
                         "created": now, "updated": now, "filename": safe_name, "size": size}
    _executor.submit(_run, job_id, dest, safe_name, size, None if assembly == "auto" else assembly)
    return JSONResponse({"job_id": job_id, "size": size}, status_code=202)


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
