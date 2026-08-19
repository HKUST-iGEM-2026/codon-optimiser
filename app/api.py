"""FastAPI surface. /api/v1 numerical behaviour is frozen; v2 gets its own module.

Sequences live in memory for the duration of the request. They are not logged,
not written to disk (permalinks store the request, on explicit demand only, see
store.py), and never sent anywhere. Error responses carry field names and
coordinates, never sequence content.
"""

from __future__ import annotations

import logging
import os
import sqlite3
import time
from pathlib import Path
from typing import Literal

from fastapi import FastAPI, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel, ConfigDict, Field

from . import engine, liabilities, store
from . import organism as org

HERE = Path(__file__).parent
templates = Jinja2Templates(directory=HERE / "templates")

# The only claim this tool is allowed to make.
CLAIM = (
    "Output matches this organism's codon statistics and removes the named sequence "
    "liabilities. It is not a prediction of expression yield."
)

# Rate limiting is on CPU-seconds per IP, not request count: one 3kb protein costs
# thousands of times more than a health check.
CPU_BUDGET_SECONDS = 20.0
CPU_WINDOW_SECONDS = 60.0
# ponytail: in-process dict, one container. Ceiling: it resets on restart and does
# not span replicas. Upgrade path = the same counters in the permalink SQLite file.
_cpu_usage: dict[str, list[float]] = {}
# Set to the proxy's client-IP header (Cloudflare sends CF-Connecting-IP) only when
# the origin cannot be reached except through that proxy. Empty, the default, means
# "trust nothing but the socket". See the deploy section of README.md.
CLIENT_IP_HEADER = os.environ.get("CLIENT_IP_HEADER", "")

log = logging.getLogger("codon-optimiser")

app = FastAPI(
    title="Codon optimiser",
    description=f"Codon optimiser for {org.ORGANISM['organism']['organism']} ({org.ACCESSION})",
    version=engine.ENGINE_VERSION,
)
app.mount("/static", StaticFiles(directory=HERE / "static"), name="static")


class OptimiseRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sequence: str = Field(min_length=1, max_length=4 * engine.INPUT_CAP_NT)
    mode: Literal["maximise", "harmonise"] | None = None
    seed: int = engine.DEFAULT_SEED
    avoid: list[str] = Field(default_factory=list, max_length=32)
    upstream: str = Field(default="", max_length=500)
    downstream: str = Field(default="", max_length=500)
    locked: list[tuple[int, int]] = Field(default_factory=list, max_length=64)
    permalink: bool = False

    def params(self) -> engine.Params:
        return engine.Params(
            mode=self.mode or org.DEFAULT_MODE,
            seed=self.seed,
            avoid=tuple(self.avoid),
            upstream=self.upstream.upper(),
            downstream=self.downstream.upper(),
            locked=tuple((a, b) for a, b in self.locked),
        )

    def stored(self) -> dict:
        payload = self.model_dump()
        payload.pop("permalink")
        return payload


# --- middleware and error handling ----------------------------------------
def _client_ip(request: Request) -> str:
    """The caller the CPU budget is charged to.

    Behind a proxy the socket peer is the proxy, so every caller would share one
    budget and any one of them could exhaust it for everybody. The forwarded header
    fixes that, but it is caller-controlled, so it is trusted only when the
    deployment names it — and it should only be named when the origin refuses
    traffic that did not come through the proxy. On an origin that is reachable
    directly, trusting the header turns the rate limit into a free reset.
    """
    if CLIENT_IP_HEADER:
        forwarded = request.headers.get(CLIENT_IP_HEADER)
        if forwarded:
            return forwarded.split(",")[0].strip()
    return request.client.host if request.client else "unknown"


def _solve(request: Request, req: OptimiseRequest) -> dict:
    """Optimise, and bill the caller for the CPU their own solve burns.

    thread_time() has to be sampled *inside* this function: FastAPI runs sync
    endpoints in anyio's worker threadpool, so the same two calls placed in the
    HTTP middleware measure the event-loop thread instead and see well under 1%
    of the real cost — which left the rate limit unreachable in practice.
    """
    started = time.thread_time()
    try:
        return engine.optimise(req.sequence, req.params())
    finally:
        ip = _client_ip(request)
        window_start, spent = _cpu_usage.get(ip, [time.monotonic(), 0.0])
        _cpu_usage[ip] = [window_start, spent + max(0.0, time.thread_time() - started)]


@app.middleware("http")
async def cpu_rate_limit(request: Request, call_next):
    """Spend the budget, then refuse. Charging happens in _solve, on the thread
    that does the work; this only opens or closes the door."""
    ip = _client_ip(request)
    now = time.monotonic()
    window_start, spent = _cpu_usage.get(ip, [now, 0.0])
    if now - window_start > CPU_WINDOW_SECONDS:
        _cpu_usage[ip] = [now, 0.0]
    elif spent >= CPU_BUDGET_SECONDS:
        retry = int(CPU_WINDOW_SECONDS - (now - window_start)) + 1
        return JSONResponse(
            {"error": "cpu_budget_exhausted", "budget_cpu_seconds": CPU_BUDGET_SECONDS,
             "retry_after_seconds": retry},
            status_code=429,
            headers={"Retry-After": str(retry)},
        )
    return await call_next(request)


@app.exception_handler(RequestValidationError)
async def validation_error(request: Request, exc: RequestValidationError):
    """Default handler echoes the offending input; that input is a sequence."""
    return JSONResponse(
        {"error": "invalid_request", "fields": [".".join(str(p) for p in e["loc"]) for e in exc.errors()]},
        status_code=422,
    )


@app.exception_handler(engine.InvalidInput)
async def invalid_input(request: Request, exc: engine.InvalidInput):
    return JSONResponse({"error": "invalid_sequence", "detail": str(exc)}, status_code=400)


@app.exception_handler(engine.Infeasible)
async def infeasible(request: Request, exc: engine.Infeasible):
    """Hard constraints are never violated: no sequence is returned here, only the
    conflicting constraint pair with coordinates."""
    return JSONResponse(
        {"error": "infeasible", "detail": str(exc), "conflict": exc.conflict}, status_code=422
    )


@app.exception_handler(engine.SolverGaveUp)
async def solver_gave_up(request: Request, exc: engine.SolverGaveUp):
    """Distinct from infeasible: the search ran out of attempts, so claiming the
    design is impossible would be a lie. Retryable with a different seed."""
    return JSONResponse(
        {"error": "solver_gave_up", "detail": str(exc), "conflict": exc.conflict,
         "hint": "retry with a different seed, or relax the avoided motifs"},
        status_code=503,
    )


@app.exception_handler(engine.InvariantViolated)
async def invariant_violated(request: Request, exc: engine.InvariantViolated):
    log.error("invariant violated on %s: %s", request.url.path, exc)  # no body, ever
    return JSONResponse({"error": "internal_invariant_failed", "detail": str(exc)}, status_code=500)


# --- routes ---------------------------------------------------------------
@app.get("/", response_class=HTMLResponse)
@app.get("/p/{permalink}", response_class=HTMLResponse)
def page(request: Request, permalink: str | None = None):
    """Server-rendered shell: organism facts and the four gate numbers come from the
    tables at import time, so the page needs no client round trip to be readable."""
    return templates.TemplateResponse(
        request,
        "index.html",
        {
            "org": org.ORGANISM,
            "gate": org.ORGANISM["gate"],
            "limits": {"input_cap_nt": engine.INPUT_CAP_NT},
            "default_seed": engine.DEFAULT_SEED,
            "five_prime_window": liabilities.FIVE_PRIME_WINDOW,
            "claim": CLAIM,
            "permalink": permalink,
        },
    )


@app.get("/healthz")
def healthz() -> dict:
    return {"status": "ok", "engine": engine.ENGINE_VERSION, "table": org.TABLE_VERSION}


@app.get("/api/v1/organism")
def get_organism() -> dict:
    """Everything the UI needs to render the header, the gate and the codon track."""
    o = org.ORGANISM
    return {
        "assembly_accession": org.ACCESSION,
        "organism": o["organism"],
        "table_version": org.TABLE_VERSION,
        "pipeline_version": o["pipeline_version"],
        "engine_version": engine.ENGINE_VERSION,
        "genes_used": o["genes_used"],
        "hygiene": o["hygiene"],
        "reference_set": {k: o["reference_set"][k] for k in ("n", "n_ribosomal", "method")},
        "gate": o["gate"],
        "gate_summary": org.gate_summary(),
        "default_mode": org.DEFAULT_MODE,
        "gc": o["gc"],
        "anti_sd": o["anti_sd"],
        "cai_weights": org.CAI_WEIGHTS,
        "relative_frequencies": org.NATIVE_FREQ,
        "reference_frequencies": org.REFERENCE_FREQ,
        "codon_families": org.FAMILIES,
        "stop_frequencies": org.STOP_FREQ,
        "limits": {
            "input_cap_nt": engine.INPUT_CAP_NT,
            "cpu_budget_seconds": CPU_BUDGET_SECONDS,
            "cpu_window_seconds": CPU_WINDOW_SECONDS,
        },
        "claim": CLAIM,
    }


@app.post("/api/v1/optimise")
def post_optimise(request: Request, req: OptimiseRequest) -> dict:
    result = _solve(request, req)
    if req.permalink:
        try:
            result["permalink"] = store.save(req.stored(), engine.ENGINE_VERSION, org.TABLE_VERSION)
        except (sqlite3.Error, OSError) as exc:
            # A read-only or full disk must not discard a correct, already-paid-for
            # solve. Report the storage failure and return the result anyway.
            log.error("permalink store failed: %s", exc)  # never the request body
            result["warnings"].append(
                "permalink could not be stored, so this result has no shareable link; "
                "the sequence and every number above are unaffected"
            )
    return result


@app.get("/api/v1/permalink/{pid}")
def get_permalink(request: Request, pid: str) -> dict:
    """Recomputed, not replayed: the row holds the request, the engine is deterministic."""
    payload = store.load(pid)
    if payload is None:
        raise HTTPException(status_code=404, detail="unknown permalink")
    req = OptimiseRequest(**payload)
    result = _solve(request, req)
    result["permalink"] = pid
    return result


@app.post("/api/v1/export.fasta", response_class=PlainTextResponse)
def export_fasta(request: Request, req: OptimiseRequest) -> str:
    """FASTA export. The manifest travels in the header line, as required."""
    result = _solve(request, req)
    m = result["manifest"]
    p = m["parameters"]
    header = (
        f">codon_optimised engine={m['engine_version']} api={m['api_version']} "
        f"table={m['table_version']} assembly={m['assembly_accession']} "
        f"transl_table={m['translation_table']} mode={p['mode']} seed={p['seed']} "
        f"avoid={','.join(p['avoid']) or 'none'} gc_global={p['gc_global_band_percent'][0]}-"
        f"{p['gc_global_band_percent'][1]} gc_window={p['gc_window_band_percent'][0]}-"
        f"{p['gc_window_band_percent'][1]}@{p['gc_window_nt']}nt cai={result['metrics']['cai']} "
        f"dnachisel={m['dnachisel_version']}"
    )
    seq = result["sequence"]
    body = "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))
    return f"{header}\n{body}\n"
