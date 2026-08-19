"""API contracts: privacy of error responses, permalink recomputation, limits."""

from __future__ import annotations

import hashlib
import json
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

from fastapi.testclient import TestClient

from app import api, engine, store
from app import organism as org

client = TestClient(api.app)
PROTEIN = "MSKIYAIAKELNVSKQTVSRWLNDESLSKEQLNLLAKELNVSPSELLGEEIKNNKVAVIYPNSNNPFWVEMKRG"
SECRET = PROTEIN + "WQWQWQ"


def test_organism_endpoint_serves_the_four_gate_numbers():
    body = client.get("/api/v1/organism").json()
    summary = body["gate_summary"]
    for key in ("ribosomal_cai_median_percentile", "s_index", "ecai", "ribosomal_self_consistent_overlap"):
        assert summary[key] is not None
    assert body["table_version"] == org.TABLE_VERSION
    assert "not a prediction of expression yield" in body["claim"]


def test_optimise_returns_manifest_and_per_codon_explanation():
    body = client.post("/api/v1/optimise", json={"sequence": PROTEIN}).json()
    assert body["manifest"]["table_version"] == org.TABLE_VERSION
    assert body["manifest"]["parameters"]["seed"] == engine.DEFAULT_SEED
    assert len(body["codons"]) == len(PROTEIN) + 1
    assert body["metrics"]["cai"] > 0.5


def test_validation_errors_never_echo_the_sequence():
    """FastAPI's default validation handler would return the input verbatim."""
    response = client.post("/api/v1/optimise", json={"sequence": SECRET, "mode": "nonsense"})
    assert response.status_code == 422
    assert SECRET not in response.text
    assert response.json()["fields"] == ["body.mode"]


def test_invalid_residues_are_reported_without_the_sequence():
    response = client.post("/api/v1/optimise", json={"sequence": SECRET + "XXX"})
    assert response.status_code == 400
    assert SECRET not in response.text


def test_unusable_avoid_motifs_are_rejected_before_the_solver_sees_them():
    """An empty motif matches everywhere: as a constraint it is unsatisfiable, and
    the solver spends minutes of CPU on its random search before saying so."""
    for motif in ("", "GAATTZ"):
        started = time.monotonic()
        response = client.post("/api/v1/optimise", json={"sequence": PROTEIN, "avoid": [motif]})
        assert response.status_code == 400, motif
        assert response.json()["error"] == "invalid_sequence"
        assert time.monotonic() - started < 5


def test_infeasible_request_returns_coordinates_and_no_sequence():
    response = client.post("/api/v1/optimise", json={"sequence": "MEEEEEEEEEK", "avoid": ["GAA", "GAG"]})
    assert response.status_code in (422, 503)
    body = response.json()
    assert "sequence" not in body
    assert body["conflict"][0]["end"] > body["conflict"][0]["start"]


def test_permalink_recomputes_the_same_result():
    first = client.post("/api/v1/optimise", json={"sequence": PROTEIN, "permalink": True}).json()
    again = client.get(f"/api/v1/permalink/{first['permalink']}").json()
    assert again["sequence"] == first["sequence"]
    assert again["manifest"] == first["manifest"]
    assert client.get("/api/v1/permalink/deadbeefdeadbeef").status_code == 404


def test_nothing_is_stored_without_an_explicit_permalink_request():
    """Permalinks are the one route that persists input, and only on demand."""
    sequence = PROTEIN + "AAK"
    payload = {"sequence": sequence, "mode": None, "seed": engine.DEFAULT_SEED,
               "avoid": [], "upstream": "", "downstream": "", "locked": []}
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    would_be = hashlib.sha256(
        f"{engine.ENGINE_VERSION}\0{org.TABLE_VERSION}\0{blob}".encode()
    ).hexdigest()[:16]
    store.delete(would_be)

    body = client.post("/api/v1/optimise", json={"sequence": sequence}).json()
    assert "permalink" not in body
    assert store.load(would_be) is None

    stored = client.post("/api/v1/optimise", json={"sequence": sequence, "permalink": True}).json()
    assert stored["permalink"] == would_be
    assert store.delete(would_be)


def test_export_fasta_carries_the_manifest_in_the_header():
    text = client.post("/api/v1/export.fasta", json={"sequence": PROTEIN}).text
    header = text.splitlines()[0]
    assert header.startswith(">codon_optimised")
    for field in ("engine=", "table=", "assembly=", "transl_table=", "mode=", "seed=", "cai="):
        assert field in header


def test_cpu_budget_returns_429_with_retry_after():
    api._cpu_usage["testclient"] = [time.monotonic(), api.CPU_BUDGET_SECONDS + 1]
    try:
        response = client.post("/api/v1/optimise", json={"sequence": PROTEIN})
        assert response.status_code == 429
        assert int(response.headers["Retry-After"]) > 0
    finally:
        api._cpu_usage.pop("testclient", None)


def test_cpu_budget_is_charged_for_the_solve_not_the_event_loop():
    """The limiter is only real if it sees the solve's CPU.

    Sampling thread_time() in the HTTP middleware measures the event-loop thread
    while the solve runs in anyio's threadpool: that billed ~0.3% of the true cost,
    so a 20 CPU-second budget was unreachable.
    """
    api._cpu_usage.pop("testclient", None)
    before = time.process_time()
    client.post("/api/v1/optimise", json={"sequence": PROTEIN})
    burned = time.process_time() - before
    billed = api._cpu_usage["testclient"][1]
    assert billed > 0.5 * burned, f"billed {billed:.4f}s of {burned:.4f}s CPU"


def test_page_renders_gate_numbers_and_config_for_the_client():
    """The shell is server-rendered: the gate must be readable with JS disabled."""
    html = client.get("/").text
    gate = org.ORGANISM["gate"]
    assert str(gate["s_index"]["value"]) in html
    assert str(gate["ecai"]["value"]) in html
    assert org.ORGANISM["organism"]["organism"] in html
    assert f'data-default-mode="{org.DEFAULT_MODE}"' in html
    assert f'data-cap="{engine.INPUT_CAP_NT}"' in html
    assert "not a prediction of expression yield" in html


def test_static_assets_and_permalink_page_are_served():
    assert client.get("/static/app.css").status_code == 200
    assert "codon" in client.get("/static/app.js").text
    page = client.get("/p/deadbeefdeadbeef")
    assert page.status_code == 200
    assert 'data-permalink="deadbeefdeadbeef"' in page.text


def test_a_failed_permalink_write_warns_but_keeps_the_result(monkeypatch):
    """The solve is the expensive part and it already succeeded: a read-only or full
    disk must degrade to 'no shareable link', not to a discarded 500."""
    def boom(*_args):
        raise sqlite3.OperationalError("attempt to write a readonly database")

    monkeypatch.setattr(api.store, "save", boom)
    response = client.post("/api/v1/optimise", json={"sequence": PROTEIN, "permalink": True})
    assert response.status_code == 200
    body = response.json()
    assert "permalink" not in body
    assert any("no shareable link" in w for w in body["warnings"])
    assert body["sequence"]


def test_missing_codon_tables_name_the_file_and_the_fix():
    """The default path is relative, so starting outside the repo root is the most
    likely deployment failure; a bare FileNotFoundError does not explain it."""
    repo_root = Path(__file__).resolve().parents[1]
    proc = subprocess.run(
        [sys.executable, "-c", "import app.organism"],
        capture_output=True, text=True, cwd=repo_root,
        env={"ORGANISM_JSON": "/nonexistent/organism.json", "PATH": ""},
    )
    assert proc.returncode != 0
    assert "cannot load the codon tables" in proc.stderr
    assert "/nonexistent/organism.json" in proc.stderr
    assert "ORGANISM_JSON" in proc.stderr
