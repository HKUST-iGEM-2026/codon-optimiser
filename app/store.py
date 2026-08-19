"""Permalink storage. SQLite, and only for permalinks.

The only place user input touches disk, and only when the caller explicitly asks
for a permalink: a result is reproducible from (input, parameters, seed, engine
version), so that tuple is all a row holds. Nothing here is read for reference
data, and no request that did not ask for a permalink is ever written.
"""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import tempfile
from pathlib import Path

DEFAULT_DB = Path("data/permalinks.sqlite")


def _default_db() -> Path:
    """`data/` beside the tables, or the system temp dir when the tree is read-only.

    Serverless runtimes mount the project read-only and give you a writable temp dir
    instead — Vercel functions, Lambda, Cloud Run with a read-only root. Detecting that
    rather than naming a platform keeps one rule for all of them, and the alternative
    (a hardcoded /tmp) would scatter permalink databases on a normal machine.
    """
    probe = DEFAULT_DB.parent if DEFAULT_DB.parent.is_dir() else Path(".")
    if os.access(probe, os.W_OK):
        return DEFAULT_DB
    return Path(tempfile.gettempdir()) / "permalinks.sqlite"


DB_PATH = Path(os.environ["PERMALINK_DB"]) if os.environ.get("PERMALINK_DB") else _default_db()

SCHEMA = """
CREATE TABLE IF NOT EXISTS permalink (
  id TEXT PRIMARY KEY,
  payload TEXT NOT NULL,
  engine_version TEXT NOT NULL,
  table_version TEXT NOT NULL
)
"""


def _connect() -> sqlite3.Connection:
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    conn.execute(SCHEMA)
    return conn


def save(payload: dict, engine_version: str, table_version: str) -> str:
    """Content-addressed: the same request always yields the same permalink."""
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    pid = hashlib.sha256(f"{engine_version}\0{table_version}\0{blob}".encode()).hexdigest()[:16]
    with _connect() as conn:
        conn.execute(
            "INSERT OR IGNORE INTO permalink (id, payload, engine_version, table_version) VALUES (?,?,?,?)",
            (pid, blob, engine_version, table_version),
        )
    return pid


def load(pid: str) -> dict | None:
    with _connect() as conn:
        row = conn.execute("SELECT payload FROM permalink WHERE id = ?", (pid,)).fetchone()
    return json.loads(row[0]) if row else None


def delete(pid: str) -> bool:
    with _connect() as conn:
        return conn.execute("DELETE FROM permalink WHERE id = ?", (pid,)).rowcount > 0
