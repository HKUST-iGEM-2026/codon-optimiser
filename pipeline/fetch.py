"""Assembly download via the NCBI Datasets CLI. Offline use and the control test."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import zipfile
from pathlib import Path

CACHE = Path(os.environ.get("ASSEMBLY_CACHE", ".work"))
INCLUDE = "cds,gbff,gff3,genome,seq-report"


def datasets_cli() -> str:
    cli = os.environ.get("DATASETS_CLI") or shutil.which("datasets") or ".tools/datasets"
    if not Path(cli).exists() and shutil.which(cli) is None:
        raise FileNotFoundError(
            "NCBI Datasets CLI not found. Install it, or set DATASETS_CLI. "
            "macOS/linux binary: https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/"
        )
    return cli


def ensure_assembly(accession: str, cache: Path = CACHE) -> Path:
    """Download <accession> (RefSeq) into the cache if absent; return its data dir."""
    out = cache / accession
    data = out / "ncbi_dataset" / "data" / accession
    if data.is_dir():
        return data
    cache.mkdir(parents=True, exist_ok=True)
    zip_path = cache / f"{accession}.zip"
    subprocess.run(
        [
            datasets_cli(), "download", "genome", "accession", accession,
            "--include", INCLUDE, "--assembly-source", "RefSeq",
            "--filename", str(zip_path),
        ],
        check=True,
    )
    with zipfile.ZipFile(zip_path) as z:
        z.extractall(out)
    if not data.is_dir():
        raise FileNotFoundError(f"{accession}: unexpected package layout under {out}")
    return data


def accession_of(asm_dir: str | Path) -> str:
    """Accession as reported by the package itself, not inferred from the path."""
    report = next(Path(asm_dir).parent.glob("assembly_data_report.jsonl"), None)
    if report is None:
        raise FileNotFoundError(f"no assembly_data_report.jsonl beside {asm_dir}")
    with report.open() as fh:
        return json.loads(fh.readline())["accession"]
