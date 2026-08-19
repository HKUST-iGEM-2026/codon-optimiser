"""The organism tables, read once at process start into a module-level dict.

No database, no pipeline import, no network. If data/organism.json is missing or
stale, the process must not start.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

ORGANISM_JSON = Path(os.environ.get("ORGANISM_JSON", "data/organism.json"))
SCHEMA = 1

try:
    ORGANISM: dict = json.loads(ORGANISM_JSON.read_text())
except (OSError, ValueError) as exc:
    # The likeliest deployment mistake by far: the default path is relative, so the
    # process must start from the repository root. Say which file and how to fix it,
    # rather than a bare FileNotFoundError on 'data/organism.json'.
    raise RuntimeError(
        f"cannot load the codon tables from {ORGANISM_JSON} "
        f"(resolved to {ORGANISM_JSON.absolute()}): {exc}. "
        "Start the service from the repository root, or set ORGANISM_JSON to the file."
    ) from exc
if ORGANISM.get("schema") != SCHEMA:
    raise RuntimeError(f"{ORGANISM_JSON}: schema {ORGANISM.get('schema')}, expected {SCHEMA}")

TABLE_ID: int = ORGANISM["organism"]["translation_table"]
TABLE_VERSION: str = ORGANISM["table_version"]
ACCESSION: str = ORGANISM["assembly_accession"]
DEFAULT_MODE: str = ORGANISM["default_mode"]

FAMILIES: dict[str, list[str]] = ORGANISM["codon_families"]  # aa -> codons, stops excluded
STOP_CODONS: list[str] = ORGANISM["stop_codons"]
CAI_WEIGHTS: dict[str, float] = ORGANISM["cai_weights"]  # reference-set derived
NATIVE_FREQ: dict[str, float] = ORGANISM["relative_frequencies"]  # genome-wide, harmonisation target
REFERENCE_FREQ: dict[str, float] = ORGANISM["reference_frequencies"]
STOP_FREQ: dict[str, float] = ORGANISM["stop_frequencies"]

CODON_TO_AA: dict[str, str] = {c: aa for aa, codons in FAMILIES.items() for c in codons}
CODON_TO_AA.update({c: "*" for c in STOP_CODONS})

NATIVE_GC_BAND: tuple[float, float] = tuple(ORGANISM["gc"]["native_band_percent"])
MANUFACTURABLE_GC_BAND: tuple[float, float] = tuple(ORGANISM["gc"]["manufacturable_band_percent"])
GC_WINDOW_NT: int = ORGANISM["gc"]["sliding_window_nt"]


def usage_table(frequencies: dict[str, float]) -> dict[str, dict[str, float]]:
    """DnaChisel codon_usage_table shape: {aa: {codon: relative frequency}}."""
    table = {aa: {c: frequencies[c] for c in codons} for aa, codons in FAMILIES.items()}
    table["*"] = dict(STOP_FREQ)
    return table


MAXIMISE_USAGE = usage_table(REFERENCE_FREQ)
HARMONISE_USAGE = usage_table(NATIVE_FREQ)


def translate(dna: str) -> str:
    """Translation under this organism's table only. Stops become '*'."""
    return "".join(CODON_TO_AA[dna[i : i + 3]] for i in range(0, len(dna), 3))


def gate_summary() -> dict:
    """The four Step-2 numbers, for /api/v1/organism and the UI."""
    g = ORGANISM["gate"]
    return {
        "ribosomal_cai_median_percentile": g["ribosomal_cai_rank"]["median_percentile"],
        "ribosomal_cai_median_rank": g["ribosomal_cai_rank"]["median_rank"],
        "s_index": g["s_index"]["value"],
        "ecai": g["ecai"]["value"],
        "ecai_observed": g["ecai"]["observed_reference_cai"],
        "ecai_passes": g["ecai"]["passes"],
        "ribosomal_self_consistent_overlap": g["overlap"]["fraction_of_ribosomal_set"],
        "ribosomal_self_consistent_shared": f"{g['overlap']['shared']}/{g['overlap']['ribosomal_set_size']}",
        "default_mode": DEFAULT_MODE,
    }
