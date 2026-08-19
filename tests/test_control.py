"""E. coli control: the only test that catches a silent hygiene-filter bug.

Runs the identical table-building code on GCF_000005845.2 and correlates CAI
against PaxDb protein abundance. If the hygiene filter, the reference set or the
CAI weights break, this correlation collapses while everything else still passes.

Needs the assembly (NCBI Datasets CLI, cached) and the PaxDb integrated dataset.
Both are cached under $ASSEMBLY_CACHE (default .work).
"""

from __future__ import annotations

import math
import os
import statistics
import urllib.request
from pathlib import Path

import pytest

from pipeline import gate, tables
from pipeline.fetch import ensure_assembly

CONTROL_ACCESSION = "GCF_000005845.2"
PAXDB_URL = "https://pax-db.org/downloads/latest/datasets/511145/511145-WHOLE_ORGANISM-integrated.txt"
CACHE = Path(os.environ.get("ASSEMBLY_CACHE", ".work"))
MIN_SPEARMAN = 0.30  # "clearly positive"; a broken filter lands near zero
MIN_S_INDEX = 0.30  # E. coli is a textbook translationally selected genome


def paxdb_abundances() -> dict[str, float]:
    """b-number -> integrated whole-organism abundance (ppm)."""
    path = CACHE / "paxdb_ecoli_integrated.txt"
    if not path.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        # pax-db.org answers 403 to urllib's default user agent
        request = urllib.request.Request(PAXDB_URL, headers={"User-Agent": "codon-optimiser/1.0"})
        with urllib.request.urlopen(request, timeout=60) as response:
            path.write_bytes(response.read())
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        _, string_id, abundance = line.split("\t")[:3]
        out[string_id.split(".", 1)[1]] = float(abundance)
    return out


@pytest.fixture(scope="module")
def control_tables():
    try:
        asm = ensure_assembly(CONTROL_ACCESSION)
    except Exception as exc:  # noqa: BLE001 - network/CLI absence is not a failure
        pytest.skip(f"control assembly unavailable: {exc}")
    records = tables.load_records(asm)
    facts = tables.facts(records)
    genes, stage, stops = tables.hygiene(records, facts["translation_table"])
    ref_idx = tables.reference_set(genes)
    ref_counts = tables.codon_counts([genes[i] for i in ref_idx])
    return {
        "facts": facts,
        "genes": genes,
        "stage": stage,
        "weights": tables.cai_weights(ref_counts, facts["translation_table"]),
        "ref_counts": ref_counts,
        "n_reference": len(ref_idx),
    }


def test_control_assembly_reads_as_expected(control_tables):
    facts = control_tables["facts"]
    assert facts["translation_table"] == 11
    assert 50.0 < facts["genomic_gc_percent"] < 51.5
    assert control_tables["stage"]["kept"] > 3000
    assert control_tables["n_reference"] >= 40


def test_cai_correlates_with_paxdb_abundance(control_tables):
    try:
        abundance = paxdb_abundances()
    except OSError as exc:
        pytest.skip(f"PaxDb unavailable: {exc}")

    weights = control_tables["weights"]
    pairs = [
        (tables.cai(g["codons"], weights), math.log10(abundance[g["locus_tag"]]))
        for g in control_tables["genes"]
        if g["locus_tag"] in abundance and abundance[g["locus_tag"]] > 0
    ]
    assert len(pairs) > 2000, "join with PaxDb b-numbers broke"
    cai_values, log_abundance = zip(*pairs)
    rho = statistics.correlation(cai_values, log_abundance, method="ranked")
    assert rho > MIN_SPEARMAN, f"CAI vs PaxDb Spearman rho = {rho:.3f}"


def test_control_s_index_shows_selection(control_tables):
    s = gate.s_index(control_tables["ref_counts"], tables.gc3(control_tables["genes"]))
    assert s["value"] > MIN_S_INDEX
    assert all(v > 0 for v in s["per_amino_acid"].values())
