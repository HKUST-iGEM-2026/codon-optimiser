"""Round-trip: re-optimising native genes of this organism.

Native ribosomal proteins are already close to the tool's target, so they must
move far less than a poorly adapted native gene. If the CAI weights were built
from the mutational background instead of the reference set, this inverts.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio import SeqIO

from app import engine as e
from app import organism as org

DATA = Path(__file__).parent / "data"
MAX_RIBOSOMAL_CHURN = 0.45
CHURN_MARGIN = 0.10


def _genes(path: Path):
    return [(r.id, str(r.seq).upper()) for r in SeqIO.parse(path, "fasta")]


RIBOSOMAL = _genes(DATA / "native_ribosomal.fna")
LOW_CAI = _genes(DATA / "native_low_cai.fna")


@pytest.mark.parametrize("locus,cds", RIBOSOMAL, ids=[locus for locus, _ in RIBOSOMAL])
def test_native_ribosomal_gene_round_trips(locus, cds):
    result = e.optimise(cds)
    metrics = result["metrics"]
    assert org.translate(result["sequence"]) == org.translate(cds)
    assert result["sequence"][:3] == cds[:3], "native initiator was not preserved"
    assert metrics["cai"] >= metrics["cai_of_input"], "re-optimising de-adapted a native gene"
    assert metrics["changed_codons_fraction"] < MAX_RIBOSOMAL_CHURN
    # the 5' convention: weights come from trimmed genes, so report both
    assert metrics["changed_codons_after_trim"] <= metrics["changed_codons"]


def test_ribosomal_genes_move_less_than_a_poorly_adapted_gene():
    churn = lambda cds: e.optimise(cds)["metrics"]["changed_codons_fraction"]  # noqa: E731
    ribosomal = [churn(cds) for _, cds in RIBOSOMAL]
    control = [churn(cds) for _, cds in LOW_CAI]
    assert max(ribosomal) < min(control) - CHURN_MARGIN


def test_reoptimising_our_own_output_is_a_fixed_point():
    """Feeding the result back in must not keep churning codons."""
    first = e.optimise(RIBOSOMAL[0][1])["sequence"]
    second = e.optimise(first)
    assert second["sequence"] == first
    assert second["metrics"]["changed_codons"] == 0
