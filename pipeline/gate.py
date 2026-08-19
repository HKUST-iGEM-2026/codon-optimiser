"""Step 2 gate metrics: does this organism show translational selection?

Four numbers decide the tool's default mode (maximise vs harmonise).
Offline only.
"""

from __future__ import annotations

import math
import random
from collections import Counter

from Bio.Data import CodonTable

from .tables import (
    REF_SET_MIN,
    cai,
    cai_weights,
    codon_counts,
    families,
    gc3,
)

GATE_SEED = 20240517  # documented, non-null default RNG seed
ECAI_REPLICATES = 1000
SELF_CONSISTENT_FRACTION = 0.05


def s_index(ref_counts: Counter, genomic_gc3: float) -> dict:
    """Sharp et al. 2005 S index, exactly as Bioconductor sscu::s_index computes it.

    Weighted mean over the four 2-fold NNT/NNC families of
    log[ (f_C / (1 - f_C)) * (1 - gc3) / gc3 ].
    """
    pairs = {"Tyr": ("TAC", "TAT"), "Asn": ("AAC", "AAT"), "Ile": ("ATC", "ATT"), "Phe": ("TTC", "TTT")}
    per_aa, num, den = {}, 0.0, 0
    for aa, (c_end, t_end) in pairs.items():
        nc, nt = ref_counts[c_end], ref_counts[t_end]
        total = nc + nt
        if total == 0 or nc == 0 or nc == total:
            raise ValueError(f"S index undefined: {aa} counts {nc}/{nt} in the reference set")
        f = nc / total
        s = math.log(f * (1 - genomic_gc3) / genomic_gc3 / (1 - f))
        per_aa[aa] = round(s, 4)
        num += s * total
        den += total
    return {"value": round(num / den, 4), "per_amino_acid": per_aa, "genomic_gc3": round(genomic_gc3, 4)}


def _null_codon_model(genes: list[dict], table_id: int) -> dict[str, list[tuple[str, float]]]:
    """Per-family codon probabilities implied by nucleotide composition alone.

    Position-specific base frequencies of the whole hygiene-passing set; a codon's
    probability inside its family is the product of its three base frequencies,
    renormalised. Matches G+C (including GC3) without carrying any synonymous
    preference, which is what a composition-matched null means.
    """
    pos = [Counter(), Counter(), Counter()]
    for g in genes:
        for c in g["codons"]:
            for i in range(3):
                pos[i][c[i]] += 1
    p = [{b: n / sum(cnt.values()) for b, n in cnt.items()} for cnt in pos]
    model = {}
    for aa, codons in families(table_id).items():
        raw = [(c, p[0].get(c[0], 0) * p[1].get(c[1], 0) * p[2].get(c[2], 0)) for c in codons]
        total = sum(w for _, w in raw) or 1.0
        model[aa] = [(c, w / total) for c, w in raw]
    return model


def ecai(
    ref_codons: list[str],
    weights: dict[str, float],
    genes: list[dict],
    table_id: int,
    replicates: int = ECAI_REPLICATES,
    seed: int = GATE_SEED,
) -> dict:
    """95th percentile CAI of composition-matched random sequences (Puigbo E-CAI).

    Random sequences keep the reference set's amino acid composition and the
    genome's nucleotide composition. Observed CAI above the 95th percentile means
    the reference set's codon bias cannot be explained by composition alone.
    """
    forward = CodonTable.unambiguous_dna_by_id[table_id].forward_table
    model = _null_codon_model(genes, table_id)
    aas = [forward[c] for c in ref_codons if c in forward]
    rng = random.Random(seed)
    scores = []
    for _ in range(replicates):
        draw = []
        for aa in aas:
            codons, probs = zip(*model[aa])
            draw.append(rng.choices(codons, weights=probs, k=1)[0])
        scores.append(cai(draw, weights))
    scores.sort()
    threshold = scores[int(0.95 * (replicates - 1))]
    observed = cai(ref_codons, weights)
    return {
        "value": round(threshold, 4),
        "observed_reference_cai": round(observed, 4),
        "passes": observed > threshold,
        "replicates": replicates,
        "seed": seed,
        "null_mean": round(sum(scores) / len(scores), 4),
    }


def self_consistent_set(genes: list[dict], table_id: int, max_iter: int = 20) -> dict:
    """Iteratively re-derive CAI weights from the top-CAI genes until the set is stable."""
    n = max(REF_SET_MIN, int(len(genes) * SELF_CONSISTENT_FRACTION))
    weights = cai_weights(codon_counts(genes), table_id)
    current: set[int] = set()
    iterations = 0
    for iterations in range(1, max_iter + 1):
        scores = [cai(g["codons"], weights) for g in genes]
        order = sorted(range(len(genes)), key=lambda i: (-scores[i], i))
        nxt = set(order[:n])
        if nxt == current:
            break
        current = nxt
        weights = cai_weights(codon_counts([genes[i] for i in order[:n]]), table_id)
    return {"members": sorted(current), "size": n, "iterations": iterations}


def ribosomal_rank(
    genes: list[dict], ribo_idx: list[int], weights: dict[str, float], weights_from: str
) -> dict:
    """Where the ribosomal proteins land when every gene is ranked by CAI.

    Run it with the reference-set weights (the CAI the tool ships) and with the
    genome-wide weights: in a mutation-biased genome the two disagree, and the
    disagreement is the point.
    """
    scores = [cai(g["codons"], weights) for g in genes]
    order = sorted(range(len(genes)), key=lambda i: (-scores[i], i))
    rank = {gene_i: r + 1 for r, gene_i in enumerate(order)}
    ranks = sorted(rank[i] for i in ribo_idx)
    median = ranks[len(ranks) // 2] if len(ranks) % 2 else (ranks[len(ranks) // 2 - 1] + ranks[len(ranks) // 2]) / 2
    return {
        "median_rank": median,
        "median_percentile": round(100 * median / len(genes), 2),
        "best_rank": ranks[0],
        "worst_rank": ranks[-1],
        "n_ranked_genes": len(genes),
        "n_ribosomal": len(ranks),
        "weights_from": weights_from,
    }


def overlap(ribo_idx: list[int], sc: dict) -> dict:
    members = set(sc["members"])
    shared = sorted(set(ribo_idx) & members)
    return {
        "shared": len(shared),
        "ribosomal_set_size": len(ribo_idx),
        "self_consistent_set_size": len(members),
        "fraction_of_ribosomal_set": round(len(shared) / len(ribo_idx), 4),
        "iterations": sc["iterations"],
    }
