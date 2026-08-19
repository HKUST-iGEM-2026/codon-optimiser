"""Build data/organism.json from a downloaded assembly. Run by hand, once.

    python -m pipeline.build GCF_000409755.1

Writes data/organism.json and stops. Not called at runtime.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from Bio.Data import CodonTable

from . import gate, tables
from .fetch import accession_of, ensure_assembly

PIPELINE_VERSION = "1.0.0"
# sscu's own documentation: S > 0.3 means C-ending codons are used significantly
# more than U-ending ones, i.e. translational selection is acting. Verified our
# s_index against their L. kunkeei example (1.5419 vs published 1.54).
S_CLEARLY_POSITIVE = 0.3


def build(asm_dir: str | Path) -> dict:
    """Assembly directory -> organism table dict. Identical code for the control."""
    records = tables.load_records(asm_dir)
    org = tables.facts(records)
    table_id = org["translation_table"]

    genes, stage, stop_counts = tables.hygiene(records, table_id)
    genome_counts = tables.codon_counts(genes)
    weights_genome = tables.cai_weights(genome_counts, table_id)

    ref_idx = tables.reference_set(genes)
    ref_genes = [genes[i] for i in ref_idx]
    ref_counts = tables.codon_counts(ref_genes)
    weights_ref = tables.cai_weights(ref_counts, table_id)

    ribo_idx = tables.ribosomal_only(genes)
    ref_codons = [c for g in ref_genes for c in g["codons"]]

    sc = gate.self_consistent_set(genes, table_id)
    rank = gate.ribosomal_rank(genes, ribo_idx, weights_ref, "reference-set CAI weights")
    rank["under_genome_weights"] = gate.ribosomal_rank(
        genes, ribo_idx, weights_genome, "genome-wide codon counts"
    )
    gate_numbers = {
        "ribosomal_cai_rank": rank,
        "s_index": gate.s_index(ref_counts, tables.gc3(genes)),
        "ecai": gate.ecai(ref_codons, weights_ref, genes, table_id),
        "overlap": gate.overlap(ribo_idx, sc),
    }
    selected = gate_numbers["s_index"]["value"] > S_CLEARLY_POSITIVE and gate_numbers["ecai"]["passes"]

    per_cds_gc = [g["gc_percent"] for g in genes]
    table = CodonTable.unambiguous_dna_by_id[table_id]
    out = {
        "schema": 1,
        "pipeline_version": PIPELINE_VERSION,
        "assembly_accession": accession_of(asm_dir),
        "organism": org,
        "hygiene": dict(sorted(stage.items())),
        "genes_used": len(genes),
        "codon_counts": {c: genome_counts[c] for c in sorted(genome_counts)},
        "stop_counts": dict(sorted(stop_counts.items())),
        "relative_frequencies": {
            c: round(v, 6) for c, v in sorted(tables.relative_frequencies(genome_counts, table_id).items())
        },
        "reference_frequencies": {
            c: round(v, 6) for c, v in sorted(tables.relative_frequencies(ref_counts, table_id).items())
        },
        "stop_frequencies": {
            c: round(n / sum(stop_counts.values()), 6) for c, n in sorted(stop_counts.items())
        },
        "cai_weights": {c: round(v, 6) for c, v in sorted(weights_ref.items())},
        "cai_weights_genome": {c: round(v, 6) for c, v in sorted(weights_genome.items())},
        "reference_set": {
            "n": len(ref_idx),
            "method": "product ~ 'ribosomal protein' (no methyl/hydroxy) + gene in tuf/tufA/tufB/tsf/fusA",
            "n_ribosomal": len(ribo_idx),
            "loci": [genes[i]["locus_tag"] for i in ref_idx],
            "genes": [genes[i]["gene"] or genes[i]["product"] for i in ref_idx],
        },
        "gate": gate_numbers,
        "default_mode": "maximise" if selected else "harmonise",
        "gc": {
            "genomic_percent": org["genomic_gc_percent"],
            "cds_mean_percent": round(sum(per_cds_gc) / len(per_cds_gc), 3),
            "native_band_percent": [
                round(tables.percentile(per_cds_gc, 5), 3),
                round(tables.percentile(per_cds_gc, 95), 3),
            ],
            "manufacturable_band_percent": [25.0, 75.0],
            "sliding_window_nt": 50,
        },
        "anti_sd": tables.anti_sd(records),
        "codon_families": tables.families(table_id),
        "start_codons": list(tables.START_CODONS),
        "stop_codons": sorted(table.stop_codons),
    }
    # table_version pins the numbers the app serves; any table change changes it.
    core = json.dumps(
        {
            k: out[k]
            for k in (
                "codon_counts", "cai_weights", "relative_frequencies", "reference_frequencies",
                "stop_frequencies", "gate", "default_mode",
            )
        },
        sort_keys=True,
    )
    out["table_version"] = f"{PIPELINE_VERSION}+{hashlib.sha256(core.encode()).hexdigest()[:12]}"
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("accession_or_dir")
    ap.add_argument("-o", "--out", default="data/organism.json")
    args = ap.parse_args()

    src = Path(args.accession_or_dir)
    asm_dir = src if src.is_dir() else ensure_assembly(args.accession_or_dir)
    out = build(asm_dir)

    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")

    g = out["gate"]
    print(f"{out['organism']['organism']}  table {out['organism']['translation_table']}  "
          f"GC {out['organism']['genomic_gc_percent']}%  CDS {out['organism']['cds_count']}")
    print(f"genes used {out['genes_used']}  reference set {out['reference_set']['n']}")
    print(f"ribosomal CAI median rank {g['ribosomal_cai_rank']['median_rank']} "
          f"({g['ribosomal_cai_rank']['median_percentile']}% of {g['ribosomal_cai_rank']['n_ranked_genes']})")
    print(f"S index {g['s_index']['value']}")
    print(f"eCAI {g['ecai']['value']}  observed {g['ecai']['observed_reference_cai']}  passes {g['ecai']['passes']}")
    print(f"overlap {g['overlap']['shared']}/{g['overlap']['ribosomal_set_size']} ribosomal in "
          f"self-consistent set of {g['overlap']['self_consistent_set_size']}")
    print(f"default mode {out['default_mode']}  table_version {out['table_version']}")
    print(f"wrote {dest}")
    for w in out["organism"]["warnings"]:
        print(w)


if __name__ == "__main__":
    main()
