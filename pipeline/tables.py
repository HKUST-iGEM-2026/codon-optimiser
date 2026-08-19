"""Codon table building from a downloaded NCBI Datasets assembly directory.

Offline only. Nothing in app/ imports this module at runtime; the app reads
data/organism.json. The E. coli control test runs this exact code on
GCF_000005845.2, which is why the assembly directory is a parameter and
nothing else is.

Hygiene filter order and every constant here are fixed by CLAUDE.md.
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable

# --- fixed conventions (CLAUDE.md) -----------------------------------------
MIN_CODONS = 80  # measured BEFORE the trim -> effective floor 68 counted codons
TRIM_CODONS = 12
START_CODONS = ("ATG", "GTG", "TTG")
PSEUDOCOUNT = 0.5
CAI_EXCLUDED_AA = ("M", "W")  # plus stops, handled separately

REF_GENE_NAMES = frozenset({"tufa", "tufb", "tuf", "tsf", "fusa"})
# "ribosomal protein" and the EcoCyc-style "ribosomal subunit protein" (RefSeq
# names E. coli K-12 that way; the plain substring finds 5 genes there and the
# control test cannot run). Match is case-insensitive.
RIBO_PRODUCT = re.compile(r"ribosomal (subunit )?protein", re.I)
RIBO_PRODUCT_EXCLUDE = ("methyl", "hydroxy")
REF_SET_MIN = 20  # below this: stop, switch to HMM search

# ponytail: annotation-keyword prophage exclusion, not a prophage caller.
# Ceiling: misses unannotated prophage regions. Upgrade path = run PHASTER/
# geNomad once offline and feed the coordinates in here.
PROPHAGE_KEYWORDS = (
    "phage",
    "capsid",
    "terminase",
    "portal protein",
    "tail fiber",
    "tail fibre",
    "tail tape measure",
    "baseplate",
    "holin",
    "endolysin",
    "head-tail",
    "head morphogenesis",
)

# Lineages that famously do not use table 11 despite bacterial annotation.
SUSPECT_LINEAGE = (
    "mycoplasm",
    "spiroplasma",
    "ureaplasma",
    "entomoplasma",
    "mesoplasma",
    "gracilibacteria",
    "absconditabacteria",
    "sr1",
)

AMBIGUOUS = re.compile(r"[^ACGT]")


# --- assembly reading ------------------------------------------------------
def gbff_path(asm_dir: str | Path) -> Path:
    """The single genomic.gbff inside an `ncbi datasets download genome` tree."""
    hits = sorted(Path(asm_dir).rglob("genomic.gbff"))
    if len(hits) != 1:
        raise FileNotFoundError(f"expected exactly one genomic.gbff under {asm_dir}, got {hits}")
    return hits[0]


def load_records(asm_dir: str | Path) -> list:
    return list(SeqIO.parse(gbff_path(asm_dir), "genbank"))


def facts(records: list) -> dict:
    """Organism identity read from the assembly itself, never a secondary source."""
    src = records[0].features[0].qualifiers
    tables = Counter(
        f.qualifiers.get("transl_table", ["missing"])[0]
        for r in records
        for f in r.features
        if f.type == "CDS"
    )
    if len(tables) != 1:
        raise ValueError(f"mixed /transl_table values in assembly: {dict(tables)}")
    table_id = int(tables.most_common(1)[0][0])  # never hardcoded, never looped 1..33
    if table_id not in CodonTable.unambiguous_dna_by_id:
        raise ValueError(f"NCBI table {table_id} is not a known table id")

    gc = at = 0
    for r in records:
        s = str(r.seq).upper()
        gc += s.count("G") + s.count("C")
        at += s.count("A") + s.count("T")

    lineage = [t.lower() for t in records[0].annotations.get("taxonomy", [])]
    organism = records[0].annotations.get("organism", "")
    warnings = []
    hay = " ".join(lineage + [organism.lower()])
    if table_id == 11 and any(k in hay for k in SUSPECT_LINEAGE):
        warnings.append(
            f"HARD WARNING: lineage looks like a non-standard-table clade but /transl_table is 11 ({organism})"
        )

    taxon = next(
        (x.split(":", 1)[1] for x in src.get("db_xref", []) if x.startswith("taxon:")), None
    )
    return {
        "organism": organism,
        "strain": src.get("strain", [None])[0],
        "taxid": int(taxon) if taxon else None,
        "translation_table": table_id,
        "genomic_gc_percent": round(100 * gc / (gc + at), 4),
        "cds_count": sum(tables.values()),
        "sequence_count": len(records),
        "genome_length": sum(len(r.seq) for r in records),
        "lineage": records[0].annotations.get("taxonomy", []),
        "warnings": warnings,
    }


# --- hygiene filter --------------------------------------------------------
def _partial(feature) -> bool:
    loc = str(feature.location)
    return "<" in loc or ">" in loc or "partial" in feature.qualifiers


def _prophage(product: str, gene: str) -> bool:
    hay = f"{product} {gene}".lower()
    return any(k in hay for k in PROPHAGE_KEYWORDS)


def hygiene(records: list, table_id: int) -> tuple[list[dict], Counter, Counter]:
    """CLAUDE.md hygiene filter, in exact order. Returns (genes, stage counts, stop counts)."""
    table = CodonTable.unambiguous_dna_by_id[table_id]
    stops = set(table.stop_codons)
    # every stage always reported, zeros included: a stage that silently never
    # fires is exactly how a broken hygiene filter hides.
    stage = Counter(
        dict.fromkeys(
            (
                "cds_total", "drop_frame", "drop_partial_pseudo", "drop_ends",
                "drop_internal_stop", "drop_short", "drop_all_ambiguous",
                "drop_duplicate", "drop_plasmid", "drop_prophage",
                "skipped_ambiguous_codons", "kept",
            ),
            0,
        )
    )
    stop_counts: Counter = Counter()
    seen: set[str] = set()
    genes: list[dict] = []

    for rec in records:
        src = rec.features[0].qualifiers if rec.features else {}
        plasmid = "plasmid" in src or "plasmid" in rec.description.lower()
        for f in rec.features:
            if f.type != "CDS":
                continue
            stage["cds_total"] += 1
            nt = str(f.extract(rec.seq)).upper()
            product = f.qualifiers.get("product", [""])[0]
            gene = f.qualifiers.get("gene", [""])[0]

            # 1. not a multiple of three
            if len(nt) % 3:
                stage["drop_frame"] += 1
                continue
            # 2. partial records and pseudogenes
            if _partial(f) or "pseudo" in f.qualifiers or "pseudogene" in f.qualifiers:
                stage["drop_partial_pseudo"] += 1
                continue
            # 3. canonical start and valid stop under this organism's table
            codons = [nt[i : i + 3] for i in range(0, len(nt), 3)]
            if codons[0] not in START_CODONS or codons[-1] not in stops:
                stage["drop_ends"] += 1
                continue
            if any(c in stops for c in codons[1:-1]):
                stage["drop_internal_stop"] += 1
                continue
            # 4. >= 80 sense codons, measured BEFORE trimming
            sense = codons[:-1]
            if len(sense) < MIN_CODONS:
                stage["drop_short"] += 1
                continue
            # 5. then trim the first 12 codons
            kept = sense[TRIM_CODONS:]
            # 6. skip codons containing ambiguity codes, keep the gene
            clean = [c for c in kept if not AMBIGUOUS.search(c)]
            stage["skipped_ambiguous_codons"] += len(kept) - len(clean)
            if not clean:
                stage["drop_all_ambiguous"] += 1
                continue
            # 8. deduplicate identical sequences (tufA/tufB are often byte-identical)
            if nt in seen:
                stage["drop_duplicate"] += 1
                continue
            seen.add(nt)
            # 9. exclude plasmid and prophage genes
            if plasmid:
                stage["drop_plasmid"] += 1
                continue
            if _prophage(product, gene):
                stage["drop_prophage"] += 1
                continue

            # 7. terminal stop tabulated separately, never in an aa family
            stop_counts[codons[-1]] += 1
            stage["kept"] += 1
            genes.append(
                {
                    "locus_tag": f.qualifiers.get("locus_tag", [""])[0],
                    "gene": gene,
                    "product": product,
                    "contig": rec.id,
                    "codons": clean,
                    "n_codons_before_trim": len(sense),
                    "gc_percent": round(100 * (nt.count("G") + nt.count("C")) / len(nt), 3),
                }
            )
    return genes, stage, stop_counts


# --- codon statistics -----------------------------------------------------
def families(table_id: int) -> dict[str, list[str]]:
    """aa -> synonymous codons, stops excluded, sorted for determinism."""
    table = CodonTable.unambiguous_dna_by_id[table_id]
    fam: dict[str, list[str]] = {}
    for codon, aa in sorted(table.forward_table.items()):
        fam.setdefault(aa, []).append(codon)
    return fam


def codon_counts(genes: list[dict]) -> Counter:
    counts: Counter = Counter()
    for g in genes:
        counts.update(g["codons"])
    return counts


def cai_weights(counts: Counter, table_id: int) -> dict[str, float]:
    """w = codon count / count of the best synonym. Pseudocount 0.5 on zeros.

    Met, Trp and stops carry no synonymous choice and are excluded from CAI.
    """
    weights: dict[str, float] = {}
    for aa, codons in families(table_id).items():
        if aa in CAI_EXCLUDED_AA or len(codons) == 1:
            continue
        adj = {c: (counts[c] if counts[c] else PSEUDOCOUNT) for c in codons}
        best = max(adj.values())
        for c, n in adj.items():
            weights[c] = n / best
    return weights


def relative_frequencies(counts: Counter, table_id: int) -> dict[str, float]:
    """Within-family relative synonymous frequency; the harmonisation target."""
    freq: dict[str, float] = {}
    for codons in families(table_id).values():
        total = sum(counts[c] for c in codons)
        for c in codons:
            freq[c] = (counts[c] / total) if total else 0.0
    return freq


def cai(codons, weights: dict[str, float]) -> float:
    """Geometric mean of weights over scorable codons (Met/Trp/stops excluded)."""
    import math

    logs = [math.log(weights[c]) for c in codons if c in weights and weights[c] > 0]
    if not logs:
        return 0.0
    return math.exp(sum(logs) / len(logs))


def gc3(genes: list[dict]) -> float:
    """GC at third codon positions over the hygiene-passing set (seqinr::GC3 convention)."""
    gc = n = 0
    for g in genes:
        for c in g["codons"]:
            n += 1
            if c[2] in "GC":
                gc += 1
    return gc / n


def percentile(values: list[float], p: float) -> float:
    """Linear-interpolated percentile; p in 0..100."""
    s = sorted(values)
    if not s:
        raise ValueError("empty")
    k = (len(s) - 1) * p / 100
    lo, hi = int(k), min(int(k) + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


# --- reference set --------------------------------------------------------
def reference_set(genes: list[dict]) -> list[int]:
    """Indices of tuf/tsf/fusA + ribosomal proteins. Duplicates already gone."""
    idx = []
    for i, g in enumerate(genes):
        if g["gene"].lower() in REF_GENE_NAMES:
            idx.append(i)
            continue
        if _is_ribosomal(g["product"]):
            idx.append(i)
    if len(idx) < REF_SET_MIN:
        raise ValueError(
            f"reference set has {len(idx)} genes (<{REF_SET_MIN}); "
            "annotation is too poor for product-name matching, switch to HMM search"
        )
    return idx


def _is_ribosomal(product: str) -> bool:
    return bool(RIBO_PRODUCT.search(product)) and not any(
        x in product.lower() for x in RIBO_PRODUCT_EXCLUDE
    )


def ribosomal_only(genes: list[dict]) -> list[int]:
    return [i for i, g in enumerate(genes) if _is_ribosomal(g["product"])]


# --- anti-Shine-Dalgarno --------------------------------------------------
def anti_sd(records: list) -> dict:
    """Helix 45 CCUCC core from the organism's own 16S. Tail is uncertain by design."""
    seqs = []
    for rec in records:
        for f in rec.features:
            if f.type == "rRNA" and "16S" in f.qualifiers.get("product", [""])[0]:
                seqs.append(str(f.extract(rec.seq)).upper())
    if not seqs:
        raise ValueError("no annotated 16S rRNA in this assembly")

    tails = []
    for s in seqs:
        pos = s.rfind("CCTCC")  # 3' half only: helix 45 sits near the 3' end
        if pos < len(s) * 0.7:
            continue
        tails.append(s[pos + 5 :].replace("T", "U"))
    if not tails:
        raise ValueError("no CCUCC core found in the 3' third of any 16S copy")
    majority, votes = Counter(tails).most_common(1)[0]
    return {
        "core": "CCUCC",
        "downstream_tail": majority,
        "copies_16S": len(seqs),
        "copies_agreeing": votes,
        "note": "tail is uncertain; never substitute the E. coli tail",
    }
