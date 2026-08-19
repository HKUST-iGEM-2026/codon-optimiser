# CLAUDE.md

## Project

A codon optimiser web app for a single bacterium: **RefSeq assembly `GCF_000409755.1`**.

Paste a protein sequence → get a DNA coding sequence, plus a per-codon explanation of why each codon was chosen. The explanation is the product; the sequence is commodity.

## Architecture

Two separate halves. Do not blur them.

**1. Offline pipeline** (`pipeline/`) — run once, by hand, on a laptop. Outputs `data/organism.json` and stops. It is not a product surface, is not called at runtime, and must not be generalised, parameterised for other organisms, or wired into CI beyond the E. coli control test.

**2. Web service** (`app/`) — reads `data/organism.json` at process start into a module-level dict. Never reads a database for reference data.

## Stack

- Python 3.12, FastAPI, one container
- Biopython (parsing/translation only), DnaChisel (constraint solving), ViennaRNA (folding, PyPI wheel)
- Jinja templates + vanilla JS. No React, no build pipeline
- SQLite, for permalink storage only
- `data/organism.json` lives in git

## UI

**Designs come from Google Stitch via MCP.** Fetch the design and implement to it. Do not invent layouts, spacing, or component structure.

One component Stitch will not provide — build it by hand: the **codon track**. A row of codons coloured on a continuous scale by relative adaptiveness, with visually distinct marks for `changed-by-constraint` / `changed-by-objective` / `unchanged`. Hovering a codon shows original codon, chosen codon, both frequencies, the exact rule that forced the change, and the local ΔG delta.

**What Stitch actually delivered (2026-08-16).** `create_design_system` worked and is the
source of the visual system: asset **`assets/17898767070783141685`**, display name
"Instrument" — LIGHT mode, seed colour `#0F766E`, Inter headline/body, `ROUND_EIGHT`, plus a
`design.md` that fixes the palette (surface `#FBFBFA`, card `#FFFFFF`, hairline `#E4E6E3`, ink
`#141614`, muted `#6B7069`, amber `#B45309` for constraint changes, indigo `#4338CA` for
objective changes, adaptiveness ramp `#F2EFE4` → `#0F5132`), the type pairing (Inter for prose,
JetBrains Mono for every nucleotide/codon/coordinate), the 8px grid with a 1200px centred
column, radii (card 12px, control 8px, cell 3px) and the component specs. `app/static/app.css`
implements exactly those tokens; the design MD is the contract.

**`generate_screen_from_text` could not be used.** Every call — both `GEMINI_3_1_PRO` and
`GEMINI_3_FLASH`, three attempts — dies on a 30 s MCP transport timeout while the tool's own
docs say generation "can take a few minutes", and the job does not survive the abort
(`get_project` reports 0 `screenInstances` after 10 × 30 s of polling). Screen mockups are
therefore not available to implement against; the layout follows the design system's documented
layout and component rules. When Stitch screens can be fetched, re-implement to them and keep
the codon track as built.

**Codon track, as built.** `app/static/app.js` + `.codon` in `app/static/app.css`: 28×22px
monospace cells, 1px gaps, background interpolated on the adaptiveness ramp, a 2px amber bottom
bar for `changed-by-constraint`, indigo for `changed-by-objective`, no bar for `unchanged`, and
a "colour by change" alternative view. Cells are real buttons: hover *and* keyboard focus open
the tooltip, Escape closes it, and each cell carries an `aria-label`. The tooltip shows codon
index, amino acid, original codon and its native frequency, chosen codon and its frequency,
relative adaptiveness, status, the exact rule with coordinates, the rule's plain-language
detail, and the local ΔG delta.

## Step 1 — before anything else

Download the assembly with the NCBI Datasets CLI (CDS + GBFF + GFF3, RefSeq). Then read from the downloaded files and record here:

- Organism name and strain: `Clostridium butyricum DSM 10702` (strain `DSM 10702`, type strain)
- Taxid: `1316931`
- Translation table (from `/transl_table` in the GBFF): `11` (unanimous across all 4,244 CDS; lineage Bacillota/Clostridia, no non-standard-table clade → no hard warning)
- Genomic GC%: `28.4992` (4,596,811 nt over 207 WGS contigs, no plasmid records)
- CDS count: `4244` (87 flagged pseudo; 3,721 survive the hygiene filter)

Do not look these up from a secondary source. Read them from the assembly.

## Step 2 — the gate

After building the codon tables, compute and record four numbers:

| Check | Record |
|---|---|
| Rank of ribosomal proteins by CAI | `median rank 132 / 3,721 = 3.55th percentile` (best 3, worst 1,162; under genome-wide weights instead: median 2,790.5 / 74.99th — the AT-mutational background and translational selection point opposite ways here) |
| *S* index (Bioconductor `sscu`) | `1.4453` (Tyr 1.61, Phe 1.69, Asn 1.35, Ile 1.27; genomic GC3 = 0.1611) |
| eCAI (95th percentile of 1,000 composition-matched random sequences) | `0.5482`, observed reference-set CAI `0.7944` → **passes** (null mean 0.5364, seed 20240517) |
| Overlap between the ribosomal set and the self-consistent set | `6 / 46 = 13.0%` (self-consistent set of 186, converged in 12 iterations). Low because iterating from genome-wide weights climbs the AT-mutational gradient, not the translational one. E. coli control run gives 33/45 = 73%, so the procedure itself is sound. |

**This decides the tool's default mode.** *S* clearly > 0 and eCAI passes → default is maximise codon adaptation. Otherwise → default is harmonisation (match the native codon distribution). Display all four numbers in the finished UI.

Recorded default mode: **maximise**. *S* = 1.4453 is far above sscu's own "translational selection is acting" threshold of 0.3, and eCAI passes by 0.25. Threshold and formula verified against sscu's L. kunkeei example: our `s_index` returns 1.5419 where the package vignette publishes 1.54.

Table version pinned by the build: `1.0.0+11359a0d03f0`.

Do not build the optimiser before this is recorded.

## Fixed conventions

Never change these without updating this file. They are the reason results are reproducible.

**Hygiene filter, in this exact order:**

1. Drop CDS not a multiple of three
2. Drop `partial` records and pseudogenes
3. Require canonical start (ATG/GTG/TTG) and valid stop under the organism's table
4. Require **≥ 80 codons, measured BEFORE trimming**
5. **Then trim the first 12 codons**
6. Skip codons containing ambiguity codes; do not drop the gene
7. Exclude terminal stop from amino-acid families; tabulate stops separately
8. Deduplicate identical sequences
9. Exclude plasmid and prophage genes

**CAI:** exclude Met, Trp and stops. Weight = codon frequency ÷ frequency of the best synonym in its family. Pseudocount 0.5 on zero counts.

**Reference set:** `tufA`, `tsf`, `fusA` + ribosomal proteins, found by case-insensitive product-name match on `ribosomal (subunit )?protein`, excluding products mentioning methyl or hydroxy. Expect ~40. Below 20, stop and switch to HMM search. (Amended 2026-08-16: RefSeq names E. coli K-12 "30S ribosomal **subunit** protein S12", so the plain substring found 5 genes there and the required control test could not run. Gene-name match covers `tuf`/`tufA`/`tufB`/`tsf`/`fusA` because this assembly annotates EF-Tu as `tuf`. Target organism: 49 genes; E. coli control: 49.)

| Parameter | Value |
|---|---|
| 5′ structure window | −30 to +40 nt from the A of the start codon, on the **assembled construct** |
| GC sliding window | 50 nt |
| Manufacturable GC band | 25–75% **global**, enforced globally |
| Native-like GC band | 5th–95th percentile of per-CDS GC in this organism = **23.98–34.266%**, enforced globally. Effective global band is the intersection: **25.0–34.266%** |
| Homopolymer limit | run of ≥10 identical A, ≥10 identical T, ≥6 identical G, or ≥6 identical C |
| Direct repeat limit | ≥20 bp; ≥8–10 bp local/tandem |
| Hairpin flag | ΔG < −15 kcal/mol, or any hairpin in the first 50 nt |
| Terminator heuristic | stem ≥5–8 bp, loop 3–8 nt, then ≥4 consecutive T within 0–5 nt |
| Input cap | 3,000 nt |
| Default RNG seed | `20240517` |

**How the 50 nt window is enforced (decided 2026-08-16).** Both GC bands above are
stated as global figures, so both are enforced globally. On the 50 nt sliding
window only the **ceiling (75%)** is a hard constraint; a window below 25% raises
an `at_rich_window` flag instead. Reason: this organism sits at 28.5% GC, its own
ribosomal genes contain 50 nt windows down to 12% GC, and a hard 25% local floor
rewrote 28–38% of every native gene and made `rplB` unsolvable — it would have
failed the round-trip test in this file. GC-rich windows are the real synthesis
liability; AT-rich ones are native here, so they are reported, not fought.

## Hard rules

- **Hard constraints are never violated.** On infeasibility, return the conflicting constraint pair with coordinates. Never return a sequence that violates one.
- **Assert round-trip translation in production**, not only in tests: output must translate to the exact input protein under the organism's table, stop included.
- **Deterministic.** Same input + same seed → byte-identical output. Ten runs, asserted in CI.
- **Every result and every export carries a manifest**: engine version, table version, assembly accession, translation table, all parameters, seed.
- **`/api/v1` numerical behaviour is frozen.** Improvements go behind `/api/v2`.
- **Never claim an expression yield or fold-change.** Permitted claim: matches this organism's codon statistics and removes named sequence liabilities.
- **Sequences are ephemeral** — in memory for the request only. Not logged, not written to disk, not sent to third parties. No session replay. Scrub request bodies before any error reporting. (FastAPI's default validation handler echoes the offending body; it is overridden to return field names only.) The single carve-out is `POST /api/v1/optimise` with `permalink: true`: because output is deterministic, the row stores the *request* and nothing else, keyed by a hash of request + engine + table version, and the result is recomputed on read. No other route writes anything.
- Rate-limit on **CPU-seconds per IP**, not request count.

## Traps

These cause silent, plausible-looking wrong answers.

- **Translation table** comes from `/transl_table`. Never hardcode 11. Hard-warn if the lineage contains Mycoplasm / Spiroplasma / Ureaplasma / Entomoplasma / Mesoplasma / Gracilibacteria / Absconditabacteria / SR1 while the table reads 11. NCBI has 27 tables with non-contiguous IDs (1–6, 9–16, 21–33) — never loop 1–33.
- **The homopolymer rule is runs of ONE base.** "10 consecutive A/T" would mean any 10-base run drawn from {A,T}, which occurs constantly and will fight the solver.
- **The ≥80-codon filter is measured before the 12-codon trim.** Effective floor is 68 counted codons.
- **`tufA` and `tufB` are often byte-identical.** Deduplicate or EF-Tu is double-weighted.
- **The round-trip metric needs a 5′ convention.** Weights are built on trimmed genes but re-optimisation rewrites the whole CDS, including the AT-rich 5′ window. Exclude the first 12 codons from the changed-codon count, or report it twice (windowed and whole-gene). Otherwise large churn appears for correct reasons.
- **5′ ΔG should be near zero (weak/open) for good expression.** Do not invert the sign. This is the opposite of mRNA-therapeutic design.
- **16S copy number is 1–15.** A single copy is normal, not a failed run. If copies disagree in the anti-SD region, use the majority and note it.
- **Anti-SD:** anchor on helix 45, take the `CCUCC` core, treat the downstream tail as uncertain. Never use E. coli's tail.
- **`Bio.SeqUtils.CodonUsage` was removed in Biopython 1.82.** Ignore any tutorial using it.
- **Do not use the Kazusa Codon Usage Database** for anything but a cross-check. It is frozen at GenBank release 160.0 (June 2007) and its tables are often built from a handful of sequences.

## Validation

- **E. coli control test (`tests/test_control.py`) — required, keep in CI.** Run the identical table-building code on `GCF_000005845.2`, correlate CAI against PaxDb protein abundance. Must be clearly positive. This is the only test that catches a silent hygiene-filter bug.
  Recorded run (2026-08-16): **Spearman ρ = 0.588** (Pearson 0.609) over 3,560 genes joined on b-numbers against the PaxDb integrated whole-organism dataset; E. coli *S* = 1.1659; 3,870 genes kept; reference set 49. Test floor is ρ > 0.30.
- **Round-trip test:** re-optimising a native ribosomal protein from this organism must change few codons.
  Recorded run: `rpsJ` 22.3%, `rplC` 35.7%, `rplB` 35.7%, `rpsC` 38.1% of codons changed, every change attributed to the codon usage objective (none constraint-forced), output CAI 0.985–1.0 vs native 0.766–0.866. A poorly adapted native control gene (`K670_RS15945`, native CAI 0.386) moves 52.4%, so the test asserts ribosomal churn < control churn − 10 points, output CAI ≥ native CAI, and churn < 45%. Re-optimising our own output is a fixed point (0 changes).
- **Determinism test:** ten identical runs, byte-identical output.
- **Invariants asserted in production:** exact translation; no internal stops; locked spans unchanged; no hard-constraint motif present in output; windowed GC in band.

**Two failure modes, never conflated.** DnaChisel raising `NoSolutionError` after its random search exhausts its attempts is *not* proof of infeasibility: that path returns `503 solver_gave_up`, while a genuinely unsatisfiable constraint set returns `422 infeasible`. Both refuse to return a sequence.

**Per-codon attribution is a whole-sequence test.** "Which rule forced this codon" is answered by putting the rejected codon back and rescanning the *entire* sequence, because a 20 bp direct repeat's partner copy can sit anywhere. A local-window scan mislabelled constraint-forced codons as objective-driven.

## Do not build

Out of scope. Do not add these, and do not suggest them.

- Machine learning models of any kind
- Any second organism, or a host/organism picker
- A job queue, Redis, Postgres, or worker processes
- User accounts, login, or payments
- Divergence maximisation, homology minimisation, k-mer/BLAST-identity avoidance, or fragmentation into sub-threshold oligos
