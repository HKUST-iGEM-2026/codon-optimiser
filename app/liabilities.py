"""Sequence liability scans. One definition of every rule, used three ways:

1. to configure the DnaChisel solver,
2. to explain which rule forced a codon (by testing the rejected codon),
3. to assert in production that the returned sequence violates nothing.

Hard rules return violations; soft rules return flags. Coordinates are 0-based
half-open on the sequence passed in.
"""

from __future__ import annotations

import re

import RNA

# --- fixed conventions (CLAUDE.md) -----------------------------------------
HOMOPOLYMER_LIMITS = {"A": 10, "T": 10, "G": 6, "C": 6}  # runs of ONE base
DIRECT_REPEAT_NT = 20
TANDEM_UNIT_NT = 8
GC_WINDOW_NT = 50
HAIRPIN_DG = -15.0  # kcal/mol
HAIRPIN_WINDOW_NT = 60
HAIRPIN_WINDOW_STEP = 20
FIVE_PRIME_WINDOW = (-30, 40)  # from the A of the start codon, assembled construct
TERMINATOR_STEM = (5, 8)
TERMINATOR_LOOP = (3, 8)
TERMINATOR_T_RUN = 4
TERMINATOR_T_GAP = 5

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[GC]",
    "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]",
    "V": "[ACG]", "N": "[ACGT]",
}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def iupac_regex(motif: str) -> re.Pattern:
    if not motif:
        # An empty pattern matches at every position, so it is a constraint no
        # sequence can satisfy: the solver grinds through its entire random-search
        # budget (minutes of CPU) before admitting defeat. Reject it at the door.
        raise ValueError("motif must not be empty")
    try:
        return re.compile("".join(IUPAC[b] for b in motif.upper()))
    except KeyError as exc:
        raise ValueError(f"{motif!r} is not a DNA/IUPAC motif") from exc


# --- hard rules ------------------------------------------------------------
def homopolymers(seq: str) -> list[dict]:
    out = []
    for base, limit in HOMOPOLYMER_LIMITS.items():
        for m in re.finditer(f"{base}{{{limit},}}", seq):
            out.append(
                {
                    "rule": "homopolymer",
                    "start": m.start(),
                    "end": m.end(),
                    "detail": f"{len(m.group())}x{base} (limit {limit})",
                }
            )
    return sorted(out, key=lambda v: v["start"])


def direct_repeats(seq: str, k: int = DIRECT_REPEAT_NT) -> list[dict]:
    """Any k-mer occurring more than once. Direct only: reverse complements are
    a different liability and are not what CLAUDE.md's repeat limit names."""
    first: dict[str, int] = {}
    out = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if kmer in first:
            out.append(
                {
                    "rule": "direct_repeat",
                    "start": i,
                    "end": i + k,
                    "detail": f"{k}bp repeat of position {first[kmer]}",
                }
            )
        else:
            first[kmer] = i
    return out


def tandem_repeats(seq: str, unit: int = TANDEM_UNIT_NT) -> list[dict]:
    out = []
    for i in range(len(seq) - 2 * unit + 1):
        if seq[i : i + unit] == seq[i + unit : i + 2 * unit]:
            out.append(
                {
                    "rule": "tandem_repeat",
                    "start": i,
                    "end": i + 2 * unit,
                    "detail": f"{unit}bp unit repeated twice",
                }
            )
    return out


def gc_percent(seq: str) -> float:
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def gc_global(seq: str, band: tuple[float, float]) -> list[dict]:
    gc = gc_percent(seq)
    if band[0] <= gc <= band[1]:
        return []
    return [
        {
            "rule": "gc_global",
            "start": 0,
            "end": len(seq),
            "detail": f"{gc:.1f}% outside {band[0]}-{band[1]}%",
        }
    ]


def at_rich_windows(seq: str, floor: float, window: int = GC_WINDOW_NT) -> list[dict]:
    """50nt windows below the manufacturable floor. A flag, not a constraint: this
    organism sits at 28.5% GC, so a 25% local floor would fight every native gene
    it is asked to reproduce. Synthesis vendors still care, hence the flag.

    Overlapping windows are collapsed into one span reporting the worst window in it;
    a 300nt AT-rich gene would otherwise raise 250 near-identical flags.
    """
    runs: list[dict] = []
    for i in range(max(1, len(seq) - window + 1)):
        gc = gc_percent(seq[i : i + window])
        if gc >= floor:
            continue
        end = min(len(seq), i + window)
        if runs and i <= runs[-1]["end"]:
            run = runs[-1]
            run["end"] = max(run["end"], end)
            run["worst"] = min(run["worst"], gc)
            run["windows"] += 1
        else:
            runs.append({"start": i, "end": end, "worst": gc, "windows": 1})
    return [
        {
            "flag": "at_rich_window",
            "start": run["start"],
            "end": run["end"],
            "detail": f"{run['windows']} overlapping {window}nt windows, worst {run['worst']:.1f}% GC, "
                      f"below the {floor}% synthesis floor",
        }
        for run in runs
    ]


def gc_windows(seq: str, band: tuple[float, float], window: int = GC_WINDOW_NT) -> list[dict]:
    if len(seq) <= window:
        return gc_global(seq, band)
    out = []
    for i in range(len(seq) - window + 1):
        gc = gc_percent(seq[i : i + window])
        if not band[0] <= gc <= band[1]:
            out.append(
                {
                    "rule": "gc_window",
                    "start": i,
                    "end": i + window,
                    "detail": f"{gc:.1f}% in {window}nt window, outside {band[0]}-{band[1]}%",
                }
            )
    return out


def motifs(seq: str, patterns: tuple[str, ...]) -> list[dict]:
    out = []
    for pattern in patterns:
        rx = iupac_regex(pattern)
        for m in rx.finditer(seq):
            out.append(
                {
                    "rule": "avoided_motif",
                    "start": m.start(),
                    "end": m.end(),
                    "detail": f"{pattern} present",
                }
            )
        for m in rx.finditer(revcomp(seq)):  # sites are double-stranded
            start = len(seq) - m.end()
            out.append(
                {
                    "rule": "avoided_motif",
                    "start": start,
                    "end": start + (m.end() - m.start()),
                    "detail": f"{pattern} present (reverse strand)",
                }
            )
    return sorted(out, key=lambda v: v["start"])


def hard_violations(seq: str, band: tuple[float, float], avoid: tuple[str, ...] = ()) -> list[dict]:
    """Every hard rule except translation, which the engine asserts separately."""
    return (
        homopolymers(seq)
        + direct_repeats(seq)
        + tandem_repeats(seq)
        + gc_global(seq, band)
        + gc_windows(seq, band)
        + motifs(seq, avoid)
    )


# --- soft flags ------------------------------------------------------------
def fold_dg(seq: str) -> float:
    """MFE in kcal/mol. Near zero means open/unstructured."""
    return RNA.fold(seq.replace("T", "U"))[1]


def five_prime_dg(construct: str, cds_start: int) -> dict:
    """-30..+40 from the A of the start codon, on the assembled construct.

    Near zero (weak, open) is what good expression wants. Do not invert the sign.
    """
    start = max(0, cds_start + FIVE_PRIME_WINDOW[0])
    end = min(len(construct), cds_start + FIVE_PRIME_WINDOW[1])
    window = construct[start:end]
    return {
        "window": [start - cds_start, end - cds_start],
        "upstream_nt_available": cds_start - start,
        "dg": round(fold_dg(window), 2),
        "structure": RNA.fold(window.replace("T", "U"))[0],
    }


def hairpin_flags(seq: str) -> list[dict]:
    """ΔG < -15 in any 60nt window, plus any base pair in the first 50nt."""
    flags = []
    for i in range(0, max(1, len(seq) - HAIRPIN_WINDOW_NT + 1), HAIRPIN_WINDOW_STEP):
        window = seq[i : i + HAIRPIN_WINDOW_NT]
        dg = fold_dg(window)
        if dg < HAIRPIN_DG:
            flags.append(
                {
                    "flag": "hairpin",
                    "start": i,
                    "end": i + len(window),
                    "detail": f"window MFE {dg:.1f} kcal/mol < {HAIRPIN_DG}",
                }
            )
    head = seq[:50]
    structure = RNA.fold(head.replace("T", "U"))[0]
    if "(" in structure:
        flags.append(
            {
                "flag": "five_prime_hairpin",
                "start": 0,
                "end": len(head),
                "detail": f"paired bases in the first 50nt: {structure}",
            }
        )
    return flags


def terminators(seq: str) -> list[dict]:
    """Stem 5-8bp, loop 3-8nt, then >=4 consecutive T within 0-5nt. Heuristic."""
    out = []
    lo_stem, hi_stem = TERMINATOR_STEM
    lo_loop, hi_loop = TERMINATOR_LOOP
    for i in range(len(seq)):
        for stem in range(hi_stem, lo_stem - 1, -1):
            for loop in range(lo_loop, hi_loop + 1):
                right = i + stem + loop
                if right + stem > len(seq):
                    continue
                if seq[i : i + stem] != revcomp(seq[right : right + stem]):
                    continue
                tail = seq[right + stem : right + stem + TERMINATOR_T_GAP + TERMINATOR_T_RUN]
                if re.search(f"T{{{TERMINATOR_T_RUN},}}", tail):
                    out.append(
                        {
                            "flag": "terminator_like",
                            "start": i,
                            "end": right + stem + len(tail),
                            "detail": f"stem {stem}bp, loop {loop}nt, T-run downstream",
                        }
                    )
                    break
            else:
                continue
            break
    return out
