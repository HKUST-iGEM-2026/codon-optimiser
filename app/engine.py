"""Frozen /api/v1 optimiser. Numerical behaviour must not change here.

Improvements go behind /api/v2 with a new engine module. Everything that can
alter a returned number is either a module constant or a request parameter that
ends up in the manifest.
"""

from __future__ import annotations

import dataclasses
import math
import random

import dnachisel as dc
import numpy

from . import liabilities as lb
from . import organism as org

ENGINE_VERSION = "1.0.0"
API_VERSION = "v1"
DEFAULT_SEED = 20240517  # fixed, documented, non-null
INPUT_CAP_NT = 3000
TRIM_CODONS = 12  # the 5' convention: weights were built on trimmed genes
LOCAL_DG_LEFT, LOCAL_DG_RIGHT = 18, 21  # per-codon folding window, 39nt
HAIRPIN_OBJECTIVE_STEM = 6
HAIRPIN_OBJECTIVE_WINDOW = 50
MODES = ("maximise", "harmonise")

# CLAUDE.md states both GC bands as global figures, so both are enforced globally:
# manufacturable 25-75% and native-like (5th-95th percentile of per-CDS GC).
GC_GLOBAL_BAND = (
    max(org.MANUFACTURABLE_GC_BAND[0], org.NATIVE_GC_BAND[0]),
    min(org.MANUFACTURABLE_GC_BAND[1], org.NATIVE_GC_BAND[1]),
)
# The 50nt sliding window is enforced on the ceiling only. A 25% floor inside a
# 28.5%-GC organism forbids sequence its own ribosomal genes contain: it rewrote
# a third of every native gene and drove rplB infeasible. GC-rich windows are the
# real synthesis liability, so the ceiling stays hard and the floor is a flag.
GC_WINDOW_BAND = (0.0, org.MANUFACTURABLE_GC_BAND[1])
GC_WINDOW_FLAG_FLOOR = org.MANUFACTURABLE_GC_BAND[0]
MAX_RANDOM_ITERS = 20000  # dnachisel's random search budget per local problem

AA_LETTERS = frozenset(org.FAMILIES) | {"*"}


class InvalidInput(ValueError):
    """Caller's fault: not a sequence we can optimise."""


class Infeasible(RuntimeError):
    """Hard constraints conflict. Never return a sequence in this case."""

    def __init__(self, message: str, conflict: list[dict]):
        super().__init__(message)
        self.conflict = conflict


class SolverGaveUp(RuntimeError):
    """The search ran out of attempts. Not the same claim as infeasible, and it
    must not be reported as one: the constraints may still be satisfiable."""

    def __init__(self, message: str, conflict: list[dict]):
        super().__init__(message)
        self.conflict = conflict


class InvariantViolated(RuntimeError):
    """A production assertion failed. Fail the request, never ship the sequence."""


@dataclasses.dataclass(frozen=True)
class Params:
    mode: str = org.DEFAULT_MODE
    seed: int = DEFAULT_SEED
    avoid: tuple[str, ...] = ()
    locked: tuple[tuple[int, int], ...] = ()  # codon index ranges, half-open

    def __post_init__(self):
        if self.mode not in MODES:
            raise InvalidInput(f"mode must be one of {MODES}")
        if not isinstance(self.seed, int) or isinstance(self.seed, bool):
            raise InvalidInput("seed must be an integer")
        for motif in self.avoid:
            # iupac_regex speaks plain ValueError; only InvalidInput maps to 400,
            # so an un-translated one would surface as a 500 on a caller mistake.
            try:
                lb.iupac_regex(motif)
            except ValueError as exc:
                raise InvalidInput(f"avoid: {exc}") from exc
        for a, b in self.locked:
            if not 0 <= a < b:
                raise InvalidInput(f"locked span {(a, b)} is not a codon range")

    def as_manifest(self) -> dict:
        return {
            "mode": self.mode,
            "seed": self.seed,
            "avoid": list(self.avoid),
            "locked_codon_spans": [list(s) for s in self.locked],
            "gc_global_band_percent": list(GC_GLOBAL_BAND),
            "gc_window_band_percent": list(GC_WINDOW_BAND),
            "gc_window_nt": lb.GC_WINDOW_NT,
            "gc_window_flag_floor_percent": GC_WINDOW_FLAG_FLOOR,
            "homopolymer_limits": lb.HOMOPOLYMER_LIMITS,
            "direct_repeat_nt": lb.DIRECT_REPEAT_NT,
            "tandem_unit_nt": lb.TANDEM_UNIT_NT,
            "input_cap_nt": INPUT_CAP_NT,
            "five_prime_window_nt": lb.FIVE_PRIME_WINDOW,
            "trim_codons_for_change_metric": TRIM_CODONS,
        }


# --- input ----------------------------------------------------------------
def parse_input(text: str) -> tuple[str, list[str] | None, str]:
    """-> (protein incl. trailing '*', original codons or None, input kind).

    A pure ACGTU string whose length is a multiple of three is read as a coding
    sequence; anything else as protein. The interpretation is echoed back to the
    caller so the guess is never silent.
    """
    cleaned = "".join(
        stripped.upper()
        for stripped in (line.strip() for line in text.splitlines())
        if not stripped.startswith(">")
    )
    cleaned = "".join(cleaned.split())
    if not cleaned:
        raise InvalidInput("empty sequence")

    if not set(cleaned) - set("ACGTU") and len(cleaned) >= 6 and len(cleaned) % 3 == 0:
        dna = cleaned.replace("U", "T")
        codons = [dna[i : i + 3] for i in range(0, len(dna), 3)]
        protein = "".join(org.CODON_TO_AA[c] for c in codons)
        if "*" in protein[:-1]:
            raise InvalidInput(f"coding sequence has an internal stop at codon {protein.index('*')}")
        if not protein.endswith("*"):
            protein += "*"
            codons.append(max(org.STOP_FREQ, key=lambda c: (org.STOP_FREQ[c], c)))
        return protein, codons, "cds"

    bad = set(cleaned) - AA_LETTERS
    if bad:
        raise InvalidInput(f"unsupported residues: {''.join(sorted(bad))}")
    if "*" in cleaned[:-1]:
        raise InvalidInput(f"protein has an internal stop at residue {cleaned.index('*')}")
    protein = cleaned if cleaned.endswith("*") else cleaned + "*"
    return protein, None, "protein"


# --- codon choice ---------------------------------------------------------
def best_codon(aa: str) -> str:
    """Most adaptive synonym; lexicographic tie-break keeps it deterministic."""
    if aa == "*":
        return max(org.STOP_FREQ, key=lambda c: (org.STOP_FREQ[c], c))
    codons = org.FAMILIES[aa]
    return max(codons, key=lambda c: (org.CAI_WEIGHTS.get(c, 1.0), c))


def _harmonised_codons(protein: str, seed: int) -> list[str]:
    """Native codon distribution, reproduced by largest-remainder allocation.

    Deterministic: the per-amino-acid quota is arithmetic, and only the order in
    which the quota is spread over that amino acid's positions is seeded.
    """
    rng = random.Random(seed)
    positions: dict[str, list[int]] = {}
    for i, aa in enumerate(protein):
        positions.setdefault(aa, []).append(i)
    out: list[str] = [""] * len(protein)
    for aa, idx in positions.items():
        codons = [best_codon(aa)] if aa == "*" or len(org.FAMILIES.get(aa, [])) == 1 else org.FAMILIES[aa]
        if len(codons) == 1:
            for i in idx:
                out[i] = codons[0]
            continue
        freq = org.NATIVE_FREQ
        exact = [(c, len(idx) * freq[c]) for c in codons]
        quota = {c: int(v) for c, v in exact}
        short = len(idx) - sum(quota.values())
        for c, _ in sorted(exact, key=lambda cv: (-(cv[1] - int(cv[1])), cv[0]))[:short]:
            quota[c] += 1
        pool = [c for c in codons for _ in range(quota[c])]
        rng.shuffle(pool)
        for i, c in zip(idx, pool):
            out[i] = c
    return out


def baseline(protein: str, params: Params) -> list[str]:
    if params.mode == "maximise":
        return [best_codon(aa) for aa in protein]
    return _harmonised_codons(protein, params.seed)


# --- solving --------------------------------------------------------------
def _constraints(protein: str, locked_nt: list[tuple[int, int]], params: Params) -> list:
    length = len(protein) * 3
    cst = [
        dc.EnforceTranslation(genetic_table="Bacterial", translation=protein, location=(0, length)),
        dc.EnforceGCContent(mini=GC_GLOBAL_BAND[0] / 100, maxi=GC_GLOBAL_BAND[1] / 100),
        dc.EnforceGCContent(
            mini=GC_WINDOW_BAND[0] / 100, maxi=GC_WINDOW_BAND[1] / 100, window=lb.GC_WINDOW_NT
        ),
        dc.UniquifyAllKmers(k=lb.DIRECT_REPEAT_NT, include_reverse_complement=False),
        dc.AvoidPattern(dc.RepeatedKmerPattern(2, lb.TANDEM_UNIT_NT)),
    ]
    cst += [
        dc.AvoidPattern(dc.HomopolymerPattern(base, n)) for base, n in lb.HOMOPOLYMER_LIMITS.items()
    ]
    cst += [dc.AvoidPattern(m) for m in params.avoid]
    cst += [dc.AvoidChanges(location=(a, b)) for a, b in locked_nt]
    return cst


def _objectives(params: Params, length: int) -> list:
    # Shallow copies: MaximizeCAI writes its cached log tables (numpy floats) back
    # into the dict it is handed, which would quietly mutate the module-level
    # reference tables that every later request reads.
    codon_objective = (
        dc.MaximizeCAI(codon_usage_table=dict(org.MAXIMISE_USAGE))
        if params.mode == "maximise"
        else dc.MatchTargetCodonUsage(codon_usage_table=dict(org.HARMONISE_USAGE))
    )
    return [
        codon_objective,
        dc.AvoidHairpins(
            stem_size=HAIRPIN_OBJECTIVE_STEM,
            hairpin_window=HAIRPIN_OBJECTIVE_WINDOW,
            location=(0, min(length, HAIRPIN_OBJECTIVE_WINDOW + 10)),
            boost=2.0,
        ),
    ]


def _conflict(seq: str, location, constraint, params: Params) -> list[dict]:
    """The conflicting pair, with coordinates: what failed, and what it fights.

    Exactly two rows. Listing every other finding that happens to overlap the window
    reads like extra causes and is not what CLAUDE.md asks for.
    """
    start = int(getattr(location, "start", 0) or 0)
    end = int(getattr(location, "end", len(seq)) or len(seq))
    return [
        {
            "constraint": str(constraint),
            "start": start,
            "end": end,
            "detail": "cannot be satisfied by any synonymous choice in this window",
        },
        {
            "constraint": "EnforceTranslation",
            "start": start - start % 3,
            "end": end + (-end % 3),
            "detail": "the amino acids in this window admit no alternative codons that clear the rule above",
        },
    ]


# --- metrics --------------------------------------------------------------
def cai(codons) -> float:
    logs = [
        math.log(org.CAI_WEIGHTS[c]) for c in codons if c in org.CAI_WEIGHTS and org.CAI_WEIGHTS[c] > 0
    ]
    return math.exp(sum(logs) / len(logs)) if logs else 0.0


def _local_dg_delta(out: str, ref: str, i: int) -> float:
    """MFE change in a 39nt window centred on codon i. Exactly 0 when the window
    is identical in both sequences, which is most of a near-native input."""
    lo = max(0, 3 * i - LOCAL_DG_LEFT)
    hi = min(len(out), 3 * i + LOCAL_DG_RIGHT)
    a, b = out[lo:hi], ref[lo:hi]
    if a == b:
        return 0.0
    return round(lb.fold_dg(a) - lb.fold_dg(b), 2)


def _explain(out_codons, ref_codons, out_seq, params: Params) -> list[dict]:
    """Per-codon record: what changed, which rule forced it, what it cost locally."""
    ref_seq = "".join(ref_codons)
    # the returned sequence violates nothing, so this is normally empty; it exists
    # so the revert test never blames a codon for a pre-existing finding.
    before = {
        (v["rule"], v["start"], v["end"])
        for v in lb.hard_violations(out_seq, GC_WINDOW_BAND, params.avoid)
    }
    records = []
    for i, (ref, chosen) in enumerate(zip(ref_codons, out_codons)):
        # Met, Trp and stops have no synonym, so their relative adaptiveness is 1
        # by definition; that is why CAI excludes them rather than scoring them.
        ref_w = org.CAI_WEIGHTS.get(ref, 1.0)
        chosen_w = org.CAI_WEIGHTS.get(chosen, 1.0)
        rule = None
        if ref == chosen:
            status = "unchanged"
        else:
            rule = _blocking_rule(out_seq, i, ref, params, before)
            status = "changed-by-constraint" if rule else "changed-by-objective"
            if rule is None:
                rule = _objective_rule(i, ref, chosen, ref_w, chosen_w, params)
        records.append(
            {
                "index": i,
                "aa": org.CODON_TO_AA[chosen],
                "original": ref,
                "chosen": chosen,
                "original_frequency": org.NATIVE_FREQ.get(ref, org.STOP_FREQ.get(ref)),
                "chosen_frequency": org.NATIVE_FREQ.get(chosen, org.STOP_FREQ.get(chosen)),
                "original_weight": ref_w,
                "chosen_weight": chosen_w,
                "relative_adaptiveness": chosen_w,
                "status": status,
                "rule": rule,
                "dg_delta": _local_dg_delta(out_seq, ref_seq, i),
            }
        )
    return records


def _objective_rule(i, ref, chosen, ref_w, chosen_w, params: Params) -> dict:
    """No hard rule forbids the original codon here, so an objective moved it.
    Name which one instead of shrugging: only two are ever active."""
    if params.mode == "maximise" and chosen_w > ref_w:
        name, detail = "codon_usage_objective", f"{chosen} is more adaptive than {ref} ({chosen_w} vs {ref_w})"
    elif params.mode == "harmonise":
        name = "codon_usage_objective"
        detail = f"{chosen} keeps the native synonymous distribution closer than {ref}"
    elif 3 * i < HAIRPIN_OBJECTIVE_WINDOW + 10:
        name, detail = "five_prime_structure_objective", f"{ref} extended a stem in the 5' window"
    else:
        name = "codon_usage_objective"
        detail = f"solver traded {ref} for {chosen} at equal or lower adaptiveness ({chosen_w} vs {ref_w})"
    return {"rule": name, "start": 3 * i, "end": 3 * i + 3, "detail": detail}


def _blocking_rule(out_seq: str, i: int, rejected: str, params: Params, before: set) -> dict | None:
    """Put the rejected codon back and see which hard rule bites.

    Scanned over the whole sequence, not a local window: a 20bp direct repeat's
    partner copy can sit anywhere, and a local scan would silently mislabel those
    constraint-forced codons as objective-driven.
    """
    variant = out_seq[: 3 * i] + rejected + out_seq[3 * i + 3 :]
    for v in lb.hard_violations(variant, GC_WINDOW_BAND, params.avoid):
        if (v["rule"], v["start"], v["end"]) in before:
            continue
        if v["start"] <= 3 * i + 2 and v["end"] > 3 * i:  # half-open: overlaps this codon
            return {
                "rule": v["rule"],
                "start": v["start"],
                "end": v["end"],
                "detail": f"{rejected} here would create {v['detail']}",
            }
    if lb.gc_global(variant, GC_GLOBAL_BAND):
        gc = lb.gc_percent(variant)
        return {
            "rule": "gc_global",
            "start": 0,
            "end": len(variant),
            "detail": f"{rejected} here would put global GC at {gc:.2f}%, outside "
                      f"{GC_GLOBAL_BAND[0]}-{GC_GLOBAL_BAND[1]}%",
        }
    return None


# --- entry point ----------------------------------------------------------
def optimise(text: str, params: Params | None = None) -> dict:
    params = params or Params()
    protein, input_codons, kind = parse_input(text)
    if len(protein) * 3 > INPUT_CAP_NT:
        raise InvalidInput(f"{len(protein) * 3}nt output exceeds the {INPUT_CAP_NT}nt cap")

    warnings: list[str] = []
    start = baseline(protein, params)
    locked = [(a, min(b, len(protein))) for a, b in params.locked]
    if any(a >= len(protein) for a, _ in locked):
        raise InvalidInput("locked span starts past the end of the protein")
    if input_codons is not None:
        for a, b in locked:
            start[a:b] = input_codons[a:b]

    # The start codon is never left to the solver: keep the native initiator for a
    # pasted CDS, otherwise ATG. GTG/TTG initiators translate as V/L here, which is
    # the honest round-trip convention (NCBI's initiator-Met rule is annotation).
    if input_codons is not None and input_codons[0] in org.ORGANISM["start_codons"]:
        start[0] = input_codons[0]
        locked.append((0, 1))
    elif protein[0] == "M":
        start[0] = "ATG"
        locked.append((0, 1))
    else:
        warnings.append(
            f"protein starts with {protein[0]}, not M: the construct has no canonical start codon"
        )

    locked_nt = sorted({(a * 3, b * 3) for a, b in locked})
    random.seed(params.seed)
    numpy.random.seed(params.seed)
    problem = dc.DnaOptimizationProblem(
        sequence="".join(start),
        constraints=_constraints(protein, locked_nt, params),
        objectives=_objectives(params, len(protein) * 3),
        logger=None,
    )
    problem.max_random_iters = MAX_RANDOM_ITERS
    try:
        problem.resolve_constraints()
        problem.optimize()
    except dc.NoSolutionError as exc:
        conflict = _conflict(problem.sequence, exc.location, exc.constraint, params)
        if "Random search did not find" in str(exc):
            raise SolverGaveUp(str(exc), conflict) from exc
        raise Infeasible(str(exc), conflict) from exc

    seq = problem.sequence
    out_codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]

    # --- production invariants; never ship a sequence that fails one ---
    if org.translate(seq) != protein:
        raise InvariantViolated("output does not translate to the input protein")
    if "*" in org.translate(seq)[:-1]:
        raise InvariantViolated("output contains an internal stop codon")
    for a, b in locked_nt:
        if seq[a:b] != "".join(start)[a:b]:
            raise InvariantViolated(f"locked span {a}-{b} was modified")
    violations = lb.hard_violations(seq, GC_WINDOW_BAND, params.avoid)
    violations += lb.gc_global(seq, GC_GLOBAL_BAND)
    if violations:
        raise Infeasible(
            "solver returned a sequence that violates a hard constraint",
            [{"constraint": v["rule"], "start": v["start"], "end": v["end"], "detail": v["detail"]} for v in violations],
        )

    reference = input_codons if input_codons is not None else start
    codon_records = _explain(out_codons, reference, seq, params)
    changed = [r for r in codon_records if r["status"] != "unchanged"]
    windows = [
        lb.gc_percent(seq[i : i + lb.GC_WINDOW_NT]) for i in range(0, max(1, len(seq) - lb.GC_WINDOW_NT + 1))
    ]

    return {
        "sequence": seq,
        "protein": protein,
        "input_kind": kind,
        "codons": codon_records,
        "metrics": {
            "length_nt": len(seq),
            "cai": round(cai(out_codons), 4),
            "cai_of_input": round(cai(input_codons), 4) if input_codons else None,
            "gc_percent": round(lb.gc_percent(seq), 2),
            "gc_window_min": round(min(windows), 2),
            "gc_window_max": round(max(windows), 2),
            "changed_codons": len(changed),
            "changed_codons_fraction": round(len(changed) / len(out_codons), 4),
            "changed_codons_after_trim": sum(1 for r in changed if r["index"] >= TRIM_CODONS),
            "changed_by_constraint": sum(1 for r in changed if r["status"] == "changed-by-constraint"),
            "five_prime": lb.five_prime_dg(seq),
        },
        "flags": (
            lb.hairpin_flags(seq)
            + lb.terminators(seq)
            + lb.at_rich_windows(seq, GC_WINDOW_FLAG_FLOOR)
        ),
        "warnings": warnings,
        "manifest": manifest(params),
    }


def manifest(params: Params) -> dict:
    """Carried by every result and every export. No timestamps: they would break
    byte-identical determinism."""
    return {
        "engine_version": ENGINE_VERSION,
        "api_version": API_VERSION,
        "table_version": org.TABLE_VERSION,
        "assembly_accession": org.ACCESSION,
        "organism": org.ORGANISM["organism"]["organism"],
        "translation_table": org.TABLE_ID,
        "dnachisel_version": dc.__version__,
        "parameters": params.as_manifest(),
    }
