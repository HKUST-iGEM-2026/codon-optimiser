"""Engine contracts: what the API promises about every returned sequence."""

from __future__ import annotations

import pytest

from app import engine as e
from app import liabilities as lb
from app import organism as org

PROTEIN = (
    "MSKIYAIAKELNVSKQTVSRWLNDESLSKEQLNLLAKELNVSPSELLGEEIKNNKVAVIYPNSNNPFWVEMKRGFEDAFKEND"
    "YNVILSNVQNDINKQVSLINSMISQGVDAIIISPNEAKELVPVVEEAMNKGIKVITLDSDAPKSKRLAYIGVDNYEMGKLAAEE"
)


def test_output_translates_to_input_protein_stop_included():
    result = e.optimise(PROTEIN)
    assert org.translate(result["sequence"]) == PROTEIN + "*"
    assert result["sequence"][-3:] in org.STOP_CODONS
    assert "*" not in org.translate(result["sequence"])[:-1]


def test_output_violates_no_hard_rule():
    result = e.optimise(PROTEIN)
    seq = result["sequence"]
    assert lb.hard_violations(seq, e.GC_WINDOW_BAND) == []
    assert lb.gc_global(seq, e.GC_GLOBAL_BAND) == []


def test_avoided_motif_is_absent_and_named_in_the_manifest():
    motif = "GAATTC"  # EcoRI
    result = e.optimise(PROTEIN, e.Params(avoid=(motif,)))
    assert motif not in result["sequence"]
    assert motif not in lb.revcomp(result["sequence"])
    assert result["manifest"]["parameters"]["avoid"] == [motif]


def test_cds_input_keeps_locked_span_and_native_start():
    native = "".join(e.baseline("MKAVILT" + PROTEIN[1:60] + "*", e.Params()))
    cds = "GTG" + native[3:]  # alternative initiator
    result = e.optimise(cds, e.Params(locked=((0, 5),)))
    assert result["sequence"][:15] == cds[:15]
    assert result["input_kind"] == "cds"


def test_every_changed_codon_names_a_rule():
    result = e.optimise(PROTEIN, e.Params(avoid=("GGATCC",)))
    changed = [c for c in result["codons"] if c["status"] != "unchanged"]
    assert changed, "nothing changed: the fixture no longer exercises the explanation"
    for codon in changed:
        assert codon["rule"] is not None
        assert codon["rule"]["rule"]
        assert codon["status"] in ("changed-by-constraint", "changed-by-objective")
    for codon in result["codons"]:
        assert codon["relative_adaptiveness"] is not None


def test_infeasible_returns_the_conflicting_pair_not_a_sequence():
    """Glu has exactly two codons; forbidding both cannot be worked around."""
    with pytest.raises((e.Infeasible, e.SolverGaveUp)) as caught:
        e.optimise("MEEEEEEEEEK", e.Params(avoid=("GAA", "GAG")))
    conflict = caught.value.conflict
    assert len(conflict) >= 2
    assert all("start" in c and "end" in c for c in conflict)
    assert any(c["constraint"] == "EnforceTranslation" for c in conflict)


def test_input_cap_is_enforced():
    with pytest.raises(e.InvalidInput):
        e.optimise("M" * (e.INPUT_CAP_NT // 3))


def test_internal_stop_is_rejected():
    with pytest.raises(e.InvalidInput):
        e.optimise("MKV*KKA")
    with pytest.raises(e.InvalidInput):
        e.optimise("ATGAAATAAAAAGCT")


def test_protein_and_cds_inputs_are_told_apart():
    protein, codons, kind = e.parse_input("ACGT")
    assert kind == "protein" and codons is None and protein == "ACGT*"
    protein, codons, kind = e.parse_input(">x\nATGAAAGCT\n")
    assert kind == "cds" and codons == ["ATG", "AAA", "GCT", "TAA"] and protein == "MKA*"


def test_harmonise_tracks_the_native_distribution_more_closely_than_maximise():
    """The two modes must not collapse into each other."""
    maximised = e.optimise(PROTEIN, e.Params(mode="maximise"))
    harmonised = e.optimise(PROTEIN, e.Params(mode="harmonise"))

    def aaa_share(result):
        lysines = [c for c in result["codons"] if c["aa"] == "K"]
        return sum(c["chosen"] == "AAA" for c in lysines) / len(lysines)

    native = org.NATIVE_FREQ["AAA"]
    assert aaa_share(maximised) > aaa_share(harmonised)
    assert abs(aaa_share(harmonised) - native) < abs(aaa_share(maximised) - native)
