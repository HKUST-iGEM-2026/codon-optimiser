"""Same input + same seed -> byte-identical output. Ten runs, both modes."""

from __future__ import annotations

import json

import pytest

from app import engine as e

PROTEIN = "MSKIYAIAKELNVSKQTVSRWLNDESLSKEQLNLLAKELNVSPSELLGEEIKNNKVAVIYPNSNNPFWVEMKRGFEDAFKEND"
CDS = "".join(e.baseline(PROTEIN + "*", e.Params()))
RUNS = 10


@pytest.mark.parametrize("mode", e.MODES)
@pytest.mark.parametrize("text", [PROTEIN, CDS], ids=["protein", "cds"])
def test_ten_runs_are_byte_identical(text, mode):
    payloads = {
        json.dumps(e.optimise(text, e.Params(mode=mode)), sort_keys=True) for _ in range(RUNS)
    }
    assert len(payloads) == 1


def test_a_different_seed_is_allowed_to_differ_but_stays_stable():
    a = [json.dumps(e.optimise(PROTEIN, e.Params(mode="harmonise", seed=7)), sort_keys=True) for _ in range(3)]
    b = json.dumps(e.optimise(PROTEIN, e.Params(mode="harmonise", seed=8)), sort_keys=True)
    assert len(set(a)) == 1
    assert a[0] != b, "the seed is not reaching the harmonisation draw"


def test_manifest_carries_no_timestamp():
    """A timestamp would silently break byte-identical determinism."""
    manifest = json.dumps(e.manifest(e.Params()))
    assert "time" not in manifest.lower() and "date" not in manifest.lower()
