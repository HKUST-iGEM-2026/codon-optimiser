"""S index checked against the reference implementation it is copied from.

Counts are the four two-fold families of sscu's own example set
(inst/sequences/L_kunkeei_highly.ffn, 40 highly expressed genes). With
gc3 = 0.308 the package vignette publishes S = 1.54.
"""

from __future__ import annotations

from collections import Counter

from pipeline.gate import s_index

SSCU_L_KUNKEEI = Counter(
    {"TAC": 138, "TAT": 29, "AAC": 235, "AAT": 66, "ATC": 204, "ATT": 203, "TTC": 118, "TTT": 58}
)
SSCU_GC3 = 0.308
SSCU_PUBLISHED_S = 1.54


def test_matches_sscu_published_value():
    assert round(s_index(SSCU_L_KUNKEEI, SSCU_GC3)["value"], 2) == SSCU_PUBLISHED_S


def test_no_bias_gives_zero():
    """C-ending frequency equal to the mutational expectation means S = 0."""
    balanced = Counter({c: 40 for c in SSCU_L_KUNKEEI})
    assert abs(s_index(balanced, 0.5)["value"]) < 1e-9


def test_sign_follows_the_direction_of_bias():
    t_biased = Counter({"TAC": 10, "TAT": 90, "AAC": 10, "AAT": 90,
                        "ATC": 10, "ATT": 90, "TTC": 10, "TTT": 90})
    assert s_index(t_biased, 0.5)["value"] < 0
