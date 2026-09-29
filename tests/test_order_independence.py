"""The redundancy filter must depend on the candidate SET, not on row order.

Equivalent motifs used to be tie-broken on reset-index position, so which of a
motif/reverse-complement pair survived depended on the order patterns were scanned in. The
single-process and multiprocessing pipelines scan in different orders, so the same candidates
could be filtered two different ways.
"""
import os

import pandas as pd
import pytest

from randseq.core import (calculate_log2fc, encode_library, filter_redundant_patterns,
                          get_lib_seq_context, get_pattern_scores)
from randseq.example_data import get_example_data_dir

LEFT, RIGHT = "GTCTAGGGCGGCGGTAAAAC", "ACTAGAGCACCAGAAGTCTA"
PATTERNS = [(6, 0, 0), (3, 4, 4), (4, 5, 4), (3, 5, 4)]


def candidates_for(sample):
    d = get_example_data_dir()
    df = pd.read_csv(os.path.join(d, "counts_panel.csv.gz")).set_index("seq")
    l2 = calculate_log2fc(df[["Control_R1", sample]], reference_column="Control_R1")[sample]
    enc = encode_library(get_lib_seq_context(l2.index.tolist(), LEFT, RIGHT))
    fc = l2.tolist()
    out = []
    for pat in PATTERNS:
        sc = get_pattern_scores(enc, fc, pat, log2FC_thr=-1)
        sc["pattern"] = str(pat)
        out.append(sc[(sc["fraction_depleted"] >= 0.5) & (sc["num_sequences"] >= 3)])
    return pd.concat(out, ignore_index=True)


@pytest.mark.parametrize("sample", ["JJ1886_R1", "LMR_503_R1", "B156_R1"])
def test_filter_is_independent_of_candidate_order(sample):
    cand = candidates_for(sample)
    assert not cand.empty, f"{sample} produced no candidates to filter"
    expected = set(filter_redundant_patterns(cand.copy())["motif"])
    for seed in range(8):
        shuffled = cand.sample(frac=1.0, random_state=seed).reset_index(drop=True)
        got = set(filter_redundant_patterns(shuffled)["motif"])
        assert got == expected, (
            f"{sample}, shuffle seed {seed}: filtering the same candidates in a different "
            f"order gave a different answer.\n"
            f"  only in unshuffled: {sorted(expected - got)}\n"
            f"  only in shuffled  : {sorted(got - expected)}"
        )


def test_equivalent_motifs_tie_break_on_the_motif_not_the_row():
    """A motif and its reverse complement always score identically, so the tie-break decides
    which survives. It must decide on content, so the same pair resolves the same way whichever
    order the two rows arrive in."""
    pair = pd.DataFrame([
        dict(motif="GGTCTC", fraction_depleted=0.9, num_sequences=50, avg_log2fc=-3.0,
             pattern="(6, 0, 0)"),
        dict(motif="GAGACC", fraction_depleted=0.9, num_sequences=50, avg_log2fc=-3.0,
             pattern="(6, 0, 0)"),
    ])
    kept_a = set(filter_redundant_patterns(pair)["motif"])
    kept_b = set(filter_redundant_patterns(pair.iloc[::-1].reset_index(drop=True))["motif"])
    assert len(kept_a) == 1, f"one of the pair should survive, got {kept_a}"
    assert kept_a == kept_b, f"row order changed the survivor: {kept_a} vs {kept_b}"
