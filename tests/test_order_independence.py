"""The redundancy filter must depend on the candidate SET, not on row order.

Equivalent motifs used to be tie-broken on reset-index position, so which of a
motif/reverse-complement pair survived depended on the order patterns were scanned in. The
single-process and multiprocessing pipelines scan in different orders, so the same candidates
could be filtered two different ways.
"""
import os

import pandas as pd
import pytest

from randseq.coding import encode_library
from randseq.core import (calculate_log2fc, filter_redundant_patterns,
                          get_lib_seq_context, get_pattern_scores)
from randseq.example_data import get_example_data_dir

LEFT, RIGHT = "GTCTAGGGCGGCGGTAAAAC", "ACTAGAGCACCAGAAGTCTA"
PATTERNS = [(6, 0, 0), (3, 4, 4), (4, 5, 4), (3, 5, 4)]


def series_for(sample):
    df = pd.read_csv(os.path.join(get_example_data_dir(), "counts_panel.csv.gz"), index_col=0)
    return calculate_log2fc(df[["Control_R1", sample]], reference_column="Control_R1")[sample]


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
    ser = series_for(sample)
    kw = dict(log2fc_thr=-1, score_thr=0.5, min_support=3)
    expected = set(filter_redundant_patterns(cand.copy(), ser, LEFT, RIGHT, **kw)["motif"])
    for seed in range(8):
        shuffled = cand.sample(frac=1.0, random_state=seed).reset_index(drop=True)
        got = set(filter_redundant_patterns(shuffled, ser, LEFT, RIGHT, **kw)["motif"])
        assert got == expected, (
            f"{sample}, shuffle seed {seed}: filtering the same candidates in a different "
            f"order gave a different answer.\n"
            f"  only in unshuffled: {sorted(expected - got)}\n"
            f"  only in shuffled  : {sorted(got - expected)}"
        )
