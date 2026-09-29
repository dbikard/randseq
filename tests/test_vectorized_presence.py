"""The vectorized motif-presence path must agree with the regex path exactly.

Presence used to be a Python regex evaluated per sequence per motif. It is now a lookup into
the window codes `_pattern_window_codes` already computes, done once per pattern. That is only
a safe substitution if the two agree on every sequence, so this checks element-for-element on
real data rather than comparing summary statistics.
"""
import os

import numpy as np
import pandas as pd
import pytest

from randseq.core import (calculate_log2fc, decode_motif_codes, encode_library, encode_motif,
                          get_lib_seq_context, get_pattern_scores, motif_presence_matrix)
from randseq.example_data import get_example_data_dir
from randseq.utils import _get_filter_for_motif

LEFT, RIGHT = "GTCTAGGGCGGCGGTAAAAC", "ACTAGAGCACCAGAAGTCTA"
PATTERNS = [(6, 0, 0), (3, 4, 4), (4, 5, 4)]


def _setup(sample):
    d = get_example_data_dir()
    df = pd.read_csv(os.path.join(d, "counts_panel.csv.gz")).set_index("seq")
    l2 = calculate_log2fc(df[["Control_R1", sample]], reference_column="Control_R1")[sample]
    enc = encode_library(get_lib_seq_context(l2.index.tolist(), LEFT, RIGHT))
    return l2, enc


@pytest.mark.parametrize("sample", ["JJ1886_R1", "B156_R1"])
def test_matrix_matches_regex_exactly(sample):
    l2, enc = _setup(sample)
    motifs, pats = [], []
    for pat in PATTERNS:
        sc = get_pattern_scores(enc, l2.tolist(), pat, log2FC_thr=-1)
        top = sc.nlargest(8, "num_sequences")
        motifs += top["motif"].tolist()
        pats += [pat] * len(top)

    matrix, handled = motif_presence_matrix(enc, motifs, pats)
    assert handled.all(), "every ACGT/N motif should take the fast path"
    for j, m in enumerate(motifs):
        expected = _get_filter_for_motif(l2, m, LEFT, RIGHT).to_numpy()
        assert np.array_equal(matrix[:, j], expected), (
            f"{sample} / {m}: vectorized presence disagrees with the regex on "
            f"{int((matrix[:, j] != expected).sum())} sequences"
        )


def test_encode_motif_round_trips():
    for pattern in [(2, 3, 2), (3, 0, 0)]:
        n = 4 ** (pattern[0] + pattern[2])
        for code in range(n):
            motif = decode_motif_codes([code], pattern)[0]
            assert encode_motif(motif, pattern) == code, (pattern, code, motif)


def test_encode_motif_declines_what_it_cannot_represent():
    """Degenerate motifs have no single code, so they must fall back rather than be encoded
    wrongly. Returning a wrong code here would silently mis-assign sequences."""
    assert encode_motif("CACNNNNGTAY", (3, 4, 4)) is None      # IUPAC ambiguity
    assert encode_motif("GGTCTC", (3, 4, 4)) is None           # wrong length
    assert encode_motif("CACGNNNGTAC", (3, 4, 4)) is None      # defined base in the spacer
    assert encode_motif("CACNNNNGTAC", (3, 4, 4)) is not None  # the real thing still works


def test_degenerate_motifs_fall_back_and_still_work():
    l2, enc = _setup("B156_R1")
    motifs = ["CACNNNNGTAY", "GGTCTC"]
    pats = [(3, 4, 4), (6, 0, 0)]
    matrix, handled = motif_presence_matrix(enc, motifs, pats)
    assert not handled[0] and handled[1], "the IUPAC motif must be left to the caller"
    assert not matrix[:, 0].any(), "an unhandled motif must be left empty, not guessed"
    expected = _get_filter_for_motif(l2, "GGTCTC", LEFT, RIGHT).to_numpy()
    assert np.array_equal(matrix[:, 1], expected)
