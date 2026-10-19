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

from randseq.coding import (_decode_motif_codes, encode_library, _encode_motif,
                            _motif_presence_matrix)
from randseq.fast import _get_pattern_scores
from randseq.utils import calculate_log2fc, get_lib_seq_context
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


def _top_motifs(l2, enc):
    """The eight most carried motifs of each pattern."""
    motifs, pats = [], []
    for pat in PATTERNS:
        sc = _get_pattern_scores(enc, l2.tolist(), pat, log2FC_thr=-1)
        top = sc.nlargest(8, "num_sequences")
        motifs += top["motif"].tolist()
        pats += [pat] * len(top)
    return motifs, pats


@pytest.mark.parametrize("sample", ["JJ1886_R1", "B156_R1"])
def test_matrix_matches_regex_exactly(sample):
    l2, enc = _setup(sample)
    motifs, pats = _top_motifs(l2, enc)

    matrix, handled = _motif_presence_matrix(enc, motifs, pats)
    assert handled.all(), "every ACGT/N motif should take the fast path"
    for j, m in enumerate(motifs):
        expected = _get_filter_for_motif(l2, m, LEFT, RIGHT).to_numpy()
        assert np.array_equal(matrix[:, j], expected), (
            f"{sample} / {m}: vectorized presence disagrees with the regex on "
            f"{int((matrix[:, j] != expected).sum())} sequences"
        )


def test_inside_matches_regex_on_inserts():
    """With `insert`, the second matrix is what the regex finds in the insert alone, read off the
    flanked windows. The panel's inserts run from 8 to 151 nt, so where the insert ends differs from
    row to row; the first matrix is unchanged."""
    l2, enc = _setup("JJ1886_R1")
    motifs, pats = _top_motifs(l2, enc)
    lengths = np.array([len(s) for s in l2.index])
    matrix, inside, handled = _motif_presence_matrix(enc, motifs, pats,
                                                     (len(LEFT), len(RIGHT), lengths))
    assert handled.all() and np.array_equal(matrix, _motif_presence_matrix(enc, motifs, pats)[0])
    assert (matrix & ~inside).any(), "some carriers should need the flanks, or the test is empty"
    for j, m in enumerate(motifs):
        expected = _get_filter_for_motif(l2, m, "", "").to_numpy()
        assert np.array_equal(inside[:, j], expected), (
            f"{m}: internal presence disagrees with the regex on the inserts on "
            f"{int((inside[:, j] != expected).sum())} sequences")


def test_encode_motif_round_trips():
    for pattern in [(2, 3, 2), (3, 0, 0)]:
        n = 4 ** (pattern[0] + pattern[2])
        for code in range(n):
            motif = _decode_motif_codes([code], pattern)[0]
            assert _encode_motif(motif, pattern) == code, (pattern, code, motif)


def test_encode_motif_declines_what_it_cannot_represent():
    """Degenerate motifs have no single code, so they must fall back rather than be encoded
    wrongly. Returning a wrong code here would silently mis-assign sequences."""
    assert _encode_motif("CACNNNNGTAY", (3, 4, 4)) is None      # IUPAC ambiguity
    assert _encode_motif("GGTCTC", (3, 4, 4)) is None           # wrong length
    assert _encode_motif("CACGNNNGTAC", (3, 4, 4)) is None      # defined base in the spacer
    assert _encode_motif("CACNNNNGTAC", (3, 4, 4)) is not None  # the real thing still works


def test_degenerate_motifs_fall_back_and_still_work():
    l2, enc = _setup("B156_R1")
    motifs = ["CACNNNNGTAY", "GGTCTC"]
    pats = [(3, 4, 4), (6, 0, 0)]
    matrix, handled = _motif_presence_matrix(enc, motifs, pats)
    assert not handled[0] and handled[1], "the IUPAC motif must be left to the caller"
    assert not matrix[:, 0].any(), "an unhandled motif must be left empty, not guessed"
    expected = _get_filter_for_motif(l2, "GGTCTC", LEFT, RIGHT).to_numpy()
    assert np.array_equal(matrix[:, 1], expected)
