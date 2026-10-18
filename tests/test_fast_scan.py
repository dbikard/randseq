"""The compiled scan (randseq.fast) must agree with the numpy path, exactly.

`randseq.fast` re-implements the window encoding, the rejection of non-ACGT windows, the strand
fold and the once-per-sequence rule as one numba kernel. It is optional, and the numpy path is
the readable statement of the algorithm, so the kernel is held to it on every shipped pattern.
"""
import os

import numpy as np
import pandas as pd
import pytest

import randseq.core as core
from randseq import fast
from randseq.coding import encode_library
from randseq.example_data import get_example_data_dir
from randseq.utils import calculate_log2fc, get_lib_seq_context, get_patterns, plasmid_background

pytestmark = pytest.mark.skipif(not fast.HAVE_NUMBA, reason="numba not installed")

LEFT, RIGHT = "GTCTAGGGCGGCGGTAAAAC", "ACTAGAGCACCAGAAGTCTA"


@pytest.fixture(scope="module")
def panel_library():
    d = get_example_data_dir()
    c = pd.read_csv(os.path.join(d, "counts_panel.csv.gz"), index_col=0)
    l2 = calculate_log2fc(c[["Control_R1", "JJ1886_R1"]], reference_column="Control_R1")["JJ1886_R1"]
    fc = l2.to_numpy(dtype=float)
    return (encode_library(get_lib_seq_context(l2.index, LEFT, RIGHT)), fc,
            plasmid_background(l2.index, fc, -1))


@pytest.mark.parametrize("pattern", get_patterns(), ids=str)
def test_totals_match_numpy_on_every_shipped_pattern(panel_library, pattern):
    enc, fc, bg = panel_library
    si, mc = core.motif_occurrences(enc, pattern)
    ref = core._motif_totals(si, mc, fc, bg, -1, 4 ** (pattern[0] + pattern[2]))
    for parallel in (True, False):
        got = fast.pattern_totals(enc, pattern, fc, bg, -1, parallel=parallel)
        assert (got[0] == ref[0]).all(), "support"
        assert (got[1] == ref[1]).all(), "depleted"
        assert np.allclose(got[2], ref[2], rtol=0, atol=1e-9), "log2FC sum"
        # plasmid_background returns multiples of 2**-20, so its sums are exact in any order
        assert (got[3] == ref[3]).all(), "background sum"
        seen = ref[0] > 0
        assert (got[4][seen] == ref[4][seen]).all(), "first sequence"


def test_get_pattern_scores_is_the_same_table_either_way(panel_library, monkeypatch):
    enc, fc, bg = panel_library
    for pattern in [(6, 0, 0), (4, 4, 3), (4, 8, 4)]:
        compiled = core.get_pattern_scores(enc, fc, pattern, background=bg)
        monkeypatch.setattr(fast, "HAVE_NUMBA", False)
        plain = core.get_pattern_scores(enc, fc, pattern, background=bg)
        monkeypatch.setattr(fast, "HAVE_NUMBA", True)
        pd.testing.assert_frame_equal(compiled, plain, check_exact=False, rtol=0, atol=1e-12)


def test_a_non_acgt_base_in_the_spacer_rejects_the_window():
    """The panel library has no N, so this rule needs its own case: an N anywhere in the window,
    the undefined spacer included, rejects it -- in both paths."""
    seqs = ["CACANNNGTAC", "CACATTTGTAC", "CACATNTGTAC"]   # (3,4,4) site; spacer clean only in #2
    enc, fc, bg = encode_library(seqs), np.array([-3.0, -3.0, -3.0]), np.full(3, 0.5)
    si, mc = core.motif_occurrences(enc, (3, 4, 4))
    ref = core._motif_totals(si, mc, fc, bg, -1, 4 ** 7)
    got = fast.pattern_totals(enc, (3, 4, 4), fc, bg, -1)
    assert (got[0] == ref[0]).all()
    assert int(got[0].sum()) == int(ref[0].sum()) == 2      # one window, both strands, seq #2 only
