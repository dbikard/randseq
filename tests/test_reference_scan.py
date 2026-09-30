"""The vectorized scan must agree with the readable reference implementation.

`randseq.core.get_pattern_scores` replaced a string-based scan with integer motif codes and
`np.bincount`. That is a large rewrite of the step every downstream number depends on, so it is
checked against the original rather than trusted: `reference_scan.py` holds the slow version,
written to be obviously correct, and these tests require the two to agree exactly.

This is the only thing standing behind the claim that optimising the scan changed no science.
"""
import os
import random
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(__file__))

from randseq.core import encode_library, get_pattern_scores
from randseq.example_data import get_example_data_dir
from randseq.utils import calculate_log2fc, get_lib_seq_context
from reference_scan import (get_fold_change_values_per_site, get_sites_in_seq,
                            get_sites_scores)

LEFT, RIGHT = "GTCCTAGGTATAATACTAGT", "GTTTTAGAGCTAGAAATAGC"


def assert_same_scores(seqs, fcs, pattern):
    """Both implementations, every column, on the same input."""
    ref = get_sites_scores(
        get_fold_change_values_per_site(get_sites_in_seq(seqs, pattern, no_ori=True), fcs),
        pattern).set_index("motif").sort_index()
    fast = get_pattern_scores(encode_library(seqs), fcs, pattern).set_index("motif").sort_index()

    assert fast.index.tolist() == ref.index.tolist(), f"{pattern}: different motifs found"
    assert fast.num_sequences.tolist() == ref.num_sequences.tolist(), f"{pattern}: support"
    assert np.allclose(fast.fraction_depleted, ref.fraction_depleted, rtol=0, atol=1e-12)
    assert np.allclose(fast.avg_log2fc, ref.avg_log2fc, rtol=0, atol=1e-9)
    return len(ref)


@pytest.fixture(scope="module")
def example_library():
    d = get_example_data_dir()
    counts = pd.read_csv(os.path.join(d, "countsTable.csv.gz"), index_col=0)
    l2 = calculate_log2fc(counts)["JJ1886_T0"]
    return get_lib_seq_context(l2.index, LEFT, RIGHT), l2.values


@pytest.mark.parametrize("pattern", [(6, 0, 0), (2, 1, 3), (4, 4, 3), (3, 6, 4)])
def test_matches_reference_on_example_library(example_library, pattern):
    """Contiguous and spaced patterns, on a real library."""
    seqs, fcs = example_library
    assert assert_same_scores(seqs, fcs, pattern) > 0


def test_matches_reference_on_ragged_and_non_acgt():
    """Sequences of different lengths, one containing an N -- the padding and rejection paths."""
    assert_same_scores(["ACGTACGTTTGACCAAAG", "GGGGCACAAAAGTATCC", "TTGACAANCCGTAC"],
                       [-2.0, 0.5, -3.0], (3, 4, 4))


def test_matches_reference_on_randomized_libraries():
    """Fuzz both implementations against each other on awkward input."""
    rng = random.Random(0)
    total = 0
    for _ in range(25):
        seqs, fcs = [], []
        for _ in range(rng.randint(3, 25)):
            alphabet = "ACGT" if rng.random() < 0.7 else "ACGTN"
            seqs.append("".join(rng.choice(alphabet) for _ in range(rng.randint(16, 45))))
            fcs.append(rng.uniform(-6, 2))
        pattern = rng.choice([(6, 0, 0), (3, 4, 4), (2, 1, 3), (3, 6, 4), (1, 2, 2), (4, 4, 3)])
        total += assert_same_scores(seqs, fcs, pattern)
    assert total > 1000, "fuzzing produced almost no motifs; the generator is probably broken"
