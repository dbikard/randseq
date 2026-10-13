"""The vectorised fixed-position scan must agree with the readable one, exactly.

Ragged inserts are the case that matters: an offset from the right-hand end is a different
absolute index in every length class, so these libraries deliberately mix lengths, include
sequences too short to reach the deeper offsets, and contain non-ACGT characters.
"""
import os
import numpy as np
import pandas as pd
import pytest

from randseq.core import identify_depleted_motifs_scanning_ends
from randseq.utils import calculate_log2fc
from randseq.example_data import get_example_data_dir
from reference_fixed_scan import reference_fixed_scan

KEY = ["motif", "anchor", "position", "length"]


def assert_same(a, b):
    a = a.sort_values(KEY).reset_index(drop=True)
    b = b.sort_values(KEY).reset_index(drop=True)
    assert len(a) == len(b), f"{len(a)} rows vs {len(b)}"
    if not len(a):
        return
    for col in KEY + ["num_sequences"]:
        assert (a[col].values == b[col].values).all(), f"{col} differs"
    for col in ["avg_log2fc", "fraction_depleted"]:
        assert np.allclose(a[col].values.astype(float), b[col].values.astype(float)), col


@pytest.mark.parametrize("strain", ["B156_R1", "JJ1886_R1", "LMR_503_R1", "EDL933_dRM_R1"])
def test_matches_reference_on_real_libraries(strain):
    d = get_example_data_dir()
    l2 = calculate_log2fc(pd.read_csv(os.path.join(d, "counts_panel.csv.gz"), index_col=0),
                          reference_column="Control_R1", count_threshold=20)
    series = l2[strain]
    assert_same(reference_fixed_scan(series), identify_depleted_motifs_scanning_ends(series))


@pytest.mark.parametrize("seed", range(10))
def test_matches_reference_on_ragged_libraries(seed):
    rng = np.random.default_rng(seed)
    seqs = ["".join(rng.choice(list("ACGTN"), rng.integers(3, 40)))
            for _ in range(int(rng.integers(200, 600)))]
    seqs = list(dict.fromkeys(seqs))
    series = pd.Series(rng.normal(-1, 1.5, len(seqs)), index=seqs)
    kw = dict(min_sequence_support=1, score_thr=0.0)
    assert_same(reference_fixed_scan(series, **kw),
                identify_depleted_motifs_scanning_ends(series, **kw))


def test_right_anchor_pools_across_insert_lengths():
    """The point of anchoring: a motif at the end of every insert is one motif, not one per
    length class."""
    rng = np.random.default_rng(0)
    seqs = ["".join(rng.choice(list("ACGT"), n)) for n in rng.integers(30, 150, 400)]
    seqs = list(dict.fromkeys(seqs))
    carries = rng.random(len(seqs)) < 0.3
    seqs = [s[:-3] + "GTC" if c else s for s, c in zip(seqs, carries)]
    carries = np.array([s.endswith("GTC") for s in seqs])
    fc = rng.normal(0, 0.3, len(seqs))
    fc[carries] = rng.normal(-5, 0.3, carries.sum())
    out = identify_depleted_motifs_scanning_ends(pd.Series(fc, index=seqs))
    hit = out[(out.motif == "GTC") & (out.anchor == "right") & (out.position == 0)]
    assert len(hit) == 1, "the planted motif must be one row, not one per insert length"
    assert int(hit.num_sequences.iloc[0]) == int(carries.sum())
