"""A readable fixed-position scan, kept outside the package as an independent check.

Nothing imports this from `randseq`; its only job is to be obviously correct so that
`tests/test_fixed_scan.py` can hold the vectorised implementation to it.
"""
import numpy as np
import pandas as pd
from randseq.utils import score


def reference_fixed_scan(log2fc_series, scan_depth_k=6, max_motif_length=4,
                         depletion_threshold=-1, score_thr=0.85, min_sequence_support=5):
    """Obviously-correct version: plain Python, one dict of lists, no arithmetic tricks."""
    acc = {}
    for seq, fc in log2fc_series.items():
        if not isinstance(seq, str) or pd.isna(fc):
            continue
        for L in range(1, max_motif_length + 1):
            for k in range(scan_depth_k + 1):
                if k + L <= len(seq):                      # from the left
                    m = seq[k:k + L]
                    if set(m) <= set("ACGT"):
                        acc.setdefault((m, "left", k, L), []).append(fc)
                if k + L <= len(seq):                      # from the right
                    end = len(seq) - k
                    m = seq[end - L:end]
                    if set(m) <= set("ACGT"):
                        acc.setdefault((m, "right", k, L), []).append(fc)
    rows = [dict(motif=m, anchor=a, position=p, length=L, num_sequences=len(v),
                 avg_log2fc=float(np.mean(v)), fraction_depleted=score(v, depletion_threshold))
            for (m, a, p, L), v in acc.items() if len(v) > min_sequence_support]
    df = pd.DataFrame(rows, columns=["motif", "anchor", "position", "length",
                                     "num_sequences", "avg_log2fc", "fraction_depleted"])
    return df[df.fraction_depleted >= score_thr] if len(df) else df
