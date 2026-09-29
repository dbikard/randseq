"""Counts -> log2FC per replicate -> mean across replicates.

log2FC is computed within each replicate first and then averaged, i.e. the geometric mean of the
fold changes. Doing it this way means each replicate is normalised against its own donor, so a
difference in sequencing depth between replicates does not leak into the result.
"""
import sys
from pathlib import Path

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

from randseq.core import calculate_log2fc  # noqa: E402

units = [tuple(u) for u in snakemake.params.units]
donor = snakemake.params.donor
replicates = snakemake.params.replicates
thr = snakemake.params.count_threshold

counts = {}
for path, (sample, rep) in zip(snakemake.input.counts, units):
    df = pd.read_csv(path, sep="\t").set_index("seq")["umis"]
    counts[(sample, rep)] = df

per_rep = []
for rep in replicates:
    cols = {s: c for (s, r), c in counts.items() if r == rep}
    if donor not in cols:
        print(f"replicate {rep}: no donor sample, skipped", file=sys.stderr)
        continue
    wide = pd.DataFrame(cols).fillna(0)
    # Row universe is the donor's: a sequence absent from a recipient is a real zero (it was
    # restricted), not missing data.
    l2 = calculate_log2fc(wide, reference_column=donor, count_threshold=thr)
    if l2 is None or l2.empty:
        print(f"replicate {rep}: empty after count_threshold={thr}", file=sys.stderr)
        continue
    per_rep.append(l2)
    print(f"replicate {rep}: {len(l2)} plasmids, {l2.shape[1]} recipients", file=sys.stderr)

if not per_rep:
    raise SystemExit(f"no replicate produced usable log2FC at count_threshold={thr}. "
                     f"Lower it, or check that the donor sample has enough depth.")

# Mean across replicates on the intersection: a plasmid must be measured in every replicate for
# its mean to mean the same thing as everyone else's.
mean = pd.concat(per_rep).groupby(level=0).mean()
common = set(per_rep[0].index)
for r in per_rep[1:]:
    common &= set(r.index)
mean = mean.loc[sorted(common)]

print(f"\n{len(mean)} plasmids measured in all {len(per_rep)} replicates", file=sys.stderr)
mean.index.name = "seq"
mean.to_csv(snakemake.output.table, sep="\t")
