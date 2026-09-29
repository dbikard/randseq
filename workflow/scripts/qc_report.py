"""One report answering: did this run work?

Per-sample checks come from randcount.qc. The run-level check that needs the whole dataset is
replicate agreement, which is what caught contamination in the published screen.
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

from randcount.qc import check_motif_coverage, check_replicates, check_sample, render_text  # noqa: E402

units = [tuple(u) for u in snakemake.params.units]
replicates = snakemake.params.replicates
r_exclude = snakemake.params.r_exclude

sections = []
n_fail = n_warn = 0

for path, (sample, rep) in zip(snakemake.input.stats, units):
    stats = json.loads(Path(path).read_text())
    checks = check_sample(stats)
    checks += check_motif_coverage(stats.get("library", {}).get("members", 0))
    n_fail += sum(c.status == "fail" for c in checks)
    n_warn += sum(c.status == "warn" for c in checks)
    sections.append(render_text(checks, f"{sample} / {rep}"))

# Replicate agreement, per sample, on normalised counts.
counts = {}
for path, (sample, rep) in zip(snakemake.input.counts, units):
    s = pd.read_csv(path, sep="\t").set_index("seq")["umis"]
    counts[(sample, rep)] = s / s.sum()

rep_checks = []
samples = sorted({s for s, _ in counts})
for sample in samples:
    present = [r for r in replicates if (sample, r) in counts]
    if len(present) < 2:
        continue
    a, b = counts[(sample, present[0])], counts[(sample, present[1])]
    joined = pd.concat([a, b], axis=1).fillna(0)
    c = check_replicates(float(joined.corr().iloc[0, 1]), r_exclude=r_exclude)
    c.name = f"replicate_correlation[{sample}]"
    rep_checks.append(c)
    n_fail += c.status == "fail"
    n_warn += c.status == "warn"

if rep_checks:
    sections.append(render_text(rep_checks, "Replicate agreement"))

verdict = ("RUN FAILED QC — do not interpret these results until the failures below are resolved"
           if n_fail else
           "Run passed QC with warnings — read them before interpreting" if n_warn else
           "Run passed QC")

header = ["RandSeq QC report", "=" * 17, "", verdict, "",
          f"{len(units)} sequencing units, {len(samples)} samples, "
          f"{n_fail} failures, {n_warn} warnings", ""]

Path(snakemake.output.report).write_text("\n".join(header) + "\n\n" +
                                         "\n\n".join(sections) + "\n")
print(verdict, file=sys.stderr)
