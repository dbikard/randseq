"""Call restriction motifs per recipient, with threshold-free confidence statistics.

Two things this does beyond calling `find_restricted_motifs`:

1. **It does not die on one bad strain.** A strain with no real motif throws up many
   near-threshold candidates and trips `TooManyCandidatesError`; an RM-deficient strain is the
   *slow* case, not the fast one, and it is exactly what a new lab runs first. Those strains are
   recorded in a skipped table and the run continues.

2. **It reports confidence in units the user already thinks in.** `num_sequences` is set by the
   count threshold, so it is not comparable between experiments. The statistics below are
   threshold-free, so calls from different runs can be compared even when their thresholds differ:

   - `median_log2fc` -- how strong, reads directly as a fold change
   - `escaper_rate`  -- how penetrant: % of motif-carrying plasmids with log2FC > -0.5
   - `n_plasmids`    -- how much evidence, reported *with* the threshold that produced it
   - `backbone_sites`-- how many copies of this motif the vector already carries
"""
import sys

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

from randseq.core import TooManyCandidatesError, find_restricted_motifs, get_patterns  # noqa: E402
from randseq.utils import _get_filter_for_motif  # noqa: E402

ESCAPER_THR = -0.5

left = snakemake.params.left
right = snakemake.params.right
recipients = snakemake.params.recipients
l2fc_thr = snakemake.params.log2fc_thr
score_thr = snakemake.params.score_thr
backbone_path = snakemake.params.backbone

backbone = None
if backbone_path:
    seq = "".join(l.strip() for l in open(backbone_path) if not l.startswith(">")).upper()
    backbone = seq + seq          # circular: catch motifs spanning the origin


def count_in_backbone(motif):
    """How many times this motif occurs in the vector itself.

    A motif already present in the backbone cannot be detected: every plasmid in the library
    carries it, so there is no motif-free baseline to deplete against. The published screen hit
    exactly this and had to change backbone to recover a masked Type I motif. Warning about it
    is the difference between a negative result and a *meaningless* negative result.
    """
    if backbone is None:
        return None
    from randseq.utils import revcomp
    import re
    from randseq.utils import IUPAC_DNA_TO_REGEX
    pat = "".join(IUPAC_DNA_TO_REGEX.get(b, b) for b in motif)
    fwd = len(re.findall(f"(?={pat})", backbone))
    rev_pat = "".join(IUPAC_DNA_TO_REGEX.get(b, b) for b in revcomp(motif))
    rev = len(re.findall(f"(?={rev_pat})", backbone))
    return (fwd + rev) // 2       # doubled sequence


table = pd.read_csv(snakemake.input.table, sep="\t").set_index("seq")

rows, skipped = [], []
for recipient in recipients:
    if recipient not in table.columns:
        skipped.append(dict(sample=recipient, reason="not in log2fc table"))
        continue
    series = table[recipient].dropna()
    print(f"\n### {recipient}: {len(series)} plasmids", flush=True, file=sys.stderr)
    try:
        _fixed, flex = find_restricted_motifs(
            series, left, right,
            flexible_motif_patterns=get_patterns(),
            flexible_motif_log2fc_thr=l2fc_thr,
            flexible_motif_score_thr=score_thr,
        )
    except TooManyCandidatesError as e:
        print(f"  skipped: {e.n_candidates} candidates > limit", file=sys.stderr)
        skipped.append(dict(sample=recipient, reason=f"too many candidates ({e.n_candidates}); "
                                                     f"no confident call. Tighten thresholds to "
                                                     f"inspect this strain."))
        continue
    except Exception as e:                      # one bad strain must not sink a 100-sample run
        print(f"  FAILED: {type(e).__name__}: {e}", file=sys.stderr)
        skipped.append(dict(sample=recipient, reason=f"{type(e).__name__}: {e}"))
        continue

    if flex is None or flex.empty:
        print("  no motif called", file=sys.stderr)
        continue

    for _, r in flex.iterrows():
        motif = r["motif"]
        mask = _get_filter_for_motif(series, motif, left, right)
        m_plus = series[mask]
        rows.append(dict(
            sample=recipient,
            motif=motif,
            pattern=r["pattern"],
            median_log2fc=round(float(m_plus.median()), 4),
            mean_log2fc=round(float(r["avg_log2fc"]), 4),
            escaper_rate=round(float((m_plus > ESCAPER_THR).mean()), 4),
            fraction_depleted=round(float(r["fraction_depleted"]), 4),
            n_plasmids=int(r["num_sequences"]),
            n_plasmids_with_motif=int(mask.sum()),
            count_threshold_note=f"n depends on the count threshold used upstream",
            backbone_sites=count_in_backbone(motif),
        ))
        print(f"  {motif}: median log2FC {m_plus.median():.2f}, "
              f"escapers {(m_plus > ESCAPER_THR).mean():.1%}", file=sys.stderr)

out = pd.DataFrame(rows)
if not out.empty and backbone is not None:
    masked = out[out["backbone_sites"] > 0]
    for _, r in masked.iterrows():
        print(f"WARNING {r['sample']}/{r['motif']}: {r['backbone_sites']} site(s) already in the "
              f"backbone -- detection is suppressed for this motif", file=sys.stderr)

out.to_csv(snakemake.output.motifs, sep="\t", index=False)
pd.DataFrame(skipped).to_csv(snakemake.output.skipped, sep="\t", index=False)
print(f"\n{len(out)} calls across {out['sample'].nunique() if not out.empty else 0} samples; "
      f"{len(skipped)} samples skipped", file=sys.stderr)
