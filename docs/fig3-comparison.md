# What the new package changes in Fig 3

Comparison of Bea's current committed table
(`randseq_paper/code/fig3/data/RandSeq_BB2_samples_debuged_v2_mp.csv`, read-only) against the
same 21-strain panel run through the current `eren` branch. Same input
(`RandSeq_BB2_log2fc_mean.csv`), same thresholds (log2FC ≤ −1, score ≥ 0.5), same patterns.

## The short version

**Nothing in the biology changes.** All 11 strains that call motifs call *exactly the same set*.
No motif gained, none lost. `avg_log2fc` — the colour of every heatmap cell — moves by at most
0.12, and mostly by under 0.03.

What changes is that the table stops containing 344 rows of junk, 17 motif labels move to the
canonical strand, and some support counts move for a reason worth understanding.

## 1. 381A stops silently poisoning the table

Her table is 374 rows. **344 of them (92%) come from strain 16223 alone**, and they are not
results: that strain exceeded `max_candidates`, so the old code returned the *unfiltered*
candidate list, which has the same columns as a real result. Those rows then flow into the IUPAC
merge and produce motifs like `SSBNNSVNN`.

The new package raises `TooManyCandidatesError` instead. Excluding 16223, both tables have
exactly **30 rows**.

Nothing published is wrong, because 16223 is dropped from the figure by hand — but the only
thing standing between those rows and a figure was remembering to drop one strain.

## 2. Seventeen of thirty rows are relabelled to the canonical strand

No other change to those rows.

| was | now |
|---|---|
| `GAGACC` | `GGTCTC` |
| `AACNNNNCTTT` | `AAAGNNNNGTT` |
| `CACNNNNGTAT` | `ATACNNNNGTG` |
| `CTANNNNNNNNTAGG` | `CCTANNNNNNNNTAG` |
| `CTANNNNNNNNTGGG` | `CCCANNNNNNNNTAG` |
| `TGGCCGC` | `GCGGCCA` |
| `TAATNNNNNNNGTCG` | `CGACNNNNNNNATTA` |
| `GACNNNNGTTC` / `GACNNNNGTTT` | `GAACNNNNGTC` / `AAACNNNNGTC` |

A motif and its reverse complement are the same site and always carry identical statistics, so
only one should be reported. Which one used to depend on scan order. This also resolves the
manuscript calling one site `GAGACC` (BsaI) in Results §2 and `GGTCTC` (Eco31I) in §3 — they are
isoschizomers.

## 3. Support counts move — and 11 of 30 do not move at all

`num_sequences` counts plasmids carrying **exactly one** candidate motif; a plasmid carrying two
cannot say which one depleted it, so it is excluded from both. That makes the count a property
of the whole candidate set, and how stable it is depends on how much competition a motif has:

| strain | motif | contains | counted | retained | change |
|---|---|---|---|---|---|
| 8097 | `CGACNNNNNNNATTA` | 138 | 138 | **100%** | unchanged |
| 16226 | `AAAGNNNNGTT` | 897 | 895 | **99.8%** | 659 → 895 |
| 15846 | `GGTCTC` | 54 | 47 | 87% | unchanged |
| 15846 | `ATACNNNNGTG` | 484 | 16 | **3.3%** | 57 → 16 |
| 15846 | `CACNNNNGTAC` | 460 | 6 | **1.3%** | 53 → 6 |

Where retention is high the count is a real measure of evidence and it barely moves. Where it is
1–3%, the count is a small residual of a large number, and any change to the candidate set swings
it. **53 and 6 are both residuals of 460; neither is more correct than the other.**

## 4. Why those particular motifs have so much competition

Not what we first assumed. Every motif competing with `CACNNNNGTAC` in LMR_503 turns out to be
*the same site at the same spacing*, over-specified:

| competitor | decomposes as |
|---|---|
| `CACNNNNGTAC` | `CAC-NNNN-GTA` + `C` |
| `TACNNNNGTGG` | `C` + `CAC-NNNN-GTA` |
| `CACTNNAGTA` | `CAC-TNNA-GTA` |
| `CACANNGGTA` | `CAC-ANNG-GTA` |
| `ACACNNNAGTA` | `A` + `CAC-NNNA-GTA` |

All spacer 4. They differ only by pinning one extra base outside the site, or bases inside the
spacer. A plasmid matching one necessarily matches several, so eight descriptions of one site
annihilate each other in the unique-hit step.

That the spacing is what matters was checked directly. Holding the half-sites fixed and varying
only the spacer:

| spacer | `CAC-Nx-GTA` median log2FC | `AAAG-Nx-GTT` median log2FC |
|---|---|---|
| 0–3 | −0.03 to −0.06 | −0.00 to −0.04 |
| **4** | **−1.81** (53% depleted) | **−3.09** (99.6% depleted) |
| 5–8 | −0.02 to −0.07 | −0.00 to −0.02 |

Baseline over all plasmids is −0.013, so every other spacing is at background — a knife-edge, as
Type I biology predicts. Controls hold: EDL933 ΔRM is flat at every spacing, and a scrambled 5′
half-site in the same strain is flat too.

**The redundancy filter cannot fix this as designed.** The natural broad motif is `CACNNNNGTA`,
but the real site is degenerate at the last position (`CACNNNNGTAY`), so `CACNNNNGTA` also matches
uncut sites and its score is roughly halved (0.52 against 1.00 for the specifics). The filter
drops a specific motif only when the general one is *not significantly worse* — and here it
genuinely is worse. The filter is behaving correctly; the design has no way to say "these eight
are one site with a degenerate position".

`find_iupac_equivalent_motifs` is exactly the right tool, and it already exists — but it runs
downstream in Bea's script, **after** the rescoring. Merging before rescoring would make these one
candidate instead of eight competitors. That is a design change with biological judgement in it
and has not been made.

## 5. The old numbers were not reproducible anyway

Bea's `fig3/README.md` records that the old search was non-deterministic across runs. Running the
old package ourselves on the same input confirms it: of 30 shared calls, **2 differ, and by ~1.5×**

| strain | motif | her run | our run of the same old code |
|---|---|---|---|
| 16171 | `GAACNNNNGTC` | 61 | 99 |
| 16171 | `AAACNNNNGTC` | 127 | 191 |

So part of the old-vs-new difference in section 3 is our changes and part is the old code
disagreeing with itself. The new pipeline is deterministic and has a test that pins it.

## 6. Runtime

The same 21-strain panel: **1 h 45 min** on the old package, minutes on the new one. An
RM-deficient strain used to be the *slow* case — nothing real to find means hundreds of
near-threshold candidates entering a quadratic step — which is exactly what a new lab runs first.

## What Bea needs to do

1. **Add `try/except TooManyCandidatesError`** to `generate_fig3_motifs.py`, or it dies on 16223
   instead of producing a table.
2. **Re-run `methylome_crossref_Fig3.R`.** It matches motif strings with reverse-complement
   expansion so it should be unaffected, but 17 labels changed and that is worth confirming rather
   than assuming.
3. Expect **IUPAC-merged names to change spelling** in `..._iupac_merged_mm2.csv` even where the
   grouping does not.
4. **Do not quote a bare `n=`.** The Results §2 placeholders need the threshold stated with them,
   and for bipartite motifs the retention (`6 of 460`) matters more than the count.
5. Correct Results §2 `GAGACC` to `GGTCTC` to match §3.

## Reproducing this

- `work/run_fig3_arm.py` — the panel under whichever `randseq` is on `sys.path`
- `work/compare_fig3.py` — the diff against her committed table
- `work/why_counts_fall2.py` — the retention analysis
- `work/flank_competition.py`, `work/spacer_scan.py`, `work/decompose.py` — section 4

All of it reads her repository and writes nothing to it.
