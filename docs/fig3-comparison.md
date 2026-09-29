# What the new package changes in Fig 3

Comparison of Bea's current committed table
(`randseq_paper/code/fig3/data/RandSeq_BB2_samples_debuged_v2_mp.csv`, read-only) against the
same 21-strain panel run through the current `eren` branch. Same input
(`RandSeq_BB2_log2fc_mean.csv`), same thresholds (log2FC ≤ −1, score ≥ 0.5), same patterns.

## The short version

**Nothing in the biology changes.** All 11 strains that call motifs call *exactly the same set*.
No motif gained, none lost. `avg_log2fc` — the colour of every heatmap cell — moves by at most
0.10, and mostly by under 0.03.

What changes:

- the table stops containing 344 rows of junk from one strain
- 17 motif labels move to the canonical strand
- support counts rise, often ten-fold, because redundant descriptions of one site were
  destroying each other's evidence
- every call gains a q-value against its own sample's noise, and a shape check

All 30 calls are significant, q = 0 to 9×10⁻⁴². One row is flagged for shape and is worth a
second look.

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

## 3. Support counts move — upwards, often ten-fold

`num_sequences` counts plasmids carrying **exactly one** candidate motif; a plasmid carrying two
cannot say which one depleted it, so it is excluded from both. Several descriptions of the same
site therefore annihilate each other, and a real motif could be left with a handful of plasmids
out of hundreds that contain it.

A filter now removes a candidate whose depletion a stronger candidate already explains (see §4),
before the counting happens. The support that was being destroyed comes back:

| strain | motif | n before | n after |
|---|---|---|---|
| 16171 | `AAACNNNNGTC` | 127 | **847** |
| 16171 | `ACANNNNNATGG` | 42 | **494** |
| 15846 | `ATACNNNNGTG` | 57 | **449** |
| 15846 | `CACNNNNGTAC` | 53 | **421** |
| 15852 | `ATACNNNNGTG` | 91 | **461** |
| 15853 | `CACNNNNGTAC` | 85 | **421** |

**Eleven of thirty are completely unchanged** — every motif in 7083, 7084, 8097, 16175 and 16169.
Those are the strains with few, well-separated motifs, which never had redundant descriptions to
remove. The change touches only what it should.

`avg_log2fc` moves by at most 0.10 throughout, and mostly by under 0.03.

## 4. Why those motifs had so much competition

Every candidate crowding out `CACNNNNGTAC` in LMR_503 is *the same site at the same spacing*,
over-specified:

| competitor | decomposes as |
|---|---|
| `CACNNNNGTAC` | `CAC-NNNN-GTA` + `C` |
| `TACNNNNGTGG` | `C` + `CAC-NNNN-GTA` |
| `CACTNNAGTA` | `CAC-TNNA-GTA` |
| `ACACNNNAGTA` | `A` + `CAC-NNNA-GTA` |

All spacer 4. They differ only by pinning an extra base outside the site or inside the spacer,
and crucially **none of them constrains the final `Y`**. The real site is `CACNNNNGTAY`; drop the
`Y` and roughly half the matches are at uncut positions, so each of these scores just over the
0.5 threshold on borrowed signal. Score tracked "fraction of its plasmids that carry the real
site" almost exactly 1:1 across all 24 candidates.

That the spacing is what matters was checked directly — holding the half-sites fixed and varying
only the spacer:

| spacer | `CAC-Nx-GTA` median log2FC | `AAAG-Nx-GTT` |
|---|---|---|
| 0–3 | −0.03 to −0.06 | −0.00 to −0.04 |
| **4** | **−1.81** | **−3.09** |
| 5–8 | −0.02 to −0.07 | −0.00 to −0.02 |

Baseline over all plasmids is −0.013, so every other spacing is at background. Controls hold:
EDL933 ΔRM is flat at every spacing, and a scrambled 5′ half-site in the same strain is flat too.

The string-based redundancy filter cannot collapse these, and is right not to: the natural
general motif `CACNNNNGTA` also matches uncut sites and scores about half as well, so the rule
"drop the specific only when the general is not significantly worse" correctly declines. The new
filter compares **plasmid sets** instead of strings, and asks whether a candidate still depletes
the plasmids that do *not* carry the stronger one. It introduces no threshold — it re-applies the
score threshold the candidate already had to pass.

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

## 7. New information: significance and shape

Every call now carries a q-value computed against **that sample's own background depletion
rate**, corrected over every motif the scan tested (~750,000). This is information that did not
exist before, and it is what makes a low-n result interpretable.

All 30 calls in the panel are significant, from q = 0 down to q = 9×10⁻⁴². Nothing in Fig 3 is
marginal.

### The background rate is the most diagnostic number about a sample

It varies 20-fold across the panel, and it changes what a *negative* result means:

| strain | background | calls |
|---|---|---|
| 8099 | 0.015 | none |
| 16225 | 0.021 | none |
| 8097 | 0.018 | 2 |
| … | | |
| 16164 | 0.108 | none |
| 16165 | 0.125 | none |
| 8098 | 0.144 | none |
| **16223 (381A)** | **0.335** | 1, plus 13 below threshold |

"No motif found" in 8099 (background 0.015) is a strong negative. The same statement for 8098
(0.144) is much weaker — nearly a seventh of that library is depleted for reasons unrelated to
any motif, so a real weak system could be hidden. Neither was distinguishable before.

For 381A the number explains what was previously a judgement call: a third of its library is
depleted, its one real call sits at q = 3×10⁻⁷⁵, and its other thirteen candidates fall between
q = 0.36 and 0.93. Its library bottleneck is visible in the statistics rather than asserted.

### Shape: one population or two

A real site shifts one population; a motif that is really a diluted version of another is a
mixture of destroyed and untouched plasmids. Sarle's bimodality coefficient separates them at the
conventional 5/9 ≈ 0.556 benchmark — shadows measured 0.647–0.667, every real call 0.19–0.47,
with BREX at 0.236 despite having the broadest distribution of any real call.

**One call in the panel is flagged:** `CAACNNNNNTCGG` in 7083, at 0.569 — marginally over the
line. Its partner `CAATNNNNNTCGG` sits at 0.468, and the two differ only at one position, so the
site is presumably `CAAYNNNNNTCGG` and both halves are real. This looks like a borderline flag
rather than a problem, but it is the one row in Fig 3 worth a second look.

Significance and shape catch different failures and neither subsumes the other: a shadow of a
real site is *highly* significant and bimodal; a chance call is neither.

## What Bea needs to do

1. **Add `try/except TooManyCandidatesError`** to `generate_fig3_motifs.py`, or it dies on 16223
   instead of producing a table.
2. **Re-run `methylome_crossref_Fig3.R`.** It matches motif strings with reverse-complement
   expansion so it should be unaffected, but 17 labels changed and that is worth confirming rather
   than assuming.
3. Expect **IUPAC-merged names to change spelling** in `..._iupac_merged_mm2.csv` even where the
   grouping does not.
4. **The `n=` placeholders in Results §2 can now be filled** — support counts are no longer a
   small residual of a large number, and every call carries a q-value to quote beside it.
5. Correct Results §2 `GAGACC` to `GGTCTC` to match §3.

## Reproducing this

- `work/run_fig3_arm.py` — the panel under whichever `randseq` is on `sys.path`
- `work/compare_fig3.py` — the diff against her committed table
- `work/why_counts_fall2.py`, `work/competitor_scores.py`, `work/explained_by.py` — §3 and §4
- `work/permutation_fdr.py`, `work/analytic_null.py` — the calibration behind §6
- `work/flank_competition.py`, `work/spacer_scan.py`, `work/decompose.py` — section 4

All of it reads her repository and writes nothing to it.
