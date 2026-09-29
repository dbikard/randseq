# Changelog

## Unreleased — for Bea, before redoing the figures

Three changes. Two of them move numbers or labels in Fig 3; the third stops a failure that has
been silently producing junk. Everything below was checked against the real BB2 panel
(`RandSeq_BB2_log2fc_mean.csv`, the same input `code/fig3/generate_fig3_motifs.py` uses), not
just the example data.

### Your scripts keep working

`core_v2` was merged into `core` (commit `b758f8f`). `randseq/core_v2.py` is now a shim that
re-exports the merged implementations, so
`from randseq.core_v2 import find_restricted_motifs_mp` still imports. It emits a
`DeprecationWarning` and will be removed later — please move to `from randseq.core import ...`
when convenient, but nothing breaks today.

### 1. Motifs are reported on a canonical strand — **labels change, numbers do not**

The library is scanned on both strands, so a motif and its reverse complement always get
identical statistics. They are one answer, and which of the two got printed was decided by
whatever order the candidates happened to be in. That is why you saw `GAGACC` in one run and
`GGTCTC` in another, and `TGGCCGC` vs `GCGGCCA`.

Now the strand is chosen by rule: **use the spelling REBASE uses**; if REBASE knows neither
strand (most novel motifs), use the alphabetically first. The REBASE Gold Standard list
(n=983, downloaded 2025-12-03) ships with the package as
`randseq/data/rebase_gold_standard_motifs.txt`.

Checked on the 21-strain panel: **all 30 motif rows still match one-to-one, `num_sequences` is
identical, and `avg_log2fc` differs by exactly 0.0.** Only labels move — 15 of 30 rows:

| strain | was | now |
|---|---|---|
| 15846, 15852, 15853, 16226 | `GAGACC` | `GGTCTC` |
| 15846, 15852, 15853, 16226 | `AACNNNNCTTT` | `AAAGNNNNGTT` |
| 15846, 15852, 15853 | `CACNNNNGTAT` | `ATACNNNNGTG` |
| 16171 | `GACNNNNGTTT` | `AAACNNNNGTC` |
| 16171 | `GACNNNNGTTC` | `GAACNNNNGTC` |
| 16172 | `CTANNNNNNNNTGGG` | `CCCANNNNNNNNTAG` |
| 16172 | `CTANNNNNNNNTAGG` | `CCTANNNNNNNNTAG` |

**This affects two things downstream.** `methylome_crossref_Fig3.R` matches motif strings after
reverse-complement expansion, so it should be unaffected — but worth re-running and re-reading
the near-miss list. And the IUPAC merge (`find_iupac_equivalent_motifs`) groups by string
similarity, so merged motif names in `..._iupac_merged_mm2.csv` will change spelling even where
the grouping does not.

**It also fixes a contradiction already in the manuscript.** Results §2 calls the JJ1886 site
"the low-abundance **BsaI** motif (**GAGACC**)"; Results §3 calls the ST131 site
"5′-**GGTCTC**-3′ ... the widespread ST131 **Eco31I**". BsaI and Eco31I are isoschizomers —
`GGTCTC(1/5)` — so that is one site written two ways in one paper. With this change the tool
says `GGTCTC` everywhere, and §2 should be corrected to match.

Fixed-position motifs are deliberately **not** canonicalised: that table records a motif *at an
offset*, and an offset only means something on one strand.

### 2. Too many candidates is now an error — **381A / 16223 will now fail instead of returning junk**

When the flexible search produces more than `max_candidates` (default 150), the unique-hit
rescoring is O(N²) and will not finish. The old behaviour was to warn and return the
**unfiltered candidate list**, which has exactly the same columns as a real result.

In your current Fig 3 table that is what happened to 16223 (381A): of the 374 rows,
**344 are unfiltered candidates from that one strain.** They then flow into the IUPAC merge and
produce motifs like `SSBNNSVNN` and `BNCCNNAGGN`. Nothing published is wrong, because 381A is
dropped from the figure by hand — but the only thing preventing it is remembering to drop it.

Now it raises `TooManyCandidatesError`, naming the strain's candidate count and suggesting
which threshold to tighten. The unfiltered candidates are still available on the exception as
`.candidates` if you want to look at them.

**Practical effect on `generate_fig3_motifs.py`:** it will now raise on 16223 unless you catch
it. The output becomes **30 rows over 20 strains** — which is the real Fig 3 content, with the
junk excluded by the tool rather than by hand. Suggested change: wrap the per-strain call in
`try/except TooManyCandidatesError`, log the strain, and carry on. `work/fig3_motifs.py` in the
`randseq` working directory does exactly this if you want a copy.

### 3. `num_sequences` is larger than in your April/September runs — the rescoring revert

This one predates these changes (it landed in `b758f8f`) but has not been reported to you yet,
and it moves numbers you may be about to quote.

The `core_v2` unique-hit rescoring change was reverted, because it lost `ATACNNNNGTG` on
JJ1886. Rescoring is now done over the full log2FC series again. `avg_log2fc` barely moves
(max |Δ| = 0.026 across the panel), but support counts go up substantially:

| strain | motif (canonical) | now | your reference |
|---|---|---|---|
| 16226 | `AAAGNNNNGTT` | 895 | 659 |
| 15852 | `AAAGNNNNGTT` | 850 | 629 |
| 15846 | `AAAGNNNNGTT` | 845 | 627 |
| 15853 | `AAAGNNNNGTT` | 845 | 626 |
| 16172 | `CCTANNNNNNNNTAG` | 449 | 404 |
| 16172 | `CCCANNNNNNNNTAG` | 432 | 398 |
| 16171 | `AAACNNNNGTC` | 341 | 127 |
| 16171 | `GAACNNNNGTC` | 99 | 61 |

`AAACNNNNGTC` at 341 vs 127 is a 2.7× difference. **Results §2 still has "n=X, n=X and n=X
plasmids for the three JJ1886 motifs" as placeholders** — please fill those from a current run,
not from an older table.

Which of the two counts is *correct* is still open. The revert was chosen because the
alternative demonstrably lost a real motif, but that is an argument against one option rather
than a positive case for the other. Settling it needs a case where the right answer is known by
construction; that is the synthetic fixture in `PLAN.md` phase 0.2, not yet built.

### 4. The search is deterministic — and now tested

Your `code/fig3/README.md` records that the search was not reproducible across runs: under
`PYTHONHASHSEED` 0/1/2 you saw the same motif on either strand and `num_sequences` swinging
61/77/80 and 85/85/33.

**That does not reproduce on this version.** We ran the full 21-strain panel at three hash
seeds: byte-identical output every time. The vectorized numpy scan from `b758f8f` replaced the
set/dict iteration the old scan depended on. `tests/test_reference_results.py` now pins this,
so it cannot regress silently.

### Regression test

`tests/test_reference_results.py` asserts the exact motif tables for both example datasets and
checks determinism across three hash seeds. Run it with the env's own pytest:

```bash
~/miniconda3/envs/randseq/bin/python -m pytest tests/test_reference_results.py -q
```

Note the two example datasets are both JJ1886 (Lab.ID 12049 *is* JJ1886) from different runs —
so they are a genuine cross-run replicate. Before canonicalisation they reported every motif on
opposite strands; they now agree, which the test asserts.
