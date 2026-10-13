# Changelog

## Unreleased — for Bea, before redoing the figures

Everything below was measured against the real BB2 panel, not just the example data.
`docs/fig3-comparison.md` has the per-strain detail and how to reproduce it.

**The short version: the biology does not change.** All 11 strains that call motifs call exactly
the same set — none gained, none lost — and `avg_log2fc`, the colour of every Fig 3 heatmap cell,
moves by at most 0.10. What changes is that the table stops containing junk, motif labels become
consistent, support counts stop being destroyed, and every call now carries a significance value.

### Your scripts keep working, but please migrate them

`core_v2` was merged into `core`. `randseq/core_v2.py` is a shim and
`find_restricted_motifs_mp` forwards `cpu_count` to `n_jobs`, both with a `DeprecationWarning`,
so nothing breaks today. Eight of your scripts use them, each needing the same one-line change:

```python
# before
from randseq.core_v2 import find_restricted_motifs_mp
... find_restricted_motifs_mp(series, left, right, cpu_count=16)
# after
from randseq.core import find_restricted_motifs
... find_restricted_motifs(series, left, right, n_jobs=16)
```

`docs/cleanup-plan.md` lists the files. The shim is removed once you confirm you have migrated —
and since your scripts have to be re-run anyway, the change is effectively free.

**Worth doing regardless:** wrap the per-strain call in `try/except TooManyCandidatesError`.
16223 no longer trips it (see 9 below), but another noisy strain could, and without the guard
`generate_fig3_motifs.py` would stop instead of producing a table.

---

### 1. 381A no longer contaminates the output

Of the 374 rows in your current Fig 3 table, **344 (92%) come from strain 16223 alone** — that
strain exceeded `max_candidates`, so the old code returned its *unfiltered candidate list*, which
has the same columns as a real result. Those rows flow into the IUPAC merge and become motifs
like `SSBNNSVNN`. Nothing published is wrong because 16223 is dropped from the figure by hand,
but the only thing preventing it was remembering to.

It now raises rather than returning junk. Excluding that strain, both tables have exactly 30 rows.
(Since 9 below, 16223 no longer reaches the limit at all.)

### 2. Motifs are reported on one strand — labels change, numbers do not

A motif and its reverse complement are the same site and always carry identical statistics, so
only one should be reported. Which one used to depend on scan order. The rule is now REBASE's
spelling, falling back to alphabetically first. **17 of 30 rows are relabelled:**

| was | now |
|---|---|
| `GAGACC` | `GGTCTC` |
| `AACNNNNCTTT` | `AAAGNNNNGTT` |
| `CACNNNNGTAT` | `ATACNNNNGTG` |
| `CTANNNNNNNNTAGG` | `CCTANNNNNNNNTAG` |
| `TGGCCGC` | `GCGGCCA` |
| `TAATNNNNNNNGTCG` | `CGACNNNNNNNATTA` |
| `GACNNNNGTTC` / `GACNNNNGTTT` | `GAACNNNNGTC` / `AAACNNNNGTC` |

This also settles a contradiction already in the manuscript: Results §2 calls the JJ1886 site
"the low-abundance **BsaI** motif (**GAGACC**)" and §3 calls the ST131 site
"5′-**GGTCTC**-3′ … the widespread ST131 **Eco31I**". Those are isoschizomers — one site, two
spellings. **§2 should be corrected to `GGTCTC`.**

`methylome_crossref_Fig3.R` matches motif strings with reverse-complement expansion so it should
be unaffected, but 17 labels changed and that is worth confirming rather than assuming. Expect
IUPAC-merged names in `..._iupac_merged_mm2.csv` to change spelling even where the grouping does
not.

### 3. Support counts go up, often ten-fold

`num_sequences` counts plasmids carrying **exactly one** candidate motif — a plasmid carrying two
cannot say which one depleted it, so it is excluded from both. Several descriptions of the *same*
site therefore destroyed each other's evidence, leaving a real motif with a handful of plasmids
out of hundreds that contain it.

A new filter removes a candidate whose depletion a stronger candidate already explains, before
the counting. The evidence comes back:

| strain | motif | before | after |
|---|---|---|---|
| 16171 | `AAACNNNNGTC` | 127 | **847** |
| 16171 | `ACANNNNNATGG` | 42 | **494** |
| 15846 | `ATACNNNNGTG` | 57 | **449** |
| 15846 | `CACNNNNGTAC` | 53 | **421** |
| 15852 | `ATACNNNNGTG` | 91 | **461** |

**Eleven of thirty are unchanged** — every motif in 7083, 7084, 8097, 16175 and 16169, the
strains with few well-separated motifs that never had redundant descriptions to remove.

The motifs doing the crowding turned out to be the same site at the same spacing, over-specified:
in LMR_503 the site is `CACNNNNGTAY`, and `TACNNNNGTGG`, `CACTNNAGTA`, `ACACNNNAGTA` and fifteen
others all decompose to `CAC-NNNN-GTA` with an extra base pinned somewhere that does not matter
and the `Y` left free. Holding the half-sites fixed and varying only the spacer gives median
log2FC −1.81 at spacer 4 against −0.03 at every other spacing, with an RM-knockout strain flat
throughout as a control.

**The `n=X` placeholders in Results §2 can now be filled**, and each call has a q-value to quote
beside it.

### 4. Every call now has a significance value

Under the null a motif's depleted count is `Binomial(n, p)` where `p` is **that sample's own
background depletion rate**, corrected across every motif the scan tested (~750,000). Calibrated
against a permutation test that strips the confident calls and reshuffles the rest: 0.60 expected
spurious calls against 0.8 measured at n=3–4 on the smaller library, and 0.00 against 0.0 on the
larger.

All 30 calls in Fig 3 are significant, q from 0 to 9×10⁻⁴². Nothing in that figure is marginal.

**The background rate is new diagnostic information and it changes what a negative result means.**
It varies 20-fold across the panel:

| strain | background | calls |
|---|---|---|
| 8099 | 0.015 | none |
| 8098 | 0.144 | none |
| **16223 (381A)** | **0.335** | 1, plus 13 below threshold |

"No motif found" in 8099 is a strong negative. The same words for 8098 are much weaker — a
seventh of that library is depleted for unrelated reasons, so a weak system could be hiding
there. Those two statements were previously indistinguishable, and nine strains in Fig 3 report
no calls.

For 381A it turns a judgement call into a measurement: a third of its library is depleted, its one
real call sits at q = 3×10⁻⁷⁵, and its other thirteen candidates fall between q = 0.36 and 0.93.

`min_support` default moves 3 → 5, set by that permutation test rather than by taste: motifs with
3–4 supporting plasmids arise by chance about once per scan on the smaller library and never at 5
or above. It cannot go higher than 6 without deleting BsaI on JJ1886, which has 6 there.

### 5. Each call is checked for shape

A real site shifts one population of plasmids. A motif that is really a diluted version of another
is a *mixture* — some plasmids destroyed, the rest untouched — and is bimodal. Sarle's coefficient
separates them at the conventional 5/9 ≈ 0.556 benchmark: shadows measured 0.647–0.667, every real
call 0.19–0.47, with BREX at 0.236 despite having the broadest distribution of any real call.

**One row in Fig 3 is flagged:** `CAACNNNNNTCGG` in 7083 at 0.569, marginally over the line. Its
partner `CAATNNNNNTCGG` is at 0.468 and they differ at one position, so the site is presumably
`CAAYNNNNNTCGG` and both halves are real. Probably a borderline flag rather than a problem, but it
is the one row worth a second look.

### 6. The search is deterministic, and much faster

Your `code/fig3/README.md` records the old search as non-deterministic across runs. Confirmed:
running the old package twice on the same input, 2 of 30 calls differ by ~1.5× (`AAACNNNNGTC`:
127 vs 191). So part of the difference in §3 is our changes and part is the old code disagreeing
with itself — the old `n` was not a stable quantity.

The new pipeline is deterministic and `tests/test_reference_results.py` pins it across three hash
seeds. The 21-strain panel went from **1 h 45 min to a few minutes**; an RM-deficient strain used
to be the *slow* case, which is exactly what a new user runs first.

---

### 7. A site and its reverse complement are one row, from the scan onwards

Previously the scan emitted `GAGACC` *and* `GGTCTC` — one physical site, two rows carrying
identical statistics — and they were reconciled later by the strand rule. They are now folded
onto a single code inside `motif_occurrences`, so the scan emits one row per plasmid per site.

**No reported result changes.** Verified by running 8 configurations (7 panel strains plus
JJ1886 from the older run) before and after: every motif, pattern, support count, depletion
fraction and `avg_log2fc` is byte-identical.

Whether the two spellings *can* be folded is a property of the pattern, not a hard-coded list: a
motif of `(d1, spacer, d2)` reverse-complements into one of `(d2, spacer, d1)`, so they share a
code space only when `spacer == 0` or `d1 == d2`. For the other patterns there was never anything
to fold — scanning `(4, 4, 3)` on both strands already finds each site exactly once, which is why
`get_patterns()` drops `(3, 4, 4)`. Folding those would merge *unrelated* motifs, so
`rc_code_table` refuses them.

Scan runtime is unchanged (8.7 s → 8.8 s over 43 patterns on 30,108 sequences), but the scored
motif space drops from 750,592 rows to 462,624 (−38%). `canonical_motif` still runs at reporting:
folding picks the smaller code, which is `GAGACC`, and the name the literature uses is `GGTCTC`.

### 8. Optional compiled scan: `pip install "randseq[fast]"`

With numba installed, the pattern scan runs as one compiled pass (`randseq.fast`) instead of the
numpy steps. **Results are identical**: on the 21-strain BB2 panel every row — motif, pattern,
support, call — matches the numpy run, and `avg_log2fc` differs by at most 5e-15 (summation
order). Total panel time 279 s -> 67 s on 8 cores; the scan alone is ~13x faster. The first call
in a process compiles for ~4 s, cached afterwards. Without numba nothing changes, and
`RANDSEQ_NO_NUMBA=1` forces the numpy path when it is installed. `tests/test_fast_scan.py` holds
the kernel to the numpy path on all 43 shipped patterns.

### 9. Step 4 tests each motif against everything stronger at once

Step 4 drops a candidate whose depletion is borrowed from stronger ones. It used to test against
each stronger candidate in turn; it now tests against the **union** of all stronger survivors. A
description that overlaps several co-occurring sites could keep a depleted residue against each
of them separately while having no plasmid of its own. Ties in strength are now broken by support
rather than mean log2FC.

On the 21-strain BB2 panel:

- **19 strains are identical**, motifs and support counts, including LMR_503's four co-occurring
  motifs.
- **15853 loses `ACGTNNNNGTTG`**, a more specific shadow, and the panel now matches your published
  Fig 3 exactly: 30 of 30 calls, none extra (it was 31 with one extra).
- **16223 (381A) now produces a result** instead of raising: 262 of its 366 candidates are
  explained away, leaving it under `max_candidates`. It calls `CCTNAGG` (584 plasmids,
  q = 4e-101), clearly real, plus six calls at q 0.02–0.04 that sit on the noise edge. Treat those
  six with suspicion.

### 10. q-values: each test counted once

`annotate_significance` added the reported motifs' p-values to the pool of all tested p-values,
which already contained them. That was not merely a slightly larger denominator: a duplicated
p-value occupies two ranks, and Benjamini-Hochberg's step-up gives both the better rank, so the
strongest motifs' q-values were up to **2x too small**. Now each test is counted once. On the
21-strain panel no call changes; q rises by up to 2x at q ~ 1e-200 and by ~0.3% near the 0.05 cut.

### 11. Defaults now match the paper; mirror patterns searched once; ~2x faster

- **`find_restricted_motifs` defaults** were 2 patterns and score >= 0.7. They are now
  `get_patterns()` -- the paper's 43 -- and score >= 0.5, the settings the paper and Fig 3 use.
  Your scripts pass both explicitly, so they are unaffected.
- **A pattern and its mirror** (`(4, 4, 3)` and `(3, 4, 4)`) find the same sites, since the scan
  reads both strands. Each pattern is now rewritten to one canonical form and searched once.
- **`get_patterns()` takes arguments** to widen the search (halves up to 5 bases, spacer up to 9,
  contiguous 4-8). Wider sets currently add redundant calls at the default thresholds, so the
  default is unchanged.
- **Speed, no result change:** the 21-strain panel runs in 31 s instead of 72 s with numba, every
  row and q-value identical. A library can be encoded once and passed to every strain screened
  with it (`encoded_library=`).

### Running the tests

```bash
~/miniconda3/envs/randseq/bin/python -m pytest tests/ -q
```

Note the two example datasets are both JJ1886 (Lab.ID 12049 *is* JJ1886) from different runs, so
they are a cross-run replicate. Before the canonical-strand rule they reported every motif on
opposite strands; the test asserts they now agree.
