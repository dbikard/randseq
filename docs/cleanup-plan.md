# Removing the legacy surface

DB, 2026-09-29: *"I eventually want to ship code that is as clean as possible and without
deprecated code paths."*

Everything deprecated in the tree exists for one reason — to let Bea's existing scripts keep
running across the handoff. So this is a sequencing problem, not a cleanup that can just be done:
remove it before she migrates and her figure generation breaks; leave it indefinitely and the
package ships with two ways to do everything.

Each item below therefore has a **removal trigger** rather than a vague "future release".

## 1. `randseq/core_v2.py` and `find_restricted_motifs_mp`

**Why they exist.** `core_v2` was merged into `core` in `b758f8f`. The shim re-exports the merged
implementations; `find_restricted_motifs_mp` forwards `cpu_count` to `n_jobs` and warns.

**Who depends on them** — 8 scripts, all needing the same one-line change:

| file | change |
|---|---|
| `randseq_paper/code/fig3/generate_fig3_motifs.py` | import, and `cpu_count=` → `n_jobs=` |
| `randseq_paper/code/fig2/generate_fig2e_grid_search.py` | same |
| `BB3/results/motif_detection/run_bb3_motif_detection.py` | same |
| `…/run_bb3_extra_strains_motif_detection.py` | same |
| `…/run_bb3_16223_motif_detection.py` | same |
| `…/run_bb3_16170_motif_detection.py` | same |
| `…/run_bb2_16223_threshold_sweep.py` | same |
| `…/run_bb3_pBB027_thr10_motif_detection.py` | same |

```python
# before
from randseq.core_v2 import find_restricted_motifs_mp
... find_restricted_motifs_mp(series, left, right, cpu_count=16)

# after
from randseq.core import find_restricted_motifs
... find_restricted_motifs(series, left, right, n_jobs=16)
```

**Removal trigger:** once Bea confirms her scripts are migrated. Not before — and note her scripts
must be re-run anyway, because the results change (see `fig3-comparison.md`), so the migration is
free in the sense that it happens during work she has to do regardless.

## 2. `old_nbs/` — done

Deleted; git history keeps it, and the path is gitignored. Nothing referenced it.

## 3. `index_files/` — keep it, and it had to be *added*

An earlier version of this plan said to delete it as regenerable build output. That was wrong.
`README.md` embeds `index_files/figure-commonmark/cell-4-output-1.png`, so it is the only image on
the GitHub front page — and it was **untracked**, which means that image had been broken on GitHub.

It is regenerable (`nbdev-readme` rewrites it from `index.ipynb`), but it has to be committed for
the README to render. Now tracked. Regenerate and re-commit it whenever `index.ipynb`'s output
changes; the current one is up to date, showing canonical-strand motif names.

## 4. The string-based reference scanner

`get_sites_in_seq`, `get_fold_change_values_per_site`, `get_sites_scores` are a second,
slower implementation of the pattern scan. They look like dead weight but they are not: the
notebook uses them to prove the vectorized scan gives identical answers on every pattern of both
example datasets. That equivalence check is the only thing standing behind the claim that
optimising the scan did not change the science.

**Decided: keep, and said so in the notebook.** The markdown above `encode_library` in
`00_core.ipynb` now states that the slow path is deliberate — it is the only independent check
that vectorising the scan changed no number, it is written to be obviously correct rather than
fast, and deleting it makes the equivalence check untestable. They stay exported so the check can
run as an ordinary notebook cell under `nbdev-test`.

## 5. Small warts from this branch

- ~~`process_single_flexible_pattern(..., return_pvalues=False)`~~ — **done.** The function now
  always returns `(kept, tested_pvalues)`. Nothing ever called it with the flag off, so this cost
  nothing; the docstring explains why both are returned.
- `randseq.egg-info/` — already gitignored and untracked; nothing to do.

## Suggested order

1. ~~Now, no risk: `old_nbs/`, `index_files/`, gitignore `randseq.egg-info/`.~~ **Done** — with one
   correction: `index_files/` was kept and committed rather than deleted (see §3).
2. ~~With the next API change: the `return_pvalues` flag.~~ **Done.**
3. After Bea migrates: `core_v2.py`, `find_restricted_motifs_mp`.
4. ~~Decide whether the reference scanner stays public or moves to `tests/`.~~ **Decided: stays.**

**Only step 3 is left, and it is blocked on Bea.** Everything else in this plan is done.
