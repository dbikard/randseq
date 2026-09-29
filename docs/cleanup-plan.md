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

## 2. `old_nbs/` — 16 MB, frozen, unused

Two notebooks (`241909_David.ipynb`, `old_analysis_David.ipynb`) kept from before the package
existed. Nothing exports from them and nothing imports them; the phase-1 cleanup confirmed the
only functions they referenced were removable. They are 16 MB of the repository.

**Removal trigger:** none needed — git history keeps them. Delete whenever.

## 3. `index_files/` — 161 KB of generated output

Rendered artifacts from an old docs build. Regenerable.

**Removal trigger:** none. Delete, and gitignore the path.

## 4. The string-based reference scanner

`get_sites_in_seq`, `get_fold_change_values_per_site`, `get_sites_scores` are a second,
slower implementation of the pattern scan. They look like dead weight but they are not: the
notebook uses them to prove the vectorized scan gives identical answers on every pattern of both
example datasets. That equivalence check is the only thing standing behind the claim that
optimising the scan did not change the science.

**Recommendation: keep**, but say so explicitly in the module docstring so the next person does
not delete them as duplication. They are a test fixture that happens to be exported.

**Alternative if the API surface matters more:** move them to `tests/` and stop exporting them.
That keeps the check and removes three names from the public interface.

## 5. Small warts from this branch

- `process_single_flexible_pattern(..., return_pvalues=False)` — a flag added so the pipeline
  can collect every tested p-value for the multiple-testing correction. Cleaner would be to
  always return both and let callers ignore one; the flag exists only to avoid changing an
  exported signature. Fold in whenever the signature can change.
- `randseq.egg-info/` is in the working tree and should be gitignored.

## Suggested order

1. Now, no risk: `old_nbs/`, `index_files/`, gitignore `randseq.egg-info/`.
2. With the next API change: the `return_pvalues` flag.
3. After Bea migrates: `core_v2.py`, `find_restricted_motifs_mp`.
4. Decide separately whether the reference scanner stays public or moves to `tests/`.

Only step 3 is blocked on anyone else.
