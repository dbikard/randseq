# Changelog

## 0.1.0 (unreleased)

A rework of motif calling. Every call now carries a significance value and an effect size, a
site is reported once, and the output is a single table. The notebooks document each step.

### Breaking changes

- `find_restricted_motifs` returns one table, not a `(fixed, flex)` tuple: the called motifs, in
  six columns (`motif`, `num_sequences`, `fraction_depleted`, `excess_depleted`, `avg_log2fc`,
  `q_value`), sorted by q. `all_motifs=True` returns every scored motif, called or not, with all
  the columns and a `reason` for each verdict (`call_reason`). The
  fixed-position scan and its `fixed_motif_*` / `core_motif_score_margin` parameters are gone:
  the flexible scan reads the vector flanks with each insert, so sites formed across the
  junction are found directly.
- `find_restricted_motifs_mp` and `randseq.core_v2` are removed; pass `n_jobs` instead.
- Defaults: each sample is centred on its median (`center=True`) and a plasmid counts as
  depleted below −0.5 (`flexible_motif_log2fc_thr`).
- `plot_motif_analysis` takes the motif table as `motifs` (was `flexible_motifs_df`).
- `max_candidates` (default 5000) is checked before step 4 and raises `TooManyCandidatesError`.
- `get_motif_filter_with_context` and `_custom_motif_to_regex` raise `ValueError` on an invalid
  motif instead of matching nothing.

```python
# before
from randseq.core_v2 import find_restricted_motifs_mp
_fixed, flex = find_restricted_motifs_mp(series, left, right, cpu_count=16)
# after
from randseq.core import find_restricted_motifs, TooManyCandidatesError
motifs = find_restricted_motifs(series, left, right, n_jobs=16)              # the calls
motifs_all = find_restricted_motifs(series, left, right, all_motifs=True)    # and the rest
```

### Calling

- Significance: each motif's depleted count is tested (binomial) against its carriers' mean
  background, where each plasmid's background is the depletion rate of inserts of similar GC
  content in the same sample (`plasmid_background`). q-values are Benjamini–Hochberg over every
  motif tested.
- A motif is `called` when q < `q_threshold` and the one-sided 95% lower bound on its depleted
  fraction reaches `flexible_motif_score_thr`. Uncalled rows are in the `all_motifs=True` table.
- One description per site: three order-independent rules (step 3) and a test of each candidate
  on the plasmids no stronger candidate explains, held to the scan's own cutoff (step 4).
- A degenerate site is reported as one IUPAC motif (e.g. `RTACNNNNGTG` rather than two variant
  rows); `merged_from` lists the rows it replaced.

### New columns

`background`, `excess_depleted` (fraction of carriers depleted beyond background),
`fraction_depleted_lower`, `p_value`, `q_value`, `called`, `bimodality`, `bimodal` (carriers
split between destroyed and untouched), `merged_from`, `reason`. All but `excess_depleted` and
`q_value` are in the `all_motifs=True` table only.

### Other

- IUPAC-aware `revcomp` and `canonical_motif`; `normalise_motif` accepts bracketed classes.
- Optional compiled scan: `pip install "randseq[fast]"` (numba); results are identical.
- `motif_presence_matrix` fills every column when two motifs share a code.
- Removed: the unused attic module.
