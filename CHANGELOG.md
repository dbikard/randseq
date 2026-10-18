# Changelog

## 0.1.0 (unreleased)

A rework of motif calling. Every call now carries a significance value, a site is reported
once, and the output is a single table. The notebooks document each step.

### Breaking changes

- `find_restricted_motifs` returns one table, not a `(fixed, flex)` tuple: the called motifs, in
  five columns (`motif`, `num_sequences`, `fraction_depleted`, `avg_log2fc`, `q_value`), sorted
  by q. `all_motifs=True` returns every scored motif, called or not, with all
  the columns and a `reason` for each verdict (`call_reason`). The
  fixed-position scan and its `fixed_motif_*` / `core_motif_score_margin` parameters are gone:
  the scan reads the vector flanks with each insert, so sites formed across the junction are
  found directly.
- `find_restricted_motifs_mp` and `randseq.core_v2` are removed; pass `n_jobs` instead.
- `find_restricted_motifs` parameters renamed, without aliases: `flexible_motif_patterns` →
  `patterns`, `flexible_motif_log2fc_thr` → `log2fc_thr`, `flexible_motif_score_thr` →
  `score_thr`, `flexible_motif_min_support` → `min_support`.
- `calculate_log2fc` median-normalises by default (`normalize="median"`): each sample's median
  log2FC is subtracted, so the median plasmid has fold change 1 (median-ratio normalisation, as in
  DESeq2). `normalize="total"` returns the old values, normalised to each sample's total only.
- `find_restricted_motifs` takes log2FC as given and no longer normalises them; it expects them
  median-normalised, as `calculate_log2fc` returns them. Values computed elsewhere must be
  normalised the same way (`series - series.median()`).
- Default depletion threshold: a plasmid counts as depleted below −0.5 (`log2fc_thr`), relative
  to the sample's median plasmid.
- `plot_motif_analysis` takes the motif table as `motifs` (was `flexible_motifs_df`).
- `max_candidates` (default 5000) is checked before step 4 and raises `TooManyCandidatesError`.
- `_custom_motif_to_regex`, and the motif lookups built on it, raise `ValueError` on an invalid
  motif instead of matching nothing.
- A smaller public API. Public: `find_restricted_motifs`, `TooManyCandidatesError` and
  `call_reason` (`randseq.core`); `calculate_log2fc`, `get_patterns`, `plasmid_background`,
  `canonical_motif`, `revcomp`, `normalise_motif`, `get_lib_seq_context`,
  `get_custom_motif_presence_in_library`, `get_log2fc_distribution_by_motif`, and for analysis
  scripts `generate_tuple_combinations`, `filter_symmetric_tuples` and
  `find_iupac_equivalent_motifs` (`randseq.utils`); `encode_library` (`randseq.coding`); the
  plotting and design functions. The pipeline's steps and helpers are private (leading
  underscore), and the scan of one pattern moved to `randseq.fast`. `score` and
  `get_motif_filter_with_context` are removed: `get_custom_motif_presence_in_library` is the
  motif lookup.

```python
# before
from randseq.core_v2 import find_restricted_motifs_mp
series = calculate_log2fc(counts, reference_column="MFDpir")["JJ1886_T0"]
_fixed, flex = find_restricted_motifs_mp(series, left, right, cpu_count=16,
                                         flexible_motif_log2fc_thr=-1, flexible_motif_score_thr=0.7)
# after
from randseq.core import find_restricted_motifs, TooManyCandidatesError
series = calculate_log2fc(counts, reference_column="MFDpir")["JJ1886_T0"]   # median-normalised
motifs = find_restricted_motifs(series, left, right, n_jobs=16,
                                log2fc_thr=-1, score_thr=0.7)                # the calls
motifs_all = find_restricted_motifs(series, left, right, all_motifs=True)    # and the rest
# log2FC from your own code: normalise them to the median first
motifs = find_restricted_motifs(my_log2fc - my_log2fc.median(), left, right)
```

### Calling

- Significance: each motif's depleted count is tested (binomial) against its carriers' mean
  background, where each plasmid's background is the depletion rate of inserts of similar GC
  content in the same sample (`plasmid_background`). q-values are Benjamini–Hochberg over every
  motif tested.
- A motif is `called` when q < `q_threshold` and the one-sided 95% lower bound on its depleted
  fraction reaches `score_thr`. Uncalled rows are in the `all_motifs=True` table.
- One description per site: three order-independent rules (step 3) and a test of each candidate
  on the plasmids no stronger candidate explains, held to the scan's own cutoff (step 4).
- A degenerate site is reported as one IUPAC motif (e.g. `RTACNNNNGTG` rather than two variant
  rows); `merged_from` lists the rows it replaced.

### New columns

`background` (the carriers' GC-matched depletion rate), `fraction_depleted_lower`, `p_value`,
`q_value`, `called`, `bimodality`, `bimodal` (carriers split between destroyed and untouched),
`merged_from`, `reason`. All but `q_value` are in the `all_motifs=True` table only.

### Other

- IUPAC-aware `revcomp` and `canonical_motif`; `normalise_motif` accepts bracketed classes.
- `randseq.design`: `plan_library` (which motif architectures a library and vector can detect,
  and the library size a site needs), `min_carriers`, and `carrier_fraction` (now exported and
  width-aware for gapped patterns).
- Optional compiled scan: `pip install "randseq[fast]"` (numba); results are identical.
- Removed: the unused attic module.
