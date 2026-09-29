# Plan: clean up randseq and rewrite the core notebook as a walkthrough

Written 2026-09-28, to be picked up in a later session.

## Where things stand

- Branch `eren` (pushed, CI green) = Bea's `bea` branch + two commits:
  - `b758f8f` Merge core_v2 into core; vectorized flexible-motif scan
  - `1854343` Upgrade to nbdev 3
- One `core` module. Kept from Bea's core_v2: the two-way `filter_redundant_patterns`, precomputed filters in
  `update_motif_scores_from_unique_hits`, and `max_candidates`. Reverted: unique-hit rescoring is done on the
  **full** log2fc series, not only on sequences without fixed-position motifs.
- Fast scan: `encode_library`, `decode_motif_codes`, `get_pattern_scores` (numpy, motif = integer over its defined
  bases, stats via `np.bincount`). Identical to the string-based reference on all 43 patterns of both example
  datasets; the 12049 search takes ~40 s and ~600 MB on one process.
- nbdev 3.3.24 in the `randseq` conda env: commands use hyphens (`nbdev-export`, `nbdev-test --n-workers 1 --save`,
  `nbdev-clean`), config lives in `pyproject.toml`.

## Working rules (important on this machine)

- WSL has 15 GB RAM and crashed twice from out-of-memory. Run anything heavy inside a memory cap:
  `systemd-run --user --scope -q -p MemoryMax=4G -p MemorySwapMax=0 <command>`
- Use the env's executables: `~/anaconda3/envs/randseq/bin/python`, `.../nbdev-test`, etc. The base Anaconda Python
  has a broken numpy/pyarrow combination.
- Run notebooks one at a time: `nbdev-test --n-workers 1 --save` executes them and saves outputs.

## Reference results (must not change in phases 1 and 3)

Default settings of the notebook examples.

JJ1886 (`countsTable.csv`, column `JJ1886_T0`, reference `MFDpir`, flanks `GTCCTAGGTATAATACTAGT` /
`GTTTTAGAGCTAGAAATAGC`, patterns `[(6,0,0),(7,0,0),(4,2,4),(4,4,3),(4,3,4),(4,4,4)]`, score threshold 0.7):

| motif | fraction_depleted | num_sequences | avg_log2fc |
|---|---|---|---|
| GAGACC | 1.0 | 6 | -3.750094 |
| GTACNNNNGTG | 1.0 | 6 | -4.922717 |
| ATACNNNNGTG | 1.0 | 146 | -4.744952 |
| AAAGNNNNGTT | 1.0 | 54 | -5.018200 |

Fixed-position motifs: GTG @0 (187 seqs), AAAG @12 (42 seqs).

12049 (`counts_12049.csv.gz`, column `12049_R1`, reference `Control`, flanks `GTCTAGGGCGGCGGTAAAAC` /
`ACTAGAGCACCAGAAGTCTA`, `get_patterns()`, score threshold 0.5):

| motif | fraction_depleted | num_sequences | avg_log2fc |
|---|---|---|---|
| GAGACC | 0.938776 | 49 | -2.768577 |
| CACNNNNGTAT | 1.000000 | 97 | -3.220028 |
| CACNNNNGTAC | 1.000000 | 78 | -3.258632 |
| AACNNNNCTTT | 0.994233 | 867 | -3.372092 |

Fixed-position motif: CTTT @4 (242 seqs).

## Findings that shape the plan

- **The fixed-position step is needed.** Fixed motifs are flexible motifs straddling the library flank (JJ1886 GTG @0
  is ATACNNNNGTG with ATAC in the left flank). Switching the step off (fixed score threshold > 1) loses
  ATACNNNNGTG and GTACNNNNGTG on JJ1886 and lowers the 12049 counts (97/78/867 → 60/56/555): many flexible motifs
  match exactly the same flank-dominated sequences and all fail the unique-hit rescoring. Removing those sequences
  before discovery, then rescoring on all sequences, recovers them. Script: see "Useful scripts" below.
- **`filter_to_core_motifs` is inverted** relative to its docstring (it keeps a longer motif only when it scores
  *worse*). No effect on either example and it cannot change which sequences are removed, but it can add useless
  rows.
- **Thresholds mix `<` and `≤`**: `score` uses `log2FC < thr`; candidate selection uses `>` for score and support;
  the final filter uses `>=` and `avg_log2fc <= thr`; notebook text says "log2FC ≤ -1".
- **`find_restricted_motifs` and `find_restricted_motifs_mp` differ**: the first filters redundancy after each
  pattern, the second once at the end; with an order-dependent filter they can disagree.
- **Order dependence**: when a motif and its reverse complement tie, which one is kept depends on row order
  (GAGACC vs GGTCTC).

## Rules for every phase

- Compare both datasets against the reference tables above after each phase.
- `nbdev-test` under the memory cap, then `nbdev-clean` + `nbdev-export`; the library must be in sync.
- One commit per phase. Phases 1 and 3 must not change results; in phase 2 every difference is reported and
  explained before committing.
- "Results" means the motif tables and fixed-position counts above, nothing else. Printed cell outputs (directory
  listings, `df.head()`, `df.columns`) do change in phase 0 and that is expected — `nbdev-test --save` rewrites
  them.

## Phase 0: Regression harness and example data (do first)

The comparison scripts from the first session are gone, and every later phase rests on "results did not change".
Build the check before touching anything, so each phase ends with one command instead of a manual read.

1. A committed check (`tests/test_reference_results.py`, or a `#| hide` cell) that runs both datasets at the
   reference settings above and asserts the motif rows and fixed-position counts, to a fixed tolerance. Run it under
   the memory cap at the end of every phase.
2. Synthetic fixture `example_data/toy_library.csv.gz`: a few hundred sequences from a fixed seed with a planted
   motif at a known depletion rate. It makes the unit tests fast and gives phase 3 the hand-checkable examples it
   asks for in sections 2, 4, 6 and 7 — where the expected answer is known by construction, not by having been
   observed once.
3. Slim the example data, 17 MB → ~1.7 MB:
   - `countsTable.csv` → `counts_jj1886.csv.gz` holding `seq,JJ1886_T0,MFDpir` (5.39 MB → 0.12 MB). Verified: no
     source cell in any notebook reads another column; the other 127 appear only inside printed `df.columns` /
     `df.head()` output. Name the sequence column `seq` to match the 12049 file (it is currently unnamed), and
     update the loading cells.
   - `counts_long.csv` and `counts_long_clean.csv`: delete (9.6 MB). Their only consumers are the `#| eval: false`
     cells that phase 1.1 removes anyway, including the one that writes `counts_long_clean.csv` back into the
     package data folder.
   - `counts_12049.csv.gz`: unchanged, already minimal.
   - Git history keeps the full table, so nothing is lost from the repo — only from the installed package.
4. Re-run the harness. The reference tables must be byte-identical after the data change; if they are not, the
   column trim was wrong.

## Phase 1: Remove dead code (no change to results)

1. Core notebook: remove `get_fold_change_values_per_site_old` + its example, the `filted_log2fc_df` cell, and the
   three `#| eval: false` leftovers (the `counts_long` cleanup that writes into the package data folder, the 16226 run,
   the duplicate JJ1886 `_mp` run).
2. Imports: drop `groupby`, `chain`, `itemgetter`, `choice`, the duplicate `re`, unused `typing` names; move
   `warnings`, `Pool`, `functools` to the top import cell.
3. Docstrings: module description (currently "Fill in a module description here"); remove "all sequences same
   length" from `get_sites_in_seq`; remove the non-existent `fc_improvement_margin` from `filter_to_core_motifs`.
4. Utils: remove functions the library never uses: `check_specific_matches_broad_iupac`,
   `create_motif_presence_matrix`, `flatten`, `get_all_sites`, `allseqs`. Checked: none is called from `core.py`
   or `plotting.py`. `flatten`, `get_all_sites` and `allseqs` are called from `old_nbs/` — confirm `old_nbs/` is
   frozen and not exported before deleting. (The `flatten` hit in `00_core.ipynb` is the word in a comment; that
   cell uses `chain.from_iterable`.)
5. Merge `get_motif_filter_with_context` and `_get_filter_for_motif` into one public function using the shared IUPAC
   table; update `plotting.py` and `update_motif_scores_from_unique_hits`.
6. Remove the unused `find_restricted_motifs` import in `plotting.py`.
7. (Moved to phase 0.3 — the data files go with the `#| eval: false` cells that use them.)

## Phase 2: Logic fixes (may change results; check each separately)

1. One pipeline function `find_restricted_motifs(..., n_jobs=1)`: scan every pattern, then run the redundancy filter
   once over all candidates. Keep `find_restricted_motifs_mp` as a thin deprecated alias for Bea's scripts.
2. Make the redundancy filter order-independent: process candidates in a fixed order (pattern, then motif); on a
   motif / reverse-complement tie keep the alphabetically first.
3. Fix the comparison in `filter_to_core_motifs` (keep a longer motif only if it beats the shorter one by more than
   the margin); add a test on a small made-up table.
4. Strict `<` everywhere (decision 1): `score` already does this; change candidate selection (`core.py:275`,
   `core.py:590`) and the final filter (`core.py:739-741`), and state the convention in the docstrings. No reference
   value sits on a boundary, so this should not move the tables — confirm with the harness rather than assume it.
5. `max_candidates` raises an error (decision 2) carrying the threshold advice from `core.py:878`, instead of
   returning unfiltered candidates that look like normal output.
6. Experiment, reported but not adopted (decision 4): repeat the unique-hit rescoring until no more motifs are
   dropped; record whether it changes either dataset, and write the answer into this file.

## Phase 3: Rewrite the core notebook as a walkthrough (no change to results)

Each section: short explanation, a small hand-checkable example, and a real-data example.

0. What RandSeq measures (random library, delivery into the strain, restriction sites depleted) + pipeline overview.
1. Counts → log2FC; library flanks and why they matter.
2. `score`: fraction of sequences depleted (toy example).
3. Step 1, fixed-position motifs: `identify_depleted_motifs_scanning_ends` + `filter_to_core_motifs` on JJ1886,
   showing that GTG @0 is ATACNNNNGTG straddling the flank.
4. Step 2, set those sequences aside, with the with/without comparison from "Findings".
5. Step 3, flexible patterns: `(d1, spacer, d2)`, `get_patterns`, `encode_library`, `get_pattern_scores`
   (toy + JJ1886).
6. Step 4, redundancy filtering: made-up table with CACNNNNGTA / CACNNNNGTAT / CACNNNNGTAC showing both directions
   of the rule.
7. Step 5, unique-hit rescoring: a motif losing its support; why it is done over all sequences.
8. Full pipeline: `find_restricted_motifs` on JJ1886 and 12049.
9. Appendix: string-based reference scanner, the fast scan's encoding, equivalence tests.

Utils and plotting keep their structure, but every exported function gets a docstring and at least one example or
test.

## Phase 4: Wrap up

- Push to `eren` and confirm CI passes.
- Short summary of result changes for David and Bea (especially anything affecting Fig3).

## Decisions (2026-09-29)

1. **Depleted is `log2FC < thr`** (strict), as `score` does today (`core.py:171`). Candidate selection
   (`core.py:275`, `core.py:590`) and the final filter (`core.py:739-741`) change to match, and the notebook text
   "log2FC ≤ -1" is corrected.
2. **`max_candidates` raises an error** carrying the threshold advice already drafted at `core.py:878`, instead of
   silently returning unfiltered candidates.
3. **Example data is slimmed** — see phase 0.3. The full 129-column counts table is archived outside the repo, so
   the package ships only what its examples use.
4. **Iterative rescoring (phase 2.6): run the experiment and report the result**, do not adopt it as the default.

## Still to tell Bea (no email sent yet — send before phase 2, not after)

The rescoring revert already landed in `b758f8f`, so the Fig3 numbers have moved whether or not the cleanup
finishes. If she is working with them now, this should not wait for phase 4.

- The unique-hit rescoring change in core_v2 was reverted (it lost ATACNNNNGTG on JJ1886); her Fig3 numbers may
  shift (AACNNNNCTTT on 12049: 640 → 867).
- The redundancy filter is order-dependent; strain 381A (> 150 candidates) should be re-checked.
- The notebook examples used `cpu_count=8`, which crashed a 16 GB machine; no longer needed.

## Useful scripts

The comparison and logic-check scripts from the first session were in the session scratchpad and are gone. Phase 0.1
replaces them with a committed check. The one-off check (fixed-position step on/off) is easy to redo: call
`find_restricted_motifs` twice with the reference settings above, the second time with `fixed_motif_score_thr=1.01`
so no fixed motif is kept.

## Verified against the code on `eren` (2026-09-29)

- `filter_to_core_motifs` inversion confirmed at `core.py:322-326`: `candidate_score >= core_score - margin` sets
  `is_subsumed = True`, so the longer motif is dropped when it scores *as well or better* and kept only when it
  scores worse by more than the margin. The docstring at `core.py:289` also documents an `fc_improvement_margin`
  parameter that does not exist.
- Threshold mixing confirmed: `<` at `core.py:171`, `>` at `core.py:275` and `core.py:590`, `>=` and `<=` at
  `core.py:739-741`.
- Example data audited; see phase 0.3 for the measured sizes.
