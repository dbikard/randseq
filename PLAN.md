# Plan: clean up randseq and rewrite the core notebook as a walkthrough

Written 2026-09-28 on a laptop, before the data was available. **Revised 2026-09-29 on the
Maestro cluster**, after reproducing Fig 3 against the real BB2 panel. The machine rules, the
reference tables and the phase order all changed; see "What the data changed" below.

**Current goal (set by DB, 2026-09-29):** get the package to a state DB is happy with, then hand
it to Bea so she redoes the manuscript figures with it. Correctness first, ship fast. Scope for
the handoff is phases 0.1 / 2.2 / 2.5 plus the strand fix — *not* the full cleanup.

## Status at a glance

| phase | item | status |
|---|---|---|
| 0.1 | regression harness | **done** — `tests/test_reference_results.py`, 6 tests |
| 0.2 | example data for the walkthrough | **done** — superseded by a curated 8-strain real panel, see below |
| 0.3 | slim example data | **done** — 17 MB → 2.8 MB, and 1 illustrative strain → 9 |
| 1 | remove dead code | **done** — 658 lines, results unchanged |
| 2.1 | one pipeline function `n_jobs=` | in progress — measuring whether the two semantics differ first |
| 2.2 | order-independence | **done** — the redundancy filter now sorts deterministically and tie-breaks on the motif, not row position |
| 2.3 | `filter_to_core_motifs` inversion | **done** — was doubly inverted; no effect on the reference tables |
| 2.4 | threshold convention | **done** — and it DID move `num_sequences`, contrary to the plan's prediction |
| 2.5 | `max_candidates` raises | **done** |
| 2.6 | iterative-rescoring experiment | not done |
| 3 | notebook walkthrough | not done |
| 4 | push + tell Bea | **partly** — `CHANGELOG.md` written; commit `7190691` not yet pushed; Bea not contacted |
| — | canonical motif strand | **done** (not in the original plan) |
| — | `core_v2` compat shim | **done** (not in the original plan) |

Unpushed commits on `eren`: 9 and counting.

### Plan items that turned out to be wrong

Recorded because the plan marked them "verified against the code", and they were not:

1. **Example data** — the plan said no cell reads a column beyond `JJ1886_T0`/`MFDpir`. The
   loading cells filter on *all 19* `_T0` columns and `01_utils` prints `min(axis=0)` across
   them, so the specified 2-column trim would have silently changed a printed output.
2. **`get_motif_filter_with_context` is live external API.** Item 1.5 said merge it away; Bea's
   Fig 5 scripts import it by name. Unified behind her name instead, after checking the two
   implementations agree on 9 real motifs.
3. **The typing imports are used.** Only `choice`, `groupby`, `chain`, `itemgetter` and a
   duplicated `re` were dead.
4. **`find_restricted_motifs` in `plotting.py`** is unused by the module but used by a notebook
   example; it moved to a non-exported cell rather than being deleted.
5. **Phase 2.4 does move results.** The plan predicted no reference value sat on a boundary.
   One did: `min_support` was compared with `>`, so `min_support=3` meant four.
6. **`get_sites_in_seq`'s docstring** claimed all sequences must be equal length and that it
   raises otherwise. Neither is true, and it matters — only ~91% of real inserts are exactly
   150 nt.

## Working rules (Maestro, not WSL)

The previous version of this file described a 15 GB WSL box and `systemd-run` memory caps.
None of that applies. On Maestro:

- **This is a SLURM cluster — never run the pipeline on the login node.** `srun` for short
  checks, `sbatch` for anything longer. The whole 21-strain panel is ~4 min on 4 cores / 32 GB;
  a single example dataset is 2–15 s. Memory has not been a constraint here, so the memory-cap
  workaround is gone.
- **`/local/scratch` is node-local.** A script a compute node must read has to live on shared
  storage — use `../work/`, not the session scratchpad.
- **Call env executables by absolute path.** `~/.local/bin` is first on PATH and wins *even
  inside `conda run -n randseq`*. `~/.local/bin/nbdev-test` is a symlink into another user's
  venv, so `conda run -n randseq nbdev-test` silently runs the wrong interpreter and fails with
  a confusing `ModuleNotFoundError`. Use `~/miniconda3/envs/randseq/bin/nbdev-test`,
  `~/miniconda3/envs/randseq/bin/python -m pytest`, etc. `conda run -n randseq python` *is*
  safe, because `.local/bin` has no `python`.
- Env is `randseq` under **miniconda** (`~/miniconda3/envs/randseq`, python 3.11), not anaconda.
  Install with `python -m pip install -e .` from the repo root — note `.[dev]` does not exist and
  silently installs nothing.
- nbdev 3.3.24, hyphenated commands, config in `pyproject.toml`. `nbdev-test --n-workers 1`.
  `nbs/llms.txt` is generated and gitignored.

## Reference results (current, post-canonicalisation)

**These replace the tables in the previous version of this file.** Canonicalisation renamed
motifs; no statistic moved. Asserted by `tests/test_reference_results.py`.

Both example datasets are **E. coli JJ1886** — Lab.ID 12049 *is* JJ1886 — from different
sequencing runs. They are a cross-run replicate of one strain, and they must now report the
*same* motifs. Before canonicalisation they reported every motif on opposite strands.

JJ1886 (`countsTable.csv`, column `JJ1886_T0`, reference `MFDpir`, flanks
`GTCCTAGGTATAATACTAGT` / `GTTTTAGAGCTAGAAATAGC`, patterns
`[(6,0,0),(7,0,0),(4,2,4),(4,4,3),(4,3,4),(4,4,4)]`, score threshold 0.7):

| motif | was called | fraction_depleted | num_sequences | avg_log2fc |
|---|---|---|---|---|
| GGTCTC | GAGACC | 1.0 | 6 | -3.750094 |
| CACNNNNGTAC | GTACNNNNGTG | 1.0 | 6 | -4.922717 |
| ATACNNNNGTG | — | 1.0 | 146 | -4.744952 |
| AAAGNNNNGTT | — | 1.0 | 54 | -5.018200 |

Fixed-position motifs: GTG @0 (187 seqs), AAAG @12 (42 seqs). Fixed-position motifs are **not**
canonicalised — an offset is only meaningful on one strand.

12049 (`counts_12049.csv.gz`, column `12049_R1`, reference `Control`, flanks
`GTCTAGGGCGGCGGTAAAAC` / `ACTAGAGCACCAGAAGTCTA`, `get_patterns()`, score threshold 0.5):

| motif | was called | fraction_depleted | num_sequences | avg_log2fc |
|---|---|---|---|---|
| GGTCTC | GAGACC | 0.938776 | 49 | -2.768577 |
| ATACNNNNGTG | CACNNNNGTAT | 1.000000 | 97 | -3.220028 |
| CACNNNNGTAC | — | 1.000000 | 78 | -3.258632 |
| AAAGNNNNGTT | AACNNNNCTTT | 0.994233 | 867 | -3.372092 |

Fixed-position motif: CTTT @4 (242 seqs).

Library sizes: JJ1886 11,629 sequences (11,503 after the count filter); 12049 34,171 (30,108).

## What the data changed

Reproducing Fig 3 against the real 21-strain panel settled three things the plan had guessed at.

1. **The non-determinism is already fixed.** The plan treated order-dependence as a phase-2
   tidiness item. Bea's `code/fig3/README.md` in the paper repo documents it as a live problem —
   under `PYTHONHASHSEED` 0/1/2 she saw `num_sequences` swing 61/77/80 and 85/85/33. **That does
   not reproduce on `eren`**: the full panel is byte-identical across three hash seeds. The
   vectorized numpy scan in `b758f8f` removed the set/dict iteration the old scan depended on.
   Phase 2.2 is therefore mostly *already done*; what it still needs is a test (now written) and,
   strictly, an order-independent filter by construction rather than by luck.
2. **The strand problem reached the manuscript.** Results §2 calls the JJ1886 site "the
   low-abundance **BsaI** motif (**GAGACC**)"; Results §3 calls the ST131 site
   "5′-**GGTCTC**-3′ ... the widespread ST131 **Eco31I**". BsaI and Eco31I are isoschizomers —
   `GGTCTC(1/5)` — so that is one site written two ways in one paper. Fixed by canonicalising to
   REBASE's spelling; **§2 needs correcting**.
3. **`max_candidates` was worse than assumed.** Bea's reference Fig 3 table is 374 rows, of which
   **344 are unfiltered candidates from strain 16223 (381A)** alone, which then flow into the
   IUPAC merge and produce motifs like `SSBNNSVNN`. Nothing published is wrong only because 381A
   is dropped by hand. Now raises; the panel output is a clean 30 rows over 20 strains.

## Findings that shape the plan

- **The fixed-position step is needed.** Fixed motifs are flexible motifs straddling the library
  flank (JJ1886 GTG @0 is ATACNNNNGTG with ATAC in the left flank). Switching the step off
  (fixed score threshold > 1) loses ATACNNNNGTG and CACNNNNGTAC on JJ1886 and lowers the 12049
  counts (97/78/867 → 60/56/555): many flexible motifs match exactly the same flank-dominated
  sequences and all fail the unique-hit rescoring. Removing those sequences before discovery,
  then rescoring on all sequences, recovers them.
- **`filter_to_core_motifs` is inverted** relative to its docstring (it keeps a longer motif only
  when it scores *worse*). No effect on either example and it cannot change which sequences are
  removed, but it can add useless rows.
- **Thresholds mix `<` and `≤`**: `score` uses `log2FC < thr`; candidate selection uses `>` for
  score and support; the final filter uses `>=` and `avg_log2fc <= thr`; notebook text says
  "log2FC ≤ -1". The same muddle exists upstream in `calculate_log2fc`, whose `>` makes the
  documented "≥ 10 UMIs" filter an effective ≥ 11 — see `../DATA.md`.
- **`find_restricted_motifs` and `find_restricted_motifs_mp` differ**: the first filters
  redundancy after each pattern, the second once at the end; with an order-dependent filter they
  can disagree.
- **Order dependence** (motif vs reverse complement): resolved at the output boundary by
  `canonical_motif`, not inside the filter. Good enough for correct results; still worth fixing
  properly in 2.1/2.2.
- **Support counts moved with the rescoring revert** and nobody has adjudicated them. Across the
  panel `avg_log2fc` shifts by at most 0.026 but `num_sequences` rises a lot — `AAACNNNNGTC` on
  16171 is 341 against Bea's 127. Which is *correct* is unknown. This is what makes 0.2 urgent.

## Rules for every phase

- Run `tests/test_reference_results.py` after every phase — it replaces the manual table
  comparison the old plan described.
- `nbdev-test --n-workers 1`, then `nbdev-clean` + `nbdev-export`; the library must be in sync.
  Use absolute env paths (see working rules).
- One commit per phase. Phases 1 and 3 must not change results; in phase 2 every difference is
  reported and explained before committing.
- "Results" means the motif tables and fixed-position counts above, nothing else.

## Phase 0: Regression harness and example data

1. ~~A committed check that runs both datasets at the reference settings and asserts the motif
   rows and fixed-position counts.~~ **Done**: `tests/test_reference_results.py`. It asserts both
   tables to 1e-6, asserts the two datasets agree on motif labels, and reruns JJ1886 in
   subprocesses at `PYTHONHASHSEED` 0/1/2 to pin determinism.
2. **Synthetic fixture `example_data/toy_library.csv.gz` — do this next.** A few hundred
   sequences from a fixed seed with a planted motif at a known depletion rate. Originally
   justified as "makes unit tests fast and gives phase 3 hand-checkable examples". It is now
   also the **only way to settle which `num_sequences` is right**: with a planted motif the
   correct support count is known by construction, so the rescoring revert can be judged on
   evidence instead of on which option broke less. Promote ahead of everything else.
3. Slim the example data, 17 MB → ~1.7 MB:
   - `countsTable.csv` → `counts_jj1886.csv.gz` holding `seq,JJ1886_T0,MFDpir` (5.4 MB →
     0.12 MB). Verified: no source cell reads another column. Name the sequence column `seq` to
     match the 12049 file (currently unnamed), and update the loading cells — including
     `tests/test_reference_results.py`, which hardcodes `Unnamed: 0`.
   - `counts_long.csv` and `counts_long_clean.csv`: delete (9.7 MB). Only consumers are the
     `#| eval: false` cells that phase 1.1 removes.
   - `counts_12049.csv.gz`: unchanged, already minimal.
   - Consider *adding* a second strain from the BB2 panel: both current examples are JJ1886, so
     the suite tests one strain twice.
4. Re-run the harness. The tables must be identical after the data change.

## Phase 1: Remove dead code (no change to results)

1. Core notebook: remove `get_fold_change_values_per_site_old` + its example (still in
   `nbs/00_core.ipynb`, not exported), the `filted_log2fc_df` cell, and the three
   `#| eval: false` leftovers.
2. Imports: drop `groupby`, `chain`, `itemgetter`, `choice`, the duplicate `re`, unused `typing`
   names; move `Pool` and `functools` to the top import cell. **`warnings` is already moved** —
   but a second `import warnings` remains at `core.py:797`; remove it.
3. Docstrings: module description (currently "Fill in a module description here"); remove "all
   sequences same length" from `get_sites_in_seq`; remove the non-existent
   `fc_improvement_margin` from `filter_to_core_motifs`.
4. Utils: remove functions the library never uses: `check_specific_matches_broad_iupac`,
   `create_motif_presence_matrix`, `flatten`, `get_all_sites`, `allseqs`. `flatten`,
   `get_all_sites` and `allseqs` are called from `old_nbs/` — confirm it is frozen first.
5. Merge `get_motif_filter_with_context` and `_get_filter_for_motif` into one public function
   using the shared IUPAC table; update `plotting.py` and `update_motif_scores_from_unique_hits`.
6. Remove the unused `find_restricted_motifs` import in `plotting.py`.

## Phase 2: Logic fixes (may change results; check each separately)

1. One pipeline function `find_restricted_motifs(..., n_jobs=1)`: scan every pattern, then run
   the redundancy filter once over all candidates. Keep `find_restricted_motifs_mp` as a thin
   deprecated alias. **Note `randseq/core_v2.py` now also exists as a compat shim for Bea's
   scripts** — fold that into the same deprecation story rather than having two.
2. Make the redundancy filter order-independent *by construction*: fixed candidate order
   (pattern, then motif). The reverse-complement tie is already handled downstream by
   `canonical_motif`, and output is empirically stable across hash seeds, so this is now about
   guaranteeing the property rather than obtaining it.
3. Fix the comparison in `filter_to_core_motifs`; add a test on a small made-up table.
4. Strict `<` everywhere (decision 1): change candidate selection and the final filter, state the
   convention in the docstrings. No reference value sits on a boundary, so the tables should not
   move — confirm with the harness.
5. ~~`max_candidates` raises an error.~~ **Done**: `TooManyCandidatesError`, with the candidate
   table attached as `.candidates`.
6. Experiment, reported but not adopted (decision 4): repeat the unique-hit rescoring until no
   more motifs are dropped; record whether it changes either dataset, and write the answer here.
   Best done after 0.2, so the fixture can say which answer is right.

## Phase 3: Rewrite the core notebook as a walkthrough (no change to results)

Unchanged from the previous version, with two additions:

- Section 0 should say what the two example datasets *are* — both JJ1886, different runs — since
  the docs currently imply they are unrelated.
- Add a short section on reporting: why a motif and its reverse complement are one answer, and
  how `canonical_motif` chooses. This is the thing most likely to confuse someone comparing new
  output against an old table.

Sections otherwise as before: 0 what RandSeq measures; 1 counts → log2FC and flanks; 2 `score`;
3 fixed-position motifs; 4 setting those sequences aside; 5 flexible patterns and the fast scan;
6 redundancy filtering; 7 unique-hit rescoring; 8 full pipeline; 9 appendix with the string-based
reference scanner and equivalence tests.

## Phase 4: Wrap up

- Push `eren` (one unpushed commit, `7190691`) and confirm CI passes.
- ~~Short summary of result changes for David and Bea.~~ Written as `CHANGELOG.md`. **Not yet
  sent to Bea** — see below.

## Decisions

1. **Depleted is `log2FC < thr`** (strict), as `score` does today. Candidate selection and the
   final filter change to match; the notebook text "log2FC ≤ -1" is corrected. *(2026-09-29)*
2. **`max_candidates` raises an error.** *(2026-09-29, done)*
3. **Example data is slimmed.** The full 129-column counts table is archived outside the repo.
   *(2026-09-29)*
4. **Iterative rescoring: run the experiment and report**, do not adopt as default. *(2026-09-29)*
5. **Motif strand: match REBASE's spelling, fall back to alphabetically-first** for motifs REBASE
   does not know. Not plain alphabetical, which would give `GAGACC` and contradict the paper.
   The Gold Standard list (n=983, downloaded 2025-12-03) ships as package data. *(2026-09-29, done)*
6. **Handoff scope: correctness only.** Threshold cleanup and absorbing the fastq→counts step are
   deferred until after Bea has a usable version. *(2026-09-29)*

## Still to tell Bea — `CHANGELOG.md` is written, nothing sent

All of this is in `CHANGELOG.md` in a form she can act on. It has not been sent.

- Motif labels change (15 of 30 rows on the panel); no statistic moves. Her
  `methylome_crossref_Fig3.R` should be unaffected since it expands reverse complements, but the
  IUPAC-merged names in `..._iupac_merged_mm2.csv` will change spelling.
- `generate_fig3_motifs.py` will now raise on 16223; it needs a `try/except
  TooManyCandidatesError`. `../work/fig3_motifs.py` does this already.
- The rescoring revert raises `num_sequences` substantially (`AAACNNNNGTC` on 16171: 341 vs 127).
  **Results §2 still has "n=X, n=X and n=X" placeholders** — fill from a current run.
- Manuscript §2 `GAGACC` should become `GGTCTC` to match §3.
- Her vendored `randseq_package/` still has `core.py` + `core_v2.py` and has been forking from
  `eren` since 2026-09-25. The shim keeps her imports working, but the two copies should be
  reconciled.

## Useful scripts

- `../work/reference_run.py` — both example datasets at the reference settings.
- `../work/fig3_motifs.py` — the 21-strain panel; handles `TooManyCandidatesError`.
- `../work/fig3_seeds.sbatch` — the same, as a 3-seed array job for the determinism check.
- `../work/patch_nb.py`, `../work/add_cells.py` — edit nbdev notebooks (string replace / insert
  cells) so the library can be regenerated with `nbdev-export`.
- One-off check of the fixed-position step: call `find_restricted_motifs` twice, the second time
  with `fixed_motif_score_thr=1.01` so no fixed motif is kept.

## Verified against the code on `eren`

- `filter_to_core_motifs` inversion confirmed: `candidate_score >= core_score - margin` sets
  `is_subsumed = True`, so the longer motif is dropped when it scores *as well or better*. Its
  docstring also documents an `fc_improvement_margin` parameter that does not exist.
- Threshold mixing confirmed (`<` in `score`, `>` in candidate selection, `>=`/`<=` in the final
  filter). Line numbers from the pre-`7190691` tree and now shifted; re-grep rather than trust them.
- Motif statistics are **exactly** strand-symmetric — measured: `GAGACC`/`GGTCTC`,
  `AACNNNNCTTT`/`AAAGNNNNGTT` and `CACNNNNGTAT`/`ATACNNNNGTG` each give identical `n`, `score`
  and mean on 12049. This is what makes canonicalisation lossless.
- Example data audited; see phase 0.3 for measured sizes.
