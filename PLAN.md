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
| 2.1 | one pipeline function `n_jobs=` | **done** — `7cdeb56`; `_mp` is a deprecated alias |
| 2.2 | order-independence | **done** — the redundancy filter now sorts deterministically and tie-breaks on the motif, not row position |
| 2.3 | `filter_to_core_motifs` inversion | **done** — was doubly inverted; no effect on the reference tables |
| 2.4 | threshold convention | **done** — and it DID move `num_sequences`, contrary to the plan's prediction |
| 2.5 | `max_candidates` raises | **done** |
| 2.6 | iterative-rescoring experiment | **done** — measured; one pass is the intended semantics (see Phase 2) |
| 3 | notebook walkthrough | **done** — `00_core` is the walkthrough; both Phase 3 additions are in |
| 4 | push + tell Bea | **partly** — `eren` is level with `origin/eren`; the 2026-10-13 work is uncommitted; Bea not contacted |
| — | canonical motif strand | **done** (not in the original plan) |
| — | `core_v2` compat shim | **done** (not in the original plan) |

Nothing on `eren` is unpushed (checked 2026-10-13); everything since `6ff2429` is uncommitted.

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
- **Prefix every absolute-path call with `PYTHONNOUSERSITE=1`.** The env sets it as a conda env
  var, which only `conda activate` / `conda run` apply. Calling `~/miniconda3/envs/randseq/bin/python`
  directly bypasses it, and numpy and scipy then load from `~/.local/lib/python3.11` (found
  2026-10-13; same versions there, so no result was affected — but it is not the env).
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
6. ~~Experiment, reported but not adopted (decision 4): repeat the unique-hit rescoring.~~
   **Done, 2026-09-29. Answer: it converges immediately but inflates support enormously, and
   must not be adopted.**

   Rescoring recomputes each motif's statistics using only plasmids that carry exactly one
   candidate, then drops motifs that no longer pass. Dropping a motif frees its plasmids to
   count as "unique" for the survivors, so the step can be iterated. Measured on four strains:

   - the **motif set is identical at every iteration** -- iterating never finds or loses a motif
   - the **statistics converge at iteration 2** and do not move again
   - but they converge to a far more permissive answer. On 12049, `CACNNNNGTAC` goes from
     n=17 to n=441 and `CACNNNNGTAT` from 31 to 483; on JJ1886 `GTACNNNNGTG` goes from 3 to 19

   That is the whole point of the step being undone. Each round re-admits plasmids that were
   excluded precisely because they carry more than one candidate, so iterating to a fixed point
   converges on something close to not having excluded them at all. One pass is the intended
   semantics. Script: `../work/iterative_rescoring.py`.

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

## Open questions raised 2026-10-02

1. **Should a motif and its reverse complement share one numerical code?** — **DONE 2026-10-13.**
   Measured, decided and implemented. Folded *within* a pattern's own code space, for the
   patterns where the two spellings live there: `spacer == 0` or `d1 == d2` (19 of the 43 in
   `get_patterns()`, but the code tests the shape per call, so it holds for any pattern a caller
   passes). The other 24 are untouched — each physical site already produced exactly one code,
   because a both-strand scan of `(4, 4, 3)` covers the whole shape class and `(3, 4, 4)` is
   dropped by `filter_symmetric_tuples`. Verified on the real library: `(3, 4, 4)` re-spelled as
   `(4, 4, 3)` is *exactly* `(4, 4, 3)`, 1,131,980 pairs, zero differences either way. Folding an
   asymmetric pattern would merge unrelated motifs, so `rc_code_table` raises for one.

   Results: identical on 8 configurations, every column. Scored motif space 750,592 -> 462,624
   rows (-38%); scan runtime unchanged. `canonical_motif` still runs at reporting, because the
   fold keeps the smaller code (`GAGACC`) and the literature name is `GGTCTC` — so this did
   *not* remove the strand rule, only the duplicate rows.

   Left deliberately in place: `_canonicalise_flexible_motifs` still collapses a duplicated pair
   defensively. The scan can no longer produce one, but a caller who passes both `(4, 4, 3)` and
   `(3, 4, 4)` explicitly still can.

2. **Step 3 is not monotonic in the pattern set.** Adding `(3, 4, 3)` to a search that already
   contains `(4, 4, 3)` deletes `ATACNNNNGTG` from JJ1886 — a call with 192 plasmids at a
   depletion fraction of 1.00 — because the shorter `CACNNNNGTA` (0.89) suppresses it under the
   "general is not significantly worse" rule at `score_margin=0.2`. Enlarging the search space
   should not remove a correct answer. Measured 2026-10-02; see the note below for what has been
   ruled out.

## Notebook tasks raised 2026-10-13 (DB, reading `00_core`)

1. ~~**"The real library" arrives ~22 cells too early, and its lead sentence is false.**~~
   **DONE 2026-10-13** — moved to just before `process_single_flexible_pattern`, its first use. Cell 14
   says *"The rest of this section works on the example library rather than three sequences by
   hand"* — but what follows is the twelve-insert `demo_lib` walkthrough (cells 17–33). `lib` and
   `fcs` are built at cell 15 and then not touched until cell 37
   (`process_single_flexible_pattern`). Move the section down to just before its first real use,
   and rewrite the lead sentence so it hands off to the demo rather than claiming to leave it.

2. ~~**The `motif_occurrences` demo exercises neither rejection rule.**~~ **DONE 2026-10-13** —
   three-insert demo (clean, `N`-bearing, too short) plus a test cell. `demo_lib` is twelve clean
   12-mers — no non-ACGT base and no short sequence — so neither the N-rejection nor the
   short-sequence/padding rejection that the docstring promises is ever shown where the reader
   meets the function. Add an `N`-bearing insert and a shorter one to the demo, show which
   windows drop out, and assert both behaviours in a test cell.

3. ~~**Say why the two strands are not collapsed at occurrence time.**~~ **DONE 2026-10-13** —
   overtaken by the fold itself (Open question 1 above). They now *are* collapsed there, and
   Step A says why, which pattern shapes it applies to, and that the surviving spelling is
   decided by the code rather than by biology. The "shape never searched" point was already made
   in the standalone-test section ("a reported pattern can be one `flexible_motif_patterns` never
   contains"). Stale fold prose fixed in the same pass: step 3's strand bullet, cells on the
   step-3 demo ("ten survive" was five), and Step B/C commentary.

## RESUME HERE — state at end of 2026-10-14 session

**Mid-change: step 3 rewritten as two order-free rules (DB approved).** Code done and exported:
`filter_redundant_patterns` = rule 1 (drop mixtures, `is_mixture`) + rule 2 (drop a one-base
extension of a rule-1 survivor, new `is_extension`); `_variants` and `report=` deleted. In the
ablation (`../work/step3_ablation.py`, mode "rule") it gave the same called sets as the old
pairwise+standalone code: 30/30, 0 extra on 43 patterns; same 11 extra on the wide test case.

**Not yet done:**
1. ~~Verification~~ **DONE**: identical row for row on the 43-pattern panel (137/137) and the
   wide test case (46/46). Slightly slower (panel 30 s -> 34 s) because rule 1 now tests every
   candidate, not only pairwise survivors -- DB asked to keep it fast, so look at this.
   Rule-2 examples on the notebook's own library (`../work/two_rules_explore.out`): with patterns
   (6,0,0),(7,0,0),(3,4,3),(3,4,4),(4,3,4),(4,4,4), rule 2 drops `CACANNNGTAC`, `CACCNNNGTAC`,
   `CACGNNNGTAC`, `CACTNNNGTAC` (the first spacer base of `CACNNNNGTAC` pinned to each base -- the
   ideal example) and `AGGTCTC`, `CGAGACC`, `GGTCTCA` (one base around `GAGACC`). With the current
   demo's three patterns rule 2 drops nothing, so the demo needs the wider pattern list.
2. **The `00_core` step-3 section is broken until rewritten**: the demo cell
   `profile = variant_profile(*_variants("CACNNNNGTA", ...))` still calls the deleted `_variants`,
   so `nbdev-test` fails on `00_core`. Rewrite the section as: intro (one site, several
   descriptions differing by one position) -> the variant profile (one tool) -> reading a profile
   (`verdict`) -> the positions next to a motif (`_undefined_positions`) -> rule 1 (mixture;
   `ATANNNNTGTG` example) -> rule 2 (extension; example from the exploration, or 7084's
   `CCCCNNNNNTTTC` from the panel) -> the filter + demo table. Ablation numbers for the prose: without
   rule 2's job (pairwise) 91 extra calls; rule 1 alone removes 16164's `GGCGNNNNNNNCGCC`.
3. Then: pytest, `nbdev-test --save` (outputs rule), CHANGELOG entry, the opposite-strand
   sentence DB queried is now in `is_extension`'s docstring.

**Also pending from today:** cfg0 still fails (`CCACNNNAGTA`, step-4 options B/C parked above);
nothing committed since `6ff2429`.

## Wide-search test case — keep for steps 3 and 4 (DB, 2026-10-13)

The default stays the paper's 43 patterns. The wider set
`get_patterns(contiguous=(4, 8), halves=(2, 5), max_spacer=9)` (68 patterns, 902/983 REBASE
motifs) is kept as a **test case**: run it, and every call it adds is an example of a redundant
description steps 3 and 4 should remove. Reproduce with `python ../work/fig3_motifs.py --wide`;
the reference output is `../work/fig3_wide68_testcase.csv`.

Against published Fig 3 it calls 29/30, loses one and adds 11:

- **ECOR64 (16169), the YGGCCR family.** Loses `GCGGCCA`; adds five 8-mers that all contain
  YGGCCR -- `GCGGCCGC`, `ACGGCCGC`, `AGTGGCCG`, `CTTGGCCA`, `GCTGGCCA`. Each is a one-base
  extension that picks out the most depleted plasmids of a weak, degenerate site (YGGCCR, mean
  log2FC -0.82, below the -1 threshold). The 8-base shape just gives that more ways to happen.
- **15853, the ACGNNNNNGTTG family.** `ACGNNNNNGTTGC`, `ACAACNNNNNCGT`, `ACGNNNNAGTTG` -- five-base
  halves and a shifted spacer around one site.
- **8098.** `CGCNGCGAG`, `CTCGCCGC`, `CTCGCGGA`.
- **381A (16223)** trips `max_candidates` again (170 > 150).

Rules any fix to steps 3-4 should satisfy on this case: remove all 11 without losing the 30, and
ideally report the degenerate site (YGGCCR) rather than any of its extensions.

## Search-space size — notes 2026-10-13 (DB: note for later, not decided)

The default was briefly widened to 68 patterns and **reverted to the paper's 43** after the
wide-search test case above (11 extra calls). Nine defined bases stays available
(`max_defined=9`: 77 patterns, 911/983) but is not the default: it scans 4.0M motifs instead of
1.6M, and in our libraries a 9-base site has few carriers (~37 in the 34k-plasmid 150 nt library,
~1 in the 11.6k-plasmid 20 nt one).

**Does widening change BH when the library cannot support the motifs?** Only motifs with
`num_sequences >= min_support` count as tests, so it should depend on support -- and it does, but
more of them clear the filter than intuition suggests. Measured, tests contributed by k=9
patterns: **2.36M on the 150 nt panel (65% of all tests -- every q roughly x3)**, and **231k on the
20 nt library (22%)**. On the short library the cause is the flanks: the scan reads
flank + insert + flank, and a window overlapping the constant vector sequence has those bases
fixed in every plasmid, so a 9-base motif with 6 bases in the flank has only 64 spellings, each on
~1/64 of the library. Junction-straddling motifs pass the support filter at any size.

**Done 2026-10-13 (no result change):** step 3 reads variant scores from step 2's per-pattern
totals (one kernel, each pattern scanned at most once), step 2 decodes only its candidates, and
the library is encoded once per run (and can be passed in, so a panel encodes it once). 21-strain
panel 72 s -> 31 s, every row and q-value identical. **Parked:** a multi-strain kernel that reads
each window once for all strains of a panel (the scan is ~0.5 s per strain now). **For the
step-4 review:** the `max_candidates` docstring and the `TooManyCandidatesError` message still
describe a "unique-hit scoring step" that no longer exists.

**Idea, parked: derive the default search space from the library.** Include a shape only if
`03_design`'s carrier formula predicts enough carriers for this library's size and insert length,
so the search follows the library rather than a fixed list (and four-base sites are searched in
short-insert libraries, where they are detectable). The flank effect above means "predicted
carriers" must be computed for the random part only, or junction motifs would justify every size.

## Settled 2026-10-13 — step 4 is cumulative, ties broken on support

Step 4 now tests each candidate against the union of every stronger survivor, and equal
depletion fractions are ordered by support rather than mean log2FC. 21-strain panel: 19 strains
identical, 15853 loses the shadow `ACGTNNNNGTTG` (now 30/30 against published Fig 3, 0 extra),
381A no longer trips `max_candidates` and calls `CCTNAGG` plus six noise-edge motifs. The support
tie-break alone changed nothing on the panel or the 8 example configs.

**Still open — `CCACNNNAGTA` on the older JJ1886 run (cfg0, the one failing test).** All 55 of its
plasmids carry one of the four real sites, yet it survives both changes. The last 2 plasmids it
keeps after removing `ATACNNNNGTG` and `AAAGNNNNGTT` carry `CACNNNNGTAC`, which has *fewer*
plasmids (20) and so ranks below it. No pass with an order fixed in advance can catch that: the order is set by each motif's full
size (55 vs 20), before anything is removed. Measured on cfg0's five motifs
(`../work/step4_orderings.py`):

- **(B) re-rank every round on what is left** (DB's idea). Each round, take the candidate with the
  most depleted plasmids not carried by anything already kept. `CCACNNNAGTA` has 2 left after
  round 1, 1 after round 2, 0 after round 3 -> dropped; the four real sites are kept. This is the
  natural statement of step 4: each motif must explain destroyed plasmids nothing chosen before
  it explains.
- **(C) audit after the pass.** Test every survivor, weakest first, against the union of all the
  *other* survivors. Same answer here. It also closes a gap (B) leaves: (B) checks a motif only
  against those taken before it, so a later pick can make an earlier one redundant.

Caveats to measure before choosing: ranking on depleted plasmids *left* favours big motifs, so a
broad half-depleted motif could be taken before a clean site and explain it away (step 3's mixture
test should remove such motifs first -- unverified); ranking on fraction left instead lets a
one-plasmid residue jump the queue. Preferred if pursued: (B) ranked on depleted plasmids left,
then (C), benchmarked on the 21-strain panel, the 8 example configs and cfg0, watching LMR_503,
381A and ECOR64. Parked 2026-10-13 while DB reviews steps 2-3.

## Settled 2026-10-02 (measured, not argued)

**Open question 2 above is resolved.** `score_margin` is gone. Step 3 now decides general-vs-
specific by scoring all four bases at the position two candidates disagree on, and calls a motif a
mixture only when that set is *bimodal* — some variant below `score_thr` **and** a spread of at
least `score_thr`. Monotonicity is now a standing benchmark check: adding `(3,4,3)` or `(2,4,4)`
to JJ1886's search no longer deletes its Type I site. The bimodality half is load-bearing: on
CFT073 three of four variant sets contain a variant below 0.5 (0.47–0.49) at spreads of 0.11–0.21,
so a threshold-only rule would delete that strain's only motif.

Each survivor is additionally tested on its own by N-padding, which catches a register-shifted
shadow (`ATANNNNTGTG` against the real `ATACNNNNGTG`) that no pairwise rule can see, because
neither motif contains the other.

Benchmark at one parameter set (log2FC < −1, score ≥ 0.5): **30/30 published Fig 3 motifs, 0
extra, 0 false positives on the 9 silent strains, panel exact, HS BREX recovered at −0.5.**
`min_support` 3 and 5 give **identical** results everywhere — it is not a sensitive parameter.

### Things tried and rejected — do not redo without new evidence

- **Dropping steps 4 and 5.** They are no-ops on 6 of 7 panel strains, which makes them look
  removable. They are not: without them sample 15853 gains `AACGNNNNNGTTG`, `ACGGNNNNGTTG` and
  `ACGTNNNNGTTG`, all more specific versions of the real `ACGNNNNNGTTG`, and 381A gains three
  motifs. Step 4 tests **independence**, not significance.
- **Tightening the call to q < 0.01.** Changes nothing, in either configuration. The extras above
  are statistically solid — their plasmids really are destroyed — so no q threshold removes them.
  Significance and independence are orthogonal; the q-value cannot substitute for step 4.
- **A q-value prefilter before step 3.** Cuts 381A's candidates 171 → 33 and touches nothing on
  the clean strains, which is the right shape for a prefilter. But it does not make steps 4 and 5
  redundant, and with step 3 now 17× faster the speed motive is largely gone. Not worth the risk
  to a weak site.
- **Caching the encoded library.** `encode_library` costs 0.12 s against a 15–45 s run. ~1%.

### Done

- **Step 3's standalone test is batched by pattern.** It was calling `_pattern_window_codes` 552
  times over 20 distinct patterns on 381A. Now the jobs are collected, `(motif, pattern)` pairs
  deduplicated, and one `motif_presence_matrix` call issued. 381A's step 3: **40.8 s → 2.4 s**;
  the test suite: 285 s → 77 s.

### Still open

- The **pairwise** comparisons in step 3 are still one `motif_presence_matrix` call each (68 of
  them on JJ1886, over 25 patterns). Batchable the same way, but now a small share of a small
  number.
- ~~`annotate_significance` double-counts the survivors' p-values in BH.~~ **FIXED 2026-10-13.**
  It was *anti*-conservative, not conservative as written here: a duplicated p-value takes two
  ranks and the step-up gives it the better one, so top q-values were up to 2x too small. No call
  changed on the panel.
- `avg_log2fc` separates signal from noise on 381A more cleanly than score or q does: the real
  `CCTNAGG` family runs −1.55 to −2.36, the noise −0.81 to −0.89. The pipeline reports it but
  never filters on it. Worth considering; not tested.
