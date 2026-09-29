# Core notebook: outline for the walkthrough rewrite

DB, 2026-09-29: *"the core notebook reading like a walkthrough and explaining all design
decisions is critical to me as it will let me review what we have done and vouch for it."*

So the test for every section is: **could DB read this and defend the choice to a reviewer?**
That means each decision needs the alternative stated, the evidence that settled it, and the
cost of getting it wrong — not just a description of what the code does.

Every claim below has been measured this session. Numbers in *italics* are the evidence to cite.

---

## 0. What RandSeq measures

Random 150 bp inserts in a mobilisable plasmid, conjugated into a recipient. A plasmid carrying
a recognition site for an active restriction system is destroyed on entry, so it is depleted
relative to the input library. Motifs are found by asking which sequence patterns are shared by
the depleted plasmids.

Then: **what the example data is.** Both `countsTable.csv.gz` (earlier run) and the `JJ1886`
columns of `counts_panel.csv.gz` (BB2) are *the same strain*, Lab.ID 12049. They are a cross-run
replicate, and the docs previously implied they were unrelated examples. The panel's other seven
strains each exist to justify a specific decision — `example_data/counts_panel_key.csv`.

## 1. Counts → log2FC

`calculate_log2fc`. Pseudocount added **before** normalisation, column sums taken over the
retained rows. Both are worth stating because the Methods currently describe the opposite order,
and adding 1 to an already-normalised count would erase the signal.

**Decision to justify — the count threshold is a power choice, not a correctness one.** It
filters on *donor abundance*, and depletion is independent of abundance (*paper, pooled Pearson
r = 0.135, Fig S4b*). So raising it costs statistical power but cannot bias which motifs are
called. The consequence: `num_sequences` is set by the threshold and is therefore **not
comparable between experiments**, while `median_log2fc` and escaper rate are. This is the single
most important thing for a reader to take away, because it defuses the temptation to compare raw
n across runs.

## 2. `score`: the fraction of plasmids depleted

Toy example. State the threshold convention here, since this function defines "depleted":
log2FC thresholds are **strict** (`< thr`); score and support thresholds are **inclusive**
(`>= thr`). This matches the paper's recommended operating thresholds, *"log2FC < −1; depletion
score ≥ 0.5"*. Note the manuscript also says "log₂FC ≤ −1" elsewhere and needs correcting.

*Cost of getting it wrong: `min_support=3` used to mean four.*

## 3. Fixed-position motifs — and why this step exists at all

`identify_depleted_motifs_scanning_ends` + `filter_to_core_motifs` on **B156**.

The worked example that makes it click: B156's fixed-position hit is `GTC` at offset 4. That is
not a position-specific site. The left flank is `GTCTAGGGCGGCGGTAAAAC`, whose last four bases are
`AAAC` — so `GTC@4` **is** `AAACNNNNGTC`, formed across the junction between vector and insert.
The flank pins it to a constant offset, which is why it shows up in a positional scan.

Same story in JJ1886: `GTG@0` is `ATACNNNNGTG` with `ATAC` from the left flank.

`filter_to_core_motifs`: keep a longer motif only if the extra base raises the depletion fraction
by more than the margin. *The comparison was inverted until this session — a longer motif was
kept only when it scored **worse**.* Show the four hand-made cases.

## 4. Setting those sequences aside

`filter_sequences_without_core_motifs`. Why discovery runs on the reduced set but **rescoring
runs on the full series**.

*Evidence, measured on B156: of `AAACNNNNGTC`'s 883 supporting plasmids, 376 carry the fixed
motif. They deplete at 0.952 against 0.966 for the rest — statistically the same, so they are
genuine support and excluding them would discard 43% of the evidence for the motif.*

The counter-example, same strain: for `AAAGNNNNGTT` (JJ1886's motif, not B156's) the 9
fixed-motif plasmids deplete at 1.000 against 0.057 for the other 888. That is what the failure
mode looks like — it just never reaches threshold here.

## 5. Flexible patterns and the fast scan

`(d1, spacer, d2)`, `get_patterns`, `encode_library`, `get_pattern_scores`. Toy first, then
JJ1886. Appendix carries the equivalence proof against the string scanner.

## 6. Redundancy filtering

**LMR_503** — four co-occurring motifs spanning −1.31 to −3.26 in one strain, which is exactly
what this filter is for.

**Decision to justify:** the filter must be a function of the candidate *set*, not of the order
patterns were scanned in. Equivalent motifs used to be tie-broken by row position, so which of a
motif/reverse-complement pair survived depended on scan order — and the single-process and
multiprocessing pipelines scan in different orders. Now sorted deterministically and tie-broken
on the motif string. *Verify by filtering shuffled copies.*

## 7. Unique-hit rescoring

Why statistics are recomputed on plasmids carrying exactly one candidate. A motif losing its
support. Report the iterative-rescoring experiment here (phase 2.6) — measured, not adopted.

## 8. Reporting: what comes out and how to read it

- **A motif and its reverse complement are one answer.** Measured: identical `n`, score and mean
  for every pair. So exactly one should be reported, and `canonical_motif` chooses by REBASE
  spelling with an alphabetical fallback.
  *Cost of getting it wrong: the manuscript currently calls the same site `GAGACC` (BsaI) in
  Results §2 and `GGTCTC` (Eco31I) in §3. They are isoschizomers.*
- **`num_sequences` is an analysis artefact** (back-reference to §1). Report it with its
  threshold or not at all.
- What to compare across experiments: `median_log2fc` (strength), escaper rate (penetrance).

## 9. Failure modes — the section a user will actually need

- **381A**: 455 candidates, `TooManyCandidatesError`. *In the published Fig 3 table this strain
  silently contributed 344 of 374 rows as unfiltered candidates, which then flowed into the
  IUPAC merge as motifs like `SSBNNSVNN`.* Raising beats returning something that looks like a
  result.
- **HS**: nothing at −1.0, `GGTAAG` at −0.5. Thresholds tuned on restriction–modification do not
  transfer to BREX. Applying the defaults and concluding "no activity" would be wrong.
- **EDL933 vs EDL933_dRM**: `CACNNNNNNNCTGG` vs nothing. Isogenic specificity control.
- **Backbone masking**: a motif already in the vector cannot be detected, because every plasmid
  carries it and there is no motif-free baseline. *The published screen hit this and had to change
  backbone.* A negative result for such a motif is meaningless, not informative.
- **An RM-free strain is the slow case**, not the fast one — nothing real to find means many
  near-threshold candidates. *Observed: >3 h on MG1655 ΔhsdR.*

## 10. Appendix

String-based reference scanner; the fast scan's integer encoding; equivalence tests; the full
pipeline on JJ1886 and the panel.

---

## Rules while writing

- Results must not change. Run `tests/` and `nbdev-test` after.
- Every number quoted must be reproducible from the shipped example data, or explicitly marked
  as coming from the BB2 panel / the manuscript.
- Prefer showing the failure to describing it.
