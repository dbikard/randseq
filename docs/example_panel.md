# The example panel

`randseq/example_data/counts_panel.csv.gz` holds eight BB2 strains and both donor replicates,
as raw UMI counts. Each strain is there because it justifies a design decision, and this page
records what each one actually does so the walkthrough's claims can be checked.

Verified 2026-09-29 with the settings the pipeline defaults to: log2FC computed per replicate
against that replicate's own donor at `count_threshold=10`, averaged, then
`find_restricted_motifs` at `score_thr=0.5`. 33,608 plasmids survive filtering.

| strain | threshold | motifs called |
|---|---|---|
| JJ1886 | −1.0 | `GGTCTC`, `ATACNNNNGTG`, `CACNNNNGTAC`, `AAAGNNNNGTT` |
| CFT073 | −1.0 | `GAGNNNNNNNGTCA` |
| HS | −1.0 | **none** |
| HS | −0.5 | `GGTAAG` |
| LMR_503 | −1.0 | `GGTCTC`, `ATACNNNNGTG`, `CACNNNNGTAC`, `AAAGNNNNGTT`, `ACGNNNNNGTTG` |
| B156 | −1.0 | `GAACNNNNGTC`, `AAACNNNNGTC`, `ACANNNNNATGG`, `ACANNNNNGTGG` |
| 381A | −1.0 | `TooManyCandidatesError` (455 candidates) |
| EDL933 | −1.0 | `CACNNNNNNNCTGG` |
| EDL933_dRM | −1.0 | **none** |

## What each one is for

**JJ1886** — the worked example. Three active RM systems, strong clean signal, and the motifs
are independently known, so it is the right place to introduce the pipeline.

**CFT073** — weak but real. A single motif that single-plasmid conjugation assays missed. It
shows why a sensitive method needs a defensible threshold rather than an obvious one.

**HS** — thresholds are not universal. BREX interference is weaker than restriction-modification,
so at the default −1.0 nothing is called; at −0.5 the cognate `GGTAAG` appears. Anyone applying
the defaults to a BREX strain and concluding "no activity" would be wrong.

**LMR_503** — co-occurring motifs of very different strength in one strain. This is the case the
redundancy filter and the unique-hit rescoring exist for: without them the strong motifs' signal
bleeds into the weak ones.

**B156** — why the fixed-position step exists. Its fixed-position hit is `GTC` at offset 4, which
is not a separate site: the left flank ends in `AAAC`, so `GTC@4` *is* `AAACNNNNGTC` formed
across the junction between vector and insert. The motif only looks position-specific because
the flank pins it there.

**381A** — what pathology looks like. It produces far more candidates than the unique-hit step
will score, and raises rather than returning them. In the published Fig 3 table this strain
silently contributed 344 of 374 rows as unfiltered candidates.

**EDL933 and EDL933_dRM** — specificity, demonstrated isogenically. The wild type yields
`CACNNNNNNNCTGG`; its RM knockout yields nothing. A negative result from an unrelated strain
proves much less than a matched pair.

## Regenerating

`work/build_example_data.py` builds the panel from
`randseq_paper/data/RandSeq_BB2_UMIcounts_EM_n10_R{1,2}.csv`;
`work/verify_panel_claims.py` reproduces the table above.
