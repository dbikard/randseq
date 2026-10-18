# The example panel

`randseq/example_data/counts_panel.csv.gz` holds eight strains and both donor replicates,
as raw UMI counts. Each strain is there because it exercises one behaviour of the pipeline, and
this page records what each one does so the walkthrough's claims can be checked.

Measured 2026-10-17 (merged motifs 2026-10-18) at the pipeline's defaults: log2FC of replicate R1 against `Control_R1`
(`calculate_log2fc`, `count_threshold=20`, 30,108 plasmids kept), then
`find_restricted_motifs(series, left, right)` — `get_patterns()`, median-normalised, depleted below −0.5,
score 0.5. The motifs called are the rows it returns; the motifs scored but not called come back
with `all_motifs=True`.

| strain | motifs called | candidates entering step 4 |
|---|---|---|
| JJ1886 | `GGTCTC`, `RTACNNNNGTG`, `AAAGNNNNGTT` | 4 |
| CFT073 | `GAGNNNNNNNGTCA` | 1 |
| HS | `GGTAAG` | 3 |
| LMR_503 | `GGTCTC`, `RTACNNNNGTG`, `AAAGNNNNGTT`, `ACGNNNNNGTTG` | 8 |
| B156 | `GACNNNNGTTY`, `CCAYNNNNNTGT` | 4 |
| 381A | `CCTNAGG` | 91 |
| EDL933 | `CACNNNNNNNCTGG` | 1 |
| EDL933_dRM | **none** | 0 |

## What each one is for

**JJ1886** — the worked example. Three active RM systems, strong clean signal, and the motifs
are independently known. The Type I site `RTACNNNNGTG` is called as two variants and reported as one motif.

**CFT073** — weak but real. A single motif that single-plasmid conjugation assays missed.

**HS** — thresholds are not universal. BREX interference is weaker than restriction-modification:
`GGTAAG` is called at the default −0.5 and missed at −1. Anyone applying a −1 threshold to a BREX
strain and concluding "no activity" would be wrong.

**LMR_503** — co-occurring motifs of very different strength in one strain. This is the case
steps 3 and 4 exist for: without them the strong motifs' signal bleeds into many shadow
descriptions.

**B156** — a site formed across the vector junction, and two degenerate sites. The left flank
ends in `AAAC`, so every insert starting `NNNNGTC` completes `AAACNNNNGTC`: 382 of its 894
carriers carry it that way. The scan reads the flanks with the insert, so the site is found like
any other. Both sites are called as two variants each and reported as `GACNNNNGTTY` (=
`RAACNNNNGTC`) and `CCAYNNNNNTGT`.

**381A** — a noisy sample: the largest step-4 input of the panel, and a third of its library
depleted for reasons unrelated to any site. Only the real `CCTNAGG` is called, against a
GC-matched background of 0.30.

**EDL933 and EDL933_dRM** — specificity, demonstrated isogenically. The wild type yields
`CACNNNNNNNCTGG`; its RM knockout yields nothing. A negative result from an unrelated strain
proves much less than a matched pair.

`00_core`'s *The example panel at a glance* runs every strain and reproduces the table above.
