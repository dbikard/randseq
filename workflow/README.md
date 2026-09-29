# RandSeq pipeline

FASTQ → molecule counts → log2FC → restriction motifs, with QC at every step.

RandSeq measures which DNA sequences a bacterium restricts, by conjugating a library of plasmids
carrying random inserts into it and sequencing what survives. Plasmids carrying a recognition
site for an active restriction system are destroyed on entry, so they are depleted relative to
the input library. The pipeline finds the motifs responsible.

## Install

```bash
conda create -n randseq python=3.11
conda activate randseq
pip install randseq randcount snakemake
```

## Check your install before touching real data

```bash
snakemake -s workflow/Snakefile --cores 4 \
    --config samplesheet=workflow/test/samplesheet.csv outdir=results_test
```

Two minutes on a laptop. It runs a synthetic library with one planted `GGTCTC` (BsaI / Eco31I)
site and must recover exactly that motif — see `test/expected.md` for the expected numbers. If
this does not pass, nothing downstream is trustworthy.

## Run your own data

Write a samplesheet, one row per FASTQ:

```csv
sample,replicate,fastq,role
Control,R1,/data/Control_R1.fastq.gz,donor
Control,R2,/data/Control_R2.fastq.gz,donor
JJ1886,R1,/data/JJ1886_R1.fastq.gz,recipient
JJ1886,R2,/data/JJ1886_R2.fastq.gz,recipient
```

`role=donor` marks the input library — the plasmid prep before conjugation. Everything is
measured relative to it, so exactly one sample must carry that role.

Then set your construct's flanking sequences in `config.yaml` and run:

```bash
snakemake -s workflow/Snakefile --profile workflow/profiles/slurm   # cluster
snakemake -s workflow/Snakefile --cores 8                           # one machine
```

## What you get

| file | contents |
|---|---|
| `qc/report.txt` | **read this first** — did the run work? |
| `motifs.tsv` | the calls, one row per motif per strain |
| `motifs_skipped.tsv` | strains that produced no confident call, and why |
| `counts/<sample>_<rep>.tsv` | molecules per library member |
| `log2fc_mean.tsv` | log2 fold change per plasmid per strain |
| `members.txt` | the library, discovered from the donor |

### Reading `motifs.tsv`

`median_log2fc` is the headline: how strongly plasmids carrying the motif were depleted, as a
log2 fold change. `-3` means roughly an 8-fold reduction in successful transfer.

`escaper_rate` is how *penetrant* the restriction is — the fraction of motif-carrying plasmids
that got through anyway. A low median log2FC with a high escaper rate means a system that blocks
most but not all plasmids, which is biologically different from one that blocks everything.

`n_plasmids` is how much evidence supports the call — **but it is set by `count_threshold`**, so
it is not comparable between experiments run at different depths unless you quote the threshold
with it. `median_log2fc` and `escaper_rate` are threshold-free and are the right things to
compare across runs.

`backbone_sites` is a warning flag: see below.

## Things that will bite you

**A motif already in your vector cannot be detected.** Every plasmid in the library carries it,
so there is no motif-free baseline to deplete against. The published screen hit exactly this and
had to change backbone to recover a masked Type I motif. Set `backbone_fasta` in the config and
the pipeline will flag any called motif that occurs in your vector — and, more importantly, tell
you that a *negative* result for such a motif is meaningless rather than informative.

**Your library size sets what you can detect.** A 6 bp motif needs roughly
`n_members × 2 × (insert_len − 5) / 4096` plasmids to carry it. Long bipartite Type I motifs
(8+ specified bases) need a much larger library than short Type II sites. The QC report computes
this for you; if `motif_coverage` fails, a negative result means "underpowered", not "absent".

**A strain with no restriction system is the *slow* case, not the fast one.** With nothing real
to find, the search wades through many near-threshold candidates. On the published data one such
strain ran for over three hours. The pipeline gives motif calling a long time limit, records
strains that exceed the candidate limit in `motifs_skipped.tsv`, and carries on rather than
losing the whole run to one sample.

**Thresholds are not universal.** The defaults (`log2FC < -1`, score `>= 0.5`) were tuned on
Enterobacteriaceae restriction–modification systems. BREX interference is weaker and needs
`-0.5`. Report whatever you used.

**Replicates catch contamination.** The published screen excluded strains whose biological
replicates correlated at r ≤ 0.75. The QC report applies the same rule.

## Configuration

Everything is in `config.yaml`, which documents each option inline. The ones you are most likely
to change:

- `left_flank` / `right_flank` — your construct's constant regions. **Check these first** if
  `flanks_found` in the QC report is low; a wrong flank is the most common setup error.
- `count_threshold` — molecules required in the donor for a plasmid to be analysed. Tune it to
  your depth. It does not bias which motifs are called (depletion is independent of plasmid
  abundance) but it does set `n_plasmids`.
- `backbone_fasta` — your vector, for the masking check above.

## Notes on the method

Library members are discovered **once**, from the deepest donor replicate, and every sample is
counted against that one list. This is what makes samples directly comparable; counting each
sample independently does not guarantee a shared set of rows.

log2FC is computed within each replicate against its own donor and then averaged across
replicates, so a difference in sequencing depth between replicates does not leak into the result.

Insert-level error correction is available (`--correct` in `randcount`) but **off by default**.
Measured end to end on real data it recovers about 6% of reads and leaves the motif calls
identical, because the gain is proportional in donor and recipient and log2FC is a ratio within
that pair. Turn it on if you have reason to think error rates differ between the samples being
compared.
