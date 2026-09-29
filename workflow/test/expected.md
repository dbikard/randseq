# Expected self-test result

The synthetic library has 6000 members of 150 nt, of which 412 carry
`GGTCTC` (or its reverse complement `GAGACC`) on either strand. In `TestStrain` those members are
depleted to 4% of their molecules; `Donor` is undepleted. 12% of reads
carry one substitution in the insert.

`results/motifs.tsv` must contain exactly one call for `TestStrain`:

| field | expected |
|---|---|
| `motif` | `GGTCTC` — reported on the REBASE strand, not `GAGACC` |
| `median_log2fc` | strongly negative (about -4) |
| `escaper_rate` | low, a few percent |
| `n_plasmids` | order 412, exact value depends on `count_threshold` |

`results/qc/report.txt` must say the run passed.

If the motif comes back as `GAGACC`, canonicalisation is not working. If nothing is called, check
the flank sequences in `workflow/config.yaml` first — that is the most common configuration
error, and `flanks_found` in the QC report will be low.
