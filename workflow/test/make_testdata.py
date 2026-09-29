"""Build the pipeline's self-test dataset.

A synthetic library with one planted restriction motif, so the expected answer is known by
construction rather than by having been observed once. `snakemake --config
samplesheet=workflow/test/samplesheet.csv` must recover exactly that motif.

Regenerate with:  python workflow/test/make_testdata.py
"""
import gzip
import random
from pathlib import Path

HERE = Path(__file__).parent
LEFT = "GTCTAGGGCGGCGGTAAAAC"
RIGHT = "ACTAGAGCACCAGAAGTCT"
BASES = "ACGT"

# The planted site. GGTCTC is BsaI / Eco31I -- a real Type IIS site, and the one the published
# screen recovers in JJ1886, so a user who knows the biology can sanity-check the answer.
MOTIF = "GGTCTC"
RC = "GAGACC"

N_MEMBERS = 6000
INSERT_LEN = 150
MOLECULES = 25             # mean; per-member abundance is drawn around this, see emit()
ABUNDANCE_SPREAD = 0.6     # lognormal sigma. A real library is far from uniform, and replicate
                           # correlation is only meaningful when abundances actually vary.
ERROR_RATE = 0.12          # per read, one substitution in the insert
DEPLETION = 0.04           # motif-carrying plasmids retain this fraction of molecules
SEED = 20260929


def rand_seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def has_motif(s):
    return MOTIF in s or RC in s


def build_library(rng):
    """Members with no site, plus members carrying exactly one, at a known rate."""
    members, with_site = [], 0
    while len(members) < N_MEMBERS:
        s = rand_seq(rng, INSERT_LEN)
        if has_motif(s):
            with_site += 1                      # keep naturally occurring sites
        members.append(s)
    # Plant extra sites so the motif has enough support to be called in a small test library.
    n_plant = 400 - with_site
    for i in rng.sample(range(len(members)), max(n_plant, 0)):
        if has_motif(members[i]):
            continue
        pos = rng.randrange(10, INSERT_LEN - len(MOTIF) - 10)
        members[i] = members[i][:pos] + MOTIF + members[i][pos + len(MOTIF):]
    return members


def member_abundance(rng_seed):
    """Per-member abundance, fixed across samples: it is a property of the library, not the run."""
    r = random.Random(rng_seed)
    return {i: max(1, round(r.lognormvariate(0, ABUNDANCE_SPREAD) * MOLECULES))
            for i in range(N_MEMBERS)}


def emit(path, members, rng, restrict, abundance):
    """Write one fastq. `restrict` applies depletion to motif-carrying members."""
    with gzip.open(path, "wt") as fh:
        n = 0
        for i, m in enumerate(members):
            mols = abundance[i]
            if restrict and has_motif(m):
                mols = sum(1 for _ in range(mols) if rng.random() < DEPLETION)
            for _ in range(mols):
                insert = m
                if rng.random() < ERROR_RATE:
                    p = rng.randrange(INSERT_LEN)
                    insert = m[:p] + rng.choice([b for b in BASES if b != m[p]]) + m[p + 1:]
                seq = rand_seq(rng, 10) + LEFT + insert + RIGHT + "GGGG"
                fh.write(f"@t{n}\n{seq}\n+\n{'I' * len(seq)}\n")
                n += 1
    return n


def main():
    rng = random.Random(SEED)
    members = build_library(rng)
    n_with = sum(1 for m in members if has_motif(m))
    print(f"library: {len(members)} members, {n_with} carrying {MOTIF}/{RC}")

    abundance = member_abundance(SEED)
    print(f"  abundance per member: min {min(abundance.values())}, "
          f"median {sorted(abundance.values())[len(abundance)//2]}, "
          f"max {max(abundance.values())}")

    rows = ["sample,replicate,fastq,role"]
    for sample, role, restrict in (("Donor", "donor", False), ("TestStrain", "recipient", True)):
        for rep in ("R1", "R2"):
            p = HERE / f"{sample}_{rep}.fastq.gz"
            n = emit(p, members, random.Random(SEED + hash((sample, rep)) % 1000), restrict,
                     abundance)
            print(f"  {p.name}: {n} reads")
            rows.append(f"{sample},{rep},workflow/test/{p.name},{role}")
    (HERE / "samplesheet.csv").write_text("\n".join(rows) + "\n")

    (HERE / "expected.md").write_text(
        f"""# Expected self-test result

The synthetic library has {len(members)} members of {INSERT_LEN} nt, of which {n_with} carry
`{MOTIF}` (or its reverse complement `{RC}`) on either strand. In `TestStrain` those members are
depleted to {DEPLETION:.0%} of their molecules; `Donor` is undepleted. {ERROR_RATE:.0%} of reads
carry one substitution in the insert.

`results/motifs.tsv` must contain exactly one call for `TestStrain`:

| field | expected |
|---|---|
| `motif` | `{MOTIF}` — reported on the REBASE strand, not `{RC}` |
| `median_log2fc` | strongly negative (about -4) |
| `escaper_rate` | low, a few percent |
| `n_plasmids` | order {n_with}, exact value depends on `count_threshold` |

`results/qc/report.txt` must say the run passed.

If the motif comes back as `{RC}`, canonicalisation is not working. If nothing is called, check
the flank sequences in `workflow/config.yaml` first — that is the most common configuration
error, and `flanks_found` in the QC report will be low.
""")
    print(f"wrote {HERE/'samplesheet.csv'} and expected.md")


if __name__ == "__main__":
    main()
