"""Reference-result regression test.

Runs the two example datasets at fixed settings and asserts the exact motif tables. These
numbers are the contract: any change to the pipeline that moves them is either a bug or a
deliberate decision that has to be recorded in CHANGELOG.md.

Both datasets are E. coli JJ1886 (Lab.ID 12049) from different sequencing runs — countsTable is
the earlier run, counts_panel is BB2 — so the two tables must report the *same motifs*. That is
the point of the canonical-strand rule, and before it they reported every motif on opposite
strands.

Run with:   python -m pytest tests/test_reference_results.py
Determinism:  PYTHONHASHSEED=<n> python -m pytest tests/test_reference_results.py
"""
import os
import subprocess
import sys

import pandas as pd
import pytest

from randseq.core import calculate_log2fc, find_restricted_motifs, get_patterns
from randseq.example_data import get_example_data_dir

TOL = 1e-6

JJ1886 = dict(
    file="countsTable.csv.gz",
    seq_col="seq",
    sample="JJ1886_T0",
    reference="MFDpir",
    left="GTCCTAGGTATAATACTAGT",
    right="GTTTTAGAGCTAGAAATAGC",
    patterns=[(6, 0, 0), (7, 0, 0), (4, 2, 4), (4, 4, 3), (4, 3, 4), (4, 4, 4)],
    score_thr=0.7,
    # motif -> (fraction_depleted, num_sequences, avg_log2fc)
    expected={
        "GGTCTC": (1.0, 6, -3.750094),
        "CACNNNNGTAC": (1.0, 3, -4.800833),
        "ATACNNNNGTG": (1.0, 146, -4.744952),
        "AAAGNNNNGTT": (1.0, 54, -5.018200),
    },
    expected_fixed={("GTG", 0): 187, ("AAAG", 12): 42},
)

S12049 = dict(
    file="counts_panel.csv.gz",
    seq_col="seq",
    sample="JJ1886_R1",
    reference="Control_R1",
    left="GTCTAGGGCGGCGGTAAAAC",
    right="ACTAGAGCACCAGAAGTCTA",
    patterns=None,  # get_patterns()
    score_thr=0.5,
    expected={
        "GGTCTC": (0.938776, 49, -2.768577),
        "ATACNNNNGTG": (1.0, 31, -3.260202),
        "CACNNNNGTAC": (1.0, 17, -3.606989),
        "AAAGNNNNGTT": (0.994213, 864, -3.369221),
    },
    expected_fixed={("CTTT", 4): 242},
)


def run_case(cfg):
    path = os.path.join(get_example_data_dir(), cfg["file"])
    df = pd.read_csv(path).set_index(cfg["seq_col"])
    log2fc = calculate_log2fc(
        df[[cfg["reference"], cfg["sample"]]], reference_column=cfg["reference"]
    )[cfg["sample"]]
    return find_restricted_motifs(
        log2fc,
        cfg["left"],
        cfg["right"],
        flexible_motif_patterns=cfg["patterns"] or get_patterns(),
        flexible_motif_score_thr=cfg["score_thr"],
    )


@pytest.mark.parametrize("cfg,name", [(JJ1886, "JJ1886"), (S12049, "12049")])
def test_reference_motifs(cfg, name):
    fixed, flex = run_case(cfg)

    assert set(flex["motif"]) == set(cfg["expected"]), (
        f"{name}: motif set changed.\n"
        f"  expected {sorted(cfg['expected'])}\n"
        f"  got      {sorted(flex['motif'])}"
    )

    for _, row in flex.iterrows():
        frac, n, fc = cfg["expected"][row["motif"]]
        assert row["num_sequences"] == n, f"{name}/{row['motif']}: num_sequences"
        assert abs(row["fraction_depleted"] - frac) < TOL, f"{name}/{row['motif']}: fraction_depleted"
        assert abs(row["avg_log2fc"] - fc) < TOL, f"{name}/{row['motif']}: avg_log2fc"

    got_fixed = {(r["motif"], r["position"]): r["num_sequences"] for _, r in fixed.iterrows()}
    assert got_fixed == cfg["expected_fixed"], f"{name}: fixed-position motifs changed"


def test_both_datasets_agree_on_motifs():
    """Both example datasets are JJ1886, so they must name the same sites the same way.

    Before canonicalisation these two tables reported every motif on opposite strands.
    """
    shared = set(JJ1886["expected"]) & set(S12049["expected"])
    assert shared == set(JJ1886["expected"]) == set(S12049["expected"]), (
        "The two JJ1886 datasets no longer agree on motif labels: "
        f"{set(JJ1886['expected']) ^ set(S12049['expected'])}"
    )


@pytest.mark.slow
@pytest.mark.parametrize("seed", ["0", "1", "2"])
def test_deterministic_across_hash_seeds(seed, tmp_path):
    """Output must not depend on PYTHONHASHSEED.

    Python randomises string hashing per process, so set and dict iteration order changes
    between runs. The search used to consume that order, which made `num_sequences` and the
    reported strand vary run to run. This pins it.
    """
    script = (
        "import pandas as pd, os\n"
        "from randseq.core import calculate_log2fc, find_restricted_motifs\n"
        "from randseq.example_data import get_example_data_dir\n"
        "df = pd.read_csv(os.path.join(get_example_data_dir(), 'countsTable.csv.gz')).set_index('seq')\n"
        "l = calculate_log2fc(df[['MFDpir','JJ1886_T0']], reference_column='MFDpir')['JJ1886_T0']\n"
        "_, flex = find_restricted_motifs(l, 'GTCCTAGGTATAATACTAGT', 'GTTTTAGAGCTAGAAATAGC',\n"
        "    flexible_motif_patterns=[(6,0,0),(7,0,0),(4,2,4),(4,4,3),(4,3,4),(4,4,4)],\n"
        "    flexible_motif_score_thr=0.7)\n"
        "print(flex[['motif','fraction_depleted','num_sequences','avg_log2fc']].to_csv(index=False))\n"
    )
    env = dict(os.environ, PYTHONHASHSEED=seed)
    out = subprocess.run([sys.executable, "-c", script], env=env, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    table = out.stdout[out.stdout.index("motif,"):].strip()

    expected = "\n".join(
        ["motif,fraction_depleted,num_sequences,avg_log2fc"]
        + [
            f"{m},{v[0]},{v[1]},{v[2]:.6f}"
            for m, v in [
                ("GGTCTC", (1.0, 6, -3.750094)),
                ("CACNNNNGTAC", (1.0, 3, -4.800833)),
                ("ATACNNNNGTG", (1.0, 146, -4.744952)),
                ("AAAGNNNNGTT", (1.0, 54, -5.018200)),
            ]
        ]
    )
    got = "\n".join(
        line if i == 0 else ",".join(line.split(",")[:3] + [f"{float(line.split(',')[3]):.6f}"])
        for i, line in enumerate(table.split("\n"))
    )
    assert got == expected, f"PYTHONHASHSEED={seed} changed the result:\n{got}\n!=\n{expected}"
