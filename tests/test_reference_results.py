"""Reference-result regression test.

Runs the two example datasets at fixed settings and asserts the exact motif tables. These
numbers are the contract: any change to the pipeline that moves them is either a bug or a
deliberate decision that has to be recorded in CHANGELOG.md.

Both datasets are E. coli JJ1886 from different sequencing runs — countsTable is
the earlier run, counts_panel a later one — so the two tables must report the *same motifs*. That is
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
        "GGTCTC": (1.0, 9, -4.313837),
        # This config restricts the pattern list and raises score_thr to 0.7. Under those
        # settings the variant test in step 3 retains a different candidate set from the
        # score-margin rule it replaced, so step 5 scores these two on fewer uniquely-carried
        # plasmids -- same motifs, same perfect depletion, smaller denominator. On the default
        # settings (get_patterns(), score_thr 0.5) both rules give identical support.
        "CACNNNNGTAC": (1.0, 18, -5.057596),
        "ATACNNNNGTG": (1.0, 138, -4.716200),
        "AAAGNNNNGTT": (1.0, 61, -5.037651),
    },
    expected_fixed={("GTG", "left", 0): 187, ("AAAG", "right", 4): 42},
)

PANEL = dict(
    file="counts_panel.csv.gz",
    seq_col="seq",
    sample="JJ1886_R1",
    reference="Control_R1",
    left="GTCTAGGGCGGCGGTAAAAC",
    right="ACTAGAGCACCAGAAGTCTA",
    patterns=None,  # get_patterns()
    score_thr=0.5,
    expected={
        # Support is now the number of plasmids that actually carry each motif -- checked
        # against a direct regex count on 2026-10-13. The previous values (52/483/441/882) were
        # undercounts left by the step-1 set-aside, which removed sequences before the scan.
        "GGTCTC": (0.945455, 55, -2.804666),
        "ATACNNNNGTG": (0.996032, 504, -3.218952),
        "CACNNNNGTAC": (1.0, 468, -3.265288),
        "AAAGNNNNGTT": (0.994524, 913, -3.378797),
    },
    expected_fixed={("CTTT", "left", 4): 242},
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


@pytest.mark.parametrize("cfg,name", [(JJ1886, "JJ1886"), (PANEL, "panel")])
def test_reference_motifs(cfg, name):
    fixed, flex = run_case(cfg)

    # Rows below the significance threshold are returned on purpose, flagged rather than
    # dropped, so the contract is about what is *called*.
    called = flex[flex["called"]] if "called" in flex.columns else flex
    assert set(called["motif"]) == set(cfg["expected"]), (
        f"{name}: called motif set changed.\n"
        f"  expected {sorted(cfg['expected'])}\n"
        f"  got      {sorted(called['motif'])}\n"
        f"  (also reported, not called: "
        f"{sorted(set(flex['motif']) - set(called['motif']))})"
    )

    for _, row in called.iterrows():
        frac, n, fc = cfg["expected"][row["motif"]]
        assert row["num_sequences"] == n, f"{name}/{row['motif']}: num_sequences"
        assert abs(row["fraction_depleted"] - frac) < TOL, f"{name}/{row['motif']}: fraction_depleted"
        assert abs(row["avg_log2fc"] - fc) < TOL, f"{name}/{row['motif']}: avg_log2fc"

    got_fixed = {(r["motif"], r["anchor"], r["position"]): r["num_sequences"]
                 for _, r in fixed.iterrows()}
    assert got_fixed == cfg["expected_fixed"], f"{name}: fixed-position motifs changed"


def test_both_datasets_agree_on_motifs():
    """Both example datasets are JJ1886, so they must name the same sites the same way.

    Before canonicalisation these two tables reported every motif on opposite strands.
    """
    shared = set(JJ1886["expected"]) & set(PANEL["expected"])
    assert shared == set(JJ1886["expected"]) == set(PANEL["expected"]), (
        "The two JJ1886 datasets no longer agree on motif labels: "
        f"{set(JJ1886['expected']) ^ set(PANEL['expected'])}"
    )


@pytest.mark.slow
def test_deterministic_across_hash_seeds():
    """Output must not depend on PYTHONHASHSEED.

    Python randomises string hashing per process, so set and dict iteration order changes
    between runs. The search used to consume that order, which made `num_sequences` and the
    reported strand vary run to run -- the paper repository's own notes record support counts
    of 61/77/80 for one motif across three seeds.

    This asserts the runs agree with each other, rather than re-pinning values that
    `test_reference_motifs` already owns.
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
        "rows = sorted((r['motif'], int(r['num_sequences']), round(float(r['avg_log2fc']), 6),\n"
        "               bool(r['called'])) for _, r in flex.iterrows())\n"
        "print('RESULT', rows)\n"
    )
    results = {}
    for seed in ("0", "1", "2"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        out = subprocess.run([sys.executable, "-c", script], env=env,
                             capture_output=True, text=True)
        assert out.returncode == 0, out.stderr
        results[seed] = out.stdout[out.stdout.index("RESULT"):].strip()

    assert results["0"] == results["1"] == results["2"], (
        "PYTHONHASHSEED changed the result:\n" +
        "\n".join(f"  seed {k}: {v}" for k, v in results.items())
    )
    assert "GGTCTC" in results["0"], "sanity: the known BsaI site should be in there"
