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

from randseq.core import call_reason, calculate_log2fc, find_restricted_motifs, get_patterns
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
        # Support is the number of plasmids that actually carry each motif, checked against a
        # direct regex count on 2026-10-15; the earlier values (138/61/18) were undercounts.
        # Rule 3 of step 3 is what removed the shadow CCACNNNAGTA, which used to be called here;
        # holding step 4's leftovers to the scan's cutoff removed ACCTNNNCTTG (2026-10-17).
        # avg_log2fc is measured from the sample's median since centring became the default.
        # The Type I site was two rows, ATACNNNNGTG (1.0, 192, -4.805313) and CACNNNNGTAC
        # (1.0, 20, -5.287198); since 2026-10-18 its called variants are reported as one IUPAC
        # motif, carried by the 212 plasmids that carry either (CHANGELOG 26).
        "GGTCTC": (1.0, 9, -4.363412),
        "RTACNNNNGTG": (1.0, 212, -4.850774),
        "AAAGNNNNGTT": (1.0, 63, -5.209253),
    },
    merged={"RTACNNNNGTG": "ATACNNNNGTG+CACNNNNGTAC"},
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
        # Support is the number of plasmids that actually carry each motif (direct regex count,
        # 2026-10-13). Since 2026-10-15 the defaults centre each series on its median and use
        # log2FC < -0.5: the same four motifs on the same plasmids, with depletion fractions and
        # avg_log2fc now measured from the sample's median. The Type I site was two rows,
        # ATACNNNNGTG (1.0, 504, -3.265818) and CACNNNNGTAC (1.0, 468, -3.312154); since
        # 2026-10-18 it is one IUPAC motif on the 963 plasmids that carry either (CHANGELOG 26).
        "GGTCTC": (0.981818, 55, -2.851532),
        "RTACNNNNGTG": (1.0, 963, -3.285517),
        "AAAGNNNNGTT": (0.997809, 913, -3.425663),
    },
    merged={"RTACNNNNGTG": "ATACNNNNGTG+CACNNNNGTAC"},
)


COMPACT = ["motif", "num_sequences", "fraction_depleted", "avg_log2fc", "q_value"]


def run_case(cfg, **kwargs):
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
        **kwargs,
    )


@pytest.mark.parametrize("cfg,name", [(JJ1886, "JJ1886"), (PANEL, "panel")])
def test_reference_motifs(cfg, name):
    motifs = run_case(cfg)

    # The default return is the called motifs alone, in five columns, sorted by q.
    assert list(motifs.columns) == COMPACT
    assert motifs["q_value"].is_monotonic_increasing
    assert set(motifs["motif"]) == set(cfg["expected"]), (
        f"{name}: called motif set changed.\n"
        f"  expected {sorted(cfg['expected'])}\n"
        f"  got      {sorted(motifs['motif'])}"
    )
    for _, row in motifs.iterrows():
        frac, n, fc = cfg["expected"][row["motif"]]
        assert row["num_sequences"] == n, f"{name}/{row['motif']}: num_sequences"
        assert abs(row["fraction_depleted"] - frac) < TOL, f"{name}/{row['motif']}: fraction_depleted"
        assert abs(row["avg_log2fc"] - fc) < TOL, f"{name}/{row['motif']}: avg_log2fc"

    # all_motifs=True: the same calls, first, with every column and a reason for each row.
    full = run_case(cfg, all_motifs=True)
    called = full[full["called"]]
    assert list(full["called"]) == sorted(full["called"], reverse=True)
    pd.testing.assert_frame_equal(called[COMPACT].sort_values(["q_value", "motif"])
                                  .reset_index(drop=True), motifs)
    assert set(called["reason"]) == {"called"}
    assert not full.loc[~full["called"], "reason"].str.startswith("called").any()
    for _, row in called.iterrows():
        assert row["merged_from"] == cfg["merged"].get(row["motif"], ""), \
            f"{name}/{row['motif']}: merged_from"


def test_reason_of_uncalled_motif():
    """At the defaults, the 20 nt JJ1886 run scores GTCGNNNNNNNCGGA as significant (q = 1.2e-3)
    but with only 7 of 12 carriers depleted: the bound on the depleted fraction fails it."""
    path = os.path.join(get_example_data_dir(), JJ1886["file"])
    df = pd.read_csv(path).set_index("seq")
    log2fc = calculate_log2fc(df[["MFDpir", "JJ1886_T0"]], reference_column="MFDpir",
                              count_threshold=20)["JJ1886_T0"]
    full = find_restricted_motifs(log2fc, JJ1886["left"], JJ1886["right"],
                                  all_motifs=True).set_index("motif")
    assert full.loc["GTCGNNNNNNNCGGA", "reason"] == \
        "most carriers not shown depleted (7/12; lower bound 0.32 < 0.5)"


def test_call_reason():
    """Each verdict in words: called, not significant, bound failed, both, bimodal note."""
    row = pd.Series({"num_sequences": 10, "fraction_depleted": 0.9, "q_value": 1e-6,
                     "fraction_depleted_lower": 0.6, "bimodal": False})
    assert call_reason(row) == "called"
    late = row.copy(); late["q_value"] = 0.21
    assert call_reason(late) == "not significant (q = 0.21)"
    weak = row.copy(); weak["fraction_depleted"], weak["fraction_depleted_lower"] = 0.5, 0.22
    assert call_reason(weak) == "most carriers not shown depleted (5/10; lower bound 0.22 < 0.5)"
    weak["q_value"] = 0.3
    assert call_reason(weak) == ("not significant (q = 0.3); "
                                 "most carriers not shown depleted (5/10; lower bound 0.22 < 0.5)")
    row["bimodal"] = True
    assert call_reason(row) == "called; carriers split between depleted and untouched"


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
        "m = find_restricted_motifs(l, 'GTCCTAGGTATAATACTAGT', 'GTTTTAGAGCTAGAAATAGC',\n"
        "    flexible_motif_patterns=[(6,0,0),(7,0,0),(4,2,4),(4,4,3),(4,3,4),(4,4,4)],\n"
        "    flexible_motif_score_thr=0.7, all_motifs=True)\n"
        "rows = sorted((r['motif'], int(r['num_sequences']), round(float(r['avg_log2fc']), 6),\n"
        "               bool(r['called'])) for _, r in m.iterrows())\n"
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
