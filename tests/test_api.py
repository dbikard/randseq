"""The functions a typical analysis needs are importable from the package itself."""
import randseq


def test_top_level_api():
    for name in ["calculate_log2fc", "get_patterns", "find_restricted_motifs", "TooManyCandidatesError",
                 "plot_motif_analysis", "plot_motif_depletion", "get_example_data_dir"]:
        assert hasattr(randseq, name), name
        assert name in randseq.__all__, name
    # the same objects as in their modules, not copies
    from randseq.core import find_restricted_motifs
    assert randseq.find_restricted_motifs is find_restricted_motifs
