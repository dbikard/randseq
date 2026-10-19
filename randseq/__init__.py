"""Infer the DNA motifs bacterial defence systems restrict, from RandSeq random-insert plasmid libraries

Modules:

- `randseq.coding`: Representing library sequences and motifs as integers, so that scoring a motif is arithmetic rather than string matching.
- `randseq.fast`: Scoring one pattern: its totals from an optional compiled kernel or from the numpy path, and the scores reported from them.
- `randseq.utils`: The helpers `core` is built from: fold changes, motif strings and strands, the search space, custom-motif lookups and significance; and helpers for analysis scripts."""

__version__ = "0.1.0"

from .utils import calculate_log2fc, get_patterns
from .core import find_restricted_motifs, TooManyCandidatesError
from .plotting import plot_motif_analysis, plot_motif_depletion
from .example_data import get_example_data_dir

__all__ = ["calculate_log2fc", "get_patterns", "find_restricted_motifs", "TooManyCandidatesError",
           "plot_motif_analysis", "plot_motif_depletion", "get_example_data_dir"]
