"""inferring restricted motifs using the  RandSeq method

Modules:

- `randseq.coding`: Representing library sequences and motifs as integers, so that scoring a motif is arithmetic rather than string matching.
- `randseq.fast`: Scoring one pattern: its totals from an optional compiled kernel or from the numpy path, and the scores reported from them.
- `randseq.utils`: The helpers `core` is built from: fold changes, motif strings and strands, the search space, custom-motif lookups and significance; and helpers for analysis scripts."""

__version__ = "0.0.2"
