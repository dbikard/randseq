import os
def get_example_data_dir():
    """
    Return the path to the example_data directory.
    """
    return os.path.join(os.path.dirname(__file__))

def get_example_backbone():
    """The pMin backbone that every plasmid in the example panel library carries.

    2,214 bp, written from the right flank round to the left flank -- i.e. reading away from the
    insert and back to it. It is not circular on its own: the insert closes it.
    """
    path = os.path.join(os.path.dirname(__file__), "pMin_backbone.fasta")
    with open(path) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">"))
