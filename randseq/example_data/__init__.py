import os
def get_example_data_dir()->str: # the directory holding the example count tables and backbone
    "The directory of the example data shipped with the package."
    return os.path.join(os.path.dirname(__file__))

def get_example_backbone()->str: # 2,214 bp, from the right flank round to the left flank
    "The pMin backbone every plasmid of the example panel carries; linear, as the insert closes it."
    path = os.path.join(os.path.dirname(__file__), "pMin_backbone.fasta")
    with open(path) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">"))
