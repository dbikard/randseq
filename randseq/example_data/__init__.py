import os
def get_example_data_dir()->str: # the directory holding the example count tables
    "The directory of the example data shipped with the package."
    return os.path.join(os.path.dirname(__file__))
