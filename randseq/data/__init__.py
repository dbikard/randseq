import os


def get_data_dir():
    """
    Return the path to the reference data directory shipped with the package.
    """
    return os.path.join(os.path.dirname(__file__))
