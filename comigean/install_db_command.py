"""
Usage:
    comigean install-db [--help] <DIR>

Required positional argument:
    DIR       Path to folder to write output

'comigean install-db' downloads the NCBI Taxonomy Database required
to download reference and outgroup genomes based on taxa ids. Note: The
directory you point to should either not exist prior to running the script or
it should be empty.
"""
from docopt import docopt


if __name__ == '__main__':
    print(docopt(__doc__))
