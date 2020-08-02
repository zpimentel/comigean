"""
Usage:
    comigean genome-stats [--help] [options] <REF_DIR> <OUTDIR>

Required positional arguments:
    OUTDIR                 Path to a directory to download genomes & proteomes
    REFDIR                 Path to database directory.


At least one the following 3 are required:
    USER_GENOMES           Path to directory containing user provided genome(s) that end with extension ".fna"

Options:
    --count_genomes        Will count the number of genomes to be downloaded.

'comigean genome-stats' profiles genomic features of species in the resultant phylogeny.
"""

from docopt import docopt


if __name__ == '__main__':
    print(docopt(__doc__))
