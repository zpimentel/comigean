"""
Usage:
    comigean find-markers [--help] [options] <PROTEOME_DIR> <REF_DIR>

Required positional arguments:
    PROTEOME_DIR                 Path to a directory to download genomes & proteomes
    REFDIR                       Path to database directory.

Options:
    --code=<11>            Specify a genetic code for gene calling. Default: 11.

"""
from docopt import docopt


if __name__ == '__main__':
    print(docopt(__doc__))
