"""
Usage:
    comigean get-genomes [--help] [options] <OUTDIR> <REF_DIR>
                         [--USER_GENOMES=<PATH>] [--REF_TAXA=<TAXA>] [--OUTGROUP_TAXA=<TAXA>]

Required positional arguments:
    OUTDIR                   Path to a directory to download genomes & proteomes
    REFDIR                   Path to database directory.

At least one the following 3 are required:
    USER_GENOMES             Path to directory containing user provided genome(s) that end with extension ".fna"
    REF_TAXA                 An NCBI Taxonomy ID for a taxa that includes your genome.
    OUTGROUP_TAXA            An NCBI Taxonomy ID for a taxa that will serve as your outgroup.

Optional positional argument:
    --refseq_category=option Download genomes of following type (all, representative,
                             reference) [default: all]
    --assembly_level=option  Download genomes of following type (all, Chromosome, Complete_Genome,
                             Scaffold, Contig) [default: Complete_Genome]

Options:
    --count                  Will count the number of genomes to be downloaded
    --code=NUM               Specify a genetic code for gene calling [default: 11]

'comigean get-genomes' takes an NCBI taxonomy ID for the taxa that will be
used as a reference in your phylogenomic analysis. Additionally, it is highly
recommended (but optional) to provide an NCBI taxonomy ID for an outgroup as
well to facilitate tree rooting.
"""

from docopt import docopt


if __name__ == '__main__':
    print(docopt(__doc__))
