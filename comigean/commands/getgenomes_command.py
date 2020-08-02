"""
Usage:
    comigean get-genomes [--help] [options] [--ref_list] [--outgroup_list] [--count] [--code=<11>]
                         <OUTDIR> <REF_DIR> [--USER_GENOMES=<PATH>] [--REF_TAXA=<TAXA>]
                         [--OUTGROUP_TAXA=<TAXA>] [--refseq_category=<all>]
                         [--assembly_level=<Complete Genome>]

Required positional arguments:
    OUTDIR                 Path to a directory to download genomes & proteomes
    REFDIR                 Path to database directory.

At least one the following 3 are required:
    USER_GENOMES           Path to directory containing user provided genome(s) that end with extension ".fna"
    REF_TAXA               An NCBI Taxonomy ID for a taxa that includes your genome. Path to list of genomes
                           can be provided if '--ref_list' flag is invoked.
    OUTGROUP_TAXA          An NCBI Taxonomy ID for a taxa that will serve as your outgroup. Path to list
                           of genomes can be provided if '--outgroup_list' flag is invoked.

Optional positional argument:
    --ref_list             If flag is invoked, a path to a file containing one
                           NCBI assembly accession per line is expected in the
                           <REF_TAXA> field.
    --outgroup_list        If flag is invoked, a path to a file containing one
                           NCBI assembly accession per line is expected in the
                           <OUTGROUP_TAXA> field.
    refseq_category        Download genomes of following type (all, representative,
                           reference). Default: all.
    assembly_level         Download genomes of following type (all, Chromosome, Complete_Genome,
                           Scaffold, Contig). Default: Complete_Genome.

Options:
    --count_genomes        Will count the number of genomes to be downloaded.
    --code=<11>            Specify a genetic code for gene calling. Default: 11.

'comigean get-genomes' takes an NCBI taxonomy ID for the taxa that will be
used as a reference in your phylogenomic analysis. Additionally, it is highly
recommended (but optional) to provide an NCBI taxonomy ID for an outgroup as
well to facilitate tree rooting. This command extends ncbi-genome-download
for genome downloading and uses prodigal for gene calling.
"""

from docopt import docopt


if __name__ == '__main__':
    print(docopt(__doc__))
