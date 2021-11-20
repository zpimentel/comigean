"""
Author: Zachary Pimentel
Institution: University of Rhode Island
"""

# coding density
# ANI
# AAI
# full taxonomy

import os
import sys

from Bio import SeqIO
import pandas as pd


def genome_stats(genome_dir):
    ''' Computes genome size and GC content  '''
    dirs = ['reference', 'outgroup', 'user']
    genome_stats_dict = {}
    #gc_content_dict = {}
    #contig_count_dict = {}

    for dir in dirs:
        if os.path.exists(genome_dir + "/" + dir):
            for filename in os.listdir(genome_dir + "/" + dir):
                if filename.endswith(".fna"):
                    basename = os.path.splitext(filename)[0]
                    genome_size = 0
                    g_count = 0
                    c_count = 0
                    contig_count = 0
                    with open(genome_dir + "/" + dir + "/" + filename, "rU") as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            genome_size += len(record.seq)
                            g_count += record.seq.count('G')
                            c_count += record.seq.count('C')
                            contig_count += 1

                    genome_size_str = str(round(genome_size/1000000, 3)) + " MB"
                    gc_content = round(((g_count + c_count) / genome_size) * 100, 2)
                    gc_content_string = str(gc_content) + "%"

                    genome_stats_dict[basename] = (genome_size_str, contig_count, gc_content_string)

    return(genome_stats_dict)


def proteome_stats(prot_dir):
    dirs = ['reference', 'outgroup', 'user']
    cds_count_dict = {}

    for dir in dirs:
        if os.path.exists(prot_dir + "/" + dir):
            for filename in os.listdir(prot_dir + "/" + dir):
                if filename.endswith(".faa"):
                    basename = os.path.splitext(filename)[0]
                    cds_count = 0
                    with open(prot_dir + "/" + dir + "/" + filename, "rU") as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            cds_count += 1
                    cds_count_dict[basename] = cds_count

    return(cds_count_dict)

def profile(db_dir, ref_dir):
    prot_dir = ref_dir + "/" + "proteomes"
    genome_dir = ref_dir + "/genomes"
    cds_count_dict = proteome_stats(prot_dir)
    genome_stats_dict = genome_stats(genome_dir)

    df = pd.DataFrame.from_dict(genome_stats_dict, orient='index')
    df.columns = ['Genome Size', '# Contigs', 'GC Content']

    print(df)
