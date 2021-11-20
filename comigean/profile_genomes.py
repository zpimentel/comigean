"""
Author: Zachary Pimentel
Institution: University of Rhode Island
"""

import os
import sys
import subprocess

from Bio import SeqIO
import pandas as pd

def make_genome_list(genome_dir):
    ''' Makes list of genomes for FastANI analysis '''
    dirs = ['reference', 'outgroup', 'user']
    file_handle = os.path.join(genome_dir, "genome_list.txt")
    f = open(file_handle, "w")

    for dir in dirs:
        if os.path.exists(genome_dir + "/" + dir):
            for filename in os.listdir(genome_dir + "/" + dir):
                if filename.endswith(".fna"):
                    genome_handle = os.path.join(genome_dir, dir, filename)
                    f.write(genome_handle + "\n")
    f.close()


def run_fastani(ref_dir):
    ''' Run the fastANI command '''
    log_file = open("logfile", 'a')
    run_fastani_cmd = subprocess.run([f'fastANI --ql {ref_dir}/genomes/genome_list.txt --rl {ref_dir}/genomes/genome_list.txt -o {ref_dir}/genomes/fast_ani_out.txt'], shell=True,
                                   stdout=log_file, stderr=subprocess.STDOUT)
    log_file.close()


def plot_fastani(ref_dir, name_dict):
    ''' Make a heatmap from the fastANI results '''
    ani_handle = os.path.join(ref_dir, "genomes/fast_ani_out.txt")

    fa_df = pd.read_csv(ani_handle, sep="\t", header=None)
    fa_df.columns = ['genome1', 'genome2', 'ani', '1', '2']

    fa_df['genome1'] = fa_df['genome1'].str.replace(f'{ref_dir}/genomes/reference/','')
    fa_df['genome2'] = fa_df['genome2'].str.replace(f'{ref_dir}/genomes/reference/','')
    fa_df['genome1'] = fa_df['genome1'].str.replace('_genomic.fna','')
    fa_df['genome2'] = fa_df['genome2'].str.replace('_genomic.fna','')

    fa_df['genome1'] = fa_df.genome1.map(name_dict)
    fa_df['genome2'] = fa_df.genome2.map(name_dict)

    fa_df = fa_df.pivot(index='genome1', columns='genome2', values='ani')
    p = fa_df.style.background_gradient(cmap='Blues')

    f=open(f"{ref_dir}/genomes/ANI_heatmap.html","w")
    f.write(p.render()) # df is the styled dataframe
    f.close()


def rename_heatmap(db_dir):
    name_dict = {}
    assembly_file = os.path.join(db_dir, "assembly_summary_refseq.txt")
    with open(assembly_file) as assem:
        next(assem)
        next(assem)
        for row in assem:
            row = row.strip().split("\t")
            accession = row[0]
            species_name = row[7]
            strain_info = row[8]
            assem_id = row[15]
            full_name = accession + "_" + assem_id
            full_name = full_name.replace(" ", "_")
            name_dict[full_name] = species_name + " " + strain_info

    return(name_dict)


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

    make_genome_list(genome_dir)
    run_fastani(ref_dir)
    name_dict = rename_heatmap(db_dir)
    #print(name_dict)
    plot_fastani(ref_dir, name_dict)

    print(df)
