"""
Author: Zachary Pimentel
Institution: University of Rhode Island
"""

# calculate genome size
# calculate the number of CDS
# coding density
# ANI
# AAI
# gene presence/absence
# full taxonomy

import os
import sys
from Bio import SeqIO


def compute_genome_size(genome_dir):
    dirs = ['reference', 'outgroup', 'user']
    genome_size_dict = {}

    for dir in dirs:
        if os.path.exists(genome_dir + "/" + dir):
            for filename in os.listdir(genome_dir + "/" + dir):
                if filename.endswith(".fna"):
                    basename = os.path.splitext(filename)[0]
                    genome_size = 0
                    with open(genome_dir + "/" + dir + "/" + filename, "rU") as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            genome_size += len(record.seq)
                    genome_size_dict[basename] = genome_size

    return(genome_size_dict)


def count_CDS(prot_dir):
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


def make_barplot(ref_dir, feature_dict, db_dir, feature):
    outdir = ref_dir + "/iTol_Viz"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    bar_base_path = sys.path[0] + "/" + "base/" + "barchart_base.txt"
    out_handle = outdir + f"/{feature}_Viz.txt"
    out_file = open(out_handle, 'w')

    with open(bar_base_path) as bar_file:
        for row in bar_file:
            if "label 1" in row:
                row = row.replace("label 1", feature)
            out_file.write(row)

    name_dict = {}
    assembly_file = db_dir + "/" + "assembly_summary_refseq.txt"
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
            strain_name = species_name + " " + strain_info
            strain_name = strain_name.replace(" ", "_").replace(',', '').replace(';', ''). \
            replace(":", "").replace("(", " ").replace(")", " "). \
            replace("/", "").replace("_", " ").replace("[", "").replace("]", ""). \
            replace("'", "").replace("  ", "")
            name_dict[full_name] = strain_name

    for genome in feature_dict:
        genome_reduced = genome.replace("_protein", "").replace("_genomic", "")
        out_file.write(f"'{name_dict[genome_reduced]}', {feature_dict[genome]}\n")


def profile(db_dir, ref_dir):
    prot_dir = ref_dir + "/" + "proteomes"
    genome_dir = ref_dir + "/genomes"
    cds_count_dict = count_CDS(prot_dir)
    genome_size_dict = compute_genome_size(genome_dir)

    make_barplot(ref_dir, cds_count_dict, db_dir, "CDS")
    make_barplot(ref_dir, genome_size_dict, db_dir, "Genome_Size")
