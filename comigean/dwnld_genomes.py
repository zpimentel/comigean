"""
Author: Zachary Pimentel
Institution: University of Rhode Island
"""

import os
import sys
import logging
import urllib
import argparse
import subprocess
import wget
from collections import defaultdict


def get_children(children, node_dict, child_list):
    """
    Recursive function that finds all species level taxa id's from a given
    taxa id.
    """
    for child in children:
        if child[1] != 'species':
            new_child = node_dict[child[0]]
            get_children(new_child, node_dict, child_list)
        else:
            child_list.append(child[0])

    return(child_list)


def parse_refseq(assembly_file):
    """ Get url of genomes to download. """
    url_dict = defaultdict(list)

    with open(assembly_file) as assem:
        next(assem)
        next(assem)
        for row in assem:
            row = row.strip().split("\t")
            accession = row[0]
            refseq_category = row[4]
            species_taxa_id = row[6]
            url = row[19]
            assembly_level = row[11]
            assem_id = row[15]

            url_dict[species_taxa_id].append((accession, url, assembly_level, assem_id, refseq_category))

    return(url_dict)


def parse_names(names_handle):
    """ Get names associated with each taxa id. """
    name_dict = {}

    names_file = open(names_handle)
    for row in names_file:
        if 'scientific name' not in row:
            continue
        row = row.strip().split('\t|\t')
        tax_id = row[0]
        name = row[1]
        name_dict[tax_id] = name

    return(name_dict)


def set_default(options, default, argument, track):
    """ Set default when argument is not provided. """
    if argument:
        if argument == "all":
            option_list = options
        else:
            option_list = []
            for cat in argument.split(","):
                if cat == "Complete_Genome":
                    cat = "Complete Genome"
                if track == 2:
                    cat = cat + " genome"
                if cat not in options:
                    raise Exception(f"{cat} is not a valid option.")
                option_list.append(cat)
    else:
        option_list = default

    return(option_list)


def mkdir_if_nonexistant(dir):
    """ Makes directory if it does not exist. """
    if not os.path.exists(dir):
        os.makedirs(dir)


def parse_nodes(node_handle):
    """ Find child of each parent node. """
    node_dict = defaultdict(list)

    nodes_file = open(node_handle)
    for row in nodes_file:
        row = row.strip().split('\t|\t')
        child = row[0]
        parent = row[1]
        tax_level = row[2]
        node_dict[parent].append((child, tax_level))

    return(node_dict)


def count_genomes(taxa_ids, url_dict, refseq_category, assembly_level):
    """ Count the genomes to be downloaded. """
    refseq_list = set_default(["reference genome", "representative genome", "na"],
                              ["reference genome", "representative genome", "na"],
                              refseq_category, 2)
    assem_list = set_default(["Chromosome", "Complete Genome", "Scaffold", "Contig"],
                             ["Complete Genome"], assembly_level, 1)

    count = 0
    for id in taxa_ids:
        for genome_data in url_dict[id]:
            if genome_data[2] in assem_list:
                if genome_data[4] in refseq_list:
                    count += 1

    return(count)


class GetGenomesClass:
    '''
    Commands related to getting genomes for all strains beneath
    specified taxa id.
    '''

    def __init__(self, taxa_id, group):
        self.taxa_id = taxa_id
        self.group = group

    def get_lineage(self, name_dict, node_dict):
        """ Get all species level taxa ids from parent taxa id. """
        child_list = []

        child = node_dict[self.taxa_id]
        child_list = get_children(child, node_dict, child_list)
        if len(child_list) < 1:
            child_list = [self.taxa_id]

        return(child_list)

    def download_sequences(self, dir, taxa_ids, url_dict, refseq_category, assembly_level):
        """ Download genomes & proteomes from NCBI RefSeq ftp. """
        protdir = f"{dir}/proteomes/{self.group}/"
        mkdir_if_nonexistant(protdir)

        genomedir = f"{dir}/genomes/{self.group}/"
        mkdir_if_nonexistant(genomedir)

        refseq_list = set_default(["reference genome", "representative genome", "na"],
                                  ["reference genome", "representative genome", "na"],
                                  refseq_category, 2)
        assem_list = set_default(["Chromosome", "Complete Genome", "Scaffold", "Contig"],
                                 ["Complete Genome"], assembly_level, 1)

        count = 0
        for id in taxa_ids:
            for genome_data in url_dict[id]:
                try:
                    if genome_data[2] in assem_list:
                        if genome_data[4] in refseq_list:
                            count += 1
                            sys.stdout.flush()
                            sys.stdout.write(f"\rDownloading {self.group} genome & proteome {genome_data[0]}")

                            proteome_url = os.path.join(genome_data[1], genome_data[0] + "_" + genome_data[3] + "_protein.faa.gz").replace(" ", "_")
                            genome_url = os.path.join(genome_data[1], genome_data[0] + "_" + genome_data[3] + "_genomic.fna.gz").replace(" ", "_")

                            proteome_out = wget.download(proteome_url, out=protdir, bar=None)
                            subprocess.run([f"gunzip {os.path.join(protdir, proteome_url.split('/')[-1])}"], shell=True)

                            genome_out = wget.download(genome_url, out=genomedir, bar=None)
                            subprocess.run([f"gunzip {os.path.join(genomedir, genome_url.split('/')[-1])}"], shell=True)
                except:
                    print(f"Download of {genome_data[0]} failed.")

        list = os.listdir(protdir)
        if len(list) == 0:
            raise Exception("While a valid NCBI taxonomy ID was provided, no genomes or proteomes were downloaded. Do you have a stable interent connection?")

        return(protdir)


    def call_genes(self, dir, genome_path, code):
        """ Call genes using prodigal. """
        outdir = os.path.join(dir, self.group)
        mkdir_if_nonexistant(outdir)

        if code:
            if str(code) not in ["{:01}".format(n) for n in range(1,34)]:
                raise Exception(f"{code} is not a valid genetic code.")
        else:
            code = 11

        print(f"Calling genes for {self.group} genomes:")
        log_file = open("logfile", 'a')
        for filename in os.listdir(genome_path):
            if self.group == "user":
                if filename.endswith(".fna"):
                    basename_fna = filename
                    basename = os.path.splitext(basename_fna)[0]
                else:
                    basename_fna = os.path.splitext(filename)[0]
                    basename = os.path.splitext(basename_fna)[0]
                    subprocess.run([f"gunzip {os.path.join(genome_path, filename)}"], shell=True)

                print(f"Calling genes for {os.path.join(genome_path, basename_fna)}")
                gene_call = subprocess.run([f'prodigal -i {os.path.join(genome_path, basename_fna)} -a {os.path.join(outdir, basename + ".faa")} -g {code}'], shell=True,
                           stdout=log_file, stderr=subprocess.STDOUT)

        log_file.close()


def get_genomes(ref_parent, out_parent, user_genomes, dir, dbdir, count, refseq_category, assembly_level, code):
    """
    Call functions and classes necessary to get all genomes within
    a particular taxa (defined by the taxonomic id). Also calls genes
    for all of the downloaded genomes.
    """
    prot_dir = f"{dir}/proteomes/"
    mkdir_if_nonexistant(prot_dir)

    if ref_parent or out_parent:
        name_dict = parse_names(f"{dbdir}/names.dmp")
        node_dict = parse_nodes(f"{dbdir}/nodes.dmp")
        url_dict = parse_refseq(f"{dbdir}/assembly_summary_refseq.txt")

    if ref_parent:
        for ref_taxa_id in ref_parent.split(","):
            ref = GetGenomesClass(ref_taxa_id, "reference")

            ref_ids = ref.get_lineage(name_dict, node_dict)

            ref_count = count_genomes(ref_ids, url_dict, refseq_category, assembly_level)

            if count:
                print(f"Number of reference genomes from {ref_taxa_id} to be downloaded: {ref_count}")
            else:
                if ref_count == 0:
                    raise Exception(f"This taxa ID ({ref_parent}) and options selected for --refseq_category and --assembly_level has resulted in no genomes being downloaded.")
                ref_genomes = ref.download_sequences(dir, ref_ids, url_dict, refseq_category, assembly_level)

    if out_parent:
        for out_taxa_id in out_parent.split(","):
            out = GetGenomesClass(out_taxa_id, "outgroup")

            out_ids = out.get_lineage(name_dict, node_dict)

            out_count = count_genomes(out_ids, url_dict, refseq_category, assembly_level)

            if count:
                print(f"Number of outgroup genomes from {out_taxa_id} to be downloaded: {out_count}")
            else:
                if out_count == 0:
                    raise Exception(f"This taxa ID ({out_parent}) and options selected for --refseq_category and --assembly_level has resulted in no genomes being downloaded.")
                out_genomes = out.download_sequences(dir, out_ids, url_dict, refseq_category, assembly_level)

    if not count:
        sys.stdout.write("\n")

    if user_genomes:
        user = GetGenomesClass(None, "user")
        user.call_genes(prot_dir, user_genomes, code)
