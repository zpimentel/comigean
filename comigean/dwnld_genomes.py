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
    # there can be many genomes associated with one taxa id
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


def set_defaults(param, param_list, default):
    """
    Set default parameters for ncbi-genome-download if None or do not match
    list of arguments.
    """

    if param:
        if param not in param_list:
            param = default
            print(f'({param} was incorrectly specified thus the default {default} was applied.)')
    else:
        param = default

    return(param)


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


    def count_genomes(self, taxa_ids, url_dict, refseq_category, assembly_level):
        """ Count the genomes to be downloaded. """
        if refseq_category:
            if refseq_category == "all":
                refseq_list = ["reference genome", "representative genome", "na"]
            else:
                refseq_list = []
                for cat in refseq_category.split(","):
                    refseq_list.append(cat + " genome")
        else:
            refseq_list = ["reference genome", "representative genome", "na"]

        if assembly_level:
            if assembly_level == "all":
                assem_list = ["Chromosome", "Complete Genome", "Scaffold", "Contig"]
            else:
                assem_list = []
                for cat in assembly_level.split(","):
                    if cat == "Complete_Genome":
                        cat = "Complete Genome"
                    assem_list.append(cat)
        else:
            assem_list = ["Complete Genome"]

        #assembly_level = [assem.replace("Complete_Genome", "Complete Genome") for assem in assembly_level]

        count = 0
        for id in taxa_ids:
            for genome_data in url_dict[id]:
                if genome_data[2] in assem_list:
                    if genome_data[4] in refseq_list:
                        count += 1

        return(count)


    def download_proteomes(self, dir, taxa_ids, url_dict, refseq_category, assembly_level):
        """ Download proteomes. """

        protdir = f"{dir}/proteomes/{self.group}/"
        if not os.path.exists(protdir):
            os.makedirs(protdir)

        genomedir = f"{dir}/genomes/{self.group}/"
        if not os.path.exists(genomedir):
            os.makedirs(genomedir)

        if refseq_category:
            if refseq_category == "all":
                refseq_list = ["reference genome", "representative genome", "na"]
            else:
                refseq_list = []
                for cat in refseq_category.split(","):
                    refseq_list.append(cat + " genome")
        else:
            refseq_list = ["reference genome", "representative genome", "na"]

        if assembly_level:
            if assembly_level == "all":
                assem_list = ["Chromosome", "Complete Genome", "Scaffold", "Contig"]
            else:
                assem_list = []
                for cat in assembly_level.split(","):
                    if cat == "Complete_Genome":
                        cat = "Complete Genome"
                    assem_list.append(cat)
        else:
            assem_list = ["Complete Genome"]

    #    assembly_level = [assem.replace("Complete_Genome", "Complete Genome") for assem in assembly_level]

        count = 0
        for id in taxa_ids:
            for genome_data in url_dict[id]:
                try:
                    if genome_data[2] in assem_list:
                        if genome_data[4] in refseq_list:
                            count += 1
                            sys.stdout.flush()
                            sys.stdout.write(f"\rDownloading {self.group} genome & proteome {genome_data[0]}")

                            proteome_url = os.path.join(genome_data[1], genome_data[0] + "_" + genome_data[3] + "_protein.faa.gz")
                            genome_url = os.path.join(genome_data[1], genome_data[0] + "_" + genome_data[3] + "_genomic.fna.gz")

                            genome_url = genome_url.replace(" ", "_")
                            proteome_url = proteome_url.replace(" ", "_")

                            proteome_out = wget.download(proteome_url, out=protdir, bar=None)
                            subprocess.run([f"gunzip {os.path.join(protdir, proteome_url.split('/')[-1])}"], shell=True)

                            genome_out = wget.download(genome_url, out=genomedir, bar=None)
                            subprocess.run([f"gunzip {os.path.join(genomedir, genome_url.split('/')[-1])}"], shell=True)

                except:
                    print(f"Download of {genome_data[0]} failed.")

        list = os.listdir(protdir)
        if len(list) == 0:
            raise Exception("Nothing was downloaded. Did you provide a valid NCBI Taxonomy ID?")

        return(protdir)


    def call_genes(self, dir, genome_path, code):
        """ Call genes using prodigal. """

        outdir = os.path.join(dir, self.group)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        if not code:
            code = 11

        print(f"Calling genes for {self.group} genomes:")
        log_file = open("logfile", 'a')
        for filename in os.listdir(genome_path):
            if self.group == "user":
                if filename.endswith(".fna"):
                    basename_fna = filename
                    basename = os.path.splitext(basename_fna)[0]
                    print(f"Calling genes for {os.path.join(genome_path, basename_fna)}")
                    gene_call = subprocess.run([f'prodigal -i {os.path.join(genome_path, basename_fna)} -a {os.path.join(outdir, basename + ".faa")} -g {code}'], shell=True,
                               stdout=log_file, stderr=subprocess.STDOUT)
            else:
                if filename.endswith(".gz"):
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
    if not os.path.exists(prot_dir):
        os.makedirs(prot_dir)

    if ref_parent or out_parent:
        name_dict = parse_names(f"{dbdir}/names.dmp")
        node_dict = parse_nodes(f"{dbdir}/nodes.dmp")
        url_dict = parse_refseq(f"{dbdir}/assembly_summary_refseq.txt")

    if ref_parent:
        for ref_taxa_id in ref_parent.split(","):
            ref = GetGenomesClass(ref_taxa_id, "reference")

            ref_ids = ref.get_lineage(name_dict, node_dict)

            if count:
                ref_count = ref.count_genomes(ref_ids, url_dict, refseq_category, assembly_level)
                print(f"Number of reference genomes from {ref_taxa_id} to be downloaded: {ref_count}")
            else:
                ref_genomes = ref.download_proteomes(dir, ref_ids, url_dict, refseq_category, assembly_level)

    if out_parent:
        for out_taxa_id in out_parent.split(","):
            out = GetGenomesClass(out_taxa_id, "outgroup")

            out_ids = out.get_lineage(name_dict, node_dict)

            if count:
                out_count = out.count_genomes(out_ids, url_dict, refseq_category, assembly_level)
                print(f"Number of outgroup genomes from {out_taxa_id} to be downloaded: {out_count}")
            else:
                out_genomes = out.download_proteomes(dir, out_ids, url_dict, refseq_category, assembly_level)

    if not count:
        sys.stdout.write("\n")

    if user_genomes:
        user = GetGenomesClass(None, "user")

        user.call_genes(prot_dir, user_genomes, code)
