"""
Author: Zachary Pimentel
Institution: University of Rhode Island
"""

from collections import defaultdict
import os
import subprocess

import wget


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

            url_dict[species_taxa_id].append((accession, url, assembly_level,
                                              assem_id, refseq_category))

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


def count_genomes(taxa_ids, url_dict, refseq_list, assem_list):
    """ Count the genomes to be downloaded. """
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

    def download_sequences(self, dir, taxa_ids, url_dict, refseq_list,
                           assem_list):
        """ Download genomes & proteomes from NCBI RefSeq ftp. """
        protdir = f"{dir}/proteomes/{self.group}/"
        mkdir_if_nonexistant(protdir)

        gendir = f"{dir}/genomes/{self.group}/"
        mkdir_if_nonexistant(gendir)

        for id in taxa_ids:
            for gen in url_dict[id]:
                if gen[2] in assem_list and gen[4] in refseq_list:
                    proturl = os.path.join(gen[1], gen[0] + "_" + gen[3] +
                                           "_protein.faa.gz").replace(" ", "_")

                    gen_url = os.path.join(gen[1], gen[0] + "_" + gen[3] +
                                           "_genomic.fna.gz").replace(" ", "_")
                    try:
                        proteome_out = wget.download(proturl, out=protdir,
                                                     bar=None)
                        genome_out = wget.download(gen_url, out=gendir,
                                                   bar=None)
                    except IOError:
                        print(f"Download of {gen[0]} failed.")

                    prot_file = os.path.join(protdir, proturl.split('/')[-1])
                    subprocess.run([f"gunzip {prot_file}"], shell=True)

                    gen_file = os.path.join(gendir, gen_url.split('/')[-1])
                    subprocess.run([f"gunzip {gen_file}"], shell=True)

        if len(os.listdir(protdir)) == 0:
            raise Exception("While a valid NCBI taxonomy ID was provided, no \
                             genomes or proteomes were downloaded. Do you  \
                             have a stable interent connection?")

        return(protdir)

    def call_genes(self, dir, gen_path, code):
        """ Call genes using prodigal. """
        dir = os.path.join(dir, self.group)
        mkdir_if_nonexistant(dir)

        if code:
            if str(code) not in ["{:01}".format(n) for n in range(1, 34)]:
                raise Exception("An invalid genetic code was provided.")
        else:
            code = 11

        log_file = open("logfile", 'a')
        for filename in os.listdir(gen_path):
            if filename.endswith(".fna"):
                bfna = filename
                base = os.path.splitext(bfna)[0]
            else:
                raise Exception("User provided genomes files must end \
                                 with the extension .fna")

            gene_call = subprocess.run([f'prodigal \
                                        -i {os.path.join(gen_path, bfna)} \
                                        -a {os.path.join(dir, base + ".faa")} \
                                        -g {code}'],
                                       shell=True, stdout=log_file,
                                       stderr=subprocess.STDOUT)

        log_file.close()


def run_group(parent, group_name, name_dict, node_dict, url_dict, rs_list,
              assem_list, count, dir):
    """
    Call functions and classes necessary to get all genomes within
    a particular taxa (defined by the taxonomic id).
    """
    for taxa_id in parent.split(","):
        group = GetGenomesClass(taxa_id, group_name)

        ids = group.get_lineage(name_dict, node_dict)

        genome_count = count_genomes(ids, url_dict, rs_list, assem_list)

        if count:
            print(f"Number of outgroup genomes from {taxa_id} to be \
                  downloaded: {genome_count}")
        else:
            if genome_count == 0:
                raise Exception(f"This taxa ID ({parent}) and options \
                                selected for --refseq_category and \
                                --assembly_level has resulted in no genomes \
                                being downloaded.")
            group.download_sequences(dir, ids, url_dict, rs_list, assem_list)


def get_genomes(ref_parent, out_parent, user_genomes, dir, dbdir, count,
                refseq_category, assembly_level, code):
    """ Prepare for genome download and gene calling. """
    prot_dir = f"{dir}/proteomes/"
    mkdir_if_nonexistant(prot_dir)

    rs_opt = ["reference genome", "representative genome", "na"]
    refseq_list = set_default(rs_opt, rs_opt, refseq_category, 2)

    assem_opt = ["Chromosome", "Complete Genome", "Scaffold", "Contig"]
    assem_list = set_default(assem_opt, ["Complete Genome"], assembly_level, 1)

    if ref_parent or out_parent:
        name_dict = parse_names(f"{dbdir}/names.dmp")
        node_dict = parse_nodes(f"{dbdir}/nodes.dmp")
        url_dict = parse_refseq(f"{dbdir}/assembly_summary_refseq.txt")

    if ref_parent:
        run_group(ref_parent, "reference", name_dict, node_dict, url_dict,
                  refseq_list, assem_list, count, dir)

    if out_parent:
        run_group(out_parent, "outgroup", name_dict, node_dict, url_dict,
                  refseq_list, assem_list, count, dir)

    if user_genomes:
        user = GetGenomesClass(None, "user")
        user.call_genes(prot_dir, user_genomes, code)
