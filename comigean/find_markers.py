"""
Author: Zachary Pimentel
Institution: University of Rhode Island
"""
import os
import glob
import collections
import subprocess
from Bio import SeqIO, AlignIO, Phylo
from pathlib import Path
from os.path import normpath, basename
# download the HMM files in first command! https://raw.githubusercontent.com/AstrobioMike/GToTree/master/hmm_sets/Bacteria.hmm
# then use hmmpress


def maketree(concat_aln, outdir):
    log_file = open("logfile", 'a')
    print(f'Running FastTree.')
    run_muscle = subprocess.run([f'FastTree {concat_aln} > {outdir}/Concat.tre'], shell=True, stdout=log_file, stderr=subprocess.STDOUT)
    log_file.close()


def rename_tree(outdir, ref_dir):

    name_dict = {}
    assembly_file = os.path.join(ref_dir, "assembly_summary_refseq.txt")
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

    tree_file = os.path.join(outdir, "Concat.tre")
    tree = Phylo.read(tree_file, 'newick')
    for node in tree.find_clades():
        if node.name:
            node.name = "_".join(str(node.name).split("_")[:-1])
            try:
                node.name = name_dict[node.name]
            except:
                print(f"Unable to rename {node.name}")
            node.name = node.name.replace(',', '').replace(';', ''). \
            replace(":", "").replace("(", " ").replace(")", " "). \
            replace("/", "").replace("_", " ").replace("[", "").replace("]", ""). \
            replace("'", "").replace("  ", "")

    Phylo.write(tree, os.path.join(outdir ,"Concat_Renamed.tre"), "newick")

def concatenate(outdir):
    genomes = []
    markers = []
    length_dict = {}
    hitdict = {}
    for file in glob.glob(outdir + '/*.aln'):
        marker = os.path.basename(os.path.splitext(file)[0])
        markers.append(marker)
        for x, seq_record in enumerate(AlignIO.read(file, "fasta")):
            genomeid = seq_record.id
            genomes.append(genomeid)
            seq = str(seq_record.seq)
            length_dict[marker] = len(seq)
            hitdict[marker, genomeid] = seq

    genomes = set(genomes)
    markers = set(markers)

    outfile_handle = f"{outdir}/Concat.aln"
    outfile = open(outfile_handle,"w")
    for genomeid in sorted(genomes):
        concatseq = ''
        for marker in sorted(markers):
            try:
                seq = hitdict[marker, genomeid].strip()
            except:
                seq = str("-" * length_dict[marker])
            concatseq += seq
        outfile.write(">{}\n{}\n".format(genomeid, concatseq))
    outfile.close()

    return(outfile_handle)


def align(outdir, hmm_count):
    log_file = open("logfile", 'a')
    for hmm in hmm_count:
        if hmm_count[hmm] > 0.99:
            basename = os.path.splitext(hmm)[0]
            print(f'Running muscle on {basename}')
            run_muscle = subprocess.run([f'muscle -in {os.path.join(outdir, basename + ".fa")} -out {os.path.join(outdir, basename + ".aln")}'], shell=True,
                               stdout=log_file, stderr=subprocess.STDOUT)
    log_file.close()


def get_sequences(prot_dir, protein2hit, outdir, hmm_hits):
    dirs = ['reference', 'outgroup', 'user']
    hmm_cnt = collections.Counter()

    proteome_count = 0
    for dir in dirs:
        if os.path.exists(os.path.join(prot_dir, dir)):
            for filename in os.listdir(os.path.join(prot_dir, dir)):
                if filename.endswith(".faa"):
                    proteome_count += 1
                    basename = os.path.splitext(filename)[0]
                    with open(os.path.join(prot_dir, dir, filename), "rU") as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            try:
                                hmm = protein2hit[(record.id, basename)]
                                hmm_hits[(basename, hmm)]
                                outfile_handle = f"{outdir}/{hmm}.fa"
                                outfile = open(outfile_handle, 'a')
                                outfile.write(f">{basename}\n{record.seq.rstrip('*')}\n")
                                outfile.close()
                                hmm_cnt[hmm] += 1
                            except:
                                pass
    for key in hmm_cnt:
        hmm_cnt[key] /= proteome_count

    return(hmm_cnt)


def parse_hmms(outdir):
    hmm_hits = {}
    protein2hit = {}
    for pfam in glob.glob(outdir + "/*.tsv"):
        pfamfile = open(pfam)
        assembly_id = os.path.splitext(pfam)[0]
        assembly_id = basename(normpath(assembly_id))

        for row in pfamfile:
            if row[0] == '#':
                continue
            if row[0] == ' ':
                continue
            row = row.strip().split(None, 22)
            domain_name = row[0]
            pfam_accession = row[1]
            seq_id = row[3]
            profile_length = float(row[2])
            query_length = float(row[5])
            evalue = float(row[12])  # i-evalue
            query_start = float(row[17])
            query_stop = float(row[18])
            profile_start = float(row[15])
            profile_stop = float(row[16])

            profile_cov = (profile_stop - profile_start + 1) / profile_length
            profile_cov = float(profile_cov)

            if profile_cov >= 0.7 and evalue <= 1E-5:
                protein2hit[(seq_id, assembly_id)] = domain_name
                try:
                    evalue2 = hmm_hits[(assembly_id, domain_name)][0] # only store the best hit
                    if evalue2 > evalue:
                        hmm_hits[(assembly_id, domain_name)] = ((evalue, seq_id))
                except:
                    hmm_hits[(assembly_id, domain_name)] = ((evalue, seq_id))

    return(hmm_hits, protein2hit)


def run_hmmer(prot_dir, ref_dir):
    dirs = ['reference', 'outgroup', 'user']
    log_file = open("logfile", 'a')

    outdir = os.path.join(str(Path(prot_dir).parents[0]), "hmm_out")
    os.makedirs(outdir)

    for dir in dirs:
        if os.path.exists(os.path.join(prot_dir, dir)):
            for filename in os.listdir(os.path.join(prot_dir, dir)):
                if filename.endswith(".faa"):
                    basename = os.path.splitext(filename)[0]
                    print(f'Running HMMER on {os.path.join(prot_dir, dir, filename)}')
                    run_hmmer = subprocess.run([f'hmmscan --domtblout {outdir}/{basename}.tsv {ref_dir}/Bacteria.hmm {os.path.join(prot_dir, dir, filename)}'], shell=True,
                               stdout=log_file, stderr=subprocess.STDOUT)

    log_file.close()

    return(outdir)


def phylogenomics(prot_dir, ref_dir):
    # switch this to outdir not protdir
    outdir = run_hmmer(prot_dir, ref_dir)

    hmm_hits, protein2hit = parse_hmms(outdir)
    hmm_count = get_sequences(prot_dir, protein2hit, outdir, hmm_hits)
    align(outdir, hmm_count)
    outfile = concatenate(outdir)

    maketree(outfile, outdir)
    rename_tree(outdir, ref_dir)
