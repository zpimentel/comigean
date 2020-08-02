"""
Author: Zachary Pimentel
Institution: University of Rhode Island
"""
import os
import wget
import subprocess

# https://github.com/AstrobioMike/GToTree/blob/master/hmm_sets/Bacteria.hmm

def unzip_db(outdir):
    cmd = ['unzip', outdir + '/taxdmp.zip', '-d', outdir]

    try:
        log_file = open("logfile", 'a')
        unzip_call = subprocess.call(cmd, shell=False,
                   stdout=log_file, stderr=subprocess.STDOUT)
        log_file.close()
    except OSError as e:
        raise Exception(f"The command {' '.join(cmd)} for database unpacking failed.")


def download_db(address, outdir, name):
    try:
        print(f"Downloading {name}:")
        out = wget.download(address, out=outdir)
        print("\n")
    except OSError as e:
        raise Exception(f"Downloading {name} failed. Note: This step requires interent access.")


def get_db(outdir):
    tax_ftp_address = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"
    refseq_ftp_address = "https://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

    download_db(tax_ftp_address, outdir, "NCBI Taxonomy Database")
    download_db(refseq_ftp_address, outdir, "NCBI RefSeq Summary")

    unzip_db(outdir)
