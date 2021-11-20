# comigean 
A platform for **co**mparitive **mi**crobial **ge**nome **an**alysis

### Installation
This package uses conda to download pre-requisities. Therefore, it is required to install comigean. Clone this git repo and change directory into it.
```
>$ conda env create -f environment.yml
>$ conda activate comigean
>$ python setup.py install
```

### 1. Download Required Reference Databases
Downloads NCBI taxonomy database, NCBI RefSeq genome summaries, and a marker gene database from GtoTree (Lee M, 2019). Just specifiy the name of a directory (that does not yet exist) and these will be downloaded.
```
>$ comigean install-db <REF_DIR>
```

### 2. Download Genomes
This function will download all of the available genomes given a NCBI Taxonomy ID from the NCBI RefSeq database. This can be done for (1) a reference taxa and (2) a potential outgroup taxa (for the construction of a phylogenetic tree) at the same time. This command can also accept the path to a directory containg genomes from the user (ending in .fna) - in this case genes will be predicted using Prodigal.  
  
Required positional arguments include the name of an output directory to put the genomes (that does not exist) and the database directory (created in Step 1). In addition, at least one of the following three are required: a reference taxa ID (--REF_TAXA), an outgroup taxa ID (--OUTGROUP_TAXA), or a path to a user (--USER_GENOMES). However, all 3 can also be specified at the same time.
```
>$ comigean get-genomes <OUTDIR> <REF_DIR> --REF_TAXA <ID1> --OUTGROUP_TAXA <ID2> --USER_GENOMES <GENOME_DIR>
```

### 3. Profile Genomes


### 4. Perform a Phylogenomic Analysis

