# comigean 
A platform for **co**mparitive **mi**crobial **ge**nome **an**alysis

### Installation
This package uses conda to download pre-requisities. Therefore, it is required to install comigean. Clone this git repo and change directory into it.
```
>$ conda env create -f environment.yml
>$ conda activate comigean
>$ python setup.py install
```

### Step 1. Download Required Reference Databases
Downloads NCBI taxonomy database, NCBI RefSeq genome summaries, and a marker gene database from GtoTree (Lee M, 2020). Just specifiy the name of a directory (that does not yet exist) and these will be downloaded.
```
>$ comigean install-db <REF_DIR>
```

### Step 2. Download Genomes
This function will download all of the available genomes given a NCBI Taxonomy ID from the NCBI RefSeq database. This can be done for (1) a reference taxa and (2) a potential outgroup taxa (for the construction of a phylogenetic tree) at the same time. This command can also accept the path to a directory containg genomes from the user (ending in .fna) - in this case genes will be predicted using Prodigal.  
  
Required positional arguments include the name of an output directory to put the genomes (that does not exist) and the database directory (created in Step 1). In addition, at least one of the following three are required: a reference taxa ID (--REF_TAXA), an outgroup taxa ID (--OUTGROUP_TAXA), or a path to a user (--USER_GENOMES). However, all 3 can also be specified at the same time.
```
>$ comigean get-genomes <OUTDIR> <REF_DIR> --REF_TAXA <ID1> --OUTGROUP_TAXA <ID2> --USER_GENOMES <GENOME_DIR>
```

### Step 3. Profile Genomes
This function will produce basic statistics for each genome including genome size, number of contigs, and GC content. This data will be printed to standard out as a table.
```
>$ comigean genome-stats <OUTDIR> <REF_DIR>
```

### Step 4. Perform a Phylogenomic Analysis
This function identifies phylogenetic marker genes with HMMER (based on a marker set from Lee M. 2020), extracts the sequences from the proteomes, indiviually aligns each marker, concatenates the alignments, constructs a phylogeny with FastTree, and renames the leaves based on the strain names of each genome.  
  
The database directory (REF_DIR) needs to be specified as does the directory with the proteomes for each strain. The proteomes can be found in <OUTDIR>/proteomes.
```
>$ comigean find-markers <PROTEOME_DIR> <REF_DIR>
``` 
