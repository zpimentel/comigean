# comigean 
A platform for **co**mparitive **mi**crobial **ge**nome **an**alysis

### Installation
This package uses conda to download pre-requisities. Therefore, it is required to install comigean. Clone this git repo and change directory into it.
```
>$ conda env create -f environment.yml
>$ conda activate comigean
>$ python setup.py install
```

### Download Required Databases
Downloads NCBI taxonomy database, NCBI RefSeq genome summaries, and a marker gene database from GtoTree (Lee M, 2019). Just specifiy the name of a directory (that does not yet exist) and these will be downloaded.
```
>$ comigean install-db <db-dir>
```

### Download Genomes


### Profile Genomes


### Perform a Phylogenomic Analysis

