# Bulk RNAseq processing pipeline
Bulk paired-end RNAseq workflow for Biowulf based on STAR read mapper.


## To build conda environment and install required tools run the following command:
```
mamba env create -f=environment.yml -n rnaseq
```
See environment.yml file for tool versions

## Generate a profile directory following the instructions here
```
https://github.com/NIH-HPC/snakemake_profile/tree/main
```

## Create the following directories within the working directory
```
mkdir data/00adapters
mkdir data/00reads
mkdir data/00ref
```

Folder data/00adapters/ should contain the fasta file with the Illumina d sequencing adapters. The name of the file should also be specified in the config.yml file.
Folder data/00reads/ should contain the gzipped fastq files
Folder data/00ref/ should contain the reference genome sequence, the annotation gtf file, a fasta file with RNA sequences to be filtered out from the RNAseq data (e.g. rRNA sequences), and will also host the STAR genome database, which will be created by the pipeline if missing. In that case, you need to specify the STAR sjdbOverhang length parameter for creating the database (usually read length - 1).   

## Configuration file
The pipeline configuration file, config.yml, is located in the config folder. Sample prefixes should be included in this file under "samples:".    

## To make a dry run of the entire snakemake pipeline
```
conda activate rnaseq
snakemake --profile slurm Snakefile -p -n 
```
Where *slurm*  is the name of the directory containing your profile 

## To run the entire snakemake pipeline
```
conda activate rnaseq
snakemake --profile slurm Snakefile -p
```

## To run a specific rule within the pipeline
First, try a dry run to check that everything works as expected
```
snakemake -R --until $MY_RULE --cores $CPUs -n
```
Then run the pipeline with the command below:
```
snakemake -R --until $MY_RULE --cores $CPUs
```
