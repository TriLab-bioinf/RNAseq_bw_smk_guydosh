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
