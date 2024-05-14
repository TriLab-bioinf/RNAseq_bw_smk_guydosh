# vim: set ft=python:

# RNAseq workflow v1.0
# Hernan Lorenzi
# hernan.lorenzi@nih.gov
# Workflow requires to have a STAR index already available within the data/00ref directory
# Also, it is necessary to configure the config.yml file accordingly to include all metadatata required.

import os
import glob

configfile: "config/config.yml"
samples = config["samples"].keys()
genome = config["reference"]["genome_file"]
annotation = config["reference"]["ensembl_gtf"]
starOverhang = config["star_db"]["sjdbOverhang"]
starIndex = config["star_db"]["star_index"]
removeDupReads = config["remove_duplicated_reads"]

# Functions
def get_fq1(wildcards):
            return [ my_files for my_files in glob.glob(f"data/00reads/{wildcards.sample}*_R1*.fastq.gz")]

def get_fq2(wildcards):
            b = [ my_files for my_files in glob.glob(f"data/00reads/{wildcards.sample}*_R1*.fastq.gz")]
            c = list(map(lambda x: str.replace(x, "_R1", "_R2"), b))
            return c

# Set what rules to run locally
localrules: all #,
            #build_abundant_db

rule all:
    input:  "results/07multiqc/multiqc_report.html"

rule trimming:
    input:  fq1 = get_fq1,
            fq2 = get_fq2
    output:
            fq1P = "results/01trim/{sample}.1P.fastq.gz",
            fq2P = "results/01trim/{sample}.2P.fastq.gz",
            fqU = "results/01trim/{sample}.U.fastq.gz"
    params: 
            "ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 overwrite=t"
    resources:
        cpus_per_task = 16,
        partition = "quick",
        time = "4:00:00"
    threads: 16
    log:    log1 = "results/01trim/{sample}.log",
            log2 = "results/01trim/{sample}.stats.log"
    benchmark:
            "benchmarks/trim/{sample}.tsv"
    shell:
        """
        bbduk.sh -Xmx1g threads={threads} \
            in1={input.fq1} in2={input.fq2} \
            out1={output.fq1P} out2={output.fq2P} outs={output.fqU} \
            ref=data/00adapters/truseq.fa.gz \
            {params} stats={log.log2} 2> {log.log1}
        """

rule build_abundant_db:
    input: config["reference"]["abundant_rna_file"]
    output: "data/00ref/abundant"
    shell:
        """
        bowtie2-build {input} {output}
        touch {output}
        """

rule rm_abundant_rnas:
    input:  fq1P = "results/01trim/{sample}.1P.fastq.gz",
            fq2P = "results/01trim/{sample}.2P.fastq.gz",
            rnas_idx = "data/00ref/abundant"
    output:
        fq1F = "results/02abundant/{sample}.fastq.1.gz",
        fq2F = "results/02abundant/{sample}.fastq.2.gz"
    
    threads: 16
    resources:
        cpus_per_task = 16,
        partition = "norm",
        time = "24:00:00",
        gres = "lscratch:20"
    log: 
        metrics = "results/02abundant/{sample}.metrics.txt",
        logs = "results/02abundant/{sample}.log"
    benchmark:
            "benchmarks/rm_abundant_rnas/{sample}.tsv"
    shell:
        """
        bowtie2 --threads {threads} -L 20 -x {input.rnas_idx} \
         --met-file {log.metrics} \
         --un-conc-gz results/02abundant/{wildcards.sample}.fastq.gz \
         -1 {input.fq1P} -2 {input.fq2P} 2> {log.logs} 1> /dev/null
        """ 
       
rule make_star_index:
    input: gen = f"{genome}", 
        ann = f"{annotation}" 
    output: db = "data/00ref/SA"
    resources: 
        mem_mb = 64000,
        cpus_per_task = 16,
        partition = "quick",
        time = "3:00:00",
        gres = "lscratch:20"
    threads: 16
    params: so = f"{starOverhang}",
            db_dir = f"{starIndex}"
    shell:
        """
            STAR --runMode genomeGenerate \
             --genomeDir {params.db_dir} \
             --genomeFastaFiles {input.gen} \
             --sjdbGTFfile {input.ann} \
             --sjdbOverhang {params.so}

            touch {output.db}  
        """

rule map_reads:
    input: fq1 = "results/02abundant/{sample}.fastq.1.gz",
           fq2 = "results/02abundant/{sample}.fastq.2.gz",
           database = "data/00ref/SA"
    output: bam = "results/03map_reads/{sample}.Aligned.sortedByCoord.out.bam",
            log = "results/03map_reads/{sample}.Log.final.out"
    threads: 16
    resources:
        cpus_per_task = 16,
        partition = "norm",
        time = "14:00:00",
        mem_mb = 32000,
        gres = "lscratch:20"
    benchmark:
        "benchmarks/map_reads/{sample}.tsv"
    params: genome_dir = f"{starIndex}",
            prefix = "results/03map_reads/{sample}."
    shell:
        """
        STAR --runMode alignReads \
                --runThreadN {threads} \
                --genomeDir {params.genome_dir} \
                --alignSJDBoverhangMin 1 \
                --alignSJoverhangMin 5 \
                --outFilterMismatchNmax 2 \
                --alignEndsType EndToEnd \
                --readFilesIn {input.fq1} {input.fq2} \
                --readFilesCommand zcat --outFileNamePrefix {params.prefix} \
                --quantMode GeneCounts \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes All
        touch {output.bam}
        touch {output.log}
        """

rule remove_duplicates:
    input: "results/03map_reads/{sample}.Aligned.sortedByCoord.out.bam"
    output: "results/04dedup/{sample}.sorted.dedup.bam"
    params: f"READ_NAME_REGEX=null REMOVE_DUPLICATES={removeDupReads}"
    log: "results/04dedup/{sample}.sorted.dedup.metrics.txt"
    benchmark:
        "benchmarks/remove_duplicates/{sample}.tsv"
    resources:
        cpus_per_task = 4,
        mem_mb = 128000,
        partition = "quick",
        time = "4:00:00",
        gres = "lscratch:40"
    shell:
        """
        picard -Xmx32g MarkDuplicates \
         I={input} \
         O={output} \
         M={log} \
         {params}
        samtools index {output}
        """
rule make_bigwig:
    input: "results/04dedup/{sample}.sorted.dedup.bam"
    output: "results/05bigwig/{sample}.bw"
    params: "--binSize 10 --normalizeUsing BPM" #  + "--filterRNAstrand [forward/reverse]" to plot strand-specific data
    threads: 8
    resources:
        cpus_per_task = 8,
        partition = "norm",
        time = "14:00:00"
    shell:
        """
        bamCoverage -p {threads} -b {input} -o {output} {params}
        """

rule counts:
    input: 
        genome = config["reference"]["genome_file"],
        gtf = config["reference"]["ensembl_gtf"],
        bam = expand("results/04dedup/{s}.sorted.dedup.bam", s=samples)
    output: counts = "results/05counts/read_counts",
            summary = "results/05counts/read_counts.summary"
    params: "-t CDS -g gene_id -O -s 2 -J -R BAM -p --ignoreDup -M --fraction"  # Current params ignore multimappers and duplicated reads
                                                                   # -p  --countReadPairs = count fragments instead of individual reads
                                                                   # -M = include multi-mapping reads -O count reads mapping overlapping features
                                                                   # --fraction = multimapped reads will be caused as a fraction 
                                                                   #              instead of 1 (1/x where x = numb alignments reported for same read)
                                                                   # -s = stranded [0 = unstranded ; 1 = forward -stranded ; 2 = reverse-stranded]
    benchmark:
        "benchmarks/counts/counts.tsv"
    threads: 16
    resources:
        cpus_per_task = 16,
        partition = "norm",
        time = "14:00:00",
        mem_mb = 32000
    shell:
        """
        featureCounts {params} -G {input.genome} -T {threads}\
         -a {input.gtf} \
         -o {output.counts} {input.bam}
        """

rule fastqc:
    input: raw1 = get_fq1,
           raw2 = get_fq2,
           trim1p = "results/01trim/{sample}.1P.fastq.gz",
           trim_u = "results/01trim/{sample}.U.fastq.gz",
           trim2p = "results/01trim/{sample}.2P.fastq.gz"
    output: 
            o1 = "results/06fastqc_raw/{sample}_R1_001_fastqc.html",
            o2 = "results/06fastqc_raw/{sample}_R2_001_fastqc.html",
            o3 = "results/06fastqc_trim/{sample}.1P_fastqc.html",
            o5 = "results/06fastqc_trim/{sample}.2P_fastqc.html",
    threads: 8
    resources:
        cpus_per_task = 8,
        partition = "quick",
        time = "4:00:00",
        mem_mb = 4000
    params:
        "--quiet"
    shell:
        """
        fastqc {params} -t {threads} -o results/06fastqc_raw {input.raw1} {input.raw2}
        fastqc {params} -t {threads} -o results/06fastqc_trim {input.trim1p} {input.trim2p}
        """
absolute_path = "/home/lorenziha/Downloads/snakemake-class/workflow/"


rule multiqc:
    input: 
           i1 = expand("results/01trim/{sample}.log", sample = samples),
           i2 = expand("results/02abundant/{sample}.log", sample = samples), 
           i3 = expand("results/03map_reads/{sample}.Log.final.out", sample = samples),
           i4 = expand("results/04dedup/{sample}.sorted.dedup.metrics.txt", sample = samples),
           i5 = "results/05counts/read_counts.summary",
           i6 = expand("results/05bigwig/{sample}.bw", sample = samples),
           i7 = expand("results/06fastqc_raw/{sample}_R{pair}_001_fastqc.html", sample = samples, pair = [1,2]),
           i8 = expand("results/06fastqc_trim/{sample}.{pair}P_fastqc.html", sample = samples, pair = [1,2])
    output: "results/07multiqc/multiqc_report.html"
    resources:
        partition = "quick",
        time = "4:00:00",
        mem_mb = 4000                    
    shell:
        """
        multiqc -f -d -o results/07multiqc results/
        """

