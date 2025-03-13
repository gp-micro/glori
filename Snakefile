configfile: "config.yaml"

RESULTS_DIR = config["results_dir"]
SAMPLE_TO_FASTQ = config["sample_to_fastq"]
SAMPLES = list(SAMPLE_TO_FASTQ.keys())

rule umitools_extract:
    input:
        FASTQ_DIR + "/{sample}.fastq.gz"
    output:
        RESULTS_DIR + "/umi_extracted/{sample}.extracted.fastq.gz"
    conda:
        "envs/umitools_env.yaml"
    log:
        "logs/umitools_extract/{sample}.log"
    shell:
        "#UMITOOLS"

rule cutadapt:
    input:
        RESULTS_DIR + "/umi_extracted/{sample}.extracted.fastq.gz"
    output:
        RESULTS_DIR + "/trimmed/{sample}.trimmed.fastq.gz"
    conda:
        "envs/cutadapt_env.yaml"
    log:
        "logs/cutadapt/{sample}.log"
    shell:
        "#CUTADAPT"

rule hisat2_mapping:
    input:
        fastq = RESULTS_DIR + "/trimmed/{sample}.trimmed.fastq.gz",
        index = HISAT2_INDEX
    output:
        RESULTS_DIR + "/bam/{sample}.bam"
    conda:
        "envs/hisat2_env.yaml"
    log:
        "logs/hisat2/{sample}.log"
    threads: 8
    shell:
        "#HISAT2"

rule samtools_sort:
    input:
        RESULTS_DIR + "/bam/{sample}.bam"
    output:
        RESULTS_DIR + "/sorted_bam/{sample}.sorted.bam"
    conda:
        "envs/samtools_env.yaml"
    threads: 4
    shell:
        "#SAMTOOLS SORT"
