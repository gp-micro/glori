configfile: "config.yaml"

RESULTS_DIR = config["results_dir"]
SAMPLE_TO_FASTQ = config["sample_to_fastq"]
SAMPLES = list(SAMPLE_TO_FASTQ.keys())

rule umitools_extract:
    input:
        lambda wildcards: SAMPLE_TO_FASTQ[wildcards.sample]
    output:
        RESULTS_DIR + "/umi_extracted/{sample}.fastq.gz"
    conda:
        "envs/umitools_env.yaml"
    log:
        "logs/umitools_extract/{sample}.log"
    shell:
        "umi_tools extract -I {input} -S {output} -p NNNNNNNNNNNN --log={log}"

rule cutadapt:
    input:
        RESULTS_DIR + "/umi_extracted/{sample}.fastq.gz"
    output:
        RESULTS_DIR + "/cutadapt_trimmed/{sample}.fastq.gz"
    conda:
        "envs/cutadapt_env.yaml"
    log:
        "logs/cutadapt/{sample}.log"
    shell:
        "cutadapt -a AGATCGGAAGAGCGTCGTG --max-n 0 --trimmed-only -e 0.1 -q 30 -m 30 --trim-n -o {output} {input} &> {log}"

rule hisat2_mapping:
    input:
        fastq = RESULTS_DIR + "/trimmed/{sample}.trimmed.fastq.gz"
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
