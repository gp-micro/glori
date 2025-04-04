configfile: "config.yaml"

RESULTS_DIR = config["results_dir"]
SAMPLE_TO_FASTQ = config["sample_to_fastq"]
SAMPLES = list(SAMPLE_TO_FASTQ.keys())
HISAT2_PATH = config["hisat2_path"]
HISAT2_INDEX_DIR = config["hisat2_index_dir"]
REFERENCE_FASTA = config["reference_fasta"]
REFERENCE_GTF = config["reference_gtf"]
DB_PATH = config["db_path"]

#Note: I ran this using 240 GB of RAM
rule hisat2_index:
    input:
        fasta=REFERENCE_FASTA,
        gtf=REFERENCE_GTF
    output:
        directory(HISAT2_INDEX_DIR)
    conda:
        "envs/python_2.7.16.yaml"
    log:
        "logs/hisat2_index.log"
    threads: 16
    shell:
        "python GLORI_pipeline/scripts/A2G_hisat2_index.py -i {input.fasta} -p {threads} -o {output} --gtf {input.gtf} --hisat2-path {HISAT2_PATH} &> {log}"

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

rule decompress:
    input:
        RESULTS_DIR + "/cutadapt_trimmed/{sample}.fastq.gz"
    output:
        temp(RESULTS_DIR + "/cutadapt_trimmed_decompressed/{sample}.fastq")
    shell:
        "zcat {input} > {output}"

rule hisat2_mapping:
    input:
        index=directory(HISAT2_INDEX_DIR),
        fastq=RESULTS_DIR + "/cutadapt_trimmed_decompressed/{sample}.fastq"
    output:
        bam=RESULTS_DIR + "/hisat2/{sample}.bam",
        multimap=RESULTS_DIR + "/hisat2/{sample}.multimappers.bam",
    conda:
        "envs/python_2.7.16.yaml"
    log:
        "logs/hisat2/{sample}.log"
    threads: 8
    shell:
        "python GLORI_pipeline/scripts/A2G_hisat2.py -F {input.fastq} -o " +
            RESULTS_DIR + "/hisat2/{wildcards.sample} -I {input.index} --index-prefix HISAT2 " +
            "--hisat2-path {HISAT2_PATH} --del-convert --del-sam >& {log}"

rule samtools_sort_index:
    input:
        RESULTS_DIR + "/hisat2/{sample}.bam",
    output:
        bam=RESULTS_DIR + "/hisat2/{sample}.sorted.bam",
        bai=RESULTS_DIR + "/hisat2/{sample}.sorted.bam.bai"
    conda:
        "envs/python_2.7.16.yaml"
    threads: 1
    shell:
        "samtools sort -m 10G -o {output.bam} {input}; samtools index {output.bam}"

rule umitools_dedup:
    input:
        RESULTS_DIR + "/hisat2/{sample}.sorted.bam"
    output:
        RESULTS_DIR + "/umi_dedup/{sample}.sorted.umi.bam"
    conda:
        "envs/umitools_env.yaml"
    log:
        "logs/umitools_dedup/{sample}.log"
    shell:
        "umi_tools dedup --stdin={input} --log={log} --method=unique > {output}"

rule samtools_sort_index_dedupped:
    input:
        RESULTS_DIR + "/umi_dedup/{sample}.sorted.umi.bam"
    output:
        bam=RESULTS_DIR + "/umi_dedup_sorted/{sample}.sorted.umi.sorted.bam",
        bai=RESULTS_DIR + "/umi_dedup_sorted/{sample}.sorted.umi.sorted.bam.bai"
    conda:
        "envs/python_2.7.16.yaml"
    threads: 1
    shell:
        "samtools sort -m 10G -o {output.bam} {input}; samtools index {output.bam}"

rule pileup:
    input:
        bam=RESULTS_DIR + "/umi_dedup_sorted/{sample}.sorted.umi.sorted.bam",
        ref=REFERENCE_FASTA
    output:
        RESULTS_DIR + "/pileups_tmp/{sample}.pileups.tmp"
    conda:
        "envs/python_2.7.16.yaml"
    threads: 1
    log: "logs/pileup/{sample}.log"
    shell:
        "python GLORI_pipeline/scripts/pileup_genome.py -P {threads} -i {input.bam} -o {output} -f {input.ref} >& {log}"

rule format_output:
    input:
        pileup=RESULTS_DIR + "/pileups_tmp/{sample}.pileups.tmp",
        db=DB_PATH + "/db"
    output:
        pileup=RESULTS_DIR + "/pileups_txt/{sample}.pileups.txt",
        pileup_CR=RESULTS_DIR + "/pileups_CR/{sample}.pileups.CR"
    log:
        "logs/format_output/{sample}.log"
    conda:
        "envs/python_2.7.16.yaml"
    threads: 1
    shell:
        "python GLORI_pipeline/scripts/format_pileups.py -i {input.pileup} -o {output.pileup} --CR {output.pileup_CR} --db {input.db} >& {log}"

rule pileups:
    input:
        expand(RESULTS_DIR + "/pileups_txt/{sample}.pileups.txt",sample=SAMPLES)

rule sample_sheet:
    input:
        expand(RESULTS_DIR + "/pileups_txt/{sample}.pileups.txt",sample=SAMPLES)
    output:
        sheet=RESULTS_DIR + "/sample_sheet.tsv"
    run:
        with open(output.sheet, "w") as out:
            for sample in SAMPLES:
                out.write(sample + "\t" + RESULTS_DIR + "/pileups_txt/" + sample + ".pileups.txt" + "\t" + "gene" + "\t" + "3" + "\n")

rule m6A_caller:
    input:
        sheet=RESULTS_DIR + "/sample_sheet.tsv"
    output:
        m6A_calls=RESULTS_DIR + "/m6A_calls.csv"
    log:
        "logs/m6A_calls.txt"
    conda:
        "envs/m6A_caller.yaml"
    threads: 1
    shell:
        "python scripts/m6A_caller_2.py -i {input.sheet} -o {output.m6A_calls} -P {threads} 2>&1 | tee {log}"

rule evaluate_calls:
    input:
        m6A_calls=RESULTS_DIR + "/m6A_calls.csv"
    output:
        m6A_calls=RESULTS_DIR + "/m6A_calls_w_reliability.csv"
    log:
        "logs/evaluate_calls.txt"
    conda:
        "envs/m6A_caller.yaml"
    threads: 1
    shell:
        "python GLORI_pipeline/scripts/evaluate_sites.py -i {input} -o {output} 2>&1 | tee {log}"

rule gtf2anno:
    input:
        REFERENCE_GTF
    output:
        DB_PATH + "/anno"
    conda:
        "envs/python_2.7.16.yaml"
    threads: 1
    shell:
        "python scripts/gtf2anno_plus_gencode.py -i {input} > {output}"

rule anno_to_base:
    input:
        DB_PATH + "/anno"
    output:
        base=DB_PATH + "/base",
        base_sorted=DB_PATH + "/base.sorted",
    conda:
        "envs/python_2.7.16.yaml"
    threads: 10
    shell:
        "python RNA-m5C/0_m5C_step-by-step_metadata/anno_to_base.py -i {input} -o {output.base}"

rule anno_to_base_remove_redundance:
    input:
        base_sorted=DB_PATH + "/base.sorted",
        anno=DB_PATH + "/anno"
    output:
        DB_PATH + "/db"
    conda:
        "envs/python_2.7.16.yaml"
    threads: 1
    shell:
        "python scripts/anno_to_base_remove_redundance_v1.1.py -i {input.base_sorted} -a {input.anno} -o {output}"
