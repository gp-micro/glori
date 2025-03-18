# GLORI Snakemake workflow

## Cloning this repository

This repository is little more than a [Snakemake](https://snakemake.readthedocs.io/en/stable/) wrapper for the following repositories (as well as some related scripts):

- [`jhfoxliu/GLORI_pipeline`](https://github.com/jhfoxliu/GLORI_pipeline)
- [`SYSU-zhanglab/RNA-m5C`](https://github.com/SYSU-zhanglab/RNA-m5C)

Those repositories are included as submodules. Here is how to clone *this* repository and include the above dependencies:

    git clone --recurse-submodules https://github.com/gp-micro/glori

## Configuration

Configuration of the pipeline takes place through the file `config.yaml`, which looks like this:

    results_dir: "results"
    sample_to_fastq:
        "test1": "GLORI_pipeline/notebook/test1.fastq.gz"
        "test2": "GLORI_pipeline/notebook/test2.fastq.gz"

    hisat2_path: "*MYSCRATCH*/GLORI/hisat2"
    hisat2_index_dir: "*MYSCRATCH*/GLORI_HISAT2_INDEXES"
    reference_fasta: "*PATH_TO_REF_FASTA*"
    reference_gtf: "*PATH_TO_REF_GTF*"
    db_path: "*MYSCRATCH*/db_will_be_created_here"

This Snakemake workflow currently *does not* perform the steps needed to generate the HISAT2 indexes. But I plan to add them.

## Software requirements

- conda
- Snakemake (I have been using version 8.16.0)

## Running

    CACHE=*MYSCRATCH*/snakemake_conda_cache
    snakemake --cores 10 --software-deployment-method conda --conda-prefix $CACHE --conda-frontend conda pileups


