# Taxonomic assignment and functional annotation pipeline for metagenomics data

This pipeline is meant to analyse whole-genome metagenomic data.
For preprocessing the pipeline follows guidelines from [methods in microbiomics](https://methods-in-microbiomics.readthedocs.io/en/latest/preprocessing/preprocessing.html) as well as for [assembly](https://methods-in-microbiomics.readthedocs.io/en/latest/assembly/metagenomic_workflows.html#id1), the [gene catalog](https://methods-in-microbiomics.readthedocs.io/en/latest/assembly/metagenomic_workflows.html#gene-catalogs) and [taxonomic profiling with mOTUs](https://methods-in-microbiomics.readthedocs.io/en/latest/taxonomic_profiling/metagenomes.html).

The pipeline is meant to handle paired-end and single-end sequencing and analyse not only bacteria, but also viruses and eukaryotes (fungi).

## Running the pipeline on obelix

To run the pipeline on obelix:

- start a screen
- get an allocation for running nextflow (`salloc --cpus-per-task=1 --partition=long`)
- activate nextflow environment `conda activate nextflow`
- run the pipeline `nextflow run ../scripts/CIDB/motus_pipeline/motus_pipeline.nf --input samples_test.csv --phanta_db /mnt/slow_storage/databases/phanta/unmasked_db_v1/ -profile apptainer -e.cpus=60 -process.executor=slurm -process.queue=short --max_time="4.h" --db_dir /mnt/slow_storage/users/njohner/CIDB/Analysis/DB/motus_pipeline/ -resume --motus_db=/mnt/slow_storage/databases/db_mOTU/ --skip_phanta=false --skip_dev=false --microeuk_db /mnt/slow_storage/databases/MicroEuk_v3/MicroEuk90.faa.gz`

## Creating the MicroEuk90 database

For Eukaryotic gene calling we use MetaEuk, which requires a database to match genes against. We propose to use the MicroEuk database from [veba](https://github.com/jolespin/veba). When downloading the database you get a database with all sequences and some files containing clusters for the database clustered at 90% and 50% sequence identity. We propose to use the 90% sequence identity clustering. For this we need to generate the database from the `MicroEuk100.faa.gz` and `MicroEuk90_clusters.tsv.gz` files. This can be easily done by getting the set of unique identifiers from the first column of `MicroEuk90_clusters.tsv.gz`, writing that to a file `MicroEuk90_representatives.csv` and reusing that to filter the database with `seqtk subseq MicroEuk100.faa.gz MicroEuk90_representatives.csv  | gzip --no-name > MicroEuk90.faa.gz`.
