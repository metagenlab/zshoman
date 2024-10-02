# Taxonomic assignment and functional annotation pipeline for metagenomics data

This pipeline is meant to analyse whole-genome metagenomic data.
For preprocessing the pipeline follows guidelines from [methods in microbiomics](https://methods-in-microbiomics.readthedocs.io/en/latest/preprocessing/preprocessing.html) as well as for [assembly](https://methods-in-microbiomics.readthedocs.io/en/latest/assembly/metagenomic_workflows.html#id1), the [gene catalog](https://methods-in-microbiomics.readthedocs.io/en/latest/assembly/metagenomic_workflows.html#gene-catalogs) and [taxonomic profiling with mOTUs](https://methods-in-microbiomics.readthedocs.io/en/latest/taxonomic_profiling/metagenomes.html).

The pipeline is meant to handle paired-end and single-end sequencing and analyse not only bacteria, but also viruses and eukaryotes (fungi).

## Running the pipeline on obelix

To run the pipeline on obelix:

- start a screen
- get an allocation for running nextflow (`salloc --cpus-per-task=16 --partition=short`)
- activate nextflow environment `conda activate nextflow`
- run the pipeline `nextflow run ../CIDB/scripts/motus_pipeline/motus_pipeline.nf --input samples_test.csv --phanta_db /mnt/slow_storage/databases/phanta/unmasked_db_v1/ -profile apptainer -e.cpus=60 -process.executor=slurm -process.queue=short --max_time="12.h" --db_dir /mnt/slow_storage/users/njohner/CIDB/Analysis/DB/motus_pipeline/ -resume --phanta_db=/mnt/slow_storage/databases/phanta/unmasked_db_v1/`
