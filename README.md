# Ze SHOtgun Metagenomics ANalysis pipeline

This pipeline is meant to analyse whole-genome metagenomic data.
For preprocessing the pipeline follows guidelines from [methods in microbiomics](https://methods-in-microbiomics.readthedocs.io/en/latest/preprocessing/preprocessing.html) as well as for [assembly](https://methods-in-microbiomics.readthedocs.io/en/latest/assembly/metagenomic_workflows.html#id1), the [gene catalog](https://methods-in-microbiomics.readthedocs.io/en/latest/assembly/metagenomic_workflows.html#gene-catalogs) and [taxonomic profiling with mOTUs](https://methods-in-microbiomics.readthedocs.io/en/latest/taxonomic_profiling/metagenomes.html).

The pipeline is meant to handle paired-end and single-end sequencing and analyse not only bacteria, but also viruses and eukaryotes (fungi).

An overview of the pipeline is shown below:

![ShoMAn] workflow](https://github.com/metagenlab/zshoman/main/assets/zdb_workflow.png)


## Running the pipeline on obelix

To run the pipeline on obelix:

- start a screen
- get an allocation for running nextflow (`salloc --cpus-per-task=1 --partition=long`)
- activate nextflow environment `conda activate nextflow`
- run the pipeline `nextflow run /mnt/slow_storage/users/njohner/scripts/CIDB/motus_pipeline/motus_pipeline.nf --input samples_test_paired.csv --db_dir  /mnt/slow_storage/databases/ -resume -c /mnt/slow_storage/metagenlab/configs/conf/metagenlab.config`

## Creating the MicroEuk90 database

For Eukaryotic gene calling we use MetaEuk, which requires a database to match genes against. We propose to use the MicroEuk database from [veba](https://github.com/jolespin/veba). When downloading the database you get a database with all sequences and some files containing clusters for the database clustered at 90% and 50% sequence identity. We propose to use the 90% sequence identity clustering. For this we need to generate the database from the `MicroEuk100.faa.gz` and `MicroEuk90_clusters.tsv.gz` files. This can be easily done by getting the set of unique identifiers from the first column of `MicroEuk90_clusters.tsv.gz`, writing that to a file `MicroEuk90_representatives.csv` and reusing that to filter the database with `seqtk subseq MicroEuk100.faa.gz MicroEuk90_representatives.csv  | gzip --no-name > MicroEuk90.faa.gz`.

## Post-processing

Post-processing scripts are python scripts run using the conda environment defined in [post_processing.yaml](https://github.com/metagenlab/zshoman/post_processing/post_processing.yaml).
To run the scripts, simply create the conda environment, activate it and run the script with python:
```
conda env create -f post_processing/post_processing.yaml
conda activate post_processing
python post_processing/annotations.py path/to/pipeline/output
```

### Eggnog annotations

The Eggnog annotation table is a bit complex to analyse as is, as it condenses many types of annotations. Notably for annotation types for which a given gene can have several annotations, the cells will contain coma-separated lists of annotations, e.g. `ko:K00336,ko:K01101`. To simplify analysis we provide a post-processing script ([annotations.py](https://github.com/metagenlab/zshoman/post_processing/annotations.py)) which will output a table for each annotation type, containing the annotation (e.g. `ko:K00336`) and its abundance in each sample.
