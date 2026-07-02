# Ze SHOtgun Metagenomics ANalysis pipeline

This pipeline is meant to analyse whole-genome metagenomic data.
For preprocessing the pipeline follows guidelines from [methods in microbiomics](https://methods-in-microbiomics.readthedocs.io/en/latest/preprocessing/preprocessing.html) as well as for [assembly](https://methods-in-microbiomics.readthedocs.io/en/latest/assembly/metagenomic_workflows.html#id1), the [gene catalog](https://methods-in-microbiomics.readthedocs.io/en/latest/assembly/metagenomic_workflows.html#gene-catalogs) and [taxonomic profiling with mOTUs](https://methods-in-microbiomics.readthedocs.io/en/latest/taxonomic_profiling/metagenomes.html).

The pipeline is meant to handle paired-end and single-end sequencing and analyse not only bacteria, but also viruses and eukaryotes (fungi).

An overview of the pipeline is shown below:

![ShoMAn workflow](https://github.com/metagenlab/zshoman/blob/main/assets/metagenomics_pipeline.drawio.svg)


## Running the pipeline

The pipeline needs nextflow 23.10.0 or later to run. You can for example install nextflow with `conda`.
You will also need to clone this repository, and you can then simply run the pipeline (`main.nf`) with nextflow.
Overall the following should work:

```

conda create -n nextflow nextflow -c bioconda
conda activate nextflow
git clone git@github.com:metagenlab/zshoman.git
cd zshoman
nextflow run main.nf --input path/to/samplesheet.csv --db_dir path/to/databases
```

To see the full list of parameters you can use `nextflow run main.nf --help`

### Input

Input is a csv file with 2 or 3 columns:

- sample: label of the samples
- fastq_1: path to forward reads
- fastq_2: optional. path to reverse reads

For single-end samples, the `fastq_2` can be omitted or left empty. The first row should contain the column headers (`sample`, `fastq_R1`, `fastq_R2`). See the [input template](https://github.com/metagenlab/zshoman/blob/main/assets/input_template.csv) or  [multilane input template](https://github.com/metagenlab/zshoman/blob/main/assets/input_template_multilane.csv) for an example.

### Databases

The pipeline requires various reference databases which are expected to be organised in subdirectories of a single database directory specified with `db_dir`. By default it expects the following organisation:

- `db_dir/phanta/uhggv2_mgv/` which can be downloaded as specified [in the phanta repository](https://github.com/bhattlab/phanta/blob/main/databases.md)
- `db_dir/db_mOTU/` for the motus database which can be downloaded using [the mOTUs software](https://github.com/motu-tool/mOTUs)
- `db_dir/MicroEuk_v3/MicroEuk90.faa.gz` for eukaryotic gene calling, which can be obtained [as explained below](#creating-the-microeuk90-database)
- `db_dir/eggnog/` for the eggnog database, which can be downloaded from the [eggnog website](http://eggnog5.embl.de/#/app/home)
- `db_dir/metagenlab_wgs_pipeline/phix174_ill.ref.fa.gz` for contaminant removal, as obtained from the [methods in microbiomics website](https://methods-in-microbiomics.readthedocs.io/en/latest/_downloads/3d3ab63ee68ea60cc2374b0690387094/Sample1_isolate.tar.gz)
- `db_dir/metagenlab_wgs_pipeline/adapters.fa` for adapter trimming, as obtained from the [methods in microbiomics website](https://methods-in-microbiomics.readthedocs.io/en/latest/_downloads/3d3ab63ee68ea60cc2374b0690387094/Sample1_isolate.tar.gz)
- `db_dir/metagenlab_wgs_pipeline/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz` for host filtering (masked human genome which can be obtained [here](https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?resourcekey=0-PsIKmg2q4EvTGWGOUjsKGQ))

If need be, you can also override the path to each database separately (see [nextflow.config](https://github.com/metagenlab/zshoman/blob/main/nextflow.config) for more details).

### Main options

The pipeline has two main branches, either making a gene catalog and doing functional annotation for that catalog, or making the functional annotation independently for each sample. The former makes comparison of genes between samples easier, but requires rerunning a large portion of the pipeline if a new sample needs to be added. These branches can be skipped with `--skip_gene_catalog` and `--skip_per_sample`. Other parts of the pipeline can also be skipped (see [nextflow.config](https://github.com/metagenlab/zshoman/blob/main/nextflow.config) for more details).


### Managing disk space usage and resuming the pipeline from output files

Disk space usage when running the pipeline can quickly become problematic because of accumulation of files in the `work` directory. One of the nextflow core developers has developed an experimental plug-in allowing to delete files during the pipeline execution (https://github.com/bentsherman/nf-boost). There are currently two issues with that plug-in:

- It breaks the `-resume` option, meaning that after a crash, the whole pipeline has to be run again.
- It is buggy and sometimes removes files that are still required for the run, leading to random crashes in the pipeline. There is already [an issue for this](https://github.com/bentsherman/nf-boost/issues/4)

We nevertheless use this plugin by default to limit disk space usage with several consequences:
1. To minimize the issue with the broken `-resume`, we have included in the pipeline ways to skip entire blocks of the pipeline if they were already run (final output files stored in the `output` directory) and allow restarting other blocks of the pipeline by loading the necessary files from the `output` folder. This is done automatically when using the `--resume_from_output` flag, allowing to skip the most time consuming steps of the pipeline when resuming.

2. To avoid having the whole pipeline stop when a single process fails due to the plugin we have set the error strategy to ignore. When running the gene catalog we do not actually want to run the catalog when any task has failed before. Current strategy (will implement a check in the workflow at some point) is to first run the pipeline with `--skip_gene_catalog`, ensure that all samples have passed, and then run again without `--skip_gene_catalog` and of course with `--resume_from_output`.


## Running the pipeline on obelix

To run the pipeline on obelix:

- start a `screen` or `tmux`, depending on your religion
- activate nextflow environment `conda activate nextflow`
- run the pipeline `nextflow run /mnt/slow_storage/metagenlab/zshoman_dev/main.nf --input samples_test.csv --db_dir  /mnt/fast_storage/databases/ -resume -c /mnt/slow_storage/metagenlab/configs/conf/metagenlab.config`


## Creating the MicroEuk90 database

For Eukaryotic gene calling we use MetaEuk, which requires a database to match genes against. We propose to use the MicroEuk database from [veba](https://github.com/jolespin/veba). When downloading the database (see https://github.com/jolespin/veba/tree/main/data/MicroEuk_v3) you get a database with all sequences and some files containing clusters for the database clustered at 90% and 50% sequence identity. We propose to use the 90% sequence identity clustering. For this we need to generate the database from the `MicroEuk100.faa.gz` and `MicroEuk90_clusters.tsv.gz` files. This can be easily done by getting the set of unique identifiers from the first column of `MicroEuk90_clusters.tsv.gz`, writing that to a file `MicroEuk90_representatives.csv` and reusing that to filter the database with `seqtk subseq MicroEuk100.faa.gz MicroEuk90_representatives.csv  | gzip --no-name > MicroEuk90.faa.gz`.

## Pre-processing, Post-processing and clean-up scripts

Pre-processing, post-processing and clean-up scripts are python scripts run using the conda environment defined in [post_processing.yaml](https://github.com/metagenlab/zshoman/post_processing/post_processing.yaml).
To run the scripts, simply create the conda environment, activate it and run the script with python:
```
conda env create -f post_processing/post_processing.yaml
conda activate post_processing
python post_processing/annotations.py path/to/pipeline/output
```

You can get help for a given script by passing `-h` then calling the command, e.g.
```
python post_processing/annotations.py -h
```

Here a quick overview of available scripts:
- pre-processing:
  - `download_files.py`: while the pipeline handles genomes hosted on online resources (i.e. the files will get downloaded when executing the pipeline), this can lead to issues as download sometimes fails and the pipeline will need to be restarted. This script allows to download the files in advance, and rewrites the input file to point to the downloaded files.
  - `filter_samples.py`: This script will remove samples for which the analysis is complete from the input file. Can be useful if the pipeline needs to be resumed but nextflow's resume is not available (e.g. because the nf-boost plugin is used to delete intermediary data from the `work` directory)
- post-processing:
  - `annotations.py`: This script will process the output from eggnog and the gene abundances and generate tables containing the annotation abundances.
  - `download_kegg_db.py`: downloads the most recent KEGG module definitions
  - `correct_module_abundances.py`: The current eggnog database uses an older version of the KEGG database, which includes modules that do not exist anymore. This script will recalculate the module abundance table from the KO abundance table produced by the `annotations.py` script using the newest module definitions downloaded using the `download_kegg_db.py` script.
  - `calculate_ko_module_completeness.py`: This script will calculate the abundance of complete KEGG modules from the corrected module abundance and kegg DB obtained with the two previous scripts.
  - `check_quality.py`: This script will extract information from the log files and prepare summary tables and plots, notably to check the quality of the data and the run.
  - `collect_output.py`: This script will gather output files from the nextflow output directory, copy (and rename them if necessary) to a different location.
  - `merge_output.py`: This script will gather the outputs (for now mOTUs) for all samples from the nextflow output directory, and merge them into a single table.
- clean-up:
  - `clean_up_assembly.py`: This script will remove the files generated by assembly in nextflow's work directory, provided that we do not need them anymore (i.e. we have the files in the output directory).
  - `clean_up_output.py`: This script will remove assembly directories from the output folder for samples for which not all necessary files are present in the assembly folder. This is to remove them for samples which were run before we stored all necessary files to skip re-running the assembly in the pipeline.
  - `clean_up_preprocessing.py`: This script will remove the files generated by the phix filtering step of preprocessing from nextflow's work directory, provided that we do not need them anymore (i.e. we have the files in the output directory).
  - `clean_up_samples.py`: This script will remove everything related to certain samples from the work directory


### Eggnog annotations

The Eggnog annotation table is a bit complex to analyse as is, as it condenses many types of annotations. Notably for annotation types for which a given gene can have several annotations, the cells will contain coma-separated lists of annotations, e.g. `ko:K00336,ko:K01101`. To simplify analysis we provide a post-processing script ([annotations.py](https://github.com/metagenlab/zshoman/post_processing/annotations.py)) which will output a table for each annotation type, containing the annotation (e.g. `ko:K00336`) and its abundance in each sample.


## Standard output analysis steps

For the purpose of this quick guide we will call `zshoman_dir` the path where zshoman is found and assume you are currently working inside the folder where the pipeline was run and that you used the defaults for the naming of output folders and such.

### Output organisation

Output of the pipeline is found in the `output` folder. It will contain one sub-folder for each sample, containing itself subfolders with the output from various steps of the pipeline (pre-processing, assembly, taxonomic profiling, gene calling, gene profiling...). It also contains a `pipeline_info` folder containing reports automatically generated by nextflow, a `logs` folder containing the log files from all processes run, a `gene_catalog` folder if you ran the pipeline with the gene catalog and a `post_processed` folder once you processed the outputs with some of the post-processing scripts as described below.

### Quality control

The pipeline collects all logs from the different steps for every sample in `output/logs`. You can extract the main information from these logs by running the `check_quality.py` post-processing script:

```
python zshoman_dir/post_processing/check_quality.py samples.csv 
```

This will create a few plots in the `analysis` folder showing how many reads/bases were removed during pre-processing, as well as write out a `statistics.csv` file containing number of reads used and output by the different processes in the pipeline.

### Taxonomic profiles

To collect taxonomic profiles into a single table for all your samples you can use the `merge_output.py` post-processing script with flags corresponding to the taxonomic profiles you generated with the pipeline (mOTUs or phanta)

```
python zshoman_dir/post_processing/merge_output.py samples.csv --motus --phanta
```

This will collect the profiles from individual samples and create a single taxonomic profiles table in the `output/post_processed` folder.

### Functional profiles

To get the functional profiles we need to combine the information from the gene profiles with the functional annotations of the genes. The gene profiles are stored in the samples subfolders (`output/sample_id/gene_counts` or `output/sample_id/gene_counts_gc` if you ran the gene catalog), while the functional annotations are either in the samples subfolders (`output/sample_id/annotations`) or in the `output/gene_catalog/annotations`. This can be done using the `annotations.py` post-processing script:

```
python zshoman_dir/post_processing/annotations.py samples.csv
```

As EggNOG uses an outdated database, some of the KEGG modules used therein do no longer exist (see release 92.0 of [KEGG](https://www.genome.jp/kegg/docs/relnote.html) and [issue in eggnog](https://github.com/eggnogdb/eggnog-mapper/issues/545)). To use up to date module definitions we need to download their definitions KEGG with the `download_kegg_db.py` script and then recompute the module abundances from these updated definitions with the `correct_module_abundances.py` script.

```
python zshoman_dir/post_processing/download_kegg_db.py
python zshoman_dir/post_processing/correct_module_abundances.py
```

### Gene profiles

If you have run the gene catalog, you can collect the profiles into a single table (as for the taxonomic profiles). Simply run

```
python zshoman_dir/post_processing/merge_output.py samples.csv --genes
```
