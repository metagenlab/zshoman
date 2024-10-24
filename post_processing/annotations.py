"""
This script will process the output from eggnog and the gene abundances
and generate tables containing the annotation abundances.
"""

import argparse
import glob
import logging
import os

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Post-processing")


def annotation_abunances(input_dir, output_dir, old_style):
    logger.info("Loading eggnog annotations")
    eggnog_file = os.path.join(input_dir, "gene_catalog", "all.emapper.annotations")
    data = pd.read_csv(eggnog_file, header=4, delimiter="\t", index_col=0)
    data.where(data != "-", np.nan, inplace=True)

    if old_style:
        count_dir = "gene_counts"
    else:
        count_dir = "gene_counts_gc"

    abundance_files = glob.glob(os.path.join(input_dir, "*", count_dir, "*_genes_per_cell.csv"))

    samples = []
    for abundance_file in abundance_files:
        sample = abundance_file.rsplit("/", 1)[-1].rstrip("_genes_per_cell.csv")
        samples.append(sample)
        logger.info(f"Merging annotations for {sample}")
        abundance = pd.read_csv(abundance_file, delimiter=",", index_col=0, names=["gene", sample])
        data = data.join(abundance)

    for colname in ['COG_category', 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway',
                    'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE',
                    'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs']:
        logger.info(f"Creating {colname} abundances table")
        annotations = data[colname][pd.notna(data[colname])].str.split(",").explode()
        annotations = annotations.to_frame().merge(
            data[samples], how="left", left_on="#query", right_index=True)
        annotations.groupby(colname).sum().to_csv(os.path.join(output_dir, f"{colname}.csv"))


if __name__ == '__main__':
    args = argparse.ArgumentParser("annotations.py")
    args.add_argument(
        "input_dir",
        help="path to the output directory of the pipeline")
    args.add_argument(
        "-o", "--output_dir",
        help="path where the post-processed tables will be written to. "
             "defaults to post_processed subdirectory in the input_dir")
    args.add_argument(
        "--old_style", action="store_true",
        help="For old versions of the pipeline where abundances for the gene "
             "catalog are stored in 'gene_catalog' instead of 'gene_catalog_gc'")

    args = args.parse_args()
    output_dir = args.output_dir or os.path.join(args.input_dir, "post_processed")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    annotation_abunances(args.input_dir, output_dir, args.old_style)
