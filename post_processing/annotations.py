"""
This script will process the output from eggnog and the gene abundances
and generate tables containing the annotation abundances.
"""

import argparse
import logging
import os
from collections import namedtuple
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Post-processing")

Sample = namedtuple("Sample", ["name", "file"])


class AnnotationAbundanceCalculator():

    cols_to_transform = [
        'COG_category', 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway',
        'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE',
        'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs']

    def __init__(self, samples_file, input_dir, output_dir, old_style):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.sample_names = self.read_samples_file(samples_file)["sample"]

        self.samples = [
            Sample(sname, Path(self.input_dir, sname, "gene_counts",
                               f"{sname}_genes_per_cell.csv"))
            for sname in self.sample_names]

    def __call__(self):
        data = self.load_annotations()
        for sample in self.samples:
            logger.info(f"Merging abundances for {sample}")
            data = data.join(self.load_abundances(sample))

        for colname in self.cols_to_transform:
            self.get_annotation_abundances(data, colname).to_csv(
                os.path.join(self.output_dir, f"{colname}.csv"))

    @staticmethod
    def read_samples_file(samples_file):
        data = pd.read_csv(args.samples_file, header=0)
        return data

    def get_annotation_abundances(self, data, colname):
        """
        There are several annotations in a cell, coma-separated,
        e.g. ko:K00336,ko:K01101. To get the counts we need to split
        these into separate rows. Note that doing this directly on the
        dataframe will not handle the other columns properly, hence we do it
        only on the annotation column and then merge with the abundances.
        Of course a given annotation can appear several times
        (i.e. for different genes) so we then need to group by annotation
        and sum the abundances.
        """
        logger.info(f"Creating {colname} abundances table")
        annotations = data[colname][pd.notna(data[colname])].str.split(",").explode()
        annotations = annotations.to_frame().merge(
            data[self.sample_names], how="left", left_on="#query", right_index=True)
        return annotations.groupby(colname).sum()

    @property
    def annotation_file(self):
        return os.path.join(self.input_dir, "gene_catalog",
                            "all.emapper.annotations")

    def load_annotations(self):
        logger.info("Loading eggnog annotations")
        data = pd.read_csv(self.annotation_file, header=4,
                           delimiter="\t", index_col=0)
        data.where(data != "-", np.nan, inplace=True)
        return data

    @staticmethod
    def load_abundances(sample):
        return pd.read_csv(sample.file, delimiter=",", index_col=0,
                           names=["gene", sample.name])


if __name__ == '__main__':
    parser = argparse.ArgumentParser("annotations.py")
    parser.add_argument("samples_file", help="path to samples csv file.")
    parser.add_argument(
        "-i", "--input_dir", default="output",
        help="path to the output directory of the pipeline")
    parser.add_argument(
        "-o", "--output_dir",
        help="path where the post-processed tables will be written to. "
             "defaults to post_processed subdirectory in the input_dir")

    args = parser.parse_args()
    output_dir = args.output_dir or os.path.join(args.input_dir, "post_processed")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    AnnotationAbundanceCalculator(args.input_dir, output_dir, args.old_style)()
