"""
This script will process the output from eggnog and the gene abundances
and generate tables containing the annotation abundances.
"""

import argparse
import logging
import os
from functools import reduce
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Post-processing")


class Sample():

    def __init__(self, name, input_dir):
        self.name = name
        self.input_dir = input_dir

    @property
    def gene_file(self):
        return Path(self.input_dir, self.name, "gene_counts",
                    f"{self.name}_genes_per_cell.csv")

    @property
    def annotation_file(self):
        return Path(self.input_dir, self.name, "annotations",
                    f"{self.name}.emapper.annotations")


class AnnotationAbundanceCalculator():

    cols_to_transform = [
        'COG_category', 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway',
        'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE',
        'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs']

    def __init__(self, samples_file, input_dir, output_dir, per_sample):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.per_sample = per_sample
        self.sample_names = self.read_samples_file(samples_file)["sample"]

        self.samples = [
            Sample(sname, self.input_dir) for sname in self.sample_names]

    def __call__(self):
        if not self.per_sample:
            annotations = self.load_annotations()
            abundances = reduce(
                lambda left, right: pd.join(left, right),
                (self.load_abundances(sample) for sample in self.samples))

            for colname in self.cols_to_transform:
                self.get_annotation_abundances(annotations, abundances, colname).to_csv(
                    os.path.join(self.output_dir, f"{colname}.csv"))
        else:
            for sample in self.samples:
                logger.info(f"Loading data for {sample.name}")
                annotations = self.load_annotations(sample)
                abundances = self.load_abundances(sample)
                for colname in self.cols_to_transform:
                    setattr(sample, colname, self.get_annotation_abundances(
                        annotations, abundances, colname))

            for colname in self.cols_to_transform:
                annotation_abundances = reduce(
                    lambda left, right: pd.merge(left, right, left_index=True,
                                                 right_index=True, how="outer"),
                    (getattr(sample, colname) for sample in self.samples))

                annotation_abundances.to_csv(
                    os.path.join(self.output_dir, f"{colname}.csv"))

    @staticmethod
    def read_samples_file(samples_file):
        data = pd.read_csv(args.samples_file, header=0)
        return data

    def get_annotation_abundances(self, annotations, abundances, colname):
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
        col = annotations[colname]
        annotations = col[pd.notna(col)].str.split(",").explode()
        annotations = annotations.to_frame().merge(
            abundances, how="left", left_on="#query", right_index=True)
        return annotations.groupby(colname).sum()

    def annotation_file(self, sample):
        if sample is None:
            return os.path.join(self.input_dir, "gene_catalog",
                                "all.emapper.annotations")
        else:
            return sample.annotation_file

    def load_annotations(self, sample=None):
        data = pd.read_csv(self.annotation_file(sample), header=4,
                           delimiter="\t", index_col=0)
        data.where(data != "-", np.nan, inplace=True)
        return data

    @staticmethod
    def load_abundances(sample):
        return pd.read_csv(sample.gene_file, delimiter=",", index_col=0,
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
    parser.add_argument(
        "--per_sample", default=False, action="store_true",
        help="Whether the pipeline was run per sample or with the gene catalog.")
    args = parser.parse_args()
    output_dir = args.output_dir or os.path.join(args.input_dir, "post_processed")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    AnnotationAbundanceCalculator(args.samples_file, args.input_dir,
                                  output_dir, args.per_sample)()
