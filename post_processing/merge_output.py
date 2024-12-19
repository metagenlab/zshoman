"""
This script will gather the outputs (for now mOTUs) for all samples from the
nextflow output directory, and merge them into a single table.
"""

import argparse
import logging
import os
from functools import reduce
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Post-processing")


class TableMerger():

    to_exclude = ["gene_catalog", "pipeline_info", "logs"]

    def __init__(self, samples_file, input_dir, output_dir, outname):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.outname = outname
        self.samples = pd.read_csv(args.samples_file, header=0)["sample"]
        logger.info(f"Found {len(self.samples)} samples.")

    def load_motus_table(self, sample):
        return pd.read_csv(
            Path(self.input_dir, sample, "motus", sample + ".motus"),
            sep="\t", header=2)

    def __call__(self, cleanup=True):
        motus_table = reduce(
            lambda left, right: pd.merge(left, right),
            (self.load_motus_table(sample) for sample in self.samples))

        if cleanup:
            motus_total = motus_table.iloc[:, 3:].sum(axis=1)
            non_zero = motus_total != 0
            motus_table = motus_table[non_zero]

        outpath = Path(self.output_dir, self.outname + ".motus")
        motus_table.to_csv(outpath, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser("annotations.py")
    parser.add_argument("samples_file", help="path to samples csv file.")
    parser.add_argument(
        "--input_dir", default="output",
        help="path to the output directory of the pipeline")
    parser.add_argument(
        "--output_dir",
        help="path where the merged tables should be written to")
    parser.add_argument(
        "--outname", default="merged",
        help="prefix for output files")
    parser.add_argument(
        "--no_cleanup", action="store_true",
        help="do not remove empty rows")

    args = parser.parse_args()
    output_dir = args.output_dir or os.path.join(args.input_dir, "post_processed")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    TableMerger(args.samples_file, args.input_dir, output_dir, args.outname)(
        not args.no_cleanup)
