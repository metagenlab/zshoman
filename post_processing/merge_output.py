"""
This script will gather the outputs (for now mOTUs or phanta) for all samples from the
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


class TableMerger:
    to_exclude = ["gene_catalog", "pipeline_info", "logs"]
    out_ext = "tsv"

    def __init__(self, samples_file, input_dir, output_dir, prefix):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.out_prefix = prefix
        self.samples = pd.read_csv(args.samples_file, header=0)["sample"]
        logger.info(f"Found {len(self.samples)} samples.")

    def load_table(self, sample):
        raise NotImplementedError()

    @property
    def outpath(self):
        return Path(
            self.output_dir, f"{self.out_prefix}_{self.out_name}.{self.out_ext}"
        )

    def filter_samples(self):
        samples = [
            sample for sample in self.samples if self.get_table_path(sample).exists()
        ]
        logger.info(f"Keeping {len(samples)} / {len(self.samples)} samples.")
        return samples

    def __call__(self, cleanup=True):
        samples = self.filter_samples()
        merged_table = reduce(
            lambda left, right: pd.merge(left, right, how="outer"),
            (self.load_table(sample) for sample in samples),
        )

        if cleanup:
            merged_total = merged_table.loc[:, samples].sum(axis=1)
            non_zero = merged_total != 0
            merged_table = merged_table[non_zero]

        merged_table.to_csv(self.outpath, index=False)


class MotusMerger(TableMerger):
    out_name = "motus"

    def get_table_path(self, sample):
        return Path(self.input_dir, sample, "motus", sample + ".motus")

    def load_table(self, sample):
        return pd.read_csv(self.get_table_path(sample), sep="\t", header=2)


class PhantaMerger(TableMerger):
    def __call__(self, table_name, cleanup=True):
        self.table_name = table_name
        super(PhantaMerger, self).__call__(cleanup)
        self.table_name = None

    @property
    def out_name(self):
        return f"phanta_{self.table_name}"

    def get_table_path(self, sample):
        return Path(
            self.input_dir,
            sample,
            "phanta",
            "final_merged_outputs",
            f"{self.table_name}.txt",
        )

    def load_table(self, sample):
        return pd.read_csv(
            self.get_table_path(sample),
            sep="\t",
            header=0,
        ).rename(columns={f"{sample}_": sample})


if __name__ == "__main__":
    parser = argparse.ArgumentParser("annotations.py")
    parser.add_argument("samples_file", help="path to samples csv file.")
    parser.add_argument(
        "--input_dir",
        default="output",
        help="path to the output directory of the pipeline",
    )
    parser.add_argument(
        "--output_dir", help="path where the merged tables should be written to"
    )
    parser.add_argument("--prefix", default="merged", help="prefix for output files")
    parser.add_argument(
        "--no_cleanup", action="store_true", help="do not remove empty rows"
    )
    parser.add_argument("--motus", action="store_true", help="merge motus table")
    parser.add_argument("--phanta", action="store_true", help="merge phanta table")

    args = parser.parse_args()
    output_dir = args.output_dir or os.path.join(args.input_dir, "post_processed")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if args.motus:
        MotusMerger(args.samples_file, args.input_dir, output_dir, args.prefix)(
            not args.no_cleanup
        )

    if args.phanta:
        merger = PhantaMerger(
            args.samples_file, args.input_dir, output_dir, args.prefix
        )
        merger("relative_taxonomic_abundance", not args.no_cleanup)
        merger("relative_read_abundance", not args.no_cleanup)
        merger("counts", not args.no_cleanup)
