"""
This script will filter out samples for which all the outputs
are present from an input file.
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("filter-samples")


class SamplesCleaner:
    def __init__(self, samples_file, output_dir, ignore_preprocessing):
        self.sample_file = Path(samples_file)
        self.output_dir = Path(output_dir)
        self.expected_subdirs = [
            "annotations",
            "assembly",
            "gene_counts",
            "motus",
            "phanta",
            "preprocessed_reads",
        ]
        if ignore_preprocessing:
            self.expected_subdirs.remove("preprocessed_reads")

    def __call__(self):
        data = pd.read_csv(args.samples_file, header=0)
        data.set_index("sample", inplace=True)
        samples = data.index
        to_keep = []
        for sample in samples:
            if not Path(self.output_dir, sample).exists():
                to_keep.append(sample)
                continue
            for subdir in self.expected_subdirs:
                if not Path(self.output_dir, sample, subdir).exists():
                    to_keep.append(sample)
                    logger.info(f"{sample}: missing {subdir}")
                    break

        logger.info(f"Keeping {len(to_keep)} samples")
        outname = Path(
            self.sample_file.parent,
            self.sample_file.stem + "_filtered" + self.sample_file.suffix,
        )
        if outname.exists():
            logger.warning(f"{outname} exists already, interrupting")
            sys.exit()
        data[data.index.isin(to_keep)].to_csv(outname, index=True)


if __name__ == "__main__":
    args = argparse.ArgumentParser("clean_up_preprocessing.py")
    args.add_argument(
        "samples_file", help="path to the input file containing the list of samples"
    )
    args.add_argument(
        "-o",
        "--output_dir",
        default="output",
        help="path to the output directory of the pipeline",
    )
    args.add_argument(
        "--ignore_preprocessing",
        action="store_true",
        help="remove input even if preprocessed_reads folder is missing",
    )

    args = args.parse_args()
    SamplesCleaner(args.samples_file, args.output_dir, args.ignore_preprocessing)()
