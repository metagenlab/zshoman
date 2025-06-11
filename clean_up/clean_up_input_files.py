"""
This script will remove the downloaded input files for samples for which
the pipeline is finished.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("remove-inputs")


class InputRemover:
    def __init__(self, samples_file, output_dir, ignore_preprocessing, dry_run):
        self.sample_file = Path(samples_file)
        self.output_dir = Path(output_dir)
        self.dry_run = dry_run
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
        to_delete = {}
        for sample in samples:
            for subdir in self.expected_subdirs:
                if not Path(self.output_dir, sample, subdir).exists():
                    to_keep.append(sample)
                    logger.info(f"{sample}: missing {subdir}")
                    break
            if sample not in to_keep:
                files = []
                for file in data.loc[sample]:
                    if pd.isna(file):
                        continue
                    file = file.strip()
                    if not file:
                        continue
                    file = Path(file.strip())
                    if file.is_file():
                        files.append(file)
                if files:
                    to_delete[sample] = files
        logger.info(f"Keeping {len(to_keep)} samples")
        logger.warning(f"Deleting files for {len(to_delete)} samples")

        print_files = input("print paths of 10 first samples? [y]/n")
        if print_files == "y":
            for sample, files in list(to_delete.items())[:10]:
                logger.info(f"{sample}: {files}")

        if self.dry_run:
            return

        input("ctl-c to cancel")
        for i, (sample, files) in enumerate(to_delete.items(), 1):
            for file in files:
                file.unlink()
            if i % 10 == 0:
                logger.info(f"Done {i}/{len(to_delete)}")


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
    args.add_argument(
        "-n",
        "--dry_run",
        action="store_true",
        help="Only list files that would get deleted.",
    )

    args = args.parse_args()
    InputRemover(
        args.samples_file, args.output_dir, args.ignore_preprocessing, args.dry_run
    )()
