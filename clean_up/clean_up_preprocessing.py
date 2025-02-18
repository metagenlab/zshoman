"""
This script will remove the preprocessed metagenomes from the output
directory for samples for which the pipeline is finished.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("remove-inputs")


class OutputRemover:
    def __init__(self, samples_file, output_dir, dry_run):
        self.sample_file = Path(samples_file)
        self.output_dir = Path(output_dir)
        self.dry_run = dry_run

    def __call__(self):
        data = pd.read_csv(args.samples_file, header=0)
        data.set_index("sample", inplace=True)
        samples = data.index
        to_keep = []
        to_delete = {}
        to_delete_nfiles = 0
        to_delete_size = 0
        for sample in samples:
            keep = False
            sample_path = Path(self.output_dir, sample)
            if not Path(sample_path, "preprocessed_reads").exists():
                continue
            for subdir in ("annotations", "assembly", "gene_counts", "motus", "phanta"):
                if not Path(sample_path, subdir).exists():
                    keep = True
                    to_keep.append(sample)
                    logger.info(f"{sample}: missing {subdir}")
                    break
            if keep:
                to_keep.append(sample)
                continue
            files = list(sample_path.glob(Path("preprocessed_reads", "*").as_posix()))
            to_delete_nfiles += len(files)
            to_delete_size += sum([file.stat().st_size for file in files])
            if files:
                to_delete[sample] = files
        logger.info(f"Keeping {len(to_keep)} samples")
        logger.warning(
            f"Deleting {to_delete_nfiles} files from {len(to_delete)} samples"
        )
        logger.info(f"This will free up {int(to_delete_size / 10**9)}GB.")

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
        "-n",
        "--dry_run",
        action="store_true",
        help="Only list files that would get deleted.",
    )

    args = args.parse_args()
    OutputRemover(args.samples_file, args.output_dir, args.dry_run)()
