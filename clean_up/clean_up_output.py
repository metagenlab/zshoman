"""
This script will remove assembly directories from the output folder
for samples for which not all necessary files are present in the
assembly folder. This is to remove them for samples which were run
before we stored all necessary files to skip re-running the
assembly in the pipeline.
"""

import argparse
import logging
import shutil
from pathlib import Path
from pprint import pprint as pp

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("clean-up")


class OutputDirCleaner():

    def __init__(self, output_dir, dry_run):
        self.output_dir = Path(output_dir)
        self.dry_run = dry_run

    def __call__(self):
        assembly_dirs = self.output_dir.glob("*/assembly")
        to_delete = set()
        for assembly_dir in assembly_dirs:
            sample = assembly_dir.parent.name
            for ending in [".scaffolds.fa.gz", ".assembly.gfa.gz", ".scaffolds.paths.gz"]:
                if not Path(assembly_dir, sample + ending).exists():
                    to_delete.add(assembly_dir)

        n_folder = len(to_delete)
        logger.info(f"Will remove {n_folder} assembly directories")

        print_files = input("print 20 first paths? [y]/n")
        if print_files == "y":
            pp(list(to_delete)[:20])

        if self.dry_run:
            return

        input("ctl-c to cancel")
        for i, folder in enumerate(to_delete, 1):
            shutil.rmtree(folder)
            if i % 10 == 0:
                logger.info(f"Done {i}/{n_folder}")


if __name__ == '__main__':
    args = argparse.ArgumentParser("clean_up_preprocessing.py")
    args.add_argument(
        "-o", "--output_dir", default="output",
        help="path to the output directory of the pipeline")
    args.add_argument(
        "-n", "--dry_run", action='store_true',
        help="Only list files that would get deleted.")

    args = args.parse_args()
    OutputDirCleaner(args.output_dir, args.dry_run)()
