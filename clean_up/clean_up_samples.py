"""
This script will remove everything related to certain samples from
the work directory
"""

import argparse
import json
import logging
import shutil
from pathlib import Path
from pprint import pprint as pp

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("clean-up")


class WorkDirCleaner():

    def __init__(self, work_dir, samples_file, dry_run):
        self.work_dir = Path(work_dir)
        self.samples_file = Path(samples_file)
        self.dry_run = dry_run

    def __call__(self):
        samples = json.load(self.samples_file.open())
        samples = set(samples)
        to_delete = []
        for sample in samples:
            process_dirs = {el.parent for el in self.work_dir.glob(f"*/*/{sample}.*")}
            to_delete.extend(process_dirs)
        assert len(to_delete) == len(set(to_delete))

        n_folders = len(to_delete)
        logger.info(f"Deleting {n_folders} directories")
        print_files = input("print 20 first folders? [y]/n")
        if print_files == "y":
            pp(to_delete[:20])
        if self.dry_run:
            return

        input("ctl-c to cancel")
        for i, folder in enumerate(to_delete, 1):
            shutil.rmtree(folder)
            if i % 100 == 0:
                logger.info(f"Done {i}/{n_folders}")


if __name__ == '__main__':
    args = argparse.ArgumentParser("clean_up_preprocessing.py")
    args.add_argument(
        "samples_file",
        help="path to a json file containing a list of samples")
    args.add_argument(
        "-w", "--work_dir", default="work",
        help="path to nextflow's work directory")
    args.add_argument(
        "-n", "--dry_run", action='store_true',
        help="Only list files that would get deleted.")

    args = args.parse_args()
    WorkDirCleaner(args.work_dir, args.samples_file, args.dry_run)()
