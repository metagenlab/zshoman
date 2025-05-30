"""
This script will remove the files generated by assembly in
nextflow's work directory, provided that we do not need them anymore
(i.e. we have the files in the output directory).
"""

import argparse
import logging
from collections import defaultdict
from pathlib import Path
from pprint import pprint as pp

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("clean-up")


class WorkDirCleaner():

    def __init__(self, work_dir, output_dir, dry_run):
        self.work_dir = Path(work_dir)
        self.output_dir = Path(output_dir)
        self.dry_run = dry_run

    def log_to_sample(self, log_file):
        filename = log_file.name
        endings = [".spades.log"]
        for end in endings:
            if filename.endswith(end):
                return filename.rstrip(end)

    def __call__(self):
        sample_process_dirs = defaultdict(list)
        spades_process_logs = self.work_dir.glob("**/*.spades.log")
        for log_file in spades_process_logs:
            sample_process_dirs[self.log_to_sample(log_file)].append(log_file.parent)
        to_handle = []
        for sample, process_dirs in sample_process_dirs.items():
            if not Path(self.output_dir,  sample, "assembly").exists():
                logger.info(f"skipping {sample}")
                continue
            to_handle.extend(process_dirs)
        to_delete = []
        for process_dir in to_handle:
            to_delete.extend([file for file in process_dir.glob("*.gz")
                              if not file.is_symlink()])
            to_delete.extend([file for file in process_dir.glob("*.fast*")
                              if not file.is_symlink()])
            to_delete.extend([file for file in process_dir.glob("*.gfa")
                              if not file.is_symlink()])

        tot_size = sum([file.stat().st_size for file in to_delete])
        n_files = len(to_delete)
        logger.info(f"Deleting {n_files} files "
                    f"from {len(to_handle)} directories.")
        logger.info(f"This will free up {int(tot_size/10**9)}GB.")
        print_files = input("print 20 first files? [y]/n")
        if print_files == "y":
            pp(to_delete[:20])
        if self.dry_run:
            return

        input("ctl-c to cancel")
        for i, file in enumerate(to_delete, 1):
            file.unlink()
            if i % 100 == 0:
                logger.info(f"Done {i}/{n_files}")


if __name__ == '__main__':
    args = argparse.ArgumentParser("clean_up_preprocessing.py")
    args.add_argument(
        "-w", "--work_dir", default="work",
        help="path to nextflow's work directory")
    args.add_argument(
        "-o", "--output_dir", default="output",
        help="path to the output directory of the pipeline")
    args.add_argument(
        "-n", "--dry_run", action='store_true',
        help="Only list files that would get deleted.")

    args = args.parse_args()
    WorkDirCleaner(args.work_dir, args.output_dir, args.dry_run)()
