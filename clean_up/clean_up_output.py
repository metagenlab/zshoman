"""
This script is used to remove directories from the output folder
which have missing files. It can remove assembly directories
for samples for which not all necessary files are present in the
assembly folder. This is to remove them for samples which were run
before we stored all necessary files to skip re-running the
assembly in the pipeline.
It can also remove preprocessing folders if some files are missing
there.
"""

import argparse
import logging
import shutil
from pathlib import Path
from pprint import pprint as pp

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("clean-up")


class OutputDirCleaner():

    def __init__(self, output_dir, dry_run, assembly, preprocessing):
        self.output_dir = Path(output_dir)
        self.dry_run = dry_run
        self.assembly = assembly
        self.preprocessing = preprocessing

    def __call__(self):
        if self.assembly:
            self.clean_assemblies()
        if self.preprocessing:
            self.clean_preprocessing()

    def get_dirs_to_delete(self, dirname, necesseray_files):
        outdirs = self.output_dir.glob(f"*/{dirname}")
        to_delete = set()
        for outdir in outdirs:
            sample = outdir.parent.name
            for ending in necesseray_files:
                if not Path(outdir, sample + ending).exists():
                    to_delete.add(outdir)
        return to_delete

    def delete_dirs(self, to_delete, label):
        n_folder = len(to_delete)
        logger.info(f"Will remove {n_folder} {label} directories")

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

    def clean_assemblies(self):
        to_delete = self.get_dirs_to_delete(
            "assembly",
            [".scaffolds.fa.gz", ".assembly.gfa.gz", ".scaffolds.paths.gz"])
        self.delete_dirs(to_delete, "assembly")

    def clean_preprocessing(self):
        paired_to_delete = self.get_dirs_to_delete(
            "preprocessed_reads",
            ["_1_unmerged.fastq.gz",
             "_2_unmerged.fastq.gz",
             "_merged.fastq.gz",
             "_host_filtered_1.fastq.gz",
             "_host_filtered_2.fastq.gz",
             "_host_filtered_singletons.fastq.gz"])
        single_to_delete = self.get_dirs_to_delete(
            "preprocessed_reads",
            ["_host_filtered.fastq.gz"])
        self.delete_dirs(paired_to_delete.intersection(single_to_delete),
                         "preprocessed_reads")


if __name__ == '__main__':
    args = argparse.ArgumentParser("clean_up_output.py")
    args.add_argument(
        "-o", "--output_dir", default="output",
        help="path to the output directory of the pipeline")
    args.add_argument(
        "-n", "--dry_run", action='store_true',
        help="Only list files that would get deleted.")
    args.add_argument(
        "--assembly", action='store_true',
        help="Clean-up assembly files")
    args.add_argument(
        "--preprocessing", action='store_true',
        help="Clean-up pre-processing files")

    args = args.parse_args()
    OutputDirCleaner(args.output_dir, args.dry_run, args.assembly, args.preprocessing)()
