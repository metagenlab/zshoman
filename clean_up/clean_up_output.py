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

import shutil
import sys
from pathlib import Path
from pprint import pprint as pp

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import logger
from utils.utils import parse_arguments


class OutputDirCleaner:
    def __init__(self, pipeline_outdir, dry_run, assembly, preprocessing):
        self.pipeline_outdir = Path(pipeline_outdir)
        self.dry_run = dry_run
        self.assembly = assembly
        self.preprocessing = preprocessing

    def __call__(self):
        if self.assembly:
            self.clean_assemblies()
        if self.preprocessing:
            self.clean_preprocessing()

    def get_dirs_to_delete(self, dirname, necesseray_files):
        outdirs = self.pipeline_outdir.glob(f"*/{dirname}")
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
            "assembly", [".scaffolds.fa.gz", ".assembly.gfa.gz", ".scaffolds.paths.gz"]
        )
        self.delete_dirs(to_delete, "assembly")

    def clean_preprocessing(self):
        paired_to_delete = self.get_dirs_to_delete(
            "preprocessed_reads",
            [
                "_1_unmerged.fastq.gz",
                "_2_unmerged.fastq.gz",
                "_merged.fastq.gz",
                "_host_filtered_1.fastq.gz",
                "_host_filtered_2.fastq.gz",
                "_host_filtered_singletons.fastq.gz",
            ],
        )
        single_to_delete = self.get_dirs_to_delete(
            "preprocessed_reads", ["_host_filtered.fastq.gz"]
        )
        self.delete_dirs(
            paired_to_delete.intersection(single_to_delete), "preprocessed_reads"
        )


if __name__ == "__main__":
    others = [
        {
            "args": ["--assembly"],
            "kwargs": {
                "action": "store_true",
                "help": "Clean-up assembly files",
            },
        },
        {
            "args": ["--preprocessing"],
            "kwargs": {
                "action": "store_true",
                "help": "Clean-up pre-processing files",
            },
        },
    ]

    args = parse_arguments(
        samples_file=False,
        pipeline_outdir=True,
        dry_run=True,
        others=others,
    )

    OutputDirCleaner(
        args.pipeline_outdir, args.dry_run, args.assembly, args.preprocessing
    )()
