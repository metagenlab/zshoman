"""
This script will remove the preprocessed metagenomes from the output
directory for samples for which the pipeline is finished.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import logger
from utils.utils import parse_arguments


class OutputRemover:
    def __init__(self, samples, pipeline_outdir, per_sample, dry_run):
        self.samples = samples
        self.pipeline_outdir = Path(pipeline_outdir)
        self.per_sample = per_sample
        self.dry_run = dry_run

    def required_subdirs(self):
        required_subdirs = ["assembly", "motus", "phanta"]
        if self.per_sample:
            required_subdirs.extend(["annotations", "gene_counts"])
        else:
            required_subdirs.extend(["gene_counts_gc"])

    def __call__(self):
        to_keep = []
        to_delete = {}
        to_delete_nfiles = 0
        to_delete_size = 0
        for sample in self.samples:
            keep = False
            sample_path = Path(self.pipeline_outdir, sample)
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
    args = parse_arguments(
        samples_file="optional",
        pipeline_outdir=True,
        per_sample=True,
        dry_run=True,
    )

    OutputRemover(args.samples, args.pipeline_outdir, args.per_sample, args.dry_run)()
