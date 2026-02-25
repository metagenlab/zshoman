"""
This script is used to copy necessary input files from a source folder to the
pipeline input folder. It will only copy over input files for samples that have
not been analysed yet. It also tries to avoid copying over files that are currently
being downloaded in the source folder.
"""

import shutil
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import logger
from utils.utils import parse_arguments


class SamplesCopier:
    def __init__(self, samples, pipeline_outdir, download_dir, ignore_preprocessing):
        self.samples = samples
        self.pipeline_outdir = pipeline_outdir
        self.download_dir = download_dir
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
        to_keep = []
        for sample in self.samples:
            if not Path(self.pipeline_outdir, sample).exists():
                to_keep.append(sample)
                continue
            for subdir in self.expected_subdirs:
                if not Path(self.pipeline_outdir, sample, subdir).exists():
                    to_keep.append(sample)
                    logger.info(f"{sample}: missing {subdir}")
                    break

        # Now among the samples that have not been analysed yet,
        # we will copy over the input if available and necessary
        to_copy = []
        for sample_name in to_keep:
            sample_data = self.samples[sample_name]
            files = sample_data["fastq1"] + sample_data["fastq2"]
            for file in files:
                if file.exists():
                    continue
                # If the file exists in the download folder and is finished downloading
                # (i.e. its modification time is not too recent), we will copy it over
                download_path = Path(self.download_dir, file.name)
                if (
                    download_path.exists()
                    and (time.time() - download_path.stat().st_mtime) > 600
                ):
                    to_copy.append((download_path, file))

        logger.info(f"Found {len(to_copy)} files to copy.")
        input("ctl-c to cancel")
        for src, dest in to_copy:
            logger.info(f"Copying {src} to {dest}")
            shutil.copy2(src, dest)


if __name__ == "__main__":
    others = [
        {
            "args": ["--src"],
            "kwargs": {
                "type": Path,
                "help": "source directory from which to copy the files",
            },
        },
        {
            "args": ["--ignore_preprocessing"],
            "kwargs": {
                "action": "store_true",
                "help": "filter out sample even if preprocessed_reads folder is missing",
            },
        },
    ]
    args = parse_arguments(
        samples_file="mandatory",
        pipeline_outdir=True,
        others=others,
    )

    args = args.parse_args()
    SamplesCopier(
        args.samples,
        args.pipeline_outdir,
        args.src,
        args.ignore_preprocessing,
    )()
