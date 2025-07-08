"""
This script is used to copy necessary input files from a source folder to the
pipeline input folder. It will only copy over input files for samples that have
not been analysed yet. It also tries to avoid copying over files that are currently
being downloaded in the source folder.
"""

import argparse
import logging
import shutil
import time
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("copy-files")


class SamplesCopier:
    def __init__(self, samples_file, output_dir, download_dir, ignore_preprocessing):
        self.sample_file = samples_file
        self.output_dir = output_dir
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

        # Now among the samples that have not been analysed yet,
        # we will copy over the input if available and necessary
        to_copy = []
        for sample in to_keep:
            for file in data.loc[sample]:
                file = file.strip()
                if not file:
                    continue
                file = Path(file)
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
    args = argparse.ArgumentParser("copy_downloaded_files.py")
    args.add_argument(
        "samples_file",
        type=Path,
        help="path to the input file containing the list of samples",
    )
    args.add_argument(
        "download_dir",
        type=Path,
        help="path to the directory where the files were downloaded",
    )
    args.add_argument(
        "-o",
        "--output_dir",
        default="output",
        type=Path,
        help="path to the output directory of the pipeline",
    )
    args.add_argument(
        "--ignore_preprocessing",
        action="store_true",
        help="filter out sample even if preprocessed_reads folder is missing",
    )

    args = args.parse_args()
    SamplesCopier(
        args.samples_file,
        args.output_dir,
        args.download_dir,
        args.ignore_preprocessing,
    )()
