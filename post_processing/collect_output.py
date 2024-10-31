"""
This script will gather output files from the nextflow output directory,
copy (and rename them if necessary) to a different location.
"""

import argparse
import logging
import os
import shutil
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Post-processing")


class FileCollector():

    to_exclude = ["gene_catalog", "pipeline_info"]

    def __init__(self, input_dir, output_dir):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.samples = {el.name for el in self.input_dir.glob("*")
                        if el.name not in self.to_exclude}
        logger.info(f"Found {len(self.samples)} samples.")

    def __call__(self, phanta=False):
        if phanta:
            for sample in self.samples:
                phanta_dir = Path(self.input_dir, sample, "phanta", "final_merged_outputs")
                if not phanta_dir.is_dir():
                    logger.warning(f"Skipping phanta for {sample}")
                    continue
                for fname in ["counts.txt", "relative_read_abundance.txt", "relative_taxonomic_abundance.txt"]:
                    shutil.copy(Path(phanta_dir, fname),
                                Path(self.output_dir, f"{sample}_{fname}"))


if __name__ == '__main__':
    args = argparse.ArgumentParser("annotations.py")
    args.add_argument(
        "input_dir",
        help="path to the output directory of the pipeline")
    args.add_argument(
        "output_dir",
        help="path to the location where the files should get copied to")
    args.add_argument(
        "--phanta", action="store_true",
        help="copy the phanta output'")

    args = args.parse_args()
    output_dir = args.output_dir or os.path.join(args.input_dir, "post_processed")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    FileCollector(args.input_dir, args.output_dir)(phanta=args.phanta)
