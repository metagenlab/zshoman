"""
This script will gather output files from the nextflow output directory,
copy (and rename them if necessary) to a different location.
"""

import argparse
import logging
import os
import subprocess
from multiprocessing import Pool
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Post-processing")


class FileDownloader():

    def __init__(self, samples_file, input_dir, output_dir):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.files = self.read_samples_file(samples_file)

        logger.info(f"Found {len(self.files)} files to download.")
        input("ctl-c to cancel")

    def read_samples_file(self, samples_file):
        files = []
        filenames = set()
        with open(samples_file) as file_handle:
            # Skip title row
            next(file_handle)
            for line in file_handle:
                res = line.split(",")
                res = [el.strip() for el in res]
                sample = res[0]
                if os.path.exists(os.path.join(self.output_dir, sample, "preprocessed_reads")):
                    logger.info(f"skipping {sample}")
                    continue
                for file in res[1:]:
                    filename = file.rsplit("/", 1)[-1]
                    if os.path.isfile(os.path.join(self.input_dir, filename)):
                        logger.info(f"skipping {sample}, {filename}")
                        continue
                    filenames.add(filename)
                    files.append(file)

        # Make sure all files are different
        assert len(filenames) == len(files)
        return files

    def download_file(self, file):
        subprocess.call(["wget", file, "-P", self.input_dir])

    def __call__(self, n):
        with Pool(n) as p:
            p.map(self.download_file, self.files)


if __name__ == '__main__':
    args = argparse.ArgumentParser("download_files.py")
    args.add_argument(
        "samples_file",
        help="path to the input file containing the list of samples")
    args.add_argument(
        "input_dir",
        help="path to the location where the files should get downloaded to")
    args.add_argument(
        "-o", "--output_dir",
        help="path to the output directory of the pipeline")
    args.add_argument(
        "-n",
        help="number of parallel processes. Default is 10.")
    args.add_argument(
        "-f", "--output_samples_file",
        help="ouput sample file")

    args = args.parse_args()
    input_dir = args.input_dir
    if not os.path.exists(input_dir):
        os.mkdir(input_dir)

    n = int(args.n or 10)
    FileDownloader(args.samples_file, args.input_dir, args.output_dir)(n)

    if args.output_samples_file:
        with open(args.samples_file) as infile_handle:
            with open(args.output_samples_file, "w") as outfile_handle:
                outfile_handle.write(next(infile_handle))
                for line in infile_handle:
                    res = [el.strip() for el in line.split(",")]
                    outfile_handle.write(", ".join(
                        [res[0]] + [os.path.join(args.input_dir, el.rsplit("/", 1)[-1])
                                    for el in res[1:]]) + "\n")
