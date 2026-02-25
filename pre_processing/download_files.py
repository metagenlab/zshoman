"""
While the pipeline handles genomes hosted on online resources
(i.e. the files will get downloaded when executing the pipeline),
this can lead to issues as download sometimes fails and the pipeline
will need to be restarted. This script allows to download the files
in advance, and rewrites the input file to point to the downloaded files.
"""

import os
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import logger
from utils.utils import parse_arguments


class FileDownloader:
    def __init__(self, samples, input_dir, output_dir):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.files = self.get_files(samples)
        logger.info(f"Found {len(self.files)} files to download.")
        input("ctl-c to cancel")

    def get_files(self, samples):
        files = []
        filenames = set()
        for sample_name, sample_data in samples.items():
            if Path(self.output_dir, sample_name, "preprocessed_reads").exists():
                logger.info(f"skipping {sample_name}")
                continue
            sample_files = sample_data["fastq1"] + sample_data["fastq2"]
            for file in sample_files:
                filename = file.name
                if Path(self.input_dir, filename).is_file():
                    logger.info(f"skipping {sample_name}, {filename}")
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


if __name__ == "__main__":
    others = [
        {
            "args": ["--postfix"],
            "kwargs": {
                "default": "downloaded",
                "help": "postfix for output samplesheet filename",
            },
        },
        {
            "args": ["--skip_download"],
            "kwargs": {
                "action": "store_true",
                "help": "Do not actually download the files, only rewrite the samplesheet.",
            },
        },
    ]

    args = parse_arguments(
        samples_file="mandatory",
        pipeline_indir=True,
        pipeline_outdir=True,
        threads=True,
        others=others,
    )

    output_samples_file = Path(
        args.samples_file.parent,
        f"{args.samples_file.stem}_{args.postfix}{args.samples_file.suffix}",
    )

    with open(args.samples_file) as infile_handle:
        with open(output_samples_file, "w") as outfile_handle:
            outfile_handle.write(next(infile_handle))
            for line in infile_handle:
                res = [el.strip() for el in line.split(",")]
                res[1:] = [
                    os.path.join(args.input_dir, el.rsplit("/", 1)[-1]) if el else ""
                    for el in res[1:]
                ]
                outfile_handle.write(", ".join(res) + "\n")

    if not args.skip_download:
        FileDownloader(args.samples, args.pipeline_indir, args.pipeline_outdir)(
            args.threads
        )
