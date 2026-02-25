"""
This script will filter out samples for which all the outputs
are present from an input file and optionally for which the input
file has not been downloaded yet.
"""

import sys
import time
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import logger
from utils.utils import parse_arguments


class SamplesCleaner:
    def __init__(
        self, samples_file, pipeline_outdir, ignore_preprocessing, only_downloaded
    ):
        self.sample_file = Path(samples_file)
        self.pipeline_outdir = Path(pipeline_outdir)
        self.only_downloaded = only_downloaded
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
            if not Path(self.pipeline_outdir, sample).exists():
                to_keep.append(sample)
                continue
            for subdir in self.expected_subdirs:
                if not Path(self.pipeline_outdir, sample, subdir).exists():
                    to_keep.append(sample)
                    logger.info(f"{sample}: missing {subdir}")
                    break

        if self.only_downloaded:
            for sample, files in data.iterrows():
                if sample not in to_keep:
                    continue
                remove = False
                for file in files:
                    file = file.strip()
                    if not file:
                        continue
                    file = Path(file)
                    if not file.exists():
                        remove = True
                        break
                    if (time.time() - file.stat().st_mtime) < 600:
                        remove = True
                        break
                if remove:
                    print(f"Removing {sample} as it is not finished downloading")
                    to_keep.remove(sample)

        logger.info(f"Keeping {len(to_keep)} samples")
        outname = Path(
            self.sample_file.parent,
            self.sample_file.stem + "_filtered" + self.sample_file.suffix,
        )
        if outname.exists():
            logger.warning(f"{outname} exists already, interrupting")
            sys.exit()
        data[data.index.isin(to_keep)].to_csv(outname, index=True)


if __name__ == "__main__":
    others = [
        {
            "args": ["--ignore_preprocessing"],
            "kwargs": {
                "action": "store_true",
                "help": "filter out sample even if preprocessed_reads folder is missing",
            },
        },
        {
            "args": ["--only_downloaded"],
            "kwargs": {
                "action": "store_true",
                "help": "filter out sample if input files have not been downloaded yet",
            },
        },
    ]

    args = parse_arguments(
        samples_file="mandatory",
        pipeline_outdir=True,
        others=others,
    )

    SamplesCleaner(
        args.samples_file,
        args.pipeline_outdir,
        args.ignore_preprocessing,
        args.only_downloaded,
    )()
