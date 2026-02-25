"""
This script will remove the downloaded input files for samples for which
the pipeline is finished.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import logger
from utils.utils import parse_arguments


class InputRemover:
    def __init__(self, samples_file, pipeline_outdir, ignore_preprocessing, dry_run):
        self.sample_file = Path(samples_file)
        self.pipeline_outdir = Path(pipeline_outdir)
        self.dry_run = dry_run
        self.expected_subdirs = [
            ".",
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
        to_delete = {}
        for sample in samples:
            files = []
            for file in data.loc[sample]:
                if pd.isna(file):
                    continue
                file = file.strip()
                if not file:
                    continue
                file = Path(file.strip())
                if file.is_file():
                    files.append(file)
            if not files:
                continue
            for subdir in self.expected_subdirs:
                if not Path(self.pipeline_outdir, sample, subdir).exists():
                    to_keep.append(sample)
                    logger.info(f"{sample}: missing {subdir}")
                    break
            if sample not in to_keep:
                to_delete[sample] = files
        logger.info(f"Keeping {len(to_keep)} samples")
        logger.warning(f"Deleting files for {len(to_delete)} samples")

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
    others = [
        {
            "args": ["--ignore_preprocessing"],
            "kwargs": {
                "action": "store_true",
                "help": "remove input even if preprocessed_reads folder is missing",
            },
        },
    ]

    args = parse_arguments(
        samples_file="mandatory",
        pipeline_outdir=True,
        dry_run=True,
        others=others,
    )

    InputRemover(
        args.samples_file, args.pipeline_outdir, args.ignore_preprocessing, args.dry_run
    )()
