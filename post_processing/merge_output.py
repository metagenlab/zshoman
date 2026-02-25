"""
This script will gather the outputs (for now mOTUs or phanta) for all samples from the
nextflow output directory, and merge them into a single table.
"""

import sys
from functools import reduce
from pathlib import Path

import pandas as pd

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import logger
from utils.utils import parse_arguments



class TableMerger:
    to_exclude = ["gene_catalog", "pipeline_info", "logs"]
    out_ext = "csv"

    def __init__(self, samples_file, pipeline_outdir, postprocessed_dir, prefix):
        self.pipeline_outdir = Path(pipeline_outdir)
        self.postprocessed_dir = Path(postprocessed_dir)
        self.out_prefix = prefix
        self.samples = pd.read_csv(args.samples_file, header=0)["sample"]
        logger.info(f"Found {len(self.samples)} samples.")

    def load_table(self, sample):
        raise NotImplementedError()

    @property
    def outpath(self):
        return Path(
            self.postprocessed_dir, f"{self.out_prefix}_{self.out_name}.{self.out_ext}"
        )

    def filter_samples(self):
        samples = [
            sample for sample in self.samples if self.get_table_path(sample).exists()
        ]
        logger.info(f"Keeping {len(samples)} / {len(self.samples)} samples.")
        return samples

    def __call__(self, cleanup=True):
        samples = self.filter_samples()
        merged_table = reduce(
            lambda left, right: pd.merge(left, right, how="outer"),
            (self.load_table(sample) for sample in samples),
        )

        if cleanup:
            merged_total = merged_table.loc[:, samples].sum(axis=1)
            non_zero = merged_total != 0
            merged_table = merged_table[non_zero]

        merged_table.to_csv(self.outpath, index=False)


class MotusMerger(TableMerger):
    out_name = "motus"

    def get_table_path(self, sample):
        return Path(self.pipeline_outdir, sample, "motus", sample + ".motus")

    def load_table(self, sample):
        return pd.read_csv(self.get_table_path(sample), sep="\t", header=2)


class PhantaMerger(TableMerger):
    def __call__(self, table_name, cleanup=True):
        self.table_name = table_name
        super(PhantaMerger, self).__call__(cleanup)
        self.table_name = None

    @property
    def out_name(self):
        return f"phanta_{self.table_name}"

    def get_table_path(self, sample):
        return Path(
            self.pipeline_outdir,
            sample,
            "phanta",
            "final_merged_outputs",
            f"{self.table_name}.txt",
        )

    def load_table(self, sample):
        return pd.read_csv(
            self.get_table_path(sample),
            sep="\t",
            header=0,
        ).rename(columns={f"{sample}_": sample})


if __name__ == "__main__":
    others = [
        {
            "args": ["--prefix"],
            "kwargs": {"default": "merged", "help": "prefix for output file names"},
        },
        {
            "args": ["--motus"],
            "kwargs": {"action": "store_true", "help": "merge motus table"},
        },
        {
            "args": ["--phanta"],
            "kwargs": {"action": "store_true", "help": "merge phanta table"},
        }
        {
            "args": ["--no_cleanup"],
            "kwargs": {"action": "store_true", "help": "do not remove empty rows from merged tables"},
        }
    ]
    args = parse_arguments(
        samples_file="mandatory",
        pipeline_outdir=True,
        postprocessed_dir=True,
        others=others
    )

    if args.motus:
        MotusMerger(args.samples_file, args.pipeline_outdir, args.postprocessed_dir, args.prefix)(
            not args.no_cleanup
        )

    if args.phanta:
        merger = PhantaMerger(
            args.samples_file, args.pipeline_outdir, args.postprocessed_dir, args.prefix
        )
        merger("relative_taxonomic_abundance", not args.no_cleanup)
        merger("relative_read_abundance", not args.no_cleanup)
        merger("counts", not args.no_cleanup)
