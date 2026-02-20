"""
This script will gather output files from the nextflow output directory,
copy (and rename them if necessary) to a different location.
"""

import shutil
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import logger
from utils.utils import parse_arguments


class FileCollector:
    to_exclude = ["gene_catalog", "pipeline_info"]

    def __init__(self, input_dir, output_dir):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.samples = {
            el.name for el in self.input_dir.glob("*") if el.name not in self.to_exclude
        }
        logger.info(f"Found {len(self.samples)} samples.")

    def __call__(self, phanta=False):
        if phanta:
            for sample in self.samples:
                phanta_dir = Path(
                    self.input_dir, sample, "phanta", "final_merged_outputs"
                )
                if not phanta_dir.is_dir():
                    logger.warning(f"Skipping phanta for {sample}")
                    continue
                for fname in [
                    "counts.txt",
                    "relative_read_abundance.txt",
                    "relative_taxonomic_abundance.txt",
                ]:
                    shutil.copy(
                        Path(phanta_dir, fname),
                        Path(self.output_dir, f"{sample}_{fname}"),
                    )


if __name__ == "__main__":
    others = [
        {
            "args": ["--phanta"],
            "kwargs": {"action": "store_true", "help": "copy the phanta output"},
        }
    ]
    args = parse_arguments(
        samples_file="mandatory",
        pipeline_outdir=True,
        postprocessed_dir=True,
        others=others,
    )

    FileCollector(args.pipeline_outdir, args.postprocessed_dir)(phanta=args.phanta)
