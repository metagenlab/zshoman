import argparse
import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("zshoman utils")


def get_argument_parser(
    samples_file="mandatory",
    pipeline_outdir=False,
    pipeline_workdir=False,
    postprocessed_dir=False,
    analysis_dir=False,
    download_dir=False,
    per_sample=False,
    dry_run=True,
):
    parser = argparse.ArgumentParser()

    if samples_file == "mandatory":
        parser.add_argument("samples_file", help="path to samples csv file.", type=Path)
    elif samples_file == "optional":
        parser.add_argument(
            "-i",
            "--samples_file",
            type=Path,
            help="path to samples csv file.",
        )

    if pipeline_outdir:
        parser.add_argument(
            "-o",
            "--pipeline_outdir",
            default="output",
            type=Path,
            help="path to the output directory of the pipeline",
        )

    if pipeline_workdir:
        parser.add_argument(
            "-w",
            "--pipeline_workdir",
            default="work",
            type=Path,
            help="path to the work directory of the pipeline",
        )

    if postprocessed_dir:
        parser.add_argument(
            "--postprocessed_dir",
            default="post_processed",
            type=Path,
            help="path to the directory where post-processed outputs will be stored",
        )

    if analysis_dir:
        parser.add_argument(
            "--analysis_dir",
            default="analysis",
            type=Path,
            help="path to the work directory of the pipeline",
        )

    if download_dir:
        parser.add_argument(
            "download_dir",
            default="input",
            type=Path,
            help="path to the location where the files should get downloaded to",
        )

    if dry_run:
        parser.add_argument(
            "-n",
            "--dry_run",
            action="store_true",
            help="Only list files that would get deleted.",
        )

    if per_sample:
        parser.add_argument(
            "--per_sample",
            default=False,
            action="store_true",
            help="Whether the pipeline was run per sample or with the gene catalog."
            "Defaults to false (i.e. gene catalog)",
        )
    return parser


def parse_arguments(
    samples_file="mandatory",
    pipeline_outdir=False,
    pipeline_workdir=False,
    postprocessed_dir=False,
    analysis_dir=False,
    download_dir=False,
    per_sample=False,
    dry_run=True,
):
    parser = get_argument_parser(
        samples_file=samples_file,
        pipeline_outdir=pipeline_outdir,
        pipeline_workdir=pipeline_workdir,
        postprocessed_dir=postprocessed_dir,
        analysis_dir=analysis_dir,
        download_dir=download_dir,
        per_sample=per_sample,
        dry_run=dry_run,
    )
    args = parser.parse_args()

    if analysis_dir:
        if not analysis_dir.exists():
            analysis_dir.mkdir()

    if postprocessed_dir:
        if not postprocessed_dir.exists():
            postprocessed_dir.mkdir()

    return args


class SamplesGetter:
    to_exclude = ["gene_catalog", "pipeline_info"]

    def __init__(self, args):
        self.samples_file = args.samples_file
        self.pipeline_outdir = args.pipeline_outdir

    def __call__(self):
        if self.samples_file:
            samples = self.from_samples_file()
        elif self.pipeline_outdir:
            samples = self.from_pipeline_outdir()
        else:
            raise RuntimeError(
                "Need to specify either samples_file or pipeline_outdir to determine samples"
            )
        logger.info(f"Found {len(samples)} samples.")
        return samples

    def from_samples_file(self):
        samples_file = self.samples_file.resolve(strict=True)
        return pd.read_csv(samples_file, header=0)["sample"].unique()

    def from_pipeline_outdir(self):
        return {
            el.name
            for el in self.pipeline_outdir.glob("*")
            if el.name not in self.to_exclude and el.is_dir()
        }
