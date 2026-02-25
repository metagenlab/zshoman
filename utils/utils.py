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
    pipeline_indir=False,
    pipeline_outdir=False,
    pipeline_workdir=False,
    postprocessed_dir=False,
    analysis_dir=False,
    db_dir=False,
    per_sample=False,
    dry_run=False,
    threads=False,
    others=None,
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

    if pipeline_indir:
        parser.add_argument(
            "-o",
            "--pipeline_indir",
            default="input",
            type=Path,
            help="path to the input directory of the pipeline where the raw reads are stored.",
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
            default="output/post_processed",
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

    if db_dir:
        parser.add_argument(
            "db_dir",
            default="db",
            type=Path,
            help="path to the database directory",
        )

    if dry_run:
        parser.add_argument(
            "-n",
            "--dry_run",
            action="store_true",
            help="Only list files that would get deleted.",
        )

    if threads:
        parser.add_argument(
            "-t",
            "--threads",
            type=int,
            default=8,
            help="Number of parallel processes to use.",
        )

    if per_sample:
        parser.add_argument(
            "--per_sample",
            default=False,
            action="store_true",
            help="Whether the pipeline was run per sample or with the gene catalog."
            "Defaults to false (i.e. gene catalog)",
        )

    if others:
        for other in others:
            parser.add_argument(*other["args"], **other["kwargs"])

    return parser


def parse_arguments(
    samples_file="mandatory",
    pipeline_indir=False,
    pipeline_outdir=False,
    pipeline_workdir=False,
    postprocessed_dir=False,
    analysis_dir=False,
    db_dir=False,
    per_sample=False,
    dry_run=False,
    threads=False,
    others=None,
):
    parser = get_argument_parser(
        samples_file=samples_file,
        pipeline_indir=pipeline_indir,
        pipeline_outdir=pipeline_outdir,
        pipeline_workdir=pipeline_workdir,
        postprocessed_dir=postprocessed_dir,
        analysis_dir=analysis_dir,
        db_dir=db_dir,
        per_sample=per_sample,
        dry_run=dry_run,
        threads=threads,
        others=others,
    )
    args = parser.parse_args()

    if pipeline_indir and not pipeline_indir.exists():
        args.pipeline_indir.mkdir()

    if analysis_dir and not args.analysis_dir.exists():
        args.analysis_dir.mkdir()

    if postprocessed_dir and not args.postprocessed_dir.exists():
        args.postprocessed_dir.mkdir()

    if db_dir and not args.db_dir.exists():
        args.db_dir.mkdir()

    if samples_file:
        args.samples = SamplesGetter(args, with_files=(samples_file == "mandatory"))()

    return args


class SamplesGetter:
    to_exclude = ["gene_catalog", "pipeline_info"]

    def __init__(self, args, with_files=False):
        self.samples_file = args.samples_file
        self.pipeline_outdir = args.pipeline_outdir
        self.with_files = with_files

    def __call__(self):
        if self.with_files and not self.samples_file:
            raise RuntimeError(
                "Cannot request samples with files when determining samples from pipeline output directory."
            )
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

    @staticmethod
    def is_paired_end(row):
        return bool(row.get("fastq_R2", "").strip())

    def from_samples_file(self):
        samples_file = self.samples_file.resolve(strict=True)
        samplesheet = pd.read_csv(samples_file, header=0)
        samplesheet["sample"] = samplesheet["sample"].astype(str)

        if not self.with_files:
            return samplesheet["sample"].unique()

        grouped_by_sample = samplesheet.groupby("sample")
        samples = {}
        for sample_name, sample_data in grouped_by_sample:
            fastq1 = []
            fastq2 = []
            paired_end = self.is_paired_end(sample_data.iloc[0])
            for i, row in sample_data.iterrows():
                if not paired_end == self.is_paired_end(row):
                    raise ValueError(f"{sample_name} has paired and single end reads")
                fastq1.append(Path(row["fastq_R1"].strip()))
                if paired_end:
                    fastq2.append(Path(row["fastq_R2"].strip()))
            samples[sample_name] = {
                "paired_end": paired_end,
                "fastq1": fastq1,
                "fastq2": fastq2,
            }
        return samples

    def from_pipeline_outdir(self):
        return {
            el.name
            for el in self.pipeline_outdir.glob("*")
            if el.name not in self.to_exclude and el.is_dir()
        }
