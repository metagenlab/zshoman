"""
This script will process the output from eggnog and the gene abundances
and generate tables containing the annotation abundances.
"""

import os
import sys
from functools import reduce
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import logger
from utils.utils import parse_arguments


class Sample:
    def __init__(self, name, pipeline_outdir):
        self.name = name
        self.pipeline_outdir = pipeline_outdir

    def has_all_files(self):
        return self.gene_file.exists() and self.annotation_file.exists()

    @property
    def gene_file(self):
        return Path(
            self.pipeline_outdir,
            self.name,
            "gene_counts",
            f"{self.name}_genes_per_cell.csv",
        )

    @property
    def annotation_file(self):
        return Path(
            self.pipeline_outdir,
            self.name,
            "annotations",
            f"{self.name}.emapper.annotations",
        )


class AnnotationAbundanceCalculator:
    cols_to_transform = [
        "COG_category",
        "GOs",
        "EC",
        "KEGG_ko",
        "KEGG_Pathway",
        "KEGG_Module",
        "KEGG_Reaction",
        "KEGG_rclass",
        "BRITE",
        "KEGG_TC",
        "CAZy",
        "BiGG_Reaction",
        "PFAMs",
    ]

    def __init__(self, samples, pipeline_outdir, output_dir, per_sample):
        self.pipeline_outdir = pipeline_outdir
        self.output_dir = output_dir
        self.per_sample = per_sample
        self.sample_names = samples

        self.samples = [
            Sample(sname, self.pipeline_outdir) for sname in self.sample_names
        ]
        if per_sample:
            to_skip = [el for el in self.samples if not el.has_all_files()]
            to_skip_names = [el.name for el in to_skip]
            logger.warning(f"Skipping {len(to_skip)} samples missing a file.")
            logger.warning(to_skip_names)
            self.samples = [el for el in self.samples if el not in to_skip]
            self.sample_names = [
                el for el in self.sample_names if el not in to_skip_names
            ]

    def __call__(self):
        if not self.per_sample:
            annotations = self.load_annotations()
            abundances = reduce(
                lambda left, right: pd.join(left, right),
                (self.load_abundances(sample) for sample in self.samples),
            )

            for colname in self.cols_to_transform:
                self.get_annotation_abundances(annotations, abundances, colname).to_csv(
                    os.path.join(self.output_dir, f"{colname}.csv")
                )
        else:
            n_samples = len(self.samples)
            for i, sample in enumerate(self.samples, 1):
                annotations = self.load_annotations(sample)
                abundances = self.load_abundances(sample)
                for colname in self.cols_to_transform:
                    setattr(
                        sample,
                        colname,
                        self.get_annotation_abundances(
                            annotations, abundances, colname
                        ),
                    )
                if i % 25 == 0:
                    logger.info(f"Loaded data for {i} / {n_samples}")

            for colname in self.cols_to_transform:
                annotation_abundances = reduce(
                    lambda left, right: pd.merge(
                        left, right, left_index=True, right_index=True, how="outer"
                    ),
                    (getattr(sample, colname) for sample in self.samples),
                )

                annotation_abundances.to_csv(
                    os.path.join(self.output_dir, f"{colname}.csv")
                )

    def get_annotation_abundances(self, annotations, abundances, colname):
        """
        There are several annotations in a cell, coma-separated,
        e.g. ko:K00336,ko:K01101. To get the counts we need to split
        these into separate rows. Note that doing this directly on the
        dataframe will not handle the other columns properly, hence we do it
        only on the annotation column and then merge with the abundances.
        Of course a given annotation can appear several times
        (i.e. for different genes) so we then need to group by annotation
        and sum the abundances.
        """
        col = annotations[colname]
        annotations = col[pd.notna(col)].str.split(",").explode()
        annotations = annotations.to_frame().merge(
            abundances, how="left", left_on="#query", right_index=True
        )
        return annotations.groupby(colname).sum()

    def annotation_file(self, sample):
        if sample is None:
            return os.path.join(
                self.pipeline_outdir, "gene_catalog", "all.emapper.annotations"
            )
        else:
            return sample.annotation_file

    def load_annotations(self, sample=None):
        data = pd.read_csv(
            self.annotation_file(sample), header=4, delimiter="\t", index_col=0
        )
        data.where(data != "-", np.nan, inplace=True)
        return data

    @staticmethod
    def load_abundances(sample):
        return pd.read_csv(
            sample.gene_file, delimiter=",", index_col=0, names=["gene", sample.name]
        )


if __name__ == "__main__":
    args = parse_arguments(
        samples_file="optional",
        pipeline_outdir=True,
        postprocessed_dir=True,
        per_sample=True,
    )

    AnnotationAbundanceCalculator(
        args.samples, args.pipeline_outdir, args.postprocessed_dir, args.per_sample
    )()
