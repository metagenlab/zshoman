#!/usr/bin/env python

"""
This script will filter hits output from blast to keep only the best hit for each
gene.
"""

import pandas as pd


def main():
    samplename = "$meta.id"
    min_seqid = float("$params.custom_annotation_seqid")

    outprefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else samplename
    hits = pd.read_csv(
        "$hits",
        sep="\t",
        header=None,
        names=["crc", "gene_id", "seqid", "leng", "evalue", "score", "qcov"],
    )

    # filter out hits with seqid too low
    hits = hits[hits["seqid"] >= min_seqid]

    # Select only the best hits
    hits = hits.sort_values(
        ["evalue", "seqid", "qcov"], ascending=[True, False, False]
    ).drop_duplicates("crc")
    hits.to_csv(outprefix + "_best_hits.csv", index=False)

    with open(f"{outprefix}_best_hits.log", "w") as logfile:
        logfile.write("sample, #hits\\n")
        logfile.write(f"{samplename}, {len(hits)}\\n")


if __name__ == "__main__":
    main()
