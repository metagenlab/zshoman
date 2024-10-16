#!/usr/bin/env python
"""
This script counts the number of cells from the motus output and uses that
as well as the length of the genes to normalize the gene counts.
"""
from collections import Counter

import pysam


def main():
    samplename = "$meta.id"
    outprefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else samplename

    with open("$motus_profile") as handle:
        for i in range(3):
            handle.readline()

        number_of_cells = 0
        for line in handle:
            number_of_cells += int(line.rsplit("\\t", 1)[-1])

    samfile = pysam.AlignmentFile("$aligned_reads", "rb")
    read_counts = Counter()
    ref_lengths = {}
    for read in samfile:
        name = read.reference_name
        read_counts[name] += 1
        if name not in ref_lengths:
            ref_lengths[name] = read.reference_length

    with open(f"{outprefix}_genes_per_cell.csv", "w") as fout:
        for key, val in read_counts.items():
            fout.write(f"{key},{val/ref_lengths[key]/number_of_cells}\\n")


if __name__ == '__main__':
    main()
