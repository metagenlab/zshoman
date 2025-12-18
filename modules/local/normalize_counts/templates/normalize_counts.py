#!/usr/bin/env python
"""
This script counts the number of cells from the motus output and uses that
as well as the length of the genes to normalize the gene counts.
"""

from statistics import fmean

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

    pysam.sort("-o", f"{samplename}.bam", "$aligned_reads")
    pysam.index(f"{samplename}.bam")

    bamfile = pysam.AlignmentFile(f"{samplename}.bam")

    stats = bamfile.get_index_statistics()

    counts = {}
    for region_stat in stats:
        if region_stat.total == 0:
            continue
        counts[region_stat.contig] = sum(
            map(fmean, bamfile.count_coverage(region_stat.contig))
        )

    with open(f"{outprefix}_genes_per_cell.csv", "w") as fout:
        for key, val in counts.items():
            fout.write(f"{key},{val / number_of_cells}\\n")


if __name__ == "__main__":
    main()
