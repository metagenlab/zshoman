#!/usr/bin/env python
"""
This script counts the number of cells from the motus output and uses that
as well as the length of the genes to normalize the gene counts.

Gene counts themselves are mean coverage with edge correction (ignoring edges of
genes when computing coverage as these are underestimated).
"""

from statistics import fmean

import pysam

min_count_length = 10
edge_correction_length = 30


def main():
    samplename = "$meta.id"
    outprefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else samplename

    with open("$motus_profile") as handle:
        for i in range(2):
            handle.readline()

        number_of_cells = 0
        for line in handle:
            number_of_cells += int(float(line.rsplit("\\t", 1)[-1]))

    pysam.sort("-o", f"{samplename}.bam", "$aligned_reads")
    pysam.index(f"{samplename}.bam")

    bamfile = pysam.AlignmentFile(f"{samplename}.bam")

    stats = bamfile.get_index_statistics()

    counts = {}
    coverage = {}
    for region_stat in stats:
        if region_stat.total == 0:
            continue

        # We ignore the counts at the edges of the contig, as they are underestimated due to edge effects.
        # We ensure to always count on at least the min_count_length amino acids.
        contig_length = bamfile.get_reference_length(region_stat.contig)
        delta = min(int((contig_length - min_count_length) / 2), edge_correction_length)
        counts = bamfile.count_coverage(
            region_stat.contig, start=delta, stop=contig_length - delta
        )
        base_coverage = tuple(map(sum, zip(*counts)))
        counts[region_stat.contig] = fmean(base_coverage)
        coverage[region_stat.contig] = 1 - (base_coverage.count(0) / len(base_coverage))

    with open(f"{outprefix}_genes_per_cell.csv", "w") as fout:
        fout.writelines(
            f"{key},{val / number_of_cells}\\n" for key, val in counts.items()
        )

    with open(f"{outprefix}_genes_coverage.csv", "w") as fout:
        fout.writelines(f"{key},{val}\\n" for key, val in coverage.items())


if __name__ == "__main__":
    main()
