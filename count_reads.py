import argparse
from collections import Counter

import pandas
import pysam


def main():
    parser = argparse.ArgumentParser(description='Count reads per gene')
    parser.add_argument('gene_counts_file', type=str, help='input file')
    parser.add_argument('motus_file', type=str, help='input file')
    parser.add_argument('output', type=str, help='input file')
    args = parser.parse_args()

    df = pandas.read_csv(args.motus_file, header=2, sep="\t")
    number_of_cells = df.iloc[:, -1].sum()

    samfile = pysam.AlignmentFile(args.gene_counts_file, "rb")
    read_counts = Counter()
    ref_lengths = {}
    for read in samfile:
        name = read.reference_name
        read_counts[name] += 1
        if name not in ref_lengths:
            ref_lengths[name] = read.reference_length

    with open(args.output, "w") as fout:
        for key, val in read_counts.items():
            fout.write(f"{key},{val/ref_lengths[key]/number_of_cells}\n")


if __name__ == '__main__':
    main()
