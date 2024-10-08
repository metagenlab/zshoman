#!/usr/bin/env python

"""
This script is adapted from
https://methods-in-microbiomics.readthedocs.io/en/latest/_downloads/a980f6f4bd1d2d49e4965aed8af17fde/scaffold_filter.py
"""
import gzip
import hashlib

import Bio.SeqIO.FastaIO as FastaIO
from Bio.Seq import Seq


def stream_fa(infile):
    if infile.endswith('fasta.gz') or infile.endswith('fa.gz'):
        with gzip.open(infile, 'rt') as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                yield (header, sequence)
    elif infile.endswith('fasta') or infile.endswith('fa'):
        with open(infile, 'rt') as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                yield (header, sequence)
    else:
        raise Exception(f'{infile} not a sequence file.')


def main():
    samplename = "$meta.id"
    seqtype = "scaffold"
    filtersize = 500

    outprefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else samplename
    sequences = []

    for cnt, (header, sequence) in enumerate(stream_fa("$scaffolds"), 1):
        sequence = sequence.upper()
        seqlen = len(sequence)
        sequence_rev = str(Seq(sequence).reverse_complement())
        md5_fw = hashlib.md5(sequence.encode()).hexdigest()
        md5_rev = hashlib.md5(sequence_rev.encode()).hexdigest()
        seqname = f'{samplename}-{seqtype}_{cnt} length={seqlen} orig={header}'
        sequences.append((seqname, sequence, md5_fw, md5_rev, seqlen))

    with open(f'{outprefix}.{seqtype}s.min{filtersize}.fasta', 'w') as handle:
        for (seqname, sequence, md5_fw, md5_rev, seqlen) in sequences:
            if seqlen >= filtersize:
                handle.write(f'>{seqname}\\n{sequence}\\n')

    with open(f'{outprefix}.{seqtype}s.hashes', 'w') as handle:
        for (seqname, sequence, md5_fw, md5_rev, seqlen) in sequences:
            if seqlen >= filtersize:
                handle.write(f'{seqname}\\t{md5_fw}\\t{md5_rev}\\t{seqlen}\\n')


if __name__ == '__main__':
    main()
