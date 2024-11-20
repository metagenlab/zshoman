#!/usr/bin/env python

"""
This script was adapted from
https://methods-in-microbiomics.readthedocs.io/en/latest/_downloads/a980f6f4bd1d2d49e4965aed8af17fde/scaffold_filter.py
and complemented with classification of contigs into prokaryotic and eukaryotic.
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


class DomainClassifier():

    def __init__(self, classification_file):
        self.eukaryotes = set()
        with open(classification_file, "r") as infile:
            for line in infile:
                name, domain = line.split(",")
                if domain.strip() == "eukarya":
                    self.eukaryotes.add(name.strip())

    def is_eukariotic(self, name):
        return name in self.eukaryotes


def main():
    samplename = "$meta.id"
    seqtype = "scaffold"
    filtersize = 500

    outprefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else samplename

    classifier = DomainClassifier("$classification")

    filtered = 0
    prefix = f'{outprefix}.{seqtype}s.min{filtersize}'
    with open(f'{prefix}.fasta', 'w') as handle_all, \
         open(f'{prefix}.hashes', 'w') as handle_hashes, \
         open(f'{prefix}_eukaryotes.fasta', 'w') as handle_euk, \
         open(f'{prefix}_prokaryotes.fasta', 'w') as handle_prok:
        for cnt, (header, sequence) in enumerate(stream_fa("$scaffolds"), 1):
            sequence = sequence.upper()
            seqlen = len(sequence)
            if seqlen < filtersize:
                filtered += 1
                continue
            sequence_rev = str(Seq(sequence).reverse_complement())
            md5_fw = hashlib.md5(sequence.encode()).hexdigest()
            md5_rev = hashlib.md5(sequence_rev.encode()).hexdigest()
            seqname = f'{samplename}-{seqtype}_{cnt} length={seqlen} orig={header}'
            handle_all.write(f'>{seqname}\\n{sequence}\\n')
            handle_hashes.write(f'{seqname}\\t{md5_fw}\\t{md5_rev}\\t{seqlen}\\n')

            if classifier.is_eukariotic(header):
                handle_euk.write(f'>{seqname}\\n{sequence}\\n')
            else:
                handle_prok.write(f'>{seqname}\\n{sequence}\\n')
    with open(f'{outprefix}.{seqtype}s.scaffold_filter.log', "w") as logfile:
        logfile.write("sample, #scaffolds, #filtered scaffolds\\n")
        logfile.write(f"{samplename}, {cnt + 1}, {filtered}\\n")


if __name__ == '__main__':
    main()
