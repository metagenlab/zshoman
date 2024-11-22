"""
This script will extract information from the log files
and prepare summary tables and plots, notably to check
the quality of the data and the run.
"""

import argparse
import logging
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Post-processing")


class ProcessWithRegex():

    def __init__(self, logfile):
        self.file = logfile

    def __call__(self):
        if not self.file.exists():
            return Counter()
        with self.file.open() as handler:
            content = handler.read()
        res = Counter()
        for pattern in self.patterns:
            for match in pattern.finditer(content):
                res.update({key: float(value)
                            for key, value in match.groupdict().items()})
        return res


class ProcessBBduk(ProcessWithRegex):

    search_terms = (
        ("KTrimmed", "ktrimmed"),
        ("QTrimmed", "qtrimmed"),
        ("Low quality discards", "lowq_discards"),
        ("Contaminants", "contaminants"),
        ("Total Removed", "removed"),
        ("Result", "result"),
    )

    search_string = r"^{0}:\s+(?P<{1}_reads>\d+) reads "\
                    r"\(.*\)\s+(?P<{1}_bases>\d+) bases \(.*\)"

    patterns = [
        re.compile(r"^Input:\s+(?P<reads>\d+) reads\s+(?P<bases>\d+) bases.",
                   re.MULTILINE)] + \
        [
            re.compile(
                r"^{0}:\s+(?P<{1}_reads>\d+) reads "
                r"\(.*\)\s+(?P<{1}_bases>\d+) bases \(.*\)".format(*terms),
                re.MULTILINE)
            for terms in search_terms
        ]


class ProcessBBMap(ProcessWithRegex):

    patterns = [
        re.compile(r"^Reads Used:\s+(?P<reads>\d+)\s+\((?P<bases>\d+) bases\)",
                   re.MULTILINE),
        re.compile(r"^mapped:\s+[0-9\.]+%\s+(?P<mapped_reads>\d+)\s+[0-9\.]+%\s+(?P<mapped_bases>\d+)",
                   re.MULTILINE),
        ]


class ProcessBBMerge(ProcessWithRegex):

    patterns = [
        re.compile(r"^Pairs:\s+(?P<reads>\d+)", re.MULTILINE),
        re.compile(r"^Joined:\s+(?P<joined_reads>\d+)", re.MULTILINE),
        ]


if __name__ == '__main__':
    parser = argparse.ArgumentParser("check_quality.py")
    parser.add_argument("samples_file", help="path to samples csv file.")
    parser.add_argument(
        "-i", "--input_dir", default="output",
        help="path to the output directory of the pipeline")
    parser.add_argument(
        "-o", "--output_dir", default="analyses/",
        help="path where the tables and plots should be saved to")

    args = parser.parse_args()
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    samples = pd.read_csv(args.samples_file, header=0)["sample"]
    log_dir = Path(args.input_dir, "logs")

    data = defaultdict(dict)
    for sample in samples:
        trim_log = Path(log_dir, f"{sample}_trimmed.bbduk.log")
        data[sample]["trim_adapters"] = ProcessBBduk(trim_log)()

        phix_log = Path(log_dir, f"{sample}_phix_filtered.bbduk.log")
        data[sample]["filter_phix"] = ProcessBBduk(phix_log)()

        qf_log = Path(log_dir, f"{sample}_quality_filtered.bbduk.log")
        data[sample]["filter_quality"] = ProcessBBduk(qf_log)()

        hf_log = Path(log_dir, f"{sample}_host_filtered.bbmap.log")
        data[sample]["filter_host"] = ProcessBBMap(hf_log)()

        hf_log = Path(log_dir, f"{sample}_host_filtered_singletons.bbmap.log")
        data[sample]["filter_host_singletons"] = ProcessBBMap(hf_log)()

        merge_log = Path(log_dir, f"{sample}.bbmerge.log")
        data[sample]["merge_reads"] = ProcessBBMerge(merge_log)()
