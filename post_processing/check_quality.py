"""
This script will extract information from the log files
and prepare summary tables and plots, notably to check
the quality of the data and the run.
"""

import argparse
import logging
import os
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Post-processing")


class ProcessBBduk():

    patterns = [
        re.compile(r"^Input:\s+(?P<reads>\d+) reads\s+(?P<bases>\d+) bases.",
                   re.MULTILINE),
        re.compile(r"^KTrimmed:\s+(?P<trimmed_reads>\d+) reads "
                   r"\(.*\)\s+(?P<trimmed_bases>\d+) bases \(.*\)",
                   re.MULTILINE),
        re.compile(r"^Contaminants:\s+(?P<contaminants_reads>\d+) reads "
                   r"\(.*\)\s+(?P<contaminants_bases>\d+) bases \(.*\)",
                   re.MULTILINE),
        re.compile(r"^Total Removed:\s+(?P<removed_reads>\d+) reads "
                   r"\(.*\)\s+(?P<removed_bases>\d+) bases \(.*\)",
                   re.MULTILINE),
        re.compile(r"^Result:\s+(?P<output_reads>\d+) reads "
                   r"\(.*\)\s+(?P<output_bases>\d+) bases \(.*\)",
                   re.MULTILINE),
    ]

    def __init__(self, logfile):
        self.file = logfile

    def __call__(self):
        with self.file.open() as handler:
            content = handler.read()
        res = {}
        for pattern in self.patterns:
            match = pattern.search(content)
            if match:
                res.update(match.groupdict())
        return res


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
    import pdb; pdb.set_trace()
