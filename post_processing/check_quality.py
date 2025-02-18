"""
This script will extract information from the log files
and prepare summary tables and plots, notably to check
the quality of the data and the run.
"""

import argparse
import logging
import os
import re
from collections import Counter
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Post-processing")

plt.style.use(["seaborn-v0_8", "seaborn-v0_8-muted"])


class ProcessWithRegex:
    def __init__(self, logfile):
        self.file = logfile

    def __call__(self):
        if not self.file.exists():
            return Counter()
        with self.file.open() as handler:
            content = handler.read()
        res = Counter()
        for pattern, match_type, transform in self.get_patterns():
            if match_type == "first":
                match = pattern.search(content)
                if match:
                    res.update(
                        {
                            key: transform(value)
                            for key, value in match.groupdict().items()
                        }
                    )
            elif match_type == "all":
                for match in pattern.finditer(content):
                    res.update(
                        {
                            key: transform(value)
                            for key, value in match.groupdict().items()
                        }
                    )
        return res

    def get_patterns(self):
        for pattern in self.patterns:
            if isinstance(pattern, tuple):
                yield pattern
            else:
                yield pattern, "all", float


class ProcessBBduk(ProcessWithRegex):
    search_terms = (
        ("KTrimmed", "ktrimmed"),
        ("QTrimmed", "qtrimmed"),
        ("Low quality discards", "lowq_discards"),
        ("Contaminants", "contaminants"),
        ("Total Removed", "removed"),
        ("Result", "result"),
    )

    search_string = (
        r"^{0}:\s+(?P<{1}_reads>\d+) reads "
        r"\(.*\)\s+(?P<{1}_bases>\d+) bases \(.*\)"
    )

    patterns = [
        re.compile(
            r"^Input:\s+(?P<reads>\d+) reads\s+(?P<bases>\d+) bases.", re.MULTILINE
        )
    ] + [
        re.compile(
            r"^{0}:\s+(?P<{1}_reads>\d+) reads "
            r"\(.*\)\s+(?P<{1}_bases>\d+) bases \(.*\)".format(*terms),
            re.MULTILINE,
        )
        for terms in search_terms
    ]


class ProcessBBMap(ProcessWithRegex):
    patterns = [
        re.compile(
            r"^Reads Used:\s+(?P<reads>\d+)\s+\((?P<bases>\d+) bases\)", re.MULTILINE
        ),
        re.compile(
            r"^mapped:\s+[0-9\.]+%\s+(?P<mapped_reads>\d+)\s+[0-9\.]+%\s+(?P<mapped_bases>\d+)",
            re.MULTILINE,
        ),
    ]


class ProcessBBMerge(ProcessWithRegex):
    patterns = [
        re.compile(r"^Pairs:\s+(?P<reads>\d+)", re.MULTILINE),
        re.compile(r"^Joined:\s+(?P<joined_reads>\d+)", re.MULTILINE),
    ]


class ProcessPhantaLog(ProcessWithRegex):
    patterns = [
        (re.compile(r"^(?P<reads>\d+) sequences", re.MULTILINE), "first", int),
        (
            re.compile(r"^\s*(?P<classified>\d+) sequences classified", re.MULTILINE),
            "first",
            int,
        ),
    ]


class ProcessMOTUsLog(ProcessWithRegex):
    patterns = [
        (
            re.compile(r"^\s*Total number of reads: (?P<reads>\d+)", re.MULTILINE),
            "all",
            int,
        ),
        (
            re.compile(
                r"^\s*Number of reads after filtering: (?P<mapped>\d+) \([0-9\.]+ percent\)",
                re.MULTILINE,
            ),
            "all",
            int,
        ),
    ]


def check_consistency(data):
    res_reads = data["trim_adapters"]["reads"]
    for el in ["trim_adapters", "filter_phix", "filter_quality"]:
        assert data[el]["reads"] - data[el]["removed_reads"] == data[el]["result_reads"]
        # Check that input number of reads is ouput of previous step
        assert data[el]["reads"] == res_reads
        res_reads = data[el]["result_reads"]
    assert data["filter_host"]["reads"] == res_reads


def plot_preprocessing_histogram(data, output_dir, type="reads"):
    initial = [el.get("trim_adapters", {}).get(type) for el in data.values()]
    initial = list(filter(None, initial))
    preprocessed = []
    for el in data.values():
        if el.get("filter_host", {}).get(type) is None:
            continue
        res = el["filter_host"][type] - el["filter_host"][f"mapped_{type}"]
        preprocessed.append(res)
    plt.figure()
    logbins = np.logspace(
        np.log10(np.min(preprocessed) - 1), np.log10(np.max(initial) + 1), 30
    )
    plt.hist([initial, preprocessed], label=["initial", "preprocessed"], bins=logbins)
    plt.xscale("log")
    plt.xlabel(f"# {type}")
    plt.ylabel("# samples")
    plt.title(f"Preprocessing {type}")
    plt.legend(loc="best")
    plt.savefig(Path(output_dir, f"{type}_preprocessing_hist.png"))
    plt.close()


def get_fractions(data, step, numerator, denominator):
    res = []
    for sample in data.values():
        if not sample.get(step):
            continue
        pp_step = sample[step]
        if numerator not in pp_step or denominator not in pp_step:
            continue
        res.append(pp_step[numerator] / pp_step[denominator])
    return np.array(res)


def plot_preprocessing_boxplot(data, output_dir, type="reads"):
    trim_adapters = get_fractions(data, "trim_adapters", f"removed_{type}", type)
    filter_phix = get_fractions(data, "filter_phix", f"removed_{type}", type)
    filter_quality = get_fractions(data, "filter_quality", f"removed_{type}", type)
    filter_host = get_fractions(data, "filter_host", f"mapped_{type}", type)
    plt.figure()
    plt.boxplot(
        [trim_adapters, filter_phix, filter_quality, filter_host],
        tick_labels=[
            "Adapter trimming",
            "PhiX filtering",
            "Quality filtering",
            "Host filtering",
        ],
    )
    plt.ylabel("Fraction removed")
    plt.yscale("log")
    plt.title(f"Fraction of {type} removed during preprocessing")
    plt.savefig(Path(output_dir, f"{type}_preprocessing_boxplot.png"))
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser("check_quality.py")
    parser.add_argument("samples_file", help="path to samples csv file.")
    parser.add_argument(
        "-i",
        "--input_dir",
        default="output",
        help="path to the output directory of the pipeline",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        default="analyses/",
        help="path where the tables and plots should be saved to",
    )

    args = parser.parse_args()
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    samples = pd.read_csv(args.samples_file, header=0)["sample"]
    logger.info(f"Found {len(samples)} samples.")
    log_dir = Path(args.input_dir, "logs")

    data = defaultdict(dict)
    for i, sample in enumerate(samples, 1):
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

        phanta_log = Path(log_dir, f"{sample}.phanta.log")
        data[sample]["phanta"] = ProcessPhantaLog(phanta_log)()

        motus_log = Path(log_dir, f"{sample}.motus.log")
        data[sample]["motus"] = ProcessMOTUsLog(motus_log)()

        check_consistency(data[sample])

        if i % 50 == 0:
            logger.info(f"Done {i}/{len(samples)}")

    plot_preprocessing_histogram(data, args.output_dir, "reads")
    plot_preprocessing_histogram(data, args.output_dir, "bases")
    plot_preprocessing_boxplot(data, args.output_dir, "reads")
    plot_preprocessing_boxplot(data, args.output_dir, "bases")

    # Write out the data as a table
    logs = [
        "trim_adapters",
        "filter_phix",
        "filter_quality",
        "filter_host",
        "filter_host_singletons",
        "merge_reads",
        "phanta",
        "motus",
    ]
    columns = [(log, key) for log in logs for key in data[samples[0]][log].keys()]
    df = pd.DataFrame(
        index=samples,
        columns=[f"{log}-{key}" for log, key in columns],
        data=[[data[sample][log][key] for log, key in columns] for sample in samples],
    )
    df.to_csv(Path(args.output_dir, "statistics.csv"))
