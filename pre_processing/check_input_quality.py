"""
This script runs a quick quality check for zshoman input FASTQ files. 
It checks sample-sheet consistency and input file integrity/availability. 
It generates a report that will be saved in the same directory as the sample sheet. 
"""

import gzip
import sys
from pathlib import Path

import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import logger
from utils.utils import parse_arguments

MAX_READS = 1000 # This is the number of reads scanned per FASTQ

def check_file(path):
    """Check that a file path is valid, exists and is not empty"""
    path= Path(path)
    return path.exists() and path.stat().st_size > 0

def sample_fastq_stats(path, max_reads=MAX_READS):
    """Check the first reads of a FASTQ and collect fast statistics"""
    n_reads = 0
    total_bases = 0

    # Get mean length for the first reads up until MAX_READS
    try:
        with gzip.open(path, "rt") as handle:
            for _, sequence, _ in FastqGeneralIterator(handle):
                if n_reads >= max_reads:
                    break

                length = len(sequence)

                n_reads += 1
                total_bases += length

    except (OSError, EOFError, ValueError):
        return None

    mean_length = total_bases / n_reads if n_reads > 0 else 0

    return {
        "reads_checked": n_reads,
        "mean_length": mean_length,
    }

def check_sample(sample, lane_index, r1, r2):
    """Run input QC checks for one paired-end sample."""

    logger.info(f"Checking {sample}, lane {lane_index}")

    #Check that both files exist
    r1_exists = check_file(r1)
    r2_exists = check_file(r2)

    result = {
        "sample": sample,
        "lane_index": lane_index,
        "R1_exists": r1_exists,
        "R2_exists": r2_exists,
        "R1_readable": None,
        "R2_readable": None,
        "R1_size_mb": None,
        "R2_size_mb": None,
        "R1_reads_checked": None,
        "R2_reads_checked": None,
        "R1_mean_length": None,
        "R2_mean_length": None,
        "status": "PASS"
    }

    # Check and register size of the files in mb
    if r1_exists:
        result["R1_size_mb"] = Path(r1).stat().st_size / 1_000_000
    if r2_exists:
        result["R2_size_mb"] = Path(r2).stat().st_size / 1_000_000

    # Sample R1
    if r1_exists:
        r1_stats = sample_fastq_stats(r1)
        result["R1_readable"] = r1_stats is not None

        if r1_stats is not None:
            result["R1_reads_checked"] = r1_stats["reads_checked"]
            result["R1_mean_length"] = r1_stats["mean_length"]

    # Sample R2
    if r2_exists:
        r2_stats = sample_fastq_stats(r2)
        result["R2_readable"] = r2_stats is not None

        if r2_stats is not None:
            result["R2_reads_checked"] = r2_stats["reads_checked"]
            result["R2_mean_length"] = r2_stats["mean_length"]

    # Assign status
    if not r1_exists or not r2_exists:
        result["status"] = "MISSING_FILE"

    elif not result["R1_readable"] or not result["R2_readable"]:
        result["status"] = "UNREADABLE_FASTQ"

    return result

def run_qc(samples_file):
    """Run input QC for all fastq pairs in the sample sheet (lane level)."""
    samples = pd.read_csv(samples_file)

    # Check no paths are duplicated for R1 or R2
    duplicate_r1 = samples["fastq_R1"].duplicated(keep=False)
    duplicate_r2 = samples["fastq_R2"].duplicated(keep=False)

    duplicate_r1_paths = (
        samples.loc[duplicate_r1, "fastq_R1"]
        .drop_duplicates()
        .tolist()
        )

    duplicate_r2_paths = (
        samples.loc[duplicate_r2,"fastq_R2"]
        .drop_duplicates()
        .tolist()
        )

    # Check whether we have the same number of lanes in all samples 
    samples["lane_index"] = samples.groupby("sample").cumcount() + 1 
    lane_counts = samples.groupby("sample").size()

    warnings = []

    if lane_counts.nunique() > 1:
        warnings.append(
            "Some sample(s) have a different number of lanes; see Sample-level QC"
        )

    # Run lane-level QC
    results = []

    for _, row in samples.iterrows():
        result = check_sample(
            row["sample"],
            row["lane_index"],
            row["fastq_R1"],
            row["fastq_R2"],
        )
        results.append(result)

    qc = pd.DataFrame(results)

    # Aggregate lane results per sample

    qc["lane_pass"] = qc["status"] == "PASS"

    sample_summary = (
        qc.groupby("sample", as_index=False)
        .agg(
            total_R1_size_mb=("R1_size_mb","sum"),
            total_R2_size_mb=("R2_size_mb","sum"),
            n_lanes=("lane_index","count"),
            all_lanes_pass=("lane_pass","all"),
        )
    )

    sample_summary["total_R1_size_mb"] = sample_summary["total_R1_size_mb"].round(1)
    sample_summary["total_R2_size_mb"] = sample_summary["total_R2_size_mb"].round(1)

    sample_summary["status"] = sample_summary["all_lanes_pass"].map(
        {True: "PASS", False: "CHECK_LANES"})

    # Reorder columns 
    sample_summary = sample_summary[
        [
            "sample",
            "n_lanes",
            "total_R1_size_mb",
            "total_R2_size_mb",
            "all_lanes_pass",
            "status",
        ]
    ]

    # Summarize for the whole sample list
    sample_sheet_summary = {
    "samples_in_sheet": samples["sample"].nunique(),
    "fastq_pairs_in_sheet": len(samples),
    "all_samples_same_n_lanes": lane_counts.nunique() == 1,
    "duplicate_R1_files": len(duplicate_r1_paths),
    "duplicate_R2_files": len(duplicate_r2_paths),
    "duplicated_R1_paths": duplicate_r1_paths,
    "duplicated_R2_paths": duplicate_r2_paths,
    "warnings": warnings,
    }

    # Add file-level problems to high-level summary
    sample_sheet_summary["missing_files"] = int(
        (~qc["R1_exists"]).sum() + (~qc["R2_exists"]).sum()
    )

    sample_sheet_summary["unreadable_fastqs"] = int(
        # Files that existed but could not be opened/parsed
        (qc["R1_readable"] == False).sum() 
        + (qc["R2_readable"] == False).sum()
    )

    sample_sheet_summary["status"] = (
        "PASS"
        if (
            sample_sheet_summary["duplicate_R1_files"] == 0
            and sample_sheet_summary["duplicate_R2_files"] ==0
            and sample_sheet_summary["missing_files"] == 0
            and sample_sheet_summary["unreadable_fastqs"] == 0
            )
        else "FAIL"
        
    )

    return sample_sheet_summary, sample_summary, qc

if __name__ == "__main__":
    args = parse_arguments(
        samples_file="mandatory",
    )

    sample_sheet_summary , sample_summary, lane_qc = run_qc(args.samples_file)

    output_file = Path(args.samples_file).parent / "input_qc_report.txt"

    with open(output_file, "w") as report:
        report.write("Sample Sheet Summary:\n")
        report.write("=====================\n")
        for key,value in sample_sheet_summary.items():
            report.write(f"{key}: {value}\n")
        
        report.write("\nSample-level QC:\n")
        report.write("=================\n")
        report.write(sample_summary.to_string(index=False))

        report.write("\n\nLane-Level QC\n")
        report.write("================\n")
        report.write(lane_qc.to_string(index=False))

    print(f"QC report written to {output_file}")