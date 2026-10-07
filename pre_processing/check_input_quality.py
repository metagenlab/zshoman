"""
Pre-run quality control for zshoman input FASTQ files. 

This script checks sample-sheet consistency, input file integrity and availability, 
and basic FASTQ statistics before running the pipeline. 
"""

import gzip
import sys
from pathlib import Path
from Bio.SeqIO.QualityIO import FastqGeneralIterator

import pandas as pd

MAX_READS = 1000

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import logger
from utils.utils import parse_arguments

def check_file(path):
    """Check that a file path is valid, exists and is not empty"""
    if pd.isna(path) or str(path).strip() == "":
        return False
	
    path= Path(path)
    return path.exists() and path.stat().st_size > 0

def get_fastq_stats(path, max_reads=MAX_READS):
    """Check the first reads of a FASTQ and collect fast statistics"""
    n_reads = 0
    total_bases = 0
    min_length = None
    max_length = 0

    # Get stats for just the first reads up until MAX_READS
    try:
        with gzip.open(path, "rt") as handle:
            for _, sequence, _ in FastqGeneralIterator(handle):
                if n_reads >= max_reads:
                    break

                length = len(sequence)

                n_reads += 1
                total_bases += length

                if min_length is None or length < min_length:
                    min_length = length

                if length > max_length:
                    max_length = length

    except (OSError, EOFError, ValueError):
        return None

    mean_length = total_bases / n_reads if n_reads > 0 else 0

    return {
        "reads_checked": n_reads,
        "mean_length": mean_length,
        "min_length": min_length,
        "max_length": max_length,
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
        "R1_gzip_ok": None,
        "R2_gzip_ok": None,
        "R1_size_gb": None,
        "R2_size_gb": None,
        "R1_reads_checked": None,
        "R2_reads_checked": None,
        "R1_mean_length": None,
        "R2_mean_length": None,
        "status": "PASS"
    }

    # Check and register size of the files
    if r1_exists:
        result["R1_size_gb"] = Path(r1).stat().st_size / (1024 ** 3)

    if r2_exists:
        result["R2_size_gb"] = Path(r2).stat().st_size / (1024 ** 3)

    # Sample R1
    if r1_exists:
        r1_stats = get_fastq_stats(r1)
        result["R1_gzip_ok"] = r1_stats is not None

        if r1_stats is not None:
            result["R1_reads_checked"] = r1_stats["reads_checked"]
            result["R1_mean_length"] = r1_stats["mean_length"]

    # Sample R2
    if r2_exists:
        r2_stats = get_fastq_stats(r2)
        result["R2_gzip_ok"] = r2_stats is not None

        if r2_stats is not None:
            result["R2_reads_checked"] = r2_stats["reads_checked"]
            result["R2_mean_length"] = r2_stats["mean_length"]

    # Assign status
    if not r1_exists or not r2_exists:
        result["status"] = "MISSING_FILE"

    elif not result["R1_gzip_ok"] or not result["R2_gzip_ok"]:
        result["status"] = "MALFORMED_FASTQ"

    return result

def run_qc(samples_file):
    """Run input QC for all fastq pairs in the sample sheet (lane level)."""
    samples = pd.read_csv(samples_file)

    # Check no paths are duplicated for R1 or R2
    duplicate_r1 = samples["fastq_R1"].duplicated(keep=False)
    duplicate_r2 = samples["fastq_R2"].duplicated(keep=False)

    duplicate_pairs = samples.duplicated(
        subset=["fastq_R1", "fastq_R2"],
        keep=False,
    )
    
    # Check whether we have the same number of lanes in all samples 
    samples["lane_index"] = samples.groupby("sample").cumcount() + 1 
    lane_counts = samples.groupby("sample").size()

    warnings = []

    if lane_counts.nunique() > 1:
        warnings.append(
            "Some sample(s) have a different number of lanes"
        )

    # Run per file - QC
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
            total_R1_size_gb=("R1_size_gb","sum"),
            total_R2_size_gb=("R2_size_gb","sum"),
            n_lanes=("lane_index","count"),
            all_lanes_pass=("lane_pass","all"),
        )
    )

    sample_summary["status"] = sample_summary["all_lanes_pass"].map(
        {True: "PASS", False: "CHECK_LANES"})

    # Summarize for the whole sample list
    sample_sheet_summary = {
    "samples_in_sheet": samples["sample"].nunique(),
    "fastq_pairs_in_sheet": len(samples),
    "duplicate_R1_paths": int(duplicate_r1.sum()),
    "duplicate_R2_paths": int(duplicate_r2.sum()),
    "duplicate_pairs": int(duplicate_pairs.sum()),
    "all_samples_same_n_lanes": lane_counts.nunique() == 1,
    "warnings": warnings,
    }

    # Add file-level problems to high-level summary
    sample_sheet_summary["missing_files"] = int(
        (~qc["R1_exists"]).sum() + (~qc["R2_exists"]).sum()
    )

    sample_sheet_summary["unreadable_fastqs"] = int(
        # Files that existed but could not be opened/parsed
        (qc["R1_gzip_ok"] == False).sum() 
        + (qc["R2_gzip_ok"] == False).sum()
    )

    sample_sheet_summary["status"] = (
        "PASS"
        if (
            sample_sheet_summary["duplicate_pairs"] == 0
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

    print("\nSample Sheet summary:")
    for key, value in sample_sheet_summary.items():
        print(f"{key}: {value}")

    print("\nSample-level QC:")
    print(sample_summary)

    print("\nLane-level QC:")
    print(lane_qc)

