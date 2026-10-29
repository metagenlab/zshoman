"""
This script runs a quick quality check for zshoman input FASTQ files. 
It checks sample-sheet consistency and input file integrity/availability. 
It generates a report that will be saved in the same directory as the sample sheet. 
"""

import gzip
import sys
from pathlib import Path
from itertools import islice
from pprint import pprint

import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator

sys.path.append(str(Path(__file__).parent.parent))

from utils.utils import logger
from utils.utils import parse_arguments

MAX_READS = 1000 # This is the number of reads scanned per FASTQ

def sample_fastq_stats(path, max_reads=MAX_READS):
    """Check the first reads of a FASTQ and collect fast statistics"""

    # Get mean length for the first reads up until MAX_READS
    try:
        path = Path(path)
        opener = gzip.open if path.suffix.lower() == ".gz" else open

        with opener(path, "rt") as handle:
            lengths = [
                len(sequence)
                for _, sequence, _ in islice(FastqGeneralIterator(handle), max_reads)
            ]
    except (TypeError, OSError, EOFError, ValueError):
        return None
    
    if not lengths:
        return None

    return {
        "reads_checked": len(lengths),
        "mean_length": sum(lengths) / len(lengths),
        "min_length": min(lengths),
        "max_length": max(lengths),
    }

def qc_file(path):
    """Check one FASTQ file and return its QC result."""
    result = {
        "path": str(path),
        "size_mb": None,
        "reads_checked": None,
        "mean_length": None,
        "min_length": None,
        "max_length": None,
        "status": None,
        "error": None,
    }

    try:
        path = Path(path)
        size_bytes = path.stat().st_size
    except (TypeError, ValueError) as exc:
        result.update(status="INVALID_PATH", error=str(exc))
        return result
    except FileNotFoundError as exc:
        result.update(status="MISSING_FILE", error=str(exc))
        return result
    except OSError as exc:
        result.update(status="UNREADABLE_FILE", error=str(exc))
        return result

    result["size_mb"] = size_bytes / 1_000_000

    if size_bytes == 0:
        result.update(status="EMPTY_FILE", error="File has 0 bytes")
        return result
    
    stats = sample_fastq_stats(path)

    if stats is None:
        result.update(
            status="UNREADABLE_FASTQ",
            error="Could not read FASTQ records from the sampled portion",
        )
        return result
    
    result.update(stats)
    result["status"] = "PASS"
    return result


def collect_file_qc(samples):
    """Run QC on every file and attach its sample information"""
    results = []
    
    for sample_name, sample_data in samples.items():
        for read_direction, paths in(
            ("R1", sample_data["fastq1"]),
            ("R2", sample_data["fastq2"]),
        ):
            for path in paths:
                logger.info(f"Checking {sample_name}, {read_direction}: {path}")

                result = qc_file(path)
                result.update(
                    sample=sample_name,
                    paired_end=sample_data["paired_end"],
                    read_direction=read_direction,
                )
                results.append(result)

    return pd.DataFrame(results)
    
def run_qc(samples):
    """Run file QC and summarize results by sample."""
    qc = collect_file_qc(samples)
    qc["file_pass"] = qc["status"] == "PASS"

    sample_summary = (
        qc.groupby("sample", as_index=False)
        .agg(
            paired_end=("paired_end", "first"),
            n_files=("path", "count"),
            total_size_mb=("size_mb", "sum"),
            all_files_pass=("file_pass", "all"),
        )
    )

    sample_summary["total_size_mb"] = sample_summary["total_size_mb"].round(1)
    sample_summary["status"] = sample_summary["all_files_pass"].map(
        {True: "PASS", False: "FAIL"}
    )

    duplicate_paths = (
        qc.loc[qc["path"].duplicated(keep=False), "path"]
        .drop_duplicates()
        .tolist()
    )

    overall_summary = {
        "n_samples": len(samples),
        "n_files": len(qc),
        "failed_files": int((~qc["file_pass"]).sum()),
        "duplicate_paths": duplicate_paths,
        "status": (
            "PASS"
            if qc["file_pass"].all() and not duplicate_paths
            else "FAIL"
        ),
    }

    return overall_summary, sample_summary, qc

if __name__ == "__main__":
    args = parse_arguments(
        samples_file="mandatory",
    )

    overall_summary, sample_summary, file_qc = run_qc(args.samples)

    output_file = Path(args.samples_file).parent / "input_qc_report.txt"

    with open(output_file, "w") as report:
        report.write("INPUT FASTQ QC REPORT\n")
        report.write("=====================\n")
        report.write(f"Sample sheet: {args.samples_file}\n")
        report.write(f"Reads sampled per file: up to {MAX_READS}\n")

        report.write("\nFILE DETAILS\n")
        report.write("------------\n")

        columns = [
            "sample",
            "read_direction",
            "path",
            "size_mb",
            "reads_checked",
            "mean_length",
            "min_length",
            "max_length",
            "status",
            "error",
        ]

        report.write(
            file_qc[columns].to_string(
                index=False,
                na_rep="-",
                formatters={
                    "size_mb": "{:.6f}".format,
                    "mean_length": "{:.1f}".format,
                },
            )
        )
        report.write("\n")

        report.write("\nSAMPLE SUMMARY\n")
        report.write("--------------\n")

        columns = [
            "sample",
            "paired_end",
            "n_files",
            "total_size_mb",
            "status",
        ]
        report.write(
            sample_summary[columns].to_string(index=False)
        )

        report.write("\n\nOVERALL QC\n")
        report.write("----------\n")
        report.write(f"Samples checked : {overall_summary['n_samples']}\n")
        report.write(f"Files checked   : {overall_summary['n_files']}\n")
        report.write(f"Failed files    : {overall_summary['failed_files']}\n")

        duplicates = overall_summary["duplicate_paths"]
        report.write(f"Duplicate paths : {len(duplicates)}\n")

        for path in duplicates:
            report.writelines(f"  - {path}\n")

        report.write(f"\nRESULT: {overall_summary['status']}\n"
        )    
        
    print(f"QC report written to {output_file}")