"""Summary I/O operations for collecting results."""

from __future__ import annotations

import glob
import os
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from rna_map_slurm.utils.logging import get_logger

log = get_logger("io.summaries")


def find_fastq_records(directory: str | Path) -> dict[str, int]:
    """Find and count FASTQ records for each barcode from job output files.

    Parses .out files looking for lines like:
    "FastQ records for barcode XXXX: N (M pairs)"

    Args:
        directory: Directory containing .out files to parse.

    Returns:
        Dictionary mapping barcode sequences to read counts.
    """
    barcode_counts: dict[str, int] = defaultdict(int)
    pattern = re.compile(r"FastQ records for barcode (\w+): \d+ \((\d+) pairs?\)")

    for root, _, files in os.walk(str(directory)):
        for file in files:
            if not file.endswith(".out"):
                continue
            counts = _parse_output_file(Path(root) / file, pattern)
            for barcode, count in counts.items():
                barcode_counts[barcode] += count

    return dict(barcode_counts)


def _parse_output_file(
    file_path: Path,
    pattern: re.Pattern[str],
) -> dict[str, int]:
    """Parse a single output file for barcode counts.

    Args:
        file_path: Path to the .out file.
        pattern: Compiled regex pattern to match.

    Returns:
        Dictionary of barcode -> count for this file.
    """
    counts: dict[str, int] = {}

    with open(file_path, encoding="utf-8") as f:
        for line in f:
            match = pattern.search(line)
            if match:
                barcode = match.group(1)
                count = int(match.group(2))
                counts[barcode] = counts.get(barcode, 0) + count

    return counts


def find_mutation_histos_files(base_directory: str | Path) -> list[pd.DataFrame]:
    """Find and load all mutation histogram JSON files.

    Args:
        base_directory: Base directory to search from.

    Returns:
        List of DataFrames loaded from mutation_histos.json files.
    """
    pattern = os.path.join(
        str(base_directory),
        "results/*/processed/*/output/BitVector_Files/mutation_histos.json",
    )

    return [pd.read_json(file_path) for file_path in glob.glob(pattern)]


def get_pop_avg_summary() -> pd.DataFrame:
    """Get summary DataFrame of all population average results.

    Returns:
        Combined DataFrame from all mutation histogram files.
        Empty DataFrame if no files found.
    """
    dfs = find_mutation_histos_files(".")

    if not dfs:
        log.warning("No mutation_histos.json files found")
        log.warning(
            "Expected path: results/*/processed/*/output/BitVector_Files/mutation_histos.json"
        )
        return pd.DataFrame()

    df = pd.concat(dfs, ignore_index=True)

    total_reads = df["num_reads"].sum()
    total_aligns = df["num_aligned"].sum()
    log.info(f"total number of reads from rna_map: {total_reads}")
    log.info(f"total number of aligns from rna_map: {total_aligns}")

    return df


def get_demultiplexing_summary() -> pd.DataFrame:
    """Get summary DataFrame of demultiplexing results.

    Returns:
        DataFrame with barcode sequences, counts, and fractions.
    """
    directory = "jobs/demultiplex/"
    barcode_counts = find_fastq_records(directory)

    df_barcodes = pd.DataFrame(
        list(barcode_counts.items()),
        columns=["sequence", "num_reads"],
    )

    total_sum = np.sum(df_barcodes["num_reads"])
    log.info(f"total number of reads from demultiplexing: {total_sum}")

    df_barcodes["fraction"] = np.round(
        df_barcodes["num_reads"] / total_sum * 100,
        3,
    )

    return df_barcodes
