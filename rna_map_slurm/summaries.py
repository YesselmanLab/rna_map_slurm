import os
import re
import glob
from collections import defaultdict
import pandas as pd
import numpy as np
from typing import Dict

from rna_map_slurm.logger import get_logger

log = get_logger("SUMMARIES")


def find_fastq_records(directory: str) -> Dict[str, int]:
    """
    Finds and counts FastQ records for each barcode in the given directory.

    Args:
        directory (str): The directory to search for FastQ record files.

    Returns:
        Dict[str, int]: A dictionary containing the sum of counts for each barcode.
                        The keys are the barcode names (str) and the values are the
                        corresponding counts (int).
    """
    # Dictionary to hold the sum of counts for each barcode
    barcode_counts = defaultdict(int)
    # Regular expression to match the desired line format
    pattern = re.compile(r"FastQ records for barcode (\w+): \d+ \((\d+) pairs?\)")
    # Traverse the directory
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".out"):
                file_path = os.path.join(root, file)
                with open(file_path, "r") as f:
                    for line in f:
                        match = pattern.search(line)
                        if match:
                            barcode = match.group(1)
                            count = int(match.group(2))
                            barcode_counts[barcode] += count

    return barcode_counts


def find_mutation_histos_files(base_directory):
    # List to hold the file information and DataFrames
    file_info = []

    # Define the glob pattern to find the mutation_histos.json files
    pattern = os.path.join(
        base_directory,
        "results/*/processed/*/output/BitVector_Files/mutation_histos.json",
    )
    # Find all files matching the pattern
    for file_path in glob.glob(pattern):
        # Load the JSON file into a Pandas DataFrame
        df = pd.read_json(file_path)
        # Append the DataFrame to the list
        file_info.append(df)
    return file_info


def get_pop_avg_summary() -> pd.DataFrame:
    dfs = find_mutation_histos_files(".")
    df = pd.concat(dfs)
    total_reads = df["num_reads"].sum()
    total_aligns = df["num_aligned"].sum()
    log.info(f"total number of reads from rna_map: {total_reads}")
    log.info(f"total number of aligns from rna_map: {total_aligns}")
    return df


def get_demultiplexing_summary() -> pd.DataFrame:
    directory = "jobs/demultiplex/"
    barcode_counts = find_fastq_records(directory)
    df_barcodes = pd.DataFrame(
        list(barcode_counts.items()), columns=["sequence", "num_reads"]
    )
    total_sum = np.sum(df_barcodes["num_reads"])
    log.info(f"total number of reads from demultiplexing : {total_sum}")
    df_barcodes["fraction"] = np.round(df_barcodes["num_reads"] / total_sum * 100, 3)
    return df_barcodes
