import os
import re
import csv
from collections import defaultdict
import pandas as pd
import numpy as np 

from rna_map_slurm.logger import get_logger

log = get_logger("SUMMARIES")


def find_fastq_records(directory):
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


def get_demultiplexing_summary():
    directory = "jobs/demultiplex/"
    barcode_counts = find_fastq_records(directory)
    df_barcodes = pd.DataFrame(list(barcode_counts.items()), columns=['sequence', 'num_reads'])
    total_sum = np.sum(df_barcodes["num_reads"])
    log.info(f"total number of reads: {total_sum}")
    df_barcodes["fraction"] = np.round(df_barcodes["num_reads"] / total_sum * 100, 3) 
    df_barcodes.to_csv("results/demultiplexing.csv", index=False)
     
    
