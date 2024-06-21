import os
import subprocess
import pandas as pd
import shutil
from typing import List, Optional
from tabulate import tabulate

from fastqsplitter import split_fastqs as fastqsplitter

from rna_map_slurm.fastq import PairedFastqFiles, FastqFile
from rna_map_slurm.demultiplex import SabreDemultiplexer
from rna_map_slurm.logger import get_logger

log = get_logger("steps")


def split_fastq_file_task(
    fastq_file: str, output_dir: str, num_chunks: int, start: int, threads: int = 1
) -> List[str]:
    """
    Splits a FASTQ file into multiple chunks and saves the chunks to the specified output directory.

    Args:
        fastq_file (str): Path to the input FASTQ file.
        output_dir (str): Directory where the split files will be saved.
        num_chunks (int): Number of chunks to split the FASTQ file into.
        start (int): Starting index for the output chunk files.
        threads (int, optional): Number of threads to use for splitting. Defaults to 1.

    Returns:
        List[str]: List of paths to the split output files.
    """
    log.info(f"Splitting {fastq_file} into {num_chunks} chunks")
    log.info(
        f"output_dir: {output_dir}, num_chunks: {num_chunks}, start: {start}, threads: {threads}"
    )
    output_file = "test_R1.fastq.gz"
    if "R2" in fastq_file:
        log.info("Detected R2 file")
        output_file = "test_R2.fastq.gz"
    else:
        log.info("Detected R1 file")
    output_files = []
    for i in range(start, num_chunks + start):
        output_files.append(f"{output_dir}/split-{i:04}/{output_file}")
    fastqsplitter(fastq_file, output_files, threads_per_file=threads)
    return output_files


def demultiplex_task(
    csv: str, r1_path: str, r2_path: str, output_dir: Optional[str] = None
) -> List[str]:
    """
    Executes the demultiplexing task using the provided CSV and FASTQ files, and returns the list of output directories.

    Args:
        csv (str): Path to the CSV file containing barcode sequences and other metadata.
        r1_path (str): Path to the R1 FASTQ file.
        r2_path (str): Path to the R2 FASTQ file.
        output_dir (Optional[str]): Directory where the output files will be saved. If None, uses the current working directory.

    Returns:
        List[str]: List of directories where the demultiplexed files are saved.
    """
    if output_dir is None:
        output_dir = os.getcwd()
    cur_dir = os.getcwd()

    log.info(f"Changing directory to {output_dir}")
    os.chdir(output_dir)

    paired_fastqs = PairedFastqFiles(FastqFile(r1_path), FastqFile(r2_path))
    df = pd.read_csv(csv)
    demultiplexer = SabreDemultiplexer()
    demultiplexer.run(df, paired_fastqs, output_dir)
    os.chdir(cur_dir)

    output_dirs = [
        os.path.join(output_dir, df_row["barcode_seq"]) for i, df_row in df.iterrows()
    ]

    return output_dirs


def join_fastq_files_task(fastq_files, joined_fastq):
    os.system(f"cat {' '.join(fastq_files)} > {joined_fastq}")
    return joined_fastq


def rna_map_task(r1_path, r2_path, fa_path, csv_path, output_dir):
    cur_dir = os.getcwd()
    os.chdir(output_dir)
    subprocess.run(
        f"rna-map -fa {fa_path} -fq1 {r1_path} -fq2 {r2_path} "
        f"--dot-bracket {csv_path} --summary-output-only --param-preset barcoded-library",
        shell=True,
    )
    shutil.rmtree("log")
    shutil.rmtree("input")
    shutil.rmtree("output/Mapping_Files")
    os.chdir(cur_dir)
