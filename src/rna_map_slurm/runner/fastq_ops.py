"""FASTQ operations for rna-map-slurm-runner CLI."""

from __future__ import annotations

import glob
import os
import subprocess

import click
import pandas as pd

from rna_map_slurm.tasks.basic import demultiplex as task_demultiplex
from rna_map_slurm.tasks.basic import split_fastq_file
from rna_map_slurm.utils.logging import get_logger, setup_logging
from rna_map_slurm.utils.timing import time_it

log = get_logger("runner.fastq")


@time_it
@click.command("split-fastqs")
@click.argument("r1_path", type=click.Path(exists=True), required=True)
@click.argument("r2_path", type=click.Path(exists=True), required=True)
@click.argument("output_dir", type=click.Path(exists=True), required=True)
@click.argument("num_chunks", type=int, required=True)
@click.option(
    "--start",
    default=0,
    show_default=True,
    help="Starting index for chunk numbering.",
)
@click.option(
    "--threads",
    default=1,
    show_default=True,
    help="Number of threads for splitting.",
)
def split_fastqs(
    r1_path: str,
    r2_path: str,
    output_dir: str,
    num_chunks: int,
    start: int,
    threads: int,
) -> None:
    """Split paired FASTQ files into multiple chunks.

    Arguments:
        r1_path: Path to R1 FASTQ file.
        r2_path: Path to R2 FASTQ file.
        output_dir: Directory for output files.
        num_chunks: Number of chunks to create.
    """
    setup_logging()
    split_fastq_file(r1_path, output_dir, num_chunks, start, threads)
    split_fastq_file(r2_path, output_dir, num_chunks, start, threads)


@time_it
@click.command()
@click.argument("csv")
@click.argument("r1_path", type=click.Path(exists=True))
@click.argument("r2_path", type=click.Path(exists=True))
@click.argument("output_dir", type=click.Path(exists=True))
def demultiplex(csv: str, r1_path: str, r2_path: str, output_dir: str) -> None:
    """Demultiplex paired FASTQ files by 3' barcodes.

    Arguments:
        csv: Path to CSV with barcode information.
        r1_path: Path to R1 FASTQ file.
        r2_path: Path to R2 FASTQ file.
        output_dir: Output directory.
    """
    setup_logging()
    task_demultiplex(csv, r1_path, r2_path, output_dir)


@time_it
@click.command("join-fastq-files")
def join_fastq_files() -> None:
    """Join demultiplexed FASTQ files by barcode."""
    setup_logging()

    os.makedirs("demultiplexed", exist_ok=True)
    df = pd.read_csv("data.csv")

    for barcode, _group in df.groupby("barcode_seq"):
        _join_barcode_files(str(barcode))


def _join_barcode_files(barcode: str) -> None:
    """Join all FASTQ files for a single barcode.

    Args:
        barcode: Barcode sequence.
    """
    os.makedirs(f"demultiplexed/{barcode}", exist_ok=True)

    r1_files = glob.glob(f"data/*/{barcode}/test_R1.fastq.gz")
    r2_files = glob.glob(f"data/*/{barcode}/test_R2.fastq.gz")

    log.info(f"joining {barcode} files")
    log.info(f"r1_files: {len(r1_files)}")
    log.info(f"outputing to: demultiplexed/{barcode}/test_R1.fastq.gz")

    _concatenate_files(r1_files, f"demultiplexed/{barcode}/test_R1.fastq.gz")

    log.info(f"r2_files: {len(r2_files)}")
    log.info(f"outputing to: demultiplexed/{barcode}/test_R2.fastq.gz")

    _concatenate_files(r2_files, f"demultiplexed/{barcode}/test_R2.fastq.gz")


def _concatenate_files(input_files: list[str], output_file: str) -> None:
    """Concatenate multiple files into one.

    Args:
        input_files: List of input file paths.
        output_file: Output file path.
    """
    files_str = " ".join(input_files)
    subprocess.run(f"cat {files_str} > {output_file}", shell=True, check=True)
