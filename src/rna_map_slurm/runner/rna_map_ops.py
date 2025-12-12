"""RNA-map operations for rna-map-slurm-runner CLI."""

from __future__ import annotations

import click
import pandas as pd

from rna_map_slurm.tasks.basic import rna_map_combine as task_rna_map_combine
from rna_map_slurm.tasks.basic import run_rna_map as task_run_rna_map
from rna_map_slurm.utils.dataframe import get_data_row
from rna_map_slurm.utils.logging import get_logger, setup_logging
from rna_map_slurm.utils.timing import time_it

log = get_logger("runner.rna_map")


@click.command("run-rna-map")
@time_it
@click.argument("fasta_path", type=click.Path(exists=True))
@click.argument("r1_path", type=click.Path(exists=True))
@click.argument("r2_path", type=click.Path(exists=True))
@click.argument("csv_path", type=click.Path(exists=True))
@click.argument("output_dir", type=click.Path(exists=True))
def run_rna_map(
    fasta_path: str,
    r1_path: str,
    r2_path: str,
    csv_path: str,
    output_dir: str,
) -> None:
    """Run RNA mapping on FASTQ files.

    Arguments:
        fasta_path: Path to FASTA reference.
        r1_path: Path to R1 FASTQ file.
        r2_path: Path to R2 FASTQ file.
        csv_path: Path to sequence CSV.
        output_dir: Output directory.
    """
    setup_logging()
    task_run_rna_map(fasta_path, r1_path, r2_path, csv_path, output_dir)


@click.command("rna-map-combine")
@time_it
@click.argument("barcode_seq")
@click.argument("construct")
def rna_map_combine(barcode_seq: str, construct: str) -> None:
    """Combine RNA mapping results from multiple chunks.

    Arguments:
        barcode_seq: Barcode sequence.
        construct: Construct/RNA name.
    """
    setup_logging()

    df = pd.read_csv("data.csv")
    row = get_data_row(df, barcode_seq, construct)

    if row is None:
        log.error(f"No barcode_seq {barcode_seq} with construct {construct} found")
        return

    task_rna_map_combine(row)
