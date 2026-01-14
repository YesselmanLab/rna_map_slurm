"""Main CLI entry point for rna-map-slurm-runner."""

from __future__ import annotations

import click

from rna_map_slurm.runner.fastq_ops import (
    demultiplex,
    join_fastq_files,
    split_fastq,
    split_fastqs,
)
from rna_map_slurm.runner.int_demultiplex_ops import (
    int_demultiplex,
    int_demultiplex_cpp,
    int_demultiplex_rna_map,
    int_demultiplex_rna_map_combine,
)
from rna_map_slurm.runner.rna_map_ops import rna_map_combine, run_rna_map


@click.group()
def cli() -> None:
    """RNA Map SLURM Runner - Execute individual workflow tasks."""
    pass


cli.add_command(split_fastq)
cli.add_command(split_fastqs)
cli.add_command(demultiplex)
cli.add_command(join_fastq_files)
cli.add_command(run_rna_map)
cli.add_command(rna_map_combine)
cli.add_command(int_demultiplex)
cli.add_command(int_demultiplex_cpp)
cli.add_command(int_demultiplex_rna_map)
cli.add_command(int_demultiplex_rna_map_combine)


if __name__ == "__main__":
    cli()
