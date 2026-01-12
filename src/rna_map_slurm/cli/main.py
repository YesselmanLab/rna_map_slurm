"""Main CLI entry point for rna-map-slurm."""

from __future__ import annotations

import click

from rna_map_slurm.cli.check_cmd import check_jobs
from rna_map_slurm.cli.clean import clean
from rna_map_slurm.cli.config_cmd import generate_example_config
from rna_map_slurm.cli.deposit import deposit_results
from rna_map_slurm.cli.run import run
from rna_map_slurm.cli.setup_cmd import get_data_csv, setup
from rna_map_slurm.cli.summaries import generate_summaries
from rna_map_slurm.cli.summary_cmd import summary
from rna_map_slurm.cli.watch import watch
from rna_map_slurm.cli.zip_cmd import zip_demultiplex_subset


@click.group()
def cli() -> None:
    """RNA Map SLURM - Process RNA sequencing data on SLURM clusters."""
    pass


cli.add_command(setup)
cli.add_command(get_data_csv)
cli.add_command(run)
cli.add_command(check_jobs)
cli.add_command(summary)
cli.add_command(generate_summaries)
cli.add_command(deposit_results)
cli.add_command(clean)
cli.add_command(zip_demultiplex_subset)
cli.add_command(generate_example_config)
cli.add_command(watch)


if __name__ == "__main__":
    cli()
