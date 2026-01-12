"""Summary generation command for rna-map-slurm CLI."""

from __future__ import annotations

import os

import click
import pandas as pd

from rna_map_slurm.io.summaries import get_demultiplexing_summary, get_pop_avg_summary
from rna_map_slurm.utils.logging import setup_logging


@click.command()
def generate_summaries() -> None:
    """Generate summary files for all runs."""
    _initialize_logging()

    df = pd.read_csv("data.csv")
    df_barcodes = get_demultiplexing_summary()

    _save_demultiplexing_summaries(df, df_barcodes)
    _save_pop_avg_summaries()


def _initialize_logging() -> None:
    """Initialize logging for summary generation."""
    if os.path.isfile("logs/generate-summaries.log"):
        os.remove("logs/generate-summaries.log")
    setup_logging(file_name="logs/generate-summaries.log")


def _save_demultiplexing_summaries(df: pd.DataFrame, df_barcodes: pd.DataFrame) -> None:
    """Save demultiplexing summaries for each run."""
    for run_name in df["run_name"].unique():
        df_sub = df.query("run_name == @run_name")
        barcode_seqs = df_sub["barcode_seq"].unique()  # noqa: F841 used in @barcode_seqs
        df_barcodes_sub = df_barcodes.query("sequence in @barcode_seqs")
        df_barcodes_sub.to_csv(
            f"results/{run_name}/summary/demultiplexing.csv",
            index=False,
        )


def _save_pop_avg_summaries() -> None:
    """Save population average summaries for each run."""
    df_summary = get_pop_avg_summary()

    if df_summary.empty:
        click.echo("No population average data to summarize (rna-map results not found)")
        return

    for run_name, group in df_summary.groupby("run_name"):
        group.to_json(
            f"results/{run_name}/summary/summary.json",
            orient="records",
        )
        group_no_data = group.drop(columns=["data"])
        group_no_data.to_csv(
            f"results/{run_name}/summary/summary.csv",
            index=False,
        )
