"""Clean command for rna-map-slurm CLI."""

from __future__ import annotations

import glob
import os
import shutil

import click

from rna_map_slurm.utils.logging import get_logger, setup_logging

log = get_logger("cli.clean")


@click.command()
@click.argument("stage")
def clean(stage: str) -> None:
    """Clean workflow directories.

    STAGE can be 'all', 'demultiplex', or 'rna_map'.
    """
    setup_logging()

    if stage == "all":
        _clean_all()
    elif stage == "demultiplex":
        _clean_demultiplex()
    elif stage == "rna_map":
        _clean_rna_map()
    else:
        log.error(f"Unknown stage: {stage}. Use 'all', 'demultiplex', or 'rna_map'")


def _clean_all() -> None:
    """Remove all workflow directories."""
    log.info("Cleaning all directories")
    for dir_name in ["jobs", "submits", "data", "inputs", "logs"]:
        shutil.rmtree(dir_name, ignore_errors=True)


def _clean_demultiplex() -> None:
    """Clean demultiplex output directories."""
    log.info("Cleaning demultiplex directories")
    dirs = glob.glob("data/split-*/[ACGT]*")
    for d in dirs:
        if os.path.isdir(d):
            shutil.rmtree(d)


def _clean_rna_map() -> None:
    """Clean RNA-map output directories."""
    log.info("Cleaning rna_map directories")
    dirs = glob.glob("data/split-*/[ACGT]*/*")
    for d in dirs:
        if os.path.isdir(d):
            shutil.rmtree(d)
