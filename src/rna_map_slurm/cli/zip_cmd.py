"""Zip demultiplexed files command for rna-map-slurm CLI."""

from __future__ import annotations

import os
import zipfile

import click
import pandas as pd

from rna_map_slurm.utils.logging import get_logger

log = get_logger("cli.zip")


@click.command()
@click.argument("csv")
def zip_demultiplex_subset(csv: str) -> None:
    """Create a zip of demultiplexed files for barcodes in CSV."""
    df = pd.read_csv(csv)
    barcodes: list[str] = df["barcode_seq"].unique().tolist()

    _create_zip(barcodes)


def _create_zip(barcodes: list[str]) -> None:
    """Create zip file containing all barcode directories."""
    base_path = "demultiplexed/"
    zip_filename = "demultiplex_subset.zip"

    with zipfile.ZipFile(zip_filename, "w", zipfile.ZIP_DEFLATED) as zipf:
        for barcode in barcodes:
            _add_barcode_to_zip(zipf, base_path, barcode)

    log.info(f"Created zip file {zip_filename}")


def _add_barcode_to_zip(
    zipf: zipfile.ZipFile,
    base_path: str,
    barcode: str,
) -> None:
    """Add all .gz files for a barcode to the zip."""
    barcode_path = os.path.join(base_path, barcode)

    if not os.path.exists(barcode_path):
        log.warning(f"Barcode path {barcode_path} does not exist, skipping.")
        return

    for root, _, files in os.walk(barcode_path):
        for file in files:
            if not file.endswith(".gz"):
                continue
            file_path = os.path.join(root, file)
            arcname = os.path.relpath(file_path, base_path)
            zipf.write(file_path, arcname)
            log.info(f"Added {file_path} to zip")
