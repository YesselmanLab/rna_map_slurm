"""Internal demultiplexing operations for rna-map-slurm-runner CLI."""

from __future__ import annotations

import click

from rna_map_slurm.tasks.int_demultiplex import (
    int_demultiplex as task_int_demultiplex,
)
from rna_map_slurm.tasks.int_demultiplex import (
    int_demultiplex_rna_map as task_int_demultiplex_rna_map,
)
from rna_map_slurm.tasks.int_demultiplex import (
    int_demultiplex_rna_map_combine as task_int_demultiplex_rna_map_combine,
)
from rna_map_slurm.utils.logging import setup_logging
from rna_map_slurm.utils.timing import time_it


@time_it
@click.command("int-demultiplex")
@click.argument("construct_barcode")
@click.argument("b1_seq")
@click.argument("b2_seq")
@click.argument("b1_min_pos", type=int)
@click.argument("b1_max_pos", type=int)
@click.argument("b2_min_pos", type=int)
@click.argument("b2_max_pos", type=int)
def int_demultiplex(
    construct_barcode: str,
    b1_seq: str,
    b2_seq: str,
    b1_min_pos: int,
    b1_max_pos: int,
    b2_min_pos: int,
    b2_max_pos: int,
) -> None:
    """Perform internal demultiplexing by barcode positions.

    Arguments:
        construct_barcode: Library barcode sequence.
        b1_seq: First internal barcode sequence.
        b2_seq: Second internal barcode sequence.
        b1_min_pos: Minimum position for barcode 1.
        b1_max_pos: Maximum position for barcode 1.
        b2_min_pos: Minimum position for barcode 2.
        b2_max_pos: Maximum position for barcode 2.
    """
    setup_logging()
    task_int_demultiplex(
        construct_barcode,
        b1_seq,
        b2_seq,
        b1_min_pos,
        b1_max_pos,
        b2_min_pos,
        b2_max_pos,
    )


@time_it
@click.command("int-demultiplex-rna-map")
@click.argument("code")
@click.argument("lib_barcode_seq")
@click.argument("construct_barcode_seq")
def int_demultiplex_rna_map(
    code: str,
    lib_barcode_seq: str,
    construct_barcode_seq: str,
) -> None:
    """Run RNA mapping on internally demultiplexed reads.

    Arguments:
        code: Construct code.
        lib_barcode_seq: Library barcode sequence.
        construct_barcode_seq: Internal construct barcode sequence.
    """
    setup_logging()
    task_int_demultiplex_rna_map(code, lib_barcode_seq, construct_barcode_seq)


@time_it
@click.command("int-demultiplex-rna-map-combine")
@click.argument("barcode_seq")
@click.argument("rna_name")
def int_demultiplex_rna_map_combine(barcode_seq: str, rna_name: str) -> None:
    """Combine internal demultiplexing RNA mapping results.

    Arguments:
        barcode_seq: Barcode sequence.
        rna_name: RNA/construct name.
    """
    setup_logging()
    task_int_demultiplex_rna_map_combine(barcode_seq, rna_name)
