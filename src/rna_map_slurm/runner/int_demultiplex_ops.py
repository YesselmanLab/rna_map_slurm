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
from rna_map_slurm.tasks.int_demultiplex_cpp import (
    int_demultiplex_batch_cpp as task_int_demultiplex_batch_cpp,
)
from rna_map_slurm.utils.logging import setup_logging
from rna_map_slurm.utils.timing import time_it


@click.command("int-demultiplex")
@time_it
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
    # Convert U to T for DNA sequences (FASTQ uses DNA, not RNA)
    b1_seq = b1_seq.replace("U", "T")
    b2_seq = b2_seq.replace("U", "T")
    task_int_demultiplex(
        construct_barcode,
        b1_seq,
        b2_seq,
        b1_min_pos,
        b1_max_pos,
        b2_min_pos,
        b2_max_pos,
    )


@click.command("int-demultiplex-cpp")
@time_it
@click.argument("lib_barcode")
@click.argument("barcode_json")
def int_demultiplex_cpp(lib_barcode: str, barcode_json: str) -> None:
    """Perform batch internal demultiplexing using C++ (218x faster).

    Processes all internal barcodes in a single pass through the FASTQ files.

    Arguments:
        lib_barcode: Library barcode sequence.
        barcode_json: Path to barcode JSON file.
    """
    setup_logging()
    r1_path = f"demultiplexed/{lib_barcode}/test_R1.fastq.gz"
    r2_path = f"demultiplexed/{lib_barcode}/test_R2.fastq.gz"
    output_dir = f"int-demultiplexed/{lib_barcode}"
    task_int_demultiplex_batch_cpp(
        r1_path=r1_path,
        r2_path=r2_path,
        output_dir=output_dir,
        barcode_json_path=barcode_json,
    )


@click.command("int-demultiplex-rna-map")
@time_it
@click.argument("code")
@click.argument("lib_barcode_seq")
@click.argument("construct_barcode_seq")
@click.option(
    "--params-file",
    type=click.Path(exists=True),
    default=None,
    help="Path to rna-map parameters YAML file. Uses bundled defaults if not specified.",
)
def int_demultiplex_rna_map(
    code: str,
    lib_barcode_seq: str,
    construct_barcode_seq: str,
    params_file: str | None,
) -> None:
    """Run RNA mapping on internally demultiplexed reads.

    Arguments:
        code: Construct code.
        lib_barcode_seq: Library barcode sequence.
        construct_barcode_seq: Internal construct barcode sequence.
        params_file: Optional path to rna-map parameters file.
    """
    setup_logging()
    task_int_demultiplex_rna_map(code, lib_barcode_seq, construct_barcode_seq, params_file=params_file)


@click.command("int-demultiplex-rna-map-combine")
@time_it
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
