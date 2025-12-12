"""Basic task implementations for the RNA mapping workflow."""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
from typing import Any

import pandas as pd
import rna_map.run
from fastqsplitter import split_fastqs as fastqsplitter
from rna_map.mutation_histogram import (
    get_mut_histos_from_pickle_file,
    merge_mut_histo_dicts,
    write_mut_histos_to_pickle_file,
)
from rna_map.parameters import get_preset_params, parse_parameters_from_file

from rna_map_slurm.io.demultiplex import SabreDemultiplexer
from rna_map_slurm.models.fastq import FastqFile, PairedFastqFiles
from rna_map_slurm.plotting.pop_avg import generate_pop_avg_plots
from rna_map_slurm.tasks.rna_map_helpers import get_mut_histo_dataframe
from rna_map_slurm.utils.logging import get_logger

log = get_logger("tasks.basic")


def split_fastq_file(
    fastq_file: str,
    output_dir: str,
    num_chunks: int,
    start: int,
    threads: int = 1,
) -> list[str]:
    """Split a FASTQ file into multiple chunks.

    Args:
        fastq_file: Path to input FASTQ file.
        output_dir: Directory for output files.
        num_chunks: Number of chunks to create.
        start: Starting index for chunk numbering.
        threads: Number of threads for splitting.

    Returns:
        List of output file paths.
    """
    log.info(f"Splitting {fastq_file} into {num_chunks} chunks")
    log.info(f"output_dir: {output_dir}, num_chunks: {num_chunks}, start: {start}")

    output_file = _get_output_filename(fastq_file)
    output_files = _build_output_paths(output_dir, output_file, num_chunks, start)

    fastqsplitter(fastq_file, output_files, threads_per_file=threads)
    return output_files


def _get_output_filename(fastq_file: str) -> str:
    """Determine output filename based on read type.

    Args:
        fastq_file: Input FASTQ file path.

    Returns:
        Output filename (test_R1.fastq.gz or test_R2.fastq.gz).
    """
    if "R2" in fastq_file:
        log.info("Detected R2 file")
        return "test_R2.fastq.gz"
    log.info("Detected R1 file")
    return "test_R1.fastq.gz"


def _build_output_paths(
    output_dir: str,
    output_file: str,
    num_chunks: int,
    start: int,
) -> list[str]:
    """Build list of output file paths.

    Args:
        output_dir: Base output directory.
        output_file: Output filename.
        num_chunks: Number of chunks.
        start: Starting index.

    Returns:
        List of output paths.
    """
    return [
        f"{output_dir}/split-{i:04}/{output_file}"
        for i in range(start, num_chunks + start)
    ]


def demultiplex(
    csv: str,
    r1_path: str,
    r2_path: str,
    output_dir: str | None = None,
) -> list[str]:
    """Demultiplex paired FASTQ files by 3' barcodes.

    Args:
        csv: Path to CSV with barcode information.
        r1_path: Path to R1 FASTQ file.
        r2_path: Path to R2 FASTQ file.
        output_dir: Output directory (defaults to current directory).

    Returns:
        List of output directories for each barcode.
    """
    if output_dir is None:
        output_dir = os.getcwd()

    cur_dir = os.getcwd()
    os.chdir(output_dir)

    try:
        paired_fastqs = PairedFastqFiles(FastqFile(r1_path), FastqFile(r2_path))
        df = pd.read_csv(csv)
        demultiplexer = SabreDemultiplexer()
        demultiplexer.run(df, paired_fastqs, output_dir)
    finally:
        os.chdir(cur_dir)

    return [os.path.join(output_dir, row["barcode_seq"]) for _, row in df.iterrows()]


def join_fastq_files(fastq_files: list[str], joined_fastq: str) -> str:
    """Concatenate multiple FASTQ files into one.

    Args:
        fastq_files: List of FASTQ file paths.
        joined_fastq: Output path for joined file.

    Returns:
        Path to joined FASTQ file.
    """
    log.info(f"Joining FASTQ files into {joined_fastq}")
    files_str = " ".join(fastq_files)

    subprocess.run(f"cat {files_str} > {joined_fastq}", shell=True, check=True)

    log.info(f"FASTQ files joined: {joined_fastq}")
    return joined_fastq


def run_rna_map(
    fa_path: str,
    r1_path: str,
    r2_path: str,
    csv_path: str,
    output_dir: str,
) -> None:
    """Run RNA mapping on FASTQ files.

    Args:
        fa_path: Path to FASTA reference.
        r1_path: Path to R1 FASTQ file.
        r2_path: Path to R2 FASTQ file.
        csv_path: Path to CSV with sequence info.
        output_dir: Output directory.
    """
    cur_dir = os.getcwd()
    log.info(f"Changing directory to {output_dir}")
    os.chdir(output_dir)

    try:
        fa_path, csv_path, params = _prepare_rna_map_inputs(fa_path, csv_path)
        _execute_rna_map(fa_path, r1_path, r2_path, csv_path, params)
        _cleanup_rna_map_outputs()
    finally:
        os.chdir(cur_dir)


def _prepare_rna_map_inputs(
    fa_path: str,
    csv_path: str,
) -> tuple[str, str, dict[str, Any]]:
    """Prepare inputs for RNA mapping, using local overrides if present.

    Args:
        fa_path: Default FASTA path.
        csv_path: Default CSV path.

    Returns:
        Tuple of (fasta_path, csv_path, params).
    """
    if os.path.isfile("input.fasta"):
        log.info("Using existing input.fasta file")
        fa_path = "input.fasta"

    if os.path.isfile("input.csv"):
        log.info("Using existing input.csv file")
        csv_path = "input.csv"

    if os.path.isfile("params.yml"):
        log.info("Using existing params.yml file")
        params = parse_parameters_from_file("params.yml")
    else:
        params = get_preset_params("barcoded-library")

    params["overwrite"] = True
    params["bit_vector"]["summary_output_only"] = True

    return fa_path, csv_path, params


def _execute_rna_map(
    fa_path: str,
    r1_path: str,
    r2_path: str,
    csv_path: str,
    params: dict[str, Any],
) -> None:
    """Execute RNA mapping.

    Args:
        fa_path: FASTA path.
        r1_path: R1 FASTQ path.
        r2_path: R2 FASTQ path.
        csv_path: CSV path.
        params: RNA-map parameters.
    """
    log.info("Starting RNA mapping")
    rna_map.run.run(fa_path, r1_path, r2_path, csv_path, params)


def _cleanup_rna_map_outputs() -> None:
    """Remove unnecessary RNA-map output files."""
    log.info("Cleaning up unnecessary files: log, input, output/Mapping_Files")
    shutil.rmtree("log", ignore_errors=True)
    shutil.rmtree("input", ignore_errors=True)
    shutil.rmtree("output/Mapping_Files", ignore_errors=True)


def rna_map_combine(row: pd.Series[Any]) -> None:
    """Combine RNA mapping results from multiple chunks.

    Args:
        row: DataFrame row with construct information.
    """
    final_path = _setup_combine_paths(row)
    merged_mut_histos = _merge_mutation_histograms(row)

    if not merged_mut_histos:
        log.warning("No mutation histogram files found to merge")
        return

    _save_combined_results(row, final_path, merged_mut_histos)


def _setup_combine_paths(row: pd.Series[Any]) -> str:
    """Set up output paths for combined results.

    Args:
        row: DataFrame row with construct info.

    Returns:
        Path for output files.
    """
    run_path = f"results/{row['run_name']}"
    dir_name = f"{row['construct']}_{row['code']}_{row['data_type']}"
    final_path = f"{run_path}/processed/{dir_name}/output/BitVector_Files/"

    log.info(f"results path: {final_path}")
    os.makedirs(run_path, exist_ok=True)
    os.makedirs(final_path, exist_ok=True)

    return final_path


def _merge_mutation_histograms(row: pd.Series[Any]) -> dict[str, Any]:
    """Merge mutation histograms from all data chunks.

    Args:
        row: DataFrame row with barcode info.

    Returns:
        Merged mutation histogram dictionary.
    """
    barcode_seq = row["barcode_seq"]
    construct = row["construct"]
    dirs = glob.glob("data/split-*")
    merged_mut_histos: dict[str, Any] = {}
    count_files = 0

    for d in dirs:
        mhs_path = f"{d}/{barcode_seq}/{construct}/output/BitVector_Files/mutation_histos.p"
        if not os.path.isfile(mhs_path):
            log.warning(f"files not found: {mhs_path}")
            continue
        merge_mut_histo_dicts(merged_mut_histos, get_mut_histos_from_pickle_file(mhs_path))
        count_files += 1

    log.info(f"merged {count_files} files")
    return merged_mut_histos


def _save_combined_results(
    row: pd.Series[Any],
    final_path: str,
    merged_mut_histos: dict[str, Any],
) -> None:
    """Save combined results to files.

    Args:
        row: DataFrame row with metadata.
        final_path: Output directory path.
        merged_mut_histos: Merged mutation histogram data.
    """
    df_results = get_mut_histo_dataframe(merged_mut_histos)
    df_results = _add_metadata_columns(df_results, row)

    df_results.to_json(f"{final_path}mutation_histos.json", orient="records")

    dir_name = f"{row['construct']}_{row['code']}_{row['data_type']}"
    generate_pop_avg_plots(df_results, row["run_name"], dir_name)

    write_mut_histos_to_pickle_file(merged_mut_histos, f"{final_path}mutation_histos.p")


def _add_metadata_columns(
    df_results: pd.DataFrame,
    row: pd.Series[Any],
) -> pd.DataFrame:
    """Add metadata columns from row to results DataFrame.

    Args:
        df_results: Results DataFrame.
        row: DataFrame row with metadata.

    Returns:
        DataFrame with added metadata columns.
    """
    cols = list(row.keys())
    for col in ["demult_cmd", "length"]:
        if col in cols:
            cols.remove(col)

    for col in cols:
        df_results[col] = row[col]

    return df_results
