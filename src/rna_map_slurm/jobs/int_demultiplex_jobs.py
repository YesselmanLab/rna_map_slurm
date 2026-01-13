"""Internal demultiplexing job generation."""

from __future__ import annotations

import os
from typing import Any

import pandas as pd

from rna_map_slurm.jobs.generator import (
    create_job_header,
    ensure_job_directories,
    generate_job_list,
    group_into_batches,
    write_job_file,
)
from rna_map_slurm.jobs.slurm import generate_submit_file
from rna_map_slurm.utils.logging import get_logger

log = get_logger("jobs.int_demultiplex")


def generate_int_demultiplex_jobs(
    params: dict[str, Any],
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Generate SLURM jobs for internal demultiplexing.

    Args:
        params: Workflow parameters.
        df: DataFrame with construct information.

    Returns:
        DataFrame with job details.
    """
    job_name = "int-demultiplex"
    job_dir = ensure_job_directories(job_name)
    os.makedirs("int-demultiplexed", exist_ok=True)

    runs_per_job = params["tasks_per_job"][job_name]
    slurm_params = params["slurm_options"][job_name]
    extra_cmds = params["slurm_options"].get("extra-header-cmds", "")

    df_barcodes = _prepare_barcode_dataframe(df)
    barcode_groups = group_into_batches(
        df_barcodes.to_dict("records"),
        runs_per_job,
    )

    job_names: list[str] = []
    for i, group in enumerate(barcode_groups):
        name = f"{job_name}-{i:04}"
        header = create_job_header(name, slurm_params, extra_cmds, str(job_dir))
        body = _build_int_demultiplex_job_body(group)
        write_job_file(job_dir, name, header + body)
        job_names.append(name)

    df_jobs = generate_job_list(job_dir, job_name, "join-fastq-files", job_names)
    generate_submit_file(f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist())
    return df_jobs


def _prepare_barcode_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Prepare combined barcode DataFrame from all constructs.

    Args:
        df: DataFrame with construct information.

    Returns:
        Combined DataFrame with barcode information.
    """
    dfs: list[pd.DataFrame] = []

    for _, row in df.iterrows():
        barcode_seq = row["barcode_seq"]
        os.makedirs(f"int-demultiplexed/{barcode_seq}", exist_ok=True)

        barcode_path = f"inputs/barcode_jsons/{row['code']}.json"
        df_barcodes = pd.read_json(barcode_path)
        df_barcodes["construct_barcode"] = barcode_seq
        dfs.append(df_barcodes)

    return pd.concat(dfs, ignore_index=True)


def _build_int_demultiplex_job_body(records: list[dict[str, Any]]) -> str:
    """Build internal demultiplexing job body.

    Args:
        records: List of barcode records.

    Returns:
        Job body as string.
    """
    lines: list[str] = []
    seen_barcodes: set[str] = set()

    for record in records:
        full_barcode = record["full_barcode"]
        if full_barcode in seen_barcodes:
            continue
        seen_barcodes.add(full_barcode)

        cmd = _build_demultiplex_command(record)
        lines.append(cmd)
        lines.append("")

    return "\n".join(lines)


def _build_demultiplex_command(record: dict[str, Any]) -> str:
    """Build single demultiplex command from record.

    Args:
        record: Barcode record dictionary.

    Returns:
        Command string.
    """
    bb1 = record["barcode_bounds"][0][0]
    bb2 = record["barcode_bounds"][0][1]

    # Map barcode location to the other read
    end_len = len(record["sequence"])
    max_len = end_len - bb2[0]
    min_len = end_len - bb2[1]
    bb2_mapped = [min_len, max_len]

    # Convert U to T for DNA sequences (FASTQ uses DNA, not RNA)
    b1_seq = record["barcodes"][0][0].replace("U", "T")
    b2_seq = record["barcodes"][0][1].replace("U", "T")

    return (
        f"rna-map-slurm-runner int-demultiplex {record['construct_barcode']} "
        f"{b1_seq} {b2_seq} {bb1[0]} {bb1[1]} {bb2_mapped[0]} {bb2_mapped[1]}"
    )


def generate_int_demultiplex_rna_map_jobs(
    params: dict[str, Any],
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Generate SLURM jobs for RNA mapping of internally demultiplexed reads.

    Args:
        params: Workflow parameters.
        df: DataFrame with construct information.

    Returns:
        DataFrame with job details.
    """
    job_name = "int-demultiplex-rna-map"
    os.makedirs("int-demultiplexed-rna-map", exist_ok=True)
    job_dir = ensure_job_directories(job_name)

    runs_per_job = params["tasks_per_job"][job_name]
    slurm_params = params["slurm_options"][job_name]
    extra_cmds = params["slurm_options"].get("extra-header-cmds", "")
    rna_map_params_file = params.get("rna_map_params_file")

    runs = _collect_int_demultiplex_runs(df)
    run_groups = group_into_batches(runs, runs_per_job)
    job_names: list[str] = []

    for i, group in enumerate(run_groups):
        name = f"{job_name}-{i:04}"
        header = create_job_header(name, slurm_params, extra_cmds, str(job_dir))
        body = _build_int_rna_map_job_body(group, rna_map_params_file)
        write_job_file(job_dir, name, header + body)
        job_names.append(name)

    df_jobs = generate_job_list(job_dir, job_name, "int-demultiplex", job_names)
    generate_submit_file(f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist())
    return df_jobs


def _collect_int_demultiplex_runs(df: pd.DataFrame) -> list[tuple[str, str, str]]:
    """Collect all internal demultiplexing run combinations.

    Args:
        df: DataFrame with construct information.

    Returns:
        List of (code, barcode_seq, full_barcode) tuples.
    """
    runs: list[tuple[str, str, str]] = []

    for _, row in df.iterrows():
        barcode_seq = row["barcode_seq"]
        os.makedirs(f"int-demultiplexed-rna-map/{barcode_seq}", exist_ok=True)

        barcode_path = f"inputs/barcode_jsons/{row['code']}.json"
        df_barcode = pd.read_json(barcode_path)

        for barcode in df_barcode["full_barcode"].unique():
            runs.append((row["code"], barcode_seq, barcode))

    return runs


def _build_int_rna_map_job_body(
    runs: list[tuple[str, str, str]],
    params_file: str | None = None,
) -> str:
    """Build internal RNA-map job body.

    Args:
        runs: List of (code, barcode_seq, full_barcode) tuples.
        params_file: Optional path to rna-map parameters file.

    Returns:
        Job body as string.
    """
    params_opt = ""
    if params_file is not None:
        params_opt = f" --params-file {params_file}"

    lines: list[str] = []
    for code, barcode_seq, full_barcode in runs:
        # Convert U to T to match file names from int-demultiplex
        full_barcode_dna = full_barcode.replace("U", "T")
        lines.append(f"rna-map-slurm-runner int-demultiplex-rna-map {code} {barcode_seq} {full_barcode_dna}{params_opt}")
        lines.append("")
    return "\n".join(lines)


def generate_int_demultiplex_rna_map_combine_jobs(
    params: dict[str, Any],
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Generate SLURM jobs for combining internal demultiplexing RNA-map results.

    Args:
        params: Workflow parameters.
        df: DataFrame with construct information.

    Returns:
        DataFrame with job details.
    """
    job_name = "int-demultiplex-rna-map-combine"
    job_dir = ensure_job_directories(job_name)

    slurm_params = params["slurm_options"][job_name]
    extra_cmds = params["slurm_options"].get("extra-header-cmds", "")
    job_names: list[str] = []

    for i, row in df.iterrows():
        name = f"{job_name}-{i:04}"
        header = create_job_header(name, slurm_params, extra_cmds, str(job_dir))
        body = f"rna-map-slurm-runner int-demultiplex-rna-map-combine {row['barcode_seq']} {row['construct']}\n"
        write_job_file(job_dir, name, header + body)
        job_names.append(name)

    df_jobs = generate_job_list(job_dir, job_name, "int-demultiplex-rna-map", job_names)
    generate_submit_file(f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist())
    return df_jobs
