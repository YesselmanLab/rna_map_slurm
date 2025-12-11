"""RNA-map job generation with large construct handling."""

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

log = get_logger("jobs.rna_map")

DEFAULT_LARGE_CONSTRUCT_THRESHOLD = 100


def count_sequences_in_csv(csv_path: str) -> int:
    """Count the number of sequences in a CSV file.

    Args:
        csv_path: Path to CSV file.

    Returns:
        Number of rows (sequences) in the CSV.
    """
    if not os.path.isfile(csv_path):
        return 0
    df = pd.read_csv(csv_path)
    return len(df)


def is_large_construct(csv_path: str, threshold: int) -> bool:
    """Check if a construct has more sequences than the threshold.

    Args:
        csv_path: Path to the construct CSV file.
        threshold: Maximum sequences before splitting.

    Returns:
        True if construct exceeds threshold.
    """
    return count_sequences_in_csv(csv_path) > threshold


def generate_rna_map_jobs(
    params: dict[str, Any],
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Generate SLURM jobs for RNA mapping.

    Large constructs (>threshold sequences) get their own dedicated jobs
    instead of being batched together.

    Args:
        params: Workflow parameters including construct_options.
        df: DataFrame with construct information.

    Returns:
        DataFrame with job details.
    """
    job_name = "rna-map"
    job_dir = ensure_job_directories(job_name)

    runs_per_job = params["tasks_per_job"][job_name]
    slurm_params = params["slurm_options"][job_name]
    extra_cmds = params["slurm_options"].get("extra-header-cmds", "")

    threshold = params.get("construct_options", {}).get(
        "large_construct_threshold",
        DEFAULT_LARGE_CONSTRUCT_THRESHOLD,
    )

    dirs = [f"data/split-{i:04}" for i in range(params["num_dirs"])]
    dir_groups = group_into_batches(dirs, runs_per_job)
    job_names: list[str] = []
    job_index = 0

    for _, row in df.iterrows():
        code = row["code"]
        fa_path = f"inputs/fastas/{code}.fasta"

        if not os.path.isfile(fa_path):
            log.warning(f"{code} does not have a FASTA file")
            continue

        csv_path = f"inputs/rnas/{code}.csv"
        construct_is_large = is_large_construct(csv_path, threshold)

        if construct_is_large:
            log.info(f"Large construct detected: {code} - creating dedicated jobs")

        job_index, new_names = _generate_jobs_for_construct(
            row=row,
            dir_groups=dir_groups,
            job_name=job_name,
            job_dir=str(job_dir),
            slurm_params=slurm_params,
            extra_cmds=extra_cmds,
            job_index=job_index,
            is_large=construct_is_large,
        )
        job_names.extend(new_names)

    df_jobs = generate_job_list(job_dir, job_name, "demultiplex", job_names)
    generate_submit_file(f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist())
    return df_jobs


def _generate_jobs_for_construct(
    row: pd.Series[Any],
    dir_groups: list[list[str]],
    job_name: str,
    job_dir: str,
    slurm_params: dict[str, Any],
    extra_cmds: str,
    job_index: int,
    is_large: bool,
) -> tuple[int, list[str]]:
    """Generate jobs for a single construct.

    Args:
        row: DataFrame row with construct info.
        dir_groups: Batched data directories.
        job_name: Base job name.
        job_dir: Job directory.
        slurm_params: SLURM parameters.
        extra_cmds: Extra SBATCH commands.
        job_index: Current job index.
        is_large: Whether this is a large construct.

    Returns:
        Tuple of (new job index, list of job names created).
    """
    job_names: list[str] = []

    for dir_group in dir_groups:
        if is_large:
            # Large constructs: one directory per job
            for single_dir in dir_group:
                name = f"{job_name}-{job_index:04}"
                header = create_job_header(name, slurm_params, extra_cmds, job_dir)
                body = _build_rna_map_job_body(row, [single_dir])
                write_job_file(job_dir, name, header + body)
                job_names.append(name)
                job_index += 1
        else:
            # Normal batching
            name = f"{job_name}-{job_index:04}"
            header = create_job_header(name, slurm_params, extra_cmds, job_dir)
            body = _build_rna_map_job_body(row, dir_group)
            write_job_file(job_dir, name, header + body)
            job_names.append(name)
            job_index += 1

    return job_index, job_names


def _build_rna_map_job_body(
    row: pd.Series[Any],
    dirs: list[str],
) -> str:
    """Build RNA-map job body.

    Args:
        row: DataFrame row with construct info.
        dirs: Data directories to process.

    Returns:
        Job body as string.
    """
    code = row["code"]
    barcode_seq = row["barcode_seq"]
    construct = row["construct"]

    fa_path = os.path.abspath(f"inputs/fastas/{code}.fasta")
    csv_path = os.path.abspath(f"inputs/rnas/{code}.csv")

    lines: list[str] = []
    for data_dir in dirs:
        output_dir = os.path.abspath(f"{data_dir}/{barcode_seq}/{construct}")
        os.makedirs(output_dir, exist_ok=True)
        fq1_path = os.path.abspath(f"{data_dir}/{barcode_seq}/test_R1.fastq.gz")
        fq2_path = os.path.abspath(f"{data_dir}/{barcode_seq}/test_R2.fastq.gz")
        lines.append(
            f"rna-map-slurm-runner run-rna-map {fa_path} {fq2_path} {fq1_path} {csv_path} {output_dir}"
        )
        lines.append("")

    return "\n".join(lines)


def generate_rna_map_combine_jobs(
    params: dict[str, Any],
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Generate SLURM jobs for combining RNA-map results.

    Args:
        params: Workflow parameters.
        df: DataFrame with construct information.

    Returns:
        DataFrame with job details.
    """
    job_name = "rna-map-combine"
    job_dir = ensure_job_directories(job_name)

    slurm_params = params["slurm_options"][job_name]
    extra_cmds = params["slurm_options"].get("extra-header-cmds", "")
    job_names: list[str] = []

    for i, row in df.iterrows():
        name = f"{job_name}-{i:04}"
        header = create_job_header(name, slurm_params, extra_cmds, str(job_dir))
        body = f"rna-map-slurm-runner rna-map-combine {row['barcode_seq']} {row['construct']}\n"
        write_job_file(job_dir, name, header + body)
        job_names.append(name)

    df_jobs = generate_job_list(job_dir, job_name, "rna-map", job_names)
    generate_submit_file(f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist())
    return df_jobs
