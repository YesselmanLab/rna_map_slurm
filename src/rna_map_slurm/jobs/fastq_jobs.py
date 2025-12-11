"""FASTQ-related job generation (split, demultiplex, join, trim)."""

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
from rna_map_slurm.models.fastq import PairedFastqFiles
from rna_map_slurm.utils.logging import get_logger

log = get_logger("jobs.fastq")


def generate_split_fastq_jobs(
    pfqs: list[PairedFastqFiles],
    params: dict[str, Any],
) -> pd.DataFrame:
    """Generate SLURM jobs for splitting FASTQ files.

    Args:
        pfqs: List of paired FASTQ files.
        params: Workflow parameters.

    Returns:
        DataFrame with job details.
    """
    job_name = "split-fastq"
    job_dir = ensure_job_directories(job_name, params["num_dirs"])

    slurm_params = params["slurm_options"][job_name]
    extra_cmds = params["slurm_options"].get("extra-header-cmds", "")
    cur_dir = os.path.abspath(os.getcwd())
    job_names: list[str] = []

    for i, pfq in enumerate(pfqs):
        name = f"{job_name}-{i:04}"
        header = create_job_header(name, slurm_params, extra_cmds, job_dir)
        body = _build_split_job_body(pfq, cur_dir, params, i)
        write_job_file(job_dir, name, header + body)
        job_names.append(name)

    df_jobs = generate_job_list(job_dir, job_name, "", job_names)
    generate_submit_file(f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist())
    return df_jobs


def _build_split_job_body(
    pfq: PairedFastqFiles,
    cur_dir: str,
    params: dict[str, Any],
    index: int,
) -> str:
    """Build the body of a split FASTQ job script.

    Args:
        pfq: Paired FASTQ files to split.
        cur_dir: Current working directory.
        params: Workflow parameters.
        index: Index of this FASTQ pair.

    Returns:
        Job body as string.
    """
    start = index * params["fastq_chunks"]
    return (
        f"rna-map-slurm-runner split-fastqs {pfq.read_1.path} {pfq.read_2.path} "
        f"{os.path.join(cur_dir, 'data')} {params['fastq_chunks']} --start {start}\n"
    )


def generate_demultiplexing_jobs(params: dict[str, Any]) -> pd.DataFrame:
    """Generate SLURM jobs for demultiplexing.

    Args:
        params: Workflow parameters.

    Returns:
        DataFrame with job details.
    """
    job_name = "demultiplex"
    job_dir = ensure_job_directories(job_name)

    runs_per_job = params["tasks_per_job"][job_name]
    slurm_params = params["slurm_options"][job_name]
    extra_cmds = params["slurm_options"].get("extra-header-cmds", "")

    dirs = [os.path.abspath(f"data/split-{i:04}") for i in range(params["num_dirs"])]
    dir_groups = group_into_batches(dirs, runs_per_job)
    csv_path = os.path.abspath("data.csv")
    job_names: list[str] = []

    for i, dir_group in enumerate(dir_groups):
        name = f"{job_name}-{i:04}"
        header = create_job_header(name, slurm_params, extra_cmds, job_dir)
        body = _build_demultiplex_job_body(dir_group, csv_path)
        write_job_file(job_dir, name, header + body)
        job_names.append(name)

    df_jobs = generate_job_list(job_dir, job_name, "split-fastq", job_names)
    generate_submit_file(f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist())
    return df_jobs


def _build_demultiplex_job_body(dirs: list[str], csv_path: str) -> str:
    """Build demultiplexing job body.

    Args:
        dirs: List of data directories.
        csv_path: Path to CSV file with barcode info.

    Returns:
        Job body as string.
    """
    lines = [
        f"rna-map-slurm-runner demultiplex {csv_path} {d}/test_R1.fastq.gz {d}/test_R2.fastq.gz {d}"
        for d in dirs
    ]
    return "\n".join(lines) + "\n"


def generate_join_fastq_files_jobs(params: dict[str, Any]) -> pd.DataFrame:
    """Generate SLURM job for joining FASTQ files.

    Args:
        params: Workflow parameters.

    Returns:
        DataFrame with job details.
    """
    job_name = "join-fastq-files"
    job_dir = ensure_job_directories(job_name)

    slurm_params = params["slurm_options"][job_name]
    extra_cmds = params["slurm_options"].get("extra-header-cmds", "")

    name = f"{job_name}-0000"
    header = create_job_header(name, slurm_params, extra_cmds, job_dir)
    body = "rna-map-slurm-runner join-fastq-files\n"
    write_job_file(job_dir, name, header + body)

    return generate_job_list(job_dir, job_name, "demultiplex", [name])


def generate_trim_galore_jobs(params: dict[str, Any]) -> pd.DataFrame:
    """Generate SLURM jobs for trim_galore processing.

    Args:
        params: Workflow parameters.

    Returns:
        DataFrame with job details.
    """
    job_name = "trim-galore"
    job_dir = ensure_job_directories(job_name)

    runs_per_job = params["tasks_per_job"][job_name]
    slurm_params = params["slurm_options"][job_name]
    extra_cmds = params["slurm_options"].get("extra-header-cmds", "")

    dirs = [os.path.abspath(f"data/split-{i:04}") for i in range(params["num_dirs"])]
    dir_groups = group_into_batches(dirs, runs_per_job)
    cur_dir = os.path.abspath(os.getcwd())
    job_names: list[str] = []

    for i, dir_group in enumerate(dir_groups):
        name = f"{job_name}-{i:04}"
        header = create_job_header(name, slurm_params, extra_cmds, job_dir)
        body = _build_trim_galore_job_body(dir_group, cur_dir)
        write_job_file(job_dir, name, header + body)
        job_names.append(name)

    df_jobs = generate_job_list(job_dir, job_name, "split-fastq", job_names)
    generate_submit_file(f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist())
    return df_jobs


def _build_trim_galore_job_body(dirs: list[str], cur_dir: str) -> str:
    """Build trim_galore job body.

    Args:
        dirs: List of data directories.
        cur_dir: Current working directory.

    Returns:
        Job body as string.
    """
    lines: list[str] = []
    for d in dirs:
        lines.extend([
            f"cd {d}",
            f"trim_galore --quality 0 --paired {d}/test_R1.fastq.gz {d}/test_R2.fastq.gz",
            f"cd {cur_dir}",
            "",
        ])
    return "\n".join(lines)
