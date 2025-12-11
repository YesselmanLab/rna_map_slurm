"""Base job generation utilities."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import pandas as pd

from rna_map_slurm.jobs.slurm import get_job_header
from rna_map_slurm.models.config import SlurmOptions
from rna_map_slurm.utils.logging import get_logger

log = get_logger("jobs.generator")


def create_job_header(
    name: str,
    slurm_params: dict[str, Any],
    extra_header_cmds: str,
    path: str | Path,
) -> str:
    """Create a SLURM job header from parameters.

    Args:
        name: Job name.
        slurm_params: Dictionary with time, mem-per-cpu, cpus-per-task.
        extra_header_cmds: Additional SBATCH directives.
        path: Directory path for job output files.

    Returns:
        Formatted SLURM job header string.
    """
    slurm_opts = SlurmOptions(
        name=name,
        time=slurm_params["time"],
        mem_per_cpu=slurm_params["mem-per-cpu"],
        cpus_per_task=slurm_params["cpus-per-task"],
        extra_header_cmds=extra_header_cmds,
    )
    return get_job_header(slurm_opts, os.path.abspath(str(path)))


def write_job_file(path: str | Path, job_name: str, job_content: str) -> None:
    """Write job content to a shell script file.

    Args:
        path: Directory to create the job file in.
        job_name: Name of the job (used for filename).
        job_content: Full job script content.
    """
    job_path = Path(path) / f"{job_name}.sh"
    with open(job_path, "w", encoding="utf-8") as f:
        f.write(job_content)


def generate_job_list(
    path: str | Path,
    job_type: str,
    requirement: str,
    job_names: list[str],
) -> pd.DataFrame:
    """Generate a DataFrame listing job details.

    Args:
        path: Directory containing job files.
        job_type: Type identifier for these jobs.
        requirement: Job dependency (job type that must complete first).
        job_names: List of job names.

    Returns:
        DataFrame with columns: job_path, job_type, job_requirement.
    """
    jobs = [[f"{path}/{name}.sh", job_type, requirement] for name in job_names]
    return pd.DataFrame(jobs, columns=["job_path", "job_type", "job_requirement"])


def group_into_batches(items: list[Any], batch_size: int) -> list[list[Any]]:
    """Split a list into batches of specified size.

    Args:
        items: List to split.
        batch_size: Maximum items per batch.

    Returns:
        List of batches.
    """
    return [items[i : i + batch_size] for i in range(0, len(items), batch_size)]


def ensure_job_directories(job_name: str, num_dirs: int | None = None) -> Path:
    """Create job directory structure.

    Args:
        job_name: Name of the job type.
        num_dirs: Optional number of split directories to create.

    Returns:
        Path to the job directory.
    """
    job_dir = Path("jobs") / job_name
    os.makedirs(job_dir, exist_ok=True)

    if num_dirs is not None:
        for i in range(num_dirs):
            os.makedirs(f"data/split-{i:04}", exist_ok=True)

    return job_dir
