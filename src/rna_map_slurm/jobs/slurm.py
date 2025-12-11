"""SLURM job header generation and utilities."""

from __future__ import annotations

import getpass
import os
import subprocess
from pathlib import Path
from typing import Any

from rna_map_slurm.models.config import SlurmOptions
from rna_map_slurm.utils.logging import get_logger

log = get_logger("jobs.slurm")


def get_job_header(args: SlurmOptions, job_dir: str | Path | None = None) -> str:
    """Generate SLURM job script header.

    Args:
        args: SLURM options configuration.
        job_dir: Optional directory for output/error files.

    Returns:
        SLURM batch script header as string.
    """
    if job_dir is not None:
        output_path = f"{job_dir}/{args.name}.out"
        error_path = f"{job_dir}/{args.name}.err"
    else:
        output_path = f"{args.name}.out"
        error_path = f"{args.name}.err"

    return f"""#!/bin/bash
#SBATCH --time={args.time}
#SBATCH --mem-per-cpu={args.mem_per_cpu}
#SBATCH --cpus-per-task={args.cpus_per_task}
#SBATCH --output={output_path}
#SBATCH --error={error_path}

{args.extra_header_cmds}
"""


def generate_submit_file(path: str | Path, jobs: list[str]) -> None:
    """Generate a submit file with sbatch commands for all jobs.

    Args:
        path: Output path for the submit file.
        jobs: List of job script paths.
    """
    with open(str(path), "w", encoding="utf-8") as f:
        for job in jobs:
            f.write(f"sbatch {job}\n")


def get_current_user() -> str:
    """Get current username from environment or system.

    Returns:
        Current username.
    """
    return os.environ.get("USER", getpass.getuser())


def get_user_jobs(user: str | None = None) -> list[dict[str, str]]:
    """Get list of SLURM jobs for a user.

    Args:
        user: Username to query. Defaults to current user.

    Returns:
        List of job dictionaries with job attributes.
    """
    if user is None:
        user = get_current_user()

    try:
        result = subprocess.run(
            [
                "squeue",
                "--user",
                user,
                "--format",
                "%i %P %j %u %t %M %l %D %R",
                "--noheader",
            ],
            capture_output=True,
            text=True,
            check=True,
        )

        lines = result.stdout.strip().split("\n")
        keys = [
            "JobID",
            "Partition",
            "Name",
            "User",
            "State",
            "Time",
            "TimeLimit",
            "Nodes",
            "NodeList",
        ]

        jobs = [dict(zip(keys, line.split(), strict=False)) for line in lines if line.strip()]

        if len(jobs) == 1 and not jobs[0]:
            log.info("No jobs found")
            return []

        return jobs

    except subprocess.CalledProcessError as e:
        log.error(f"Error querying jobs: {e.stderr}")
        return []
    except Exception as e:
        log.error(f"Unexpected error: {e}")
        return []


def is_job_type_completed(job_type: str, jobs: list[dict[str, Any]]) -> bool:
    """Check if all jobs of a specific type have completed.

    Args:
        job_type: Job type substring to match against job names.
        jobs: List of current job dictionaries.

    Returns:
        True if no matching jobs are currently running.
    """
    return not any(job_type in job.get("Name", "") for job in jobs)
