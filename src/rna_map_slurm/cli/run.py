"""Run command for rna-map-slurm CLI."""

from __future__ import annotations

import os
import time
from typing import Any

import click
import pandas as pd

from rna_map_slurm.cli.utils import submit_jobs
from rna_map_slurm.jobs.slurm import get_current_user, get_user_jobs, is_job_type_completed
from rna_map_slurm.utils.logging import get_logger, setup_logging
from rna_map_slurm.utils.timing import time_it

log = get_logger("cli.run")

MAX_CONCURRENT_JOBS = 999
POLL_INTERVAL_SECONDS = 60


@click.command()
@time_it
def run() -> None:
    """Run the SLURM workflow, submitting and monitoring jobs."""
    _initialize_run()
    df = _load_jobs()

    completed_types: list[str] = []
    submitted_types = _submit_initial_jobs(df)

    _run_job_loop(df, submitted_types, completed_types)
    log.info("All jobs are completed")


def _initialize_run() -> None:
    """Initialize run logging."""
    if os.path.isfile("logs/run.log"):
        os.remove("logs/run.log")
    setup_logging(file_name="logs/run.log")


def _load_jobs() -> pd.DataFrame:
    """Load jobs DataFrame and initialize status column."""
    df = pd.read_csv("jobs.csv")
    df["status"] = "not_started"
    return df


def _submit_initial_jobs(df: pd.DataFrame) -> list[str]:
    """Submit jobs with no requirements and return their types."""
    df_can_run = df[df["job_requirement"].isna()]
    df.loc[df["job_requirement"].isna(), "status"] = "run"
    submit_jobs(df_can_run)
    return df_can_run["job_type"].tolist()


def _run_job_loop(
    df: pd.DataFrame,
    submitted_types: list[str],
    completed_types: list[str],
) -> None:
    """Main job monitoring and submission loop."""
    user = get_current_user()

    while True:
        time.sleep(POLL_INTERVAL_SECONDS)

        jobs = get_user_jobs(user)
        _log_job_status(jobs, submitted_types, completed_types)

        _update_completed_types(submitted_types, completed_types, jobs)
        df_not_run = _get_pending_jobs(df, completed_types)

        if df_not_run.empty:
            log.info("All jobs are submitted")
            break

        _submit_pending_jobs(df, df_not_run, completed_types, submitted_types, jobs)


def _log_job_status(
    jobs: list[dict[str, Any]],
    submitted_types: list[str],
    completed_types: list[str],
) -> None:
    """Log current job status."""
    log.info(f"num_jobs: {len(jobs)}")
    log.info(f"submitted_types: {submitted_types}")
    log.info(f"completed_types: {completed_types}")


def _update_completed_types(
    submitted_types: list[str],
    completed_types: list[str],
    jobs: list[dict[str, Any]],
) -> None:
    """Update list of completed job types."""
    newly_completed = [
        job_type
        for job_type in submitted_types
        if is_job_type_completed(job_type, jobs)
    ]

    for job_type in newly_completed:
        log.info(f"Job type {job_type} is completed")
        completed_types.append(job_type)
        submitted_types.remove(job_type)


def _get_pending_jobs(df: pd.DataFrame, completed_types: list[str]) -> pd.DataFrame:
    """Get jobs that haven't been run yet and aren't completed."""
    return df[
        (~df["job_type"].isin(completed_types)) & (df["status"] == "not_started")
    ]


def _submit_pending_jobs(
    df: pd.DataFrame,
    df_not_run: pd.DataFrame,
    completed_types: list[str],
    submitted_types: list[str],
    jobs: list[dict[str, Any]],
) -> None:
    """Submit pending jobs whose requirements are met."""
    for job_type, group in df_not_run.groupby("job_type"):
        requirement = group["job_requirement"].iloc[0]
        if requirement not in completed_types:
            continue

        count = _submit_job_batch(df, group, jobs, submitted_types, str(job_type))
        if count > 0:
            log.info(f"Submitted {count} jobs for {job_type}")
            break


def _submit_job_batch(
    df: pd.DataFrame,
    group: pd.DataFrame,
    jobs: list[dict[str, Any]],
    submitted_types: list[str],
    job_type: str,
) -> int:
    """Submit a batch of jobs up to the max limit.

    Returns:
        Number of jobs submitted.
    """
    job_num = len(jobs)
    count = 0

    if job_type not in submitted_types:
        submitted_types.append(job_type)

    log.info(f"Submitting jobs for {job_type}")

    for idx, row in group.iterrows():
        if job_num >= MAX_CONCURRENT_JOBS:
            break
        os.system(f"sbatch {row['job_path']}")
        df.loc[idx, "status"] = "run"  # type: ignore[index]
        job_num += 1
        count += 1

    return count
