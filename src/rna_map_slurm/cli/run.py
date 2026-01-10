"""Run command for rna-map-slurm CLI."""

from __future__ import annotations

import os
import time
from dataclasses import dataclass, field
from pathlib import Path

import click
import pandas as pd

from rna_map_slurm.jobs.checker import CheckReport, JobChecker
from rna_map_slurm.jobs.slurm import submit_job
from rna_map_slurm.utils.logging import get_logger, setup_logging
from rna_map_slurm.utils.timing import time_it

log = get_logger("cli.run")


@dataclass
class SubmitResult:
    """Result of batch job submission."""

    submitted: int = 0
    failed: int = 0
    job_ids: list[str] = field(default_factory=list)
    failed_jobs: list[str] = field(default_factory=list)
    elapsed_seconds: float = 0.0


@click.command()
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show jobs that would be submitted without actually submitting.",
)
@click.option(
    "--verify",
    is_flag=True,
    help="Check job outputs for errors after submission.",
)
@click.option(
    "--job-dir",
    default="jobs",
    type=click.Path(exists=False),
    help="Directory containing job scripts and outputs.",
)
@click.option(
    "--jobs-csv",
    default="jobs.csv",
    type=click.Path(exists=False),
    help="Path to jobs.csv file listing jobs to submit.",
)
@time_it
def run(dry_run: bool, verify: bool, job_dir: str, jobs_csv: str) -> None:
    """Submit SLURM jobs and optionally verify completion.

    Submits all jobs defined in jobs.csv to the SLURM scheduler.
    Unlike a polling approach, this submits all jobs immediately
    and returns.

    \b
    Examples:
        rna-map-slurm run                    # Submit all jobs
        rna-map-slurm run --dry-run          # Show what would be submitted
        rna-map-slurm run --verify           # Submit and check outputs
    """
    _initialize_run()

    # Load and validate jobs
    jobs_csv_path = Path(jobs_csv)
    if not jobs_csv_path.exists():
        log.error(f"Jobs file not found: {jobs_csv}")
        click.echo(f"Error: {jobs_csv} not found. Run 'rna-map-slurm setup' first.")
        raise SystemExit(1)

    df = _load_jobs(jobs_csv_path)
    _log_job_summary(df)

    # Validate job scripts exist
    missing = _validate_job_scripts(df)
    if missing:
        log.error(f"Missing {len(missing)} job scripts")
        for path in missing[:5]:
            log.error(f"  - {path}")
        if len(missing) > 5:
            log.error(f"  ... and {len(missing) - 5} more")
        raise SystemExit(1)

    if dry_run:
        _print_dry_run_summary(df)
        return

    # Submit all jobs
    result = _submit_all_jobs(df)
    _log_submit_result(result)

    # Save submitted job IDs for reference
    if result.job_ids:
        _save_submitted_jobs(result, jobs_csv_path.parent)

    if result.failed > 0:
        log.warning(f"{result.failed} jobs failed to submit")

    # Verify job outputs if requested
    if verify:
        click.echo("\nVerifying job outputs...")
        report = _verify_jobs(Path(job_dir), jobs_csv_path)
        _log_verification_result(report)

        if report.failure_count > 0:
            raise SystemExit(1)

    log.info("Done")


def _initialize_run() -> None:
    """Initialize run logging."""
    os.makedirs("logs", exist_ok=True)
    if os.path.isfile("logs/run.log"):
        os.remove("logs/run.log")
    setup_logging(file_name="logs/run.log")


def _load_jobs(jobs_csv: Path) -> pd.DataFrame:
    """Load jobs DataFrame.

    Args:
        jobs_csv: Path to jobs.csv file.

    Returns:
        DataFrame with job information.
    """
    df = pd.read_csv(jobs_csv)
    log.info(f"Loaded {len(df)} jobs from {jobs_csv}")
    return df


def _log_job_summary(df: pd.DataFrame) -> None:
    """Log summary of jobs to submit."""
    job_types = df["job_type"].unique().tolist()
    log.info(f"Job types: {job_types}")

    for job_type in job_types:
        count = len(df[df["job_type"] == job_type])
        log.info(f"  {job_type}: {count} jobs")


def _validate_job_scripts(df: pd.DataFrame) -> list[str]:
    """Validate that all job scripts exist.

    Args:
        df: DataFrame with job_path column.

    Returns:
        List of missing job script paths.
    """
    missing = []
    for job_path in df["job_path"]:
        if not Path(job_path).exists():
            missing.append(job_path)
    return missing


def _print_dry_run_summary(df: pd.DataFrame) -> None:
    """Print what would be submitted in dry-run mode."""
    click.echo("\n=== DRY RUN - No jobs will be submitted ===\n")

    for job_type in df["job_type"].unique():
        group = df[df["job_type"] == job_type]
        click.echo(f"{job_type}: {len(group)} jobs")
        for _, row in group.head(3).iterrows():
            click.echo(f"  - {row['job_path']}")
        if len(group) > 3:
            click.echo(f"  ... and {len(group) - 3} more")
        click.echo()

    click.echo(f"Total: {len(df)} jobs would be submitted")


def _submit_all_jobs(df: pd.DataFrame) -> SubmitResult:
    """Submit all jobs to SLURM.

    Args:
        df: DataFrame with job_path column.

    Returns:
        SubmitResult with submission statistics.
    """
    start = time.time()
    job_ids: list[str] = []
    failed_jobs: list[str] = []

    total = len(df)
    submitted_count = 0

    for job_type, group in df.groupby("job_type"):
        log.info(f"Submitting {len(group)} jobs for: {job_type}")

        for _, row in group.iterrows():
            job_path = row["job_path"]
            success, job_id = submit_job(job_path)

            if success:
                job_ids.append(job_id or "unknown")
                submitted_count += 1
                if submitted_count % 50 == 0:
                    log.info(f"  Progress: {submitted_count}/{total} submitted")
            else:
                failed_jobs.append(job_path)
                log.warning(f"  Failed: {job_path}")

    elapsed = time.time() - start

    return SubmitResult(
        submitted=len(job_ids),
        failed=len(failed_jobs),
        job_ids=job_ids,
        failed_jobs=failed_jobs,
        elapsed_seconds=elapsed,
    )


def _log_submit_result(result: SubmitResult) -> None:
    """Log submission results."""
    log.info(f"Submission complete in {result.elapsed_seconds:.1f}s")
    log.info(f"  Submitted: {result.submitted}")
    log.info(f"  Failed:    {result.failed}")


def _save_submitted_jobs(result: SubmitResult, output_dir: Path) -> None:
    """Save submitted job IDs to file for reference.

    Args:
        result: SubmitResult with job IDs.
        output_dir: Directory to save the file.
    """
    output_path = output_dir / "submitted_jobs.txt"
    with open(output_path, "w") as f:
        for job_id in result.job_ids:
            f.write(f"{job_id}\n")
    log.info(f"Saved {len(result.job_ids)} job IDs to {output_path}")


def _verify_jobs(job_dir: Path, jobs_csv: Path) -> CheckReport:
    """Verify job outputs using JobChecker.

    Args:
        job_dir: Directory containing job output files.
        jobs_csv: Path to jobs.csv for expected job list.

    Returns:
        CheckReport with verification results.
    """
    checker = JobChecker(
        job_dir=job_dir,
        jobs_csv=jobs_csv if jobs_csv.exists() else None,
    )
    return checker.check_all()


def _log_verification_result(report: CheckReport) -> None:
    """Log verification results."""
    log.info("Job verification complete:")
    log.info(f"  Succeeded: {report.success_count}/{report.total}")
    log.info(f"  Failed:    {report.failure_count}")
    log.info(f"  Missing:   {report.missing_count}")

    if report.failed:
        log.error("Failed jobs:")
        for job in report.failed:
            log.error(f"  {job.job_name}: {job.reason}")

    if report.missing:
        log.warning("Missing output files:")
        for job in report.missing[:5]:
            log.warning(f"  {job.job_name}")
        if len(report.missing) > 5:
            log.warning(f"  ... and {len(report.missing) - 5} more")
