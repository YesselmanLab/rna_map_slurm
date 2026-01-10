"""Run command for rna-map-slurm CLI."""

from __future__ import annotations

import os
from pathlib import Path

import click
import pandas as pd

from rna_map_slurm.jobs.checker import CheckReport, JobChecker
from rna_map_slurm.jobs.executor import JobExecutor, SubmitResult
from rna_map_slurm.jobs.pre_validator import PreRunValidator
from rna_map_slurm.jobs.slurm import submit_job
from rna_map_slurm.models.config import SlurmOptions
from rna_map_slurm.utils.logging import get_logger, setup_logging
from rna_map_slurm.utils.timing import time_it

log = get_logger("cli.run")


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
    "--validate/--no-validate",
    default=True,
    help="Run pre-submission validation checks.",
)
@click.option(
    "--use-arrays/--no-arrays",
    default=False,
    help="Submit jobs as SLURM arrays for faster submission.",
)
@click.option(
    "--max-concurrent",
    default=50,
    type=int,
    help="Maximum concurrent array tasks (only with --use-arrays).",
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
def run(
    dry_run: bool,
    verify: bool,
    validate: bool,
    use_arrays: bool,
    max_concurrent: int,
    job_dir: str,
    jobs_csv: str,
) -> None:
    """Submit SLURM jobs and optionally verify completion.

    Submits all jobs defined in jobs.csv to the SLURM scheduler.
    Supports both individual job submission and SLURM job arrays.

    \b
    Examples:
        rna-map-slurm run                    # Submit all jobs
        rna-map-slurm run --dry-run          # Show what would be submitted
        rna-map-slurm run --use-arrays       # Submit as job arrays (faster)
        rna-map-slurm run --verify           # Submit and check outputs
        rna-map-slurm run --no-validate      # Skip pre-run validation
    """
    _initialize_run()

    jobs_csv_path = Path(jobs_csv)
    job_dir_path = Path(job_dir)

    # Run pre-submission validation
    if validate:
        log.info("Running pre-submission validation...")
        validator = PreRunValidator()
        report = validator.validate_all(jobs_csv=jobs_csv_path)

        for error in report.errors:
            log.error(f"  {error}")
        for warning in report.warnings:
            log.warning(f"  {warning}")

        if not report.valid:
            log.error("Pre-run validation failed. Fix errors before submitting.")
            raise SystemExit(1)

        log.info(f"Validation passed: {report.checks_passed} checks OK")

    # Load jobs
    if not jobs_csv_path.exists():
        log.error(f"Jobs file not found: {jobs_csv}")
        click.echo(f"Error: {jobs_csv} not found. Run 'rna-map-slurm setup' first.")
        raise SystemExit(1)

    df = _load_jobs(jobs_csv_path)
    _log_job_summary(df)

    if dry_run:
        _print_dry_run_summary(df, use_arrays)
        return

    # Submit jobs
    if use_arrays:
        result = _submit_as_arrays(df, job_dir_path, max_concurrent)
    else:
        result = _submit_individual_jobs(df)

    _log_submit_result(result)

    # Save submitted job IDs for reference
    if result.job_ids:
        _save_submitted_jobs(result, jobs_csv_path.parent)

    # Track if any errors occurred
    has_errors = False

    if result.failed > 0:
        log.error(f"{result.failed} jobs failed to submit")
        has_errors = True

    # Verify job outputs if requested
    if verify:
        click.echo("\nVerifying job outputs...")
        check_report = _verify_jobs(job_dir_path, jobs_csv_path)
        _log_verification_result(check_report)

        if check_report.failure_count > 0 or check_report.missing_count > 0:
            has_errors = True

    if has_errors:
        log.error("Completed with errors")
        raise SystemExit(1)

    log.info("Done")


def _initialize_run() -> None:
    """Initialize run logging."""
    os.makedirs("logs", exist_ok=True)
    if os.path.isfile("logs/run.log"):
        os.remove("logs/run.log")
    setup_logging(file_name="logs/run.log")


def _load_jobs(jobs_csv: Path) -> pd.DataFrame:
    """Load jobs DataFrame."""
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


def _print_dry_run_summary(df: pd.DataFrame, use_arrays: bool) -> None:
    """Print what would be submitted in dry-run mode."""
    mode = "as SLURM arrays" if use_arrays else "individually"
    click.echo(f"\n=== DRY RUN - Jobs would be submitted {mode} ===\n")

    for job_type in df["job_type"].unique():
        group = df[df["job_type"] == job_type]
        if use_arrays:
            click.echo(f"{job_type}: {len(group)} jobs (1 array job)")
        else:
            click.echo(f"{job_type}: {len(group)} jobs")

        for _, row in group.head(3).iterrows():
            click.echo(f"  - {row['job_path']}")
        if len(group) > 3:
            click.echo(f"  ... and {len(group) - 3} more")
        click.echo()

    if use_arrays:
        num_arrays = len(df["job_type"].unique())
        click.echo(f"Total: {len(df)} jobs in {num_arrays} array(s)")
    else:
        click.echo(f"Total: {len(df)} individual jobs")


def _submit_individual_jobs(df: pd.DataFrame) -> SubmitResult:
    """Submit jobs individually using sbatch."""
    import time

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


def _submit_as_arrays(
    df: pd.DataFrame,
    job_dir: Path,
    max_concurrent: int,
) -> SubmitResult:
    """Submit jobs as SLURM arrays."""
    executor = JobExecutor(max_concurrent=max_concurrent)

    # Build slurm options map (use defaults for now)
    slurm_options_map: dict[str, SlurmOptions] = {}
    for job_type in df["job_type"].unique():
        slurm_options_map[job_type] = SlurmOptions(name=job_type)

    total_result = SubmitResult()

    for job_type, group in df.groupby("job_type"):
        scripts = [Path(p) for p in group["job_path"]]
        options = slurm_options_map.get(str(job_type), SlurmOptions(name=str(job_type)))

        type_job_dir = job_dir / str(job_type)
        type_job_dir.mkdir(parents=True, exist_ok=True)

        result = executor.submit_scripts_as_array(
            job_type=str(job_type),
            scripts=scripts,
            slurm_options=options,
            job_dir=type_job_dir,
        )

        total_result.submitted += result.submitted
        total_result.failed += result.failed
        total_result.job_ids.extend(result.job_ids)
        total_result.failed_jobs.extend(result.failed_jobs)
        total_result.elapsed_seconds += result.elapsed_seconds

    return total_result


def _log_submit_result(result: SubmitResult) -> None:
    """Log submission results."""
    log.info(f"Submission complete in {result.elapsed_seconds:.1f}s")
    log.info(f"  Submitted: {result.submitted}")
    log.info(f"  Failed:    {result.failed}")


def _save_submitted_jobs(result: SubmitResult, output_dir: Path) -> None:
    """Save submitted job IDs to file for reference."""
    output_path = output_dir / "submitted_jobs.txt"
    with open(output_path, "w") as f:
        for job_id in result.job_ids:
            f.write(f"{job_id}\n")
    log.info(f"Saved {len(result.job_ids)} job IDs to {output_path}")


def _verify_jobs(job_dir: Path, jobs_csv: Path) -> CheckReport:
    """Verify job outputs using JobChecker."""
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
