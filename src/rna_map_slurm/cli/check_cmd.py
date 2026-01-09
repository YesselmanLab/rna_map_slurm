"""CLI command for checking SLURM job outputs."""

from __future__ import annotations

from pathlib import Path

import click
from tabulate import tabulate

from rna_map_slurm.jobs.checker import CheckReport, JobChecker


def _format_report(report: CheckReport, verbose: bool = False) -> str:
    """Format the check report for display.

    Args:
        report: CheckReport from JobChecker
        verbose: Whether to show detailed information

    Returns:
        Formatted string for display
    """
    lines: list[str] = []

    # Summary table
    summary_data = [
        ["Total Jobs", report.total],
        ["Succeeded", report.success_count],
        ["Failed", report.failure_count],
        ["Missing", report.missing_count],
    ]
    lines.append(tabulate(summary_data, headers=["Status", "Count"], tablefmt="simple"))

    # Failed jobs details
    if report.failed:
        lines.append("\n--- Failed Jobs ---")
        failed_data = []
        for job in report.failed:
            row = [job.job_name, job.reason or "Unknown"]
            if verbose and job.error_line:
                # Truncate long error lines
                error_preview = job.error_line[:80] + "..." if len(job.error_line) > 80 else job.error_line
                row.append(error_preview)
            failed_data.append(row)

        headers = ["Job Name", "Reason"]
        if verbose:
            headers.append("Error Line")
        lines.append(tabulate(failed_data, headers=headers, tablefmt="simple"))

    # Missing jobs details (only in verbose mode)
    if report.missing and verbose:
        lines.append("\n--- Missing Jobs ---")
        missing_data = [[job.job_name] for job in report.missing]
        lines.append(tabulate(missing_data, headers=["Job Name"], tablefmt="simple"))

    return "\n".join(lines)


@click.command(name="check-jobs")
@click.option(
    "--job-dir",
    default="jobs",
    type=click.Path(exists=False),
    help="Directory containing job output files (.out files).",
)
@click.option(
    "--jobs-csv",
    default="jobs.csv",
    type=click.Path(exists=False),
    help="Path to jobs.csv file listing expected jobs.",
)
@click.option(
    "--fix",
    is_flag=True,
    help="Show suggested fixes for failed jobs.",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Show detailed output including error lines.",
)
def check_jobs(job_dir: str, jobs_csv: str, fix: bool, verbose: bool) -> None:
    """Check SLURM job outputs for errors and completion status.

    Scans job output files for common error patterns like time limit
    exceeded, out of memory, and other SLURM failures. Reports which
    jobs succeeded, failed, or are missing output files.

    \b
    Examples:
        rna-map-slurm check-jobs
        rna-map-slurm check-jobs --job-dir logs --verbose
        rna-map-slurm check-jobs --fix
    """
    job_dir_path = Path(job_dir)
    jobs_csv_path = Path(jobs_csv)

    # Validate paths
    if not job_dir_path.exists():
        click.echo(f"Error: Job directory not found: {job_dir_path}", err=True)
        click.echo("Run 'rna-map-slurm setup' first to generate jobs.", err=True)
        raise SystemExit(1)

    # Initialize checker
    checker = JobChecker(
        job_dir=job_dir_path,
        jobs_csv=jobs_csv_path if jobs_csv_path.exists() else None,
    )

    # Run checks
    report = checker.check_all()

    if report.total == 0:
        click.echo("No jobs found to check.")
        click.echo(f"Looked in: {job_dir_path}")
        return

    # Display report
    click.echo(_format_report(report, verbose=verbose))

    # Show fix suggestions if requested
    if fix and (report.failed or report.missing):
        suggestions = checker.suggest_fixes(report)
        if suggestions:
            click.echo("\n--- Suggested Fixes ---")
            for suggestion in suggestions:
                click.echo(f"  - {suggestion}")

    # Exit with error code if any failures
    if report.failure_count > 0 or report.missing_count > 0:
        raise SystemExit(1)
