"""CLI command for generating pipeline summary reports."""

from __future__ import annotations

from pathlib import Path

import click

from rna_map_slurm.jobs.summary import SummaryBuilder
from rna_map_slurm.jobs.timing import get_timing_report
from rna_map_slurm.utils.logging import get_logger

log = get_logger("cli.summary")


@click.command(name="summary")
@click.option(
    "--job-dir",
    default="jobs",
    type=click.Path(exists=False),
    help="Directory containing job output files.",
)
@click.option(
    "--jobs-csv",
    default="jobs.csv",
    type=click.Path(exists=False),
    help="Path to jobs.csv file.",
)
@click.option(
    "--validate-outputs",
    is_flag=True,
    help="Validate that expected output files were created.",
)
@click.option(
    "--no-timing",
    is_flag=True,
    help="Skip SLURM timing analysis (faster if sacct is slow).",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["table", "json"]),
    default="table",
    help="Output format.",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    help="Write output to file instead of stdout.",
)
def summary(
    job_dir: str,
    jobs_csv: str,
    validate_outputs: bool,
    no_timing: bool,
    output_format: str,
    output: str | None,
) -> None:
    """Generate a summary report of pipeline job status and timing.

    Shows success/failure counts, timing statistics, and parallelization
    metrics for each job type.

    \b
    Examples:
        rna-map-slurm summary
        rna-map-slurm summary --validate-outputs
        rna-map-slurm summary --no-timing
        rna-map-slurm summary --format json -o report.json
    """
    job_dir_path = Path(job_dir)
    jobs_csv_path = Path(jobs_csv)

    if not job_dir_path.exists():
        click.echo(f"Error: Job directory not found: {job_dir}", err=True)
        raise SystemExit(1)

    # Build job status summary
    builder = SummaryBuilder(
        job_dir=job_dir_path,
        jobs_csv=jobs_csv_path if jobs_csv_path.exists() else None,
    )

    try:
        pipeline_summary = builder.build(validate_outputs=validate_outputs)
    except Exception as e:
        click.echo(f"Error generating summary: {e}", err=True)
        raise SystemExit(1) from None

    # Get timing report if available
    timing_report = None
    if not no_timing:
        submitted_jobs_file = Path("submitted_jobs.txt")
        if submitted_jobs_file.exists():
            timing_report = get_timing_report(submitted_jobs_file)

    # Format output
    if output_format == "json":
        result = pipeline_summary.to_json()
    else:
        result = _format_table_output(
            pipeline_summary, timing_report, validate_outputs
        )

    # Write output
    if output:
        Path(output).write_text(result)
        click.echo(f"Summary written to: {output}")
    else:
        click.echo(result)

    # Exit with error if there were failures
    if pipeline_summary.total_failed > 0 or pipeline_summary.total_missing > 0:
        raise SystemExit(1)


def _format_table_output(
    pipeline_summary,
    timing_report,
    validate_outputs: bool,
) -> str:
    """Format the complete summary output as a table.

    Args:
        pipeline_summary: Job status summary.
        timing_report: Optional timing report from sacct.
        validate_outputs: Whether output validation was enabled.

    Returns:
        Formatted string output.
    """
    sections = []

    # Job status section
    status_section = "Pipeline Summary\n" + "=" * 70 + "\n\n"
    status_section += pipeline_summary.to_table()
    if validate_outputs:
        status_section += "\n\n(Output validation enabled)"
    sections.append(status_section)

    # Timing section
    if timing_report:
        sections.append("\n\n" + timing_report.to_table())
    else:
        sections.append(
            "\n\n(Timing data not available - "
            "run on SLURM cluster with submitted_jobs.txt)"
        )

    return "".join(sections)
