"""Watch command for monitoring SLURM job status."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from rna_map_slurm.jobs.watcher import JobWatcher, TimingParser, WatcherConfig, WatcherStatus
from rna_map_slurm.utils.logging import get_logger, setup_logging

log = get_logger("cli.watch")


def _print_status_update(status: WatcherStatus) -> None:
    """Print status update to terminal."""
    # Clear screen and print status table
    click.clear()
    click.echo("RNA Map SLURM - Job Watcher")
    click.echo("=" * 50)
    click.echo()
    click.echo(status.to_table())
    click.echo()
    click.echo(status.to_summary_line())
    click.echo()
    click.echo("Press Ctrl+C to stop watching")


def _print_status_line(status: WatcherStatus) -> None:
    """Print single-line status update."""
    click.echo(status.to_summary_line())


@click.command("watch")
@click.option(
    "--poll-interval",
    "-p",
    default=60,
    type=int,
    help="Seconds between status checks (default: 60).",
)
@click.option(
    "--timeout",
    "-t",
    default=48.0,
    type=float,
    help="Hours before timeout (default: 48).",
)
@click.option(
    "--job-dir",
    default="jobs",
    type=click.Path(exists=False),
    help="Directory containing job outputs.",
)
@click.option(
    "--stop-on-failure",
    is_flag=True,
    help="Stop watching when any job fails.",
)
@click.option(
    "--once",
    is_flag=True,
    help="Check status once and exit (no watching).",
)
@click.option(
    "--compact",
    is_flag=True,
    help="Use compact single-line output instead of table.",
)
@click.option(
    "--timing",
    is_flag=True,
    help="Show timing summary for completed jobs.",
)
def watch(
    poll_interval: int,
    timeout: float,
    job_dir: str,
    stop_on_failure: bool,
    once: bool,
    compact: bool,
    timing: bool,
) -> None:
    """Watch SLURM jobs and display status.

    Monitors job status by querying SLURM and checking output files.
    Displays a live-updating status table showing job progress.

    \b
    Examples:
        rna-map-slurm watch                    # Watch with defaults
        rna-map-slurm watch --poll-interval 30 # Check every 30 seconds
        rna-map-slurm watch --once             # Check once and exit
        rna-map-slurm watch --timing           # Show timing summary
        rna-map-slurm watch --compact          # Single-line output
    """
    setup_logging()
    job_dir_path = Path(job_dir)

    if not job_dir_path.exists():
        click.echo(f"Error: Job directory '{job_dir}' not found.", err=True)
        click.echo("Run 'rna-map-slurm setup' first to create jobs.", err=True)
        sys.exit(1)

    # Show timing summary if requested
    if timing:
        parser = TimingParser(job_dir_path)
        click.echo("\nTiming Summary")
        click.echo("=" * 50)
        click.echo(parser.get_summary())
        click.echo()
        if once:
            return

    config = WatcherConfig(
        poll_interval_seconds=poll_interval,
        timeout_hours=timeout,
        job_dir=job_dir_path,
    )
    watcher = JobWatcher(config)

    # Discover jobs
    jobs = watcher.discover_jobs()
    if not jobs:
        click.echo("No jobs found in job directory.", err=True)
        sys.exit(1)

    log.info(f"Found {len(jobs)} jobs to watch")

    if once:
        # Single check mode
        status = watcher.get_status()
        if compact:
            _print_status_line(status)
        else:
            click.echo()
            click.echo(status.to_table())
            click.echo()
            click.echo(status.to_summary_line())

        if status.has_failures:
            sys.exit(1)
        return

    # Watch mode
    try:
        callback = _print_status_line if compact else _print_status_update
        final_status = watcher.watch_until_complete(
            callback=callback,
            stop_on_failure=stop_on_failure,
        )

        # Print final summary
        click.echo()
        click.echo("=" * 50)
        click.echo("Final Status:")
        click.echo(final_status.to_table())
        click.echo()

        if timing:
            parser = TimingParser(job_dir_path)
            click.echo("\nTiming Summary")
            click.echo("=" * 50)
            click.echo(parser.get_summary())

        if final_status.has_failures:
            click.echo(
                f"\nWarning: {final_status.total_failed + final_status.total_dependency_failed} "
                "jobs failed",
                err=True,
            )
            sys.exit(1)

        click.echo("\nAll jobs completed successfully!")

    except KeyboardInterrupt:
        click.echo("\n\nStopped watching.")
        sys.exit(0)
