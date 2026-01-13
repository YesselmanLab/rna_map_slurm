"""SLURM job timing analysis using sacct data."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from rna_map_slurm.models.job_types import PIPELINE_ORDER
from rna_map_slurm.utils.logging import get_logger

log = get_logger("jobs.timing")


@dataclass
class JobTimingRecord:
    """Timing data for a single job."""

    job_id: str
    job_name: str
    job_type: str
    state: str
    elapsed_seconds: int
    submit_time: datetime | None = None
    start_time: datetime | None = None
    end_time: datetime | None = None

    @property
    def wait_seconds(self) -> int:
        """Time spent waiting in queue."""
        if self.submit_time and self.start_time:
            return int((self.start_time - self.submit_time).total_seconds())
        return 0


@dataclass
class JobTypeTimingStats:
    """Aggregated timing statistics for a job type."""

    job_type: str
    count: int = 0
    total_seconds: int = 0
    min_seconds: int = 0
    max_seconds: int = 0
    avg_seconds: float = 0.0
    total_wait_seconds: int = 0
    completed: int = 0
    failed: int = 0

    def add_job(self, record: JobTimingRecord) -> None:
        """Add a job's timing to the stats."""
        self.count += 1
        self.total_seconds += record.elapsed_seconds
        self.total_wait_seconds += record.wait_seconds

        if self.count == 1:
            self.min_seconds = record.elapsed_seconds
            self.max_seconds = record.elapsed_seconds
        else:
            self.min_seconds = min(self.min_seconds, record.elapsed_seconds)
            self.max_seconds = max(self.max_seconds, record.elapsed_seconds)

        self.avg_seconds = self.total_seconds / self.count

        if record.state == "COMPLETED":
            self.completed += 1
        elif record.state == "FAILED":
            self.failed += 1


@dataclass
class TimingReport:
    """Complete timing report for a pipeline run."""

    job_types: dict[str, JobTypeTimingStats] = field(default_factory=dict)
    first_start: datetime | None = None
    last_end: datetime | None = None
    total_jobs: int = 0
    peak_concurrent: int = 0

    @property
    def wall_clock_seconds(self) -> int:
        """Total wall clock time from first job start to last job end."""
        if self.first_start and self.last_end:
            return int((self.last_end - self.first_start).total_seconds())
        return 0

    @property
    def total_compute_seconds(self) -> int:
        """Sum of all job runtimes."""
        return sum(jt.total_seconds for jt in self.job_types.values())

    @property
    def speedup_factor(self) -> float:
        """Parallelization speedup (compute time / wall clock time)."""
        if self.wall_clock_seconds > 0:
            return self.total_compute_seconds / self.wall_clock_seconds
        return 0.0

    @property
    def avg_concurrency(self) -> float:
        """Average number of concurrent jobs."""
        return self.speedup_factor

    def to_table(self) -> str:
        """Format report as markdown."""
        from tabulate import tabulate

        lines = []
        lines.append("## Timing Report")
        lines.append("")

        # Summary metrics as markdown table
        summary_rows = [
            ["Wall clock time", _format_duration(self.wall_clock_seconds)],
            ["Total compute time", _format_duration(self.total_compute_seconds)],
            ["Speedup factor", f"{self.speedup_factor:.1f}x"],
            ["Peak concurrent", f"{self.peak_concurrent} jobs"],
            ["Avg concurrency", f"{self.avg_concurrency:.1f} jobs"],
            ["Total jobs", str(self.total_jobs)],
        ]
        lines.append(tabulate(summary_rows, headers=["Metric", "Value"], tablefmt="github"))
        lines.append("")

        # Per-job-type table
        lines.append("### Job Timing by Type")
        lines.append("")
        headers = ["Job Type", "Count", "Total", "Avg", "Min", "Max", "Done", "Fail"]
        rows = []

        sorted_types = sorted(
            self.job_types.values(),
            key=lambda x: -x.total_seconds,
        )

        for jt in sorted_types:
            rows.append([
                jt.job_type,
                jt.count,
                _format_duration(jt.total_seconds),
                _format_duration(int(jt.avg_seconds)),
                _format_duration(jt.min_seconds),
                _format_duration(jt.max_seconds),
                jt.completed,
                jt.failed,
            ])

        lines.append(tabulate(rows, headers=headers, tablefmt="github"))
        lines.append("")

        # Queue wait times
        lines.append("### Queue Wait Times")
        lines.append("")
        wait_headers = ["Job Type", "Avg Wait", "Jobs"]
        wait_rows = []

        for jt in sorted_types:
            if jt.count > 0:
                avg_wait = jt.total_wait_seconds / jt.count
                wait_rows.append([
                    jt.job_type,
                    _format_duration(int(avg_wait)),
                    jt.count,
                ])

        lines.append(tabulate(wait_rows, headers=wait_headers, tablefmt="github"))

        return "\n".join(lines)


def _format_duration(seconds: int) -> str:
    """Format seconds as human-readable duration."""
    if seconds < 60:
        return f"{seconds}s"
    elif seconds < 3600:
        mins = seconds // 60
        secs = seconds % 60
        return f"{mins}m {secs}s"
    else:
        hours = seconds // 3600
        mins = (seconds % 3600) // 60
        return f"{hours}h {mins}m"


def _parse_elapsed(elapsed_str: str) -> int:
    """Parse SLURM elapsed time string to seconds."""
    # Format: HH:MM:SS or D-HH:MM:SS
    try:
        if "-" in elapsed_str:
            days, time_part = elapsed_str.split("-")
            days = int(days)
        else:
            days = 0
            time_part = elapsed_str

        parts = time_part.split(":")
        if len(parts) == 3:
            hours, mins, secs = int(parts[0]), int(parts[1]), int(parts[2])
            return days * 86400 + hours * 3600 + mins * 60 + secs
    except (ValueError, IndexError):
        pass
    return 0


def _parse_datetime(dt_str: str) -> datetime | None:
    """Parse SLURM datetime string."""
    if not dt_str or dt_str == "Unknown":
        return None
    try:
        return datetime.strptime(dt_str, "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        return None


def _extract_job_type(job_name: str) -> str:
    """Extract job type from job name (e.g., 'rna-map-0001.sh' -> 'rna-map')."""
    name = job_name.replace(".sh", "")
    parts = name.rsplit("-", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return name


def _get_pipeline_order(job_type: str) -> int:
    """Get sort order for job type based on pipeline order."""
    for i, jt in enumerate(PIPELINE_ORDER):
        if jt.value == job_type:
            return i
    return 999


def get_timing_report(submitted_jobs_file: Path) -> TimingReport | None:
    """Generate timing report from submitted jobs file.

    Args:
        submitted_jobs_file: Path to submitted_jobs.txt with job IDs.

    Returns:
        TimingReport or None if unable to query SLURM.
    """
    if not submitted_jobs_file.exists():
        log.warning(f"Submitted jobs file not found: {submitted_jobs_file}")
        return None

    # Read job IDs
    job_ids = submitted_jobs_file.read_text().strip().split("\n")
    if not job_ids:
        log.warning("No job IDs found in submitted_jobs.txt")
        return None

    # Query sacct for timing data
    records = _query_sacct(job_ids)
    if not records:
        log.warning("No timing data returned from sacct")
        return None

    # Build report
    return _build_report(records)


def _query_sacct(job_ids: list[str]) -> list[JobTimingRecord]:
    """Query sacct for job timing data.

    Supports both regular job IDs and array job IDs (format: JOBID_TASKID).

    Args:
        job_ids: List of SLURM job IDs.

    Returns:
        List of JobTimingRecord objects.
    """
    records = []
    job_set = set(job_ids)

    # Separate array jobs from regular jobs
    array_base_ids = set()
    regular_ids = set()
    for jid in job_ids:
        if "_" in jid:
            # Array job: extract base ID (e.g., "12345_0" -> "12345")
            array_base_ids.add(jid.split("_")[0])
        else:
            regular_ids.add(jid)

    # Build comma-separated job list for sacct
    # For arrays, we query the base ID which returns all tasks
    all_query_ids = list(regular_ids | array_base_ids)

    # Query in batches to avoid command line length limits
    batch_size = 100
    for i in range(0, len(all_query_ids), batch_size):
        batch = all_query_ids[i : i + batch_size]
        batch_records = _query_sacct_batch(batch, job_set)
        records.extend(batch_records)

    return records


def _parse_sacct_line(line: str, valid_ids: set[str]) -> JobTimingRecord | None:
    """Parse a single sacct output line into a JobTimingRecord.

    Args:
        line: Pipe-separated sacct output line.
        valid_ids: Set of valid job IDs to filter results.

    Returns:
        JobTimingRecord or None if line should be skipped.
    """
    parts = line.split("|")
    if len(parts) < 7:
        return None

    job_id, job_name, state, elapsed, submit, start, end = parts[:7]

    # Skip batch/extern steps (e.g., "12345.batch", "12345.extern")
    if "." in job_id:
        return None

    # Only include jobs from our submission
    if job_id not in valid_ids:
        return None

    return JobTimingRecord(
        job_id=job_id,
        job_name=job_name,
        job_type=_extract_job_type(job_name),
        state=state,
        elapsed_seconds=_parse_elapsed(elapsed),
        submit_time=_parse_datetime(submit),
        start_time=_parse_datetime(start),
        end_time=_parse_datetime(end),
    )


def _query_sacct_batch(
    query_ids: list[str], valid_ids: set[str]
) -> list[JobTimingRecord]:
    """Query sacct for a batch of job IDs.

    Args:
        query_ids: Job IDs to query.
        valid_ids: Set of valid job IDs to filter results.

    Returns:
        List of JobTimingRecord objects.
    """
    cmd = [
        "sacct",
        "-j",
        ",".join(query_ids),
        "--format=JobID,JobName%50,State,Elapsed,Submit,Start,End",
        "-P",
        "--noheader",
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120
        )
        if result.returncode != 0:
            log.warning(f"sacct failed: {result.stderr}")
            return []

        records = []
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            record = _parse_sacct_line(line, valid_ids)
            if record:
                records.append(record)
        return records

    except subprocess.TimeoutExpired:
        log.warning("sacct query timed out")
    except FileNotFoundError:
        log.warning("sacct command not found - not on a SLURM cluster?")
    except Exception as e:
        log.warning(f"Error querying sacct: {e}")

    return []


def _build_report(records: list[JobTimingRecord]) -> TimingReport:
    """Build timing report from job records.

    Args:
        records: List of JobTimingRecord objects.

    Returns:
        Populated TimingReport.
    """
    report = TimingReport()
    report.total_jobs = len(records)

    # Aggregate by job type
    for record in records:
        if record.job_type not in report.job_types:
            report.job_types[record.job_type] = JobTypeTimingStats(
                job_type=record.job_type
            )
        report.job_types[record.job_type].add_job(record)

        # Track overall time bounds
        if record.start_time and (
            report.first_start is None or record.start_time < report.first_start
        ):
            report.first_start = record.start_time
        if record.end_time and (
            report.last_end is None or record.end_time > report.last_end
        ):
            report.last_end = record.end_time

    # Calculate peak concurrency
    report.peak_concurrent = _calculate_peak_concurrent(records)

    return report


def _calculate_peak_concurrent(records: list[JobTimingRecord]) -> int:
    """Calculate peak number of concurrent jobs.

    Args:
        records: List of JobTimingRecord objects.

    Returns:
        Maximum number of jobs running simultaneously.
    """
    events: list[tuple[datetime, int]] = []

    for record in records:
        if record.start_time and record.end_time:
            events.append((record.start_time, 1))
            events.append((record.end_time, -1))

    if not events:
        return 0

    events.sort(key=lambda x: x[0])

    current = 0
    peak = 0
    for _, delta in events:
        current += delta
        peak = max(peak, current)

    return peak
