"""Lightweight job watcher for monitoring SLURM job status."""

from __future__ import annotations

import re
import time
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable

from rna_map_slurm.jobs.slurm import get_user_jobs
from rna_map_slurm.models.job_types import PIPELINE_ORDER, JobType
from rna_map_slurm.utils.logging import get_logger

log = get_logger("jobs.watcher")


@dataclass
class WatcherConfig:
    """Configuration for the job watcher.

    Attributes:
        poll_interval_seconds: Seconds between status checks (default 60).
        timeout_hours: Maximum hours to watch before timing out (default 48).
        job_dir: Directory containing job outputs.
    """

    poll_interval_seconds: int = 60
    timeout_hours: float = 48.0
    job_dir: Path = field(default_factory=lambda: Path("jobs"))


@dataclass
class JobTypeStatus:
    """Status counts for a single job type."""

    job_type: str
    pending: int = 0
    running: int = 0
    completed: int = 0
    failed: int = 0
    dependency_failed: int = 0

    @property
    def total(self) -> int:
        return self.pending + self.running + self.completed + self.failed + self.dependency_failed

    @property
    def is_done(self) -> bool:
        """True if no jobs are pending or running."""
        return self.pending == 0 and self.running == 0


@dataclass
class WatcherStatus:
    """Current status of all watched jobs.

    Attributes:
        job_types: Status per job type.
        elapsed_seconds: Time since watching started.
        is_complete: True if all jobs are done.
        has_failures: True if any jobs failed.
    """

    job_types: dict[str, JobTypeStatus] = field(default_factory=dict)
    elapsed_seconds: float = 0.0
    timestamp: datetime = field(default_factory=datetime.now)

    @property
    def total_pending(self) -> int:
        return sum(jt.pending for jt in self.job_types.values())

    @property
    def total_running(self) -> int:
        return sum(jt.running for jt in self.job_types.values())

    @property
    def total_completed(self) -> int:
        return sum(jt.completed for jt in self.job_types.values())

    @property
    def total_failed(self) -> int:
        return sum(jt.failed for jt in self.job_types.values())

    @property
    def total_dependency_failed(self) -> int:
        return sum(jt.dependency_failed for jt in self.job_types.values())

    @property
    def total_jobs(self) -> int:
        return sum(jt.total for jt in self.job_types.values())

    @property
    def is_complete(self) -> bool:
        """True if no jobs are pending or running."""
        return self.total_pending == 0 and self.total_running == 0

    @property
    def has_failures(self) -> bool:
        """True if any jobs failed."""
        return self.total_failed > 0 or self.total_dependency_failed > 0

    def to_summary_line(self) -> str:
        """Single-line summary of status."""
        elapsed = _format_duration(self.elapsed_seconds)
        return (
            f"[{elapsed}] "
            f"Pending: {self.total_pending} | "
            f"Running: {self.total_running} | "
            f"Done: {self.total_completed} | "
            f"Failed: {self.total_failed + self.total_dependency_failed}"
        )

    def to_table(self) -> str:
        """Format status as a table."""
        from tabulate import tabulate

        headers = ["Job Type", "Pending", "Running", "Done", "Failed", "Dep.Fail"]
        rows = []

        # Sort by pipeline order
        sorted_types = sorted(
            self.job_types.values(),
            key=lambda x: _get_pipeline_order(x.job_type),
        )

        for jt in sorted_types:
            rows.append([
                jt.job_type,
                jt.pending,
                jt.running,
                jt.completed,
                jt.failed,
                jt.dependency_failed,
            ])

        # Add totals row
        rows.append([
            "TOTAL",
            self.total_pending,
            self.total_running,
            self.total_completed,
            self.total_failed,
            self.total_dependency_failed,
        ])

        return tabulate(rows, headers=headers, tablefmt="simple")


def _get_pipeline_order(job_type: str) -> int:
    """Get the order index for a job type."""
    for i, jt in enumerate(PIPELINE_ORDER):
        if jt.value == job_type:
            return i
    return 999


def _format_duration(seconds: float) -> str:
    """Format duration as human-readable string."""
    if seconds < 60:
        return f"{int(seconds)}s"
    elif seconds < 3600:
        mins = int(seconds // 60)
        secs = int(seconds % 60)
        return f"{mins}m {secs}s"
    else:
        hours = int(seconds // 3600)
        mins = int((seconds % 3600) // 60)
        return f"{hours}h {mins}m"


def _extract_job_type(job_name: str) -> str:
    """Extract job type from job name (e.g., 'rna-map-0001' -> 'rna-map')."""
    parts = job_name.rsplit("-", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return job_name


class JobWatcher:
    """Lightweight watcher for monitoring SLURM job status.

    Designed for speed over completeness - uses fast file existence checks
    and minimal SLURM queries.
    """

    def __init__(self, config: WatcherConfig | None = None) -> None:
        """Initialize the watcher.

        Args:
            config: Watcher configuration. Uses defaults if not provided.
        """
        self.config = config or WatcherConfig()
        self._start_time: float | None = None
        self._known_jobs: set[str] = set()

    def discover_jobs(self) -> set[str]:
        """Discover jobs from .sh files in job directory.

        Returns:
            Set of job names found.
        """
        jobs: set[str] = set()
        if not self.config.job_dir.exists():
            return jobs

        for sh_file in self.config.job_dir.glob("**/*.sh"):
            jobs.add(sh_file.stem)

        self._known_jobs = jobs
        return jobs

    def get_status(self) -> WatcherStatus:
        """Get current status with a single poll.

        This is the main entry point for getting job status.
        Fast and lightweight - designed for frequent polling.

        Returns:
            Current WatcherStatus.
        """
        if self._start_time is None:
            self._start_time = time.time()

        # Discover jobs if not done yet
        if not self._known_jobs:
            self.discover_jobs()

        # Query SLURM for running/pending jobs
        slurm_jobs = self._query_slurm()

        # Check outputs for completion
        completed_jobs = self._check_outputs_fast()

        # Build status per job type
        status = self._build_status(slurm_jobs, completed_jobs)
        status.elapsed_seconds = time.time() - self._start_time

        return status

    def _query_slurm(self) -> dict[str, str]:
        """Query SLURM for job states.

        Returns:
            Dict mapping job name to state (PD, R, etc.)
        """
        states: dict[str, str] = {}

        try:
            jobs = get_user_jobs()
            for job in jobs:
                name = job.get("Name", "")
                state = job.get("State", "")
                if name:
                    states[name] = state
        except Exception as e:
            log.warning(f"Failed to query SLURM: {e}")

        return states

    def _check_outputs_fast(self) -> dict[str, bool]:
        """Fast check for completed jobs by output file existence.

        Returns:
            Dict mapping job name to completion status (True = has output).
        """
        completed: dict[str, bool] = {}

        if not self.config.job_dir.exists():
            return completed

        for out_file in self.config.job_dir.glob("**/*.out"):
            job_name = out_file.stem
            # Just check existence - don't parse content for speed
            completed[job_name] = True

        return completed

    def _check_job_failed(self, job_name: str) -> tuple[bool, bool]:
        """Quick check if a job failed.

        Returns:
            Tuple of (is_failed, is_dependency_failed).
        """
        # Check for .err file with content
        job_type = _extract_job_type(job_name)
        err_file = self.config.job_dir / job_type / f"{job_name}.err"
        out_file = self.config.job_dir / job_type / f"{job_name}.out"

        if not out_file.exists():
            return False, False

        # Quick scan for common error patterns
        try:
            with open(out_file, "r", errors="replace") as f:
                content = f.read(8192)  # Only read first 8KB for speed

            # Check for dependency failure
            if "DependencyNeverSatisfied" in content or "Dependency" in content:
                return False, True

            # Check for other failures (quick patterns)
            fail_patterns = ["CANCELLED", "oom-kill", "OUT OF MEMORY", "Segmentation fault"]
            for pattern in fail_patterns:
                if pattern in content:
                    return True, False

        except OSError:
            pass

        return False, False

    def _build_status(
        self,
        slurm_jobs: dict[str, str],
        completed_jobs: dict[str, bool],
    ) -> WatcherStatus:
        """Build status from SLURM and output data.

        Args:
            slurm_jobs: Job name -> SLURM state.
            completed_jobs: Job name -> has output file.

        Returns:
            WatcherStatus with aggregated counts.
        """
        type_status: dict[str, JobTypeStatus] = {}

        for job_name in self._known_jobs:
            job_type = _extract_job_type(job_name)

            if job_type not in type_status:
                type_status[job_type] = JobTypeStatus(job_type=job_type)

            jts = type_status[job_type]

            # Check SLURM state first
            slurm_state = slurm_jobs.get(job_name, "")

            if slurm_state in ("PD", "PENDING"):
                jts.pending += 1
            elif slurm_state in ("R", "RUNNING"):
                jts.running += 1
            elif slurm_state in ("CG", "COMPLETING"):
                jts.running += 1  # Count as running
            elif job_name in completed_jobs:
                # Job finished - check if it failed
                failed, dep_failed = self._check_job_failed(job_name)
                if dep_failed:
                    jts.dependency_failed += 1
                elif failed:
                    jts.failed += 1
                else:
                    jts.completed += 1
            else:
                # No output yet, not in SLURM - might be queued or not submitted
                jts.pending += 1

        return WatcherStatus(job_types=type_status)

    def watch_until_complete(
        self,
        callback: Callable[[WatcherStatus], None] | None = None,
        stop_on_failure: bool = False,
    ) -> WatcherStatus:
        """Watch jobs until all complete or timeout.

        Args:
            callback: Optional function called after each poll with current status.
            stop_on_failure: If True, stop watching when any job fails.

        Returns:
            Final WatcherStatus.
        """
        self._start_time = time.time()
        timeout_seconds = self.config.timeout_hours * 3600

        while True:
            status = self.get_status()

            if callback:
                callback(status)

            # Check termination conditions
            if status.is_complete:
                log.info("All jobs completed")
                return status

            if stop_on_failure and status.has_failures:
                log.warning("Stopping due to job failures")
                return status

            if status.elapsed_seconds >= timeout_seconds:
                log.warning(f"Timeout after {self.config.timeout_hours} hours")
                return status

            # Wait before next poll
            time.sleep(self.config.poll_interval_seconds)

        return status


@dataclass
class JobTiming:
    """Timing information for a single job."""

    job_name: str
    job_type: str
    total_runtime_seconds: float = 0.0
    task_timings: list[tuple[str, float]] = field(default_factory=list)

    def to_human_readable(self) -> str:
        """Format timing as human-readable string."""
        total = _format_duration(self.total_runtime_seconds)
        return f"{self.job_name}: {total}"


class TimingParser:
    """Parse timing information from job output files."""

    # Pattern for time_it decorator output
    TIMING_PATTERN = re.compile(
        r"Function '(\w+)' executed in ([\d.]+) seconds"
    )

    def __init__(self, job_dir: Path) -> None:
        """Initialize the parser.

        Args:
            job_dir: Directory containing job outputs.
        """
        self.job_dir = job_dir

    def parse_job(self, job_name: str) -> JobTiming | None:
        """Parse timing from a job's output file.

        Args:
            job_name: Name of the job.

        Returns:
            JobTiming or None if output not found.
        """
        job_type = _extract_job_type(job_name)
        out_file = self.job_dir / job_type / f"{job_name}.out"

        if not out_file.exists():
            return None

        timing = JobTiming(job_name=job_name, job_type=job_type)
        task_timings: list[tuple[str, float]] = []

        try:
            with open(out_file, "r", errors="replace") as f:
                for line in f:
                    match = self.TIMING_PATTERN.search(line)
                    if match:
                        func_name = match.group(1)
                        seconds = float(match.group(2))
                        task_timings.append((func_name, seconds))
                        timing.total_runtime_seconds += seconds

            timing.task_timings = task_timings
        except OSError:
            return None

        return timing

    def parse_all(self) -> dict[str, list[JobTiming]]:
        """Parse timing for all jobs.

        Returns:
            Dict mapping job type to list of JobTiming.
        """
        timings: dict[str, list[JobTiming]] = defaultdict(list)

        for out_file in self.job_dir.glob("**/*.out"):
            job_name = out_file.stem
            timing = self.parse_job(job_name)
            if timing:
                timings[timing.job_type].append(timing)

        return dict(timings)

    def get_summary(self) -> str:
        """Get human-readable timing summary.

        Returns:
            Formatted timing summary.
        """
        from tabulate import tabulate

        all_timings = self.parse_all()

        if not all_timings:
            return "No timing data found"

        headers = ["Job Type", "Jobs", "Total Time", "Avg Time", "Min", "Max"]
        rows = []

        for job_type in sorted(all_timings.keys(), key=_get_pipeline_order):
            timings = all_timings[job_type]
            if not timings:
                continue

            runtimes = [t.total_runtime_seconds for t in timings]
            total = sum(runtimes)
            avg = total / len(runtimes)
            min_time = min(runtimes)
            max_time = max(runtimes)

            rows.append([
                job_type,
                len(timings),
                _format_duration(total),
                _format_duration(avg),
                _format_duration(min_time),
                _format_duration(max_time),
            ])

        return tabulate(rows, headers=headers, tablefmt="simple")
