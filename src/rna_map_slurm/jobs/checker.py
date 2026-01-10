"""Job output checking and error detection for SLURM jobs."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import pandas as pd


@dataclass
class JobStatus:
    """Status of a single SLURM job."""

    job_name: str
    status: Literal["success", "failed", "missing"]
    reason: str | None = None
    output_file: Path | None = None
    error_line: str | None = None


@dataclass
class CheckReport:
    """Summary report of job status checks."""

    total: int
    succeeded: list[JobStatus] = field(default_factory=list)
    failed: list[JobStatus] = field(default_factory=list)
    missing: list[JobStatus] = field(default_factory=list)

    @property
    def success_count(self) -> int:
        """Number of successfully completed jobs."""
        return len(self.succeeded)

    @property
    def failure_count(self) -> int:
        """Number of failed jobs."""
        return len(self.failed)

    @property
    def missing_count(self) -> int:
        """Number of jobs with no output file."""
        return len(self.missing)

    def summary(self) -> str:
        """Return a text summary of the report."""
        lines = [
            f"Total jobs: {self.total}",
            f"Succeeded:  {self.success_count}",
            f"Failed:     {self.failure_count}",
            f"Missing:    {self.missing_count}",
        ]
        return "\n".join(lines)


class JobChecker:
    """Check SLURM job output files for errors and completion status."""

    # Error patterns to search for in job output files
    # Maps regex pattern -> (error_type, human-readable description)
    SLURM_ERROR_PATTERNS: dict[str, tuple[str, str]] = {
        r"CANCELLED.*DUE TO TIME LIMIT": (
            "time_limit",
            "Job exceeded wall time limit",
        ),
        r"slurmstepd:.*CANCELLED.*TIME LIMIT": (
            "time_limit",
            "Job exceeded wall time limit",
        ),
        r"oom-kill|OUT OF MEMORY|Killed.*memory|Cannot allocate memory": (
            "out_of_memory",
            "Job ran out of memory",
        ),
        r"slurmstepd: error: Detected \d+ oom-kill": (
            "out_of_memory",
            "Job killed by OOM killer",
        ),
        r"DUE TO NODE FAILURE": (
            "node_failure",
            "Job failed due to node failure",
        ),
        r"Segmentation fault|SIGSEGV": (
            "segfault",
            "Job crashed with segmentation fault",
        ),
        r"srun: error: .* task \d+: Exited with exit code": (
            "exit_error",
            "Job exited with non-zero exit code",
        ),
    }

    # Error patterns for stderr files (Python exceptions, etc.)
    STDERR_ERROR_PATTERNS: dict[str, tuple[str, str]] = {
        r"^Traceback \(most recent call last\):": (
            "python_exception",
            "Python exception occurred",
        ),
        r"^Error:|^ERROR:": (
            "error_message",
            "Error message in stderr",
        ),
    }

    # Suggestions for fixing different error types
    FIX_SUGGESTIONS: dict[str, str] = {
        "time_limit": (
            "Increase time limit in slurm_options.{job_type}.time or reduce tasks_per_job"
        ),
        "out_of_memory": (
            "Increase memory in slurm_options.{job_type}.mem-per-cpu or reduce tasks_per_job"
        ),
        "node_failure": "This is a cluster issue - resubmit the failed jobs",
        "segfault": "Check input data for corruption or report bug to developers",
        "exit_error": "Check the job output for specific error messages",
        "python_exception": "Check the .err file for the full traceback",
        "error_message": "Check the .err file for error details",
    }

    def __init__(
        self,
        job_dir: Path,
        jobs_csv: Path | None = None,
        max_lines_to_check: int = 1000,
    ) -> None:
        """Initialize the job checker.

        Args:
            job_dir: Directory containing job output files (.out files)
            jobs_csv: Path to jobs.csv file listing expected jobs (optional)
            max_lines_to_check: Maximum lines to scan from end of each file
        """
        self.job_dir = Path(job_dir)
        self.jobs_csv = Path(jobs_csv) if jobs_csv else None
        self.max_lines_to_check = max_lines_to_check

    def _load_expected_jobs(self) -> set[str]:
        """Load expected job names from jobs.csv."""
        if self.jobs_csv is None or not self.jobs_csv.exists():
            return set()

        df = pd.read_csv(self.jobs_csv)
        if "job_path" not in df.columns:
            return set()

        # Extract job names from paths (e.g., "jobs/split-fastq-001.sh" -> "split-fastq-001")
        job_names = set()
        for job_path in df["job_path"]:
            job_name = Path(job_path).stem
            job_names.add(job_name)
        return job_names

    def _find_output_files(self) -> dict[str, Path]:
        """Find all .out files in the job directory and its subdirectories.

        Returns:
            Dict mapping job name to output file path
        """
        output_files: dict[str, Path] = {}

        if not self.job_dir.exists():
            return output_files

        # Search recursively for .out files in subdirectories
        for out_file in self.job_dir.glob("**/*.out"):
            job_name = out_file.stem
            output_files[job_name] = out_file

        return output_files

    def _parse_output_file(self, path: Path) -> tuple[bool, str | None, str | None]:
        """Scan an output file for error patterns.

        Args:
            path: Path to the job output file

        Returns:
            Tuple of (success, error_type, matching_line)
        """
        if not path.exists():
            return False, "file_not_found", None

        try:
            with open(path, "r", errors="replace") as f:
                lines = f.readlines()
        except OSError as e:
            return False, f"read_error: {e}", None

        # Check the last N lines for errors (errors usually appear at the end)
        lines_to_check = lines[-self.max_lines_to_check :]

        for line in lines_to_check:
            for pattern, (error_type, _) in self.SLURM_ERROR_PATTERNS.items():
                if re.search(pattern, line, re.IGNORECASE):
                    return False, error_type, line.strip()

        return True, None, None

    def _parse_stderr_file(self, path: Path) -> tuple[bool, str | None, str | None]:
        """Scan a stderr file for error patterns.

        Args:
            path: Path to the job stderr file (.err)

        Returns:
            Tuple of (success, error_type, matching_line)
        """
        if not path.exists():
            return True, None, None  # No stderr file is OK

        try:
            with open(path, "r", errors="replace") as f:
                lines = f.readlines()
        except OSError:
            return True, None, None  # Can't read is OK for stderr

        # Check all lines for error patterns
        for line in lines:
            for pattern, (error_type, _) in self.STDERR_ERROR_PATTERNS.items():
                if re.search(pattern, line, re.MULTILINE):
                    return False, error_type, line.strip()

        return True, None, None

    def _check_single_job(
        self, job_name: str, output_file: Path | None
    ) -> JobStatus:
        """Check the status of a single job.

        Args:
            job_name: Name of the job
            output_file: Path to output file, or None if not found

        Returns:
            JobStatus with the check results
        """
        if output_file is None:
            return JobStatus(
                job_name=job_name,
                status="missing",
                reason="No output file found",
                output_file=None,
            )

        # Check stdout file for SLURM errors
        success, error_type, error_line = self._parse_output_file(output_file)

        if not success:
            # Get human-readable description for the error
            reason = error_type
            for pattern, (etype, description) in self.SLURM_ERROR_PATTERNS.items():
                if etype == error_type:
                    reason = description
                    break

            return JobStatus(
                job_name=job_name,
                status="failed",
                reason=reason,
                output_file=output_file,
                error_line=error_line,
            )

        # Check stderr file for Python exceptions and other errors
        stderr_file = output_file.with_suffix(".err")
        stderr_success, stderr_error_type, stderr_error_line = self._parse_stderr_file(
            stderr_file
        )

        if not stderr_success:
            # Get human-readable description for stderr error
            reason = stderr_error_type
            for pattern, (etype, description) in self.STDERR_ERROR_PATTERNS.items():
                if etype == stderr_error_type:
                    reason = description
                    break

            return JobStatus(
                job_name=job_name,
                status="failed",
                reason=reason,
                output_file=output_file,
                error_line=stderr_error_line,
            )

        return JobStatus(
            job_name=job_name,
            status="success",
            output_file=output_file,
        )

    def check_all(self) -> CheckReport:
        """Check all jobs and return a report.

        Returns:
            CheckReport with results for all jobs
        """
        output_files = self._find_output_files()
        expected_jobs = self._load_expected_jobs()

        # If we have a jobs.csv, use it to determine expected jobs
        # Otherwise, just check whatever output files exist
        if expected_jobs:
            all_jobs = expected_jobs
        else:
            all_jobs = set(output_files.keys())

        succeeded: list[JobStatus] = []
        failed: list[JobStatus] = []
        missing: list[JobStatus] = []

        for job_name in sorted(all_jobs):
            output_file = output_files.get(job_name)
            status = self._check_single_job(job_name, output_file)

            if status.status == "success":
                succeeded.append(status)
            elif status.status == "failed":
                failed.append(status)
            else:
                missing.append(status)

        return CheckReport(
            total=len(all_jobs),
            succeeded=succeeded,
            failed=failed,
            missing=missing,
        )

    def suggest_fixes(self, report: CheckReport) -> list[str]:
        """Generate fix suggestions based on the report.

        Args:
            report: CheckReport from check_all()

        Returns:
            List of suggested fixes
        """
        suggestions: list[str] = []
        error_types_seen: dict[str, list[str]] = {}

        # Combine all error patterns for lookup
        all_patterns = {**self.SLURM_ERROR_PATTERNS, **self.STDERR_ERROR_PATTERNS}

        # Group failures by error type
        for job in report.failed:
            if job.reason:
                # Find the error type for this reason
                error_type = None
                for pattern, (etype, description) in all_patterns.items():
                    if description == job.reason:
                        error_type = etype
                        break

                if error_type:
                    if error_type not in error_types_seen:
                        error_types_seen[error_type] = []
                    error_types_seen[error_type].append(job.job_name)

        # Generate suggestions for each error type
        for error_type, job_names in error_types_seen.items():
            if error_type in self.FIX_SUGGESTIONS:
                # Extract job type from first job name (e.g., "split-fastq-001" -> "split-fastq")
                job_type = "-".join(job_names[0].rsplit("-", 1)[0].split("-"))
                suggestion = self.FIX_SUGGESTIONS[error_type].format(job_type=job_type)
                count = len(job_names)
                suggestions.append(f"{count} job(s) with {error_type}: {suggestion}")

        if report.missing:
            suggestions.append(
                f"{len(report.missing)} job(s) have no output file - "
                "they may not have been submitted or are still running"
            )

        return suggestions
