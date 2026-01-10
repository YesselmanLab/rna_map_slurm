"""Pipeline summary and reporting for job status."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
from tabulate import tabulate

from rna_map_slurm.jobs.checker import CheckReport, JobChecker
from rna_map_slurm.jobs.validation import OutputValidator, ValidationResult
from rna_map_slurm.utils.logging import get_logger

log = get_logger("jobs.summary")


@dataclass
class JobTypeSummary:
    """Summary statistics for a single job type.

    Attributes:
        job_type: Name of the job type.
        total: Total number of jobs.
        succeeded: Number of successful jobs.
        failed: Number of failed jobs.
        missing: Number of jobs with missing output.
        success_rate: Percentage of successful jobs.
        validated: Number of jobs with validated outputs.
        validation_failed: Number of jobs that failed validation.
    """

    job_type: str
    total: int = 0
    succeeded: int = 0
    failed: int = 0
    missing: int = 0
    validated: int = 0
    validation_failed: int = 0

    @property
    def success_rate(self) -> float:
        """Calculate success rate as percentage."""
        if self.total == 0:
            return 0.0
        return (self.succeeded / self.total) * 100


@dataclass
class PipelineSummary:
    """Complete pipeline summary across all job types.

    Attributes:
        job_types: List of per-type summaries.
        total_jobs: Total number of jobs across all types.
        total_succeeded: Total successful jobs.
        total_failed: Total failed jobs.
        total_missing: Total missing jobs.
        overall_success_rate: Overall success percentage.
        validation_enabled: Whether output validation was run.
    """

    job_types: list[JobTypeSummary] = field(default_factory=list)
    validation_enabled: bool = False

    @property
    def total_jobs(self) -> int:
        return sum(jt.total for jt in self.job_types)

    @property
    def total_succeeded(self) -> int:
        return sum(jt.succeeded for jt in self.job_types)

    @property
    def total_failed(self) -> int:
        return sum(jt.failed for jt in self.job_types)

    @property
    def total_missing(self) -> int:
        return sum(jt.missing for jt in self.job_types)

    @property
    def overall_success_rate(self) -> float:
        if self.total_jobs == 0:
            return 0.0
        return (self.total_succeeded / self.total_jobs) * 100

    def to_table(self) -> str:
        """Generate a formatted table of the summary."""
        headers = ["Job Type", "Total", "Succeeded", "Failed", "Missing", "Success %"]
        rows = []

        for jt in self.job_types:
            rows.append([
                jt.job_type,
                jt.total,
                jt.succeeded,
                jt.failed,
                jt.missing,
                f"{jt.success_rate:.1f}%",
            ])

        # Add totals row
        rows.append([
            "TOTAL",
            self.total_jobs,
            self.total_succeeded,
            self.total_failed,
            self.total_missing,
            f"{self.overall_success_rate:.1f}%",
        ])

        return tabulate(rows, headers=headers, tablefmt="simple")

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "job_types": [
                {
                    "job_type": jt.job_type,
                    "total": jt.total,
                    "succeeded": jt.succeeded,
                    "failed": jt.failed,
                    "missing": jt.missing,
                    "success_rate": jt.success_rate,
                    "validated": jt.validated,
                    "validation_failed": jt.validation_failed,
                }
                for jt in self.job_types
            ],
            "totals": {
                "total_jobs": self.total_jobs,
                "total_succeeded": self.total_succeeded,
                "total_failed": self.total_failed,
                "total_missing": self.total_missing,
                "overall_success_rate": self.overall_success_rate,
            },
            "validation_enabled": self.validation_enabled,
        }

    def to_json(self) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(), indent=2)


class SummaryBuilder:
    """Build pipeline summary from job check results."""

    def __init__(
        self,
        job_dir: Path = Path("jobs"),
        jobs_csv: Path | None = None,
    ) -> None:
        """Initialize the summary builder.

        Args:
            job_dir: Directory containing job outputs.
            jobs_csv: Path to jobs.csv for job metadata.
        """
        self.job_dir = job_dir
        self.jobs_csv = jobs_csv
        self.checker = JobChecker(job_dir=job_dir, jobs_csv=jobs_csv)
        self.validator = OutputValidator()

    def build(self, validate_outputs: bool = False) -> PipelineSummary:
        """Build the pipeline summary.

        Args:
            validate_outputs: Whether to validate job outputs.

        Returns:
            PipelineSummary with all statistics.
        """
        # Get check report
        check_report = self.checker.check_all()

        # Group by job type
        job_type_stats = self._group_by_job_type(check_report)

        # Build summaries
        summaries = []
        for job_type, stats in job_type_stats.items():
            summary = JobTypeSummary(
                job_type=job_type,
                total=stats["total"],
                succeeded=stats["succeeded"],
                failed=stats["failed"],
                missing=stats["missing"],
            )

            # Run output validation if requested
            if validate_outputs:
                validated, failed = self._validate_job_type(
                    job_type, stats["succeeded_jobs"]
                )
                summary.validated = validated
                summary.validation_failed = failed

            summaries.append(summary)

        # Sort by job type name
        summaries.sort(key=lambda x: x.job_type)

        return PipelineSummary(
            job_types=summaries,
            validation_enabled=validate_outputs,
        )

    def _group_by_job_type(
        self, report: CheckReport
    ) -> dict[str, dict[str, Any]]:
        """Group check results by job type.

        Args:
            report: CheckReport from JobChecker.

        Returns:
            Dict mapping job type to statistics.
        """
        stats: dict[str, dict[str, Any]] = {}

        def get_job_type(job_name: str) -> str:
            """Extract job type from job name (e.g., 'rna-map-0001' -> 'rna-map')."""
            parts = job_name.rsplit("-", 1)
            if len(parts) == 2 and parts[1].isdigit():
                return parts[0]
            return job_name

        # Process succeeded jobs
        for job in report.succeeded:
            job_type = get_job_type(job.job_name)
            if job_type not in stats:
                stats[job_type] = {
                    "total": 0,
                    "succeeded": 0,
                    "failed": 0,
                    "missing": 0,
                    "succeeded_jobs": [],
                }
            stats[job_type]["total"] += 1
            stats[job_type]["succeeded"] += 1
            stats[job_type]["succeeded_jobs"].append(job)

        # Process failed jobs
        for job in report.failed:
            job_type = get_job_type(job.job_name)
            if job_type not in stats:
                stats[job_type] = {
                    "total": 0,
                    "succeeded": 0,
                    "failed": 0,
                    "missing": 0,
                    "succeeded_jobs": [],
                }
            stats[job_type]["total"] += 1
            stats[job_type]["failed"] += 1

        # Process missing jobs
        for job in report.missing:
            job_type = get_job_type(job.job_name)
            if job_type not in stats:
                stats[job_type] = {
                    "total": 0,
                    "succeeded": 0,
                    "failed": 0,
                    "missing": 0,
                    "succeeded_jobs": [],
                }
            stats[job_type]["total"] += 1
            stats[job_type]["missing"] += 1

        return stats

    def _validate_job_type(
        self,
        job_type: str,
        succeeded_jobs: list[Any],
    ) -> tuple[int, int]:
        """Validate outputs for succeeded jobs of a type.

        Args:
            job_type: The job type to validate.
            succeeded_jobs: List of JobStatus objects that succeeded.

        Returns:
            Tuple of (validated_count, failed_count).
        """
        validated = 0
        failed = 0

        for job in succeeded_jobs:
            if job.output_file:
                output_dir = job.output_file.parent
                result = self.validator.validate_job(
                    job_name=job.job_name,
                    job_type=job_type,
                    output_dir=output_dir,
                )
                if result.valid:
                    validated += 1
                else:
                    failed += 1

        return validated, failed


def generate_summary(
    job_dir: str | Path = "jobs",
    jobs_csv: str | Path | None = None,
    validate_outputs: bool = False,
    output_format: str = "table",
) -> str:
    """Generate a pipeline summary report.

    Args:
        job_dir: Directory containing job outputs.
        jobs_csv: Path to jobs.csv file.
        validate_outputs: Whether to validate job outputs.
        output_format: Output format ("table", "json").

    Returns:
        Formatted summary string.
    """
    builder = SummaryBuilder(
        job_dir=Path(job_dir),
        jobs_csv=Path(jobs_csv) if jobs_csv else None,
    )

    summary = builder.build(validate_outputs=validate_outputs)

    if output_format == "json":
        return summary.to_json()
    else:
        return summary.to_table()
