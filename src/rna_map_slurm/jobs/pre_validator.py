"""Pre-run validation to check inputs before job submission."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
from pydantic import ValidationError

from rna_map_slurm.models.config import WorkflowConfig
from rna_map_slurm.utils.logging import get_logger

log = get_logger("jobs.pre_validator")


@dataclass
class PreRunReport:
    """Report from pre-run validation.

    Attributes:
        valid: Whether all required checks passed.
        errors: List of error messages (blocking issues).
        warnings: List of warning messages (non-blocking issues).
        checks_passed: Number of checks that passed.
        checks_failed: Number of checks that failed.
    """

    valid: bool = True
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    checks_passed: int = 0
    checks_failed: int = 0

    def add_error(self, message: str) -> None:
        """Add an error and mark as invalid."""
        self.errors.append(message)
        self.valid = False
        self.checks_failed += 1

    def add_warning(self, message: str) -> None:
        """Add a warning (doesn't affect validity)."""
        self.warnings.append(message)

    def add_pass(self) -> None:
        """Record a passed check."""
        self.checks_passed += 1

    def summary(self) -> str:
        """Get a summary string."""
        status = "PASSED" if self.valid else "FAILED"
        return (
            f"Pre-run validation: {status}\n"
            f"  Checks passed: {self.checks_passed}\n"
            f"  Checks failed: {self.checks_failed}\n"
            f"  Warnings: {len(self.warnings)}"
        )


class PreRunValidator:
    """Validate inputs and configuration before job submission."""

    def __init__(self, base_dir: Path | None = None) -> None:
        """Initialize the validator.

        Args:
            base_dir: Base directory for the workflow. Defaults to cwd.
        """
        self.base_dir = Path(base_dir) if base_dir else Path.cwd()

    def validate_all(
        self,
        jobs_csv: Path | None = None,
        config: dict[str, Any] | None = None,
    ) -> PreRunReport:
        """Run all pre-run validation checks.

        Args:
            jobs_csv: Path to jobs.csv file.
            config: Configuration dictionary (optional).

        Returns:
            PreRunReport with all check results.
        """
        report = PreRunReport()

        # Check jobs.csv
        jobs_csv_path = jobs_csv or self.base_dir / "jobs.csv"
        self._check_jobs_csv(jobs_csv_path, report)

        # Check config if provided
        if config:
            self._check_config(config, report)

        # Check directory structure
        self._check_directories(report)

        # Check job scripts exist
        if jobs_csv_path.exists():
            self._check_job_scripts(jobs_csv_path, report)

        return report

    def _check_jobs_csv(self, path: Path, report: PreRunReport) -> None:
        """Check that jobs.csv exists and has required columns."""
        if not path.exists():
            report.add_error(f"jobs.csv not found: {path}")
            return

        try:
            df = pd.read_csv(path)
            report.add_pass()

            required_cols = ["job_type", "job_path"]
            missing_cols = [c for c in required_cols if c not in df.columns]

            if missing_cols:
                report.add_error(
                    f"jobs.csv missing required columns: {missing_cols}"
                )
            else:
                report.add_pass()

            if df.empty:
                report.add_error("jobs.csv is empty - no jobs to submit")
            else:
                report.add_pass()
                log.info(f"Found {len(df)} jobs in jobs.csv")

        except Exception as e:
            report.add_error(f"Failed to read jobs.csv: {e}")

    def _check_config(self, config: dict[str, Any], report: PreRunReport) -> None:
        """Validate configuration using Pydantic models."""
        try:
            WorkflowConfig.model_validate(config)
            report.add_pass()
            log.info("Configuration validation passed")
        except ValidationError as e:
            for error in e.errors():
                loc = ".".join(str(x) for x in error["loc"])
                msg = error["msg"]
                report.add_error(f"Config error at {loc}: {msg}")

    def _check_directories(self, report: PreRunReport) -> None:
        """Check that required directories exist."""
        required_dirs = ["jobs"]
        optional_dirs = ["logs", "submits", "inputs"]

        for dir_name in required_dirs:
            dir_path = self.base_dir / dir_name
            if dir_path.exists():
                report.add_pass()
            else:
                report.add_error(f"Required directory not found: {dir_name}/")

        for dir_name in optional_dirs:
            dir_path = self.base_dir / dir_name
            if not dir_path.exists():
                report.add_warning(f"Optional directory not found: {dir_name}/")

    def _check_job_scripts(self, jobs_csv: Path, report: PreRunReport) -> None:
        """Check that all job scripts referenced in jobs.csv exist."""
        try:
            df = pd.read_csv(jobs_csv)

            if "job_path" not in df.columns:
                return

            missing_scripts = []
            for job_path in df["job_path"]:
                script_path = Path(job_path)
                if not script_path.is_absolute():
                    script_path = self.base_dir / script_path

                if not script_path.exists():
                    missing_scripts.append(job_path)

            if missing_scripts:
                report.add_error(
                    f"{len(missing_scripts)} job scripts not found. "
                    f"First missing: {missing_scripts[0]}"
                )
                if len(missing_scripts) > 1:
                    report.add_warning(
                        f"Additional missing scripts: {missing_scripts[1:5]}"
                    )
            else:
                report.add_pass()
                log.info(f"All {len(df)} job scripts exist")

        except Exception as e:
            report.add_error(f"Error checking job scripts: {e}")

    def validate_inputs_for_job_type(
        self,
        job_type: str,
        report: PreRunReport,
    ) -> None:
        """Validate inputs specific to a job type.

        Args:
            job_type: Type of job to validate inputs for.
            report: Report to add results to.
        """
        # Job-type specific input validation
        input_checks: dict[str, list[tuple[str, str]]] = {
            "split-fastq": [
                ("inputs/fastqs", "FASTQ input directory"),
            ],
            "demultiplex": [
                ("inputs/barcodes", "Barcode files"),
            ],
            "rna-map": [
                ("inputs/fastas", "FASTA reference files"),
                ("inputs/rnas", "RNA CSV files"),
            ],
        }

        if job_type not in input_checks:
            return

        for path_pattern, description in input_checks[job_type]:
            path = self.base_dir / path_pattern
            if path.exists():
                report.add_pass()
            else:
                report.add_warning(f"{description} not found: {path_pattern}")


def validate_before_run(
    jobs_csv: Path | str = "jobs.csv",
    config: dict[str, Any] | None = None,
    base_dir: Path | str | None = None,
) -> PreRunReport:
    """Convenience function for pre-run validation.

    Args:
        jobs_csv: Path to jobs.csv file.
        config: Optional configuration dictionary.
        base_dir: Base directory for the workflow.

    Returns:
        PreRunReport with validation results.
    """
    validator = PreRunValidator(base_dir=Path(base_dir) if base_dir else None)
    return validator.validate_all(
        jobs_csv=Path(jobs_csv),
        config=config,
    )
