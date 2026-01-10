"""Output validation for completed SLURM jobs."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from rna_map_slurm.utils.logging import get_logger

log = get_logger("jobs.validation")


@dataclass
class ExpectedOutput:
    """Definition of an expected output file.

    Attributes:
        pattern: Glob or regex pattern for the file.
        required: Whether this output is required for success.
        check_non_empty: Whether to verify the file is non-empty.
        description: Human-readable description of the output.
    """

    pattern: str
    required: bool = True
    check_non_empty: bool = True
    description: str = ""


@dataclass
class OutputCheck:
    """Result of checking a single expected output.

    Attributes:
        pattern: The pattern that was checked.
        found: Whether matching files were found.
        files: List of matching file paths.
        empty_files: List of files that were empty.
        valid: Whether this output check passed.
        message: Description of the result.
    """

    pattern: str
    found: bool = False
    files: list[Path] = field(default_factory=list)
    empty_files: list[Path] = field(default_factory=list)
    valid: bool = False
    message: str = ""


@dataclass
class ValidationResult:
    """Result of validating all outputs for a job.

    Attributes:
        job_name: Name of the job that was validated.
        job_type: Type of job.
        output_checks: Results for each expected output.
        valid: Whether all required outputs were found.
        error_summary: Summary of validation errors.
    """

    job_name: str
    job_type: str
    output_checks: list[OutputCheck] = field(default_factory=list)
    valid: bool = True
    error_summary: str = ""

    @property
    def missing_outputs(self) -> list[str]:
        """Get list of missing required output patterns."""
        return [c.pattern for c in self.output_checks if not c.valid]


# Expected outputs for each job type
# Patterns use {var} for variables that will be substituted
OUTPUT_SPECS: dict[str, list[ExpectedOutput]] = {
    "split-fastq": [
        ExpectedOutput(
            pattern="test_R1.fastq.gz",
            description="Split R1 FASTQ file",
        ),
        ExpectedOutput(
            pattern="test_R2.fastq.gz",
            description="Split R2 FASTQ file",
        ),
    ],
    "demultiplex": [
        ExpectedOutput(
            pattern="*/test_R1.fastq.gz",
            description="Demultiplexed R1 FASTQ files",
        ),
        ExpectedOutput(
            pattern="*/test_R2.fastq.gz",
            description="Demultiplexed R2 FASTQ files",
        ),
    ],
    "join-fastq-files": [
        ExpectedOutput(
            pattern="test_R1.fastq.gz",
            description="Joined R1 FASTQ file",
        ),
        ExpectedOutput(
            pattern="test_R2.fastq.gz",
            description="Joined R2 FASTQ file",
        ),
    ],
    "rna-map": [
        ExpectedOutput(
            pattern="output/BitVector_Files/mutation_histos.p",
            description="Mutation histogram pickle file",
        ),
    ],
    "rna-map-combine": [
        ExpectedOutput(
            pattern="output/BitVector_Files/mutation_histos.json",
            description="Combined mutation histogram JSON",
        ),
        ExpectedOutput(
            pattern="output/BitVector_Files/mutation_histos.p",
            description="Combined mutation histogram pickle",
            required=False,
        ),
    ],
    "int-demultiplex": [
        ExpectedOutput(
            pattern="*_mate1.fastq.gz",
            description="Internal demultiplex mate1 files",
        ),
        ExpectedOutput(
            pattern="*_mate2.fastq.gz",
            description="Internal demultiplex mate2 files",
        ),
    ],
    "int-demultiplex-rna-map": [
        ExpectedOutput(
            pattern="mutation_histos_*.p",
            description="Internal demultiplex RNA-map histograms",
        ),
    ],
    "int-demultiplex-rna-map-combine": [
        ExpectedOutput(
            pattern="output/BitVector_Files/mutation_histos.json",
            description="Combined internal demultiplex histograms",
        ),
    ],
}


class OutputValidator:
    """Validate expected outputs for completed jobs."""

    def __init__(self, base_dir: Path | None = None) -> None:
        """Initialize the validator.

        Args:
            base_dir: Base directory for output files. Defaults to cwd.
        """
        self.base_dir = Path(base_dir) if base_dir else Path.cwd()

    def validate_job(
        self,
        job_name: str,
        job_type: str,
        output_dir: Path | None = None,
    ) -> ValidationResult:
        """Validate outputs for a single job.

        Args:
            job_name: Name of the job.
            job_type: Type of job (must be in OUTPUT_SPECS).
            output_dir: Directory containing job outputs.

        Returns:
            ValidationResult with check details.
        """
        if job_type not in OUTPUT_SPECS:
            return ValidationResult(
                job_name=job_name,
                job_type=job_type,
                valid=True,  # Unknown job types pass by default
                error_summary=f"No output specs defined for job type: {job_type}",
            )

        specs = OUTPUT_SPECS[job_type]
        search_dir = output_dir or self.base_dir

        result = ValidationResult(job_name=job_name, job_type=job_type)
        errors = []

        for spec in specs:
            check = self._check_output(search_dir, spec)
            result.output_checks.append(check)

            if spec.required and not check.valid:
                errors.append(f"Missing: {spec.description} ({spec.pattern})")
                result.valid = False

        result.error_summary = "; ".join(errors) if errors else ""
        return result

    def _check_output(self, search_dir: Path, spec: ExpectedOutput) -> OutputCheck:
        """Check for a single expected output.

        Args:
            search_dir: Directory to search in.
            spec: Expected output specification.

        Returns:
            OutputCheck with results.
        """
        check = OutputCheck(pattern=spec.pattern)

        try:
            # Use glob to find matching files
            matches = list(search_dir.glob(spec.pattern))

            if not matches:
                # Try recursive glob
                matches = list(search_dir.glob(f"**/{spec.pattern}"))

            check.files = matches
            check.found = len(matches) > 0

            if not check.found:
                check.message = "No matching files found"
                check.valid = not spec.required
                return check

            # Check for empty files if required
            if spec.check_non_empty:
                for f in matches:
                    if f.is_file() and f.stat().st_size == 0:
                        check.empty_files.append(f)

                if check.empty_files:
                    check.message = f"{len(check.empty_files)} empty file(s)"
                    check.valid = not spec.required
                    return check

            check.valid = True
            check.message = f"Found {len(matches)} file(s)"

        except Exception as e:
            check.message = f"Error checking: {e}"
            check.valid = False

        return check

    def validate_all_jobs(
        self,
        jobs: list[tuple[str, str, Path]],
    ) -> list[ValidationResult]:
        """Validate outputs for multiple jobs.

        Args:
            jobs: List of (job_name, job_type, output_dir) tuples.

        Returns:
            List of ValidationResult objects.
        """
        results = []
        for job_name, job_type, output_dir in jobs:
            result = self.validate_job(job_name, job_type, output_dir)
            results.append(result)
        return results

    def get_failed_validations(
        self, results: list[ValidationResult]
    ) -> list[ValidationResult]:
        """Filter to only failed validations.

        Args:
            results: List of validation results.

        Returns:
            List of failed ValidationResult objects.
        """
        return [r for r in results if not r.valid]


def validate_job_outputs(
    job_type: str,
    output_dir: Path,
    job_name: str = "",
) -> ValidationResult:
    """Convenience function to validate a single job's outputs.

    Args:
        job_type: Type of job.
        output_dir: Directory containing outputs.
        job_name: Optional job name for reporting.

    Returns:
        ValidationResult.
    """
    validator = OutputValidator()
    return validator.validate_job(
        job_name=job_name or output_dir.name,
        job_type=job_type,
        output_dir=output_dir,
    )
