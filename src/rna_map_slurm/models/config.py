"""Configuration data models for SLURM jobs using Pydantic for validation."""

from __future__ import annotations

import re
from typing import Any

from pydantic import BaseModel, Field, field_validator, model_validator


class SlurmOptions(BaseModel):
    """SLURM job configuration options with validation.

    Attributes:
        name: Job name for SLURM scheduler.
        time: Maximum wall time (e.g., "12:00:00").
        mem_per_cpu: Memory per CPU (e.g., "2GB", "4000MB").
        cpus_per_task: Number of CPUs per task.
        extra_header_cmds: Additional SBATCH directives.
    """

    name: str = ""
    time: str = "12:00:00"
    mem_per_cpu: str = "2GB"
    cpus_per_task: int = Field(default=1, ge=1)
    extra_header_cmds: str = ""

    model_config = {"frozen": True}

    @field_validator("time")
    @classmethod
    def validate_time_format(cls, v: str) -> str:
        """Validate SLURM time format (HH:MM:SS or D-HH:MM:SS)."""
        patterns = [
            r"^\d{1,2}:\d{2}:\d{2}$",  # HH:MM:SS
            r"^\d+-\d{1,2}:\d{2}:\d{2}$",  # D-HH:MM:SS
            r"^\d+$",  # minutes
        ]
        if not any(re.match(p, v) for p in patterns):
            raise ValueError(
                f"Invalid time format: {v}. Use HH:MM:SS, D-HH:MM:SS, or minutes"
            )
        return v

    @field_validator("mem_per_cpu")
    @classmethod
    def validate_memory_format(cls, v: str) -> str:
        """Validate memory format (e.g., 2GB, 4000MB, 4000)."""
        if not re.match(r"^\d+(\.\d+)?(GB|MB|G|M)?$", v, re.IGNORECASE):
            raise ValueError(
                f"Invalid memory format: {v}. Use format like 2GB, 4000MB, or 4000"
            )
        return v


class PathConfig(BaseModel):
    """Path configuration for the workflow.

    Attributes:
        log: Directory for log files.
        jobs: Directory for job scripts.
        submits: Directory for submit files.
        inputs: Directory for input files.
        tmp: Temporary directory for scratch space.
        seq_path: Path to sequence data.
    """

    log: str = "logs"
    jobs: str = "jobs"
    submits: str = "submits"
    inputs: str = "inputs"
    tmp: str = "/scratch"
    seq_path: str = ""


class ConstructOptions(BaseModel):
    """Options for construct handling.

    Attributes:
        large_construct_threshold: Constructs with more sequences than this
            get their own dedicated job instead of being batched.
    """

    large_construct_threshold: int = Field(default=100, ge=1)


class WorkflowConfig(BaseModel):
    """Complete workflow configuration with validation.

    Attributes:
        fastq_chunks: Number of chunks to split FASTQ files into.
        num_dirs: Total number of data directories.
        paths: Path configuration.
        construct_options: Construct handling options.
        tasks_per_job: Mapping of job type to tasks per job.
        slurm_options: Mapping of job type to SLURM options.
        use_cpp_demultiplex: Use C++ batch mode for internal demultiplexing (218x faster).
    """

    fastq_chunks: int = Field(default=100, ge=1)
    num_dirs: int = Field(default=0, ge=0)
    paths: PathConfig = Field(default_factory=PathConfig)
    construct_options: ConstructOptions = Field(default_factory=ConstructOptions)
    tasks_per_job: dict[str, int] = Field(default_factory=dict)
    slurm_options: dict[str, dict[str, Any]] = Field(default_factory=dict)
    use_cpp_demultiplex: bool = Field(default=False)

    @field_validator("tasks_per_job")
    @classmethod
    def validate_tasks_per_job(cls, v: dict[str, int]) -> dict[str, int]:
        """Validate tasks_per_job values are positive."""
        for job_type, count in v.items():
            if count < 1:
                raise ValueError(
                    f"tasks_per_job[{job_type}] must be >= 1, got {count}"
                )
        return v

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> WorkflowConfig:
        """Create WorkflowConfig from a dictionary.

        Args:
            data: Configuration dictionary (typically from YAML).

        Returns:
            WorkflowConfig instance.

        Raises:
            ValidationError: If configuration is invalid.
        """
        return cls.model_validate(data)

    def get_slurm_options(self, job_type: str) -> SlurmOptions:
        """Get SLURM options for a specific job type.

        Args:
            job_type: The job type to get options for.

        Returns:
            SlurmOptions for the job type, with defaults if not specified.
        """
        opts = self.slurm_options.get(job_type, {})
        return SlurmOptions(name=job_type, **opts)


# Keep backward compatibility with dataclass-style usage
JobConfig = None  # Removed - use WorkflowConfig.get_slurm_options() instead
