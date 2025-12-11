"""Configuration data models for SLURM jobs."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class SlurmOptions:
    """SLURM job configuration options.

    Attributes:
        name: Job name for SLURM scheduler.
        time: Maximum wall time (e.g., "12:00:00").
        mem_per_cpu: Memory per CPU (e.g., "2GB").
        cpus_per_task: Number of CPUs per task.
        extra_header_cmds: Additional SBATCH directives.
    """

    name: str
    time: str = "12:00:00"
    mem_per_cpu: str = "2GB"
    cpus_per_task: int = 1
    extra_header_cmds: str = ""


@dataclass
class JobConfig:
    """Configuration for a batch of jobs.

    Attributes:
        job_type: Type of job (e.g., "rna-map", "demultiplex").
        runs_per_job: Number of tasks to batch per SLURM job.
        slurm_options: SLURM-specific configuration.
        requirement: Job dependency (job type that must complete first).
    """

    job_type: str
    runs_per_job: int
    slurm_options: SlurmOptions
    requirement: str = ""


@dataclass
class PathConfig:
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


@dataclass
class ConstructOptions:
    """Options for construct handling.

    Attributes:
        large_construct_threshold: Constructs with more sequences than this
            get their own dedicated job instead of being batched.
    """

    large_construct_threshold: int = 100


@dataclass
class WorkflowConfig:
    """Complete workflow configuration.

    Attributes:
        fastq_chunks: Number of chunks to split FASTQ files into.
        num_dirs: Total number of data directories.
        paths: Path configuration.
        construct_options: Construct handling options.
        tasks_per_job: Mapping of job type to tasks per job.
        slurm_options: Mapping of job type to SLURM options.
    """

    fastq_chunks: int = 100
    num_dirs: int = 0
    paths: PathConfig = field(default_factory=PathConfig)
    construct_options: ConstructOptions = field(default_factory=ConstructOptions)
    tasks_per_job: dict[str, int] = field(default_factory=dict)
    slurm_options: dict[str, dict[str, Any]] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> WorkflowConfig:
        """Create WorkflowConfig from a dictionary.

        Args:
            data: Configuration dictionary (typically from YAML).

        Returns:
            WorkflowConfig instance.
        """
        paths = PathConfig(**data.get("paths", {}))
        construct_options = ConstructOptions(**data.get("construct_options", {}))

        return cls(
            fastq_chunks=data.get("fastq_chunks", 100),
            num_dirs=data.get("num_dirs", 0),
            paths=paths,
            construct_options=construct_options,
            tasks_per_job=data.get("tasks_per_job", {}),
            slurm_options=data.get("slurm_options", {}),
        )
