"""Job type definitions for RNA Map SLURM pipeline."""

from __future__ import annotations

from enum import Enum


class JobType(str, Enum):
    """All job types in the RNA Map SLURM pipeline.

    Each job type corresponds to a stage in the processing pipeline.
    The string values match the keys used in default.yml configuration.
    """

    SPLIT_FASTQ = "split-fastq"
    DEMULTIPLEX = "demultiplex"
    TRIM_GALORE = "trim-galore"
    JOIN_FASTQ_FILES = "join-fastq-files"
    RNA_MAP = "rna-map"
    RNA_MAP_COMBINE = "rna-map-combine"
    INT_DEMULTIPLEX = "int-demultiplex"
    INT_DEMULTIPLEX_RNA_MAP = "int-demultiplex-rna-map"
    INT_DEMULTIPLEX_RNA_MAP_COMBINE = "int-demultiplex-rna-map-combine"

    @classmethod
    def from_string(cls, value: str) -> JobType | None:
        """Get JobType from string value.

        Args:
            value: Job type string (e.g., "rna-map")

        Returns:
            JobType enum member or None if not found.
        """
        for member in cls:
            if member.value == value:
                return member
        return None

    @classmethod
    def all_values(cls) -> list[str]:
        """Get all job type string values."""
        return [member.value for member in cls]

    def get_dependency(self) -> JobType | None:
        """Get the job type this job depends on.

        Returns:
            The prerequisite JobType or None if this is the first stage.
        """
        dependencies = {
            JobType.SPLIT_FASTQ: None,
            JobType.DEMULTIPLEX: JobType.SPLIT_FASTQ,
            JobType.TRIM_GALORE: JobType.DEMULTIPLEX,
            JobType.JOIN_FASTQ_FILES: JobType.DEMULTIPLEX,
            JobType.RNA_MAP: JobType.JOIN_FASTQ_FILES,
            JobType.RNA_MAP_COMBINE: JobType.RNA_MAP,
            JobType.INT_DEMULTIPLEX: JobType.JOIN_FASTQ_FILES,
            JobType.INT_DEMULTIPLEX_RNA_MAP: JobType.INT_DEMULTIPLEX,
            JobType.INT_DEMULTIPLEX_RNA_MAP_COMBINE: JobType.INT_DEMULTIPLEX_RNA_MAP,
        }
        return dependencies.get(self)


# Pipeline order for display purposes
PIPELINE_ORDER = [
    JobType.SPLIT_FASTQ,
    JobType.DEMULTIPLEX,
    JobType.TRIM_GALORE,
    JobType.JOIN_FASTQ_FILES,
    JobType.RNA_MAP,
    JobType.RNA_MAP_COMBINE,
    JobType.INT_DEMULTIPLEX,
    JobType.INT_DEMULTIPLEX_RNA_MAP,
    JobType.INT_DEMULTIPLEX_RNA_MAP_COMBINE,
]
