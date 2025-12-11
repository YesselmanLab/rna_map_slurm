"""Data models and dataclasses."""

from rna_map_slurm.models.config import JobConfig, SlurmOptions
from rna_map_slurm.models.fastq import FastqFile, PairedFastqFiles

__all__ = ["FastqFile", "JobConfig", "PairedFastqFiles", "SlurmOptions"]
