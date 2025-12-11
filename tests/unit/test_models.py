"""Tests for data models."""

from __future__ import annotations

from rna_map_slurm.models.config import (
    SlurmOptions,
    WorkflowConfig,
)
from rna_map_slurm.models.fastq import FastqFile, PairedFastqFiles


class TestFastqFile:
    """Tests for FastqFile model."""

    def test_is_compressed_gz(self) -> None:
        """Test detection of gzipped files."""
        fq = FastqFile("/path/to/file.fastq.gz")
        assert fq.is_compressed() is True

    def test_is_compressed_uncompressed(self) -> None:
        """Test detection of uncompressed files."""
        fq = FastqFile("/path/to/file.fastq")
        assert fq.is_compressed() is False

    def test_is_r1(self) -> None:
        """Test R1 file detection."""
        fq = FastqFile("/path/to/sample_R1_001.fastq.gz")
        assert fq.is_r1() is True
        assert fq.is_r2() is False

    def test_is_r2(self) -> None:
        """Test R2 file detection."""
        fq = FastqFile("/path/to/sample_R2_001.fastq.gz")
        assert fq.is_r1() is False
        assert fq.is_r2() is True

    def test_neither_r1_nor_r2(self) -> None:
        """Test file that is neither R1 nor R2."""
        fq = FastqFile("/path/to/sample.fastq.gz")
        assert fq.is_r1() is False
        assert fq.is_r2() is False


class TestPairedFastqFiles:
    """Tests for PairedFastqFiles model."""

    def test_is_compressed_both_compressed(self) -> None:
        """Test paired files are both compressed."""
        pair = PairedFastqFiles(
            FastqFile("/path/R1.fastq.gz"),
            FastqFile("/path/R2.fastq.gz"),
        )
        assert pair.is_compressed() is True

    def test_is_compressed_one_uncompressed(self) -> None:
        """Test paired files where one is uncompressed."""
        pair = PairedFastqFiles(
            FastqFile("/path/R1.fastq.gz"),
            FastqFile("/path/R2.fastq"),
        )
        assert pair.is_compressed() is False


class TestSlurmOptions:
    """Tests for SlurmOptions model."""

    def test_defaults(self) -> None:
        """Test default values."""
        opts = SlurmOptions(name="test_job")
        assert opts.name == "test_job"
        assert opts.time == "12:00:00"
        assert opts.mem_per_cpu == "2GB"
        assert opts.cpus_per_task == 1
        assert opts.extra_header_cmds == ""

    def test_custom_values(self) -> None:
        """Test custom values."""
        opts = SlurmOptions(
            name="custom",
            time="06:00:00",
            mem_per_cpu="8GB",
            cpus_per_task=4,
            extra_header_cmds="#SBATCH --partition=gpu",
        )
        assert opts.time == "06:00:00"
        assert opts.mem_per_cpu == "8GB"
        assert opts.cpus_per_task == 4


class TestWorkflowConfig:
    """Tests for WorkflowConfig model."""

    def test_from_dict(self) -> None:
        """Test creating config from dictionary."""
        data = {
            "fastq_chunks": 50,
            "paths": {"tmp": "/scratch"},
            "construct_options": {"large_construct_threshold": 200},
        }
        config = WorkflowConfig.from_dict(data)
        assert config.fastq_chunks == 50
        assert config.paths.tmp == "/scratch"
        assert config.construct_options.large_construct_threshold == 200

    def test_defaults(self) -> None:
        """Test default values."""
        config = WorkflowConfig()
        assert config.fastq_chunks == 100
        assert config.construct_options.large_construct_threshold == 100
