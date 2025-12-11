"""Tests for FASTQ I/O operations."""

from __future__ import annotations

from pathlib import Path

import pytest

from rna_map_slurm.io.fastq import (
    find_fastq_files,
    get_paired_fastqs,
    pair_fastq_files,
    validate_directory,
)


class TestValidateDirectory:
    """Tests for validate_directory function."""

    def test_valid_directory(self, temp_dir: Path) -> None:
        """Test with valid directory."""
        # Should not raise
        validate_directory(temp_dir)

    def test_invalid_directory(self) -> None:
        """Test with non-existent directory."""
        with pytest.raises(NotADirectoryError):
            validate_directory("/nonexistent/path")

    def test_file_not_directory(self, temp_dir: Path) -> None:
        """Test with file instead of directory."""
        file_path = temp_dir / "test.txt"
        file_path.write_text("test")
        with pytest.raises(NotADirectoryError):
            validate_directory(file_path)


class TestFindFastqFiles:
    """Tests for find_fastq_files function."""

    def test_find_fastq_files(self, temp_dir: Path) -> None:
        """Test finding FASTQ files."""
        # Create test files
        (temp_dir / "sample_R1_001.fastq").touch()
        (temp_dir / "sample_R2_001.fastq").touch()
        (temp_dir / "sample_R1_001.fastq.gz").touch()
        (temp_dir / "sample_R2_001.fastq.gz").touch()

        r1_paths, r2_paths = find_fastq_files(temp_dir)

        assert len(r1_paths) == 2
        assert len(r2_paths) == 2

    def test_no_fastq_files(self, temp_dir: Path) -> None:
        """Test empty directory."""
        r1_paths, r2_paths = find_fastq_files(temp_dir)
        assert len(r1_paths) == 0
        assert len(r2_paths) == 0


class TestPairFastqFiles:
    """Tests for pair_fastq_files function."""

    def test_pair_matching_files(self) -> None:
        """Test pairing matching R1/R2 files."""
        r1_paths = [
            "/path/sample_R1_001.fastq",
            "/path/sample_R1_002.fastq",
        ]
        r2_paths = [
            "/path/sample_R2_001.fastq",
            "/path/sample_R2_002.fastq",
        ]

        pairs = pair_fastq_files(r1_paths, r2_paths)

        assert len(pairs) == 2
        assert pairs[0].read_1.is_r1()
        assert pairs[0].read_2.is_r2()

    def test_no_matching_pairs(self) -> None:
        """Test with no matching pairs."""
        r1_paths = ["/path/A_R1_001.fastq"]
        r2_paths = ["/path/B_R2_001.fastq"]

        pairs = pair_fastq_files(r1_paths, r2_paths)

        assert len(pairs) == 0


class TestGetPairedFastqs:
    """Tests for get_paired_fastqs function."""

    def test_get_paired_fastqs(self, temp_dir: Path) -> None:
        """Test getting paired FASTQ files."""
        (temp_dir / "sample_R1_001.fastq").touch()
        (temp_dir / "sample_R2_001.fastq").touch()

        pairs = get_paired_fastqs(temp_dir)

        assert len(pairs) == 1
        assert pairs[0].read_1.is_r1()
        assert pairs[0].read_2.is_r2()

    def test_invalid_directory(self) -> None:
        """Test with invalid directory."""
        with pytest.raises(NotADirectoryError):
            get_paired_fastqs("/nonexistent")

    def test_no_fastq_files(self, temp_dir: Path) -> None:
        """Test with no FASTQ files."""
        with pytest.raises(FileNotFoundError):
            get_paired_fastqs(temp_dir)
