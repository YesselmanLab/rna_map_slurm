"""Tests for job generation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from rna_map_slurm.jobs.generator import (
    generate_job_list,
    group_into_batches,
    write_job_file,
)
from rna_map_slurm.jobs.rna_map_jobs import count_sequences_in_csv, is_large_construct
from rna_map_slurm.jobs.slurm import get_job_header, is_job_type_completed
from rna_map_slurm.models.config import SlurmOptions


class TestGetJobHeader:
    """Tests for get_job_header function."""

    def test_basic_header(self) -> None:
        """Test basic job header generation."""
        opts = SlurmOptions(
            name="test_job",
            time="01:00:00",
            mem_per_cpu="4GB",
            cpus_per_task=2,
        )

        header = get_job_header(opts)

        assert "#!/bin/bash" in header
        assert "#SBATCH --time=01:00:00" in header
        assert "#SBATCH --mem-per-cpu=4GB" in header
        assert "#SBATCH --cpus-per-task=2" in header

    def test_header_with_job_dir(self, temp_dir: Path) -> None:
        """Test header with output directory specified."""
        opts = SlurmOptions(name="test_job")

        header = get_job_header(opts, str(temp_dir))

        assert f"{temp_dir}/test_job.out" in header
        assert f"{temp_dir}/test_job.err" in header


class TestIsJobTypeCompleted:
    """Tests for is_job_type_completed function."""

    def test_completed_when_no_jobs(self) -> None:
        """Test that job type is completed when no jobs running."""
        jobs: list[dict[str, str]] = []
        assert is_job_type_completed("rna-map", jobs) is True

    def test_not_completed_when_running(self) -> None:
        """Test that job type is not completed when jobs running."""
        jobs = [{"Name": "rna-map-0001"}]
        assert is_job_type_completed("rna-map", jobs) is False

    def test_completed_when_other_type_running(self) -> None:
        """Test that job type is completed when other types running."""
        jobs = [{"Name": "demultiplex-0001"}]
        assert is_job_type_completed("rna-map", jobs) is True


class TestGroupIntoBatches:
    """Tests for group_into_batches function."""

    def test_even_batches(self) -> None:
        """Test grouping into even batches."""
        items = list(range(10))

        result = group_into_batches(items, 2)

        assert len(result) == 5
        assert all(len(batch) == 2 for batch in result)

    def test_uneven_batches(self) -> None:
        """Test grouping with remainder."""
        items = list(range(10))

        result = group_into_batches(items, 3)

        assert len(result) == 4
        assert len(result[-1]) == 1  # Last batch has 1 item


class TestWriteJobFile:
    """Tests for write_job_file function."""

    def test_write_job_file(self, temp_dir: Path) -> None:
        """Test writing job file."""
        content = "#!/bin/bash\necho hello"

        write_job_file(temp_dir, "test_job", content)

        job_path = temp_dir / "test_job.sh"
        assert job_path.exists()
        assert job_path.read_text() == content


class TestGenerateJobList:
    """Tests for generate_job_list function."""

    def test_generate_job_list(self) -> None:
        """Test generating job list DataFrame."""
        df = generate_job_list(
            path="jobs/rna-map",
            job_type="rna-map",
            requirement="demultiplex",
            job_names=["rna-map-0001", "rna-map-0002"],
        )

        assert len(df) == 2
        assert "job_path" in df.columns
        assert "job_type" in df.columns
        assert "job_requirement" in df.columns
        assert df["job_type"].iloc[0] == "rna-map"
        assert df["job_requirement"].iloc[0] == "demultiplex"


class TestCountSequencesInCsv:
    """Tests for count_sequences_in_csv function."""

    def test_count_sequences(self, temp_dir: Path) -> None:
        """Test counting sequences in CSV."""
        csv_path = temp_dir / "test.csv"
        df = pd.DataFrame({"name": ["seq1", "seq2", "seq3"]})
        df.to_csv(csv_path, index=False)

        count = count_sequences_in_csv(str(csv_path))

        assert count == 3

    def test_nonexistent_file(self) -> None:
        """Test with nonexistent file."""
        count = count_sequences_in_csv("/nonexistent.csv")
        assert count == 0


class TestIsLargeConstruct:
    """Tests for is_large_construct function."""

    def test_large_construct(self, temp_dir: Path) -> None:
        """Test detection of large construct."""
        csv_path = temp_dir / "large.csv"
        df = pd.DataFrame({"name": [f"seq{i}" for i in range(150)]})
        df.to_csv(csv_path, index=False)

        assert is_large_construct(str(csv_path), threshold=100) is True

    def test_small_construct(self, temp_dir: Path) -> None:
        """Test detection of small construct."""
        csv_path = temp_dir / "small.csv"
        df = pd.DataFrame({"name": [f"seq{i}" for i in range(50)]})
        df.to_csv(csv_path, index=False)

        assert is_large_construct(str(csv_path), threshold=100) is False
