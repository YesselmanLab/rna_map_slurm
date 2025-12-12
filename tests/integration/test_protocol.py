"""Integration tests for the full RNA map SLURM protocol using C0098 example data."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pandas as pd
import pytest

# Check if optional dependencies are available
try:
    from rna_map_slurm.tasks.basic import demultiplex, split_fastq_file

    HAS_FULL_DEPENDENCIES = True
except ImportError:
    HAS_FULL_DEPENDENCIES = False
    demultiplex = None  # type: ignore[assignment, misc]
    split_fastq_file = None  # type: ignore[assignment, misc]

requires_full_deps = pytest.mark.skipif(
    not HAS_FULL_DEPENDENCIES,
    reason="Requires full dependencies (rna_map, fastqsplitter, etc.)",
)


@requires_full_deps
class TestSplitFastqWithRealData:
    """Test FASTQ splitting with real C0098 example data."""

    def test_split_fastq_creates_chunks(
        self,
        c0098_fastq_r1: Path,
        c0098_fastq_r2: Path,
        temp_dir: Path,
    ) -> None:
        """Test that split_fastq_file creates the expected number of chunks."""
        num_chunks = 2

        # Create output directories
        for i in range(num_chunks):
            (temp_dir / f"split-{i:04}").mkdir()

        # Split R1
        r1_outputs = split_fastq_file(
            str(c0098_fastq_r1),
            str(temp_dir),
            num_chunks,
            start=0,
            threads=1,
        )

        # Split R2
        r2_outputs = split_fastq_file(
            str(c0098_fastq_r2),
            str(temp_dir),
            num_chunks,
            start=0,
            threads=1,
        )

        assert len(r1_outputs) == num_chunks
        assert len(r2_outputs) == num_chunks

        # Verify files were created
        for i in range(num_chunks):
            r1_path = temp_dir / f"split-{i:04}" / "test_R1.fastq.gz"
            r2_path = temp_dir / f"split-{i:04}" / "test_R2.fastq.gz"
            assert r1_path.exists(), f"R1 chunk {i} not created"
            assert r2_path.exists(), f"R2 chunk {i} not created"
            assert r1_path.stat().st_size > 0, f"R1 chunk {i} is empty"
            assert r2_path.stat().st_size > 0, f"R2 chunk {i} is empty"


@requires_full_deps
class TestDemultiplexWithRealData:
    """Test demultiplexing with real C0098 example data."""

    @pytest.fixture
    def demultiplex_csv(self, temp_dir: Path) -> Path:
        """Create a sample CSV for demultiplexing."""
        csv_content = """barcode,barcode_seq,construct
BC01,ACGTACGT,test_construct
"""
        csv_path = temp_dir / "barcodes.csv"
        csv_path.write_text(csv_content)
        return csv_path

    def test_demultiplex_creates_barcode_dirs(
        self,
        c0098_fastq_r1: Path,
        c0098_fastq_r2: Path,
        demultiplex_csv: Path,
        temp_dir: Path,
    ) -> None:
        """Test that demultiplex creates directories for each barcode."""
        output_dir = temp_dir / "demux_output"
        output_dir.mkdir()

        # Run demultiplexing
        output_dirs = demultiplex(
            str(demultiplex_csv),
            str(c0098_fastq_r1),
            str(c0098_fastq_r2),
            str(output_dir),
        )

        assert len(output_dirs) == 1
        assert "ACGTACGT" in output_dirs[0]


class TestProtocolSetup:
    """Test the setup phase of the protocol."""

    def test_can_read_c0098_data_files(
        self,
        c0098_test_case_dir: Path,
        c0098_csv: Path,
        c0098_fasta: Path,
    ) -> None:
        """Test that C0098 test data files are readable."""
        # Check CSV is readable
        df = pd.read_csv(c0098_csv)
        assert len(df) > 0
        assert "name" in df.columns
        assert "sequence" in df.columns

        # Check FASTA is readable
        fasta_content = c0098_fasta.read_text()
        assert fasta_content.startswith(">")
        assert "CAUGG_CCUAAA" in fasta_content

    def test_c0098_fastq_files_exist_and_are_valid(
        self,
        c0098_fastq_r1: Path,
        c0098_fastq_r2: Path,
    ) -> None:
        """Test that C0098 FASTQ files exist and have content."""
        assert c0098_fastq_r1.exists()
        assert c0098_fastq_r2.exists()
        assert c0098_fastq_r1.stat().st_size > 0
        assert c0098_fastq_r2.stat().st_size > 0

    def test_data_csv_has_required_columns(
        self,
        c0098_test_case_dir: Path,
    ) -> None:
        """Test that data.csv has the required columns for the protocol."""
        data_csv = c0098_test_case_dir / "data.csv"
        df = pd.read_csv(data_csv)

        required_columns = [
            "run_name",
            "exp_name",
            "code",
            "construct",
            "barcode_seq",
            "exp_type",
            "data_type",
        ]

        for col in required_columns:
            assert col in df.columns, f"Missing required column: {col}"
