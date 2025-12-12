"""Shared pytest fixtures for rna_map_slurm tests."""

from __future__ import annotations

import tempfile
from collections.abc import Generator
from pathlib import Path

import pytest


@pytest.fixture
def test_resources_dir() -> Path:
    """Get path to test resources directory."""
    return Path(__file__).parent / "resources"


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_fastq_dir(test_resources_dir: Path) -> Path:
    """Get path to sample FASTQ files."""
    return test_resources_dir / "test_fastqs"


@pytest.fixture
def sample_fastq_gzipped_dir(test_resources_dir: Path) -> Path:
    """Get path to gzipped sample FASTQ files."""
    return test_resources_dir / "test_fastqs_gziped"


@pytest.fixture
def c0098_test_case_dir(test_resources_dir: Path) -> Path:
    """Get path to C0098 test case directory with real FASTQ files."""
    return test_resources_dir / "test_cases" / "C0098"


@pytest.fixture
def c0098_fastq_r1(c0098_test_case_dir: Path) -> Path:
    """Get path to C0098 R1 FASTQ file."""
    return c0098_test_case_dir / "R1.sub.fastq.gz"


@pytest.fixture
def c0098_fastq_r2(c0098_test_case_dir: Path) -> Path:
    """Get path to C0098 R2 FASTQ file."""
    return c0098_test_case_dir / "R2.sub.fastq.gz"


@pytest.fixture
def c0098_fasta(c0098_test_case_dir: Path) -> Path:
    """Get path to C0098 FASTA file."""
    return c0098_test_case_dir / "C0098.fasta"


@pytest.fixture
def c0098_csv(c0098_test_case_dir: Path) -> Path:
    """Get path to C0098 CSV file."""
    return c0098_test_case_dir / "C0098.csv"


@pytest.fixture
def c0098_barcodes_json(c0098_test_case_dir: Path) -> Path:
    """Get path to C0098 barcodes JSON file."""
    return c0098_test_case_dir / "C0098_barcodes.json"


@pytest.fixture
def c0098_data_csv(c0098_test_case_dir: Path) -> Path:
    """Get path to C0098 data CSV file for workflow setup."""
    return c0098_test_case_dir / "data.csv"


@pytest.fixture
def sample_data_csv(test_resources_dir: Path, temp_dir: Path) -> Path:
    """Create a sample data.csv file for testing."""
    csv_content = """barcode,barcode_seq,construct,code,run_name,exp_name,exp_type,data_type
BC01,ACGTACGT,RNA1,C0001,run1,test_exp,DMS,reactivity
BC02,TGCATGCA,RNA2,C0002,run1,test_exp,DMS,reactivity
"""
    csv_path = temp_dir / "data.csv"
    csv_path.write_text(csv_content)
    return csv_path


@pytest.fixture
def sample_workflow_params() -> dict[str, object]:
    """Create sample workflow parameters for testing."""
    return {
        "fastq_chunks": 10,
        "num_dirs": 10,
        "paths": {
            "log": "logs",
            "jobs": "jobs",
            "submits": "submits",
            "inputs": "inputs",
            "tmp": "/tmp",
        },
        "construct_options": {
            "large_construct_threshold": 100,
        },
        "tasks_per_job": {
            "split-fastq": 1,
            "demultiplex": 5,
            "rna-map": 10,
            "rna-map-combine": 1,
        },
        "slurm_options": {
            "extra-header-cmds": "",
            "split-fastq": {
                "time": "04:00:00",
                "cpus-per-task": 4,
                "mem-per-cpu": "8GB",
            },
            "demultiplex": {
                "time": "02:00:00",
                "cpus-per-task": 1,
                "mem-per-cpu": "2GB",
            },
            "rna-map": {
                "time": "06:00:00",
                "cpus-per-task": 1,
                "mem-per-cpu": "2GB",
            },
            "rna-map-combine": {
                "time": "02:00:00",
                "cpus-per-task": 1,
                "mem-per-cpu": "2GB",
            },
        },
    }
