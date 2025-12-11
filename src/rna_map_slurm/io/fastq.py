"""FASTQ file I/O operations."""

from __future__ import annotations

import glob
import os
import re
from pathlib import Path

from rna_map_slurm.models.fastq import FastqFile, PairedFastqFiles
from rna_map_slurm.utils.logging import get_logger

log = get_logger("io.fastq")

FASTQ_EXTENSIONS = re.compile(r"\.(fastq|fq)(\.gz|\.bz2|\.xz)?$")


def validate_directory(dir_path: str | Path) -> None:
    """Validate that a path exists and is a directory.

    Args:
        dir_path: Path to validate.

    Raises:
        NotADirectoryError: If path doesn't exist or isn't a directory.
    """
    if not os.path.isdir(str(dir_path)):
        raise NotADirectoryError(f"{dir_path} is not a valid directory.")


def find_fastq_files(dir_path: str | Path) -> tuple[list[str], list[str]]:
    """Find R1 and R2 FASTQ files in a directory.

    Searches for files with recognized FASTQ extensions (.fastq, .fq)
    and optional compression extensions (.gz, .bz2, .xz).

    Args:
        dir_path: Directory to search.

    Returns:
        Tuple of (R1 file paths, R2 file paths).
    """
    dir_str = str(dir_path)

    r1_paths = [
        f
        for f in glob.glob(os.path.join(dir_str, "*_R1*"))
        if FASTQ_EXTENSIONS.search(f)
    ]

    r2_paths = [
        f
        for f in glob.glob(os.path.join(dir_str, "*_R2*"))
        if FASTQ_EXTENSIONS.search(f)
    ]

    return r1_paths, r2_paths


def _extract_base_name(path: str, read_marker: str) -> str:
    """Extract base name by removing read marker and everything after.

    Args:
        path: File path.
        read_marker: Read marker to remove (e.g., "_R1_", "_R2_").

    Returns:
        Base name without read marker.
    """
    return re.sub(rf"{read_marker}.*", "", path)


def pair_fastq_files(
    r1_paths: list[str],
    r2_paths: list[str],
) -> list[PairedFastqFiles]:
    """Pair R1 and R2 FASTQ files based on matching base names.

    Args:
        r1_paths: List of R1 file paths.
        r2_paths: List of R2 file paths.

    Returns:
        List of PairedFastqFiles objects.
    """
    paired_files: list[PairedFastqFiles] = []

    for r1 in r1_paths:
        base_name = _extract_base_name(r1, "_R1_")
        matching_r2 = next(
            (r2 for r2 in r2_paths if _extract_base_name(r2, "_R2_") == base_name),
            None,
        )
        if matching_r2:
            paired_files.append(
                PairedFastqFiles(FastqFile(r1), FastqFile(matching_r2))
            )

    return paired_files


def get_paired_fastqs(dir_path: str | Path) -> list[PairedFastqFiles]:
    """Get all paired FASTQ files from a directory.

    Args:
        dir_path: Directory containing FASTQ files.

    Returns:
        List of PairedFastqFiles objects.

    Raises:
        NotADirectoryError: If dir_path is not a valid directory.
        FileNotFoundError: If no FASTQ files are found.
        ValueError: If no matching pairs are found.
    """
    validate_directory(dir_path)

    r1_paths, r2_paths = find_fastq_files(dir_path)

    if not r1_paths or not r2_paths:
        raise FileNotFoundError(f"Could not find paired fastq files in {dir_path}")

    paired_files = pair_fastq_files(r1_paths, r2_paths)

    if not paired_files:
        raise ValueError(f"Could not find matching pairs of fastq files in {dir_path}")

    return paired_files
