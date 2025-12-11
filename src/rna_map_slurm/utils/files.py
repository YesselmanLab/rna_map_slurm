"""File manipulation utilities."""

from __future__ import annotations

import gzip
import os
import random
import shutil
import string
import zipfile
from pathlib import Path

from rna_map_slurm.utils.logging import get_logger

log = get_logger("files")


def get_file_size(file_path: str | Path) -> int:
    """Get the size of a file in bytes.

    Args:
        file_path: Path to the file.

    Returns:
        File size in bytes.
    """
    resolved_path = os.path.realpath(str(file_path))
    return os.path.getsize(resolved_path)


def random_string(length: int) -> str:
    """Generate a random string of ASCII letters.

    Args:
        length: Desired length of the random string.

    Returns:
        Random string of specified length.
    """
    return "".join(random.choices(string.ascii_letters, k=length))


def gzip_files(directory: str | Path) -> None:
    """Compress all uncompressed files in a directory with gzip.

    Args:
        directory: Directory containing files to compress.

    Note:
        Original files are removed after compression.
        Files already ending in .gz are skipped.
    """
    for root, _dirs, files in os.walk(str(directory)):
        for file in files:
            if file.endswith(".gz"):
                continue
            file_path = Path(root) / file
            compressed_path = Path(f"{file_path}.gz")
            with (
                open(file_path, "rb") as f_in,
                gzip.open(compressed_path, "wb") as f_out,
            ):
                shutil.copyfileobj(f_in, f_out)
            os.remove(file_path)


def flatten_and_zip_directory(input_directory: str | Path, output_zip: str | Path) -> None:
    """Create a flat zip archive from a directory tree.

    All files are placed at the root level of the zip, regardless of
    their original directory structure.

    Args:
        input_directory: Directory to compress.
        output_zip: Path for output zip file.
    """
    with zipfile.ZipFile(str(output_zip), "w") as zip_ref:
        for root, _, files in os.walk(str(input_directory)):
            for file in files:
                file_path = Path(root) / file
                zip_ref.write(file_path, file)
