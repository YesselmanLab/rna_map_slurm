"""Data models for FASTQ file handling."""

from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class FastqFile:
    """Represents a FASTQ file with path and metadata.

    Attributes:
        path: Absolute path to the FASTQ file.
    """

    path: str

    def is_compressed(self) -> bool:
        """Check if the file is gzip compressed.

        Returns:
            True if file ends with .gz, False otherwise.
        """
        return self.path.endswith(".gz")

    def is_r1(self) -> bool:
        """Check if this is an R1 (forward read) file.

        Returns:
            True if filename contains _R1_, False otherwise.
        """
        return bool(re.search(r"_R1_", self.path))

    def is_r2(self) -> bool:
        """Check if this is an R2 (reverse read) file.

        Returns:
            True if filename contains _R2_, False otherwise.
        """
        return bool(re.search(r"_R2_", self.path))


@dataclass(frozen=True, order=True)
class PairedFastqFiles:
    """Represents a pair of R1 and R2 FASTQ files.

    Attributes:
        read_1: The R1 (forward read) FASTQ file.
        read_2: The R2 (reverse read) FASTQ file.
    """

    read_1: FastqFile
    read_2: FastqFile

    def is_compressed(self) -> bool:
        """Check if both files are gzip compressed.

        Returns:
            True if both R1 and R2 files are compressed, False otherwise.
        """
        return self.read_1.is_compressed() and self.read_2.is_compressed()
