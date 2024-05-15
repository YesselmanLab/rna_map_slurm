import os
import glob
import re
from dataclasses import dataclass
from typing import List, Tuple
import logging

from rna_map_slurm.logger import get_logger

log = get_logger(__name__)


@dataclass(frozen=True, order=True)
class FastqFile:
    """
    Holds the path to a fastq file
    """

    path: str

    def is_compressed(self) -> bool:
        """
        Check if the file is compressed
        """
        return self.path.endswith(".gz")

    def is_r1(self) -> bool:
        """
        Check if the file is R1
        """
        return bool(re.search(r"_R1_", self.path))

    def is_r2(self) -> bool:
        """
        Check if the file is R2
        """
        return bool(re.search(r"_R2_", self.path))


@dataclass(frozen=True, order=True)
class PairedFastqFiles:
    """
    Holds the paths to paired fastq files
    """

    read_1: FastqFile
    read_2: FastqFile

    def is_compressed(self) -> bool:
        """
        Check if the files are compressed
        """
        return self.read_1.is_compressed() and self.read_2.is_compressed()


def validate_directory(dir_path: str) -> None:
    """
    Validate if the provided directory exists and is a directory
    :param dir_path: path to directory
    """
    if not os.path.isdir(dir_path):
        raise NotADirectoryError(f"{dir_path} is not a valid directory.")


def find_fastq_files(dir_path: str) -> Tuple[List[str], List[str]]:
    """
    Find R1 and R2 fastq files in the given directory with common recognized extensions.
    :param dir_path: path to directory
    :return: tuple containing lists of R1 and R2 file paths
    """
    # Regular expression to match recognized FASTQ file extensions
    recognized_extensions = re.compile(r"\.(fastq|fq)(\.gz|\.bz2|\.xz)?$")

    f1_paths = [
        f
        for f in glob.glob(os.path.join(dir_path, "*_R1_*"))
        if recognized_extensions.search(f)
    ]

    f2_paths = [
        f
        for f in glob.glob(os.path.join(dir_path, "*_R2_*"))
        if recognized_extensions.search(f)
    ]

    return f1_paths, f2_paths


def pair_fastq_files(
    f1_paths: List[str], f2_paths: List[str]
) -> List[PairedFastqFiles]:
    """
    Pair R1 and R2 fastq files based on their common base names
    :param f1_paths: list of R1 file paths
    :param f2_paths: list of R2 file paths
    :return: list of PairedFastqFiles objects
    """
    paired_files = []
    for f1 in f1_paths:
        base_name = re.sub(r"_R1_.*", "", f1)
        matching_f2 = next(
            (f2 for f2 in f2_paths if re.sub(r"_R2_.*", "", f2) == base_name), None
        )
        if matching_f2:
            paired_files.append(PairedFastqFiles(FastqFile(f1), FastqFile(matching_f2)))
    return paired_files


def get_paired_fastqs(dir_path: str) -> List[PairedFastqFiles]:
    """
    Get the paired fastq files from a directory
    :param dir_path: path to directory
    :return: list of paired fastq files
    """
    validate_directory(dir_path)

    f1_paths, f2_paths = find_fastq_files(dir_path)

    if not f1_paths or not f2_paths:
        raise FileNotFoundError(f"Could not find paired fastq files in {dir_path}")

    paired_files = pair_fastq_files(f1_paths, f2_paths)

    if not paired_files:
        raise ValueError(f"Could not find matching pairs of fastq files in {dir_path}")

    return paired_files
