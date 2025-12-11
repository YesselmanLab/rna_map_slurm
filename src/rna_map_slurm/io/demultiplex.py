"""Demultiplexing I/O operations using sabre."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd
from tabulate import tabulate

from rna_map_slurm.models.fastq import PairedFastqFiles
from rna_map_slurm.utils.logging import get_logger

log = get_logger("io.demultiplex")

REQUIRED_COLUMNS = ["barcode", "barcode_seq", "construct"]


class SabreDemultiplexer:
    """Demultiplexer using the sabre tool for paired-end reads."""

    def run(
        self,
        df: pd.DataFrame,
        paired_fqs: PairedFastqFiles,
        demultiplex_path: str | Path,
    ) -> None:
        """Run the demultiplexing process.

        Args:
            df: DataFrame containing barcode information.
            paired_fqs: Paired FASTQ files to demultiplex.
            demultiplex_path: Output directory for demultiplexed files.
        """
        if not os.path.isdir(str(demultiplex_path)):
            log.error(f"{demultiplex_path} does not exist")
            return

        log.info("Preparing barcodes.txt file for demultiplexing")
        self._generate_barcode_file(df)
        self._run_sabre(paired_fqs)
        self._compress_outputs(df)

    def _generate_barcode_file(
        self,
        df: pd.DataFrame,
        fname: str = "barcode.txt",
    ) -> None:
        """Generate barcode file for sabre demultiplexing.

        Args:
            df: DataFrame containing barcode information.
            fname: Output filename for barcode file.

        Raises:
            ValueError: If required columns are missing.
        """
        self._validate_columns(df)
        self._log_constructs(df)
        lines = self._build_barcode_lines(df)
        self._write_barcode_file(fname, lines)

    def _validate_columns(self, df: pd.DataFrame) -> None:
        """Validate required columns exist in DataFrame.

        Args:
            df: DataFrame to validate.

        Raises:
            ValueError: If required columns are missing.
        """
        missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")

    def _log_constructs(self, df: pd.DataFrame) -> None:
        """Log construct information in table format.

        Args:
            df: DataFrame with construct information.
        """
        table = tabulate(
            df[REQUIRED_COLUMNS].values.tolist(),
            REQUIRED_COLUMNS,
            tablefmt="github",
            showindex=False,
        )
        log.info(f"Constructs:\n\n{table}\n")

    def _build_barcode_lines(self, df: pd.DataFrame) -> list[str]:
        """Build barcode file lines and create output directories.

        Args:
            df: DataFrame with barcode information.

        Returns:
            List of barcode file lines.
        """
        seen: set[str] = set()
        lines: list[str] = []

        for _, row in df.iterrows():
            barcode = str(row["barcode"])
            barcode_seq = str(row["barcode_seq"])

            if barcode in seen:
                log.warning(f"{barcode} has been used more than once; this may be an issue.")
                continue

            line = f"{barcode_seq}\t{barcode_seq}/test_R1.fastq\t{barcode_seq}/test_R2.fastq"
            lines.append(line)
            os.makedirs(barcode_seq, exist_ok=True)
            seen.add(barcode)

        os.makedirs("NC", exist_ok=True)
        log.info(f"{len(seen)} unique barcodes found in the CSV file.")

        return lines

    def _write_barcode_file(self, fname: str, lines: list[str]) -> None:
        """Write barcode lines to file.

        Args:
            fname: Output filename.
            lines: Lines to write.
        """
        with open(fname, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    def _run_sabre(self, paired_fqs: PairedFastqFiles) -> None:
        """Execute sabre demultiplexing command.

        Args:
            paired_fqs: Paired FASTQ files to demultiplex.
        """
        r1_path = paired_fqs.read_1.path
        r2_path = paired_fqs.read_2.path
        command = (
            f"sabre pe -f {r1_path} -r {r2_path} -b barcode.txt "
            f"-u NC/test_R1.fastq -w NC/test_R2.fastq -m 2"
        )
        log.info(f"Running sabre with command: {command}")

        try:
            result = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
            )
            log.info(f"Output from sabre:\n{result.stdout}")
        except subprocess.CalledProcessError as e:
            log.error(f"Error running sabre: {e.stderr}")

    def _compress_outputs(self, df: pd.DataFrame) -> None:
        """Compress demultiplexed FASTQ files.

        Args:
            df: DataFrame with barcode information.
        """
        for _, row in df.iterrows():
            self._gzip_barcode_files(str(row["barcode_seq"]))

    def _gzip_barcode_files(self, barcode_seq: str) -> None:
        """Gzip demultiplexed files for a barcode.

        Args:
            barcode_seq: Barcode sequence identifier.
        """
        try:
            subprocess.run(
                f"gzip {barcode_seq}/test_R1.fastq",
                shell=True,
                check=True,
            )
            subprocess.run(
                f"gzip {barcode_seq}/test_R2.fastq",
                shell=True,
                check=True,
            )
            log.info(f"Gzipped files for barcode sequence: {barcode_seq}")
        except subprocess.CalledProcessError as e:
            log.error(f"Error gzipping files for barcode sequence {barcode_seq}: {e}")
