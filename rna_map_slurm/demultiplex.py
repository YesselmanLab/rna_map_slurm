import os
import subprocess
import pandas as pd
from typing import List
from tabulate import tabulate

from rna_map_slurm.fastq import PairedFastqFiles
from rna_map_slurm.logger import get_logger

log = get_logger(__name__)


class SabreDemultiplexer:
    def run(
        self, df: pd.DataFrame, paired_fqs: PairedFastqFiles, demultiplex_path: str
    ) -> None:
        """
        Run the demultiplexing process.

        Args:
            df (pd.DataFrame): DataFrame containing barcode information.
            paired_fqs (PairedFastqFiles): PairedFastqFiles object containing paths to the FASTQ files.
            demultiplex_path (str): Path to the directory where demultiplexing should be performed.
        """
        if not os.path.isdir(demultiplex_path):
            log.error(f"{demultiplex_path} does not exist")
            return

        log.info("Preparing barcodes.txt file for demultiplexing")
        self._generate_barcode_file(df)

        r1_path = paired_fqs.read_1.path
        r2_path = paired_fqs.read_2.path
        command = (
            f"sabre pe -f {r1_path} -r {r2_path} -b barcode.txt "
            f"-u NC/test_R1.fastq -w NC/test_R2.fastq -m 4"
        )
        log.info(f"Running sabre with command: {command}")

        try:
            output = subprocess.check_output(command, shell=True)
            output = output.decode("UTF-8")
            log.info(f"Output from sabre:\n{output}")
        except subprocess.CalledProcessError as e:
            log.error(f"Error running sabre: {e.output.decode('UTF-8')}")
            return

        for _, row in df.iterrows():
            self._gzip_files(row["barcode_seq"])

    def _generate_barcode_file(
        self, df: pd.DataFrame, fname: str = "barcode.txt"
    ) -> None:
        """
        Generate barcode file for sabre demultiplexing.

        Args:
            df (pd.DataFrame): DataFrame containing barcode information.
            fname (str): Filename for the barcode file.
        """
        expects = ["barcode", "barcode_seq", "construct"]
        self._check_if_columns_exist(df, expects)

        seen = set()
        warning = False
        log.info(
            "Constructs:\n\n"
            + tabulate(df[expects], expects, tablefmt="github", showindex=False)
            + "\n"
        )

        lines = []
        for _, row in df.iterrows():
            barcode = row["barcode"]
            barcode_seq = row["barcode_seq"]
            if barcode in seen:
                log.warning(
                    f"{barcode} has been used more than once; this may be an issue."
                )
                warning = True
                continue

            line = f"{barcode_seq}\t{barcode_seq}/test_R1.fastq\t{barcode_seq}/test_R2.fastq"
            lines.append(line)
            os.makedirs(barcode_seq, exist_ok=True)
            seen.add(barcode)

        os.makedirs("NC", exist_ok=True)
        log.info(f"{len(seen)} unique barcodes found in the CSV file.")
        if not warning:
            log.info("No barcode conflicts detected.")

        with open(fname, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    def _gzip_files(self, barcode_seq: str) -> None:
        """
        Gzip the demultiplexed files for a given barcode sequence.

        Args:
            barcode_seq (str): Barcode sequence for which to gzip files.
        """
        try:
            subprocess.check_call(f"gzip {barcode_seq}/test_R1.fastq", shell=True)
            subprocess.check_call(f"gzip {barcode_seq}/test_R2.fastq", shell=True)
            log.info(f"Gzipped files for barcode sequence: {barcode_seq}")
        except subprocess.CalledProcessError as e:
            log.error(f"Error gzipping files for barcode sequence {barcode_seq}: {e}")

    def _check_if_columns_exist(self, df: pd.DataFrame, columns: List[str]) -> None:
        """
        Check if required columns exist in the DataFrame.

        Args:
            df (pd.DataFrame): DataFrame to check.
            columns (List[str]): List of required column names.

        Raises:
            ValueError: If any required columns are missing from the DataFrame.
        """
        missing_columns = [col for col in columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")
