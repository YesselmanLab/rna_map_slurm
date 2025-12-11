import os
import subprocess
import pandas as pd
import numpy as np
from typing import List, Tuple
from tabulate import tabulate

from rna_library import SecStruct

from rna_map_slurm.fastq import PairedFastqFiles
from rna_map_slurm.logger import get_logger

log = get_logger("demultiplex")


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
            f"-u NC/test_R1.fastq -w NC/test_R2.fastq -m 2"
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


# int demultiplexing ##############################################


def find_helix_barcodes(
    df: pd.DataFrame, helices: List[Tuple[int, int, int]]
) -> pd.DataFrame:
    """
    Finds the sequence and bounds of helix barcodes in a dataframe of sequences and structures.

    Args:
        df: A dataframe with columns "sequence" and "structure".
        helices: A list of tuples of the form (helix_index, start_pos, end_pos).

    Returns:
        A dataframe with the same columns as the input, plus the columns barcodes,
        barcode_bounds, and full_barcode.
    """
    if len(helices) > 1:
        raise ValueError("only supporting one helix now! Will be changed soon!!")

    def __get_subsection(h1: str, h2: str, pos1: int, pos2: int) -> List[str]:
        """
        Get a subsection of two strings based on the given positions.

        Args:
            h1 (str): The first string.
            h2 (str): The second string.
            pos1 (int): The starting position of the subsection.
            pos2 (int): The ending position of the subsection.

        Returns:
            List[str]: A list containing the subsections of h1 and h2.

        """
        h1_new = h1[pos1 : pos2 + 1]
        h2_new = h2[::-1][pos1 : pos2 + 1][::-1]
        return [h1_new, h2_new]

    df["barcodes"] = [[] for _ in range(len(df))]
    df["barcode_bounds"] = [[] for _ in range(len(df))]
    df["full_barcode"] = ""
    for i, row in df.iterrows():
        s = SecStruct(row["sequence"].replace("U", "T"), row["structure"])
        row_helices = list(s.get_helices())
        all_barcodes = []
        all_bounds = []
        for j, h in enumerate(helices):
            row_h = row_helices[h[0]]
            seqs = row_h.sequence.split("&")
            strands = row_h.strands
            b_seq = __get_subsection(seqs[0], seqs[1], h[1], h[2])
            b_strands = __get_subsection(strands[0], strands[1], h[1], h[2])
            b_bounds = [
                [min(b_strands[0]), max(b_strands[0])],
                [min(b_strands[1]), max(b_strands[1])],
            ]
            all_barcodes.append(b_seq)
            all_bounds.append(b_bounds)
        full_barcode = "_".join(np.concatenate(all_barcodes).flat)
        df.at[i, "barcodes"] = all_barcodes
        df.at[i, "barcode_bounds"] = all_bounds
        df.at[i, "full_barcode"] = full_barcode
    return df
