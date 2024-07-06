import os
import pandas as pd
import shutil
import subprocess
import glob
from typing import List, Optional

from fastqsplitter import split_fastqs as fastqsplitter

from seq_tools.sequence import get_reverse_complement
from rna_map.parameters import get_preset_params, parse_parameters_from_file
import rna_map
from rna_map.mutation_histogram import (
    get_dataframe,
    get_mut_histos_from_pickle_file,
    merge_mut_histo_dicts,
    write_mut_histos_to_pickle_file,
)

from rna_map_slurm.demultiplex import SabreDemultiplexer
from rna_map_slurm.fastq import PairedFastqFiles, FastqFile
from rna_map_slurm.logger import get_logger
from rna_map_slurm.rna_map_funcs import generate_pop_avg_plots, get_mut_histo_dataframe
from rna_map_slurm.util import get_data_row, get_file_size, random_string

log = get_logger("TASKS")


class BasicTasks:
    @staticmethod
    def split_fastq_file(
        fastq_file: str, output_dir: str, num_chunks: int, start: int, threads: int = 1
    ) -> List[str]:
        """
        Splits a FASTQ file into multiple chunks and saves the chunks to the specified output directory.

        Args:
            fastq_file (str): Path to the input FASTQ file.
            output_dir (str): Directory where the split files will be saved.
            num_chunks (int): Number of chunks to split the FASTQ file into.
            start (int): Starting index for the output chunk files.
            threads (int, optional): Number of threads to use for splitting. Defaults to 1.

        Returns:
            List[str]: List of paths to the split output files.
        """
        log.info(f"Splitting {fastq_file} into {num_chunks} chunks")
        log.info(
            f"output_dir: {output_dir}, num_chunks: {num_chunks}, start: {start}, threads: {threads}"
        )
        output_file = "test_R1.fastq.gz"
        if "R2" in fastq_file:
            log.info("Detected R2 file")
            output_file = "test_R2.fastq.gz"
        else:
            log.info("Detected R1 file")
        output_files = []
        for i in range(start, num_chunks + start):
            output_files.append(f"{output_dir}/split-{i:04}/{output_file}")
        fastqsplitter(fastq_file, output_files, threads_per_file=threads)
        return output_files

    @staticmethod
    def trim_galore(
        r1_path: str, r2_path: str, output_dir: Optional[str] = None
    ) -> None:
        if output_dir is None:
            output_dir = os.getcwd()
        cur_dir = os.getcwd()
        subprocess.call(
            f"trim_galore --quality 0 --paired {r1_path} {r2_path}", shell=True
        )
        os.chdir(cur_dir)

    @staticmethod
    def demultiplex(
        csv: str, r1_path: str, r2_path: str, output_dir: Optional[str] = None
    ) -> List[str]:
        """
        Executes the demultiplexing task using the provided CSV and FASTQ files, and returns the list of output directories.

        Args:
            csv (str): Path to the CSV file containing barcode sequences and other metadata.
            r1_path (str): Path to the R1 FASTQ file.
            r2_path (str): Path to the R2 FASTQ file.
            output_dir (Optional[str]): Directory where the output files will be saved. If None, uses the current working directory.

        Returns:
            List[str]: List of directories where the demultiplexed files are saved.
        """
        if output_dir is None:
            output_dir = os.getcwd()
        cur_dir = os.getcwd()
        os.chdir(output_dir)
        paired_fastqs = PairedFastqFiles(FastqFile(r1_path), FastqFile(r2_path))
        df = pd.read_csv(csv)
        demultiplexer = SabreDemultiplexer()
        demultiplexer.run(df, paired_fastqs, output_dir)
        os.chdir(cur_dir)
        output_dirs = [
            os.path.join(output_dir, df_row["barcode_seq"])
            for i, df_row in df.iterrows()
        ]
        return output_dirs

    @staticmethod
    def join_fastq_files(fastq_files: List[str], joined_fastq: str) -> str:
        """
        Concatenates multiple FASTQ files into a single file.

        Args:
            fastq_files (List[str]): List of paths to FASTQ files to concatenate.
            joined_fastq (str): Path to the output joined FASTQ file.

        Returns:
            str: Path to the joined FASTQ file.
        """
        log.info(f"Joining FASTQ files into {joined_fastq}")
        os.system(f"cat {' '.join(fastq_files)} > {joined_fastq}")
        log.info(f"FASTQ files joined: {joined_fastq}")
        return joined_fastq

    @staticmethod
    def rna_map(
        fa_path: str, r1_path: str, r2_path: str, csv_path: str, output_dir: str
    ) -> None:
        """
        Runs the RNA mapping task with the provided FASTA, FASTQ, and CSV files.

        Args:
            fa_path (str): Path to the FASTA file.
            r1_path (str): Path to the R1 FASTQ file.
            r2_path (str): Path to the R2 FASTQ file.
            csv_path (str): Path to the CSV file with mapping parameters.
            output_dir (str): Directory to store the output files.

        Returns:
            None
        """
        cur_dir = os.getcwd()
        log.info(f"Changing directory to {output_dir}")
        os.chdir(output_dir)
        # TODO find a way to resolve this hack
        # a hack now to deal with adding additional constructs to be processed
        if os.path.isfile("input.fasta"):
            log.info("Using existing input.fa file")
            fa_path = "input.fasta"
        if os.path.isfile("input.csv"):
            log.info("Using existing input.csv file")
            csv_path = "input.csv"
        if os.path.isfile("params.yml"):
            log.info("Using existing parameter fileparams.yml file")
            params = parse_parameters_from_file("params.yml")
        else:
            params = get_preset_params("barcoded-library")
        params["overwrite"] = True
        # TODO do I want to do this?
        params["bit_vector"]["summary_output_only"] = True
        log.info("Starting RNA mapping")
        rna_map.run.run(fa_path, r1_path, r2_path, csv_path, params)

        log.info("Cleaning up unnecessary files: log, input, output/Mapping_Files")
        shutil.rmtree("log")
        shutil.rmtree("input")
        shutil.rmtree("output/Mapping_Files")
        os.chdir(cur_dir)

    @staticmethod
    def rna_map_combine(row) -> None:
        barcode_seq = row["barcode_seq"]
        construct = row["construct"]
        run_path = "results/" + row["run_name"]
        dir_name = row["construct"] + "_" + row["code"] + "_" + row["data_type"]
        final_path = f"{run_path}/processed/{dir_name}/output/BitVector_Files/"
        log.info(f"results path: {final_path}")
        os.makedirs(run_path, exist_ok=True)
        os.makedirs(final_path, exist_ok=True)
        dirs = glob.glob("data/split-*")
        merged_mut_histos = None
        count_files = 0
        merged_mut_histos = {}
        for d in dirs:
            mhs_path = f"{d}/{barcode_seq}/{construct}/output/BitVector_Files/mutation_histos.p"
            if not os.path.isfile(mhs_path):
                log.warning("files not found:" + mhs_path)
                continue
            merge_mut_histo_dicts(
                merged_mut_histos, get_mut_histos_from_pickle_file(mhs_path)
            )
            count_files += 1
        log.info(f"merged {count_files} files")
        df_results = get_mut_histo_dataframe(merged_mut_histos)
        cols = list(row.keys())
        for c in "demult_cmd,length".split(","):
            cols.remove(c)
        for c in cols:
            df_results[c] = row[c]
        df_results.to_json(final_path + "mutation_histos.json", orient="records")
        generate_pop_avg_plots(df_results, row["run_name"], dir_name)
        write_mut_histos_to_pickle_file(
            merged_mut_histos, final_path + "mutation_histos.p"
        )


class IntDemultiplexTasks:
    """ """

    @staticmethod
    def int_demultiplex(
        construct_barcode,
        b1_seq,
        b2_seq,
        b1_min_pos,
        b1_max_pos,
        b2_min_pos,
        b2_max_pos,
    ) -> None:
        tmp_dir = "/scratch/" + random_string(10)
        log.info(f"tmp_dir: {tmp_dir}")
        os.makedirs(tmp_dir, exist_ok=True)
        r2_path = f"demultiplexed/{construct_barcode}/test_R2.fastq.gz"
        r1_path = f"demultiplexed/{construct_barcode}/test_R1.fastq.gz"
        b2_seq_rc = get_reverse_complement(b2_seq)
        os.system(
            f'seqkit grep -s -p "{b1_seq}" -P -R {b1_min_pos-2}:{b1_max_pos+2} {r2_path} -o {tmp_dir}/test_R2.fastq.gz'
        )
        os.system(
            f'seqkit grep -s -p "{b2_seq_rc}" -P -R {b2_min_pos-2}:{b2_max_pos+2} {r1_path} -o {tmp_dir}/test_R1.fastq.gz'
        )
        os.system(f"seqkit seq -n {tmp_dir}/test_R2.fastq.gz > {tmp_dir}/R2_names.txt")
        os.system(f"seqkit seq -n {tmp_dir}/test_R1.fastq.gz > {tmp_dir}/R1_names.txt")
        os.system(f"sort {tmp_dir}/R1_names.txt > {tmp_dir}/R1_names_sorted.txt")
        os.system(f"sort {tmp_dir}/R2_names.txt > {tmp_dir}/R2_names_sorted.txt")
        cmd = (
            "awk 'NR==FNR{a[$1]; next} $1 in a' "
            + f"{tmp_dir}/R1_names_sorted.txt {tmp_dir}/R2_names_sorted.txt "
            + "| awk '{print($1)}' >"
            + f"{tmp_dir}/common_names.txt"
        )
        os.system(cmd)
        os.system(
            f"seqkit grep -f {tmp_dir}/common_names.txt {tmp_dir}/test_R2.fastq.gz -o int-demultiplexed/{construct_barcode}/{b1_seq}_{b2_seq}_mate1.fastq.gz"
        )
        os.system(
            f"seqkit grep -f {tmp_dir}/common_names.txt {tmp_dir}/test_R1.fastq.gz -o int-demultiplexed/{construct_barcode}/{b1_seq}_{b2_seq}_mate2.fastq.gz"
        )

    @staticmethod
    def int_demultiplex_rna_map():
        pass
