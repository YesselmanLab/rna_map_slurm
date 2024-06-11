import os
import subprocess
import pandas as pd
import shutil
from tabulate import tabulate

from fastqsplitter import split_fastqs as fastqsplitter

from rna_map_slurm.fastq import PairedFastqFiles, FastqFile
from rna_map_slurm.demultiplex import SabreDemultiplexer
from rna_map_slurm.logger import get_logger

log = get_logger("steps")


def split_fastq_file(fastq_file, output_dir, num_chunks, start, threads=1):
    # determine if fastq_file is R1 or R2
    output_file = "test_R1.fastq.gz"
    if "R2" in fastq_file:
        output_file = "test_R2.fastq.gz"
    output_files = []
    for i in range(start, num_chunks + start):
        output_files.append(f"{output_dir}/split-{i:04}/{output_file}")
    fastqsplitter(fastq_file, output_files, threads_per_file=threads)
    return output_files


def demultiplex(csv, r1_path, r2_path, output_dir):
    if output_dir is None:
        output_dir = os.getcwd()
    paired_fastqs = PairedFastqFiles(FastqFile(r1_path), FastqFile(r2_path))
    df = pd.read_csv(csv)
    demultiplexer = SabreDemultiplexer()
    demultiplexer.run(df, paired_fastqs, output_dir)
    output_dirs = []
    for df_row in df.iterrows():
        output_dirs.append(f"{output_dir}/{df_row['barcode_seq']}")
    return output_dirs


def rna_map(r1_path, r2_path, fa_path, csv_path, output_dir):
    cur_dir = os.getcwd()
    os.chdir(output_dir)
    subprocess.run(
        f"rna-map -fa {fa_path} -fq1 {r1_path} -fq2 {r2_path} "
        f"--dot-bracket {csv_path} --summary-output-only --param-preset barcoded-library",
        shell=True,
    )
    shutil.rmtree("log")
    shutil.rmtree("input")
    shutil.rmtree("output/Mapping_Files")
    os.chdir(cur_dir)
