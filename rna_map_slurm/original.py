import pandas as pd
import numpy as np
import click
import glob
import subprocess
import logging
import sys
import os
import yaml
import json
import shutil
import gzip
import random
import string
import zipfile
from pathlib import Path
from tabulate import tabulate
from dataclasses import dataclass, field
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

from seq_tools.dataframe import to_fasta, to_dna
from barcode_demultiplex.demultiplex import demultiplex as barcode_demultiplex
from barcode_demultiplex.logger import (
    setup_applevel_logger as setup_applevel_logger_barcode_demultiplex,
)

import rna_map
from rna_map.mutation_histogram import (
    merge_mut_histo_dicts,
    get_mut_histos_from_pickle_file,
    write_mut_histos_to_pickle_file,
)
from rna_map.logger import setup_applevel_logger as setup_applevel_logger_rna_map
from rna_map.run import validate_fasta_file
from rna_map.parameters import get_preset_params

from rna_map_slurm.logger import get_logger, setup_logging

from fastqsplitter import fastqsplitter


# logging ############################################################################


log = get_logger("MAIN")

# UTILS ##############################################################################


# CONVERTED
@dataclass(frozen=True, order=True)
class FastqFile:
    """
    Holds the path to a fastq file
    """

    path: str

    def __post_init__(self):
        """
        Check that the file exists
        """
        if not os.path.exists(self.path):
            raise FileNotFoundError(f"File {self.path} does not exist")
        if not self.is_r1() and not self.is_r2():
            raise ValueError(
                f"File {self.path} is not R1 or R2 please check the file name "
                f"must have either _R1_ or _R2_ in the name"
            )

    def is_compressed(self):
        """
        Check if the file is compressed
        """
        if self.path.endswith(".gz"):
            return True
        return False

    def is_r1(self):
        """
        Check if the file is R1
        """
        if "_R1" in self.path:
            return True
        return False

    def is_r2(self):
        """
        Check if the file is R2
        """
        if "_R2" in self.path:
            return True
        return False


# CONVERTED
@dataclass(frozen=True, order=True)
class PairedFastqFiles:
    """
    Holds the paths to paired fastq files
    """

    read_1: FastqFile
    read_2: FastqFile

    def is_compressed(self):
        """
        Check if the files are compressed
        """
        if self.read_1.is_compressed() and self.read_2.is_compressed():
            return True
        return False


# CONVERTED
def get_paired_fastqs(dir_path: str) -> PairedFastqFiles:
    """
    Get the paired fastq files from a directory
    :dir_path: path to directory
    :return: list of paired fastq files
    """
    if os.path.isdir(dir_path):
        f1_paths = glob.glob(os.path.join(dir_path, "*_R1*"))
        f2_paths = glob.glob(os.path.join(dir_path, "*_R2*"))
    else:
        f1_paths = glob.glob(dir_path + "*_R1*")
        f2_paths = glob.glob(dir_path + "*_R2*")
    if len(f1_paths) == 0 or len(f2_paths) == 0:
        raise FileNotFoundError(f"Could not find paired fastq files in {dir_path}")
    if len(f1_paths) > 1 or len(f2_paths) > 1:
        raise ValueError(f"Found more than one fastq file in {dir_path}")
    return PairedFastqFiles(FastqFile(f1_paths[0]), FastqFile(f2_paths[0]))


def check_if_columns_exist(df: pd.DataFrame, cols) -> None:
    """
    validates the data in the dataframe
    :param df: a dataframe with data
    :return: None
    """
    for e in cols:
        if e not in df:
            raise ValueError(
                f"{e} column is required in the dataframe! colums are "
                f"{list(df.columns)} "
            )


# CONVERTED
def gzip_files(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if not file.endswith(".gz"):  # Ignore already compressed files
                file_path = os.path.join(root, file)
                compressed_file_path = f"{file_path}.gz"
                with open(file_path, "rb") as f_in:
                    with gzip.open(compressed_file_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)

                os.remove(file_path)  # Remove the original file


# CONVERTED
def random_string(length):
    return "".join(random.choices(string.ascii_letters, k=length))


def flatten_and_zip_directory(input_directory, output_zip):
    with zipfile.ZipFile(output_zip, "w") as zip_ref:
        for root, _, files in os.walk(input_directory):
            for file in files:
                file_path = os.path.join(root, file)
                # Save the file in the zip with only its base name
                zip_ref.write(file_path, os.path.basename(file))


def combine_gzipped_fastq(input_files, output_file):
    with gzip.open(output_file, "wb") as output_gz:
        for input_file in input_files:
            with gzip.open(input_file, "rb") as input_gz:
                for line in input_gz:
                    output_gz.write(line)


def process_pair(name_pair, barcode, zip_files, outdir, tmp_dir):
    construct_barcode, pair = name_pair
    count = 0
    for zip_file in zip_files:
        try:
            with zipfile.ZipFile(zip_file, "r") as zip_ref:
                zip_ref.extract(pair[0], f"{tmp_dir}/{count}")
                zip_ref.extract(pair[1], f"{tmp_dir}/{count}")
            count += 1
        except:
            continue
    mate1_files = glob.glob(f"{tmp_dir}/*/{pair[0]}")
    mate2_files = glob.glob(f"{tmp_dir}/*/{pair[1]}")
    # if len(mate1_files) < 1:
    print(construct_barcode, len(mate1_files), len(mate2_files))
    combine_gzipped_fastq(mate1_files, f"{outdir}/{pair[0]}")
    combine_gzipped_fastq(mate2_files, f"{outdir}/{pair[1]}")
    subprocess.call(f"rm -r {tmp_dir}/*/{pair[0]}", shell=True)
    subprocess.call(f"rm -r {tmp_dir}/*/{pair[1]}", shell=True)


def get_file_size(file_path):
    file_path = os.path.realpath(file_path)
    return os.path.getsize(file_path)


# DEMULTIPLEXING ####################################################################


# CONVERTED
class SabreDemultiplexer(object):
    def setup(self, params) -> None:
        """
        setup demultiplexer
        :param params: setup parameters which are in the format of a dictionary
        :return: None
        """
        self._params = params

    def run(
        self,
        df: pd.DataFrame,
        paired_fqs: PairedFastqFiles,
        demultiplex_path,
    ):
        if not os.path.isdir(demultiplex_path):
            log.error(f"{demultiplex_path} does not exist")
            exit()
        os.chdir(demultiplex_path)
        log.info("preparing barcodes.txt file for demultiplexing")
        self.__generate_barcode_file(df)
        r1_path = paired_fqs.read_1.path
        r2_path = paired_fqs.read_2.path
        output = subprocess.check_output(
            f" sabre pe -f {r1_path} -r {r2_path} -b barcode.txt "
            f"-u NC/test_R1.fastq -w NC/test_R2.fastq -m 4",
            shell=True,
        )
        output = output.decode("UTF-8")
        log.info(f"output from sabre:\n{output}")
        for _, row in df.iterrows():
            gzip_files(row["barcode_seq"])

    def __generate_barcode_file(self, df, fname="barcode.txt"):
        expects = ["barcode", "barcode_seq", "construct"]
        check_if_columns_exist(df, expects)
        seen = []
        warning = False
        log.info(
            "constructs:\n\n"
            + tabulate(
                df[expects],
                expects,
                tablefmt="github",
                showindex=False,
            )
            + "\n"
        )
        s = ""
        for i, row in df.iterrows():
            if row["barcode"] in seen:
                log.warning(
                    f"{row['barcode']} has been used more than once this may "
                    f"be an issue"
                )
                warning = True
                continue
            s += f"{row['barcode_seq']}\t"
            s += f"{row['barcode_seq']}/test_R1.fastq\t"
            s += f"{row['barcode_seq']}/test_R2.fastq\n"
            os.makedirs(f"{row['barcode_seq']}", exist_ok=True)
            seen.append(row["barcode"])
        os.makedirs("NC", exist_ok=True)
        log.info(f"{len(seen)} unique barcodes found from csv file")
        if not warning:
            log.info("no barcode conflicts detected")
        with open(fname, "w", encoding="utf8") as f:
            f.write(s)


# SLURM JOBS ########################################################################


# CONVERTED
@dataclass(frozen=True, order=True)
class JobArguments:
    time: str = "12:00:00"
    mem_per_job: str = "2GB"
    cpus_per_task: int = 1


# CONVERTED
def get_job_header(args: JobArguments):
    return f"""#!/bin/bash
#SBATCH --time={args.time}          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu={args.mem_per_job}       # Maximum memory required per CPU (in megabytes)
#SBATCH --cpus-per-task={args.cpus_per_task}    # Number of CPUs per task
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

module load bowtie/2.3 trim_galore fastqc anaconda
conda activate
conda activate rna-map-nextflow

"""


def get_num_of_active_jobs(cmd):
    output = subprocess.check_output(cmd, shell=True)
    output = output.decode("utf8")
    lines = output.split("\n")
    count = 0
    for l in lines:
        if len(l) == 0:
            continue
        if l[0] == "[":
            continue
        count += 1
    return count - 1


# CONVERTED
def generate_submit_file(path, jobs):
    f = open(path, "w")
    for j in jobs:
        f.write(f"sbatch {j[0]}\n")
    f.close()


# JOB GENERATION ####################################################################


# CONVERTED
def generate_fastq_spliting_jobs(params):
    os.makedirs("jobs/fastq_spliting", exist_ok=True)
    dirs = params["fastq_dirs"].split(",")
    # generate jobs
    cur_dir = os.path.abspath(os.getcwd())
    job_args = JobArguments(time="04:00:00", mem_per_job="64GB", cpus_per_task=1)
    job_header = get_job_header(job_args)
    jobs = []
    for i, d in enumerate(dirs):
        try:
            fq = get_paired_fastqs(d)
        except:
            log.error(f"Could not find R1 and R2 files in {d}")
            continue
        job_body = (
            f"python {cur_dir}/run.py split-fastqs {fq.read_1.path} {fq.read_2.path} "
            f"{params['output_dir']} {params['split_fastq_chunks']} --start "
            f"{i*params['split_fastq_chunks']}\n"
        )
        job = job_header + job_body
        f = open(f"jobs/fastq_spliting/fastq-splitting-{i:04}.sh", "w")
        f.write(job)
        f.close()
        jobs.append(
            [f"jobs/fastq_spliting/fastq-splitting-{i:04}.sh", "FASTQ_SPLITTING", ""]
        )
    generate_submit_file("submits/README_FASTQ_SPLITTING", jobs)
    return jobs


# CONVERTED
def generate_demultiplexing_jobs(params):
    os.makedirs("jobs/demultiplex", exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())
    job_args = JobArguments(time="04:00:00", mem_per_job="2GB", cpus_per_task=1)
    job_header = get_job_header(job_args)
    csv_path = f"{params['input_dir']}/data.csv"
    runs_per_job = params["demultiplex_runs_per_job"]
    num_dirs = params["num_dirs"]
    dirs = [f"data/split-{i:04}" for i in range(0, num_dirs)]
    dir_groups = [dirs[i : i + runs_per_job] for i in range(0, len(dirs), runs_per_job)]
    jobs = []
    for i, dg in enumerate(dir_groups):
        job_body = ""
        for d in dg:
            r1 = f"test_R1.fastq.gz"
            r2 = f"test_R2.fastq.gz"
            job_body += (
                f"cd {cur_dir}/{d}\n"
                f"python {cur_dir}/run.py demultiplex {csv_path} {r1} {r2}\n"
            )
        job = job_header + job_body
        f = open(f"jobs/demultiplex/demultiplex-{i:04}.sh", "w")
        f.write(job)
        f.close()
        jobs.append(
            [
                f"jobs/demultiplex/demultiplex-{i:04}.sh",
                "DEMULTIPLEXING",
                "FASTQ_SPLITTING",
            ]
        )
    generate_submit_file("submits/README_DEMULTIPLEXING", jobs)
    return jobs


def generate_internal_demultiplex_jobs(params):
    os.makedirs("jobs/internal_demultiplex", exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())
    df = pd.read_csv(f"{params['input_dir']}/data.csv")
    job_args = JobArguments(time="48:00:00", mem_per_job="2GB", cpus_per_task=1)
    job_header = get_job_header(job_args)
    num_dirs = params["num_dirs"]
    runs_per_job = params["int_demultiplex_runs_per_job"]
    dirs = [f"data/split-{i:04}" for i in range(0, num_dirs)]
    jobs = []
    i = 0
    dir_groups = [dirs[i : i + runs_per_job] for i in range(0, len(dirs), runs_per_job)]
    if "demult_cmd" not in df.columns:
        df["demult_cmd"] = np.nan
    for _, row in df.iterrows():
        if pd.isnull(row["demult_cmd"]):
            continue
        for dg in dir_groups:
            job_body = ""
            for dir in dg:
                job_body += (
                    f"cd {cur_dir}\n"
                    f"python {cur_dir}/run.py int-demultiplex {cur_dir} {dir}/{row['barcode_seq']} "
                    f"--output_dir /scratch/\n\n"
                )
            job = job_header + job_body
            f = open(f"jobs/internal_demultiplex/internal-demultiplex-{i:04}.sh", "w")
            f.write(job)
            f.close()
            jobs.append(
                [
                    f"jobs/internal_demultiplex/internal-demultiplex-{i:04}.sh",
                    "INTERNAL_DEMULTIPLEXING",
                    "DEMULTIPLEXING",
                ]
            )
            i += 1
    generate_submit_file("submits/README_INTERNAL_DEMULTIPLEXING", jobs)
    return jobs


def generate_join_int_demultiplex_jobs(params):
    cur_dir = os.path.abspath(os.getcwd())
    threads = 8
    job_args = JobArguments(time="24:00:00", mem_per_job="16GB", cpus_per_task=threads)
    job_header = get_job_header(job_args)
    df = pd.read_csv(f"{params['input_dir']}/data.csv")
    if "demult_cmd" not in df.columns:
        df["demult_cmd"] = np.nan
    fsum = open("submits/README_JOIN_INT_DEMULTIPLEX", "w")
    for i, row in df.iterrows():
        if pd.isnull(row["demult_cmd"]):
            continue
        job_body = (
            f"python {cur_dir}/run.py join-int-demultiplex {cur_dir} "
            f"{row['barcode_seq']} --tmp_dir /scratch/ --threads {threads}\n"
        )
        job = job_header + job_body
        f = open(f"jobs/join_int_demultiplex/join-int-demultiplex-{i:04}.sh", "w")
        f.write(job)
        f.close()
        fsum.write(f"sbatch jobs/join_int_demultiplex/join-int-demultiplex-{i:04}.sh\n")
    fsum.close()


def generate_fastq_joining_jobs(params):
    os.makedirs("jobs/fastq_joining", exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())
    job_args = JobArguments(time="24:00:00", mem_per_job="2GB", cpus_per_task=1)
    job_header = get_job_header(job_args)
    job_body = f"python {cur_dir}/run.py fastq-concat {params['input_dir']}/data.csv\n"
    i = 0
    job = job_header + job_body
    f = open(f"jobs/fastq_joining/fastq-joining-{i:04}.sh", "w")
    f.write(job)
    f.close()
    jobs = [
        [
            f"jobs/fastq_joining/fastq-joining-{i:04}.sh",
            "FASTQ_JOINING",
            "FASTQ_SPLITTING",
        ]
    ]
    generate_submit_file("submits/README_FASTQ_JOINING", jobs)


# just normal rna-map
def generate_rna_map_jobs(params):
    """
    These are jobs that do not do the additional internal demultiplexing required
    of some large libraries
    """
    os.makedirs("jobs/rna_map", exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())
    df = pd.read_csv(f"{params['input_dir']}/data.csv")
    job_args = JobArguments(time="6:00:00", mem_per_job="2GB", cpus_per_task=1)
    job_header = get_job_header(job_args)
    runs_per_job = params["rna_map_runs_per_job"]
    num_dirs = params["num_dirs"]
    dirs = [f"data/split-{i:04}" for i in range(0, num_dirs)]
    dir_groups = [dirs[i : i + runs_per_job] for i in range(0, len(dirs), runs_per_job)]
    count = 0
    jobs = []
    if "demult_cmd" not in df.columns:
        df["demult_cmd"] = np.nan
    for _, row in df.iterrows():
        if not pd.isnull(row["demult_cmd"]):
            continue
        if not os.path.isfile(f"{params['input_dir']}/fastas/{row['code']}.fasta"):
            continue
        for dg in dir_groups:
            job_body = ""
            for dir in dg:
                os.makedirs(
                    f"{cur_dir}/{dir}/{row['barcode_seq']}/{row['construct']}",
                    exist_ok=True,
                )
                job_body += (
                    f"cd {cur_dir}/{dir}/{row['barcode_seq']}/{row['construct']}\n"
                    f"rna-map -fa {cur_dir}/inputs/fastas/{row['code']}.fasta "
                    f"-fq1 ../test_R1.fastq.gz -fq2 ../test_R2.fastq.gz "
                    f"--dot-bracket {cur_dir}/inputs/rnas/{row['code']}.csv "
                    f"--summary-output-only --param-preset barcoded-library\n"
                    f"rm -rf log input output/Mapping_Files\n"
                )
            job = job_header + job_body
            f = open(f"jobs/rna_map/rna-map-{count:04}.sh", "w")
            f.write(job)
            f.close()
            jobs.append(
                [f"jobs/rna_map/rna-map-{count:04}.sh", "RNA_MAP", "DEMULTIPLEXING"]
            )
            count += 1
    generate_submit_file("submits/README_RNA_MAP", jobs)
    return jobs


def generate_rna_map_combine_jobs(params):
    os.makedirs("jobs/rna_map_combine", exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())
    df = pd.read_csv(f"{params['input_dir']}/data.csv")
    job_args = JobArguments(time="2:00:00", mem_per_job="2GB", cpus_per_task=1)
    job_header = get_job_header(job_args)
    jobs = []
    if "demult_cmd" not in df.columns:
        df["demult_cmd"] = np.nan
    fsum = open("submits/README_RNA_MAP_COMBINE", "w")
    for i, row in df.iterrows():
        if not pd.isnull(row["demult_cmd"]):
            continue
        job_body = f"python {cur_dir}/run.py combine-rna-map {cur_dir} {row['barcode_seq']} {row['construct']}\n"
        job = job_header + job_body
        f = open(f"jobs/rna_map_combine/rna-map-combine-{i:04}.sh", "w")
        f.write(job)
        f.close()
        jobs.append(
            [f"jobs/rna_map/rna-map-combine-{i:04}.sh", "RNA_MAP_COMBINE", "RNA_MAP"]
        )
        fsum.write(f"sbatch jobs/rna_map_combine/rna-map-combine-{i:04}.sh\n")
    fsum.close()
    return jobs


def generate_rna_map_single_barcode_jobs(params):
    cur_dir = os.path.abspath(os.getcwd())
    df = pd.read_csv(f"{params['input_dir']}/data.csv")
    job_args = JobArguments(time="24:00:00", mem_per_job="2GB", cpus_per_task=1)
    job_header = get_job_header(job_args)
    runs = []
    if "demult_cmd" not in df.columns:
        df["demult_cmd"] = np.nan
    for _, row in df.iterrows():
        if pd.isnull(row["demult_cmd"]):
            continue
        os.makedirs(f"rna_map/{row['barcode_seq']}", exist_ok=True)
        df_barcode = pd.read_json(
            f"{params['input_dir']}/barcode_jsons/{row['code']}.json"
        )
        unique_barcodes = df_barcode["full_barcode"].unique()
        for barcode in unique_barcodes:
            runs.append([row["barcode_seq"], barcode])
    runs_per_job = params["rna_map_single_barcode_runs_per_job"]
    run_groups = [runs[i : i + runs_per_job] for i in range(0, len(runs), runs_per_job)]
    fsum = open("submits/README_RNA_MAP_SB", "w")
    count = 0
    for group in run_groups:
        job_body = ""
        for run in group:
            job_body += (
                f"python {cur_dir}/run.py rna-map-single-barcode {cur_dir} "
                f"{run[0]} {run[1]} --tmp_dir /scratch/\n"
            )
        job = job_header + job_body
        f = open(f"jobs/rna_map_sb/rna-map-{count:04}.sh", "w")
        f.write(job)
        f.close()
        fsum.write(f"sbatch jobs/rna_map_sb/rna-map-{count:04}.sh\n")
        count += 1
    fsum.close()


def generate_rna_map_single_barcode_combine_jobs(params):
    cur_dir = os.path.abspath(os.getcwd())
    df = pd.read_csv(f"{params['input_dir']}/data.csv")
    job_args = JobArguments(time="2:00:00", mem_per_job="2GB", cpus_per_task=1)
    job_header = get_job_header(job_args)
    fsum = open("submits/README_RNA_MAP_COMBINE", "w")
    if "demult_cmd" not in df.columns:
        df["demult_cmd"] = np.nan
    for i, row in df.iterrows():
        if pd.isnull(row["demult_cmd"]):
            continue
        job_body = (
            f"python {cur_dir}/run.py combine-rna-map-single-barcode "
            f"{cur_dir} {row['barcode_seq']}\n"
        )
        job = job_header + job_body
        f = open(f"jobs/rna_map_combine_sb/rna-map-combine-sb-{i:04}.sh", "w")
        f.write(job)
        f.close()
        fsum.write(f"sbatch jobs/rna_map_combine_sb/rna-map-combine-sb-{i:04}.sh\n")
    fsum.close()


# CLI ################################################################################


@click.group()
def cli():
    pass


@cli.command()
def test():
    setup_logging()
    log.info("testing")
    validate_fasta_file("test.fa")


# CONVERTED
@cli.command()
@click.argument("input_file", type=click.Path(exists=True))
def setup(input_file):
    setup_applevel_logger()
    # TODO add more info to print like whats in data.csv and what the params are
    # df = pd.read_csv(f"inputs/data.csv")
    log.info("creating jobs/ all jobs required to run will go here")
    log.info("creating data/ all data will be plot into split directories here")
    # make directories for jobs and data
    os.makedirs("jobs", exist_ok=True)
    os.makedirs("jobs/join_int_demultiplex", exist_ok=True)
    os.makedirs("jobs/rna_map_sb", exist_ok=True)
    os.makedirs("jobs/rna_map_combine_sb", exist_ok=True)
    os.makedirs("data", exist_ok=True)
    os.makedirs("submits", exist_ok=True)
    # log.info(f"\n params: {json.dumps(params, indent=4)}")
    params = yaml.safe_load(open(input_file))
    params["output_dir"] = os.path.abspath("data")
    num_dirs = params["split_fastq_chunks"] * len(params["fastq_dirs"].split(","))
    params["num_dirs"] = num_dirs
    for i in range(0, num_dirs):
        os.makedirs(f"data/split-{i:04}", exist_ok=True)
    # generate all jobs
    jobs = []
    jobs += generate_fastq_spliting_jobs(params)
    jobs += generate_demultiplexing_jobs(params)
    jobs += generate_rna_map_jobs(params)
    jobs += generate_rna_map_combine_jobs(params)
    jobs += generate_internal_demultiplex_jobs(params)
    generate_join_int_demultiplex_jobs(params)
    generate_rna_map_single_barcode_jobs(params)
    # generate_rna_map_single_barcode_combine_jobs(params)
    # generate_rna_map_combine_jobs(params)


@cli.command()
def job_manager():
    setup_applevel_logger()
    log.info("submitting jobs")
    jobs = []

    log.info("starting submit loop")


@cli.command()
@click.argument("input_file", type=click.Path(exists=True))
def status(input_file):
    setup_applevel_logger()
    log.info("checking status of data")
    params = yaml.safe_load(open(input_file))
    num_dirs = params["split_fastq_chunks"] * len(params["fastq_dirs"].split(","))


######################################################################################
# WORKER SCRIPTS                                                                     #
######################################################################################

# split fastqs #######################################################################

# CONVERTED


@cli.command()
@click.argument("r1_path", type=click.Path(exists=True))
@click.argument("r2_path", type=click.Path(exists=True))
@click.argument("output_dir", type=click.Path(exists=True))
@click.argument("num_chunks", type=int)
@click.option("--start", default=0)
@click.option("--debug", is_flag=True)
def split_fastqs(r1_path, r2_path, output_dir, num_chunks, start, debug):
    setup_applevel_logger()
    r1_output_files = ""
    r2_output_files = ""
    for i in range(start, num_chunks + start):
        r1_output_files += f"-o {output_dir}/split-{i:04}/test_R1.fastq.gz "
        r2_output_files += f"-o {output_dir}/split-{i:04}/test_R2.fastq.gz "
    if debug:
        log.info(f"fastqsplitter -i {r1_path} {r1_output_files}")
        log.info(f"fastqsplitter -i {r2_path} {r2_output_files}")
    subprocess.call(f"fastqsplitter -i {r1_path} {r1_output_files}", shell=True)
    subprocess.call(f"fastqsplitter -i {r2_path} {r2_output_files}", shell=True)


# demultiplex #######################################################################

# CONVERTED


@cli.command()
@click.argument("csv")
@click.argument("r1_path", type=click.Path(exists=True))
@click.argument("r2_path", type=click.Path(exists=True))
@click.option("--output_dir", type=click.Path(exists=True))
@click.option("--debug", is_flag=True)
def demultiplex(csv, r1_path, r2_path, output_dir, debug):
    """
    demultiplexes paired fastq files given 3' end barcodes
    """
    setup_applevel_logger(file_name="demultiplex.log", is_debug=debug)
    if output_dir is None:
        output_dir = os.getcwd()
    paired_fastqs = PairedFastqFiles(FastqFile(r1_path), FastqFile(r2_path))
    df = pd.read_csv(csv)
    demultiplexer = SabreDemultiplexer()
    demultiplexer.setup({})
    demultiplexer.run(df, paired_fastqs, output_dir)


# CONVERTED


# combine fastqs #######################################################################
@cli.command()
@click.argument("csv_file", type=click.Path(exists=True))
def fastq_concat(csv_file):
    os.makedirs(f"demultiplexed", exist_ok=True)
    df = pd.read_csv(csv_file)
    seq_path = os.environ["SEQPATH"]
    for barcode, g in df.groupby("barcode_seq"):
        rna_count = 0
        for row in g.iterrows():
            try:
                df_rna = pd.read_csv(f"{seq_path}/rna/{row['code']}.csv")
                if len(df_rna) < 1000:
                    rna_count += 1
            except:
                rna_count += 1
                continue
        if rna_count == 0:
            continue
        print(barcode, rna_count)
        os.makedirs(f"demultiplexed/{barcode}", exist_ok=True)
        r1_files = glob.glob(f"data/*/{barcode}/test_R1.fastq.gz")
        r2_files = glob.glob(f"data/*/{barcode}/test_R2.fastq.gz")
        os.system(
            f"cat {' '.join(r1_files)} > demultiplexed/{barcode}/test_R1.fastq.gz"
        )
        os.system(
            f"cat {' '.join(r2_files)} > demultiplexed/{barcode}/test_R2.fastq.gz"
        )


# CONVERTED


# combine rna-map results ############################################################
# takes regular rna-map results and combines them into a mutation_histos.p
@cli.command()
@click.argument("home_dir", type=click.Path(exists=True))
@click.argument("barcode_seq")
@click.argument("rna_name")
def combine_rna_map(home_dir, barcode_seq, rna_name):
    df = pd.read_csv(f"{home_dir}/inputs/data.csv")
    row = df[df["barcode_seq"] == barcode_seq].iloc[0]
    os.makedirs("results", exist_ok=True)
    run_path = "results/" + row["run_name"]
    dir_name = row["construct"] + "_" + row["code"] + "_" + row["data_type"]
    final_path = f"{run_path}/processed/{dir_name}/output/BitVector_Files/"
    os.makedirs(run_path, exist_ok=True)
    os.makedirs(run_path, exist_ok=True)
    os.makedirs(final_path, exist_ok=True)
    dirs = glob.glob("data/split-*")
    merged_mut_histos = None
    for d in dirs:
        mhs_path = (
            f"{d}/{barcode_seq}/{rna_name}/output/BitVector_Files/mutation_histos.p"
        )
        if not os.path.isfile(mhs_path):
            print(mhs_path)
            continue
        if merged_mut_histos is None:
            try:
                merged_mut_histos = get_mut_histos_from_pickle_file(mhs_path)
            except:
                print(f"could not open file: {mhs_path}")
        else:
            try:
                merge_mut_histo_dicts(
                    merged_mut_histos, get_mut_histos_from_pickle_file(mhs_path)
                )
            except:
                print(f"could not open file: {mhs_path}")
    write_mut_histos_to_pickle_file(merged_mut_histos, final_path + "mutation_histos.p")


# internal demultiplex ###############################################################
@cli.command()
@click.argument("home_dir", type=click.Path(exists=True))
@click.argument("fastq_dir", type=click.Path(exists=True))
@click.option("--output_dir", default=None)
def int_demultiplex(home_dir, fastq_dir, output_dir):
    if output_dir is None:
        output_dir = os.getcwd()
    df = pd.read_csv(f"{home_dir}/inputs/data.csv")
    fastq1_path = glob.glob(os.path.join(fastq_dir, "*R1*.fastq.gz"))[0]
    fastq2_path = glob.glob(os.path.join(fastq_dir, "*R2*.fastq.gz"))[0]
    barcode_seq = Path(fastq_dir).stem
    df_sub = df.loc[df["barcode_seq"] == barcode_seq]
    if df_sub.empty:
        raise ValueError("No barcode found in csv file")
    row = df_sub.iloc[0]
    # get helices from commandline
    helices = []
    args = row["demult_cmd"].split()
    for i in range(0, len(args)):
        if args[i] == "--helix" or args[i] == "-helix":
            helices.append([int(args[i + 1]), int(args[i + 2]), int(args[i + 3])])
    unique_code = random_string(10)
    data_path = f"{output_dir}/{unique_code}"
    df_seq = pd.read_csv(f"{home_dir}/inputs/rnas/{row['code']}.csv")
    barcode_demultiplex(
        df_seq, Path(fastq2_path), Path(fastq1_path), helices, data_path
    )
    zip_path = f"{fastq_dir}/int_demultiplexed.zip"
    flatten_and_zip_directory(data_path, zip_path)
    shutil.rmtree(data_path)


# join internal demultiplex ##########################################################
@cli.command()
@click.argument("home_dir", type=click.Path(exists=True))
@click.argument("barcode_seq")
@click.option("--threads", default=1)
@click.option("--tmp_dir", default=None)
def join_int_demultiplex(home_dir, barcode_seq, threads, tmp_dir):
    setup_applevel_logger()
    setup_applevel_logger_barcode_demultiplex()
    if tmp_dir is None:
        tmp_dir = os.getcwd()
    # get all the zip files
    dirs = glob.glob(f"{home_dir}/data/split-*")
    zip_files = []
    count = 0
    for d in dirs:
        zip_file = f"{d}/{barcode_seq}/int_demultiplexed.zip"
        if not os.path.isfile(zip_file):
            log.warning(f"{zip_file} does not exist")
            continue
        zip_files.append(zip_file)
    # get the name pairs
    df = pd.read_csv(f"{home_dir}/inputs/data.csv")
    row = df[df["barcode_seq"] == barcode_seq].iloc[0]
    df_barcode = pd.read_json(f"{home_dir}/inputs/barcode_jsons/{row['code']}.json")
    pairs = {}
    for _, r in df_barcode.iterrows():
        full_barcode = r["full_barcode"]
        pairs[full_barcode] = [
            f"{full_barcode}_mate1.fastq.gz",
            f"{full_barcode}_mate2.fastq.gz",
        ]

    outdir = f"{home_dir}/joined_fastqs/{barcode_seq}"
    os.makedirs(outdir, exist_ok=True)
    tmp_dir = tmp_dir + "/" + random_string(10)
    os.makedirs(tmp_dir, exist_ok=True)

    with ProcessPoolExecutor(max_workers=threads) as executor:
        executor.map(
            process_pair,
            pairs.items(),
            [barcode_seq] * len(pairs),
            [zip_files] * len(pairs),
            [outdir] * len(pairs),
            [tmp_dir] * len(pairs),
        )

    count = 0
    for zip_file in zip_files:
        shutil.rmtree(f"{tmp_dir}/{count}")
        count += 1


# rna-map on int demultiplexed #######################################################


@cli.command()
@click.argument("home_dir", type=click.Path(exists=True))
@click.argument("lib_barcode_seq")
@click.argument("construct_barcode_seq")
@click.option("--tmp_dir", default=None)
def rna_map_single_barcode(home_dir, lib_barcode_seq, construct_barcode_seq, tmp_dir):
    if tmp_dir is None:
        tmp_dir = os.getcwd()
    # get fastq files
    fastq_dir = f"{home_dir}/joined_fastqs/{lib_barcode_seq}"
    print(fastq_dir)
    print(f"{fastq_dir}/{construct_barcode_seq}_mate1.fastq.gz")
    mate_1_path = glob.glob(f"{fastq_dir}/{construct_barcode_seq}_mate1.fastq.gz")[0]
    mate_2_path = glob.glob(f"{fastq_dir}/{construct_barcode_seq}_mate2.fastq.gz")[0]
    # check to make sure files actually have stuff in them
    fsize_1 = get_file_size(mate_1_path)
    fsize_2 = get_file_size(mate_2_path)
    if fsize_1 < 100 or fsize_2 < 100:
        print(f"skipping {construct_barcode_seq} because file size is too small")
        return
    # get sequences to make input files
    df = pd.read_csv(f"{home_dir}/inputs/data.csv")
    row = df[df["barcode_seq"] == lib_barcode_seq].iloc[0]
    df_barcode = pd.read_json(f"{home_dir}/inputs/barcode_jsons/{row['code']}.json")
    df_barcode = df_barcode[df_barcode["full_barcode"] == construct_barcode_seq]
    tmp_dir = tmp_dir + "/" + random_string(10)
    os.makedirs(tmp_dir, exist_ok=True)
    to_fasta(to_dna(df_barcode), f"{tmp_dir}/test.fasta")
    df_barcode = df_barcode[["name", "sequence", "structure"]]
    df_barcode.to_csv(f"{tmp_dir}/test.csv", index=False)
    # update params
    params = get_preset_params("barcoded-library")
    params["dirs"]["log"] = f"{tmp_dir}/log"
    params["dirs"]["input"] = f"{tmp_dir}/input"
    params["dirs"]["output"] = f"{tmp_dir}/output"
    setup_applevel_logger_rna_map()
    # run rna-map
    rna_map.run.run(
        f"{tmp_dir}/test.fasta", mate_1_path, mate_2_path, f"{tmp_dir}/test.csv", params
    )
    # clean up unncessary files
    shutil.rmtree(f"{tmp_dir}/log")
    os.makedirs(f"{tmp_dir}/input", exist_ok=True)
    # move to unique name to avoid collisions
    shutil.move(
        f"{tmp_dir}/output/BitVector_Files/mutation_histos.p",
        f"{home_dir}/rna_map/{lib_barcode_seq}/mutation_histos_{construct_barcode_seq}.p",
    )
    shutil.rmtree(tmp_dir)


# combine rna-map single barcode results ##############################################
@cli.command()
@click.argument("home_dir", type=click.Path(exists=True))
@click.argument("barcode_seq")
def combine_rna_map_single_barcode(home_dir, barcode_seq):
    df = pd.read_csv(f"{home_dir}/inputs/data.csv")
    row = df[df["barcode_seq"] == barcode_seq].iloc[0]
    os.makedirs("results", exist_ok=True)
    run_path = "results/" + row["run_name"]
    dir_name = row["construct"] + "_" + row["code"] + "_" + row["data_type"]
    final_path = f"{run_path}/processed/{dir_name}/output/BitVector_Files/"
    os.makedirs(run_path, exist_ok=True)
    os.makedirs(run_path, exist_ok=True)
    os.makedirs(final_path, exist_ok=True)
    merged_mut_histos = None
    mhs_files = glob.glob(f"{home_dir}/rna_map/{barcode_seq}/mutation_histos_*")
    print(len(mhs_files))
    for mhs_file in mhs_files:
        if merged_mut_histos is None:
            merged_mut_histos = get_mut_histos_from_pickle_file(mhs_file)
        else:
            merge_mut_histo_dicts(
                merged_mut_histos, get_mut_histos_from_pickle_file(mhs_file)
            )
    write_mut_histos_to_pickle_file(merged_mut_histos, final_path + "mutation_histos.p")


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    cli()
