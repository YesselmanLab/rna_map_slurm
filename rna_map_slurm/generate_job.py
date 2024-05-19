import os
import pandas as pd
import numpy as np
from typing import List, Dict

from rna_map_slurm.logger import get_logger
from rna_map_slurm.jobs import SlurmOptions, get_job_header, generate_submit_file
from rna_map_slurm.fastq import get_paired_fastqs, PairedFastqFiles

log = get_logger(__name__)


################################################################################
############################ Use everytime function ############################
################################################################################


# TODO split into two jobs for each pair of fastq files
def generate_split_fastq_jobs(
    pfqs: List[PairedFastqFiles], params: Dict
) -> pd.DataFrame:
    os.makedirs("jobs/split-fastq", exist_ok=True)
    # generate jobs
    cur_dir = os.path.abspath(os.getcwd())
    slurm_params = params["slurm_options"]["split_fastq"]
    jobs = []
    for i, pfq in enumerate(pfqs):
        name = f"split-fastq-{i:04}"
        slurm_opts = SlurmOptions(
            name,
            slurm_params["time"],
            slurm_params["mem-per-cpu"],
            slurm_params["cpus-per-task"],
            params["slurm_options"]["extra_header_cmds"],
        )
        job_header = get_job_header(slurm_opts, os.path.abspath("jobs/split-fastq/"))
        job_body = (
            f"rna-map-slurm-runner split-fastqs {pfq.read_1.path} {pfq.read_2.path} "
            f"{cur_dir + '/data'} {params['fastq_chunks']} --start "
            f"{i*params['fastq_chunks']}\n"
        )
        job = job_header + job_body
        f = open(f"jobs/split-fastq/{name}.sh", "w")
        f.write(job)
        f.close()
        jobs.append([f"jobs/split-fastq/{name}.sh", "SPLIT_FASTQ", ""])
    df_jobs = pd.DataFrame(jobs, columns=["job", "type", "status"])
    generate_submit_file("submits/README_SPLIT_FASTQ", df_jobs["job"].tolist())
    return df_jobs


def generate_demultiplexing_jobs(params, num_dirs):
    os.makedirs("jobs/demultiplex", exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())
    csv_path = os.path.abspath("data.csv")
    runs_per_job = params["tasks_per_job"]["demultiplex"]
    dirs = [f"data/split-{i:04}" for i in range(0, num_dirs)]
    slurm_params = params["slurm_options"]["demultiplex"]
    dir_groups = [dirs[i : i + runs_per_job] for i in range(0, len(dirs), runs_per_job)]
    jobs = []
    for i, dg in enumerate(dir_groups):
        name = f"demultiplex-{i:04}"
        slurm_opts = SlurmOptions(
            name,
            slurm_params["time"],
            slurm_params["mem-per-cpu"],
            slurm_params["cpus-per-task"],
            params["slurm_options"]["extra_header_cmds"],
        )
        job_header = get_job_header(slurm_opts, os.path.abspath("jobs/demultiplex/"))
        job_body = ""
        for d in dg:
            job_body += (
                f"cd {cur_dir}/{d}\n"
                f"rna-map-slurm-runner demultiplex {csv_path} test_R1.fastq.gz test_R2.fastq.gz\n"
            )
        job = job_header + job_body
        f = open(f"jobs/demultiplex/{name}.sh", "w")
        f.write(job)
        f.close()
        jobs.append(
            [
                f"jobs/demultiplex/demultiplex-{i:04}.sh",
                "DEMULTIPLEXING",
                "FASTQ_SPLITTING",
            ]
        )
    df_jobs = pd.DataFrame(jobs, columns=["job", "type", "status"])
    generate_submit_file("submits/README_DEMULTIPLEXING", df_jobs["job"].tolist())
    return df_jobs


def generate_rna_map_jobs(params, num_dirs):
    """
    These are jobs that do not do the additional internal demultiplexing required
    of some large libraries
    """
    os.makedirs("jobs/rna-map", exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())
    csv_path = os.path.abspath("data.csv")
    df = pd.read_csv(csv_path)
    runs_per_job = params["tasks_per_job"]["rna_map"]
    slurm_params = params["slurm_options"]["rna_map"]
    dirs = [f"data/split-{i:04}" for i in range(0, num_dirs)]
    dir_groups = [dirs[i : i + runs_per_job] for i in range(0, len(dirs), runs_per_job)]
    count = 0
    jobs = []
    # should not be necessary anymore but worth leaving
    if "demult_cmd" not in df.columns:
        df["demult_cmd"] = np.nan
    for _, row in df.iterrows():
        if not pd.isnull(row["demult_cmd"]):
            continue
        if row["exp_name"].lower().startswith("eich"):
            continue
        if not os.path.isfile(f"inputs/fastas/{row['code']}.fasta"):
            log.warning(f"{row['code']} does not have a RNA CSV file!!!")
            continue
        for dg in dir_groups:
            name = f"rna-map-{count:04}"
            slurm_opts = SlurmOptions(
                name,
                slurm_params["time"],
                slurm_params["mem-per-cpu"],
                slurm_params["cpus-per-task"],
                params["slurm_options"]["extra_header_cmds"],
            )
            job_header = get_job_header(slurm_opts, os.path.abspath("jobs/rna-map/"))
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
            f = open(f"jobs/rna-map/{name}.sh", "w")
            f.write(job)
            f.close()
            jobs.append([f"jobs/rna-map/{name}.sh", "RNA_MAP", "DEMULTIPLEXING"])
            count += 1
    df_jobs = pd.DataFrame(jobs, columns=["job", "type", "status"])
    generate_submit_file("submits/README_RNA_MAP", df_jobs["job"].tolist())
    return df_jobs


def generate_rna_map_combine_jobs(params):
    os.makedirs("jobs/rna-map-combine", exist_ok=True)
    df = pd.read_csv("data.csv")
    slurm_params = params["slurm_options"]["rna_map_combine"]
    jobs = []
    if "demult_cmd" not in df.columns:
        df["demult_cmd"] = np.nan
    for i, row in df.iterrows():
        name = f"rna-map-combine-{i:04}"
        slurm_opts = SlurmOptions(
            name,
            slurm_params["time"],
            slurm_params["mem-per-cpu"],
            slurm_params["cpus-per-task"],
            params["slurm_options"]["extra_header_cmds"],
        )
        job_header = get_job_header(
            slurm_opts, os.path.abspath("jobs/rna-map-combine/")
        )
        if not pd.isnull(row["demult_cmd"]):
            continue
        if row["exp_name"].lower().startswith("eich"):
            continue
        job_body = f"rna-map-slurm-runner combine-rna-map {row['barcode_seq']} {row['construct']}\n"
        job = job_header + job_body
        f = open(f"jobs/rna-map-combine/{name}.sh", "w")
        f.write(job)
        f.close()
        jobs.append(
            [
                f"jobs/rna-map-combine/{name}.sh",
                "RNA_MAP_COMBINE",
                "RNA_MAP",
            ]
        )
    df_jobs = pd.DataFrame(jobs, columns=["job", "type", "status"])
    generate_submit_file("submits/README_RNA_MAP_COMBINE", df_jobs["job"].tolist())
    return df_jobs


################################################################################
############################ With internal barcodes ############################
################################################################################


def generate_internal_demultiplex_jobs(params, num_dirs):
    os.makedirs("jobs/internal-demultiplex", exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())
    df = pd.read_csv("data.csv")
    slurm_params = params["slurm_options"]["internal_demultiplex"]
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
        name = f"int-demultiplex-{i:04}"
        slurm_opts = SlurmOptions(
            name,
            slurm_params["time"],
            slurm_params["mem-per-cpu"],
            slurm_params["cpus-per-task"],
            params["slurm_options"]["extra_header_cmds"],
        )
        job_header = get_job_header(
            slurm_opts, os.path.abspath("jobs/int-demultiplex/")
        )
        for dg in dir_groups:
            job_body = ""
            for dir in dg:
                job_body += (
                    f"cd {cur_dir}\n"
                    f"rna-map-slurm-runner int-demultiplex {dir}/{row['barcode_seq']} "
                    f"--output_dir /scratch/\n\n"
                )
            job = job_header + job_body
            f = open(f"jobs/internal_demultiplex/{name}.sh", "w")
            f.write(job)
            f.close()
            jobs.append(
                [
                    f"jobs/internal_demultiplex/{name}.sh",
                    "INTERNAL_DEMULTIPLEXING",
                    "DEMULTIPLEXING",
                ]
            )
            i += 1
    generate_submit_file("submits/README_INTERNAL_DEMULTIPLEXING", jobs)
    return jobs
