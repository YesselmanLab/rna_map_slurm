import os
import pandas as pd
import numpy as np
from typing import List, Dict

from seq_tools.sequence import get_reverse_complement

from rna_map_slurm.logger import get_logger
from rna_map_slurm.jobs import SlurmOptions, get_job_header, generate_submit_file
from rna_map_slurm.fastq import get_paired_fastqs, PairedFastqFiles

log = get_logger(__name__)


def split_into_n(df, n):
    """
    Split a list into n sublists with sizes as close to equal as possible.
    
    Args:
    - lst (pd.dataframe): The list to be split.
    - n (int): The number of sublists.
    
    Returns:
    - list of dataframe: The split sublists.
    """
    avg = len(df) // n
    remainder = len(df) % n
    result = []
    idx = 0
    for i in range(n):
        size = avg + (1 if i < remainder else 0)
        result.append(df[idx:idx+size].copy())
        idx += size
    return result

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
        jobs.append([f"jobs/split-fastq/{name}.sh", "split-fastq", ""])
    df_jobs = pd.DataFrame(jobs, columns=["job_path", "job_type", "job_requirement"])
    generate_submit_file("submits/README_SPLIT_FASTQ", df_jobs["job_path"].tolist())
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
                "demultiplex",
                "split-fastq",
            ]
        )
    df_jobs = pd.DataFrame(jobs, columns=["job_path", "job_type", "job_requirement"])
    generate_submit_file("submits/README_DEMULTIPLEXING", df_jobs["job_path"].tolist())
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
            jobs.append([f"jobs/rna-map/{name}.sh", "rna-map", "demultiplex"])
            count += 1
    df_jobs = pd.DataFrame(jobs, columns=["job_path", "job_type", "job_requirement"])
    generate_submit_file("submits/README_RNA_MAP", df_jobs["job_path"].tolist())
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
                "rna-map-combine",
                "rna-map",
            ]
        )
    df_jobs = pd.DataFrame(jobs, columns=["job_path", "job_type", "job_requirement"])
    generate_submit_file("submits/README_RNA_MAP_COMBINE", df_jobs["job_path"].tolist())
    return df_jobs


################################################################################
############################ With internal barcodes ############################
################################################################################


def generate_internal_demultiplex_jobs(params, num_dirs):
    os.makedirs("jobs/int-demultiplex", exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())
    df = pd.read_csv("data.csv")
    slurm_params = params["slurm_options"]["internal_demultiplex"]
    runs_per_job = params["tasks_per_job"]["internal_demultiplex"]
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
            job_body = ""
            for dir in dg:
                job_body += (
                    f"cd {cur_dir}\n"
                    f"rna-map-slurm-runner int-demultiplex {dir}/{row['barcode_seq']} "
                    f"--output_dir /scratch/\n\n"
                )
            job = job_header + job_body
            f = open(f"jobs/int-demultiplex/{name}.sh", "w")
            f.write(job)
            f.close()
            jobs.append(
                [
                    f"jobs/int-demultiplex/{name}.sh",
                    "int-demultiplex",
                    "demultiplex",
                ]
            )
            i += 1
    df_jobs = pd.DataFrame(jobs, columns=["job_path", "job_type", "job_requirement"])
    generate_submit_file(
        "submits/README_INTERNAL_DEMULTIPLEXING", df_jobs["job_path"].tolist()
    )
    return df_jobs

def generate_internal_demultiplex_single_barcode(params):
    os.makedirs("jobs/int-demultiplex-sb", exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())
    df = pd.read_csv("data.csv")
    slurm_params = SlurmOptions("int-demultiplex-sb") #TODO fill in later
    runs_per_job = params["tasks_per_job"]["internal_demultiplex"] 
    add_dfs = []
    for _, row in df.iterrows():
        if pd.isnull(row["demult_cmd"]):
            continue
        df_barcodes = pd.read_json(f"inputs/barcode_jsons/{row['code']}.json")
        df_barcodes["construct_barcode"] = row["barcode_seq"]
        add_dfs.append(df_barcodes)
    df_barcodes = pd.concat(add_dfs)
    jobs = []
    dfs = split_into_n(df_barcodes, 500)
    for i, df in enumerate(dfs):
        name = f"int-demultiplex-sb-{i:04}"
        slurm_params = SlurmOptions(name)
        job_header = get_job_header(
            slurm_params, os.path.abspath("jobs/int-demultiplex-sb/")
        )
        job_body = ""
        for full_barcode, group in df.groupby("full_barcode"):
            row = group.iloc[0]
            bb1 = row["barcode_bounds"][0][0]
            bb2 = row["barcode_bounds"][0][1]
            end_len = len(row['sequence'])
            max_len = end_len - bb2[0]
            min_len = end_len - bb2[1]
            bb2 = [min_len, max_len]
            b1_seq = row['barcodes'][0][0]
            b2_seq = get_reverse_complement(row['barcodes'][0][1])
            job_body += (
                f"rna-map-slurm-runner int-demultiplex-single-barcode {row['construct_barcode']} "
                f"{b1_seq} {b2_seq} {bb1[0]} {bb1[1]} {bb2[0]} {bb2[1]}\n\n"
            )
        job = job_header + job_body
        f = open(f"jobs/int-demultiplex-sb/{name}.sh", "w")
        f.write(job)
        f.close()
        jobs.append(f"jobs/int-demultiplex-sb/{name}.sh")
    generate_submit_file("submits/README_INTERNAL_DEMULTIPLEXING_SB", jobs)


def generate_join_int_demultiplex_jobs(params):
    os.makedirs("jobs/join-int-demultiplex", exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())
    df = pd.read_csv("data.csv")
    slurm_params = params["slurm_options"]["join_internal_demultiplex"]
    if "demult_cmd" not in df.columns:
        df["demult_cmd"] = np.nan
    fsum = open("submits/README_JOIN_INT_DEMULTIPLEX", "w")
    for i, row in df.iterrows():
        name = f"join-int-demultiplex-{i:04}"
        slurm_opts = SlurmOptions(
            name,
            slurm_params["time"],
            slurm_params["mem-per-cpu"],
            slurm_params["cpus-per-task"],
            params["slurm_options"]["extra_header_cmds"],
        )
        job_header = get_job_header(
            slurm_opts, os.path.abspath("jobs/join-int-demultiplex/")
        )
        if pd.isnull(row["demult_cmd"]):
            continue
        job_body = (
            f"rna-map-slurm-runner join-int-demultiplex {cur_dir} "
            f"{row['barcode_seq']} --tmp_dir /scratch/ --threads {slurm_params['cpus-per-task']}\n"
        )
        job = job_header + job_body
        f = open(f"jobs/join-int-demultiplex/join-int-demultiplex-{i:04}.sh", "w")
        f.write(job)
        f.close()
        fsum.write(f"sbatch jobs/join-int-demultiplex/join-int-demultiplex-{i:04}.sh\n")
    fsum.close()


def generate_rna_map_single_barcode_jobs(params):
    cur_dir = os.path.abspath(os.getcwd())
    df = pd.read_csv("data.csv")
    slurm_params = params["slurm_options"]["rna_map_single_barcode"]
    runs_per_job = params["tasks_per_job"]["rna_map_single_barcode"]
    runs = []
    if "demult_cmd" not in df.columns:
        df["demult_cmd"] = np.nan
    for i, row in df.iterrows():
        if pd.isnull(row["demult_cmd"]):
            continue
        os.makedirs(f"rna_map/{row['barcode_seq']}", exist_ok=True)
        df_barcode = pd.read_json(
            f"{params['input_dir']}/barcode_jsons/{row['code']}.json"
        )
        unique_barcodes = df_barcode["full_barcode"].unique()
        for barcode in unique_barcodes:
            runs.append([row["barcode_seq"], barcode])
    run_groups = [runs[i : i + runs_per_job] for i in range(0, len(runs), runs_per_job)]
    fsum = open("submits/README_RNA_MAP_SB", "w")
    count = 0
    for i, group in enumerate(run_groups):
        name = f"rna-map-single-barcode{i:04}"
        slurm_opts = SlurmOptions(
            name,
            slurm_params["time"],
            slurm_params["mem-per-cpu"],
            slurm_params["cpus-per-task"],
            params["slurm_options"]["extra_header_cmds"],
        )
        job_header = get_job_header(
            slurm_opts, os.path.abspath("jobs/rna-map-single-barcode/")
        )
        job_body = ""
        for run in group:
            job_body += (
                f"rna-map-slurm-runner rna-map-single-barcode {cur_dir} "
                f"{run[0]} {run[1]} --tmp_dir /scratch/\n"
            )
        job = job_header + job_body
        f = open(
            f"jobs/rna-map-single-barcode/rna-map-single-barcode-{count:04}.sh", "w"
        )
        f.write(job)
        f.close()
        fsum.write(
            f"sbatch jobs/rna-map-single-barcode/rna-map-single-barcode-{count:04}.sh\n"
        )
        count += 1
    fsum.close()
