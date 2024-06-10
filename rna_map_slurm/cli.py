import click
import os
import pandas as pd
import numpy as np
import shutil
import json
import glob
import zipfile
import time
import subprocess

from gsheets.sheet import get_sequence_run_info_sheet, get_sequence_sheet

from barcode_demultiplex.demultiplex import find_helix_barcodes

from rna_map_slurm.fastq import get_paired_fastqs
from rna_map_slurm.logger import setup_applevel_logger, setup_logging, get_logger
from rna_map_slurm.parameters import get_parameters_from_file, get_default_parameters
from rna_map_slurm.generate_job import (
    generate_split_fastq_jobs,
    generate_demultiplexing_jobs,
    generate_rna_map_jobs,
    generate_rna_map_combine_jobs,
    generate_internal_demultiplex_jobs,
    generate_join_int_demultiplex_jobs,
)
from rna_map_slurm.jobs import get_user_jobs

log = get_logger(__name__)


def replace_spaces_warn(df, column_name):
    """
    Replaces spaces with underscores in a specified column of a DataFrame and logs a warning
      for each change.

    Parameters:
    df (pd.DataFrame): The DataFrame containing the data.
    column_name (str): The name of the column in which to replace spaces with underscores.

    Returns:
    pd.DataFrame: The DataFrame with spaces replaced by underscores in the specified column.
    """
    index = 0
    for value in df[column_name]:
        if " " in str(value):
            log.warning(
                f"Replacing spaces with underscores in row {index} for column '{column_name}'. Original value: '{value}'"
            )
            df.at[index, column_name] = value.replace(" ", "_")
        index += 1
    return df


def format_sequencing_run_info(df: pd.DataFrame):
    """
    Formats the sequencing run information.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the sequencing run information.

    Returns:
        pandas.DataFrame: The formatted DataFrame containing the sequencing run information.
    """
    if len(df) == 0:
        log.error("No sequencing run information found")
        exit()
    log.info("Formatting sequencing run information")
    df = replace_spaces_warn(df, "exp_name")
    df = replace_spaces_warn(df, "construct")
    df = replace_spaces_warn(df, "exp_type")
    df_seq = get_sequence_sheet()
    demult_cmds = []
    for _, row in df.iterrows():
        seq = df_seq[df_seq["code"] == row["code"]]
        if len(seq) == 0:
            if not row["exp_name"].lower().startswith("eich"):
                log.warning(f"No sequence information found for {row['code']}")
            demult_cmds.append(None)
            continue
        seq = seq.iloc[0]
        demult_cmds.append(seq["demultiplex"])
        if not pd.isnull(seq["demultiplex"]):
            log.info(
                f"Found a demultiplexing command for {row['construct']}: {seq['demultiplex']}"
            )
    df["demult_cmd"] = demult_cmds
    return df


# TODO check if valid path?
def get_seq_path(params) -> str:
    """
    Gets the path where sequence information is stored. First check the environ
    variable SEQPATH is that is not set check the params file params["paths"]["seq_path"]

    Args:
        params (dict): The parameters dictionary.

    Returns:
        str: The path where sequence information is stored.
    """
    seq_path = ""
    try:
        seq_path = os.environ["SEQPATH"]
        log.info(f"setting seq_path from environment variable: {seq_path}")
    except:
        log.info("SEQPATH not set")
    if seq_path == "":
        seq_path = params["paths"]["seq_path"]
        log.info(f"setting seq_path from params file: {seq_path}")
    if seq_path == "":
        log.error("SEQPATH not set in environment variable or params file")
        exit()
    return seq_path


def submit_jobs(df):
    """
    Submits jobs to the SLURM scheduler.

    Args:
        df (pandas.DataFrame): The DataFrame containing the jobs to submit.
    """
    log.info("Submitting jobs")
    for _, row in df.iterrows():
        os.system(f"sbatch {row['job_path']}")


# TODO need to update this so can use sacct
# need to record which just has which id 
def is_job_type_completed(job_type, jobs):
    """
    Check if a specific job type is completed. the job_type should be a substring of the job name
    in the jobs Name dictionary

    Args:
        job_type (str): The type of job to check.
        jobs (list): A list of job objects, each is a dictionary.

    Returns:
        bool: True if all jobs are completed i.e. not running anymore.
    """
    for job in jobs:
        if job_type in job["Name"]:
            return False
    return True
              
@click.group()
def cli():
    pass

# TODO need a way to check to see if existing jobs have been run
@cli.command()
def run():
    setup_logging(file_name="run.log")
    #user_jobs = get_user_jobs("jyesselm")
    df = pd.read_csv("jobs.csv")
    df["status"] = "not_started"
    df_can_run = df[df['job_requirement'].isna()]
    df.loc[ df['job_requirement'].isna(), "status"] = "run"
    completed_types = [] 
    submitted_types = df_can_run["job_type"].to_list()
    submit_jobs(df_can_run)
    while True:        
        # Wait for the submitted jobs to finish
        time.sleep(60)
        jobs = get_user_jobs("jyesselm")
        log.info("num_jobs: " + str(len(jobs)))
        log.info("submitted_types: " + str(submitted_types))
        log.info("completed_types: " + str(completed_types))
        for job_type in submitted_types:
            if is_job_type_completed(job_type, jobs):
                log.info(f"Job type {job_type} is completed")
                completed_types.append(job_type)
        for job_type in completed_types:
            if job_type in submitted_types:
                submitted_types.remove(job_type)
        
        df_not_run = df[(~df["job_type"].isin(completed_types)) & (df["status"] == "not_started")]
        submitted = False
        for job_type, g in df_not_run.groupby("job_type"):
            requirement = g["job_requirement"].iloc[0]
            if requirement not in completed_types:
                continue
            job_num = len(jobs)
            count = 0
            log.info(f"Submitting jobs for {job_type}")
            if job_type not in submitted_types:
                submitted_types.append(job_type)
            for _, row in g.iterrows():
                if job_num > 999:
                    break
                os.system(f"sbatch {row['job_path']}")
                job_num += 1
                df.loc[row.name, "status"] = "run"
                submitted = True
                count += 1
            log.info(f"Submitted {count} jobs for {job_type}")
            if submitted:
                break

            
@cli.command()
def setup_celery():
    output = subprocess.run("hostname -I", shell=True, capture_output=True, text=True)
    print(output.stdout)


@cli.command()
@click.argument("run_name")
def get_data_csv(run_name):
    setup_logging(file_name="get_data_csv.log")
    df = get_sequence_run_info_sheet()
    # setup data.csv file ##########################################
    df = get_sequence_run_info_sheet()
    df = df[df["run_name"] == run_name]
    df = format_sequencing_run_info(df)
    df.to_csv("data.csv", index=False)


@cli.command()
@click.argument("data_csv")
@click.argument("data_dirs", nargs=-1)
@click.option("--param-file", type=click.Path(exists=True), default=None)
def setup(data_csv, data_dirs, param_file):
    setup_logging(file_name="setup.log")
    df = pd.read_csv(data_csv)
    if param_file is not None:
        log.info(f"Reading param file: {param_file}")
        params = get_parameters_from_file(param_file)
    else:
        log.info("Using default parameters")
        params = get_default_parameters()
    log.info("\n" + json.dumps(params, indent=4))
    log.info(f"data.csv has {len(df)} constructs")
    log.info("Setting up run")
    # setup directories ###########################################
    log.info("Setting up directories")
    log.info("Creating directories: jobs, submits, data, inputs")
    os.makedirs("jobs", exist_ok=True)
    os.makedirs("submits", exist_ok=True)
    os.makedirs("data", exist_ok=True)
    os.makedirs("inputs", exist_ok=True)
    # setup input directories ######################################
    os.makedirs("inputs/barcode_jsons", exist_ok=True)
    os.makedirs("inputs/fastas", exist_ok=True)
    os.makedirs("inputs/rnas", exist_ok=True)
    # setup input files ############################################
    seq_path = get_seq_path(params)
    if "demult_cmd" not in df.columns:
        df["demult_cmd"] = np.nan
    # generate data dirs ###########################################
    num_int_demult = 0
    for i, row in df.iterrows():
        if row["exp_name"].lower().startswith("eich"):
            continue
        if not os.path.isfile(f"{seq_path}/rna/{row['code']}.csv"):
            log.warning(f"{row['code']} does not have a RNA CSV file!!!")
            continue
        df_seq = pd.read_csv(f"{seq_path}/rna/{row['code']}.csv")
        shutil.copy(f"{seq_path}/fastas/{row['code']}.fasta", "inputs/fastas/")
        df_seq.to_csv(f"inputs/rnas/{row['code']}.csv", index=False)
        helices = []
        if pd.isnull(row["demult_cmd"]):
            continue
        args = row["demult_cmd"].split()
        if len(args) > 1:
            num_int_demult += 1
        for i in range(0, len(args)):
            if args[i] == "--helix" or args[i] == "-helix":
                helices.append([int(args[i + 1]), int(args[i + 2]), int(args[i + 3])])
        df_barcodes = find_helix_barcodes(df_seq, helices)
        df_barcodes.to_json(
            f"inputs/barcode_jsons/{row['code']}.json", orient="records"
        )
    all_pfqs = []
    for d in data_dirs:
        d = os.path.abspath(d)
        all_pfqs.extend(get_paired_fastqs(d))
    num_dirs = params["fastq_chunks"] * len(all_pfqs)
    params["num_dirs"] = num_dirs
    for i in range(0, num_dirs):
        os.makedirs(f"data/split-{i:04}", exist_ok=True)
    # generate all jobs
    df_jobs = []
    df_jobs.append(generate_split_fastq_jobs(all_pfqs, params))
    df_jobs.append(generate_demultiplexing_jobs(params, num_dirs))
    df_jobs.append(generate_rna_map_jobs(params, num_dirs))
    df_jobs.append(generate_rna_map_combine_jobs(params))
    if num_int_demult > 0:
        df_jobs.append(generate_internal_demultiplex_jobs(params, num_dirs))
        df_jobs.append(generate_join_int_demultiplex_jobs(params))
    df_job = pd.concat(df_jobs)
    df_job.to_csv("jobs.csv", index=False)


@cli.command()
@click.argument("stage")
def clean(stage):
    setup_logging()
    if stage == "all":
        log.info("Cleaning all directories")
        shutil.rmtree("jobs", ignore_errors=True)
        shutil.rmtree("submits", ignore_errors=True)
        shutil.rmtree("data", ignore_errors=True)
        shutil.rmtree("inputs", ignore_errors=True)
    elif stage == "demultiplex":
        log.info("Cleaning demultiplex directories")
        dirs = glob.glob("data/split-*/[ACGT]*")
    elif stage == "rna_map":
        log.info("Cleaning rna_map directories")
        dirs = glob.glob("data/split-*/[ACGT]*/*")
        for d in dirs:
            if not os.path.isdir(d):
                continue
            os.system(f"rm -rf {d}/*")


@cli.command()
@click.argument("csv")
def zip_demultiplex_subset(csv):
    df = pd.read_csv(csv)
    base_path = "demultiplexed/"
    # Get the list of barcodes from the DataFrame
    barcodes = df["barcode_seq"].unique()

    # Create a zip file for all barcodes
    zip_filename = "demultiplex_subset.zip"
    with zipfile.ZipFile(zip_filename, "w", zipfile.ZIP_DEFLATED) as zipf:
        for barcode in barcodes:
            # Create the path to the barcode directory
            barcode_path = os.path.join(base_path, barcode)
            if not os.path.exists(barcode_path):
                print(f"Barcode path {barcode_path} does not exist, skipping.")
                continue

            # Find all .gz files in the barcode directory
            for root, _, files in os.walk(barcode_path):
                for file in files:
                    if file.endswith(".gz"):
                        file_path = os.path.join(root, file)
                        # Add the file to the zip file, preserving the directory structure
                        zipf.write(file_path, os.path.relpath(file_path, base_path))
                        print(f"Added {file_path} to {zip_filename}")

    print(f"Created zip file {zip_filename}")


if __name__ == "__main__":
    cli()
