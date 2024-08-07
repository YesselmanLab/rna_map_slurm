# Standard library imports
import datetime
import glob
import json
import os
import re
import shutil
import sys
import time
import yaml
import zipfile
from typing import Any, Callable, Dict

# Third-party library imports
import click
import pandas as pd

# Local application imports
from barcode_demultiplex.demultiplex import find_helix_barcodes
from gsheets.sheet import get_sequence_run_info_sheet, get_sequence_sheet

# current application
from rna_map_slurm.fastq import get_paired_fastqs
from rna_map_slurm.generate_job import (
    generate_demultiplexing_jobs,
    generate_int_demultiplex_jobs,
    generate_int_demultiplex_rna_map_combine_jobs,
    generate_int_demultiplex_rna_map_jobs,
    generate_join_fastq_files_jobs,
    generate_rna_map_combine_jobs,
    generate_rna_map_jobs,
    generate_split_fastq_jobs,
    generate_trim_galore_jobs,
)
from rna_map_slurm.jobs import get_user_jobs
from rna_map_slurm.logger import get_logger, setup_logging
from rna_map_slurm.parameters import get_default_parameters, get_parameters_from_file
from rna_map_slurm.run_dask import dask_runner
from rna_map_slurm.summaries import get_demultiplexing_summary, get_pop_avg_summary

log = get_logger("cli")


def time_it(func: Callable) -> Callable:
    """
    Decorator to measure the execution time of a function.

    Args:
        func (Callable): The function to be timed.

    Returns:
        Callable: The wrapped function with timing.
    """

    def wrapper(*args, **kwargs) -> Any:
        start_time = datetime.datetime.now()
        result = func(*args, **kwargs)
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        log.info(
            f"Function '{func.__name__}' executed in {elapsed_time.total_seconds():.4f} seconds"
        )
        return result

    return wrapper


# helper functions for get_data_csv ##################################################


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
    for index, row in df.iterrows():
        value = row[column_name]
        if " " in str(value):
            log.warning(
                f"Replacing spaces with underscores in row {index} for column '{column_name}'. Original value: '{value}'"
            )
            df.at[index, column_name] = value.replace(" ", "_")
    return df


def format_sequencing_run_info(df: pd.DataFrame) -> pd.DataFrame:
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


# helper functions for setup ##########################################################


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


def setup_directories(df):
    """
    Sets up the directories for the run.
    """
    log.info("Setting up directories")
    log.info("Creating directories: jobs, submits, data, inputs, results")
    os.makedirs("jobs", exist_ok=True)
    os.makedirs("submits", exist_ok=True)
    os.makedirs("data", exist_ok=True)
    os.makedirs("inputs", exist_ok=True)
    # setup results directories ######################################
    os.makedirs("results", exist_ok=True)
    run_names = df["run_name"].unique()
    for run_name in run_names:
        os.makedirs(f"results/{run_name}", exist_ok=True)
        os.makedirs(f"results/{run_name}/processed", exist_ok=True)
        os.makedirs(f"results/{run_name}/summary", exist_ok=True)
        os.makedirs(f"results/{run_name}/plots", exist_ok=True)
        os.makedirs(f"results/{run_name}/plots/pop_avg_pngs", exist_ok=True)
        os.makedirs(f"results/{run_name}/plots/pop_avg_pngs_0_10", exist_ok=True)
        os.makedirs(f"results/{run_name}/plots/pop_avg_pngs_0_05", exist_ok=True)
    # setup input directories ######################################
    os.makedirs("inputs/barcode_jsons", exist_ok=True)
    os.makedirs("inputs/fastas", exist_ok=True)
    os.makedirs("inputs/rnas", exist_ok=True)


def setup_input_files(df, seq_path):
    """
    Sets up input files for RNA mapping.

    Args:
        df (pandas.DataFrame): The input DataFrame containing information about the RNA samples.

    Returns:
        None
    """
    # generate data dirs ###########################################
    for i, row in df.iterrows():
        df_seq = pd.read_csv(f"{seq_path}/rna/{row['code']}.csv")
        shutil.copy(f"{seq_path}/fastas/{row['code']}.fasta", "inputs/fastas/")
        df_seq.to_csv(f"inputs/rnas/{row['code']}.csv", index=False)
        helices = []
        if pd.isnull(row["demult_cmd"]):
            continue
        args = row["demult_cmd"].split()
        for i in range(0, len(args)):
            if args[i] == "--helix" or args[i] == "-helix":
                helices.append([int(args[i + 1]), int(args[i + 2]), int(args[i + 3])])
        if len(helices) == 0:
            log.error(
                f"trying to setup demultiplexing on {row['code']} but no helices are supplied!!"
            )
            continue
        df_barcodes = find_helix_barcodes(df_seq, helices)
        df_barcodes.to_json(
            f"inputs/barcode_jsons/{row['code']}.json", orient="records"
        )


def generate_jobs(df, params, all_pfqs):
    # generate all jobs
    single_df = df.query("demult_cmd.isnull()")
    single_df.to_csv("csvs/data-single.csv", index=False)
    int_mult_df = df.query("not demult_cmd.isnull()")
    df_jobs = []
    df_jobs.append(generate_split_fastq_jobs(all_pfqs, params))
    # df_jobs.append(generate_trim_galore_jobs(params))
    df_jobs.append(generate_demultiplexing_jobs(params))
    df_jobs.append(generate_join_fastq_files_jobs(params))
    df_jobs.append(generate_rna_map_jobs(params, single_df))
    df_jobs.append(generate_rna_map_combine_jobs(params, single_df))
    if len(int_mult_df) > 0:
        int_mult_df.to_csv("csvs/data-int_multiplex.csv", index=False)
        df_jobs.append(generate_int_demultiplex_jobs(params, int_mult_df))
        df_jobs.append(generate_int_demultiplex_rna_map_jobs(params, int_mult_df))
        df_jobs.append(
            generate_int_demultiplex_rna_map_combine_jobs(params, int_mult_df)
        )
    df_job = pd.concat(df_jobs)
    df_job.to_csv("jobs.csv", index=False)


# helper functions for run ############################################################


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
@time_it
@cli.command()
def run():
    start_time = time.time()
    if os.path.isfile("logs/run.log"):
        os.remove("logs/run.log")
    setup_logging(file_name="logs/run.log")
    df = pd.read_csv("jobs.csv")
    df["status"] = "not_started"
    df_can_run = df[df["job_requirement"].isna()]
    df.loc[df["job_requirement"].isna(), "status"] = "run"
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

        df_not_run = df[
            (~df["job_type"].isin(completed_types)) & (df["status"] == "not_started")
        ]
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
        if len(df_not_run) == 0:
            log.info("All jobs are submitted")
            break
    log.info("All jobs are completed")
    end_time = time.time()
    elapsed_time = end_time - start_time
    log.info(f"Elapsed time: {elapsed_time:.2f} seconds")


@time_it
@cli.command()
@click.argument("data_dirs", nargs=-1)
@click.option(
    "--num-workers", default=100, help="Number of workers to use.", required=True
)
@click.option(
    "--num-splits", default=20, help="Number of splits to use.", required=True
)
@click.option(
    "--start-step", default=0, type=int, help="Start at a specific job.", required=False
)
@click.option("--debug", is_flag=True, help="Run in debug mode.")
def run_dask(data_dirs, num_workers, num_splits, start_step, debug):
    setup_logging()
    dask_runner(data_dirs, num_workers, num_splits, start_step, debug)


@cli.command()
@click.argument("run_name")
def get_data_csv(run_name):
    os.makedirs("logs", exist_ok=True)
    if os.path.isfile("logs/get_data_csv.log"):
        os.remove("logs/get_data_csv.log")
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
    os.makedirs("logs", exist_ok=True)
    if os.path.isfile("logs/setup.log"):
        os.remove("logs/setup.log")
    setup_logging(file_name="logs/setup.log")
    df = pd.read_csv(data_csv)
    if param_file is not None:
        log.info(f"Reading param file: {param_file}")
        params = get_parameters_from_file(param_file)
    else:
        log.info("Using default parameters")
        params = get_default_parameters()
    # Log all Click parameters
    ctx = click.get_current_context()
    log.info("Ran at commandline as: %s", " ".join(sys.argv))
    log.info("Command line arguments and options:")
    for param, value in ctx.params.items():
        log.info(f"{param}: {value}")
    yaml.dump(params, open("logs/params.yaml", "w"))
    # setup data files #############################################
    seq_path = get_seq_path(params)
    os.makedirs("csvs", exist_ok=True)
    rm_df = df.query("exp_name.str.lower().str.startswith('eich')")
    rm_df.to_csv("csvs/data-eichhorn-constructs.csv", index=False)
    sub_df = df.query("not exp_name.str.lower().str.startswith('sub')")
    sub_df.to_csv("csvs/data-yesselman-constructs.csv", index=False)
    setup_directories(sub_df)
    keep = []
    for i, row in sub_df.iterrows():
        if not os.path.isfile(f"{seq_path}/rna/{row['code']}.csv"):
            log.warning(f"{row['code']} does not have a RNA CSV file!!!")
            continue
        keep.append(i)
    sub_df = sub_df.iloc[keep]
    setup_input_files(sub_df, seq_path)
    log.info("\n" + json.dumps(params, indent=4))
    log.info("saved params to logs/params.yaml")
    log.info(f"data.csv has {len(df)} constructs")
    all_pfqs = []
    for d in data_dirs:
        d = os.path.abspath(d)
        all_pfqs.extend(get_paired_fastqs(d))
    params["num_dirs"] = params["fastq_chunks"] * len(all_pfqs)
    generate_jobs(df, params, all_pfqs)


@cli.command()
def generate_summaries():
    if os.path.isfile("logs/generate-summaries.log"):
        os.remove("logs/generate-summaries.log")
    setup_logging(file_name="logs/generate-summaries.log")
    df = pd.read_csv("data.csv")
    df_barcodes = get_demultiplexing_summary()
    run_names = df["run_name"].unique()
    for run_name in run_names:
        df_sub = df.query("run_name == @run_name")
        barcodes = df_sub["barcode_seq"].unique()
        df_barcodes_sub = df_barcodes.query("sequence in @barcodes")
        df_barcodes_sub.to_csv(
            f"results/{run_name}/summary/demultiplexing.csv", index=False
        )

    df_summary = get_pop_avg_summary()
    for run_name, g in df_summary.groupby("run_name"):
        g.to_json(f"results/{run_name}/summary/summary.json", orient="records")
        g.drop(columns=["data"], inplace=True)
        g.to_csv(f"results/{run_name}/summary/summary.csv", index=False)


@cli.command()
@click.option("-v", "--version", default="v1", type=str)
@click.option("-p", "--path", default=None)
@click.option("--overwrite", type=bool)
def deposit_results(version, path, overwrite):
    setup_logging()
    if not version.startswith("v"):
        log.error("version needs to start with v")
        exit(1)
    if path is None:
        path = os.environ.get("NRDSTORSHARED")
    log_file_path = "logs/setup.log"
    data_dirs = None
    with open(log_file_path, "r") as log_file:
        for line in log_file:
            match = re.match(r"^rna-map-slurm\.cli - INFO - data_dirs: \((.*)\)", line)
            if match:
                data_dirs = match.group(1).split(",")[:1]
                break
    # remove the extra quotes
    data_dirs = [x[1:-1] for x in data_dirs]
    if len(data_dirs) != 1:
        log.error(
            "Have not implemented dealing with multiple data dirs - need to handle this manually!"
        )
    df = pd.read_csv("data.csv")
    run_names = df["run_name"].unique()
    for run_name in run_names:
        run_save_path = f"{path}/{run_name}"
        if not os.path.isdir(run_save_path):
            log.info(f"{run_save_path} does not exist, creating")
            os.makedirs(run_save_path)
        df_sub = df.query("run_name == @run_name")
        df_sub.to_csv(f"{run_save_path}/data.csv", index=False)
        log.info(f"copying data.csv to {run_save_path}")
        os.system(f"cp -r logs {run_save_path}")
        os.system(f"cp -r csvs {run_save_path}")
        log.info(f"copying logs and csvs to {run_save_path}")
        raw_path = f"{run_save_path}/raw"
        if not os.path.isdir(raw_path) and len(data_dirs) == 1:
            log.info(f"copying {data_dirs[0]} to {raw_path}")
            os.system(f"cp -r {data_dirs[0]} {raw_path}")
        demultiplex_path = f"{run_save_path}/demultiplexed"
        if not os.path.isdir(demultiplex_path):
            log.info(f"copying demultiplexed to {demultiplex_path}")
            os.system(f"cp -r demultiplexed {run_save_path}")
        analysis_path = f"{run_save_path}/analysis"
        if not os.path.isdir(analysis_path):
            log.info(f"{analysis_path} does not exist, creating")
            os.makedirs(analysis_path)
        v_path = f"{run_save_path}/analysis/{version}"
        if os.path.isdir(v_path) and not overwrite:
            log.error(
                f"{v_path} exists please pick another version number or supply --overwrite"
            )
            return
        log.info(f"copying results/{run_name} to {v_path}")
        os.system(f"cp -r results/{run_name} {v_path}")


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
        shutil.rmtree("logs", ignore_errors=True)
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
