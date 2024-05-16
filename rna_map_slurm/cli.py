import click
import os
import pandas as pd
import numpy as np
import shutil

from gsheets.sheet import get_sequence_run_info_sheet, get_sequence_sheet

from barcode_demultiplex.demultiplex import find_helix_barcodes

from rna_map_slurm.fastq import get_paired_fastqs
from rna_map_slurm.logger import setup_applevel_logger, setup_logging, get_logger
from rna_map_slurm.parameters import get_parameters_from_file


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
    for index, value in df[column_name].iteritems():
        if " " in str(value):
            log.warning(
                f"Replacing spaces with underscores in row {index} for column '{column_name}'. Original value: '{value}'"
            )
            df.at[index, column_name] = value.replace(" ", "_")
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
    return seq_path


@click.group()
def cli():
    pass


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
@click.argument("data_dir")
@click.option("--param-file", type=click.Path(exists=True), default=None)
def setup(data_csv, data_dir, param_file):
    setup_logging(file_name="setup.log")
    df = pd.read_csv(data_csv)
    params = get_parameters_from_file(param_file)
    log.info(f"data.csv has {len(df)} constructs")
    log.info("Setting up run")
    log.info(f"Reading param file: {param_file}")
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
        for i in range(0, len(args)):
            if args[i] == "--helix" or args[i] == "-helix":
                helices.append([int(args[i + 1]), int(args[i + 2]), int(args[i + 3])])
        df_barcodes = find_helix_barcodes(df_seq, helices)
        df_barcodes.to_json(
            f"inputs/barcode_jsons/{row['code']}.json", orient="records"
        )


if __name__ == "__main__":
    cli()
