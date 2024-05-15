import click
import pandas as pd
from fastqsplitter import split_fastqs as fastqsplitter

from gsheets.sheet import get_sequence_run_info_sheet, get_sequence_sheet

from rna_map_slurm.fastq import get_paired_fastqs
from rna_map_slurm.logger import setup_applevel_logger, setup_logging, get_logger


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
        None
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
            continue
        seq = seq.iloc[0]
        demult_cmds.append(seq["demultiplex"])
        if not pd.isnull(seq["demultiplex"]):
            log.info(
                f"Found a demultiplexing command for {row['construct']}: {seq['demultiplex']}"
            )


@click.group()
def cli():
    pass


# TODO add schema valiadation of params!
# TODO compute a way to figure out the optimial number of chunks to split into
# TODO hide job output when its working!
@cli.command()
@click.argument("run_name")
@click.argument("input_file", type=click.Path(exists=True))
def setup(run_name, input_file):
    setup_logging(file_name="rna_map_slurm.log")
    log.info(f"Setting up run: {run_name}")
    log.info(f"Reading input file: {input_file}")
    # df = get_sequence_run_info_sheet()
    # df = df[df["run_name"] == run_name]
    # df.to_csv("data.csv", index=False)
    df = pd.read_csv("data.csv")
    df = format_sequencing_run_info(df)


@cli.command()
@click.argument("r1_path", type=click.Path(exists=True))
@click.argument("r2_path", type=click.Path(exists=True))
@click.argument("output_dir", type=click.Path(exists=True))
@click.argument("num_chunks", type=int)
@click.option("--start", default=0)
@click.option("--threads", default=1)
def split_fastqs(r1_path, r2_path, output_dir, num_chunks, start, threads):
    r1_output_files = []
    r2_output_files = []
    for i in range(start, num_chunks + start):
        r1_output_files.append(f"{output_dir}/split-{i:04}/test_R1.fastq.gz")
        r2_output_files.append(f"{output_dir}/split-{i:04}/test_R2.fastq.gz")
    fastqsplitter(r1_path, r1_output_files, threads_per_file=threads)
    fastqsplitter(r2_path, r2_output_files, threads_per_file=threads)


if __name__ == "__main__":
    cli()
