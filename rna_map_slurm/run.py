# Standard library imports
import datetime
import glob
import os
import shutil
from typing import Any, Callable

# Third-party library imports
import click
import matplotlib.pyplot as plt
import pandas as pd

# Local application imports
from barcode_demultiplex.demultiplex import demultiplex as barcode_demultiplex
from fastqsplitter import split_fastqs as fastqsplitter
from seq_tools.dataframe import to_dna, to_fasta
from seq_tools.sequence import get_reverse_complement
import rna_map
from rna_map.mutation_histogram import (
    get_dataframe,
    get_mut_histos_from_pickle_file,
    merge_mut_histo_dicts,
    write_mut_histos_to_pickle_file,
)
from rna_map.parameters import get_preset_params

# Current application imports
from rna_map_slurm.logger import get_logger, setup_logging
from rna_map_slurm.plotting import plot_pop_avg_from_row
from rna_map_slurm.tasks import BasicTasks, IntDemultiplexTasks
from rna_map_slurm.util import get_data_row, get_file_size, random_string

log = get_logger("CLI")


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


# helper functions #################################################################


# cli ##############################################################################


@click.group()
def cli():
    pass


################################################################################
############################ Use everytime function ############################
################################################################################


@time_it
@cli.command()
@click.argument("r1_path", type=click.Path(exists=True), required=True)
@click.argument("r2_path", type=click.Path(exists=True), required=True)
@click.argument("output_dir", type=click.Path(exists=True), required=True)
@click.argument("num_chunks", type=int, required=True)
@click.option(
    "--start",
    default=0,
    show_default=True,
    help="Starting index for the chunk numbering.",
)
@click.option(
    "--threads",
    default=1,
    show_default=True,
    help="Number of threads to use for splitting the FASTQ files.",
)
def split_fastqs(r1_path, r2_path, output_dir, num_chunks, start, threads):
    """
    Split the input FASTQ files into multiple chunks.

    Arguments:
        r1_path (str): Path to the input R1 FASTQ file.
        r2_path (str): Path to the input R2 FASTQ file.
        output_dir (str): Path to the output directory where the split FASTQ files will be saved.
        num_chunks (int): Number of chunks to split the FASTQ files into.
    """
    setup_logging()
    BasicTasks.split_fastq_file(r1_path, output_dir, num_chunks, start, threads)
    BasicTasks.split_fastq_file(r2_path, output_dir, num_chunks, start, threads)


@time_it
@cli.command()
@click.argument("csv")
@click.argument("r1_path", type=click.Path(exists=True))
@click.argument("r2_path", type=click.Path(exists=True))
@click.argument("output_dir", type=click.Path(exists=True))
def demultiplex(csv, r1_path, r2_path, output_dir):
    """
    demultiplexes paired fastq files given 3' end barcodes
    """
    setup_logging()
    BasicTasks.demultiplex(csv, r1_path, r2_path, output_dir)


@time_it
@cli.command()
@click.argument("fasta_path", type=click.Path(exists=True))
@click.argument("r1_path", type=click.Path(exists=True))
@click.argument("r2_path", type=click.Path(exists=True))
@click.argument("csv_path", type=click.Path(exists=True))
@click.argument("output_dir", type=click.Path(exists=True))
def run_rna_map(fasta_path, r1_path, r2_path, csv_path, output_dir):
    setup_logging()
    BasicTasks.rna_map(fasta_path, r1_path, r2_path, csv_path, output_dir)


@time_it
@cli.command()
@click.argument("barcode_seq")
@click.argument("construct")
def rna_map_combine(barcode_seq, construct):
    setup_logging()
    df = pd.read_csv("data.csv")
    row = get_data_row(df, barcode_seq, construct)
    if row is None:
        log.error(f"no barcode_seq {barcode_seq} with construct {construct} found")
        return
    BasicTasks.rna_map_combine(row)


@time_it
@cli.command()
def join_fastq_files():
    setup_logging()
    os.makedirs(f"demultiplexed", exist_ok=True)
    df = pd.read_csv("data.csv")
    for barcode, g in df.groupby("barcode_seq"):
        os.makedirs(f"demultiplexed/{barcode}", exist_ok=True)
        r1_files = glob.glob(f"data/*/{barcode}/test_R1.fastq.gz")
        r2_files = glob.glob(f"data/*/{barcode}/test_R2.fastq.gz")
        log.info(f"joining {barcode} files")
        log.info(f"r1_files: {len(r1_files)}")
        log.info(f"outputing to: demultiplexed/{barcode}/test_R1.fastq.gz")
        os.system(
            f"cat {' '.join(r1_files)} > demultiplexed/{barcode}/test_R1.fastq.gz"
        )
        log.info(f"r2_files: {len(r2_files)}")
        log.info(f"outputing to: demultiplexed/{barcode}/test_R2.fastq.gz")
        os.system(
            f"cat {' '.join(r2_files)} > demultiplexed/{barcode}/test_R2.fastq.gz"
        )


################################################################################
############################ With internal barcodes ############################
################################################################################


@time_it
@cli.command()
@click.argument("construct_barcode")
@click.argument("b1_seq")
@click.argument("b2_seq")
@click.argument("b1_min_pos", type=int)
@click.argument("b1_max_pos", type=int)
@click.argument("b2_min_pos", type=int)
@click.argument("b2_max_pos", type=int)
def int_demultiplex(
    construct_barcode, b1_seq, b2_seq, b1_min_pos, b1_max_pos, b2_min_pos, b2_max_pos
):
    setup_logging()
    IntDemultiplexTasks.int_demultiplex(
        construct_barcode,
        b1_seq,
        b2_seq,
        b1_min_pos,
        b1_max_pos,
        b2_min_pos,
        b2_max_pos,
    )


@time_it
@cli.command()
@click.argument("code")
@click.argument("lib_barcode_seq")
@click.argument("construct_barcode_seq")
def int_demultiplex_rna_map(code, lib_barcode_seq, construct_barcode_seq):
    setup_logging()
    # get fastq files
    fastq_dir = f"int-demultiplexed/{lib_barcode_seq}"
    log.info(fastq_dir)
    log.info(f"{fastq_dir}/{construct_barcode_seq}_mate1.fastq.gz")
    mate_1_path = glob.glob(f"{fastq_dir}/{construct_barcode_seq}_mate1.fastq.gz")[0]
    mate_2_path = glob.glob(f"{fastq_dir}/{construct_barcode_seq}_mate2.fastq.gz")[0]
    mate_1_path = os.path.abspath(mate_1_path)
    mate_2_path = os.path.abspath(mate_2_path)
    # check to make sure files actually have stuff in them
    fsize_1 = get_file_size(mate_1_path)
    fsize_2 = get_file_size(mate_2_path)
    if fsize_1 < 100 or fsize_2 < 100:
        log.warning(f"skipping {construct_barcode_seq} because file size is too small")
        return
    # get sequences to make input files
    df = pd.read_csv("data.csv")
    df_sub = df.query("barcode_seq == @lib_barcode_seq and code == @code")
    if len(df_sub) == 0:
        log.error(f"no barcode_seq {lib_barcode_seq} with code {code} found")
        return
    row = df_sub.iloc[0]
    df_barcode = pd.read_json(f"inputs/barcode_jsons/{row['code']}.json")
    df_barcode = df_barcode[df_barcode["full_barcode"] == construct_barcode_seq]
    tmp_dir = "/scratch/" + random_string(10)
    os.makedirs(tmp_dir, exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())
    to_fasta(to_dna(df_barcode), f"{tmp_dir}/test.fasta")
    df_barcode = df_barcode[["name", "sequence", "structure"]]
    df_barcode.to_csv(f"{tmp_dir}/test.csv", index=False)
    os.chdir(tmp_dir)
    # update params
    params = get_preset_params("barcoded-library")
    # run rna-map
    try:
        rna_map.run.run(f"test.fasta", mate_1_path, mate_2_path, "test.csv", params)
    except:
        log.error(f"rna-map failed for {construct_barcode_seq}")
        return
    # move to unique name to avoid collisions
    shutil.move(
        f"output/BitVector_Files/mutation_histos.p",
        f"{cur_dir}/int-demultiplexed-rna-map/{lib_barcode_seq}/mutation_histos_{construct_barcode_seq}.p",
    )
    os.chdir(cur_dir)
    shutil.rmtree(tmp_dir)


@time_it
@cli.command()
@click.argument("barcode_seq")
@click.argument("rna_name")
def int_demultiplex_rna_map_combine(barcode_seq, rna_name):
    setup_logging()
    df = pd.read_csv("data.csv")
    row = get_data_row(df, barcode_seq, rna_name)
    if row is None:
        return
    run_path = "results/" + row["run_name"]
    dir_name = row["construct"] + "_" + row["code"] + "_" + row["data_type"]
    final_path = f"{run_path}/processed/{dir_name}/output/BitVector_Files/"
    log.info(f"results path: {final_path}")
    os.makedirs(final_path, exist_ok=True)
    mut_histo_files = glob.glob(f"int-demultiplexed-rna-map/{barcode_seq}/*p")
    log.info(f"found {len(mut_histo_files)} files")
    merged_mut_histos = {}
    i = 0
    for mut_hist_file in mut_histo_files:
        if i % 100 == 0:
            log.info(f"merged {i} mut histos")
        merge_mut_histo_dicts(
            merged_mut_histos, get_mut_histos_from_pickle_file(mut_hist_file)
        )
        i += 1
    df_results = get_mut_histo_dataframe(merged_mut_histos)
    df_results["rna_name"] = rna_name
    cols = list(row.keys())
    for c in "demult_cmd,length".split(","):
        cols.remove(c)
    for c in cols:
        df_results[c] = row[c]
    df_results.to_json(final_path + "mutation_histos.json", orient="records")
    df_results = df_results.sort_values("num_aligned", ascending=False)
    generate_pop_avg_plots(df_results, row["run_name"], dir_name)
    write_mut_histos_to_pickle_file(merged_mut_histos, final_path + "mutation_histos.p")


# single ###########################################################################


@time_it
@cli.command()
@click.argument("run_name")
def single_run(run_name):
    setup_logging(file_name="demultiplex.log")
    pass


################################################################################
############################## Summary functions ###############################
################################################################################


if __name__ == "__main__":
    cli()
