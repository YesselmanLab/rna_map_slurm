import os
import glob
import shutil
import click
import pandas as pd
import matplotlib.pyplot as plt
from typing import Callable, Any
import datetime

from fastqsplitter import split_fastqs as fastqsplitter

from seq_tools.dataframe import to_fasta, to_dna
from seq_tools.sequence import get_reverse_complement

from barcode_demultiplex.demultiplex import demultiplex as barcode_demultiplex

import rna_map
from rna_map.mutation_histogram import (
    merge_mut_histo_dicts,
    get_mut_histos_from_pickle_file,
    write_mut_histos_to_pickle_file,
    get_dataframe,
)
from rna_map.parameters import get_preset_params

from rna_map_slurm.logger import setup_logging, get_logger
from rna_map_slurm.plotting import plot_pop_avg_from_row
from rna_map_slurm.util import random_string
from rna_map_slurm.tasks import BasicTasks

log = get_logger("cli")


def get_file_size(file_path):
    file_path = os.path.realpath(file_path)
    return os.path.getsize(file_path)


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
    end_time = datetime.datetime.now()


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


# combine rna-map results ############################################################
# takes regular rna-map results and combines them into a mutation_histos.p


@time_it
@cli.command()
@click.argument("barcode_seq")
@click.argument("rna_name")
def rna_map_combine(barcode_seq, rna_name):
    setup_logging()
    df = pd.read_csv("data.csv")
    df_sub = df.query("barcode_seq == @barcode_seq and construct == @rna_name")
    if len(df_sub) == 0:
        log.error(
            f"barcode_seq {barcode_seq} with rna {rna_name} not found in data.csv"
        )
        return
    if len(df_sub) > 1:
        log.warning(
            f"barcode_seq {barcode_seq} with rna {rna_name} has multiple entries in data.csv"
        )
    row = df_sub.iloc[0]
    run_path = "results/" + row["run_name"]
    dir_name = row["construct"] + "_" + row["code"] + "_" + row["data_type"]
    final_path = f"{run_path}/processed/{dir_name}/output/BitVector_Files/"
    log.info(f"results path: {final_path}")
    os.makedirs(run_path, exist_ok=True)
    os.makedirs(run_path, exist_ok=True)
    os.makedirs(final_path, exist_ok=True)
    dirs = glob.glob("data/split-*")
    merged_mut_histos = None
    count_files = 0
    for d in dirs:
        mhs_path = (
            f"{d}/{barcode_seq}/{rna_name}/output/BitVector_Files/mutation_histos.p"
        )
        if not os.path.isfile(mhs_path):
            log.warning("files not found:" + mhs_path)
            continue
        if merged_mut_histos is None:
            try:
                merged_mut_histos = get_mut_histos_from_pickle_file(mhs_path)
                count_files += 1
            except:
                log.warning(f"could not open file: {mhs_path}")
        else:
            try:
                merge_mut_histo_dicts(
                    merged_mut_histos, get_mut_histos_from_pickle_file(mhs_path)
                )
                count_files += 1
            except:
                log.warning(f"could not open file: {mhs_path}")
    log.info(f"merged {count_files} files")
    cols = [
        "name",
        "sequence",
        "structure",
        "pop_avg",
        "sn",
        "num_reads",
        "num_aligned",
        "no_mut",
        "1_mut",
        "2_mut",
        "3_mut",
        "3plus_mut",
    ]
    df_results = get_dataframe(merged_mut_histos, cols)
    df_results.rename(columns={"pop_avg": "data"}, inplace=True)
    df_results.to_json(final_path + "mutation_histos.json", orient="records")
    i = 0
    for _, row in df_results.iterrows():
        plot_pop_avg_from_row(row)
        plt.title("num_aligned: " + str(row["num_aligned"]) + "\tsn: " + str(row["sn"]))
        plt.savefig(final_path + f"{row['name']}.png")
        plt.close()
        i += 1
        if i > 100:
            break
        shutil.copy(
            final_path + f"{row['name']}.png",
            f"results/pop_avg_pngs/{dir_name}_{row['name']}.png",
        )
    write_mut_histos_to_pickle_file(merged_mut_histos, final_path + "mutation_histos.p")


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
    tmp_dir = "/scratch/" + random_string(10)
    print(tmp_dir)
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
        f"seqkit grep -f {tmp_dir}/common_names.txt {tmp_dir}/test_R2.fastq.gz -o int_demultiplexed/{construct_barcode}/{b1_seq}_{b2_seq}_mate1.fastq.gz"
    )
    os.system(
        f"seqkit grep -f {tmp_dir}/common_names.txt {tmp_dir}/test_R1.fastq.gz -o int_demultiplexed/{construct_barcode}/{b1_seq}_{b2_seq}_mate2.fastq.gz"
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
        f"{tmp_dir}/output/BitVector_Files/mutation_histos.p",
        f"int_demultiplexed_rna_map/{lib_barcode_seq}/mutation_histos_{construct_barcode_seq}.p",
    )
    os.chdir(cur_dir)
    shutil.rmtree(tmp_dir)


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
