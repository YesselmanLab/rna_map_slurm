import os
import glob
import shutil
import click
import gzip
import zipfile
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from fastqsplitter import split_fastqs as fastqsplitter

from seq_tools.dataframe import to_fasta, to_dna

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
from rna_map_slurm.fastq import PairedFastqFiles, FastqFile, get_paired_fastqs
from rna_map_slurm.demultiplex import SabreDemultiplexer
from rna_map_slurm.plotting import plot_pop_avg_from_row
from rna_map_slurm.util import random_string, gzip_files, flatten_and_zip_directory

log = get_logger(__name__)


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
    if len(mate1_files) == 0:
        log.warning(f"no files found for {construct_barcode}")
        return
    if len(mate1_files) != len(mate2_files):
        log.warning(f"mismatched files for {construct_barcode}")
        return
    log.info(f"{construct_barcode} {len(mate1_files)} {len(mate2_files)}")
    combine_gzipped_fastq(mate1_files, f"{outdir}/{pair[0]}")
    combine_gzipped_fastq(mate2_files, f"{outdir}/{pair[1]}")
    subprocess.call(f"rm -r {tmp_dir}/*/{pair[0]}", shell=True)
    subprocess.call(f"rm -r {tmp_dir}/*/{pair[1]}", shell=True)


def get_file_size(file_path):
    file_path = os.path.realpath(file_path)
    return os.path.getsize(file_path)


@click.group()
def cli():
    pass


################################################################################
############################ Use everytime function ############################
################################################################################


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
    r1_output_files = []
    r2_output_files = []
    for i in range(start, num_chunks + start):
        r1_output_files.append(f"{output_dir}/split-{i:04}/test_R1.fastq.gz")
        r2_output_files.append(f"{output_dir}/split-{i:04}/test_R2.fastq.gz")
    fastqsplitter(r1_path, r1_output_files, threads_per_file=threads)
    fastqsplitter(r2_path, r2_output_files, threads_per_file=threads)


@cli.command()
@click.argument("csv")
@click.argument("r1_path", type=click.Path(exists=True))
@click.argument("r2_path", type=click.Path(exists=True))
@click.option("--output_dir", type=click.Path(exists=True))
def demultiplex(csv, r1_path, r2_path, output_dir):
    """
    demultiplexes paired fastq files given 3' end barcodes
    """
    setup_logging(file_name="demultiplex.log")
    if output_dir is None:
        output_dir = os.getcwd()
    paired_fastqs = PairedFastqFiles(FastqFile(r1_path), FastqFile(r2_path))
    df = pd.read_csv(csv)
    demultiplexer = SabreDemultiplexer()
    demultiplexer.run(df, paired_fastqs, output_dir)


# combine rna-map results ############################################################
# takes regular rna-map results and combines them into a mutation_histos.p
@cli.command()
@click.argument("barcode_seq")
@click.argument("rna_name")
def combine_rna_map(barcode_seq, rna_name):
    setup_logging()
    df = pd.read_csv("data.csv")
    df_sub = df.query(f"barcode_seq == {barcode_seq} and rna_name == {rna_name}")
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
    os.makedirs("results", exist_ok=True)
    os.makedirs("results/pop_avg_pngs", exist_ok=True)
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
        plt.savefig(final_path + f"{row['name']}.png")
        plt.close()
        i += 1
        if i > 100:
            break
    shutil.copy(
        final_path + f"{row['name']}.png",
        f"results/pop_avg_pngs/{row['rna_name']}_{row['name']}.png",
    )
    write_mut_histos_to_pickle_file(merged_mut_histos, final_path + "mutation_histos.p")


################################################################################
############################ With internal barcodes ############################
################################################################################


@cli.command()
@click.argument("fastq_dir", type=click.Path(exists=True))
@click.option("--output_dir", default=None)
def int_demultiplex(fastq_dir, output_dir):
    setup_logging(file_name=f"{fastq_dir}/int_demultiplex.log")
    if output_dir is None:
        output_dir = os.getcwd()
    df = pd.read_csv(f"data.csv")
    read_1_path = fastq_dir + "/test_R1.fastq.gz"
    read_2_path = fastq_dir + "/test_R2.fastq.gz"
    barcode_seq = Path(fastq_dir).stem
    df_sub = df.loc[df["barcode_seq"] == barcode_seq]
    if df_sub.empty:
        log.error(f"barcode_seq {barcode_seq} not found in data.csv")
        return
    row = df_sub.iloc[0]
    # get helices from commandline
    helices = []
    args = row["demult_cmd"].split()
    for i in range(0, len(args)):
        if args[i] == "--helix" or args[i] == "-helix":
            helices.append([int(args[i + 1]), int(args[i + 2]), int(args[i + 3])])
    unique_code = random_string(10)
    data_path = f"{output_dir}/{unique_code}"
    df_seq = pd.read_csv(f"inputs/rnas/{row['code']}.csv")
    barcode_demultiplex(
        df_seq, Path(read_2_path), Path(read_1_path), helices, data_path
    )
    zip_path = f"{fastq_dir}/int_demultiplexed.zip"
    flatten_and_zip_directory(data_path, zip_path)
    shutil.rmtree(data_path)


@cli.command()
@click.argument("home_dir", type=click.Path(exists=True))
@click.argument("barcode_seq")
@click.option("--threads", default=1)
@click.option("--tmp_dir", default=None)
def join_int_demultiplex(home_dir, barcode_seq, threads, tmp_dir):
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


@cli.command()
@click.argument("home_dir", type=click.Path(exists=True))
@click.argument("lib_barcode_seq")
@click.argument("construct_barcode_seq")
@click.option("--tmp_dir", default=None)
def rna_map_single_barcode(home_dir, lib_barcode_seq, construct_barcode_seq, tmp_dir):
    setup_logging()
    if tmp_dir is None:
        tmp_dir = os.getcwd()
    # get fastq files
    fastq_dir = f"{home_dir}/joined_fastqs/{lib_barcode_seq}"
    log.info(fastq_dir)
    log.info(f"{fastq_dir}/{construct_barcode_seq}_mate1.fastq.gz")
    mate_1_path = glob.glob(f"{fastq_dir}/{construct_barcode_seq}_mate1.fastq.gz")[0]
    mate_2_path = glob.glob(f"{fastq_dir}/{construct_barcode_seq}_mate2.fastq.gz")[0]
    # check to make sure files actually have stuff in them
    fsize_1 = get_file_size(mate_1_path)
    fsize_2 = get_file_size(mate_2_path)
    if fsize_1 < 100 or fsize_2 < 100:
        log.warning(f"skipping {construct_barcode_seq} because file size is too small")
        return
    # get sequences to make input files
    df = pd.read_csv("data.csv")
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


# single ###########################################################################


@cli.command()
@click.argument("run_name")
def single_run(run_name):
    setup_logging(file_name="demultiplex.log")

    paired_fastqs = PairedFastqFiles(FastqFile(r1_path), FastqFile(r2_path))
    df = pd.read_csv(csv)
    demultiplexer = SabreDemultiplexer()
    demultiplexer.run(df, paired_fastqs, "demultiplexed")


################################################################################
############################## Summary functions ###############################
################################################################################


# TODO needs be refactored to grab data from that transfered in and doesnt do its own copying
@cli.command()
@click.option("--combine-all", is_flag=True)
def fastq_concat(combine_all):
    setup_logging()
    os.makedirs(f"demultiplexed", exist_ok=True)
    df = pd.read_csv("data.csv")
    seq_path = os.environ["SEQPATH"]
    for barcode, g in df.groupby("barcode_seq"):
        rna_count = 0
        for row in g.iterrows():
            try:
                df_rna = pd.read_csv(f"{seq_path}/rna/{row['code']}.csv")
                if len(df_rna) < 1000 and not combine_all:
                    rna_count += 1
            except:
                rna_count += 1
                continue
        if rna_count == 0:
            continue
        log.info(f"in {barcode} {rna_count} rnas are  present")
        os.makedirs(f"demultiplexed/{barcode}", exist_ok=True)
        r1_files = glob.glob(f"data/*/{barcode}/test_R1.fastq.gz")
        r2_files = glob.glob(f"data/*/{barcode}/test_R2.fastq.gz")
        os.system(
            f"cat {' '.join(r1_files)} > demultiplexed/{barcode}/test_R1.fastq.gz"
        )
        os.system(
            f"cat {' '.join(r2_files)} > demultiplexed/{barcode}/test_R2.fastq.gz"
        )


if __name__ == "__main__":
    cli()
