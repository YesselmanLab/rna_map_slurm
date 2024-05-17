import os
import glob
import click
import pandas as pd
import matplotlib.pyplot as plt

from fastqsplitter import split_fastqs as fastqsplitter

from rna_map.mutation_histogram import (
    merge_mut_histo_dicts,
    get_mut_histos_from_pickle_file,
    write_mut_histos_to_pickle_file,
    get_dataframe,
)

from rna_map_slurm.logger import setup_logging, get_logger
from rna_map_slurm.fastq import PairedFastqFiles, FastqFile
from rna_map_slurm.demultiplex import SabreDemultiplexer
from rna_map_slurm.plotting import plot_pop_avg_from_row
import click

log = get_logger(__name__)


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
    df_sub = df[df["barcode_seq"] == barcode_seq]
    if len(df_sub) == 0:
        log.error(f"barcode_seq {barcode_seq} not found in data.csv")
        return
    row = df_sub.iloc[0]
    os.makedirs("results", exist_ok=True)
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
    write_mut_histos_to_pickle_file(merged_mut_histos, final_path + "mutation_histos.p")


################################################################################
############################ With internal barcodes ############################
################################################################################


################################################################################
############################## Summary functions ###############################
################################################################################


# TODO needs be refactored to grab data from that transfered in and doesnt do its own copying
@cli.command()
@click.option("--combine-all", is_flag=True)
def fastq_concat(combine_all):
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
