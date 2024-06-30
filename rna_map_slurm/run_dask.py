import os
import time
import glob
import shutil
import pandas as pd
import click
from dask.distributed import Client
from dask_jobqueue import SLURMCluster

from barcode_demultiplex.demultiplex import find_helix_barcodes

from rna_map_slurm.logger import setup_logging, get_logger
from rna_map_slurm.tasks import BasicTasks
from rna_map_slurm.fastq import get_paired_fastqs


log = get_logger("RUN-DASK")

# dask tasks ##################################################################


def split_fastq_file_task(fastq_file, output_dir, num_chunks, start):
    result = BasicTasks.split_fastq_file(fastq_file, output_dir, num_chunks, start)
    time.sleep(2)  # add small time buffer
    return result


def demultiplex_task(data_csv, r1_fastq, r2_fastq, output_dir):
    result = BasicTasks.demultiplex(data_csv, r1_fastq, r2_fastq, output_dir)
    time.sleep(2)  # add small time buffer
    return result


def join_fastq_files_task(fastq_files, joined_fastq):
    os.system(f"cat {' '.join(fastq_files)} > {joined_fastq}")
    time.sleep(1)
    return joined_fastq


def rna_map_task(fa, r1_fastq, r2_fastq, csv, output_dir):
    result = None
    try:
        result = BasicTasks.rna_map(fa, r1_fastq, r2_fastq, csv, output_dir)
    except Exception as e:
        log.error("rna_map_task raised an exception")
        log.error("An error occurred", exc_info=True)
    time.sleep(1)
    return result


# setup tasks ################################################################


def setup_split_fastq_file_tasks(all_pfqs, params):
    split_tasks = []
    for i, pfq in enumerate(all_pfqs):
        split_tasks.append(
            [
                pfq.read_1.path,
                os.path.abspath("data"),
                params["num_splits"],
                i * params["num_splits"],
            ]
        )
        split_tasks.append(
            [
                pfq.read_2.path,
                os.path.abspath("data"),
                params["num_splits"],
                i * params["num_splits"],
            ]
        )
    return split_tasks


def setup_demultiplex_tasks(params):
    demultiplex_tasks = []
    split_dirs = [
        os.path.abspath(f"data/split-{i:04}") for i in range(params["num_dirs"])
    ]
    for split_dir in split_dirs:
        demultiplex_tasks.append(
            [
                os.path.abspath("data.csv"),
                f"{split_dir}/test_R1.fastq.gz",
                f"{split_dir}/test_R2.fastq.gz",
                split_dir,
            ]
        )
    return demultiplex_tasks


def setup_join_fastq_files_tasks(barcodes):
    os.makedirs("demultiplexed", exist_ok=True)
    fastq_join_tasks = []
    for barcode in barcodes:
        os.makedirs(f"demultiplexed/{barcode}", exist_ok=True)
        r1_files = glob.glob(f"data/*/{barcode}/test_R1.fastq.gz")
        r2_files = glob.glob(f"data/*/{barcode}/test_R2.fastq.gz")
        fastq_join_tasks.append([r1_files, f"demultiplexed/{barcode}/test_R1.fastq.gz"])
        fastq_join_tasks.append([r2_files, f"demultiplexed/{barcode}/test_R2.fastq.gz"])
    return fastq_join_tasks


def setup_rna_map_tasks(df, params):
    rna_map_tasks = []
    split_dirs = [f"data/split-{i:04}" for i in range(0, params["num_dirs"])]
    for i, row in df.iterrows():
        fa_path = os.path.abspath(f"inputs/fastas/{row['code']}.fasta")
        dot_bracket_path = os.path.abspath(f"inputs/rnas/{row['code']}.csv")
        for split_dir in split_dirs:
            fastq_dir = os.path.abspath(f"{split_dir}/{row['barcode_seq']}")
            output_dir = f"{fastq_dir}/{row['construct']}"
            os.makedirs(output_dir, exist_ok=True)
            rna_map_tasks.append(
                [
                    fa_path,
                    f"{fastq_dir}/test_R1.fastq.gz",
                    f"{fastq_dir}/test_R2.fastq.gz",
                    dot_bracket_path,
                    output_dir,
                ]
            )
    return rna_map_tasks


# main func ########################################################################


def dask_runner(data_dirs, num_workers, num_splits, debug):
    """
    main function for script
    """ 
    log.info("removing worker output files")
    os.system("rm -rf slurm-*")
    cluster = SLURMCluster(
        processes=1,
        cores=1,
        memory="4GB",
        walltime="5:00:00",
        job_script_prologue=[
            "#!/bin/bash",
            "module load bowtie/2.3 trim_galore fastqc anaconda",
            "conda activate",
            "conda activate dask",
        ],
    )
    params = {
        "data_dirs": data_dirs,
        "num_splits": num_splits,
    }
    cluster.scale(num_workers)  # this may take a few seconds to launch
    # start client
    client = Client(cluster)
    all_pfqs = []
    for d in data_dirs:
        d = os.path.abspath(d)
        all_pfqs.extend(get_paired_fastqs(d))
    params["num_dirs"] = params["num_splits"] * len(all_pfqs)
    # split fastq files ###############################################################
    split_tasks = setup_split_fastq_file_tasks(all_pfqs, params)
    log.info(f"currently {len(split_tasks)} split tasks")
    futures = client.map(lambda args: split_fastq_file_task(*args), split_tasks)
    client.gather(futures)
    log.info("finished with splitting")
    # demultiplex  #####################################################################
    demultiplex_tasks = setup_demultiplex_tasks(params)
    log.info(f"currently {len(demultiplex_tasks)} demultiplex tasks")
    futures = client.map(lambda args: demultiplex_task(*args), demultiplex_tasks)
    client.gather(futures)
    log.info("finished with demultiplexing")
    # recombine fastq files  ###########################################################
    df_single = pd.read_csv("csvs/data-single.csv")
    barcodes = df_single["barcode_seq"].unique()
    fastq_join_tasks = setup_join_fastq_files_tasks(barcodes)
    print(f"currently {len(fastq_join_tasks)} fastq_join tasks")
    futures = client.map(lambda args: join_fastq_files_task(*args), fastq_join_tasks)
    client.gather(futures)
    log.info("finished with joining fastq files")
    # run rna map ######################################################################
    rna_map_tasks = setup_rna_map_tasks(df_single, params)
    print(f"currently {len(rna_map_tasks)} rna_map tasks")
    print(rna_map_tasks[0])
    futures = client.map(lambda args: rna_map_task(*args), rna_map_tasks)
    client.gather(futures)
    print("finished with rna_map tasks")
    print("DONE #####################################################################")
