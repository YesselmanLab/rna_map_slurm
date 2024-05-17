import os
import click
import pandas as pd

from fastqsplitter import split_fastqs as fastqsplitter

from rna_map_slurm.logger import setup_logging, get_logger
from rna_map_slurm.fastq import PairedFastqFiles, FastqFile
from rna_map_slurm.demultiplex import SabreDemultiplexer

log = get_logger(__name__)


@click.group()
def cli():
    pass


@cli.command()
@click.argument("r1_path", type=click.Path(exists=True))
@click.argument("r2_path", type=click.Path(exists=True))
@click.argument("output_dir", type=click.Path(exists=True))
@click.argument("num_chunks", type=int)
@click.option("--start", default=0)
@click.option("--threads", default=1)
def split_fastqs(r1_path, r2_path, output_dir, num_chunks, start, threads):
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


if __name__ == "__main__":
    cli()
