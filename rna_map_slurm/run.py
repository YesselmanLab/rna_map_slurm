import click

from fastqsplitter import split_fastqs as fastqsplitter


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
    r1_output_files = []
    r2_output_files = []
    for i in range(start, num_chunks + start):
        r1_output_files.append(f"{output_dir}/split-{i:04}/test_R1.fastq.gz")
        r2_output_files.append(f"{output_dir}/split-{i:04}/test_R2.fastq.gz")
    fastqsplitter(r1_path, r1_output_files, threads_per_file=threads)
    fastqsplitter(r2_path, r2_output_files, threads_per_file=threads)


if __name__ == "__main__":
    cli()
