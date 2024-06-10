from celery import Celery
from fastqsplitter import split_fastqs as fastqsplitter

app = Celery('rna_map_slurm')

@app.task
def split_fastq_file(fastq_file, output_dir, num_chunks, start):
    # determine if fastq_file is R1 or R2
    output_file = "test_R1.fastq.gz"
    if "R2" in fastq_file:
        output_file = "test_R2.fastq.gz"
    output_files = []
    for i in range(start, num_chunks + start):
        output_files.append(f"{output_dir}/split-{i:04}/{output_file}")
    fastqsplitter(fastq_file, output_files, threads_per_file=1)
