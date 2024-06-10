from celery import Celery, group

from rna_map_slurm import app

from fastqsplitter import split_fastqs as fastqsplitter

import time

@app.task
def master_task(tasks):
    tasks = [
        ["/work/yesselmanlab/jyesselm/data/2024_06_03_Nextseq_Run9/2024_06_03_Nextseq_Run9_S1_L001_R1_001.fastq.gz", "data", 20, 0],
        ["/work/yesselmanlab/jyesselm/data/2024_06_03_Nextseq_Run9/2024_06_03_Nextseq_Run9_S1_L001_R2_001.fastq.gz", "data", 20, 0]
    ]
    result = group(tasks).apply_async()
    return result

@app.task
def split_fastq_file(fastq_file, output_dir, num_chunks, start):
    # determine if fastq_file is R1 or R2
    output_file = "test_R1.fastq.gz"
    if "R2" in fastq_file:
        output_file = "test_R2.fastq.gz"
    output_files = []
    for i in range(start, num_chunks + start):
        output_files.append(f"{output_dir}/split-{i:04}/{output_file}")
    #fastqsplitter(fastq_file, output_files, threads_per_file=1)
    time.sleep(10)
    return output_files
