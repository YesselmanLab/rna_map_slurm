import os
import pandas as pd
from typing import List, Dict

from rna_map_slurm.logger import get_logger
from rna_map_slurm.jobs import SlurmOptions, get_job_header, generate_submit_file
from rna_map_slurm.fastq import get_paired_fastqs, PairedFastqFiles

log = get_logger(__name__)


# TODO split into two jobs for each pair of fastq files
def generate_split_fastq_jobs(
    pfqs: List[PairedFastqFiles], params: Dict
) -> pd.DataFrame:
    os.makedirs("jobs/split-fastq", exist_ok=True)
    # generate jobs
    cur_dir = os.path.abspath(os.getcwd())
    slurm_params = params["slurm_options"]["split_fastq"]
    jobs = []
    for i, pfq in enumerate(pfqs):
        name = f"split-fastq-{i:04}"
        slurm_opts = SlurmOptions(
            name,
            slurm_params["time"],
            slurm_params["mem-per-cpu"],
            slurm_params["cpus-per-task"],
            params["slurm_options"]["extra_header_cmds"],
        )
        job_header = get_job_header(slurm_opts, os.path.abspath("jobs/split-fastq/"))
        job_body = (
            f"rna-map-slurm-runner split-fastqs {pfq.read_1.path} {pfq.read_2.path} "
            f"{cur_dir + '/data'} {params['fastq_chunks']} --start "
            f"{i*params['fastq_chunks']}\n"
        )
        job = job_header + job_body
        f = open(f"jobs/split-fastq/{name}.sh", "w")
        f.write(job)
        f.close()
        jobs.append([f"jobs/split_fastq/{name}.sh", "SPLIT_FASTQ", ""])
    df_jobs = pd.DataFrame(jobs, columns=["job", "type", "status"])
    generate_submit_file("submits/README_SPLIT_FASTQ", df_jobs["job"].tolist())
    return df_jobs
