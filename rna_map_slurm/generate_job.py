import os
import pandas as pd
import numpy as np
from typing import List, Dict


from rna_map_slurm.logger import get_logger
from rna_map_slurm.jobs import SlurmOptions, get_job_header, generate_submit_file
from rna_map_slurm.fastq import get_paired_fastqs, PairedFastqFiles

log = get_logger("generate-jobs")


def split_into_n(df, n):
    """
    Split a list into n sublists with sizes as close to equal as possible.

    Args:
    - lst (pd.dataframe): The list to be split.
    - n (int): The number of sublists.

    Returns:
    - list of dataframe: The split sublists.
    """
    avg = len(df) // n
    remainder = len(df) % n
    result = []
    idx = 0
    for i in range(n):
        size = avg + (1 if i < remainder else 0)
        result.append(df[idx : idx + size].copy())
        idx += size
    return result


def create_job_header(
    name: str, slurm_params: Dict, extra_header_cmds: str, path: str
) -> str:
    """
    Generates the SLURM job header using the provided parameters.

    Args:
        name (str): The name of the SLURM job.
        slurm_params (Dict): Dictionary containing SLURM parameters such as time, memory
            per CPU, and CPUs per task.
        extra_header_cmds (str): Additional header commands to be included in the job header.
        path (str): The directory path where the job header will be used.

    Returns:
        str: The formatted SLURM job header.

    Raises:
        None

    Examples:
        >>> slurm_params = {
        ...     "time": "00:10:00",
        ...     "mem-per-cpu": "2G",
        ...     "cpus-per-task": 4
        ... }
        >>> extra_header_cmds = "#SBATCH --constraint=skylake"
        >>> path = "/path/to/job"
        >>> create_job_header("my_job", slurm_params, extra_header_cmds, path)
        '#SBATCH --job-name=my_job\\n#SBATCH --time=00:10:00\\n#SBATCH --mem-per-cpu=2G\\n#SBATCH --cpus-per-task=4\\n#SBATCH --constraint=skylake\\n'
    """
    slurm_opts = SlurmOptions(
        name,
        slurm_params["time"],
        slurm_params["mem-per-cpu"],
        slurm_params["cpus-per-task"],
        extra_header_cmds,
    )
    return get_job_header(slurm_opts, os.path.abspath(path))


def write_job_file(path: str, job_name: str, job_content: str) -> None:
    """
    Writes the job content to a specified file.

    Parameters:
    path (str): The directory path where the job file will be created.
    job_name (str): The name of the job file.
    job_content (str): The content of the job script to be written into the file.

    Returns:
    None
    """
    with open(f"{path}/{job_name}.sh", "w") as f:
        f.write(job_content)


def generate_job_list(
    path: str, job_type: str, requirement: str, job_names: List[str]
) -> pd.DataFrame:
    """
    Generates a DataFrame containing the job list details.

    Parameters:
    path (str): The directory path where the job files are located.
    job_type (str): The type of the job.
    requirement (str): The job requirement (dependency).
    job_names (List[str]): List of job file names.

    Returns:
    pd.DataFrame: DataFrame containing job details with columns
      ['job_path', 'job_type', 'job_requirement'].
    """
    jobs = [[f"{path}/{name}.sh", job_type, requirement] for name in job_names]
    return pd.DataFrame(jobs, columns=["job_path", "job_type", "job_requirement"])


################################################################################
############################ Use everytime function ############################
################################################################################


def generate_split_fastq_jobs(
    pfqs: List[PairedFastqFiles], params: Dict
) -> pd.DataFrame:
    """
    Generates SLURM job files for splitting FASTQ files and returns a DataFrame of job details.

    Args:
        pfqs (List[PairedFastqFiles]): List of paired FASTQ files.
        params (Dict): Parameters for the SLURM job and FASTQ splitting.

    Returns:
        pd.DataFrame: DataFrame containing job details.
    """
    job_name = "split-fastq"
    os.makedirs(f"jobs/{job_name}", exist_ok=True)
    for i in range(0, params["num_dirs"]):
        os.makedirs(f"data/split-{i:04}", exist_ok=True)
    slurm_params = params["slurm_options"][job_name]
    extra_header_cmds = params["slurm_options"]["extra-header-cmds"]
    cur_dir = os.path.abspath(os.getcwd())
    job_names = []
    for i, pfq in enumerate(pfqs):
        name = f"{job_name}-{i:04}"
        job_header = create_job_header(
            name, slurm_params, extra_header_cmds, f"jobs/{job_name}"
        )
        job_body = (
            f"rna-map-slurm-runner split-fastqs {pfq.read_1.path} {pfq.read_2.path} "
            f"{os.path.join(cur_dir, 'data')} {params['fastq_chunks']} --start {i * params['fastq_chunks']}\n"
        )
        write_job_file(f"jobs/{job_name}", name, job_header + job_body)
        job_names.append(name)

    df_jobs = generate_job_list(f"jobs/{job_name}", job_name, "", job_names)
    generate_submit_file(
        f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist()
    )
    return df_jobs


def generate_trim_galore_jobs(params: Dict) -> pd.DataFrame:
    job_name = "trim-galore"
    os.makedirs(f"jobs/{job_name}", exist_ok=True)
    runs_per_job = params["tasks_per_job"][job_name]
    slurm_params = params["slurm_options"][job_name]
    extra_header_cmds = params["slurm_options"]["extra-header-cmds"]
    dirs = [os.path.abspath(f"data/split-{i:04}") for i in range(params["num_dirs"])]
    dir_groups = [dirs[i : i + runs_per_job] for i in range(0, len(dirs), runs_per_job)]
    cur_dir = os.path.abspath(os.getcwd())
    job_names = []
    for i, dg in enumerate(dir_groups):
        name = f"{job_name}-{i:04}"
        job_header = create_job_header(
            name, slurm_params, extra_header_cmds, f"jobs/{job_name}"
        )
        job_body = ""
        for d in dg:
            job_body += f"cd {d}\n"
            job_body += f"trim_galore --quality 0 --paired {d}/test_R1.fastq.gz {d}/test_R2.fastq.gz\n"
            job_body += f"cd {cur_dir}\n\n"
        write_job_file(f"jobs/{job_name}", name, job_header + job_body)
        job_names.append(name)
    df_jobs = generate_job_list(f"jobs/{job_name}", job_name, "split-fastq", job_names)
    generate_submit_file(
        f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist()
    )
    return df_jobs


def generate_demultiplexing_jobs(params: Dict) -> pd.DataFrame:
    job_name = "demultiplex"
    os.makedirs(f"jobs/{job_name}", exist_ok=True)
    runs_per_job = params["tasks_per_job"][job_name]
    slurm_params = params["slurm_options"][job_name]
    extra_header_cmds = params["slurm_options"]["extra-header-cmds"]
    dirs = [os.path.abspath(f"data/split-{i:04}") for i in range(params["num_dirs"])]
    dir_groups = [dirs[i : i + runs_per_job] for i in range(0, len(dirs), runs_per_job)]
    csv_path = os.path.abspath("data.csv")
    job_names = []
    for i, dg in enumerate(dir_groups):
        name = f"{job_name}-{i:04}"
        job_header = create_job_header(
            name, slurm_params, extra_header_cmds, f"jobs/{job_name}"
        )
        job_body = "\n".join(
            [
                f"rna-map-slurm-runner demultiplex {csv_path} {d}/test_R1.fastq.gz {d}/test_R2.fastq.gz {d}\n"
                for d in dg
            ]
        )
        write_job_file(f"jobs/{job_name}", name, job_header + job_body)
        job_names.append(name)
    df_jobs = generate_job_list(f"jobs/{job_name}", job_name, "split-fastq", job_names)
    generate_submit_file(
        f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist()
    )
    return df_jobs


def generate_rna_map_jobs(params: Dict, df: pd.DataFrame) -> pd.DataFrame:
    """
    These are jobs that do not do the additional internal demultiplexing required
    of some large libraries
    """
    job_name = "rna-map"
    os.makedirs(f"jobs/{job_name}", exist_ok=True)
    runs_per_job = params["tasks_per_job"][job_name]
    slurm_params = params["slurm_options"][job_name]
    extra_header_cmds = params["slurm_options"]["extra-header-cmds"]
    dirs = [f"data/split-{i:04}" for i in range(0, params["num_dirs"])]
    dir_groups = [dirs[i : i + runs_per_job] for i in range(0, len(dirs), runs_per_job)]
    job_names = []
    i = 0
    for _, row in df.iterrows():
        if not os.path.isfile(f"inputs/fastas/{row['code']}.fasta"):
            log.warning(f"{row['code']} does not have a RNA CSV file!!!")
            continue
        for dg in dir_groups:
            name = f"{job_name}-{i:04}"
            job_header = create_job_header(
                name, slurm_params, extra_header_cmds, f"jobs/{job_name}"
            )
            job_body = ""
            fa_path = os.path.abspath(f"inputs/fastas/{row['code']}.fasta")
            dot_bracket_path = os.path.abspath(f"inputs/rnas/{row['code']}.csv")
            for dir in dg:
                output_dir = os.path.abspath(
                    f"{dir}/{row['barcode_seq']}/{row['construct']}"
                )
                os.makedirs(output_dir, exist_ok=True)
                fq1_path = os.path.abspath(
                    f"{dir}/{row['barcode_seq']}/test_R1.fastq.gz"
                )
                fq2_path = os.path.abspath(
                    f"{dir}/{row['barcode_seq']}/test_R2.fastq.gz"
                )
                job_body += f"rna-map-slurm-runner run-rna-map {fa_path} {fq2_path} {fq1_path} {dot_bracket_path} {output_dir}\n\n"
            write_job_file(f"jobs/{job_name}", name, job_header + job_body)
            job_names.append(name)
            i += 1
    df_jobs = generate_job_list(f"jobs/{job_name}", job_name, "demultiplex", job_names)
    generate_submit_file(
        f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist()
    )
    return df_jobs


def generate_rna_map_combine_jobs(params: Dict, df: pd.DataFrame) -> pd.DataFrame:
    job_name = "rna-map-combine"
    os.makedirs(f"jobs/{job_name}", exist_ok=True)
    slurm_params = params["slurm_options"][job_name]
    extra_header_cmds = params["slurm_options"]["extra-header-cmds"]
    job_names = []
    for i, row in df.iterrows():
        name = f"rna-map-combine-{i:04}"
        job_header = create_job_header(
            name, slurm_params, extra_header_cmds, f"jobs/{job_name}"
        )
        job_body = f"rna-map-slurm-runner rna-map-combine {row['barcode_seq']} {row['construct']}\n"
        write_job_file(f"jobs/{job_name}", name, job_header + job_body)
        job_names.append(name)
    df_jobs = generate_job_list(f"jobs/{job_name}", job_name, "rna-map", job_names)
    generate_submit_file(
        f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist()
    )
    return df_jobs


# maybe need to generalize this but now just 1 job
def generate_join_fastq_files_jobs(params: Dict) -> pd.DataFrame:
    job_name = "join-fastq-files"
    os.makedirs(f"jobs/{job_name}", exist_ok=True)
    slurm_params = params["slurm_options"][job_name]
    extra_header_cmds = params["slurm_options"]["extra-header-cmds"]
    name = f"{job_name}-0000"
    job_header = create_job_header(
        name, slurm_params, extra_header_cmds, f"jobs/{job_name}"
    )
    job_body = f"rna-map-slurm-runner join-fastq-files\n"
    write_job_file(f"jobs/{job_name}", name, job_header + job_body)
    df_jobs = generate_job_list(f"jobs/{job_name}", job_name, "demultiplex", [name])
    return df_jobs


################################################################################
############################ With internal barcodes ############################
################################################################################


def generate_int_demultiplex_jobs(params, df):
    job_name = "int-demultiplex"
    os.makedirs(f"jobs/{job_name}", exist_ok=True)
    os.makedirs("int-demultiplexed", exist_ok=True)
    runs_per_job = params["tasks_per_job"][job_name]
    slurm_params = params["slurm_options"][job_name]
    extra_header_cmds = params["slurm_options"]["extra-header-cmds"]
    add_dfs = []
    for _, row in df.iterrows():
        os.makedirs(f"int-demultiplexed/{row['barcode_seq']}", exist_ok=True)
        df_barcodes = pd.read_json(f"inputs/barcode_jsons/{row['code']}.json")
        df_barcodes["construct_barcode"] = row["barcode_seq"]
        add_dfs.append(df_barcodes)
    df_barcodes = pd.concat(add_dfs)
    job_names = []
    dfs = [
        df_barcodes[i : i + runs_per_job]
        for i in range(0, len(df_barcodes), runs_per_job)
    ]
    for i, df in enumerate(dfs):
        name = f"{job_name}-{i:04}"
        job_header = create_job_header(
            name, slurm_params, extra_header_cmds, f"jobs/{job_name}"
        )
        job_body = ""
        for _, group in df.groupby("full_barcode"):
            row = group.iloc[0]
            # TODO check to make sure this is true but I think it works out
            # since the barcodes are at a constant position from the 3' end
            bb1 = row["barcode_bounds"][0][0]
            bb2 = row["barcode_bounds"][0][1]
            # map barcode location to the other read as it starts from the end
            end_len = len(row["sequence"])
            max_len = end_len - bb2[0]
            min_len = end_len - bb2[1]
            bb2 = [min_len, max_len]
            b1_seq = row["barcodes"][0][0]
            b2_seq = row["barcodes"][0][1]
            job_body += (
                f"rna-map-slurm-runner int-demultiplex {row['construct_barcode']} "
                f"{b1_seq} {b2_seq} {bb1[0]} {bb1[1]} {bb2[0]} {bb2[1]}\n\n"
            )
        write_job_file(f"jobs/{job_name}", name, job_header + job_body)
        job_names.append(name)
    df_jobs = generate_job_list(
        f"jobs/{job_name}", job_name, "join-fastq-files", job_names
    )
    generate_submit_file(
        f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist()
    )
    return df_jobs


def generate_int_demultiplex_rna_map_jobs(params, df):
    job_name = "int-demultiplex-rna-map"
    os.makedirs("int-demultiplexed-rna-map", exist_ok=True)
    os.makedirs(f"jobs/{job_name}", exist_ok=True)
    runs_per_job = params["tasks_per_job"][job_name]
    slurm_params = params["slurm_options"][job_name]
    extra_header_cmds = params["slurm_options"]["extra-header-cmds"]
    runs = []
    for i, row in df.iterrows():
        os.makedirs(f"int-demultiplexed-rna-map/{row['barcode_seq']}", exist_ok=True)
        df_barcode = pd.read_json(f"inputs/barcode_jsons/{row['code']}.json")
        unique_barcodes = df_barcode["full_barcode"].unique()
        for barcode in unique_barcodes:
            runs.append([row["code"], row["barcode_seq"], barcode])
    run_groups = [runs[i : i + runs_per_job] for i in range(0, len(runs), runs_per_job)]
    job_names = []
    for i, group in enumerate(run_groups):
        name = f"{job_name}-{i:04}"
        job_header = create_job_header(
            name, slurm_params, extra_header_cmds, f"jobs/{job_name}"
        )
        job_body = ""
        for run in group:
            job_body += f"rna-map-slurm-runner int-demultiplex-rna-map {run[0]} {run[1]} {run[2]}\n\n"
        write_job_file(f"jobs/{job_name}", name, job_header + job_body)
        job_names.append(name)
    df_jobs = generate_job_list(
        f"jobs/{job_name}", job_name, "int-demultiplex", job_names
    )
    generate_submit_file(
        f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist()
    )
    return df_jobs


def generate_int_demultiplex_rna_map_combine_jobs(params, df):
    job_name = "int-demultiplex-rna-map-combine"
    os.makedirs(f"jobs/{job_name}", exist_ok=True)
    slurm_params = params["slurm_options"][job_name]
    extra_header_cmds = params["slurm_options"]["extra-header-cmds"]
    job_names = []
    for i, row in df.iterrows():
        name = f"{job_name}-{i:04}"
        job_header = create_job_header(
            name, slurm_params, extra_header_cmds, f"jobs/{job_name}"
        )
        job_body = f"rna-map-slurm-runner int-demultiplex-rna-map-combine {row['barcode_seq']} {row['construct']}\n"
        write_job_file(f"jobs/{job_name}", name, job_header + job_body)
        job_names.append(name)
    df_jobs = generate_job_list(
        f"jobs/{job_name}", job_name, "int-demultiplexed-rna-map", job_names
    )
    generate_submit_file(
        f"submits/README-{job_name.upper()}", df_jobs["job_path"].tolist()
    )
    return df_jobs
