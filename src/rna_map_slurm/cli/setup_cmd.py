"""Setup command for rna-map-slurm CLI."""

from __future__ import annotations

import json
import os
import shutil
import sys
from typing import Any, cast

import click
import pandas as pd
from barcode_demultiplex.demultiplex import find_helix_barcodes
from ylab_gdrive import get_sequence_run_info_df, get_sequences_df

from rna_map_slurm.cli.utils import (
    get_seq_path,
    replace_spaces_warn,
    save_params_to_yaml,
)
from rna_map_slurm.config.parameters import get_default_parameters, get_parameters_from_file
from rna_map_slurm.io.fastq import get_paired_fastqs
from rna_map_slurm.jobs.fastq_jobs import (
    generate_demultiplexing_jobs,
    generate_join_fastq_files_jobs,
    generate_split_fastq_jobs,
)
from rna_map_slurm.jobs.int_demultiplex_jobs import (
    generate_int_demultiplex_jobs,
    generate_int_demultiplex_rna_map_combine_jobs,
    generate_int_demultiplex_rna_map_jobs,
)
from rna_map_slurm.jobs.rna_map_jobs import generate_rna_map_combine_jobs, generate_rna_map_jobs
from rna_map_slurm.utils.logging import get_logger, setup_logging

log = get_logger("cli.setup")


@click.command()
@click.argument("run_name")
def get_data_csv(run_name: str) -> None:
    """Fetch data CSV from Google Sheets for a run."""
    os.makedirs("logs", exist_ok=True)
    _remove_old_log("logs/get_data_csv.log")
    setup_logging(file_name="get_data_csv.log")

    df = get_sequence_run_info_df()
    df = df[df["run_name"] == run_name].reset_index(drop=True)
    df = _format_sequencing_run_info(df)
    df.to_csv("data.csv", index=False)


@click.command()
@click.argument("data_csv")
@click.argument("data_dirs", nargs=-1)
@click.option("--param-file", type=click.Path(exists=True), default=None)
@click.option(
    "--rna-map-params",
    type=click.Path(exists=True),
    default=None,
    help="Path to rna-map parameters YAML file. Copied to inputs/ for use by jobs.",
)
def setup(
    data_csv: str,
    data_dirs: tuple[str, ...],
    param_file: str | None,
    rna_map_params: str | None,
) -> None:
    """Set up workflow directories and generate job files."""
    os.makedirs("logs", exist_ok=True)
    _remove_old_log("logs/setup.log")
    setup_logging(file_name="logs/setup.log")

    df = pd.read_csv(data_csv)
    params = _load_params(param_file)
    _log_setup_info(params)

    seq_path = get_seq_path(params)
    _setup_data_csvs(df)

    sub_df = _filter_constructs(df)
    sub_df = _validate_and_filter_codes(sub_df, seq_path)

    _setup_directories(sub_df)
    _setup_input_files(sub_df, seq_path)
    _setup_rna_map_params(rna_map_params)

    all_pfqs = _collect_paired_fastqs(data_dirs)
    fastq_chunks = cast(int, params["fastq_chunks"])
    params["num_dirs"] = fastq_chunks * len(all_pfqs)
    params["rna_map_params_file"] = _get_rna_map_params_path(rna_map_params)

    _generate_all_jobs(df, params, all_pfqs)


def _remove_old_log(path: str) -> None:
    """Remove old log file if it exists."""
    if os.path.isfile(path):
        os.remove(path)


def _load_params(param_file: str | None) -> dict[str, object]:
    """Load parameters from file or use defaults."""
    if param_file is not None:
        log.info(f"Reading param file: {param_file}")
        return get_parameters_from_file(param_file)
    log.info("Using default parameters")
    return get_default_parameters()


def _log_setup_info(params: dict[str, object]) -> None:
    """Log setup information."""
    ctx = click.get_current_context()
    log.info("Ran at commandline as: %s", " ".join(sys.argv))
    log.info("Command line arguments and options:")
    for param, value in ctx.params.items():
        log.info(f"{param}: {value}")
    save_params_to_yaml(params, "logs/params.yaml")
    log.info("\n" + json.dumps(params, indent=4, default=str))
    log.info("saved params to logs/params.yaml")


def _format_sequencing_run_info(df: pd.DataFrame) -> pd.DataFrame:
    """Format sequencing run information from Google Sheets."""
    if len(df) == 0:
        log.error("No sequencing run information found")
        sys.exit(1)

    log.info("Formatting sequencing run information")
    df = replace_spaces_warn(df, "exp_name")
    df = replace_spaces_warn(df, "construct")
    df = replace_spaces_warn(df, "exp_type")

    df_seq = get_sequences_df()
    demult_cmds = _get_demultiplex_commands(df, df_seq)
    df["demult_cmd"] = demult_cmds

    return df


def _get_demultiplex_commands(df: pd.DataFrame, df_seq: pd.DataFrame) -> list[str | None]:
    """Get demultiplex commands for each row."""
    demult_cmds: list[str | None] = []

    for _, row in df.iterrows():
        seq = df_seq[df_seq["code"] == row["code"]]
        if len(seq) == 0:
            if not str(row["exp_name"]).lower().startswith("eich"):
                log.warning(f"No sequence information found for {row['code']}")
            demult_cmds.append(None)
            continue

        seq_row = seq.iloc[0]
        demult_cmds.append(seq_row.get("demultiplex"))
        if pd.notna(seq_row.get("demultiplex")):
            log.info(f"Found demultiplexing command for {row['construct']}: {seq_row['demultiplex']}")

    return demult_cmds


def _setup_data_csvs(df: pd.DataFrame) -> None:
    """Set up data CSV files."""
    os.makedirs("csvs", exist_ok=True)
    rm_df = df.query("exp_name.str.lower().str.startswith('eich')")
    rm_df.to_csv("csvs/data-eichhorn-constructs.csv", index=False)


def _filter_constructs(df: pd.DataFrame) -> pd.DataFrame:
    """Filter to Yesselman lab constructs."""
    sub_df = df.query("not exp_name.str.lower().str.startswith('sub')")
    sub_df.to_csv("csvs/data-yesselman-constructs.csv", index=False)
    return sub_df


def _validate_and_filter_codes(df: pd.DataFrame, seq_path: str) -> pd.DataFrame:
    """Validate codes have required files and filter invalid ones."""
    keep_indices = []
    for i, row in df.iterrows():
        csv_path = f"{seq_path}/rna/{row['code']}.csv"
        if not os.path.isfile(csv_path):
            log.warning(f"{row['code']} does not have a RNA CSV file")
            continue
        keep_indices.append(i)
    return df.loc[keep_indices]


def _setup_directories(df: pd.DataFrame) -> None:
    """Create directory structure for the workflow."""
    log.info("Setting up directories")
    log.info("Creating directories: jobs, submits, data, inputs, results")

    for dir_name in ["jobs", "submits", "data", "inputs", "results"]:
        os.makedirs(dir_name, exist_ok=True)

    for run_name in df["run_name"].unique():
        _create_run_directories(run_name)

    os.makedirs("inputs/barcode_jsons", exist_ok=True)
    os.makedirs("inputs/fastas", exist_ok=True)
    os.makedirs("inputs/rnas", exist_ok=True)


def _setup_rna_map_params(rna_map_params: str | None) -> None:
    """Copy rna-map parameters file to inputs directory if provided."""
    if rna_map_params is None:
        log.info("No custom rna-map params file provided, will use bundled defaults")
        return

    dest_path = "inputs/rna-map-params.yml"
    log.info(f"Copying rna-map params file to {dest_path}")
    shutil.copy(rna_map_params, dest_path)


def _get_rna_map_params_path(rna_map_params: str | None) -> str | None:
    """Get the path to the rna-map params file for job generation."""
    if rna_map_params is None:
        return None
    return os.path.abspath("inputs/rna-map-params.yml")


def _create_run_directories(run_name: str) -> None:
    """Create directories for a specific run."""
    base = f"results/{run_name}"
    for subdir in [
        "",
        "processed",
        "summary",
        "plots",
        "plots/pop_avg_pngs",
        "plots/pop_avg_pngs_0_10",
        "plots/pop_avg_pngs_0_05",
    ]:
        os.makedirs(f"{base}/{subdir}", exist_ok=True)


def _setup_input_files(df: pd.DataFrame, seq_path: str) -> None:
    """Copy input files and generate barcode JSONs."""
    for _, row in df.iterrows():
        _copy_sequence_files(row, seq_path)
        _generate_barcode_json(row, seq_path)


def _copy_sequence_files(row: pd.Series[Any], seq_path: str) -> None:
    """Copy FASTA and CSV files for a construct."""
    code = row["code"]
    df_seq = pd.read_csv(f"{seq_path}/rna/{code}.csv")
    shutil.copy(f"{seq_path}/fastas/{code}.fasta", "inputs/fastas/")
    df_seq.to_csv(f"inputs/rnas/{code}.csv", index=False)


def _generate_barcode_json(row: pd.Series[Any], seq_path: str) -> None:
    """Generate barcode JSON for internal demultiplexing."""
    if pd.isna(row["demult_cmd"]):
        return

    helices = _parse_helix_arguments(str(row["demult_cmd"]))
    if not helices:
        log.error(f"Trying to setup demultiplexing on {row['code']} but no helices supplied")
        return

    df_seq = pd.read_csv(f"{seq_path}/rna/{row['code']}.csv")
    df_barcodes = find_helix_barcodes(df_seq, helices)
    df_barcodes.to_json(f"inputs/barcode_jsons/{row['code']}.json", orient="records")


def _parse_helix_arguments(demult_cmd: str) -> list[list[int]]:
    """Parse helix arguments from demultiplex command."""
    args = demult_cmd.split()
    helices: list[list[int]] = []

    for i, arg in enumerate(args):
        if arg in ("--helix", "-helix"):
            helices.append([int(args[i + 1]), int(args[i + 2]), int(args[i + 3])])

    return helices


def _collect_paired_fastqs(data_dirs: tuple[str, ...]) -> list[Any]:
    """Collect all paired FASTQ files from data directories."""
    all_pfqs: list[Any] = []
    for d in data_dirs:
        d = os.path.abspath(d)
        all_pfqs.extend(get_paired_fastqs(d))
    return all_pfqs


def _generate_all_jobs(
    df: pd.DataFrame,
    params: dict[str, object],
    all_pfqs: list[Any],
) -> None:
    """Generate all job files for the workflow."""
    single_df = df.query("demult_cmd.isnull()")
    single_df.to_csv("csvs/data-single.csv", index=False)

    int_mult_df = df.query("not demult_cmd.isnull()")
    params_any = cast(dict[str, Any], params)

    df_jobs_list = [
        generate_split_fastq_jobs(all_pfqs, params_any),
        generate_demultiplexing_jobs(params_any),
        generate_join_fastq_files_jobs(params_any),
        generate_rna_map_jobs(params_any, single_df),
        generate_rna_map_combine_jobs(params_any, single_df),
    ]

    if len(int_mult_df) > 0:
        int_mult_df.to_csv("csvs/data-int_multiplex.csv", index=False)
        df_jobs_list.extend([
            generate_int_demultiplex_jobs(params_any, int_mult_df),
            generate_int_demultiplex_rna_map_jobs(params_any, int_mult_df),
            generate_int_demultiplex_rna_map_combine_jobs(params_any, int_mult_df),
        ])

    df_job = pd.concat(df_jobs_list, ignore_index=True)
    df_job.to_csv("jobs.csv", index=False)
    log.info(f"Generated {len(df_job)} jobs")
