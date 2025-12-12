"""Deposit results command for rna-map-slurm CLI."""

from __future__ import annotations

import os
import re
import sys

import click
import pandas as pd

from rna_map_slurm.utils.logging import get_logger, setup_logging

log = get_logger("cli.deposit")


@click.command()
@click.option("-v", "--version", default="v1", type=str, help="Version string (e.g., v1)")
@click.option("-p", "--path", default=None, help="Base path for depositing results")
@click.option("--overwrite", is_flag=True, help="Overwrite existing version")
def deposit_results(version: str, path: str | None, overwrite: bool) -> None:
    """Deposit analysis results to storage location."""
    setup_logging()

    if not version.startswith("v"):
        log.error("version needs to start with 'v'")
        sys.exit(1)

    if path is None:
        path = os.environ.get("NRDSTORSHARED", "")

    data_dirs = _parse_data_dirs_from_log()
    if data_dirs is None:
        return

    df = pd.read_csv("data.csv")
    _deposit_all_runs(df, path, version, data_dirs, overwrite)


def _parse_data_dirs_from_log() -> list[str] | None:
    """Parse data directories from setup log."""
    log_file_path = "logs/setup.log"

    if not os.path.isfile(log_file_path):
        log.error(f"Log file not found: {log_file_path}")
        return None

    with open(log_file_path, encoding="utf-8") as log_file:
        for line in log_file:
            match = re.match(r"^rna-map-slurm\.cli\.setup - INFO - data_dirs: \((.*)\)", line)
            if match:
                data_dirs = match.group(1).split(",")[:1]
                return [x.strip().strip("'\"") for x in data_dirs]

    log.error("Could not find data_dirs in log file")
    return None


def _deposit_all_runs(
    df: pd.DataFrame,
    path: str,
    version: str,
    data_dirs: list[str],
    overwrite: bool,
) -> None:
    """Deposit results for all runs."""
    for run_name in df["run_name"].unique():
        _deposit_single_run(df, path, run_name, version, data_dirs, overwrite)


def _deposit_single_run(
    df: pd.DataFrame,
    path: str,
    run_name: str,
    version: str,
    data_dirs: list[str],
    overwrite: bool,
) -> None:
    """Deposit results for a single run."""
    run_save_path = f"{path}/{run_name}"

    _ensure_run_directory(run_save_path)
    _save_run_data(df, run_name, run_save_path)
    _copy_logs_and_csvs(run_save_path)
    _copy_raw_data(run_save_path, data_dirs)
    _copy_demultiplexed(run_save_path)
    _copy_analysis(run_save_path, run_name, version, overwrite)


def _ensure_run_directory(run_save_path: str) -> None:
    """Ensure run directory exists."""
    if not os.path.isdir(run_save_path):
        log.info(f"{run_save_path} does not exist, creating")
        os.makedirs(run_save_path)


def _save_run_data(df: pd.DataFrame, run_name: str, run_save_path: str) -> None:
    """Save run-specific data CSV."""
    df_sub = df.query("run_name == @run_name")
    df_sub.to_csv(f"{run_save_path}/data.csv", index=False)
    log.info(f"copying data.csv to {run_save_path}")


def _copy_logs_and_csvs(run_save_path: str) -> None:
    """Copy logs and CSVs to run directory."""
    os.system(f"cp -r logs {run_save_path}")
    os.system(f"cp -r csvs {run_save_path}")
    log.info(f"copying logs and csvs to {run_save_path}")


def _copy_raw_data(run_save_path: str, data_dirs: list[str]) -> None:
    """Copy raw data if not already present."""
    raw_path = f"{run_save_path}/raw"

    if os.path.isdir(raw_path):
        return

    if len(data_dirs) != 1:
        log.error("Multiple data dirs not yet supported for raw data copy")
        return

    log.info(f"copying {data_dirs[0]} to {raw_path}")
    os.system(f"cp -r {data_dirs[0]} {raw_path}")


def _copy_demultiplexed(run_save_path: str) -> None:
    """Copy demultiplexed data if not already present."""
    demultiplex_path = f"{run_save_path}/demultiplexed"

    if os.path.isdir(demultiplex_path):
        return

    log.info(f"copying demultiplexed to {demultiplex_path}")
    os.system(f"cp -r demultiplexed {run_save_path}")


def _copy_analysis(
    run_save_path: str,
    run_name: str,
    version: str,
    overwrite: bool,
) -> None:
    """Copy analysis results to versioned directory."""
    analysis_path = f"{run_save_path}/analysis"
    if not os.path.isdir(analysis_path):
        log.info(f"{analysis_path} does not exist, creating")
        os.makedirs(analysis_path)

    v_path = f"{analysis_path}/{version}"
    if os.path.isdir(v_path) and not overwrite:
        log.error(f"{v_path} exists. Use another version or --overwrite")
        return

    log.info(f"copying results/{run_name} to {v_path}")
    os.system(f"cp -r results/{run_name} {v_path}")
