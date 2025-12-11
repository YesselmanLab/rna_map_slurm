"""Shared CLI utilities."""

from __future__ import annotations

import os
import sys

import pandas as pd
import yaml

from rna_map_slurm.utils.logging import get_logger

log = get_logger("cli.utils")


def replace_spaces_warn(df: pd.DataFrame, column_name: str) -> pd.DataFrame:
    """Replace spaces with underscores in a DataFrame column.

    Args:
        df: DataFrame to modify.
        column_name: Column to process.

    Returns:
        Modified DataFrame.
    """
    for index, row in df.iterrows():
        value = row[column_name]
        if " " in str(value):
            log.warning(
                f"Replacing spaces with underscores in row {index} "
                f"for column '{column_name}'. Original value: '{value}'"
            )
            df.loc[index, column_name] = str(value).replace(" ", "_")  # type: ignore[index]
    return df


def get_seq_path(params: dict[str, object]) -> str:
    """Get path to sequence files from environment or params.

    Args:
        params: Parameters dictionary.

    Returns:
        Path to sequence files.
    """
    seq_path = os.environ.get("SEQPATH", "")
    if seq_path:
        log.info(f"Setting seq_path from environment variable: {seq_path}")
        return seq_path

    paths = params.get("paths", {})
    if isinstance(paths, dict):
        seq_path = str(paths.get("seq_path", ""))
        if seq_path:
            log.info(f"Setting seq_path from params file: {seq_path}")
            return seq_path

    log.error("SEQPATH not set in environment variable or params file")
    sys.exit(1)


def save_params_to_yaml(params: dict[str, object], path: str) -> None:
    """Save parameters to a YAML file.

    Args:
        params: Parameters to save.
        path: Output file path.
    """
    with open(path, "w", encoding="utf-8") as f:
        yaml.dump(params, f)


def submit_jobs(df: pd.DataFrame) -> None:
    """Submit jobs to SLURM scheduler.

    Args:
        df: DataFrame with job information.
    """
    log.info("Submitting jobs")
    for _, row in df.iterrows():
        os.system(f"sbatch {row['job_path']}")
