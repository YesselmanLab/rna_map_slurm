"""DataFrame utility functions."""

from __future__ import annotations

from typing import Any

import pandas as pd

from rna_map_slurm.utils.logging import get_logger

log = get_logger("dataframe")


def get_data_row(
    df: pd.DataFrame,
    barcode_seq: str,
    rna_name: str,
) -> pd.Series[Any] | None:
    """Get a single row from DataFrame matching barcode and RNA name.

    Args:
        df: DataFrame containing the data.
        barcode_seq: Barcode sequence to match.
        rna_name: RNA/construct name to match.

    Returns:
        Matching row as a Series, or None if not found or multiple matches.
    """
    df_sub = df.query("barcode_seq == @barcode_seq and construct == @rna_name")

    if len(df_sub) == 0:
        log.error(f"barcode_seq {barcode_seq} with rna {rna_name} not found in data.csv")
        return None

    if len(df_sub) > 1:
        log.warning(
            f"barcode_seq {barcode_seq} with rna {rna_name} has multiple entries in data.csv"
        )
        return None

    result: pd.Series[Any] = df_sub.iloc[0]
    return result


def split_dataframe_into_n(df: pd.DataFrame, n: int) -> list[pd.DataFrame]:
    """Split a DataFrame into n roughly equal parts.

    Args:
        df: DataFrame to split.
        n: Number of parts to split into.

    Returns:
        List of DataFrames with sizes as close to equal as possible.
    """
    if n <= 0:
        return [df.copy()]

    avg = len(df) // n
    remainder = len(df) % n
    result: list[pd.DataFrame] = []
    idx = 0

    for i in range(n):
        size = avg + (1 if i < remainder else 0)
        result.append(df.iloc[idx : idx + size].copy())
        idx += size

    return result
