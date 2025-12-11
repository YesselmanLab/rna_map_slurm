"""Helper functions for RNA-map result processing."""

from __future__ import annotations

from typing import Any

import pandas as pd
from rna_map.mutation_histogram import get_dataframe


def get_mut_histo_dataframe(mut_histos: dict[str, Any]) -> pd.DataFrame:
    """Convert mutation histogram dictionary to DataFrame.

    Args:
        mut_histos: Dictionary of mutation histogram data.

    Returns:
        DataFrame with processed mutation histogram data.
    """
    cols = [
        "name",
        "sequence",
        "structure",
        "pop_avg",
        "sn",
        "num_reads",
        "num_aligned",
        "no_mut",
        "1_mut",
        "2_mut",
        "3_mut",
        "3plus_mut",
    ]
    df_results: pd.DataFrame = get_dataframe(mut_histos, cols)
    df_results.rename(columns={"pop_avg": "data"}, inplace=True)
    return df_results
