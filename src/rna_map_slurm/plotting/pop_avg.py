"""Population average plotting utilities."""

from __future__ import annotations

import shutil
from typing import TYPE_CHECKING, Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rna_map_slurm.utils.logging import get_logger

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

log = get_logger("plotting")

MAX_PLOTS = 50


def colors_for_sequence(seq: str) -> list[str]:
    """Get colors for each nucleotide in a sequence.

    Args:
        seq: RNA/DNA sequence.

    Returns:
        List of color strings for each nucleotide.
    """
    color_map = {
        "A": "red",
        "C": "blue",
        "G": "orange",
        "T": "green",
        "U": "green",
    }
    return [color_map.get(nt, "gray") for nt in seq]


def plot_pop_avg(
    seq: str,
    ss: str,
    reactivities: list[float],
    ax: Axes | None = None,
) -> Axes:
    """Plot population average reactivity data.

    Args:
        seq: RNA sequence.
        ss: Secondary structure in dot-bracket notation.
        reactivities: Reactivity values for each position.
        ax: Optional matplotlib axes to plot on.

    Returns:
        Matplotlib axes object.
    """
    colors = colors_for_sequence(seq)
    x = list(range(len(seq)))

    if ax is None:
        _, ax = plt.subplots(1, figsize=(20, 4))

    ax.bar(range(len(reactivities)), reactivities, color=colors)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}\n{nt}" for s, nt in zip(seq, ss, strict=True)])

    return ax


def plot_pop_avg_from_row(
    row: pd.Series[Any],
    data_col: str = "data",
    ax: Axes | None = None,
) -> Axes:
    """Plot population average from a DataFrame row.

    Args:
        row: DataFrame row with sequence, structure, and data.
        data_col: Column name for reactivity data.
        ax: Optional matplotlib axes.

    Returns:
        Matplotlib axes object.
    """
    return plot_pop_avg(row["sequence"], row["structure"], row[data_col], ax)


def plot_pop_avg_diff_from_rows(
    row1: pd.Series[Any],
    row2: pd.Series[Any],
    data_col: str = "data",
    **kwargs: Any,
) -> Figure:
    """Plot comparison of two population averages with difference.

    Args:
        row1: First DataFrame row.
        row2: Second DataFrame row.
        data_col: Column name for reactivity data.
        **kwargs: Additional arguments for plt.subplots.

    Returns:
        Matplotlib figure object.
    """
    fig, axes = plt.subplots(3, 1, **kwargs)
    plot_pop_avg_from_row(row1, data_col, axes[0])
    plot_pop_avg_from_row(row2, data_col, axes[1])

    diff = {
        "sequence": row1["sequence"],
        "structure": row1["structure"],
        data_col: np.array(row1[data_col]) - np.array(row2[data_col]),
    }
    plot_pop_avg_from_row(pd.Series(diff), data_col, axes[2])

    result: Figure = fig
    return result


def plot_pop_avg_all(
    df: pd.DataFrame,
    data_col: str = "data",
    **kwargs: Any,
) -> Figure:
    """Plot all population averages in a DataFrame.

    Args:
        df: DataFrame with sequence data.
        data_col: Column name for reactivity data.
        **kwargs: Additional arguments for plt.subplots.

    Returns:
        Matplotlib figure object.
    """
    fig, axes = plt.subplots(len(df), 1, **kwargs)

    for j, (_, row) in enumerate(df.iterrows()):
        colors = colors_for_sequence(row["sequence"])
        axes[j].bar(range(len(row[data_col])), row[data_col], color=colors)
        axes[j].set_title(row.get("rna_name", row.get("name", f"Sample {j}")))

    result: Figure = fig
    return result


def generate_pop_avg_plots(
    df_results: pd.DataFrame,
    run_name: str,
    dir_name: str,
) -> None:
    """Generate and save population average plots at multiple scales.

    Args:
        df_results: DataFrame with reactivity results.
        run_name: Run name for path organization.
        dir_name: Directory name for the construct.
    """
    final_path = f"results/{run_name}/processed/{dir_name}/output/BitVector_Files/"

    for i, (_, row) in enumerate(df_results.iterrows()):
        if i >= MAX_PLOTS:
            break

        log.info(f"plotting: {row['name']}")
        _save_plot_at_scale(row, run_name, dir_name, final_path, None)
        _save_plot_at_scale(row, run_name, dir_name, final_path, 0.10)
        _save_plot_at_scale(row, run_name, dir_name, final_path, 0.05)


def _save_plot_at_scale(
    row: pd.Series[Any],
    run_name: str,
    dir_name: str,
    final_path: str,
    y_max: float | None,
) -> None:
    """Save a plot at a specific y-axis scale.

    Args:
        row: DataFrame row with data.
        run_name: Run name.
        dir_name: Directory name.
        final_path: Base output path.
        y_max: Maximum y-axis value, or None for auto-scale.
    """
    ax = plot_pop_avg_from_row(row)
    title = f"num_aligned: {row['num_aligned']}   sn: {row['sn']}"
    plt.title(title)

    if y_max is not None:
        ax.set_ylim(0, y_max)
        scale_suffix = f"_0_{int(y_max * 100):02d}"
        dest_dir = f"results/{run_name}/plots/pop_avg_pngs{scale_suffix}"
    else:
        dest_dir = f"results/{run_name}/plots/pop_avg_pngs"
        # Save to final_path for full-scale version
        plt.savefig(f"{final_path}{row['name']}.png")
        shutil.copy(
            f"{final_path}{row['name']}.png",
            f"{dest_dir}/{dir_name}_{row['name']}.png",
        )
        plt.close()
        return

    plt.savefig(f"{dest_dir}/{dir_name}_{row['name']}.png")
    plt.close()
