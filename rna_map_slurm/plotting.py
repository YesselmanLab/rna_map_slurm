import matplotlib.pyplot as plt
import numpy as np


def colors_for_sequence(seq: str):
    colors = []
    for e in seq:
        if e == "A":
            colors.append("red")
        elif e == "C":
            colors.append("blue")
        elif e == "G":
            colors.append("orange")
        else:
            colors.append("green")
    return colors


# sub plots ###################################################################


def plot_pop_avg(seq, ss, reactivities, ax=None):
    colors = colors_for_sequence(seq)
    x = list(range(len(seq)))
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(20, 4))
    ax.bar(range(0, len(reactivities)), reactivities, color=colors)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}\n{nt}" for s, nt in zip(seq, ss)])
    return ax


def plot_pop_avg_from_row(row, data_col="data", ax=None):
    return plot_pop_avg(row["sequence"], row["structure"], row[data_col], ax)


# full plots ###################################################################


def plot_pop_avg_diff_from_rows(row1, row2, data_col="data", **kwargs):
    fig, axes = plt.subplots(3, 1, **kwargs)
    plot_pop_avg_from_row(row1, data_col, axes[0])
    plot_pop_avg_from_row(row2, data_col, axes[1])
    diff = {
        "sequence": row1["sequence"],
        "structure": row1["structure"],
        data_col: np.array(row1[data_col]) - np.array(row2[data_col]),
    }
    plot_pop_avg_from_row(diff, data_col, axes[2])
    return fig


def plot_pop_avg_all(df, data_col="data", **kwargs):
    fig, axes = plt.subplots(len(df), 1, **kwargs)
    j = 0
    for i, row in df.iterrows():
        colors = colors_for_sequence(row["sequence"])
        axes[j].bar(range(0, len(row[data_col])), row[data_col], color=colors)
        axes[j].set_title(row["rna_name"])
        j += 1
    plot_pop_avg_from_row(df.iloc[-1], ax=axes[-1])
    return fig


def plot_pop_avg_traces_all(df, **kwargs):
    fig, ax = plt.subplots(1, 1, **kwargs)
    for i, row in df.iterrows():
        plt.plot(row["data"], label=row["rna_name"])
    fig.legend(loc="upper left")
    return fig
