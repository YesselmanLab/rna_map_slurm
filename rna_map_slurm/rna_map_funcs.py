import shutil
import matplotlib.pyplot as plt
import pandas as pd

from rna_map.mutation_histogram import (
    get_dataframe,
    get_mut_histos_from_pickle_file,
    merge_mut_histo_dicts,
    write_mut_histos_to_pickle_file,
)

from rna_map_slurm.logger import get_logger
from rna_map_slurm.plotting import plot_pop_avg_from_row

log = get_logger("RNA-MAP-FUNCS")


def get_mut_histo_dataframe(mut_histos):
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
    df_results = get_dataframe(mut_histos, cols)
    df_results.rename(columns={"pop_avg": "data"}, inplace=True)
    return df_results


def generate_pop_avg_plots(
    df_results: pd.DataFrame, run_name: str, dir_name: str
) -> None:
    """
    Generates population average plots from the DataFrame rows, saves them,
    and copies them to specified directories.

    Args:
        df_results (pd.DataFrame): DataFrame containing the results.
        run_name (str): Name of the run for organizing the saved plots.
        dir_name (str): Directory name for organizing the saved plots.
    """
    i = 0
    for _, row in df_results.iterrows():
        log.info("plotting : " + row["name"])
        plot_pop_avg_from_row(row)
        final_path = f"results/{run_name}/processed/{dir_name}/output/BitVector_Files/"
        plt.title(
            "num_aligned: " + str(row["num_aligned"]) + "   sn: " + str(row["sn"])
        )
        plt.savefig(final_path + f"{row['name']}.png")
        plt.close()
        shutil.copy(
            final_path + f"{row['name']}.png",
            f"results/{run_name}/plots/pop_avg_pngs/{dir_name}_{row['name']}.png",
        )
        ax = plot_pop_avg_from_row(row)
        plt.title(
            "num_aligned: " + str(row["num_aligned"]) + "   sn: " + str(row["sn"])
        )
        ax.set_ylim(0, 0.10)
        plt.savefig(
            f"results/{run_name}/plots/pop_avg_pngs_0_10/{dir_name}_{row['name']}.png"
        )
        plt.close()
        ax = plot_pop_avg_from_row(row)
        plt.title(
            "num_aligned: " + str(row["num_aligned"]) + "   sn: " + str(row["sn"])
        )
        ax.set_ylim(0, 0.05)
        plt.savefig(
            f"results/{run_name}/plots/pop_avg_pngs_0_05/{dir_name}_{row['name']}.png"
        )
        plt.close()
        i += 1
        if i > 50:
            break
