"""Internal demultiplexing task implementations."""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
from typing import Any

import pandas as pd
import rna_map
from rna_map.mutation_histogram import (
    get_mut_histos_from_pickle_file,
    merge_mut_histo_dicts,
    write_mut_histos_to_pickle_file,
)
from rna_map.parameters import parse_parameters_from_file
from seq_tools.dataframe import to_dna, to_fasta
from seq_tools.sequence import get_reverse_complement

from rna_map_slurm.config.paths import get_rna_map_defaults_path
from rna_map_slurm.plotting.pop_avg import generate_pop_avg_plots
from rna_map_slurm.tasks.rna_map_helpers import get_mut_histo_dataframe
from rna_map_slurm.utils.files import get_file_size, random_string
from rna_map_slurm.utils.logging import get_logger

log = get_logger("tasks.int_demultiplex")


def get_tmp_dir(params: dict[str, Any] | None = None) -> str:
    """Get temporary directory path from params or environment.

    Args:
        params: Optional workflow parameters.

    Returns:
        Temporary directory path.
    """
    if params is not None:
        paths = params.get("paths", {})
        if isinstance(paths, dict):
            return str(paths.get("tmp", "/scratch"))
        return "/scratch"
    return os.environ.get("SCRATCH", "/scratch")


def int_demultiplex(
    construct_barcode: str,
    b1_seq: str,
    b2_seq: str,
    b1_min_pos: int,
    b1_max_pos: int,
    b2_min_pos: int,
    b2_max_pos: int,
    params: dict[str, Any] | None = None,
) -> None:
    """Perform internal demultiplexing using seqkit.

    Args:
        construct_barcode: Barcode sequence for the construct.
        b1_seq: First barcode sequence.
        b2_seq: Second barcode sequence.
        b1_min_pos: Minimum position for barcode 1.
        b1_max_pos: Maximum position for barcode 1.
        b2_min_pos: Minimum position for barcode 2.
        b2_max_pos: Maximum position for barcode 2.
        params: Optional workflow parameters.
    """
    tmp_base = get_tmp_dir(params)
    tmp_dir = f"{tmp_base}/{random_string(10)}"
    log.info(f"tmp_dir: {tmp_dir}")
    os.makedirs(tmp_dir, exist_ok=True)

    try:
        _run_seqkit_filtering(construct_barcode, b1_seq, b2_seq, b1_min_pos, b1_max_pos, b2_min_pos, b2_max_pos, tmp_dir)
        _find_common_reads(tmp_dir)
        _extract_common_reads(construct_barcode, b1_seq, b2_seq, tmp_dir)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _run_seqkit_filtering(
    construct_barcode: str,
    b1_seq: str,
    b2_seq: str,
    b1_min_pos: int,
    b1_max_pos: int,
    b2_min_pos: int,
    b2_max_pos: int,
    tmp_dir: str,
) -> None:
    """Run seqkit grep to filter reads by barcode sequences.

    Args:
        construct_barcode: Construct barcode.
        b1_seq: First barcode sequence.
        b2_seq: Second barcode sequence.
        b1_min_pos: Min position for barcode 1.
        b1_max_pos: Max position for barcode 1.
        b2_min_pos: Min position for barcode 2.
        b2_max_pos: Max position for barcode 2.
        tmp_dir: Temporary directory.
    """
    r2_path = f"demultiplexed/{construct_barcode}/test_R2.fastq.gz"
    r1_path = f"demultiplexed/{construct_barcode}/test_R1.fastq.gz"
    # Convert U to T for barcode search (FASTQ uses DNA, not RNA)
    b1_seq_dna = b1_seq.replace("U", "T")
    b2_seq_dna = b2_seq.replace("U", "T")
    b2_seq_rc = get_reverse_complement(b2_seq_dna)

    _run_command(
        f'seqkit grep -s -p "{b1_seq_dna}" -P -R {b1_min_pos - 2}:{b1_max_pos + 2} '
        f'{r2_path} -o {tmp_dir}/test_R2.fastq.gz'
    )
    _run_command(
        f'seqkit grep -s -p "{b2_seq_rc}" -P -R {b2_min_pos - 2}:{b2_max_pos + 2} '
        f'{r1_path} -o {tmp_dir}/test_R1.fastq.gz'
    )


def _find_common_reads(tmp_dir: str) -> None:
    """Find read names common to both R1 and R2 filtered files.

    Args:
        tmp_dir: Temporary directory with filtered files.
    """
    _run_command(f"seqkit seq -n {tmp_dir}/test_R2.fastq.gz > {tmp_dir}/R2_names.txt")
    _run_command(f"seqkit seq -n {tmp_dir}/test_R1.fastq.gz > {tmp_dir}/R1_names.txt")
    _run_command(f"sort {tmp_dir}/R1_names.txt > {tmp_dir}/R1_names_sorted.txt")
    _run_command(f"sort {tmp_dir}/R2_names.txt > {tmp_dir}/R2_names_sorted.txt")

    awk_cmd = (
        f"awk 'NR==FNR{{a[$1]; next}} $1 in a' "
        f"{tmp_dir}/R1_names_sorted.txt {tmp_dir}/R2_names_sorted.txt "
        f"| awk '{{print($1)}}' > {tmp_dir}/common_names.txt"
    )
    _run_command(awk_cmd)


def _extract_common_reads(
    construct_barcode: str,
    b1_seq: str,
    b2_seq: str,
    tmp_dir: str,
) -> None:
    """Extract reads with common names to final output.

    Args:
        construct_barcode: Construct barcode.
        b1_seq: First barcode sequence.
        b2_seq: Second barcode sequence.
        tmp_dir: Temporary directory.
    """
    output_dir = f"int-demultiplexed/{construct_barcode}"
    os.makedirs(output_dir, exist_ok=True)

    _run_command(
        f"seqkit grep -f {tmp_dir}/common_names.txt {tmp_dir}/test_R2.fastq.gz "
        f"-o {output_dir}/{b1_seq}_{b2_seq}_mate1.fastq.gz"
    )
    _run_command(
        f"seqkit grep -f {tmp_dir}/common_names.txt {tmp_dir}/test_R1.fastq.gz "
        f"-o {output_dir}/{b1_seq}_{b2_seq}_mate2.fastq.gz"
    )


def _run_command(cmd: str) -> None:
    """Run a shell command.

    Args:
        cmd: Command to run.
    """
    subprocess.run(cmd, shell=True, check=True)


def int_demultiplex_rna_map(
    code: str,
    lib_barcode_seq: str,
    construct_barcode_seq: str,
    params: dict[str, Any] | None = None,
    params_file: str | None = None,
) -> None:
    """Run RNA mapping on internally demultiplexed reads.

    Args:
        code: Construct code.
        lib_barcode_seq: Library barcode sequence.
        construct_barcode_seq: Construct barcode sequence.
        params: Optional workflow parameters.
        params_file: Optional path to rna-map parameters file.
    """
    fastq_paths = _find_demultiplexed_fastqs(lib_barcode_seq, construct_barcode_seq)
    if fastq_paths is None:
        return

    mate_1_path, mate_2_path = fastq_paths
    if _files_too_small(mate_1_path, mate_2_path, construct_barcode_seq):
        return

    _run_int_demultiplex_rna_map(
        code, lib_barcode_seq, construct_barcode_seq, mate_1_path, mate_2_path, params, params_file
    )


def _find_demultiplexed_fastqs(
    lib_barcode_seq: str,
    construct_barcode_seq: str,
) -> tuple[str, str] | None:
    """Find demultiplexed FASTQ files.

    Args:
        lib_barcode_seq: Library barcode.
        construct_barcode_seq: Construct barcode.

    Returns:
        Tuple of (mate1_path, mate2_path) or None if not found.
    """
    fastq_dir = f"int-demultiplexed/{lib_barcode_seq}"
    log.info(f"Looking in {fastq_dir}")

    mate1_files = glob.glob(f"{fastq_dir}/{construct_barcode_seq}_mate1.fastq.gz")
    mate2_files = glob.glob(f"{fastq_dir}/{construct_barcode_seq}_mate2.fastq.gz")

    if not mate1_files or not mate2_files:
        log.error(f"Could not find FASTQ files for {construct_barcode_seq}")
        return None

    return os.path.abspath(mate1_files[0]), os.path.abspath(mate2_files[0])


def _files_too_small(
    mate_1_path: str,
    mate_2_path: str,
    construct_barcode_seq: str,
    min_size: int = 100,
) -> bool:
    """Check if files are too small to process.

    Args:
        mate_1_path: Path to mate 1 file.
        mate_2_path: Path to mate 2 file.
        construct_barcode_seq: Construct barcode for logging.
        min_size: Minimum file size in bytes.

    Returns:
        True if files are too small.
    """
    if get_file_size(mate_1_path) < min_size or get_file_size(mate_2_path) < min_size:
        log.warning(f"Skipping {construct_barcode_seq} because file size is too small")
        return True
    return False


def _run_int_demultiplex_rna_map(
    code: str,
    lib_barcode_seq: str,
    construct_barcode_seq: str,
    mate_1_path: str,
    mate_2_path: str,
    params: dict[str, Any] | None = None,
    params_file: str | None = None,
) -> None:
    """Execute RNA mapping for internally demultiplexed reads.

    Args:
        code: Construct code.
        lib_barcode_seq: Library barcode.
        construct_barcode_seq: Construct barcode.
        mate_1_path: Path to mate 1 FASTQ.
        mate_2_path: Path to mate 2 FASTQ.
        params: Optional workflow parameters.
        params_file: Optional path to rna-map parameters file.
    """
    df_barcode = _prepare_barcode_dataframe(code, lib_barcode_seq, construct_barcode_seq)
    if df_barcode is None:
        return

    tmp_base = get_tmp_dir(params)
    tmp_dir = f"{tmp_base}/{random_string(10)}"
    os.makedirs(tmp_dir, exist_ok=True)
    cur_dir = os.path.abspath(os.getcwd())

    try:
        to_fasta(to_dna(df_barcode), f"{tmp_dir}/test.fasta")
        df_barcode[["name", "sequence", "structure"]].to_csv(f"{tmp_dir}/test.csv", index=False)

        os.chdir(tmp_dir)
        rna_map_params = _load_rna_map_params(params_file)
        rna_map.run.run("test.fasta", mate_1_path, mate_2_path, "test.csv", rna_map_params)

        output_path = f"{cur_dir}/int-demultiplexed-rna-map/{lib_barcode_seq}"
        shutil.move(
            "output/BitVector_Files/mutation_histos.p",
            f"{output_path}/mutation_histos_{construct_barcode_seq}.p",
        )
    except Exception as e:
        log.error(f"rna-map failed for {construct_barcode_seq}: {e}")
    finally:
        os.chdir(cur_dir)
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _load_rna_map_params(params_file: str | None = None) -> dict[str, Any]:
    """Load rna-map parameters from file or defaults.

    Args:
        params_file: Optional path to rna-map parameters file.

    Returns:
        Dictionary of rna-map parameters.
    """
    if params_file is not None:
        log.info(f"Using provided params file: {params_file}")
        return parse_parameters_from_file(params_file)

    default_params_file = get_rna_map_defaults_path()
    log.info(f"Using bundled rna-map defaults: {default_params_file}")
    return parse_parameters_from_file(default_params_file)


def _prepare_barcode_dataframe(
    code: str,
    lib_barcode_seq: str,
    construct_barcode_seq: str,
) -> pd.DataFrame | None:
    """Prepare barcode DataFrame for RNA mapping.

    Args:
        code: Construct code.
        lib_barcode_seq: Library barcode.
        construct_barcode_seq: Construct barcode.

    Returns:
        Filtered barcode DataFrame or None if not found.
    """
    df = pd.read_csv("data.csv")
    df_sub = df.query("barcode_seq == @lib_barcode_seq and code == @code")

    if len(df_sub) == 0:
        log.error(f"No barcode_seq {lib_barcode_seq} with code {code} found")
        return None

    row = df_sub.iloc[0]
    df_barcode = pd.read_json(f"inputs/barcode_jsons/{row['code']}.json")
    return df_barcode[df_barcode["full_barcode"] == construct_barcode_seq]


def int_demultiplex_rna_map_combine(
    barcode_seq: str,
    rna_name: str,
) -> None:
    """Combine RNA mapping results from internal demultiplexing.

    Args:
        barcode_seq: Barcode sequence.
        rna_name: RNA/construct name.
    """
    from rna_map_slurm.utils.dataframe import get_data_row

    df = pd.read_csv("data.csv")
    row = get_data_row(df, barcode_seq, rna_name)
    if row is None:
        return

    final_path = _setup_int_combine_paths(row)
    merged_mut_histos = _merge_int_mutation_histograms(barcode_seq)

    if not merged_mut_histos:
        log.warning("No mutation histogram files found to merge")
        return

    _save_int_combined_results(row, final_path, merged_mut_histos, rna_name)


def _setup_int_combine_paths(row: pd.Series[Any]) -> str:
    """Set up paths for internal demultiplexing combine results.

    Args:
        row: DataFrame row with construct info.

    Returns:
        Final output path.
    """
    run_path = f"results/{row['run_name']}"
    dir_name = f"{row['construct']}_{row['code']}_{row['data_type']}"
    final_path = f"{run_path}/processed/{dir_name}/output/BitVector_Files/"

    log.info(f"results path: {final_path}")
    os.makedirs(final_path, exist_ok=True)

    return final_path


def _merge_int_mutation_histograms(barcode_seq: str) -> dict[str, Any]:
    """Merge mutation histograms from internal demultiplexing.

    Args:
        barcode_seq: Barcode sequence.

    Returns:
        Merged mutation histogram dictionary.
    """
    mut_histo_files = glob.glob(f"int-demultiplexed-rna-map/{barcode_seq}/*p")
    log.info(f"found {len(mut_histo_files)} files")

    merged_mut_histos: dict[str, Any] = {}
    for i, mut_hist_file in enumerate(mut_histo_files):
        if i % 100 == 0:
            log.info(f"merged {i} mut histos")
        merge_mut_histo_dicts(merged_mut_histos, get_mut_histos_from_pickle_file(mut_hist_file))

    return merged_mut_histos


def _save_int_combined_results(
    row: pd.Series[Any],
    final_path: str,
    merged_mut_histos: dict[str, Any],
    rna_name: str,
) -> None:
    """Save combined internal demultiplexing results.

    Args:
        row: DataFrame row with metadata.
        final_path: Output directory.
        merged_mut_histos: Merged histogram data.
        rna_name: RNA name.
    """
    df_results = get_mut_histo_dataframe(merged_mut_histos)
    df_results["rna_name"] = rna_name

    cols = list(row.keys())
    for col in ["demult_cmd", "length"]:
        if col in cols:
            cols.remove(col)
    for col in cols:
        df_results[col] = row[col]

    df_results.to_json(f"{final_path}mutation_histos.json", orient="records")
    df_results = df_results.sort_values("num_aligned", ascending=False)

    dir_name = f"{row['construct']}_{row['code']}_{row['data_type']}"
    generate_pop_avg_plots(df_results, row["run_name"], dir_name)
    write_mut_histos_to_pickle_file(merged_mut_histos, f"{final_path}mutation_histos.p")
