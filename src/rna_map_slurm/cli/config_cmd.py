"""CLI command for generating example configuration files."""

from __future__ import annotations

from pathlib import Path

import click
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.scalarstring import DoubleQuotedScalarString as DQS


def _build_commented_config() -> CommentedMap:
    """Build a commented configuration map with all default values and descriptions."""
    config = CommentedMap()

    # Header comment
    config.yaml_set_start_comment(
        "RNA Map SLURM Configuration File\n"
        "Generated example with all available options and their descriptions.\n"
        "Copy this file and modify values as needed for your workflow.\n"
    )

    # fastq_chunks
    config["fastq_chunks"] = 100
    config.yaml_set_comment_before_after_key(
        "fastq_chunks",
        before=(
            "Number of chunks to split each FASTQ file into for parallel processing.\n"
            "Higher values = more parallelism but more overhead. Typical range: 50-200."
        ),
    )

    # paths section
    paths = CommentedMap()
    paths["log"] = "logs"
    paths["jobs"] = "jobs"
    paths["submits"] = "submits"
    paths["inputs"] = "inputs"
    paths["tmp"] = "/scratch"
    paths["seq_path"] = ""

    paths.yaml_set_comment_before_after_key(
        "log", before="Directory where log files will be written"
    )
    paths.yaml_set_comment_before_after_key(
        "jobs", before="Directory where SLURM job scripts are generated"
    )
    paths.yaml_set_comment_before_after_key(
        "submits", before="Directory where submit reference files are stored"
    )
    paths.yaml_set_comment_before_after_key(
        "inputs", before="Directory for input data files"
    )
    paths.yaml_set_comment_before_after_key(
        "tmp",
        before=(
            "Scratch space for temporary files. Use cluster-specific path.\n"
            "Common values: /scratch, /tmp, $TMPDIR"
        ),
    )
    paths.yaml_set_comment_before_after_key(
        "seq_path",
        before=(
            "Base path to sequence reference data.\n"
            "Can also be set via SEQPATH environment variable."
        ),
    )

    config["paths"] = paths
    config.yaml_set_comment_before_after_key(
        "paths", before="\nPath configuration for workflow directories"
    )

    # construct_options section
    construct_options = CommentedMap()
    construct_options["large_construct_threshold"] = 100
    construct_options.yaml_set_comment_before_after_key(
        "large_construct_threshold",
        before=(
            "Constructs with more sequences than this threshold get dedicated SLURM jobs\n"
            "instead of being batched together. Prevents single large constructs from\n"
            "causing timeouts for batched jobs."
        ),
    )

    config["construct_options"] = construct_options
    config.yaml_set_comment_before_after_key(
        "construct_options", before="\nOptions for handling different construct sizes"
    )

    # tasks_per_job section
    tasks_per_job = CommentedMap()
    tasks_per_job["default"] = 1
    tasks_per_job["split-fastq"] = 1
    tasks_per_job["demultiplex"] = 10
    tasks_per_job["trim-galore"] = 10
    tasks_per_job["join-fastq-files"] = 1
    tasks_per_job["rna-map"] = 25
    tasks_per_job["int-demultiplex"] = 10
    tasks_per_job["int-demultiplex-rna-map"] = 25
    tasks_per_job["int-demultiplex-rna-map-combine"] = 1

    tasks_per_job.yaml_set_comment_before_after_key(
        "default", before="Default number of tasks per job if not specified"
    )
    tasks_per_job.yaml_set_comment_before_after_key(
        "split-fastq", before="FASTQ splitting (keep at 1, memory-intensive)"
    )
    tasks_per_job.yaml_set_comment_before_after_key(
        "demultiplex", before="Barcode demultiplexing tasks per job"
    )
    tasks_per_job.yaml_set_comment_before_after_key(
        "rna-map", before="RNA mapping tasks per job (higher = more efficient)"
    )

    config["tasks_per_job"] = tasks_per_job
    config.yaml_set_comment_before_after_key(
        "tasks_per_job",
        before=(
            "\nNumber of tasks to batch per SLURM job.\n"
            "Higher values reduce job overhead but increase per-job runtime."
        ),
    )

    # slurm_options section
    slurm_options = CommentedMap()

    slurm_options["extra-header-cmds"] = ""
    slurm_options.yaml_set_comment_before_after_key(
        "extra-header-cmds",
        before=(
            "Additional SBATCH directives to include in all jobs.\n"
            "Example: '#SBATCH --partition=gpu\\n#SBATCH --gres=gpu:1'"
        ),
    )

    # default slurm options
    default_opts = CommentedMap()
    default_opts["time"] = DQS("6:00:00")
    default_opts["cpus-per-task"] = 1
    default_opts["mem-per-cpu"] = "2GB"
    default_opts.yaml_set_comment_before_after_key(
        "time", before="Wall time limit in HH:MM:SS format"
    )
    default_opts.yaml_set_comment_before_after_key(
        "cpus-per-task", before="Number of CPU cores per task"
    )
    default_opts.yaml_set_comment_before_after_key(
        "mem-per-cpu", before="Memory allocation per CPU core"
    )
    slurm_options["default"] = default_opts
    slurm_options.yaml_set_comment_before_after_key(
        "default", before="\nDefault SLURM options used if job type not specified"
    )

    # split-fastq options
    split_fastq_opts = CommentedMap()
    split_fastq_opts["time"] = DQS("04:00:00")
    split_fastq_opts["cpus-per-task"] = 8
    split_fastq_opts["mem-per-cpu"] = "64GB"
    slurm_options["split-fastq"] = split_fastq_opts
    slurm_options.yaml_set_comment_before_after_key(
        "split-fastq", before="\nFASTQ splitting (high memory for large files)"
    )

    # trim-galore options
    trim_galore_opts = CommentedMap()
    trim_galore_opts["time"] = DQS("04:00:00")
    trim_galore_opts["cpus-per-task"] = 1
    trim_galore_opts["mem-per-cpu"] = "2GB"
    slurm_options["trim-galore"] = trim_galore_opts
    slurm_options.yaml_set_comment_before_after_key(
        "trim-galore", before="\nTrim Galore adapter trimming"
    )

    # demultiplex options
    demultiplex_opts = CommentedMap()
    demultiplex_opts["time"] = DQS("04:00:00")
    demultiplex_opts["cpus-per-task"] = 1
    demultiplex_opts["mem-per-cpu"] = "2GB"
    slurm_options["demultiplex"] = demultiplex_opts
    slurm_options.yaml_set_comment_before_after_key(
        "demultiplex", before="\nBarcode demultiplexing"
    )

    # int-demultiplex options
    int_demultiplex_opts = CommentedMap()
    int_demultiplex_opts["time"] = DQS("12:00:00")
    int_demultiplex_opts["cpus-per-task"] = 1
    int_demultiplex_opts["mem-per-cpu"] = "2GB"
    slurm_options["int-demultiplex"] = int_demultiplex_opts
    slurm_options.yaml_set_comment_before_after_key(
        "int-demultiplex", before="\nInternal barcode demultiplexing (longer runtime)"
    )

    # join-fastq-files options
    join_fastq_opts = CommentedMap()
    join_fastq_opts["time"] = DQS("12:00:00")
    join_fastq_opts["cpus-per-task"] = 1
    join_fastq_opts["mem-per-cpu"] = "4GB"
    slurm_options["join-fastq-files"] = join_fastq_opts
    slurm_options.yaml_set_comment_before_after_key(
        "join-fastq-files", before="\nJoin FASTQ files across chunks"
    )

    # rna-map options
    rna_map_opts = CommentedMap()
    rna_map_opts["time"] = DQS("6:00:00")
    rna_map_opts["cpus-per-task"] = 1
    rna_map_opts["mem-per-cpu"] = "2GB"
    slurm_options["rna-map"] = rna_map_opts
    slurm_options.yaml_set_comment_before_after_key(
        "rna-map", before="\nRNA mapping analysis"
    )

    # rna-map-combine options
    rna_map_combine_opts = CommentedMap()
    rna_map_combine_opts["time"] = DQS("6:00:00")
    rna_map_combine_opts["cpus-per-task"] = 1
    rna_map_combine_opts["mem-per-cpu"] = "2GB"
    slurm_options["rna-map-combine"] = rna_map_combine_opts
    slurm_options.yaml_set_comment_before_after_key(
        "rna-map-combine", before="\nCombine RNA mapping results from chunks"
    )

    # int-demultiplex-rna-map options
    int_dm_rna_map_opts = CommentedMap()
    int_dm_rna_map_opts["time"] = DQS("2:00:00")
    int_dm_rna_map_opts["cpus-per-task"] = 1
    int_dm_rna_map_opts["mem-per-cpu"] = "2GB"
    slurm_options["int-demultiplex-rna-map"] = int_dm_rna_map_opts
    slurm_options.yaml_set_comment_before_after_key(
        "int-demultiplex-rna-map", before="\nRNA mapping for internally demultiplexed reads"
    )

    # int-demultiplex-rna-map-combine options
    int_dm_rna_map_combine_opts = CommentedMap()
    int_dm_rna_map_combine_opts["time"] = DQS("2:00:00")
    int_dm_rna_map_combine_opts["cpus-per-task"] = 1
    int_dm_rna_map_combine_opts["mem-per-cpu"] = "2GB"
    slurm_options["int-demultiplex-rna-map-combine"] = int_dm_rna_map_combine_opts
    slurm_options.yaml_set_comment_before_after_key(
        "int-demultiplex-rna-map-combine",
        before="\nCombine internal demultiplex RNA mapping results",
    )

    config["slurm_options"] = slurm_options
    config.yaml_set_comment_before_after_key(
        "slurm_options",
        before=(
            "\nSLURM job resource configuration.\n"
            "Each job type can have its own time, CPU, and memory settings.\n"
            "If a job type is not specified here, 'default' values are used."
        ),
    )

    return config


@click.command(name="generate-example-config")
@click.option(
    "--output",
    "-o",
    default="example_config.yml",
    type=click.Path(),
    help="Output path for the example configuration file.",
)
def generate_example_config(output: str) -> None:
    """Generate an example YAML configuration file with comments.

    Creates a fully documented configuration file showing all available
    options and their default values. Use this as a starting point for
    customizing your workflow.
    """
    output_path = Path(output)

    config = _build_commented_config()

    yaml = YAML()
    yaml.default_flow_style = False
    yaml.indent(mapping=2, sequence=4, offset=2)

    with open(output_path, "w") as f:
        yaml.dump(config, f)

    click.echo(f"Generated example configuration file: {output_path}")
    click.echo("Edit this file and use with: rna-map-slurm setup --params <config.yml>")
