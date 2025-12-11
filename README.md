# rna_map_slurm

[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

A tool that takes care of all processes required to run rna_map on a SLURM cluster.

## Installation

### Option 1: Conda Environment (Recommended)

This will install all required tools (seqkit, sabre, trim-galore) and Python packages:

```shell
# Clone the repository
git clone https://github.com/jyesselm/rna_map_slurm.git
cd rna_map_slurm

# Create and activate the conda environment
conda env create -f environment.yml
conda activate rna_map_slurm
```

### Option 2: Pip Only

If you already have the bioinformatics tools installed:

```shell
# Basic installation
pip install git+https://github.com/jyesselm/rna_map_slurm

# With lab-specific dependencies
pip install "rna_map_slurm[full] @ git+https://github.com/jyesselm/rna_map_slurm"

# For development
pip install -e ".[dev]"
```

## Configuration

The workflow can be configured via a YAML file. Key options include:

```yaml
# Number of chunks to split FASTQ files into
fastq_chunks: 100

# Paths configuration
paths:
  tmp: "/scratch"  # Temporary directory (cluster-specific)

# Construct handling - constructs with more sequences than this
# threshold will be processed in dedicated jobs
construct_options:
  large_construct_threshold: 100
```

## Usage

```shell
# Fetch data from Google Sheets
rna-map-slurm get-data-csv <run_name>

# Set up the workflow
rna-map-slurm setup data.csv /path/to/fastq/dir

# Run the workflow
rna-map-slurm run

# Generate summaries
rna-map-slurm generate-summaries
```

## Development

```shell
# Install dev dependencies
pip install -e ".[dev]"

# Run linting
ruff check src/ tests/

# Run type checking
mypy src/rna_map_slurm

# Run tests
pytest tests/ -v
```

## Requirements

- Python >= 3.10
- seqkit
- sabre
- trim-galore
- SLURM cluster access
