"""Path utilities for locating package resources."""

from __future__ import annotations

import os
from pathlib import Path


def get_lib_path() -> str:
    """Get the path to the library's src directory.

    Returns:
        Path to the src directory containing rna_map_slurm.
    """
    this_file = Path(__file__)
    # Go up from config/paths.py to src/rna_map_slurm/config -> src/rna_map_slurm -> src
    return str(this_file.parent.parent.parent)


def get_resources_path() -> str:
    """Get the path to the resources directory.

    Returns:
        Path to the resources directory.
    """
    return os.path.join(get_lib_path(), "rna_map_slurm", "resources")
